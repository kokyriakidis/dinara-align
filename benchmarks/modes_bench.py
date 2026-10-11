#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Times dinara-align's modes beyond global alignment against the aligners that offer the same one: free
ends at unit and gap-affine costs, seed extension, two-piece gaps, substitution tables, and local, overlap
and scored infix alignment.

    pixi run bench-modes                # builds and times the rivals the first time, a few minutes; then seconds
    pixi run bench-modes --only local-long table-local
    pixi run bench-modes --remeasure    # times the rivals again rather than replaying their kept rows

| workload | dinara-align | rivals |
| :-- | :-- | :-- |
| infix-edit, infix-edit-long | `Costs.edit()`, `Mode.INFIX` | Edlib's HW, WFA2-lib's ends-free edit |
| prefix-edit | `Costs.edit()`, `Mode.PREFIX` | Edlib's SHW, WFA2-lib's ends-free edit |
| infix-affine | `Costs.affine(4, 6, 2)`, `Mode.INFIX` | WFA2-lib's ends-free gap-affine |
| extension | `Costs.affine(4, 6, 2)`, `Mode.extension(2)` | KSW2's `extz2`, extension only, no Z-drop |
| extension-bonus | the same with `end_bonus=50` | KSW2's `extz2` with an end bonus of 50 |
| two-piece | `Costs.two_piece(4, 6, 2, 24, 1)`, global | KSW2's `extd2`, WFA2-lib's two-piece gap-affine |
| table-global, table-local | a `Scoring` with transitions apart from transversions | parasail, and SSW locally |
| local-short, local-window, local-long | `Costs.affine(4, 6, 2)`, `Mode.local(2)` | SSW, parasail's `sw` |
| overlap-reads | `Costs.affine(4, 6, 2)`, `Mode.overlap(2)` | parasail's `sg` |
| infix-reads | `Costs.affine(4, 6, 2)`, `Mode.INFIX.with_match_score(2)` | parasail's `sg_dx` |

Every tool aligns every pair with its CIGAR, on one thread, and its time is the faster of two passes over
a workload, its mean per pair. Its answer, the sum and position-weighted sum of its costs or scores, must
equal every other tool's on the same workload, or the run fails. KSW2's kernels are SSE only, so it runs on
x86 alone. Z-drop is left out: KSW2 gauges it by anti-diagonal and dinara-align as WFA2-lib does, a cost at
a time, so the two stop at different places and their times compare different work.
"""

import argparse
import random
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import harness
from harness import DATA, HERE, RESULTS, build_note, c_driver, kept_rows, mojo_binary, mutate, set_cpu

MODES = DATA / "modes"
WORKLOADS = [
    "infix-edit",
    "infix-edit-long",
    "prefix-edit",
    "infix-affine",
    "extension",
    "extension-bonus",
    "two-piece",
    "table-global",
    "table-local",
    "local-short",
    "local-window",
    "local-long",
    "overlap-reads",
    "infix-reads",
]
TOOLS = ["dinara-align", "Edlib", "WFA2-lib", "KSW2", "parasail", "SSW"]


def generate() -> None:
    """Writes each workload's pairs, `name, reference, query` a row, from fixed seeds."""
    MODES.mkdir(parents=True, exist_ok=True)
    rng = random.Random(11)

    def sequence(length: int) -> str:
        """`length` random bases."""
        return "".join(rng.choice("ACGT") for _ in range(length))

    def write(name: str, pairs: list) -> None:
        """One workload's file."""
        with open(MODES / f"{name}.tsv", "w") as out:
            for reference, query in pairs:
                out.write(f"{name}\t{reference}\t{query}\n")

    def placed(count: int, read: int, window: int, rate: float) -> list:
        """`count` reads of `read` bases at `rate` edits a base, each in a `window` of its reference."""
        pairs = []
        for _ in range(count):
            reference = sequence(window)
            start = rng.randint(0, window - read)
            pairs.append((reference, mutate(reference[start : start + read], rate, "ACGT", rng)))
        return pairs

    infix = placed(100, 1_000, 3_000, 0.1)
    write("infix-edit", infix)
    write("infix-affine", infix)
    write("infix-edit-long", placed(20, 10_000, 30_000, 0.05))
    # A read's start against a reference running on past its end.
    write("prefix-edit", [(reference, mutate(reference[:1_000], 0.1, "ACGT", rng)) for reference in (sequence(2_000) for _ in range(100))])
    # A seed's extension: the read follows the reference for a while at 5%, then turns to noise.
    extension = []
    for _ in range(200):
        reference = sequence(2_000)
        follows = rng.randint(300, 1_400)
        extension.append((reference, mutate(reference[:follows], 0.05, "ACGT", rng) + sequence(1_500 - follows)))
    write("extension", extension)
    # Long gaps, where the second piece is the cheaper: 5 kbp at 5% with three indels of 100 to 400 bases.
    two_piece = []
    for _ in range(50):
        reference = sequence(5_000)
        query = mutate(reference, 0.05, "ACGT", rng)
        for _ in range(3):
            at = rng.randint(0, len(query) - 500)
            length = rng.randint(100, 400)
            query = query[:at] + query[at + length :] if rng.random() < 0.5 else query[:at] + sequence(length) + query[at:]
        two_piece.append((reference, query))
    write("two-piece", two_piece)
    write("table-global", [(reference, mutate(reference, 0.1, "ACGT", rng)) for reference in (sequence(1_000) for _ in range(100))])
    write("table-local", placed(20, 1_000, 10_000, 0.1))
    # Reads that follow the reference at 5% and end in a short stretch of noise, which an end bonus may or
    # may not align through.
    bonus = []
    for _ in range(200):
        reference = sequence(2_000)
        follows = rng.randint(300, 1_400)
        bonus.append((reference, mutate(reference[:follows], 0.05, "ACGT", rng) + sequence(rng.randint(0, 80))))
    write("extension-bonus", bonus)

    # The local workloads, from a seed of their own: short noisy pairs with planted cores, a 1 kbp read in a
    # 10 kbp window, 10 kbp against 12 kbp, two 2 kbp reads overlapping by 0.5 to 1.5 kbp, and a 1 kbp read
    # placed whole in a 3 kbp window.
    rng = random.Random(9)
    short = []
    for trial in range(600):
        core = sequence(rng.randint(20, 300))
        reference = sequence(rng.randint(0, 200)) + core + sequence(rng.randint(0, 200))
        query = sequence(rng.randint(0, 50)) + mutate(core, [0.02, 0.08, 0.2, 0.35][trial % 4], "ACGT", rng)
        query += sequence(rng.randint(0, 50))
        if trial % 10 == 0:
            query = sequence(rng.randint(5, 100))
        short.append((reference, query))
    write("local-short", short)
    window = []
    for _ in range(20):
        reference = sequence(10_000)
        start = rng.randint(0, 9_000)
        window.append((reference, mutate(reference[start : start + 1000], 0.1, "ACGT", rng)))
    write("local-window", window)
    long = []
    for _ in range(5):
        core = sequence(10_000)
        long.append((sequence(1000) + core + sequence(1000), mutate(core, 0.05, "ACGT", rng)))
    write("local-long", long)
    overlap = []
    for _ in range(100):
        genome = sequence(3_500)
        shift = rng.randint(500, 1_500)
        first = mutate(genome[shift : shift + 2_000], 0.05, "ACGT", rng)
        second = mutate(genome[shift + 2_000 - rng.randint(500, 1_500) :][:2_000], 0.05, "ACGT", rng)
        overlap.append((first, second))
    write("overlap-reads", overlap)
    placed_reads = []
    for _ in range(100):
        reference = sequence(3_000)
        start = rng.randint(0, 2_000)
        placed_reads.append((reference, mutate(reference[start : start + 1000], 0.1, "ACGT", rng)))
    write("infix-reads", placed_reads)


def main() -> None:
    """Generates the workloads, runs every tool on each, fails unless their answers agree, and writes the
    table of times to `results/modes.md`."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cpu", default="native", help="the CPU every tool is built for (default: the host's)")
    parser.add_argument("--only", nargs="*", help="run these workloads alone")
    parser.add_argument("--remeasure", action="store_true", help="run the rivals again, not their kept rows")
    options = parser.parse_args()
    set_cpu(options.cpu)
    harness.REMEASURE = options.remeasure
    generate()
    ours = mojo_binary(HERE / "modes.mojo", "dinara-align-modes")
    # The rivals' driver is built only when some workload has no kept rows for it (see `harness.kept_rows`).
    built = {}

    def rows_of(command: list) -> list:
        """One runner's rows on one workload."""
        print(f"{Path(command[-1]).stem}: {Path(command[0]).name} ...", file=sys.stderr, flush=True)
        output = subprocess.run([str(part) for part in command], check=True, capture_output=True, text=True).stdout
        return [line.split("\t") for line in output.strip().splitlines()]

    def rivals() -> Path:
        """The driver over Edlib, WFA2-lib, KSW2, parasail and SSW, built once."""
        if "rivals" not in built:
            built["rivals"] = c_driver(HERE / "modes" / "rivals.cpp", "mode-rivals")
        return built["rivals"]

    rows = []
    for workload in options.only or WORKLOADS:
        path = MODES / f"{workload}.tsv"
        rows.extend(rows_of([ours, path]))
        rows.extend(
            kept_rows(
                f"Edlib, WFA2-lib, KSW2, parasail and SSW on {workload}",
                ["edlib", "WFA2-lib", "ksw2", "parasail", "SSW"],
                [HERE / "modes" / "rivals.cpp", path],
                lambda: rows_of([rivals(), path]),
            )
        )
    answers = defaultdict(dict)
    times = defaultdict(dict)
    for tool, workload, seconds, answer in rows:
        answers[workload][tool] = answer
        times[workload][tool] = float(seconds)
    disagreeing = {workload: found for workload, found in answers.items() if len(set(found.values())) > 1}
    if disagreeing:
        sys.exit(f"the tools' answers disagree: {disagreeing}")
    lines = ["| workload | " + " | ".join(TOOLS) + " |", "| :-- | " + " | ".join("--:" for _ in TOOLS) + " |"]
    for workload, measured in times.items():
        fastest = min(measured.values())
        cells = []
        for tool in TOOLS:
            seconds = measured.get(tool)
            if seconds is None:
                cells.append("—")
                continue
            text = harness.duration(seconds)
            cells.append(f"**{text}**" if seconds == fastest else text)
        lines.append(f"| {workload} | " + " | ".join(cells) + " |")
    table = build_note() + "\n".join(lines)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "modes.md").write_text(table + "\n")
    print(table)


if __name__ == "__main__":
    main()
