#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Times batches of short-read global scores on every thread, as Accelign's short-read case study does
(Kallenborn et al., BMC Bioinformatics 2026), against the CPU aligners it compares.

    pixi run bench-batch     # builds the rivals the first time; then a minute or two

Accelign's study scores 10 million Illumina reads of 148 bp against the reference sections BWA placed
them in, 103 to 209 bp long, each pair's global alignment score alone, with a match 0, a mismatch -1 and
a gap of `k` letters `-(3 + (k - 1))`: dinara-align's `Costs.affine(1, 2, 1)`. Those reads are not ours
to ship, so the pairs here are drawn to the same shape from a fixed seed: 500,000 reads of 148 bp at 1%
error against sections of their source, most of them the read's own 148 bp, a fifth with up to 22 bases
fewer or 30 more at either end.

| workload | dinara-align | rivals |
| :-- | :-- | :-- |
| illumina-affine | `distances(..., Costs.affine(1, 2, 1))` | WFA2-lib, exact, score only; KSW2's `extz2`, score only, no band; parasail's striped `nw` |
| illumina-edit | `distances(...)` at unit costs | Edlib's NW, distance only, its own band |

Every tool runs on every thread, the rivals through OpenMP, one aligner a thread, and its time is the
faster of two passes over the whole batch. The costs, summed and position-weighted, must agree between
every tool on a workload, or the run fails. Each rival is pinned (see `run.RIVALS`), built for the same
CPU as everything else, and timed once (see `run.kept_rows`).
"""

import random
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import mode_bench
import run
from run import CACHE, DATA, HERE, RESULTS, ROOT, build_note, kept_rows, mojo_cpu_args, mutate

BATCHES = DATA / "batches"
BUILD = CACHE / "batches"
WORKLOADS = ["illumina-affine", "illumina-edit"]
PAIRS = 500_000


def generate() -> None:
    """The short-read pairs, `name, reference, query` a row, written once for both workloads."""
    BATCHES.mkdir(parents=True, exist_ok=True)
    path = BATCHES / "illumina.tsv"
    if path.exists():
        return
    rng = random.Random(13)
    rows = []
    for _ in range(PAIRS):
        # A source with room either side, the read from its middle, the section around where it lies.
        source = "".join(rng.choices("ACGT", k=260))
        read = mutate(source[56 : 56 + 148], 0.01, "ACGT", rng)
        before, after = (rng.randint(-22, 30), rng.randint(-22, 30)) if rng.random() < 0.2 else (0, 0)
        rows.append(f"illumina\t{source[56 - before : 56 + 148 + after]}\t{read}\n")
    path.write_text("".join(rows))


def build_ours() -> Path:
    """dinara-align's runner, `batches.mojo`, built for the chosen CPU."""
    binary = BUILD / f"dinara-align-batches-{run.CPU}"
    BUILD.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["mojo", "build", "-I", str(ROOT), str(HERE / "batches.mojo"), "-o", str(binary), *mojo_cpu_args()], check=True
    )
    return binary


def main() -> None:
    """Generates the pairs, runs every tool on both workloads, fails unless their answers agree, and writes
    the table of times to `results/batch-results.md`."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cpu", default="native", help="the CPU every tool is built for (default: the host's)")
    parser.add_argument("--remeasure", action="store_true", help="run the rivals again, not their kept rows")
    arguments = parser.parse_args()
    run.set_cpu(arguments.cpu)
    run.REMEASURE = arguments.remeasure
    generate()
    path = BATCHES / "illumina.tsv"
    ours = build_ours()
    built = {}

    def rows_of(command: list) -> list:
        """One runner's rows on the pairs."""
        print(f"{command[-1]}: {Path(command[0]).name} ...", file=sys.stderr, flush=True)
        output = subprocess.run([str(part) for part in command], check=True, capture_output=True, text=True).stdout
        return [line.split("\t") for line in output.strip().splitlines()]

    def rivals() -> Path:
        """The OpenMP driver over WFA2-lib, KSW2, parasail and Edlib, built once."""
        if "rivals" not in built:
            built["rivals"] = mode_bench.build_rivals(HERE / "batches" / "rivals.cpp", "batch-rivals", ["-fopenmp"])
        return built["rivals"]

    rows = []
    for workload in WORKLOADS:
        rows.extend(rows_of([ours, path, workload]))
        rows.extend(
            kept_rows(
                f"WFA2-lib, KSW2, parasail and Edlib on {workload}",
                ["WFA2-lib", "ksw2", "parasail", "edlib"],
                [HERE / "batches" / "rivals.cpp", path],
                lambda: rows_of([rivals(), path, workload]),
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
    tools = ["dinara-align", "WFA2-lib", "KSW2", "parasail", "Edlib"]
    lines = ["| workload | " + " | ".join(tools) + " |", "| :-- | " + " | ".join("--:" for _ in tools) + " |"]
    for workload, measured in times.items():
        fastest = min(measured.values())
        cells = []
        for tool in tools:
            seconds = measured.get(tool)
            if seconds is None:
                cells.append("—")
                continue
            text = f"{seconds * 1e3:.0f} ms" if seconds < 10 else f"{seconds:.1f} s"
            cells.append(f"**{text}**" if seconds == fastest else text)
        lines.append(f"| {workload} | " + " | ".join(cells) + " |")
    table = build_note() + "\n".join(lines)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "batch-results.md").write_text(table + "\n")
    print(table)


if __name__ == "__main__":
    main()
