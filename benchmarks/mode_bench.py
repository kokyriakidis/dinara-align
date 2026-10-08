#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Times dinara-align's other modes against the aligners that offer the same one: free ends at unit and
gap-affine costs, seed extension, two-piece gaps, and substitution tables with more than one mismatch
score.

    pixi run bench-modes     # builds and times the rivals the first time, a few minutes; then seconds
    pixi run bench-modes --remeasure   # times the rivals again rather than replaying their kept rows

| workload | dinara-align | rivals |
| :-- | :-- | :-- |
| infix-edit, infix-edit-long | `Costs.edit()`, `Mode.INFIX` | Edlib's HW, WFA2-lib's ends-free edit |
| prefix-edit | `Costs.edit()`, `Mode.PREFIX` | Edlib's SHW, WFA2-lib's ends-free edit |
| infix-affine | `Costs.affine(4, 6, 2)`, `Mode.INFIX` | WFA2-lib's ends-free gap-affine |
| extension | `Costs.affine(4, 6, 2)`, `Mode.extension(2)` | KSW2's `extz2`, extension only, no Z-drop |
| extension-bonus | the same with `end_bonus=50` | KSW2's `extz2` with an end bonus of 50 |
| two-piece | `Costs.two_piece(4, 6, 2, 24, 1)`, global | KSW2's `extd2`, WFA2-lib's two-piece gap-affine |
| table-global, table-local | a `Scoring` with transitions apart from transversions | parasail, and SSW locally |

Every tool aligns every pair with its CIGAR, on one thread, and its time is the faster of two passes over
a workload, its mean per pair. Its answer, the sum and position-weighted sum of its costs or scores, must
equal every other tool's on the same workload, or the run fails. Each rival is cloned at a pinned commit
(see `run.RIVALS`) and built in `benchmarks/.cache/` for the same CPU as everything else. KSW2's kernels are
SSE only, so it runs on x86 alone. Z-drop is left out: KSW2 gauges it by anti-diagonal and dinara-align
as WFA2-lib does, a cost at a time, so the two stop at different places and their times compare different
work. Deletions priced apart from insertions are left out too: no rival here offers them.
"""

import random
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import run
from run import CACHE, DATA, HERE, RESULTS, ROOT, build_note, c_cpu_flag, fetch, kept_rows, mojo_cpu_args, mutate

MODES = DATA / "modes"
BUILD = CACHE / "modes"
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
]


def generate() -> None:
    """Writes each workload's pairs, `name, reference, query` a row, from a fixed seed."""
    MODES.mkdir(parents=True, exist_ok=True)
    rng = random.Random(11)

    def sequence(length: int, alphabet: str = "ACGT") -> str:
        """`length` random letters of `alphabet`."""
        return "".join(rng.choice(alphabet) for _ in range(length))

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
    prefix = []
    for _ in range(100):
        reference = sequence(2_000)
        prefix.append((reference, mutate(reference[:1_000], 0.1, "ACGT", rng)))
    write("prefix-edit", prefix)
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

    # Reads that follow the reference at 5% and end in a short stretch of noise, which an end bonus may
    # or may not align through: last, so every earlier workload's pairs stay as they were.
    bonus = []
    for _ in range(200):
        reference = sequence(2_000)
        follows = rng.randint(300, 1_400)
        bonus.append((reference, mutate(reference[:follows], 0.05, "ACGT", rng) + sequence(rng.randint(0, 80))))
    write("extension-bonus", bonus)


def build_rivals(driver: Path = HERE / "modes" / "rivals.cpp", name: str = "rivals", extra: list = []) -> Path:
    """A C++ driver, `driver`, over Edlib, WFA2-lib, parasail and SSW, and on x86 KSW2, built as `name` for the
    chosen CPU with `extra` flags besides: by default this benchmark's own."""
    edlib = fetch("edlib")
    wfa = fetch("WFA2-lib")
    parasail = fetch("parasail")
    ssw = fetch("SSW")
    # The same per-CPU parasail build `local_bench.py` makes, built here if that has not run.
    parasail_build = parasail / f"build-{run.CPU}"
    if not (parasail_build / "libparasail.a").exists():
        parasail_build.mkdir(exist_ok=True)
        subprocess.run(
            ["cmake", "-DCMAKE_POLICY_VERSION_MINIMUM=3.5", "-DCMAKE_BUILD_TYPE=Release", "-DBUILD_SHARED_LIBS=OFF",
             f"-DCMAKE_C_FLAGS={c_cpu_flag()}", ".."],
            cwd=parasail_build, check=True, capture_output=True,
        )
        subprocess.run(["make", "-j4", "parasail"], cwd=parasail_build, check=True, capture_output=True)
    # WFA2-lib's Makefile builds for the host's own CPU; given the chosen one instead, kept per CPU.
    wfa_library = wfa / "lib" / f"libwfa-{run.CPU}.a"
    if not wfa_library.exists():
        subprocess.run(["make", "clean"], cwd=wfa, check=True, capture_output=True)
        subprocess.run(["make", "setup", "lib_wfa", f"CC_FLAGS=-Wall -fPIE -O3 {c_cpu_flag()}"], cwd=wfa, check=True, capture_output=True)
        (wfa / "lib" / "libwfa.a").rename(wfa_library)
    BUILD.mkdir(parents=True, exist_ok=True)
    binary = BUILD / f"{name}-{run.CPU}"
    flags = ["-O3", c_cpu_flag(), *extra]
    sources = [str(ssw / "src" / "ssw.c")]
    defines = []
    if not run.arm():
        ksw2 = fetch("ksw2")
        sources += [str(ksw2 / "ksw2_extz2_sse.c"), str(ksw2 / "ksw2_extd2_sse.c"), str(ksw2 / "kalloc.c")]
        defines += ["-DWITH_KSW2", f"-I{ksw2}"]
    objects = []
    for source in sources:
        target = BUILD / (Path(source).stem + f"-{run.CPU}.o")
        subprocess.run(["cc", "-c", *flags, *defines, f"-I{ssw / 'src'}", source, "-o", str(target)], check=True)
        objects.append(str(target))
    subprocess.run(
        [
            "c++", "-std=c++17", *flags, *defines,
            str(driver), str(edlib / "edlib" / "src" / "edlib.cpp"), *objects,
            f"-I{edlib / 'edlib' / 'include'}", f"-I{wfa}", f"-I{ssw / 'src'}", f"-I{parasail}", f"-I{parasail_build}",
            str(wfa_library), str(parasail_build / "libparasail.a"),
            "-lz", "-lm", "-lpthread", "-o", str(binary),
        ],
        check=True,
    )
    return binary


def build_ours() -> Path:
    """dinara-align's runner, `modes.mojo`, built for the chosen CPU."""
    binary = BUILD / f"dinara-align-modes-{run.CPU}"
    BUILD.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["mojo", "build", "-I", str(ROOT), str(HERE / "modes.mojo"), "-o", str(binary), *mojo_cpu_args()], check=True
    )
    return binary


def main() -> None:
    """Generates the workloads, runs every tool on each, fails unless their answers agree, and writes the
    table of times to `results/mode-results.md`."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cpu", default="native", help="the CPU every tool is built for (default: the host's)")
    parser.add_argument("--only", nargs="*", help="run these workloads alone")
    parser.add_argument("--remeasure", action="store_true", help="run the rivals again, not their kept rows")
    arguments = parser.parse_args()
    run.set_cpu(arguments.cpu)
    run.REMEASURE = arguments.remeasure
    generate()
    ours = build_ours()
    # The rivals' driver is built only when some workload has no kept rows for it (see `run.kept_rows`).
    built = {}

    def rows_of(command: list) -> list:
        """One runner's rows on one workload."""
        print(f"{Path(command[-1]).stem}: {Path(command[0]).name} ...", file=sys.stderr, flush=True)
        output = subprocess.run([str(part) for part in command], check=True, capture_output=True, text=True).stdout
        return [line.split("\t") for line in output.strip().splitlines()]

    def rivals() -> Path:
        """The driver over Edlib, WFA2-lib, KSW2, parasail and SSW, built once."""
        if "rivals" not in built:
            built["rivals"] = build_rivals()
        return built["rivals"]

    rows = []
    for workload in arguments.only or WORKLOADS:
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
    tools = ["dinara-align", "Edlib", "WFA2-lib", "KSW2", "parasail", "SSW"]
    lines = ["| workload | " + " | ".join(tools) + " |", "| :-- | " + " | ".join("--:" for _ in tools) + " |"]
    for workload, measured in times.items():
        fastest = min(measured.values())
        cells = []
        for tool in tools:
            seconds = measured.get(tool)
            if seconds is None:
                cells.append("—")
                continue
            text = f"{seconds * 1e6:.0f} µs" if seconds < 1e-3 else f"{seconds * 1e3:.2f} ms"
            cells.append(f"**{text}**" if seconds == fastest else text)
        lines.append(f"| {workload} | " + " | ".join(cells) + " |")
    table = build_note() + "\n".join(lines)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "mode-results.md").write_text(table + "\n")
    print(table)


if __name__ == "__main__":
    main()
