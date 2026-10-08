#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Times dinara-align's local, overlap and scored infix alignment against SSW, parasail, abPOA and hyalite.

    pixi run bench-local     # builds and times the rivals the first time, a few minutes; then seconds
    pixi run bench-local --remeasure   # times the rivals again rather than replaying their kept rows

Every tool aligns every pair with its CIGAR at the same scores: a match 2, a mismatch -4, a gap of
`k` letters `6 + 2k` (dinara-align's `Mode.local(2)` and `Mode.overlap(2)` under
`Costs.affine(4, 6, 2)`). Each rival is cloned at a pinned commit (see `run.RIVALS`) and built in
`benchmarks/.cache/`. A tool's time is the faster of two passes over a workload, its mean per pair, on
one thread; its answer, the sum and position-weighted sum of its scores, must equal every other's on
the same workload, or the run fails. abPOA aligns to a graph, so its time includes adding the
reference to one, which any pairwise use of it pays. SSW and abPOA have no overlap or infix mode; the
infix rows score with a reward, `Mode.INFIX.with_match_score(2)`, as parasail's `sg_dx` and hyalite's
HW do.
"""

import platform
import random
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import run
from run import CACHE, DATA, HERE, RESULTS, ROOT, build_note, c_cpu_flag, cargo_runner, fetch, kept_rows, mojo_cpu_args, mutate

LOCAL = DATA / "local"
BUILD = CACHE / "local"
WORKLOADS = ["local-short", "local-window", "local-long", "overlap-reads", "infix-reads"]
"""Short noisy pairs with planted cores, a 1 kbp read in a 10 kbp window, 10 kbp against 12 kbp, two
2 kbp reads overlapping by 0.5 to 1.5 kbp, and a 1 kbp read placed whole in a 3 kbp window."""


def generate() -> None:
    """Writes the workload files, `name, reference, query` a row, from a fixed seed."""
    LOCAL.mkdir(parents=True, exist_ok=True)
    rng = random.Random(9)

    def sequence(length: int) -> str:
        """`length` random bases."""
        return "".join(rng.choice("ACGT") for _ in range(length))

    with open(LOCAL / "local-short.tsv", "w") as out:
        for trial in range(600):
            core = sequence(rng.randint(20, 300))
            reference = sequence(rng.randint(0, 200)) + core + sequence(rng.randint(0, 200))
            query = sequence(rng.randint(0, 50)) + mutate(core, [0.02, 0.08, 0.2, 0.35][trial % 4], "ACGT", rng)
            query += sequence(rng.randint(0, 50))
            if trial % 10 == 0:
                query = sequence(rng.randint(5, 100))
            out.write(f"local-short\t{reference}\t{query}\n")
    with open(LOCAL / "local-window.tsv", "w") as out:
        for _ in range(20):
            reference = sequence(10_000)
            start = rng.randint(0, 9_000)
            out.write(f"local-window\t{reference}\t{mutate(reference[start : start + 1000], 0.1, 'ACGT', rng)}\n")
    with open(LOCAL / "local-long.tsv", "w") as out:
        for _ in range(5):
            core = sequence(10_000)
            out.write(f"local-long\t{sequence(1000) + core + sequence(1000)}\t{mutate(core, 0.05, 'ACGT', rng)}\n")
    with open(LOCAL / "overlap-reads.tsv", "w") as out:
        for _ in range(100):
            genome = sequence(3_500)
            shift = rng.randint(500, 1_500)
            first = mutate(genome[shift : shift + 2_000], 0.05, "ACGT", rng)
            second = mutate(genome[shift + 2_000 - rng.randint(500, 1_500) :][:2_000], 0.05, "ACGT", rng)
            out.write(f"overlap-reads\t{first}\t{second}\n")
    with open(LOCAL / "infix-reads.tsv", "w") as out:
        for _ in range(100):
            reference = sequence(3_000)
            start = rng.randint(0, 2_000)
            out.write(f"infix-reads\t{reference}\t{mutate(reference[start : start + 1000], 0.1, 'ACGT', rng)}\n")


def build_rivals() -> Path:
    """parasail and abPOA's libraries, and the C driver over them and SSW's source."""
    parasail = fetch("parasail")
    abpoa = fetch("abPOA")
    ssw = fetch("SSW")
    # Each CPU its own build: parasail dispatches its kernels at run time either way, but its own C
    # follows the flag too.
    build = parasail / f"build-{run.CPU}"
    if not (build / "libparasail.a").exists():
        build.mkdir(exist_ok=True)
        subprocess.run(
            [
                "cmake",
                "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
                "-DCMAKE_BUILD_TYPE=Release",
                "-DBUILD_SHARED_LIBS=OFF",
                f"-DCMAKE_C_FLAGS={c_cpu_flag()}",
                "..",
            ],
            cwd=build,
            check=True,
            capture_output=True,
        )
        subprocess.run(["make", "-j4", "parasail"], cwd=build, check=True, capture_output=True)
    # abPOA's Makefile builds for the host's own CPU on x86 and for the M1 on Apple's ARM by default, its
    # `SIMD_FLAG`; given the CPU every tool is built for instead, kept per CPU.
    abpoa_library = abpoa / "lib" / f"libabpoa-{run.CPU}.a"
    if not abpoa_library.exists():
        subprocess.run(["git", "submodule", "update", "--init", "--depth", "1"], cwd=abpoa, check=True)
        subprocess.run(["make", "clean"], cwd=abpoa, check=True, capture_output=True)
        simd = c_cpu_flag() + (" -D__AVX2__" if run.arm() else "")
        subprocess.run(["make", "libabpoa", f"SIMD_FLAG={simd}"], cwd=abpoa, check=True, capture_output=True)
        shutil.copy(abpoa / "lib" / "libabpoa.a", abpoa_library)
    BUILD.mkdir(parents=True, exist_ok=True)
    binary = BUILD / f"rivals-{run.CPU}"
    subprocess.run(
        [
            "cc",
            "-O3",
            c_cpu_flag(),
            str(HERE / "local" / "rivals.c"),
            str(ssw / "src" / "ssw.c"),
            f"-I{ssw / 'src'}",
            f"-I{parasail}",
            f"-I{build}",
            f"-I{abpoa / 'include'}",
            str(build / "libparasail.a"),
            str(abpoa_library),
            "-lz",
            "-lm",
            "-lpthread",
            "-o",
            str(binary),
        ],
        check=True,
    )
    return binary


def build_ours() -> Path:
    """dinara-align's runner, `local.mojo`, built for the chosen CPU."""
    binary = BUILD / f"dinara-align-local-{run.CPU}"
    BUILD.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["mojo", "build", "-I", str(ROOT), str(HERE / "local.mojo"), "-o", str(binary), *mojo_cpu_args()], check=True
    )
    return binary


def main() -> None:
    """Generates the workloads, runs every tool on each, fails unless their answers agree, and writes the
    table of times to `results/local-results.md`."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cpu", default="native", help="the CPU every tool is built for (default: the host's)")
    parser.add_argument("--remeasure", action="store_true", help="run the rivals again, not their kept rows")
    options = parser.parse_args()
    run.set_cpu(options.cpu)
    run.REMEASURE = options.remeasure
    generate()
    ours = build_ours()
    # The rivals are built only when some workload has no kept rows for them (see `run.kept_rows`).
    built = {}

    def rows_of(command: list) -> list:
        """One runner's rows on one workload."""
        print(f"{Path(command[-1]).stem}: {Path(command[0]).name} ...", file=sys.stderr, flush=True)
        output = subprocess.run([str(part) for part in command], check=True, capture_output=True, text=True).stdout
        return [line.split("\t") for line in output.strip().splitlines()]

    def rivals() -> Path:
        """The C driver over SSW, parasail and abPOA, built once."""
        if "rivals" not in built:
            built["rivals"] = build_rivals()
        return built["rivals"]

    def hyalite() -> Path:
        """hyalite's runner, built once."""
        if "hyalite" not in built:
            fetch("hyalite")
            built["hyalite"] = cargo_runner("hyalite", dict(__import__("os").environ))
        return built["hyalite"]

    rows = []
    for workload in WORKLOADS:
        path = LOCAL / f"{workload}.tsv"
        rows.extend(rows_of([ours, path]))
        rows.extend(
            kept_rows(f"SSW, parasail and abPOA on {workload}", ["SSW", "parasail", "abPOA"],
                      [HERE / "local" / "rivals.c", path], lambda: rows_of([rivals(), path]))
        )
        if shutil.which("cargo"):
            rows.extend(
                kept_rows(f"hyalite on {workload}", ["hyalite"], [HERE / "hyalite", path],
                          lambda: rows_of([hyalite(), "local", path]))
            )
    answers = defaultdict(set)
    times = defaultdict(dict)
    for tool, workload, task, seconds, answer in rows:
        answers[(workload, task)].add(answer)
        times[(workload, task)][tool] = float(seconds)
    disagreeing = [key for key, values in answers.items() if len(values) > 1]
    if disagreeing:
        sys.exit(f"the tools' scores disagree on {disagreeing}")
    tools = ["dinara-align", "SSW", "parasail", "abPOA", "hyalite"]
    lines = ["| workload | task | " + " | ".join(tools) + " |", "| :-- | :-- | " + " | ".join("--:" for _ in tools) + " |"]
    for (workload, task), measured in times.items():
        cells = []
        for tool in tools:
            seconds = measured.get(tool)
            cells.append("—" if seconds is None else (f"{seconds * 1e6:.0f} µs" if seconds < 1e-3 else f"{seconds * 1e3:.2f} ms"))
        lines.append(f"| {workload} | {task} | " + " | ".join(cells) + " |")
    table = build_note() + "\n".join(lines)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "local-results.md").write_text(table + "\n")
    print(table)


if __name__ == "__main__":
    main()
