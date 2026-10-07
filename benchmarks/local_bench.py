#!/usr/bin/env python3
"""Times dinara-align's local and overlap alignment against SSW, parasail, abPOA and hyalite.

    pixi run bench-local     # builds the rivals the first time, a few minutes; then about a minute

Every tool aligns every pair with its CIGAR at the same scores: a match 2, a mismatch -4, a gap of
`k` letters `6 + 2k` (dinara-align's `Mode.local(2)` and `Mode.overlap(2)` under
`Costs.affine(4, 6, 2)`). Each rival is cloned at a pinned commit (see `run.RIVALS`) and built in
`benchmarks/.cache/`. A tool's time is the faster of two passes over a workload, its mean per pair, on
one thread; its answer, the sum and position-weighted sum of its scores, must equal every other's on
the same workload, or the run fails. abPOA aligns to a graph, so its time includes adding the
reference to one, which any pairwise use of it pays. SSW and abPOA have no overlap mode.
"""

import platform
import random
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from run import CACHE, DATA, HERE, RESULTS, ROOT, cargo_runner, fetch, mutate

LOCAL = DATA / "local"
BUILD = CACHE / "local"
WORKLOADS = ["local-short", "local-window", "local-long", "overlap-reads"]
"""Short noisy pairs with planted cores, a 1 kbp read in a 10 kbp window, 10 kbp against 12 kbp, and
two 2 kbp reads overlapping by 0.5 to 1.5 kbp."""


def generate() -> None:
    """Writes the workload files, `name, reference, query` a row, from a fixed seed."""
    LOCAL.mkdir(parents=True, exist_ok=True)
    rng = random.Random(9)

    def sequence(length: int) -> str:
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


def build_rivals() -> Path:
    """parasail and abPOA's libraries, and the C driver over them and SSW's source."""
    parasail = fetch("parasail")
    abpoa = fetch("abPOA")
    ssw = fetch("SSW")
    build = parasail / "build"
    if not (build / "libparasail.a").exists():
        build.mkdir(exist_ok=True)
        subprocess.run(
            [
                "cmake",
                "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
                "-DCMAKE_BUILD_TYPE=Release",
                "-DBUILD_SHARED_LIBS=OFF",
                "..",
            ],
            cwd=build,
            check=True,
            capture_output=True,
        )
        subprocess.run(["make", "-j4", "parasail"], cwd=build, check=True, capture_output=True)
    if not (abpoa / "lib" / "libabpoa.a").exists():
        subprocess.run(["git", "submodule", "update", "--init", "--depth", "1"], cwd=abpoa, check=True)
        subprocess.run(["make", "libabpoa"], cwd=abpoa, check=True, capture_output=True)
    BUILD.mkdir(parents=True, exist_ok=True)
    binary = BUILD / "rivals"
    native = "-mcpu=native" if platform.machine().lower() in ("arm64", "aarch64") else "-march=native"
    subprocess.run(
        [
            "cc",
            "-O3",
            native,
            str(HERE / "local" / "rivals.c"),
            str(ssw / "src" / "ssw.c"),
            f"-I{ssw / 'src'}",
            f"-I{parasail}",
            f"-I{build}",
            f"-I{abpoa / 'include'}",
            str(build / "libparasail.a"),
            str(abpoa / "lib" / "libabpoa.a"),
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
    binary = BUILD / "dinara-align-local"
    BUILD.mkdir(parents=True, exist_ok=True)
    subprocess.run(["mojo", "build", "-I", str(ROOT), str(HERE / "local.mojo"), "-o", str(binary)], check=True)
    return binary


def main() -> None:
    generate()
    ours = build_ours()
    rivals = build_rivals()
    fetch("hyalite")
    hyalite = cargo_runner("hyalite", dict(__import__("os").environ)) if shutil.which("cargo") else None
    rows = []
    for workload in WORKLOADS:
        path = LOCAL / f"{workload}.tsv"
        commands = [[str(ours), str(path)], [str(rivals), str(path)]]
        if hyalite:
            commands.append([str(hyalite), "local", str(path)])
        for command in commands:
            print(f"{workload}: {Path(command[0]).name} ...", file=sys.stderr, flush=True)
            output = subprocess.run(command, check=True, capture_output=True, text=True).stdout
            rows.extend(line.split("\t") for line in output.strip().splitlines())
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
    table = "\n".join(lines)
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "local-results.md").write_text(table + "\n")
    print(table)


if __name__ == "__main__":
    main()
