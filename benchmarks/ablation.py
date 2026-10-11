#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Measures what each technique contributes: dinara-align built with one switched off at a time (see
`dinara_align/ablation.mojo`) against the build with none, everything else alike.

    pixi run bench-ablation                 # three rounds, every configuration in turn each round
    pixi run bench-ablation --rounds 5

Each configuration aligns samples of the main benchmark's datasets with traceback (`ablation.mojo`, the
samples drawn as `bench.py` draws them, larger for the 500 kbp reads), and the batch switch scores the
short-read batch of `batch_bench.py`. The configurations alternate within each round, so a drift in the machine falls on
all of them, and the table gives each one's median against the baseline's. Every configuration's costs
must equal the baseline's, or the run fails: a switch may change the time and nothing else. On Linux
each run is pinned to core 2 by `taskset`.

On x86 the alignment builds keep their jumps clear of 32-byte boundaries, which Skylake-era cores
penalise (Intel's JCC erratum fix): otherwise where a rebuild happens to place a loop moves the time by
several percent, as much as the smaller switches change it. `mojo build` offers no such option, so the
builds take Mojo's own assembly through the system's `clang -mbranches-within-32B-boundaries` and link
it as `mojo build` does (see `mitigated_build`). A baseline built so from another source file ran within
0.7% of this one on every sample. The batch harness carries GPU kernels the assembly leaves out, so it
is built as usual; its one switch changes the time threefold.
"""

import argparse
import os
import platform
import re
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import batch_bench
import bench
import harness
from harness import CACHE, DATA, HERE, ROOT, mojo_cpu_args

SWITCHES = ["ABLATE_REGROUP", "ABLATE_PAIRED", "ABLATE_GATHER", "ABLATE_AGREEMENT", "ABLATE_NEIGHBOURS"]
BATCH_SWITCHES = ["ABLATE_CERTIFIED"]
DATASETS = {
    "ont-500k": 20_000_000,
    "ont-500k-genvar": 20_000_000,
    "ont-1k": 2_000_000,
    "ont-10k": 2_000_000,
    "sars-cov-2": 2_000_000,
    "l100000.n1000.e1": 2_000_000,
    "l100000.n1000.e5": 2_000_000,
    "l100000.n1000.e10": 2_000_000,
    "l100000.n1000.e20": 2_000_000,
}
"""Each dataset the ablation samples, and the sequence its sample holds: the long reads take ten times the
default, which holds four of them, so their time is not one read's."""
BUILD = CACHE / "ablation"


def build(source: Path, switch: str, mitigated: bool = False) -> Path:
    """`source` built with `switch` defined, or with none for `NONE`; with `mitigated`, its jumps kept
    clear of 32-byte boundaries where the machine and its tools allow (see `mitigated_build`)."""
    binary = BUILD / f"{source.stem}-{switch}-{harness.CPU}"
    defines = [] if switch == "NONE" else ["-D", f"{switch}=1"]
    if mitigated and platform.machine() in ("x86_64", "AMD64") and shutil.which("clang") and shutil.which("cc"):
        mitigated_build(source, binary, defines)
        return binary
    subprocess.run(["mojo", "build", *defines, "-I", str(ROOT), str(source), "-o", str(binary), *mojo_cpu_args()], check=True)
    return binary


def mitigated_build(source: Path, binary: Path, defines: list) -> None:
    """Mojo's assembly for `source`, assembled with `clang -mbranches-within-32B-boundaries` and linked as
    `mojo build` links: the object, the Mojo runtime's two libraries, `--gc-sections` and `-lm`."""
    assembly = binary.with_suffix(".s")
    subprocess.run(
        ["mojo", "build", *defines, "-I", str(ROOT), "--emit", "asm", str(source), "-o", str(assembly), *mojo_cpu_args()],
        check=True,
    )
    # Mojo's LLVM writes `.prefalign`, a preferred alignment older assemblers predate: the same power as `.p2align`.
    text = re.sub(r"^\s*\.prefalign\s+(\d+),.*$", r"\t.p2align\t\1, 0x90", assembly.read_text(), flags=re.M)
    assembly.write_text(text)
    objects = binary.with_suffix(".o")
    subprocess.run(["clang", "-c", "-mbranches-within-32B-boundaries", str(assembly), "-o", str(objects)], check=True)
    library = Path(os.environ["CONDA_PREFIX"]) / "lib"
    subprocess.run(
        [
            "cc", str(objects), str(library / "libKGENCompilerRTShared.so"), "-Xlinker", "-rpath", "-Xlinker",
            str(library), str(library / "libAsyncRTMojoBindings.so"), "-o", str(binary), "-Wl,--gc-sections", "-lm",
        ],
        check=True,
    )


def pinned(command: list) -> list:
    """`command` pinned to one core where `taskset` exists."""
    return ["taskset", "-c", "2", *command] if shutil.which("taskset") else command


def main() -> None:
    """Builds every configuration, runs the rounds, fails on any cost that differs, prints the table."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--rounds", type=int, default=3)
    arguments = parser.parse_args()
    datasets = dict(bench.download())
    published = bench.published_costs()
    samples = [bench.sample(name, datasets[name], bases, published)[0] for name, bases in DATASETS.items()]
    names = [path.stem for path in samples]
    batch_bench.generate()
    pairs = DATA / "batches" / "illumina.tsv"
    BUILD.mkdir(parents=True, exist_ok=True)
    aligners = {switch: build(HERE / "ablation.mojo", switch, mitigated=True) for switch in ["NONE", *SWITCHES]}
    batchers = {switch: build(HERE / "batches.mojo", switch) for switch in ["NONE", *BATCH_SWITCHES]}
    times = defaultdict(list)
    answers = defaultdict(set)
    files = [str(path) for path in samples]
    for _ in range(arguments.rounds):
        for switch, binary in aligners.items():
            output = subprocess.run(pinned([str(binary), *files]), check=True, capture_output=True, text=True).stdout
            for line in output.splitlines():
                name, _, milliseconds, costs = line.split()
                name = name.removesuffix(".seq")
                times[(switch, name)].append(float(milliseconds))
                answers[name].add(costs)
        for switch, binary in batchers.items():
            for workload in batch_bench.WORKLOADS:
                output = subprocess.run(
                    pinned([str(binary), str(pairs), workload, "1"]), check=True, capture_output=True, text=True
                ).stdout
                _, _, seconds, checksum = output.split()
                times[(switch, workload)].append(float(seconds) * 1000)
                answers[workload].add(checksum)
    disagreeing = [name for name, found in answers.items() if len(found) != 1]
    if disagreeing:
        sys.exit(f"the configurations' costs differ on {disagreeing}")
    print(f"Median time, ms, and each switch's change against it, over {arguments.rounds} rounds ({harness.build_note()}):")
    columns = SWITCHES + BATCH_SWITCHES
    print(f"| dataset | baseline | {' | '.join(switch.removeprefix('ABLATE_').lower() for switch in columns)} |")
    print("| :-- | --: |" + " --: |" * len(columns))
    for name in names + batch_bench.WORKLOADS:
        base = statistics.median(times[("NONE", name)])
        cells = []
        for switch in columns:
            found = times.get((switch, name))
            cells.append(f"{100 * (statistics.median(found) / base - 1):+.1f}%" if found else "")
        print(f"| {name} | {base:.4g} | {' | '.join(cells)} |")


if __name__ == "__main__":
    main()
