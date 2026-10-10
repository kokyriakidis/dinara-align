#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Measures what each technique contributes: dinara-align built with one switched off at a time (see
`dinara_align/ablation.mojo`) against the build with none, everything else alike.

    pixi run bench-ablation                 # three rounds, every configuration in turn each round
    pixi run bench-ablation --rounds 5

Each configuration aligns the A*PA2 benchmark's samples with traceback (`ablation.mojo`; the samples
come from `pixi run bench-astarpa2`), and the batch switch scores the short-read batch of
`batch_bench.py`. The configurations alternate within each round, so a drift in the machine falls on
all of them, and the table gives each one's median against the baseline's. Every configuration's costs
must equal the baseline's, or the run fails: a switch may change the time and nothing else. On Linux
each run is pinned to core 2 by `taskset`.
"""

import argparse
import shutil
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import batch_bench
import run
from run import CACHE, DATA, HERE, ROOT, mojo_cpu_args

SWITCHES = ["ABLATE_REGROUP", "ABLATE_PAIRED", "ABLATE_GATHER", "ABLATE_AGREEMENT", "ABLATE_NEIGHBOURS"]
BATCH_SWITCHES = ["ABLATE_CERTIFIED"]
SAMPLES = DATA / "pa-bench" / "samples"
SETS = [
    "ont-500k-b20000000-p4",
    "ont-500k-genvar-b20000000-p4",
    "ont-1k-b2000000-p4",
    "ont-10k-b2000000-p4",
    "sars-cov-2-b2000000-p4",
    "Uniform-t1000000-n100000-e0.01-b2000000-p4",
    "Uniform-t1000000-n100000-e0.05-b2000000-p4",
    "Uniform-t1000000-n100000-e0.1-b2000000-p4",
    "Uniform-t1000000-n100000-e0.15-b2000000-p4",
]
BUILD = CACHE / "ablation"


def build(source: Path, switch: str) -> Path:
    """`source` built with `switch` defined, or with none for `NONE`."""
    binary = BUILD / f"{source.stem}-{switch}-{run.CPU}"
    defines = [] if switch == "NONE" else ["-D", f"{switch}=1"]
    subprocess.run(["mojo", "build", *defines, "-I", str(ROOT), str(source), "-o", str(binary), *mojo_cpu_args()], check=True)
    return binary


def pinned(command: list) -> list:
    """`command` pinned to one core where `taskset` exists."""
    return ["taskset", "-c", "2", *command] if shutil.which("taskset") else command


def main() -> None:
    """Builds every configuration, runs the rounds, fails on any cost that differs, prints the table."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--rounds", type=int, default=3)
    arguments = parser.parse_args()
    missing = [name for name in SETS if not (SAMPLES / f"{name}.seq").exists()]
    if missing:
        sys.exit(f"missing samples {missing}: run `pixi run bench-astarpa2` first")
    batch_bench.generate()
    pairs = DATA / "batches" / "illumina.tsv"
    BUILD.mkdir(parents=True, exist_ok=True)
    aligners = {switch: build(HERE / "ablation.mojo", switch) for switch in ["NONE", *SWITCHES]}
    batchers = {switch: build(HERE / "batches.mojo", switch) for switch in ["NONE", *BATCH_SWITCHES]}
    times = defaultdict(list)
    answers = defaultdict(set)
    files = [str(SAMPLES / f"{name}.seq") for name in SETS]
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
    print(f"Median time, ms, and each switch's change against it, over {arguments.rounds} rounds ({run.build_note()}):")
    columns = SWITCHES + BATCH_SWITCHES
    print(f"| dataset | baseline | {' | '.join(switch.removeprefix('ABLATE_').lower() for switch in columns)} |")
    print("| :-- | --: |" + " --: |" * len(columns))
    for name in SETS + batch_bench.WORKLOADS:
        base = statistics.median(times[("NONE", name)])
        cells = []
        for switch in columns:
            found = times.get((switch, name))
            cells.append(f"{100 * (statistics.median(found) / base - 1):+.1f}%" if found else "")
        print(f"| {name} | {base:.4g} | {' | '.join(cells)} |")


if __name__ == "__main__":
    main()
