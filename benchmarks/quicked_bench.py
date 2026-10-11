#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Times dinara-align on QuickEd's datasets as QuickEd's paper measures them (Doblas et al.,
Bioinformatics 2025, Supplementary Table S3): each dataset whole, its total time and peak memory.

    pixi run bench-quicked                  # every dataset, 300 s per tool per dataset
    pixi run bench-quicked --budget 1200
    pixi run bench-quicked --only l100000 miniion

The datasets are QuickEd's:

- its simulated pairs, downloaded from Zenodo (https://doi.org/10.5281/zenodo.14526793): random
  sequences of 10 kbp, 100 kbp, 500 kbp and 1 Mbp with 1, 5, 10 and 20% uniformly spread errors, 100
  Mbp of sequence a file, 3.2 GB in all;
- ONT MiniION-a and MiniION-b, which are A*PA2's `ont-500k` and `ont-500k-genvar` (pa-bench's release,
  as `pa_bench.py` downloads them), each aligned whole here where `pa_bench.py` samples them.

QuickEd's other real datasets, its Illumina, PacBio HiFi, ONT UltraLong and PromethION pairs, were not
published as pairs, so they are left out.

Each tool aligns every pair of a dataset once with its traceback, single-threaded and, on Linux, pinned
by the harness's runners as `pa_bench.py` runs them (see `pa_bench.run_tool`), until its budget runs
out. The time table gives the total seconds over the dataset, as QuickEd's does; a tool stopped by its
budget shows its pairs' mean times the dataset's pairs, marked with the pairs it finished, an estimate.
QuickEd's memory figure is each process's peak, its runners reading a pair at a time. The harness's
runners hold their dataset whole, up to 200 MB, so the first memory table gives instead what aligning a
pair adds to the most the process held before it, the largest over the dataset's pairs (A*PA2's
measure; see `pa_bench.growth`), which leaves the input and runtime out alike; the second gives the
whole processes' peaks (GNU time's on Linux), input included, a runner stopped mid-pair having reached
at least that. Every cost is checked against
every other tool's on the pairs both aligned; a disagreement fails the run. The rivals' runs are kept
as `pa_bench.py` keeps them, so a rerun times dinara-align alone.
"""

import argparse
import json
import os
import sys
from pathlib import Path
from urllib.request import urlretrieve

import pa_bench
import run
from pa_bench import KEPT, dinara, load_kept, measure, megabytes, quicked_runner
from run import CACHE, RESULTS, build_note, cargo_runner, fetch, mojo_runner, nightly_environment, set_cpu

DATA = CACHE / "data" / "quicked"
# The record of the DOI's one version, which holds the files.
ZENODO = "https://zenodo.org/api/records/14526794/files"
SIMULATED = [
    (length, pairs, rate)
    for length, pairs in [(10_000, 10_000), (100_000, 1_000), (500_000, 200), (1_000_000, 100)]
    for rate in (1, 5, 10, 20)
]
"""QuickEd's simulated datasets: length, pairs and percent error of each Zenodo file."""
REAL = {"miniion-a": "ont-500k", "miniion-b": "ont-500k-genvar"}
"""QuickEd's ONT MiniION-a and MiniION-b, as the pa-bench datasets they are."""


def datasets() -> list[tuple[str, Path]]:
    """Every dataset as one `.seq` file, downloading or gathering what is missing."""
    DATA.mkdir(parents=True, exist_ok=True)
    found = []
    for length, pairs, rate in SIMULATED:
        name = f"l{length}.n{pairs}.e{rate}.seq"
        path = DATA / name
        if not path.exists():
            print(f"downloading {name} (200 MB) ...", file=sys.stderr, flush=True)
            partial = path.with_suffix(".part")
            urlretrieve(f"{ZENODO}/{name}/content", partial)
            partial.rename(path)
        found.append((path.stem, path))
    pa_bench.download()
    for name, directory in REAL.items():
        path = DATA / f"{name}.seq"
        if not path.exists():
            path.write_text("".join(file.read_text().rstrip("\n") + "\n" for file in sorted((pa_bench.DATA / directory).glob("*.seq"))))
        found.append((name, path))
    return found


def seconds(value: float) -> str:
    """A total time's cell: QuickEd's tables give seconds with one decimal."""
    return f"{value:.1f}" if value >= 0.1 else f"{value:.3f}"


def main() -> None:
    """Runs every tool on every dataset whole, checks their costs agree, and writes the time and memory
    tables; exits nonzero when any disagree."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--budget", type=float, default=300.0, help="seconds per tool per dataset (default 300)")
    parser.add_argument("--only", nargs="*", default=[], help="run only datasets whose names contain one of these")
    parser.add_argument("--install-rust", action="store_true", help="install A*PA's nightly under the cache")
    parser.add_argument("--cpu", default="native", help="the CPU every tool is built for (default: the host's)")
    options = parser.parse_args()
    set_cpu(options.cpu)

    chosen_sets = datasets()
    if options.only:
        chosen_sets = [(name, path) for name, path in chosen_sets if any(part in name for part in options.only)]

    ours = mojo_runner()
    fetch("pa-bench")
    wrapper = cargo_runner("pa-wrapper", dict(os.environ))
    nightly = nightly_environment(options.install_rust)
    if nightly is None:
        sys.exit("A*PA needs `rustup` for its pinned nightly; pass --install-rust")
    fetch("astar-pairwise-aligner")
    astarpa = cargo_runner("astarpa", nightly)
    quicked = quicked_runner()
    kept = load_kept()

    def columns(name: str) -> list[tuple[str, Path, str]]:
        """Every tool QuickEd's table compares, dinara-align's pairs one at a time, as the others run; A*PA
        with the settings A*PA2's evaluation gave its datasets, the ONT ones for the MiniIONs."""
        named = REAL.get(name, name)
        return [column for column in pa_bench.tools(named, ours, astarpa, wrapper, quicked) if "batch" not in column[0]]

    names = [column for column, _, _ in columns("")]
    time_lines = ["| dataset | pairs | " + " | ".join(names) + " | agree |", "| :-- | --: | " + " | ".join("--:" for _ in names) + " | :-: |"]
    memory_lines = ["| dataset | " + " | ".join(names) + " |", "| :-- | " + " | ".join("--:" for _ in names) + " |"]
    peak_lines = memory_lines[:2]
    failed = []
    for name, path in chosen_sets:
        pairs = path.read_bytes().count(b"\n") // 2
        cells, added, peaks, seen = [], [], [], []
        for column, binary, tool in columns(name):
            print(f"{name}: {column} ...", file=sys.stderr, flush=True)
            rows, peak, stopped = measure(binary, tool, path, options.budget, None if binary == ours else kept)
            peaks.append(megabytes(peak, stopped))
            grown = [row[2] for row in rows if len(row) > 2 and row[2] >= 0]
            added.append(f"{max(grown) / 2**20:.1f} MB" if grown else "—")
            seen.append((column, [row[1] for row in rows]))
            if not rows:
                cells.append(f"> {options.budget:g} s, none")
                continue
            total = sum(row[0] for row in rows)
            cells.append(seconds(total) if len(rows) == pairs else f"≈ {seconds(total / len(rows) * pairs)} ({len(rows)}/{pairs})")
        agree = True
        for index, (column, costs) in enumerate(seen):
            for other, other_costs in seen[index + 1 :]:
                common = min(len(costs), len(other_costs))
                if costs[:common] != other_costs[:common]:
                    agree = False
                    print(f"DISAGREEMENT on {name}: {column} and {other}", file=sys.stderr)
        if not agree:
            failed.append(name)
        time_lines.append(f"| {name} | {pairs} | " + " | ".join(cells) + f" | {'✓' if agree else '✗'} |")
        memory_lines.append(f"| {name} | " + " | ".join(added) + " |")
        peak_lines.append(f"| {name} | " + " | ".join(peaks) + " |")
        print(time_lines[-1], file=sys.stderr, flush=True)
        KEPT.write_text(json.dumps(kept))
    table = (
        build_note()
        + f"Total time over each dataset, seconds, at most {options.budget:g} s a tool a dataset; ≈ marks a tool its budget\n"
        + "stopped, its mean time per pair times the dataset's pairs, with the pairs it finished:\n\n"
        + "\n".join(time_lines)
        + "\n\nMemory an alignment adds, the largest over each dataset's pairs, input and runtime left out:\n\n"
        + "\n".join(memory_lines)
        + "\n\nPeak resident memory of each runner, its copy of the dataset included:\n\n"
        + "\n".join(peak_lines)
        + "\n"
    )
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "quicked.md").write_text(table)
    print(table)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
