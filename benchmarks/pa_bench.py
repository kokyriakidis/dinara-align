#!/usr/bin/env python3
"""Times dinara-align head to head with A*PA2's evaluation, on its own datasets.

    pixi run bench-astarpa2                 # real datasets and synthetic pairs up to 1 Mbp
    pixi run bench-astarpa2 --full          # adds the 3 and 10 Mbp synthetic pairs
    pixi run bench-astarpa2 --budget 60     # seconds each tool may spend on each dataset

The datasets are A*PA2's (https://curiouscoding.nl/posts/astarpa2/), as pa-bench defines them:

- the real ones, Oxford Nanopore reads and SARS-CoV-2 genomes, downloaded from pa-bench's `datasets`
  release;
- the synthetic ones regenerated exactly, by pa-generate at the commit pa-bench locked when the
  evaluation ran: seed 31415, uniform errors at 5% and 15%, 10 Mbp of pairs per length.

Each tool aligns every pair with its traceback, once, as pa-bench times them, and the table gives the
average time per alignment. A tool that has not finished a dataset when its budget runs out stops
there, and its average covers the pairs it aligned, which the table counts. Every cost is checked
against every other tool's on the pairs both aligned, and against the costs A*PA2's published results
recorded for the same files; a disagreement fails the run.

The aligners are the evaluation's exact ones: Edlib, BiWFA, A*PA, A*PA2-simple and A*PA2-full, with
its parameters. Like it, the times here are wall-clock on one thread; dinara-align also runs on all.
"""

import argparse
import json
import os
import subprocess
import sys
import zipfile
from pathlib import Path
from urllib.request import urlretrieve

from run import CACHE, HERE, RESULTS, cargo_runner, duration, fetch, mojo_runner, nightly_environment

DATA = CACHE / "data" / "pa-bench"
PUBLISHED = CACHE / "astarpa2-evals"
RELEASE = "https://github.com/pairwise-alignment/pa-bench/releases/download"

REAL = ["ont-1k", "ont-10k", "ont-50k", "ont-500k", "ont-500k-genvar", "sars-cov-2"]
"""The evaluation's real datasets, each a directory of `.seq` files."""

SEED = 31415
TOTAL_SIZE = 10_000_000
RATES = ["0.05", "0.15"]
LENGTHS = [3_000, 10_000, 30_000, 100_000, 300_000, 1_000_000]
FULL_LENGTHS = [3_000_000, 10_000_000]
"""The evaluation's synthetic datasets (`scaling-n.yaml`): one file per rate and length."""

EVALUATION_LOCK = "78c1cf7"
"""The pa-bench commit whose `Cargo.lock` the evaluation's data was generated with."""


def dinara(threads: str) -> str:
    return f"dinara-align (bit-parallel, {threads})"


def astarpa_settings(dataset: str) -> str:
    """A*PA as the evaluation ran it on each dataset: GCSH with diagonal transition, `r` and pruning per data."""
    if dataset.startswith("ont"):
        return "a*pa r=2 prune=start"
    if dataset == "sars-cov-2":
        return "a*pa r=1 prune=start"
    return "a*pa r=1 prune=both" if "-e0.05" in dataset else "a*pa r=2 prune=both"


# region Inputs


def download() -> None:
    """The real datasets, from pa-bench's release, and A*PA2's published results to check against."""
    DATA.mkdir(parents=True, exist_ok=True)
    for name in REAL:
        directory = DATA / name
        if directory.exists() and any(directory.glob("*.seq")):
            continue
        archive = DATA / f"{name}.zip"
        print(f"downloading {name} ...", file=sys.stderr, flush=True)
        urlretrieve(f"{RELEASE}/datasets/{name}.zip", archive)
        with zipfile.ZipFile(archive) as opened:
            opened.extractall(directory)
        archive.unlink()
    if not (PUBLISHED / "results").exists():
        PUBLISHED.mkdir(parents=True, exist_ok=True)
        archive = PUBLISHED / "results.zip"
        print("downloading A*PA2's published results ...", file=sys.stderr, flush=True)
        urlretrieve(f"{RELEASE}/astarpa2-evals/results.zip", archive)
        with zipfile.ZipFile(archive) as opened:
            opened.extractall(PUBLISHED)


def generate(lengths: list[int]) -> list[Path]:
    """The synthetic files, by the evaluation's own generator, under the names pa-bench gives them."""
    pa_bench = fetch("pa-bench")
    crate = HERE / "pa-generate"
    lock = subprocess.run(
        ["git", "-C", str(pa_bench), "show", f"{EVALUATION_LOCK}:Cargo.lock"], check=True, capture_output=True
    ).stdout
    (crate / "Cargo.lock").write_bytes(lock)
    generator = cargo_runner("pa-generate", dict(os.environ))
    directory = DATA / "generated"
    subprocess.run(
        [str(generator), str(directory), str(SEED), str(TOTAL_SIZE), *RATES, "--", *map(str, lengths)], check=True
    )
    return [directory / f"Uniform-t{TOTAL_SIZE}-n{length}-e{rate}.seq" for rate in RATES for length in lengths]


def published_costs() -> dict[str, list[int]]:
    """Every file's costs as A*PA2's published results recorded them, by dataset and file name."""
    costs = {}
    for path in (PUBLISHED / "results").glob("*.json"):
        if path.name == "cache.json":
            continue
        for result in json.loads(path.read_text()):
            output = result.get("output", {})
            if "Ok" not in output or result["job"]["costs"] != {"sub": 1, "open": 0, "extend": 1}:
                continue
            # Only an exact aligner's costs are the optimum: not Block Aligner's, nor WFA-adaptive's.
            algorithm = result["job"]["algo"]
            if not any(name in algorithm for name in ("Edlib", "AstarPa", "AstarPa2", "Wfa")):
                continue
            if "Wfa" in algorithm and str(algorithm["Wfa"].get("heuristic", "None")) != "None":
                continue
            if not output["Ok"]["costs"]:
                continue
            dataset = result["job"]["dataset"]
            if "File" in dataset:
                file = Path(dataset["File"])
                key = f"{file.parent.name}/{file.name}"
            elif "Generated" in dataset:
                generated = dataset["Generated"]
                key = (
                    f"generated/{generated['error_model']}-t{generated['total_size']}"
                    f"-n{generated['length']}-e{generated['error_rate']}.seq"
                )
            else:
                continue
            costs.setdefault(key, output["Ok"]["costs"])
    return costs


# endregion Inputs

# region Runs


def run_tool(binary: Path, tool: str, files: list[Path], budget: float) -> dict[str, tuple[int, float, list[int]]]:
    """One tool over one dataset's files: per file, the pairs it aligned, their time, and their costs.

    A single pair can outlast the budget, so the runner is also stopped once it overruns the budget
    by a margin; the files it finished by then still count.
    """
    limit = 3 * budget + 60
    try:
        output = subprocess.run(
            [str(binary), "seq", tool, str(budget), *map(str, files)], capture_output=True, text=True, timeout=limit
        ).stdout
    except subprocess.TimeoutExpired as stopped:
        output = stopped.stdout.decode() if isinstance(stopped.stdout, bytes) else (stopped.stdout or "")
    rows = {}
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) != 5 or fields[0] != tool:
            continue
        costs = [int(cost) for cost in fields[4].split(",")] if fields[4] else []
        rows[f"{Path(fields[1]).parent.name}/{Path(fields[1]).name}"] = (int(fields[2]), float(fields[3]), costs)
    return rows


def pair_count(files: list[Path]) -> tuple[int, float]:
    """A dataset's pairs and their mean length, over both sequences."""
    pairs = 0
    letters = 0
    for file in files:
        with open(file) as opened:
            for line in opened:
                letters += len(line) - 2
                pairs += line.startswith(">")
    return pairs, letters / max(2 * pairs, 1)


# endregion Runs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--full", action="store_true", help="add the 3 and 10 Mbp synthetic pairs")
    parser.add_argument("--budget", type=float, default=20.0, help="seconds per tool per dataset (default 20)")
    parser.add_argument("--install-rust", action="store_true", help="install A*PA's nightly under the cache")
    options = parser.parse_args()

    download()
    datasets = [(name, sorted((DATA / name).glob("*.seq"))) for name in REAL]
    datasets += [(file.stem, [file]) for file in generate(LENGTHS + (FULL_LENGTHS if options.full else []))]
    published = published_costs()

    ours = mojo_runner()
    fetch("pa-bench")
    wrapper = cargo_runner("pa-wrapper", dict(os.environ))
    nightly = nightly_environment(options.install_rust)
    if nightly is None:
        sys.exit("A*PA needs `rustup` for its pinned nightly; pass --install-rust")
    fetch("astar-pairwise-aligner")
    astarpa = cargo_runner("astarpa", nightly)

    columns = [dinara("1 thread"), dinara("8 threads"), "a*pa2-full", "a*pa2-simple", "a*pa", "edlib", "biwfa"]
    lines = [
        "| dataset | pairs | mean length | " + " | ".join(columns) + " | agree |",
        "| :-- | --: | --: | " + " | ".join("--:" for _ in columns) + " | :-: |",
    ]
    failed = []
    for name, files in datasets:
        pairs, mean_length = pair_count(files)
        tools = [
            (dinara("1 thread"), ours, dinara("1 thread")),
            (dinara("8 threads"), ours, dinara("8 threads")),
            ("a*pa2-full", astarpa, "a*pa2-full"),
            ("a*pa2-simple", astarpa, "a*pa2-simple"),
            ("a*pa", astarpa, astarpa_settings(name)),
            ("edlib", wrapper, "edlib"),
            ("biwfa", wrapper, "biwfa"),
        ]
        results = {}
        for column, binary, tool in tools:
            print(f"{name}: {column} ...", file=sys.stderr, flush=True)
            results[column] = run_tool(binary, tool, files, options.budget)

        # Every tool's costs against every other's and the published ones, on the pairs both have.
        agree = True
        for file in files:
            key = f"{file.parent.name}/{file.name}"
            reference = published.get(key)
            seen = [(column, rows[key][2]) for column, rows in results.items() if key in rows]
            if reference is not None:
                seen.append(("published", reference))
            for index, (column, costs) in enumerate(seen):
                for other, other_costs in seen[index + 1 :]:
                    common = min(len(costs), len(other_costs))
                    if costs[:common] != other_costs[:common]:
                        agree = False
                        print(f"DISAGREEMENT on {key}: {column} and {other}", file=sys.stderr)
        if not agree:
            failed.append(name)

        cells = []
        for column in columns:
            rows = results[column].values()
            done = sum(row[0] for row in rows)
            seconds = sum(row[1] for row in rows)
            if done == 0:
                cells.append("—")
                continue
            cell = duration(seconds / done)
            cells.append(cell if done == pairs else f"{cell} ({done}/{pairs})")
        lines.append(
            f"| {name} | {pairs} | {mean_length / 1000:.3g} kbp | " + " | ".join(cells) + f" | {'✓' if agree else '✗'} |"
        )
        print(lines[-1], file=sys.stderr, flush=True)

    table = "\n".join(lines) + "\n"
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "astarpa2.md").write_text(table)
    print(table)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
