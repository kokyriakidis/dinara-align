#!/usr/bin/env python3
"""Times dinara-align head to head with A*PA2's evaluation, on its own datasets.

    pixi run bench-astarpa2                 # real datasets and synthetic pairs up to 1 Mbp
    pixi run bench-astarpa2 --fresh         # the rivals too, rather than their kept results
    pixi run bench-astarpa2 --full          # adds the 3 and 10 Mbp synthetic pairs
    pixi run bench-astarpa2 --bases 2e7     # a larger sample of each dataset
    pixi run bench-astarpa2 --budget 20     # seconds each tool may spend on each dataset

The datasets are A*PA2's (https://curiouscoding.nl/posts/astarpa2/), as pa-bench defines them:

- the real ones, Oxford Nanopore reads and SARS-CoV-2 genomes, downloaded from pa-bench's `datasets`
  release;
- the synthetic ones regenerated exactly, by pa-generate at the commit pa-bench locked when the
  evaluation ran: seed 31415, uniform errors at 5% and 15%, 10 Mbp of pairs per length.

Each dataset is sampled: its pairs in a fixed random order, as many as fit 2 Mbp of sequence and at
least four, the same sample for every tool. Each tool aligns the sample's pairs with their traceback,
once each, as pa-bench times them, and the table gives the average time per alignment. A tool that
has not finished when its budget runs out stops there, its average covering the pairs it aligned,
which the table counts; a pair still running a second past the budget is stopped, and a tool that
finished none reads as more than the budget. Every cost is checked against every other tool's on the
pairs both aligned, and against the costs A*PA2's published results recorded for the same pairs; a
disagreement fails the run.

The aligners are the evaluation's exact ones: Edlib, BiWFA, A*PA, A*PA2-simple and A*PA2-full, with
its parameters. Like it, the times here are wall-clock on one thread; dinara-align also runs on all.

The rivals are pinned, so each one's results are kept and reused while its binary, the sample and
the budget stay the same; only dinara-align's columns run every time, unless `--fresh` asks for all.

A second table gives each tool's peak resident memory over the same run, from the operating system's
account of the process. Its first row, one 8 bp pair, is what each runner holds before any aligner's
work; a run stopped mid-pair had reached at least what it shows.
"""

import argparse
import json
import os
import random
import subprocess
import sys
import threading
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

SAMPLE_BASES = 2_000_000
"""The sequence a dataset's sample holds by default, both sides of every pair counted."""

SAMPLE_PAIRS = 4
"""The fewest pairs a sample holds, however long: one 500 kbp read says little of a dataset whose
pairs take from 60 ms to 1.5 s."""

GRACE = 1.0
"""The seconds past its budget a runner may take, starting up included, before it is stopped."""

KEPT = CACHE / "pa-bench-rivals.json"
"""The rivals' results from earlier runs, by tool, sample, budget and binary."""


def dinara(threads: str) -> str:
    return f"dinara-align (bit-parallel, {threads})"


def tools(dataset: str, ours: Path, astarpa: Path, wrapper: Path) -> list[tuple[str, Path, str]]:
    """Each column's runner and the tool name it is given, as the evaluation ran them on a dataset."""
    return [
        (dinara("1 thread"), ours, dinara("1 thread")),
        (dinara("8 threads"), ours, dinara("8 threads")),
        ("a*pa2-full", astarpa, "a*pa2-full"),
        ("a*pa2-simple", astarpa, "a*pa2-simple"),
        ("a*pa", astarpa, astarpa_settings(dataset)),
        ("edlib", wrapper, "edlib"),
        ("biwfa", wrapper, "biwfa"),
    ]


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


def sample(name: str, files: list[Path], bases: int, published: dict[str, list[int]]) -> tuple[Path, dict]:
    """A dataset's sample as one file, and its pair count, mean sequence length and published costs.

    pa-bench's files often run from close pairs to divergent ones, so the pairs are shuffled first:
    a sample is then unbiased, and so are the first pairs of it a tool stopped by its budget aligned.
    The costs are in the sample's order, or None unless every pair has one. A sample is kept beside
    its dataset and made again only when the dataset or the size asked for changes.
    """
    directory = DATA / "samples"
    path = directory / f"{name}-b{bases}-p{SAMPLE_PAIRS}.seq"
    facts = path.with_suffix(".json")
    newest = max(file.stat().st_mtime for file in files)
    if path.exists() and facts.exists() and path.stat().st_mtime > newest:
        return path, json.loads(facts.read_text())
    pairs = []
    for file in files:
        key = f"{file.parent.name}/{file.name}"
        lines = file.read_text().splitlines()
        for index in range(0, len(lines) - 1, 2):
            pairs.append((key, index // 2, lines[index], lines[index + 1]))
    random.Random(SEED).shuffle(pairs)
    chosen, letters = [], 0
    for pair in pairs:
        if len(chosen) >= SAMPLE_PAIRS and letters + len(pair[2]) + len(pair[3]) - 2 > bases:
            break
        chosen.append(pair)
        letters += len(pair[2]) + len(pair[3]) - 2
    reference = None
    if all(key in published for key, _, _, _ in chosen):
        reference = [published[key][index] for key, index, _, _ in chosen]
    directory.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{first}\n{second}\n" for _, _, first, second in chosen))
    found = {"pairs": len(chosen), "of": len(pairs), "mean_length": letters / (2 * len(chosen)), "reference": reference}
    facts.write_text(json.dumps(found))
    return path, found


# endregion Inputs

# region Runs


def run_tool(binary: Path, tool: str, path: Path, budget: float) -> tuple[list[tuple[float, int]], int, bool]:
    """One tool over one sample: each pair's time and cost it finished, its peak resident bytes, and
    whether it was stopped mid-pair.

    The runner prints a row as each pair finishes and checks its budget only between pairs, so one
    still busy a grace past the budget is stopped, keeping the rows it printed. `wait4` reports the
    process's peak memory, even when stopped, where `subprocess` would discard it.
    """
    process = subprocess.Popen(
        [str(binary), "seq", tool, str(budget), str(path)], stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True
    )
    stopped = threading.Event()

    def stop() -> None:
        stopped.set()
        process.kill()

    timer = threading.Timer(budget + GRACE, stop)
    timer.start()
    output = process.stdout.read()
    process.stdout.close()
    _, status, usage = os.wait4(process.pid, 0)
    timer.cancel()
    process.returncode = os.waitstatus_to_exitcode(status)
    # macOS counts peak memory in bytes, Linux in kilobytes.
    peak = usage.ru_maxrss if sys.platform == "darwin" else usage.ru_maxrss * 1024
    rows = []
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) == 4 and fields[0] == tool:
            rows.append((float(fields[2]), int(fields[3])))
    return rows, peak, stopped.is_set()


def measure(
    binary: Path, tool: str, path: Path, budget: float, kept: dict | None
) -> tuple[list[tuple[float, int]], int, bool]:
    """`run_tool`, or a rival's result from an earlier run when `kept` holds one for the same binary,
    sample and budget: the rivals are pinned, so only dinara-align's columns change between runs."""
    if kept is None:
        return run_tool(binary, tool, path, budget)
    key = "\t".join([tool, path.name, f"{budget:g}", str(binary.stat().st_mtime_ns), str(path.stat().st_mtime_ns)])
    if key not in kept:
        kept[key] = run_tool(binary, tool, path, budget)
    rows, peak, stopped = kept[key]
    return [tuple(row) for row in rows], peak, stopped


def megabytes(peak: int, stopped: bool) -> str:
    return f"{'≥ ' if stopped else ''}{peak / 2**20:.0f} MB"


# endregion Runs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--full", action="store_true", help="add the 3 and 10 Mbp synthetic pairs")
    parser.add_argument("--budget", type=float, default=5.0, help="seconds per tool per dataset (default 5)")
    parser.add_argument(
        "--bases", type=float, default=SAMPLE_BASES, help=f"sequence per dataset's sample (default {SAMPLE_BASES:.0e})"
    )
    parser.add_argument("--install-rust", action="store_true", help="install A*PA's nightly under the cache")
    parser.add_argument("--only", nargs="*", default=[], help="run only datasets whose names contain one of these")
    parser.add_argument("--fresh", action="store_true", help="run the rivals again rather than reuse their results")
    options = parser.parse_args()
    bases = int(options.bases)

    download()
    datasets = [(name, sorted((DATA / name).glob("*.seq"))) for name in REAL]
    datasets += [(file.stem, [file]) for file in generate(LENGTHS + (FULL_LENGTHS if options.full else []))]
    if options.only:
        datasets = [(name, files) for name, files in datasets if any(part in name for part in options.only)]
    published = published_costs()

    ours = mojo_runner()
    fetch("pa-bench")
    wrapper = cargo_runner("pa-wrapper", dict(os.environ))
    nightly = nightly_environment(options.install_rust)
    if nightly is None:
        sys.exit("A*PA needs `rustup` for its pinned nightly; pass --install-rust")
    fetch("astar-pairwise-aligner")
    astarpa = cargo_runner("astarpa", nightly)

    kept = json.loads(KEPT.read_text()) if KEPT.exists() and not options.fresh else {}

    def run(binary: Path, tool: str, path: Path) -> tuple[list[tuple[float, int]], int, bool]:
        return measure(binary, tool, path, options.budget, None if binary == ours else kept)

    columns = [column for column, _, _ in tools("", ours, astarpa, wrapper)]
    lines = [
        "| dataset | pairs | mean length | " + " | ".join(columns) + " | agree |",
        "| :-- | --: | --: | " + " | ".join("--:" for _ in columns) + " | :-: |",
    ]
    memory_lines = [
        "| dataset | " + " | ".join(columns) + " |",
        "| :-- | " + " | ".join("--:" for _ in columns) + " |",
    ]

    # What each runner holds aligning almost nothing: its runtime, before any aligner's work.
    tiny = DATA / "samples" / "baseline.seq"
    tiny.parent.mkdir(parents=True, exist_ok=True)
    if not tiny.exists():
        tiny.write_text(">ACGTACGT\n<ACGAACGT\n")
    peaks = [megabytes(*run(binary, tool, tiny)[1:]) for _, binary, tool in tools("ont", ours, astarpa, wrapper)]
    memory_lines.append("| one 8 bp pair | " + " | ".join(peaks) + " |")

    failed = []
    for name, files in datasets:
        path, facts = sample(name, files, bases, published)
        pairs = facts["pairs"]
        cells, peaks, seen = [], [], []
        for column, binary, tool in tools(name, ours, astarpa, wrapper):
            print(f"{name}: {column} ...", file=sys.stderr, flush=True)
            rows, peak, stopped = run(binary, tool, path)
            peaks.append(megabytes(peak, stopped))
            seen.append((column, [cost for _, cost in rows]))
            if not rows:
                cells.append(f"> {options.budget:g} s")
                continue
            cell = duration(sum(seconds for seconds, _ in rows) / len(rows))
            cells.append(cell if len(rows) == pairs else f"{cell} ({len(rows)}/{pairs})")

        # Every tool's costs against every other's and the published ones, on the pairs both have.
        agree = True
        if facts["reference"] is not None:
            seen.append(("published", facts["reference"]))
        for index, (column, costs) in enumerate(seen):
            for other, other_costs in seen[index + 1 :]:
                common = min(len(costs), len(other_costs))
                if costs[:common] != other_costs[:common]:
                    agree = False
                    print(f"DISAGREEMENT on {name}: {column} and {other}", file=sys.stderr)
        if not agree:
            failed.append(name)

        lines.append(
            f"| {name} | {pairs} of {facts['of']} | {facts['mean_length'] / 1000:.3g} kbp | "
            + " | ".join(cells)
            + f" | {'✓' if agree else '✗'} |"
        )
        memory_lines.append(f"| {name} | " + " | ".join(peaks) + " |")
        print(lines[-1], file=sys.stderr, flush=True)
        KEPT.write_text(json.dumps(kept))

    table = "\n".join(lines) + "\n\nPeak resident memory over the same runs:\n\n" + "\n".join(memory_lines) + "\n"
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "astarpa2.md").write_text(table)
    print(table)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
