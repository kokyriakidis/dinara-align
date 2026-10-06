#!/usr/bin/env python3
"""The results of `benchmarks/README.md`, after A*PA2's own (https://curiouscoding.nl/posts/astarpa2/#results).

    pixi run results-astarpa2 collect     # measure: forty minutes the first time, then dinara-align alone, minutes
    pixi run results-astarpa2 tables      # the measurements as Markdown tables, in .cache/results/astarpa2-results.md

Every aligner runs as A*PA2's evaluation ran them: one single-threaded job at a time, each pair aligned
once with its traceback, the time the average wall clock per alignment. One table for each of A*PA2's
figures:

- the real datasets, each aligner's mean and median time per alignment;
- 100 kbp uniform pairs from 0 to 15% divergence, the mean time at each divergence;
- uniform pairs from 3 kbp to 1 Mbp at 5% and 15%, the mean time at each length;
- dinara-align's history, the same long reads aligned by the package as of the commit that added each method.

The rivals and dinara-align's older commits never change, so each runs on a sample once and is kept
for good (see `pa_bench.measure`); a later `collect` times today's dinara-align alone. Every cost is
checked against every other aligner's on the pairs both aligned, and the run fails on a disagreement. A rival that has not finished a pair when its budget runs out stops there: its numbers
cover the pairs it finished, which the table counts, and a rival that finished none shows a dash.
"""

import json
import math
import platform
import subprocess
import sys
import tarfile
import io
from pathlib import Path

from pa_bench import (
    DATA,
    KEPT,
    REAL,
    SEED,
    astarpa_settings,
    dinara,
    download,
    generate,
    load_kept,
    measure,
    published_costs,
    run_tool,
    sample,
)
from run import CACHE, HERE, RESULTS, cargo_runner, duration, fetch, mojo_runner, nightly_environment

import os

MEASURED = RESULTS / "astarpa2-results.json"
BUDGET = 20.0
"""Seconds each aligner may spend on each sample."""

REAL_BASES = {"ont-500k": 20_000_000, "ont-500k-genvar": 20_000_000}
"""A larger sample of the long reads, sixteen pairs or so, so their boxes have something to show."""

DIVERGENCE_LENGTH = 100_000
DIVERGENCE_RATES = [f"{rate / 100:g}" for rate in range(0, 16)]
"""As pa-generate names its files: `0`, `0.01`, ..., `0.1`."""
DIVERGENCE_TOTAL = 1_000_000
"""Ten 100 kbp pairs a divergence, as A*PA2's divergence plot."""

LENGTH_RATES = ["0.05", "0.15"]
LENGTHS = [3_000, 10_000, 30_000, 100_000, 300_000, 1_000_000]

HISTORY = [
    ("641819a", "port"),
    ("5a6ee58", "+ diagonal transition"),
    ("3851282", "+ seed heuristic"),
    ("28c1241", "+ local pruning"),
    ("c53edc6", "+ two-ended search"),
    ("b60ece0", "+ real-read fixes"),
    ("501a6bb", "+ retries aimed"),
    ("ef0e73a", "+ inexact seeds"),
]
"""The commits that added each method, cumulatively, from the port of A*PA2-simple's band doubling."""

HISTORY_DATASETS = ["ont-500k", "ont-500k-genvar"]

DRIVER = """
from std.sys import argv
from std.time import perf_counter_ns

from dinara_align import edit_alignment


def main() raises:
    var tool = String(argv()[2])
    var budget = Float64(String(argv()[3]))
    var spun = perf_counter_ns()
    while perf_counter_ns() - spun < 50_000_000:
        pass
    var spent = 0.0
    for argument in range(4, len(argv())):
        var path = String(argv()[argument])
        var lines = open(path, "r").read().split("\\n")
        var index = 0
        while index + 1 < len(lines) and spent < budget:
            var first = String(lines[index][byte=1:])
            var second = String(lines[index + 1][byte=1:])
            var started = perf_counter_ns()
            var aligned = edit_alignment(first, second, 1)
            var seconds = Float64(perf_counter_ns() - started) / 1e9
            spent += seconds
            print(tool, path, seconds, Int(aligned.score), sep="\\t", flush=True)
            index += 2
"""
"""A runner for an older commit: only `edit_alignment`, which every commit has, on one thread."""

ALIGNERS = [dinara("1 thread"), "a*pa2-full", "a*pa2-simple", "a*pa", "edlib", "biwfa", "wfa"]
"""The single-threaded exact aligners, in the figures' fixed order and colours."""

APPROXIMATE = ["wfa-adaptive", "block-aligner"]
"""The approximate aligners A*PA2's evaluation sets beside the exact ones on the real datasets, with
its parameters. They may return a worse alignment than the optimum, so their costs are held against
the exact distance, never checked for agreement, and the share they get right goes beside their time."""

LABELS = {
    dinara("1 thread"): "dinara-align",
    "a*pa2-full": "A*PA2-full",
    "a*pa2-simple": "A*PA2-simple",
    "a*pa": "A*PA",
    "edlib": "Edlib",
    "biwfa": "BiWFA",
    "wfa": "WFA",
    "wfa-adaptive": "WFA-adaptive",
    "block-aligner": "Block Aligner",
}


# region Measuring


def runners() -> dict:
    ours = mojo_runner()
    fetch("pa-bench")
    wrapper = cargo_runner("pa-wrapper", dict(os.environ))
    nightly = nightly_environment(False)
    if nightly is None:
        sys.exit("A*PA needs `rustup` for its pinned nightly; run `pixi run bench-astarpa2 --install-rust` once")
    fetch("astar-pairwise-aligner")
    astarpa = cargo_runner("astarpa", nightly)
    return {"ours": ours, "astarpa": astarpa, "wrapper": wrapper}


def aligners(dataset: str, binaries: dict) -> list[tuple[str, Path, str]]:
    return [
        (dinara("1 thread"), binaries["ours"], dinara("1 thread")),
        ("a*pa2-full", binaries["astarpa"], "a*pa2-full"),
        ("a*pa2-simple", binaries["astarpa"], "a*pa2-simple"),
        ("a*pa", binaries["astarpa"], astarpa_settings(dataset)),
        ("edlib", binaries["wrapper"], "edlib"),
        ("biwfa", binaries["wrapper"], "biwfa"),
        ("wfa", binaries["wrapper"], "wfa"),
    ]


def times(
    name: str,
    path: Path,
    reference,
    binaries: dict,
    kept: dict,
    astarpa_choices: list[str] | None = None,
    correct: dict | None = None,
) -> dict[str, list[float]]:
    """Every aligner's seconds per pair on one sample, its costs checked against the others' and the published ones.

    With `astarpa_choices`, A*PA runs with each of those settings and keeps the faster: on a sweep the
    evaluation's per-dataset settings do not reach, the fairest reading of A*PA.

    With `correct`, the approximate aligners run too, and `correct` takes, for each, how many of the
    pairs it finished it aligned at the exact distance, and how many it finished. A cost below the
    exact distance would be an error in the harness or the aligner, and stops the run.
    """
    found, seen = {}, []
    for column, binary, tool in aligners(name, binaries):
        print(f"{name}: {LABELS[column]} ...", file=sys.stderr, flush=True)
        choices = astarpa_choices if column == "a*pa" and astarpa_choices else [tool]
        best = None
        for choice in choices:
            rows, _, _ = measure(binary, choice, path, BUDGET, None if binary == binaries["ours"] else kept)
            KEPT.write_text(json.dumps(kept))
            mean = sum(t for t, _ in rows) / len(rows) if rows else math.inf
            if best is None or (len(rows), -mean) > (len(best[0]), -best[1]):
                best = (rows, mean)
        rows = best[0]
        found[column] = [seconds for seconds, _ in rows]
        seen.append((column, [cost for _, cost in rows]))
    if reference is not None:
        seen.append(("published", reference))
    for index, (column, costs) in enumerate(seen):
        for other, other_costs in seen[index + 1 :]:
            common = min(len(costs), len(other_costs))
            if costs[:common] != other_costs[:common]:
                sys.exit(f"DISAGREEMENT on {name}: {column} and {other}")
    if correct is not None:
        exact = max((costs for _, costs in seen), key=len)
        for column in APPROXIMATE:
            print(f"{name}: {LABELS[column]} ...", file=sys.stderr, flush=True)
            rows, _, _ = measure(binaries["wrapper"], column, path, BUDGET, kept)
            KEPT.write_text(json.dumps(kept))
            found[column] = [seconds for seconds, _ in rows]
            costs = [cost for _, cost in rows]
            if any(cost < best for cost, best in zip(costs, exact)):
                sys.exit(f"{LABELS[column]} reports a cost below the exact distance on {name}")
            correct[column] = [sum(cost == best for cost, best in zip(costs, exact)), min(len(costs), len(exact))]
    return found


def generated(total: int, rates: list[str], lengths: list[int]) -> list[Path]:
    """pa-generate's files for any rates, lengths and total, under pa-bench's names (see `pa_bench.generate`)."""
    generate([lengths[0]])  # pins the generator's lock and builds it
    generator = cargo_runner("pa-generate", dict(os.environ))
    directory = DATA / "generated"
    subprocess.run([str(generator), str(directory), str(SEED), str(total), *rates, "--", *map(str, lengths)], check=True)
    return [directory / f"Uniform-t{total}-n{length}-e{rate}.seq" for rate in rates for length in lengths]


def historic_runner(commit: str) -> Path:
    """The package as of `commit`, built with `DRIVER` into a runner of its own, once.

    Commits before 66a0b82 name a NEON register in `opaque`'s inline assembly, which x86 cannot
    allocate: on x86 their copy takes that commit's fix, `x` for `w`, and nothing else. A runner that
    cannot align one short pair is refused rather than timed.
    """
    root = CACHE / "history" / commit
    binary = root / "runner"
    if binary.exists():
        return binary
    root.mkdir(parents=True, exist_ok=True)
    archive = subprocess.run(
        ["git", "-C", str(HERE.parent), "archive", commit, "dinara_align"], check=True, capture_output=True
    ).stdout
    with tarfile.open(fileobj=io.BytesIO(archive)) as opened:
        opened.extractall(root)
    if platform.machine().lower() in ("x86_64", "amd64"):
        source = root / "dinara_align" / "edit_distance.mojo"
        source.write_text(source.read_text().replace('constraints="=w,0"', 'constraints="=x,0"'))
    (root / "driver.mojo").write_text(DRIVER)
    print(f"building the package as of {commit} ...", file=sys.stderr, flush=True)
    subprocess.run(["mojo", "build", "-I", str(root), str(root / "driver.mojo"), "-o", str(binary)], check=True)
    probe = root / "probe.seq"
    probe.write_text(">ACGTACGTAC\n<ACGAACGTAC\n")
    checked = subprocess.run([str(binary), "seq", "history", "10", str(probe)], capture_output=True, text=True, timeout=60)
    if not checked.stdout.strip().endswith("\t1"):
        binary.unlink()
        sys.exit(f"the package as of {commit} built a runner that cannot align one pair")
    return binary


def collect() -> None:
    download()
    published = published_costs()
    binaries = runners()
    kept = load_kept()
    measured = {"machine": f"{platform.machine()}, {platform.system()} {platform.release()}"}

    real = {}
    for name in REAL:
        files = sorted((DATA / name).glob("*.seq"))
        path, facts = sample(name, files, REAL_BASES.get(name, 2_000_000), published)
        correct = {}
        real[name] = {
            "pairs": facts["pairs"],
            "mean_length": facts["mean_length"],
            "times": times(name, path, facts["reference"], binaries, kept, correct=correct),
            "correct": correct,
        }
    measured["real"] = real

    divergence = {}
    for file in generated(DIVERGENCE_TOTAL, DIVERGENCE_RATES, [DIVERGENCE_LENGTH]):
        path, facts = sample(file.stem, [file], 2_000_000, published)
        rate = file.stem.rsplit("-e", 1)[1]
        divergence[rate] = {
            "pairs": facts["pairs"],
            "times": times(file.stem, path, None, binaries, kept, ["a*pa r=1 prune=both", "a*pa r=2 prune=both"]),
        }
    measured["divergence"] = divergence

    length = {rate: {} for rate in LENGTH_RATES}
    for file in generate(LENGTHS):
        path, facts = sample(file.stem, [file], 2_000_000, published)
        rate = file.stem.rsplit("-e", 1)[1]
        size = int(file.stem.split("-n")[1].split("-")[0])
        length[rate][str(size)] = {
            "pairs": facts["pairs"],
            "times": times(file.stem, path, facts["reference"], binaries, kept),
        }
    measured["length"] = length

    history = {}
    for name in HISTORY_DATASETS:
        files = sorted((DATA / name).glob("*.seq"))
        path, _ = sample(name, files, REAL_BASES[name], published)
        steps = {}
        for commit, label in HISTORY:
            # An old commit never changes, so its times are kept for good beside the rivals'.
            key = "\t".join(["history", path.name, commit])
            if key not in kept:
                print(f"{name}: dinara-align as of {commit} ({label}) ...", file=sys.stderr, flush=True)
                rows, peak, stopped = run_tool(historic_runner(commit), "history", path, 3 * BUDGET)
                kept[key] = {"budget": 3 * BUDGET, "rows": rows, "peak": peak, "stopped": stopped}
                KEPT.write_text(json.dumps(kept))
            steps[commit] = [seconds for seconds, _ in kept[key]["rows"]]
        history[name] = {"steps": steps, "current": real[name]["times"][dinara("1 thread")]}
    measured["history"] = history

    RESULTS.mkdir(parents=True, exist_ok=True)
    MEASURED.write_text(json.dumps(measured, indent=1))
    print(f"measured into {MEASURED}", file=sys.stderr)


# endregion Measuring

# region Tables


def cell(values: list[float], pairs: int) -> str:
    """An aligner's mean and median time per alignment, and how many pairs it finished if not all."""
    if not values:
        return "—"
    ordered = sorted(values)
    middle = len(ordered) // 2
    median = ordered[middle] if len(ordered) % 2 else (ordered[middle - 1] + ordered[middle]) / 2
    text = f"{duration(sum(values) / len(values))} ({duration(median)})"
    return text if len(values) == pairs else f"{text}, {len(values)}/{pairs}"


def mean_cell(values: list[float], pairs: int) -> str:
    if not values:
        return "—"
    text = duration(sum(values) / len(values))
    return text if len(values) == pairs else f"{text} ({len(values)}/{pairs})"


def header(first: str, columns: list[str]) -> list[str]:
    return [
        f"| {first} | " + " | ".join(columns) + " |",
        "| :-- | " + " | ".join("--:" for _ in columns) + " |",
    ]


def tables(measured: dict) -> str:
    """The results as Markdown tables, one a figure of A*PA2's."""
    names = [LABELS[column] for column in ALIGNERS]
    out = ["### Real datasets", "", "Mean time per alignment, median in brackets; for the approximate aligners, marked with an", "asterisk, also the share of the pairs they finished that they aligned optimally:", ""]
    out += header("dataset", ["pairs", "mean length"] + names + [LABELS[column] + "*" for column in APPROXIMATE])
    for name, data in measured["real"].items():
        cells = [cell(data["times"].get(column, []), data["pairs"]) for column in ALIGNERS]
        for column in APPROXIMATE:
            text = cell(data["times"].get(column, []), data["pairs"])
            right, finished = data.get("correct", {}).get(column, [0, 0])
            cells.append(f"{text}, {100 * right / finished:.0f}% optimal" if finished else text)
        out.append(f"| {name} | {data['pairs']} | {data['mean_length'] / 1000:.3g} kbp | " + " | ".join(cells) + " |")

    out += ["", f"### Divergence, {DIVERGENCE_LENGTH // 1000} kbp pairs", "", "Mean time per alignment:", ""]
    out += header("divergence", names)
    for rate in sorted(measured["divergence"], key=float):
        data = measured["divergence"][rate]
        cells = [mean_cell(data["times"].get(column, []), data["pairs"]) for column in ALIGNERS]
        out.append(f"| {float(rate) * 100:.0f}% | " + " | ".join(cells) + " |")

    for rate in LENGTH_RATES:
        out += ["", f"### Length, {float(rate) * 100:.0f}% divergence", "", "Mean time per alignment:", ""]
        out += header("length", names)
        for size in sorted(measured["length"][rate], key=int):
            data = measured["length"][rate][size]
            cells = [mean_cell(data["times"].get(column, []), data["pairs"]) for column in ALIGNERS]
            size_label = f"{int(size) / 1e6:g} Mbp" if int(size) >= 1_000_000 else f"{int(size) / 1000:g} kbp"
            out.append(f"| {size_label} | " + " | ".join(cells) + " |")

    out += ["", "### dinara-align's history on the long reads", "", "Mean time per alignment, median in brackets:", ""]
    out += header("commit", [f"{name}" for name in HISTORY_DATASETS])
    for commit, label in HISTORY:
        cells = []
        for name in HISTORY_DATASETS:
            data = measured["history"][name]
            cells.append(cell(data["steps"].get(commit, []), len(data["current"])))
        out.append(f"| `{commit}` {label} | " + " | ".join(cells) + " |")
    current = [cell(measured["history"][name]["current"], len(measured["history"][name]["current"])) for name in HISTORY_DATASETS]
    out.append("| today | " + " | ".join(current) + " |")
    return "\n".join(out) + "\n"


def write_tables() -> None:
    measured = json.loads(MEASURED.read_text())
    text = tables(measured)
    (RESULTS / "astarpa2-results.md").write_text(text)
    print(text)


# endregion Tables


if __name__ == "__main__":
    step = sys.argv[1] if len(sys.argv) > 1 else ""
    if step == "collect":
        collect()
    elif step == "tables":
        write_tables()
    else:
        sys.exit(__doc__)
