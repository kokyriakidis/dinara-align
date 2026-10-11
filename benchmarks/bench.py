#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Times dinara-align against every exact global DNA aligner on the datasets of the two latest evaluations,
A*PA2's and QuickEd's, and checks that every answer agrees.

    pixi run bench                       # unit costs, every dataset, 5 s per tool per dataset
    pixi run bench --affine 4,6,2        # affine costs, as WFA counts them: mismatch, opening, extension
    pixi run bench --only ont l10000.    # the datasets whose names hold one of these
    pixi run bench --budget 60 --bases 2e7
    pixi run bench --fresh               # the rivals too, rather than their kept results

The datasets:

- A*PA2's real ones (Groot Koerkamp, Bioinformatics 2025), as pa-bench releases them: ONT reads of about
  1, 10, 50 and 500 kbp, the 500 kbp ones also against a reference with genetic variation, and SARS-CoV-2
  genomes. QuickEd's MiniION-a and MiniION-b are `ont-500k` and `ont-500k-genvar`.
- QuickEd's simulated ones (Doblas et al., Bioinformatics 2025), from Zenodo: random sequences of 10 kbp,
  100 kbp, 500 kbp and 1 Mbp with 1, 5, 10 and 20% uniformly spread errors, 100 Mbp a file, 3.2 GB in all,
  downloaded the first time. They cover the lengths and divergences A*PA2's synthetic sets did.

Each dataset is sampled: its pairs in a fixed random order, as many as fit `--bases` of sequence, 2 Mbp by
default, and at least four, the same sample for every tool. Each tool aligns the sample's pairs with their
traceback, once each, single-threaded, until its budget runs out, a pair still running a second past it
stopped. One run gives both evaluations' views: A*PA2's mean time per alignment and memory growth per pair,
and QuickEd's total over the whole dataset, the mean times the dataset's pairs, exact where the sample is
the dataset and marked ≈ where it is not. Every cost is checked against every other tool's on the pairs
both aligned, and on A*PA2's datasets against the costs its published results recorded; a disagreement
fails the run.

The aligners are every exact one either evaluation compared: Edlib, BiWFA, WFA2-lib keeping every front,
A*PA with the settings A*PA2's evaluation gave each kind of data, A*PA2-simple and A*PA2-full; on x86-64
KSW2, QuickEd's bound-and-align and BSAlign's striped kernels with their band off, so exact. WFA, KSW2
and BSAlign, whose memory grows with the matrix or the square of the distance, run under an 8 GB cap.
At affine costs, the ones that take them: WFA, BiWFA, and on x86-64 KSW2 and BSAlign. dinara-align runs
pair by pair, as the others do, and as a batch, every pair of the sample in one call on the same thread,
whose column is a throughput rather than a latency.

The rivals are pinned, so each runs on a sample once and its results are kept (see `harness.measure`); a
rerun times dinara-align alone.
"""

import argparse
import json
import math
import os
import random
import statistics
import subprocess
import sys
import zipfile
from pathlib import Path
from urllib.request import urlretrieve

import harness
from harness import (
    CACHE,
    DATA,
    HERE,
    KEPT,
    RESULTS,
    build_note,
    c_cpu_flag,
    cargo_runner,
    duration,
    fetch,
    growth,
    load_kept,
    measure,
    megabytes,
    mojo_runner,
    nightly_environment,
    set_cpu,
    x86,
)

PA_BENCH = DATA / "pa-bench"
PUBLISHED = CACHE / "astarpa2-evals"
RELEASE = "https://github.com/pairwise-alignment/pa-bench/releases/download"
REAL = ["ont-1k", "ont-10k", "ont-50k", "ont-500k", "ont-500k-genvar", "sars-cov-2"]
"""A*PA2's real datasets, each a directory of `.seq` files."""

QUICKED = DATA / "quicked"
ZENODO = "https://zenodo.org/api/records/14526794/files"
"""The one version of QuickEd's datasets' record (https://doi.org/10.5281/zenodo.14526793)."""
SIMULATED = [
    (length, pairs, rate)
    for length, pairs in [(10_000, 10_000), (100_000, 1_000), (500_000, 200), (1_000_000, 100)]
    for rate in (1, 5, 10, 20)
]
"""QuickEd's simulated datasets: each file's length, pairs and percent error."""

SEED = 31415
SAMPLE_BASES = 2_000_000
"""The sequence a dataset's sample holds by default, both sides of every pair counted."""
SAMPLE_PAIRS = 4
"""The fewest pairs a sample holds, however long."""


# region Inputs


def download() -> list[tuple[str, list[Path]]]:
    """Every dataset by name and its `.seq` files, downloading what is missing, and A*PA2's published results."""
    PA_BENCH.mkdir(parents=True, exist_ok=True)
    found = []
    for name in REAL:
        directory = PA_BENCH / name
        if not (directory.exists() and any(directory.glob("*.seq"))):
            archive = PA_BENCH / f"{name}.zip"
            print(f"downloading {name} ...", file=sys.stderr, flush=True)
            urlretrieve(f"{RELEASE}/datasets/{name}.zip", archive)
            with zipfile.ZipFile(archive) as opened:
                opened.extractall(directory)
            archive.unlink()
        found.append((name, sorted(directory.glob("*.seq"))))
    if not (PUBLISHED / "results").exists():
        PUBLISHED.mkdir(parents=True, exist_ok=True)
        archive = PUBLISHED / "results.zip"
        print("downloading A*PA2's published results ...", file=sys.stderr, flush=True)
        urlretrieve(f"{RELEASE}/astarpa2-evals/results.zip", archive)
        with zipfile.ZipFile(archive) as opened:
            opened.extractall(PUBLISHED)
    QUICKED.mkdir(parents=True, exist_ok=True)
    for length, pairs, rate in SIMULATED:
        name = f"l{length}.n{pairs}.e{rate}"
        path = QUICKED / f"{name}.seq"
        if not path.exists():
            print(f"downloading {name} (200 MB) ...", file=sys.stderr, flush=True)
            partial = path.with_suffix(".part")
            urlretrieve(f"{ZENODO}/{name}.seq/content", partial)
            partial.rename(path)
        found.append((name, [path]))
    return found


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
            if not output["Ok"]["costs"] or "File" not in result["job"]["dataset"]:
                continue
            file = Path(result["job"]["dataset"]["File"])
            costs.setdefault(f"{file.parent.name}/{file.name}", output["Ok"]["costs"])
    return costs


def sample(name: str, files: list[Path], bases: int, published: dict[str, list[int]]) -> tuple[Path, dict]:
    """A dataset's sample as one file, and its pair count, the dataset's, the mean sequence length and the
    published costs.

    Files often run from close pairs to divergent ones, so the pairs are shuffled first: a sample is then
    unbiased, and so are the first pairs of it a tool stopped by its budget aligned. The costs are in the
    sample's order, or None unless every pair has one. A sample is kept beside its dataset and made again
    only when the dataset or the size asked for changes.
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
    found = {
        "pairs": len(chosen),
        "of": len(pairs),
        "mean_length": letters / (2 * len(chosen)),
        "reference": reference,
        "lengths": [max(len(first), len(second)) - 1 for _, _, first, second in chosen],
    }
    facts.write_text(json.dumps(found))
    return path, found


# endregion Inputs

# region Tools


def dinara(threads: str) -> str:
    """The bit-parallel column's tool name, `threads` saying how it runs."""
    return f"dinara-align (bit-parallel, {threads})"


def astarpa_settings(dataset: str) -> str:
    """A*PA as A*PA2's evaluation ran it: GCSH with diagonal transition, its `r` and pruning per kind of data,
    its synthetic settings for QuickEd's simulated sets by their error rate."""
    if dataset.startswith("ont"):
        return "a*pa r=2 prune=start"
    if dataset == "sars-cov-2":
        return "a*pa r=1 prune=start"
    rate = int(dataset.rsplit(".e", 1)[1])
    return "a*pa r=1 prune=both" if rate <= 5 else "a*pa r=2 prune=both"


def quicked_runner() -> Path | None:
    """QuickEd's runner (`quicked/runner.c`) against its library at the pinned commit, or None where it does
    not build: its sources take x86 vector extensions."""
    if not x86():
        return None
    source = fetch("QuickEd")
    build = CACHE / "target-quicked"
    configure = ["cmake", "-S", str(source), "-B", str(build), "-DCMAKE_BUILD_TYPE=Release"]
    # Its bundled Edlib still asks for a CMake older than 3.5, which current CMake refuses without this.
    configure.append("-DCMAKE_POLICY_VERSION_MINIMUM=3.5")
    if harness.CPU != "native":
        configure.append("-DQUICKED_NONATIVE=ON")
    subprocess.run(configure, check=True, capture_output=True)
    subprocess.run(["cmake", "--build", str(build), "-j", "4", "--target", "quicked"], check=True, capture_output=True)
    binary = build / "quicked-runner"
    subprocess.run(
        ["cc", "-O3", c_cpu_flag(), f"-I{source}", f"-I{source / 'quicked'}", str(HERE / "quicked" / "runner.c"),
         str(source / "lib" / "libquicked.a"), "-lm", "-o", str(binary)],
        check=True,
    )
    return binary


def bsalign_runners() -> tuple[Path, Path] | None:
    """BSAlign's runners (`bsalign/runner.c`) at the pinned commit, its edit distance's and its affine
    alignment's, or None where they do not build: its kernels are x86 vector code. Its affine kernel stops on
    AVX2, so that runner keeps to SSE4.2, as BSAlign's own build does by default."""
    if not x86():
        return None
    source = fetch("bsalign")
    build = CACHE / "target-bsalign"
    build.mkdir(exist_ok=True)
    common = ["cc", "-O3", c_cpu_flag(), "-mpopcnt", "-D_GNU_SOURCE", f"-I{source}", str(HERE / "bsalign" / "runner.c")]
    libraries = ["-lm", "-lz", "-lpthread"]
    edit, affine = build / "edit-runner", build / "affine-runner"
    subprocess.run([*common, *libraries, "-o", str(edit)], check=True)
    subprocess.run([*common, "-mno-avx", *libraries, "-o", str(affine)], check=True)
    return edit, affine


def unit_tools(dataset: str, runners: dict) -> list[tuple[str, Path, str]]:
    """The columns at unit costs: each column's name, runner and the tool name it passes."""
    chosen = [
        (dinara("1 thread"), runners["ours"], dinara("1 thread")),
        (dinara("batch, 1 thread"), runners["ours"], dinara("batch, 1 thread")),
        ("a*pa2-full", runners["astarpa"], "a*pa2-full"),
        ("a*pa2-simple", runners["astarpa"], "a*pa2-simple"),
        ("a*pa", runners["astarpa"], astarpa_settings(dataset)),
        ("edlib", runners["wrapper"], "edlib"),
        ("biwfa", runners["wrapper"], "biwfa"),
        ("wfa", runners["wrapper"], "wfa"),
    ]
    if x86():
        # KSW2 at unit costs, a mismatch 1 and a gap of `k` letters `0 + k`, as QuickEd's evaluation runs it.
        chosen.append(("ksw2", runners["wrapper"], "ksw2:1,0,1"))
    if runners["quicked"] is not None:
        chosen.append(("quicked", runners["quicked"], "quicked"))
    if runners["bsalign"] is not None:
        chosen.append(("bsalign", runners["bsalign"][0], "bsalign-edit"))
    return chosen


def affine_tools(costs: str, runners: dict) -> list[tuple[str, Path, str]]:
    """The columns at affine costs `x,o,e`, as WFA counts them: the exact aligners that take them."""
    batch = "dinara-align (batch, 1 thread)"
    chosen = [
        ("dinara-align (1 thread)", runners["ours"], f"dinara-align:{costs}"),
        (batch, runners["ours"], f"{batch}:{costs}"),
        ("wfa", runners["wrapper"], f"wfa:{costs}"),
        ("biwfa", runners["wrapper"], f"biwfa:{costs}"),
    ]
    if x86():
        chosen.append(("ksw2", runners["wrapper"], f"ksw2:{costs}"))
    if runners["bsalign"] is not None:
        chosen.append(("bsalign", runners["bsalign"][1], f"bsalign:{costs}"))
    return chosen


# endregion Tools


def lengths_of(path: Path, facts: dict) -> list[int]:
    """Each pair's longer sequence in a sample, from its facts or, for a sample made before they kept them,
    from the file."""
    if "lengths" in facts:
        return facts["lengths"]
    lines = path.read_text().splitlines()
    return [max(len(lines[index]), len(lines[index + 1])) - 1 for index in range(0, len(lines) - 1, 2)]


def spread(values: list[float], unit: str = "") -> str:
    """A minimum, mean and maximum, as QuickEd's dataset table gives them."""
    if not values:
        return "—"
    return f"{min(values):,.0f} / {sum(values) / len(values):,.0f} / {max(values):,.0f}{unit}"


def setup_note(runners: dict) -> str:
    """What every table rests on beside the CPU: each rival's pinned commit and the compilers' versions, as
    QuickEd's supplement (S1.1) lists them."""

    def version(command: list) -> str:
        try:
            return subprocess.run(command, capture_output=True, text=True, check=True).stdout.splitlines()[0].strip()
        except (OSError, subprocess.CalledProcessError, IndexError):
            return "unknown"

    pins = ", ".join(f"{name} `{commit[:7]}`" for name, (_, commit) in harness.RIVALS.items())
    return (
        f"Rivals at {pins}. Compilers: {version(['cc', '--version'])}; {version(['rustc', '--version'])}; "
        f"{version(['mojo', '--version'])}. One thread each, pinned on Linux by the caller.\n\n"
    )


def main() -> None:
    """Runs every tool on a sample of every dataset, checks their costs agree with each other's and the
    published ones, and writes the tables; exits nonzero when any disagree."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--budget", type=float, default=5.0, help="seconds per tool per dataset (default 5)")
    parser.add_argument("--bases", type=float, default=SAMPLE_BASES, help=f"sequence per sample (default {SAMPLE_BASES:.0e})")
    parser.add_argument("--only", nargs="*", default=[], help="run only datasets whose names hold one of these")
    parser.add_argument("--affine", metavar="X,O,E", help="align at affine costs: WFA's mismatch, opening and extension")
    parser.add_argument("--fresh", action="store_true", help="run the rivals again rather than reuse their results")
    parser.add_argument("--install-rust", action="store_true", help="install A*PA's nightly under the cache")
    parser.add_argument("--cpu", default="native", help="the CPU every tool is built for (default: the host's)")
    options = parser.parse_args()
    set_cpu(options.cpu)
    harness.REMEASURE = options.fresh
    bases = int(options.bases)

    datasets = download()
    if options.only:
        datasets = [(name, files) for name, files in datasets if any(part in name for part in options.only)]
    published = published_costs()

    fetch("pa-bench")
    nightly = nightly_environment(options.install_rust)
    if nightly is None:
        sys.exit("A*PA needs `rustup` for its pinned nightly; pass --install-rust")
    fetch("astar-pairwise-aligner")
    runners = {
        "ours": mojo_runner(),
        "wrapper": cargo_runner("pa-wrapper", dict(os.environ)),
        "astarpa": cargo_runner("astarpa", nightly),
        "quicked": quicked_runner(),
        "bsalign": bsalign_runners(),
    }
    kept = load_kept()

    def chosen(dataset: str) -> list[tuple[str, Path, str]]:
        """The columns for `dataset`: the affine ones when `--affine` asks, else the unit-cost ones."""
        return affine_tools(options.affine, runners) if options.affine else unit_tools(dataset, runners)

    columns = [column for column, _, _ in chosen("ont-1k")]
    ours = columns[0]
    rivals = columns[2:]
    head = " | ".join(columns)
    rule = " | ".join("--:" for _ in columns)
    data_lines = [
        "| dataset | pairs | length, min / mean / max (bp) | divergence, min / mean / max (%) |",
        "| :-- | --: | --: | --: |",
    ]
    time_lines = [f"| dataset | pairs | {head} | agree |", f"| :-- | --: | {rule} | :-: |"]
    median_lines = [f"| dataset | pairs | {head} |", f"| :-- | --: | {rule} |"]
    total_lines = [f"| dataset | pairs | {head} |", f"| :-- | --: | {rule} |"]
    speed_lines = [f"| dataset | {' | '.join(rivals)} |", f"| :-- | {' | '.join('--:' for _ in rivals)} |"]
    memory_lines = [f"| dataset | {head} |", f"| :-- | {rule} |"]
    peak_lines = memory_lines[:2]
    speedups = {rival: [] for rival in rivals}
    scaling = []
    failed = []
    for name, files in datasets:
        path, facts = sample(name, files, bases, published)
        pairs, whole = facts["pairs"], facts["of"]
        times, medians, totals, grown, peaks, seen = [], [], [], [], [], []
        measured = {}
        for column, binary, tool in chosen(name):
            print(f"{name}: {column} ...", file=sys.stderr, flush=True)
            rows, peak, status = measure(binary, tool, path, options.budget, None if binary == runners["ours"] else kept)
            measured[column] = (rows, status)
            peaks.append(megabytes(peak, status))
            grown.append(growth(rows))
            seen.append((column, [row[1] for row in rows]))
            if not rows:
                marker = {harness.MEMORY: "o.o.m.", harness.FAILED: "failed"}.get(status, f"t.o. ({options.budget:g} s)")
                times.append(marker)
                medians.append(marker)
                totals.append(marker)
                continue
            mean = sum(row[0] for row in rows) / len(rows)
            median = statistics.median(row[0] for row in rows)
            if len(rows) == pairs:
                times.append(duration(mean))
                medians.append(duration(median))
            else:
                why = "o.o.m. after " if status == harness.MEMORY else ""
                times.append(f"{duration(mean)} ({why}{len(rows)}/{pairs})")
                medians.append(f"{duration(median)} ({why}{len(rows)}/{pairs})")
            if status == harness.MEMORY:
                totals.append("o.o.m.")
            else:
                totals.append(duration(mean * whole) if len(rows) == whole else f"≈ {duration(mean * whole)}")
        # Every tool's costs against every other's, on the pairs both have, and the published ones at unit costs.
        agree = True
        if facts["reference"] is not None and not options.affine:
            seen.append(("published", facts["reference"]))
        for index, (column, costs) in enumerate(seen):
            for other, other_costs in seen[index + 1 :]:
                common = min(len(costs), len(other_costs))
                if costs[:common] != other_costs[:common]:
                    agree = False
                    print(f"DISAGREEMENT on {name}: {column} and {other}", file=sys.stderr)
        if not agree:
            failed.append(name)
        # The dataset as QuickEd's Table S1 describes one: its pairs' lengths and, at unit costs, their
        # divergence, the agreed distance over the longer sequence.
        lengths = lengths_of(path, facts)
        known = max((costs for _, costs in seen), key=len)
        divergence = [100 * cost / max(length, 1) for cost, length in zip(known, lengths)] if not options.affine else []
        data_lines.append(f"| {name} | {pairs} of {whole} | {spread(lengths)} | {spread(divergence)} |")
        # dinara-align against each rival on the pairs both aligned, the rival's time over dinara-align's.
        ours_rows = measured[ours][0]
        cells = []
        for rival in rivals:
            rival_rows = measured[rival][0]
            common = min(len(ours_rows), len(rival_rows))
            if common == 0:
                cells.append("—")
                continue
            ratio = sum(row[0] for row in rival_rows[:common]) / sum(row[0] for row in ours_rows[:common])
            speedups[rival].append(ratio)
            cells.append(f"{ratio:.1f}×" if ratio < 100 else f"{ratio:.0f}×")
        speed_lines.append(f"| {name} | " + " | ".join(cells) + " |")
        time_lines.append(f"| {name} | {pairs} of {whole} | " + " | ".join(times) + f" | {'✓' if agree else '✗'} |")
        median_lines.append(f"| {name} | {pairs} of {whole} | " + " | ".join(medians) + " |")
        total_lines.append(f"| {name} | {whole} | " + " | ".join(totals) + " |")
        memory_lines.append(f"| {name} | " + " | ".join(grown) + " |")
        peak_lines.append(f"| {name} | " + " | ".join(peaks) + " |")
        if name.startswith("l") and ".e" in name:
            length, rate = int(name[1:].split(".")[0]), int(name.rsplit(".e", 1)[1])
            scaling.append({
                "length": length,
                "error": rate,
                "pairs": pairs,
                "seconds_per_pair": {
                    column: sum(row[0] for row in rows) / len(rows) for column, (rows, _) in measured.items() if rows
                },
                "finished": {column: len(rows) for column, (rows, _) in measured.items()},
            })
        print(time_lines[-1], file=sys.stderr, flush=True)
        KEPT.write_text(json.dumps(kept))

    def geometric(values: list[float]) -> float:
        return math.exp(sum(math.log(value) for value in values) / len(values))

    summary = [
        f"| {rival} | {min(values):.1f}× | {geometric(values):.1f}× | {max(values):.1f}× | {len(values)} |"
        for rival, values in speedups.items()
        if values
    ]
    table = (
        build_note()
        + setup_note(runners)
        + "Datasets, over each sample (QuickEd's Table S1):\n\n"
        + "\n".join(data_lines)
        + f"\n\nMean time per alignment over each dataset's sample, at most {options.budget:g} s a tool a dataset, the pairs\n"
        + "a tool aligned before it was stopped counted; t.o. a tool its budget stopped before a pair, o.o.m. one the\n"
        + "8 GB memory cap stopped (A*PA2's view):\n\n"
        + "\n".join(time_lines)
        + "\n\nMedian time per alignment over the same pairs (A*PA2's supplement):\n\n"
        + "\n".join(median_lines)
        + "\n\nTotal over each whole dataset, its mean time per pair times its pairs, ≈ where the sample is not the\n"
        + "whole (QuickEd's Table 2 and S3):\n\n"
        + "\n".join(total_lines)
        + f"\n\nEach rival's time over {ours}'s, on the pairs both aligned:\n\n"
        + "\n".join(speed_lines)
        + "\n\nThe same across datasets, as QuickEd's text summarizes its own:\n\n"
        + "| rival | least | geometric mean | most | datasets |\n| :-- | --: | --: | --: | --: |\n"
        + "\n".join(summary)
        + "\n\nGrowth of peak memory aligning a pair, median / largest, input and runtime left out (A*PA2's Table 10):\n\n"
        + "\n".join(memory_lines)
        + "\n\nPeak resident memory of each runner, its sample included, a run stopped mid-pair at least this\n"
        + "(QuickEd's Table 2 and S3):\n\n"
        + "\n".join(peak_lines)
        + "\n"
    )
    RESULTS.mkdir(parents=True, exist_ok=True)
    stem = f"bench-affine-{options.affine.replace(',', '-')}" if options.affine else "bench"
    (RESULTS / f"{stem}.md").write_text(table)
    # Time per pair against length at each error rate, and against error at each length, over QuickEd's
    # simulated sets (its Figures S6 and S7), for the paper's charts.
    (RESULTS / f"{stem}-scaling.json").write_text(json.dumps(scaling, indent=1))
    print(table)
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
