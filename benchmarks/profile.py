#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Pinpoints where dinara-align spends its time and where a rival beats it, dataset by dataset and pair by
pair, on the samples `bench.py` times.

    pixi run profile duel                    # every pair against every rival: where dinara-align loses
    pixi run profile duel --affine 4,6,2
    pixi run profile probes --only ont-10k   # the engines each pair took and the work each did
    pixi run profile perf --only ont-10k     # Linux: hardware counters and the hottest functions and lines
    pixi run profile compare main~1          # this checkout against a commit, pair by pair, with intervals
    pixi run profile compare a1b2c3d main    # one commit against another
    pixi run profile all                     # duel, probes and perf, one report

`duel` times dinara-align over `--rounds` rounds and sets each pair's median beside every rival's time for
the same pair, the rivals' kept results (`harness.measure`) or a fresh run of them. A dataset's mean can
win while single pairs lose: the report lists, for each dataset, the margin over the fastest rival, the
pairs some rival aligned faster, those it aligned faster than dinara-align's fastest round (beyond its
own noise), the time those losses cost, and the worst of them with their length, cost and divergence, the
place to look first.

`probes` builds `probed.mojo` with `-D DINARA_PROBES` (see `dinara_align/probes.mojo`), which records per
pair the engine that settled it, the time each phase took and the work it did: diagonal-transition steps,
band rounds and the bounds they tried against the distance found, tiles swept and recomputed, seeds built
and pruned, wavefront costs and sweep cells. Around each pair it also reads the processor's counters
(instructions, cycles, branch, cache and TLB misses, through `perf_event_open` on Linux), and with
`probes/counting.c` preloaded the allocations and the bytes copied, moved and set through the C library,
and the page faults. A normal build compiles none of it. The report gives each dataset's time by phase and
its pairs by engine, the work against what the pair needed, the processor's and the allocator's view per
pair, and flags the pairs whose probes point at a choice that cost time (see `probes_report`).

`perf` (Linux, `perf` installed, `perf_event_paranoid` at most 2) runs the runner built with line tables
under `perf stat` and `perf record` on each sample: cycles, instructions per cycle, branch and cache
misses, and the share of samples by function, by source file (each file one engine: `diagonal.mojo` the
diagonal transition, `bit_parallel.mojo` the band, and so on) and by source line, inlined code charged
to the line it came from. Samples in the C library's copy, set and allocation routines are summed, and,
where the processor records its last branches (Intel's LBR), each is shown with the code that called it.

`compare` builds the runner from two commits, or a commit and this checkout as it stands, from the same
runner source, and times both over `--rounds` rounds that alternate (A B, then B A), every pair's costs
held equal. Each dataset gets the ratio of the two mean times with a 95% bootstrap interval that resamples
the pairs and, for each pair, a round of each side, so it covers the rounds' noise and the pairs' spread;
a ratio whose interval excludes 1 is reported as faster or slower, any other as within noise. On x86 both
builds keep their jumps clear of 32-byte boundaries (see `ablation.mitigated_build`), without which a
rebuild alone moves times by several percent.

Every run pins itself to one core on Linux (`--core`, default 2) and refuses to start while that core, or
the machine, is busy (`--force` to override), since a shared core corrupts every timing silently. Reports
go to `benchmarks/.cache/results/profile-*.md`, with each pair's figures in a JSON beside it.
"""

import argparse
import json
import math
import os
import random
import re
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

import bench
import harness
from harness import CACHE, HERE, RESULTS, ROOT, duration, mojo_cpu_args, x86

BUILDS = CACHE / "profile"
"""The runners this script builds: with probes, with line tables, and one per commit compared."""
WORKTREES = CACHE / "worktrees"
"""A detached worktree per commit `compare` builds."""

ENGINE_FILES = {
    "diagonal.mojo": "diagonal transition",
    "slides.mojo": "diagonal transition",
    "edit_search.mojo": "diagonal transition",
    "search.mojo": "diagonal transition",
    "bit_parallel.mojo": "band sweep",
    "band.mojo": "band rounds",
    "seeds.mojo": "seed heuristic",
    "traceback.mojo": "traceback",
    "cigar.mojo": "CIGAR",
    "gap_affine.mojo": "gap-affine wavefront",
    "differences.mojo": "difference sweep",
    "anti_diagonals.mojo": "anti-diagonal sweep",
    "gotoh.mojo": "Gotoh",
    "scoring.mojo": "scoring",
    "modes.mojo": "modes",
    "api.mojo": "API",
    "common.mojo": "profile and buffers",
}
"""The engine each source file of the library holds, for `perf`'s view by file."""


# region The machine


def pin(core: int) -> None:
    """Pins this process, and so every runner it starts, to `core` on Linux; elsewhere does nothing."""
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, {core})


def core_busy(core: int, seconds: float = 1.0) -> float | None:
    """The share of `seconds` that `core` spent off idle, from `/proc/stat`, or None off Linux."""
    stat = Path("/proc/stat")
    if not stat.exists():
        return None

    def sample() -> tuple[int, int]:
        for line in stat.read_text().splitlines():
            fields = line.split()
            if fields[0] == f"cpu{core}":
                values = [int(value) for value in fields[1:]]
                # idle and iowait
                return sum(values), values[3] + values[4]
        raise RuntimeError(f"no cpu{core} in /proc/stat")

    total0, idle0 = sample()
    time.sleep(seconds)
    total1, idle1 = sample()
    return 1 - (idle1 - idle0) / max(total1 - total0, 1)


def check_quiet(core: int, force: bool) -> None:
    """Stops unless `core` is idle and the machine's load is low: two processes on a pinned core take turns,
    and every time measured then is wrong without any error. `force` warns instead."""
    problems = []
    busy = core_busy(core)
    if busy is not None and busy > 0.05:
        problems.append(f"core {core} was {100 * busy:.0f}% busy over a second")
    if hasattr(os, "getloadavg"):
        load = os.getloadavg()[0]
        cores = os.cpu_count() or 1
        if load > max(1.5, cores / 4):
            problems.append(f"load average {load:.2f}")
    if problems:
        message = "the machine is not quiet: " + "; ".join(problems)
        if not force:
            sys.exit(message + " (pass --force to time anyway)")
        print("WARNING: " + message, file=sys.stderr)


# endregion The machine

# region Builds


def mojo_build(root: Path, source: Path, binary: Path, extra: list = (), mitigated: bool = False) -> Path:
    """`source` built against the library at `root` into `binary`, kept while the command and every source are
    unchanged; with `mitigated` on x86, jumps kept clear of 32-byte boundaries (see `ablation`)."""
    import ablation

    binary.parent.mkdir(parents=True, exist_ok=True)
    command = ["mojo", "build", *extra, "-I", str(root), str(source), "-o", str(binary), *mojo_cpu_args()]
    stamp = binary.with_suffix(".command")
    sources = [source, *(root / "dinara_align").rglob("*.mojo")]
    signature = " ".join(command) + (" mitigated" if mitigated else "")
    if (
        binary.exists()
        and stamp.exists()
        and stamp.read_text() == signature
        and all(item.stat().st_mtime < binary.stat().st_mtime for item in sources)
    ):
        return binary
    if mitigated and x86() and shutil.which("clang") and shutil.which("cc"):
        saved = ablation.ROOT
        ablation.ROOT = root
        try:
            ablation.mitigated_build(source, binary, list(extra))
        finally:
            ablation.ROOT = saved
    else:
        subprocess.run(command, check=True)
    stamp.write_text(signature)
    return binary


def checkout(reference: str) -> tuple[Path, str]:
    """The library at `reference`, a commit or `.` for this checkout as it stands, and a label for it: a
    detached worktree under the cache, made once per commit."""
    if reference == ".":
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "dinara_align"],
                               capture_output=True, text=True, check=True).stdout.strip()
        head = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
        return ROOT, f"{head}{'+changes' if dirty else ''}"
    sha = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--verify", f"{reference}^{{commit}}"],
                         capture_output=True, text=True, check=True).stdout.strip()
    tree = WORKTREES / sha[:12]
    if not (tree / "dinara_align").exists():
        WORKTREES.mkdir(parents=True, exist_ok=True)
        subprocess.run(["git", "-C", str(ROOT), "worktree", "add", "--detach", str(tree), sha], check=True,
                       capture_output=True)
    return tree, sha[:9]


# endregion Builds

# region Samples


def samples(options) -> list[tuple[str, Path, dict]]:
    """The samples `bench.py` times, those `--only` names, with their facts."""
    datasets = bench.download(options.only)
    if options.only:
        datasets = [(name, files) for name, files in datasets if any(part in name for part in options.only)]
    published = bench.published_costs()
    found = []
    for name, files in datasets:
        path, facts = bench.sample(name, files, int(options.bases), published)
        facts = dict(facts, lengths=bench.lengths_of(path, facts))
        found.append((name, path, facts))
    return found


def ours_tool(options) -> str:
    """The name the runner takes for one pair at a time, at the costs asked for."""
    return f"dinara-align:{options.affine}" if options.affine else bench.dinara("1 thread")


def rounds_of(binary: Path, tool: str, path: Path, options) -> list[list[tuple]]:
    """`--rounds` runs of `binary` over the sample at `path`, each its rows (seconds, cost, growth[, more])."""
    runs = []
    for _ in range(options.rounds):
        rows, _, status = harness.run_tool(binary, tool, path, options.budget)
        if status == harness.FAILED:
            sys.exit(f"{binary.name} failed on {path.name}")
        runs.append(rows)
    return runs


def per_pair(runs: list[list[tuple]]) -> list[dict]:
    """The pairs every run finished, each with its median, fastest and slowest time and its cost, which must
    be the same in every run."""
    finished = min(len(run) for run in runs)
    pairs = []
    for index in range(finished):
        times = [run[index][0] for run in runs]
        costs = {run[index][1] for run in runs}
        if len(costs) != 1:
            sys.exit(f"pair {index}: the rounds found different costs {sorted(costs)}")
        pairs.append({"median": statistics.median(times), "fastest": min(times), "slowest": max(times),
                      "cost": runs[0][index][1]})
    return pairs


# endregion Samples

# region duel


def rival_tools(name: str, options, runners: dict) -> list[tuple[str, Path, str]]:
    """The rivals `bench.py` would run on dataset `name`, at the costs asked for."""
    chosen = bench.affine_tools(options.affine, runners) if options.affine else bench.unit_tools(name, runners)
    return [tool for tool in chosen if not tool[0].startswith("dinara-align")]


def rival_runners(options) -> dict:
    """Every rival's runner, as `bench.py` builds them."""
    bench.fetch("pa-bench")
    nightly = harness.nightly_environment(False)
    if nightly is None:
        sys.exit("A*PA needs `rustup` for its pinned nightly; run `pixi run bench --install-rust` once")
    bench.fetch("astar-pairwise-aligner")
    return {
        "ours": harness.mojo_runner(),
        "wrapper": harness.cargo_runner("pa-wrapper", dict(os.environ)),
        "astarpa": harness.cargo_runner("astarpa", nightly),
        "quicked": bench.quicked_runner(),
        "bsalign": bench.bsalign_runners(),
    }


def duel(options) -> str:
    """dinara-align against every rival, pair by pair, on every sample; returns the report."""
    harness.REMEASURE = options.fresh
    runners = rival_runners(options)
    kept = harness.load_kept()
    tool = ours_tool(options)
    lines = [
        f"# Duel: dinara-align against every rival, pair by pair{', at ' + options.affine if options.affine else ''}",
        "",
        harness.build_note().strip(),
        "",
        f"dinara-align's time per pair is its median over {options.rounds} rounds; each rival's is its kept time"
        " (one run). *Margin*: the fastest rival's mean over dinara-align's on the pairs both finished. *Lost*: pairs"
        " some rival aligned faster; *beyond noise*: faster than dinara-align's fastest round. *Cost of losses*: the"
        " time dinara-align would save by matching the fastest rival on every pair it lost, as a share of its total.",
        "",
        "| dataset | pairs | dinara-align mean | fastest rival | margin | lost | beyond noise | cost of losses |",
        "| :-- | --: | --: | :-- | --: | --: | --: | --: |",
    ]
    details, worst_lines, everything = [], [], {}
    for name, path, facts in samples(options):
        print(f"{name}: dinara-align, {options.rounds} rounds ...", file=sys.stderr, flush=True)
        ours = per_pair(rounds_of(runners["ours"], tool, path, options))
        rivals = {}
        for column, binary, rival in rival_tools(name, options, runners):
            if options.kept_only:
                key = "\t".join([rival, path.name, harness.fingerprint(path), harness.identity(binary)])
                if key not in kept:
                    continue
            print(f"{name}: {column} ...", file=sys.stderr, flush=True)
            rows, _, status = harness.measure(binary, rival, path, options.budget, kept)
            rivals[column] = (rows, status)
        harness.KEPT.write_text(json.dumps(kept))
        everything[name] = {"ours": ours, "rivals": {column: rows for column, (rows, _) in rivals.items()},
                            "lengths": facts["lengths"]}
        if not ours:
            lines.append(f"| {name} | {facts['pairs']} | t.o. | | | | | |")
            continue
        # The fastest rival by mean over the pairs both it and dinara-align finished.
        best, margin = None, math.inf
        for column, (rows, _) in rivals.items():
            common = min(len(rows), len(ours))
            if common == 0:
                continue
            ratio = sum(row[0] for row in rows[:common]) / sum(pair["median"] for pair in ours[:common])
            if ratio < margin:
                best, margin = column, ratio
        losses = []
        for index, pair in enumerate(ours):
            times = [(rows[index][0], column) for column, (rows, _) in rivals.items() if index < len(rows)]
            if not times:
                continue
            fastest, column = min(times)
            if fastest < pair["median"]:
                length = facts["lengths"][index] if index < len(facts["lengths"]) else 0
                losses.append({
                    "pair": index, "rival": column, "rival_time": fastest, "ours": pair["median"],
                    "ours_fastest": pair["fastest"], "cost": pair["cost"], "length": length,
                    "lost": pair["median"] - fastest, "beyond_noise": fastest < pair["fastest"],
                })
        total = sum(pair["median"] for pair in ours)
        lost = sum(loss["lost"] for loss in losses)
        beyond = sum(1 for loss in losses if loss["beyond_noise"])
        mean = total / len(ours)
        finished = f"{len(ours)}" if len(ours) == facts["pairs"] else f"{len(ours)}/{facts['pairs']}"
        lines.append(
            f"| {name} | {finished} | {duration(mean)} | {best or '—'} | "
            f"{f'{margin:.2f}×' if best else '—'} | {len(losses)} | {beyond} | {100 * lost / total:.1f}% |"
        )
        for loss in sorted(losses, key=lambda loss: -loss["lost"])[: options.worst]:
            divergence = 100 * loss["cost"] / max(loss["length"], 1)
            worst_lines.append(
                f"| {name} | {loss['pair']} | {loss['length']:,} | {loss['cost']:,} | {divergence:.1f}% | "
                f"{duration(loss['ours'])} | {duration(loss['rival_time'])} ({loss['rival']}) | "
                f"{loss['ours'] / loss['rival_time']:.2f}× | {'yes' if loss['beyond_noise'] else 'no'} |"
            )
        details.append((name, margin if best else math.inf))
    closest = sorted(details, key=lambda item: item[1])[:5]
    lines += [
        "",
        "Closest margins first: " + ", ".join(f"{name} {margin:.2f}×" for name, margin in closest if margin < math.inf),
        "",
        f"## The pairs lost by the most time (at most {options.worst} a dataset)",
        "",
        "Divergence is the cost over the longer sequence. Ratio: dinara-align's time over the rival's.",
        "",
        "| dataset | pair | length | cost | divergence | dinara-align | fastest rival | ratio | beyond noise |",
        "| :-- | --: | --: | --: | --: | --: | --: | --: | :-: |",
        *worst_lines,
    ]
    stem = f"profile-duel-{options.affine.replace(',', '-')}" if options.affine else "profile-duel"
    write(stem, "\n".join(lines) + "\n", everything)
    return "\n".join(lines) + "\n"


# endregion duel

# region compare


def bootstrap(a: list[list[float]], b: list[list[float]], draws: int = 2000) -> tuple[float, float, float]:
    """The ratio of B's mean time to A's, and its 95% interval: each draw resamples the pairs and, for each
    pair drawn, one round of each side. `a[i]` and `b[i]` are pair `i`'s times over the rounds."""
    point = sum(statistics.median(times) for times in b) / sum(statistics.median(times) for times in a)
    generator = random.Random(20261012)
    count = len(a)
    ratios = []
    for _ in range(draws):
        total_a = total_b = 0.0
        for _ in range(count):
            index = generator.randrange(count)
            total_a += generator.choice(a[index])
            total_b += generator.choice(b[index])
        ratios.append(total_b / total_a)
    ratios.sort()
    return point, ratios[int(0.025 * draws)], ratios[int(0.975 * draws) - 1]


def compare(options) -> str:
    """Two builds over the same samples, alternating rounds; returns the report."""
    base, base_label = checkout(options.base)
    head, head_label = checkout(options.head)
    source = HERE / "ours.mojo"
    mitigated = not options.unmitigated
    binaries = {}
    for label, root in ((base_label, base), (head_label, head)):
        binary = BUILDS / f"runner-{label}-{harness.CPU}"
        try:
            binaries[label] = mojo_build(root, source, binary, mitigated=mitigated)
        except subprocess.CalledProcessError:
            # An older library may not take this runner; its own runner then, which the report notes.
            print(f"{label}: this runner does not build against it; using its own", file=sys.stderr)
            binaries[label] = mojo_build(root, root / "benchmarks" / "ours.mojo", binary, mitigated=mitigated)
    tool = ours_tool(options)
    lines = [
        f"# Compare: {head_label} against {base_label}{', at ' + options.affine if options.affine else ''}",
        "",
        harness.build_note().strip(),
        "",
        f"Mean time per pair over {options.rounds} alternating rounds, and the ratio {head_label} / {base_label} with"
        " its 95% bootstrap interval over pairs and rounds. Every pair's costs are equal in both.",
        "",
        f"| dataset | pairs | {base_label} | {head_label} | ratio | 95% interval | verdict |",
        "| :-- | --: | --: | --: | --: | :-: | :-- |",
    ]
    everything = {}
    for name, path, facts in samples(options):
        runs = {base_label: [], head_label: []}
        for round in range(options.rounds):
            order = [base_label, head_label] if round % 2 == 0 else [head_label, base_label]
            for label in order:
                print(f"{name}: {label}, round {round + 1} of {options.rounds} ...", file=sys.stderr, flush=True)
                rows, _, status = harness.run_tool(binaries[label], tool, path, options.budget)
                if status == harness.FAILED:
                    sys.exit(f"{label} failed on {name}")
                runs[label].append(rows)
        finished = min(len(rows) for side in runs.values() for rows in side)
        if finished == 0:
            lines.append(f"| {name} | 0/{facts['pairs']} | t.o. | t.o. | | | |")
            continue
        costs = {label: [rows[index][1] for index in range(finished) for rows in side[:1]]
                 for label, side in runs.items()}
        if costs[base_label] != costs[head_label]:
            sys.exit(f"{name}: the two builds' costs differ")
        times = {label: [[rows[index][0] for rows in side] for index in range(finished)]
                 for label, side in runs.items()}
        ratio, low, high = bootstrap(times[base_label], times[head_label])
        verdict = "faster" if high < 1 else "slower" if low > 1 else "within noise"
        means = {label: sum(statistics.median(pair) for pair in times[label]) / finished for label in times}
        pairs = f"{finished}" if finished == facts["pairs"] else f"{finished}/{facts['pairs']}"
        lines.append(
            f"| {name} | {pairs} | {duration(means[base_label])} | {duration(means[head_label])} | {ratio:.3f} | "
            f"{low:.3f}–{high:.3f} | {verdict} |"
        )
        everything[name] = {"times": times, "costs": costs[base_label], "lengths": facts["lengths"][:finished]}
    write(f"profile-compare-{base_label}-{head_label}", "\n".join(lines) + "\n", everything)
    return "\n".join(lines) + "\n"


# endregion compare

# region perf

COUNTERS = ["cycles", "instructions", "branches", "branch-misses", "L1-dcache-load-misses", "LLC-load-misses"]


def perf_ready() -> str | None:
    """Why `perf` cannot run here, or None when it can."""
    if not sys.platform.startswith("linux"):
        return "perf runs on Linux only"
    if not shutil.which("perf"):
        return "perf is not installed"
    paranoid = Path("/proc/sys/kernel/perf_event_paranoid")
    if paranoid.exists() and int(paranoid.read_text()) > 2:
        return f"perf_event_paranoid is {paranoid.read_text().strip()}; at most 2 is needed"
    return None


def perf_stat(command: list) -> dict:
    """`perf stat`'s counts for `command`, by event."""
    output = subprocess.run(["perf", "stat", "-x", ",", "-e", ",".join(COUNTERS), "--", *command],
                            capture_output=True, text=True)
    counts = {}
    for line in output.stderr.splitlines():
        fields = line.split(",")
        if len(fields) > 2 and fields[0].strip().replace(".", "").isdigit():
            counts[fields[2].split(":")[0]] = float(fields[0])
    return counts


def perf_top(data: Path, sort: str, limit: int) -> list[tuple[float, str]]:
    """The `limit` largest shares of samples in a `perf record` file, by `sort` (`sym`, `srcfile`, `srcline`)."""
    output = subprocess.run(
        ["perf", "report", "-i", str(data), "--stdio", "--no-children", "--sort", sort, "--percent-limit", "0.3"],
        capture_output=True, text=True,
    ).stdout
    found = []
    for line in output.splitlines():
        match = re.match(r"\s*([\d.]+)%\s+(.*\S)", line)
        if match and not line.lstrip().startswith("#"):
            found.append((float(match.group(1)), match.group(2)))
    return found[:limit]


COPYING = re.compile(r"\b(__)?(mem(cpy|move|set)|malloc|calloc|realloc|free|cfree|_int_malloc|_int_free|"
                     r"posix_memalign|aligned_alloc|memalign|malloc_consolidate|KGEN_CompilerRT_\w*Alloc\w*)")
"""The C library's copy, set and allocation routines, and the Mojo runtime's allocator, by symbol."""


def copies_and_allocations(data: Path, traced: bool) -> list[str]:
    """The share of samples in copy, set and allocation routines, and with call graphs their callers."""
    found = [(share, symbol) for share, symbol in perf_top(data, "sym", 400) if COPYING.search(symbol)]
    if not found:
        return ["Copies and allocations: no samples in the C library's copy, set or allocation routines.", ""]
    lines = [f"Copies and allocations: {sum(share for share, _ in found):.1f}% of samples, "
             + ", ".join(f"`{symbol.split()[-1]}` {share:.1f}%" for share, symbol in found[:6]), ""]
    if traced:
        for _, symbol in found[:3]:
            name = symbol.split()[-1]
            output = subprocess.run(
                ["perf", "report", "-i", str(data), "--stdio", "--no-children", "--symbols", name,
                 "-g", "caller,0.5,callee,function,percent", "--percent-limit", "0.3"],
                capture_output=True, text=True,
            ).stdout
            callers = [line.rstrip() for line in output.splitlines() if line.strip() and not line.startswith("#")]
            lines += [f"Who calls `{name}`:", "", "```", *callers[:30], "```", ""]
    return lines


def perf(options) -> str:
    """Hardware counters and sampled hot spots of dinara-align on every sample; returns the report."""
    reason = perf_ready()
    if reason:
        return f"# Perf\n\nSkipped: {reason}.\n"
    binary = mojo_build(ROOT, HERE / "ours.mojo", BUILDS / f"runner-lines-{harness.CPU}", extra=["-g1"])
    tool = ours_tool(options)
    lines = [
        f"# Perf: where dinara-align's cycles go{', at ' + options.affine if options.affine else ''}",
        "",
        harness.build_note().strip(),
        "",
        "The whole runner on each sample, its start and reading included (small beside the pairs except on the"
        " shortest samples). Misses are per thousand instructions.",
        "",
        "| dataset | cycles | IPC | branch misses / branch | L1d misses / kinst | LLC misses / kinst |",
        "| :-- | --: | --: | --: | --: | --: |",
    ]
    hot = []
    for name, path, _ in samples(options):
        command = [str(binary), "seq", tool, str(options.budget), str(path)]
        print(f"{name}: perf stat ...", file=sys.stderr, flush=True)
        counts = perf_stat(command)
        instructions = counts.get("instructions", 0) or math.nan
        cycles = counts.get("cycles", 0) or math.nan

        def per_kilo(event: str) -> str:
            return f"{1000 * counts[event] / instructions:.2f}" if event in counts else "—"

        branch = (f"{100 * counts['branch-misses'] / counts['branches']:.2f}%"
                  if counts.get("branches") else "—")
        lines.append(f"| {name} | {cycles:.3g} | {instructions / cycles:.2f} | {branch} | "
                     f"{per_kilo('L1-dcache-load-misses')} | {per_kilo('LLC-load-misses')} |")
        print(f"{name}: perf record ...", file=sys.stderr, flush=True)
        data = BUILDS / f"perf-{name}.data"
        record = ["perf", "record", "-q", "-F", str(options.frequency), "-o", str(data)]
        # Callers through the last-branch record where the processor keeps one, which needs no frame pointers.
        traced = subprocess.run([*record, "--call-graph", "lbr", "--", *command], capture_output=True).returncode == 0
        if not traced:
            subprocess.run([*record, "--", *command], capture_output=True, check=True)
        engines = {}
        for share, file in perf_top(data, "srcfile", 200):
            engine = ENGINE_FILES.get(Path(file.split(":")[0]).name, Path(file).name or "other")
            engines[engine] = engines.get(engine, 0.0) + share
        hot += [
            f"## {name}",
            "",
            "By engine (source file): " + ", ".join(f"{engine} {share:.1f}%" for engine, share in
                                                     sorted(engines.items(), key=lambda item: -item[1])[:8]),
            "",
            "| share | function |",
            "| --: | :-- |",
            *[f"| {share:.1f}% | `{symbol[:110]}` |" for share, symbol in perf_top(data, "sym", options.top)],
            "",
            "| share | source line |",
            "| --: | :-- |",
            *[f"| {share:.1f}% | `{line}` |" for share, line in perf_top(data, "srcline", options.top)],
            "",
            *copies_and_allocations(data, traced),
        ]
    report = "\n".join(lines + [""] + hot) + "\n"
    stem = f"profile-perf-{options.affine.replace(',', '-')}" if options.affine else "profile-perf"
    write(stem, report, None)
    return report


# endregion perf

# region probes


def probes(options) -> str:
    """The engines and work of every pair, from a runner built with probes; returns the report."""
    import probes_report

    binary = mojo_build(ROOT, HERE / "probed.mojo", BUILDS / f"runner-probes-{harness.CPU}",
                        extra=["-D", "DINARA_PROBES=1"])
    tool = ours_tool(options)
    found = {}
    for name, path, facts in samples(options):
        print(f"{name}: probes ...", file=sys.stderr, flush=True)
        found[name] = (probes_report.run(binary, tool, path, options.budget), facts)
    report = probes_report.report(found, options)
    stem = f"profile-probes-{options.affine.replace(',', '-')}" if options.affine else "profile-probes"
    write(stem, report, {name: records for name, (records, _) in found.items()})
    return report


# endregion probes


def write(stem: str, report: str, data) -> None:
    """The report as `stem.md` and its per-pair figures as `stem.json` under the results directory."""
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f"{stem}.md").write_text(report)
    if data is not None:
        (RESULTS / f"{stem}.json").write_text(json.dumps(data))
    print(f"wrote {RESULTS / (stem + '.md')}", file=sys.stderr)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=["duel", "probes", "perf", "compare", "all"])
    parser.add_argument("base", nargs="?", help="compare: the commit to compare against")
    parser.add_argument("head", nargs="?", default=".", help="compare: the commit compared, or . (default) for this checkout")
    parser.add_argument("--only", nargs="*", default=[], help="only datasets whose names hold one of these")
    parser.add_argument("--affine", metavar="X,O,E", help="at affine costs: WFA's mismatch, opening and extension")
    parser.add_argument("--bases", type=float, default=bench.SAMPLE_BASES, help="sequence per sample, as bench.py")
    parser.add_argument("--budget", type=float, default=30.0, help="seconds per run per sample (default 30)")
    parser.add_argument("--rounds", type=int, default=3, help="rounds of dinara-align (default 3)")
    parser.add_argument("--fresh", action="store_true", help="duel: run the rivals again rather than reuse them")
    parser.add_argument("--kept-only", action="store_true", help="duel: only rivals with kept results; never run one")
    parser.add_argument("--worst", type=int, default=5, help="duel: the losing pairs listed per dataset (default 5)")
    parser.add_argument("--top", type=int, default=12, help="perf: functions and lines listed per dataset")
    parser.add_argument("--frequency", type=int, default=4999, help="perf: samples a second (default 4999)")
    parser.add_argument("--unmitigated", action="store_true", help="compare: plain builds, jumps not kept clear")
    parser.add_argument("--core", type=int, default=2, help="the core every run is pinned to on Linux (default 2)")
    parser.add_argument("--force", action="store_true", help="time even when the machine is not quiet")
    parser.add_argument("--cpu", default="native", help="the CPU every tool is built for (default: the host's)")
    options = parser.parse_args()
    if options.command == "compare" and not options.base:
        parser.error("compare needs the commit to compare against")
    harness.set_cpu(options.cpu)
    pin(options.core)
    check_quiet(options.core, options.force)
    commands = ["duel", "probes", "perf"] if options.command == "all" else [options.command]
    reports = [globals()[command](options) for command in commands]
    if options.command == "all":
        write("profile-all" + (f"-{options.affine.replace(',', '-')}" if options.affine else ""),
              "\n\n".join(reports), None)
    print("\n\n".join(reports))


if __name__ == "__main__":
    main()
