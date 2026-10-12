#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Pinpoints where dinara-align spends its time, where it could spend less, and where a rival beats it,
dataset by dataset and pair by pair, on the samples `bench.py` times.

    pixi run profile duel [--affine 4,6,2]   # every pair against every rival: where dinara-align loses
    pixi run profile rivals                  # the lost pairs, ours and the winner's hardware counts side by side
    pixi run profile regret --only ont-10k   # every pair under every engine policy: what each choice cost
    pixi run profile probes --only ont-10k   # per pair: engine, phases, work, counters, allocations, copies
    pixi run profile peak                    # each kernel's rate alone, which probes measure phases against
    pixi run profile perf [--target batch]   # Linux: counters, top-down, hot functions, lines, instructions
    pixi run profile modes                   # the probes of the other modes' workloads
    pixi run profile scaling                 # the batch over threads, false sharing, the GPU's breakdown
    pixi run profile noise                   # what moving the code alone does: the floor for comparisons
    pixi run profile compare main~1          # this checkout against a commit, pair by pair, with intervals
    pixi run profile all                     # every view but compare, one report

`duel` times dinara-align over `--rounds` rounds and sets each pair's median beside every rival's time for
the same pair, the rivals' kept results (`harness.measure`) or a fresh run of them. A dataset's mean can
win while single pairs lose: the report lists, for each dataset, the margin over the fastest rival, the
pairs some rival aligned faster, those it aligned faster than dinara-align's fastest round (beyond its
own noise), the time those losses cost, and the worst of them with their length, cost and divergence.
`rivals` then counts each such pair alone, for dinara-align and for the rival that won it (the pair run
once and many times, the difference per copy): more instructions mean more work, the same instructions at
a lower IPC slower work.

`regret` builds the runner once with the library's own choices and once with each policy forced
(`-D DINARA_FORCE`, see `dinara_align/probes.mojo`: band without diagonal transition, diagonal transition
to its limit, the whole matrix, each kind of seeds; at gap-affine costs the difference sweep from the
start, never, or both searches from the start), times every pair under each, holds each policy's costs to
the default's, and reports the time the default would save by taking, pair by pair, the policy that beat
it beyond noise: the cost of every choice the library makes, and the pairs and features where it chose
wrong.

`probes` builds `probed.mojo` with `-D DINARA_PROBES`, which records per pair the engine that settled it,
the time each phase took and the instructions and cycles it retired (`rdpmc` on Linux x86), and the work
it did: diagonal-transition steps, band rounds and the bounds they tried against the distance found,
tiles swept and recomputed, seeds built and pruned, wavefront costs and sweep cells. Around each pair it
reads the processor's counters (instructions, cycles, branch, cache and TLB misses, `perf_event_open`),
and with `probes/counting.c` preloaded the allocations and the bytes copied, moved and set through the C
library, and the page faults. A normal build compiles none of it. The report gives each dataset's time by
engine and phase, each phase's IPC and its rate against the kernel's peak (`peak`), the work against what
the pair needed, the processor's and the allocator's view per pair, and the pairs whose probes point at a
choice that cost time (see `probes_report`). `modes` does the same for the other modes' workloads.

`peak` runs each kernel alone, the engine forced, on the workload that keeps it busiest: the band's
bit-parallel sweep, diagonal transition, the gap-affine wavefront and the difference sweep, each a rate in
its own units, against which `probes` measures every phase doing that kernel's work.

`perf` (Linux, `perf` installed, `perf_event_paranoid` at most 2) runs the runner built with line tables,
or with `--target` the batch or the modes drivers, under `perf stat` and `perf record`: cycles, IPC, branch
and cache misses, Intel's top-down shares (retiring, bad speculation, front-end and back-end bound, and
level 2), the share of samples by function, by source file (each file one engine) and by source line,
inlined code charged to its own line, the hottest functions down to their instructions with the share on
vector instructions and on copies and sets the compiler inlined, and the samples in the C library's copy,
set and allocation routines, with their callers where the processor records its last branches (LBR).

`scaling` times the short-read batch on 1, 2, 4 … every thread, records the cache lines threads write
under each other (`perf c2c`), and breaks the GPU's call into kernels, copies and host work where Nsight
Systems is installed.

`noise` builds the same code with its functions moved by 0, 16, 32 and 48 bytes and times them alike: how
far layout alone moves each dataset's mean, the floor below which `compare` calls a difference noise.
`compare` builds the runner from two commits, or a commit and this checkout as it stands, times both over
alternating rounds, holds every pair's costs equal, and gives each dataset's ratio of mean times with a
95% bootstrap interval over pairs and rounds: faster or slower only when the interval excludes 1 and the
difference exceeds the layout floor. On x86 every build keeps its jumps clear of 32-byte boundaries.

Every run but `scaling` pins itself to one core on Linux (`--core`, default 2), and every run refuses to
start while that core, or the machine, is busy (`--force` to override), since a shared core corrupts every
timing silently. Reports go to `benchmarks/.cache/results/profile-*.md`, each pair's figures to a JSON
beside it.
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


def on_core(core: int) -> str:
    """The busiest processes last seen on `core`, for a refusal to name what keeps it busy."""
    if not shutil.which("ps"):
        return ""
    output = subprocess.run(["ps", "-eo", "psr,pcpu,comm", "--sort=-pcpu"], capture_output=True, text=True).stdout
    found = [" ".join(line.split()[1:]) for line in output.splitlines()[1:]
             if line.split() and line.split()[0] == str(core)]
    return ", ".join(found[:4])


def check_quiet(core: int, force: bool, patience: int = 60) -> None:
    """Waits, up to `patience` seconds, for `core` to stay idle through a one-second sample and for the
    machine's load to be low, and stops if it never does: two processes on a pinned core take turns, and
    every time measured then is wrong without any error. A busy moment, this script's own start among them,
    only delays the run. `force` warns instead of stopping."""
    cores = os.cpu_count() or 1
    problems = []
    for _ in range(patience):
        problems = []
        busy = core_busy(core)
        if busy is not None and busy > 0.05:
            problems.append(f"core {core} was {100 * busy:.0f}% busy over a second ({on_core(core)})")
        if hasattr(os, "getloadavg") and os.getloadavg()[0] > max(1.5, cores / 4):
            problems.append(f"load average {os.getloadavg()[0]:.2f}")
        if not problems or busy is None and not hasattr(os, "getloadavg"):
            return
    message = f"the machine stayed busy for {patience} s: " + "; ".join(problems)
    if not force:
        sys.exit(message + " (pass --force to time anyway)")
    print("WARNING: " + message, file=sys.stderr)


# endregion The machine

# region Builds


def mojo_build(root: Path, source: Path, binary: Path, extra: list = (), mitigated: bool = False, pad: int = 0) -> Path:
    """`source` built against the library at `root` into `binary`, kept while the command and every source are
    unchanged. With `mitigated` on x86, its jumps kept clear of 32-byte boundaries (see `ablation`), and with
    `pad`, its code moved `pad` bytes along, which `noise` uses to measure what layout alone does."""
    binary.parent.mkdir(parents=True, exist_ok=True)
    command = ["mojo", "build", *extra, "-I", str(root), str(source), "-o", str(binary), *mojo_cpu_args()]
    stamp = binary.with_name(binary.name + ".command")
    sources = [*source.parent.glob("*.mojo"), *(root / "dinara_align").rglob("*.mojo")]
    assembled = (mitigated or pad) and x86() and shutil.which("clang") and shutil.which("cc")
    signature = " ".join(command) + (f" assembled mitigated={mitigated} pad={pad}" if assembled else "")
    if (
        binary.exists()
        and stamp.exists()
        and stamp.read_text() == signature
        and all(item.stat().st_mtime < binary.stat().st_mtime for item in sources)
    ):
        return binary
    # A build times nothing, so it takes every core, and the pin the timings need comes back after it.
    pinned = os.sched_getaffinity(0) if hasattr(os, "sched_getaffinity") else None
    if pinned is not None:
        os.sched_setaffinity(0, set(range(os.cpu_count() or 1)))
    try:
        if assembled:
            assembled_build(root, source, binary, list(extra), mitigated, pad)
        else:
            subprocess.run(command, check=True)
    finally:
        if pinned is not None:
            os.sched_setaffinity(0, pinned)
    stamp.write_text(signature)
    return binary


def assembled_build(root: Path, source: Path, binary: Path, defines: list, mitigated: bool, pad: int) -> None:
    """Mojo's assembly for `source`, `pad` bytes of padding ahead of its code, assembled by clang, with
    `-mbranches-within-32B-boundaries` when `mitigated`, and linked as `mojo build` links (see
    `ablation.mitigated_build`)."""
    assembly = binary.with_name(binary.name + ".s")
    subprocess.run(["mojo", "build", *defines, "-I", str(root), "--emit", "asm", str(source), "-o", str(assembly),
                    *mojo_cpu_args()], check=True)
    # Mojo's LLVM writes `.prefalign`, a preferred alignment older assemblers predate: the same power as `.p2align`.
    text = re.sub(r"^\s*\.prefalign\s+(\d+),.*$", r"\t.p2align\t\1, 0x90", assembly.read_text(), flags=re.M)
    if pad:
        # A block of no-ops ahead of the first function: every function after it moves `pad` bytes along,
        # modulo the alignment each asks for.
        text = f"\t.text\n\t.p2align\t6, 0x90\n\t.space\t{pad}, 0x90\n" + text
    assembly.write_text(text)
    objects = binary.with_name(binary.name + ".o")
    flags = ["-mbranches-within-32B-boundaries"] if mitigated else []
    subprocess.run(["clang", "-c", *flags, str(assembly), "-o", str(objects)], check=True)
    library = Path(os.environ["CONDA_PREFIX"]) / "lib"
    subprocess.run(["cc", str(objects), str(library / "libKGENCompilerRTShared.so"), "-Xlinker", "-rpath", "-Xlinker",
                    str(library), str(library / "libAsyncRTMojoBindings.so"), "-o", str(binary), "-Wl,--gc-sections",
                    "-lm"], check=True)


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
    noise_file = RESULTS / "profile-noise.json"
    floors = json.loads(noise_file.read_text()) if noise_file.exists() else {}
    lines = [
        f"# Compare: {head_label} against {base_label}{', at ' + options.affine if options.affine else ''}",
        "",
        harness.build_note().strip(),
        "",
        f"Mean time per pair over {options.rounds} alternating rounds, and the ratio {head_label} / {base_label} with"
        " its 95% bootstrap interval over pairs and rounds. Every pair's costs are equal in both."
        + (" A difference no larger than what moving the code alone does (`profile noise`) is within layout noise."
           if floors else " Run `profile noise` once to also mark differences that layout alone could make."),
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
        # A difference the code's layout alone could make (see `noise`) is not the code's.
        floor = floors.get(name)
        if floor is not None and verdict != "within noise" and abs(ratio - 1) <= floor:
            verdict = f"within layout noise ({100 * floor:.1f}%)"
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

# region regret

UNIT_POLICIES = ["", "band", "diagonal", "full", "seeds-none", "seeds-exact", "seeds-inexact"]
"""The policies `regret` times at unit costs, "" the library's own choices (see `dinara_align.probes.FORCE`)."""
AFFINE_POLICIES = ["", "sweep", "no-sweep", "bidirectional"]
"""The policies it times at gap-affine costs."""


def policy_runner(policy: str, options) -> Path:
    """The plain runner built with `policy` forced, or with none."""
    extra = ["-D", f"DINARA_FORCE={policy}"] if policy else []
    return mojo_build(ROOT, HERE / "ours.mojo", BUILDS / f"runner-force-{policy or 'default'}-{harness.CPU}", extra,
                      mitigated=not options.unmitigated)


def probe_engines(options) -> dict:
    """Each pair's engine by dataset, from the last `probes` report at the same costs, if any."""
    stem = f"profile-probes-{options.affine.replace(',', '-')}" if options.affine else "profile-probes"
    path = RESULTS / f"{stem}.json"
    if not path.exists():
        return {}
    return {name: [record.get("engine", "") for record in records] for name, records in json.loads(path.read_text()).items()}


def regret(options) -> str:
    """Every pair under every policy the library could have followed; returns the report."""
    policies = AFFINE_POLICIES if options.affine else UNIT_POLICIES
    binaries = {policy: policy_runner(policy, options) for policy in policies}
    tool = ours_tool(options)
    engines = probe_engines(options)
    label = {policy: policy or "default" for policy in policies}
    lines = [
        f"# Regret: what each choice of engine cost{', at ' + options.affine if options.affine else ''}",
        "",
        harness.build_note().strip(),
        "",
        f"Every pair under the library's own choices and under each policy forced in turn ("
        + ", ".join(f"`{policy}`" for policy in policies if policy)
        + f"; see `dinara_align/probes.mojo`), {options.rounds} alternating rounds, every policy's costs held to "
        "the default's. A policy *wins* a pair when its slowest round beat the default's fastest. *Regret*: the "
        "time the default would save by taking, on each pair, the policy that won it, as a share of its total; "
        "*oracle*: that total over the default's. A policy stopped by its budget is left out of the pairs it "
        "did not reach.",
        "",
        "| dataset | pairs | default | regret | oracle | " + " | ".join(f"{label[p]} wins (saves)" for p in policies if p)
        + " |",
        "| :-- | --: | --: | --: | --: |" + " --: |" * (len(policies) - 1),
    ]
    worst_lines, everything = [], {}
    for name, path, facts in samples(options):
        runs = {policy: [] for policy in policies}
        for round in range(options.rounds):
            order = policies[round % len(policies):] + policies[: round % len(policies)]
            for policy in order:
                print(f"{name}: {label[policy]}, round {round + 1} of {options.rounds} ...", file=sys.stderr, flush=True)
                rows, _, status = harness.run_tool(binaries[policy], tool, path, options.budget)
                if status == harness.FAILED:
                    sys.exit(f"{label[policy]} failed on {name}")
                runs[policy].append(rows)
        timed = {}
        for policy, rounds in runs.items():
            finished = min(len(rows) for rows in rounds)
            timed[policy] = [
                {"median": statistics.median(rows[index][0] for rows in rounds),
                 "fastest": min(rows[index][0] for rows in rounds),
                 "slowest": max(rows[index][0] for rows in rounds),
                 "cost": rounds[0][index][1]}
                for index in range(finished)
            ]
        default = timed[""]
        for policy, pairs in timed.items():
            for index, pair in enumerate(pairs[: len(default)]):
                if pair["cost"] != default[index]["cost"]:
                    sys.exit(f"{name}: policy {label[policy]} found cost {pair['cost']} on pair {index}, the default "
                             f"{default[index]['cost']}: a forced policy must stay exact")
        if not default:
            lines.append(f"| {name} | 0/{facts['pairs']} | t.o. | | |" + " |" * (len(policies) - 1))
            continue
        total = sum(pair["median"] for pair in default)
        oracle, regret_total = 0.0, 0.0
        wins = {policy: [0, 0.0] for policy in policies if policy}
        per_pair = []
        for index, pair in enumerate(default):
            best, best_policy = pair["median"], ""
            for policy in policies[1:]:
                if index < len(timed[policy]):
                    other = timed[policy][index]
                    if other["slowest"] < pair["fastest"]:
                        wins[policy][0] += 1
                        wins[policy][1] += pair["median"] - other["median"]
                        if other["median"] < best:
                            best, best_policy = other["median"], policy
            oracle += best
            regret_total += pair["median"] - best
            per_pair.append({"pair": index, "default": pair["median"], "best": best, "policy": best_policy,
                             "cost": pair["cost"]})
        pairs = f"{len(default)}" if len(default) == facts["pairs"] else f"{len(default)}/{facts['pairs']}"
        cells = [f"{count} ({share_of(saved, total)})" for count, saved in wins.values()]
        lines.append(f"| {name} | {pairs} | {duration(total / len(default))} | {share_of(regret_total, total)} | "
                     f"{oracle / total:.3f} | " + " | ".join(cells) + " |")
        for item in sorted(per_pair, key=lambda item: item["best"] - item["default"])[: options.worst]:
            if not item["policy"]:
                continue
            index = item["pair"]
            length = facts["lengths"][index] if index < len(facts["lengths"]) else 0
            engine = engines.get(name, [])[index] if index < len(engines.get(name, [])) else ""
            worst_lines.append(
                f"| {name} | {index} | {length:,} | {item['cost']:,} | {100 * item['cost'] / max(length, 1):.1f}% | "
                f"{engine or '—'} | {duration(item['default'])} | {label[item['policy']]} {duration(item['best'])} | "
                f"{item['default'] / item['best']:.2f}× |"
            )
        everything[name] = {"timed": {label[p]: v for p, v in timed.items()}, "lengths": facts["lengths"]}
    lines += [
        "",
        f"## The pairs the default's choice cost most (at most {options.worst} a dataset)",
        "",
        "Engine: the one the default took, from the last `probes` report, where there is one.",
        "",
        "| dataset | pair | length | cost | divergence | engine | default | best policy | ratio |",
        "| :-- | --: | --: | --: | --: | :-- | --: | --: | --: |",
        *worst_lines,
    ]
    stem = f"profile-regret-{options.affine.replace(',', '-')}" if options.affine else "profile-regret"
    write(stem, "\n".join(lines) + "\n", everything)
    return "\n".join(lines) + "\n"


def share_of(part: float, whole: float) -> str:
    """`part` as a percentage of `whole`."""
    return f"{100 * part / whole:.1f}%" if whole else "—"


# endregion regret

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


def perf_stat(command: list, delay: int = 0) -> dict:
    """`perf stat`'s counts for `command`, by event, the first `delay` milliseconds left out."""
    output = subprocess.run(["perf", "stat", "-x", ",", "-D", str(delay), "-e", ",".join(COUNTERS), "--", *command],
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


TOPDOWN = re.compile(r"([\d.]+)\s*%\s+(tma_\w+)")
"""A top-down metric as `perf stat -M Topdown...` prints it: its share and its name."""
ANNOTATED = re.compile(r"^\s*([\d.]+)\s*:\s*([0-9a-f]+):\s*(\S.*)$")
"""A line of `perf annotate --stdio`: the share of the symbol's samples, the address and the instruction."""
BULK = re.compile(r"\brep(z|nz|e|ne)?\s+(movs|stos)|\b(movs|stos)[bwlq]\b")
"""String instructions, which a compiler emits for copies and sets it inlines rather than call."""
VECTOR = re.compile(r"%[xyz]mm\d|%k[0-7]\b")
"""An instruction on vector or mask registers."""


def topdown(command: list, delay: int = 0) -> dict:
    """Intel's top-down shares for `command`, level 1 and level 2, by metric name; empty where the processor or
    perf has none."""
    found = {}
    for group in ("TopdownL1", "TopdownL2"):
        output = subprocess.run(["perf", "stat", "-D", str(delay), "-M", group, "--", *command], capture_output=True,
                                text=True)
        for share, name in TOPDOWN.findall(output.stderr):
            found[name] = float(share)
    return found


def annotate(data: Path, symbol: str, top: int) -> tuple[list[str], float, float, float]:
    """The hottest instructions of `symbol` in a `perf record` file, the shares of its samples on string copies
    and sets and on vector instructions, and the share annotated."""
    output = subprocess.run(["perf", "annotate", "-i", str(data), "--stdio", "--no-source", "-s", symbol],
                            capture_output=True, text=True).stdout
    instructions = []
    for line in output.splitlines():
        match = ANNOTATED.match(line)
        if match:
            instructions.append((float(match.group(1)), match.group(3).strip()))
    bulk = sum(share for share, text in instructions if BULK.search(text))
    vector = sum(share for share, text in instructions if VECTOR.search(text))
    hottest = [f"| {share:.1f}% | `{text[:100]}` |" for share, text in sorted(instructions, key=lambda item: -item[0])[:top]
               if share >= 0.5]
    return hottest, bulk, vector, sum(share for share, _ in instructions)


START_DELAY_MS = 100
"""What `perf` leaves out at the start of a runner: its 50 ms spin onto a fast core (see `ours.mojo`) and its
reading, which a pair's time leaves out too."""


def perf_targets(options) -> list[tuple[str, list]]:
    """What `perf` runs: each sample through the runner, or with `--target batch` each short-read workload on
    one thread, or with `--target modes` each modes workload; every binary built with line tables."""
    if options.target == "batch":
        import batch_bench

        batch_bench.generate()
        binary = mojo_build(ROOT, HERE / "batches.mojo", BUILDS / f"batches-lines-{harness.CPU}", extra=["-g1"])
        pairs = harness.DATA / "batches" / "illumina.tsv"
        return [(workload, [str(binary), str(pairs), workload, "1"]) for workload in batch_bench.WORKLOADS]
    if options.target == "modes":
        import modes_bench

        modes_bench.generate()
        binary = mojo_build(ROOT, HERE / "modes.mojo", BUILDS / f"modes-lines-{harness.CPU}", extra=["-g1"])
        chosen = [workload for workload in modes_bench.WORKLOADS
                  if not options.only or any(part in workload for part in options.only)]
        return [(workload, [str(binary), str(modes_bench.MODES / f"{workload}.tsv")]) for workload in chosen]
    binary = mojo_build(ROOT, HERE / "ours.mojo", BUILDS / f"runner-lines-{harness.CPU}", extra=["-g1"])
    tool = ours_tool(options)
    return [(name, [str(binary), "seq", tool, str(options.budget), str(path)]) for name, path, _ in samples(options)]


def perf(options) -> str:
    """Hardware counters, top-down shares and sampled hot spots, down to instructions; returns the report."""
    reason = perf_ready()
    if reason:
        return f"# Perf\n\nSkipped: {reason}.\n"
    targets = perf_targets(options)
    what = {"pairs": "each sample through the runner", "batch": "each short-read workload on one thread",
            "modes": "each modes workload"}[options.target]
    lines = [
        f"# Perf: where dinara-align's cycles go{', at ' + options.affine if options.affine else ''}",
        "",
        harness.build_note().strip(),
        "",
        f"The whole program on {what}, its start and reading included (small beside the work except on the"
        " shortest). Misses are per thousand instructions. Top-down: Intel's shares of the pipeline's issue"
        " slots, retiring (useful work), bad speculation, front-end bound (instructions not delivered) and back-end"
        " bound (waiting on memory or on execution units), and level 2 beneath them.",
        "",
        "| target | cycles | IPC | branch misses / branch | L1d misses / kinst | LLC misses / kinst | retiring | bad"
        " speculation | front-end bound | back-end bound |",
        "| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: |",
    ]
    hot = []
    for name, command in targets:
        print(f"{name}: perf stat ...", file=sys.stderr, flush=True)
        # A runner's spin and reading left out; the batch and modes drivers time after reading too, but their
        # reading is long and varies, so they are taken whole.
        delay = START_DELAY_MS if options.target == "pairs" else 0
        counts = perf_stat(command, delay)
        shares = topdown(command, delay)
        instructions = counts.get("instructions", 0) or math.nan
        cycles = counts.get("cycles", 0) or math.nan

        def per_kilo(event: str) -> str:
            return f"{1000 * counts[event] / instructions:.2f}" if event in counts else "—"

        def level(name: str) -> str:
            return f"{shares[name]:.1f}%" if name in shares else "—"

        branch = (f"{100 * counts['branch-misses'] / counts['branches']:.2f}%"
                  if counts.get("branches") else "—")
        lines.append(f"| {name} | {cycles:.3g} | {instructions / cycles:.2f} | {branch} | "
                     f"{per_kilo('L1-dcache-load-misses')} | {per_kilo('LLC-load-misses')} | {level('tma_retiring')} | "
                     f"{level('tma_bad_speculation')} | {level('tma_frontend_bound')} | {level('tma_backend_bound')} |")
        print(f"{name}: perf record ...", file=sys.stderr, flush=True)
        data = BUILDS / f"perf-{options.target}-{name}.data"
        record = ["perf", "record", "-q", "-D", str(delay), "-F", str(options.frequency), "-o", str(data)]
        # Callers through the last-branch record where the processor keeps one, which needs no frame pointers.
        traced = subprocess.run([*record, "--call-graph", "lbr", "--", *command], capture_output=True).returncode == 0
        if not traced:
            subprocess.run([*record, "--", *command], capture_output=True, check=True)
        engines = {}
        for share, file in perf_top(data, "srcfile", 200):
            engine = ENGINE_FILES.get(Path(file.split(":")[0]).name, Path(file).name or "other")
            engines[engine] = engines.get(engine, 0.0) + share
        level_two = sorted(((share, metric) for metric, share in shares.items()
                            if metric not in ("tma_retiring", "tma_bad_speculation", "tma_frontend_bound",
                                              "tma_backend_bound")), reverse=True)
        symbols = perf_top(data, "sym", options.top)
        hot += [
            f"## {name}",
            "",
            "By engine (source file): " + ", ".join(f"{engine} {share:.1f}%" for engine, share in
                                                     sorted(engines.items(), key=lambda item: -item[1])[:8]),
            "",
        ]
        if level_two:
            hot += ["Top-down, level 2: " + ", ".join(f"{metric.removeprefix('tma_')} {share:.1f}%"
                                                     for share, metric in level_two[:8]), ""]
        hot += [
            "| share | function |",
            "| --: | :-- |",
            *[f"| {share:.1f}% | `{symbol[:110]}` |" for share, symbol in symbols],
            "",
            "| share | source line |",
            "| --: | :-- |",
            *[f"| {share:.1f}% | `{line}` |" for share, line in perf_top(data, "srcline", options.top)],
            "",
            *copies_and_allocations(data, traced),
        ]
        # The hottest functions down to their instructions: where inside them the samples land, how much on
        # copies and sets the compiler inlined, and how much on vector instructions.
        for share, symbol in symbols[: options.annotate]:
            name_only = symbol.split()[-1]
            hottest, bulk, vector, covered = annotate(data, name_only, options.top)
            if not hottest:
                continue
            hot += [
                f"`{name_only[:100]}`, {share:.1f}% of samples: {vector:.0f}% of its samples on vector instructions, "
                f"{bulk:.1f}% on inlined copies and sets (`rep movs`, `rep stos`).",
                "",
                "| share of the function | instruction |",
                "| --: | :-- |",
                *hottest,
                "",
            ]
    report = "\n".join(lines + [""] + hot) + "\n"
    suffix = "" if options.target == "pairs" else f"-{options.target}"
    stem = (f"profile-perf{suffix}-{options.affine.replace(',', '-')}" if options.affine else f"profile-perf{suffix}")
    write(stem, report, None)
    return report


# endregion perf

# region peak


def synthetic(name: str, pairs: int, length: int, rate: float, unrelated: bool = False) -> Path:
    """A `.seq` file of `pairs` random pairs of `length` bases, the second the first at `rate` errors spread
    uniformly, or unrelated; made once from a fixed seed."""
    path = BUILDS / "peak" / f"{name}.seq"
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    generator = random.Random(20261012)
    lines = []
    for _ in range(pairs):
        first = "".join(generator.choices("ACGT", k=length))
        second = ("".join(generator.choices("ACGT", k=length)) if unrelated
                  else harness.mutate(first, rate, "ACGT", generator))
        lines += [f">{first}", f"<{second}"]
    path.write_text("\n".join(lines) + "\n")
    return path


PEAKS = [
    # name, policy, tool kind, workload, the probe counting the work, the probe timing it, the unit
    ("band_cells_per_ns", "full", "distance", ("sweep-30k", 6, 30_000, 0.3, False), "full_word_columns", "full_ns", 64),
    ("diagonal_steps_per_ns", "diagonal", "distance", ("near-100k", 6, 100_000, 0.004, False), "two_ended_steps",
     "two_ended_ns", 1),
    ("wavefront_diagonals_per_ns", "no-sweep", "affine", ("wave-10k", 6, 10_000, 0.05, False), "wavefront_diagonals",
     "wavefront_ns", 1),
    ("sweep_cells_per_ns", "sweep", "affine", ("unrelated-5k", 6, 5_000, 0.0, True), "sweep_cells", "sweep_ns", 1),
]
"""Each kernel's peak: the policy that makes it alone settle a pair, the workload that keeps it busiest, and
the probes whose ratio is its rate."""


def peak(options) -> str:
    """Each kernel's rate in isolation, the most of three runs, for `probes` to measure each phase against;
    returns the report."""
    import probes_report

    costs = options.affine or "4,6,2"
    rates = {}
    lines = [
        "# Peak: each kernel's rate alone",
        "",
        harness.build_note().strip(),
        "",
        "Each kernel alone on the workload that keeps it busiest, the engine forced (see"
        " `dinara_align/probes.mojo`), the most of three runs: the rate `probes` measures each phase against.",
        "",
        "| kernel | workload | rate |",
        "| :-- | :-- | --: |",
    ]
    for name, policy, kind, (workload, pairs, length, rate, unrelated), work, spent, unit in PEAKS:
        binary = mojo_build(ROOT, HERE / "probed.mojo", BUILDS / f"runner-probes-{policy}-{harness.CPU}",
                            extra=["-D", "DINARA_PROBES=1", "-D", f"DINARA_FORCE={policy}"])
        tool = f"dinara-align:{costs}" if kind == "affine" else "dinara-align (bit-parallel, 1 thread) distance"
        path = synthetic(workload, pairs, length, rate, unrelated)
        best = 0.0
        for _ in range(3):
            records = probes_report.run(binary, tool, path, options.budget)
            done = sum(record.get(work, 0) for record in records[1:]) * unit
            taken = sum(record.get(spent, 0) for record in records[1:])
            if taken:
                best = max(best, done / taken)
        rates[name] = best
        lines.append(f"| {name.replace('_', ' ')} | {pairs} pairs of {length:,} bp"
                     f"{', unrelated' if unrelated else f' at {100 * rate:g}%'} | {best:.3g} |")
    rates["machine"] = harness.processor()
    rates["costs"] = costs
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "profile-peak.json").write_text(json.dumps(rates, indent=1))
    report = "\n".join(lines) + "\n"
    write("profile-peak", report, None)
    return report


# endregion peak

# region rivals


def stat_per_pair(command_of, pair: tuple[str, str], copies: int) -> dict:
    """The hardware counts of one pair, alone: `command_of(path)` run on the pair once and `copies` times, the
    difference over `copies - 1`, so the program's start and reading cancel."""
    folder = BUILDS / "rivals"
    folder.mkdir(parents=True, exist_ok=True)
    one, many = folder / "one.seq", folder / "many.seq"
    one.write_text(f"{pair[0]}\n{pair[1]}\n")
    many.write_text(f"{pair[0]}\n{pair[1]}\n" * copies)
    once = perf_stat(command_of(one))
    repeated = perf_stat(command_of(many))
    return {event: (repeated[event] - once.get(event, 0)) / (copies - 1) for event in repeated}


def rivals(options) -> str:
    """dinara-align and the rival that beat it, side by side in hardware counts, on the pairs `duel` found it
    lost; returns the report."""
    reason = perf_ready()
    if reason:
        return f"# Rivals\n\nSkipped: {reason}.\n"
    stem = f"profile-duel-{options.affine.replace(',', '-')}" if options.affine else "profile-duel"
    found = RESULTS / f"{stem}.json"
    if not found.exists():
        return f"# Rivals\n\nSkipped: no `{stem}.json`; run `profile duel` first.\n"
    duels = json.loads(found.read_text())
    runners = rival_runners(options)
    ours = ours_tool(options)
    lines = [
        f"# Rivals: why a rival wins the pairs it wins{', at ' + options.affine if options.affine else ''}",
        "",
        harness.build_note().strip(),
        "",
        "Each pair `duel` found dinara-align losing beyond its noise, counted alone for dinara-align and for the"
        " rival that won it (the pair run once and many times, the difference per copy, so start-up cancels)."
        " More instructions mean more work done; the same instructions at a lower IPC mean the work runs slower"
        " (misses, dependencies, mispredictions).",
        "",
        "| dataset | pair | rival | instructions, ours / rival | IPC, ours / rival | branch misses / kinst, ours / rival"
        " | L1d misses / kinst, ours / rival | LLC misses / kinst, ours / rival | reading |",
        "| :-- | --: | :-- | --: | --: | --: | --: | --: | :-- |",
    ]
    for name, path, _ in samples(options):
        if name not in duels:
            continue
        entry = duels[name]
        ours_pairs, rival_rows = entry["ours"], entry["rivals"]
        losses = []
        for index, pair in enumerate(ours_pairs):
            best = min(((rows[index][0], column) for column, rows in rival_rows.items() if index < len(rows)),
                       default=None)
            if best and best[0] < pair["fastest"]:
                losses.append((pair["median"] - best[0], index, best[1]))
        text = path.read_text().splitlines()
        tools = {column: (binary, tool) for column, binary, tool in rival_tools(name, options, runners)}
        for _, index, column in sorted(losses, reverse=True)[: options.worst]:
            pair = (text[2 * index], text[2 * index + 1])
            copies = int(min(max(0.3 / max(ours_pairs[index]["median"], 1e-7), 3), 2000))
            print(f"{name}: pair {index} against {column} ...", file=sys.stderr, flush=True)
            mine = stat_per_pair(lambda file: [str(runners["ours"]), "seq", ours, "1e9", str(file)], pair, copies)
            binary, tool = tools[column]
            theirs = stat_per_pair(lambda file: [str(binary), "seq", tool, "1e9", str(file)], pair, copies)

            def ratio(event: str) -> str:
                return f"{mine.get(event, 0):.3g} / {theirs.get(event, 0):.3g}"

            def ipc(counts: dict) -> float:
                return counts.get("instructions", 0) / max(counts.get("cycles", 1), 1)

            def kilo(counts: dict, event: str) -> float:
                return 1000 * counts.get(event, 0) / max(counts.get("instructions", 1), 1)

            more = mine.get("instructions", 0) / max(theirs.get("instructions", 1), 1)
            reading = ("more work" if more > 1.2 else "slower work" if ipc(mine) < 0.83 * ipc(theirs) else "both close")
            lines.append(
                f"| {name} | {index} | {column} | {ratio('instructions')} ({more:.2f}×) | {ipc(mine):.2f} / "
                f"{ipc(theirs):.2f} | {kilo(mine, 'branch-misses'):.2f} / {kilo(theirs, 'branch-misses'):.2f} | "
                f"{kilo(mine, 'L1-dcache-load-misses'):.2f} / {kilo(theirs, 'L1-dcache-load-misses'):.2f} | "
                f"{kilo(mine, 'LLC-load-misses'):.2f} / {kilo(theirs, 'LLC-load-misses'):.2f} | {reading} |"
            )
    report = "\n".join(lines) + "\n"
    write(f"profile-rivals{'-' + options.affine.replace(',', '-') if options.affine else ''}", report, None)
    return report


# endregion rivals

# region scaling


def scaling(options) -> str:
    """The short-read batch on one thread and on more, the false sharing between threads, and the GPU's share;
    returns the report."""
    import batch_bench

    batch_bench.generate()
    binary = mojo_build(ROOT, HERE / "batches.mojo", BUILDS / f"batches-{harness.CPU}")
    pairs = harness.DATA / "batches" / "illumina.tsv"
    # Every core this process may use, not the one the timings are pinned to.
    cores = sorted(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else []
    if hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, set(range(os.cpu_count() or 1)))
    most = os.cpu_count() or 1
    counts = [count for count in (1, 2, 4, 8, 16, 32, 64) if count < most] + [most]
    lines = [
        "# Scaling: the short-read batch over threads",
        "",
        harness.build_note().strip(),
        "",
        "Each workload's 500,000 pairs on every thread count, the faster of the driver's two passes."
        " Efficiency: one thread's time over `threads` times this one's; 100% is perfect scaling.",
        "",
        "| workload | threads | time | speed-up | efficiency |",
        "| :-- | --: | --: | --: | --: |",
    ]
    gpu_lines = []
    for workload in batch_bench.WORKLOADS:
        single = None
        for count in counts:
            print(f"{workload}: {count} threads ...", file=sys.stderr, flush=True)
            output = subprocess.run([str(binary), str(pairs), workload, str(count)], capture_output=True, text=True,
                                    check=True).stdout
            for line in output.splitlines():
                fields = line.split("\t")
                if fields[0] == "dinara-align":
                    seconds = float(fields[2])
                    single = single or seconds
                    lines.append(f"| {workload} | {count} | {duration(seconds)} | {single / seconds:.2f}× | "
                                 f"{100 * single / (seconds * count):.0f}% |")
                elif fields[0] == "dinara-align GPU" and count == most:
                    gpu_lines.append(f"| {workload} | {duration(float(fields[2]))} |")
    if gpu_lines:
        lines += ["", "On the GPU, the whole call, packing included:", "", "| workload | time |", "| :-- | --: |",
                  *gpu_lines]
    lines += ["", *false_sharing(binary, pairs, most), *gpu_breakdown(binary, pairs, most)]
    if cores and hasattr(os, "sched_setaffinity"):
        os.sched_setaffinity(0, set(cores))
    report = "\n".join(lines) + "\n"
    write("profile-scaling", report, None)
    return report


def false_sharing(binary: Path, pairs: Path, threads: int) -> list[str]:
    """`perf c2c`'s cache lines that threads wrote under each other (HITM), on the affine batch on every
    thread: false or true sharing worth a look. Empty where `perf c2c` cannot run."""
    if perf_ready() or threads < 2:
        return []
    data = BUILDS / "c2c.data"
    recorded = subprocess.run(["perf", "c2c", "record", "-o", str(data), "--", str(binary), str(pairs),
                               "illumina-affine", str(threads)], capture_output=True, text=True)
    if recorded.returncode != 0:
        return [f"False sharing: `perf c2c` could not record ({recorded.stderr.strip().splitlines()[-1:]}).", ""]
    output = subprocess.run(["perf", "c2c", "report", "-i", str(data), "--stdio", "--stats"], capture_output=True,
                            text=True).stdout
    keep = [line for line in output.splitlines() if "HITM" in line or "Load" in line and ":" in line][:20]
    shared = subprocess.run(["perf", "c2c", "report", "-i", str(data), "--stdio", "-N", "--coalesce", "iaddr,dso,sym",
                             "--full-symbols"], capture_output=True, text=True).stdout
    table = [line for line in shared.splitlines() if line.strip() and not line.startswith("#")][:30]
    return ["False sharing (`perf c2c`, the affine batch on every thread):", "", "```", *keep, "", *table, "```", ""]


def gpu_breakdown(binary: Path, pairs: Path, threads: int) -> list[str]:
    """Where the GPU call's time goes, kernels against copies and the host, from Nsight Systems where it is
    installed; empty elsewhere."""
    if not shutil.which("nsys"):
        return ["GPU breakdown: Nsight Systems (`nsys`) is not installed; the GPU call is timed whole above.", ""]
    report = BUILDS / "gpu"
    subprocess.run(["nsys", "profile", "-f", "true", "-o", str(report), "--stats=false", str(binary), str(pairs),
                    "illumina-affine", str(threads)], capture_output=True)
    lines = ["GPU breakdown (Nsight Systems, the affine batch):", ""]
    for name in ("cuda_gpu_kern_sum", "cuda_gpu_mem_time_sum", "cuda_api_sum"):
        output = subprocess.run(["nsys", "stats", "-r", name, "-f", "csv", f"{report}.nsys-rep"], capture_output=True,
                                text=True).stdout
        rows = [line for line in output.splitlines() if line and not line.startswith("Processing")][:8]
        lines += [f"{name}:", "", "```", *rows, "```", ""]
    return lines


# endregion scaling

# region modes


def modes(options) -> str:
    """The probes of every pair of the modes workloads; returns the report."""
    import modes_bench
    import probes_report

    modes_bench.generate()
    binary = mojo_build(ROOT, HERE / "probed_modes.mojo", BUILDS / f"modes-probes-{harness.CPU}",
                        extra=["-D", "DINARA_PROBES=1"])
    environment = dict(os.environ)
    preload = probes_report.shim()
    if preload is not None:
        environment["LD_PRELOAD"] = str(preload)
    found = {}
    for workload in modes_bench.WORKLOADS:
        if options.only and not any(part in workload for part in options.only):
            continue
        print(f"{workload}: probes ...", file=sys.stderr, flush=True)
        output = subprocess.run([str(binary), str(modes_bench.MODES / f"{workload}.tsv")], capture_output=True,
                                text=True, env=environment, check=True).stdout
        records = []
        for line in output.splitlines():
            fields = line.split("\t")
            if len(fields) != 6:
                continue
            record = json.loads(fields[5])
            record.update(seconds=float(fields[1]), cost=int(fields[2]), columns=int(fields[3]), rows=int(fields[4]))
            records.append(record)
        found[workload] = (records, {})
    report = probes_report.report(found, options).replace("# Probes: what each pair took and did",
                                                          "# Modes: what each pair took and did", 1)
    write("profile-modes", report, {name: records for name, (records, _) in found.items()})
    return report


# endregion modes

# region noise

PADS = [0, 16, 32, 48]
"""The bytes `noise` moves the code by: every offset of a 64-byte line in steps of 16, the alignment the
compiler gives a function."""


def noise(options) -> str:
    """How much a rebuild that moves the code, and nothing else, moves the times: the floor below which a
    difference between two builds says nothing; returns the report."""
    assembled = x86() and shutil.which("clang") and shutil.which("cc")
    pads = PADS if assembled else [0]
    binaries = {pad: mojo_build(ROOT, HERE / "ours.mojo", BUILDS / f"runner-pad{pad}-{harness.CPU}", pad=pad,
                                mitigated=not options.unmitigated)
                for pad in pads}
    tool = ours_tool(options)
    floors = {}
    lines = [
        "# Noise: what layout alone does to the times",
        "",
        harness.build_note().strip(),
        "",
        ("The same code built with its functions moved by " + ", ".join(f"{pad}" for pad in pads)
         + " bytes, alternating rounds. *Layout*: the slowest build's mean over the fastest's, less 1, the floor"
         " below which two builds' difference is layout, not code; `compare` reports a difference below it as"
         " within noise. *Rounds*: the median spread of one build's rounds on a pair."
         if assembled else "Moving the code needs clang and x86; only the rounds' spread is measured here."),
        "",
        "| dataset | pairs | layout | rounds |",
        "| :-- | --: | --: | --: |",
    ]
    for name, path, facts in samples(options):
        runs = {pad: [] for pad in pads}
        for round in range(options.rounds):
            order = pads[round % len(pads):] + pads[: round % len(pads)]
            for pad in order:
                print(f"{name}: pad {pad}, round {round + 1} of {options.rounds} ...", file=sys.stderr, flush=True)
                rows, _, _ = harness.run_tool(binaries[pad], tool, path, options.budget)
                runs[pad].append(rows)
        finished = min(len(rows) for side in runs.values() for rows in side)
        if finished == 0:
            continue
        means = [sum(statistics.median(rows[index][0] for rows in runs[pad]) for index in range(finished))
                 for pad in pads]
        spreads = [(max(rows[index][0] for rows in runs[pad]) - min(rows[index][0] for rows in runs[pad]))
                   / max(statistics.median(rows[index][0] for rows in runs[pad]), 1e-12)
                   for pad in pads for index in range(finished)]
        layout = max(means) / min(means) - 1
        floors[name] = layout
        lines.append(f"| {name} | {finished} | {100 * layout:.1f}% | {100 * statistics.median(spreads):.1f}% |")
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / "profile-noise.json").write_text(json.dumps(floors, indent=1))
    report = "\n".join(lines) + "\n"
    write("profile-noise", report, None)
    return report


# endregion noise

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
    parser.add_argument("command", choices=["duel", "regret", "probes", "peak", "perf", "rivals", "modes", "scaling",
                                            "noise", "compare", "all"])
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
    parser.add_argument("--target", choices=["pairs", "batch", "modes"], default="pairs",
                        help="perf: the samples one pair at a time (default), the short-read batch, or the modes")
    parser.add_argument("--annotate", type=int, default=3, help="perf: hottest functions shown by instruction")
    parser.add_argument("--unmitigated", action="store_true", help="compare: plain builds, jumps not kept clear")
    parser.add_argument("--core", type=int, default=2, help="the core every run is pinned to on Linux (default 2)")
    parser.add_argument("--force", action="store_true", help="time even when the machine is not quiet")
    parser.add_argument("--cpu", default="native", help="the CPU every tool is built for (default: the host's)")
    options = parser.parse_args()
    if options.command == "compare" and not options.base:
        parser.error("compare needs the commit to compare against")
    harness.set_cpu(options.cpu)
    if options.command != "scaling":
        pin(options.core)
    check_quiet(options.core, options.force)
    # `peak` first, so that `probes` measures each phase against it, and `duel` before `rivals` and `regret`,
    # so that they find the pairs it lost and the engines `probes` saw.
    commands = (["peak", "noise", "duel", "rivals", "probes", "regret", "perf", "modes", "scaling"]
                if options.command == "all" else [options.command])
    reports = [globals()[command](options) for command in commands]
    if options.command == "all":
        write("profile-all" + (f"-{options.affine.replace(',', '-')}" if options.affine else ""),
              "\n\n".join(reports), None)
    print("\n\n".join(reports))


if __name__ == "__main__":
    main()
