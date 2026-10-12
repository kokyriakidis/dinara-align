# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""The `probes` view of `profile.py`: runs the runner built with probes over a sample, the copy-counting
shim preloaded where it builds, and turns each pair's probes and counters into a report.

The report, per dataset: the pairs and time by engine; the time by phase; the work each engine did against
what the pair needed (band cells per base, the share swept by rounds that failed, the first and last bound
against the distance, diagonal-transition steps spent on pairs it handed over); the processor's view
(instructions per pair and per cell, IPC, branch and cache misses); allocations, bytes allocated, copied
and set, and page faults per pair, after the first pair, whose setup a warm aligner does not repeat; and
the findings, each pair whose probes match a pattern that costs time, with the slowest pairs in full.
The findings' tests are diagnostic, not tuned: each names a pattern worth a look, nothing more.
"""

import json
import os
import shutil
import statistics
import subprocess
import sys
from pathlib import Path

import harness
from harness import CACHE, HERE, duration

SHIM_SOURCE = HERE / "probes" / "counting.c"
SHIM = CACHE / "profile" / "counting.so"

PHASES = [
    ("encode_ns", "encode"),
    ("diagonal_ns", "diagonal transition, one end"),
    ("two_ended_ns", "diagonal transition, two ends"),
    ("grow_ns", "fronts grown to the distance"),
    ("seeds_ns", "seed heuristic"),
    ("band_ns", "band rounds"),
    ("traceback_ns", "band traceback"),
    ("cigar_ns", "CIGAR"),
    ("full_ns", "full-matrix sweep"),
    ("wavefront_ns", "gap-affine wavefront"),
    ("sweep_ns", "difference sweep"),
    ("search_ns", "end search (prefix, suffix, infix)"),
    ("scored_ns", "scored modes (local, extension, tables)"),
]

EFFICIENCY = [
    # phase, its work in units of the peak, the phase's time, the peak it is measured against
    ("band rounds", lambda r: 64 * r.get("word_columns", 0), "band_ns", "band_cells_per_ns"),
    ("full-matrix sweep", lambda r: 64 * r.get("full_word_columns", 0), "full_ns", "band_cells_per_ns"),
    ("diagonal transition, two ends", lambda r: r.get("two_ended_steps", 0), "two_ended_ns", "diagonal_steps_per_ns"),
    ("diagonal transition, one end", lambda r: r.get("diagonal_steps", 0), "diagonal_ns", "diagonal_steps_per_ns"),
    ("gap-affine wavefront", lambda r: r.get("wavefront_diagonals", 0), "wavefront_ns", "wavefront_diagonals_per_ns"),
    ("difference sweep", lambda r: r.get("sweep_cells", 0), "sweep_ns", "sweep_cells_per_ns"),
]
"""Each phase whose work has a kernel's unit: band cells, diagonal-transition steps, wavefront diagonals,
sweep cells; its rate over the kernel's peak (see `profile.py peak`) is how much of the kernel it gets."""


def peaks() -> dict:
    """The kernels' peaks from the last `profile.py peak`, if any."""
    path = harness.RESULTS / "profile-peak.json"
    return json.loads(path.read_text()) if path.exists() else {}


def shim() -> Path | None:
    """The copy-counting shim, built once where the system takes `LD_PRELOAD`; None elsewhere."""
    if not sys.platform.startswith("linux") or not shutil.which("cc"):
        return None
    if not SHIM.exists() or SHIM.stat().st_mtime < SHIM_SOURCE.stat().st_mtime:
        SHIM.parent.mkdir(parents=True, exist_ok=True)
        # No builtins, so the fallback loops are never turned into calls of the functions they stand in for.
        subprocess.run(["cc", "-O2", "-shared", "-fPIC", "-fno-builtin", "-fno-tree-loop-distribute-patterns",
                        str(SHIM_SOURCE), "-o", str(SHIM), "-ldl"], check=True)
    return SHIM


def run(binary: Path, tool: str, path: Path, budget: float) -> list[dict]:
    """Every pair the runner finished on the sample, as its time, cost, memory growth, probes and counters."""
    environment = dict(os.environ)
    preload = shim()
    if preload is not None:
        environment["LD_PRELOAD"] = str(preload)
    output = subprocess.run([str(binary), "seq", tool, str(budget + harness.GRACE), str(path)],
                            capture_output=True, text=True, env=environment, timeout=budget * 4 + 60).stdout
    records = []
    lengths = None
    for line in output.splitlines():
        fields = line.split("\t")
        if len(fields) != 6 or fields[0] != tool:
            continue
        record = json.loads(fields[5])
        record.update(seconds=float(fields[2]), cost=int(fields[3]), growth=int(fields[4]))
        records.append(record)
    if lengths is None:
        lines = path.read_text().splitlines()
        lengths = [(len(lines[index]) - 1, len(lines[index + 1]) - 1) for index in range(0, len(lines) - 1, 2)]
    for record, (columns, rows) in zip(records, lengths):
        record.update(columns=columns, rows=rows)
    return records


def share(part: float, whole: float) -> str:
    return f"{100 * part / whole:.1f}%" if whole else "—"


def median(values: list[float]) -> float:
    return statistics.median(values) if values else 0.0


def number(value: float) -> str:
    """A count in three significant figures with a thousands separator."""
    if value == 0:
        return "0"
    if abs(value) >= 1000:
        return f"{value:,.0f}"
    return f"{value:.3g}"


def findings(records: list[dict]) -> list[str]:
    """The patterns worth a look, each with the pairs that show it and the time those pairs took."""
    total = sum(record["seconds"] for record in records) or 1.0
    tests = [
        ("bands that needed three rounds or more", lambda r: r.get("rounds", 0) >= 3),
        ("first bound below the distance (a round wasted)",
         lambda r: r.get("rounds", 0) > 1 and r.get("first_bound", 0) < r.get("distance", 0)),
        ("last bound more than 1.25× the distance (a band wider than needed)",
         lambda r: r.get("distance", 0) > 0 and r.get("last_bound", 0) > 1.25 * r["distance"]),
        ("rounds that failed swept more than the round that succeeded",
         lambda r: r.get("last_round_word_columns", 0) and
         r.get("word_columns", 0) - r["last_round_word_columns"] > r["last_round_word_columns"]),
        ("seed construction took longer than the band rounds", lambda r: r.get("seeds_ns", 0) > r.get("band_ns", 1)),
        ("exact seeds built, then rebuilt inexact", lambda r: r.get("seed_rebuilt", 0)),
        ("too many matches: seeds abandoned", lambda r: r.get("seed_gave_up", 0)),
        ("diagonal transition took more than a third of a pair it handed to the band",
         lambda r: r.get("engine") in ("band", "wide full", "short full")
         and r.get("diagonal_ns", 0) + r.get("two_ended_ns", 0) > r["seconds"] * 1e9 / 3),
        ("a band too wide, swept whole instead", lambda r: r.get("engine") == "wide full"),
        ("tiles swept again for the traceback (edges beyond the memory limit)", lambda r: r.get("replayed_tiles", 0)),
        ("tiles recomputed in the traceback (wavefront over its limit)", lambda r: r.get("recomputed_segments", 0)),
        ("gap-affine pair split for the memory limit", lambda r: r.get("splits", 0)),
        ("gap-affine pair swept in differences", lambda r: r.get("sweep_cells", 0)),
        ("instructions per cycle below 1", lambda r: r.get("cycles", 0) and r.get("instructions", 0) < r["cycles"]),
        ("allocations during the pair (after the first)", lambda r: r.get("allocations", 0) > 0),
        ("more bytes copied, moved or set than 4× the sequences",
         lambda r: r.get("copied_bytes", 0) + r.get("moved_bytes", 0) + r.get("set_bytes", 0)
         > 4 * (r["columns"] + r["rows"])),
    ]
    lines = []
    for label, test in tests:
        hits = [index for index, record in enumerate(records) if index > 0 and test(record)]
        if not hits:
            continue
        spent = sum(records[index]["seconds"] for index in hits)
        sample = ", ".join(str(index) for index in hits[:8]) + (" …" if len(hits) > 8 else "")
        lines.append(f"| {label} | {len(hits)} | {share(spent, total)} | {sample} |")
    return lines


def dataset_report(name: str, records: list[dict], facts: dict, top: int) -> list[str]:
    """One dataset's section."""
    if not records:
        return [f"## {name}", "", "No pair finished.", ""]
    total = sum(record["seconds"] for record in records)
    warm = records[1:] or records
    lines = [f"## {name}", "", f"{len(records)} pairs, {duration(total / len(records))} a pair in this build.", ""]

    # Engines.
    engines = {}
    for record in records:
        engine = record.get("engine", "none")
        count, seconds = engines.get(engine, (0, 0.0))
        engines[engine] = (count + 1, seconds + record["seconds"])
    lines += ["| engine | pairs | time |", "| :-- | --: | --: |"]
    lines += [f"| {engine} | {count} | {share(seconds, total)} |"
              for engine, (count, seconds) in sorted(engines.items(), key=lambda item: -item[1][1])]
    lines.append("")

    # Phases, with the instructions and cycles each retired where the counters could be read.
    measured = 0.0
    rows = []
    for key, label in PHASES:
        seconds = sum(record.get(key, 0) for record in records) / 1e9
        if not seconds:
            continue
        measured += seconds
        stem = key.removesuffix("_ns")
        instructions = sum(record.get(f"{stem}_instructions", 0) for record in records)
        cycles = sum(record.get(f"{stem}_cycles", 0) for record in records)
        counted = f"{number(instructions)} | {instructions / cycles:.2f}" if cycles else "— | —"
        rows.append(f"| {label} | {duration(seconds)} | {share(seconds, total)} | {counted} |")
    lines += ["| phase | time | share | instructions | IPC |", "| :-- | --: | --: | --: | --: |", *rows]
    lines.append(f"| outside the phases | {duration(max(total - measured, 0))} | {share(max(total - measured, 0), total)}"
                 " | | |")
    lines.append("")
    # Each kernel-shaped phase's rate against the kernel alone.
    known = peaks()
    rates = []
    for label, work, key, peak in EFFICIENCY:
        done = sum(work(record) for record in records)
        taken = sum(record.get(key, 0) for record in records)
        if done and taken and known.get(peak):
            rate = done / taken
            rates.append(f"| {label} | {rate:.3g} | {known[peak]:.3g} | {100 * rate / known[peak]:.0f}% |")
    if rates:
        lines += ["| phase | its rate, units a ns | the kernel's peak | share of the peak |",
                  "| :-- | --: | --: | --: |", *rates, ""]

    # Work against need.
    work = []
    banded = [record for record in records if record.get("word_columns")]
    if banded:
        cells_per_base = median([64 * record["word_columns"] / max(record["columns"], 1) for record in banded])
        failed = sum(record["word_columns"] - record.get("last_round_word_columns", record["word_columns"])
                     for record in banded) / max(sum(record["word_columns"] for record in banded), 1)
        first = median([record["first_bound"] / max(record.get("distance", 1), 1) for record in banded
                        if record.get("distance")])
        last = median([record["last_bound"] / max(record.get("distance", 1), 1) for record in banded
                       if record.get("distance")])
        rounds = median([record.get("rounds", 0) for record in banded])
        work += [
            f"| band cells per base (median) | {number(cells_per_base)} |",
            f"| band cells swept by rounds that failed | {100 * failed:.1f}% |",
            f"| rounds a band took (median) | {number(rounds)} |",
            f"| first bound over the distance (median) | {first:.2f} |",
            f"| last bound over the distance (median) | {last:.2f} |",
        ]
        origin = [record["heuristic_origin"] / record["distance"] for record in banded
                  if record.get("seeds") and record.get("distance")]
        if origin:
            work.append(f"| seed heuristic at the origin over the distance (median) | {median(origin):.2f} |")
    handed = [record for record in records if record.get("engine") in ("band", "wide full", "short full")]
    if handed:
        spent = sum(record.get("diagonal_ns", 0) + record.get("two_ended_ns", 0) for record in handed) / 1e9
        work.append(f"| diagonal transition spent on pairs it handed over | {share(spent, total)} of all time |")
    affine = [record for record in records if record.get("wavefront_diagonals")]
    if affine:
        steps = [record.get("forward_steps", 0) + record.get("backward_steps", 0) for record in affine]
        work += [
            f"| wavefront diagonals per cost step (median) | "
            f"{number(median([r['wavefront_diagonals'] / max(s, 1) for r, s in zip(affine, steps)]))} |",
            f"| meeting checks per cost step (median) | "
            f"{number(median([r.get('meets', 0) / max(s, 1) for r, s in zip(affine, steps)]))} |",
        ]
    if work:
        lines += ["| work | value |", "| :-- | --: |", *work, ""]

    # The processor and memory, per pair after the first.
    counted = [record for record in warm if "instructions" in record]
    machine = []
    if counted:
        instructions = sum(record["instructions"] for record in counted)
        cycles = sum(record.get("cycles", 0) for record in counted)
        machine += [
            f"| instructions per pair (median) | {number(median([r['instructions'] for r in counted]))} |",
            f"| instructions per base (median) | "
            f"{number(median([r['instructions'] / max(r['columns'] + r['rows'], 1) for r in counted]))} |",
            f"| instructions per cycle | {instructions / max(cycles, 1):.2f} |",
        ]
        if any("branch_misses" in record for record in counted):
            misses = sum(record.get("branch_misses", 0) for record in counted)
            machine.append(f"| branch misses per thousand instructions | {1000 * misses / instructions:.2f} |")
        for key, label in (("l1d_misses", "L1d"), ("llc_misses", "LLC")):
            if any(key in record for record in counted):
                missed = sum(record.get(key, 0) for record in counted)
                machine.append(f"| {label} load misses per thousand instructions | {1000 * missed / instructions:.2f} |")
        banded_counted = [r for r in counted if r.get("word_columns")]
        if banded_counted:
            machine.append(f"| instructions per band cell (median) | "
                           f"{median([r['instructions'] / (64 * r['word_columns']) for r in banded_counted]):.3f} |")
    shimmed = [record for record in warm if "allocations" in record]
    if shimmed:
        machine += [
            f"| allocations per pair (median / most) | {number(median([r['allocations'] for r in shimmed]))} / "
            f"{number(max(r['allocations'] for r in shimmed))} |",
            f"| bytes allocated per pair (median) | {number(median([r['allocated_bytes'] for r in shimmed]))} |",
            f"| bytes copied per pair (median) | {number(median([r.get('copied_bytes', 0) for r in shimmed]))} "
            f"in {number(median([r.get('copies', 0) for r in shimmed]))} calls |",
            f"| bytes moved per pair (median) | {number(median([r.get('moved_bytes', 0) for r in shimmed]))} "
            f"in {number(median([r.get('moves', 0) for r in shimmed]))} calls |",
            f"| bytes set per pair (median) | {number(median([r.get('set_bytes', 0) for r in shimmed]))} "
            f"in {number(median([r.get('sets', 0) for r in shimmed]))} calls |",
        ]
    machine.append(f"| minor page faults per pair (median) | {number(median([r.get('minor_faults', 0) for r in warm]))} |")
    if machine:
        lines += ["| per pair, after the first | value |", "| :-- | --: |", *machine, ""]

    found = findings(records)
    if found:
        lines += ["Findings (the first pair, cold, left out):", "",
                  "| pattern | pairs | their time | pairs (index in the sample) |", "| :-- | --: | --: | :-- |",
                  *found, ""]

    slowest = sorted(range(len(records)), key=lambda index: -records[index]["seconds"])[:top]
    lines += [f"The {len(slowest)} slowest pairs, every probe:", ""]
    for index in slowest:
        record = records[index]
        rest = {key: value for key, value in record.items() if key not in ("seconds", "columns", "rows")}
        lines.append(f"- pair {index}, {record['columns']:,} × {record['rows']:,}, {duration(record['seconds'])}: "
                     f"`{json.dumps(rest)}`")
    lines.append("")
    return lines


def report(found: dict, options) -> str:
    """Every dataset's section, after a summary of where the time goes."""
    note = ("hardware counters" if any("instructions" in r for records, _ in found.values() for r in records)
            else "no hardware counters (perf_event_open unavailable)")
    shimmed = any("allocations" in r for records, _ in found.values() for r in records)
    lines = [
        f"# Probes: what each pair took and did{', at ' + options.affine if options.affine else ''}",
        "",
        harness.build_note().strip(),
        "",
        f"From a build with `-D DINARA_PROBES`, {note}, "
        + ("allocations and copies counted by `probes/counting.c`." if shimmed
           else "no allocation counts (the shim preloads on Linux only)."),
        " The probes' clocks add about 20 ns a phase, so times here run slightly above the plain build's.",
        "",
        "| dataset | pairs | mean | top engine | top phase |",
        "| :-- | --: | --: | :-- | :-- |",
    ]
    sections = []
    for name, (records, facts) in found.items():
        if records:
            total = sum(r["seconds"] for r in records)
            engines = {}
            for r in records:
                engines[r.get("engine", "none")] = engines.get(r.get("engine", "none"), 0.0) + r["seconds"]
            engine = max(engines.items(), key=lambda item: item[1])
            phases = [(label, sum(r.get(key, 0) for r in records) / 1e9) for key, label in PHASES]
            phase = max(phases, key=lambda item: item[1])
            lines.append(f"| {name} | {len(records)} | {duration(total / len(records))} | {engine[0]} "
                         f"({share(engine[1], total)}) | {phase[0]} ({share(phase[1], total)}) |")
        sections += dataset_report(name, records, facts, options.worst)
    return "\n".join(lines + [""] + sections) + "\n"
