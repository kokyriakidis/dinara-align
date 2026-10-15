"""Rewrites the paper's benchmark tables and figure data from bench.py's exact JSON and its tables, then prints
the figures the text quotes; render.py then turns them into paper/.

    python3 tools/paper/tables.py benchmarks/results/2026-10-15   # the run the paper reports (see merge_runs.py)
"""

import json
import re
import sys
from pathlib import Path

results, src = Path(sys.argv[1]), Path(__file__).resolve().parent
unit = json.loads((results / "bench-exact.json").read_text())
affine = json.loads((results / "bench-affine-4-6-2-exact.json").read_text())
unit_md = (results / "bench.md").read_text()

OURS_UNIT = "dinara-align (bit-parallel, 1 thread)"
OURS_AFFINE = "dinara-align (1 thread)"
SIZES = {"10000": "10 kbp", "100000": "100 kbp", "500000": "500 kbp", "1000000": "1 Mbp"}
REAL = ["ont-1k", "ont-10k", "ont-50k", "ont-promethion", "ont-500k", "ont-500k-genvar", "sars-cov-2", "illumina-250", "pacbio-hifi"]
REAL_LABELS = {"sars-cov-2": "SARS-CoV-2", "illumina-250": "Illumina", "pacbio-hifi": "PacBio HiFi"}
SIMULATED = [f"l{length}.n{count}.e{rate}" for length, count in
             [("10000", "10000"), ("100000", "1000"), ("500000", "200"), ("1000000", "100")] for rate in (1, 5, 10, 20)]


def label(name: str) -> str:
    if name in REAL_LABELS:
        return REAL_LABELS[name]
    if name.startswith("l") and ".e" in name:
        return f"{SIZES[name[1:].split('.')[0]]}, {name.rsplit('.e', 1)[1]}%"
    return name


def number(ms: float) -> str:
    """Milliseconds as the tables print them: three significant figures, thousands separated."""
    if ms >= 999.5:
        return f"{float(f'{ms:.3g}'):,.0f}"
    return f"{ms:#.3g}".rstrip(".")


def best(record: dict, columns: list[str]) -> str | None:
    """The column to set in bold: the lowest mean of the tools that finished the most pairs."""
    tools = record["tools"]
    shown = [column for column in columns if tools[column]["finished"] > 0]
    # A tool alone in completing any pair is compared with nothing.
    if len(shown) < 2:
        return None
    most = max(tools[column]["finished"] for column in shown)
    return min((column for column in shown if tools[column]["finished"] == most), key=lambda c: tools[c]["mean"])


def cell(record: dict, column: str, none: str, median: bool = False, bold: bool = False) -> str:
    tool, pairs = record["tools"][column], record["pairs"]
    if tool["finished"] == 0:
        return ended(tool, none, column)
    text = number(1e3 * tool["mean"])
    if tool["finished"] < pairs:
        text += f" [{tool['finished']}/{pairs}]"
    if median:
        text += f" ({number(1e3 * tool['median'])})"
    return f"**{text}**" if bold else text


# Tools whose failed unit-cost runs were checked without the memory limit and crash there too: BSAlign's edit
# distance segfaults on every pair of about 150 kbp or more (ont-promethion's 153 kbp eighth pair, the 500 kbp
# sets), which the harness, seeing a capped runner die, records as stopped by the limit.
CRASHES = {"bsalign"}
UNIT = False


def ended(tool: dict, none: str, column: str = "") -> str:
    """The cell of a tool that completed no pair: stopped by the memory limit, crashed (a run that ended on its
    own, which bench.py's kept results from before statuses record as done), or stopped by the budget."""
    if UNIT and column in CRASHES and tool["status"] in ("memory", "failed", "done"):
        return "crash"
    if tool["status"] == "memory":
        return "o.o.m."
    if tool["status"] in ("failed", "done"):
        return "crash"
    return none


def time_rows(data: dict, names: list[str], columns: list[str], none: str, median: bool = False, labels=label):
    rows = []
    for name in names:
        record = data[name]
        top = best(record, columns)
        rows.append([labels(name)] + [cell(record, column, none, median, column == top) for column in columns])
    return rows


def peak(record: dict, column: str, none: str) -> str:
    tool = record["tools"][column]
    if tool["finished"] == 0:
        return ended(tool, none, column)
    if UNIT and column in CRASHES and tool["status"] in ("memory", "failed"):
        return "crash"
    return f"{tool['peak'] / 2**20:.0f}" if tool["peak"] >= 0 else "—"


def md_table(text: str, heading: str) -> tuple[list[str], dict[str, list[str]]]:
    """The markdown table after `heading`: its column names and each dataset's cells."""
    lines = text[text.index(heading):].splitlines()
    start = next(index for index, line in enumerate(lines) if line.startswith("| dataset"))
    header = [cell.strip() for cell in lines[start].strip("|").split("|")]
    rows = {}
    for line in lines[start + 2:]:
        if not line.startswith("|"):
            break
        cells = [cell.strip() for cell in line.strip("|").split("|")]
        rows[cells[0]] = cells
    return header, rows


def replace_rows(path: Path, table: str, rows: list[list[str]]):
    """Replaces the `rows` list of the table keyed `table` in a paper source."""
    text = path.read_text()
    at = text.index(f'"{table}": {{')
    begin = text.index('"rows": [', at) + len('"rows": [')
    end = text.index("\n  ]", begin)
    body = "".join(f"\n    {json.dumps(row, ensure_ascii=False)}," for row in rows)
    path.write_text(text[:begin] + body + text[end:])


main, supplement = src / "main_text.py", src / "supplement.py"

# Table 2: the real datasets at unit costs, time and peak memory.
UNIT_MAIN = [OURS_UNIT, "a*pa2-full", "a*pa2-simple", "a*pa", "quicked", "edlib", "biwfa", "wfa", "bsalign", "ksw2"]
_, facts = md_table(unit_md, "Datasets, over each sample")


def kbp(name: str) -> str:
    mean = int(facts[name][2].split("/")[1].strip().replace(",", ""))
    size = f"{mean} bp" if mean < 1000 else (f"{mean / 1000:.1f} kbp" if mean < 10_000 else f"{mean / 1000:.0f} kbp")
    return f"{label(name)} ({size})"


blank = [""] * len(UNIT_MAIN)
UNIT = True
rows = [["*Mean time per alignment*"] + blank] + time_rows(unit, REAL, UNIT_MAIN, "t.o.", labels=kbp)
rows += [["*Peak resident memory (MiB)*"] + blank]
rows += [[label(name)] + [peak(unit[name], column, "t.o.") for column in UNIT_MAIN] for name in REAL]
replace_rows(main, "unit", rows)

REAL_HEAD = "*Real data*"
SIMULATED_HEAD = "*Simulated: random sequences, uniform errors (QuickEd)*"


def grouped(rows: list, real: int) -> list:
    """Rows under a heading for the real datasets, the first `real` of them, and one for the simulated ones."""
    width = len(rows[0]) - 1
    out = [[REAL_HEAD] + [""] * width] + rows[:real]
    if len(rows) > real:
        out += [[SIMULATED_HEAD] + [""] * width] + rows[real:]
    return out


UNIT = False
# Table 3: gap-affine, the real datasets and the 10 and 100 kbp ones its rivals finish.
AFFINE_MAIN = [OURS_AFFINE, "wfa", "biwfa", "bsalign", "ksw2"]
replace_rows(main, "affine", time_rows(affine, REAL, AFFINE_MAIN, "t.o."))

# S3 and S4: every dataset.
replace_rows(supplement, "affine_full",
             grouped(time_rows(affine, REAL + SIMULATED, [OURS_AFFINE, "wfa", "biwfa", "bsalign", "ksw2"], "t.o."), len(REAL)))
UNIT = True
UNIT_FULL = [OURS_UNIT, "a*pa2-full", "a*pa2-simple", "a*pa", "edlib", "biwfa", "wfa", "quicked"]
replace_rows(supplement, "unit_full", grouped(time_rows(unit, REAL + SIMULATED, UNIT_FULL, "t.o.", median=True), len(REAL)))

# S5: the samples.
replace_rows(supplement, "datasets", grouped([[label(name)] + facts[name][1:] for name in REAL + SIMULATED], len(REAL)))

# S6: memory, the largest growth and the run's peak.
MEMORY = [OURS_UNIT, "a*pa2-full", "a*pa2-simple", "edlib", "biwfa", "quicked"]
header, growth = md_table(unit_md, "Growth of peak memory aligning a pair")
memory_rows = []
for name in ["ont-500k", "ont-500k-genvar", "l100000.n1000.e20", "l500000.n200.e20", "l1000000.n100.e20"]:
    record, row = unit[name], [label(name)]
    for column in MEMORY:
        tool = record["tools"][column]
        largest = growth[name][header.index(column)].split("/")[-1].replace("MB", "").strip()
        if tool["finished"] == 0:
            row.append(ended(tool, "t.o."))
        elif tool["finished"] < record["pairs"]:
            row.append(f"{largest} [{tool['finished']}/{record['pairs']}]")
        else:
            row.append(f"{largest} ({tool['peak'] / 2**20:.0f})")
    memory_rows.append(row)
replace_rows(supplement, "memory", grouped(memory_rows, 2))

# S3 and S8: the other modes and the batches, the median of the runs bench_logs keeps (modes-1.md, ...), each
# run the faster of two passes.
def suite(kind: str, keys: int) -> dict:
    """Each row's cells, by tool, as (median, least, most) in ms over the runs of modes_bench.py or batch_bench.py."""
    def ms(text: str):
        found = re.match(r"([\d.,]+)\s*(µs|ms|s)$", text.replace("*", "").strip())
        if not found:
            return None
        value = float(found.group(1).replace(",", ""))
        return value / 1000 if found.group(2) == "µs" else value * 1000 if found.group(2) == "s" else value

    runs = []
    for path in sorted(results.glob(f"{kind}-*.md")):
        lines = [line for line in path.read_text().splitlines() if line.startswith("|")]
        head = [cell.strip() for cell in lines[0].strip("|").split("|")]
        rows = {}
        for line in lines[2:]:
            cells = [cell.strip() for cell in line.strip("|").split("|")]
            rows["|".join(cells[:keys])] = {name: ms(cell) for name, cell in zip(head[keys:], cells[keys:])}
        runs.append(rows)
    combined = {}
    for key in runs[0] if runs else []:
        combined[key] = {}
        for tool in runs[0][key]:
            values = [run[key][tool] for run in runs if run[key][tool] is not None]
            combined[key][tool] = sorted(values)[len(values) // 2] if len(values) == len(runs) else None
    return combined


def replace_rows_by(path: Path, table: str, make):
    """`replace_rows` with each new row made from the table's current one, its labels kept."""
    text = path.read_text()
    at = text.index(f'"{table}": {{')
    begin = text.index('"rows": [', at) + len('"rows": [')
    end = text.index("\n  ]", begin)
    rows = json.loads("[" + text[begin:end].rstrip().rstrip(",") + "]")
    replace_rows(path, table, [make(row, index) for index, row in enumerate(rows)])


modes = suite("modes", 1)
if modes:
    order = list(modes)

    def mode_row(row, index):
        tools = modes[order[index]]
        rivals = {name: value for name, value in tools.items() if name != "dinara-align" and value}
        best = min(rivals, key=rivals.get)
        return [row[0], f"**{number(tools['dinara-align'])}**", f"{best} {number(rivals[best])}"]

    replace_rows_by(supplement, "modes", mode_row)
batches = suite("batch", 2)
if batches:
    BATCH_KEYS = ["illumina-affine|one core", "illumina-edit|one core", "illumina-affine|10 cores", "illumina-edit|10 cores"]

    def batch_row(row, index):
        tools = batches[BATCH_KEYS[index]]
        cell = lambda name: number(tools[name]) if tools.get(name) else "–"
        return row[:2] + [f"**{cell('dinara-align')}**", cell("WFA2-lib"), cell("KSW2"), cell("parasail"), cell("Edlib")]

    replace_rows_by(supplement, "batch", batch_row)

# Figure 2's data.
NAMES = {OURS_UNIT: "dinara-align", "a*pa2-full": "A*PA2-full", "a*pa2-simple": "A*PA2-simple", "a*pa": "A*PA",
         "edlib": "Edlib", "biwfa": "BiWFA", "wfa": "WFA", "quicked": "QuickEd"}


def point(name: str, column: str):
    tool, pairs = unit[name]["tools"][column], unit[name]["pairs"]
    if tool["finished"] == 0:
        return None
    value = float(f"{1e3 * tool['mean']:.3g}")
    return [value] if tool["finished"] == pairs else [value, 1]


def series(names: list[str], indent: str) -> str:
    lines = [f'{indent}"{NAMES[column]}": {json.dumps([point(name, column) for name in names]).replace("null", "None")}' for column in NAMES]
    return ",\n".join(lines)


charts = src / "charts_data.py"
text = charts.read_text()
divergence = [f"l100000.n1000.e{rate}" for rate in (1, 5, 10, 20)]
lengths = {rate: [f"l{length}.n{count}.e{rate}" for length, count in
                  [("10000", "10000"), ("100000", "1000"), ("500000", "200"), ("1000000", "100")]] for rate in (5, 20)}
text = re.sub(r"DIVERGENCE = \{\n.*?\}\n", "DIVERGENCE = {\n" + series(divergence, "    ") + "}\n", text, flags=re.S)
text = re.sub(r"    5: \{\n.*?\},\n", "    5: {\n" + series(lengths[5], "        ") + "},\n", text, flags=re.S)
text = re.sub(r"    20: \{\n.*?\},\n", "    20: {\n" + series(lengths[20], "         ") + "},\n", text, flags=re.S)
charts.write_text(text)

# The figures the text quotes.
def ratio_line(data, ours, names, rival):
    out = []
    for name in names:
        a, b = data[name]["tools"][ours], data[name]["tools"][rival]
        if a["finished"] and b["finished"] == data[name]["pairs"] and a["finished"] == data[name]["pairs"]:
            out.append((name, b["mean"] / a["mean"]))
    return out


for name in REAL + SIMULATED:
    record = unit[name]
    ours = record["tools"][OURS_UNIT]
    rivals = {c: t for c, t in record["tools"].items() if not c.startswith("dinara") and t["finished"] == record["pairs"]}
    if rivals:
        fastest = min(rivals, key=lambda c: rivals[c]["mean"])
        print(f"unit {name}: ours {1e3 * ours['mean']:.4g} ms, fastest rival {fastest} "
              f"{1e3 * rivals[fastest]['mean']:.4g} ms, x{rivals[fastest]['mean'] / ours['mean']:.2f}")
for name in REAL + SIMULATED:
    record = affine[name]
    ours = record["tools"][OURS_AFFINE]
    for rival in ("wfa", "biwfa"):
        tool = record["tools"][rival]
        if ours["finished"] == record["pairs"] and tool["finished"] == record["pairs"]:
            print(f"affine {name}: {rival} x{tool['mean'] / ours['mean']:.2f}")
quick = [r for _, r in ratio_line(unit, OURS_UNIT, REAL + SIMULATED, "quicked")]
import math
print(f"quicked: {min(quick):.1f} to {max(quick):.1f}, geometric {math.exp(sum(map(math.log, quick)) / len(quick)):.1f}, "
      f"{len(quick)} datasets")


# How far apart the rounds were: each cell's largest round mean over its smallest, over every completed cell.
def spreads(data, label):
    worst = []
    for name, record in data.items():
        for column, tool in record["tools"].items():
            means = tool.get("round_means") or []
            if len(means) > 1 and min(means) > 0:
                worst.append((max(means) / min(means) - 1, name, column))
    worst.sort(reverse=True)
    if worst:
        median = sorted(w[0] for w in worst)[len(worst) // 2]
        print(f"{label} rounds: median spread {100 * median:.1f}%, largest {[(round(100 * w, 1), n, c) for w, n, c in worst[:5]]}")


spreads(unit, "unit")
spreads(affine, "affine")
