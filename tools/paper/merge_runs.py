"""Merges a run that timed dinara-align alone, the rivals replayed from their kept results, with the full run
before it: each rival's figures, the median of its three alternating rounds, from the full run, and dinara-align's
from the new one; every other file, the samples and memory tables among them, from the new run.

    python3 tools/paper/merge_runs.py benchmarks/results/2026-10-11 <new run> benchmarks/results/2026-10-12
"""
import json, shutil, sys
from pathlib import Path
full, new, out = map(Path, sys.argv[1:4])
out.mkdir(parents=True, exist_ok=True)
for item in new.iterdir():
    shutil.copy(item, out / item.name)
for name in ["bench-exact.json", "bench-affine-4-6-2-exact.json"]:
    old = json.loads((full / name).read_text())
    fresh = json.loads((new / name).read_text())
    for dataset, record in fresh.items():
        if dataset not in old:
            continue
        assert record["pairs"] == old[dataset]["pairs"], dataset
        for tool in list(record["tools"]):
            if not tool.startswith("dinara-align") and tool in old[dataset]["tools"]:
                record["tools"][tool] = old[dataset]["tools"][tool]
    (out / name).write_text(json.dumps(fresh))
print("merged into", out)
