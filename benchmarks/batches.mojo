# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The side of `batch_bench.py` dinara-align runs: a whole batch of global costs on every thread.

    mojo build -I . benchmarks/batches.mojo -o <binary> && <binary> <pairs file> <workload>

`illumina-affine` takes `Costs.affine(1, 2, 1)`, Accelign's case study's scores; `illumina-edit` unit
costs. The time is the faster of two passes over the batch, and the answer the sum and position-weighted
sum of the costs, as `batches/rivals.cpp` prints them.
"""

from std.sys import argv
from std.time import perf_counter_ns

from dinara_align import Costs, distances


def main() raises:
    """Scores the file's pairs as the workload asks and prints one row: tool, workload, the batch's
    seconds and the costs' checksum."""
    var references = List[String]()
    var queries = List[String]()
    for line in open(String(argv()[1]), "r").read().split("\n"):
        if line.byte_length() == 0:
            continue
        var fields = line.split("\t")
        references.append(String(fields[1]))
        queries.append(String(fields[2]))
    var workload = String(argv()[2])
    var costs = Costs.affine(1, 2, 1) if workload == "illumina-affine" else Costs.edit()
    var best = Float64.MAX
    var total = 0
    var weighted = 0
    for _ in range(2):
        var started = perf_counter_ns()
        var found = distances(references, queries, costs)
        best = min(best, Float64(perf_counter_ns() - started) / 1e9)
        total = 0
        weighted = 0
        for index in range(len(found)):
            total += found[index]
            weighted += (index + 1) * found[index]
    print("dinara-align", workload, best, String(total, ":", weighted), sep="\t")
