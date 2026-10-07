"""
dinara-align's side of `local_bench.py`: local alignment, `Mode.local(2)`, overlap alignment,
`Mode.overlap(2)`, or a query placed whole in a reference, `Mode.INFIX.with_match_score(2)`, under
`Costs.affine(4, 6, 2)`, each with its CIGAR, one pair at a time on one
thread.

    mojo build -I . benchmarks/local.mojo -o <binary> && <binary> <workload file>

Each pair of the file, `name<TAB>reference<TAB>query`, is aligned; the time is the faster of two
passes over the file, its mean per pair, and the answer the sum and position-weighted sum of the
scores, as `rivals.c` prints them.
"""

from std.sys import argv
from std.time import perf_counter_ns

from dinara_align import Costs, Mode, align


def main() raises:
    var names = List[String]()
    var references = List[String]()
    var queries = List[String]()
    for line in open(String(argv()[1]), "r").read().split("\n"):
        if line.byte_length() == 0:
            continue
        var fields = line.split("\t")
        names.append(String(fields[0]))
        references.append(String(fields[1]))
        queries.append(String(fields[2]))
    var costs = Costs.affine(4, 6, 2)
    var overlap = "overlap" in names[0]
    var infix = "infix" in names[0]
    var mode = Mode.overlap(2) if overlap else (Mode.INFIX.with_match_score(2) if infix else Mode.local(2))
    var best = Float64.MAX
    var total = 0
    var weighted = 0
    for _ in range(2):
        total = 0
        weighted = 0
        var started = perf_counter_ns()
        for index in range(len(references)):
            var score = align(references[index], queries[index], costs, mode).score
            total += score
            weighted += (index + 1) * score
        best = min(best, Float64(perf_counter_ns() - started) / 1e9)
    var task = "overlap" if overlap else ("infix" if infix else "local")
    print("dinara-align", names[0], task, best / Float64(len(references)), String(total, ":", weighted), sep="\t")
