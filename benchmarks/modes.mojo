# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The side of `mode_bench.py` dinara-align runs: each workload's mode, with its CIGAR, one pair at a time on one
thread.

    mojo build -I . benchmarks/modes.mojo -o <binary> && <binary> <workload file>

The workload's name picks the call (see `mode_bench.py`); the time is the faster of two passes over the
file, its mean per pair, and the answer the sum and position-weighted sum of each pair's cost, or of its
score where the mode rewards matches, as `modes/rivals.cpp` prints them.
"""

from std.sys import argv
from std.time import perf_counter_ns

from dinara_align import Costs, Mode, Scoring, align


def dna_table() raises -> Scoring:
    """A match 2, a transition (A and G, C and T) -2, a transversion -4, a gap of `k` letters `-(4 + 2k)`."""
    var cells = List[Int8]()
    comptime BASES = "ACGT"
    for row in range(4):
        for column in range(4):
            if row == column:
                cells.append(2)
            elif (row + column) % 2 == 0:
                # A with G and C with T: their codes differ by two.
                cells.append(-2)
            else:
                cells.append(-4)
    return Scoring.tabulated(BASES, cells^, -4, -2)


def answer(reference: String, query: String, workload: String, scoring: Scoring) raises -> Int:
    """The pair's cost, or its score where the workload's mode rewards matches."""
    if workload.startswith("infix-edit"):
        return align(reference, query, Costs.edit(), Mode.INFIX).cost
    if workload == "prefix-edit":
        return align(reference, query, Costs.edit(), Mode.PREFIX).cost
    if workload == "infix-affine":
        return align(reference, query, Costs.affine(4, 6, 2), Mode.INFIX).cost
    if workload == "extension":
        return align(reference, query, Costs.affine(4, 6, 2), Mode.extension(2)).score
    if workload == "extension-bonus":
        return align(reference, query, Costs.affine(4, 6, 2), Mode.extension(2, end_bonus=50)).score
    if workload == "two-piece":
        return align(reference, query, Costs.two_piece(4, 6, 2, 24, 1)).cost
    if workload == "table-global":
        return align(reference, query, scoring, Mode.GLOBAL).score
    return align(reference, query, scoring, Mode.LOCAL).score


def main() raises:
    """Aligns the workload file's pairs and prints one row: tool, workload, the mean seconds a pair and
    the answers' checksum."""
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
    var workload = names[0]
    var scoring = dna_table()
    var best = Float64.MAX
    var total = 0
    var weighted = 0
    for _ in range(2):
        total = 0
        weighted = 0
        var started = perf_counter_ns()
        for index in range(len(references)):
            var found = answer(references[index], queries[index], workload, scoring)
            total += found
            weighted += (index + 1) * found
        best = min(best, Float64(perf_counter_ns() - started) / 1e9)
    print("dinara-align", workload, best / Float64(len(references)), String(total, ":", weighted), sep="\t")
