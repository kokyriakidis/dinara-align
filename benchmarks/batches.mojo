# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The side of `batch_bench.py` dinara-align runs: a whole batch of global costs, on one thread or on
`threads` of them.

    mojo build -I . benchmarks/batches.mojo -o <binary> && <binary> <pairs file> <workload> [threads]

The library runs on its caller's thread, as a library should: an application spreads its own calls.
So on one thread the batch is one call, and on more this program spreads it, as the rivals' OpenMP
loop spreads theirs: cut beforehand into a few pieces a thread, each piece one call, the threads taking
the next piece as each finishes.

`illumina-affine` takes `Costs.affine(1, 2, 1)`, Accelign's case study's scores; `illumina-edit` unit
costs. The time is the faster of two passes over the batch, and the answer the sum and position-weighted
sum of the costs, as `batches/rivals.cpp` prints them.

Where a GPU answers, the affine workload is scored there too, by `scores` at the same scores as a
match 0, a mismatch -1 and a gap of `k` letters `-(2 + k)`, each pair's cost its score negated: a row
of its own, `dinara-align GPU`, the faster of five passes after three that bring the device's clocks
up from idle. Its time is the call's: packing the batch, copying it over, scoring, copying back.
"""

from std.sys import argv
from std.time import perf_counter_ns

from std.sys import has_accelerator

from std.atomic import Atomic
from std.os import abort
from max.algorithm import parallelize

from dinara_align import Costs, Mode, Placement, Scoring, distances, scores
from dinara_align.common import next_share


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
    var threads = Int(String(argv()[3])) if len(argv()) > 3 else 1
    var costs = Costs.affine(1, 2, 1) if workload == "illumina-affine" else Costs.edit()
    var pairs = len(references)
    # The pieces a thread takes, cut before the clock starts, as the rivals' pairs are read beforehand.
    var pieces = 1 if threads <= 1 else 4 * threads
    var piece_references = List[List[String]]()
    var piece_queries = List[List[String]]()
    for piece in range(pieces):
        piece_references.append(List[String]())
        piece_queries.append(List[String]())
        for index in range(pairs * piece // pieces, pairs * (piece + 1) // pieces):
            piece_references[piece].append(references[index])
            piece_queries[piece].append(queries[index])
    var best = Float64.MAX
    var total = 0
    var weighted = 0
    for _ in range(2):
        var found = List[Int](length=pairs, fill=0)
        var out = found.unsafe_ptr()
        var taken = Atomic[Int64](0)
        var started = perf_counter_ns()

        def run_pieces(worker: Int) {mut taken, imm piece_references, imm piece_queries, imm costs, imm out, imm pairs, imm pieces, imm threads}:
            """Takes the next piece until none is left, each one call."""
            var last = 0
            while True:
                var share = next_share(taken, pieces, threads, last)
                if share[0] >= pieces:
                    return
                for piece in range(share[0], share[1]):
                    try:
                        var costs_found = distances(piece_references[piece], piece_queries[piece], costs)
                        var first = pairs * piece // pieces
                        for index in range(len(costs_found)):
                            out[unsafe_offset=first + index] = costs_found[index]
                    except error:
                        abort(String("a piece failed: ", error))

        if threads <= 1:
            run_pieces(0)
        else:
            parallelize(run_pieces, threads, threads)
        best = min(best, Float64(perf_counter_ns() - started) / 1e9)
        total = 0
        weighted = 0
        for index in range(len(found)):
            total += found[index]
            weighted += (index + 1) * found[index]
    print("dinara-align", workload, best, String(total, ":", weighted), sep="\t")
    comptime if has_accelerator():
        if workload == "illumina-affine" and threads > 1:
            var scoring = Scoring.uniform(0, -1, -2, -1)
            # The host threads that pack the batch for the device, which this program asks for.
            var device = Placement.on_gpu(0, threads)
            var fastest = Float64.MAX
            for attempt in range(8):
                var started = perf_counter_ns()
                var found = scores(references, queries, scoring, Mode.GLOBAL, placement=device)
                var seconds = Float64(perf_counter_ns() - started) / 1e9
                if attempt >= 3:
                    fastest = min(fastest, seconds)
                total = 0
                weighted = 0
                for index in range(len(found)):
                    total -= found[index]
                    weighted -= (index + 1) * found[index]
            print("dinara-align GPU", workload, fastest, String(total, ":", weighted), sep="\t")
