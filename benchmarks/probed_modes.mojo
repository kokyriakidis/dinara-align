# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The modes runner of `profile.py modes`: each pair of a `modes_bench.py` workload aligned once as
`modes.mojo` aligns it, through one `Aligner`, with the pair's probes and counters as a JSON object (see
`probed.mojo`). The substitution-table workloads take the free `align`, which keeps no tally, so their rows
carry the counters alone.

    mojo build -D DINARA_PROBES -I . benchmarks/probed_modes.mojo -o <binary> && <binary> <workload.tsv>

Each row: workload, seconds, the answer, columns, rows and the JSON object.
"""

from std.sys import argv
from std.time import perf_counter_ns

from dinara_align import Aligner, Costs, Mode, align
from dinara_align.probes import PROBES

from counters import Counters
from modes import dna_table


def main() raises:
    """Aligns the workload file's pairs, one row each."""
    comptime assert PROBES, "build with -D DINARA_PROBES"
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
    var aligner = Aligner()
    var counters = Counters()
    for index in range(len(references)):
        ref reference = references[index]
        ref query = queries[index]
        counters.start()
        var started = perf_counter_ns()
        var found: Int
        if workload.startswith("local"):
            found = aligner.align(reference, query, Costs.affine(4, 6, 2), Mode.local(2)).score
        elif workload == "overlap-reads":
            found = aligner.align(reference, query, Costs.affine(4, 6, 2), Mode.overlap(2)).score
        elif workload == "infix-reads":
            found = aligner.align(reference, query, Costs.affine(4, 6, 2), Mode.INFIX.with_match_score(2)).score
        elif workload.startswith("infix-edit"):
            found = aligner.align(reference, query, Costs.edit(), Mode.INFIX).cost
        elif workload == "prefix-edit":
            found = aligner.align(reference, query, Costs.edit(), Mode.PREFIX).cost
        elif workload == "infix-affine":
            found = aligner.align(reference, query, Costs.affine(4, 6, 2), Mode.INFIX).cost
        elif workload == "extension":
            found = aligner.align(reference, query, Costs.affine(4, 6, 2), Mode.extension(2)).score
        elif workload == "extension-bonus":
            found = aligner.align(reference, query, Costs.affine(4, 6, 2), Mode.extension(2, end_bonus=50)).score
        elif workload == "two-piece":
            found = aligner.align(reference, query, Costs.two_piece(4, 6, 2, 24, 1)).cost
        elif workload == "table-global":
            found = align(reference, query, scoring, Mode.GLOBAL).score
        else:
            found = align(reference, query, scoring, Mode.local()).score
        var seconds = Float64(perf_counter_ns() - started) / 1e9
        counters.stop()
        ref tally = aligner.space.edit.tally
        var probes = tally.json()
        var counted = counters.json()
        var comma = ", " if probes.byte_length() > 0 and counted.byte_length() > 0 else ""
        print(
            workload,
            seconds,
            found,
            reference.byte_length(),
            query.byte_length(),
            "{" + probes + comma + counted + "}",
            sep="\t",
        )
        tally.clear()
