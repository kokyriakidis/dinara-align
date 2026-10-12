# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The runner of `profile.py probes`: `ours.mojo`'s `seq` protocol for one pair at a time, built with
`-D DINARA_PROBES`, each row carrying a sixth field, a JSON object of the pair's probes (see
`dinara_align/probes.mojo`) and of the counters read around it (see `counters.mojo`). A runner apart from
`ours.mojo`, which must still build against libraries from before the probes (see `profile.py compare`).

    mojo build -D DINARA_PROBES -I . benchmarks/probed.mojo -o <binary>
    <binary> seq <tool> <budget seconds> <file.seq>...
"""

from std.sys import argv
from std.time import perf_counter_ns

from dinara_align import Aligner, Costs
from dinara_align.probes import PROBES

from counters import Counters
from ours import peak_resident, SeqFile


def main() raises:
    """`seq <tool> <budget seconds> <file>...`: each pair aligned once with its traceback, as `ours.mojo`
    aligns it, at affine costs for a tool named `name:x,o,e`, or its distance alone for one named with
    `distance`."""
    if len(argv()) < 5 or String(argv()[1]) != "seq":
        raise Error("usage: probed seq <tool> <budget seconds> <file.seq>...")
    comptime assert PROBES, "build with -D DINARA_PROBES"
    var tool = String(argv()[2])
    var affine = ":" in tool
    var costs = Costs.edit()
    if affine:
        var parts = tool.split(":")[1].split(",")
        costs = Costs.affine(Int(String(parts[0])), Int(String(parts[1])), Int(String(parts[2])))
    var budget = Float64(String(argv()[3]))
    # A short spin first, so the scheduler has moved this process onto a fast core.
    var spun = perf_counter_ns()
    while perf_counter_ns() - spun < 50_000_000:
        pass
    var spent = 0.0
    for argument in range(4, len(argv())):
        if spent >= budget:
            break
        var path = String(argv()[argument])
        var pairs = SeqFile(path)
        var aligner = Aligner()
        var counters = Counters()
        aligner.space.edit.tally.use_counters(counters.instruction_page, counters.cycle_page)
        var first = String()
        var second = String()
        while spent < budget and pairs.next(first, second):
            var before = peak_resident()
            counters.start()
            var started = perf_counter_ns()
            var cost: Int
            if affine:
                cost = aligner.align(first, second, costs).cost
            elif "distance" in tool:
                cost = aligner.distance(first, second)
            else:
                cost = aligner.align(first, second).cost
            var seconds = Float64(perf_counter_ns() - started) / 1e9
            counters.stop()
            var growth = peak_resident() - before
            spent += seconds
            ref tally = aligner.space.edit.tally
            var probes = tally.json()
            var counted = counters.json()
            var comma = ", " if probes.byte_length() > 0 and counted.byte_length() > 0 else ""
            print(tool, path, seconds, cost, growth, "{" + probes + comma + counted + "}", sep="\t", flush=True)
            tally.clear()
