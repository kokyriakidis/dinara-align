# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The side of `ablation.py` dinara-align runs: every pair of each sample file aligned with its traceback,
the mean time per alignment over the faster of two passes, and the costs' sum, which every build of the
ablation must agree on.

    mojo build -I . [-D ABLATE_...=1] benchmarks/ablation.mojo -o <binary> && <binary> <file.seq>...

A sample file holds a pair as a `>` line and a `<` line, as pa-bench writes them.
"""

from std.sys import argv
from std.time import perf_counter_ns

from dinara_align import align


def read_pairs(path: String, mut firsts: List[String], mut seconds: List[String]) raises:
    """The pairs of a pa-bench sample file."""
    with open(path, "r") as handle:
        for line in handle.read().splitlines():
            if line.startswith(">"):
                firsts.append(String(line[byte=1:]))
            elif line.startswith("<"):
                seconds.append(String(line[byte=1:]))


def main() raises:
    """One row a file: its name, the pairs, the mean milliseconds per alignment and the costs' sum."""
    var arguments = argv()
    for index in range(1, len(arguments)):
        var path = String(arguments[index])
        var firsts = List[String]()
        var seconds = List[String]()
        read_pairs(path, firsts, seconds)
        var best = Int.MAX
        var costs = 0
        for _ in range(2):
            costs = 0
            var start = perf_counter_ns()
            for pair in range(len(firsts)):
                costs += align(firsts[pair], seconds[pair]).cost
            best = min(best, Int(perf_counter_ns() - start))
        var parts = path.split("/")
        print(parts[len(parts) - 1], len(firsts), Float64(best) / 1e6 / Float64(len(firsts)), costs)
