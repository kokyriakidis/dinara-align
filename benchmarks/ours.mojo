# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The dinara-align runner of `bench.py`, speaking the `seq` protocol every runner there speaks (see
`harness.run_tool`).

    mojo build -I . benchmarks/ours.mojo -o <binary> && <binary> seq <tool> <budget seconds> <file.seq>...
"""

from std.sys import argv
from std.ffi import external_call
from std.os.path import getsize
from std.sys.info import CompilationTarget
from std.time import perf_counter_ns

from dinara_align import Aligner, Costs, alignments


struct SeqFile:
    """A `.seq` file's pairs, a `>` line and a `<` line each, read in one call of exactly the file's size
    as the Rust runners read theirs; each line becomes a `String`, its mark left off, only as its pair
    comes up, so the file and one pair are all that is held. Read through a growing buffer and split
    into lines first, the file left several times its size resident at the peak."""

    var data: List[UInt8]
    """The whole file's bytes."""
    var at: Int
    """Where the next line starts."""

    def __init__(out self, path: String) raises:
        """The file at `path`, read whole, its first line next."""
        var handle = open(path, "r")
        self.data = handle.read_bytes(getsize(path))
        handle.close()
        self.at = 0

    def line(mut self) -> String:
        """The next line, without its first byte, the mark; empty past the last."""
        comptime CHUNK = 16
        # Sixteen bytes at a time while that many are left, then a byte at a time to the newline.
        var length = len(self.data)
        var start = self.at
        var end = start
        var bytes = self.data.unsafe_ptr()
        while end + CHUNK <= length:
            var newline = bytes.unsafe_offset(end).unsafe_load[width=CHUNK]().eq(UInt8(ord("\n")))
            if newline.reduce_or():
                break
            end += CHUNK
        while end < length and bytes[unsafe_offset=end] != UInt8(ord("\n")):
            end += 1
        self.at = min(end + 1, length)
        if end - start < 1:
            return String()
        return String(StringSlice(unsafe_from_utf8=Span(self.data)[start + 1 : end]))

    def next(mut self, mut first: String, mut second: String) -> Bool:
        """The next pair into `first` and `second`, false past the last."""
        if self.at >= len(self.data):
            return False
        first = self.line()
        if self.at >= len(self.data):
            return False
        second = self.line()
        return True


def peak_resident() -> Int:
    """The process's peak resident memory so far, in bytes: `getrusage`'s `ru_maxrss`, after the two
    16-byte times that open `struct rusage`, counts kilobytes on Linux and bytes on macOS. Read before
    and after each alignment, its growth is A*PA2's memory measure."""
    var usage = List[Int64](length=20, fill=0)
    _ = external_call["getrusage", Int32](Int32(0), usage.unsafe_ptr())
    var peak = Int(usage[4])
    comptime if CompilationTarget.is_macos():
        return peak
    return peak * 1024


def seq_mode() raises:
    """A*PA2's evaluation datasets: `seq <tool> <budget seconds> <file>...`, every pair aligned once.

    pa-bench's `.seq` files hold pairs as a `>` line and a `<` line. Each pair is aligned with its
    traceback, once, as pa-bench times every aligner, until the budget is spent, through one `Aligner`
    a file as every rival's runner builds one aligner a file and keeps it; one row per pair,
    flushed as it finishes, gives its file, its time and its cost, so a run stopped mid-pair still
    reports the pairs before it.

    A batch tool instead aligns each file's pairs in one call on this one thread, as the library runs,
    and gives each pair the batch's time shared evenly: a throughput, not one pair's latency.

    A tool named `name:x,o,e` aligns at affine costs instead, as WFA counts them: a mismatch `x` and a
    gap of `k` letters `o + k e`, on one thread, with a CIGAR as the rivals' traceback hands back.

    A tool named with `distance` gives the edit distance alone, no traceback, for the rivals that give
    no alignment.
    """
    var tool = String(argv()[2])
    var affine = False
    var mismatch = 0
    var opening = 0
    var extension = 0
    if ":" in tool:
        var costs = tool.split(":")[1].split(",")
        affine = True
        mismatch = Int(String(costs[0]))
        opening = Int(String(costs[1]))
        extension = Int(String(costs[2]))
    var budget = Float64(String(argv()[3]))
    # A short spin first, so the scheduler has moved this process onto a fast core.
    var spun = perf_counter_ns()
    while perf_counter_ns() - spun < 50_000_000:
        pass
    if "batch" in tool:
        for argument in range(4, len(argv())):
            var path = String(argv()[argument])
            var pairs = SeqFile(path)
            var firsts = List[String]()
            var seconds = List[String]()
            var first = String()
            var second = String()
            while pairs.next(first, second):
                firsts.append(first)
                seconds.append(second)
            var started = perf_counter_ns()
            var costs = List[Int]()
            if affine:
                var model = Costs.affine(mismatch, opening, extension)
                for found in alignments(firsts, seconds, model):
                    costs.append(found.cost)
            else:
                for found in alignments(firsts, seconds):
                    costs.append(found.cost)
            var share = Float64(perf_counter_ns() - started) / 1e9 / Float64(max(len(firsts), 1))
            for cost in costs:
                print(tool, path, share, cost, sep="\t")
        return
    var spent = 0.0
    for argument in range(4, len(argv())):
        if spent >= budget:
            break
        var path = String(argv()[argument])
        var pairs = SeqFile(path)
        var aligner = Aligner()
        var first = String()
        var second = String()
        while spent < budget and pairs.next(first, second):
            var before = peak_resident()
            var started = perf_counter_ns()
            var cost: Int
            if affine:
                cost = aligner.align(first, second, Costs.affine(mismatch, opening, extension)).cost
            elif "distance" in tool:
                cost = aligner.distance(first, second)
            else:
                # A CIGAR, as every rival's traceback hands back, rather than the two gapped rows.
                cost = aligner.align(first, second).cost
            var seconds = Float64(perf_counter_ns() - started) / 1e9
            var growth = peak_resident() - before
            spent += seconds
            print(tool, path, seconds, cost, growth, sep="\t", flush=True)


def main() raises:
    """The `seq` protocol (see `seq_mode`)."""
    if len(argv()) < 5 or String(argv()[1]) != "seq":
        raise Error("usage: ours seq <tool> <budget seconds> <file.seq>...")
    seq_mode()
