"""
The dinara-align side of the comparison: every workload in the data directory, on the host and,
where an accelerator answers, on the device.

    mojo build -I . benchmarks/ours.mojo -o <binary> && <binary> <data directory>

Prints one tab-separated row per measurement, in the shape `run.py` reads from every tool:
`tool, workload, task, device, seconds, answer`. Every measurement is taken warm and in-process
(see `measure`), the way every runner takes its own.
"""

from std.sys import argv
from std.ffi import external_call
from std.os.path import getsize
from std.sys.info import CompilationTarget
from std.time import perf_counter_ns

from dinara_align import (
    AlignmentMode,
    Placement,
    Scoring,
    align,
    alignments,
    edit_alignment,
    edit_alignments,
    edit_cigar,
    edit_distance,
    hardware_threads,
    score,
    scores,
)

comptime GLOBAL = AlignmentMode.GLOBAL


struct Pairs(Movable):
    """One workload file: a workload name per pair, and both sequences."""

    var names: List[String]
    var firsts: List[String]
    var seconds: List[String]

    def __init__(out self, path: String) raises:
        self.names = List[String]()
        self.firsts = List[String]()
        self.seconds = List[String]()
        # `splitlines` also breaks on tabs in Mojo 1.1, so lines are split on the newline alone.
        for line in open(path, "r").read().split("\n"):
            if line.byte_length() == 0:
                continue
            var fields = line.split(chr(9))
            self.names.append(String(fields[0]))
            self.firsts.append(String(fields[1]))
            self.seconds.append(String(fields[2]))


def checksum(values: List[Int]) -> String:
    """The sum and a position-weighted sum, so a reordered or shifted batch cannot pass."""
    var total = 0
    var weighted = 0
    for index in range(len(values)):
        total += values[index]
        weighted += (index + 1) * values[index]
    return String(total, ":", weighted)


comptime BATCH_SCORE = 0
"""The tasks `call` runs: a whole file as one batch, scored or aligned."""
comptime BATCH_ALIGNMENT = 1
comptime PAIR_SCORE = 2
"""One pair, scored or aligned, by the affine-gap kernels."""
comptime PAIR_ALIGNMENT = 3
comptime EDIT_DISTANCE = 4
"""One pair by the bit-parallel edit distance, distance alone or with its alignment."""
comptime EDIT_ALIGNMENT = 5

comptime BATCHES = 20
"""Batches a short call is repeated in; the fastest batch's average is its time."""
comptime BATCH_SECONDS = 0.01
"""About how long one batch runs."""
comptime ONCE_SECONDS = 0.1
"""A call at least this long is timed once: noise is small next to it."""


def call[task: Int](data: Pairs, index: Int, scoring: Scoring, placement: Placement) raises -> String:
    """Runs one task once and returns its answer as the checksum `run.py` compares."""
    comptime if task == BATCH_SCORE:
        var values = List[Int]()
        for value in scores[GLOBAL](data.firsts, data.seconds, scoring, placement):
            values.append(Int(value))
        return checksum(values)
    elif task == BATCH_ALIGNMENT:
        var values = List[Int]()
        for result in alignments[GLOBAL](data.firsts, data.seconds, scoring, placement):
            values.append(Int(result.score))
        return checksum(values)
    elif task == PAIR_SCORE:
        return checksum([Int(score[GLOBAL](data.firsts[index], data.seconds[index], scoring, placement))])
    elif task == PAIR_ALIGNMENT:
        return checksum([Int(align[GLOBAL](data.firsts[index], data.seconds[index], scoring, placement).score)])
    elif task == EDIT_DISTANCE:
        return checksum([edit_distance(data.firsts[index], data.seconds[index])])
    else:
        return checksum([Int(edit_alignment(data.firsts[index], data.seconds[index]).score)])


def measure[
    task: Int
](data: Pairs, index: Int, scoring: Scoring, placement: Placement) raises -> Tuple[Float64, String]:
    """The time of one call, in seconds, and its answer, measured as every runner measures its own.

    A first call sizes the batches; a call shorter than `ONCE_SECONDS` is then repeated in `BATCHES`
    batches of about `BATCH_SECONDS` each, and the fastest batch's average is the time, since
    anything slowing a batch down comes from outside the aligner. A longer call is timed once.
    """
    var started = perf_counter_ns()
    var answer = call[task](data, index, scoring, placement)
    var once = Float64(perf_counter_ns() - started) / 1e9
    if once >= ONCE_SECONDS:
        return (once, answer)
    var size = max(Int(BATCH_SECONDS / max(once, 1e-9)), 1)
    var best = Float64.MAX
    for _ in range(BATCHES):
        started = perf_counter_ns()
        for _ in range(size):
            answer = call[task](data, index, scoring, placement)
        best = min(best, Float64(perf_counter_ns() - started) / 1e9 / Float64(size))
    return (best, answer)


def emit(tool: String, workload: String, task: String, device: String, timed: Tuple[Float64, String]):
    print(tool, workload, task, device, timed[0], timed[1], sep="\t")


def gpu_answers(scoring: Scoring) raises -> Bool:
    """Whether an accelerator serves a real alignment, which also warms the device before timing."""
    try:
        _ = align[GLOBAL]("ACGT", "AGT", scoring, Placement.on_gpu(0, hardware_threads()))
        return True
    except:
        return False


def run_batch(data: Pairs, scoring: Scoring, device: String, placement: Placement) raises:
    """A whole file as one batch, which is how a caller with many pairs would use the package."""
    emit("dinara-align", data.names[0], "score", device, measure[BATCH_SCORE](data, 0, scoring, placement))
    emit("dinara-align", data.names[0], "alignment", device, measure[BATCH_ALIGNMENT](data, 0, scoring, placement))


def run_pairs(data: Pairs, scoring: Scoring, device: String, placement: Placement) raises:
    """Each pair as its own workload, so long pairs are timed one at a time."""
    for index in range(len(data.names)):
        emit("dinara-align", data.names[index], "score", device, measure[PAIR_SCORE](data, index, scoring, placement))
        emit(
            "dinara-align",
            data.names[index],
            "alignment",
            device,
            measure[PAIR_ALIGNMENT](data, index, scoring, placement),
        )


def run_bit_parallel(data: Pairs, scoring: Scoring) raises:
    """The bit-parallel edit distance and alignment, one pair on one thread, as A*PA runs them."""
    var placement = Placement.on_cpu(1)
    var tool = String("dinara-align (bit-parallel, 1 thread)")
    for index in range(len(data.names)):
        emit(tool, data.names[index], "score", "cpu", measure[EDIT_DISTANCE](data, index, scoring, placement))
        emit(tool, data.names[index], "alignment", "cpu", measure[EDIT_ALIGNMENT](data, index, scoring, placement))


def write_scoring(directory: String, scoring: Scoring) raises:
    """The DNA default, written out so every other runner scores with exactly these numbers.

    One line each for the alphabet, the gap opening, the gap extension, then the table row-major.
    """
    var text = String(scoring.alphabet, "\n", scoring.gaps.open, "\n", scoring.gaps.extend, "\n")
    for value in scoring.substitutions:
        text += String(Int(value)) + "\n"
    with open(directory + "/dna_scoring.txt", "w") as out:
        out.write(text)


struct SeqFile:
    """A `.seq` file's pairs, a `>` line and a `<` line each, read in one call of exactly the file's size
    as the Rust runners read theirs; each line becomes a `String`, its mark left off, only as its pair
    comes up, so the file and one pair are all that is held. Read through a growing buffer and split
    into lines first, the file left several times its size resident at the peak."""

    var data: List[UInt8]
    var at: Int

    def __init__(out self, path: String) raises:
        var handle = open(path, "r")
        self.data = handle.read_bytes(getsize(path))
        handle.close()
        self.at = 0

    def line(mut self) -> String:
        """The next line, without its first byte, the mark; empty past the last."""
        comptime CHUNK = 16
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
    traceback, once, as pa-bench times every aligner, until the budget is spent; one row per pair,
    flushed as it finishes, gives its file, its time and its cost, so a run stopped mid-pair still
    reports the pairs before it.

    A batch tool instead aligns each file's pairs in one call across every thread, and gives each
    pair the batch's wall-clock time shared evenly: a throughput, not one pair's latency.

    A tool named `name:x,o,e` aligns at affine costs instead, as WFA counts them: a mismatch `x` and a
    gap of `k` letters `o + k e`, the global Gotoh alignment at scores `0`, `-x`, `-(o + e)` and `-e`,
    on one thread; its cost is the score negated.
    """
    var tool = String(argv()[2])
    var affine = Optional[Scoring]()
    if ":" in tool:
        var costs = tool.split(":")[1].split(",")
        var mismatch = Int(String(costs[0]))
        var opening = Int(String(costs[1]))
        var extension = Int(String(costs[2]))
        affine = Scoring.uniform(0, -mismatch, -(opening + extension), -extension)
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
            var aligned = edit_alignments(firsts, seconds, hardware_threads())
            var share = Float64(perf_counter_ns() - started) / 1e9 / Float64(max(len(firsts), 1))
            for pair in range(len(aligned)):
                print(tool, path, share, Int(aligned[pair].score), sep="\t")
        return
    var spent = 0.0
    for argument in range(4, len(argv())):
        if spent >= budget:
            break
        var path = String(argv()[argument])
        var pairs = SeqFile(path)
        var first = String()
        var second = String()
        while spent < budget and pairs.next(first, second):
            var before = peak_resident()
            var started = perf_counter_ns()
            var cost: Int
            if affine:
                cost = -Int(align[GLOBAL](first, second, affine.value(), Placement.on_cpu(1)).score)
            else:
                # A CIGAR, as every rival's traceback hands back, rather than the two gapped rows.
                cost = edit_cigar(first, second).distance
            var seconds = Float64(perf_counter_ns() - started) / 1e9
            var growth = peak_resident() - before
            spent += seconds
            print(tool, path, seconds, cost, growth, sep="\t", flush=True)


def main() raises:
    if String(argv()[1]) == "seq":
        seq_mode()
        return
    var directory = String(argv()[1])
    var dna = Scoring.dna()
    write_scoring(directory, dna)
    # A short spin first, so the scheduler has moved this process onto a fast core.
    var spun = perf_counter_ns()
    while perf_counter_ns() - spun < 200_000_000:
        pass

    var devices = List[String]()
    devices.append("cpu")
    if gpu_answers(dna):
        devices.append("gpu")

    var reads = Pairs(directory + "/dna_reads.tsv")
    var kilobase = Pairs(directory + "/dna_kilobase.tsv")
    var affine = Pairs(directory + "/dna_affine.tsv")
    var edit = Pairs(directory + "/dna_edit.tsv")
    for device in devices:
        var threads = hardware_threads()
        var placement = Placement.on_gpu(0, threads) if device == "gpu" else Placement.on_cpu(threads)
        run_batch(reads, dna, device, placement)
        run_batch(kilobase, dna, device, placement)
        run_pairs(affine, dna, device, placement)
    # The edit-distance pairs go to the bit-parallel path alone: the affine-gap kernels would answer
    # them with a full quadratic sweep, a minute and more per 100 kbp pair, comparing nothing new.
    run_bit_parallel(edit, dna)
