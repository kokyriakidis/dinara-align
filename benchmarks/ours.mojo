"""
The dinara-align side of the comparison: every workload in the data directory, on the host and,
where an accelerator answers, on the device.

    mojo build -I . benchmarks/ours.mojo -o <binary> && <binary> <data directory>

Prints one tab-separated row per measurement, in the shape `run.py` reads from every tool:
`tool, workload, task, device, seconds, answer`.
"""

from std.sys import argv
from std.time import perf_counter_ns

from dinara_align import (
    AlignmentMode,
    Placement,
    Scoring,
    align,
    alignments,
    edit_alignment,
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


def emit(workload: String, task: String, device: String, started: Int, values: List[Int]):
    var seconds = Float64(perf_counter_ns() - started) / 1e9
    print("dinara-align", workload, task, device, seconds, checksum(values), sep="\t")


def gpu_answers(scoring: Scoring) raises -> Bool:
    """Whether an accelerator serves a real alignment, which also warms the device before timing."""
    try:
        _ = align[GLOBAL]("ACGT", "AGT", scoring, Placement.on_gpu(0, hardware_threads()))
        return True
    except:
        return False


def run_batch(data: Pairs, scoring: Scoring, device: String, placement: Placement) raises:
    """A whole file as one batch, which is how a caller with many pairs would use the package."""
    var started = perf_counter_ns()
    var scored = scores[GLOBAL](data.firsts, data.seconds, scoring, placement)
    var values = List[Int]()
    for value in scored:
        values.append(Int(value))
    emit(data.names[0], "score", device, started, values)

    started = perf_counter_ns()
    var aligned = alignments[GLOBAL](data.firsts, data.seconds, scoring, placement)
    values = List[Int]()
    for result in aligned:
        values.append(Int(result.score))
    emit(data.names[0], "alignment", device, started, values)


def run_pairs(data: Pairs, scoring: Scoring, negate: Bool, device: String, placement: Placement) raises:
    """Each pair as its own workload, so long pairs are timed one at a time."""
    var sign = -1 if negate else 1
    for index in range(len(data.names)):
        var started = perf_counter_ns()
        var value = score[GLOBAL](data.firsts[index], data.seconds[index], scoring, placement)
        emit(data.names[index], "score", device, started, [sign * Int(value)])
        started = perf_counter_ns()
        var aligned = align[GLOBAL](data.firsts[index], data.seconds[index], scoring, placement)
        emit(data.names[index], "alignment", device, started, [sign * Int(aligned.score)])


comptime BATCHES = 20
"""Batches a short bit-parallel call is repeated in; the fastest batch's average is its time."""
comptime BATCH_SECONDS = 0.01
"""About how long one batch runs."""
comptime ONCE_SECONDS = 0.1
"""A call at least this long is timed once: noise is small next to it."""


def call[align: Bool](first: String, second: String, threads: Int) raises -> Int:
    comptime if align:
        return Int(edit_alignment(first, second, threads).score)
    else:
        return edit_distance(first, second, threads)


def measure[align: Bool](first: String, second: String, threads: Int) raises -> Tuple[Float64, Int]:
    """The time of one call, in seconds, and its answer, measured as A*PA's runner measures its own.

    A first call sizes the batches; a short call is then repeated in `BATCHES` batches of about
    `BATCH_SECONDS` each, and the fastest batch's average is the time.
    """
    var started = perf_counter_ns()
    var value = call[align](first, second, threads)
    var once = Float64(perf_counter_ns() - started) / 1e9
    if once >= ONCE_SECONDS:
        return (once, value)
    var size = max(Int(BATCH_SECONDS / max(once, 1e-9)), 1)
    var best = Float64.MAX
    for _ in range(BATCHES):
        started = perf_counter_ns()
        for _ in range(size):
            value = call[align](first, second, threads)
        best = min(best, Float64(perf_counter_ns() - started) / 1e9 / Float64(size))
    return (best, value)


def run_bit_parallel(data: Pairs) raises:
    """The bit-parallel edit distance and alignment, on one thread and on all of them, each its own column.

    One thread is the like-for-like against A*PA, which never forks; all threads is what a caller
    gets. Timed warm and in-process, as A*PA's runner times its own aligners, since a single cold
    call of a few microseconds measures the allocator and the scheduler as much as the aligner.
    """
    # A short spin first, so the scheduler has moved this process onto a fast core.
    var spun = perf_counter_ns()
    while perf_counter_ns() - spun < 200_000_000:
        pass
    var widths: List[Int] = [1, hardware_threads()]
    for threads in widths:
        var tool = String("dinara-align (bit-parallel, ", threads, " thread", "" if threads == 1 else "s", ")")
        for index in range(len(data.names)):
            var scored = measure[False](data.firsts[index], data.seconds[index], threads)
            print(tool, data.names[index], "score", "cpu", scored[0], checksum([scored[1]]), sep="\t")
            var aligned = measure[True](data.firsts[index], data.seconds[index], threads)
            print(tool, data.names[index], "alignment", "cpu", aligned[0], checksum([aligned[1]]), sep="\t")


def write_scoring(directory: String, scoring: Scoring) raises:
    """The DNA default, written out so every other runner scores with exactly these numbers.

    One line each for the alphabet, the gap opening, the gap extension, then the table row-major.
    """
    var text = String(scoring.alphabet, "\n", scoring.gaps.open, "\n", scoring.gaps.extend, "\n")
    for value in scoring.substitutions:
        text += String(Int(value)) + "\n"
    with open(directory + "/dna_scoring.txt", "w") as out:
        out.write(text)


def main() raises:
    var directory = String(argv()[1])
    var dna = Scoring.dna()
    var unit = Scoring.edit_distance()
    write_scoring(directory, dna)

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
        run_pairs(affine, dna, False, device, placement)
        # Edit distance is the unit-cost global recurrence, negated.
        run_pairs(edit, unit, True, device, placement)
    run_bit_parallel(edit)
