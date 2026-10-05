"""
The dinara-align side of the comparison: every workload in the data directory, on the host and,
where an accelerator answers, on the device.

    mojo build -I . benchmarks/ours.mojo -o <binary> && <binary> <data directory>

Prints one tab-separated row per measurement, in the shape `run.py` reads from every tool:
`tool, workload, task, device, seconds, answer`. Every measurement is taken warm and in-process
(see `measure`), the way every runner takes its own.
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


def call[
    task: Int
](data: Pairs, index: Int, scoring: Scoring, placement: Placement, threads: Int) raises -> String:
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
        return checksum([edit_distance(data.firsts[index], data.seconds[index], threads)])
    else:
        return checksum([Int(edit_alignment(data.firsts[index], data.seconds[index], threads).score)])


def measure[
    task: Int
](data: Pairs, index: Int, scoring: Scoring, placement: Placement, threads: Int) raises -> Tuple[Float64, String]:
    """The time of one call, in seconds, and its answer, measured as every runner measures its own.

    A first call sizes the batches; a call shorter than `ONCE_SECONDS` is then repeated in `BATCHES`
    batches of about `BATCH_SECONDS` each, and the fastest batch's average is the time, since
    anything slowing a batch down comes from outside the aligner. A longer call is timed once.
    """
    var started = perf_counter_ns()
    var answer = call[task](data, index, scoring, placement, threads)
    var once = Float64(perf_counter_ns() - started) / 1e9
    if once >= ONCE_SECONDS:
        return (once, answer)
    var size = max(Int(BATCH_SECONDS / max(once, 1e-9)), 1)
    var best = Float64.MAX
    for _ in range(BATCHES):
        started = perf_counter_ns()
        for _ in range(size):
            answer = call[task](data, index, scoring, placement, threads)
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
    emit("dinara-align", data.names[0], "score", device, measure[BATCH_SCORE](data, 0, scoring, placement, 1))
    emit("dinara-align", data.names[0], "alignment", device, measure[BATCH_ALIGNMENT](data, 0, scoring, placement, 1))


def run_pairs(data: Pairs, scoring: Scoring, device: String, placement: Placement) raises:
    """Each pair as its own workload, so long pairs are timed one at a time."""
    for index in range(len(data.names)):
        emit("dinara-align", data.names[index], "score", device, measure[PAIR_SCORE](data, index, scoring, placement, 1))
        emit(
            "dinara-align",
            data.names[index],
            "alignment",
            device,
            measure[PAIR_ALIGNMENT](data, index, scoring, placement, 1),
        )


def run_bit_parallel(data: Pairs, scoring: Scoring) raises:
    """The bit-parallel edit distance and alignment, on one thread and on all of them, each its own column.

    One thread is the like-for-like against A*PA, which never forks; all threads is what a caller gets.
    """
    var placement = Placement.on_cpu(1)
    var widths: List[Int] = [1, hardware_threads()]
    for threads in widths:
        var tool = String("dinara-align (bit-parallel, ", threads, " thread", "" if threads == 1 else "s", ")")
        for index in range(len(data.names)):
            emit(tool, data.names[index], "score", "cpu", measure[EDIT_DISTANCE](data, index, scoring, placement, threads))
            emit(
                tool,
                data.names[index],
                "alignment",
                "cpu",
                measure[EDIT_ALIGNMENT](data, index, scoring, placement, threads),
            )


def write_scoring(directory: String, scoring: Scoring) raises:
    """The DNA default, written out so every other runner scores with exactly these numbers.

    One line each for the alphabet, the gap opening, the gap extension, then the table row-major.
    """
    var text = String(scoring.alphabet, "\n", scoring.gaps.open, "\n", scoring.gaps.extend, "\n")
    for value in scoring.substitutions:
        text += String(Int(value)) + "\n"
    with open(directory + "/dna_scoring.txt", "w") as out:
        out.write(text)


def seq_mode() raises:
    """A*PA2's evaluation datasets: `seq <tool> <budget seconds> <file>...`, every pair aligned once.

    pa-bench's `.seq` files hold pairs as a `>` line and a `<` line. Each pair is aligned with its
    traceback, once, as pa-bench times every aligner, until the budget is spent; one row per pair,
    flushed as it finishes, gives its file, its time and its cost, so a run stopped mid-pair still
    reports the pairs before it.
    """
    var tool = String(argv()[2])
    var threads = 1 if tool.endswith("1 thread)") else hardware_threads()
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
        var lines = open(path, "r").read().split("\n")
        var index = 0
        while index + 1 < len(lines) and spent < budget:
            var first = String(lines[index][byte=1:])
            var second = String(lines[index + 1][byte=1:])
            var started = perf_counter_ns()
            var aligned = edit_alignment(first, second, threads)
            var seconds = Float64(perf_counter_ns() - started) / 1e9
            spent += seconds
            print(tool, path, seconds, Int(aligned.score), sep="\t", flush=True)
            index += 2


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
