# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""What every engine shares: the element types and markers, how a call's work is spread over threads,
where a call runs, and the device's memory, which every launch takes through one place that refuses
what the device cannot hold.
"""

from max.algorithm import parallelize
from std.atomic import Atomic
from std.ffi import c_int, c_size_t, external_call
from std.memory import stack_allocation
from std.sys.info import CompilationTarget, has_apple_gpu_accelerator, num_logical_cores, size_of

from max.gpu import WARP_SIZE
from max.gpu.host import DeviceAttribute, DeviceBuffer, DeviceContext

from .errors import AlignmentError, ErrorKind

# region Types and markers

comptime ScoreDType = DType.int32
"""A score or a cost in a `Scoring`'s sweeps, on the host and the device alike."""
comptime SymbolDType = DType.uint8
"""A letter, as its place in a `Scoring`'s alphabet."""
comptime SubstitutionDType = DType.int8
"""A substitution table's entry: what aligning one letter with another scores."""
comptime OffsetDType = DType.uint64
"""A place on a batch's tape, every sequence of the batch end to end, so up to their summed lengths."""

comptime GAP_BYTE = Byte(ord("-"))
"""What a gapped row holds where its sequence has no letter."""
comptime UNKNOWN_SYMBOL = UInt8(255)
"""The code of a byte an alphabet lacks: no alphabet holds 255 letters (see `MAX_ALPHABET_SIZE`)."""

comptime MAX_ALPHABET_SIZE = 32
"""The most letters a `Scoring` takes: its table, 32 by 32 bytes, sits in each GPU block's shared memory.
The IUPAC nucleotide codes in either case fit."""

comptime THREADS_PER_BLOCK = 256
"""Threads of a GPU block whose threads combine their results; the sweeps' blocks are one warp."""
comptime WARPS_PER_BLOCK = THREADS_PER_BLOCK // WARP_SIZE
"""The warps of such a block, each handing on one partial result."""

comptime NEGATIVE_INFINITY = Int32.MIN // 4
"""A score no alignment has, with room below it for a few costs added before it is compared."""

comptime UNREACHED = Int32(-(1 << 28))
"""A diagonal or a cell no path reaches: far enough below zero that a few more columns or gap costs keep it
there."""

comptime FIRST_SENTINEL = UInt8(0xFE)
"""Past the first sequence's last letter: never a letter, UTF-8 never holding it, and unequal to
`SECOND_SENTINEL`."""

comptime SECOND_SENTINEL = UInt8(0xFF)
"""Past the second sequence's last letter."""

# endregion Types and markers

# region Threads


@always_inline
def next_share(mut taken: Atomic[Int64], count: Int, workers: Int, mut last: Int) -> Tuple[Int, Int]:
    """The next items a worker of `workers` takes from `count`, `taken` of them already handed out, and
    `last` the size of its last share, 0 before its first.

    A share is half of what is left divided among the workers, as OpenMP's guided schedule deals them,
    so the shares shrink as the work runs out and the last ones balance the workers; but never more
    than twice the worker's last, starting from one. The batch's longest pairs come first and can cost
    several times the rest, and a first share of a twentieth of the batch left one worker holding them
    while the others finished: on the Skylake-X half again the batch's time. Doubling from one spreads
    them over every worker, and a worker reaches the guided size within a few shares, each one taking
    the counter's cache line once."""
    var left = count - Int(taken.load())
    var share = max(min(left // (2 * max(workers, 1)), 2 * last), 1)
    last = share
    var first = Int(taken.fetch_add(Int64(share)))
    return (first, min(first + share, count))


def spread[F: def(Int) -> None](work: F, items: Int, workers: Int):
    """Runs `work` on each of `items` items over `workers` threads; on the caller's own thread, starting
    none, for one worker or one item. The library starts threads only when a caller asks for more than
    one: an application that calls it from threads of its own knows best how to spread its work."""
    if workers <= 1 or items <= 1:
        for item in range(items):
            work(item)
        return
    parallelize(work, items, workers)


def thread_count(asked: Int, items: Int) -> Int:
    """Threads for `items` items when `asked` for so many: at least one, and no more than there are items nor
    than the machine runs at once, past which they would only take turns, each holding its own memory."""
    return max(min(asked, items, hardware_threads()), 1)


def hardware_threads() -> Int:
    """The threads this process may run on at once. On Linux its affinity mask, which a scheduler or a
    container narrows, counts them; elsewhere the logical cores, `sched_getaffinity` being Linux's alone."""
    comptime if CompilationTarget.is_linux():
        comptime WORDS = 16
        var mask = stack_allocation[WORDS, UInt64]()
        for index in range(WORDS):
            mask[unsafe_offset=index] = 0
        if Int(external_call["sched_getaffinity", c_int](c_int(0), c_size_t(WORDS * 8), mask)) == 0:
            var total = 0
            for index in range(WORDS):
                total += Int(mask[unsafe_offset=index].reduce_bit_count())
            return max(total, 1)
    return max(Int(num_logical_cores()), 1)


# endregion Threads

# region Placement


@fieldwise_init
struct Device(Equatable, ImplicitlyCopyable, TrivialRegisterPassable):
    """The host's cores or a GPU."""

    var kind: UInt8

    comptime CPU = Self(0)
    comptime GPU = Self(1)


struct Placement(ImplicitlyCopyable, TrivialRegisterPassable):
    """Where a call runs, which GPU, and how many of the host's threads it may take. Built by `on_cpu` or
    `on_gpu`, so a call on the host never carries a GPU's number."""

    var device: Device
    var gpu_id: Int
    """The GPU's number, zero on the host."""
    var threads: Int
    """Host threads the call may take, at least one and no more than the process can run."""

    def __init__(out self, device: Device, gpu_id: Int, threads: Int):
        """`device`, the GPU numbered `gpu_id` on a GPU, and `threads` held to what the process runs."""
        self.device = device
        self.gpu_id = max(gpu_id, 0) if device == Device.GPU else 0
        self.threads = thread_count(threads, Int.MAX)

    @staticmethod
    def on_cpu(threads: Int) -> Self:
        """The host, over up to `threads` threads."""
        return Self(Device.CPU, 0, threads)

    @staticmethod
    def on_gpu(gpu_id: Int, threads: Int) -> Self:
        """The GPU numbered `gpu_id`, the host's part of the work, its packing, over up to `threads` threads."""
        return Self(Device.GPU, gpu_id, threads)

    @staticmethod
    def default() -> Self:
        """The host, on the caller's own thread: an application spreads its calls over its threads itself,
        and asks for more here only when it wants this call spread too."""
        return Self.on_cpu(1)


# endregion Placement

# region Device memory


@fieldwise_init
struct GpuSpecs(ImplicitlyCopyable, TrivialRegisterPassable):
    """What a launch is sized by, as the GPU reports it."""

    var shared_memory_per_multiprocessor: Int
    """Shared memory on one multiprocessor, in bytes."""
    var reserved_memory_per_block: Int
    """Of that, what no block may claim for itself."""
    var largest_allocation: Int
    """The largest buffer the device allocates at once, in bytes."""
    var streaming_multiprocessors: Int
    var max_blocks_per_multiprocessor: Int
    """Blocks a multiprocessor runs at once."""

    @staticmethod
    def of(context: DeviceContext) raises -> Self:
        """The specs `context`'s GPU reports. Metal reports a block's shared memory alone, every block given
        the same, and not how many blocks run at once: a block's threads over a warp stand in for that,
        which only sets how finely a launch is cut."""
        var largest = Int(context.max_single_alloc_size())
        var multiprocessors = Int(context.get_attribute(DeviceAttribute.MULTIPROCESSOR_COUNT))
        comptime if has_apple_gpu_accelerator():
            var per_block = Int(context.get_attribute(DeviceAttribute.MAX_SHARED_MEMORY_PER_BLOCK))
            var resident = Int(context.get_attribute(DeviceAttribute.MAX_THREADS_PER_BLOCK)) // WARP_SIZE
            return Self(per_block, 0, largest, multiprocessors, resident)
        var shared = Int(context.get_attribute(DeviceAttribute.MAX_SHARED_MEMORY_PER_MULTIPROCESSOR))
        var claimable = Int(context.get_attribute(DeviceAttribute.MAX_SHARED_MEMORY_PER_BLOCK_OPTIN))
        var resident = Int(context.get_attribute(DeviceAttribute.MAX_BLOCKS_PER_MULTIPROCESSOR))
        return Self(shared, shared - claimable, largest, multiprocessors, resident)


struct DeviceScope(Copyable, Movable):
    """A GPU opened for a call, with its specs, read once."""

    var context: DeviceContext
    var specs: GpuSpecs

    def __init__(out self, gpu_id: Int) raises:
        """The GPU numbered `gpu_id`."""
        self.context = DeviceContext(device_id=gpu_id)
        self.specs = GpuSpecs.of(self.context)


def allocate[dtype: DType](scope: DeviceScope, count: Int) raises -> DeviceBuffer[dtype]:
    """A device buffer of `count` elements, one at least so an empty launch needs no case of its own;
    refused, `SEQUENCE_TOO_LONG`, past the largest buffer the device allocates."""
    var length = max(count, 1)
    var bytes = length * size_of[Scalar[dtype]]()
    if bytes > scope.specs.largest_allocation:
        raise AlignmentError(ErrorKind.SEQUENCE_TOO_LONG, String(bytes, " bytes over ", scope.specs.largest_allocation))
    return scope.context.enqueue_create_buffer[dtype](length)


def upload[dtype: DType](scope: DeviceScope, values: ImmSpan[Scalar[dtype], _]) raises -> DeviceBuffer[dtype]:
    """`values` copied to a new device buffer."""
    var buffer = allocate[dtype](scope, len(values))
    if len(values) > 0:
        scope.context.enqueue_copy(buffer, values)
    return buffer^


def filled[dtype: DType](scope: DeviceScope, count: Int, value: Scalar[dtype]) raises -> DeviceBuffer[dtype]:
    """A new device buffer of `count` elements, each `value`."""
    var buffer = allocate[dtype](scope, count)
    scope.context.enqueue_memset(buffer, value)
    return buffer^


def zeroed[dtype: DType](scope: DeviceScope, count: Int) raises -> DeviceBuffer[dtype]:
    """A new device buffer of `count` zeros."""
    return filled[dtype](scope, count, Scalar[dtype](0))


# endregion Device memory

# region Letters


def table_entry(name: StaticString, score: Int) raises AlignmentError -> Scalar[SubstitutionDType]:
    """`score` as a substitution table's entry; refused, `INVALID_SCORING`, past what an entry holds."""
    if score < Int(Scalar[SubstitutionDType].MIN) or score > Int(Scalar[SubstitutionDType].MAX):
        raise AlignmentError(ErrorKind.INVALID_SCORING, String(name, " ", score))
    return Scalar[SubstitutionDType](score)


def uniform_matrix(
    alphabet_size: Int, match_score: Int, mismatch_score: Int
) raises AlignmentError -> List[Scalar[SubstitutionDType]]:
    """The table of `alphabet_size` letters scoring `match_score` for two alike and `mismatch_score` for two
    that differ, row by row."""
    var hit = table_entry("match", match_score)
    var table = List[Scalar[SubstitutionDType]](
        length=alphabet_size * alphabet_size, fill=table_entry("mismatch", mismatch_score)
    )
    for letter in range(alphabet_size):
        table[letter * (alphabet_size + 1)] = hit
    return table^


def code_table(alphabet: String) -> Array[UInt8, 256]:
    """Every byte's place in `alphabet`, `UNKNOWN_SYMBOL` where it lacks the byte."""
    var codes = Array[UInt8, 256](fill=UNKNOWN_SYMBOL)
    var place = 0
    for letter in alphabet.as_bytes():
        codes[Int(letter)] = UInt8(place)
        place += 1
    return codes^


def translate(text: String, alphabet: String) raises AlignmentError -> List[Scalar[SymbolDType]]:
    """`text` as its letters' places in `alphabet`; refused, `UNKNOWN_SYMBOL`, naming the text, when it holds
    a letter the alphabet lacks."""
    var codes = code_table(alphabet)
    var out = List[Scalar[SymbolDType]](capacity=text.byte_length())
    for letter in text.as_bytes():
        var code = codes[Int(letter)]
        if code == UNKNOWN_SYMBOL:
            raise AlignmentError(ErrorKind.UNKNOWN_SYMBOL, text)
        out.append(Scalar[SymbolDType](code))
    return out^


def raise_unknown(first: String, second: String, alphabet: String) raises AlignmentError:
    """Raises what translating `first`, then `second`, raises: for a pair a batch's packing, spread over
    threads, found a letter in that `alphabet` lacks, raised in the batch's order as a serial loop would."""
    _ = translate(first, alphabet)
    _ = translate(second, alphabet)


# endregion Letters
