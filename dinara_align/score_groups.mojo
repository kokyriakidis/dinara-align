# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
A batch of short pairs scored on the GPU, several pairs to a warp: the shape Accelign's paper
(Kallenborn et al., BMC Bioinformatics 2026) describes for short reads, a group of a warp's lanes
to each pair, every lane holding a run of the pair's columns in registers and handing its last
column to the next lane by shuffle, the group sweeping the pair's matrix anti-diagonal by
anti-diagonal of lanes, a row a lane a step.

A pair takes `lanes` lanes, each `columns_per_lane` columns, so a warp scores `WARP_SIZE // lanes`
pairs at once. The shape is the batch's: of the shapes compiled (see `SHAPES`), the one whose lanes
cover the batch's longest second sequence and update the fewest cells to do it, each pair's rows
plus the steps that fill the group's pipeline, times its columns. One warp of one pair, as
`alignment.strip_pair_kernel` scores, sweeps 256 columns and fills 32 lanes for every pair: a read
of 148 bases used half its updates. A pair wider than the widest shape goes to that kernel.

The groups of a warp step together, as far as the warp's longest pair needs, so every shuffle has
every lane: a lane whose pair has ended, or whose warp has no pair for it, computes nothing that is
read. Scores only, no traceback.

The pairs travel packed: each letter in as few bits as the alphabet needs, two for DNA, sixteen to a
32-bit word, each sequence from a word of its own, so the threads packing a batch never share a word;
and each pair three numbers, the word its first sequence starts at and the two lengths, its second
sequence from the next word after its first. They are packed straight into page-locked memory, which
the device copies from at the bus's full rate: 500,000 short reads went over as 150 MB of bytes from
pageable memory, 20 ms, a fifth of the batch's time.
"""

from std.math import ceildiv, clamp
from std.memory import stack_allocation
from std.memory.pointer import AddressSpace
from max.gpu import WARP_SIZE, barrier, block_idx, thread_idx
from max.gpu.host import DeviceContext
from max.gpu.primitives.warp import shuffle_up, shuffle_xor

from .alignment import AffineGapCosts, AlignmentMode
from .bit_parallel import base_codes, not_bases
from .common import (
    spread,
    DeviceScope,
    MAX_ALPHABET_SIZE,
    ScoreDType,
    SubstitutionDType,
    THREADS_PER_BLOCK,
    UNKNOWN_SYMBOL,
    code_table,
    raise_unknown,
    allocate,
    translate,
    upload,
    zeroed,
)

comptime SHAPE_LANES: Array[Int, 18] = [8, 8, 8, 8, 8, 8, 16, 16, 16, 16, 16, 16, 32, 32, 32, 32, 32, 32]
"""Each compiled shape's lanes a pair."""
comptime SHAPE_COLUMNS: Array[Int, 18] = [6, 8, 10, 12, 14, 16, 6, 8, 10, 12, 14, 16, 6, 8, 10, 12, 14, 16]
"""Each compiled shape's columns a lane: a pair spans 48 to 512 columns, in steps of an eighth or less, so
a read is swept little wider than it is."""
comptime SHAPES = 18


def unpack_kernel[bits: Int](codes: Pointer[UInt32, MutAnyOrigin], letters: Pointer[UInt8, MutAnyOrigin], words: Int32):
    """Each word of the packed tape as its letters, a byte each, one thread a word: the tape crosses
    the bus packed, and the scoring kernel reads a byte where unpacking a letter itself every step cost
    it a fifth of its time."""
    comptime per_word = 32 // bits
    var index = Int(block_idx.x) * THREADS_PER_BLOCK + Int(thread_idx.x)
    if index >= Int(words):
        return
    var word = codes[unsafe_offset=index]
    var unpacked = SIMD[DType.uint8, per_word]()
    comptime for slot in range(per_word):
        unpacked[slot] = UInt8((word >> UInt32(slot * bits)) & UInt32((1 << bits) - 1))
    # A word's letters in one store: a byte a store left a warp's writes sixteen bytes apart.
    letters.unsafe_offset(index * per_word).unsafe_store[alignment=per_word](unpacked)


def group_score_kernel[
    mode: AlignmentMode, lanes: Int, columns_per_lane: Int, per_word: Int
](
    letters: Pointer[UInt8, MutAnyOrigin],
    shapes: Pointer[UInt32, MutAnyOrigin],
    substitutions: Pointer[Scalar[SubstitutionDType], MutAnyOrigin],
    results: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    pairs: Int32,
    alphabet_size: Int32,
    open: Int32,
    extend: Int32,
):
    """`group_scores` over `pairs` pairs."""
    group_scores[mode, lanes, columns_per_lane, per_word](
        letters, shapes, substitutions, results, Int(pairs), alphabet_size, open, extend
    )


@always_inline
def group_scores[
    mode: AlignmentMode, lanes: Int, columns_per_lane: Int, per_word: Int
](
    letters: Pointer[UInt8, MutAnyOrigin],
    shapes: Pointer[UInt32, MutAnyOrigin],
    substitutions: Pointer[Scalar[SubstitutionDType], MutAnyOrigin],
    results: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    pairs: Int,
    alphabet_size: Int32,
    open: Int32,
    extend: Int32,
):
    """Scores `THREADS_PER_BLOCK // lanes` pairs a block, `lanes` lanes a pair (see the module's
    notes), a letter a byte of `letters`, the tape unpacked. Pair `p`'s first sequence starts at word
    `shapes[3 p]`, `per_word` letters a word, and is `shapes[3 p + 1]` letters long, the rows;
    its second, `shapes[3 p + 2]` letters, the columns, starts at the next word after the first's last.

    A global alignment's cells are held shifted by their row plus their column times `extend`, so Gotoh's
    recurrence takes three additions a cell where it takes five: a gap's layer becomes the larger of the
    cell before plus `open - extend` and the layer's own value there, and a cell the largest of the one
    diagonally before plus its substitution less `2 extend`, kept so in the table, and the two layers.
    The matrix's borders all become `open - extend`, and the corner is shifted back. What a GPU with
    fused add-and-max instructions does in hardware, this does by algebra on any. A local alignment's
    cells are held as they are, since its floor of zero would move with every cell."""
    comptime pairs_per_block = THREADS_PER_BLOCK // lanes
    comptime shifted = mode == AlignmentMode.GLOBAL
    var width = Int(alphabet_size)
    var table = stack_allocation[MAX_ALPHABET_SIZE * MAX_ALPHABET_SIZE, Int32, address_space=AddressSpace.SHARED]()
    for index in range(Int(thread_idx.x), width * width, THREADS_PER_BLOCK):
        var substitution = Int32(substitutions[unsafe_offset=index])
        comptime if shifted:
            substitution -= 2 * extend
        table[unsafe_offset=index] = substitution
    barrier()

    var thread = Int(thread_idx.x)
    var member = thread % lanes
    var pair = Int(block_idx.x) * pairs_per_block + thread // lanes
    if Int(block_idx.x) * pairs_per_block >= pairs:
        # A block past the last pair, launched for as many as there might have been.
        return
    var live = pair < pairs
    var first_start = 0
    var second_start = 0
    var rows = 0
    var columns = 0
    if live:
        first_start = Int(shapes[unsafe_offset=3 * pair]) * per_word
        rows = Int(shapes[unsafe_offset=3 * pair + 1])
        columns = Int(shapes[unsafe_offset=3 * pair + 2])
        second_start = first_start + ceildiv(rows, per_word) * per_word

    # The warp steps as far as its longest pair needs, every lane in every step.
    var steps = Int32(rows + lanes)
    var span = UInt32(1)
    while span < UInt32(WARP_SIZE):
        steps = max(steps, shuffle_xor(steps, span))
        span *= 2

    var first_column = member * columns_per_lane
    var owned = clamp(columns - first_column, 0, columns_per_lane)
    var symbols = Array[Int32, columns_per_lane](fill=0)
    # A pair with no columns has no letters of its own to read: the last pair's would lie past the tape.
    if columns > 0:
        comptime for slot in range(columns_per_lane):
            var column = clamp(first_column + slot, 0, columns - 1)
            symbols[slot] = Int32(letters[unsafe_offset=second_start + column])

    # The border: a global alignment's, shifted, `open - extend` but at the origin, which is 0; a local
    # one's 0. A gap layer at the border is a gap's cost more than the cell, so a path never takes it.
    var opening = open - extend
    var border = Int32(0)
    comptime if shifted:
        border = opening
    var scores = Array[Int32, columns_per_lane](fill=border)
    var deletions = Array[Int32, columns_per_lane](fill=border + open + extend)
    var edge_score = border
    var edge_insertion = border + open + extend
    var above_left_carry = border if first_column > 0 else Int32(0)
    var reported = Int32(0)
    var best = Int32(0)
    for step in range(1, Int(steps) + 1):
        var row = step - member
        # The lane before hands over its last column of this row, which it computed a step ago; the
        # group's first lane reads the matrix's left edge instead.
        var left_score = shuffle_up(edge_score, 1)
        var left_insertion = shuffle_up(edge_insertion, 1)
        if member == 0:
            left_score = border if row > 0 else Int32(0)
            left_insertion = border + open + extend
        var above_left = above_left_carry
        above_left_carry = left_score
        if row >= 1 and row <= rows:
            var row_base = Int(letters[unsafe_offset=first_start + row - 1]) * width
            var running_score = left_score
            var running_insertion = left_insertion
            comptime for slot in range(columns_per_lane):
                var diagonal = above_left + table[unsafe_offset=row_base + Int(symbols[slot])]
                var deletion: Int32
                var insertion: Int32
                comptime if shifted:
                    deletion = max(scores[slot] + opening, deletions[slot])
                    insertion = max(running_score + opening, running_insertion)
                else:
                    deletion = max(scores[slot] + open, deletions[slot] + extend)
                    insertion = max(running_score + open, running_insertion + extend)
                var score = max(max(diagonal, deletion), insertion)
                comptime if mode == AlignmentMode.LOCAL:
                    score = max(score, 0)
                    if slot < owned:
                        best = max(best, score)
                above_left = scores[slot]
                scores[slot] = score
                deletions[slot] = deletion
                running_score = score
                running_insertion = insertion
            edge_score = running_score
            edge_insertion = running_insertion
            comptime if shifted:
                if row == rows and owned > 0 and first_column + owned == columns:
                    reported = scores[owned - 1]

    comptime if mode == AlignmentMode.LOCAL:
        # The group's best, every lane of the warp taking part, the groups' halves never mixing.
        var reach = UInt32(1)
        while reach < UInt32(lanes):
            best = max(best, shuffle_xor(best, reach))
            reach *= 2
        if live and member == 0:
            results[unsafe_offset=pair] = best
    else:
        # The lane holding the last column holds the corner, shifted back; with no columns, or no rows,
        # the corner is the border's, a gap over the other side or nothing.
        var letters_total = Int32(rows + columns)
        var corner = Int32(0) if letters_total == 0 else open + (letters_total - 1) * extend
        if rows > 0 and columns > 0:
            corner = reported + letters_total * extend
        var holder = member == 0 if columns == 0 or rows == 0 else (owned > 0 and first_column + owned == columns)
        if live and holder:
            results[unsafe_offset=pair] = corner


def shape_for(longest_rows: Int, longest_columns: Int) -> Int:
    """The compiled shape that spans `longest_columns` with the fewest cell updates a pair: its rows
    and the steps that fill its group, times its columns. -1 when none spans them."""
    var chosen = -1
    var least = Int.MAX
    comptime for shape in range(SHAPES):
        comptime lanes = SHAPE_LANES[shape]
        comptime spanned = lanes * SHAPE_COLUMNS[shape]
        if spanned >= longest_columns:
            var updates = (longest_rows + lanes) * spanned
            if updates < least:
                least = updates
                chosen = shape
    return chosen


struct Chunks(Movable):
    """A batch's pairs in the chunks it crosses to the device in, each packed while the device scores the
    one before: the first as many pairs as the device holds at once in the narrowest pairs a block, each
    after twice the last, so only the first chunk's packing waits and the device, slower a pair than the
    packing, never does. Each chunk is split into a stretch a thread."""

    var starts: List[Int]
    """Each chunk's first stretch, then the number of stretches."""
    var bounds: List[Int]
    """Each stretch's first pair, then the number of pairs."""

    def __init__(out self, scope: DeviceScope, pairs: Int, workers: Int):
        """The chunks of `pairs` pairs packed over `workers` threads for the device of `scope`."""
        var holds = (
            scope.specs.streaming_multiprocessors
            * scope.specs.max_blocks_per_multiprocessor
            * (THREADS_PER_BLOCK // WARP_SIZE)
        )
        self.starts = List[Int]()
        self.bounds = List[Int]()
        var chunk_size = max(holds, 1)
        var first = 0
        while first < pairs:
            var last = min(first + chunk_size, pairs)
            self.starts.append(len(self.bounds))
            var parts = min(workers, last - first)
            for part in range(parts):
                self.bounds.append(first + (last - first) * part // parts)
            first = last
            chunk_size *= 2
        self.starts.append(len(self.bounds))
        self.bounds.append(pairs)


def copy_stream(gpu_id: Int) raises -> Optional[DeviceContext]:
    """A second context on device `gpu_id`, a second stream for the copies: each chunk crosses the bus
    there while the kernels score the chunk before. None on a device whose driver has no streams to
    offer (Metal's), whose copies go on the kernels' own stream, in order."""
    var copier = DeviceContext(device_id=gpu_id)
    try:
        var probe = copier.create_event()
        copier.stream().record_event(probe)
    except:
        return None
    return copier^


def grouped_scores[
    mode: AlignmentMode
](
    scope: DeviceScope,
    firsts: List[String],
    seconds: List[String],
    indices: List[Int],
    alphabet: String,
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    scoring: AffineGapCosts,
    threads: Int,
    gpu_id: Int,
) raises -> Optional[List[Int32]]:
    """The named pairs' scores, several pairs a warp (see the module's notes), or None when the widest
    second sequence among them is wider than any shape spans. A letter outside the alphabet raises as
    `translate` does, for the first pair in order holding one."""
    var pairs = len(indices)
    if pairs == 0:
        return None
    var size = alphabet.byte_length()
    var bits = 2 if size <= 4 else (4 if size <= 16 else 8)
    var per_word = 32 // bits
    var bases = alphabet == "ACGT"
    var workers = max(threads, 1)

    var chunks = Chunks(scope, pairs, workers)
    ref chunk_starts = chunks.starts
    ref bounds = chunks.bounds
    var stretches = len(bounds) - 1

    # Every stretch measured over the threads asked for: its words on the tape, its longest first and second
    # sequences.
    var measures = List[Int](length=3 * stretches, fill=0)
    var measure_ptr = measures.unsafe_ptr()
    var bound_ptr = bounds.unsafe_ptr()

    def measure(
        stretch: Int,
    ) {imm firsts, imm seconds, imm indices, imm per_word, imm measure_ptr, imm bound_ptr}:
        """Stretch `stretch`'s words and longest sides."""
        var words = 0
        var rows = 0
        var columns = 0
        for slot in range(bound_ptr[unsafe_offset=stretch], bound_ptr[unsafe_offset=stretch + 1]):
            var first = firsts[indices[slot]].byte_length()
            var second = seconds[indices[slot]].byte_length()
            words += ceildiv(first, per_word) + ceildiv(second, per_word)
            rows = max(rows, first)
            columns = max(columns, second)
        measure_ptr[unsafe_offset=3 * stretch] = words
        measure_ptr[unsafe_offset=3 * stretch + 1] = rows
        measure_ptr[unsafe_offset=3 * stretch + 2] = columns

    spread(measure, stretches, min(workers, stretches))
    var longest_rows = 0
    var longest_columns = 0
    var words = 0
    var begins = List[Int](capacity=stretches + 1)
    for stretch in range(stretches):
        begins.append(words)
        words += measures[3 * stretch]
        longest_rows = max(longest_rows, measures[3 * stretch + 1])
        longest_columns = max(longest_columns, measures[3 * stretch + 2])
    begins.append(words)
    var chosen = shape_for(longest_rows, longest_columns)
    if chosen < 0:
        return None

    var shapes = scope.context.enqueue_create_host_buffer[DType.uint32](3 * pairs)
    var codes = scope.context.enqueue_create_host_buffer[DType.uint32](max(words, 1))
    var codes_buffer = allocate[DType.uint32](scope, max(words, 1))
    var letters_buffer = allocate[DType.uint8](scope, max(words, 1) * per_word)
    var shapes_buffer = allocate[DType.uint32](scope, 3 * pairs)
    var substitutions_buffer = upload(scope, substitutions)
    var results_buffer = zeroed[ScoreDType](scope, pairs)
    scope.context.synchronize()
    # The kernels wait only for their own chunk's copies.
    var copies = copy_stream(gpu_id)
    var streams = Bool(copies)
    var copier = copies.take() if streams else DeviceContext(device_id=gpu_id)
    var codes_by_byte = code_table(alphabet)
    var failed = List[Bool](length=pairs, fill=False)
    var tape = codes.unsafe_ptr()
    var places = shapes.unsafe_ptr()
    var flags = failed.unsafe_ptr()
    var begin_ptr = begins.unsafe_ptr()

    for chunk in range(len(chunk_starts) - 1):
        var first_stretch = chunk_starts[chunk]
        var stop_stretch = chunk_starts[chunk + 1]

        def pack(
            part: Int,
        ) {
            imm firsts,
            imm seconds,
            imm indices,
            imm codes_by_byte,
            imm bits,
            imm per_word,
            imm bases,
            imm tape,
            imm places,
            imm flags,
            imm begin_ptr,
            imm bound_ptr,
            imm first_stretch,
        }:
            """Stretch `first_stretch + part`'s pairs onto the tape and their places, flagging a pair
            holding a letter outside the alphabet."""
            var stretch = first_stretch + part
            var start = begin_ptr[unsafe_offset=stretch]
            for slot in range(bound_ptr[unsafe_offset=stretch], bound_ptr[unsafe_offset=stretch + 1]):
                var index = indices[slot]
                var first = firsts[index].byte_length()
                var second = seconds[index].byte_length()
                places[unsafe_offset=3 * slot] = UInt32(start)
                places[unsafe_offset=3 * slot + 1] = UInt32(first)
                places[unsafe_offset=3 * slot + 2] = UInt32(second)
                var known = packed_into(firsts[index], codes_by_byte, bits, per_word, bases, tape.unsafe_offset(start))
                start += ceildiv(first, per_word)
                known = (
                    packed_into(seconds[index], codes_by_byte, bits, per_word, bases, tape.unsafe_offset(start))
                    and known
                )
                start += ceildiv(second, per_word)
                if not known:
                    flags[unsafe_offset=slot] = True

        spread(pack, stop_stretch - first_stretch, stop_stretch - first_stretch)

        # The chunk's tape and places over, unpacked, scored: all enqueued, the host on to the next.
        var first_pair = bounds[first_stretch]
        var chunk_pairs = bounds[stop_stretch] - first_pair
        var first_word = begins[first_stretch]
        var chunk_words = begins[stop_stretch] - first_word
        ref copies = copier if streams else scope.context
        if chunk_words > 0:
            copies.enqueue_copy(
                codes_buffer.create_sub_buffer[DType.uint32](first_word, chunk_words),
                codes.create_sub_buffer[DType.uint32](first_word, chunk_words),
            )
        copies.enqueue_copy(
            shapes_buffer.create_sub_buffer[DType.uint32](3 * first_pair, 3 * chunk_pairs),
            shapes.create_sub_buffer[DType.uint32](3 * first_pair, 3 * chunk_pairs),
        )
        if streams:
            var landed = copier.create_event()
            copier.stream().record_event(landed)
            scope.context.stream().enqueue_wait_for(landed)
        comptime for packing in range(3):
            comptime letter_bits = [2, 4, 8][packing]
            if bits == letter_bits:
                if chunk_words > 0:
                    scope.context.enqueue_function[unpack_kernel[letter_bits]](
                        codes_buffer.create_sub_buffer[DType.uint32](first_word, chunk_words).unsafe_ptr(),
                        letters_buffer.create_sub_buffer[DType.uint8](
                            first_word * per_word, chunk_words * per_word
                        ).unsafe_ptr(),
                        Int32(chunk_words),
                        grid_dim=ceildiv(chunk_words, THREADS_PER_BLOCK),
                        block_dim=THREADS_PER_BLOCK,
                    )
                comptime for shape in range(SHAPES):
                    if shape == chosen:
                        comptime lanes = SHAPE_LANES[shape]
                        comptime kernel = group_score_kernel[mode, lanes, SHAPE_COLUMNS[shape], 32 // letter_bits]
                        scope.context.enqueue_function[kernel](
                            letters_buffer.unsafe_ptr(),
                            shapes_buffer.create_sub_buffer[DType.uint32](3 * first_pair, 3 * chunk_pairs).unsafe_ptr(),
                            substitutions_buffer.unsafe_ptr(),
                            results_buffer.create_sub_buffer[ScoreDType](first_pair, chunk_pairs).unsafe_ptr(),
                            Int32(chunk_pairs),
                            Int32(size),
                            scoring.open,
                            scoring.extend,
                            grid_dim=ceildiv(chunk_pairs, THREADS_PER_BLOCK // lanes),
                            block_dim=THREADS_PER_BLOCK,
                        )

    scope.context.synchronize()
    copier.synchronize()
    for slot in range(pairs):
        if failed[slot]:
            raise_unknown(firsts[indices[slot]], seconds[indices[slot]], alphabet)
    # Back in one copy to page-locked memory, and from there in one more.
    var landed = scope.context.enqueue_create_host_buffer[ScoreDType](pairs)
    scope.context.enqueue_copy(landed, results_buffer)
    scope.context.synchronize()
    return List[Int32](Span(unsafe_ptr=landed.unsafe_ptr(), length=pairs))


comptime PAIR_SHIFTS = SIMD[DType.uint32, 16](0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30)
"""Where each of sixteen two-bit codes sits in its word."""


@always_inline
def packed_into(
    text: String,
    codes_by_byte: Array[UInt8, 256],
    bits: Int,
    per_word: Int,
    bases: Bool,
    target: MutPointer[UInt32, _],
) -> Bool:
    """Writes `text`'s letters `bits` bits each, `per_word` to a word, from `target` on, and whether
    every letter is in the alphabet. Over `ACGT` itself, sixteen letters a word at once: their codes
    are those `bit_parallel.base_codes` reads off their ASCII bits, in the alphabet's own order."""
    var bytes = text.unsafe_ptr()
    var length = text.byte_length()
    var unknown = False
    var position = 0
    var word_index = 0
    if bases:
        var others = SIMD[DType.bool, 16](fill=False)
        while position + 16 <= length:
            var chunk = bytes.unsafe_offset(position).unsafe_load[width=16]()
            others |= not_bases[16](chunk)
            target[unsafe_offset=word_index] = (base_codes[16](chunk).cast[DType.uint32]() << PAIR_SHIFTS).reduce_or()
            word_index += 1
            position += 16
        unknown = others.reduce_or()
    while position < length:
        var word = UInt32(0)
        var stop = min(position + per_word, length)
        for at in range(position, stop):
            var code = codes_by_byte[Int(bytes[unsafe_offset=at])]
            unknown = unknown or code == UNKNOWN_SYMBOL
            word |= UInt32(code) << UInt32((at - position) * bits)
        target[unsafe_offset=word_index] = word
        word_index += 1
        position = stop
    return not unknown
