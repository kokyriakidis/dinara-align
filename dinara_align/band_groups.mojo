# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
A batch of global scores on the GPU over a band of diagonals each pair's own cost proves, as the CPU's
lanes sweep one (see `lanes`): a thread a pair, holding `BAND` diagonals of it in registers, so a read
of 150 letters costs 2,400 cells, not the 22,500 of its whole matrix.

A cell on row `r` and diagonal `d`, column minus row, lies on anti-diagonal `t = 2 r + d`, and its three
sources on the ones before: the diagonal step on its own diagonal two back, the step down from the
diagonal above one back, and the step across from the diagonal below one back. So a thread sweeps its
pair anti-diagonal by anti-diagonal, its even-numbered slots on even `t` and its odd-numbered ones on
odd, every slot on every second step, from a band whose lowest diagonal is even.

Cells are held shifted by their anti-diagonal times the gap extension, as `score_groups` holds them by
row plus column, which is the same number: three additions a cell, and every border `open - extend`.

A pair's band is its end diagonal, the main one and those between, one more either side, and the
slack shared either side. Its score is proven when no path off the band could score as well: in costs,
a global score being one once the table's best is folded in (see `scoring.tabled_scores`), a path
visiting a diagonal outside runs a gap out to it and one back. The thread weighs that itself, and gives
up as soon as the band cannot prove it: a path's cost only grows, so once every cell of the anti-diagonal
it is on costs more than the bound, so does the corner.

A pair the band does not prove, or whose end lies too far off the main diagonal for it, goes on its
chunk's list of leftovers, and the chunk's leftovers then go whole, several a warp, by `score_groups`
from the same tape, all on the device while the host packs the next chunk; the order a list takes is
the threads', the scores each pair's own. Any batch under costs that fold to a cheaper gap the more it
is split goes to `score_groups` whole, and leftovers wider than any of its shapes go back unscored.

On an RTX 2070, 500,000 short reads scored in 30 ms where whole matrices took 44, the host's packing
of them most of that; reads too noisy or too long for the band cost a few percent more, packed for a band
that then proves little. A wider band than one thread's registers hold, a group of lanes to a pair
handing its edges on by shuffle, proved slower than whole matrices: 33,000 noisy reads over 128
diagonals, eight lanes each, took 10 ms, several times what their whole matrices take.
"""

from std.math import ceildiv
from std.memory import stack_allocation
from std.memory.pointer import AddressSpace
from std.atomic import Atomic
from max.gpu import barrier, block_idx, thread_idx
from max.gpu.host import DeviceContext

from .alignment import AffineGapCosts, AlignmentMode
from .common import (
    MAX_ALPHABET_SIZE,
    NEGATIVE_INFINITY,
    DeviceScope,
    ScoreDType,
    SubstitutionDType,
    THREADS_PER_BLOCK,
    UNKNOWN_SYMBOL,
    code_table,
    raise_unknown,
    allocate,
    spread,
    translate,
    upload,
    zeroed,
)
from .score_groups import (
    SHAPES,
    SHAPE_COLUMNS,
    SHAPE_LANES,
    Chunks,
    copy_stream,
    group_scores,
    packed_into,
    shape_for,
    unpack_kernel,
)

comptime BAND = 16
"""A pair's diagonals, a register each for its scores, deletions and insertions."""


comptime SHAPE_FIELDS = 4
"""A pair's numbers for the kernel: the word its first sequence starts at on the packed tape, its second
from the next word after its first, its rows and columns, and its band's lowest diagonal."""


def band_score_kernel[
    bits: Int
](
    codes: Pointer[UInt32, MutAnyOrigin],
    shapes: Pointer[Int32, MutAnyOrigin],
    substitutions: Pointer[Scalar[SubstitutionDType], MutAnyOrigin],
    results: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    left_count: Pointer[UInt32, MutAnyOrigin],
    left_places: Pointer[UInt32, MutAnyOrigin],
    left_pairs: Pointer[UInt32, MutAnyOrigin],
    first_pair: Int32,
    pairs: Int32,
    alphabet_size: Int32,
    open: Int32,
    extend: Int32,
    best: Int32,
):
    """Scores pairs `first_pair` on, `pairs` of them, over their bands (see the module's notes), from the
    packed tape, `bits` a letter, pair `p`'s numbers at `shapes[SHAPE_FIELDS p ...]`, its score to
    `results[p]` where no path off the band could score as well, the table's best entry `best`. A pair
    the band does not prove goes on the leftovers' list, at the next of `left_count`'s slots: its index to
    `left_pairs` and its numbers to `left_places` as `score_groups` reads them. A thread reads a letter at
    a time off a word it keeps, loading the next every `32 / bits` steps."""
    comptime per_word = 32 // bits
    comptime mask = UInt32((1 << bits) - 1)
    comptime half = BAND // 2
    var width = Int(alphabet_size)
    var table = stack_allocation[MAX_ALPHABET_SIZE * MAX_ALPHABET_SIZE, Int32, address_space=AddressSpace.SHARED]()
    for index in range(Int(thread_idx.x), width * width, THREADS_PER_BLOCK):
        table[unsafe_offset=index] = Int32(substitutions[unsafe_offset=index]) - 2 * extend
    barrier()

    var pair = Int(block_idx.x) * THREADS_PER_BLOCK + Int(thread_idx.x)
    if pair >= Int(pairs):
        return
    pair += Int(first_pair)
    var first_word = Int(shapes[unsafe_offset=SHAPE_FIELDS * pair])
    var rows = Int(shapes[unsafe_offset=SHAPE_FIELDS * pair + 1])
    var columns = Int(shapes[unsafe_offset=SHAPE_FIELDS * pair + 2])
    var low = Int(shapes[unsafe_offset=SHAPE_FIELDS * pair + 3])
    var second_word = first_word + ceildiv(rows, per_word)
    var end = columns - rows

    @always_inline
    def leave() {imm left_count, imm left_places, imm left_pairs, imm pair, imm first_word, imm rows, imm columns}:
        """The pair on the leftovers' list."""
        var slot = Int(Atomic.fetch_add(left_count, UInt32(1)))
        left_pairs[unsafe_offset=slot] = UInt32(pair)
        left_places[unsafe_offset=3 * slot] = UInt32(first_word)
        left_places[unsafe_offset=3 * slot + 1] = UInt32(rows)
        left_places[unsafe_offset=3 * slot + 2] = UInt32(columns)

    if end < low or end >= low + BAND:
        leave()
        return

    var opening = open - extend
    var scores = Array[Int32, BAND](fill=NEGATIVE_INFINITY)
    var deletions = Array[Int32, BAND](fill=NEGATIVE_INFINITY)
    var insertions = Array[Int32, BAND](fill=NEGATIVE_INFINITY)
    var reported = NEGATIVE_INFINITY

    # On even anti-diagonal `t` slot `2 j` holds row `top - j` and column `left + j`, and on the odd one
    # after slot `2 j + 1` the same row and column `left + j + 1`, `top` and `left` half of `t` less and
    # more the band's lowest diagonal. So a window of the rows' letters and one of the columns', each
    # moving one on a step, serve every slot: two letters a step, not two a cell.
    var row_entries = Array[Int32, half](fill=0)
    var column_codes = Array[Int32, half + 1](fill=0)
    var row_word_index = -1
    var row_word = UInt32(0)
    var column_word_index = -1
    var column_word = UInt32(0)

    @always_inline
    def row_entry(row: Int) {mut row_word_index, mut row_word, imm codes, imm first_word, imm rows, imm width} -> Int32:
        """Row `row`'s letter as the first half of a table entry, or nothing past the rows."""
        if row < 1 or row > rows:
            return 0
        var index = first_word + (row - 1) // per_word
        if index != row_word_index:
            row_word_index = index
            row_word = codes[unsafe_offset=index]
        return Int32((row_word >> UInt32(((row - 1) % per_word) * bits)) & mask) * Int32(width)

    @always_inline
    def column_code(
        column: Int,
    ) {mut column_word_index, mut column_word, imm codes, imm second_word, imm columns} -> Int32:
        """Column `column`'s letter, or nothing past the columns."""
        if column < 1 or column > columns:
            return 0
        var index = second_word + (column - 1) // per_word
        if index != column_word_index:
            column_word_index = index
            column_word = codes[unsafe_offset=index]
        return Int32((column_word >> UInt32(((column - 1) % per_word) * bits)) & mask)

    var top = -(low >> 1)
    var left = low >> 1
    comptime for j in range(half):
        row_entries[j] = row_entry(top - j)
    comptime for k in range(half + 1):
        column_codes[k] = column_code(left + k)

    @always_inline
    def cell[
        edged: Bool
    ](
        slot: Int,
        row: Int,
        column: Int,
        entry: Int32,
        up_score: Int32,
        up_deletion: Int32,
        left_score: Int32,
        left_insertion: Int32,
    ) {
        mut scores,
        mut deletions,
        mut insertions,
        mut reported,
        imm rows,
        imm columns,
        imm end,
        imm opening,
        imm table,
        imm low,
    }:
        """Slot `slot`'s cell on `row` and `column`, its pair's table entry `entry`, from its own last,
        two steps back, and its neighbours', one back; `edged` where the cell may lie off the matrix, on
        its borders or on its last row."""
        comptime if edged:
            if row < 0 or row > rows or column < 0 or column > columns:
                scores[slot] = NEGATIVE_INFINITY
                deletions[slot] = NEGATIVE_INFINITY
                insertions[slot] = NEGATIVE_INFINITY
                return
            if row == 0 or column == 0:
                # A border: the origin, or a gap along the first row or column, shifted.
                scores[slot] = Int32(0) if row == column else opening
                deletions[slot] = NEGATIVE_INFINITY
                insertions[slot] = NEGATIVE_INFINITY
                return
        var deletion = max(up_score + opening, up_deletion)
        var insertion = max(left_score + opening, left_insertion)
        scores[slot] = max(max(scores[slot] + table[unsafe_offset=Int(entry)], deletion), insertion)
        deletions[slot] = deletion
        insertions[slot] = insertion
        comptime if edged:
            if row == rows and low + slot == end:
                reported = scores[slot]

    # The proof's bound: every path through a diagonal just off the band, a gap out to it and one back,
    # costs at least this much, so a band no dearer proves its score.
    var fold_opening = 2 * Int(extend - open)
    var fold_extension = Int(best) - 2 * Int(extend)

    @always_inline
    def gap(letters: Int) {imm fold_opening, imm fold_extension} -> Int:
        """A gap of `letters` letters as a cost, nothing for none."""
        return fold_opening + letters * fold_extension if letters > 0 else 0

    var high = low + BAND - 1
    var off_high = gap(high + 1) + gap(high + 1 - end) if high < columns else Int.MAX
    var off_low = gap(1 - low) + gap(end - low + 1) if low > -rows else Int.MAX
    var bound = min(off_high, off_low)

    var step = 0
    while step <= rows + columns:
        # Inside, every slot's row and column within the matrix, past its borders and short of its last
        # row: most of a pair's steps.
        var inside = top - half + 1 >= 1 and top < rows and left >= 1 and left + half <= columns
        comptime for edge in range(2):
            comptime edged = edge == 1
            if inside != edged:
                # Even slots on this even anti-diagonal: the one above each is its odd neighbour, a step
                # back.
                comptime for j in range(half):
                    comptime slot = 2 * j
                    cell[edged](
                        slot,
                        top - j,
                        left + j,
                        row_entries[j] + column_codes[j],
                        scores[slot + 1],
                        deletions[slot + 1],
                        NEGATIVE_INFINITY if slot == 0 else scores[slot - 1],
                        NEGATIVE_INFINITY if slot == 0 else insertions[slot - 1],
                    )
                comptime for j in range(half):
                    comptime slot = 2 * j + 1
                    cell[edged](
                        slot,
                        top - j,
                        left + j + 1,
                        row_entries[j] + column_codes[j + 1],
                        NEGATIVE_INFINITY if slot == BAND - 1 else scores[slot + 1],
                        NEGATIVE_INFINITY if slot == BAND - 1 else deletions[slot + 1],
                        scores[slot - 1],
                        insertions[slot - 1],
                    )
        # A cell on anti-diagonal `t` held `h` costs `t (best - 2 extend) - 2 h`, and a path's cost only
        # grows: once every cell of these two costs more than the bound, the band cannot prove the pair.
        var even = scores[0]
        var odd = scores[1]
        comptime for j in range(1, half):
            even = max(even, scores[2 * j])
            odd = max(odd, scores[2 * j + 1])
        if min(step * fold_extension - 2 * Int(even), (step + 1) * fold_extension - 2 * Int(odd)) > bound:
            leave()
            return
        # Both windows one on: every slot a row further and a column further.
        top += 1
        left += 1
        comptime for j in range(half - 1, 0, -1):
            row_entries[j] = row_entries[j - 1]
        row_entries[0] = row_entry(top)
        comptime for k in range(half):
            column_codes[k] = column_codes[k + 1]
        column_codes[half] = column_code(left + half)
        step += 2

    # The corner shifted back by its anti-diagonal.
    if reported == NEGATIVE_INFINITY:
        leave()
        return
    var score = reported + Int32(rows + columns) * extend
    results[unsafe_offset=pair] = score
    if Int(best) * (rows + columns) - 2 * Int(score) > bound:
        leave()


def leftover_kernel[
    lanes: Int, columns_per_lane: Int, per_word: Int
](
    letters: Pointer[UInt8, MutAnyOrigin],
    places: Pointer[UInt32, MutAnyOrigin],
    substitutions: Pointer[Scalar[SubstitutionDType], MutAnyOrigin],
    scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    left_count: Pointer[UInt32, MutAnyOrigin],
    alphabet_size: Int32,
    open: Int32,
    extend: Int32,
):
    """A chunk's leftovers whole, as many as `left_count` holds, by `score_groups.group_scores`: launched
    for as many as the chunk has pairs, before the host knows how many the bands left."""
    group_scores[AlignmentMode.GLOBAL, lanes, columns_per_lane, per_word](
        letters, places, substitutions, scores, Int(left_count[]), alphabet_size, open, extend
    )


def scatter_kernel(
    scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    left_pairs: Pointer[UInt32, MutAnyOrigin],
    results: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    left_count: Pointer[UInt32, MutAnyOrigin],
):
    """Each of a chunk's leftovers' scores, as many as `left_count` holds, to its pair's place among the
    results."""
    var slot = Int(block_idx.x) * THREADS_PER_BLOCK + Int(thread_idx.x)
    if slot < Int(left_count[]):
        results[unsafe_offset=Int(left_pairs[unsafe_offset=slot])] = scores[unsafe_offset=slot]


def band_low(rows: Int, columns: Int) -> Int:
    """The lowest diagonal of a pair's band: its end diagonal, the main one and those between and one more
    either side, the slack shared either side, made even. Past the end where those do not fit, so the
    kernel leaves the pair unproven at once."""
    var end = columns - rows
    var needed_low = min(0, end) - 1
    var width = max(0, end) + 1 - needed_low + 1
    if width >= BAND:
        return end + 1
    var low = needed_low - (BAND - width) // 2
    return low - 1 if low % 2 != 0 else low


def banded_scores(
    scope: DeviceScope,
    firsts: List[String],
    seconds: List[String],
    alphabet: String,
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    scoring: AffineGapCosts,
    threads: Int,
    gpu_id: Int,
) raises -> Optional[Tuple[List[Int32], List[Bool]]]:
    """Every pair's global score over its band, or over its whole matrix where the band did not prove it
    (see the module's notes), and which pairs are scored: all but those wider than any of `score_groups`'
    shapes spans, when the band did not prove them. None when the costs fold to gaps the proof cannot
    bound: a gap letter worth no less than the table's best pair, or a gap cheaper split than whole. The
    batch crosses in chunks, as `score_groups` sends its own. A letter outside the alphabet raises as
    `translate` does, for the first pair in order holding one."""
    var size = alphabet.byte_length()
    var best = Int(Int32.MIN)
    for cell in range(size * size):
        best = max(best, Int(substitutions[cell]))
    if scoring.extend < scoring.open or best - 2 * Int(scoring.extend) <= 0:
        return None
    var pairs = len(firsts)
    if pairs == 0:
        return (List[Int32](), List[Bool]())
    var workers = max(threads, 1)
    var bits = 2 if size <= 4 else (4 if size <= 16 else 8)
    var per_word = 32 // bits
    var bases = alphabet == "ACGT"
    var chunks = Chunks(scope, pairs, workers)
    var stretches = len(chunks.bounds) - 1
    var bound_ptr = chunks.bounds.unsafe_ptr()

    # Every stretch's words on the tape and its longest sides, measured over the threads asked for.
    var begins = List[Int](length=stretches + 1, fill=0)
    var longest = List[Int](length=2 * stretches, fill=0)
    var begin_ptr = begins.unsafe_ptr()
    var longest_ptr = longest.unsafe_ptr()

    def measure(
        stretch: Int,
    ) {imm firsts, imm seconds, imm per_word, imm begin_ptr, imm longest_ptr, imm bound_ptr}:
        """Stretch `stretch`'s words and longest sides."""
        var words = 0
        var rows = 0
        var columns = 0
        for index in range(bound_ptr[unsafe_offset=stretch], bound_ptr[unsafe_offset=stretch + 1]):
            var first = firsts[index].byte_length()
            var second = seconds[index].byte_length()
            words += ceildiv(first, per_word) + ceildiv(second, per_word)
            rows = max(rows, first)
            columns = max(columns, second)
        begin_ptr[unsafe_offset=stretch + 1] = words
        longest_ptr[unsafe_offset=2 * stretch] = rows
        longest_ptr[unsafe_offset=2 * stretch + 1] = columns

    spread(measure, stretches, min(workers, stretches))
    var longest_rows = 0
    var longest_columns = 0
    for stretch in range(stretches):
        begins[stretch + 1] += begins[stretch]
        longest_rows = max(longest_rows, longest[2 * stretch])
        longest_columns = max(longest_columns, longest[2 * stretch + 1])
    var words = begins[stretches]

    var codes = scope.context.enqueue_create_host_buffer[DType.uint32](max(words, 1))
    var shapes = scope.context.enqueue_create_host_buffer[DType.int32](max(SHAPE_FIELDS * pairs, 1))
    var codes_buffer = allocate[DType.uint32](scope, max(words, 1))
    var shapes_buffer = allocate[DType.int32](scope, max(SHAPE_FIELDS * pairs, 1))
    var substitutions_buffer = upload(scope, substitutions)
    var results_buffer = allocate[ScoreDType](scope, max(pairs, 1))
    # Each chunk's leftovers in the chunk's own stretch of these, its count its own.
    var chunk_count = len(chunks.starts) - 1
    var left_counts = zeroed[DType.uint32](scope, max(chunk_count, 1))
    var left_places = allocate[DType.uint32](scope, max(3 * pairs, 1))
    var left_pairs = allocate[DType.uint32](scope, max(pairs, 1))
    var left_scores = allocate[ScoreDType](scope, max(pairs, 1))
    # The leftovers go whole in the shape that spans the batch, if one does, from the tape unpacked.
    var chosen = shape_for(longest_rows, longest_columns)
    var letters_buffer = allocate[DType.uint8](scope, max(words, 1) * per_word if chosen >= 0 else 1)
    scope.context.synchronize()
    var copies = copy_stream(gpu_id)
    var streams = Bool(copies)
    var copier = copies.take() if streams else DeviceContext(device_id=gpu_id)
    var codes_by_byte = code_table(alphabet)
    var failed = List[Bool](length=pairs, fill=False)
    var tape = codes.unsafe_ptr()
    var fields = shapes.unsafe_ptr()
    var flags = failed.unsafe_ptr()

    for chunk in range(chunk_count):
        var first_stretch = chunks.starts[chunk]
        var stop_stretch = chunks.starts[chunk + 1]

        def pack(
            part: Int,
        ) {
            imm firsts,
            imm seconds,
            imm codes_by_byte,
            imm bits,
            imm per_word,
            imm bases,
            imm tape,
            imm fields,
            imm flags,
            imm begin_ptr,
            imm bound_ptr,
            imm first_stretch,
        }:
            """Stretch `first_stretch + part`'s pairs onto the tape and their numbers, flagging a pair
            holding a letter outside the alphabet."""
            var stretch = first_stretch + part
            var start = begin_ptr[unsafe_offset=stretch]
            for index in range(bound_ptr[unsafe_offset=stretch], bound_ptr[unsafe_offset=stretch + 1]):
                var rows = firsts[index].byte_length()
                var columns = seconds[index].byte_length()
                fields[unsafe_offset=SHAPE_FIELDS * index] = Int32(start)
                fields[unsafe_offset=SHAPE_FIELDS * index + 1] = Int32(rows)
                fields[unsafe_offset=SHAPE_FIELDS * index + 2] = Int32(columns)
                fields[unsafe_offset=SHAPE_FIELDS * index + 3] = Int32(band_low(rows, columns))
                var known = packed_into(firsts[index], codes_by_byte, bits, per_word, bases, tape.unsafe_offset(start))
                start += ceildiv(rows, per_word)
                known = (
                    packed_into(seconds[index], codes_by_byte, bits, per_word, bases, tape.unsafe_offset(start))
                    and known
                )
                start += ceildiv(columns, per_word)
                if not known:
                    flags[unsafe_offset=index] = True

        spread(pack, stop_stretch - first_stretch, stop_stretch - first_stretch)

        # The chunk's tape and numbers over and scored: all enqueued, the host on to the next.
        var first_pair = chunks.bounds[first_stretch]
        var chunk_pairs = chunks.bounds[stop_stretch] - first_pair
        var first_word = begins[first_stretch]
        var chunk_words = begins[stop_stretch] - first_word
        ref copying = copier if streams else scope.context
        if chunk_words > 0:
            copying.enqueue_copy(
                codes_buffer.create_sub_buffer[DType.uint32](first_word, chunk_words),
                codes.create_sub_buffer[DType.uint32](first_word, chunk_words),
            )
        copying.enqueue_copy(
            shapes_buffer.create_sub_buffer[DType.int32](SHAPE_FIELDS * first_pair, SHAPE_FIELDS * chunk_pairs),
            shapes.create_sub_buffer[DType.int32](SHAPE_FIELDS * first_pair, SHAPE_FIELDS * chunk_pairs),
        )
        if streams:
            var landed = copier.create_event()
            copier.stream().record_event(landed)
            scope.context.stream().enqueue_wait_for(landed)
        # The chunk's bands, then its leftovers whole and their scores to their places.
        var count = left_counts.create_sub_buffer[DType.uint32](chunk, 1).unsafe_ptr()
        var places = left_places.create_sub_buffer[DType.uint32](3 * first_pair, 3 * chunk_pairs).unsafe_ptr()
        var left = left_pairs.create_sub_buffer[DType.uint32](first_pair, chunk_pairs).unsafe_ptr()
        var scores = left_scores.create_sub_buffer[ScoreDType](first_pair, chunk_pairs).unsafe_ptr()
        comptime for packing in range(3):
            comptime letter_bits = [2, 4, 8][packing]
            if bits == letter_bits:
                scope.context.enqueue_function[band_score_kernel[letter_bits]](
                    codes_buffer.unsafe_ptr(),
                    shapes_buffer.unsafe_ptr(),
                    substitutions_buffer.unsafe_ptr(),
                    results_buffer.unsafe_ptr(),
                    count,
                    places,
                    left,
                    Int32(first_pair),
                    Int32(chunk_pairs),
                    Int32(size),
                    scoring.open,
                    scoring.extend,
                    Int32(best),
                    grid_dim=ceildiv(chunk_pairs, THREADS_PER_BLOCK),
                    block_dim=THREADS_PER_BLOCK,
                )
                if chosen >= 0 and chunk_words > 0:
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
                        comptime kernel = leftover_kernel[lanes, SHAPE_COLUMNS[shape], 32 // letter_bits]
                        scope.context.enqueue_function[kernel](
                            letters_buffer.unsafe_ptr(),
                            places,
                            substitutions_buffer.unsafe_ptr(),
                            scores,
                            count,
                            Int32(size),
                            scoring.open,
                            scoring.extend,
                            grid_dim=ceildiv(chunk_pairs, THREADS_PER_BLOCK // lanes),
                            block_dim=THREADS_PER_BLOCK,
                        )
        if chosen >= 0:
            scope.context.enqueue_function[scatter_kernel](
                scores,
                left,
                results_buffer.unsafe_ptr(),
                count,
                grid_dim=ceildiv(chunk_pairs, THREADS_PER_BLOCK),
                block_dim=THREADS_PER_BLOCK,
            )

    var scored = scope.context.enqueue_create_host_buffer[ScoreDType](max(pairs, 1))
    scope.context.enqueue_copy(scored, results_buffer)
    scope.context.synchronize()
    copier.synchronize()
    for index in range(pairs):
        if failed[index]:
            raise_unknown(firsts[index], seconds[index], alphabet)
    # With no shape to score them whole, the leftovers go back unscored.
    var settled = List[Bool](length=pairs, fill=True)
    if chosen < 0:
        var counted = scope.context.enqueue_create_host_buffer[DType.uint32](max(chunk_count, 1))
        var unscored = scope.context.enqueue_create_host_buffer[DType.uint32](max(pairs, 1))
        scope.context.enqueue_copy(counted, left_counts)
        scope.context.enqueue_copy(unscored, left_pairs)
        scope.context.synchronize()
        var counts = List[UInt32](Span(unsafe_ptr=counted.unsafe_ptr(), length=chunk_count))
        var lefts = List[UInt32](Span(unsafe_ptr=unscored.unsafe_ptr(), length=pairs))
        for chunk in range(chunk_count):
            var first_pair = chunks.bounds[chunks.starts[chunk]]
            for slot in range(Int(counts[chunk])):
                settled[Int(lefts[first_pair + slot])] = False
    return (List[Int32](Span(unsafe_ptr=scored.unsafe_ptr(), length=pairs)), settled^)
