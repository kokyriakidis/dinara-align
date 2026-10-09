# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Gotoh's affine-gap alignment under a substitution table, on the GPU, with the same recurrence, decisions
and tie rules as the host's (see `gotoh`).

A warp sweeps the matrix in strips `STRIP_WIDTH` columns wide: lane `t` holds `STRIP_COLUMNS` columns in
registers and works one row behind lane `t - 1`, so the cell to its left arrives by one shuffle, and the
strip runs to its end without a barrier. Two engines are built on that:

- A batch of pairs, a warp each (`device_scores`, `device_alignments`), the column left of each strip
  carried in shared memory, so a pair's first sequence is as long as that memory holds (see
  `band_length`). Scoring keeps the last row alone; aligning packs every cell's decision, four bits,
  eight to a word, and walks them back on the device.
- One pair of any size (`device_score`, `device_align`): `gotoh.linear_path` with every half's sweep on
  the device, cut into tiles a warp each; all the frames of a level of the split sweep together, one
  launch per anti-diagonal of tiles, the launches standing in for the grid-wide barrier Mojo lacks.
"""

from std.math import ceildiv, clamp
from std.memory import stack_allocation
from std.memory.pointer import AddressSpace
from std.sys.info import size_of

from max.gpu import WARP_SIZE, barrier, block_dim, block_idx, lane_id, thread_idx
from max.gpu.primitives.warp import shuffle_down, shuffle_up, shuffle_xor
from max.gpu.host import DeviceBuffer, FuncAttribute
from max.gpu.memory import external_memory

from .common import (
    DeviceScope,
    GAP_BYTE,
    GpuSpecs,
    MAX_ALPHABET_SIZE,
    NEGATIVE_INFINITY,
    OffsetDType,
    Placement,
    ScoreDType,
    SubstitutionDType,
    SymbolDType,
    THREADS_PER_BLOCK,
    WARPS_PER_BLOCK,
    allocate,
    spread,
    upload,
    zeroed,
)
from .errors import AlignmentError, ErrorKind
from .gotoh import (
    AffineGapCosts,
    AlignmentMode,
    Cell,
    Cut,
    Decision,
    Frame,
    GapRun,
    GappedAlignment,
    Layer,
    PathView,
    Rectangle,
    RowPath,
    SweepHalf,
    advance,
    beats,
    decide,
    gotoh_cell,
    solve_outright,
)

# region Sizes

comptime ChangeDType = DType.uint32
"""A word of packed decisions."""
comptime DECISION_BITS = 4
comptime STRIP_COLUMNS = size_of[Scalar[ChangeDType]]() * 8 // DECISION_BITS
"""Columns a lane holds: as many as one word packs decisions for."""
comptime STRIP_LANES = WARP_SIZE
comptime STRIP_WIDTH = STRIP_COLUMNS * STRIP_LANES
comptime TILE_SIDE = STRIP_WIDTH
"""A tile's width: one strip."""
comptime MIN_TILE_HEIGHT = STRIP_LANES
"""A tile's least height: a tile of `h` rows takes `h + STRIP_LANES` steps, the last lane starting that late."""

comptime CARRY_LAYERS = 2
"""A batch strip carries the score and insertion layers of the column to its left, a row each."""
comptime STATIC_SHARED_BYTES = MAX_ALPHABET_SIZE * MAX_ALPHABET_SIZE + 16
"""A block's shared memory beside its carry: the substitution table and one score's aligned slot."""

comptime DEFAULT_LEAF_CELLS = 4096
"""Cells from which a linear-space frame is split rather than solved over its whole matrix."""
comptime PARALLEL_LEAF_FLOOR = 256
"""Leaves of a level from which the host solves them over its threads: fewer cost less than a fork."""
comptime DEVICE_STORED_CELLS = 1_000_000
"""Cells past which one pair goes to the linear-space engine rather than a batch of one: one warp sweeps a
batch pair, the whole device a linear-space one, and they took as long at about a million cells."""


@fieldwise_init
struct Space(Equatable, ImplicitlyCopyable, TrivialRegisterPassable):
    """Which engine a pair's height allows: a batch strip, its carry in shared memory, or tiles."""

    var kind: UInt8

    comptime BANDED = Self(0)
    comptime TILED = Self(1)


def band_length(specs: GpuSpecs) -> Int:
    """The longest first sequence a batch strip's carry holds on this device."""
    var usable = specs.shared_memory_per_multiprocessor - specs.reserved_memory_per_block
    return (usable - STATIC_SHARED_BYTES) // (CARRY_LAYERS * size_of[Scalar[ScoreDType]]()) - 1


def serving_space(rows: Int, band: Int) -> Space:
    """The engine a pair of `rows` first-sequence letters takes, `band` the longest a strip carries."""
    return Space.BANDED if rows <= band else Space.TILED


def target_tiles(specs: GpuSpecs) -> Int:
    """The warps the device runs at once, which a level's tiles aim to fill."""
    return max(specs.streaming_multiprocessors * specs.max_blocks_per_multiprocessor, 1)


def launch_bytes(rows: Int, columns: Int) -> Int:
    """What a `rows` by `columns` pair's decisions take in `device_alignments`."""
    return ceildiv(columns, STRIP_WIDTH) * (rows + STRIP_LANES) * STRIP_LANES * size_of[Scalar[ChangeDType]]()


# endregion Sizes

# region Shared steps


@always_inline
def ahead(score: Int32, place: Int64, best: Int32, best_place: Int64) -> Bool:
    """Whether a cut scoring `score` at `place` displaces the best so far: a higher score, or an equal one at
    an earlier place, which is `gotoh.best_cut`'s first maximum."""
    return score > best or (score == best and place < best_place)


@always_inline
def keeps[every_tie: Bool](score: Int32, place: Int64, best: Int32, best_place: Int64) -> Bool:
    """Whether `(score, place)` displaces the best so far: by `ahead` when places break `every_tie`, a cut's
    rule, else by `beats`, a local alignment's end's."""
    comptime if every_tie:
        return ahead(score, place, best, best_place)
    return beats(score, place, best, best_place)


@always_inline
def warp_best[every_tie: Bool](mut best: Int32, mut place: Int64):
    """The warp's best of every lane's `(best, place)`, by `keeps`, into every lane."""
    var span = UInt32(1)
    while span < UInt32(WARP_SIZE):
        var theirs = shuffle_xor(best, span)
        var their_place = shuffle_xor(place, span)
        if keeps[every_tie](theirs, their_place, best, place):
            best = theirs
            place = their_place
        span *= 2


@always_inline
def block_best[
    every_tie: Bool
](
    scores: Pointer[Scalar[ScoreDType], MutUntrackedOrigin, address_space=AddressSpace.SHARED],
    places: Pointer[Scalar[DType.int64], MutUntrackedOrigin, address_space=AddressSpace.SHARED],
    best: Int32,
    place: Int64,
):
    """The block's best of every thread's `(best, place)`, by `keeps`, into `scores[0]` and `places[0]`: each
    warp's best through registers, then thread zero's pass over the warps'."""
    var winner = best
    var winner_place = place
    warp_best[every_tie](winner, winner_place)
    var warp = Int(thread_idx.x) // WARP_SIZE
    if Int(lane_id()) == 0:
        scores[unsafe_offset=warp] = winner
        places[unsafe_offset=warp] = winner_place
    barrier()
    if thread_idx.x != 0:
        return
    for other in range(1, WARPS_PER_BLOCK):
        if keeps[every_tie](scores[unsafe_offset=other], places[unsafe_offset=other], winner, winner_place):
            winner = scores[unsafe_offset=other]
            winner_place = places[unsafe_offset=other]
    scores[unsafe_offset=0] = winner
    places[unsafe_offset=0] = winner_place


@always_inline
def staged_table(
    substitutions: Pointer[Scalar[SubstitutionDType], MutAnyOrigin], alphabet_size: Int
) -> Pointer[Scalar[SubstitutionDType], MutAnyOrigin, address_space=AddressSpace.SHARED]:
    """The substitution table copied into shared memory by the whole block; a barrier must follow before it
    is read."""
    var table = stack_allocation[
        MAX_ALPHABET_SIZE * MAX_ALPHABET_SIZE, Scalar[SubstitutionDType], address_space=AddressSpace.SHARED
    ]()
    for index in range(Int(thread_idx.x), alphabet_size * alphabet_size, Int(block_dim.x)):
        table[unsafe_offset=index] = substitutions[unsafe_offset=index]
    return table.unsafe_origin_cast[MutAnyOrigin]()


@fieldwise_init
struct StripRow(ImplicitlyCopyable, TrivialRegisterPassable):
    """What a lane's pass over its columns of a row leaves: the score and insertion layer of its last column,
    which the lane to its right reads next, and with recording its cells' decisions, packed."""

    var score: Int32
    var insertion: Int32
    var packed: UInt32


@always_inline
def strip_row[
    mode: AlignmentMode, record: Bool
](
    table: Pointer[Scalar[SubstitutionDType], MutAnyOrigin, address_space=AddressSpace.SHARED],
    table_row: Int,
    symbols: Array[Int32, STRIP_COLUMNS],
    mut scores: Array[Int32, STRIP_COLUMNS],
    mut deletions: Array[Int32, STRIP_COLUMNS],
    above_left: Int32,
    left_score: Int32,
    left_insertion: Int32,
    gaps: AffineGapCosts,
    owned: Int,
    first_place: Int64,
    mut best: Int32,
    mut best_place: Int64,
) -> StripRow:
    """A lane's `STRIP_COLUMNS` cells of one row, left to right: `scores` and `deletions` hold the row above
    and are left holding this one; `table_row` is the row letter's offset in `table`; `first_place` is the
    first cell's place, row-major, for a local alignment's best of the `owned` cells that lie in the matrix."""
    var diagonal = above_left
    var running_score = left_score
    var running_insertion = left_insertion
    var packed = UInt32(0)
    comptime for slot in range(STRIP_COLUMNS):
        var pair = Int32(table[unsafe_offset=table_row + Int(symbols[slot])])
        var cell = gotoh_cell[mode](
            diagonal, scores[slot], deletions[slot], running_score, running_insertion, pair, gaps
        )
        comptime if record:
            var decision = decide[mode](
                cell,
                diagonal + pair,
                Cell(scores[slot], deletions[slot], 0),
                Cell(running_score, 0, running_insertion),
                gaps,
            )
            packed |= UInt32(decision.code) << UInt32(DECISION_BITS * slot)
        comptime if mode == AlignmentMode.LOCAL:
            if slot < owned and beats(cell.score, first_place + Int64(slot), best, best_place):
                best = cell.score
                best_place = first_place + Int64(slot)
        diagonal = scores[slot]
        scores[slot] = cell.score
        deletions[slot] = cell.deletion
        running_score = cell.score
        running_insertion = cell.insertion
    return StripRow(running_score, running_insertion, packed)


@always_inline
def global_border(length: Int, gaps: AffineGapCosts, entering: GapRun = GapRun.OPENS) -> Int32:
    """A global sweep's border cell `length` letters from the corner: a gap, opened already above when a
    deletion run enters open."""
    return Int32(length) * gaps.extend if entering == GapRun.EXTENDS else gaps.run(length)


# endregion Shared steps

# region Batches


def batch_kernel[
    mode: AlignmentMode, record: Bool
](
    sequences: Pointer[Scalar[SymbolDType], MutAnyOrigin],
    offsets: Pointer[Scalar[OffsetDType], MutAnyOrigin],
    substitutions: Pointer[Scalar[SubstitutionDType], MutAnyOrigin],
    letters: Pointer[Scalar[SymbolDType], MutAnyOrigin],
    changes: Pointer[Scalar[ChangeDType], MutAnyOrigin],
    change_offsets: Pointer[Scalar[DType.int64], MutAnyOrigin],
    results: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    first_gapped: Pointer[Scalar[SymbolDType], MutAnyOrigin],
    second_gapped: Pointer[Scalar[SymbolDType], MutAnyOrigin],
    gapped_lengths: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    gapped_stride: Int32,
    carry_stride: Int32,
    alphabet_size: Int32,
    open: Int32,
    extend: Int32,
):
    """Pair `block_idx.x` of a batch, its score into `results`; with `record`, its decisions packed into
    `changes` and walked back by lane zero into its gapped rows, reversed.

    A word of decisions is stored at its strip's step rather than its row, so a warp's 32 stores of one step
    land side by side; the walk finds a cell's word from its row and column in closed form."""
    var pair = Int(block_idx.x)
    var first_start = Int(offsets[unsafe_offset=2 * pair])
    var second_start = Int(offsets[unsafe_offset=2 * pair + 1])
    var rows = second_start - first_start
    var columns = Int(offsets[unsafe_offset=2 * pair + 2]) - second_start
    var size = Int(alphabet_size)
    var gaps = AffineGapCosts(open, extend)
    var lane = Int(lane_id())
    var stride = Int64(columns + 1)
    var strip_words = (rows + STRIP_LANES) * STRIP_LANES
    var decisions_at = 0
    comptime if record:
        decisions_at = Int(change_offsets[unsafe_offset=pair])

    var table = staged_table(substitutions, size)
    # The last cell's lane is not lane zero, which reports, so the score crosses the warp in shared memory.
    var corner = stack_allocation[1, Scalar[ScoreDType], address_space=AddressSpace.SHARED]()
    if thread_idx.x == 0:
        corner[unsafe_offset=0] = gaps.run(rows + columns) if mode == AlignmentMode.GLOBAL else 0
    var carry = external_memory[Scalar[ScoreDType], address_space=AddressSpace.SHARED, alignment=16, name="carry"]()
    var carry_scores = carry
    var carry_insertions = carry.unsafe_offset(Int(carry_stride))
    for row in range(Int(thread_idx.x), rows + 1, Int(block_dim.x)):
        var border = global_border(row, gaps) if mode == AlignmentMode.GLOBAL else Int32(0)
        carry_scores[unsafe_offset=row] = border
        carry_insertions[unsafe_offset=row] = border + gaps.open + gaps.extend
    barrier()

    var best = Int32(0)
    var best_place = Int64(0)
    for strip in range(ceildiv(columns, STRIP_WIDTH)):
        var first_column = strip * STRIP_WIDTH + lane * STRIP_COLUMNS
        var owned = clamp(columns - first_column, 0, STRIP_COLUMNS)
        var symbols = Array[Int32, STRIP_COLUMNS](fill=0)
        var scores = Array[Int32, STRIP_COLUMNS](fill=0)
        var deletions = Array[Int32, STRIP_COLUMNS](fill=0)
        comptime for slot in range(STRIP_COLUMNS):
            # A lane past the last column reads the last one's letter, its cells never counted.
            symbols[slot] = Int32(sequences[unsafe_offset=second_start + max(min(first_column + slot, columns - 1), 0)])
            comptime if mode == AlignmentMode.GLOBAL:
                scores[slot] = gaps.run(first_column + slot + 1)
            deletions[slot] = scores[slot] + gaps.open + gaps.extend
        var edge_score = gaps.run(first_column + STRIP_COLUMNS) if mode == AlignmentMode.GLOBAL else Int32(0)
        var edge_insertion = edge_score + gaps.open + gaps.extend
        var above_left = gaps.run(first_column) if mode == AlignmentMode.GLOBAL else Int32(0)

        for step in range(1, rows + STRIP_LANES + 1):
            var row = step - lane
            var left_score = shuffle_up(edge_score, 1)
            var left_insertion = shuffle_up(edge_insertion, 1)
            if lane == 0:
                left_score = carry_scores[unsafe_offset=clamp(row, 0, rows)]
                left_insertion = carry_insertions[unsafe_offset=clamp(row, 0, rows)]
            var diagonal = above_left
            above_left = left_score
            if row < 1 or row > rows:
                continue
            var done = strip_row[mode, record](
                table,
                Int(sequences[unsafe_offset=first_start + row - 1]) * size,
                symbols,
                scores,
                deletions,
                diagonal,
                left_score,
                left_insertion,
                gaps,
                owned,
                Int64(row) * stride + Int64(first_column + 1),
                best,
                best_place,
            )
            comptime if record:
                changes[unsafe_offset=decisions_at + strip * strip_words + step * STRIP_LANES + lane] = done.packed
            edge_score = done.score
            edge_insertion = done.insertion
            if lane == STRIP_LANES - 1 and owned > 0:
                carry_scores[unsafe_offset=row] = scores[STRIP_COLUMNS - 1]
                carry_insertions[unsafe_offset=row] = done.insertion
            comptime if mode == AlignmentMode.GLOBAL:
                if row == rows and owned > 0 and first_column + owned == columns:
                    corner[unsafe_offset=0] = scores[owned - 1]
        barrier()

    comptime if mode == AlignmentMode.LOCAL:
        warp_best[False](best, best_place)
    if thread_idx.x != 0:
        return

    var score = corner[unsafe_offset=0]
    var row = rows
    var column = columns
    comptime if mode == AlignmentMode.LOCAL:
        score = best
        row = Int(best_place // stride) if best != 0 else 0
        column = Int(best_place % stride) if best != 0 else 0
    results[unsafe_offset=pair] = score
    comptime if not record:
        return

    var gap = Scalar[SymbolDType](GAP_BYTE)
    var written = 0
    var at = pair * Int(gapped_stride)

    @always_inline
    def emit(
        top: Scalar[SymbolDType], bottom: Scalar[SymbolDType]
    ) {mut written, imm at, imm first_gapped, imm second_gapped}:
        first_gapped[unsafe_offset=at + written] = top
        second_gapped[unsafe_offset=at + written] = bottom
        written += 1

    @always_inline
    def letter(start: Int, index: Int) {imm letters, imm sequences} -> Scalar[SymbolDType]:
        return letters[unsafe_offset=Int(sequences[unsafe_offset=start + index])]

    var state = Layer.ALIGNING
    while row > 0 and column > 0:
        var offset = column - 1
        var owner = (offset % STRIP_WIDTH) // STRIP_COLUMNS
        var word = changes[
            unsafe_offset=decisions_at + (offset // STRIP_WIDTH) * strip_words + (row + owner) * STRIP_LANES + owner
        ]
        var decision = Decision(UInt8((word >> UInt32(DECISION_BITS * (offset % STRIP_COLUMNS))) & 0x0F))
        if mode == AlignmentMode.LOCAL and state == Layer.ALIGNING and decision.ends():
            break
        var step = advance(state, decision)
        emit(
            letter(first_start, row - 1) if step.row_advance != 0 else gap,
            letter(second_start, column - 1) if step.column_advance != 0 else gap,
        )
        row += step.row_advance
        column += step.column_advance
        state = step.lands_in
    comptime if mode == AlignmentMode.GLOBAL:
        for rest in range(row, 0, -1):
            emit(letter(first_start, rest - 1), gap)
        for rest in range(column, 0, -1):
            emit(gap, letter(second_start, rest - 1))
    gapped_lengths[unsafe_offset=pair] = Int32(written)


struct BatchLaunch(Movable):
    """A batch's sequences and table on the device, and what its launch is sized by."""

    var pairs: Int
    var sequences: DeviceBuffer[SymbolDType]
    var offsets: DeviceBuffer[OffsetDType]
    var substitutions: DeviceBuffer[SubstitutionDType]
    var carry_stride: Int
    """A carry layer's length: one more than the batch's longest first sequence."""

    def __init__(
        out self,
        scope: DeviceScope,
        sequences: ImmSpan[Scalar[SymbolDType], _],
        offsets: List[Scalar[OffsetDType]],
        substitutions: ImmSpan[Scalar[SubstitutionDType], _],
        alphabet_size: Int,
    ) raises:
        """The batch whose pair `i` is `sequences[offsets[2i] ..< offsets[2i + 1]]` against what follows up to
        `offsets[2i + 2]`, copied to the device; refused past `MAX_ALPHABET_SIZE` letters."""
        if alphabet_size > MAX_ALPHABET_SIZE:
            raise AlignmentError(ErrorKind.ALPHABET_TOO_LARGE, "substitution table")
        self.pairs = (len(offsets) - 1) // 2
        var longest = 0
        for pair in range(self.pairs):
            longest = max(longest, Int(offsets[2 * pair + 1]) - Int(offsets[2 * pair]))
        self.carry_stride = longest + 1
        self.sequences = upload(scope, sequences)
        self.offsets = upload(scope, offsets)
        self.substitutions = upload(scope, substitutions)

    def shared_bytes(self) -> Int:
        """A block's carry, in bytes."""
        return CARRY_LAYERS * self.carry_stride * size_of[Scalar[ScoreDType]]()


def device_scores[
    mode: AlignmentMode
](
    scope: DeviceScope,
    sequences: ImmSpan[Scalar[SymbolDType], _],
    offsets: List[Scalar[OffsetDType]],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    gaps: AffineGapCosts,
) raises -> List[Int32]:
    """Every pair's best score, a warp each (see `BatchLaunch` for the batch's layout)."""
    var batch = BatchLaunch(scope, sequences, offsets, substitutions, alphabet_size)
    var results = zeroed[ScoreDType](scope, batch.pairs)
    # A scoring launch reads none of the recording's buffers; one placeholder stands in for each.
    var symbols = allocate[SymbolDType](scope, 1)
    var nowhere = symbols.unsafe_ptr().unsafe_origin_cast[MutAnyOrigin]()
    var changes = allocate[ChangeDType](scope, 1)
    var change_offsets = zeroed[DType.int64](scope, 1)
    var lengths = zeroed[ScoreDType](scope, 1)
    scope.context.enqueue_function[batch_kernel[mode, False]](
        batch.sequences.unsafe_ptr(),
        batch.offsets.unsafe_ptr(),
        batch.substitutions.unsafe_ptr(),
        nowhere,
        changes.unsafe_ptr(),
        change_offsets.unsafe_ptr(),
        results.unsafe_ptr(),
        nowhere,
        nowhere,
        lengths.unsafe_ptr(),
        Int32(0),
        Int32(batch.carry_stride),
        Int32(alphabet_size),
        gaps.open,
        gaps.extend,
        grid_dim=batch.pairs,
        block_dim=STRIP_LANES,
        shared_mem_bytes=batch.shared_bytes(),
        func_attribute=FuncAttribute.MAX_DYNAMIC_SHARED_SIZE_BYTES(UInt32(batch.shared_bytes())),
    )
    scope.context.synchronize()
    var out = List[Int32](capacity=batch.pairs)
    with results.map_to_host() as host:
        for pair in range(batch.pairs):
            out.append(host[pair])
    return out^


def device_alignments[
    mode: AlignmentMode
](
    scope: DeviceScope,
    sequences: ImmSpan[Scalar[SymbolDType], _],
    offsets: List[Scalar[OffsetDType]],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet: String,
    gaps: AffineGapCosts,
) raises -> List[GappedAlignment]:
    """Every pair's optimal alignment, a warp each, its decisions recorded and walked back on the device."""
    var alphabet_bytes = alphabet.as_bytes()
    var batch = BatchLaunch(scope, sequences, offsets, substitutions, len(alphabet_bytes))
    var change_offsets = List[Int64](capacity=batch.pairs + 1)
    var words = Int64(0)
    var widest = 1
    for pair in range(batch.pairs):
        var rows = Int(offsets[2 * pair + 1]) - Int(offsets[2 * pair])
        var columns = Int(offsets[2 * pair + 2]) - Int(offsets[2 * pair + 1])
        change_offsets.append(words)
        words += Int64(launch_bytes(rows, columns) // size_of[Scalar[ChangeDType]]())
        widest = max(widest, rows + columns)
    change_offsets.append(words)

    var letters = List[Scalar[SymbolDType]](capacity=len(alphabet_bytes))
    for byte in alphabet_bytes:
        letters.append(Scalar[SymbolDType](byte))
    var letters_buffer = upload(scope, letters)
    var changes = allocate[ChangeDType](scope, Int(words))
    var change_offsets_buffer = upload(scope, change_offsets)
    var results = zeroed[ScoreDType](scope, batch.pairs)
    var first_rows = allocate[SymbolDType](scope, batch.pairs * widest)
    var second_rows = allocate[SymbolDType](scope, batch.pairs * widest)
    var lengths = zeroed[ScoreDType](scope, batch.pairs)
    scope.context.enqueue_function[batch_kernel[mode, True]](
        batch.sequences.unsafe_ptr(),
        batch.offsets.unsafe_ptr(),
        batch.substitutions.unsafe_ptr(),
        letters_buffer.unsafe_ptr(),
        changes.unsafe_ptr(),
        change_offsets_buffer.unsafe_ptr(),
        results.unsafe_ptr(),
        first_rows.unsafe_ptr(),
        second_rows.unsafe_ptr(),
        lengths.unsafe_ptr(),
        Int32(widest),
        Int32(batch.carry_stride),
        Int32(len(alphabet_bytes)),
        gaps.open,
        gaps.extend,
        grid_dim=batch.pairs,
        block_dim=STRIP_LANES,
        shared_mem_bytes=batch.shared_bytes(),
        func_attribute=FuncAttribute.MAX_DYNAMIC_SHARED_SIZE_BYTES(UInt32(batch.shared_bytes())),
    )
    scope.context.synchronize()

    var aligned = List[GappedAlignment](capacity=batch.pairs)
    with results.map_to_host() as scores, lengths.map_to_host() as written:
        with first_rows.map_to_host() as first_host, second_rows.map_to_host() as second_host:
            for pair in range(batch.pairs):
                # The walk wrote each pair's columns last first.
                var top = List[UInt8](capacity=Int(written[pair]))
                var bottom = List[UInt8](capacity=Int(written[pair]))
                for column in range(pair * widest + Int(written[pair]) - 1, pair * widest - 1, -1):
                    top.append(UInt8(first_host[column]))
                    bottom.append(UInt8(second_host[column]))
                aligned.append(
                    GappedAlignment(scores[pair], String(unsafe_from_utf8=top), String(unsafe_from_utf8=bottom))
                )
    return aligned^


# endregion Batches

# region Tiles


@fieldwise_init
struct SweepPlan(ImplicitlyCopyable, TrivialRegisterPassable):
    """A half's sweep over a rectangle, as the tile kernel reads it: the rectangle, where its scratch lies,
    how it is tiled, the deletion run entering it, and which way it reads its letters. It crosses to the
    device as raw words."""

    var rows: Int64
    var columns: Int64
    var row_from: Int64
    """Its first row's place in the first sequence."""
    var column_from: Int64
    """Its first column's place on the batch's tape: the second sequence follows the first."""
    var left_base: Int64
    """Where its left-edge carry starts, a row each."""
    var top_base: Int64
    """Where its frontier starts, the last row's scores and deletions, a column each."""
    var corner_base: Int64
    """Where its tile corners start: three rotating rows of one per tile row."""
    var tile_rows: Int64
    var tile_columns: Int64
    var tile_height: Int64
    var entering: GapRun
    var half: SweepHalf


comptime PLAN_WORDS = ceildiv(size_of[SweepPlan](), size_of[Int64]())


def plan_sweep(
    rows: Int, columns: Int, row_from: Int, column_from: Int, entering: GapRun, half: SweepHalf
) -> SweepPlan:
    """A sweep of a `rows` by `columns` rectangle, its scratch and tiling still to be placed (see `plan_level`)."""
    return SweepPlan(Int64(rows), Int64(columns), Int64(row_from), Int64(column_from), 0, 0, 0, 1, 1, 1, entering, half)


def plan_level(mut sweeps: List[SweepPlan], specs: GpuSpecs) -> Tuple[Int, Int, Int, Int, Int]:
    """Places a level's sweeps side by side in the scratch and tiles each: tiles as tall as still fill the
    device, the level's sweeps together, and no taller than square in count once a sweep's columns alone
    set its widest anti-diagonal. Returns the widest anti-diagonal of tiles, the most anti-diagonals any
    sweep takes, and the corner, left and frontier entries the level uses."""
    var want = max(target_tiles(specs) // max(len(sweeps), 1), 1)
    var widest = 1
    var deepest = 0
    var corners = 0
    var lefts = 0
    var tops = 0
    for ref sweep in sweeps:
        var tile_columns = max(ceildiv(Int(sweep.columns), TILE_SIDE), 1)
        var height = clamp(ceildiv(Int(sweep.rows), min(want, tile_columns)), MIN_TILE_HEIGHT, TILE_SIDE)
        var tile_rows = max(ceildiv(Int(sweep.rows), height), 1)
        sweep.tile_columns = Int64(tile_columns)
        sweep.tile_height = Int64(height)
        sweep.tile_rows = Int64(tile_rows)
        sweep.corner_base = Int64(corners)
        sweep.left_base = Int64(lefts)
        sweep.top_base = Int64(tops)
        corners += 3 * (tile_rows + 1)
        lefts += Int(sweep.rows) + 2
        tops += Int(sweep.columns) + 2
        widest = max(widest, min(tile_rows, tile_columns))
        deepest = max(deepest, tile_rows + tile_columns - 1)
    return (widest, deepest, corners, lefts, tops)


def tile_kernel[
    mode: AlignmentMode
](
    sequences: Pointer[Scalar[SymbolDType], MutAnyOrigin],
    substitutions: Pointer[Scalar[SubstitutionDType], MutAnyOrigin],
    plans: Pointer[Scalar[DType.int64], MutAnyOrigin],
    block_scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    block_places: Pointer[Scalar[DType.int64], MutAnyOrigin],
    forward_scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    forward_deletes: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    reverse_scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    reverse_deletes: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    left_scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    left_inserts: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    corner_scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    tile_diagonal: Int32,
    widest: Int32,
    alphabet_size: Int32,
    open: Int32,
    extend: Int32,
):
    """The tile of anti-diagonal `tile_diagonal` that block `block_idx.x` holds: block `b` is tile
    `b % widest` of the diagonal in sweep `b // widest`, one flat index, a level holding more sweeps than a
    grid's second dimension allows.

    A tile reads its left column from the tile to its left (`left_*`), its top row from the one above (the
    half's frontier), and its corner from the one above and to the left (`corner_scores`, written two
    launches before); it writes its own right column, last row and bottom-right corner in turn, which the
    sweep's last row leaves as the half's frontier. A local sweep also keeps its best cell, block by block."""
    var widest_tiles = Int(widest)
    var block = Int(block_idx.x)
    var plan = plans.unsafe_bitcast[SweepPlan]()[unsafe_offset=block // widest_tiles]
    var rows = Int(plan.rows)
    var columns = Int(plan.columns)
    var tile_rows = Int(plan.tile_rows)
    var tile_columns = Int(plan.tile_columns)
    var diagonal = Int(tile_diagonal)
    var tile_row = diagonal - min(diagonal, tile_columns - 1) + block % widest_tiles
    if tile_row > min(diagonal, tile_rows - 1):
        return
    var tile_height = Int(plan.tile_height)
    var row_begin = tile_row * tile_height
    var column_begin = (diagonal - tile_row) * TILE_SIDE
    var height = min(tile_height, rows - row_begin)
    var span = min(TILE_SIDE, columns - column_begin)
    if height <= 0 or span <= 0:
        return

    comptime local = mode == AlignmentMode.LOCAL
    var gaps = AffineGapCosts(open, extend)
    var backward = plan.half == SweepHalf.REVERSE
    var entering = plan.entering
    var row_from = Int(plan.row_from)
    var column_from = Int(plan.column_from)
    var left_base = Int(plan.left_base)
    var top_base = Int(plan.top_base)
    var corners = Int(plan.corner_base)
    var corner_row = tile_rows + 1
    var top_scores = reverse_scores if backward else forward_scores
    var top_deletes = reverse_deletes if backward else forward_deletes
    var lane = Int(lane_id())
    var first_column = lane * STRIP_COLUMNS
    var owned = clamp(span - first_column, 0, STRIP_COLUMNS)
    var size = Int(alphabet_size)

    @always_inline
    def border(row: Int) {imm gaps, imm entering} -> Int32:
        """The left border's score at `row`: zero for a local sweep."""
        comptime if local:
            return 0
        return global_border(row, gaps, entering)

    @always_inline
    def top_border(column: Int) {imm gaps} -> Int32:
        """The top border's score at `column`: zero for a local sweep."""
        comptime if local:
            return 0
        return gaps.run(column)

    var table = staged_table(substitutions, size)
    # The left column, staged by the whole warp: lane zero reads one entry a step, on the path every shuffle
    # waits for, and this tile's own right column, written as it sweeps, shares its slots with it.
    var edge_scores = stack_allocation[TILE_SIDE, Scalar[ScoreDType], address_space=AddressSpace.SHARED]()
    var edge_inserts = stack_allocation[TILE_SIDE, Scalar[ScoreDType], address_space=AddressSpace.SHARED]()
    for index in range(Int(thread_idx.x), height, Int(block_dim.x)):
        var row = row_begin + index + 1
        if column_begin == 0:
            edge_scores[unsafe_offset=index] = border(row)
            edge_inserts[unsafe_offset=index] = border(row) + gaps.open + gaps.extend
        else:
            edge_scores[unsafe_offset=index] = left_scores[unsafe_offset=left_base + row]
            edge_inserts[unsafe_offset=index] = left_inserts[unsafe_offset=left_base + row]
    barrier()

    # The row above, for this lane's columns: the border on the first tile row, else the frontier.
    var symbols = Array[Int32, STRIP_COLUMNS](fill=0)
    var scores = Array[Int32, STRIP_COLUMNS](fill=0)
    var deletions = Array[Int32, STRIP_COLUMNS](fill=0)
    comptime for slot in range(STRIP_COLUMNS):
        # A lane past the tile's last column repeats it, its cells never counted.
        var column = column_begin + min(first_column + slot, span - 1) + 1
        symbols[slot] = Int32(sequences[unsafe_offset=column_from + (columns - column if backward else column - 1)])
        if row_begin == 0:
            scores[slot] = top_border(column)
            deletions[slot] = scores[slot] + gaps.open + gaps.extend
        else:
            scores[slot] = top_scores[unsafe_offset=top_base + column]
            deletions[slot] = top_deletes[unsafe_offset=top_base + column]
    var edge_score = scores[STRIP_COLUMNS - 1]
    var edge_insertion = edge_score + gaps.open + gaps.extend

    # The cell above and left of this lane's first: lane zero's is the tile's corner, every other lane's a
    # cell of the row above.
    var above_left: Int32
    if lane != 0:
        var column = column_begin + min(first_column, span - 1)
        above_left = top_border(column) if row_begin == 0 else top_scores[unsafe_offset=top_base + column]
    elif row_begin == 0:
        above_left = top_border(column_begin)
    elif column_begin == 0:
        above_left = border(row_begin)
    else:
        above_left = corner_scores[unsafe_offset=corners + (diagonal % 3) * corner_row + tile_row]

    var best = Int32(0)
    var best_place = Int64(0)
    for step in range(1, height + STRIP_LANES + 1):
        var local_row = step - lane
        var left_score = shuffle_up(edge_score, 1)
        var left_insertion = shuffle_up(edge_insertion, 1)
        if lane == 0:
            left_score = edge_scores[unsafe_offset=clamp(local_row, 1, height) - 1]
            left_insertion = edge_inserts[unsafe_offset=clamp(local_row, 1, height) - 1]
        var diagonal_score = above_left
        above_left = left_score
        if local_row < 1 or local_row > height:
            continue
        var row = row_begin + local_row
        var letter = Int(sequences[unsafe_offset=row_from + (rows - row if backward else row - 1)])
        var done = strip_row[mode, False](
            table,
            letter * size,
            symbols,
            scores,
            deletions,
            diagonal_score,
            left_score,
            left_insertion,
            gaps,
            owned,
            Int64(row) * Int64(columns + 1) + Int64(column_begin + first_column + 1),
            best,
            best_place,
        )
        edge_score = done.score
        edge_insertion = done.insertion
        # A tile short of the full side lies against the matrix's right edge: no tile reads its right column.
        var full = span == TILE_SIDE
        if lane == STRIP_LANES - 1 and full:
            left_scores[unsafe_offset=left_base + row] = done.score
            left_inserts[unsafe_offset=left_base + row] = done.insertion
        if local_row != height:
            continue
        # The tile's last row onto the frontier; column zero, the border, which the cut reads too, from the
        # tile on the left edge alone.
        if lane == 0 and column_begin == 0:
            top_scores[unsafe_offset=top_base] = border(row)
            top_deletes[unsafe_offset=top_base] = border(row) if not local else gaps.open + gaps.extend
        comptime for slot in range(STRIP_COLUMNS):
            if slot < owned:
                var column = column_begin + first_column + slot + 1
                top_scores[unsafe_offset=top_base + column] = scores[slot]
                top_deletes[unsafe_offset=top_base + column] = deletions[slot]
        if lane == STRIP_LANES - 1 and full:
            corner_scores[unsafe_offset=corners + ((diagonal + 2) % 3) * corner_row + tile_row + 1] = scores[
                STRIP_COLUMNS - 1
            ]

    comptime if local:
        warp_best[False](best, best_place)
        if thread_idx.x == 0 and beats(
            best, best_place, block_scores[unsafe_offset=block], block_places[unsafe_offset=block]
        ):
            block_scores[unsafe_offset=block] = best
            block_places[unsafe_offset=block] = best_place


struct TileScratch(Movable):
    """What one alignment's tiled sweeps share on the device, sized for its widest level: the sequences and
    table, each half's frontier, the left-edge carries, the tile corners, each block's best cell for a
    local scan, and each split's cut."""

    var sequences: DeviceBuffer[SymbolDType]
    var substitutions: DeviceBuffer[SubstitutionDType]
    var forward_scores: DeviceBuffer[ScoreDType]
    var forward_deletes: DeviceBuffer[ScoreDType]
    var reverse_scores: DeviceBuffer[ScoreDType]
    var reverse_deletes: DeviceBuffer[ScoreDType]
    var left_scores: DeviceBuffer[ScoreDType]
    var left_inserts: DeviceBuffer[ScoreDType]
    var corner_scores: DeviceBuffer[ScoreDType]
    var block_scores: DeviceBuffer[ScoreDType]
    var block_places: DeviceBuffer[DType.int64]
    var cuts: DeviceBuffer[ScoreDType]
    var frontier_span: Int
    var left_span: Int
    var corner_span: Int
    var block_slots: Int

    def __init__(
        out self,
        scope: DeviceScope,
        first: ImmSpan[Scalar[SymbolDType], _],
        second: ImmSpan[Scalar[SymbolDType], _],
        substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    ) raises:
        """Scratch for aligning `first` against `second`. A level's frames split the rows and each frame's two
        halves span its columns, so a level holds the rows once, the columns twice, and a few entries more a
        sweep, never more than two sweeps a row."""
        # The kernel stages the whole table into shared memory unchecked, so it is refused here.
        if len(substitutions) > MAX_ALPHABET_SIZE * MAX_ALPHABET_SIZE:
            raise AlignmentError(ErrorKind.ALPHABET_TOO_LARGE, "staged substitution table")
        var rows = len(first)
        var columns = len(second)
        self.frontier_span = 2 * columns + 4 * rows + 32
        self.left_span = 5 * rows + 32
        self.corner_span = 3 * ceildiv(rows, MIN_TILE_HEIGHT) + 6 * rows + 32
        # A local scan is one sweep, its grid one anti-diagonal of tiles wide.
        self.block_slots = min(ceildiv(rows, MIN_TILE_HEIGHT), ceildiv(columns, TILE_SIDE) + 2) + 4
        var tape = List[Scalar[SymbolDType]](capacity=rows + columns)
        tape.extend(first)
        tape.extend(second)
        self.sequences = upload(scope, Span(tape))
        self.substitutions = upload(scope, substitutions)
        self.forward_scores = allocate[ScoreDType](scope, self.frontier_span)
        self.forward_deletes = allocate[ScoreDType](scope, self.frontier_span)
        self.reverse_scores = allocate[ScoreDType](scope, self.frontier_span)
        self.reverse_deletes = allocate[ScoreDType](scope, self.frontier_span)
        self.left_scores = allocate[ScoreDType](scope, self.left_span)
        self.left_inserts = allocate[ScoreDType](scope, self.left_span)
        self.corner_scores = allocate[ScoreDType](scope, self.corner_span)
        self.block_scores = zeroed[ScoreDType](scope, self.block_slots)
        self.block_places = zeroed[DType.int64](scope, self.block_slots)
        self.cuts = allocate[ScoreDType](scope, 4 * max(rows, 1))


def sweep_level[
    mode: AlignmentMode
](
    scope: DeviceScope, mut scratch: TileScratch, mut sweeps: List[SweepPlan], alphabet_size: Int, gaps: AffineGapCosts
) raises:
    """Sweeps a level's rectangles together, each half's last row left on its frontier at `top_base`."""
    if len(sweeps) == 0:
        return
    var widest, deepest, corners, lefts, tops = plan_level(sweeps, scope.specs)
    if corners > scratch.corner_span or lefts > scratch.left_span or tops > scratch.frontier_span:
        raise AlignmentError(
            ErrorKind.SCRATCH_TOO_SMALL,
            String("corners ", corners, " of ", scratch.corner_span, ", left ", lefts, " of ", scratch.left_span),
        )
    var words = List[Int64](length=len(sweeps) * PLAN_WORDS, fill=0)
    for index in range(len(sweeps)):
        words.unsafe_ptr().unsafe_offset(index * PLAN_WORDS).unsafe_bitcast[SweepPlan]()[] = sweeps[index]
    var plans = upload(scope, Span(words))
    for tile_diagonal in range(deepest):
        scope.context.enqueue_function[tile_kernel[mode]](
            scratch.sequences.unsafe_ptr(),
            scratch.substitutions.unsafe_ptr(),
            plans.unsafe_ptr(),
            scratch.block_scores.unsafe_ptr(),
            scratch.block_places.unsafe_ptr(),
            scratch.forward_scores.unsafe_ptr(),
            scratch.forward_deletes.unsafe_ptr(),
            scratch.reverse_scores.unsafe_ptr(),
            scratch.reverse_deletes.unsafe_ptr(),
            scratch.left_scores.unsafe_ptr(),
            scratch.left_inserts.unsafe_ptr(),
            scratch.corner_scores.unsafe_ptr(),
            Int32(tile_diagonal),
            Int32(widest),
            Int32(alphabet_size),
            gaps.open,
            gaps.extend,
            grid_dim=widest * len(sweeps),
            block_dim=STRIP_LANES,
        )
    scope.context.synchronize()


def best_end[
    half: SweepHalf
](
    scope: DeviceScope,
    mut scratch: TileScratch,
    rows: Int,
    columns: Int,
    tape_rows: Int,
    alphabet_size: Int,
    gaps: AffineGapCosts,
) raises -> Tuple[Int, Int, Int32]:
    """The best local alignment's end in the first `rows` rows and `columns` columns, read forward, or
    backward for the `REVERSE` half, which from an end found forward gives how far back it starts: its row,
    its column and its score, the first by row then column on a tie. `tape_rows` is the first sequence's
    length, where the second starts on the tape."""
    var sweeps: List[SweepPlan] = [plan_sweep(rows, columns, 0, tape_rows, GapRun.OPENS, half)]
    # A tile outside the sweep writes nothing, so an earlier scan's bests are cleared first.
    scope.context.enqueue_memset(scratch.block_scores, Scalar[ScoreDType](0))
    scope.context.enqueue_memset(scratch.block_places, Scalar[DType.int64](0))
    sweep_level[AlignmentMode.LOCAL](scope, scratch, sweeps, alphabet_size, gaps)
    var best = Int32(0)
    var best_place = Int64(0)
    var blocks = min(Int(sweeps[0].tile_rows), Int(sweeps[0].tile_columns), scratch.block_slots)
    with scratch.block_scores.map_to_host() as scores, scratch.block_places.map_to_host() as places:
        for block in range(blocks):
            if beats(scores[block], places[block], best, best_place):
                best = scores[block]
                best_place = places[block]
    var stride = Int64(columns + 1)
    return (Int(best_place // stride), Int(best_place % stride), best)


def cut_kernel(
    forward_scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    forward_deletes: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    reverse_scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    reverse_deletes: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    joins: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    cuts: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    refund: Int32,
):
    """Split `block_idx.x`'s cut, as `gotoh.best_cut` finds it, from its halves' frontiers at `joins[3s]` and
    `joins[3s + 1]`, its width `joins[3s + 2]`: the aligning layer's best and its column, then the straddling
    deletion's, into `cuts[4s ..< 4s + 4]`. Four numbers a split cross back, not two frontiers."""
    var split = Int(block_idx.x)
    var forward_base = Int(joins[unsafe_offset=3 * split])
    var reverse_base = Int(joins[unsafe_offset=3 * split + 1])
    var width = Int(joins[unsafe_offset=3 * split + 2])
    var aligned = NEGATIVE_INFINITY
    var aligned_at = Int64(0)
    var straddling = NEGATIVE_INFINITY
    var straddling_at = Int64(0)
    for offset in range(Int(thread_idx.x), width + 1, Int(block_dim.x)):
        var near = forward_base + offset
        var far = reverse_base + width - offset
        var through = forward_scores[unsafe_offset=near] + reverse_scores[unsafe_offset=far]
        if through > aligned:
            aligned = through
            aligned_at = Int64(offset)
        var across = forward_deletes[unsafe_offset=near] + reverse_deletes[unsafe_offset=far] + refund
        if across > straddling:
            straddling = across
            straddling_at = Int64(offset)
    var scores = stack_allocation[THREADS_PER_BLOCK, Scalar[ScoreDType], address_space=AddressSpace.SHARED]()
    var places = stack_allocation[THREADS_PER_BLOCK, Scalar[DType.int64], address_space=AddressSpace.SHARED]()
    block_best[True](scores, places, aligned, aligned_at)
    if thread_idx.x == 0:
        cuts[unsafe_offset=4 * split] = scores[unsafe_offset=0]
        cuts[unsafe_offset=4 * split + 1] = Scalar[ScoreDType](places[unsafe_offset=0])
    barrier()
    block_best[True](scores, places, straddling, straddling_at)
    if thread_idx.x == 0:
        cuts[unsafe_offset=4 * split + 2] = scores[unsafe_offset=0]
        cuts[unsafe_offset=4 * split + 3] = Scalar[ScoreDType](places[unsafe_offset=0])


def linear_path_on_device(
    scope: DeviceScope,
    mut scratch: TileScratch,
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    window: Rectangle,
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    gaps: AffineGapCosts,
    leaf_cells: Int,
    placement: Placement,
    path: PathView,
) raises:
    """`gotoh.linear_path` with every half's sweep on the device, a level of frames at a time: the frames of
    a level share no row or column, so all their halves sweep together and their leaves are solved
    together, on the host's threads."""
    var tape_rows = len(first)
    path.leaves(window.row_to, window.column_to)
    var level: List[Frame] = [Frame(window, GapRun.OPENS, GapRun.OPENS)]
    while len(level) > 0:
        var splitting = List[Frame]()
        var leaves = List[Frame]()
        for frame in level:
            if frame.area.height() == 0:
                continue
            if frame.solved_outright(leaf_cells):
                leaves.append(frame)
            else:
                splitting.append(frame)

        def solve_leaf_on_host(index: Int) {imm}:
            solve_outright(first, second, leaves[index], substitutions, alphabet_size, gaps, path)

        spread(solve_leaf_on_host, len(leaves), placement.threads if len(leaves) >= PARALLEL_LEAF_FLOOR else 1)
        if len(splitting) == 0:
            return

        var sweeps = List[SweepPlan](capacity=2 * len(splitting))
        for frame in splitting:
            var area = frame.area
            var middle = area.middle()
            var columns_at = tape_rows + area.column_from
            sweeps.append(
                plan_sweep(
                    middle - area.row_from, area.width(), area.row_from, columns_at, frame.top, SweepHalf.FORWARD
                )
            )
            sweeps.append(
                plan_sweep(area.row_to - middle, area.width(), middle, columns_at, frame.bottom, SweepHalf.REVERSE)
            )
        sweep_level[AlignmentMode.GLOBAL](scope, scratch, sweeps, alphabet_size, gaps)

        var joins = List[Scalar[ScoreDType]](capacity=3 * len(splitting))
        for index in range(len(splitting)):
            joins.append(Scalar[ScoreDType](sweeps[2 * index].top_base))
            joins.append(Scalar[ScoreDType](sweeps[2 * index + 1].top_base))
            joins.append(Scalar[ScoreDType](splitting[index].area.width()))
        var joins_buffer = upload(scope, Span(joins))
        scope.context.enqueue_function[cut_kernel](
            scratch.forward_scores.unsafe_ptr(),
            scratch.forward_deletes.unsafe_ptr(),
            scratch.reverse_scores.unsafe_ptr(),
            scratch.reverse_deletes.unsafe_ptr(),
            joins_buffer.unsafe_ptr(),
            scratch.cuts.unsafe_ptr(),
            gaps.extend - gaps.open,
            grid_dim=len(splitting),
            block_dim=THREADS_PER_BLOCK,
        )
        scope.context.synchronize()

        var next_level = List[Frame](capacity=2 * len(splitting))
        with scratch.cuts.map_to_host() as cuts:
            for index in range(len(splitting)):
                var cut = Cut.choosing(
                    cuts[4 * index], Int(cuts[4 * index + 1]), cuts[4 * index + 2], Int(cuts[4 * index + 3])
                )
                var halves = splitting[index].halves(cut)
                next_level.append(halves[0])
                next_level.append(halves[1])
        level = next_level^


def device_score[
    mode: AlignmentMode
](
    scope: DeviceScope,
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    gaps: AffineGapCosts,
) raises -> Int32:
    """One pair's best score, its sweep tiled over the device, a pair of any size."""
    var rows = len(first)
    var columns = len(second)
    # No tile of an empty rectangle writes a frontier.
    if rows == 0 or columns == 0:
        return Int32(0) if mode == AlignmentMode.LOCAL else gaps.run(rows + columns)
    var scratch = TileScratch(scope, first, second, substitutions)
    comptime if mode == AlignmentMode.LOCAL:
        return best_end[SweepHalf.FORWARD](scope, scratch, rows, columns, rows, alphabet_size, gaps)[2]
    var sweeps: List[SweepPlan] = [plan_sweep(rows, columns, 0, rows, GapRun.OPENS, SweepHalf.FORWARD)]
    sweep_level[AlignmentMode.GLOBAL](scope, scratch, sweeps, alphabet_size, gaps)
    with scratch.forward_scores.map_to_host() as frontier:
        return frontier[Int(sweeps[0].top_base) + columns]


def device_align[
    mode: AlignmentMode
](
    scope: DeviceScope,
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    gaps: AffineGapCosts,
    alphabet: String,
    leaf_cells: Int,
    placement: Placement,
) raises -> GappedAlignment:
    """One pair's optimal alignment in linear space, its sweeps tiled over the device, a pair of any size: a
    global one over the whole matrix; a local one over the span its two scans find, its end forward and its
    start backward from there."""
    var rows = len(first)
    var columns = len(second)
    var path = RowPath(rows)
    var scratch = TileScratch(scope, first, second, substitutions)
    comptime if mode == AlignmentMode.GLOBAL:
        linear_path_on_device(
            scope,
            scratch,
            first,
            second,
            Rectangle(0, rows, 0, columns),
            substitutions,
            alphabet_size,
            gaps,
            leaf_cells,
            placement,
            path.view(),
        )
        var whole = path.gapped(first, second, alphabet, mode, 0, rows)
        return GappedAlignment(path.score(first, second, substitutions, alphabet_size, gaps), whole[0], whole[1])

    var end_row, end_column, score = best_end[SweepHalf.FORWARD](
        scope, scratch, rows, columns, rows, alphabet_size, gaps
    )
    # A best of zero is the empty alignment, at the forward scan's end.
    var start_row = end_row
    if score > 0:
        var back_rows, back_columns, _ = best_end[SweepHalf.REVERSE](
            scope, scratch, end_row, end_column, rows, alphabet_size, gaps
        )
        start_row = end_row - back_rows
        linear_path_on_device(
            scope,
            scratch,
            first,
            second,
            Rectangle(start_row, end_row, end_column - back_columns, end_column),
            substitutions,
            alphabet_size,
            gaps,
            leaf_cells,
            placement,
            path.view(),
        )
    var span = path.gapped(first, second, alphabet, mode, start_row, end_row)
    return GappedAlignment(score, span[0], span[1])


# endregion Tiles
