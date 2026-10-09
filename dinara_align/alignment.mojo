# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Derived from AffineGaps (https://github.com/unum-science/AffineGaps), Copyright Ash Vardanian, under the
# Apache License, Version 2.0, and changed since: see LICENSES/Apache-2.0.txt and NOTICE.
"""
Gotoh affine-gap alignment on the GPU, with the reconstruction itself on the device; the host's is `gotoh`.

A row is the wrong sweep axis for a GPU, because the insertion term of a cell reads the
insertion term of its left neighbour. This module sweeps anti-diagonals instead, where every
cell of `d = i + j` reads only `d - 1` and `d - 2` and the whole above_left is independent.

Traceback is Hirschberg with a Myers-Miller affine join, splitting on rows rather than
anti-diagonals: a substitution step advances `i + j` by two and can skip an above-left cell entirely,
while the row index advances by zero or one per anti_diagonal. The recursion bottoms out in a direct
traceback over a stored decision tile.

The traceback walks all three layers — match, deletion and insertion — so every path realizes
the score reported alongside it, which is not automatic for affine gaps.
"""

from std.math import ceildiv, clamp
from std.memory import stack_allocation
from std.memory.pointer import AddressSpace
from std.sys.info import size_of

from max.gpu import WARP_SIZE, barrier, block_dim, block_idx, grid_dim, lane_id, thread_idx
from max.gpu.primitives.warp import shuffle_down, shuffle_up, shuffle_xor
from max.gpu.host import DeviceBuffer, FuncAttribute
from max.gpu.memory import external_memory

from .cigar import CigarWriter
from .errors import AlignmentError, ErrorKind
from .common import (
    spread,
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
    upload,
    zeroed,
)
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

# region Scoring

comptime ChangeDType = DType.uint32

comptime DECISION_BITS = 4
"""Bits one cell's traceback decision occupies inside a change word, which `CellDecision.nibble` writes."""

comptime CORNER_BYTES = 16
"""One aligned slot for the corner score the walk reads, which is all the strip stages beyond its table."""
comptime STATIC_SHARED_USED = MAX_ALPHABET_SIZE * MAX_ALPHABET_SIZE + CORNER_BYTES

comptime CARRY_LAYERS = 2
"""
A strip carries the score and insertion layers of the column to its left, one entry per row, which is what bounds how
long a first sequence one block can take.
"""


@fieldwise_init
struct Space(Equatable, ImplicitlyCopyable, TrivialRegisterPassable):
    """Which sweep serves a pair, which its height decides."""

    var identifier: UInt8
    """Which case this names."""
    comptime BANDED = Self(0)
    """One block per pair, carrying the column to its left in shared memory."""
    comptime TILED = Self(1)
    """Tiles over global-memory bands, bounded by nothing."""


def band_length(specs: GpuSpecs) -> Int:
    """How many rows one block's carry can index on this device.

    Derived from what the card reported rather than from the architecture the kernels were built
    for, so one artifact serves every target it is run on.
    """
    var usable = specs.shared_memory_per_multiprocessor - specs.reserved_memory_per_block
    return (usable - STATIC_SHARED_USED) // (CARRY_LAYERS * size_of[Scalar[ScoreDType]]()) - 1


def target_tiles_for(specs: GpuSpecs) -> Int:
    """Warps a level aims to put in flight, which is what this machine holds resident at once.

    Past roughly this many the strip stops gaining, so further splitting only pays the skew ramp
    again — and the number is the card's, not one card's product written down.
    """
    return max(specs.streaming_multiprocessors * specs.max_blocks_per_multiprocessor, 1)


def serving_space(rows: Int, band: Int) -> Space:
    """The one place a pair is measured against what a block can carry.

    The strip kernels index their carry by the first sequence, so height alone decides; the
    stored-traceback crossover is a throughput choice layered on top of this one.
    """
    return Space.BANDED if rows <= band else Space.TILED


# endregion Scoring


# region GPU Wavefront


@inline(.always)
def block_argmax(
    scores: Pointer[Scalar[ScoreDType], MutUntrackedOrigin, address_space=AddressSpace.SHARED],
    places: Pointer[Scalar[DType.int64], MutUntrackedOrigin, address_space=AddressSpace.SHARED],
    best: Int32,
    best_place: Int64,
):
    """Block-wide maximum keeping the earliest cell in row-major order on a tie.

    The stdlib block reductions carry a scalar, and this one has to carry the place alongside the
    score to break ties the way the reference scan's strict `>` does, so the tree stays here.
    Thread zero holds the winner afterwards.

    Two stages: each warp collapses through registers, then one barrier and a walk over the per-warp
    winners. The shuffle travels downward rather than by butterfly because the tie rule keeps the
    lower index, which only a directional reduction reproduces.
    """
    var lane = Int(thread_idx.x) % WARP_SIZE
    var winner = best
    var winner_place = best_place
    var reach = UInt32(WARP_SIZE // 2)
    while reach > 0:
        var theirs = shuffle_down(winner, reach)
        var their_place = shuffle_down(winner_place, reach)
        if beats(theirs, their_place, winner, winner_place):
            winner = theirs
            winner_place = their_place
        reach //= 2
    if lane == 0:
        scores[unsafe_offset=Int(thread_idx.x) // WARP_SIZE] = winner
        places[unsafe_offset=Int(thread_idx.x) // WARP_SIZE] = winner_place
    barrier()
    if Int(thread_idx.x) == 0:
        for index in range(1, WARPS_PER_BLOCK):
            var theirs = scores[unsafe_offset=index]
            var their_place = places[unsafe_offset=index]
            if beats(theirs, their_place, winner, winner_place):
                winner = theirs
                winner_place = their_place
        scores[unsafe_offset=0] = winner
        places[unsafe_offset=0] = winner_place


def device_scores[
    mode: AlignmentMode
](
    scope: DeviceScope,
    sequences: ImmSpan[Scalar[SymbolDType], _],
    offsets: List[Scalar[OffsetDType]],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    scoring: AffineGapCosts,
) raises -> List[Int32]:
    """Scores every pair in the batch, one thread block each."""
    var pairs = (len(offsets) - 1) // 2
    if alphabet_size > MAX_ALPHABET_SIZE:
        raise AlignmentError(ErrorKind.ALPHABET_TOO_LARGE, "substitution table")

    # The bands are indexed by row, so one launch only needs the longest first sequence it carries.
    var longest_first = 0
    for pair in range(pairs):
        longest_first = max(longest_first, Int(offsets[2 * pair + 1]) - Int(offsets[2 * pair]))
    var band_stride = longest_first + 1
    var dynamic_bytes = CARRY_LAYERS * band_stride * size_of[Scalar[ScoreDType]]()

    var sequences_buffer = upload(scope, sequences)
    var offsets_buffer = upload(scope, offsets)
    var substitutions_buffer = upload(scope, substitutions)
    var results_buffer = zeroed[ScoreDType](scope, pairs)
    # A discarding sweep never reads or writes these, but the one kernel still names them.
    var unused_symbols = allocate[SymbolDType](scope, 1)
    var unused_changes = allocate[ChangeDType](scope, 1)
    var unused_offsets = zeroed[DType.int64](scope, 1)
    var unused_lengths = zeroed[ScoreDType](scope, 1)
    # One placeholder fills every symbol slot. The origin cast is what lets it appear more than once in a launch, and it
    # is sound because a discarding sweep reads none of them.
    var nowhere = unused_symbols.unsafe_ptr().unsafe_origin_cast[MutAnyOrigin]()

    scope.context.enqueue_function[strip_pair_kernel[mode, Recording.DISCARDED]](
        sequences_buffer.unsafe_ptr(),
        offsets_buffer.unsafe_ptr(),
        substitutions_buffer.unsafe_ptr(),
        nowhere,
        unused_changes.unsafe_ptr(),
        unused_offsets.unsafe_ptr(),
        results_buffer.unsafe_ptr(),
        nowhere,
        nowhere,
        unused_lengths.unsafe_ptr(),
        Int32(0),
        Int32(band_stride),
        Int32(alphabet_size),
        scoring.open,
        scoring.extend,
        grid_dim=pairs,
        block_dim=STRIP_LANES,
        shared_mem_bytes=dynamic_bytes,
        func_attribute=FuncAttribute.MAX_DYNAMIC_SHARED_SIZE_BYTES(UInt32(dynamic_bytes)),
    )
    scope.context.synchronize()

    var results = List[Int32](capacity=pairs)
    with results_buffer.map_to_host() as host:
        for index in range(pairs):
            results.append(host[index])
    return results^


def strip_pair_kernel[
    mode: AlignmentMode, recording: Recording
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
    """A recording strip: every decision is packed as the cell is computed, then walked back.

    A lane packs its `STRIP_COLUMNS` decisions into one word and stores it under `anti_diagonal` rather
    than `row`. The 32 lanes of one anti-diagonal sit on 32 different rows but share the anti-diagonal, so the
    An anti-diagonal-major address makes the warp lay down 128 contiguous bytes where a row-major address
    would scatter the same warp over 32 sectors. The walk inverts it in closed form.

    Four bits hold a cell because `advance` reads exactly the source layer and the two run flags,
    and a clamped local cell borrows `Layer`'s spare source code instead of a fifth bit.
    """
    var pair = Int(block_idx.x)
    var first_start = Int(offsets[unsafe_offset=2 * pair])
    var second_start = Int(offsets[unsafe_offset=2 * pair + 1])
    var rows = second_start - first_start
    var columns = Int(offsets[unsafe_offset=2 * pair + 2]) - second_start
    var width = Int(alphabet_size)
    var scoring = AffineGapCosts(open, extend)
    var lane = Int(lane_id())
    var stride = columns + 1
    var tile = 0
    comptime if recording == Recording.TO_GLOBAL_MEMORY:
        tile = Int(change_offsets[unsafe_offset=pair])
    var strip_span = (rows + STRIP_LANES) * STRIP_LANES

    var table = stack_allocation[
        MAX_ALPHABET_SIZE * MAX_ALPHABET_SIZE, Scalar[SubstitutionDType], address_space=AddressSpace.SHARED
    ]()
    for index in range(Int(thread_idx.x), width * width, Int(block_dim.x)):
        table[unsafe_offset=index] = substitutions[unsafe_offset=index]

    # The bottom-right cell belongs to whichever lane owns the last column, and the walk runs on lane zero, so the
    # global answer crosses the warp through one shared word.
    var reported = stack_allocation[1, Scalar[ScoreDType], address_space=AddressSpace.SHARED]()
    if thread_idx.x == 0:
        var span = rows + columns
        comptime if mode == AlignmentMode.GLOBAL:
            reported[unsafe_offset=0] = scoring.run(span)
        else:
            reported[unsafe_offset=0] = 0

    var carry = external_memory[Scalar[ScoreDType], address_space=AddressSpace.SHARED, alignment=16, name="carry"]()
    var carry_scores = carry
    var carry_insertions = carry.unsafe_offset(Int(carry_stride))
    for index in range(Int(thread_idx.x), rows + 1, Int(block_dim.x)):
        var border = Int32(0)
        comptime if mode == AlignmentMode.GLOBAL:
            border = scoring.run(index)
        carry_scores[unsafe_offset=index] = border
        carry_insertions[unsafe_offset=index] = border + scoring.open + scoring.extend
    barrier()

    var best = Int32(0)
    var best_place = Int64(0)
    var strips = ceildiv(columns, STRIP_WIDTH)

    for strip in range(strips):
        var first_column = strip * STRIP_WIDTH + lane * STRIP_COLUMNS
        var owned = clamp(columns - first_column, 0, STRIP_COLUMNS)

        var symbols = Array[Int32, STRIP_COLUMNS](fill=0)
        comptime for column_slot in range(STRIP_COLUMNS):
            var column = min(first_column + column_slot, columns - 1)
            symbols[column_slot] = Int32(sequences[unsafe_offset=second_start + max(column, 0)])

        var scores = Array[Int32, STRIP_COLUMNS](fill=0)
        var deletions = Array[Int32, STRIP_COLUMNS](fill=0)
        comptime for column_slot in range(STRIP_COLUMNS):
            comptime if mode == AlignmentMode.GLOBAL:
                scores[column_slot] = scoring.run(first_column + column_slot + 1)
            deletions[column_slot] = scores[column_slot] + scoring.open + scoring.extend

        var past = first_column + STRIP_COLUMNS
        var edge_score = Int32(0)
        comptime if mode == AlignmentMode.GLOBAL:
            edge_score = scoring.run(past)
        var edge_insertion = edge_score + scoring.open + scoring.extend

        var above_left_carry = Int32(0)
        comptime if mode == AlignmentMode.GLOBAL:
            above_left_carry = scoring.run(first_column)

        for anti_diagonal in range(1, rows + STRIP_LANES + 1):
            var row = anti_diagonal - lane
            var left_score = shuffle_up(edge_score, 1)
            var left_insertion = shuffle_up(edge_insertion, 1)
            if lane == 0:
                var here = clamp(row, 0, rows)
                left_score = carry_scores[unsafe_offset=here]
                left_insertion = carry_insertions[unsafe_offset=here]
            var above_left_score = above_left_carry
            above_left_carry = left_score

            if row >= 1 and row <= rows:
                var symbol = Int(sequences[unsafe_offset=first_start + row - 1])
                var above_left = above_left_score
                var running_score = left_score
                var running_insertion = left_insertion
                var packed = UInt32(0)
                comptime for column_slot in range(STRIP_COLUMNS):
                    var substitution = Int32(table[unsafe_offset=symbol * width + Int(symbols[column_slot])])
                    var computed = gotoh_cell[mode](
                        above_left,
                        scores[column_slot],
                        deletions[column_slot],
                        running_score,
                        running_insertion,
                        substitution,
                        scoring,
                    )
                    comptime if recording == Recording.TO_GLOBAL_MEMORY:
                        # Both neighbours the decision needs are still in registers, unread.
                        var decision = decide[mode](
                            computed,
                            above_left + substitution,
                            Cell(scores[column_slot], deletions[column_slot], 0),
                            Cell(running_score, 0, running_insertion),
                            scoring,
                        )
                        packed |= UInt32(decision.code) << UInt32(4 * column_slot)
                    above_left = scores[column_slot]
                    scores[column_slot] = computed.score
                    deletions[column_slot] = computed.deletion
                    running_score = computed.score
                    running_insertion = computed.insertion
                    comptime if mode == AlignmentMode.LOCAL:
                        if column_slot < owned:
                            var place = Int64(row) * Int64(stride) + Int64(first_column + column_slot + 1)
                            if computed.score > best:
                                best = computed.score
                                best_place = place
                            elif computed.score == best and computed.score != 0 and place < best_place:
                                best_place = place
                comptime if recording == Recording.TO_GLOBAL_MEMORY:
                    changes[unsafe_offset=tile + strip * strip_span + anti_diagonal * STRIP_LANES + lane] = packed
                edge_score = running_score
                edge_insertion = running_insertion
                if lane == STRIP_LANES - 1 and owned > 0:
                    carry_scores[unsafe_offset=row] = scores[STRIP_COLUMNS - 1]
                    carry_insertions[unsafe_offset=row] = running_insertion
                comptime if mode == AlignmentMode.GLOBAL:
                    if row == rows and owned > 0 and first_column + owned == columns:
                        reported[unsafe_offset=0] = scores[owned - 1]
        barrier()

    # Warp-wide, keeping the earliest cell in row-major order on a tie, which is what the
    # reference scan's strict `>` picks. Every lane holds the winner afterwards.
    comptime if mode == AlignmentMode.LOCAL:
        var span = UInt32(1)
        while span < UInt32(STRIP_LANES):
            var theirs = shuffle_xor(best, span)
            var their_place = shuffle_xor(best_place, span)
            if beats(theirs, their_place, best, best_place):
                best = theirs
                best_place = their_place
            span *= 2

    if thread_idx.x != 0:
        return

    var start_row = rows
    var start_column = columns
    var final_score = reported[unsafe_offset=0]
    comptime if mode == AlignmentMode.LOCAL:
        final_score = best
        start_row = Int(best_place // Int64(stride))
        start_column = Int(best_place % Int64(stride))
        if final_score == 0:
            start_row = 0
            start_column = 0

    results[unsafe_offset=pair] = final_score
    comptime if recording == Recording.DISCARDED:
        return

    var row = start_row
    var column = start_column
    var produced = 0
    var base = pair * Int(gapped_stride)
    var state = Layer.ALIGNING

    while row > 0 and column > 0:
        var offset = column - 1
        var walk_lane = (offset % STRIP_WIDTH) // STRIP_COLUMNS
        var word = changes[
            unsafe_offset=tile + (offset // STRIP_WIDTH) * strip_span + (row + walk_lane) * STRIP_LANES + walk_lane
        ]
        var decision = Decision(UInt8((word >> UInt32(4 * (offset % STRIP_COLUMNS))) & 0x0F))
        if mode == AlignmentMode.LOCAL and state == Layer.ALIGNING and decision.ends():
            break
        var anti_diagonal = advance(state, decision)
        var gap = Scalar[SymbolDType](GAP_BYTE)
        first_gapped[unsafe_offset=base + produced] = (
            letters[unsafe_offset=Int(sequences[unsafe_offset=first_start + row - 1])] if anti_diagonal.row_advance
            != 0 else gap
        )
        second_gapped[unsafe_offset=base + produced] = (
            letters[
                unsafe_offset=Int(sequences[unsafe_offset=second_start + column - 1])
            ] if anti_diagonal.column_advance
            != 0 else gap
        )
        row += anti_diagonal.row_advance
        column += anti_diagonal.column_advance
        produced += 1
        state = anti_diagonal.lands_in

    # Only a global path is required to reach the origin; see the host reconstruction.
    if mode == AlignmentMode.GLOBAL:
        while row > 0:
            first_gapped[unsafe_offset=base + produced] = letters[
                unsafe_offset=Int(sequences[unsafe_offset=first_start + row - 1])
            ]
            second_gapped[unsafe_offset=base + produced] = Scalar[SymbolDType](GAP_BYTE)
            row -= 1
            produced += 1
        while column > 0:
            first_gapped[unsafe_offset=base + produced] = Scalar[SymbolDType](GAP_BYTE)
            second_gapped[unsafe_offset=base + produced] = letters[
                unsafe_offset=Int(sequences[unsafe_offset=second_start + column - 1])
            ]
            column -= 1
            produced += 1

    gapped_lengths[unsafe_offset=pair] = Int32(produced)


def launch_bytes(rows: Int, columns: Int) -> Int:
    """What a `rows` by `columns` pair's recorded decisions take in `device_alignments`, its strips' changes."""
    return ceildiv(columns, STRIP_WIDTH) * (rows + STRIP_LANES) * STRIP_LANES * size_of[Scalar[ChangeDType]]()


def device_alignments[
    mode: AlignmentMode
](
    scope: DeviceScope,
    sequences: ImmSpan[Scalar[SymbolDType], _],
    offsets: List[Scalar[OffsetDType]],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet: String,
    scoring: AffineGapCosts,
) raises -> List[GappedAlignment]:
    """Aligns a batch on the GPU by recording every decision, one thread block per pair."""
    var alphabet_bytes = alphabet.as_bytes()
    var alphabet_size = len(alphabet_bytes)
    if alphabet_size > MAX_ALPHABET_SIZE:
        raise AlignmentError(ErrorKind.ALPHABET_TOO_LARGE, "substitution table")
    var pairs = (len(offsets) - 1) // 2

    var change_offsets = List[Int64](capacity=pairs + 1)
    var running = Int64(0)
    var widest = 1
    var longest_first = 0
    for pair in range(pairs):
        var rows = Int(offsets[2 * pair + 1]) - Int(offsets[2 * pair])
        longest_first = max(longest_first, rows)
        var columns = Int(offsets[2 * pair + 2]) - Int(offsets[2 * pair + 1])
        change_offsets.append(running)
        var strips = ceildiv(columns, STRIP_WIDTH)
        running += Int64(strips) * Int64(rows + STRIP_LANES) * Int64(STRIP_LANES)
        widest = max(widest, rows + columns)
    change_offsets.append(running)

    var sequences_buffer = upload(scope, sequences)
    var offsets_buffer = upload(scope, offsets)
    var substitutions_buffer = upload(scope, substitutions)
    var letters = List[Scalar[SymbolDType]](capacity=alphabet_size)
    for index in range(alphabet_size):
        letters.append(Scalar[SymbolDType](alphabet_bytes[index]))
    var letters_buffer = upload(scope, letters)
    var changes_buffer = allocate[ChangeDType](scope, Int(max(running, Int64(1))))
    var change_offsets_buffer = upload(scope, change_offsets)
    var results_buffer = zeroed[ScoreDType](scope, pairs)
    var first_buffer = allocate[SymbolDType](scope, max(pairs * widest, 1))
    var second_buffer = allocate[SymbolDType](scope, max(pairs * widest, 1))
    var lengths_buffer = zeroed[ScoreDType](scope, pairs)

    scope.context.enqueue_function[strip_pair_kernel[mode, Recording.TO_GLOBAL_MEMORY]](
        sequences_buffer.unsafe_ptr(),
        offsets_buffer.unsafe_ptr(),
        substitutions_buffer.unsafe_ptr(),
        letters_buffer.unsafe_ptr(),
        changes_buffer.unsafe_ptr(),
        change_offsets_buffer.unsafe_ptr(),
        results_buffer.unsafe_ptr(),
        first_buffer.unsafe_ptr(),
        second_buffer.unsafe_ptr(),
        lengths_buffer.unsafe_ptr(),
        Int32(widest),
        Int32(longest_first + 1),
        Int32(alphabet_size),
        scoring.open,
        scoring.extend,
        grid_dim=pairs,
        block_dim=STRIP_LANES,
        shared_mem_bytes=CARRY_LAYERS * (longest_first + 1) * size_of[Scalar[ScoreDType]](),
        func_attribute=FuncAttribute.MAX_DYNAMIC_SHARED_SIZE_BYTES(
            UInt32(CARRY_LAYERS * (longest_first + 1) * size_of[Scalar[ScoreDType]]())
        ),
    )
    scope.context.synchronize()

    var aligned = List[GappedAlignment](capacity=pairs)
    with results_buffer.map_to_host() as scores_host, lengths_buffer.map_to_host() as lengths_host:
        with first_buffer.map_to_host() as first_host, second_buffer.map_to_host() as second_host:
            for pair in range(pairs):
                var produced = Int(lengths_host[pair])
                var base = pair * widest
                var first_gapped = List[UInt8](capacity=produced + 1)
                var second_gapped = List[UInt8](capacity=produced + 1)
                for index in range(produced - 1, -1, -1):
                    first_gapped.append(UInt8(first_host[base + index]))
                    second_gapped.append(UInt8(second_host[base + index]))
                aligned.append(
                    GappedAlignment(
                        scores_host[pair], String(unsafe_from_utf8=first_gapped), String(unsafe_from_utf8=second_gapped)
                    )
                )
    return aligned^


def device_hirschberg(
    scope: DeviceScope,
    mut buffers: SweepBuffers,
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    window: Rectangle,
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    scoring: AffineGapCosts,
    leaf_cells: Int,
    placement: Placement,
    path: PathView,
) raises:
    """`gotoh.linear_path` with every half's sweep on the device: a level of frames at a time, all of a
    level's sweeps launched together, as frames of a level share no row or column. Each launch is a
    barrier, which stands in for the grid-wide synchronization Mojo does not expose."""
    # The second sequence is staged after the first, so its indices carry that offset.
    var rows = len(first)
    path.leaves(window.row_to, window.column_to)
    var frames: List[Frame] = [Frame(window, GapRun.OPENS, GapRun.OPENS)]
    while len(frames) > 0:
        var splitting = List[Frame]()
        var leaves = List[Frame]()
        for frame in frames:
            if frame.area.height() == 0:
                continue
            if frame.solved_outright(leaf_cells):
                leaves.append(frame)
            else:
                splitting.append(frame)

        # A fork-join costs milliseconds, so a level of few leaves solves them in place.
        if len(leaves) >= PARALLEL_LEAF_FLOOR:

            def solve_leaf(slot: Int) {imm}:
                solve_outright(first, second, leaves[slot], substitutions, alphabet_size, scoring, path)

            spread(solve_leaf, len(leaves), placement.threads)
        else:
            for leaf in leaves:
                solve_outright(first, second, leaf, substitutions, alphabet_size, scoring, path)
        if len(splitting) == 0:
            break

        var sweeps = List[Sweep]()
        var joins = List[Scalar[ScoreDType]](length=len(splitting) * 3, fill=0)
        for index in range(len(splitting)):
            var frame = splitting[index]
            var area = frame.area
            var middle = area.middle()
            sweeps.append(
                Sweep(
                    middle - area.row_from,
                    area.width(),
                    area.row_from,
                    rows + area.column_from,
                    area.row_from,
                    area.column_from,
                    frame.top,
                    SweepHalf.FORWARD,
                )
            )
            sweeps.append(
                Sweep(
                    area.row_to - middle,
                    area.width(),
                    middle,
                    rows + area.column_from,
                    middle,
                    area.column_from,
                    frame.bottom,
                    SweepHalf.REVERSE,
                )
            )
            joins[index * 3 + 2] = Scalar[ScoreDType](area.width())

        device_sweep_level[AlignmentMode.GLOBAL](scope, buffers, sweeps, alphabet_size, scoring)

        for index in range(len(splitting)):
            joins[index * 3] = Scalar[ScoreDType](sweeps[index * 2].column_base)
            joins[index * 3 + 1] = Scalar[ScoreDType](sweeps[index * 2 + 1].column_base)
        var joins_buffer = upload(scope, Span(joins))
        scope.context.enqueue_function[crossing_kernel](
            buffers.top_scores.unsafe_ptr(),
            buffers.top_deletes.unsafe_ptr(),
            buffers.reverse_scores.unsafe_ptr(),
            buffers.reverse_deletes.unsafe_ptr(),
            joins_buffer.unsafe_ptr(),
            buffers.crossing.unsafe_ptr(),
            scoring.extend - scoring.open,
            grid_dim=len(splitting),
            block_dim=THREADS_PER_BLOCK,
        )
        scope.context.synchronize()

        var children = List[Frame]()
        with buffers.crossing.map_to_host() as host:
            for index in range(len(splitting)):
                var cut = Cut.choosing(
                    host[index * 4], Int(host[index * 4 + 1]), host[index * 4 + 2], Int(host[index * 4 + 3])
                )
                var halves = splitting[index].halves(cut)
                children.append(halves[0])
                children.append(halves[1])
        frames = children^


comptime DEFAULT_LEAF_CELLS = 4096
"""
Cells at or below which a Hirschberg frame is solved outright rather than split again. Raising it past a whole matrix
collapses the recursion to one direct traceback, which is what makes the two paths comparable under test.
"""

comptime PARALLEL_LEAF_FLOOR = 256
"""
A fork-join costs a few milliseconds, and the shallow levels of the recursion hold only a handful of leaves, so below
this many the dispatch costs more than the work it spreads.
"""

comptime DEVICE_STORED_CELLS = 1_000_000
"""
One pair on the stored device path gets a single warp, while the linear recursion spreads the same matrix over the
whole machine, so the device crossover sits far below what device memory would allow. Measured here: level near a
million cells, and the recursion is three times faster by sixteen million. A batch inverts the argument, since the
recursion takes its pairs in turn, and keeps to the caller's budget.
"""

comptime STRIP_COLUMNS = size_of[Scalar[ChangeDType]]() * 8 // DECISION_BITS
"""
Columns one lane owns, which is however many decisions a change word holds. A decision is a nibble, so the two are
locked together and widening `ChangeDType` widens the strip rather than silently truncating the traceback.
"""
comptime STRIP_LANES = WARP_SIZE
"""One lane per thread of a warp, so the skew shuffles reach their neighbour without crossing a warp boundary."""
comptime STRIP_WIDTH = STRIP_COLUMNS * STRIP_LANES


comptime TILE_SIDE = STRIP_WIDTH

comptime MIN_TILE_HEIGHT = STRIP_LANES
"""
Below this a tile spends more steps ramping the skew in and out than sweeping rows: a tile of `height` rows runs
`height + STRIP_LANES - 1` steps, so the ramp is the warp width rather than a number of its own.
"""


comptime PlanDType = DType.int64
"""
`rows`, `columns`, `first_from`, `second_from`, `row_base`, `column_base`, `corner_base`, `entering_run`, `half`,
`tile_rows_count`, `tile_columns_count`, `tile_height`. A device buffer is typed by `DType`, so the plan travels as
words and is read back
through this struct on both sides. Naming the fields once is what keeps the host writer and the two device readers
from disagreeing about a position seven hundred lines apart.
"""


@fieldwise_init
struct Recording(Equatable, ImplicitlyCopyable, TrivialRegisterPassable):
    """Where a sweep puts the decision it makes for each cell, if it records one at all.

    A sweep that only needs a score throws every decision away; one that has to reconstruct keeps
    them, and where they go is what separates a batch pair from a recursion leaf.
    """

    var identifier: UInt8
    """Which case this names."""
    comptime DISCARDED = Self(0)
    """The sweep keeps only scores, because the caller will recurse instead."""
    comptime TO_GLOBAL_MEMORY = Self(1)
    """The sweep stores one decision nibble per cell, eight to a word, for a direct traceback."""


@fieldwise_init
struct SweepPlan(ImplicitlyCopyable, TrivialRegisterPassable):
    """One sub-rectangle's assignment: what to sweep, and where its edges live."""

    var rows: Int64
    """How many rows this sub-rectangle covers."""
    var columns: Int64
    """How many columns it covers."""
    var first_from: Int64
    """First row of the sub-rectangle, in the full sequence."""
    var second_from: Int64
    """First column of the sub-rectangle, in the full sequence."""
    var row_base: Int64
    """Where its row-indexed scratch begins."""
    var column_base: Int64
    """Where its column-indexed scratch begins."""
    var corner_base: Int64
    """Where its tile-corner scratch begins."""
    var tile_rows_count: Int64
    """Tiles down."""
    var tile_columns_count: Int64
    """Tiles across."""
    var tile_height: Int64
    """Rows per tile."""
    var entering_run: GapRun
    """Whether a gap run is already open where the sweep starts."""
    var half: SweepHalf
    """Which direction the sweep runs."""


comptime PLAN_WORDS = size_of[SweepPlan]() // size_of[Int64]()


def tiled_sweep_kernel[
    mode: AlignmentMode
](
    sequences: Pointer[Scalar[SymbolDType], MutAnyOrigin],
    substitutions: Pointer[Scalar[SubstitutionDType], MutAnyOrigin],
    plans: Pointer[Scalar[PlanDType], MutAnyOrigin],
    block_best: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    block_place: Pointer[Scalar[DType.int64], MutAnyOrigin],
    forward_scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    forward_deletes: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    reverse_scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    reverse_deletes: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    left_scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    left_inserts: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    corner_scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    tile_anti_diagonal: Int32,
    widest_tiles: Int32,
    alphabet_size: Int32,
    open: Int32,
    extend: Int32,
):
    """One warp per tile, one tile-anti-diagonal per launch, every sweep of a level at once.

    Tiles sharing an anti-diagonal are independent, so a launch boundary supplies the barrier that
    Mojo 1.0 has no grid-wide primitive for. The flat block index carries both the sweep and the
    tile, which is what lets independent subproblems share a launch instead of taking turns. It is
    flat rather than two-dimensional because a deep recursion level holds more sweeps than the
    65535 a grid's second dimension allows.

    Inside a tile the sweep is a skewed register strip: lane `t` owns `STRIP_COLUMNS` columns and
    sits on row `s - t`, so the left-neighbour dependency travels by one `shuffle_up` and the tile
    runs to completion without a single barrier.
    """
    var lanes_wide = Int(widest_tiles)
    var sweep_index = Int(block_idx.x) // lanes_wide
    var tile_slot = Int(block_idx.x) % lanes_wide
    var sweep = plans.unsafe_bitcast[SweepPlan]()[unsafe_offset=sweep_index]
    var rows = Int(sweep.rows)
    var columns = Int(sweep.columns)
    var first_from = Int(sweep.first_from)
    var second_from = Int(sweep.second_from)
    var row_base = Int(sweep.row_base)
    var column_base = Int(sweep.column_base)
    var corner_base = Int(sweep.corner_base)
    var extends = sweep.entering_run == GapRun.EXTENDS
    var reversed_order = sweep.half == SweepHalf.REVERSE
    var tile_rows_count = Int(sweep.tile_rows_count)
    var tile_columns_count = Int(sweep.tile_columns_count)
    comptime local = mode == AlignmentMode.LOCAL
    var tile_height = Int(sweep.tile_height)
    # A corner is written on one tile-anti-diagonal and read on the next but one, by exactly one tile, so three rotating
    # buffers of one entry per tile row hold every live corner at once.
    var corner_rows = tile_rows_count + 1
    var slot = Int(block_idx.x)

    var tile_row_low = Int(tile_anti_diagonal) - min(Int(tile_anti_diagonal), tile_columns_count - 1)
    var tile_row_high = min(Int(tile_anti_diagonal), tile_rows_count - 1)
    var tile_row = tile_row_low + tile_slot
    if tile_row > tile_row_high:
        return
    var tile_column = Int(tile_anti_diagonal) - tile_row
    var row_begin = tile_row * tile_height
    var column_begin = tile_column * TILE_SIDE
    var height = min(tile_height, rows - row_begin)
    var width_span = min(TILE_SIDE, columns - column_begin)
    if height <= 0 or width_span <= 0:
        return

    var width = Int(alphabet_size)
    var scoring = AffineGapCosts(open, extend)
    var lane = Int(lane_id())
    var first_column = lane * STRIP_COLUMNS
    var owned = clamp(width_span - first_column, 0, STRIP_COLUMNS)
    var top_scores = reverse_scores if reversed_order else forward_scores
    var top_deletes = reverse_deletes if reversed_order else forward_deletes

    var table = stack_allocation[
        MAX_ALPHABET_SIZE * MAX_ALPHABET_SIZE, Scalar[SubstitutionDType], address_space=AddressSpace.SHARED
    ]()
    for index in range(Int(thread_idx.x), width * width, Int(block_dim.x)):
        table[unsafe_offset=index] = substitutions[unsafe_offset=index]

    # The tile's left column, staged once by the whole warp. Lane zero consumes one entry per anti-diagonal, and a
    # global load there would sit on the dependency chain that feeds every shuffle. Staging also decouples the read of
    # the neighbour's frontier from this tile's write of its own, which land in the same slots.
    var edge_scores = stack_allocation[TILE_SIDE, Scalar[ScoreDType], address_space=AddressSpace.SHARED]()
    var edge_inserts = stack_allocation[TILE_SIDE, Scalar[ScoreDType], address_space=AddressSpace.SHARED]()
    for index in range(Int(thread_idx.x), height, Int(block_dim.x)):
        var row = row_begin + index + 1
        if column_begin == 0:
            var border = Int32(0) if local else (Int32(row) * extend if extends else open + Int32(row - 1) * extend)
            edge_scores[unsafe_offset=index] = border
            edge_inserts[unsafe_offset=index] = border + open + extend
        else:
            edge_scores[unsafe_offset=index] = left_scores[unsafe_offset=row_base + row]
            edge_inserts[unsafe_offset=index] = left_inserts[unsafe_offset=row_base + row]
    barrier()

    # The tile's top row, for the columns this lane owns. Outside the matrix it is the affine ramp, or zero for a local
    # sweep; inside it is whatever the tile above left behind.
    var symbols = Array[Int32, STRIP_COLUMNS](fill=0)
    var scores = Array[Int32, STRIP_COLUMNS](fill=0)
    var deletions = Array[Int32, STRIP_COLUMNS](fill=0)
    comptime for column_slot in range(STRIP_COLUMNS):
        var column = column_begin + min(first_column + column_slot, width_span - 1) + 1
        var right_index = second_from + columns - column if reversed_order else second_from + column - 1
        symbols[column_slot] = Int32(sequences[unsafe_offset=right_index])
        if local and row_begin == 0:
            scores[column_slot] = 0
            deletions[column_slot] = open + extend
        elif row_begin == 0:
            scores[column_slot] = open + Int32(column - 1) * extend
            deletions[column_slot] = scores[column_slot] + open + extend
        else:
            scores[column_slot] = top_scores[unsafe_offset=column_base + column]
            deletions[column_slot] = top_deletes[unsafe_offset=column_base + column]

    var edge_score = scores[STRIP_COLUMNS - 1]
    var edge_insertion = edge_score + open + extend

    # The cell above and to the left of this lane's first one. For lane zero that is the tile's own corner; for every
    # other lane it is a cell of the top row.
    var above_left_carry = Int32(0)
    if lane == 0:
        if local and (row_begin == 0 or column_begin == 0):
            above_left_carry = 0
        elif row_begin == 0 and column_begin == 0:
            above_left_carry = 0
        elif row_begin == 0:
            above_left_carry = open + Int32(column_begin - 1) * extend
        elif column_begin == 0:
            above_left_carry = Int32(row_begin) * extend if extends else open + Int32(row_begin - 1) * extend
        else:
            above_left_carry = corner_scores[
                unsafe_offset=corner_base + (Int(tile_anti_diagonal) % 3) * corner_rows + tile_row
            ]
    else:
        # Clamped as the frontier reads above are: a lane owning no columns must not read past them.
        var column = column_begin + min(first_column, width_span - 1)
        if local and row_begin == 0:
            above_left_carry = 0
        elif row_begin == 0:
            above_left_carry = open + Int32(column - 1) * extend
        else:
            above_left_carry = top_scores[unsafe_offset=column_base + column]

    var best_cell = Int32(0)
    var best_at = Int64(0)

    for anti_diagonal in range(1, height + STRIP_LANES + 1):
        var local_row = anti_diagonal - lane
        var left_score = shuffle_up(edge_score, 1)
        var left_insertion = shuffle_up(edge_insertion, 1)
        if lane == 0:
            var index = min(max(local_row, 1), height) - 1
            left_score = edge_scores[unsafe_offset=index]
            left_insertion = edge_inserts[unsafe_offset=index]
        var above_left_score = above_left_carry
        above_left_carry = left_score

        if local_row >= 1 and local_row <= height:
            var row = row_begin + local_row
            var left_index = first_from + rows - row if reversed_order else first_from + row - 1
            var symbol = Int(sequences[unsafe_offset=left_index])
            var above_left = above_left_score
            var running_score = left_score
            var running_insertion = left_insertion
            comptime for column_slot in range(STRIP_COLUMNS):
                var substitution = Int32(table[unsafe_offset=symbol * width + Int(symbols[column_slot])])
                var computed = gotoh_cell[mode](
                    above_left,
                    scores[column_slot],
                    deletions[column_slot],
                    running_score,
                    running_insertion,
                    substitution,
                    scoring,
                )
                var cell = computed.score
                if local and column_slot < owned and cell > best_cell:
                    best_cell = cell
                    best_at = Int64(row) * Int64(columns + 1) + Int64(column_begin + first_column + column_slot + 1)
                above_left = scores[column_slot]
                scores[column_slot] = cell
                deletions[column_slot] = computed.deletion
                running_score = cell
                running_insertion = computed.insertion
            edge_score = running_score
            edge_insertion = running_insertion

            # A tile narrower than the full side is against the matrix's right edge, and nobody
            # reads its right column, so the gate is exact rather than conservative.
            if lane == STRIP_LANES - 1 and width_span == TILE_SIDE:
                left_scores[unsafe_offset=row_base + row] = running_score
                left_inserts[unsafe_offset=row_base + row] = running_insertion
            if local_row == height:
                # Column zero is the matrix border, so no lane computes it, but the crossing
                # reduction reads it as a candidate cut. Only the leftmost tile may publish it;
                # elsewhere that slot belongs to the tile on the left.
                if lane == 0 and column_begin == 0:
                    var border = Int32(0)
                    var border_delete = open + extend
                    if not local:
                        border = Int32(row) * extend if extends else open + Int32(row - 1) * extend
                        border_delete = border
                    top_scores[unsafe_offset=column_base] = border
                    top_deletes[unsafe_offset=column_base] = border_delete
                comptime for column_slot in range(STRIP_COLUMNS):
                    if column_slot < owned:
                        var column = column_begin + first_column + column_slot + 1
                        top_scores[unsafe_offset=column_base + column] = scores[column_slot]
                        top_deletes[unsafe_offset=column_base + column] = deletions[column_slot]
                if lane == STRIP_LANES - 1 and width_span == TILE_SIDE:
                    corner_scores[
                        unsafe_offset=corner_base + ((Int(tile_anti_diagonal) + 2) % 3) * corner_rows + tile_row + 1
                    ] = scores[STRIP_COLUMNS - 1]

    # A local sweep reports the best cell this warp saw, for the host to reduce across blocks.
    if local:
        var span = UInt32(1)
        while span < UInt32(STRIP_LANES):
            var theirs = shuffle_xor(best_cell, span)
            var their_place = shuffle_xor(best_at, span)
            if beats(theirs, their_place, best_cell, best_at):
                best_cell = theirs
                best_at = their_place
            span *= 2
        if thread_idx.x == 0:
            # One launch per tile-anti-diagonal, so this slot accumulates rather than replaces.
            var held = block_best[unsafe_offset=slot]
            if beats(best_cell, best_at, held, block_place[unsafe_offset=slot]):
                block_best[unsafe_offset=slot] = best_cell
                block_place[unsafe_offset=slot] = best_at


@fieldwise_init
struct SweepBuffers(Movable):
    """Device scratch reused across every sweep of one alignment."""

    var sequences: DeviceBuffer[SymbolDType]
    """Both sequences, concatenated."""
    var substitutions: DeviceBuffer[SubstitutionDType]
    """The substitution table, staged once per alignment."""
    var top_scores: DeviceBuffer[ScoreDType]
    """Score frontier along the top edge of each frame."""
    var top_deletes: DeviceBuffer[ScoreDType]
    """Deletion frontier along the same edge."""
    var reverse_scores: DeviceBuffer[ScoreDType]
    """Score frontier of the reverse half."""
    var reverse_deletes: DeviceBuffer[ScoreDType]
    """Deletion frontier of the reverse half."""
    var crossing: DeviceBuffer[ScoreDType]
    """Where the two halves meet, which is what the join reads."""
    var left_scores: DeviceBuffer[ScoreDType]
    """Score carry down the left edge of a strip."""
    var left_inserts: DeviceBuffer[ScoreDType]
    """Insertion carry down the same edge."""
    var corner_scores: DeviceBuffer[ScoreDType]
    """Tile corners, which is how one tile hands off to the next."""
    var corner_span: Int
    """How much of each shared array a whole recursion level may claim."""
    var left_span: Int
    """How many entries the left carry holds."""
    var frontier_span: Int
    """How many entries a frontier holds."""
    var block_best: DeviceBuffer[ScoreDType]
    """Best score each block found, for the local-alignment scan."""
    var block_place: DeviceBuffer[DType.int64]
    """Where each block found it, so ties break on position."""
    var block_slots: Int
    """How many blocks the scan reduces over."""


def device_sweep_buffers(
    scope: DeviceScope,
    rows: Int,
    columns: Int,
    sequences: ImmSpan[Scalar[SymbolDType], _],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
) raises -> SweepBuffers:
    """Device scratch for one alignment, sized so a whole recursion level fits side by side.

    Every sweep of a level claims its own slice, and a deep level is many tiny sweeps, so each
    array carries the real extent plus a constant per sweep. A level never holds more than two
    sweeps per row.

    The column extent is doubled because a frame's forward and reverse halves split its rows but
    both span all of its columns, so one level's sweeps cover the column axis twice.
    """
    # The sweep stages the whole table into shared memory unchecked, so the bound belongs here,
    # where every tiled launch passes, rather than in the kernel that cannot raise.
    if len(substitutions) > MAX_ALPHABET_SIZE * MAX_ALPHABET_SIZE:
        raise AlignmentError(ErrorKind.ALPHABET_TOO_LARGE, "staged substitution table")

    var frontier_span = 2 * columns + 4 * rows + 32
    var left_span = 5 * rows + 32
    var tile_columns_count = ceildiv(columns, TILE_SIDE) + 2
    # Square-in-count tiling makes a sweep's tile grid as wide as it is tall, and a level's frames partition the
    # columns, so the corner arrays of one level are bounded by twice the square of the full column count. Three
    # rotating buffers per sweep, each one entry per tile row. Linear in sequence length, where the full tile grid was
    # quadratic and overflowed its own int32 plan field.
    var corner_span = 3 * ceildiv(rows, MIN_TILE_HEIGHT) + 6 * rows + 32
    # Only a local scan writes here, and that is one sweep, so the grid is one tile-diagonal wide.
    var block_slots = min(ceildiv(rows, MIN_TILE_HEIGHT), tile_columns_count) + 4
    var buffers = SweepBuffers(
        allocate[SymbolDType](scope, max(rows + columns, 1)),
        allocate[SubstitutionDType](scope, len(substitutions)),
        allocate[ScoreDType](scope, frontier_span),
        allocate[ScoreDType](scope, frontier_span),
        allocate[ScoreDType](scope, frontier_span),
        allocate[ScoreDType](scope, frontier_span),
        allocate[ScoreDType](scope, 4 * max(rows, 1)),
        allocate[ScoreDType](scope, left_span),
        allocate[ScoreDType](scope, left_span),
        allocate[ScoreDType](scope, corner_span),
        corner_span,
        left_span,
        frontier_span,
        zeroed[ScoreDType](scope, block_slots),
        zeroed[DType.int64](scope, block_slots),
        block_slots,
    )
    scope.context.enqueue_copy(buffers.sequences, sequences)
    scope.context.enqueue_copy(buffers.substitutions, substitutions)
    return buffers^


def crossing_kernel(
    forward_scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    forward_deletes: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    reverse_scores: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    reverse_deletes: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    joins: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    crossing: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    refund: Int32,
):
    """Picks where each split of a level crosses its cut, one block per split.

    Reading the frontiers back per split costs more than the sweeps once the recursion is
    thousands of nodes deep, so the reduction happens here and four numbers per split travel.
    """
    var split = Int(block_idx.x)
    var forward_base = Int(joins[unsafe_offset=split * 3])
    var reverse_base = Int(joins[unsafe_offset=split * 3 + 1])
    var span = Int(joins[unsafe_offset=split * 3 + 2])

    var plain_scores = stack_allocation[THREADS_PER_BLOCK, Scalar[ScoreDType], address_space=AddressSpace.SHARED]()
    var plain_places = stack_allocation[THREADS_PER_BLOCK, Scalar[DType.int64], address_space=AddressSpace.SHARED]()
    var gapped_scores = stack_allocation[THREADS_PER_BLOCK, Scalar[ScoreDType], address_space=AddressSpace.SHARED]()
    var gapped_places = stack_allocation[THREADS_PER_BLOCK, Scalar[DType.int64], address_space=AddressSpace.SHARED]()

    var best_plain = NEGATIVE_INFINITY
    var best_plain_column = Int64(0)
    var best_gapped = NEGATIVE_INFINITY
    var best_gapped_column = Int64(0)
    for offset in range(Int(thread_idx.x), span + 1, Int(block_dim.x)):
        var near = forward_base + offset
        var far = reverse_base + span - offset
        var plain = forward_scores[unsafe_offset=near] + reverse_scores[unsafe_offset=far]
        if plain > best_plain:
            best_plain = plain
            best_plain_column = Int64(offset)
        var gapped = forward_deletes[unsafe_offset=near] + reverse_deletes[unsafe_offset=far] + refund
        if gapped > best_gapped:
            best_gapped = gapped
            best_gapped_column = Int64(offset)

    block_argmax(plain_scores, plain_places, best_plain, best_plain_column)
    block_argmax(gapped_scores, gapped_places, best_gapped, best_gapped_column)
    if thread_idx.x == 0:
        crossing[unsafe_offset=split * 4] = plain_scores[unsafe_offset=0]
        crossing[unsafe_offset=split * 4 + 1] = Scalar[ScoreDType](plain_places[unsafe_offset=0])
        crossing[unsafe_offset=split * 4 + 2] = gapped_scores[unsafe_offset=0]
        crossing[unsafe_offset=split * 4 + 3] = Scalar[ScoreDType](gapped_places[unsafe_offset=0])


@fieldwise_init
struct Sweep(ImplicitlyCopyable, TrivialRegisterPassable):
    """One independent sub-rectangle sweep, as the plan-driven kernel needs to see it."""

    var rows: Int
    """How many rows this sub-rectangle covers."""
    var columns: Int
    """How many columns it covers."""
    var first_from: Int
    """First row of the sub-rectangle."""
    var second_from: Int
    """First column of the sub-rectangle."""
    var row_base: Int
    """Where its row-indexed scratch begins."""
    var column_base: Int
    """Where its column-indexed scratch begins."""
    var entering_run: GapRun
    """Whether a gap run is already open where the sweep starts."""
    var half: SweepHalf
    """Which direction the sweep runs."""

    def tile_height(self, sweeps_in_level: Int, target_tiles: Int) -> Int:
        """The tallest tile that still fills the machine, given how many sweeps share the level.

        Concurrency is `sweeps_in_level * min(tile_rows, tile_columns)`, so shorter tiles buy
        parallelism only until the machine is full; past that they cost skew, since a tile `h`
        rows tall runs `h + 31` steps. Once `tile_columns` alone caps the anti-diagonal, the best
        a height can do is make the grid square in count.
        """
        var wide = self.tile_columns()
        var want = max(target_tiles // max(sweeps_in_level, 1), 1)
        var enough = min(want, wide)
        var fair = ceildiv(self.rows, enough)
        return clamp(fair, MIN_TILE_HEIGHT, TILE_SIDE)

    def tile_rows(self, sweeps_in_level: Int, target_tiles: Int) -> Int:
        """How many tiles of `tile_height` rows cover the sweep, at least one."""
        var height = self.tile_height(sweeps_in_level, target_tiles)
        return max(ceildiv(self.rows, height), 1)

    def tile_columns(self) -> Int:
        """How many `TILE_SIDE`-wide tiles cover the sweep, at least one."""
        return max(ceildiv(self.columns, TILE_SIDE), 1)

    def tile_anti_diagonals(self, sweeps_in_level: Int, target_tiles: Int) -> Int:
        """How many anti-diagonals of tiles the sweep takes, one launch each."""
        return self.tile_rows(sweeps_in_level, target_tiles) + self.tile_columns() - 1


def device_local_extremum[
    half: SweepHalf
](
    scope: DeviceScope,
    mut buffers: SweepBuffers,
    rows: Int,
    first_to: Int,
    second_to: Int,
    alphabet_size: Int,
    scoring: AffineGapCosts,
) raises -> Tuple[Int, Int, Int32]:
    """Finds where the best local alignment ends, tiled across the device rather than one block.

    Local alignment has unknown endpoints, so this runs twice: forwards it says where the best
    alignment ends, and backwards over those prefixes how far the same alignment reaches back.
    Ties resolve to the earliest cell in row-major order, matching the reference scan.
    """
    var sweeps = List[Sweep]()
    sweeps.append(Sweep(first_to, second_to, 0, rows, 0, 0, GapRun.OPENS, half))
    # Blocks that fall outside their sweep return without writing, and the buffers outlive the
    # scan, so anything left from an earlier one would be read as a candidate.
    scope.context.enqueue_memset(buffers.block_best, Scalar[ScoreDType](0))
    scope.context.enqueue_memset(buffers.block_place, Scalar[DType.int64](0))
    device_sweep_level[AlignmentMode.LOCAL](scope, buffers, sweeps, alphabet_size, scoring)

    var best = Int32(0)
    var best_place = Int64(0)
    var slots = min(sweeps[0].tile_rows(len(sweeps), target_tiles_for(scope.specs)), sweeps[0].tile_columns())
    with buffers.block_best.map_to_host() as scores, buffers.block_place.map_to_host() as places:
        for slot in range(min(slots, buffers.block_slots)):
            var candidate = scores[slot]
            if beats(candidate, places[slot], best, best_place):
                best = candidate
                best_place = places[slot]
    var stride = Int64(second_to + 1)
    return (Int(best_place // stride), Int(best_place % stride), best)


def device_sweep_level[
    mode: AlignmentMode
](
    scope: DeviceScope, mut buffers: SweepBuffers, mut sweeps: List[Sweep], alphabet_size: Int, scoring: AffineGapCosts
) raises:
    """Sweeps every independent sub-rectangle of one recursion level together.

    Depth-first recursion offers the parallelism of a single frame, which halves as the recursion
    deepens while its work halves too, so utilization falls as fast as the work does. Frames at a
    level partition both axes, so they share the frontier arrays without touching, and one flat
    block index carries both the sweep and its tile.
    """
    if len(sweeps) == 0:
        return

    var target_tiles = target_tiles_for(scope.specs)
    var plan = List[SweepPlan]()
    var widest_tiles = 1
    var deepest = 0
    var corner_base = 0
    var top_base = 0
    var left_base = 0
    for index in range(len(sweeps)):
        var sweep = sweeps[index]
        var tile_rows = sweep.tile_rows(len(sweeps), target_tiles)
        var tile_columns = sweep.tile_columns()
        plan.append(
            SweepPlan(
                Int64(sweep.rows),
                Int64(sweep.columns),
                Int64(sweep.first_from),
                Int64(sweep.second_from),
                Int64(left_base),
                Int64(top_base),
                Int64(corner_base),
                Int64(tile_rows),
                Int64(tile_columns),
                Int64(sweep.tile_height(len(sweeps), target_tiles)),
                sweep.entering_run,
                sweep.half,
            )
        )
        corner_base += 3 * (tile_rows + 1)
        sweeps[index].row_base = left_base
        sweeps[index].column_base = top_base
        left_base += sweep.rows + 2
        top_base += sweep.columns + 2
        widest_tiles = max(widest_tiles, min(tile_rows, tile_columns))
        deepest = max(deepest, sweep.tile_anti_diagonals(len(sweeps), target_tiles))

    if corner_base > buffers.corner_span or left_base > buffers.left_span or top_base > buffers.frontier_span:
        raise AlignmentError(
            ErrorKind.SCRATCH_TOO_SMALL,
            String("corner {}/{} left {}/{} top {}/{}").format(
                corner_base, buffers.corner_span, left_base, buffers.left_span, top_base, buffers.frontier_span
            ),
        )

    # The plan crosses to the device as raw words, because a kernel takes a pointer, not a `List`.
    #
    # Every word below is written, so filling first would be a memset immediately overwritten.
    var words = List[Scalar[PlanDType]](unsafe_uninit_length=len(plan) * PLAN_WORDS)
    var source = plan.unsafe_ptr().unsafe_bitcast[Scalar[PlanDType]]()
    for index in range(len(words)):
        words[index] = source[unsafe_offset=index]
    var plan_buffer = upload(scope, Span(words))
    for tile_anti_diagonal in range(deepest):
        scope.context.enqueue_function[tiled_sweep_kernel[mode]](
            buffers.sequences.unsafe_ptr(),
            buffers.substitutions.unsafe_ptr(),
            plan_buffer.unsafe_ptr(),
            buffers.block_best.unsafe_ptr(),
            buffers.block_place.unsafe_ptr(),
            buffers.top_scores.unsafe_ptr(),
            buffers.top_deletes.unsafe_ptr(),
            buffers.reverse_scores.unsafe_ptr(),
            buffers.reverse_deletes.unsafe_ptr(),
            buffers.left_scores.unsafe_ptr(),
            buffers.left_inserts.unsafe_ptr(),
            buffers.corner_scores.unsafe_ptr(),
            Int32(tile_anti_diagonal),
            Int32(widest_tiles),
            Int32(alphabet_size),
            scoring.open,
            scoring.extend,
            grid_dim=widest_tiles * len(sweeps),
            block_dim=STRIP_LANES,
        )
    scope.context.synchronize()


def border_score[mode: AlignmentMode](rows: Int, columns: Int, scoring: AffineGapCosts) -> Int32:
    """The score of a rectangle with no interior, which no sweep writes a frontier for."""
    comptime if mode == AlignmentMode.LOCAL:
        return Int32(0)
    return scoring.run(rows + columns)


def device_score[
    mode: AlignmentMode
](
    scope: DeviceScope,
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    scoring: AffineGapCosts,
) raises -> Int32:
    """One pair scored on the device in linear space, over global-memory bands and without a band.

    The counterpart of `device_align` for a caller that wants only the number, so a pair too tall
    for one block's carry is scored rather than refused.
    """
    var rows = len(first)
    var columns = len(second)
    # Every tile of an empty rectangle returns before writing, and the frontier is never zeroed,
    # so the borders are computed here rather than read back.
    if rows == 0 or columns == 0:
        return border_score[mode](rows, columns, scoring)

    var sequences = List[Scalar[SymbolDType]](capacity=rows + columns)
    sequences.extend(first)
    sequences.extend(second)
    var buffers = device_sweep_buffers(scope, rows, columns, Span(sequences), substitutions)

    comptime if mode == AlignmentMode.LOCAL:
        var last_row, last_column, best = device_local_extremum[SweepHalf.FORWARD](
            scope, buffers, rows, rows, columns, alphabet_size, scoring
        )
        _ = last_row
        _ = last_column
        return best

    var sweeps = List[Sweep]()
    sweeps.append(Sweep(rows, columns, 0, rows, 0, 0, GapRun.OPENS, SweepHalf.FORWARD))
    device_sweep_level[AlignmentMode.GLOBAL](scope, buffers, sweeps, alphabet_size, scoring)
    with buffers.top_scores.map_to_host() as frontier:
        return frontier[sweeps[0].column_base + columns]


def device_align[
    mode: AlignmentMode
](
    scope: DeviceScope,
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    scoring: AffineGapCosts,
    alphabet: String,
    leaf_cells: Int,
    placement: Placement,
) raises -> GappedAlignment:
    """One pair aligned with every sweep on the device, in linear space.

    The device counterpart of `serial_align`, so a caller chooses where the work runs without
    also choosing how the traceback stores its state.
    """
    var rows = len(first)
    var columns = len(second)
    var path = RowPath(rows)

    var sequences = List[Scalar[SymbolDType]](capacity=rows + columns)
    sequences.extend(first)
    sequences.extend(second)
    # One scratch set for the whole alignment: the local scans and the recursion take it in turn.
    var buffers = device_sweep_buffers(scope, rows, columns, Span(sequences), substitutions)

    comptime if mode == AlignmentMode.GLOBAL:
        device_hirschberg(
            scope,
            buffers,
            first,
            second,
            Rectangle(0, rows, 0, columns),
            substitutions,
            alphabet_size,
            scoring,
            leaf_cells,
            placement,
            path.view(),
        )
        var whole = path.gapped(first, second, alphabet, mode, 0, rows)
        return GappedAlignment(path.score(first, second, substitutions, alphabet_size, scoring), whole[0], whole[1])

    var last_row, last_column, score = device_local_extremum[SweepHalf.FORWARD](
        scope, buffers, rows, rows, columns, alphabet_size, scoring
    )

    # A non-positive best means the local alignment is empty, so the walk never runs and the two
    # ends stay where the forward scan left them.
    var first_row = last_row
    if score > 0:
        var back_rows, back_columns, _ = device_local_extremum[SweepHalf.REVERSE](
            scope, buffers, rows, last_row, last_column, alphabet_size, scoring
        )
        first_row = last_row - back_rows
        var first_column = last_column - back_columns
        device_hirschberg(
            scope,
            buffers,
            first,
            second,
            Rectangle(first_row, last_row, first_column, last_column),
            substitutions,
            alphabet_size,
            scoring,
            leaf_cells,
            placement,
            path.view(),
        )

    var window = path.gapped(first, second, alphabet, mode, first_row, last_row)
    return GappedAlignment(score, window[0], window[1])


# endregion GPU Wavefront
