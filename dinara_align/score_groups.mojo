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
"""

from std.math import ceildiv, clamp
from std.memory import stack_allocation
from std.memory.pointer import AddressSpace
from max.gpu import WARP_SIZE, barrier, block_idx, thread_idx
from max.gpu.primitives.warp import shuffle_up, shuffle_xor

from .alignment import AffineGapCosts, AlignmentMode, gotoh_cell
from .common import (
    DeviceScope,
    MAX_ALPHABET_SIZE,
    OffsetDType,
    ScoreDType,
    SubstitutionDType,
    SymbolDType,
    THREADS_PER_BLOCK,
    upload,
    zeroed,
)

comptime SHAPE_LANES: Array[Int, 12] = [4, 4, 4, 8, 8, 8, 16, 16, 16, 32, 32, 32]
"""Each compiled shape's lanes a pair."""
comptime SHAPE_COLUMNS: Array[Int, 12] = [8, 12, 16, 8, 12, 16, 8, 12, 16, 8, 12, 16]
"""Each compiled shape's columns a lane: a pair spans 32 to 512 columns in steps of about a third."""
comptime SHAPES = 12


def group_score_kernel[
    mode: AlignmentMode, lanes: Int, columns_per_lane: Int
](
    sequences: Pointer[Scalar[SymbolDType], MutAnyOrigin],
    offsets: Pointer[Scalar[OffsetDType], MutAnyOrigin],
    substitutions: Pointer[Scalar[SubstitutionDType], MutAnyOrigin],
    results: Pointer[Scalar[ScoreDType], MutAnyOrigin],
    pairs: Int32,
    alphabet_size: Int32,
    open: Int32,
    extend: Int32,
):
    """Scores `THREADS_PER_BLOCK // lanes` pairs a block, `lanes` lanes a pair (see the module's
    notes). Pair `p`'s first sequence runs from `offsets[2 p]`, its second from `offsets[2 p + 1]` to
    `offsets[2 p + 2]`; the first's letters are the rows."""
    comptime pairs_per_block = THREADS_PER_BLOCK // lanes
    var width = Int(alphabet_size)
    var table = stack_allocation[
        MAX_ALPHABET_SIZE * MAX_ALPHABET_SIZE, Scalar[SubstitutionDType], address_space=AddressSpace.SHARED
    ]()
    for index in range(Int(thread_idx.x), width * width, THREADS_PER_BLOCK):
        table[unsafe_offset=index] = substitutions[unsafe_offset=index]
    barrier()

    var scoring = AffineGapCosts(open, extend)
    var thread = Int(thread_idx.x)
    var member = thread % lanes
    var pair = Int(block_idx.x) * pairs_per_block + thread // lanes
    var live = pair < Int(pairs)
    var first_start = 0
    var second_start = 0
    var rows = 0
    var columns = 0
    if live:
        first_start = Int(offsets[unsafe_offset=2 * pair])
        second_start = Int(offsets[unsafe_offset=2 * pair + 1])
        rows = second_start - first_start
        columns = Int(offsets[unsafe_offset=2 * pair + 2]) - second_start

    # The warp steps as far as its longest pair needs, every lane in every step.
    var steps = Int32(rows + lanes)
    var span = UInt32(1)
    while span < UInt32(WARP_SIZE):
        steps = max(steps, shuffle_xor(steps, span))
        span *= 2

    var first_column = member * columns_per_lane
    var owned = clamp(columns - first_column, 0, columns_per_lane)
    var symbols = Array[Int32, columns_per_lane](fill=0)
    comptime for slot in range(columns_per_lane):
        var column = clamp(first_column + slot, 0, max(columns - 1, 0))
        symbols[slot] = Int32(sequences[unsafe_offset=second_start + column])

    # Row zero: a global alignment pays a gap for every column it starts past; a local one nothing.
    var scores = Array[Int32, columns_per_lane](fill=0)
    var deletions = Array[Int32, columns_per_lane](fill=0)
    comptime for slot in range(columns_per_lane):
        comptime if mode == AlignmentMode.GLOBAL:
            scores[slot] = scoring.open + Int32(first_column + slot) * scoring.extend
        deletions[slot] = scores[slot] + scoring.open + scoring.extend
    var edge_score = Int32(0)
    comptime if mode == AlignmentMode.GLOBAL:
        edge_score = scoring.open + Int32(first_column + columns_per_lane - 1) * scoring.extend
    var edge_insertion = edge_score + scoring.open + scoring.extend
    var above_left_carry = Int32(0)
    comptime if mode == AlignmentMode.GLOBAL:
        above_left_carry = 0 if first_column == 0 else scoring.open + Int32(first_column - 1) * scoring.extend

    var reported = Int32(0)
    comptime if mode == AlignmentMode.GLOBAL:
        var letters = rows + columns
        reported = 0 if letters == 0 else scoring.open + Int32(letters - 1) * scoring.extend
    var best = Int32(0)

    for step in range(1, Int(steps) + 1):
        var row = step - member
        # The lane before hands over its last column of this row, which it computed a step ago; the
        # group's first lane reads the matrix's left edge instead.
        var left_score = shuffle_up(edge_score, 1)
        var left_insertion = shuffle_up(edge_insertion, 1)
        if member == 0:
            var border = Int32(0)
            comptime if mode == AlignmentMode.GLOBAL:
                var here = clamp(row, 0, rows)
                border = 0 if here == 0 else scoring.open + Int32(here - 1) * scoring.extend
            left_score = border
            left_insertion = border + scoring.open + scoring.extend
        var above_left = above_left_carry
        above_left_carry = left_score
        if row >= 1 and row <= rows:
            var symbol = Int(sequences[unsafe_offset=first_start + row - 1])
            var running_score = left_score
            var running_insertion = left_insertion
            comptime for slot in range(columns_per_lane):
                var substitution = Int32(table[unsafe_offset=symbol * width + Int(symbols[slot])])
                var computed = gotoh_cell[mode](
                    above_left,
                    scores[slot],
                    deletions[slot],
                    running_score,
                    running_insertion,
                    substitution,
                    scoring,
                )
                above_left = scores[slot]
                scores[slot] = computed.score
                deletions[slot] = computed.deletion
                running_score = computed.score
                running_insertion = computed.insertion
                comptime if mode == AlignmentMode.LOCAL:
                    if slot < owned:
                        best = max(best, computed.score)
            edge_score = running_score
            edge_insertion = running_insertion
            comptime if mode == AlignmentMode.GLOBAL:
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
        # The lane holding the last column holds the corner; with no columns, the first lane.
        var corner = member == 0 if columns == 0 else (owned > 0 and first_column + owned == columns)
        if live and corner:
            results[unsafe_offset=pair] = reported


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


def grouped_scores[
    mode: AlignmentMode
](
    scope: DeviceScope,
    sequences: ImmSpan[Scalar[SymbolDType], _],
    offsets: List[Scalar[OffsetDType]],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    scoring: AffineGapCosts,
) raises -> Optional[List[Int32]]:
    """Every pair's score, several pairs a warp (see the module's notes), or None when the batch's
    widest second sequence is wider than any shape spans."""
    var pairs = (len(offsets) - 1) // 2
    var longest_rows = 0
    var longest_columns = 0
    for pair in range(pairs):
        longest_rows = max(longest_rows, Int(offsets[2 * pair + 1]) - Int(offsets[2 * pair]))
        longest_columns = max(longest_columns, Int(offsets[2 * pair + 2]) - Int(offsets[2 * pair + 1]))
    var chosen = shape_for(longest_rows, longest_columns)
    if chosen < 0 or pairs == 0:
        return None
    var sequences_buffer = upload(scope, sequences)
    var offsets_buffer = upload(scope, offsets)
    var substitutions_buffer = upload(scope, substitutions)
    var results_buffer = zeroed[ScoreDType](scope, pairs)
    comptime for shape in range(SHAPES):
        if shape == chosen:
            comptime lanes = SHAPE_LANES[shape]
            comptime kernel = group_score_kernel[mode, lanes, SHAPE_COLUMNS[shape]]
            scope.context.enqueue_function[kernel](
                sequences_buffer.unsafe_ptr(),
                offsets_buffer.unsafe_ptr(),
                substitutions_buffer.unsafe_ptr(),
                results_buffer.unsafe_ptr(),
                Int32(pairs),
                Int32(alphabet_size),
                scoring.open,
                scoring.extend,
                grid_dim=ceildiv(pairs, THREADS_PER_BLOCK // lanes),
                block_dim=THREADS_PER_BLOCK,
            )
    scope.context.synchronize()
    var results = List[Int32](capacity=pairs)
    with results_buffer.map_to_host() as host:
        for index in range(pairs):
            results.append(host[index])
    return results^
