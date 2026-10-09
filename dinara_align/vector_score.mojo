# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The host's affine-gap sweeps under any substitution table that find more than a score.

The same Gotoh recurrence and borders as `alignment.serial_align`, so the same scores, swept by anti-diagonal
a vector of cells a step (see `anti_diagonals`): the cell, swept back from a local alignment's end, where
it earns its score, which is where it starts (`reach_back`); and a global alignment's three layers, stored
in the band of diagonals its score bounds, to be traced (`vector_align`). A score alone is the tabled
sweep's (see `scoring.swept_score`).
"""

from .common import UNREACHED
from .anti_diagonals import AntiDiagonals, GapLanes, column_letters, gotoh_lanes, row_letters
from .gotoh import (
    AffineGapCosts,
    AlignmentMode,
    AntiDiagonalMajor,
    BAND_PADDING,
    GappedAlignment,
    reconstruct,
    serial_align,
)
from .substitutions import SubstitutionLookup

comptime WIDTH = 16
"""Cells a step computes at once."""

comptime Lanes = SIMD[DType.int32, WIDTH]


def reach_back(
    first: List[UInt8],
    second: List[UInt8],
    end_row: Int,
    end_column: Int,
    lookup: SubstitutionLookup,
    gaps: AffineGapCosts,
    target: Int32,
    narrow: Bool,
) -> Tuple[Int, Int]:
    """Where a local alignment ending at `(end_row, end_column)` and scoring `target`, the best there is,
    starts: the first cell, by anti-diagonal and then row, that a global alignment from it to the end
    scores `target` from. In lanes of 16 bits with `narrow`, which the caller says the prefixes' scores
    fit (see `anti_diagonals.fits_16_bits`).

    The sweep runs back from the end over both prefixes, as a global one anchored there: a cell's score
    is the best global alignment of the letters between it and the end. None can pass `target`, which no
    local alignment passes, and the start of one that earns it reaches it, so the first cell that does
    starts an optimal alignment, and the global alignment between it and the end is one. The sweep stops
    on the diagonal it finds it on, so it covers about the alignment's own cells, not the prefixes'.
    """
    if narrow:
        return reached_from[DType.int16, 2 * WIDTH](first, second, end_row, end_column, lookup, gaps, target)
    return reached_from[DType.int32, WIDTH](first, second, end_row, end_column, lookup, gaps, target)


def reached_from[
    dtype: DType, width: Int
](
    first: List[UInt8],
    second: List[UInt8],
    end_row: Int,
    end_column: Int,
    lookup: SubstitutionLookup,
    gaps: AffineGapCosts,
    target: Int32,
) -> Tuple[Int, Int]:
    """`reach_back` in lanes of `dtype`, `width` of them."""
    comptime Value = Scalar[dtype]
    comptime Lanes = SIMD[dtype, width]
    var rows = end_row
    var columns = end_column
    var opened = Int(gaps.open + gaps.extend)

    @inline(.always)
    def border(length: Int) {imm gaps} -> Int:
        """The score of a gap run of `length` along the anchored border."""
        return Int(gaps.run(length))

    # Row `i` reads the `i`-th letter back from the end, and so does column `j`.
    var sweep = AntiDiagonals[dtype, width](Span(first)[:end_row], Span(second)[:end_column], 0, True, True)
    var cells = sweep.cells()
    cells.one_back[unsafe_offset=0] = Value(border(1))
    cells.deletes_back[unsafe_offset=0] = Value(border(1) + opened)
    cells.one_back[unsafe_offset=1] = Value(border(1))
    cells.inserts_back[unsafe_offset=1] = Value(border(1) + opened)

    var substitute = lookup.lanes[width, dtype]()
    var steps = GapLanes[dtype, width].symmetric(Int(gaps.open), Int(gaps.extend))
    var wanted = Lanes(Value(target))
    var lane_rows = Lanes()
    comptime for lane in range(width):
        lane_rows[lane] = Value(lane)
    for diagonal in range(2, rows + columns + 1):
        var high = min(rows, diagonal - 1)
        var lag = columns - diagonal
        var row = max(1, diagonal - columns)
        while row <= high:
            var score = cells.step(row, lag, substitute, steps)
            var reached = score.ge(wanted) & lane_rows.lt(Lanes(Value(high - row + 1)))
            if reached.reduce_or():
                # The lanes run by row, so the first that reached it is the diagonal's earliest.
                for lane in range(width):
                    if reached[lane]:
                        var back_rows = row + lane
                        return (end_row - back_rows, end_column - (diagonal - back_rows))
            row += width
        if diagonal <= columns:
            cells.current[unsafe_offset=0] = Value(border(diagonal))
            cells.deletes[unsafe_offset=0] = Value(border(diagonal) + opened)
        if diagonal <= rows:
            cells.current[unsafe_offset=diagonal] = Value(border(diagonal))
            cells.inserts[unsafe_offset=diagonal] = Value(border(diagonal) + opened)
        cells.advance()
    # The whole prefixes, which a caller with a score above zero never needs.
    return (0, 0)


def vector_cells(rows: Int, columns: Int, band_cells: Int) -> Int:
    """What `vector_align` holds for a `rows` by `columns` pair whose band has `band_cells` cells, in cells of
    its three 32-bit layers: every anti-diagonal's cells and their padding, a step's spare lanes, and the
    three indexes of the diagonals, two cells' bytes a diagonal."""
    var diagonals = rows + columns + 1
    return band_cells + diagonals * (2 * BAND_PADDING + 2) + WIDTH


def vector_align(
    first: List[UInt8],
    second: List[UInt8],
    lookup: SubstitutionLookup,
    gaps: AffineGapCosts,
    substitutions: List[Scalar[DType.int8]],
    alphabet_size: Int,
    alphabet: String,
    low_diagonal: Int = Int.MIN,
    high_diagonal: Int = Int.MAX,
) -> GappedAlignment:
    """`serial_align`'s global alignment, its three layers swept sixteen cells at a time.

    Only the cells on diagonals `low_diagonal ..= high_diagonal`, column minus row, are computed and
    stored, the rest reading as unreachable; given a band holding every optimal path (see
    `optimal_band`), those cells hold the same values as `serial_align`'s, so `reconstruct` walks the
    same path, in memory that grows with the band rather than the matrix. They are stored
    anti-diagonal by anti-diagonal, each diagonal's band contiguous by row, so a step loads its
    neighbours from the two diagonals before as vectors (see `AntiDiagonalMajor`).
    """
    var rows = len(first)
    var columns = len(second)
    if rows == 0 or columns == 0:
        return serial_align[AlignmentMode.GLOBAL](first, second, substitutions, alphabet_size, gaps, alphabet)
    var open = gaps.open
    var extend = gaps.extend
    var low_band = max(low_diagonal, -rows)
    var high_band = min(high_diagonal, columns)

    @inline(.always)
    def border(length: Int) {imm gaps} -> Int32:
        """The score of a gap run of `length` along the global border."""
        return gaps.run(length)

    # Diagonal `d`'s cells lie on rows `lows[d] ..= highs[d]`: inside the matrix, and with
    # `column - row` inside the band, column being `d - row`.
    var diagonals = rows + columns + 1
    var lows = List[Int](length=diagonals, fill=0)
    var highs = List[Int](length=diagonals, fill=0)
    var starts = List[Int](length=diagonals + 1, fill=0)
    for diagonal in range(diagonals):
        var low = max(0, diagonal - columns, -((high_band - diagonal) // 2))
        var high = min(rows, diagonal, (diagonal - low_band) // 2)
        lows[diagonal] = low
        highs[diagonal] = high
        starts[diagonal + 1] = starts[diagonal] + max(high - low + 1, 0) + 2 * BAND_PADDING
    var cells = starts[diagonals]
    var scores = List[Int32](unsafe_uninit_length=cells + WIDTH)
    var deletes = List[Int32](unsafe_uninit_length=cells + WIDTH)
    var inserts = List[Int32](unsafe_uninit_length=cells + WIDTH)
    var letters = row_letters[WIDTH](Span(first))
    var reversed = column_letters[WIDTH](Span(second))

    var score_cells = scores.unsafe_ptr()
    var delete_cells = deletes.unsafe_ptr()
    var insert_cells = inserts.unsafe_ptr()

    @inline(.always)
    def unreachable(index: Int) {imm score_cells, imm delete_cells, imm insert_cells}:
        """Marks stored cell `index` unreachable in all three layers, as band padding reads."""
        score_cells[unsafe_offset=index] = UNREACHED
        delete_cells[unsafe_offset=index] = UNREACHED
        insert_cells[unsafe_offset=index] = UNREACHED

    # Diagonal zero: the origin, which every band holds.
    for index in range(BAND_PADDING):
        unreachable(index)
        unreachable(BAND_PADDING + 1 + index)
    score_cells[unsafe_offset=BAND_PADDING] = 0
    delete_cells[unsafe_offset=BAND_PADDING] = 0
    insert_cells[unsafe_offset=BAND_PADDING] = 0
    var substitute = lookup.lanes[WIDTH]()
    var steps = GapLanes[DType.int32, WIDTH].symmetric(Int(open), Int(extend))
    for diagonal in range(1, diagonals):
        var first_row = lows[diagonal]
        var last_row = highs[diagonal]
        # Diagonal `d`'s cell on row `i` sits at `here + i`; its neighbours on the two before likewise.
        var here = starts[diagonal] + BAND_PADDING - first_row
        var one_back = starts[diagonal - 1] + BAND_PADDING - lows[diagonal - 1]
        var two_back = starts[diagonal - 2] + BAND_PADDING - lows[diagonal - 2] if diagonal >= 2 else 0
        for index in range(BAND_PADDING):
            unreachable(starts[diagonal] + index)
        var low = max(1, first_row)
        var high = min(last_row, diagonal - 1)
        var lag = columns - diagonal
        var row = low
        while row <= high:
            var above = score_cells.unsafe_offset(one_back + row - 1).unsafe_load[width=WIDTH]()
            var above_delete = delete_cells.unsafe_offset(one_back + row - 1).unsafe_load[width=WIDTH]()
            var left = score_cells.unsafe_offset(one_back + row).unsafe_load[width=WIDTH]()
            var left_insert = insert_cells.unsafe_offset(one_back + row).unsafe_load[width=WIDTH]()
            var above_left = score_cells.unsafe_offset(two_back + row - 1).unsafe_load[width=WIDTH]()
            var mine = letters.unsafe_ptr().unsafe_offset(row).unsafe_load[width=WIDTH]()
            var theirs = reversed.unsafe_ptr().unsafe_offset(lag + row).unsafe_load[width=WIDTH]()
            var cell = gotoh_lanes[DType.int32, WIDTH](
                above_left,
                above,
                above_delete,
                left,
                left_insert,
                substitute(mine, theirs),
                steps.down_first,
                steps.down_further,
                steps.across_first,
                steps.across_further,
            )
            score_cells.unsafe_offset(here + row).unsafe_store(cell[0])
            delete_cells.unsafe_offset(here + row).unsafe_store(cell[1])
            insert_cells.unsafe_offset(here + row).unsafe_store(cell[2])
            row += WIDTH
        # The border cells inside the band, then the padding after it, which the lanes may have run over.
        if first_row == 0:
            score_cells[unsafe_offset=here] = border(diagonal)
            delete_cells[unsafe_offset=here] = score_cells[unsafe_offset=here] + open + extend
            insert_cells[unsafe_offset=here] = 0
        if last_row == diagonal:
            score_cells[unsafe_offset=here + diagonal] = border(diagonal)
            delete_cells[unsafe_offset=here + diagonal] = 0
            insert_cells[unsafe_offset=here + diagonal] = score_cells[unsafe_offset=here + diagonal] + open + extend
        for index in range(BAND_PADDING):
            unreachable(here + last_row + 1 + index)

    var layout = AntiDiagonalMajor(
        starts.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        lows.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
    )
    var reconstruction = reconstruct[AlignmentMode.GLOBAL](
        scores,
        deletes,
        inserts,
        layout,
        first,
        second,
        substitutions,
        alphabet_size,
        rows,
        columns,
        alphabet,
        gaps,
    )
    var final_score = scores[layout.index(rows, columns)]
    # The layout reads `starts` and `lows` through pointers, so both must outlive every use of it.
    _ = len(starts)
    _ = len(lows)
    _ = len(lookup.cells)
    return GappedAlignment(final_score, reconstruction[0], reconstruction[1])


def optimal_band(rows: Int, columns: Int, reward: Int, gaps: AffineGapCosts, score: Int) -> Tuple[Int, Int]:
    """The diagonals, column minus row, every optimal global alignment of score `score` stays on, under a
    table whose largest score is `reward`.

    No substitution earns more than `reward`, so an alignment of `n` and `m` letters with `g` of them in
    gaps scores at most `reward (n + m - g) / 2 + e g`, `e` the extension's score: each gapped letter
    gives up at least `2 e + reward` of the most `reward (n + m)` it could score twice over. A path
    reaching diagonal `t` past the start's 0 and the end's `m - n` gaps at least `|m - n|` letters plus
    twice the detour, so an optimal one never strays further than what it gave up allows. With nothing
    to spend on gaps the band is the two diagonals and the run between them.
    """
    var difference = columns - rows
    var per_gap = 2 * -Int(gaps.extend) + reward
    var cost = reward * (rows + columns) - 2 * score
    var reach = max(0, (cost - per_gap * abs(difference)) // (2 * per_gap)) if per_gap > 0 else rows + columns
    return (min(0, difference) - reach, max(0, difference) + reach)
