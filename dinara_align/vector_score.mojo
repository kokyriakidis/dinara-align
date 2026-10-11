# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The host's affine-gap sweeps under any substitution table that find more than a score.

The same Gotoh recurrence and borders as `gotoh.serial_align`, so the same scores, swept by anti-diagonal
a vector of cells a step (see `anti_diagonals`): the cell, swept back from a local alignment's end, where
it earns its score, which is where it starts (`reach_back`); and a global alignment's decisions, stored in
the band of diagonals its score bounds, to be traced (`vector_align`). A score alone is the tabled
sweep's (see `scoring.swept_score`).
"""

from .common import GAP_BYTE
from .anti_diagonals import (
    AntiDiagonals,
    GapLanes,
    best_move,
    column_letters,
    fits_16_bits,
    gap_layer,
    row_letters,
    straight_deficit,
)
from .gotoh import AffineGapCosts, AlignmentMode, Decision, GapRun, GappedAlignment, Layer, advance, serial_align
from .substitutions import SubstitutionLookup, table_extremes

comptime WIDTH = 16
"""Cells a step computes at once."""


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


def vector_bytes(rows: Int, columns: Int, band_cells: Int) -> Int:
    """What `vector_align` holds for a `rows` by `columns` pair whose band has `band_cells` cells, in bytes: a
    decision a cell and a step's spare lanes, the three indexes of every diagonal, and three diagonals of
    each of its three layers, at 32 bits."""
    var diagonals = rows + columns + 1
    return band_cells + 2 * 2 * WIDTH + diagonals * 3 * 8 + 9 * (rows + 2 + 4 * WIDTH) * 4


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
    """`serial_align`'s global alignment, swept by anti-diagonal a vector of cells a step, in lanes of 16 bits
    while the scores fit (see `anti_diagonals.fits_16_bits`), else 32.

    Only the cells on diagonals `low_diagonal ..= high_diagonal`, column minus row, are computed, the rest
    reading as unreachable; given a band holding every optimal path (see `optimal_band`), those cells hold
    the same values as `serial_align`'s. Each keeps only its decision, the byte `gotoh.decide` would give it
    from those values, stored anti-diagonal by anti-diagonal, each diagonal's band contiguous by row; the
    scores keep three diagonals. The walk back reads the decisions as `gotoh.walk` reads the scores, so it
    takes the same path, in a byte a cell where the scores took twelve.
    """
    var rows = len(first)
    var columns = len(second)
    if rows == 0 or columns == 0:
        return serial_align[AlignmentMode.GLOBAL](first, second, substitutions, alphabet_size, gaps, alphabet)
    var extremes = table_extremes(substitutions, alphabet_size)
    var substitution = -min(extremes[1], 0)
    var deficit = straight_deficit(substitution, -Int(gaps.open), -Int(gaps.extend), rows, columns)
    if fits_16_bits[False](max(extremes[0], 0), max(substitution, -Int(gaps.open)), deficit, rows, columns):
        return traced_band[DType.int16, 2 * WIDTH](first, second, lookup, gaps, alphabet, low_diagonal, high_diagonal)
    return traced_band[DType.int32, WIDTH](first, second, lookup, gaps, alphabet, low_diagonal, high_diagonal)


def traced_band[
    dtype: DType, width: Int
](
    first: List[UInt8],
    second: List[UInt8],
    lookup: SubstitutionLookup,
    gaps: AffineGapCosts,
    alphabet: String,
    low_diagonal: Int,
    high_diagonal: Int,
) -> GappedAlignment:
    """`vector_align` in `width` lanes of `dtype`, which hold every score of the band."""
    comptime Lanes = SIMD[dtype, width]
    comptime Value = Scalar[dtype]
    # Unreachable: below every score the lanes hold, with a move to spare (see `fits_16_bits`).
    comptime LOW = Value.MIN // 4
    var rows = len(first)
    var columns = len(second)
    var open = Int(gaps.open)
    var extend = Int(gaps.extend)
    var low_band = max(low_diagonal, -rows)
    var high_band = min(high_diagonal, columns)

    # Diagonal `d`'s cells lie on rows `lows[d] ..= highs[d]`: inside the matrix, and with `column - row`
    # inside the band, column being `d - row`. Its decisions start at `starts[d]`.
    var diagonals = rows + columns + 1
    var lows = List[Int](length=diagonals, fill=0)
    var highs = List[Int](length=diagonals, fill=0)
    var starts = List[Int](length=diagonals + 1, fill=0)
    for diagonal in range(diagonals):
        var low = max(0, diagonal - columns, -((high_band - diagonal) // 2))
        var high = min(rows, diagonal, (diagonal - low_band) // 2)
        lows[diagonal] = low
        highs[diagonal] = high
        starts[diagonal + 1] = starts[diagonal] + max(high - low + 1, 0)
    # A step's spare lanes run past a diagonal's last decision, the last diagonal's past the end.
    var decisions = List[UInt8](unsafe_uninit_length=starts[diagonals] + 2 * width)
    var letters = row_letters[width](Span(first))
    var reversed = column_letters[width](Span(second))

    # Three diagonals of each layer, a cell on row `i` at `i` within its diagonal's slot; the slots turn as
    # the sweep moves on.
    var span = rows + 2 + 2 * width
    var scores = List[Value](length=3 * span, fill=LOW)
    var deletes = List[Value](length=3 * span, fill=LOW)
    var inserts = List[Value](length=3 * span, fill=LOW)
    var score_cells = scores.unsafe_ptr()
    var delete_cells = deletes.unsafe_ptr()
    var insert_cells = inserts.unsafe_ptr()
    var decision_cells = decisions.unsafe_ptr()

    # Diagonal zero: the origin, which every band holds.
    score_cells[unsafe_offset=0] = 0
    delete_cells[unsafe_offset=0] = 0
    insert_cells[unsafe_offset=0] = 0
    var substitute = lookup.lanes[width, dtype]()
    var steps = GapLanes[dtype, width].symmetric(open, extend)
    var source_deleting = Lanes(Value(Int(Layer.DELETING.kind)))
    var source_inserting = Lanes(Value(Int(Layer.INSERTING.kind)))
    var deletion_extends = Lanes(Value(Int(GapRun.EXTENDS.kind) << 2))
    var insertion_extends = Lanes(Value(Int(GapRun.EXTENDS.kind) << 3))
    var aligning = Lanes(Value(Int(Layer.ALIGNING.kind)))
    var nothing = Lanes(0)
    var final_score = Value(0)
    for diagonal in range(1, diagonals):
        var first_row = lows[diagonal]
        var last_row = highs[diagonal]
        var here = (diagonal % 3) * span
        var one_back = ((diagonal - 1) % 3) * span
        var two_back = ((diagonal - 2) % 3) * span if diagonal >= 2 else 0
        var low = max(1, first_row)
        var high = min(last_row, diagonal - 1)
        var lag = columns - diagonal
        var decided = starts[diagonal] - first_row
        var row = low
        while row <= high:
            var above = score_cells.unsafe_offset(one_back + row - 1).unsafe_load[width=width]()
            var above_delete = delete_cells.unsafe_offset(one_back + row - 1).unsafe_load[width=width]()
            var left = score_cells.unsafe_offset(one_back + row).unsafe_load[width=width]()
            var left_insert = insert_cells.unsafe_offset(one_back + row).unsafe_load[width=width]()
            var above_left = score_cells.unsafe_offset(two_back + row - 1).unsafe_load[width=width]()
            var mine = letters.unsafe_ptr().unsafe_offset(row).unsafe_load[width=width]()
            var theirs = reversed.unsafe_ptr().unsafe_offset(lag + row).unsafe_load[width=width]()
            var aligned = above_left + substitute(mine, theirs)
            var deletion = gap_layer(above, above_delete, steps.down_first, steps.down_further)
            var insertion = gap_layer(left, left_insert, steps.across_first, steps.across_further)
            var score = best_move[False](aligned, deletion, insertion)
            # `gotoh.decide`'s byte: the source, the aligned move before a deletion before an insertion, and
            # each run extending only where extending scores strictly more than opening.
            var source = aligned.eq(score).select(
                aligning, deletion.eq(score).select(source_deleting, source_inserting)
            )
            var deleting = (above_delete + steps.down_further).gt(above + steps.down_first)
            var inserting = (left_insert + steps.across_further).gt(left + steps.across_first)
            var code = (
                source | deleting.select(deletion_extends, nothing) | inserting.select(insertion_extends, nothing)
            )
            score_cells.unsafe_offset(here + row).unsafe_store(score)
            delete_cells.unsafe_offset(here + row).unsafe_store(deletion)
            insert_cells.unsafe_offset(here + row).unsafe_store(insertion)
            decision_cells.unsafe_offset(decided + row).unsafe_store(code.cast[DType.uint8]())
            row += width
        # The border cells inside the band, then the cells either side of it, which the lanes may have run
        # over and the next two diagonals read as unreachable.
        if first_row == 0:
            var border = Value(Int(gaps.run(diagonal)))
            score_cells[unsafe_offset=here] = border
            delete_cells[unsafe_offset=here] = border + Value(open + extend)
            insert_cells[unsafe_offset=here] = 0
        else:
            score_cells[unsafe_offset=here + first_row - 1] = LOW
            delete_cells[unsafe_offset=here + first_row - 1] = LOW
            insert_cells[unsafe_offset=here + first_row - 1] = LOW
        if last_row == diagonal:
            var border = Value(Int(gaps.run(diagonal)))
            score_cells[unsafe_offset=here + diagonal] = border
            delete_cells[unsafe_offset=here + diagonal] = 0
            insert_cells[unsafe_offset=here + diagonal] = border + Value(open + extend)
        score_cells[unsafe_offset=here + last_row + 1] = LOW
        delete_cells[unsafe_offset=here + last_row + 1] = LOW
        insert_cells[unsafe_offset=here + last_row + 1] = LOW
        if diagonal == diagonals - 1:
            final_score = score_cells[unsafe_offset=here + rows]

    # The walk back, `gotoh.walk`'s over the decisions, to the first row or column; a global alignment's rest
    # gaps against the letters left.
    var symbols = alphabet.as_bytes()
    var top = List[UInt8]()
    var bottom = List[UInt8]()
    var row = rows
    var column = columns
    var layer = Layer.ALIGNING
    while row > 0 and column > 0:
        var diagonal = row + column
        var decision = Decision(decisions[starts[diagonal] + row - lows[diagonal]])
        var step = advance(layer, decision)
        top.append(symbols[Int(first[row - 1])] if step.row_advance != 0 else GAP_BYTE)
        bottom.append(symbols[Int(second[column - 1])] if step.column_advance != 0 else GAP_BYTE)
        row += step.row_advance
        column += step.column_advance
        layer = step.lands_in
    for letter in range(row, 0, -1):
        top.append(symbols[Int(first[letter - 1])])
        bottom.append(GAP_BYTE)
    for letter in range(column, 0, -1):
        top.append(GAP_BYTE)
        bottom.append(symbols[Int(second[letter - 1])])
    top.reverse()
    bottom.reverse()
    return GappedAlignment(Int32(final_score), String(unsafe_from_utf8=top), String(unsafe_from_utf8=bottom))


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
