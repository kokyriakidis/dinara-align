# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The host's affine-gap sweeps under any substitution table, sixteen cells at a time.

The same Gotoh recurrence and borders as `alignment.serial_align`, so the same scores, swept by anti-diagonal:
every cell of `d = row + column` reads only diagonals `d - 1` and `d - 2`, so a whole diagonal is
independent and fills vector lanes. Indexed by row, a diagonal's cells read the first sequence
forward and the second backward, so the second is stored reversed and both load contiguously; each
step scores its sixteen pairs at once (see `substitutions`).

Two sweeps share the recurrence: the cell, swept back from a local alignment's end, where it earns its
score, which is where it starts (`reach_back`); and a global alignment's three layers, stored in the band of
diagonals its score bounds, to be traced (`vector_align`). A score alone is the tabled sweep's (see
`scoring.swept_score`).
"""

from .alignment import (
    AffineGapCosts,
    AlignmentMode,
    GappedAlignment,
    AntiDiagonalMajor,
    BAND_PADDING,
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
) -> Tuple[Int, Int]:
    """Where a local alignment ending at `(end_row, end_column)` and scoring `target`, the best there is,
    starts: the first cell, by anti-diagonal and then row, that a global alignment from it to the end
    scores `target` from.

    The sweep runs back from the end over both prefixes, as a global one anchored there: a cell's score
    is the best global alignment of the letters between it and the end. None can pass `target`, which no
    local alignment passes, and the start of one that earns it reaches it, so the first cell that does
    starts an optimal alignment, and the global alignment between it and the end is one. The sweep stops
    on the diagonal it finds it on, so it covers about the alignment's own cells, not the prefixes'.
    """
    var rows = end_row
    var columns = end_column
    var open = gaps.open
    var extend = gaps.extend

    @inline(.always)
    def border(length: Int) {imm open, imm extend} -> Int32:
        """The score of a gap run of `length` along the anchored border."""
        if length == 0:
            return 0
        return open + Int32(length - 1) * extend

    # Row `i` reads the `i`-th letter back from the end; a diagonal reads the second sequence's prefix
    # backward from the end, which stored for loads running forward is the prefix as it stands.
    var letters = List[UInt8](length=rows + 1 + WIDTH, fill=0xFE)
    for index in range(1, rows + 1):
        letters[index] = first[end_row - index]
    var reversed = List[UInt8](length=columns + WIDTH, fill=0xFF)
    for index in range(columns):
        reversed[index] = second[index]
    var size = rows + 1 + WIDTH
    var two_back = List[Int32](length=size, fill=0)
    var one_back = List[Int32](length=size, fill=0)
    var current = List[Int32](length=size, fill=0)
    var deletes_back = List[Int32](length=size, fill=0)
    var deletes = List[Int32](length=size, fill=0)
    var inserts_back = List[Int32](length=size, fill=0)
    var inserts = List[Int32](length=size, fill=0)
    one_back[0] = border(1)
    deletes_back[0] = one_back[0] + open + extend
    one_back[1] = border(1)
    inserts_back[1] = one_back[1] + open + extend

    var substitute = lookup.lanes[WIDTH]()
    var opening = Lanes(open)
    var extension = Lanes(extend)
    var wanted = Lanes(target)
    var lane_rows = Lanes()
    comptime for lane in range(WIDTH):
        lane_rows[lane] = Int32(lane)
    for diagonal in range(2, rows + columns + 1):
        var low = max(1, diagonal - columns)
        var high = min(rows, diagonal - 1)
        var lag = columns - diagonal
        var row = low
        while row <= high:
            var above = one_back.unsafe_ptr().unsafe_offset(row - 1).unsafe_load[width=WIDTH]()
            var above_delete = deletes_back.unsafe_ptr().unsafe_offset(row - 1).unsafe_load[width=WIDTH]()
            var left = one_back.unsafe_ptr().unsafe_offset(row).unsafe_load[width=WIDTH]()
            var left_insert = inserts_back.unsafe_ptr().unsafe_offset(row).unsafe_load[width=WIDTH]()
            var above_left = two_back.unsafe_ptr().unsafe_offset(row - 1).unsafe_load[width=WIDTH]()
            var mine = letters.unsafe_ptr().unsafe_offset(row).unsafe_load[width=WIDTH]()
            var theirs = reversed.unsafe_ptr().unsafe_offset(lag + row).unsafe_load[width=WIDTH]()
            var deletion = max(above + opening, above_delete + extension)
            var insertion = max(left + opening, left_insert + extension)
            var score = max(above_left + substitute(mine, theirs), max(deletion, insertion))
            var reached = score.ge(wanted) & (lane_rows + Int32(row)).le(Lanes(Int32(high)))
            if reached.reduce_or():
                # The lanes run by row, so the first that reached it is the diagonal's earliest.
                for lane in range(WIDTH):
                    if reached[lane]:
                        var back_rows = row + lane
                        return (end_row - back_rows, end_column - (diagonal - back_rows))
            current.unsafe_ptr().unsafe_offset(row).unsafe_store(score)
            deletes.unsafe_ptr().unsafe_offset(row).unsafe_store(deletion)
            inserts.unsafe_ptr().unsafe_offset(row).unsafe_store(insertion)
            row += WIDTH
        if diagonal <= columns:
            current[0] = border(diagonal)
            deletes[0] = current[0] + open + extend
        if diagonal <= rows:
            current[diagonal] = border(diagonal)
            inserts[diagonal] = current[diagonal] + open + extend
        swap(two_back, one_back)
        swap(one_back, current)
        swap(deletes_back, deletes)
        swap(inserts_back, inserts)
    # The whole prefixes, which a caller with a score above zero never needs.
    return (0, 0)


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
    def border(length: Int) {imm open, imm extend} -> Int32:
        """The score of a gap run of `length` along the global border."""
        if length == 0:
            return 0
        return open + Int32(length - 1) * extend

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
    var letters = List[UInt8](length=rows + 1 + WIDTH, fill=0xFE)
    for index in range(rows):
        letters[index + 1] = first[index]
    var reversed = List[UInt8](length=columns + WIDTH, fill=0xFF)
    for index in range(columns):
        reversed[index] = second[columns - 1 - index]

    var score_cells = scores.unsafe_ptr()
    var delete_cells = deletes.unsafe_ptr()
    var insert_cells = inserts.unsafe_ptr()

    @inline(.always)
    def unreachable(index: Int) {imm score_cells, imm delete_cells, imm insert_cells}:
        """Marks stored cell `index` unreachable in all three layers, as band padding reads."""
        score_cells[unsafe_offset=index] = UNREACHABLE
        delete_cells[unsafe_offset=index] = UNREACHABLE
        insert_cells[unsafe_offset=index] = UNREACHABLE

    # Diagonal zero: the origin, which every band holds.
    for index in range(BAND_PADDING):
        unreachable(index)
        unreachable(BAND_PADDING + 1 + index)
    score_cells[unsafe_offset=BAND_PADDING] = 0
    delete_cells[unsafe_offset=BAND_PADDING] = 0
    insert_cells[unsafe_offset=BAND_PADDING] = 0
    var substitute = lookup.lanes[WIDTH]()
    var opening = Lanes(open)
    var extension = Lanes(extend)
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
            var deletion = max(above + opening, above_delete + extension)
            var insertion = max(left + opening, left_insert + extension)
            var score = max(above_left + substitute(mine, theirs), max(deletion, insertion))
            score_cells.unsafe_offset(here + row).unsafe_store(score)
            delete_cells.unsafe_offset(here + row).unsafe_store(deletion)
            insert_cells.unsafe_offset(here + row).unsafe_store(insertion)
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
    var reconstruction = reconstruct(
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
        AlignmentMode.GLOBAL,
    )
    var final_score = scores[layout.index(rows, columns)]
    # The layout reads `starts` and `lows` through pointers, so both must outlive every use of it.
    _ = len(starts)
    _ = len(lows)
    _ = len(lookup.cells)
    return GappedAlignment(final_score, reconstruction[0], reconstruction[1])


comptime UNREACHABLE = Int32(-(1 << 28))
"""A cell outside the band: far enough below any score that a few gap costs keep it there."""


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
    var per_gap = 2 * Int(-gaps.extend) + reward
    var cost = reward * (rows + columns) - 2 * score
    var reach = max(0, (cost - per_gap * abs(difference)) // (2 * per_gap)) if per_gap > 0 else rows + columns
    return (min(0, difference) - reach, max(0, difference) + reach)
