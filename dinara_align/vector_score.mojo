"""
The host's affine-gap score by full sweep, sixteen cells at a time, for a table of one match and
one mismatch score.

The same Gotoh recurrence and borders as `serial_score`, so the same score, swept by anti-diagonal:
every cell of `d = row + column` reads only diagonals `d - 1` and `d - 2`, so a whole diagonal is
independent and fills vector lanes. Indexed by row, a diagonal's cells read the first sequence
forward and the second backward, so the second is stored reversed and both load contiguously; a
uniform table turns each substitution into one comparison.
"""

from .alignment import AffineGapCosts, AlignmentMode, AlignmentResult, AntiDiagonalMajor, reconstruct, serial_align

comptime WIDTH = 16
"""Cells a step computes at once."""

comptime Lanes = SIMD[DType.int32, WIDTH]


def uniform_table(substitutions: List[Scalar[DType.int8]], alphabet_size: Int) -> Optional[Tuple[Int, Int]]:
    """The match and mismatch scores of a table with one of each, or None."""
    if alphabet_size < 1:
        return None
    var reward = Int(substitutions[0])
    var mismatch = Int(substitutions[1]) if alphabet_size > 1 else reward - 1
    for row in range(alphabet_size):
        for column in range(alphabet_size):
            if Int(substitutions[row * alphabet_size + column]) != (reward if row == column else mismatch):
                return None
    return (reward, mismatch)


def vector_score[
    mode: AlignmentMode
](first: List[UInt8], second: List[UInt8], reward: Int, mismatch: Int, gaps: AffineGapCosts) -> Int32:
    """The optimal score under a uniform table, equal to `serial_score`'s."""
    var rows = len(first)
    var columns = len(second)
    var open = gaps.open
    var extend = gaps.extend

    @always_inline
    def border(length: Int) {imm open, imm extend} -> Int32:
        """The score of a gap run of `length` along the global border; the local border is zero."""
        comptime if mode == AlignmentMode.LOCAL:
            return 0
        if length == 0:
            return 0
        return open + Int32(length - 1) * extend

    if rows == 0 or columns == 0:
        return border(rows + columns)

    # Row `i` reads letter `i - 1` of the first sequence, so it is stored one place on; the second
    # is stored back to front. Both carry a step's width of padding past their ends.
    var letters = List[UInt8](length=rows + 1 + WIDTH, fill=0xFE)
    for index in range(rows):
        letters[index + 1] = first[index]
    var reversed = List[UInt8](length=columns + WIDTH, fill=0xFF)
    for index in range(columns):
        reversed[index] = second[columns - 1 - index]
    var size = rows + 1 + WIDTH
    var two_back = List[Int32](length=size, fill=0)
    var one_back = List[Int32](length=size, fill=0)
    var current = List[Int32](length=size, fill=0)
    var deletes_back = List[Int32](length=size, fill=0)
    var deletes = List[Int32](length=size, fill=0)
    var inserts_back = List[Int32](length=size, fill=0)
    var inserts = List[Int32](length=size, fill=0)
    # Diagonal one: the border cells beside the origin.
    one_back[0] = border(1)
    deletes_back[0] = one_back[0] + open + extend
    one_back[1] = border(1)
    inserts_back[1] = one_back[1] + open + extend

    var opening = Lanes(open)
    var extension = Lanes(extend)
    var matched = Lanes(Int32(reward))
    var mismatched = Lanes(Int32(mismatch))
    var zero = Lanes(0)
    var lane_rows = Lanes()
    comptime for lane in range(WIDTH):
        lane_rows[lane] = Int32(lane)
    var best = Lanes(0)
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
            var substitution = mine.eq(theirs).select(matched, mismatched)
            var deletion = max(above + opening, above_delete + extension)
            var insertion = max(left + opening, left_insert + extension)
            var score = max(above_left + substitution, max(deletion, insertion))
            comptime if mode == AlignmentMode.LOCAL:
                score = max(score, zero)
                # Lanes past the diagonal's last row hold no cell.
                best = max(best, (lane_rows + Int32(row)).le(Lanes(Int32(high))).select(score, zero))
            current.unsafe_ptr().unsafe_offset(row).unsafe_store(score)
            deletes.unsafe_ptr().unsafe_offset(row).unsafe_store(deletion)
            inserts.unsafe_ptr().unsafe_offset(row).unsafe_store(insertion)
            row += WIDTH
        # The border cells of this diagonal, written after the lanes that may have run over them.
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
    comptime if mode == AlignmentMode.LOCAL:
        return best.reduce_max()
    return one_back[rows]


def vector_align[
    mode: AlignmentMode
](
    first: List[UInt8],
    second: List[UInt8],
    reward: Int,
    mismatch: Int,
    gaps: AffineGapCosts,
    substitutions: List[Scalar[DType.int8]],
    alphabet_size: Int,
    alphabet: String,
) -> AlignmentResult:
    """`serial_align`'s alignment under a uniform table, its three layers swept sixteen cells at a time.

    The layers hold the same values at every cell, borders included, so `reconstruct` walks the same
    path; they are only stored anti-diagonal by anti-diagonal, each diagonal's cells contiguous by
    row, so a step loads its neighbours from the two diagonals before as vectors (see
    `AntiDiagonalMajor`). A local alignment starts where `serial_align`'s row-major scan first meets
    the best score: the cell earliest by row, then column, among those holding it.
    """
    var rows = len(first)
    var columns = len(second)
    if rows == 0 or columns == 0:
        return serial_align[mode](first, second, substitutions, alphabet_size, gaps, alphabet)
    var open = gaps.open
    var extend = gaps.extend

    @always_inline
    def border(length: Int) {imm open, imm extend} -> Int32:
        comptime if mode == AlignmentMode.LOCAL:
            return 0
        if length == 0:
            return 0
        return open + Int32(length - 1) * extend

    var diagonals = rows + columns + 1
    var starts = List[Int](length=diagonals + 1, fill=0)
    for diagonal in range(diagonals):
        var length = min(rows, diagonal) - max(0, diagonal - columns) + 1
        starts[diagonal + 1] = starts[diagonal] + length
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
    score_cells[0] = 0
    delete_cells[0] = 0
    insert_cells[0] = 0
    var opening = Lanes(open)
    var extension = Lanes(extend)
    var matched = Lanes(Int32(reward))
    var mismatched = Lanes(Int32(mismatch))
    var zero = Lanes(0)
    var lane_rows = Lanes()
    comptime for lane in range(WIDTH):
        lane_rows[lane] = Int32(lane)
    var best = Int32(0)
    var best_row = 0
    var best_column = 0
    for diagonal in range(1, diagonals):
        var first_row = max(0, diagonal - columns)
        var here = starts[diagonal] - first_row
        # Diagonal `d`'s cell on row `i` sits at `here + i`; its neighbours on the two before likewise.
        var one_back = starts[diagonal - 1] - max(0, diagonal - 1 - columns)
        var two_back = starts[diagonal - 2] - max(0, diagonal - 2 - columns) if diagonal >= 2 else 0
        var low = max(1, diagonal - columns)
        var high = min(rows, diagonal - 1)
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
            var substitution = mine.eq(theirs).select(matched, mismatched)
            var deletion = max(above + opening, above_delete + extension)
            var insertion = max(left + opening, left_insert + extension)
            var score = max(above_left + substitution, max(deletion, insertion))
            comptime if mode == AlignmentMode.LOCAL:
                score = max(score, zero)
                var counted = (lane_rows + Int32(row)).le(Lanes(Int32(high))).select(score, zero)
                var top = counted.reduce_max()
                if top > best or (top == best and top != 0):
                    # The earliest lane holding it is the earliest cell of this diagonal's in row order.
                    for lane in range(WIDTH):
                        if counted[lane] == top:
                            var cell_row = row + lane
                            var cell_column = diagonal - cell_row
                            if (
                                top > best
                                or cell_row < best_row
                                or (cell_row == best_row and cell_column < best_column)
                            ):
                                best = top
                                best_row = cell_row
                                best_column = cell_column
                            break
            score_cells.unsafe_offset(here + row).unsafe_store(score)
            delete_cells.unsafe_offset(here + row).unsafe_store(deletion)
            insert_cells.unsafe_offset(here + row).unsafe_store(insertion)
            row += WIDTH
        # The border cells, written after the lanes that may have run over the last of them.
        if diagonal <= columns:
            score_cells[here] = border(diagonal)
            delete_cells[here] = score_cells[here] + open + extend
            insert_cells[here] = 0
        if diagonal <= rows:
            score_cells[here + diagonal] = border(diagonal)
            delete_cells[here + diagonal] = 0
            insert_cells[here + diagonal] = score_cells[here + diagonal] + open + extend

    var start_row = rows
    var start_column = columns
    comptime if mode == AlignmentMode.LOCAL:
        start_row = best_row
        start_column = best_column
    var layout = AntiDiagonalMajor(starts.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](), columns)
    var reconstruction = reconstruct(
        scores,
        deletes,
        inserts,
        layout,
        first,
        second,
        substitutions,
        alphabet_size,
        start_row,
        start_column,
        alphabet,
        gaps,
        mode,
    )
    var final_score = scores[layout.index(start_row, start_column)]
    # The layout reads `starts` through a pointer, so the list must outlive every use of it.
    _ = len(starts)
    return AlignmentResult(final_score, reconstruction[0], reconstruction[1])
