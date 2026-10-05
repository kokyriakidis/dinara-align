"""
The host's affine-gap score by full sweep, sixteen cells at a time, for a table of one match and
one mismatch score.

The same Gotoh recurrence and borders as `serial_score`, so the same score, swept by anti-diagonal:
every cell of `d = row + column` reads only diagonals `d - 1` and `d - 2`, so a whole diagonal is
independent and fills vector lanes. Indexed by row, a diagonal's cells read the first sequence
forward and the second backward, so the second is stored reversed and both load contiguously; a
uniform table turns each substitution into one comparison.
"""

from .alignment import AffineGapCosts, AlignmentMode

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
