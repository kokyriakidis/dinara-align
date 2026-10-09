# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The host's affine-gap sweep by anti-diagonal, a vector of cells a step, which every sweep under a score
runs.

Every cell of `d = row + column` reads only diagonals `d - 1` and `d - 2`, so a whole anti-diagonal is
independent and fills vector lanes. Indexed by row, a diagonal's cells read the rows' sequence forward and
the columns' backward, so the columns' is stored back to front and both load contiguously; each step
scores its pairs at once (see `substitutions`), in lanes of 16 bits while the scores fit (see
`fits_16_bits`). The sweeps differ in their borders, in what they watch each step for and in what they
keep, which stay theirs: `scored.swept_cells` a best end, from local or free starts; `vector_score.reach_back`
where a local alignment starts; `alignment.vector_sweep_bands` a half's last row. `vector_score.vector_align`
keeps every diagonal of its band rather than the last two, so it shares the letters and the recurrence.
"""

from .substitutions import SubstitutionLanes


@always_inline
def gotoh_lanes[
    dtype: DType, width: Int, local: Bool = False
](
    above_left: SIMD[dtype, width],
    above: SIMD[dtype, width],
    above_delete: SIMD[dtype, width],
    left: SIMD[dtype, width],
    left_insert: SIMD[dtype, width],
    substitution: SIMD[dtype, width],
    deletion_opening: SIMD[dtype, width],
    deletion_extension: SIMD[dtype, width],
    insertion_opening: SIMD[dtype, width],
    insertion_extension: SIMD[dtype, width],
) -> Tuple[SIMD[dtype, width], SIMD[dtype, width], SIMD[dtype, width]]:
    """`width` interior cells of the Gotoh recurrence, a gap's first letter scoring its `opening` and each
    further its `extension`, with the local clamp folded in at comptime: the score, the deletion layer and
    the insertion layer.

    This is the single transcription of the recurrence that AffineGaps' NumPy reference holds as the oracle;
    every sweep of it on the host and on the device goes through it, a cell or an anti-diagonal's lanes at
    a time. The device's kernels holding cells shifted by their anti-diagonal (see `score_groups`) take
    the same recurrence in three additions instead."""
    var deletion = max(above + deletion_opening, above_delete + deletion_extension)
    var insertion = max(left + insertion_opening, left_insert + insertion_extension)
    # The substitution, the slowest to come, waits on one maximum alone.
    var score = max(above_left + substitution, max(deletion, insertion))
    comptime if local:
        score = max(score, SIMD[dtype, width](0))
    return (score, deletion, insertion)


def fits_16_bits[floored: Bool](reward: Int, dearest: Int, rows: Int, columns: Int) -> Bool:
    """Whether every score of a sweep fits 16 bits, a pair earning at most `reward` and a move costing at
    most `dearest`. None passes the reward of the shorter sequence matched throughout. A `floored` sweep's,
    a local alignment's, none falls further below zero than the dearest single move, as each cell takes the
    best of its moves from cells of zero or more; any other's none falls below every letter of both paying
    the dearest move. Nor does a gap's sentinel, a quarter of the way down, fall further than one extension
    below it."""
    var fits = reward * (min(rows, columns) + 1) < 32000 and dearest < 4000
    comptime if not floored:
        fits = fits and dearest * (rows + columns + 1) < 8000
    return fits


@fieldwise_init
struct GapLanes[dtype: DType, width: Int](ImplicitlyCopyable):
    """What a gap scores in a sweep's lanes, its first letter and each further one, down the lanes, a gap
    of the rows' letters, and across them, of the columns', for each of two pieces."""

    var down_first: SIMD[Self.dtype, Self.width]
    var down_further: SIMD[Self.dtype, Self.width]
    var across_first: SIMD[Self.dtype, Self.width]
    var across_further: SIMD[Self.dtype, Self.width]
    var down_first2: SIMD[Self.dtype, Self.width]
    var down_further2: SIMD[Self.dtype, Self.width]
    var across_first2: SIMD[Self.dtype, Self.width]
    var across_further2: SIMD[Self.dtype, Self.width]

    @staticmethod
    def symmetric(first: Int, further: Int) -> Self:
        """One piece, the same both ways: a gap's first letter scoring `first`, each further `further`."""
        comptime Lanes = SIMD[Self.dtype, Self.width]
        return Self(
            Lanes(Scalar[Self.dtype](first)),
            Lanes(Scalar[Self.dtype](further)),
            Lanes(Scalar[Self.dtype](first)),
            Lanes(Scalar[Self.dtype](further)),
            Lanes(0),
            Lanes(0),
            Lanes(0),
            Lanes(0),
        )


def row_letters[width: Int](down: ImmSpan[UInt8, _], backward: Bool = False) -> List[UInt8]:
    """The rows' letters as a sweep of `width` lanes loads them, row `i`'s at `i`, from one, so the cells of a
    diagonal from row `i` on read them from `i`; padded past the end by code zero, a letter of every table,
    whose cells nothing reads. A `backward` sequence is read from its end, its last letter the first row."""
    var rows = len(down)
    var letters = List[UInt8](length=rows + 1 + width, fill=0)
    var to = letters.unsafe_ptr()
    var source = down.unsafe_ptr()
    # One loop each way, so each copies as a vector.
    if backward:
        for row in range(1, rows + 1):
            to[unsafe_offset=row] = source[unsafe_offset=rows - row]
    else:
        for row in range(1, rows + 1):
            to[unsafe_offset=row] = source[unsafe_offset=row - 1]
    return letters^


def column_letters[width: Int](across: ImmSpan[UInt8, _], backward: Bool = False) -> List[UInt8]:
    """The columns' letters as a sweep of `width` lanes loads them, column `j`'s at `len(across) - j`, so the
    cells of diagonal `d` from row `i` on read them from `len(across) - d + i`; padded as `row_letters` pads.
    A `backward` sequence is read from its end, its last letter the first column."""
    var columns = len(across)
    var letters = List[UInt8](length=columns + width, fill=0)
    var to = letters.unsafe_ptr()
    var source = across.unsafe_ptr()
    if backward:
        for index in range(columns):
            to[unsafe_offset=index] = source[unsafe_offset=index]
    else:
        for index in range(columns):
            to[unsafe_offset=index] = source[unsafe_offset=columns - 1 - index]
    return letters^


struct AntiDiagonals[dtype: DType, width: Int, pieces: Int = 1](Movable):
    """A sweep's letters (see `row_letters`, `column_letters`) and its cells on the diagonal it fills and the
    two before, each layer's cell on row `i` at `i`, with `width` places past the last row for a step's spare
    lanes.

    A sweep reads and writes them through `cells`: it fills a diagonal a step at a time (see
    `DiagonalCells.step`), writes its border cells after, which a step's spare lanes may have run over, and
    moves on (see `DiagonalCells.advance`)."""

    comptime Value = Scalar[Self.dtype]

    var letters: List[UInt8]
    var others: List[UInt8]
    var two_back: List[Self.Value]
    var one_back: List[Self.Value]
    var current: List[Self.Value]
    var deletes_back: List[Self.Value]
    """A gap of the rows' letters grows down the lanes, one of the columns' across diagonals."""
    var deletes: List[Self.Value]
    var inserts_back: List[Self.Value]
    var inserts: List[Self.Value]
    var deletes2_back: List[Self.Value]
    """The second piece's layers, empty with one."""
    var deletes2: List[Self.Value]
    var inserts2_back: List[Self.Value]
    var inserts2: List[Self.Value]

    def __init__(
        out self,
        down: ImmSpan[UInt8, _],
        across: ImmSpan[UInt8, _],
        gaps_fill: Self.Value,
        down_backward: Bool = False,
        across_backward: Bool = False,
    ):
        """A sweep of `down`'s letters as rows and `across`'s as columns, every cell zero and every gap
        layer's `gaps_fill`."""
        self.letters = row_letters[Self.width](down, down_backward)
        self.others = column_letters[Self.width](across, across_backward)
        var size = len(down) + 1 + Self.width
        comptime two = Self.pieces == 2
        self.two_back = List[Self.Value](length=size, fill=0)
        self.one_back = List[Self.Value](length=size, fill=0)
        self.current = List[Self.Value](length=size, fill=0)
        self.deletes_back = List[Self.Value](length=size, fill=gaps_fill)
        self.deletes = List[Self.Value](length=size, fill=gaps_fill)
        self.inserts_back = List[Self.Value](length=size, fill=gaps_fill)
        self.inserts = List[Self.Value](length=size, fill=gaps_fill)
        self.deletes2_back = List[Self.Value](length=size if two else 0, fill=gaps_fill)
        self.deletes2 = List[Self.Value](length=size if two else 0, fill=gaps_fill)
        self.inserts2_back = List[Self.Value](length=size if two else 0, fill=gaps_fill)
        self.inserts2 = List[Self.Value](length=size if two else 0, fill=gaps_fill)

    def cells(mut self) -> DiagonalCells[Self.dtype, Self.width, Self.pieces, origin_of(self)]:
        """The diagonals as a sweep's steps read and write them, held in registers, the sweep's buffers kept
        alive for as long as they are used (see `DiagonalCells`)."""
        return DiagonalCells[Self.dtype, Self.width, Self.pieces, origin_of(self)](self)


struct DiagonalCells[dtype: DType, width: Int, pieces: Int, origin: MutOrigin](
    ImplicitlyCopyable, TrivialRegisterPassable
):
    """An `AntiDiagonals`' letters and cells as its sweep reads and writes them, in registers rather than
    read through the sweep at every step: each layer's cell on row `i` of the diagonal being filled, and
    of the two before, `i` places on (see `AntiDiagonals.cells`)."""

    comptime Lanes = SIMD[Self.dtype, Self.width]
    comptime Cells = MutPointer[Scalar[Self.dtype], Self.origin]

    var letters: MutPointer[UInt8, Self.origin]
    var others: MutPointer[UInt8, Self.origin]
    var two_back: Self.Cells
    var one_back: Self.Cells
    var current: Self.Cells
    var deletes_back: Self.Cells
    var deletes: Self.Cells
    var inserts_back: Self.Cells
    var inserts: Self.Cells
    var deletes2_back: Self.Cells
    var deletes2: Self.Cells
    var inserts2_back: Self.Cells
    var inserts2: Self.Cells

    def __init__(out self, mut sweep: AntiDiagonals[Self.dtype, Self.width, Self.pieces]):
        """`sweep`'s buffers as they stand, its diagonal one the one before the first it fills."""
        self.letters = sweep.letters.unsafe_ptr().unsafe_origin_cast[Self.origin]()
        self.others = sweep.others.unsafe_ptr().unsafe_origin_cast[Self.origin]()
        self.two_back = sweep.two_back.unsafe_ptr().unsafe_origin_cast[Self.origin]()
        self.one_back = sweep.one_back.unsafe_ptr().unsafe_origin_cast[Self.origin]()
        self.current = sweep.current.unsafe_ptr().unsafe_origin_cast[Self.origin]()
        self.deletes_back = sweep.deletes_back.unsafe_ptr().unsafe_origin_cast[Self.origin]()
        self.deletes = sweep.deletes.unsafe_ptr().unsafe_origin_cast[Self.origin]()
        self.inserts_back = sweep.inserts_back.unsafe_ptr().unsafe_origin_cast[Self.origin]()
        self.inserts = sweep.inserts.unsafe_ptr().unsafe_origin_cast[Self.origin]()
        self.deletes2_back = sweep.deletes2_back.unsafe_ptr().unsafe_origin_cast[Self.origin]()
        self.deletes2 = sweep.deletes2.unsafe_ptr().unsafe_origin_cast[Self.origin]()
        self.inserts2_back = sweep.inserts2_back.unsafe_ptr().unsafe_origin_cast[Self.origin]()
        self.inserts2 = sweep.inserts2.unsafe_ptr().unsafe_origin_cast[Self.origin]()

    @always_inline
    def step[
        local: Bool = False, transposed: Bool = False
    ](
        self,
        row: Int,
        lag: Int,
        substitute: SubstitutionLanes[Self.width, Self.dtype, _],
        gaps: GapLanes[Self.dtype, Self.width],
    ) -> Self.Lanes:
        """The cells on rows `row ..< row + width` of the diagonal `columns - lag`, stored, their scores
        returned: a local alignment's floored at zero, and with `transposed` the table's row the columns'
        letter."""
        var above = self.one_back.unsafe_offset(row - 1).unsafe_load[width=Self.width]()
        var left = self.one_back.unsafe_offset(row).unsafe_load[width=Self.width]()
        var above_left = self.two_back.unsafe_offset(row - 1).unsafe_load[width=Self.width]()
        var above_delete = self.deletes_back.unsafe_offset(row - 1).unsafe_load[width=Self.width]()
        var left_insert = self.inserts_back.unsafe_offset(row).unsafe_load[width=Self.width]()
        var mine = self.letters.unsafe_offset(row).unsafe_load[width=Self.width]()
        var theirs = self.others.unsafe_offset(lag + row).unsafe_load[width=Self.width]()
        var pair = substitute(theirs, mine) if transposed else substitute(mine, theirs)
        var cell = gotoh_lanes[Self.dtype, Self.width, local and Self.pieces == 1](
            above_left,
            above,
            above_delete,
            left,
            left_insert,
            pair,
            gaps.down_first,
            gaps.down_further,
            gaps.across_first,
            gaps.across_further,
        )
        var score = cell[0]
        self.deletes.unsafe_offset(row).unsafe_store(cell[1])
        self.inserts.unsafe_offset(row).unsafe_store(cell[2])
        comptime if Self.pieces == 2:
            var deletion2 = max(
                above + gaps.down_first2,
                self.deletes2_back.unsafe_offset(row - 1).unsafe_load[width=Self.width]() + gaps.down_further2,
            )
            var insertion2 = max(
                left + gaps.across_first2,
                self.inserts2_back.unsafe_offset(row).unsafe_load[width=Self.width]() + gaps.across_further2,
            )
            score = max(score, max(deletion2, insertion2))
            comptime if local:
                score = max(score, Self.Lanes(0))
            self.deletes2.unsafe_offset(row).unsafe_store(deletion2)
            self.inserts2.unsafe_offset(row).unsafe_store(insertion2)
        self.current.unsafe_offset(row).unsafe_store(score)
        return score

    @always_inline
    def filled(self, row: Int) -> Scalar[Self.dtype]:
        """The cell on `row` of the diagonal filled last, which past its rows is anything."""
        return self.one_back[unsafe_offset=row]

    @always_inline
    def advance(mut self):
        """On to the next diagonal: the one filled is the one before it, and the oldest's buffers hold the
        next."""
        var oldest = self.two_back
        self.two_back = self.one_back
        self.one_back = self.current
        self.current = oldest
        var deletes = self.deletes_back
        self.deletes_back = self.deletes
        self.deletes = deletes
        var inserts = self.inserts_back
        self.inserts_back = self.inserts
        self.inserts = inserts
        comptime if Self.pieces == 2:
            var deletes2 = self.deletes2_back
            self.deletes2_back = self.deletes2
            self.deletes2 = deletes2
            var inserts2 = self.inserts2_back
            self.inserts2_back = self.inserts2
            self.inserts2 = inserts2
