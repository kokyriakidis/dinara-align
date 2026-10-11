# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Gotoh's affine-gap alignment under a substitution table, on the host.

Each cell of the matrix, `row` letters of the first sequence against `column` of the second, keeps three
scores (Gotoh 1982): the best alignment of those prefixes, the best ending in a deletion, a run of the
first sequence's letters against gaps, and the best ending in an insertion, a run of the second's. A gap
of `k` letters scores `open + (k - 1) extend`, both negative.

`serial_align` stores the three layers of the whole matrix and walks back through them. `linear_path`
keeps two rows: it splits a rectangle at its middle row, sweeps each half toward the cut, joins them
where an optimal path crosses (Hirschberg 1975), in the aligning layer or inside a deletion that
straddles the cut (Myers and Miller 1988), and solves rectangles small enough outright.

Of equally good moves into a cell the walk takes aligning, then deleting, then inserting, and a gap run
opening before one extending: every walk here, and the device's, decides a cell by `decide` and moves by
`advance`, so they agree. A path walks all three layers, so the score it reports is the one it realizes.
"""

from .anti_diagonals import AntiDiagonals, GapLanes, fits_16_bits, gotoh_lanes, straight_deficit
from .cigar import CigarWriter
from .common import GAP_BYTE, NEGATIVE_INFINITY, SubstitutionDType, SymbolDType
from .errors import AlignmentError, ErrorKind
from .substitutions import SubstitutionLookup, table_extremes

# region The recurrence


@fieldwise_init
struct AlignmentMode(Equatable, ImplicitlyCopyable, TrivialRegisterPassable):
    """Global, both sequences end to end, or local, any part of each, every cell's score floored at zero."""

    var kind: UInt8

    comptime GLOBAL = Self(0)
    comptime LOCAL = Self(1)


@fieldwise_init
struct AffineGapCosts(ImplicitlyCopyable, TrivialRegisterPassable):
    """A gap's scores: `open` for its first letter, `extend` for each further one, neither above zero."""

    var open: Int32
    var extend: Int32

    @staticmethod
    def checked(open: Int32, extend: Int32) raises AlignmentError -> Self:
        """These scores, refused when a gap's first letter would score more than a further one, or a letter
        of a gap would score more than nothing."""
        if open > extend:
            raise AlignmentError(ErrorKind.INVALID_SCORING, "a gap's first letter scoring above a further one")
        if extend > 0:
            raise AlignmentError(ErrorKind.INVALID_SCORING, "a gap that earns")
        return Self(open, extend)

    @always_inline
    def run(self, length: Int) -> Int32:
        """A gap of `length` letters, nothing for none."""
        return 0 if length == 0 else self.open + Int32(length - 1) * self.extend


@fieldwise_init
struct Layer(Equatable, ImplicitlyCopyable, TrivialRegisterPassable):
    """One of a cell's three scores, or the one a walk back is in: aligning a pair of letters, deleting
    (inside a run of the first sequence's letters) or inserting (the second's)."""

    var kind: UInt8

    comptime ALIGNING = Self(0)
    comptime DELETING = Self(1)
    comptime INSERTING = Self(2)


@fieldwise_init
struct GapRun(Equatable, ImplicitlyCopyable, TrivialRegisterPassable):
    """Whether a gap run is opened at a cell or extends one already open: a run opened on a tie."""

    var kind: UInt8

    comptime OPENS = Self(0)
    comptime EXTENDS = Self(1)


@fieldwise_init
struct Cell(ImplicitlyCopyable, TrivialRegisterPassable):
    """A cell's three scores: its best, its best ending in a deletion, and in an insertion."""

    var score: Int32
    var deletion: Int32
    var insertion: Int32


@fieldwise_init
struct Step(ImplicitlyCopyable, TrivialRegisterPassable):
    """A step back: `-1` or `0` along each sequence, a `0` a gap in it, and the layer it lands in."""

    var row_advance: Int
    var column_advance: Int
    var lands_in: Layer


@always_inline
def gotoh_cell[
    mode: AlignmentMode
](
    above_left: Int32,
    above: Int32,
    above_delete: Int32,
    left: Int32,
    left_insert: Int32,
    substitution: Int32,
    gaps: AffineGapCosts,
) -> Cell:
    """A cell from its neighbours' scores and its pair's `substitution` (see `anti_diagonals.gotoh_lanes`)."""
    var found = gotoh_lanes[DType.int32, 1, mode == AlignmentMode.LOCAL](
        above_left, above, above_delete, left, left_insert, substitution, gaps.open, gaps.extend, gaps.open, gaps.extend
    )
    return Cell(found[0], found[1], found[2])


@always_inline
def beats(score: Int32, place: Int64, best: Int32, best_place: Int64) -> Bool:
    """Whether a local alignment's end scoring `score` at `place` displaces the best so far: a higher score,
    or an equal one above zero at an earlier place. The host's scans and the device's reductions all ask
    this, so they pick the same end."""
    return score > best or (score == best and score != 0 and place < best_place)


@fieldwise_init
struct Decision(ImplicitlyCopyable, TrivialRegisterPassable):
    """What a walk back needs of a cell, in four bits, as the device packs them eight to a word: the layer
    its score came from in the low two, `ENDS` where a local path starts, and whether its deletion and its
    insertion extend a run, bits 2 and 3."""

    var code: UInt8

    comptime ENDS = UInt8(3)
    """The low two bits of a local cell whose score fell to zero, where its path starts."""

    @staticmethod
    @always_inline
    def of(source: UInt8, deletion: GapRun, insertion: GapRun) -> Self:
        """A decision from its source's two bits and its two runs."""
        return Self(source | (deletion.kind << 2) | (insertion.kind << 3))

    @always_inline
    def ends(self) -> Bool:
        """Whether a local path starts here, its score fallen to zero."""
        return (self.code & 3) == Self.ENDS

    @always_inline
    def source(self) -> Layer:
        """The layer the cell's score came from, which a walk asks only of a cell where no path starts."""
        return Layer(self.code & 3)

    @always_inline
    def deletion(self) -> GapRun:
        return GapRun((self.code >> 2) & 1)

    @always_inline
    def insertion(self) -> GapRun:
        return GapRun((self.code >> 3) & 1)


@always_inline
def decide[mode: AlignmentMode](cell: Cell, aligned: Int32, above: Cell, left: Cell, gaps: AffineGapCosts) -> Decision:
    """`cell`'s decision, from its scores, `aligned` (the cell above and left plus its pair's score), and the
    cells above and to the left. A run extends only where extending scores strictly more than opening."""
    var deletion = GapRun.EXTENDS if above.deletion + gaps.extend > above.score + gaps.open else GapRun.OPENS
    var insertion = GapRun.EXTENDS if left.insertion + gaps.extend > left.score + gaps.open else GapRun.OPENS
    var source: UInt8
    if mode == AlignmentMode.LOCAL and cell.score <= 0:
        source = Decision.ENDS
    elif cell.score == aligned:
        source = Layer.ALIGNING.kind
    elif cell.score == cell.deletion:
        source = Layer.DELETING.kind
    else:
        source = Layer.INSERTING.kind
    return Decision.of(source, deletion, insertion)


@always_inline
def advance(state: Layer, decision: Decision) -> Step:
    """The step back from a cell whose decision is `decision`, in layer `state`. From the aligning layer it
    follows the cell's source, entering a gap run and taking its first step at once."""
    var layer = decision.source() if state == Layer.ALIGNING else state
    if layer == Layer.DELETING:
        return Step(-1, 0, Layer.DELETING if decision.deletion() == GapRun.EXTENDS else Layer.ALIGNING)
    if layer == Layer.INSERTING:
        return Step(0, -1, Layer.INSERTING if decision.insertion() == GapRun.EXTENDS else Layer.ALIGNING)
    return Step(-1, -1, Layer.ALIGNING)


@fieldwise_init
struct GappedAlignment(Copyable, Movable):
    """An alignment as its score and its two gapped rows, one column a letter or a gap of each."""

    var score: Int32
    var first_gapped: String
    var second_gapped: String

    def cigar(self, eqx: Bool = True) -> String:
        """The rows as a CIGAR, the first sequence the reference: `=` a match and `X` a substitution, or `M`
        for both without `eqx`, `D` a letter of the first alone and `I` one of the second."""
        var top = self.first_gapped.as_bytes()
        var bottom = self.second_gapped.as_bytes()
        # A run of `n` columns writes its digits and a letter, no more than `2 n` bytes.
        var writer = CigarWriter(2 * len(top) + 1)
        for column in range(len(top)):
            var letter = UInt8(ord("M"))
            if top[column] == GAP_BYTE:
                letter = UInt8(ord("I"))
            elif bottom[column] == GAP_BYTE:
                letter = UInt8(ord("D"))
            elif eqx:
                letter = UInt8(ord("=")) if top[column] == bottom[column] else UInt8(ord("X"))
            writer.add(letter, 1)
        return writer^.finish()


# endregion The recurrence

# region The whole matrix


trait CellLayout(ImplicitlyCopyable):
    """Where cell `(row, column)` of stored layers sits in each layer."""

    def index(self, row: Int, column: Int) -> Int:
        ...


@fieldwise_init
struct RowMajor(CellLayout, TrivialRegisterPassable):
    """Row after row, `stride` cells each."""

    var stride: Int

    @always_inline
    def index(self, row: Int, column: Int) -> Int:
        return row * self.stride + column


@fieldwise_init
struct RollingRows(CellLayout, TrivialRegisterPassable):
    """Two rows of `stride` cells taking turns, a linear-space sweep's: the row being filled and the one above."""

    var stride: Int

    @always_inline
    def index(self, row: Int, column: Int) -> Int:
        return (row & 1) * self.stride + column


comptime BAND_PADDING = 2
"""Cells kept either side of each anti-diagonal's band in `AntiDiagonalMajor`, read as unreachable, so a step
and a walk read past a band's edge unchecked."""


@fieldwise_init
struct AntiDiagonalMajor(CellLayout, TrivialRegisterPassable):
    """Anti-diagonal after anti-diagonal, each a band of rows from `lows[d]`, its padding from `starts[d]`."""

    var starts: MutPointer[Int, MutUntrackedOrigin]
    var lows: MutPointer[Int, MutUntrackedOrigin]

    @always_inline
    def index(self, row: Int, column: Int) -> Int:
        var diagonal = row + column
        return self.starts[unsafe_offset=diagonal] + BAND_PADDING + row - self.lows[unsafe_offset=diagonal]


@fieldwise_init
struct SweepHalf(Equatable, ImplicitlyCopyable, TrivialRegisterPassable):
    """Which way a sweep reads its rectangle: from its top left, or from its bottom right with both sequences
    read backward, which is how a linear-space split sweeps its lower half toward the cut."""

    var kind: UInt8

    comptime FORWARD = Self(0)
    comptime REVERSE = Self(1)


def fill_rows[
    mode: AlignmentMode, half: SweepHalf, Layout: CellLayout
](
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    gaps: AffineGapCosts,
    top: GapRun,
    layout: Layout,
    scores: MutSpan[Int32, _],
    deletes: MutSpan[Int32, _],
    inserts: MutSpan[Int32, _],
) -> Tuple[Int32, Int, Int]:
    """The three layers of `first`'s letters (the rows) against `second`'s (the columns), row by row where
    `layout` keeps them, read from their ends for the `REVERSE` half; for a local alignment the best cell,
    the first by row then column, as its score, row and column.

    A global alignment's borders are gaps from the corner, the first column's opened above it already when
    `top` extends a run into it. Its first row's deletion layer reads one more opening, so a deletion from
    there pays one; nothing reads the first column's deletion layer or the first row's insertion layer.
    A local alignment's borders are zero."""
    var rows = len(first)
    var columns = len(second)
    comptime local = mode == AlignmentMode.LOCAL
    comptime backward = half == SweepHalf.REVERSE
    var opening = gaps.open + gaps.extend
    # Every index below lies inside the layers `layout` spans: unchecked.
    var score_at = scores.unsafe_ptr()
    var delete_at = deletes.unsafe_ptr()
    var insert_at = inserts.unsafe_ptr()
    var table = substitutions.unsafe_ptr()
    var row_letters = first.unsafe_ptr()
    var column_letters = second.unsafe_ptr()

    for column in range(columns + 1):
        var cell = layout.index(0, column)
        var border = Int32(0) if local else gaps.run(column)
        score_at[unsafe_offset=cell] = border
        delete_at[unsafe_offset=cell] = border + opening
        insert_at[unsafe_offset=cell] = 0

    var best = Int32(0)
    var best_row = 0
    var best_column = 0
    for row in range(1, rows + 1):
        var border = Int32(0)
        comptime if not local:
            border = Int32(row) * gaps.extend if top == GapRun.EXTENDS else gaps.run(row)
        var first_cell = layout.index(row, 0)
        score_at[unsafe_offset=first_cell] = border
        delete_at[unsafe_offset=first_cell] = border
        insert_at[unsafe_offset=first_cell] = border + opening
        var table_row = table.unsafe_offset(
            Int(row_letters[unsafe_offset=rows - row if backward else row - 1]) * alphabet_size
        )
        for column in range(1, columns + 1):
            var cell = layout.index(row, column)
            var above = layout.index(row - 1, column)
            var letter = Int(column_letters[unsafe_offset=columns - column if backward else column - 1])
            var filled = gotoh_cell[mode](
                score_at[unsafe_offset=above - 1],
                score_at[unsafe_offset=above],
                delete_at[unsafe_offset=above],
                score_at[unsafe_offset=cell - 1],
                insert_at[unsafe_offset=cell - 1],
                Int32(table_row[unsafe_offset=letter]),
                gaps,
            )
            score_at[unsafe_offset=cell] = filled.score
            delete_at[unsafe_offset=cell] = filled.deletion
            insert_at[unsafe_offset=cell] = filled.insertion
            comptime if local:
                if filled.score > best:
                    best = filled.score
                    best_row = row
                    best_column = column
    return (best, best_row, best_column)


def walk[
    mode: AlignmentMode, Layout: CellLayout, F: def(Int, Int, Step) -> None
](
    scores: ImmSpan[Int32, _],
    deletes: ImmSpan[Int32, _],
    inserts: ImmSpan[Int32, _],
    layout: Layout,
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    gaps: AffineGapCosts,
    from_row: Int,
    from_column: Int,
    state: Layer,
    taken: F,
) -> Tuple[Int, Int]:
    """Walks back over stored layers from `(from_row, from_column)` in `state`, handing each step to `taken`
    with the cell it leaves, to the first row or column, or for a local alignment to the cell its path
    starts at; returns where it stopped."""
    var row = from_row
    var column = from_column
    var layer = state

    @always_inline
    def cell_at(index: Int) {imm scores, imm deletes, imm inserts} -> Cell:
        return Cell(scores[index], deletes[index], inserts[index])

    while row > 0 and column > 0:
        var here = layout.index(row, column)
        var pair = Int32(substitutions[Int(first[row - 1]) * alphabet_size + Int(second[column - 1])])
        var decision = decide[mode](
            cell_at(here),
            scores[layout.index(row - 1, column - 1)] + pair,
            cell_at(layout.index(row - 1, column)),
            cell_at(layout.index(row, column - 1)),
            gaps,
        )
        comptime if mode == AlignmentMode.LOCAL:
            if layer == Layer.ALIGNING and decision.ends():
                break
        var step = advance(layer, decision)
        taken(row, column, step)
        row += step.row_advance
        column += step.column_advance
        layer = step.lands_in
    return (row, column)


def reconstruct[
    mode: AlignmentMode, Layout: CellLayout
](
    scores: ImmSpan[Int32, _],
    deletes: ImmSpan[Int32, _],
    inserts: ImmSpan[Int32, _],
    layout: Layout,
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    from_row: Int,
    from_column: Int,
    alphabet: String,
    gaps: AffineGapCosts,
) -> Tuple[String, String]:
    """The gapped rows of the path `walk` takes back from `(from_row, from_column)`; a global one's rest, to
    the corner, gaps against the letters left."""
    var letters = alphabet.as_bytes()
    var top = List[UInt8]()
    var bottom = List[UInt8]()

    def taken(row: Int, column: Int, step: Step) {mut top, mut bottom, imm letters, imm first, imm second}:
        top.append(letters[Int(first[row - 1])] if step.row_advance != 0 else GAP_BYTE)
        bottom.append(letters[Int(second[column - 1])] if step.column_advance != 0 else GAP_BYTE)

    var stop = walk[mode](
        scores,
        deletes,
        inserts,
        layout,
        first,
        second,
        substitutions,
        alphabet_size,
        gaps,
        from_row,
        from_column,
        Layer.ALIGNING,
        taken,
    )
    comptime if mode == AlignmentMode.GLOBAL:
        for row in range(stop[0], 0, -1):
            top.append(letters[Int(first[row - 1])])
            bottom.append(GAP_BYTE)
        for column in range(stop[1], 0, -1):
            top.append(GAP_BYTE)
            bottom.append(letters[Int(second[column - 1])])
    top.reverse()
    bottom.reverse()
    return (String(unsafe_from_utf8=top), String(unsafe_from_utf8=bottom))


def serial_align[
    mode: AlignmentMode
](
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    gaps: AffineGapCosts,
    alphabet: String,
) -> GappedAlignment:
    """An optimal alignment over the whole matrix's three layers: a global one from the corner, a local one
    from its best cell."""
    var layout = RowMajor(len(second) + 1)
    var cells = (len(first) + 1) * (len(second) + 1)
    # The fill writes every cell, borders included.
    var scores = List[Int32](unsafe_uninit_length=cells)
    var deletes = List[Int32](unsafe_uninit_length=cells)
    var inserts = List[Int32](unsafe_uninit_length=cells)
    var best = fill_rows[mode, SweepHalf.FORWARD](
        first, second, substitutions, alphabet_size, gaps, GapRun.OPENS, layout, scores, deletes, inserts
    )
    var end_row = best[1] if mode == AlignmentMode.LOCAL else len(first)
    var end_column = best[2] if mode == AlignmentMode.LOCAL else len(second)
    var rows = reconstruct[mode](
        scores,
        deletes,
        inserts,
        layout,
        first,
        second,
        substitutions,
        alphabet_size,
        end_row,
        end_column,
        alphabet,
        gaps,
    )
    return GappedAlignment(scores[layout.index(end_row, end_column)], rows[0], rows[1])


# endregion The whole matrix

# region Linear space


@fieldwise_init
struct Rectangle(ImplicitlyCopyable, TrivialRegisterPassable):
    """Rows `row_from ..< row_to` of the first sequence against columns `column_from ..< column_to` of the second."""

    var row_from: Int
    var row_to: Int
    var column_from: Int
    var column_to: Int

    @always_inline
    def height(self) -> Int:
        return self.row_to - self.row_from

    @always_inline
    def width(self) -> Int:
        return self.column_to - self.column_from

    @always_inline
    def middle(self) -> Int:
        """The row a split cuts at."""
        return (self.row_from + self.row_to) // 2


@fieldwise_init
struct Frame(ImplicitlyCopyable, TrivialRegisterPassable):
    """A rectangle a linear-space traceback still has to solve, and whether a deletion run is already open
    across its top edge, and across its bottom edge, where a split cut inside one."""

    var area: Rectangle
    var top: GapRun
    var bottom: GapRun

    @always_inline
    def solved_outright(self, leaf_cells: Int) -> Bool:
        """Whether this frame is solved over its whole matrix rather than split: too narrow or too short to
        split, or holding no more than `leaf_cells` cells."""
        var area = self.area
        return area.width() == 0 or area.height() <= 2 or (area.height() + 1) * (area.width() + 1) <= leaf_cells

    def halves(self, cut: Cut) -> Tuple[Frame, Frame]:
        """The two frames on either side of `cut` at the middle row, the upper first."""
        var area = self.area
        var middle = area.middle()
        var column = area.column_from + cut.offset
        var shared = GapRun.EXTENDS if cut.in_deletion else GapRun.OPENS
        return (
            Frame(Rectangle(area.row_from, middle, area.column_from, column), self.top, shared),
            Frame(Rectangle(middle, area.row_to, column, area.column_to), shared, self.bottom),
        )


@fieldwise_init
struct Cut(ImplicitlyCopyable, TrivialRegisterPassable):
    """Where an optimal path crosses a split's middle row: `offset` columns into its frame, and whether
    inside a deletion that straddles it."""

    var offset: Int
    var in_deletion: Bool

    @staticmethod
    @always_inline
    def choosing(aligned: Int32, aligned_offset: Int, straddling: Int32, straddling_offset: Int) -> Self:
        """The better of the two crossings, the aligning layer's on a tie."""
        if straddling > aligned:
            return Self(straddling_offset, True)
        return Self(aligned_offset, False)


def best_cut(
    forward_scores: ImmSpan[Int32, _],
    forward_deletes: ImmSpan[Int32, _],
    reverse_scores: ImmSpan[Int32, _],
    reverse_deletes: ImmSpan[Int32, _],
    width: Int,
    gaps: AffineGapCosts,
) -> Cut:
    """The cut joining the two halves' last rows best, the first column on a tie: either the path crosses in
    the aligning layer, or a deletion straddles the cut, whose opening both halves charged, one of them
    refunded."""
    var aligned = NEGATIVE_INFINITY
    var aligned_offset = 0
    var straddling = NEGATIVE_INFINITY
    var straddling_offset = 0
    var refund = gaps.extend - gaps.open
    for offset in range(width + 1):
        var through = forward_scores[offset] + reverse_scores[width - offset]
        if through > aligned:
            aligned = through
            aligned_offset = offset
        var across = forward_deletes[offset] + reverse_deletes[width - offset] + refund
        if across > straddling:
            straddling = across
            straddling_offset = offset
    return Cut.choosing(aligned, aligned_offset, straddling, straddling_offset)


struct RowPath(Movable):
    """An alignment as the column it leaves each row at, `columns[row]` for `row` in `0 ... rows`, and the layer
    it enters each row in, `layers[row]` for `row` in `1 ... rows`: aligning, or deleting the row's letter."""

    var columns: List[Int32]
    var layers: List[Layer]

    def __init__(out self, rows: Int):
        """A path over `rows` rows, every row's place still to be written."""
        self.columns = List[Int32](length=rows + 1, fill=0)
        self.layers = List[Layer](length=rows + 1, fill=Layer.ALIGNING)

    def view(mut self) -> PathView:
        """Where frames write their rows' places, each its own rows, from any thread."""
        return PathView(
            self.columns.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
            self.layers.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin](),
        )

    def score(
        self,
        first: ImmSpan[Scalar[SymbolDType], _],
        second: ImmSpan[Scalar[SymbolDType], _],
        substitutions: ImmSpan[Scalar[SubstitutionDType], _],
        alphabet_size: Int,
        gaps: AffineGapCosts,
    ) -> Int32:
        """What the path over all of `first`'s rows scores, read off it in one pass."""
        var total = Int32(0)
        var deleting = False
        var inserting = False

        @always_inline
        def insert(mut total: Int32, mut deleting: Bool, mut inserting: Bool) {imm gaps}:
            total += gaps.extend if inserting else gaps.open
            inserting = True
            deleting = False

        for _ in range(Int(self.columns[0])):
            insert(total, deleting, inserting)
        for row in range(1, len(first) + 1):
            var column = Int(self.columns[row - 1])
            if self.layers[row] == Layer.DELETING:
                total += gaps.extend if deleting else gaps.open
                deleting = True
                inserting = False
            else:
                total += Int32(substitutions[Int(first[row - 1]) * alphabet_size + Int(second[column])])
                column += 1
                deleting = False
                inserting = False
            for _ in range(column, Int(self.columns[row])):
                insert(total, deleting, inserting)
        return total

    def gapped(
        self,
        first: ImmSpan[Scalar[SymbolDType], _],
        second: ImmSpan[Scalar[SymbolDType], _],
        alphabet: String,
        mode: AlignmentMode,
        from_row: Int,
        to_row: Int,
    ) -> Tuple[String, String]:
        """The gapped rows of the path over rows `from_row + 1 ... to_row`; a global path's first starts with
        the insertions before its first row."""
        var letters = alphabet.as_bytes()
        var top = List[UInt8]()
        var bottom = List[UInt8]()
        if mode == AlignmentMode.GLOBAL:
            for column in range(Int(self.columns[0])):
                top.append(GAP_BYTE)
                bottom.append(letters[Int(second[column])])
        for row in range(from_row + 1, to_row + 1):
            var column = Int(self.columns[row - 1])
            top.append(letters[Int(first[row - 1])])
            if self.layers[row] == Layer.DELETING:
                bottom.append(GAP_BYTE)
            else:
                bottom.append(letters[Int(second[column])])
                column += 1
            for inserted in range(column, Int(self.columns[row])):
                top.append(GAP_BYTE)
                bottom.append(letters[Int(second[inserted])])
        return (String(unsafe_from_utf8=top), String(unsafe_from_utf8=bottom))


@fieldwise_init
struct PathView(ImplicitlyCopyable, TrivialRegisterPassable):
    """A `RowPath`'s rows, written in place."""

    var columns: MutPointer[Int32, MutUntrackedOrigin]
    var layers: MutPointer[Layer, MutUntrackedOrigin]

    @always_inline
    def leaves(self, row: Int, column: Int):
        """The path leaves `row` at `column`."""
        self.columns[unsafe_offset=row] = Int32(column)

    @always_inline
    def enters(self, row: Int, layer: Layer):
        """The path enters `row` in `layer`."""
        self.layers[unsafe_offset=row] = layer


def solve_outright(
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    frame: Frame,
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    gaps: AffineGapCosts,
    path: PathView,
):
    """Solves `frame` over its whole matrix and writes its rows' places on `path`: where it leaves each of its
    rows but the last, and the layer it enters each but the first, so frames side by side share no row."""
    var area = frame.area
    var layout = RowMajor(area.width() + 1)
    var cells = (area.height() + 1) * (area.width() + 1)
    var scores = List[Int32](unsafe_uninit_length=cells)
    var deletes = List[Int32](unsafe_uninit_length=cells)
    var inserts = List[Int32](unsafe_uninit_length=cells)
    var rows = first[area.row_from : area.row_to]
    var columns = second[area.column_from : area.column_to]
    _ = fill_rows[AlignmentMode.GLOBAL, SweepHalf.FORWARD](
        rows, columns, substitutions, alphabet_size, gaps, frame.top, layout, scores, deletes, inserts
    )

    # A step along a row writes nothing: only a step down a row does.
    def record(row: Int, column: Int, step: Step) {imm path, imm area}:
        if step.row_advance == 0:
            return
        path.enters(area.row_from + row, Layer.ALIGNING if step.column_advance != 0 else Layer.DELETING)
        path.leaves(area.row_from + row - 1, area.column_from + column + step.column_advance)

    # A deletion open across the bottom edge was the split's finding: the walk leaves deleting.
    var stop = walk[AlignmentMode.GLOBAL](
        scores,
        deletes,
        inserts,
        layout,
        rows,
        columns,
        substitutions,
        alphabet_size,
        gaps,
        area.height(),
        area.width(),
        Layer.DELETING if frame.bottom == GapRun.EXTENDS else Layer.ALIGNING,
        record,
    )
    for row in range(stop[0], 0, -1):
        path.enters(area.row_from + row, Layer.DELETING)
        path.leaves(area.row_from + row - 1, area.column_from + stop[1])


struct RollingBands(Movable):
    """Two rows of each of the three layers, which a scalar sweep fills in turn (see `RollingRows`)."""

    var scores: List[Int32]
    var deletes: List[Int32]
    var inserts: List[Int32]

    def __init__(out self, columns: Int):
        """Room for frames up to `columns` wide; every entry is written before it is read."""
        self.scores = List[Int32](unsafe_uninit_length=2 * (columns + 1))
        self.deletes = List[Int32](unsafe_uninit_length=2 * (columns + 1))
        self.inserts = List[Int32](unsafe_uninit_length=2 * (columns + 1))


comptime SWEEP_LANES = 16
"""Cells a vector sweep fills at once in 32-bit lanes; twice that in 16-bit ones."""

comptime VECTOR_SWEEP_ROWS = 2 * SWEEP_LANES
"""Rows from which a half is swept by anti-diagonal in lanes: fewer would leave lanes idle."""


def vector_sweep_bands[
    half: SweepHalf, dtype: DType, width: Int
](
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    area: Rectangle,
    entering: GapRun,
    lookup: SubstitutionLookup,
    gaps: AffineGapCosts,
    last_scores: MutSpan[Int32, _],
    last_deletes: MutSpan[Int32, _],
):
    """The last row of `area` toward its cut, its scores and deletion layer, as `fill_rows` gives them,
    swept by anti-diagonal `width` cells at a time in lanes of `dtype` (see `anti_diagonals`): each
    diagonal's cell on the last row taken as the diagonal passes it. `entering` opens a deletion run above
    the first row."""
    comptime Value = Scalar[dtype]
    comptime backward = half == SweepHalf.REVERSE
    var rows = area.height()
    var columns = area.width()
    var open = Int(gaps.open)
    var extend = Int(gaps.extend)
    var sweep = AntiDiagonals[dtype, width](
        first[area.row_from : area.row_to], second[area.column_from : area.column_to], 0, backward, backward
    )

    @always_inline
    def top_border(column: Int) {imm open, imm extend} -> Int:
        return 0 if column == 0 else open + (column - 1) * extend

    @always_inline
    def left_border(row: Int) {imm open, imm extend, imm entering} -> Int:
        return row * extend if entering == GapRun.EXTENDS else open + (row - 1) * extend

    var cells = sweep.cells()
    # The corner, and the two border cells beside it, diagonal one.
    cells.two_back[unsafe_offset=0] = 0
    cells.one_back[unsafe_offset=0] = Value(top_border(1))
    cells.deletes_back[unsafe_offset=0] = Value(top_border(1) + open + extend)
    cells.one_back[unsafe_offset=1] = Value(left_border(1))
    cells.deletes_back[unsafe_offset=1] = Value(left_border(1))
    cells.inserts_back[unsafe_offset=1] = Value(left_border(1) + open + extend)
    if rows == 1:
        last_scores[0] = Int32(cells.one_back[unsafe_offset=1])
        last_deletes[0] = Int32(cells.deletes_back[unsafe_offset=1])

    var substitute = lookup.lanes[width, dtype]()
    var lanes = GapLanes[dtype, width].symmetric(open, extend)
    for diagonal in range(2, rows + columns + 1):
        var row = max(1, diagonal - columns)
        while row <= min(rows, diagonal - 1):
            _ = cells.step(row, columns - diagonal, substitute, lanes)
            row += width
        # The diagonal's border cells, after the lanes that may have run over them.
        if diagonal <= columns:
            cells.current[unsafe_offset=0] = Value(top_border(diagonal))
            cells.deletes[unsafe_offset=0] = Value(top_border(diagonal) + open + extend)
        if diagonal <= rows:
            cells.current[unsafe_offset=diagonal] = Value(left_border(diagonal))
            cells.deletes[unsafe_offset=diagonal] = Value(left_border(diagonal))
            cells.inserts[unsafe_offset=diagonal] = Value(left_border(diagonal) + open + extend)
        if diagonal >= rows:
            last_scores[diagonal - rows] = Int32(cells.current[unsafe_offset=rows])
            last_deletes[diagonal - rows] = Int32(cells.deletes[unsafe_offset=rows])
        cells.advance()


def sweep_half[
    half: SweepHalf
](
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    area: Rectangle,
    entering: GapRun,
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    lookup: SubstitutionLookup,
    gaps: AffineGapCosts,
    in_lanes: Bool,
    fits_16: Bool,
    mut rolling: RollingBands,
    last_scores: MutSpan[Int32, _],
    last_deletes: MutSpan[Int32, _],
):
    """`area`'s last row toward its cut, its scores and deletion layer: `in_lanes`, 16-bit ones when `fits_16`,
    else a row at a time through `rolling`."""
    if in_lanes:
        if fits_16:
            vector_sweep_bands[half, DType.int16, 2 * SWEEP_LANES](
                first, second, area, entering, lookup, gaps, last_scores, last_deletes
            )
        else:
            vector_sweep_bands[half, DType.int32, SWEEP_LANES](
                first, second, area, entering, lookup, gaps, last_scores, last_deletes
            )
        return
    var layout = RollingRows(area.width() + 1)
    _ = fill_rows[AlignmentMode.GLOBAL, half](
        first[area.row_from : area.row_to],
        second[area.column_from : area.column_to],
        substitutions,
        alphabet_size,
        gaps,
        entering,
        layout,
        rolling.scores,
        rolling.deletes,
        rolling.inserts,
    )
    for column in range(area.width() + 1):
        last_scores[column] = rolling.scores[layout.index(area.height(), column)]
        last_deletes[column] = rolling.deletes[layout.index(area.height(), column)]


def linear_path(
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    window: Rectangle,
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet_size: Int,
    gaps: AffineGapCosts,
    leaf_cells: Int,
    mut path: RowPath,
):
    """A global alignment of `window` in linear space, its rows written on `path`: frames split at their
    middle row and joined at the best cut (see `best_cut`), the upper solved first, until each is solved
    outright (see `Frame.solved_outright`). A half's sweep runs in 16-bit lanes when its scores fit them."""
    var view = path.view()
    view.leaves(window.row_to, window.column_to)
    # A split's two last rows, each half's scores and deletion layer; no frame is wider than the window.
    var forward_scores = List[Int32](unsafe_uninit_length=window.width() + 1)
    var forward_deletes = List[Int32](unsafe_uninit_length=window.width() + 1)
    var reverse_scores = List[Int32](unsafe_uninit_length=window.width() + 1)
    var reverse_deletes = List[Int32](unsafe_uninit_length=window.width() + 1)
    var rolling = RollingBands(window.width())
    var lookup = SubstitutionLookup(substitutions, alphabet_size)
    # Whether a half's scores fit 16 bits follows from what a pair earns and what a move costs at most.
    var extremes = table_extremes(substitutions, alphabet_size)
    var reward = max(extremes[0], 0)
    var substitution = -min(extremes[1], 0)
    var dearest = max(substitution, -Int(gaps.open))

    var pending: List[Frame] = [Frame(window, GapRun.OPENS, GapRun.OPENS)]
    while len(pending) > 0:
        var frame = pending.pop()
        var area = frame.area
        if area.height() == 0:
            continue
        if frame.solved_outright(leaf_cells):
            solve_outright(first, second, frame, substitutions, alphabet_size, gaps, view)
            continue
        var middle = area.middle()
        # Both halves by the lower, the taller: its lanes are not idle, and its fit holds for both.
        var in_lanes = area.row_to - middle >= VECTOR_SWEEP_ROWS
        var half_rows = area.row_to - middle
        var deficit = straight_deficit(substitution, -Int(gaps.open), -Int(gaps.extend), half_rows, area.width())
        var fits_16 = fits_16_bits[False](reward, dearest, deficit, half_rows, area.width())
        var upper = Rectangle(area.row_from, middle, area.column_from, area.column_to)
        var lower = Rectangle(middle, area.row_to, area.column_from, area.column_to)
        sweep_half[SweepHalf.FORWARD](
            first,
            second,
            upper,
            frame.top,
            substitutions,
            alphabet_size,
            lookup,
            gaps,
            in_lanes,
            fits_16,
            rolling,
            forward_scores,
            forward_deletes,
        )
        sweep_half[SweepHalf.REVERSE](
            first,
            second,
            lower,
            frame.bottom,
            substitutions,
            alphabet_size,
            lookup,
            gaps,
            in_lanes,
            fits_16,
            rolling,
            reverse_scores,
            reverse_deletes,
        )
        var cut = best_cut(forward_scores, forward_deletes, reverse_scores, reverse_deletes, area.width(), gaps)
        var halves = frame.halves(cut)
        # The lower half pushed first, so the upper is solved first and the rows fill in order.
        pending.append(halves[1])
        pending.append(halves[0])


# endregion Linear space
