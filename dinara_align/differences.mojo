# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
A global alignment at gap costs over the whole matrix, swept by anti-diagonal in difference recurrences
(Suzuki and Kasahara, BMC Bioinformatics 2018), for the pairs too divergent for the wavefront, whose work
grows with the cost times the length where the matrix's is fixed.

Each cell keeps its cost less its neighbours' rather than its cost: `u` its cost less the cell above's, `v`
less the cell to its left's, and for each gap layer what a run through it costs the next cell, less this
cell's cost. Then no value depends on the length: each lies within five times the dearest single move
of zero (see `difference_bits`), so DNA's costs sweep in 8-bit lanes, 64 to a 64-byte register, at any
length, where the costs themselves would need 32 bits past a few kilobases. The pair's cost is the first
column's last cell plus the last row's `v`s.

A cell's cost less the one above and to its left, `z`, is the least of its pair's cost, the deletion
layer's run less that cell, the cell above's own `v` added, and the insertion layer's, the cell to the
left's `u` added; its `u` is then `z` less the cell above's `v`, its `v` `z` less the left one's `u`. A
deletion run into a cell costs the cheaper of the cell above's run continued, its value, and one opened
there, an opening and an extension; on a tie it continues, as the wavefront's backtrace takes it.

The traceback keeps every few anti-diagonals and sweeps each stretch between two again as the walk
reaches it (see `traced`), and each diagonal is swept in place (see `DifferenceSweep.run`), so the memory
stays a few cache-sized buffers however long the pair.

With a traceback each cell keeps a flag, the one the lanes keep (see `lanes.swept`): the source of its
cost by the wavefront's priority, a substitution, then a deletion, an insertion, the second piece before
the first, and only then a match; and a bit for each gap layer that extends into the cell. The whole
matrix holds every path, so the walk over the flags takes the path the wavefront's backtrace takes:
`Ties.RIGHT` forward, `Ties.LEFT` over both sequences reversed, as the lanes trace it.
"""

from std.math import ceildiv, sqrt
from std.memory import unsafe_memcpy
from std.sys import size_of

from .cigar import reverse_bytes, reversed_list
from .gap_affine import ALIGNED, ENTRY_MASK, FIRST_GAP, SECOND_GAP, Penalties, gap_layer, layer_bit
from .modes import Ties

comptime SECOND_DELETION = gap_layer(1, True)
"""The second piece's deletion layer."""
comptime SECOND_INSERTION = gap_layer(1, False)
"""The second piece's insertion layer."""


def difference_bits[pieces: Int](penalties: Penalties) -> Int:
    """The narrowest lanes, 8 or 16 bits, that hold every value of a sweep under `penalties`, or 0 for none.

    With `F` the dearest single move (see `Penalties.window`), a cell's `u` and `v` lie within `F` of zero:
    a path into the cell from the row above crosses it at some column, and the cell above is reached from
    there along its row for no more than the crossing's own path pays plus a gap's first letter. A gap
    layer's value is at least its extension and at most twice `F` and an extension, and `z` within twice
    `F`; so every sum and difference the recurrence forms, junk in a step's spare lanes aside, lies within
    five times `F`."""
    var dearest = penalties.window[pieces]()
    if 5 * dearest < 127:
        return 8
    if 5 * dearest < 32767:
        return 16
    return 0


def path_bytes(rows: Int, columns: Int) -> Int:
    """What `swept_path` keeps for a `rows` by `columns` pair at most, its lanes 16 bits and two gap pieces: the
    anti-diagonals it keeps and one stretch's flags (see `stretch_diagonals`), and the letters and two
    diagonals of every layer."""
    var stretch = stretch_diagonals(rows + columns, MAX_LAYERS * 2)
    var kept = (ceildiv(rows + columns, stretch) + 1) * MAX_LAYERS * 2 * (rows + 1)
    return kept + (stretch + 2) * (rows + 128) + 16 * (rows + columns + 256)


comptime MAX_LAYERS = 6
"""A cell's values with two gap pieces: `u`, `v`, and four gap layers."""


def stretch_diagonals(diagonals: Int, cell_bytes: Int) -> Int:
    """The anti-diagonals between two the traceback keeps, of a sweep over `diagonals` of them keeping
    `cell_bytes` a cell: the square root of their product, which keeps the least memory, the kept diagonals
    and one stretch's flags, a byte a cell, weighing about the same (as `lanes.stretch_rows` keeps a band's
    rows)."""
    return max(Int(sqrt(Float64(diagonals * cell_bytes))), 1)


def gap_cost[pieces: Int](penalties: Penalties, letters: Int, deleted: Bool) -> Int:
    """A run of `letters` letters of the first sequence alone with `deleted`, else of the second, at its cheaper
    piece; nothing for none."""
    if letters == 0:
        return 0
    var opening = penalties.deletion_opening if deleted else penalties.opening
    var extension = penalties.deletion_extension if deleted else penalties.extension
    var cost = opening + letters * extension
    comptime if pieces == 2:
        var opening2 = penalties.deletion_opening2 if deleted else penalties.opening2
        var extension2 = penalties.deletion_extension2 if deleted else penalties.extension2
        cost = min(cost, opening2 + letters * extension2)
    return cost


def swept_cost[pieces: Int](first: Span[UInt8, _], second: Span[UInt8, _], penalties: Penalties) -> Int:
    """The least cost of a global alignment of `first` and `second` under `penalties`, both holding a letter,
    which `difference_bits` must allow."""
    if difference_bits[pieces](penalties) == 8:
        return swept_cost_in[DType.int8, pieces](first, second, penalties)
    return swept_cost_in[DType.int16, pieces](first, second, penalties)


def swept_cost_in[
    value: DType, pieces: Int
](first: Span[UInt8, _], second: Span[UInt8, _], penalties: Penalties) -> Int:
    """`swept_cost` in lanes of `value`."""
    var sweep = DifferenceSweep[value, pieces](first, second, penalties)
    var unused = List[UInt8]()
    sweep.run[False](1, len(first) + len(second), unused)
    return sweep.total


def swept_path[
    pieces: Int
](first: Span[UInt8, _], second: Span[UInt8, _], penalties: Penalties, ties: Ties, mut moves: List[UInt8]) -> Int:
    """The least cost of a global alignment of `first` and `second` under `penalties`, both holding a letter,
    which `difference_bits` must allow, with the moves of the path `ties` picks appended right to left, as
    `gap_affine.solve` appends them."""
    if ties == Ties.RIGHT:
        return traced_by[pieces](first, second, penalties, moves)
    # Over both sequences reversed the walk runs from the origin on, left to right: turned around after.
    var reversed_first = reversed_list(first)
    var reversed_second = reversed_list(second)
    var appended = len(moves)
    var cost = traced_by[pieces](Span(reversed_first), Span(reversed_second), penalties, moves)
    reverse_bytes(moves.unsafe_ptr().unsafe_offset(appended), len(moves) - appended)
    return cost


def traced_by[
    pieces: Int
](first: Span[UInt8, _], second: Span[UInt8, _], penalties: Penalties, mut moves: List[UInt8]) -> Int:
    """`traced` in the lanes `difference_bits` picks."""
    if difference_bits[pieces](penalties) == 8:
        return traced[DType.int8, pieces](first, second, penalties, moves)
    return traced[DType.int16, pieces](first, second, penalties, moves)


def traced[
    value: DType, pieces: Int
](first: Span[UInt8, _], second: Span[UInt8, _], penalties: Penalties, mut moves: List[UInt8]) -> Int:
    """The least cost of a global alignment of `first` and `second` in lanes of `value`, with the moves of the
    path the wavefront's backtrace takes from the corner appended right to left.

    The flags of the whole matrix would take a byte a cell, megabytes of fresh pages for every pair of a
    few kilobases, which cost more to fault in than the sweep takes. So the sweep keeps every
    `stretch_diagonals`-th anti-diagonal's values, and the walk sweeps each stretch between two again from
    the one kept before it, keeping its flags alone, from the last stretch back; the last one's flags are
    kept as the first sweep passes."""
    comptime Sweep = DifferenceSweep[value, pieces]
    var rows = len(first)
    var diagonals = rows + len(second)
    var sweep = Sweep(first, second, penalties)
    var state = Sweep.LAYERS * (rows + 1)
    var stretch = stretch_diagonals(diagonals, Sweep.LAYERS * size_of[value]())
    var stretches = ceildiv(diagonals, stretch)
    var kept = List[Scalar[value]](unsafe_uninit_length=stretches * state)
    var flags = List[UInt8]()
    sweep.save(kept, 0)
    for index in range(stretches):
        var last = min((index + 1) * stretch, diagonals)
        if index == stretches - 1:
            sweep.run[True](index * stretch + 1, last, flags)
        else:
            sweep.run[False](index * stretch + 1, last, flags)
            sweep.save(kept, index + 1)
    # Sweeping a stretch again adds its last-row cells again: the cost is the first sweep's.
    var cost = sweep.total
    var walk = Walk(rows, len(second))
    for index in reversed(range(stretches)):
        if walk.row == 0 or walk.column == 0:
            break
        var start = index * stretch
        if index < stretches - 1:
            sweep.restore(kept, index)
            sweep.run[True](start + 1, min((index + 1) * stretch, diagonals), flags)
        walk.through(flags, start, interior_before(start + 1, rows, len(second)), moves)
    walk.finish(moves)
    return cost


struct DifferenceSweep[value: DType, pieces: Int](Movable):
    """A pair's sweep in differences (see the module's notes): its letters, laid out for a diagonal's lanes,
    the costs as lanes, and one anti-diagonal of every layer, a cell at its row, which each diagonal's sweep
    overwrites as it goes (see `run`); and the first column's last cell plus the last row's `v`s swept so
    far, which once every diagonal is swept is the pair's cost."""

    comptime WIDTH = 64 // size_of[Self.value]()
    comptime LAYERS = 4 if Self.pieces == 1 else 6
    comptime Value = Scalar[Self.value]
    comptime Lanes = SIMD[Self.value, Self.WIDTH]
    comptime UP = 0
    """A cell's cost less the one above's."""
    comptime LEFT = 1
    """A cell's cost less the one to its left's."""
    comptime DELETES = 2
    """What a deletion run through the cell costs the one below, less this cell's cost; then the insertions,
    then the second piece's two."""

    var rows: Int
    var columns: Int
    var span: Int
    var down: List[UInt8]
    var across: List[UInt8]
    var cells: List[Self.Value]
    """Layer `l`'s cells of the last diagonal swept, at `l span`."""
    var penalties: Penalties
    var total: Int

    def __init__(out self, first: Span[UInt8, _], second: Span[UInt8, _], penalties: Penalties):
        """The sweep of `first` down the rows and `second` across, before its first diagonal."""
        self.rows = len(first)
        self.columns = len(second)
        self.span = self.rows + 2 + 2 * Self.WIDTH
        # Row `i`'s letter at `i`; column `j`'s at `columns - j`, so a diagonal's cells load both contiguously.
        self.down = List[UInt8](length=self.span, fill=0)
        self.across = List[UInt8](length=self.columns + self.span, fill=0)
        for row in range(self.rows):
            self.down[row + 1] = first[row]
        for column in range(self.columns):
            self.across[self.columns - 1 - column] = second[column]
        self.cells = List[Self.Value](length=Self.LAYERS * self.span, fill=0)
        self.penalties = penalties
        self.total = gap_cost[Self.pieces](penalties, self.rows, True)

    def save(self, mut kept: List[Self.Value], index: Int):
        """Keeps the last diagonal swept in place `index` of `kept`."""
        var state = Self.LAYERS * (self.rows + 1)
        for layer in range(Self.LAYERS):
            unsafe_memcpy(
                dest=kept.unsafe_ptr().unsafe_offset(index * state + layer * (self.rows + 1)),
                src=self.cells.unsafe_ptr().unsafe_offset(layer * self.span),
                count=self.rows + 1,
            )

    def restore(mut self, kept: List[Self.Value], index: Int):
        """Puts back the diagonal kept in place `index` as the last swept."""
        var state = Self.LAYERS * (self.rows + 1)
        for layer in range(Self.LAYERS):
            unsafe_memcpy(
                dest=self.cells.unsafe_ptr().unsafe_offset(layer * self.span),
                src=kept.unsafe_ptr().unsafe_offset(index * state + layer * (self.rows + 1)),
                count=self.rows + 1,
            )

    def run[record: Bool](mut self, start: Int, end: Int, mut flags: List[UInt8]):
        """Sweeps diagonals `start ..= end`, the last swept `start - 1`; with `record`, each interior cell's flag
        in `flags`, anti-diagonal after anti-diagonal, each by row, from diagonal `start`'s first.

        A diagonal overwrites the last in place, its steps from the last row up: a cell reads the last diagonal
        on its own row and the row above, which no step below it has reached yet and which the step above it
        read before this one writes. One diagonal of each layer, not two, stays in the core's first cache to
        twice the length: two took 0.36 ns a cell on 5 kbp pairs and 0.45 on 10 kbp, where 1 to 2 kbp took
        0.21 (Skylake-X)."""
        comptime WIDTH = Self.WIDTH
        comptime Lanes = Self.Lanes
        comptime Value = Self.Value
        comptime Flags = SIMD[DType.uint8, WIDTH]
        comptime two = Self.pieces == 2
        var rows = self.rows
        var columns = self.columns
        var span = self.span
        comptime if record:
            var count = interior_before(end + 1, rows, columns) - interior_before(start, rows, columns)
            flags.resize(unsafe_uninit_length=count + 2 * WIDTH)
        var penalties = self.penalties
        var delete_open = penalties.deletion_opening + penalties.deletion_extension
        var insert_open = penalties.opening + penalties.extension
        var delete_open2 = penalties.deletion_opening2 + penalties.deletion_extension2
        var insert_open2 = penalties.opening2 + penalties.extension2
        var mismatch = Lanes(Value(penalties.mismatch))
        var opened_deletion = Lanes(Value(delete_open))
        var opened_insertion = Lanes(Value(insert_open))
        var deletion_extension = Lanes(Value(penalties.deletion_extension))
        var insertion_extension = Lanes(Value(penalties.extension))
        var opened_deletion2 = Lanes(Value(delete_open2 if two else 0))
        var opened_insertion2 = Lanes(Value(insert_open2 if two else 0))
        var deletion_extension2 = Lanes(Value(penalties.deletion_extension2 if two else 0))
        var insertion_extension2 = Lanes(Value(penalties.extension2 if two else 0))
        var nothing = Lanes(0)

        var down_letters = self.down.unsafe_ptr()
        var across_letters = self.across.unsafe_ptr()
        var cells = self.cells.unsafe_ptr()
        var up_cells = cells.unsafe_offset(Self.UP * span)
        var left_cells = cells.unsafe_offset(Self.LEFT * span)
        var delete_cells = cells.unsafe_offset(Self.DELETES * span)
        var insert_cells = cells.unsafe_offset((Self.DELETES + 1) * span)
        var delete2_cells = cells.unsafe_offset((Self.DELETES + 2 if two else Self.DELETES) * span)
        var insert2_cells = cells.unsafe_offset((Self.DELETES + 3 if two else Self.DELETES) * span)
        var flag_cells = flags.unsafe_ptr()
        var written = 0
        for diagonal in range(start, end + 1):
            var low = max(1, diagonal - columns)
            var high = min(rows, diagonal - 1)
            var lag = columns - diagonal
            var steps = ceildiv(high - low + 1, WIDTH) if high >= low else 0
            for step in reversed(range(steps)):
                var row = low + step * WIDTH
                # The cell above is the last diagonal's on the row before, the one to the left its on this row.
                var above = row - 1
                var left = row
                var mine = down_letters.unsafe_offset(row).unsafe_load[width=WIDTH]()
                var theirs = across_letters.unsafe_offset(lag + row).unsafe_load[width=WIDTH]()
                var equal = mine.eq(theirs)
                var substituted = equal.select(nothing, mismatch)
                # The cell above's cost and the left one's, less the cell above and to the left's.
                var above_over = left_cells.unsafe_offset(above).unsafe_load[width=WIDTH]()
                var left_over = up_cells.unsafe_offset(left).unsafe_load[width=WIDTH]()
                var extended_deletion = delete_cells.unsafe_offset(above).unsafe_load[width=WIDTH]()
                var extended_insertion = insert_cells.unsafe_offset(left).unsafe_load[width=WIDTH]()
                var deletion = min(opened_deletion, extended_deletion) + above_over
                var insertion = min(opened_insertion, extended_insertion) + left_over
                var cost = min(substituted, min(deletion, insertion))
                var deletion2 = nothing
                var insertion2 = nothing
                var extended_deletion2 = nothing
                var extended_insertion2 = nothing
                comptime if two:
                    extended_deletion2 = delete2_cells.unsafe_offset(above).unsafe_load[width=WIDTH]()
                    extended_insertion2 = insert2_cells.unsafe_offset(left).unsafe_load[width=WIDTH]()
                    deletion2 = min(opened_deletion2, extended_deletion2) + above_over
                    insertion2 = min(opened_insertion2, extended_insertion2) + left_over
                    cost = min(cost, min(deletion2, insertion2))
                    delete2_cells.unsafe_offset(row).unsafe_store(deletion2 - cost + deletion_extension2)
                    insert2_cells.unsafe_offset(row).unsafe_store(insertion2 - cost + insertion_extension2)
                up_cells.unsafe_offset(row).unsafe_store(cost - above_over)
                left_cells.unsafe_offset(row).unsafe_store(cost - left_over)
                delete_cells.unsafe_offset(row).unsafe_store(deletion - cost + deletion_extension)
                insert_cells.unsafe_offset(row).unsafe_store(insertion - cost + insertion_extension)
                comptime if record:
                    # The lowest priority first, each source as cheap as the cell taking over from the last.
                    var flag = insertion.eq(cost).select(Flags(UInt8(SECOND_GAP)), Flags(UInt8(ALIGNED)))
                    comptime if two:
                        flag = insertion2.eq(cost).select(Flags(UInt8(SECOND_INSERTION)), flag)
                    flag = deletion.eq(cost).select(Flags(UInt8(FIRST_GAP)), flag)
                    comptime if two:
                        flag = deletion2.eq(cost).select(Flags(UInt8(SECOND_DELETION)), flag)
                    flag = (~equal & substituted.eq(cost)).select(Flags(UInt8(ALIGNED)), flag)
                    flag |= extended_deletion.le(opened_deletion).select(Flags(layer_bit(FIRST_GAP)), Flags(0))
                    flag |= extended_insertion.le(opened_insertion).select(Flags(layer_bit(SECOND_GAP)), Flags(0))
                    comptime if two:
                        flag |= extended_deletion2.le(opened_deletion2).select(
                            Flags(layer_bit(SECOND_DELETION)), Flags(0)
                        )
                        flag |= extended_insertion2.le(opened_insertion2).select(
                            Flags(layer_bit(SECOND_INSERTION)), Flags(0)
                        )
                    flag_cells.unsafe_offset(written + row - low).unsafe_store(flag)
            written += max(high - low + 1, 0)
            # The borders, written after the step's spare lanes may have run over them. Row zero's cell is at
            # index zero: an insertion of every column so far, no deletion run reaching it. Column zero's is at
            # its row: a deletion of every row so far, no insertion run reaching it.
            if diagonal <= columns:
                left_cells[unsafe_offset=0] = Value(
                    gap_cost[Self.pieces](penalties, diagonal, False)
                    - gap_cost[Self.pieces](penalties, diagonal - 1, False)
                )
                delete_cells[unsafe_offset=0] = Value(delete_open + 1)
                comptime if two:
                    delete2_cells[unsafe_offset=0] = Value(delete_open2 + 1)
            if diagonal <= rows:
                up_cells[unsafe_offset=diagonal] = Value(
                    gap_cost[Self.pieces](penalties, diagonal, True)
                    - gap_cost[Self.pieces](penalties, diagonal - 1, True)
                )
                insert_cells[unsafe_offset=diagonal] = Value(insert_open + 1)
                comptime if two:
                    insert2_cells[unsafe_offset=diagonal] = Value(insert_open2 + 1)
            # The last row's cell on this diagonal, `(rows, diagonal - rows)`.
            if diagonal > rows:
                self.total += Int(left_cells[unsafe_offset=rows])


struct Walk(Movable):
    """The walk back over a sweep's flags (see `traced`), as `lanes.walked` takes one: the cell it stands at
    and the layer it is in."""

    var row: Int
    var column: Int
    var columns: Int
    var rows: Int
    var layer: Int

    def __init__(out self, rows: Int, columns: Int):
        """A walk from the corner of a `rows` by `columns` matrix, in the alignment layer."""
        self.row = rows
        self.column = columns
        self.rows = rows
        self.columns = columns
        self.layer = ALIGNED

    def through(mut self, flags: List[UInt8], before: Int, first_flag: Int, mut moves: List[UInt8]):
        """Appends the moves the flags of a stretch name, the anti-diagonals after `before` whose flags start at
        interior cell `first_flag`, until the walk leaves the stretch or reaches the first row or column: at
        each cell the source its flag names, a gap layer's run until it opens."""
        while self.row > 0 and self.column > 0 and self.row + self.column > before:
            var diagonal = self.row + self.column
            var at = interior_before(diagonal, self.rows, self.columns) - first_flag
            var flag = flags[at + self.row - max(1, diagonal - self.columns)]
            if self.layer == ALIGNED:
                self.layer = Int(flag & ENTRY_MASK)
                if self.layer == ALIGNED:
                    moves.append(UInt8(ALIGNED))
                    self.row -= 1
                    self.column -= 1
                continue
            var extended = flag & layer_bit(self.layer) != 0
            if self.layer % 2 == 1:
                moves.append(UInt8(FIRST_GAP))
                self.row -= 1
            else:
                moves.append(UInt8(SECOND_GAP))
                self.column -= 1
            if not extended:
                self.layer = ALIGNED

    def finish(mut self, mut moves: List[UInt8]):
        """Appends the one way left from the first row or column, a gap along it."""
        for _ in range(self.row):
            moves.append(UInt8(FIRST_GAP))
        for _ in range(self.column):
            moves.append(UInt8(SECOND_GAP))
        self.row = 0
        self.column = 0


@always_inline
def interior_before(diagonal: Int, rows: Int, columns: Int) -> Int:
    """The interior cells, those off the first row and column, on the anti-diagonals before `diagonal`."""
    # Diagonal `d` holds rows `max(1, d - columns) ..= min(rows, d - 1)`: summed in closed form, the
    # diagonals up to `d - 1` minus those past the columns' end and past the rows' end.
    var last = diagonal - 1
    var count = last * (last - 1) // 2
    var past_columns = max(last - columns - 1, 0)
    var past_rows = max(last - rows - 1, 0)
    return count - past_columns * (past_columns + 1) // 2 - past_rows * (past_rows + 1) // 2
