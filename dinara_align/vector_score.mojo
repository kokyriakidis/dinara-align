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

from std.math import ceildiv, sqrt
from std.memory import unsafe_memcpy
from std.sys import size_of

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
    strip_rows,
    unreachable,
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

    var substitute = lookup.lanes[width, dtype]()
    var steps = GapLanes[dtype, width].symmetric(Int(gaps.open), Int(gaps.extend))
    var wanted = Lanes(Value(target))
    var lane_rows = Lanes()
    comptime for lane in range(width):
        lane_rows[lane] = Value(lane)
    # The rows go in strips that stay in the core's first cache (see `anti_diagonals.strip_rows`), each swept
    # from the row above it, which the strip before leaves in `top`, scores and deletion runs alike. The start
    # is the first cell that reaches the target by anti-diagonal, then by row: a strip's own first, and the
    # strips below it, whose rows come later, only on an earlier anti-diagonal, so each stops once it reaches
    # the best one found.
    var strip = strip_rows[dtype, width, 1]()
    var top = List[Value](length=columns + 1 + width, fill=0)
    var top_deletes = List[Value](length=columns + 1 + width, fill=0)
    var top_cells = top.unsafe_ptr()
    var top_delete_cells = top_deletes.unsafe_ptr()
    for column in range(1, columns + 1):
        top_cells[unsafe_offset=column] = Value(border(column))
        top_delete_cells[unsafe_offset=column] = Value(border(column) + opened)
    var best_diagonal = Int.MAX
    var best_row = 0
    for strip_start in range(0, rows, strip):
        # A strip's first interior cell lies on anti-diagonal `strip_start + 2`.
        if strip_start + 2 >= best_diagonal:
            break
        var height = min(strip, rows - strip_start)
        var strip_end = strip_start + height
        var below = strip_end < rows
        # Row `i` reads the `i`-th letter back from the end, and so does column `j`.
        var sweep = AntiDiagonals[dtype, width](
            Span(first)[end_row - strip_end : end_row - strip_start], Span(second)[:end_column], 0, True, True
        )
        var cells = sweep.cells()
        cells.two_back[unsafe_offset=0] = Value(border(strip_start))
        cells.one_back[unsafe_offset=0] = top_cells[unsafe_offset=1]
        cells.deletes_back[unsafe_offset=0] = top_delete_cells[unsafe_offset=1]
        cells.one_back[unsafe_offset=1] = Value(border(strip_start + 1))
        cells.inserts_back[unsafe_offset=1] = Value(border(strip_start + 1) + opened)
        var found = False
        for local_diagonal in range(2, height + columns + 1):
            var diagonal = strip_start + local_diagonal
            if diagonal >= best_diagonal:
                break
            var high = min(height, local_diagonal - 1)
            var lag = columns - local_diagonal
            var row = max(1, local_diagonal - columns)
            while row <= high:
                var score = cells.step(row, lag, substitute, steps)
                var reached = score.ge(wanted) & lane_rows.lt(Lanes(Value(high - row + 1)))
                if reached.reduce_or():
                    # The lanes run by row, so the first that reached it is the diagonal's earliest.
                    for lane in range(width):
                        if reached[lane]:
                            best_row = strip_start + row + lane
                            best_diagonal = diagonal
                            break
                    found = True
                    break
                row += width
            if found:
                break
            # The strip's last row, on this diagonal past its first column, for the strip below, if any.
            var bottom_column = local_diagonal - height
            if below and bottom_column >= 1 and bottom_column <= columns:
                top_cells[unsafe_offset=bottom_column] = cells.current[unsafe_offset=height]
                top_delete_cells[unsafe_offset=bottom_column] = cells.deletes[unsafe_offset=height]
            if local_diagonal <= columns:
                cells.current[unsafe_offset=0] = top_cells[unsafe_offset=local_diagonal]
                cells.deletes[unsafe_offset=0] = top_delete_cells[unsafe_offset=local_diagonal]
            if local_diagonal <= height:
                cells.current[unsafe_offset=local_diagonal] = Value(border(diagonal))
                cells.inserts[unsafe_offset=local_diagonal] = Value(border(diagonal) + opened)
            cells.advance()
    if best_diagonal < Int.MAX:
        return (end_row - best_row, end_column - (best_diagonal - best_row))
    # The whole prefixes, which a caller with a score above zero never needs.
    return (0, 0)


def vector_bytes(rows: Int, columns: Int, band_cells: Int) -> Int:
    """What `vector_align` holds for a `rows` by `columns` pair whose band has `band_cells` cells, in bytes, at
    32 bits: every cell's decision in one go, else, past the memory allowed, a few kept anti-diagonals and one
    stretch's decisions (see `traced_band`), whichever it would take."""
    return min(whole_band_bytes(rows, columns, band_cells), stretched_band_bytes(rows, columns, band_cells))


def whole_band_bytes(rows: Int, columns: Int, band_cells: Int) -> Int:
    """`vector_align`'s bytes keeping every cell's decision: those, a step's spare lanes, the three indexes of
    every diagonal, and one diagonal of each layer and two of scores."""
    var diagonals = rows + columns + 1
    return band_cells + 4 * WIDTH + diagonals * 3 * 8 + 4 * (rows + 2 + 4 * WIDTH) * 4


def stretched_band_bytes(rows: Int, columns: Int, band_cells: Int) -> Int:
    """`vector_align`'s bytes keeping a few anti-diagonals and one stretch's decisions (see `band_stretch`)."""
    var diagonals = rows + columns + 1
    var stretch = band_stretch(diagonals, rows, band_cells, 4)
    var kept = (ceildiv(diagonals, stretch) + 1) * 4 * (rows + 1) * 4
    return kept + stretch * (rows + 1 + 4 * WIDTH) + diagonals * 3 * 8 + 4 * (rows + 2 + 4 * WIDTH) * 4


def band_stretch(diagonals: Int, rows: Int, band_cells: Int, value_bytes: Int) -> Int:
    """The anti-diagonals between two a band's traceback keeps: the square root of the diagonals times a kept
    diagonal's bytes, four layers of `rows + 1` cells of `value_bytes`, over a diagonal's decisions, which
    keeps the least memory, the kept diagonals and one stretch's decisions weighing about the same."""
    var decisions = max(band_cells // max(diagonals, 1), 1)
    return max(Int(sqrt(Float64(diagonals) * Float64(4 * (rows + 1) * value_bytes) / Float64(decisions))), 1)


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
    budget: Int = Int.MAX,
) -> GappedAlignment:
    """`serial_align`'s global alignment, swept by anti-diagonal a vector of cells a step, in lanes of 16 bits
    while the scores fit (see `anti_diagonals.fits_16_bits`), else 32.

    Only the cells on diagonals `low_diagonal ..= high_diagonal`, column minus row, are computed, the rest
    reading as unreachable; given a band holding every optimal path (see `optimal_band`), those cells hold
    the same values as `serial_align`'s. Each keeps only its decision, the byte `gotoh.decide` would give it
    from those values, stored anti-diagonal by anti-diagonal, each diagonal's band contiguous by row. The
    walk back reads the decisions as `gotoh.walk` reads the scores, so it takes the same path, in a byte a
    cell where the scores took twelve. A band whose decisions would pass `budget` bytes keeps a few
    anti-diagonals instead and sweeps each stretch between two again for its decisions (see `traced_band`).
    """
    var rows = len(first)
    var columns = len(second)
    if rows == 0 or columns == 0:
        return serial_align[AlignmentMode.GLOBAL](first, second, substitutions, alphabet_size, gaps, alphabet)
    var extremes = table_extremes(substitutions, alphabet_size)
    var substitution = -min(extremes[1], 0)
    var deficit = straight_deficit(substitution, -Int(gaps.open), -Int(gaps.extend), rows, columns)
    if fits_16_bits[False](max(extremes[0], 0), max(substitution, -Int(gaps.open)), deficit, rows, columns):
        return traced_band[DType.int16, 2 * WIDTH](
            first, second, lookup, gaps, alphabet, low_diagonal, high_diagonal, budget
        )
    return traced_band[DType.int32, WIDTH](first, second, lookup, gaps, alphabet, low_diagonal, high_diagonal, budget)


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
    budget: Int,
) -> GappedAlignment:
    """`vector_align` in `width` lanes of `dtype`, which hold every score of the band.

    A band whose decisions fit `budget` keeps them all from one sweep. Past it, the decisions of the whole
    band would fault in fresh pages a byte a cell, so the sweep keeps every `band_stretch`-th anti-diagonal's
    cells instead, and the walk sweeps each stretch between two again from the one kept before it, keeping
    its decisions alone, from the last stretch back."""
    comptime Sweep = BandSweep[dtype, width]
    var rows = len(first)
    var columns = len(second)
    var sweep = Sweep(first, second, lookup, gaps, low_diagonal, high_diagonal)
    var diagonals = rows + columns + 1
    var last = diagonals - 1
    var band_cells = sweep.starts[diagonals]
    var decisions = List[UInt8]()
    var walk = BandWalk(rows, columns)
    if whole_band_bytes(rows, columns, band_cells) <= budget:
        sweep.run[True](1, last, decisions)
        walk.through(sweep, first, second, alphabet, decisions, 0)
    else:
        var stretch = band_stretch(diagonals, rows, band_cells, size_of[dtype]())
        var stretches = ceildiv(last, stretch)
        var kept = List[Scalar[dtype]](unsafe_uninit_length=stretches * Sweep.KEPT * (rows + 1))
        sweep.save(kept, 0)
        for index in range(stretches):
            var end = min((index + 1) * stretch, last)
            if index == stretches - 1:
                sweep.run[True](index * stretch + 1, end, decisions)
            else:
                sweep.run[False](index * stretch + 1, end, decisions)
                sweep.save(kept, index + 1)
        var score = sweep.final_score
        for index in reversed(range(stretches)):
            if walk.row == 0 or walk.column == 0:
                break
            var start = index * stretch
            if index < stretches - 1:
                sweep.restore(kept, index)
                sweep.run[True](start + 1, min((index + 1) * stretch, last), decisions)
            walk.through(sweep, first, second, alphabet, decisions, start)
        sweep.final_score = score
    walk.finish(first, second, alphabet)
    return GappedAlignment(Int32(sweep.final_score), walk.top_text(), walk.bottom_text())


struct BandSweep[dtype: DType, width: Int](Movable):
    """A band's sweep (see `vector_align`): its rows on each anti-diagonal, its letters laid out for a diagonal's
    lanes, and the last two anti-diagonals of scores and the last of each gap layer, a cell at its row, which
    each diagonal's sweep overwrites as it goes (see `run`)."""

    comptime KEPT = 4
    """The layers a kept anti-diagonal holds: the last two of scores and the last of each gap layer."""
    comptime Value = Scalar[Self.dtype]

    var rows: Int
    var columns: Int
    var lows: List[Int]
    """Diagonal `d`'s cells lie on rows `lows[d] ..= highs[d]`: inside the matrix, and with `column - row`
    inside the band, column being `d - row`."""
    var highs: List[Int]
    var starts: List[Int]
    """Where diagonal `d`'s decisions start among the band's, row by row."""
    var letters: List[UInt8]
    var reversed: List[UInt8]
    var lookup: SubstitutionLookup
    var gaps: AffineGapCosts
    var scores: List[Self.Value]
    """Two anti-diagonals of scores, the last swept in slot `newest`, at `slot span`."""
    var deletes: List[Self.Value]
    var inserts: List[Self.Value]
    var newest: Int
    var span: Int
    var final_score: Self.Value

    def __init__(
        out self,
        first: List[UInt8],
        second: List[UInt8],
        lookup: SubstitutionLookup,
        gaps: AffineGapCosts,
        low_diagonal: Int,
        high_diagonal: Int,
    ):
        """The sweep of a band of diagonals `low_diagonal ..= high_diagonal`, before its first anti-diagonal."""
        comptime LOW = unreachable[Self.dtype]()
        self.rows = len(first)
        self.columns = len(second)
        var low_band = max(low_diagonal, -self.rows)
        var high_band = min(high_diagonal, self.columns)
        var diagonals = self.rows + self.columns + 1
        self.lows = List[Int](length=diagonals, fill=0)
        self.highs = List[Int](length=diagonals, fill=0)
        self.starts = List[Int](length=diagonals + 1, fill=0)
        for diagonal in range(diagonals):
            var low = max(0, diagonal - self.columns, -((high_band - diagonal) // 2))
            var high = min(self.rows, diagonal, (diagonal - low_band) // 2)
            self.lows[diagonal] = low
            self.highs[diagonal] = high
            self.starts[diagonal + 1] = self.starts[diagonal] + max(high - low + 1, 0)
        self.letters = row_letters[Self.width](Span(first))
        self.reversed = column_letters[Self.width](Span(second))
        self.lookup = lookup.copy()
        self.gaps = gaps
        self.span = self.rows + 2 + 2 * Self.width
        self.scores = List[Self.Value](length=2 * self.span, fill=LOW)
        self.deletes = List[Self.Value](length=self.span, fill=LOW)
        self.inserts = List[Self.Value](length=self.span, fill=LOW)
        # Diagonal zero: the origin, which every band holds.
        self.newest = 0
        self.scores[0] = 0
        self.deletes[0] = 0
        self.inserts[0] = 0
        self.final_score = 0

    def save(self, mut kept: List[Self.Value], index: Int):
        """Keeps the last anti-diagonal swept, and the one before it, in place `index` of `kept`."""
        var length = self.rows + 1
        var at = kept.unsafe_ptr().unsafe_offset(index * Self.KEPT * length)
        unsafe_memcpy(dest=at, src=self.scores.unsafe_ptr().unsafe_offset(self.newest * self.span), count=length)
        unsafe_memcpy(
            dest=at.unsafe_offset(length),
            src=self.scores.unsafe_ptr().unsafe_offset((1 - self.newest) * self.span),
            count=length,
        )
        unsafe_memcpy(dest=at.unsafe_offset(2 * length), src=self.deletes.unsafe_ptr(), count=length)
        unsafe_memcpy(dest=at.unsafe_offset(3 * length), src=self.inserts.unsafe_ptr(), count=length)

    def restore(mut self, kept: List[Self.Value], index: Int):
        """Puts back the anti-diagonals kept in place `index` as the last two swept."""
        var length = self.rows + 1
        var at = kept.unsafe_ptr().unsafe_offset(index * Self.KEPT * length)
        self.newest = 0
        unsafe_memcpy(dest=self.scores.unsafe_ptr(), src=at, count=length)
        unsafe_memcpy(
            dest=self.scores.unsafe_ptr().unsafe_offset(self.span), src=at.unsafe_offset(length), count=length
        )
        unsafe_memcpy(dest=self.deletes.unsafe_ptr(), src=at.unsafe_offset(2 * length), count=length)
        unsafe_memcpy(dest=self.inserts.unsafe_ptr(), src=at.unsafe_offset(3 * length), count=length)

    def run[record: Bool](mut self, start: Int, end: Int, mut decisions: List[UInt8]):
        """Sweeps anti-diagonals `start ..= end`, the last swept `start - 1`; with `record`, each interior cell's
        decision in `decisions`, from diagonal `start`'s first.

        A diagonal overwrites the one two back in place, its steps from the last row up: a cell reads the one
        two back on the row above, which no step below it has reached yet, and the gap layers on its own row
        and the row above, which the step above it read before this one writes. Two diagonals of scores and
        one of each gap layer, not three of each."""
        comptime Lanes = SIMD[Self.dtype, Self.width]
        comptime Value = Self.Value
        comptime LOW = unreachable[Self.dtype]()
        var rows = self.rows
        var columns = self.columns
        var span = self.span
        var open = Int(self.gaps.open)
        var extend = Int(self.gaps.extend)
        var base = self.starts[start]
        comptime if record:
            decisions.resize(unsafe_uninit_length=self.starts[end + 1] - base + 2 * Self.width)
        var decision_cells = decisions.unsafe_ptr()
        var letters = self.letters.unsafe_ptr()
        var across = self.reversed.unsafe_ptr()
        var delete_cells = self.deletes.unsafe_ptr()
        var insert_cells = self.inserts.unsafe_ptr()
        var substitute = self.lookup.lanes[Self.width, Self.dtype]()
        var steps = GapLanes[Self.dtype, Self.width].symmetric(open, extend)
        var source_deleting = Lanes(Value(Int(Layer.DELETING.kind)))
        var source_inserting = Lanes(Value(Int(Layer.INSERTING.kind)))
        var deletion_extends = Lanes(Value(Int(GapRun.EXTENDS.kind) << 2))
        var insertion_extends = Lanes(Value(Int(GapRun.EXTENDS.kind) << 3))
        var aligning = Lanes(Value(Int(Layer.ALIGNING.kind)))
        var nothing = Lanes(0)
        var lows = self.lows.unsafe_ptr()
        var highs = self.highs.unsafe_ptr()
        var starts = self.starts.unsafe_ptr()
        for diagonal in range(start, end + 1):
            var first_row = lows[unsafe_offset=diagonal]
            var last_row = highs[unsafe_offset=diagonal]
            # The last diagonal's scores, and the one two back's, which this one overwrites.
            var last_cells = self.scores.unsafe_ptr().unsafe_offset(self.newest * span)
            var score_cells = self.scores.unsafe_ptr().unsafe_offset((1 - self.newest) * span)
            var low = max(1, first_row)
            var high = min(last_row, diagonal - 1)
            var lag = columns - diagonal
            var decided = starts[unsafe_offset=diagonal] - base - first_row
            var row = low + (ceildiv(high - low + 1, Self.width) - 1) * Self.width if high >= low else low - Self.width
            while row >= low:
                var above = last_cells.unsafe_offset(row - 1).unsafe_load[width=Self.width]()
                var above_delete = delete_cells.unsafe_offset(row - 1).unsafe_load[width=Self.width]()
                var left = last_cells.unsafe_offset(row).unsafe_load[width=Self.width]()
                var left_insert = insert_cells.unsafe_offset(row).unsafe_load[width=Self.width]()
                var above_left = score_cells.unsafe_offset(row - 1).unsafe_load[width=Self.width]()
                var mine = letters.unsafe_offset(row).unsafe_load[width=Self.width]()
                var theirs = across.unsafe_offset(lag + row).unsafe_load[width=Self.width]()
                var aligned = above_left + substitute(mine, theirs)
                var deletion = gap_layer(above, above_delete, steps.down_first, steps.down_further)
                var insertion = gap_layer(left, left_insert, steps.across_first, steps.across_further)
                var score = best_move[False](aligned, deletion, insertion)
                score_cells.unsafe_offset(row).unsafe_store(score)
                delete_cells.unsafe_offset(row).unsafe_store(deletion)
                insert_cells.unsafe_offset(row).unsafe_store(insertion)
                comptime if record:
                    # `gotoh.decide`'s byte: the source, the aligned move before a deletion before an insertion,
                    # and each run extending only where extending scores strictly more than opening.
                    var source = aligned.eq(score).select(
                        aligning, deletion.eq(score).select(source_deleting, source_inserting)
                    )
                    var deleting = (above_delete + steps.down_further).gt(above + steps.down_first)
                    var inserting = (left_insert + steps.across_further).gt(left + steps.across_first)
                    var code = (
                        source
                        | deleting.select(deletion_extends, nothing)
                        | inserting.select(insertion_extends, nothing)
                    )
                    decision_cells.unsafe_offset(decided + row).unsafe_store(code.cast[DType.uint8]())
                row -= Self.width
            # The border cells inside the band, then the cells either side of it, which the lanes may have run
            # over and the next two diagonals read as unreachable.
            if first_row == 0:
                var border = Value(Int(self.gaps.run(diagonal)))
                score_cells[unsafe_offset=0] = border
                delete_cells[unsafe_offset=0] = border + Value(open + extend)
                insert_cells[unsafe_offset=0] = 0
            else:
                score_cells[unsafe_offset=first_row - 1] = LOW
                delete_cells[unsafe_offset=first_row - 1] = LOW
                insert_cells[unsafe_offset=first_row - 1] = LOW
            if last_row == diagonal:
                var border = Value(Int(self.gaps.run(diagonal)))
                score_cells[unsafe_offset=diagonal] = border
                delete_cells[unsafe_offset=diagonal] = 0
                insert_cells[unsafe_offset=diagonal] = border + Value(open + extend)
            score_cells[unsafe_offset=last_row + 1] = LOW
            delete_cells[unsafe_offset=last_row + 1] = LOW
            insert_cells[unsafe_offset=last_row + 1] = LOW
            self.newest = 1 - self.newest
            if diagonal == rows + columns:
                self.final_score = score_cells[unsafe_offset=rows]


struct BandWalk(Movable):
    """The walk back over a band's decisions (see `traced_band`), `gotoh.walk`'s over the decisions: the cell it
    stands at, the layer it is in, and the gapped rows so far, right to left."""

    var row: Int
    var column: Int
    var layer: Layer
    var top: List[UInt8]
    var bottom: List[UInt8]

    def __init__(out self, rows: Int, columns: Int):
        """A walk from the corner, in the aligning layer."""
        self.row = rows
        self.column = columns
        self.layer = Layer.ALIGNING
        self.top = List[UInt8]()
        self.bottom = List[UInt8]()

    def through[
        dtype: DType, width: Int
    ](
        mut self,
        sweep: BandSweep[dtype, width],
        first: List[UInt8],
        second: List[UInt8],
        alphabet: String,
        decisions: List[UInt8],
        before: Int,
    ):
        """Takes the walk up through the decisions of anti-diagonals `before + 1` on, which start there, until it
        leaves them or reaches the first row or column."""
        var symbols = alphabet.as_bytes()
        var base = sweep.starts[before + 1]
        while self.row > 0 and self.column > 0 and self.row + self.column > before:
            var diagonal = self.row + self.column
            var decision = Decision(decisions[sweep.starts[diagonal] - base + self.row - sweep.lows[diagonal]])
            var step = advance(self.layer, decision)
            self.top.append(symbols[Int(first[self.row - 1])] if step.row_advance != 0 else GAP_BYTE)
            self.bottom.append(symbols[Int(second[self.column - 1])] if step.column_advance != 0 else GAP_BYTE)
            self.row += step.row_advance
            self.column += step.column_advance
            self.layer = step.lands_in

    def finish(mut self, first: List[UInt8], second: List[UInt8], alphabet: String):
        """A global alignment's rest from the first row or column: gaps against the letters left."""
        var symbols = alphabet.as_bytes()
        for letter in range(self.row, 0, -1):
            self.top.append(symbols[Int(first[letter - 1])])
            self.bottom.append(GAP_BYTE)
        for letter in range(self.column, 0, -1):
            self.top.append(GAP_BYTE)
            self.bottom.append(symbols[Int(second[letter - 1])])
        self.row = 0
        self.column = 0
        self.top.reverse()
        self.bottom.reverse()

    def top_text(self) -> String:
        """The first sequence's gapped row, once finished."""
        return String(unsafe_from_utf8=self.top)

    def bottom_text(self) -> String:
        """The second sequence's gapped row, once finished."""
        return String(unsafe_from_utf8=self.bottom)


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
