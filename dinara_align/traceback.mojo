# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Ported from `pa-bitpacking` in A*PA (https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner,
# commit bf2e14e), by Ragnar Groot Koerkamp and Pesho Ivanov, itself translated from Edlib.
"""
An optimal alignment from a band's recorded tile edges, retraced a tile at a time from the corner by
WFA2-lib's rule for ties (see `trace_back`), and written out as two gapped rows or as a CIGAR.
"""

from std.bit import count_trailing_zeros
from std.math import ceildiv, clamp

from .diagonal import best_source, DiagonalFronts, slide_forward
from .slides import GATHERED_SLIDES, LANES, gathered_slides, slide
from .bit_parallel import (
    advance,
    ALL_ONES,
    append_diagonals,
    COLUMN_PADDING,
    DIAGONAL,
    Edge,
    LEFT,
    Profile,
    Trail,
    UP,
    WORD_BITS,
    word_value,
)


comptime WAVEFRONT_FLOOR = 32
"""The fewest costs a tile's forward search may climb before the tile is recomputed instead."""

comptime TILE_WORK = 1 << 16
"""Diagonals a tile's forward search may step before it hands the tile to the recompute."""


struct TileFronts(Movable):
    """A tile's forward search, kept across tiles so a tile allocates nothing: the left-edge rows that
    may start a path, and for each cost from the least of theirs, the furthest column of every
    diagonal a path to the traced cell may use, one row of `columns` a cost."""

    var columns: List[Int32]
    var lows: List[Int]
    """The diagonals each cost's row grew, `lows[level] ..= highs[level]`; the rest read unreached."""
    var highs: List[Int]
    var entry_diagonals: List[Int]
    var entry_costs: List[Int]
    var starts: List[Int]
    var ordered: List[Int]
    """The left-edge rows' diagonals by cost, each cost's from `starts` at its level."""

    def __init__(out self):
        """Empty buffers, with room for a typical tile."""
        self.columns = List[Int32](capacity=4096)
        self.lows = List[Int](capacity=64)
        self.highs = List[Int](capacity=64)
        self.entry_diagonals = List[Int](capacity=256)
        self.entry_costs = List[Int](capacity=256)
        self.starts = List[Int](capacity=64)
        self.ordered = List[Int](capacity=256)


@inline(.never)
def forward_segment(
    profile: Profile,
    edge: Edge,
    first_column: Int,
    end_column: Int,
    end_row: Int,
    score: Int,
    limit: Int,
    mut fronts: TileFronts,
    mut moves: List[UInt8],
) -> Int:
    """Traces from `(end_column, end_row)`, scoring `score`, back to the tile's left edge by WFA2-lib's
    rule, a wavefront grown forward across the tile from its left edge; appends the moves right to
    left and returns the left-edge row, or -1 when the search would cost more than `limit` costs or
    `TILE_WORK` diagonals, and the tile is recomputed instead.

    Each left-edge row starts at its recorded score, which is exact on every optimal path and never
    below the truth elsewhere; so a front, the furthest column reached at no more than its cost,
    never reaches past the truth, and reaches it wherever an optimal path passes. A diagonal is kept
    only while its cost plus how far it lies from the traced cell's leaves room to get there, and
    only rows that could start such a path start at all: walking out from the traced cell's
    diagonal, a row's score plus its distance never falls, so the first row too dear ends the walk.
    The backtrace then takes at each cell the furthest source an optimal path passes, as a
    wavefront's backtrace does, and stops where its matches reach the left edge.

    Every diagonal a path may use lies within the costs to climb of the traced cell's, so the rows
    share one window of diagonals, laid out once per tile. Each row is written over the diagonals it
    grew and two either side, which read unreached; a row steps eight diagonals at a time, and on
    AVX-512 slides them together, as the diagonal transition does (see `step_front`).
    """
    var first = profile.column_codes.unsafe_ptr()
    var second = profile.row_codes.unsafe_ptr()
    var target = end_column - end_row
    comptime UNREACHED = Int32(-(1 << 30))
    comptime PAD = 2
    var low_row = edge.low_row
    var high_row = min(edge.high_row, end_row)
    if low_row > high_row:
        return -1
    # The left-edge rows a path to the traced cell within `score` could start from.
    fronts.entry_diagonals.clear()
    fronts.entry_costs.clear()
    var base = Int.MAX
    var middle = clamp(first_column - target, low_row, high_row)
    var row = middle
    while row >= low_row:
        var cost = edge.score(row)
        if cost + abs(first_column - row - target) > score:
            break
        fronts.entry_diagonals.append(first_column - row)
        fronts.entry_costs.append(cost)
        base = min(base, cost)
        row -= 1
    row = middle + 1
    while row <= high_row:
        var cost = edge.score(row)
        if cost + abs(first_column - row - target) > score:
            break
        fronts.entry_diagonals.append(first_column - row)
        fronts.entry_costs.append(cost)
        base = min(base, cost)
        row += 1
    var count = len(fronts.entry_costs)
    var levels = score - base + 1
    if count == 0 or levels > limit:
        return -1
    # One window of diagonals for every row: no path strays further from the traced cell's than the
    # costs it has left to climb.
    var reach = levels - 1
    var window_low = target - reach
    var stride = 2 * reach + 1 + 2 * PAD
    if levels * stride > TILE_WORK:
        return -1
    # Room past the last row for a lane group's stores past its end.
    fronts.columns.resize(unsafe_uninit_length=levels * stride + LANES)
    var grid = fronts.columns.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
    fronts.lows.resize(levels, 0)
    fronts.highs.resize(levels, 0)
    var lows = fronts.lows.unsafe_ptr()
    var highs = fronts.highs.unsafe_ptr()
    var entry_diagonals = fronts.entry_diagonals.unsafe_ptr()
    var entry_costs = fronts.entry_costs.unsafe_ptr()
    # The rows by cost, `ordered[starts[level] ..< starts[level + 1]]`, so a cost reads only its own.
    fronts.starts.clear()
    fronts.starts.resize(levels + 1, 0)
    fronts.ordered.resize(unsafe_uninit_length=count)
    var starts = fronts.starts.unsafe_ptr()
    var ordered = fronts.ordered.unsafe_ptr()
    for index in range(count):
        starts[unsafe_offset=entry_costs[unsafe_offset=index] - base] += 1
    var running = 0
    for level in range(levels + 1):
        var here = starts[unsafe_offset=level]
        starts[unsafe_offset=level] = running
        running += here
    for index in range(count):
        var level = entry_costs[unsafe_offset=index] - base
        ordered[unsafe_offset=starts[unsafe_offset=level]] = entry_diagonals[unsafe_offset=index]
        starts[unsafe_offset=level] += 1
    # Placing the rows moved each cost's start on to the next cost's, so each is taken back one.
    for level in range(levels, 0, -1):
        starts[unsafe_offset=level] = starts[unsafe_offset=level - 1]
    starts[unsafe_offset=0] = 0

    @inline(.always)
    def row_of(level: Int) {imm grid, imm stride, imm window_low} -> MutPointer[Int32, MutUntrackedOrigin]:
        """A cost's row, indexed by diagonal, two diagonals either side of what it grew reading unreached."""
        return grid.unsafe_offset(level * stride + PAD - window_low)

    comptime Lanes = SIMD[DType.int32, LANES]
    var lane_diagonals = Lanes()
    comptime for lane in range(LANES):
        lane_diagonals[lane] = Int32(lane)
    var column_limit = Lanes(Int32(end_column))
    var row_limit = Lanes(Int32(end_row))
    var unreached = Lanes(UNREACHED)
    var low = Int.MAX
    var high = Int.MIN
    for level in range(levels):
        var cost = base + level
        var new_low = low - 1
        var new_high = high + 1
        for index in range(starts[unsafe_offset=level], starts[unsafe_offset=level + 1]):
            new_low = min(new_low, ordered[unsafe_offset=index])
            new_high = max(new_high, ordered[unsafe_offset=index])
        # Room to reach the traced cell: no more diagonals off than costs left.
        new_low = max(new_low, target - (score - cost))
        new_high = min(new_high, target + (score - cost))
        lows[unsafe_offset=level] = new_low
        highs[unsafe_offset=level] = new_high
        if new_low > new_high:
            low = Int.MAX
            high = Int.MIN
            continue
        var current = row_of(level)
        # The diagonals the cost before reaches: the furthest column at no more than it, then one more
        # edit from it. A substitution past the tile's end leaves the front where it was, as `same`
        # keeps it. The rest of the row reads unreached, the lanes stored past the last included.
        var from_low = max(new_low, low - 1)
        var from_high = min(new_high, high + 1)
        # When the cost before reaches none of them, an empty range past the row, so all of it reads
        # unreached below but for the left-edge rows starting at this cost.
        if from_low > from_high:
            from_low = new_high + PAD + 1
            from_high = new_high + PAD
        var previous = row_of(level - 1)
        var diagonal = from_low
        while diagonal <= from_high:
            var diagonals = lane_diagonals + Int32(diagonal)
            var stop = min(column_limit, row_limit + diagonals)
            var same = previous.unsafe_offset(diagonal).unsafe_load[width=LANES]()
            var below = previous.unsafe_offset(diagonal - 1).unsafe_load[width=LANES]()
            var above = previous.unsafe_offset(diagonal + 1).unsafe_load[width=LANES]()
            var deleted = below.lt(column_limit).select(below + 1, unreached)
            var inserted = above.le(stop).select(above, unreached)
            var entry = max(min(same + 1, stop), max(deleted, inserted))
            comptime if GATHERED_SLIDES:
                var inside = diagonals.le(Int32(from_high))
                var slid = gathered_slides(first, second, inside.select(entry, unreached), diagonals)
                current.unsafe_offset(diagonal).unsafe_store(min(slid, stop))
            else:
                current.unsafe_offset(diagonal).unsafe_store(entry)
            diagonal += LANES
        for diagonal in range(new_low - PAD, from_low):
            current[unsafe_offset=diagonal] = UNREACHED
        for diagonal in range(from_high + 1, new_high + PAD + 1):
            current[unsafe_offset=diagonal] = UNREACHED
        comptime if not GATHERED_SLIDES:
            for diagonal in range(from_low, from_high + 1):
                var column = Int(current[unsafe_offset=diagonal])
                if column >= 0:
                    var stop = min(end_column, end_row + diagonal)
                    current[unsafe_offset=diagonal] = Int32(min(slide(first, second, column, diagonal), stop))
        for index in range(starts[unsafe_offset=level], starts[unsafe_offset=level + 1]):
            var diagonal = ordered[unsafe_offset=index]
            if diagonal < new_low or diagonal > new_high:
                continue
            var stop = min(end_column, end_row + diagonal)
            var slid = min(slide_forward(first, second, first_column, first_column - diagonal), stop)
            current[unsafe_offset=diagonal] = max(current[unsafe_offset=diagonal], Int32(slid))
        low = new_low
        high = new_high

    @inline(.always)
    def at(level: Int, diagonal: Int) {imm row_of, imm lows, imm highs} -> Int:
        """The furthest column on `diagonal` at cost level `level`, far below zero where that cost's row
        did not grow it."""
        if level < 0 or diagonal < lows[unsafe_offset=level] or diagonal > highs[unsafe_offset=level]:
            return Int(UNREACHED)
        return Int(row_of(level)[unsafe_offset=diagonal])

    if at(levels - 1, target) != end_column:
        return -1

    # The backtrace: at each cell's own cost, the furthest source, a substitution before a base of
    # the first sequence alone before one of the second; matches back to the left edge end the tile.
    var block = len(moves)
    var diagonal = target
    var column = end_column
    var level = levels - 1
    while True:
        while level > 0 and at(level - 1, diagonal) >= column:
            level -= 1
        var same = at(level - 1, diagonal)
        var below = at(level - 1, diagonal - 1)
        var above = at(level - 1, diagonal + 1)
        var substituted = same + 1 if same >= first_column and same < end_column else -1
        var deleted = below + 1 if below >= first_column and below < end_column else -1
        var inserted = above if above >= first_column else -1
        var entry = min(column, max(substituted, max(deleted, inserted)))
        if entry <= first_column:
            append_diagonals(moves, column - first_column)
            return first_column - diagonal
        append_diagonals(moves, column - entry)
        if substituted >= entry:
            moves.append(DIAGONAL)
            column = entry - 1
        elif deleted >= entry:
            moves.append(LEFT)
            column = entry - 1
            diagonal -= 1
        else:
            moves.append(UP)
            column = entry
            diagonal += 1
        level -= 1
        if level < 0:
            moves.resize(block, 0)
            return -1


comptime RECOMPUTE_WORDS = 4
"""Words a tile's recompute first takes above the point it traces from, doubling while that falls short."""


struct Recompute(Movable):
    """A tile recompute's buffers, kept across tiles so a recompute neither allocates nor zeroes: every
    slot it reads it has written first."""

    var plus: List[UInt64]
    var minus: List[UInt64]
    var bases: List[Int]
    var low_plus: List[UInt64]
    """The same window swept again from a top row that falls by one a column, a floor under every score."""
    var low_minus: List[UInt64]
    var low_bases: List[Int]

    def __init__(out self):
        """Empty buffers, which the first recompute grows."""
        self.plus = List[UInt64]()
        self.minus = List[UInt64]()
        self.bases = List[Int]()
        self.low_plus = List[UInt64]()
        self.low_minus = List[UInt64]()
        self.low_bases = List[Int]()


def recomputed_segment(
    profile: Profile,
    trail: Trail,
    tile: Int,
    end_column: Int,
    end_row: Int,
    score: Int,
    mut buffers: Recompute,
    mut moves: List[UInt8],
) -> Int:
    """Traces a tile back cell by cell, after sweeping it again and keeping every column's differences.

    The fallback for a tile with more edits than the wavefront search is allowed. As A*PA2 does, it
    first recomputes only a window of rows above the point it traces from, `RECOMPUTE_WORDS` words
    and doubling, rather than the band's whole height (see `window_segment`). Appends the moves
    right to left and returns the left-edge row.
    """
    var top = trail.tops[tile]
    # A tile is recorded only with a word in it, so its end is past its top.
    var end_word = clamp(ceildiv(end_row, WORD_BITS), top + 1, trail.ends[tile])
    var window = RECOMPUTE_WORDS
    while True:
        var first_word = max(top, end_word - window)
        var left: Int
        if profile.extended:
            left = window_segment[True](
                profile, trail, tile, end_column, end_row, score, first_word, end_word, buffers, moves
            )
        else:
            left = window_segment[False](
                profile, trail, tile, end_column, end_row, score, first_word, end_word, buffers, moves
            )
        if left >= 0 or first_word == top:
            return left
        window *= 2


def window_segment[
    extended: Bool
](
    profile: Profile,
    trail: Trail,
    tile: Int,
    end_column: Int,
    end_row: Int,
    score: Int,
    first_word: Int,
    end_word: Int,
    mut buffers: Recompute,
    mut moves: List[UInt8],
) -> Int:
    """`recomputed_segment` over words `[first_word, end_word)` alone, or -1 when they cannot settle it.

    The path is WFA2-lib's, as `Ties.RIGHT` names it: from the traced cell back, an edit into the cell
    whenever one is optimal, a substitution before a letter of the first sequence alone before one of
    the second, and a match only when none is; which is the furthest-reaching source a wavefront's
    backtrace takes, cell by cell. So every candidate's own score is needed, not only the path's.

    The window's left edge comes from the recorded one. Its top reads `+1` a column from above, as a
    band's top does, so every score is a real path's and never below the true one; a second sweep
    from a top falling by one a column, as no score can fall faster along a row, gives a floor. A
    candidate is optimal when its score is the one the path needs, and is not when even its floor
    exceeds that; in between, a taller window must tell. Over the band's whole height every optimal
    cell lies inside the band and reads true, so its scores alone decide.
    """
    var first_column = trail.first_columns[tile]
    var top = first_word
    var count = end_word - first_word
    var width = end_column - first_column
    var whole = first_word == trail.tops[tile]
    var offset = trail.offsets[tile] + first_word - trail.tops[tile]
    var sweeps = 1 if whole else 2
    buffers.plus.resize(unsafe_uninit_length=(width + 1) * count)
    buffers.minus.resize(unsafe_uninit_length=(width + 1) * count)
    buffers.bases.resize(unsafe_uninit_length=(width + 1) * (count + 1))
    buffers.low_plus.resize(unsafe_uninit_length=(width + 1) * count)
    buffers.low_minus.resize(unsafe_uninit_length=(width + 1) * count)
    buffers.low_bases.resize(unsafe_uninit_length=(width + 1) * (count + 1))
    # The left edge's score at the window's top, carried down from the band's.
    var anchor = trail.anchors[tile]
    for word in range(trail.offsets[tile], offset):
        anchor += word_value(trail.edge_plus[word], trail.edge_minus[word])
    var column_low = profile.column_low.unsafe_ptr().unsafe_offset(COLUMN_PADDING + first_column - 1)
    var column_high = profile.column_high.unsafe_ptr().unsafe_offset(COLUMN_PADDING + first_column - 1)
    var row_low = profile.row_low.unsafe_ptr().unsafe_offset(top)
    var row_high = profile.row_high.unsafe_ptr().unsafe_offset(top)
    # The third plane exists only for symbols past `ACGT` (see `Profile.extended`).
    var column_extra = profile.column_extra.unsafe_ptr().unsafe_offset(COLUMN_PADDING + first_column - 1)
    var row_extra = profile.row_extra.unsafe_ptr().unsafe_offset(top)
    for sweep in range(sweeps):
        var floor = sweep == 1
        var plus = (buffers.low_plus if floor else buffers.plus).unsafe_ptr()
        var minus = (buffers.low_minus if floor else buffers.minus).unsafe_ptr()
        var bases = (buffers.low_bases if floor else buffers.bases).unsafe_ptr()
        for word in range(count):
            plus[unsafe_offset=word] = trail.edge_plus[offset + word]
            minus[unsafe_offset=word] = trail.edge_minus[offset + word]
        for step in range(1, width + 1):
            # The top row's score rises by one a column, or for the floor falls by one.
            var horizontal_plus = UInt64(0) if floor else UInt64(1)
            var horizontal_minus = UInt64(1) if floor else UInt64(0)
            var low = column_low[unsafe_offset=step]
            var high = column_high[unsafe_offset=step]
            var extra = column_extra[unsafe_offset=step] if extended else UInt64(0)
            for word in range(count):
                var vertical_plus = plus[unsafe_offset=(step - 1) * count + word]
                var vertical_minus = minus[unsafe_offset=(step - 1) * count + word]
                var matches = (low ^ row_low[unsafe_offset=word]) & (high ^ row_high[unsafe_offset=word])
                comptime if extended:
                    matches &= extra ^ row_extra[unsafe_offset=word]
                advance[1](horizontal_plus, horizontal_minus, vertical_plus, vertical_minus, matches)
                plus[unsafe_offset=step * count + word] = vertical_plus
                minus[unsafe_offset=step * count + word] = vertical_minus
        # Scores at each word's top on every column, from the window's top row.
        for step in range(width + 1):
            var running = anchor - step if floor else anchor + step
            bases[unsafe_offset=step * (count + 1)] = running
            for word in range(count):
                running += word_value(plus[unsafe_offset=step * count + word], minus[unsafe_offset=step * count + word])
                bases[unsafe_offset=step * (count + 1) + word + 1] = running

    @inline(.always)
    def score_in(
        plus: ImmPointer[UInt64, _], minus: ImmPointer[UInt64, _], bases: ImmPointer[Int, _], step: Int, row: Int
    ) {imm count, imm top} -> Int:
        """The score at `row` on the window's column `step`, from one sweep's differences and its scores at
        each word's top."""
        var word = (row - 1) // WORD_BITS - top if row > top * WORD_BITS else 0
        if row == top * WORD_BITS:
            return bases[unsafe_offset=step * (count + 1)]
        var bits = row - (top + word) * WORD_BITS
        var kept = ALL_ONES if bits == WORD_BITS else (UInt64(1) << UInt64(bits)) - 1
        return bases[unsafe_offset=step * (count + 1) + word] + word_value(
            plus[unsafe_offset=step * count + word] & kept, minus[unsafe_offset=step * count + word] & kept
        )

    var plus = buffers.plus.unsafe_ptr()
    var minus = buffers.minus.unsafe_ptr()
    var bases = buffers.bases.unsafe_ptr()
    var low_plus = buffers.low_plus.unsafe_ptr()
    var low_minus = buffers.low_minus.unsafe_ptr()
    var low_bases = buffers.low_bases.unsafe_ptr()

    @inline(.always)
    def optimal(
        step: Int, row: Int, needed: Int
    ) {imm plus, imm minus, imm bases, imm low_plus, imm low_minus, imm low_bases, imm whole, imm score_in} -> Int:
        """1 when the cell's score is `needed`, so a path of it joins the trace optimally; 0 when it
        cannot be; -1 when this window cannot tell."""
        if score_in(plus, minus, bases, step, row) == needed:
            return 1
        if whole or score_in(low_plus, low_minus, low_bases, step, row) > needed:
            return 0
        return -1

    if not whole and score_in(plus, minus, bases, width, end_row) != score:
        return -1
    var lowest = top * WORD_BITS
    var start = len(moves)
    var step = width
    var row = end_row
    var current = score

    @inline(.always)
    def give_up(mut moves: List[UInt8]) {imm start}:
        """Drops the moves this window appended, so a taller one starts where it did."""
        moves.resize(start, 0)

    while step > 0:
        var column = first_column + step - 1
        var inside = row > lowest
        var differs = inside and profile.column_codes[column] != profile.row_codes[row - 1]
        if differs:
            var found = optimal(step - 1, row - 1, current - 1)
            if found < 0:
                give_up(moves)
                return -1
            if found == 1:
                moves.append(DIAGONAL)
                current -= 1
                step -= 1
                row -= 1
                continue
        var found = optimal(step - 1, row, current - 1)
        if found < 0:
            give_up(moves)
            return -1
        if found == 1:
            moves.append(LEFT)
            current -= 1
            step -= 1
            continue
        if not inside:
            if not whole:
                # A gap through the window's top may be the one to take: a taller window must tell.
                give_up(moves)
                return -1
            # The band's top: the path keeps to the band, so it goes on up the gap.
            moves.append(UP)
            current -= 1
            row -= 1
            continue
        found = optimal(step, row - 1, current - 1)
        if found < 0:
            give_up(moves)
            return -1
        if found == 1:
            moves.append(UP)
            current -= 1
            row -= 1
            continue
        # No edit enters the cell optimally, so it is a match the path came along.
        if differs or optimal(step - 1, row - 1, current) != 1:
            give_up(moves)
            return -1
        moves.append(DIAGONAL)
        step -= 1
        row -= 1
    return row


def trace_back(profile: Profile, trail: Trail, start_column: Int, start_row: Int, score: Int, mut moves: List[UInt8]):
    """An optimal path from `(start_column, start_row)`, scoring `score`, back to the origin, as moves right to left.

    One tile at a time from the last: each segment ends at a left-edge row whose recorded score is
    exact, because it plus the segment's cost is the exact score it started from, so the next tile
    starts from a known score. Each tile is traced by WFA2-lib's rule, by a forward wavefront across it
    (see `forward_segment`), or, when that would cost too much, swept again and traced cell by cell
    (see `window_segment`): the same path either way, whichever band found the distance.
    """
    var edge = Edge()
    var fronts = TileFronts()
    var buffers = Recompute()
    var column = start_column
    var row = start_row
    var current = score
    for tile in range(len(trail.first_columns) - 1, -1, -1):
        var first_column = trail.first_columns[tile]
        edge.load(trail, tile, profile.rows)
        # The tile's share of the distance, by its share of the columns: its wavefront may climb three
        # times that, or `WAVEFRONT_FLOOR`, before the tile is recomputed instead.
        var share = ceildiv(score * (column - first_column), max(start_column, 1))
        var limit = max(WAVEFRONT_FLOOR, 3 * share)
        var left = forward_segment(profile, edge, first_column, column, row, current, limit, fronts, moves)
        if left < 0:
            left = recomputed_segment(profile, trail, tile, column, row, current, buffers, moves)
        # Exact, since it plus the segment's cost is the exact score the segment started from.
        current = edge.score(left)
        column = first_column
        row = left
    # The first column is the border: the rest of the way up is gaps against the second sequence.
    for _ in range(row):
        moves.append(UP)


@fieldwise_init
struct EditPath(Movable):
    """An optimal alignment as traceback moves, as `edit_cigar` finds it: `prefix` from the origin to
    `(middle_column, middle_row)` right to left, as the traceback appends them, and `suffix` from
    there to the corner left to right, so neither is reversed first."""

    var prefix: List[UInt8]
    var suffix: List[UInt8]
    var middle_column: Int
    var middle_row: Int
    var distance: Int


@inline(.always)
def copy_bytes(destination: MutPointer[UInt8, _], source: ImmPointer[UInt8, _], count: Int):
    """`count` bytes, sixteen at a time while that many are left."""
    var index = 0
    while index + 16 <= count:
        destination.unsafe_offset(index).unsafe_store(source.unsafe_offset(index).unsafe_load[width=16]())
        index += 16
    while index < count:
        destination[unsafe_offset=index] = source[unsafe_offset=index]
        index += 1


@inline(.always)
def diagonal_run(moves: ImmPointer[UInt8, _], start: Int, end: Int) -> Int:
    """How many moves from `start` on, before `end`, are `DIAGONAL`, eight at a time."""
    var index = start
    while index + 8 <= end:
        var eight = moves.unsafe_offset(index).unsafe_bitcast[UInt64]().unsafe_load()
        if eight != 0:
            return index - start + Int(count_trailing_zeros(eight)) // 8
        index += 8
    while index < end and moves[unsafe_offset=index] == DIAGONAL:
        index += 1
    return index - start


struct CigarWriter:
    """A CIGAR string written a run at a time into bytes reserved once, a run of the same letter as the
    last one joining it; each length's digits are written by hand, as formatting one through a `String`,
    or growing the bytes a run at a time, took longer than the gapped rows' whole copy on short reads."""

    var text: List[UInt8]
    var used: Int
    var letter: UInt8
    """The letter of the run still being added to, zero before the first."""
    var length: Int
    """The length of that run so far."""

    def __init__(out self, capacity: Int):
        """Room for `capacity` bytes, which the caller bounds: nothing past it is checked."""
        self.text = List[UInt8](capacity=capacity)
        self.text.resize(unsafe_uninit_length=capacity)
        self.used = 0
        self.letter = 0
        self.length = 0

    @inline(.always)
    def add(mut self, letter: UInt8, length: Int):
        """`length` more of `letter`, joining the run being added to when it has the same letter."""
        if letter != self.letter:
            self.flush()
            self.letter = letter
        self.length += length

    def flush(mut self):
        """Writes the run being added to, if any, as its length's digits then its letter."""
        if self.length == 0:
            return
        var digits = 1
        var power = 10
        while power <= self.length:
            digits += 1
            power *= 10
        var at = self.used
        self.used += digits + 1
        var out = self.text.unsafe_ptr()
        var rest = self.length
        for place in range(digits - 1, -1, -1):
            out[unsafe_offset=at + place] = UInt8(ord("0") + rest % 10)
            rest //= 10
        out[unsafe_offset=at + digits] = self.letter
        self.length = 0

    def finish(var self) -> String:
        """The CIGAR string, its last run written."""
        self.flush()
        self.text.resize(self.used, 0)
        return String(unsafe_from_utf8=self.text)


@inline(.always)
def equal_run(first: ImmPointer[UInt8, _], second: ImmPointer[UInt8, _], column: Int, row: Int, limit: Int) -> Int:
    """How many of the next `limit` bases along the diagonal from `(column, row)` are equal, eight at a time."""
    var length = 0
    while length + 8 <= limit:
        var differing = (
            first.unsafe_offset(column + length).unsafe_bitcast[UInt64]().unsafe_load()
            ^ second.unsafe_offset(row + length).unsafe_bitcast[UInt64]().unsafe_load()
        )
        if differing != 0:
            return length + Int(count_trailing_zeros(differing)) // 8
        length += 8
    while length < limit and first[unsafe_offset=column + length] == second[unsafe_offset=row + length]:
        length += 1
    return length


def cigar_string(first: String, second: String, path: EditPath, eqx: Bool) -> String:
    """`path` as a CIGAR string (see `EditCigar`): its moves put left to right, the prefix reversed, and
    each run of one move written as one entry, a diagonal run split into its matches and substitutions
    by comparing the bases eight at a time unless `M` stands for both."""
    var before = len(path.prefix)
    var count = before + len(path.suffix)
    var moves = List[UInt8](capacity=count)
    moves.resize(unsafe_uninit_length=count)
    var ordered = moves.unsafe_ptr()
    var prefix = path.prefix.unsafe_ptr()
    comptime CHUNK = 16
    var index = 0
    while index + CHUNK <= before:
        var chunk = prefix.unsafe_offset(before - index - CHUNK).unsafe_load[width=CHUNK]()
        ordered.unsafe_offset(index).unsafe_store(chunk.reversed())
        index += CHUNK
    while index < before:
        ordered[unsafe_offset=index] = prefix[unsafe_offset=before - 1 - index]
        index += 1
    copy_bytes(ordered.unsafe_offset(before), path.suffix.unsafe_ptr(), len(path.suffix))

    var first_bytes = first.unsafe_ptr()
    var second_bytes = second.unsafe_ptr()
    # At most two runs an edit and one more, each its length's digits and a letter.
    var digits = 1
    var power = 10
    while power <= max(first.byte_length(), second.byte_length()):
        digits += 1
        power *= 10
    var writer = CigarWriter((digits + 1) * (2 * path.distance + 2))
    var column = 0
    var row = 0
    index = 0
    while index < count:
        var move = ordered[unsafe_offset=index]
        if move == DIAGONAL:
            var run = diagonal_run(ordered, index, count)
            index += run
            if not eqx:
                writer.add(UInt8(ord("M")), run)
                column += run
                row += run
                continue
            var end = column + run
            while column < end:
                var same = equal_run(first_bytes, second_bytes, column, row, end - column)
                if same > 0:
                    writer.add(UInt8(ord("=")), same)
                    column += same
                    row += same
                var differ = 0
                while (
                    column + differ < end
                    and first_bytes[unsafe_offset=column + differ] != second_bytes[unsafe_offset=row + differ]
                ):
                    differ += 1
                if differ > 0:
                    writer.add(UInt8(ord("X")), differ)
                    column += differ
                    row += differ
            continue
        var start = index
        while index < count and ordered[unsafe_offset=index] == move:
            index += 1
        if move == LEFT:
            writer.add(UInt8(ord("D")), index - start)
            column += index - start
        else:
            writer.add(UInt8(ord("I")), index - start)
            row += index - start
    return writer^.finish()


def diagonal_cigar(profile: Profile, fronts: DiagonalFronts, distance: Int, reversed: Bool, eqx: Bool) -> String:
    """The CIGAR of the path `trace_diagonals` traces through `fronts`, written straight from them: each
    score undoes the matches its front slid over, `=`, and the edit that reached its start, a
    substitution `X`, as one after the furthest point is a mismatch, or a gap; `M` for both kinds of
    pair without `eqx`. The runs come right to left, unless the fronts are the reversed pair's,
    `reversed`, whose right to left is the pair's left to right."""
    var columns = profile.columns
    var rows = profile.rows
    var matched = UInt8(ord("=")) if eqx else UInt8(ord("M"))
    var substituted = UInt8(ord("X")) if eqx else UInt8(ord("M"))
    var letters = List[UInt8](capacity=2 * distance + 2)
    var lengths = List[Int](capacity=2 * distance + 2)
    var diagonal = columns - rows
    var column = columns
    var score = distance
    while score > 0:
        var source = best_source(fronts, score, diagonal, columns, rows)
        var best = source[0]
        if column > best:
            letters.append(matched)
            lengths.append(column - best)
        var move = source[1]
        if move == DIAGONAL:
            letters.append(substituted)
            column = best - 1
        elif move == LEFT:
            letters.append(UInt8(ord("D")))
            column = best - 1
            diagonal -= 1
        else:
            letters.append(UInt8(ord("I")))
            column = best
            diagonal += 1
        lengths.append(1)
        score -= 1
    if column > 0:
        letters.append(matched)
        lengths.append(column)
    var digits = 1
    var power = 10
    while power <= max(columns, rows):
        digits += 1
        power *= 10
    var writer = CigarWriter((digits + 1) * (2 * distance + 2))
    var count = len(letters)
    for index in range(count):
        var at = index if reversed else count - 1 - index
        writer.add(letters[at], lengths[at])
    return writer^.finish()
