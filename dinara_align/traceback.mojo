# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Ported from `pa-bitpacking` in A*PA (https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner,
# commit bf2e14e), by Ragnar Groot Koerkamp and Pesho Ivanov, itself translated from Edlib.
"""
An optimal alignment from a band's recorded tile edges, retraced a tile at a time from the corner (see
`trace_back`), and written out as two gapped rows or as a CIGAR.
"""

from std.bit import count_leading_zeros, count_trailing_zeros
from std.math import ceildiv, clamp

from .alignment import AlignmentResult
from .bit_parallel import (
    advance,
    ALL_ONES,
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
from .seeds import differing


comptime WAVEFRONT_FLOOR = 32
"""The fewest edits a tile's wavefront search may spend before the tile is recomputed instead."""


comptime TRACE_PADDING = 2
"""Unreached diagonals stored either side of a traceback front, so the next reads its neighbours unchecked."""


comptime FRONT_DROP = 20
"""
How far behind the furthest front, in column plus row, a front may fall before the search drops it;
about ten diagonal steps, A*PA2's `fr_drop`.
"""


struct Wavefronts(Movable):
    """The wavefront search's buffers, kept across tiles so a tile allocates nothing.

    Cost `s` keeps its diagonals from offset `lows[s]` from the start's, padding included, at position
    `starts[s]` of the flat buffers on.
    """

    var landed: List[Int]
    var reached: List[Int]
    var how: List[UInt8]
    var starts: List[Int]
    var lows: List[Int]

    def __init__(out self):
        # Room for a typical tile up front: growing from empty cost more than the search itself.
        comptime DIAGONALS = 1024
        comptime COSTS = 64
        self.landed = List[Int](capacity=DIAGONALS)
        self.reached = List[Int](capacity=DIAGONALS)
        self.how = List[UInt8](capacity=DIAGONALS)
        self.starts = List[Int](capacity=COSTS)
        self.lows = List[Int](capacity=COSTS)


@inline(.always)
def slide(
    first: ImmPointer[UInt8, _], second: ImmPointer[UInt8, _], column: Int, diagonal: Int, first_column: Int
) -> Int:
    """How far back along a diagonal matches carry `column`, eight bases at a time while both have eight.

    The eight bases ending just before the position are loaded from each sequence as one word; their
    XOR is zero in every matching byte, so the matching run is the count of zero bytes from the top.
    """
    var at = column
    var row = at + diagonal
    while at - first_column >= 8 and row >= 8:
        var mismatches = (
            first.unsafe_offset(at - 8).unsafe_bitcast[UInt64]().unsafe_load()
            ^ second.unsafe_offset(row - 8).unsafe_bitcast[UInt64]().unsafe_load()
        )
        if mismatches != 0:
            var run = Int(count_leading_zeros(mismatches)) // 8
            return at - run
        at -= 8
        row -= 8
    while at > first_column and row > 0 and first[unsafe_offset=at - 1] == second[unsafe_offset=row - 1]:
        at -= 1
        row -= 1
    return at


def wavefront_segment(
    profile: Profile,
    edge: Edge,
    first_column: Int,
    end_column: Int,
    end_row: Int,
    score: Int,
    limit: Int,
    mut fronts: Wavefronts,
    mut moves: List[UInt8],
) -> Int:
    """Traces from `(end_column, end_row)`, scoring `score`, back to the tile's left edge by wavefront.

    Diagonal transition run backwards, as A*PA2's trace does: for each cost `s`, the furthest a path
    of `s` edits reaches back along every live diagonal, sliding over matches for free. A path is
    accepted only when it reaches the left edge at a row whose score plus `s` is `score`, so it is
    optimal and joins the trace exactly.

    As in A*PA2, the live diagonals shrink from both ends while the outermost front has fallen more
    than `FRONT_DROP` behind the furthest or already stands on the left edge, and the search gives up
    when half its budget has not carried any front halfway across. Each only ever turns a success
    into a miss, and a miss falls back to recomputing the tile. Appends the segment's moves right to
    left and returns its left-edge row, or -1 on a miss.
    """
    var first = profile.column_codes.unsafe_ptr()
    var second = profile.row_codes.unsafe_ptr()
    var home = end_row - end_column
    # Unreached reads as a column far left of any tile, so every move from it fails its own test
    # unchecked; each stored row carries `TRACE_PADDING` of it either side, so the next row reads its
    # neighbours without range checks.
    comptime UNREACHED = -(1 << 40)
    comptime NONE = Int.MAX

    fronts.landed.clear()
    fronts.reached.clear()
    fronts.how.clear()
    fronts.starts.clear()
    fronts.lows.clear()

    @inline(.always)
    def ends_here(column: Int, offset: Int, cost: Int) {imm edge, imm first_column, imm score, imm home} -> Bool:
        if column != first_column:
            return False
        var row = column + home + offset
        return row >= edge.low_row and row <= edge.high_row and edge.score(row) + cost == score

    @inline(.always)
    def row(mut fronts: Wavefronts, low: Int, high: Int) -> Int:
        """Lays out the next cost's row over diagonals `low ..= high`, its padding unreached, and
        returns where the row's diagonal zero would sit, so diagonal `k` is that plus `k`."""
        var start = len(fronts.reached)
        var width = high - low + 1 + 2 * TRACE_PADDING
        fronts.starts.append(start)
        fronts.lows.append(low - TRACE_PADDING)
        fronts.reached.resize(unsafe_uninit_length=start + width)
        fronts.landed.resize(unsafe_uninit_length=start + width)
        fronts.how.resize(unsafe_uninit_length=start + width)
        var zero = start + TRACE_PADDING - low
        var reached = fronts.reached.unsafe_ptr()
        for index in range(1, TRACE_PADDING + 1):
            reached[unsafe_offset=zero + low - index] = UNREACHED
            reached[unsafe_offset=zero + high + index] = UNREACHED
        return zero

    var finish_cost = -1
    var finish_offset = 0
    var zero = row(fronts, 0, 0)
    fronts.landed[zero] = end_column
    fronts.reached[zero] = slide(first, second, end_column, home, first_column)
    fronts.how[zero] = DIAGONAL
    if ends_here(fronts.reached[zero], 0, 0):
        finish_cost = 0
    var low = 0
    var high = 0
    var cost = 1
    var halfway = (first_column + end_column) // 2
    while finish_cost < 0 and cost <= limit:
        low -= 1
        high += 1
        var current_zero = row(fronts, low, high)
        # The previous row is laid out from its own first stored diagonal, unreached outside its live ones.
        var previous = fronts.reached.unsafe_ptr().unsafe_offset(zero)
        var current = fronts.reached.unsafe_ptr().unsafe_offset(current_zero)
        var landed = fronts.landed.unsafe_ptr().unsafe_offset(current_zero)
        var how = fronts.how.unsafe_ptr().unsafe_offset(current_zero)
        zero = current_zero
        var furthest = NONE
        var nearest_column = NONE
        for offset in range(low, high + 1):
            var same = previous[unsafe_offset=offset]
            var below = previous[unsafe_offset=offset - 1]
            var above = previous[unsafe_offset=offset + 1]
            var best = NONE
            var move = DIAGONAL
            var diagonal = home + offset
            # Substitution, along the same diagonal.
            if same > first_column and same + diagonal > 0:
                best = same - 1
            # A base of the first sequence against a gap, from the diagonal below.
            if below > first_column and below - 1 < best:
                best = below - 1
                move = LEFT
            # A base of the second sequence against a gap, from the diagonal above.
            if above + diagonal + 1 > 0 and above < best:
                best = above
                move = UP
            how[unsafe_offset=offset] = move
            if best == NONE:
                landed[unsafe_offset=offset] = UNREACHED
                current[unsafe_offset=offset] = UNREACHED
                continue
            landed[unsafe_offset=offset] = best
            var slid = slide(first, second, best, diagonal, first_column)
            current[unsafe_offset=offset] = slid
            furthest = min(furthest, 2 * slid - offset)
            nearest_column = min(nearest_column, slid)
            if finish_cost < 0 and ends_here(slid, offset, cost):
                finish_cost = cost
                finish_offset = offset
        if finish_cost >= 0:
            break
        if furthest == NONE:
            return -1
        if 2 * cost >= limit and nearest_column > halfway:
            return -1

        @inline(.always)
        def dropped(offset: Int) {imm current, imm furthest, imm first_column} -> Bool:
            var at = current[unsafe_offset=offset]
            return at <= first_column or 2 * at - offset > furthest + FRONT_DROP

        # Shrink the live diagonals from both ends, marking the dropped ones unreached in place.
        var new_low = low
        var new_high = high
        while new_low < new_high and dropped(new_low):
            current[unsafe_offset=new_low] = UNREACHED
            new_low += 1
        while new_high > new_low and dropped(new_high):
            current[unsafe_offset=new_high] = UNREACHED
            new_high -= 1
        if dropped(new_low):
            return -1
        low = new_low
        high = new_high
        cost += 1
    if finish_cost < 0:
        return -1

    # Walk the costs back down, left to right: once to size the segment, an edit per cost plus the
    # matches each slid over, then again writing its moves from the block's far end, since they go
    # after the moves already traced right to left.
    var length = finish_cost
    var offset = finish_offset
    var at_cost = finish_cost
    while True:
        var index = fronts.starts[at_cost] + offset - fronts.lows[at_cost]
        length += fronts.landed[index] - fronts.reached[index]
        if at_cost == 0:
            break
        var move = fronts.how[index]
        if move == LEFT:
            offset -= 1
        elif move == UP:
            offset += 1
        at_cost -= 1
    var block = len(moves)
    moves.resize(unsafe_uninit_length=block + length)
    var out = moves.unsafe_ptr().unsafe_offset(block + length)
    offset = finish_offset
    at_cost = finish_cost
    while True:
        var index = fronts.starts[at_cost] + offset - fronts.lows[at_cost]
        var run = fronts.landed[index] - fronts.reached[index]
        out = out.unsafe_offset(-run)
        for step in range(run):
            out[unsafe_offset=step] = DIAGONAL
        if at_cost == 0:
            break
        var move = fronts.how[index]
        out = out.unsafe_offset(-1)
        out[] = move
        if move == LEFT:
            offset -= 1
        elif move == UP:
            offset += 1
        at_cost -= 1
    return first_column + home + finish_offset


comptime RECOMPUTE_WORDS = 4
"""Words a tile's recompute first takes above the point it traces from, doubling while that falls short."""


struct Recompute(Movable):
    """A tile recompute's buffers, kept across tiles so a recompute neither allocates nor zeroes: every
    slot it reads it has written first."""

    var plus: List[UInt64]
    var minus: List[UInt64]
    var bases: List[Int]

    def __init__(out self):
        self.plus = List[UInt64]()
        self.minus = List[UInt64]()
        self.bases = List[Int]()


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
    """`recomputed_segment` over words `[first_word, end_word)` alone, or -1 when they do not hold the path.

    The window's left edge comes from the recorded one, and its top reads `+1` from above, as a
    band's top does, so every score is a real path's and never below the true one. When the
    recomputed score at the traced cell is its exact `score`, a path back to the left edge along
    those scores is optimal; the window falls short when the scores differ, or when the path would
    leave through its top. Below the band's own top, which is the full recompute, neither can happen.
    """
    var first_column = trail.first_columns[tile]
    var top = first_word
    var count = end_word - first_word
    var width = end_column - first_column
    var offset = trail.offsets[tile] + first_word - trail.tops[tile]
    buffers.plus.resize(unsafe_uninit_length=(width + 1) * count)
    buffers.minus.resize(unsafe_uninit_length=(width + 1) * count)
    buffers.bases.resize(unsafe_uninit_length=(width + 1) * (count + 1))
    var plus = buffers.plus.unsafe_ptr()
    var minus = buffers.minus.unsafe_ptr()
    var bases = buffers.bases.unsafe_ptr()
    # The left edge's score at the window's top, carried down from the band's.
    var anchor = trail.anchors[tile]
    for word in range(trail.offsets[tile], offset):
        anchor += word_value(trail.edge_plus[word], trail.edge_minus[word])
    for word in range(count):
        plus[unsafe_offset=word] = trail.edge_plus[offset + word]
        minus[unsafe_offset=word] = trail.edge_minus[offset + word]
    var column_low = profile.column_low.unsafe_ptr().unsafe_offset(COLUMN_PADDING + first_column - 1)
    var column_high = profile.column_high.unsafe_ptr().unsafe_offset(COLUMN_PADDING + first_column - 1)
    var row_low = profile.row_low.unsafe_ptr().unsafe_offset(top)
    var row_high = profile.row_high.unsafe_ptr().unsafe_offset(top)
    # The third plane exists only for symbols past `ACGT` (see `Profile.extended`).
    var column_extra = profile.column_extra.unsafe_ptr().unsafe_offset(COLUMN_PADDING + first_column - 1)
    var row_extra = profile.row_extra.unsafe_ptr().unsafe_offset(top)
    for step in range(1, width + 1):
        var horizontal_plus = UInt64(1)
        var horizontal_minus = UInt64(0)
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

    # Scores at each word's top on every column: the window's top scores the anchor plus one per column.
    for step in range(width + 1):
        var running = anchor + step
        bases[unsafe_offset=step * (count + 1)] = running
        for word in range(count):
            running += word_value(plus[unsafe_offset=step * count + word], minus[unsafe_offset=step * count + word])
            bases[unsafe_offset=step * (count + 1) + word + 1] = running

    @inline(.always)
    def score_at(step: Int, row: Int) {imm bases, imm plus, imm minus, imm count, imm top} -> Int:
        var word = (row - 1) // WORD_BITS - top if row > top * WORD_BITS else 0
        if row == top * WORD_BITS:
            return bases[unsafe_offset=step * (count + 1)]
        var bits = row - (top + word) * WORD_BITS
        var kept = ALL_ONES if bits == WORD_BITS else (UInt64(1) << UInt64(bits)) - 1
        return bases[unsafe_offset=step * (count + 1) + word] + word_value(
            plus[unsafe_offset=step * count + word] & kept, minus[unsafe_offset=step * count + word] & kept
        )

    var whole = first_word == trail.tops[tile]
    if not whole and score_at(width, end_row) != score:
        return -1
    var lowest = top * WORD_BITS
    var start = len(moves)
    var step = width
    var row = end_row
    var current = score
    while step > 0:
        var column = first_column + step - 1
        if row > lowest:
            var change = 0 if profile.column_codes[column] == profile.row_codes[row - 1] else 1
            if score_at(step - 1, row - 1) + change == current:
                moves.append(DIAGONAL)
                current -= change
                step -= 1
                row -= 1
                continue
        if score_at(step - 1, row) + 1 == current:
            moves.append(LEFT)
            current -= 1
            step -= 1
            continue
        if not whole and row == lowest:
            # The path leaves through the window's top: a taller window must hold it.
            moves.resize(start, 0)
            return -1
        moves.append(UP)
        current -= 1
        row -= 1
    return row


def trace_back(profile: Profile, trail: Trail, start_column: Int, start_row: Int, score: Int, mut moves: List[UInt8]):
    """An optimal path from `(start_column, start_row)`, scoring `score`, back to the origin, as moves right to left.

    One tile at a time from the last: each segment ends at a left-edge row whose recorded score is
    exact, because it plus the segment's cost is the exact score it started from, so the next tile
    starts from a known score. The wavefront search allows a few times the tile's share of the
    score in edits before handing the tile to the recompute.
    """
    var edge = Edge()
    var fronts = Wavefronts()
    var buffers = Recompute()
    var column = start_column
    var row = start_row
    var current = score
    for tile in range(len(trail.first_columns) - 1, -1, -1):
        var first_column = trail.first_columns[tile]
        edge.load(trail, tile, profile.rows)
        var share = ceildiv(score * (column - first_column), max(start_column, 1))
        var limit = max(WAVEFRONT_FLOOR, 3 * share)
        var left = wavefront_segment(profile, edge, first_column, column, row, current, limit, fronts, moves)
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
    """An optimal alignment as traceback moves, as `edit_path` finds it: `prefix` from the origin to
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


def gapped_rows(
    first: String,
    second: String,
    prefix: List[UInt8],
    middle: Int,
    meeting_row: Int,
    suffix: List[UInt8],
    distance: Int,
) -> AlignmentResult:
    """An alignment written out as the two gapped rows.

    `prefix` holds the moves from the origin to `(middle, meeting_row)` right to left, as the
    traceback appends them, and `suffix` the moves from there to the corner left to right, so
    neither is reversed first. A run of diagonal moves, most of any alignment, copies both
    sequences' bases a block at a time; a gap move writes one base against `-`.
    """
    comptime GAP = UInt8(ord("-"))
    var length = len(prefix) + len(suffix)
    var top_row = List[UInt8](capacity=length)
    var bottom_row = List[UInt8](capacity=length)
    top_row.resize(unsafe_uninit_length=length)
    bottom_row.resize(unsafe_uninit_length=length)
    var first_bytes = first.unsafe_ptr()
    var second_bytes = second.unsafe_ptr()
    var top = top_row.unsafe_ptr()
    var bottom = bottom_row.unsafe_ptr()

    # The prefix, from its last move back to its first, fills the rows from `len(prefix)` down.
    var moves = prefix.unsafe_ptr()
    var count = len(prefix)
    var column = middle
    var row = meeting_row
    var at = count
    var index = 0
    while index < count:
        var run = diagonal_run(moves, index, count)
        if run > 0:
            at -= run
            column -= run
            row -= run
            copy_bytes(top.unsafe_offset(at), first_bytes.unsafe_offset(column), run)
            copy_bytes(bottom.unsafe_offset(at), second_bytes.unsafe_offset(row), run)
            index += run
            continue
        var move = moves[unsafe_offset=index]
        at -= 1
        if move == LEFT:
            column -= 1
            top[unsafe_offset=at] = first_bytes[unsafe_offset=column]
            bottom[unsafe_offset=at] = GAP
        else:
            row -= 1
            top[unsafe_offset=at] = GAP
            bottom[unsafe_offset=at] = second_bytes[unsafe_offset=row]
        index += 1

    # The suffix, first move first, fills the rest.
    var later = suffix.unsafe_ptr()
    count = len(suffix)
    column = middle
    row = meeting_row
    at = len(prefix)
    index = 0
    while index < count:
        var run = diagonal_run(later, index, count)
        if run > 0:
            copy_bytes(top.unsafe_offset(at), first_bytes.unsafe_offset(column), run)
            copy_bytes(bottom.unsafe_offset(at), second_bytes.unsafe_offset(row), run)
            at += run
            column += run
            row += run
            index += run
            continue
        if later[unsafe_offset=index] == LEFT:
            top[unsafe_offset=at] = first_bytes[unsafe_offset=column]
            bottom[unsafe_offset=at] = GAP
            column += 1
        else:
            top[unsafe_offset=at] = GAP
            bottom[unsafe_offset=at] = second_bytes[unsafe_offset=row]
            row += 1
        at += 1
        index += 1
    return AlignmentResult(Int32(distance), String(unsafe_from_utf8=top_row), String(unsafe_from_utf8=bottom_row))


struct CigarWriter:
    """A CIGAR string written a run at a time into bytes reserved once, a run of the same letter as the
    last one joining it; each length's digits are written by hand, as formatting one through a `String`,
    or growing the bytes a run at a time, took longer than the gapped rows' whole copy on short reads."""

    var text: List[UInt8]
    var used: Int
    var letter: UInt8
    var length: Int

    def __init__(out self, capacity: Int):
        """Room for `capacity` bytes, which the caller bounds: nothing past it is checked."""
        self.text = List[UInt8](capacity=capacity)
        self.text.resize(unsafe_uninit_length=capacity)
        self.used = 0
        self.letter = 0
        self.length = 0

    @inline(.always)
    def add(mut self, letter: UInt8, length: Int):
        if letter != self.letter:
            self.flush()
            self.letter = letter
        self.length += length

    def flush(mut self):
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


def cigar_string(first: String, second: String, path: EditPath, extended: Bool) -> String:
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
            if not extended:
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
