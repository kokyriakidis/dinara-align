"""
The global affine-gap score by wavefront, as WFA computes it, for a table with one match score and
one mismatch score.

The Gotoh recurrence maximizes a score with a reward for every match, which a wavefront cannot
search: it grows by cost, and a match must cost nothing. For a global alignment the reward folds
away (Eizenga and Lindquist, 2022). A path over sequences of `n` and `m` letters with `M`
matches, `X` mismatches, `g` gap runs and `G` gapped letters spends `2M + 2X + G = n + m`, so with
match `a`, mismatch `b`, and a run of `k` costing `o + k e` under the gap penalties,

    a (n + m) - 2 score = 2 (a - b) X + 2 o g + (2 e + a) G,

a sum of non-negative penalties that are zero for a match. Minimizing it is maximizing the score:
mismatch `x = 2 (a - b)`, opening `2 o` and extension `2 e + a`, all divided by their common factor.

The wavefront then runs as in WFA (Marco-Sola et al., 2021): for each cost `s`, three fronts hold
the furthest column every diagonal reaches at exactly that cost, ending in an alignment of two
letters or a gap either way, each built from the fronts `x`, `o + e` and `e` costs back, and the
alignment front slides over matches for free. The work grows with the square of the cost, not with
the matrix, so a close pair costs a small share of a full sweep.

It runs from both ends at once, as BiWFA does (Marco-Sola et al., 2023): a search from the origin
and one from the corner over the reversed sequences, a cost at a time each, until where they
overlap proves the optimum (see `bidirectional`). Each grows to about half the cost, so the two
together step about half the diagonals one search would.

For an alignment both searches keep, of every cost, the alignment front's columns and a byte of
which source won each layer (see `History`), five bytes a diagonal where WFA's high-memory mode
keeps twelve, and the path is traced from where they met back to each end (see `trace`). A pair
whose fronts would grow past `HISTORY_LIMIT` is split where an optimal path crosses, which the two
searches find keeping only their last few costs, and each piece is aligned the same way (see
`solve`), so the memory stays bounded whatever the pair.
"""

from std.bit import count_trailing_zeros

from .errors import AlignmentError, ErrorKind
from .traceback import cigar_string, EditPath

comptime Slot = MutPointer[Int32, MutUntrackedOrigin]
"""A front read or written by `step`, indexed by diagonal: a source and its destination may share a row."""

comptime UNREACHED = Int32(-(1 << 28))
"""A diagonal no path of the cost reaches: far enough below zero that a few more columns stay negative."""

comptime PADDING = 16
"""Sentinels after each sequence's codes, so a slide reads eight letters at a time past the end."""

comptime FIRST_SENTINEL = UInt8(0xFE)
"""Ends the first sequence's codes; never a letter, and never the second's sentinel."""

comptime SECOND_SENTINEL = UInt8(0xFF)
"""Ends the second sequence's codes."""

comptime LANES = 8
"""Diagonals a front step computes at once."""

comptime CHECK_STRIDE = 64
"""Costs between the search's checks of its projected work."""

comptime CHECK_START = 128
"""The cost from which the search judges its projection, so it has a few edits to go on."""

comptime HISTORY_LIMIT = 1 << 24
"""Diagonals the two searches' kept fronts may hold together, 80 MB, past which the pair is split
(see `solve`)."""

comptime ALIGNED = 0
"""The front ending in two letters aligned, a match or a substitution."""
comptime FIRST_GAP = 1
"""The front ending in a letter of the first sequence against a gap."""
comptime SECOND_GAP = 2
"""The front ending in a letter of the second sequence against a gap."""

comptime CELLS_PER_STEP = 4
"""Cells of the vectorized full sweep (see `vector_score`) one diagonal step of the three fronts
costs about as much as: about 1.5 against 0.4 ns."""

comptime FREE_START = 0
"""A search's origin as a whole alignment's: any first move, each at its own cost."""
comptime IN_FIRST_GAP = 1
"""The origin inside a gap of the first sequence's letters, which a first move of one continues for
an extension alone: the piece after a split inside that gap."""
comptime IN_SECOND_GAP = 2
"""The origin inside a gap of the second sequence's letters."""
comptime OPENING_FIRST_GAP = 3
"""The first move a letter of the first sequence against a gap, opened there: the backward search of
the piece before a split inside that gap, which must end in it."""
comptime OPENING_SECOND_GAP = 4
"""The first move a letter of the second sequence against a gap, opened there."""


@fieldwise_init
struct Penalties(ImplicitlyCopyable, TrivialRegisterPassable):
    """The wavefront's costs for one scoring, and what turns a cost back into a score."""

    var mismatch: Int
    var opening: Int
    """Charged once a gap run, on top of `extension` for its first letter."""
    var extension: Int
    var scale: Int
    """The common factor the costs were divided by."""
    var reward: Int
    """The match score, which every letter of both sequences earns half of before the costs."""

    def score(self, cost: Int, letters: Int) -> Int:
        """The Gotoh score of an alignment over `letters` letters in all, costing `cost`."""
        return (self.reward * letters - cost * self.scale) // 2


def greatest_common_divisor(first: Int, second: Int) -> Int:
    var a = first
    var b = second
    while b != 0:
        var rest = a % b
        a = b
        b = rest
    return a


def wavefront_penalties(
    substitutions: List[Scalar[DType.int8]], alphabet_size: Int, open: Int, extend: Int
) -> Optional[Penalties]:
    """The wavefront costs for a table of one match and one mismatch score, if it is one.

    `open` and `extend` are the Gotoh gap scores, a run of `k` scoring `open + (k - 1) extend`.
    None when the table is not uniform, or when the costs would not all be positive: a free
    mismatch or extension would let a front grow without its cost growing.
    """
    if alphabet_size < 2:
        return None
    var reward = Int(substitutions[0])
    var mismatch = Int(substitutions[1])
    for row in range(alphabet_size):
        for column in range(alphabet_size):
            var expected = reward if row == column else mismatch
            if Int(substitutions[row * alphabet_size + column]) != expected:
                return None
    # A run of `k` scores `open + (k - 1) extend`, so it costs `(extend - open) - k extend`.
    var gap_opening = extend - open
    var gap_extension = -extend
    var x = 2 * (reward - mismatch)
    var o = 2 * gap_opening
    var e = 2 * gap_extension + reward
    if x <= 0 or e <= 0 or o < 0:
        return None
    var scale = greatest_common_divisor(greatest_common_divisor(x, e), o)
    return Penalties(x // scale, o // scale, e // scale, scale, reward)


def padded(codes: Span[UInt8, _], sentinel: UInt8, reverse: Bool) -> List[UInt8]:
    """The codes, back to front with `reverse`, and `PADDING` sentinels after them."""
    var count = len(codes)
    var out = List[UInt8](length=count + PADDING, fill=sentinel)
    var source = codes.unsafe_ptr()
    var destination = out.unsafe_ptr()
    if not reverse:
        Span(unsafe_ptr=destination, length=count).copy_from(codes)
        return out^
    comptime CHUNK = 16
    var index = 0
    while index + CHUNK <= count:
        var chunk = source.unsafe_offset(count - index - CHUNK).unsafe_load[width=CHUNK]()
        destination.unsafe_offset(index).unsafe_store(chunk.reversed())
        index += CHUNK
    while index < count:
        destination[unsafe_offset=index] = source[unsafe_offset=count - 1 - index]
        index += 1
    return out^


@inline(.always)
def slide(first: ImmPointer[UInt8, _], second: ImmPointer[UInt8, _], start: Int, diagonal: Int) -> Int:
    """How far matches carry column `start` of `diagonal`, eight letters at a time, to the sentinels."""
    var column = start
    var lag = second.unsafe_offset(-diagonal)
    var mismatches = (
        first.unsafe_offset(column).unsafe_bitcast[UInt64]().unsafe_load()
        ^ lag.unsafe_offset(column).unsafe_bitcast[UInt64]().unsafe_load()
    )
    while mismatches == 0:
        column += 8
        mismatches = (
            first.unsafe_offset(column).unsafe_bitcast[UInt64]().unsafe_load()
            ^ lag.unsafe_offset(column).unsafe_bitcast[UInt64]().unsafe_load()
        )
    return column + (Int(count_trailing_zeros(mismatches)) >> 3)


comptime FROM_FIRST_GAP = UInt8(1)
"""A flag's low two bits: the aligned front's entry came from the first sequence's gap layer; zero
for a substitution, `FROM_SECOND_GAP` for the second's."""
comptime FROM_SECOND_GAP = UInt8(2)
comptime FIRST_OPENED = UInt8(4)
"""A flag bit: the first sequence's gap layer came from an opening, not an extension."""
comptime SECOND_OPENED = UInt8(8)
"""A flag bit: the second sequence's gap layer came from an opening, not an extension."""


struct History(Movable):
    """What the traceback needs of every cost's fronts: the alignment front's column on each diagonal
    kept, and a flag of which source won each layer there (see `FROM_FIRST_GAP` and the bits after
    it). Cost `s` holds diagonals `lows[s] ..= highs[s]`, its columns from `starts[s]` and its flags
    from `flag_starts[s]`; elsewhere it reads as unreached."""

    var starts: List[Int]
    var flag_starts: List[Int]
    var lows: List[Int]
    var highs: List[Int]
    var aligned: List[Int32]
    var flags: List[UInt8]

    def __init__(out self):
        self.starts = List[Int]()
        self.flag_starts = List[Int]()
        self.lows = List[Int]()
        self.highs = List[Int]()
        self.aligned = List[Int32]()
        self.flags = List[UInt8]()

    def flags_for(mut self, low: Int, high: Int) -> MutPointer[UInt8, MutUntrackedOrigin]:
        """Room for the next cost's flags on `low ..= high` and a lane group past, indexed by diagonal."""
        var start = len(self.flags)
        # The step writes every flag of the range before any is read.
        self.flags.resize(unsafe_uninit_length=start + high - low + 1 + LANES)
        return self.flags.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]().unsafe_offset(start - low)

    def record(mut self, low: Int, high: Int, flagged_low: Int, flagged_high: Int, aligned: Slot):
        """The next cost's alignment front on `low ..= high`, none when `low > high`, whose flags
        `flags_for` placed for `flagged_low ..= flagged_high`."""
        var flagged_start = len(self.flags) - LANES - (flagged_high - flagged_low + 1)
        self.starts.append(len(self.aligned))
        self.lows.append(low)
        self.highs.append(high)
        var count = high - low + 1
        if count <= 0:
            self.flags.shrink(flagged_start)
            self.flag_starts.append(flagged_start)
            return
        self.flag_starts.append(flagged_start + low - flagged_low)
        self.flags.shrink(flagged_start + high - flagged_low + 1)
        var at = len(self.aligned)
        self.aligned.resize(unsafe_uninit_length=at + count)
        Span(unsafe_ptr=self.aligned.unsafe_ptr().unsafe_offset(at), length=count).copy_from(
            Span(unsafe_ptr=aligned.unsafe_offset(low), length=count)
        )

    def skip(mut self):
        """The next cost reached nothing."""
        self.starts.append(len(self.aligned))
        self.flag_starts.append(len(self.flags))
        self.lows.append(1)
        self.highs.append(0)

    @inline(.always)
    def column(self, cost: Int, diagonal: Int) -> Int:
        """The alignment front's furthest column on `diagonal` at `cost`, or `UNREACHED`."""
        if cost < 0 or cost >= len(self.lows) or diagonal < self.lows[cost] or diagonal > self.highs[cost]:
            return Int(UNREACHED)
        return Int(self.aligned[self.starts[cost] + diagonal - self.lows[cost]])

    @inline(.always)
    def flag(self, cost: Int, diagonal: Int) -> UInt8:
        """The flag of `diagonal` at `cost`, which the traceback reads only where a front reached."""
        return self.flags[self.flag_starts[cost] + diagonal - self.lows[cost]]


struct Fronts(Movable):
    """The last `slots` costs' three fronts, in rings indexed by cost, each over the diagonals reached.

    Diagonal `k` holds the cells whose column minus row is `k`, from `-rows` to `columns`, stored
    at `k + base` of each of the `3 slots` rows of one buffer, a slot's layers side by side. The rows
    cover only the diagonals the fronts have reached, about `2 s / e` of them at cost `s`, and double
    towards whichever side they outgrow, so the memory grows with the score rather than the
    sequences. A slot holds its cost's values on `lows[slot] ..= highs[slot]` and reads as unreached
    everywhere else, so a front reads its sources' neighbours unchecked: the rows start unreached,
    and a slot reused for a later cost has the old cost's diagonals outside the new range cleared.
    """

    var slots: Int
    var current: Int
    """The slot of the last cost grown."""
    var width: Int
    var base: Int
    var buffer: List[Int32]
    var lows: List[Int]
    var highs: List[Int]
    var reach: List[Int]
    """Each slot's furthest anti-diagonal, `2 column - diagonal`, over its alignment front."""
    var least: Int
    """The diagonals a front may ever read: the matrix's and two either side, and a lane group past."""
    var most: Int

    def __init__(out self, slots: Int, columns: Int, rows: Int):
        self.slots = slots
        self.current = 0
        self.least = -rows - 2
        self.most = columns + LANES + 2
        self.width = 4 * LANES
        self.base = self.width // 2
        self.buffer = List[Int32](length=3 * slots * self.width, fill=UNREACHED)
        self.lows = List[Int](length=slots, fill=1)
        self.highs = List[Int](length=slots, fill=0)
        self.reach = List[Int](length=slots, fill=Int.MIN // 2)

    @inline(.always)
    def back(self, lag: Int) -> Int:
        """The slot of the cost `lag` before the last one grown, `lag < slots`; before the first cost,
        one no cost has taken yet."""
        var slot = self.current - lag
        return slot + self.slots if slot < 0 else slot

    @inline(.always)
    def row(mut self, slot: Int, layer: Int) -> Slot:
        """A slot's front of one layer, indexed by diagonal."""
        return (
            self.buffer.unsafe_ptr()
            .unsafe_origin_cast[MutUntrackedOrigin]()
            .unsafe_offset((3 * slot + layer) * self.width + self.base)
        )

    def ready(mut self, low: Int, high: Int):
        """Makes every row cover `low ..= high`, within the diagonals a front may read."""
        if max(low, self.least) + self.base < 0 or min(high, self.most) + self.base >= self.width:
            self.grow(max(low, self.least), min(high, self.most))

    def grow(mut self, low: Int, high: Int):
        """Widens every row to cover `low ..= high`, doubling towards the side outgrown, and moves
        what they hold with it; the rest reads as unreached."""
        var size = self.width
        var first = -self.base
        var last = size - 1 - self.base
        var new_first = first
        var new_last = last
        if low < first:
            new_first = max(min(low, first - size), self.least)
        if high > last:
            new_last = min(max(high, last + size), self.most)
        var new_size = new_last - new_first + 1
        var shift = first - new_first
        var rows = 3 * self.slots
        var wider = List[Int32](length=rows * new_size, fill=UNREACHED)
        var source = self.buffer.unsafe_ptr()
        var destination = wider.unsafe_ptr()
        for row in range(rows):
            Span(unsafe_ptr=destination.unsafe_offset(row * new_size + shift), length=size).copy_from(
                Span(unsafe_ptr=source.unsafe_offset(row * size), length=size)
            )
        self.buffer = wider^
        self.width = new_size
        self.base = -new_first

    def claim(mut self, slot: Int, low: Int, high: Int):
        """Hands `slot` to a new cost over `low ..= high`, clearing what its old cost left outside."""
        var old_low = self.lows[slot]
        var old_high = self.highs[slot]
        if old_low <= old_high:
            # Only the old diagonals outside the new range: the new cost writes over the rest.
            var below_end = min(old_high, low - 1) if low <= high else old_high
            var above_start = max(old_low, high + 1) if low <= high else old_high + 1
            comptime for layer in range(3):
                var values = self.row(slot, layer)
                for diagonal in range(old_low, below_end + 1):
                    values[unsafe_offset=diagonal] = UNREACHED
                for diagonal in range(above_start, old_high + 1):
                    values[unsafe_offset=diagonal] = UNREACHED
        self.lows[slot] = low
        self.highs[slot] = high


@inline(.never)
def step[
    record: Bool
](
    mismatched: Slot,
    opening: Slot,
    first_gaps: Slot,
    second_gaps: Slot,
    aligned: Slot,
    opened_first: Slot,
    opened_second: Slot,
    flags: MutPointer[UInt8, MutUntrackedOrigin],
    first: ImmPointer[UInt8, _],
    second: ImmPointer[UInt8, _],
    low: Int,
    high: Int,
    columns: Int,
    rows: Int,
) -> Int:
    """One cost's three fronts on diagonals `low ..= high`, every pointer indexed by diagonal, and with
    `record` each diagonal's flag (see `FROM_FIRST_GAP`).

    `mismatched` is the alignment front a mismatch back, `opening` the one an opened gap back, and
    the gap fronts are their own layers an extension back. A letter of the first sequence against
    a gap comes from the diagonal below and moves one column; one of the second, from the diagonal
    above, stays in its column and moves one row. Each only where it stays inside the matrix.

    Each lane group's alignment front then slides over its matches straight away, the eight slides
    independent of each other, and the furthest anti-diagonal, `2 column - diagonal`, comes back.
    """
    comptime Lanes = SIMD[DType.int32, LANES]
    var lane_diagonals = Lanes()
    comptime for lane in range(LANES):
        lane_diagonals[lane] = Int32(lane)
    var column_limit = Lanes(Int32(columns))
    var row_limit = Lanes(Int32(rows))
    var unreached = Lanes(UNREACHED)
    var reach = Lanes(Int32.MIN // 2)
    var diagonal = low
    while diagonal <= high:
        var diagonals = lane_diagonals + Int32(diagonal)
        var same = mismatched.unsafe_offset(diagonal).unsafe_load[width=LANES]()
        var opened_below = opening.unsafe_offset(diagonal - 1).unsafe_load[width=LANES]()
        var extended_below = first_gaps.unsafe_offset(diagonal - 1).unsafe_load[width=LANES]()
        var opened_above = opening.unsafe_offset(diagonal + 1).unsafe_load[width=LANES]()
        var extended_above = second_gaps.unsafe_offset(diagonal + 1).unsafe_load[width=LANES]()
        var below = max(opened_below, extended_below)
        var above = max(opened_above, extended_above)
        var first_gap = below.lt(column_limit).select(below + 1, unreached)
        var second_gap = (above - diagonals).le(row_limit).select(above, unreached)
        var substituted = (same.lt(column_limit) & (same - diagonals).lt(row_limit)).select(same + 1, unreached)
        var gapped = max(first_gap, second_gap)
        opened_first.unsafe_offset(diagonal).unsafe_store(first_gap)
        opened_second.unsafe_offset(diagonal).unsafe_store(second_gap)
        var entry = max(substituted, gapped)
        comptime for lane in range(LANES):
            var column = Int(entry[lane])
            if column >= 0:
                entry[lane] = Int32(slide(first, second, column, diagonal + lane))
        aligned.unsafe_offset(diagonal).unsafe_store(entry)
        reach = max(reach, entry + entry - diagonals)
        comptime if record:
            # Worked out in the fronts' own lanes and narrowed once.
            var entry = substituted.ge(gapped).select(
                Lanes(0), first_gap.ge(second_gap).select(Lanes(Int32(FROM_FIRST_GAP)), Lanes(Int32(FROM_SECOND_GAP)))
            )
            var opened = opened_below.ge(extended_below).select(Lanes(Int32(FIRST_OPENED)), Lanes(0)) | opened_above.ge(
                extended_above
            ).select(Lanes(Int32(SECOND_OPENED)), Lanes(0))
            flags.unsafe_offset(diagonal).unsafe_store((entry | opened).cast[DType.uint8]())
        diagonal += LANES
    return Int(reach.reduce_max())


struct Wavefront(Movable):
    """One direction's search over two encoded sequences: the three fronts of each cost in turn, from
    cost zero at the origin, and, for the costs grown recording, what the traceback needs of each
    (see `History`).

    The backward search is the same search over both sequences reversed, whose origin is the corner.
    """

    var first: List[UInt8]
    var second: List[UInt8]
    var columns: Int
    var rows: Int
    var penalties: Penalties
    var fronts: Fronts
    var history: History
    var cost: Int
    """The last cost whose fronts are grown."""
    var furthest: Int
    """The furthest anti-diagonal, `2 column - diagonal`, any alignment front has reached."""
    var work: Int
    """Diagonals stepped so far."""
    var origin: Int
    """What the origin allows, `FREE_START` or one of the constants after it."""

    def __init__(
        out self,
        first: Span[UInt8, _],
        second: Span[UInt8, _],
        penalties: Penalties,
        origin: Int,
        record: Bool,
        reverse: Bool,
    ):
        """The fronts of cost zero over two encoded sequences, both back to front with `reverse`: the
        matches from the origin, unless its first move must open a gap, kept with `record`. Both
        sequences must hold a letter."""
        self.columns = len(first)
        self.rows = len(second)
        self.first = padded(first, FIRST_SENTINEL, reverse)
        self.second = padded(second, SECOND_SENTINEL, reverse)
        self.penalties = penalties
        self.origin = origin
        # Every source a cost reads lies at most this far back, and a slot is reused after as many.
        self.fronts = Fronts(
            max(penalties.mismatch, penalties.opening + penalties.extension) + 1, self.columns, self.rows
        )
        self.history = History()
        self.cost = 0
        self.work = 0
        self.furthest = Int.MIN // 2
        self.fronts.ready(-1 - LANES, 1 + LANES)
        var front = self.fronts.row(0, ALIGNED)
        if origin >= OPENING_FIRST_GAP:
            # Nothing at cost zero: the opening gap enters at its own cost (see `advance`). The kept
            # fronts still hold the origin at column zero, where that gap's walk back ends.
            if record:
                front[unsafe_offset=0] = 0
                _ = self.history.flags_for(0, 0)
                self.history.record(0, 0, 0, 0, front)
                front[unsafe_offset=0] = UNREACHED
            return
        self.fronts.claim(0, 0, 0)
        var start = slide(self.first.unsafe_ptr(), self.second.unsafe_ptr(), 0, 0)
        front[unsafe_offset=0] = Int32(start)
        # Inside a gap, its layer holds the origin too, which extends without a second opening.
        if origin == IN_FIRST_GAP:
            self.fronts.row(0, FIRST_GAP)[unsafe_offset=0] = 0
        elif origin == IN_SECOND_GAP:
            self.fronts.row(0, SECOND_GAP)[unsafe_offset=0] = 0
        self.fronts.reach[0] = 2 * start
        self.furthest = 2 * start
        if record:
            # The origin's flag is never read: the traceback stops at cost zero.
            _ = self.history.flags_for(0, 0)
            self.history.record(0, 0, 0, 0, front)

    def advance[record: Bool](mut self):
        """Grows the next cost's three fronts from the ring, and with `record` keeps what the traceback
        needs of them."""
        self.cost += 1
        var cost = self.cost
        var x = self.penalties.mismatch
        var o = self.penalties.opening
        var e = self.penalties.extension
        var columns = self.columns
        var rows = self.rows
        self.fronts.current = self.fronts.back(self.fronts.slots - 1)
        var slot = self.fronts.current
        # A source before the first cost reads a slot no cost has taken yet, which reads as unreached.
        var mismatch_slot = self.fronts.back(x)
        var opening_slot = self.fronts.back(o + e)
        var extension_slot = self.fronts.back(e)
        var low = Int.MAX
        var high = Int.MIN
        if cost >= x and self.fronts.lows[mismatch_slot] <= self.fronts.highs[mismatch_slot]:
            low = min(low, self.fronts.lows[mismatch_slot])
            high = max(high, self.fronts.highs[mismatch_slot])
        if cost >= o + e and self.fronts.lows[opening_slot] <= self.fronts.highs[opening_slot]:
            low = min(low, self.fronts.lows[opening_slot] - 1)
            high = max(high, self.fronts.highs[opening_slot] + 1)
        if cost >= e and self.fronts.lows[extension_slot] <= self.fronts.highs[extension_slot]:
            low = min(low, self.fronts.lows[extension_slot] - 1)
            high = max(high, self.fronts.highs[extension_slot] + 1)
        # An origin that must open a gap reaches the gap's first letter at the opening's cost.
        var opened = 0
        if self.origin >= OPENING_FIRST_GAP and cost == o + e:
            opened = 1 if self.origin == OPENING_FIRST_GAP else -1
            low = min(low, opened)
            high = max(high, opened)
        low = max(low, -rows)
        high = min(high, columns)
        if low > high:
            self.empty[record](slot)
            return
        # The step reads a lane group past `high` and a diagonal either side of the range.
        self.fronts.ready(low - 1, high + LANES + 1)
        self.fronts.claim(slot, low, high)
        var front = self.fronts.row(slot, ALIGNED)
        var first_gaps = self.fronts.row(slot, FIRST_GAP)
        var second_gaps = self.fronts.row(slot, SECOND_GAP)
        var flags = front.unsafe_bitcast[UInt8]()
        comptime if record:
            flags = self.history.flags_for(low, high)
        var reach = step[record](
            self.fronts.row(mismatch_slot, ALIGNED),
            self.fronts.row(opening_slot, ALIGNED),
            self.fronts.row(extension_slot, FIRST_GAP),
            self.fronts.row(extension_slot, SECOND_GAP),
            front,
            first_gaps,
            second_gaps,
            flags,
            self.first.unsafe_ptr(),
            self.second.unsafe_ptr(),
            low,
            high,
            columns,
            rows,
        )
        # The step wrote whole lane groups; past `high` the slot must read unreached again.
        for diagonal in range(high + 1, high + LANES):
            front[unsafe_offset=diagonal] = UNREACHED
            first_gaps[unsafe_offset=diagonal] = UNREACHED
            second_gaps[unsafe_offset=diagonal] = UNREACHED
        # The gap an origin must open enters past the step, and slides as the step's columns did.
        if opened != 0:
            var column = 1 if opened == 1 else 0
            if opened == 1:
                first_gaps[unsafe_offset=1] = 1
                comptime if record:
                    flags[unsafe_offset=1] = FROM_FIRST_GAP | FIRST_OPENED
            else:
                second_gaps[unsafe_offset=-1] = 0
                comptime if record:
                    flags[unsafe_offset=-1] = FROM_SECOND_GAP | SECOND_OPENED
            if Int(front[unsafe_offset=opened]) < column:
                column = slide(self.first.unsafe_ptr(), self.second.unsafe_ptr(), column, opened)
                front[unsafe_offset=opened] = Int32(column)
                reach = max(reach, 2 * column - opened)
        self.fronts.reach[slot] = reach
        self.furthest = max(self.furthest, reach)
        self.work += high - low + 1

        # Diagonals no layer reached at either end are dropped, as WFA trims them, so the range
        # tracks the paths alive rather than every diagonal the gap costs allow.
        @inline(.always)
        def dead(diagonal: Int) {imm front, imm first_gaps, imm second_gaps} -> Bool:
            return (
                front[unsafe_offset=diagonal] < 0
                and first_gaps[unsafe_offset=diagonal] < 0
                and second_gaps[unsafe_offset=diagonal] < 0
            )

        var kept_low = low
        var kept_high = high
        while kept_low <= kept_high and dead(kept_low):
            kept_low += 1
        while kept_high >= kept_low and dead(kept_high):
            kept_high -= 1
        if kept_low > kept_high:
            self.fronts.claim(slot, 1, 0)
            self.fronts.reach[slot] = Int.MIN // 2
            comptime if record:
                self.history.record(1, 0, low, high, front)
            return
        self.fronts.lows[slot] = kept_low
        self.fronts.highs[slot] = kept_high
        comptime if record:
            self.history.record(kept_low, kept_high, low, high, front)

    def empty[record: Bool](mut self, slot: Int):
        """Leaves the cost just begun with no diagonal reached."""
        self.fronts.claim(slot, 1, 0)
        self.fronts.reach[slot] = Int.MIN // 2
        comptime if record:
            self.history.skip()


@fieldwise_init
struct Meeting(ImplicitlyCopyable, TrivialRegisterPassable):
    """A cell an optimal path passes, in one layer, and what it costs from each end: the forward front
    of `forward_cost` reaches it on `diagonal`, and the backward front of `backward_cost` comes back
    to it, in the same layer for a gap, which then pays its opening once for both halves."""

    var cost: Int
    var layer: Int
    var diagonal: Int
    var column: Int
    var forward_cost: Int
    var backward_cost: Int

    @staticmethod
    def none() -> Meeting:
        """No meeting yet, dearer than any."""
        return Meeting(Int.MAX, ALIGNED, 0, 0, 0, 0)


def meet(mut forward: Wavefront, forward_cost: Int, mut backward: Wavefront, backward_cost: Int, mut best: Meeting):
    """Lowers `best` to where the forward fronts of one cost and the backward fronts of another
    overlap, if that costs less.

    The backward front's diagonal `target - k` mirrors the forward's `k`, and its column `c` the
    forward's `columns - c`. Two alignment fronts overlap where their columns add up to `columns` or
    more: the forward path reaches a cell the backward one comes back past, and edit costs never fall
    along a diagonal, so the cell costs at most the sum. Two gap fronts of one layer overlap the same
    way, at a cell inside the gap both may hold, and join into one gap of one opening.
    """
    var columns = forward.columns
    var rows = forward.rows
    var o = forward.penalties.opening
    var total = forward_cost + backward_cost
    if total - o >= best.cost:
        return
    var ahead = forward.fronts.back(forward.cost - forward_cost)
    var behind = backward.fronts.back(backward.cost - backward_cost)
    # A front's furthest anti-diagonal bounds all its layers': the two must add up to the matrix's.
    if forward.fronts.reach[ahead] + backward.fronts.reach[behind] < columns + rows:
        return
    var target = columns - rows
    var low = max(forward.fronts.lows[ahead], target - backward.fronts.highs[behind])
    var high = min(forward.fronts.highs[ahead], target - backward.fronts.lows[behind])
    var aligned = forward.fronts.row(ahead, ALIGNED)
    var first_gaps = forward.fronts.row(ahead, FIRST_GAP)
    var second_gaps = forward.fronts.row(ahead, SECOND_GAP)
    var back_aligned = backward.fronts.row(behind, ALIGNED)
    var back_first_gaps = backward.fronts.row(behind, FIRST_GAP)
    var back_second_gaps = backward.fronts.row(behind, SECOND_GAP)

    @inline(.always)
    def check(
        diagonal: Int,
    ) {
        mut best,
        imm aligned,
        imm first_gaps,
        imm second_gaps,
        imm back_aligned,
        imm back_first_gaps,
        imm back_second_gaps,
        imm target,
        imm columns,
        imm rows,
        imm total,
        imm o,
        imm forward_cost,
        imm backward_cost,
    }:
        var mirrored = target - diagonal
        var column = Int(aligned[unsafe_offset=diagonal])
        if total < best.cost and column + Int(back_aligned[unsafe_offset=mirrored]) >= columns:
            best = Meeting(total, ALIGNED, diagonal, column, forward_cost, backward_cost)
        if total - o >= best.cost:
            return
        # A cell inside a gap of the first sequence's letters has consumed one of them either way.
        var gap = Int(first_gaps[unsafe_offset=diagonal])
        var back_gap = Int(back_first_gaps[unsafe_offset=mirrored])
        if gap + back_gap >= columns:
            var cell = min(gap, min(columns - 1, rows + diagonal))
            if cell >= max(max(1, diagonal), columns - back_gap):
                best = Meeting(total - o, FIRST_GAP, diagonal, cell, forward_cost, backward_cost)
        # And inside a gap of the second's, one of its rows either way.
        gap = Int(second_gaps[unsafe_offset=diagonal])
        back_gap = Int(back_second_gaps[unsafe_offset=mirrored])
        if gap + back_gap >= columns:
            var cell = min(gap, rows + diagonal - 1)
            if cell >= max(diagonal + 1, columns - back_gap):
                best = Meeting(total - o, SECOND_GAP, diagonal, cell, forward_cost, backward_cost)

    # A lane group at a time, and the cells of a group that might meet one by one. An alignment front
    # reaches as far as either gap front of its cost, so where no two alignment fronts overlap, no
    # two gap fronts do.
    comptime Lanes = SIMD[DType.int32, LANES]
    var needed = Lanes(Int32(columns))
    var diagonal = low
    while diagonal + LANES - 1 <= high:
        var mirrored = target - diagonal - LANES + 1
        var meets = (
            aligned.unsafe_offset(diagonal).unsafe_load[width=LANES]()
            + back_aligned.unsafe_offset(mirrored).unsafe_load[width=LANES]().reversed()
        ).ge(needed)
        if meets.reduce_or():
            for lane in range(LANES):
                check(diagonal + lane)
        diagonal += LANES
    while diagonal <= high:
        check(diagonal)
        diagonal += 1


def bidirectional[
    record: Bool
](mut forward: Wavefront, mut backward: Wavefront, mut best: Meeting, give_up: Bool, limit: Int) -> Bool:
    """Grows the two searches a cost at a time each in turn, lowering `best` to the cheapest place an
    optimal path splits between them, until no cheaper one is left unchecked; False once a full sweep
    would be cheaper, with `give_up`, or, with `record`, once the kept fronts would pass `limit`
    entries. Either way the searches may resume from where they stopped.

    Each new front is checked against the other side's last `window` costs, which its ring holds. Every
    optimal path has a cell, between two moves or inside a gap, whose costs from the two ends differ by
    at most `window`, the dearest single move: the difference rises from `-cost` to `cost` in steps of
    at most twice that. So once both sides have passed half of `best`, plus a gap opening the two
    halves both paid and `window`, no cheaper path is left unchecked. Each side then grows to about
    half the optimum, and the two together about half the diagonals one search grows alone.
    """
    var columns = forward.columns
    var rows = forward.rows
    var letters = columns + rows
    var o = forward.penalties.opening
    var window = max(forward.penalties.mismatch, o + forward.penalties.extension)
    if forward.cost == 0 and backward.cost == 0:
        meet(forward, 0, backward, 0, best)
    var budget = columns * rows // CELLS_PER_STEP
    var next_check = CHECK_START
    while True:
        if best.cost != Int.MAX and 2 * min(forward.cost, backward.cost) >= best.cost - 1 + o + window:
            return True
        if forward.cost <= backward.cost:
            forward.advance[record]()
            var ahead = forward.cost
            for behind in range(max(0, ahead - window, backward.cost - window), backward.cost + 1):
                meet(forward, ahead, backward, behind, best)
        else:
            backward.advance[record]()
            var behind = backward.cost
            for ahead in range(max(0, behind - window, forward.cost - window), forward.cost + 1):
                meet(forward, ahead, backward, behind, best)
        comptime if record:
            if len(forward.history.aligned) + len(backward.history.aligned) > limit:
                return False
        var spent = forward.cost + backward.cost
        if give_up and spent >= next_check:
            next_check = spent + CHECK_STRIDE
            # The fronts widen with the cost, so the work grows with its square; projected from how
            # far across the two furthest anti-diagonals have come together.
            var projected = spent * letters // max(forward.furthest + backward.furthest, 1)
            var ratio = Float64(projected) / Float64(spent)
            var work = forward.work + backward.work
            if Float64(work) * (ratio * ratio - 1.0) > Float64(budget - work):
                return False


def wavefront_score(
    first: List[UInt8], second: List[UInt8], penalties: Penalties, give_up: Bool = True
) -> Optional[Int]:
    """The optimal global score of two encoded sequences, or None once a full sweep would be cheaper.

    With `give_up` off, the search runs to the end whatever it costs.
    """
    var letters = len(first) + len(second)
    if len(first) == 0 or len(second) == 0:
        var cost = 0 if letters == 0 else penalties.opening + penalties.extension * letters
        return penalties.score(cost, letters)
    var forward = Wavefront(Span(first), Span(second), penalties, FREE_START, False, False)
    var backward = Wavefront(Span(first), Span(second), penalties, FREE_START, False, True)
    var best = Meeting.none()
    if not bidirectional[False](forward, backward, best, give_up, Int.MAX):
        return None
    return penalties.score(best.cost, letters)


def trace(
    history: History,
    penalties: Penalties,
    layer: Int,
    cost: Int,
    diagonal: Int,
    column: Int,
    mut moves: List[UInt8],
):
    """An optimal path from the origin to a cell in `layer` that the front of `cost` in that layer
    reaches on `diagonal`, as moves right to left, through every cost's kept fronts.

    The cell need not be the front's furthest, as where the two searches met: every cell short of a
    front on its diagonal costs no more than the front, edit costs never falling along a diagonal. So
    the walk holds a layer, a cost and a cell its front of that cost reaches, and each step follows
    the source that won the front there, which reaches at least as far as any other. In the
    alignment layer that is a substitution from the front `x` back, or a gap layer at the same cost,
    entered after the matches back to its column when the cell lies past it. A gap cell came by a
    letter of its sequence from the diagonal beside it, opened from the alignment front `o + e` back
    or extended from its own layer `e` back. A cell on the first row or column has one path left, a
    gap along it.
    """
    var x = penalties.mismatch
    var o = penalties.opening
    var e = penalties.extension
    var spent = cost
    var current = layer
    var k = diagonal
    var c = column
    while True:
        if current == ALIGNED:
            var row = c - k
            if c == 0:
                for _ in range(row):
                    moves.append(UInt8(SECOND_GAP))
                return
            if row == 0:
                for _ in range(c):
                    moves.append(UInt8(FIRST_GAP))
                return
            if spent == 0:
                for _ in range(c):
                    moves.append(UInt8(ALIGNED))
                return
            var source = history.flag(spent, k) & 3
            # The column the winning source brought the front to, before its matches.
            var entry: Int
            if source == 0:
                entry = history.column(spent - x, k) + 1
            else:
                # A gap layer's column is where its run opened, and one letter on per extension.
                var gap_cost = spent
                var gap_diagonal = k
                var letters = 0
                # A run that reaches cost zero continues a gap the origin was inside, at column zero.
                if source == FROM_FIRST_GAP:
                    while gap_cost > 0 and history.flag(gap_cost, gap_diagonal) & FIRST_OPENED == 0:
                        gap_cost -= e
                        gap_diagonal -= 1
                        letters += 1
                    if gap_cost == 0:
                        entry = letters
                    else:
                        entry = history.column(gap_cost - o - e, gap_diagonal - 1) + letters + 1
                else:
                    while gap_cost > 0 and history.flag(gap_cost, gap_diagonal) & SECOND_OPENED == 0:
                        gap_cost -= e
                        gap_diagonal += 1
                    entry = 0 if gap_cost == 0 else history.column(gap_cost - o - e, gap_diagonal + 1)
            if entry <= c:
                for _ in range(c - entry):
                    moves.append(UInt8(ALIGNED))
                c = entry
            if source == 0:
                moves.append(UInt8(ALIGNED))
                c -= 1
                spent -= x
            elif source == FROM_FIRST_GAP:
                current = FIRST_GAP
            else:
                current = SECOND_GAP
        elif current == FIRST_GAP:
            var opened = history.flag(spent, k) & FIRST_OPENED != 0
            moves.append(UInt8(FIRST_GAP))
            c -= 1
            k -= 1
            if c == 0:
                for _ in range(-k):
                    moves.append(UInt8(SECOND_GAP))
                return
            if opened:
                spent -= o + e
                current = ALIGNED
            else:
                spent -= e
        else:
            var opened = history.flag(spent, k) & SECOND_OPENED != 0
            moves.append(UInt8(SECOND_GAP))
            k += 1
            if c - k == 0:
                for _ in range(c):
                    moves.append(UInt8(FIRST_GAP))
                return
            if opened:
                spent -= o + e
                current = ALIGNED
            else:
                spent -= e


def solve(
    first: Span[UInt8, _],
    second: Span[UInt8, _],
    penalties: Penalties,
    start: Int,
    finish: Int,
    limit: Int,
    mut moves: List[UInt8],
    keep: Bool = True,
) -> Int:
    """Appends an optimal path's moves right to left and returns its cost, from an origin as `start`
    allows and to a corner the backward search's origin `finish` allows (see `FREE_START`).

    With `keep`, both searches keep every cost's fronts while they stay within `limit` entries, and
    the path is traced from where they met: back to the origin through the forward fronts, and to the
    corner through the backward ones. A pair too large is split instead where an optimal path
    crosses, which the two searches find keeping only their rings, as BiWFA does; a crossing inside a
    gap leaves the piece before it to end in that gap and the piece after it to begin there, the
    opening paid once. A piece is about a quarter of the pair, its two searches about half of the
    diagonals the pair's search grew from its end, so one that cannot fit skips keeping at once.
    """
    var columns = len(first)
    var rows = len(second)
    if columns == 0 or rows == 0:
        for _ in range(columns):
            moves.append(UInt8(FIRST_GAP))
        for _ in range(rows):
            moves.append(UInt8(SECOND_GAP))
        var letters = columns + rows
        var continued = (start == IN_FIRST_GAP and rows == 0) or (start == IN_SECOND_GAP and columns == 0)
        return 0 if letters == 0 else penalties.extension * letters + (0 if continued else penalties.opening)
    var forward = Wavefront(first, second, penalties, start, keep, False)
    var backward = Wavefront(first, second, penalties, finish, keep, True)
    var best = Meeting.none()
    if keep:
        if bidirectional[True](forward, backward, best, False, limit):
            # The backward walk runs from the meeting to the corner, left to right as the forward path goes.
            var behind = List[UInt8](capacity=columns + rows)
            trace(
                backward.history,
                penalties,
                best.layer,
                best.backward_cost,
                columns - rows - best.diagonal,
                columns - best.column,
                behind,
            )
            for index in range(len(behind) - 1, -1, -1):
                moves.append(behind[index])
            trace(forward.history, penalties, best.layer, best.forward_cost, best.diagonal, best.column, moves)
            return best.cost
        # Too large to keep: the searches go on from where they stopped keeping only their rings.
        forward.history = History()
        backward.history = History()
    _ = bidirectional[False](forward, backward, best, False, Int.MAX)
    var column = best.column
    var row = column - best.diagonal
    # A crossing at either end splits nothing: such a pair costs too little for its fronts not to fit.
    if (column == 0 and row == 0) or (column == columns and row == rows):
        return solve(first, second, penalties, start, finish, Int.MAX, moves)
    var before_finish = FREE_START
    var after_start = FREE_START
    if best.layer == FIRST_GAP:
        before_finish = OPENING_FIRST_GAP
        after_start = IN_FIRST_GAP
    elif best.layer == SECOND_GAP:
        before_finish = OPENING_SECOND_GAP
        after_start = IN_SECOND_GAP
    _ = solve(first[column:], second[row:], penalties, after_start, finish, limit, moves, backward.work // 2 <= limit)
    _ = solve(first[:column], second[:row], penalties, start, before_finish, limit, moves, forward.work // 2 <= limit)
    return best.cost


def wavefront_align(
    first: List[UInt8], second: List[UInt8], penalties: Penalties, alphabet: String, limit: Int = HISTORY_LIMIT
) -> Tuple[Int, String, String]:
    """The optimal global score of two encoded sequences and the gapped rows of an alignment that earns
    it, keeping at most `limit` entries of fronts at once (see `solve`)."""
    var letters = len(first) + len(second)
    var moves = List[UInt8](capacity=letters)
    var cost = solve(Span(first), Span(second), penalties, FREE_START, FREE_START, limit, moves)
    comptime GAP = UInt8(ord("-"))
    var symbols = alphabet.as_bytes()
    var top = List[UInt8](capacity=len(moves))
    var bottom = List[UInt8](capacity=len(moves))
    var column = 0
    var row = 0
    for index in range(len(moves) - 1, -1, -1):
        var move = Int(moves[index])
        if move == ALIGNED:
            top.append(symbols[Int(first[column])])
            bottom.append(symbols[Int(second[row])])
            column += 1
            row += 1
        elif move == FIRST_GAP:
            top.append(symbols[Int(first[column])])
            bottom.append(GAP)
            column += 1
        else:
            top.append(GAP)
            bottom.append(symbols[Int(second[row])])
            row += 1
    return (penalties.score(cost, letters), String(unsafe_from_utf8=top), String(unsafe_from_utf8=bottom))


@fieldwise_init
struct AffineCigar(Copyable, Movable, Writable):
    """The least cost of a global alignment of two sequences under gap-affine costs and an optimal
    alignment as a CIGAR string, written as `EditCigar`'s is."""

    var cost: Int
    var cigar: String


def affine_cigar(
    first: String, second: String, mismatch: Int, opening: Int, extension: Int, extended: Bool = True
) raises AlignmentError -> AffineCigar:
    """The least cost of a global alignment of two sequences under gap-affine costs as WFA counts them,
    a substitution `mismatch` and a gap of `k` letters `opening + k extension`, and an optimal
    alignment as a CIGAR string, `=` and `X`, or with `extended` false `M` for both (see `EditCigar`).

    Every byte is a symbol matching only itself, so DNA in either case, or any other text, needs no
    alphabet. The two-ended wavefront finds it (see the module): its work grows with the square of
    the cost rather than with the matrix, and its memory stays bounded.
    """
    if mismatch <= 0 or extension <= 0 or opening < 0:
        raise AlignmentError(
            ErrorKind.INVALID_SCORING,
            String("costs ", mismatch, ", ", opening, ", ", extension, ": a mismatch and an extension must cost"),
        )
    var scale = greatest_common_divisor(greatest_common_divisor(mismatch, extension), opening)
    var penalties = Penalties(mismatch // scale, opening // scale, extension // scale, scale, 0)
    var columns = first.byte_length()
    var rows = second.byte_length()
    # UTF-8 never holds the sentinels' bytes, so the text's own bytes serve as codes.
    var moves = List[UInt8](capacity=columns + rows)
    var cost = solve(first.as_bytes(), second.as_bytes(), penalties, FREE_START, FREE_START, HISTORY_LIMIT, moves)
    # The CIGAR's room is bounded by the edits: every gapped letter, and a substitution per mismatch cost.
    var gapped = 0
    for move in moves:
        if move != UInt8(ALIGNED):
            gapped += 1
    var path = EditPath(moves^, List[UInt8](), columns, rows, gapped + cost // penalties.mismatch)
    return AffineCigar(cost * scale, cigar_string(first, second, path, extended))
