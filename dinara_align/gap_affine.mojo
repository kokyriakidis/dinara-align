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
alignment front slides over matches for free. The first cost whose front reaches the corner is the
optimum. The work grows with the square of the cost, not with the matrix, so a close pair costs a
small share of a full sweep; a pair whose projected work would pass a share of the sweep is handed
back, and the dynamic programming answers it.
"""

from std.bit import count_trailing_zeros

comptime Front = List[Int32]

comptime Slot = MutPointer[Int32, MutUntrackedOrigin]
"""A front read or written by `step`: a source and its destination may share a ring slot's list."""

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

comptime CELLS_PER_STEP = 4
"""Cells of the vectorized full sweep (see `vector_score`) one diagonal step of the three fronts
costs about as much as: about 1.5 against 0.4 ns."""


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


def padded(codes: List[UInt8], sentinel: UInt8) -> List[UInt8]:
    """The codes with `PADDING` sentinels after them."""
    var out = List[UInt8](capacity=len(codes) + PADDING)
    out.extend(codes.copy())
    for _ in range(PADDING):
        out.append(sentinel)
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


struct Fronts(Movable):
    """The last `slots` costs' three fronts, in rings indexed by cost, each over the diagonals reached.

    Diagonal `k` holds the cells whose column minus row is `k`, from `-rows` to `columns`, stored
    at `k + base`. The buffers cover only the diagonals the fronts have reached, about `2 s / e`
    of them at cost `s`, and double towards whichever side they outgrow, so the memory grows with
    the score rather than the sequences. A slot holds its cost's values on `lows[slot] ..= highs[slot]` and reads as
    unreached everywhere else, so a front reads its sources' neighbours unchecked: the diagonals
    the rings have ever been asked for are kept initialized, and a slot reused for a later cost has
    the old cost's diagonals outside the new range cleared.
    """

    var slots: Int
    var base: Int
    var aligned: List[Front]
    var opened_first: List[Front]
    """Fronts ending in a letter of the first sequence against a gap."""
    var opened_second: List[Front]
    """Fronts ending in a letter of the second sequence against a gap."""
    var lows: List[Int]
    var highs: List[Int]
    var ready_low: Int
    """The diagonals every slot holds a value for, unreached or not."""
    var ready_high: Int
    var least: Int
    """The diagonals a front may ever read: the matrix's and two either side, and a lane group past."""
    var most: Int

    def __init__(out self, slots: Int, columns: Int, rows: Int):
        self.slots = slots
        self.least = -rows - 2
        self.most = columns + LANES + 2
        var width = 4 * LANES
        self.base = width // 2
        self.aligned = List[Front](capacity=slots)
        self.opened_first = List[Front](capacity=slots)
        self.opened_second = List[Front](capacity=slots)
        for _ in range(slots):
            var empty = Front(capacity=width)
            empty.resize(unsafe_uninit_length=width)
            self.aligned.append(empty.copy())
            self.opened_first.append(empty.copy())
            self.opened_second.append(empty^)
        self.lows = List[Int](length=slots, fill=1)
        self.highs = List[Int](length=slots, fill=0)
        self.ready_low = 1
        self.ready_high = 0

    def ready(mut self, wanted_low: Int, wanted_high: Int):
        """Makes every slot hold unreached on `low ..= high` wherever it held nothing yet."""
        var low = max(wanted_low, self.least)
        var high = min(wanted_high, self.most)
        if low + self.base < 0 or high + self.base >= len(self.aligned[0]):
            self.grow(low, high)
        if self.ready_low > self.ready_high:
            self.ready_low = low
            self.ready_high = low - 1
        for diagonal in range(low, self.ready_low):
            self.clear(diagonal)
        for diagonal in range(self.ready_high + 1, high + 1):
            self.clear(diagonal)
        self.ready_low = min(self.ready_low, low)
        self.ready_high = max(self.ready_high, high)

    def grow(mut self, low: Int, high: Int):
        """Widens every buffer to cover `low ..= high`, doubling towards the side outgrown, and
        moves what they hold with it."""
        var size = len(self.aligned[0])
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

        @inline(.always)
        def moved(old: Front) {imm new_size, imm shift, imm size} -> Front:
            var wider = Front(capacity=new_size)
            wider.resize(unsafe_uninit_length=new_size)
            for index in range(size):
                wider[index + shift] = old[index]
            return wider^

        for slot in range(self.slots):
            self.aligned[slot] = moved(self.aligned[slot])
            self.opened_first[slot] = moved(self.opened_first[slot])
            self.opened_second[slot] = moved(self.opened_second[slot])
        self.base = -new_first

    @inline(.always)
    def clear(mut self, diagonal: Int):
        for slot in range(self.slots):
            self.aligned[slot][diagonal + self.base] = UNREACHED
            self.opened_first[slot][diagonal + self.base] = UNREACHED
            self.opened_second[slot][diagonal + self.base] = UNREACHED

    @inline(.always)
    def unset(mut self, slot: Int, diagonal: Int):
        self.aligned[slot][diagonal + self.base] = UNREACHED
        self.opened_first[slot][diagonal + self.base] = UNREACHED
        self.opened_second[slot][diagonal + self.base] = UNREACHED

    def claim(mut self, slot: Int, low: Int, high: Int):
        """Hands `slot` to a new cost over `low ..= high`, clearing what its old cost left outside."""
        var old_low = self.lows[slot]
        var old_high = self.highs[slot]
        if old_low <= old_high:
            # Only the old diagonals outside the new range: the new cost writes over the rest.
            var below_end = min(old_high, low - 1) if low <= high else old_high
            for diagonal in range(old_low, below_end + 1):
                self.unset(slot, diagonal)
            if low <= high:
                for diagonal in range(max(old_low, high + 1), old_high + 1):
                    self.unset(slot, diagonal)
        self.lows[slot] = low
        self.highs[slot] = high


@inline(.never)
def step(
    mismatched: Slot,
    opening: Slot,
    first_gaps: Slot,
    second_gaps: Slot,
    aligned: Slot,
    opened_first: Slot,
    opened_second: Slot,
    low: Int,
    high: Int,
    columns: Int,
    rows: Int,
):
    """One cost's three fronts on diagonals `low ..= high`, every pointer indexed by diagonal.

    `mismatched` is the alignment front a mismatch back, `opening` the one an opened gap back, and
    the gap fronts are their own layers an extension back. A letter of the first sequence against
    a gap comes from the diagonal below and moves one column; one of the second, from the diagonal
    above, stays in its column and moves one row. Each only where it stays inside the matrix.
    """
    comptime Lanes = SIMD[DType.int32, LANES]
    var lane_diagonals = Lanes()
    comptime for lane in range(LANES):
        lane_diagonals[lane] = Int32(lane)
    var column_limit = Lanes(Int32(columns))
    var row_limit = Lanes(Int32(rows))
    var unreached = Lanes(UNREACHED)
    var diagonal = low
    while diagonal <= high:
        var diagonals = lane_diagonals + Int32(diagonal)
        var same = mismatched.unsafe_offset(diagonal).unsafe_load[width=LANES]()
        var below = max(
            opening.unsafe_offset(diagonal - 1).unsafe_load[width=LANES](),
            first_gaps.unsafe_offset(diagonal - 1).unsafe_load[width=LANES](),
        )
        var above = max(
            opening.unsafe_offset(diagonal + 1).unsafe_load[width=LANES](),
            second_gaps.unsafe_offset(diagonal + 1).unsafe_load[width=LANES](),
        )
        var first_gap = below.lt(column_limit).select(below + 1, unreached)
        var second_gap = (above - diagonals).le(row_limit).select(above, unreached)
        var substituted = (same.lt(column_limit) & (same - diagonals).lt(row_limit)).select(same + 1, unreached)
        opened_first.unsafe_offset(diagonal).unsafe_store(first_gap)
        opened_second.unsafe_offset(diagonal).unsafe_store(second_gap)
        aligned.unsafe_offset(diagonal).unsafe_store(max(substituted, max(first_gap, second_gap)))
        diagonal += LANES


def wavefront_score(
    first: List[UInt8], second: List[UInt8], penalties: Penalties, give_up: Bool = True
) -> Optional[Int]:
    """The optimal global score of two encoded sequences, or None once a full sweep would be cheaper.

    With `give_up` off, the search runs to the end whatever it costs.
    """
    var columns = len(first)
    var rows = len(second)
    var letters = columns + rows
    var x = penalties.mismatch
    var o = penalties.opening
    var e = penalties.extension
    if columns == 0 or rows == 0:
        if letters == 0:
            return 0
        return penalties.score(o + e * letters, letters)
    var first_codes = padded(first, FIRST_SENTINEL)
    var second_codes = padded(second, SECOND_SENTINEL)
    var first_pointer = first_codes.unsafe_ptr()
    var second_pointer = second_codes.unsafe_ptr()
    var target = columns - rows
    # Every source a cost reads lies at most this far back, and a slot is reused after as many.
    var slots = max(x, o + e) + 1
    var fronts = Fronts(slots, columns, rows)
    fronts.ready(-1 - LANES, 1 + LANES)
    var base = fronts.base
    fronts.claim(0, 0, 0)
    var start = slide(first_pointer, second_pointer, 0, 0)
    fronts.aligned[0][base] = Int32(start)
    if target == 0 and start >= columns:
        return penalties.score(0, letters)
    var budget = columns * rows // CELLS_PER_STEP
    var work = 0
    var furthest = 2 * start
    var next_check = CHECK_START
    var cost = 0
    while True:
        cost += 1
        var slot = cost % slots

        @inline(.always)
        def source(lag: Int) {imm cost, imm slots} -> Int:
            """The slot of the cost `lag` back, or -1 before the first."""
            return (cost - lag) % slots if cost >= lag else -1

        var from_mismatch = source(x)
        var from_opening = source(o + e)
        var from_extension = source(e)
        var low = Int.MAX
        var high = Int.MIN
        if from_mismatch >= 0 and fronts.lows[from_mismatch] <= fronts.highs[from_mismatch]:
            low = min(low, fronts.lows[from_mismatch])
            high = max(high, fronts.highs[from_mismatch])
        if from_opening >= 0 and fronts.lows[from_opening] <= fronts.highs[from_opening]:
            low = min(low, fronts.lows[from_opening] - 1)
            high = max(high, fronts.highs[from_opening] + 1)
        if from_extension >= 0 and fronts.lows[from_extension] <= fronts.highs[from_extension]:
            low = min(low, fronts.lows[from_extension] - 1)
            high = max(high, fronts.highs[from_extension] + 1)
        low = max(low, -rows)
        high = min(high, columns)
        if low > high:
            fronts.claim(slot, 1, 0)
            continue
        # The step reads a lane group past `high` and a diagonal either side of the range.
        fronts.ready(low - 1, high + LANES + 1)
        # Growing moved the buffers, and the diagonals with them.
        base = fronts.base
        fronts.claim(slot, low, high)
        # A source before the first cost reads the slot it will take, which no cost has written yet.
        var mismatch_slot = from_mismatch if from_mismatch >= 0 else (cost + slots - x) % slots
        var opening_slot = from_opening if from_opening >= 0 else (cost + slots - o - e) % slots
        var extension_slot = from_extension if from_extension >= 0 else (cost + slots - e) % slots

        @inline(.always)
        def at(mut front: Front) {imm base} -> Slot:
            return front.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]().unsafe_offset(base)

        var front = at(fronts.aligned[slot])
        step(
            at(fronts.aligned[mismatch_slot]),
            at(fronts.aligned[opening_slot]),
            at(fronts.opened_first[extension_slot]),
            at(fronts.opened_second[extension_slot]),
            front,
            at(fronts.opened_first[slot]),
            at(fronts.opened_second[slot]),
            low,
            high,
            columns,
            rows,
        )
        # The step wrote whole lane groups; past `high` the slot must read unreached again.
        for diagonal in range(high + 1, high + LANES):
            fronts.aligned[slot][diagonal + base] = UNREACHED
            fronts.opened_first[slot][diagonal + base] = UNREACHED
            fronts.opened_second[slot][diagonal + base] = UNREACHED
        # Then each alignment front slides over its matches.
        for diagonal in range(low, high + 1):
            var column = Int(front[unsafe_offset=diagonal])
            if column >= 0:
                column = slide(first_pointer, second_pointer, column, diagonal)
                front[unsafe_offset=diagonal] = Int32(column)
                furthest = max(furthest, 2 * column - diagonal)
        if target >= low and target <= high and Int(front[unsafe_offset=target]) >= columns:
            return penalties.score(cost, letters)
        work += high - low + 1
        # Diagonals no layer reached at either end are dropped, as WFA trims them, so the range
        # tracks the paths alive rather than every diagonal the gap costs allow.
        var first_gaps = at(fronts.opened_first[slot])
        var second_gaps = at(fronts.opened_second[slot])

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
            fronts.claim(slot, 1, 0)
        elif kept_low != low or kept_high != high:
            fronts.lows[slot] = kept_low
            fronts.highs[slot] = kept_high
        if give_up and cost >= next_check:
            next_check = cost + CHECK_STRIDE
            # The fronts widen with the cost, so the work grows with its square; projected from how
            # far across the furthest diagonal has come.
            var projected = cost * letters // max(furthest, 1)
            var ratio = Float64(projected) / Float64(cost)
            if Float64(work) * (ratio * ratio - 1.0) > Float64(budget - work):
                return None
