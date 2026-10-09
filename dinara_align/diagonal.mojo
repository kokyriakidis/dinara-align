# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Ported from `pa-bitpacking` in A*PA (https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner,
# commit bf2e14e), by Ragnar Groot Koerkamp and Pesho Ivanov, itself translated from Edlib.
"""
Diagonal transition, as WFA runs it, before any band: from the start alone, or from both ends at
once, for as long as it stays cheaper than a band would be. Near-identical pairs, which the post on
A*PA2 names as its weak spot, finish here, traceback included; any other pair leaves with a projection
of its distance, which becomes the band's first bound (see `band.band_doubling`).
"""

from std.bit import count_trailing_zeros
from std.math import sqrt
from std.sys import simd_width_of

from .common import UNREACHED
from .cigar import reversed_into
from .slides import GATHERED_SLIDES, gathered_slides, slide
from .bit_parallel import (
    advance,
    append_diagonals,
    CODE_PADDING,
    DIAGONAL,
    FIRST_SENTINEL,
    LEFT,
    Profile,
    SECOND_SENTINEL,
    UP,
)


comptime FRONT_PADDING = 2
"""Unreached diagonals stored either side of a front, so the next front reads its neighbours unchecked."""


comptime FRONT_LANES = 8
"""Diagonals a front step computes at once."""


comptime PROBE_STRIDE = 4
"""Scores between the diagonal transition's checks of its projection."""


comptime PROBE_START = 8
"""The score from which the diagonal transition judges whether to go on, so the projection has a few
edits to go on."""


comptime STEP_TENTHS_DISTANCE = 60
"""Diagonal steps, in tenths, a distance's diagonal transition may still take per column before a band
would be cheaper, at no distance: a step costs about 1.3 ns, a band column about 6 ns."""


comptime STEP_TENTHS_ALIGNMENT = 50
"""`STEP_TENTHS_DISTANCE` when an alignment is wanted: the diagonal transition then keeps every front
and the band records and retraces, which roughly cancel, so the budget sits a little under the
distance's."""


comptime TWO_ENDED_PERCENT = 57
"""The two-ended search's cost per square edit of distance, in percent of one front's step."""


comptime ALIGNED_TWO_ENDED_PERCENT = 73
"""The two-ended alignment's cost per square edit of distance, in percent of one front's step: the
search's own `TWO_ENDED_PERCENT` and, since a tie rule needs the path traced from the corner, the
forward fronts grown on to the distance (see `grow_to`)."""


comptime TWO_ENDED_SETUP = 2600
"""What starting the two-ended alignment costs over one front, in one front's steps, besides a tenth
of a step a base (see `two_ended_setup`): its second front, its histories, and growing the forward
fronts on. An alignment's search runs one front while what it has left costs less than that. Fitted
on 1,300 pairs from 1 to 30 kbp: one front wins every pair below about 100 edits, both ends nearly
every pair above."""


@always_inline
def two_ended_setup(columns: Int, rows: Int) -> Int:
    """What starting the two-ended alignment costs a pair over one front, in one front's steps."""
    return TWO_ENDED_SETUP + (columns + rows) // 10


comptime NO_DIAGONAL = Int.MIN
"""No diagonal: the fronts have not met."""


comptime GIVE_UP_SHARE = 24
"""A diagonal transition gives up on its projection only once it has spent at least one part in this
many of what the cheapest band could cost, the band whose bound is the score already searched, since
the distance is at least that.

The projection from a front's first edits assumes them spread along the pair. Real reads and genomes
gather theirs at the ends, a primer or a poly-A tail apart, and a projection from those alone can run
hundreds of times over a distance of a few dozen; giving up on it hands a pair the search would finish
in microseconds to a band over the whole matrix. Spending this share first bounds what a search that
was going to give up anyway wastes to it, of a band's cost."""


comptime PROBE_SPREAD = 1.5
"""Spreads of a projection's noise the search allows for before it gives up (see `noisy_budget`)."""


comptime PROJECTION_ONLY = -(1 << 40)
"""A budget no search fits: the diagonal transition stops at its first check, with its projection."""


comptime EDITS_PER_STEP = 83
"""A band's column grows about 1.3 ns, a step, per this many edits of projected distance, its band
growing taller with the distance."""


comptime PROBE_MARGIN = 16
"""How far past the projected distance the band's first bound reaches, on top of an eighth of it."""


comptime SHORT_COLUMNS = 4096
"""Pairs up to this many columns aim the band's first bound further past the projection.

On a short pair the projection from a handful of edits ran as much as 1.7 times under the distance,
and a first round that fails there costs nearly a whole round more, where aiming wide costs a band
some rows taller; a long pair's projection landed within about a sixth either way."""


comptime SHORT_BAND_COLUMNS = 2048
"""Pairs up to this many columns whose band costs more than its columns: aimed wide (see `SHORT_AIM`),
its first bound lifts it over the whole matrix, which it sweeps and then traces back through
full-height tiles. On ONT reads of about a kilobase that came to some 40 µs a read, several times what
the columns count, so a diagonal transition the per-column budget gave up on still beat it."""


comptime SHORT_BAND_SETUP = 30_000
"""What a band costs a pair of up to `SHORT_BAND_COLUMNS` columns beyond its columns, in diagonal
steps: ONT reads of about a kilobase aligned fastest from 20,000 to 40,000, a quarter faster than with
none, and longer pairs, whose bands stay narrow, are left out."""


comptime PROBE_CEILING = 2048
"""The highest score the diagonal transition reaches, which bounds its memory to a few megabytes."""


comptime TWO_ENDED_LIMIT = PROBE_CEILING // 2 + 1
"""Scores each direction of the two-ended search keeps room for, the two together the ceiling."""


comptime FRONT_CENTER = TWO_ENDED_LIMIT + FRONT_PADDING
"""Where diagonal zero sits in a two-ended front's slot: room on the left for its lowest diagonal and
padding, on the right also for a step's last vector, which reads `FRONT_LANES` past the front."""


comptime FRONT_WIDTH = 2 * FRONT_CENTER + FRONT_LANES + 1
"""A two-ended front's slot."""


struct DiagonalFronts(Movable):
    """Every score's wavefront: for score `s`, the furthest column each diagonal reaches.

    Diagonal `k` holds the cells whose column minus row is `k`. Score `s` keeps diagonals
    `lows[s] ..= highs[s]`, with `FRONT_PADDING` unreached diagonals either side, from position
    `starts[s]` of `offsets`; all of them stay, since the traceback reads them back.
    """

    var offsets: List[Int32]
    var starts: List[Int]
    var lows: List[Int]
    var highs: List[Int]

    def __init__(out self, reserve: Bool = True):
        """Room for a typical search up front, unless `reserve` is false and none is kept."""
        var offsets = 4096 if reserve else 0
        var scores = 64 if reserve else 0
        self.offsets = List[Int32](capacity=offsets)
        self.starts = List[Int](capacity=scores)
        self.lows = List[Int](capacity=scores)
        self.highs = List[Int](capacity=scores)

    @inline(.always)
    def at(self, score: Int, diagonal: Int) -> Int:
        """The furthest column of `diagonal` at `score`, or -1 when no path of that score reaches it."""
        var low = Int(self.lows.unsafe_ptr()[unsafe_offset=score])
        if diagonal < low or diagonal > Int(self.highs.unsafe_ptr()[unsafe_offset=score]):
            return -1
        var start = Int(self.starts.unsafe_ptr()[unsafe_offset=score])
        var column = Int(self.offsets.unsafe_ptr()[unsafe_offset=start + FRONT_PADDING + diagonal - low])
        return column if column >= 0 else -1


@fieldwise_init
struct Probe(ImplicitlyCopyable, TrivialRegisterPassable):
    """What the diagonal transition learned: the distance, or a projection of it and a floor under it."""

    var distance: Int
    """The edit distance, or -1 when the search stopped first."""
    var estimate: Int
    """Where the distance was heading when the search stopped: its score scaled by the progress made."""
    var floor: Int
    """Every score up to this was searched, so the distance exceeds it."""


@inline(.always)
def slide_forward(first: ImmPointer[UInt8, _], second: ImmPointer[UInt8, _], column: Int, row: Int) -> Int:
    """How far matches carry `(column, row)` along its diagonal, sixteen bases at a time.

    The sequences' sentinels differ from every base and from each other, so the run stops at the
    matrix's edge by itself (see `CODE_PADDING`).
    """
    var at = column
    var down = row
    while True:
        var low = (
            first.unsafe_offset(at).unsafe_bitcast[UInt64]().unsafe_load()
            ^ second.unsafe_offset(down).unsafe_bitcast[UInt64]().unsafe_load()
        )
        if low != 0:
            return at + Int(count_trailing_zeros(low)) // 8
        var high = (
            first.unsafe_offset(at + 8).unsafe_bitcast[UInt64]().unsafe_load()
            ^ second.unsafe_offset(down + 8).unsafe_bitcast[UInt64]().unsafe_load()
        )
        if high != 0:
            return at + 8 + Int(count_trailing_zeros(high)) // 8
        at += 16
        down += 16


@inline(.always)
def extend(
    first: ImmPointer[UInt8, _], second: ImmPointer[UInt8, _], column: Int, row: Int, columns: Int, rows: Int
) -> Int:
    """How far matches carry `(column, row)` along its diagonal, eight bases at a time while both have eight."""
    var at = column
    var down = row
    while at + 8 <= columns and down + 8 <= rows:
        var mismatches = (
            first.unsafe_offset(at).unsafe_bitcast[UInt64]().unsafe_load()
            ^ second.unsafe_offset(down).unsafe_bitcast[UInt64]().unsafe_load()
        )
        if mismatches != 0:
            return at + Int(count_trailing_zeros(mismatches)) // 8
        at += 8
        down += 8
    while at < columns and down < rows and first[unsafe_offset=at] == second[unsafe_offset=down]:
        at += 1
        down += 1
    return at


@inline(.always)
def best_source(fronts: DiagonalFronts, score: Int, diagonal: Int, columns: Int, rows: Int) -> Tuple[Int, UInt8]:
    """The furthest column a path of `score` reaches on `diagonal` before sliding over matches, and its last move.

    One more edit after the furthest point of score `score - 1` on this diagonal or a neighbour,
    each only where that edit stays inside the matrix; -1 when none does.
    """
    var best = -1
    var move = DIAGONAL
    var same = fronts.at(score - 1, diagonal)
    if same >= 0 and same < columns and same - diagonal < rows:
        best = same + 1
    # A base of the first sequence against a gap, from the diagonal below.
    var left = fronts.at(score - 1, diagonal - 1)
    if left >= 0 and left < columns and left + 1 > best:
        best = left + 1
        move = LEFT
    # A base of the second sequence against a gap, from the diagonal above.
    var up = fronts.at(score - 1, diagonal + 1)
    if up >= 0 and up - diagonal - 1 < rows and up > best:
        best = up
        move = UP
    return (best, move)


@inline(.always)
def noisy_budget(budget: Int, edits: Int) -> Int:
    """`budget` widened for a projection from `edits` edits, about one in `sqrt(edits)` off.

    A projection from a few edits can land far above the distance, and giving up on it then hands a
    pair the band would align more slowly; so an alignment's search goes on until the projected work
    passes the budget by more than the projection's noise, `PROBE_SPREAD` of its spreads. A distance
    has a cheaper band to fall back on, and gives up on the plain budget.
    """
    if budget < 0:
        return budget
    return Int(Float64(budget) * (1.0 + PROBE_SPREAD / sqrt(Float64(max(edits, 1)))))


@inline(.always)
def step_budget(columns: Int, step_tenths: Int, estimate: Int) -> Int:
    """How many diagonal steps cost what a band over `columns` columns would, at a projected distance.

    A band's column costs a fixed part plus a part growing with the distance, its band taller; in
    steps that is `step_tenths / 10 + estimate / EDITS_PER_STEP` a column. A short pair's band costs
    `SHORT_BAND_SETUP` more (see `SHORT_BAND_COLUMNS`).
    """
    var setup = SHORT_BAND_SETUP if columns <= SHORT_BAND_COLUMNS else 0
    return setup + columns * step_tenths // 10 + columns * estimate // EDITS_PER_STEP


@inline(.never)
def step_front[
    measure: Bool
](
    previous: ImmPointer[Int32, _],
    current: MutPointer[Int32, _],
    low: Int,
    high: Int,
    columns: Int,
    rows: Int,
    first: ImmPointer[UInt8, _],
    second: ImmPointer[UInt8, _],
) -> Int:
    """One score's front from the last: every diagonal `low ..= high` of `current`, both indexed by diagonal.

    A function of its own, so its handful of values stay in registers. With `measure`, returns
    the furthest anti-diagonal, column plus row, any diagonal reached; otherwise zero.
    """
    comptime Lanes = SIMD[DType.int32, FRONT_LANES]
    # Every diagonal's start, eight at a time: one more edit after the previous front on the same
    # diagonal or either neighbour, each only where that edit stays inside the matrix. The previous
    # front's padding reads as unreached, far below any column.
    var lane_diagonals = Lanes()
    comptime for lane in range(FRONT_LANES):
        lane_diagonals[lane] = Int32(lane)
    var column_limit = Lanes(Int32(columns))
    var row_limit = Lanes(Int32(rows))
    var unreached = Lanes(UNREACHED)
    var reaches = Lanes(0)
    var diagonal = low
    while diagonal <= high:
        var diagonals = lane_diagonals + Int32(diagonal)
        var same = previous.unsafe_offset(diagonal).unsafe_load[width=FRONT_LANES]()
        var below = previous.unsafe_offset(diagonal - 1).unsafe_load[width=FRONT_LANES]()
        var above = previous.unsafe_offset(diagonal + 1).unsafe_load[width=FRONT_LANES]()
        var substituted = (same.lt(column_limit) & (same - diagonals).lt(row_limit)).select(same + 1, unreached)
        # A base of the first sequence against a gap, from the diagonal below.
        var deleted = below.lt(column_limit).select(below + 1, unreached)
        # A base of the second sequence against a gap, from the diagonal above.
        var inserted = (above - diagonals).le(row_limit).select(above, unreached)
        var entry = max(substituted, max(deleted, inserted))
        comptime if GATHERED_SLIDES:
            # Each lane group slides as it is computed, both sequences' next eight bases gathered per lane.
            # A lane past `high` reads the previous front past its padding, an older front's values,
            # so it slides nothing: the gather would read wherever such a column points.
            var inside = diagonals.le(Int32(high))
            var slid = gathered_slides(first, second, inside.select(entry, unreached), diagonals)
            current.unsafe_offset(diagonal).unsafe_store(slid)
            comptime if measure:
                reaches = max(reaches, inside.select(slid + slid - diagonals, Lanes(0)))
        else:
            current.unsafe_offset(diagonal).unsafe_store(entry)
        diagonal += FRONT_LANES
    for index in range(FRONT_PADDING):
        current[unsafe_offset=high + 1 + index] = UNREACHED
    comptime if GATHERED_SLIDES:
        return Int(reaches.reduce_max())

    # Then each slides over its matches, eight bases at a time; the sentinels stop it at the edge.
    # Two diagonals a turn, so their loads overlap and they share the loop's bookkeeping.
    var furthest = 0
    diagonal = low
    while diagonal + 1 <= high:
        var left = Int(current[unsafe_offset=diagonal])
        var right = Int(current[unsafe_offset=diagonal + 1])
        if left >= 0:
            left = slide(first, second, left, diagonal)
            current[unsafe_offset=diagonal] = Int32(left)
            comptime if measure:
                furthest = max(furthest, 2 * left - diagonal)
        if right >= 0:
            right = slide(first, second, right, diagonal + 1)
            current[unsafe_offset=diagonal + 1] = Int32(right)
            comptime if measure:
                furthest = max(furthest, 2 * right - diagonal - 1)
        diagonal += 2
    if diagonal <= high:
        var column = Int(current[unsafe_offset=diagonal])
        if column >= 0:
            column = slide(first, second, column, diagonal)
            current[unsafe_offset=diagonal] = Int32(column)
            comptime if measure:
                furthest = max(furthest, 2 * column - diagonal)
    return furthest


def diagonal_transition(
    profile: Profile,
    step_tenths: Int,
    mut fronts: DiagonalFronts,
    switch_setup: Int = -1,
) -> Probe:
    """The edit distance by diagonal transition, as WFA computes it, while it stays cheaper than a band.

    Score `s` reaches, on every diagonal, the furthest cell some path of `s` edits reaches; matches
    are free, so each front slides as far as they carry it. The distance is the first score whose
    front reaches the corner, after about `d²` diagonals and the matches along the way, where the
    band's sweep pays for every column whatever the distance. So near-identical pairs, the case the
    post on A*PA2 names as its weak spot, finish here.

    From `PROBE_START` on, the furthest anti-diagonal reached projects the distance, and once the
    search still to do, the projection squared less the score squared, passes what a band would cost
    (see `step_budget`), it stops and hands the projection on as the band's first bound. With
    `switch_setup`, it also stops once the two-ended search would finish the rest more cheaply.
    """
    var columns = profile.columns
    var rows = profile.rows
    var first = profile.column_codes.unsafe_ptr()
    var second = profile.row_codes.unsafe_ptr()
    var target = columns - rows
    fronts.offsets.clear()
    fronts.starts.clear()
    fronts.lows.clear()
    fronts.highs.clear()

    var start = slide_forward(first, second, 0, 0)
    fronts.starts.append(0)
    fronts.lows.append(0)
    fronts.highs.append(0)
    for _ in range(FRONT_PADDING):
        fronts.offsets.append(UNREACHED)
    fronts.offsets.append(Int32(start))
    for _ in range(FRONT_PADDING):
        fronts.offsets.append(UNREACHED)
    if target == 0 and start == columns:
        return Probe(0, 0, 0)
    var score = 0
    while score < PROBE_CEILING:
        score += 1
        var previous_low = fronts.lows[score - 1]
        var previous_start = fronts.starts[score - 1]
        var low = max(-score, -rows)
        var high = min(score, columns)
        var count = high - low + 1
        var row_start = len(fronts.offsets)
        fronts.starts.append(row_start)
        fronts.lows.append(low)
        fronts.highs.append(high)
        # Room for the front, its padding, and a last vector's spill past the end.
        fronts.offsets.resize(unsafe_uninit_length=row_start + count + 2 * FRONT_PADDING + FRONT_LANES)
        var offsets = fronts.offsets.unsafe_ptr()
        for index in range(FRONT_PADDING):
            offsets[unsafe_offset=row_start + index] = UNREACHED
        # The previous front and the new one, both indexed by diagonal.
        var previous = offsets.unsafe_origin_cast[MutUntrackedOrigin]().unsafe_offset(
            previous_start + FRONT_PADDING - previous_low
        )
        var current = offsets.unsafe_origin_cast[MutUntrackedOrigin]().unsafe_offset(row_start + FRONT_PADDING - low)

        # The projection is checked every `PROBE_STRIDE` scores, so most fronts skip measuring it.
        var checking = score >= PROBE_START and score % PROBE_STRIDE == 0
        var furthest: Int
        if checking:
            furthest = step_front[True](previous, current, low, high, columns, rows, first, second)
        else:
            furthest = step_front[False](previous, current, low, high, columns, rows, first, second)
        if target >= low and target <= high and Int(current[unsafe_offset=target]) == columns:
            fronts.offsets.resize(unsafe_uninit_length=row_start + count + 2 * FRONT_PADDING)
            return Probe(score, score, score - 1)
        fronts.offsets.resize(unsafe_uninit_length=row_start + count + 2 * FRONT_PADDING)
        if checking:
            var estimate = score * (columns + rows) // max(furthest, 1)
            # With `switch_setup`, give way to the two-ended search once what this front has left costs
            # more than that search would from scratch, setup included; the work done is spent.
            if switch_setup >= 0 and estimate * estimate - score * score > (
                estimate * estimate * ALIGNED_TWO_ENDED_PERCENT // 100 + switch_setup
            ):
                return Probe(-1, max(estimate, score + 1), score)
            # What is left to search, about `estimate² - score²` diagonals, against what a band
            # would cost; the work already done is spent either way. Only once the search has spent
            # its share of the cheapest band there could be (see `GIVE_UP_SHARE`).
            if estimate * estimate - score * score > noisy_budget(
                step_budget(columns, step_tenths, estimate), score
            ) and score * score * GIVE_UP_SHARE >= step_budget(columns, step_tenths, score):
                return Probe(-1, max(estimate, score + 1), score)
    return Probe(-1, PROBE_CEILING + 1, PROBE_CEILING)


def reversed_codes(codes: List[UInt8], count: Int, sentinel: UInt8) -> List[UInt8]:
    """The first `count` codes back to front, with `CODE_PADDING` sentinels after them."""
    var flipped = List[UInt8]()
    reversed_codes_into(flipped, codes, count, sentinel)
    return flipped^


def reversed_codes_into(mut flipped: List[UInt8], codes: List[UInt8], count: Int, sentinel: UInt8):
    """`reversed_codes` written over `flipped`, whose memory is kept."""
    flipped.resize(unsafe_uninit_length=count + CODE_PADDING)
    var target = flipped.unsafe_ptr()
    reversed_into(target, codes.unsafe_ptr(), count)
    target.unsafe_offset(count).unsafe_store(SIMD[DType.uint8, CODE_PADDING](sentinel))


comptime FRONT_RING = 2
"""Fronts a direction keeps in turn: the latest, and the one before, which an overlap check reads; the
next grows over the one before."""


struct FrontPair(Movable):
    """One direction's latest fronts, in a ring of `FRONT_RING` slots of `FRONT_WIDTH`, each indexed
    by diagonal, for up to `TWO_ENDED_LIMIT` scores."""

    var buffers: List[Int32]
    var low: Int
    var high: Int
    var previous_low: Int
    """The previous front's diagonals, which the ring still holds."""
    var previous_high: Int
    var score: Int
    var slot: Int
    """Which of the ring holds the latest front."""
    var furthest: Int
    """The furthest anti-diagonal, column plus row, the last measured front reached."""
    var history: DiagonalFronts
    """Every front so far, when `record`, laid out as `diagonal_transition` keeps them, for a traceback."""
    var record: Bool

    def __init__(out self, record: Bool = False):
        """A ring holding only score zero's front, still to be set, and a history when `record`."""
        # Every front writes its own diagonals and the padding either side before the next reads it,
        # so only the first front's surroundings need setting.
        self.buffers = List[Int32](capacity=FRONT_RING * FRONT_WIDTH)
        self.buffers.resize(unsafe_uninit_length=FRONT_RING * FRONT_WIDTH)
        self.low = 0
        self.high = 0
        self.previous_low = 0
        self.previous_high = -1
        self.score = 0
        self.slot = 0
        self.furthest = 0
        self.history = DiagonalFronts(reserve=record)
        self.record = record
        var first = self.front_mut(0)
        for diagonal in range(-FRONT_PADDING, FRONT_PADDING + 1):
            first[unsafe_offset=diagonal] = UNREACHED

    def reset(mut self):
        """As new, recording nothing, the ring's memory kept for a batch's next pair."""
        self.low = 0
        self.high = 0
        self.previous_low = 0
        self.previous_high = -1
        self.score = 0
        self.slot = 0
        self.furthest = 0
        self.history = DiagonalFronts(reserve=False)
        self.record = False
        var first = self.front_mut(0)
        for diagonal in range(-FRONT_PADDING, FRONT_PADDING + 1):
            first[unsafe_offset=diagonal] = UNREACHED

    def take_history(deinit self) -> DiagonalFronts:
        """The fronts kept, the rest given up."""
        return self.history^

    @always_inline
    def keep(mut self):
        """Copies the latest front, with its padding, onto the history, when recording: a check in place
        for a distance, which records nothing and once paid a call a step for it."""
        if self.record:
            self.kept()

    def kept(mut self):
        """`keep`'s copy."""
        var source = self.front(self.slot)
        var start = len(self.history.offsets)
        self.history.starts.append(start)
        self.history.lows.append(self.low)
        self.history.highs.append(self.high)
        var count = self.high - self.low + 1 + 2 * FRONT_PADDING
        self.history.offsets.resize(unsafe_uninit_length=start + count)
        var target = self.history.offsets.unsafe_ptr().unsafe_offset(start)
        var origin = source.unsafe_offset(self.low - FRONT_PADDING)
        # Sixteen at a time: the compiler keeps a plain loop here scalar.
        var index = 0
        while index + 16 <= count:
            target.unsafe_offset(index).unsafe_store(origin.unsafe_offset(index).unsafe_load[width=16]())
            index += 16
        while index < count:
            target[unsafe_offset=index] = origin[unsafe_offset=index]
            index += 1

    @inline(.always)
    def previous(self) -> Int:
        """The ring slot of the front before the latest."""
        return self.slot - 1 if self.slot > 0 else FRONT_RING - 1

    @inline(.always)
    def advance(
        mut self, codes: ImmPointer[UInt8, _], others: ImmPointer[UInt8, _], measure: Bool, columns: Int, rows: Int
    ):
        """One more score: the next front in the ring, grown from the latest."""
        var previous = self.front(self.slot)
        var next_slot = self.slot + 1 if self.slot + 1 < FRONT_RING else 0
        var current = self.front_mut(next_slot)
        self.score += 1
        var low = max(-self.score, -rows)
        var high = min(self.score, columns)
        for index in range(1, FRONT_PADDING + 1):
            current[unsafe_offset=low - index] = UNREACHED
        if measure:
            self.furthest = step_front[True](previous, current, low, high, columns, rows, codes, others)
        else:
            _ = step_front[False](previous, current, low, high, columns, rows, codes, others)
        self.previous_low = self.low
        self.previous_high = self.high
        self.low = low
        self.high = high
        self.slot = next_slot
        self.keep()

    @inline(.always)
    def front(self, which: Int) -> ImmPointer[Int32, ImmUntrackedOrigin]:
        """Front `which` of the ring, to read, indexed by diagonal."""
        return (
            self.buffers.unsafe_ptr()
            .unsafe_origin_cast[ImmUntrackedOrigin]()
            .unsafe_offset(which * FRONT_WIDTH + FRONT_CENTER)
        )

    @inline(.always)
    def front_mut(mut self, which: Int) -> MutPointer[Int32, MutUntrackedOrigin]:
        """Front `which` of the ring, to write, indexed by diagonal. Untracked, as a front is written
        while `front` reads the one before it from the same buffer."""
        return (
            self.buffers.unsafe_ptr()
            .unsafe_origin_cast[MutUntrackedOrigin]()
            .unsafe_offset(which * FRONT_WIDTH + FRONT_CENTER)
        )


@inline(.always)
def overlap(
    forward: ImmPointer[Int32, _],
    forward_low: Int,
    forward_high: Int,
    backward: ImmPointer[Int32, _],
    back_low: Int,
    back_high: Int,
    target: Int,
    columns: Int,
) -> Int:
    """A diagonal whose forward front reaches past where the backward front comes back to, or
    `NO_DIAGONAL`; the backward front's diagonal `target - k` mirrors the forward's `k`."""
    var low = max(forward_low, target - back_high)
    var high = min(forward_high, target - back_low)
    # An unreached diagonal holds a value far below zero, so the sum alone rules it out.
    var needed = SIMD[DType.int32, FRONT_LANES](Int32(columns))
    var diagonal = low
    while diagonal + FRONT_LANES - 1 <= high:
        var reached = forward.unsafe_offset(diagonal).unsafe_load[width=FRONT_LANES]()
        var back = (
            backward.unsafe_offset(target - diagonal - FRONT_LANES + 1).unsafe_load[width=FRONT_LANES]().reversed()
        )
        if (reached + back).ge(needed).reduce_or():
            break
        diagonal += FRONT_LANES
    while diagonal <= high:
        if Int(forward[unsafe_offset=diagonal]) + Int(backward[unsafe_offset=target - diagonal]) >= columns:
            return diagonal
        diagonal += 1
    return NO_DIAGONAL


@fieldwise_init
struct Meeting(ImplicitlyCopyable, TrivialRegisterPassable):
    """Where the two-ended search's fronts met: the forward front's furthest cell on a diagonal both
    reached, and each side's score there; or, when it gave up, `probe` alone."""

    var probe: Probe
    var diagonal: Int
    var column: Int
    var forward_score: Int
    var backward_score: Int


def two_ended_distance(profile: Profile, step_tenths: Int, ceiling: Int = Int.MAX) -> Probe:
    """The edit distance by `two_ended`, keeping no history, giving up past `ceiling` steps too."""
    var first_back = List[UInt8]()
    var second_back = List[UInt8]()
    var ahead = FrontPair()
    var behind = FrontPair()
    return two_ended_distance(profile, step_tenths, ceiling, first_back, second_back, ahead, behind)


def two_ended_distance(
    profile: Profile,
    step_tenths: Int,
    ceiling: Int,
    mut first_back: List[UInt8],
    mut second_back: List[UInt8],
    mut ahead: FrontPair,
    mut behind: FrontPair,
) -> Probe:
    """`two_ended_distance` over memory a batch's worker keeps from pair to pair: the reversed codes
    and the two fronts' rings, begun afresh."""
    reversed_codes_into(first_back, profile.column_codes, profile.columns, FIRST_SENTINEL)
    reversed_codes_into(second_back, profile.row_codes, profile.rows, SECOND_SENTINEL)
    ahead.reset()
    behind.reset()
    return two_ended(profile, first_back, second_back, step_tenths, ahead, behind, ceiling).probe


def two_ended(
    profile: Profile,
    first_back: List[UInt8],
    second_back: List[UInt8],
    step_tenths: Int,
    mut ahead: FrontPair,
    mut behind: FrontPair,
    ceiling: Int = Int.MAX,
) -> Meeting:
    """The edit distance by diagonal transition from both ends at once, as BiWFA scores, while cheap.

    One front grows from the start and one from the end, over the reversed sequences, a score at a
    time each in turn, and the distance is the first total score at which they overlap: on some
    diagonal the forward front reaches at least as far as the backward one comes back. An overlap
    joins a real path from the start to one to the end, so the total is at least the distance; and
    every cost level an optimal path passes splits it in two whose halves the fronts have reached
    by the time their scores add up to the distance, edit costs never falling along a diagonal. So
    the first overlap is exact, after about half the diagonals of one front grown alone, keeping
    only the last two fronts each way.

    It stops as `diagonal_transition` does, once the search still to do passes the budget (see
    `step_budget`), or `ceiling` steps where that is less, projecting the distance from both fronts'
    progress.
    """
    var columns = profile.columns
    var rows = profile.rows
    var target = columns - rows
    var first = profile.column_codes.unsafe_ptr()
    var second = profile.row_codes.unsafe_ptr()
    var first_reversed = first_back.unsafe_ptr()
    var second_reversed = second_back.unsafe_ptr()
    var limit = TWO_ENDED_LIMIT
    ahead.front_mut(0)[unsafe_offset=0] = Int32(slide_forward(first, second, 0, 0))
    behind.front_mut(0)[unsafe_offset=0] = Int32(slide_forward(first_reversed, second_reversed, 0, 0))
    ahead.keep()
    behind.keep()

    @inline(.always)
    def overlapping(
        ahead: FrontPair, behind: FrontPair, forward_earlier: Bool, backward_earlier: Bool
    ) {imm target, imm columns} -> Int:
        """`overlap` of the forward front, or the one before it, and the backward front, or the one
        before it."""
        return overlap(
            ahead.front(ahead.previous() if forward_earlier else ahead.slot),
            ahead.previous_low if forward_earlier else ahead.low,
            ahead.previous_high if forward_earlier else ahead.high,
            behind.front(behind.previous() if backward_earlier else behind.slot),
            behind.previous_low if backward_earlier else behind.low,
            behind.previous_high if backward_earlier else behind.high,
            target,
            columns,
        )

    @inline(.always)
    def met(
        ahead: FrontPair, behind: FrontPair, diagonal: Int, forward_earlier: Bool, backward_earlier: Bool
    ) -> Meeting:
        """Where the fronts `overlapping` paired met on `diagonal`: the forward front's column there, and
        each side's score."""
        var column = Int(ahead.front(ahead.previous() if forward_earlier else ahead.slot)[unsafe_offset=diagonal])
        var forward_score = ahead.score - 1 if forward_earlier else ahead.score
        var backward_score = behind.score - 1 if backward_earlier else behind.score
        var total = forward_score + backward_score
        return Meeting(Probe(total, total, total - 1), diagonal, column, forward_score, backward_score)

    var meeting = overlapping(ahead, behind, False, False)
    if meeting != NO_DIAGONAL:
        return met(ahead, behind, meeting, False, False)
    var total = 0
    while ahead.score < limit and behind.score < limit:
        total += 1
        # Each side measures its progress on its step just before a check, every `PROBE_STRIDE`, from
        # when each has taken `PROBE_START` edits, so the projection has as many to go on as one front's.
        var checking = total >= 2 * PROBE_START and total % PROBE_STRIDE == 0
        var measuring = total + 1 >= 2 * PROBE_START and (total + 1) % PROBE_STRIDE == 0
        if total % 2 == 1:
            ahead.advance(first, second, checking or measuring, columns, rows)
        else:
            behind.advance(first_reversed, second_reversed, checking or measuring, columns, rows)
        # Checked after the backward front's steps alone: a first overlap one score sooner is found
        # then, against the backward front's previous score, which the ring still holds.
        if total % 2 == 0:
            meeting = overlapping(ahead, behind, False, False)
            if meeting != NO_DIAGONAL:
                var sooner = overlapping(ahead, behind, False, True)
                if sooner != NO_DIAGONAL:
                    return met(ahead, behind, sooner, False, True)
                return met(ahead, behind, meeting, False, False)
        if checking:
            var estimate = two_ended_gives_up(
                total, ahead.furthest + behind.furthest, columns, rows, step_tenths, ceiling
            )
            if estimate >= 0:
                return Meeting(Probe(-1, estimate, total), 0, 0, 0, 0)
    return Meeting(Probe(-1, 2 * limit + 1, total), 0, 0, 0, 0)


comptime AGREEMENT = 3
"""Two projections agree when neither exceeds the other by more than one part in this many."""


def trusted_projection(search: Probe, projected: Probe) -> Bool:
    """Whether one front's projection from its first edits agrees with the search's from many more.

    On pairs whose errors spread along them, as mutated sequences' do, the two land within a few
    percent of each other and of the distance, and the band aims straight at them. Real reads gather
    errors at their ends, a primer or a tail, and there the two disagree several times over; neither
    is then worth aiming at.
    """
    if search.floor < 2 * PROBE_START or search.estimate <= search.floor:
        return True
    var first = projected.estimate
    var later = search.estimate
    return first * AGREEMENT <= later * (AGREEMENT + 1) and later * AGREEMENT <= first * (AGREEMENT + 1)


def two_ended_gives_up(
    total: Int, reached: Int, columns: Int, rows: Int, step_tenths: Int, ceiling: Int = Int.MAX
) -> Int:
    """The two-ended search's projected distance once what is left of it passes the budget, or the
    `ceiling` where that is less, or -1."""
    var estimate = total * (columns + rows) // max(reached, 1)
    # A ceiling the search's own fallback costs, widened by the projection's noise (see `noisy_budget`),
    # so a projection that overshoots does not give up a search that would have finished cheaper.
    var cap = noisy_budget(ceiling, total) if ceiling != Int.MAX else ceiling
    # What is left, about half of `estimate² - total²` diagonals, against what a band costs.
    # Both fronts' steps, about 0.57 ns per square edit with the overlap check, in one-front steps.
    if (estimate * estimate - total * total) * TWO_ENDED_PERCENT // 100 > min(
        step_budget(columns, step_tenths, estimate), cap
    ) and total * total * TWO_ENDED_PERCENT // 100 * GIVE_UP_SHARE >= min(
        step_budget(columns, step_tenths, total), cap
    ):
        return max(estimate, total + 1)
    return -1


comptime FULL_WORD_HUNDREDTHS = 35 if simd_width_of[DType.uint64]() >= 8 else (
    45 if simd_width_of[DType.uint64]() >= 4 else 77
)
"""What the whole matrix's sweep costs a word, a column of 64 rows, in hundredths of a diagonal step,
by vector width: measured on 1 kbp pairs at 15%, where both run to the end, as 0.64 ns a word against
1.83 ns a step on AVX-512 (a Skylake-X at 3.3 GHz), 1.05 against 2.34 on AVX2 (the same machine
built for Haswell) and 1.03 against 1.33 on NEON (an M2): the wider the vectors, the cheaper the
sweep beside the diagonal transition, which steps a diagonal at a time."""


def full_matrix_steps(columns: Int, rows: Int) -> Int:
    """What sweeping a pair's whole matrix costs, in diagonal steps (see `FULL_WORD_HUNDREDTHS`)."""
    return columns * ((rows + 63) // 64) * FULL_WORD_HUNDREDTHS // 100


def grow_to(profile: Profile, mut ahead: DiagonalFronts, behind: DiagonalFronts, distance: Int):
    """Grows the forward fronts `ahead` kept on to `distance`, which the two-ended search proved, each
    score kept only on the diagonals an optimal path passes: where the backward fronts `behind` kept,
    at the rest of the distance, come back as far. So `trace_diagonals` then traces the path WFA2-lib's
    rule picks, whose every step takes the furthest source, which an optimal path passes too and the
    pruning keeps, while each score grows only the few diagonals left. The last front the search kept
    is pruned the same way first, so the first score grown is as narrow as the rest."""
    var columns = profile.columns
    var rows = profile.rows
    var target = columns - rows
    var first = profile.column_codes.unsafe_ptr()
    var second = profile.row_codes.unsafe_ptr()
    var last = len(behind.lows) - 1

    var top = len(ahead.lows) - 1
    if top >= distance:
        return
    # The last kept front, pruned in place: what it loses reads unreached, and its range narrows. Only
    # diagonals the backward front of the rest also holds can meet it, eight at a time.
    var low = ahead.lows[top]
    var high = ahead.highs[top]
    var row = ahead.offsets.unsafe_ptr().unsafe_offset(ahead.starts[top] + FRONT_PADDING - low)
    var kept_low = high + 1
    var kept_high = low - 1
    var rest = distance - top
    if rest <= last:
        var back_low = behind.lows[rest]
        var back_high = behind.highs[rest]
        var back = behind.offsets.unsafe_ptr().unsafe_offset(behind.starts[rest] + FRONT_PADDING - back_low)
        var from_diagonal = max(low, target - back_high)
        var to_diagonal = min(high, target - back_low)
        comptime Lanes = SIMD[DType.int32, FRONT_LANES]
        var needed = Lanes(Int32(columns))
        var diagonal = from_diagonal
        while diagonal <= to_diagonal:
            if diagonal + FRONT_LANES - 1 <= to_diagonal:
                var reached = row.unsafe_offset(diagonal).unsafe_load[width=FRONT_LANES]()
                var behind_reached = (
                    back.unsafe_offset(target - diagonal - FRONT_LANES + 1).unsafe_load[width=FRONT_LANES]().reversed()
                )
                var live = reached.ge(Lanes(0)) & behind_reached.ge(Lanes(0)) & (reached + behind_reached).ge(needed)
                if live.reduce_or():
                    for lane in range(FRONT_LANES):
                        if live[lane]:
                            kept_low = min(kept_low, diagonal + lane)
                            kept_high = diagonal + lane
                row.unsafe_offset(diagonal).unsafe_store(live.select(reached, Lanes(UNREACHED)))
                diagonal += FRONT_LANES
                continue
            var column = Int(row[unsafe_offset=diagonal])
            var behind_column = Int(back[unsafe_offset=target - diagonal])
            if column >= 0 and behind_column >= 0 and column + behind_column >= columns:
                kept_low = min(kept_low, diagonal)
                kept_high = diagonal
            else:
                row[unsafe_offset=diagonal] = UNREACHED
            diagonal += 1
        for dead in range(low, min(from_diagonal, high + 1)):
            row[unsafe_offset=dead] = UNREACHED
        for dead in range(max(to_diagonal + 1, low), high + 1):
            row[unsafe_offset=dead] = UNREACHED
    else:
        kept_low = low
        kept_high = high
    if kept_low > kept_high:
        kept_low = low
        kept_high = low - 1
    ahead.starts[top] += kept_low - low
    ahead.lows[top] = kept_low
    ahead.highs[top] = kept_high
    # Room for every score to come at once, so the rows grow without copying what is kept.
    var scores = distance - top
    ahead.offsets.reserve(len(ahead.offsets) + scores * (8 + 2 * FRONT_PADDING) + kept_high - kept_low + 1)
    ahead.starts.reserve(len(ahead.starts) + scores)
    ahead.lows.reserve(len(ahead.lows) + scores)
    ahead.highs.reserve(len(ahead.highs) + scores)
    for score in range(top + 1, distance + 1):
        var previous_low = ahead.lows[score - 1]
        var previous_high = ahead.highs[score - 1]
        var previous = ahead.offsets.unsafe_ptr().unsafe_offset(ahead.starts[score - 1] + FRONT_PADDING - previous_low)
        var new_low = max(previous_low - 1, -rows)
        var new_high = min(previous_high + 1, columns)
        var row_start = len(ahead.offsets)
        # Room for the row and its padding, every slot written below; the row is narrowed to what
        # lives once it is grown.
        ahead.offsets.resize(unsafe_uninit_length=row_start + max(new_high - new_low + 1, 0) + 2 * FRONT_PADDING)
        var current = ahead.offsets.unsafe_ptr().unsafe_offset(row_start + FRONT_PADDING - new_low)
        for pad in range(1, FRONT_PADDING + 1):
            current[unsafe_offset=new_low - pad] = UNREACHED
            current[unsafe_offset=new_high + pad] = UNREACHED
        rest = distance - score
        # The backward front of the rest, as a row of the forward diagonals it mirrors.
        var back_low = behind.lows[rest]
        var back_high = behind.highs[rest]
        var back = behind.offsets.unsafe_ptr().unsafe_offset(behind.starts[rest] + FRONT_PADDING - back_low)
        kept_low = new_high + 1
        kept_high = new_low - 1
        for diagonal in range(new_low, new_high + 1):
            # A diagonal the backward front of the rest does not reach meets nothing there.
            var mirrored = target - diagonal
            if mirrored < back_low or mirrored > back_high or back[unsafe_offset=mirrored] < 0:
                current[unsafe_offset=diagonal] = UNREACHED
                continue
            # One more edit after the previous front, as `best_source` takes it; the previous row's
            # padding reads unreached either side.
            var best = -1
            var same = Int(previous[unsafe_offset=diagonal])
            if same >= 0 and same < columns and same - diagonal < rows:
                best = same + 1
            var left = Int(previous[unsafe_offset=diagonal - 1])
            if left >= 0 and left < columns and left + 1 > best:
                best = left + 1
            var up = Int(previous[unsafe_offset=diagonal + 1])
            if up >= 0 and up - diagonal - 1 < rows and up > best:
                best = up
            var column = slide_forward(first, second, best, best - diagonal) if best >= 0 else -1
            if column >= 0 and column + Int(back[unsafe_offset=mirrored]) >= columns:
                current[unsafe_offset=diagonal] = Int32(column)
                kept_low = min(kept_low, diagonal)
                kept_high = diagonal
            else:
                current[unsafe_offset=diagonal] = UNREACHED
        if kept_low > kept_high:
            kept_low = new_low
            kept_high = new_low - 1
        ahead.starts.append(row_start + kept_low - new_low)
        ahead.lows.append(kept_low)
        ahead.highs.append(kept_high)


def trace_diagonals(profile: Profile, fronts: DiagonalFronts, distance: Int, mut moves: List[UInt8]):
    """The optimal path `diagonal_transition` found, as moves right to left, like `trace_back`.

    From the corner, each score's front is undone: the matches it slid over, then the edit that
    reached its start from the furthest front of one score less, which `best_source` finds again.
    """
    var columns = profile.columns
    var rows = profile.rows
    var diagonal = columns - rows
    var column = columns
    var score = distance
    while score > 0:
        var source = best_source(fronts, score, diagonal, columns, rows)
        var best = source[0]
        append_diagonals(moves, column - best)
        var move = source[1]
        moves.append(move)
        if move == DIAGONAL:
            column = best - 1
        elif move == LEFT:
            column = best - 1
            diagonal -= 1
        else:
            column = best
            diagonal += 1
        score -= 1
    append_diagonals(moves, column)
