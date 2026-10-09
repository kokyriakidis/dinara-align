# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Ported from `pa-bitpacking` in A*PA (https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner,
# commit bf2e14e), by Ragnar Groot Koerkamp and Pesho Ivanov, itself translated from Edlib.
"""
A*PA2-full's gap-chaining seed heuristic, with A*PA's inexact seeds and local pruning: a lower bound on
the cost from any cell to the end, which narrows the band (see `SeedHeuristic`).
"""

from std.bit import count_leading_zeros, count_trailing_zeros, pop_count
from std.math import ceildiv
from std.sys import simd_width_of
from std.sys.intrinsics import likely, prefetch, PrefetchOptions

from .bit_parallel import folded, PAIRED_GROUPS, Profile
from .diagonal import extend


comptime SEED_LENGTH = 12
"""Bases per seed: A*PA2-full's `k`, long enough that a random match is rare, short enough that most
seeds of a moderately divergent pair still match exactly."""


comptime INEXACT_LENGTH = 16
"""Bases per seed when seeds may match with one edit: A*PA's inexact seeds, long enough that a random
one-edit match is still rare, each seed then charging two edits to a path that matches it nowhere."""


comptime INEXACT_COLUMNS = 65_536
"""Columns from which the seeds may match with one edit. Their setup grows with the length and the band
they save with its square: on real reads of 16 to 70 kbp exact seeds win however divergent the pair,
and from about 100 kbp, and on uniform errors of one in seven, inexact ones do."""


comptime INEXACT_DIVERGENCE = 10
"""Columns per projected edit at or below which the seeds may match with one edit: past about one edit
in ten bases most exact seeds are broken, while most still match within one edit."""


comptime INEXACT_CHAINED = 40 if simd_width_of[DType.uint64]() <= 2 else (
    30 if simd_width_of[DType.uint64]() <= 4 else 20
)
"""Percent of the exact seeds chained from the origin below which the seeds are rebuilt to match within
one edit. Where that pays depends on the machine: the band is vector work and the seeds' setup scalar,
memory-bound lookups, so the wider the vectors beside the scalar core, the less the band inexact seeds
narrow is worth beside the setup they cost. On an M2, two 64-bit lanes to a NEON register, 40 percent,
about one edit in fifteen bases on real reads; on a Skylake-X at a fixed 3.3 GHz, eight to an AVX-512
register, its seeds' setup twice as slow as the M2's and its band no slower, 20 percent, short of which
the long reads ran faster on exact seeds. Four lanes, AVX2, are between the two and unmeasured."""


comptime INEXACT_CHAINED_SPREAD = 20
"""`INEXACT_CHAINED` for a pair whose two projections agree (see `trusted_projection`), its errors
spread along it as mutated sequences' are, not gathered at an end as real reads' are. Spread errors
break fewer exact seeds' chains than gathered ones for the same narrowing of the band, so exact seeds
still pay down to this share: on the M2, 100 kbp pairs at 6 to 8% aligned in 3.7 to 4.3 ms on exact
seeds against 5.7 ms rebuilt, where the real reads still wanted `INEXACT_CHAINED`'s 40 percent."""


comptime HALF_BITS = INEXACT_LENGTH
"""Bits in half an inexact seed's two-bit code: a one-edit match matches one half exactly, so each half
indexes a table of this many bits."""


comptime ENTRY_BITS = 32
"""Bits below an inexact seed's index in its half tables' entries, its two-bit code there: one load
gives a lookup both."""


comptime SCAN_BATCH = 16
"""Rows whose half tables' buckets are looked up together before any is searched (see
`inexact_scan`): enough misses in flight at once to hide most of their wait."""


comptime SEED_EDITS = 1500
"""Projected edits from which a band prunes with the seed heuristic. Its setup costs about ten
nanoseconds a column; what it saves grows with the distance, the band otherwise sweeping rows in
proportion to it, and passes the setup about here."""


comptime SEED_DIVERGENCE = 7
"""Columns per projected edit below which the seeds are left out: past about one edit in seven bases,
almost no seed survives local pruning, the heuristic is little more than an edit a seed, and its
setup is not repaid."""


comptime SEED_COLUMNS = 16_384 if not PAIRED_GROUPS else 86_016
"""Columns from which a band prunes with the seeds whatever the projection: its setup is a small share
of any band this long, and on real reads, whose errors gather at the ends, the projection that gates
shorter pairs can put a divergence of one edit in ten at one in two.

On AVX-512, whose band sweeps two groups at a time (see `PAIRED_GROUPS`), the band runs fast enough
beside the seeds' setup that they pay only from between 80 and 90 kbp: below, the Skylake-X aligned
real reads of 16 to 64 kbp 14% faster without them, uniform 30 kbp pairs at 5% 30% faster, and at
100 kbp seeds won by 7 to 24%. With AVX2 they broke even on those reads and won from 64 kbp, as on
the M2, which keeps the shorter gates too (see `SHORT_SEEDS`)."""


comptime SHORT_SEEDS = not PAIRED_GROUPS
"""Whether a pair shorter than `SEED_COLUMNS` may still take seeds on its projection (see `SEED_EDITS`):
not on AVX-512, where seeds lost on every such pair tried."""


comptime SEEDED_GROWTH = 4
"""How many times a seeded band's margin over the heuristic at the origin grows after a round that
died early."""


comptime SEEDED_TRUST_SHARE = 8
"""A seeded round must cross one part in this many of the columns before its death projects the next
bound."""


comptime CHAINED_SHARE = 32
"""With seeds, the band's first bound starts at the origin's when at least one seed in this many is
chained there; with fewer, the bound is little more than an edit a seed, and the projection leads."""


comptime SEED_SLACK = 64
"""With seeds, the band's first bound reaches at least this far past the heuristic at the origin."""


comptime LOOKAHEAD_SEEDS = 14
"""Seeds a match's local pruning looks ahead, its own included: A*PA2-full's `p`."""


comptime LAYER_SLOTS = 8
"""Starts a layer holds in place before the rest spill into a staircase of its own."""


@inline(.always)
def differing(first: UInt64, second: UInt64) -> UInt64:
    """The bases two two-bit codes disagree on, each as the low bit of its pair."""
    var mismatched = first ^ second
    return (mismatched | (mismatched >> 1)) & UInt64(0x5555555555555555)


@inline(.always)
def covers(prefix: UInt64, suffix: UInt64) -> Bool:
    """Whether the leading bases one comparison agrees on and the trailing bases another agrees on,
    both of the same bases, cover all of them between them.

    The leading run reaches past the trailing run's start when every base the first comparison
    differs on lies below the lowest the second differs on, an unsigned comparison of the first's
    differing bits against the second's lowest, with no count of leading or trailing zeros.
    """
    return suffix == 0 or prefix < (suffix & (0 - suffix))


@inline(.always)
def ends_within_one_edit(shorter: UInt64, level: UInt64, longer: UInt64, code: UInt64) -> Bool:
    """`within_one_edit` for the windows ending a seed's length on from `level`'s start, which start
    at `shorter`, `level` and `longer`, less any whose left half matches the seed's too: that one is
    the left half's lookup's to take."""
    var left = code >> UInt64(HALF_BITS)
    return within_one_edit(
        shorter >> 2,
        level,
        longer,
        level,
        code,
        left != shorter >> UInt64(HALF_BITS),
        left != level >> UInt64(HALF_BITS),
        left != longer >> UInt64(HALF_BITS),
    )


@inline(.always)
def first_reachable(entries: ImmPointer[UInt64, _], start: Int, end: Int, seed: Int) -> Int:
    """The first of a bucket's entries, `[start, end)` in seed order, of `seed` or a later seed: a
    step at a time in a small bucket, by bisection in a large one, a repeat's."""
    comptime LINEAR = 8
    var low = start
    var high = end
    while high - low > LINEAR:
        var middle = (low + high) // 2
        if Int(entries[unsafe_offset=middle] >> UInt64(ENTRY_BITS)) < seed:
            low = middle + 1
        else:
            high = middle
    while low < high and Int(entries[unsafe_offset=low] >> UInt64(ENTRY_BITS)) < seed:
        low += 1
    return low


@inline(.always)
def within_one_edit(
    shorter: UInt64,
    level: UInt64,
    head: UInt64,
    tail: UInt64,
    code: UInt64,
    shorter_open: Bool = True,
    level_open: Bool = True,
    longer_open: Bool = True,
) -> Bool:
    """Whether a seed's code is within one edit of any of three windows `try_windows` would keep:
    `shorter`, one base shorter than the seed, `level`, as long, or the window one longer, whose
    first and last `INEXACT_LENGTH` bases are `head` and `tail`, each only while open.

    The three tests folded together, with no branch to mispredict: an equal-length window within
    one substitution differs on at most one base, a shorter one is the seed less a base when the
    seed's first and last bases but one cover it between them, and a longer one the seed plus a base
    when its first and last `INEXACT_LENGTH` bases cover the seed. An exact match leaves out the
    shorter and longer windows, as `try_windows` does.
    """
    comptime K = INEXACT_LENGTH
    comptime SHORTER = (UInt64(1) << UInt64(2 * K - 2)) - 1
    var whole = differing(level, code)
    var substituted = (whole & (whole - 1)) == 0
    var deleted = covers(differing(shorter, code >> 2), differing(shorter, code & SHORTER))
    var inserted = covers(differing(head, code), differing(tail, code))
    return (substituted & level_open) | (((deleted & shorter_open) | (inserted & longer_open)) & (whole != 0))


@always_inline
def first_at_least(values: List[Int32], wanted: Int) -> Int:
    """The first place in the rising `values` holding `wanted` or more, `len(values)` past them all."""
    var low = 0
    var high = len(values)
    var at = values.unsafe_ptr()
    while low < high:
        var middle = (low + high) // 2
        if Int(at[unsafe_offset=middle]) < wanted:
            low = middle + 1
        else:
            high = middle
    return low


struct SeedHeuristic(Movable):
    """A*PA2-full's gap-chaining seed heuristic, without match pruning.

    The first sequence is cut into disjoint seeds of `SEED_LENGTH` bases, and every exact occurrence
    of a seed in the second sequence is a match. A path that misses a seed must spend an edit on it,
    so the seeds still ahead of a cell, its potential `P`, bound the cost to the end, less one for
    each match a path can still chain. Chaining charges the gap between matches: in the coordinates
    `T(i, j) = (i - j - P(i), j - i - P(i))` one match can follow another exactly when the second's
    start lies above and right of the first's end, so the longest chain from a cell is a dominance
    query, answered from layers of match starts, layer `v` holding the starts from which `v` matches
    chain (see `score`).

    With inexact seeds, as A*PA's `r = 2`, seeds are `INEXACT_LENGTH` bases, a match may cost one
    edit, and a seed matched nowhere costs a path two: the potential counts two a seed, a match
    scores two less its cost, and the transform is the same in those units. A match scoring two
    starts in two layers, so the layers still nest.

    The heuristic never overestimates the cost to the end, and changes by at most `cost` a row down a
    column or up it (see `climb`), so the band's pruning stays exact and its jumps stay valid. With no
    seeds it is the plain gap to the end's diagonal.
    """

    var columns: Int
    var rows: Int
    var length: Int
    """Bases per seed."""
    var cost: Int
    """What a path pays crossing a seed it matches nowhere: one for exact seeds, two for inexact."""
    var seeds: Int
    var counted: Int
    """The seeds the potential counts: all of them, but for those holding a symbol past `ACGT`."""
    var remaining: List[Int32]
    """Per seed, the counted seeds from it to the end, one past the last; empty when every seed counts.

    A seed holding a symbol past `ACGT`, an `N`, has no matches and costs a path nothing: the heuristic
    is then the same on the other seeds alone, which bound the cost as all of them do. Matching it with
    `N` folded into a base instead would flood a run of `N` in both sequences with matches, a run of
    `A` against a run of `A`."""
    var slot_x: List[Int32]
    """Per layer, `LAYER_SLOTS` slots for the transformed starts of the matches from which that many
    matches chain; a layer rarely holds more than three, and any past the slots spill over."""
    var slot_y: List[Int32]
    var counts: List[Int32]
    """Per layer, how many of its slots hold a start; only a layer whose slots are full can have spilled."""
    var spill_x: List[List[Int32]]
    """Per layer, the starts past its slots as a staircase: none at or above and right of another, so by
    `x` rising `y` falls, and whether one lies at or above and right of a point is one binary search. In
    a repeat a layer holds thousands; a chain walked a start at a time took poly-A of 32,000 bases
    against 16,000 seventeen seconds."""
    var spill_y: List[List[Int32]]
    var hint: Int
    """The layer the last query ended on: neighbouring queries land near it."""

    def __init__(out self, columns: Int, rows: Int):
        """No seeds: the gap heuristic."""
        self.columns = columns
        self.rows = rows
        self.length = SEED_LENGTH
        self.cost = 1
        self.seeds = 0
        self.counted = 0
        self.remaining = List[Int32]()
        self.slot_x = List[Int32]()
        self.slot_y = List[Int32]()
        self.counts = List[Int32]()
        self.spill_x = List[List[Int32]]()
        self.spill_y = List[List[Int32]]()
        self.hint = 0

    def __init__(
        out self,
        profile: Profile,
        inexact: Bool = False,
        choose: Bool = False,
        cutoff: Int = INEXACT_CHAINED,
    ):
        """Seeds of the profile's first sequence, matched exactly in its second, or with `inexact`
        within one edit. With `choose`, exact seeds, rebuilt inexact when fewer than `cutoff` percent of
        them chain from the origin.

        A match is kept only if a path from its start crosses the next `LOOKAHEAD_SEEDS` seeds for
        less than they would cost unmatched (see `worth_keeping`).
        """
        self = Self(profile.columns, profile.rows)
        if profile.extended:
            var first = folded(profile.column_codes)
            var second = folded(profile.row_codes)
            self.build_choosing(first, second, profile.column_codes, inexact, choose, cutoff)
        else:
            self.build_choosing(profile.column_codes, profile.row_codes, List[UInt8](), inexact, choose, cutoff)

    def build_choosing(
        mut self,
        first_codes: List[UInt8],
        second_codes: List[UInt8],
        symbols: List[UInt8],
        inexact: Bool,
        choose: Bool,
        cutoff: Int,
    ):
        """The seeds of `first_codes` matched in `second_codes`, both two bits a base, rebuilt inexact
        as `__init__` says; `symbols` as `build` takes them."""
        self.build(first_codes, second_codes, symbols, inexact)
        if choose and not inexact and self.seeds > 0:
            var chained = self.counted - self.h(0, 0)
            if chained * 100 < cutoff * self.counted:
                self.build(first_codes, second_codes, symbols, True)

    def build(mut self, first_codes: List[UInt8], second_codes: List[UInt8], symbols: List[UInt8], inexact: Bool):
        """The seeds, their matches, and the layers, from scratch. `symbols` are the first sequence's
        own codes when some lie past `ACGT`, folded in `first_codes`, else empty: the seeds holding
        one go uncounted (see `remaining`)."""
        self.length = INEXACT_LENGTH if inexact else SEED_LENGTH
        self.cost = 2 if inexact else 1
        self.seeds = self.columns // self.length
        self.counted = self.seeds
        self.remaining.clear()
        if len(symbols) > 0:
            self.count_seeds(symbols)
            if self.counted == 0:
                # No seed left to count: the gap heuristic.
                self.seeds = 0
                self.remaining.clear()
        self.slot_x.clear()
        self.slot_y.clear()
        self.counts.clear()
        self.spill_x.clear()
        self.spill_y.clear()
        self.slot_x.reserve(LAYER_SLOTS * (self.seeds + 1))
        self.slot_y.reserve(LAYER_SLOTS * (self.seeds + 1))
        self.counts.reserve(self.seeds + 1)
        self.hint = 0
        self.add_sentinel()
        if self.seeds == 0 or self.rows < self.length + 1:
            return
        var first = first_codes.unsafe_ptr()
        var second = second_codes.unsafe_ptr()
        # Each match a seed and its start row, and for inexact seeds its end row and cost; an exact
        # match ends a seed's length further down, at no cost.
        var found_seed = List[Int32]()
        var found_row = List[Int32]()
        var found_end = List[Int32]()
        var found_cost = List[Int32]()
        if inexact:
            self.inexact_matches(first, second, found_seed, found_row, found_end, found_cost)
        else:
            self.exact_matches(first, second, found_seed, found_row)

        # Bucketed by seed: a start can dominate another match's end only from a later seed, so taking
        # the seeds last first is an order the layers can be built in, no sort.
        var firsts = List[Int32](length=self.seeds + 1, fill=0)
        for index in range(len(found_seed)):
            firsts[Int(found_seed[index]) + 1] += 1
        for seed in range(self.seeds):
            firsts[seed + 1] += firsts[seed]
        var rows_by_seed = List[Int32](length=len(found_seed), fill=0)
        var ends_by_seed = List[Int32](length=len(found_end), fill=0)
        var costs_by_seed = List[Int32](length=len(found_cost), fill=0)
        var filled = firsts.copy()
        for index in range(len(found_seed)):
            var seed = Int(found_seed[index])
            var slot = Int(filled[seed])
            rows_by_seed[slot] = found_row[index]
            if inexact:
                ends_by_seed[slot] = found_end[index]
                costs_by_seed[slot] = found_cost[index]
            filled[seed] += 1

        # A match starts `cost` less its own cost layers above the best its end can chain on to. Right
        # to left on every diagonal, as the seeds go last first, so local pruning sees the matches
        # kept after it.
        var leftmost = List[Int32](length=self.columns + self.rows + 1, fill=Int32.MAX)
        var fronts = List[Int](length=4 * self.cost * LOOKAHEAD_SEEDS + 3, fill=0)
        for seed in range(self.seeds - 1, -1, -1):
            var column = seed * self.length
            var potential = self.potential(column)
            var end_potential = potential - self.cost
            for slot in range(Int(firsts[seed]), Int(firsts[seed + 1])):
                var start_row = Int(rows_by_seed[slot])
                var end_row = Int(ends_by_seed[slot]) if inexact else start_row + SEED_LENGTH
                var match_cost = Int(costs_by_seed[slot]) if inexact else 0
                if not self.worth_keeping(first, second, seed, start_row, end_row, match_cost, leftmost, fronts):
                    continue
                leftmost[column - start_row + self.rows] = Int32(column)
                var x = column - start_row - potential
                var y = start_row - column - potential
                var end_x = column + self.length - end_row - end_potential
                var end_y = end_row - column - self.length - end_potential
                var score = self.cost - match_cost
                var layer = self.score(end_x, end_y) + score
                while layer >= len(self.counts):
                    self.add_layer()
                for below in range(score):
                    # A start some start of the layer already dominates answers no query differently. In a
                    # repeat, the next seed's match a few rows on dominates nearly every match, and kept,
                    # they made each layer a long chain every query walked: poly-A of 16,384 bases
                    # against 8,192 took two seconds, 32,000 against 16,000 seventeen.
                    if not self.contains(layer - below, x, y):
                        self.add_point(layer - below, x, y)

    def count_seeds(mut self, symbols: List[UInt8]):
        """`remaining` and `counted`, a seed uncounted when it holds a code past `ACGT`."""
        self.remaining = List[Int32](length=self.seeds + 1, fill=0)
        var codes = symbols.unsafe_ptr()
        for seed in range(self.seeds - 1, -1, -1):
            var plain = True
            for offset in range(self.length):
                if codes[unsafe_offset=seed * self.length + offset] >= 4:
                    plain = False
                    break
            self.remaining[seed] = self.remaining[seed + 1] + Int32(1 if plain else 0)
        self.counted = Int(self.remaining[0])

    @inline(.always)
    def is_counted(self, seed: Int) -> Bool:
        """Whether the potential counts `seed` (see `remaining`)."""
        return len(self.remaining) == 0 or self.remaining[seed] != self.remaining[seed + 1]

    # Out of line: inlined into `build`, it slows the pruning loop there by a few percent.
    @inline(.never)
    def exact_matches(
        self,
        first: ImmPointer[UInt8, _],
        second: ImmPointer[UInt8, _],
        mut found_seed: List[Int32],
        mut found_row: List[Int32],
    ):
        """Every exact occurrence of a seed whose chain can still reach the end.

        Every seed's two-bit code is hashed by open addressing on the multiply's top bits; a slot
        holds its code and the first seed with it in one word, and seeds sharing a code chain on.
        The table is half full, so a window that matches nothing, nearly every one, would land on a
        taken slot half the time and probe on; a filter of `FILTER_BITS` a seed, the same multiply's
        top bits, turns almost all of them away first, on a branch that rarely goes the other way.
        """
        comptime MASK = (1 << (2 * SEED_LENGTH)) - 1
        comptime EMPTY = Int64(-1)
        comptime FILTER_BITS = 32
        var bits = 1
        while (1 << bits) < 2 * self.seeds:
            bits += 1
        var size = 1 << bits
        var filter_bits = 6
        while (1 << filter_bits) < FILTER_BITS * self.seeds:
            filter_bits += 1
        var filtered = List[UInt64](length=1 << (filter_bits - 6), fill=0)
        var filter = filtered.unsafe_ptr()
        var table = List[Int64](length=size, fill=EMPTY)
        var chained = List[Int32](length=self.seeds, fill=-1)
        var slots = table.unsafe_ptr()
        var chain = chained.unsafe_ptr()
        for seed in range(self.seeds):
            if not self.is_counted(seed):
                continue
            var code = 0
            for offset in range(SEED_LENGTH):
                code = (code << 2) | Int(first[unsafe_offset=seed * SEED_LENGTH + offset])
            var hash = UInt64(code) * 0x9E3779B97F4A7C15
            var bit = Int(hash >> UInt64(64 - filter_bits))
            filter[unsafe_offset=bit >> 6] |= UInt64(1) << UInt64(bit & 63)
            var slot = Int(hash >> UInt64(64 - bits))
            while slots[unsafe_offset=slot] != EMPTY and Int(slots[unsafe_offset=slot] >> 32) != code:
                slot = (slot + 1) & (size - 1)
            var held = slots[unsafe_offset=slot]
            chain[unsafe_offset=seed] = Int32(held & 0xFFFFFFFF) if held != EMPTY else -1
            slots[unsafe_offset=slot] = (Int64(code) << 32) | Int64(seed)

        # Every window of the second sequence looked up. A match's end, one seed on, is `(x + 1, y + 1)`
        # in transformed coordinates, and must lie at or below and left of the end's.
        var target_x = self.columns - self.rows
        var target_y = self.rows - self.columns
        var code = 0
        for row in range(self.rows):
            code = ((code << 2) | Int(second[unsafe_offset=row])) & MASK
            if row + 1 < SEED_LENGTH:
                continue
            var hash = UInt64(code) * 0x9E3779B97F4A7C15
            var bit = Int(hash >> UInt64(64 - filter_bits))
            if (filter[unsafe_offset=bit >> 6] >> UInt64(bit & 63)) & 1 == 0:
                continue
            var slot = Int(hash >> UInt64(64 - bits))
            var held = slots[unsafe_offset=slot]
            while held != EMPTY and Int(held >> 32) != code:
                slot = (slot + 1) & (size - 1)
                held = slots[unsafe_offset=slot]
            if held == EMPTY:
                continue
            var start_row = row + 1 - SEED_LENGTH
            var seed = Int(held & 0xFFFFFFFF)
            while seed >= 0:
                var column = seed * SEED_LENGTH
                # The potential with every seed counted, the most it can be: with some uncounted this
                # keeps more matches than it must, which only weaken the bound.
                var potential = self.seeds - seed
                if column - start_row - potential + 1 <= target_x and start_row - column - potential + 1 <= target_y:
                    found_seed.append(Int32(seed))
                    found_row.append(Int32(start_row))
                seed = Int(chain[unsafe_offset=seed])

    # Out of line: inlined into `build`, it slows the pruning loop there by a few percent.
    @inline(.never)
    def inexact_matches(
        self,
        first: ImmPointer[UInt8, _],
        second: ImmPointer[UInt8, _],
        mut found_seed: List[Int32],
        mut found_row: List[Int32],
        mut found_end: List[Int32],
        mut found_cost: List[Int32],
    ):
        """Every occurrence of a seed within one edit whose chain can still reach the end.

        One edit leaves one half of the seed matching exactly, at the window's start or its end. Each
        half's code indexes a table of the seeds with it, and every window of the second sequence is
        looked up as either half: a left half found at a row tries the windows starting there, one
        base shorter, as long, or one longer than the seed; a right half found tries those ending
        where it ends. A window both halves find is taken from the left half alone.

        Within a bucket the seeds go in order, and only a range of them can be found at any row with
        their chain reaching the end (see `reachable_seeds`), so each lookup tries that range alone.
        """
        comptime HALF = INEXACT_LENGTH // 2
        comptime HALF_MASK = (1 << HALF_BITS) - 1
        comptime BUCKETS = 1 << HALF_BITS
        var codes = List[UInt64](length=self.seeds, fill=0)
        # Each half's seeds bucketed by its code, contiguous, in order, each with its code below it in
        # one word (see `ENTRY_BITS`).
        var left_start = List[Int32](length=BUCKETS + 1, fill=0)
        var right_start = List[Int32](length=BUCKETS + 1, fill=0)
        for seed in range(self.seeds):
            var code = UInt64(0)
            for offset in range(INEXACT_LENGTH):
                code = (code << 2) | UInt64(first[unsafe_offset=seed * INEXACT_LENGTH + offset])
            codes[seed] = code
            if not self.is_counted(seed):
                continue
            left_start[Int(code >> UInt64(HALF_BITS)) + 1] += 1
            right_start[(Int(code) & HALF_MASK) + 1] += 1
        for bucket in range(BUCKETS):
            left_start[bucket + 1] += left_start[bucket]
            right_start[bucket + 1] += right_start[bucket]
        var left_entries = List[UInt64](length=self.seeds, fill=0)
        var right_entries = List[UInt64](length=self.seeds, fill=0)
        var left_fill = left_start.copy()
        var right_fill = right_start.copy()
        for seed in range(self.seeds):
            if not self.is_counted(seed):
                continue
            var code = codes[seed]
            var entry = (UInt64(seed) << UInt64(ENTRY_BITS)) | code
            left_entries[Int(left_fill[Int(code >> UInt64(HALF_BITS))])] = entry
            left_fill[Int(code >> UInt64(HALF_BITS))] += 1
            right_entries[Int(right_fill[Int(code) & HALF_MASK])] = entry
            right_fill[Int(code) & HALF_MASK] += 1

        # Every window of `INEXACT_LENGTH` bases of the second sequence as one code, the first base in
        # the high bits; past the end, the bases read as zero, and no window reaching there is tried.
        var windows = List[UInt64](length=self.rows + 1, fill=0)
        var rolling = UInt64(0)
        for row in range(self.rows + INEXACT_LENGTH - 1, -1, -1):
            var base = UInt64(second[unsafe_offset=row]) if row < self.rows else UInt64(0)
            rolling = (rolling >> 2) | (base << UInt64(2 * INEXACT_LENGTH - 2))
            if row <= self.rows:
                windows[row] = rolling
        self.inexact_scan(
            left_start.unsafe_ptr(),
            right_start.unsafe_ptr(),
            left_entries.unsafe_ptr(),
            right_entries.unsafe_ptr(),
            windows.unsafe_ptr(),
            found_seed,
            found_row,
            found_end,
            found_cost,
        )

    def inexact_scan(
        self,
        left_start: ImmPointer[Int32, _],
        right_start: ImmPointer[Int32, _],
        left_entries: ImmPointer[UInt64, _],
        right_entries: ImmPointer[UInt64, _],
        window: ImmPointer[UInt64, _],
        mut found_seed: List[Int32],
        mut found_row: List[Int32],
        mut found_end: List[Int32],
        mut found_cost: List[Int32],
    ):
        """The inexact matches every row of the second sequence finds, from the half tables (each
        half's bucket bounds and entries) and the windows.

        Every seed a half finds is tested against the three windows it could match without a branch
        (see `within_one_edit`), most failing all three, and only one that passes is tried in full.

        The tables are larger than a core's own cache on some machines, and every lookup goes from
        a bucket's bounds to its entries, each load a miss waiting on the last. So the rows go in
        batches of `SCAN_BATCH`: first every row's bounds, the loads independent and overlapping,
        their entries fetched ahead, and then the tests, the entries on their way or arrived.
        """
        comptime HALF = INEXACT_LENGTH // 2
        comptime K = INEXACT_LENGTH
        comptime CODE = (UInt64(1) << UInt64(ENTRY_BITS)) - 1
        comptime AHEAD = PrefetchOptions().for_read().high_locality()
        var batch = List[Int32](length=4 * SCAN_BATCH, fill=0)
        var bounds = batch.unsafe_ptr()
        # A left half needs half a seed after it.
        var end_row = self.rows - HALF + 1
        for batch_row in range(0, end_row, SCAN_BATCH):
            var batch_end = min(batch_row + SCAN_BATCH, end_row)
            for row in range(batch_row, batch_end):
                var half = Int(window[unsafe_offset=row] >> UInt64(HALF_BITS))
                var at = 4 * (row - batch_row)
                bounds[unsafe_offset=at] = left_start[unsafe_offset=half]
                bounds[unsafe_offset=at + 1] = left_start[unsafe_offset=half + 1]
                bounds[unsafe_offset=at + 2] = right_start[unsafe_offset=half]
                bounds[unsafe_offset=at + 3] = right_start[unsafe_offset=half + 1]
                prefetch[AHEAD](left_entries.unsafe_offset(Int(bounds[unsafe_offset=at])))
                prefetch[AHEAD](right_entries.unsafe_offset(Int(bounds[unsafe_offset=at + 2])))
            for row in range(batch_row, batch_end):
                var at = 4 * (row - batch_row)
                # A left half found here: the windows start at `row`, one base shorter than the seed,
                # as long, and one longer, whose last `K` bases start a row on.
                var head = window[unsafe_offset=row]
                var tail = window[unsafe_offset=row + 1]
                var starting = self.reachable_seeds(row + K - 1, row + K + 1)
                var after = Int(bounds[unsafe_offset=at + 1])
                for slot in range(
                    first_reachable(left_entries, Int(bounds[unsafe_offset=at]), after, starting[0]), after
                ):
                    var seed = Int(left_entries[unsafe_offset=slot] >> UInt64(ENTRY_BITS))
                    if seed > starting[1]:
                        break
                    var code = left_entries[unsafe_offset=slot] & CODE
                    if within_one_edit(head >> 2, head, head, tail, code):
                        self.try_windows(code, seed, row, -1, window, found_seed, found_row, found_end, found_cost)
                # A right half ending at `end`: the windows end there, so start a base later, at the
                # same row, or a base earlier than the seed's length back.
                var end = row + HALF
                var shorter = window[unsafe_offset=max(end - K + 1, 0)]
                var level = window[unsafe_offset=max(end - K, 0)]
                var longer = window[unsafe_offset=max(end - K - 1, 0)]
                var ending = self.reachable_seeds(end, end)
                after = Int(bounds[unsafe_offset=at + 3])
                for slot in range(
                    first_reachable(right_entries, Int(bounds[unsafe_offset=at + 2]), after, ending[0]), after
                ):
                    var seed = Int(right_entries[unsafe_offset=slot] >> UInt64(ENTRY_BITS))
                    if seed > ending[1]:
                        break
                    var code = right_entries[unsafe_offset=slot] & CODE
                    if ends_within_one_edit(shorter, level, longer, code):
                        self.try_windows(code, seed, -1, end, window, found_seed, found_row, found_end, found_cost)

    @inline(.always)
    def try_windows(
        self,
        code: UInt64,
        seed: Int,
        start: Int,
        end: Int,
        window: ImmPointer[UInt64, _],
        mut found_seed: List[Int32],
        mut found_row: List[Int32],
        mut found_end: List[Int32],
        mut found_cost: List[Int32],
    ):
        """The windows within one edit of a seed that start at `start`, or end at `end` when `start`
        is -1: one base shorter than the seed, as long, and one longer.

        Two codes of equal length agree base by base where their XOR has neither bit of a pair set;
        the leading agreeing bases are the common prefix, the trailing ones the common suffix. A
        window one base shorter is the seed less one base when the two cover it between them, and
        one longer is the seed plus one base when they cover the seed.

        An exact match's windows one base shorter or longer, sharing its start or its end, are left
        out: a path crossing the seed along one costs at least one edit there, and taking the exact
        match in its place moves the chain's gaps by one edit at most, so the heuristic stays a
        lower bound without them.
        """
        comptime K = INEXACT_LENGTH
        var exact = False
        for step in range(3):
            # As long as the seed first, so an exact match is known before its neighbours.
            var extra = 0 if step == 0 else (-1 if step == 1 else 1)
            if exact and extra != 0:
                return
            var size = K + extra
            var low = start if start >= 0 else end - size
            var high = low + size
            if low < 0 or high > self.rows:
                continue
            # A window whose left half matches too is the left half's lookup's to take.
            var taken = start < 0 and (window[unsafe_offset=low] >> UInt64(HALF_BITS)) == (code >> UInt64(HALF_BITS))
            var cost = 1
            if extra == 0:
                var bases = Int(pop_count(differing(window[unsafe_offset=low], code)))
                exact = bases == 0
                if bases > 1 or taken:
                    continue
                cost = bases
            elif taken:
                continue
            elif extra == -1:
                # The window's bases, `K - 1` of them, against the seed's first and last `K - 1`.
                var shorter = window[unsafe_offset=low] >> 2
                var prefix = self.agreeing_prefix(shorter, code >> 2, K - 1)
                var suffix = self.agreeing_suffix(shorter, code & ((UInt64(1) << UInt64(2 * K - 2)) - 1), K - 1)
                if prefix + suffix < K - 1:
                    continue
            else:
                # The window's first and last `K` bases against the seed.
                var head = window[unsafe_offset=low]
                var tail = window[unsafe_offset=low + 1]
                var prefix = self.agreeing_prefix(head, code, K)
                var suffix = self.agreeing_suffix(tail, code, K)
                if prefix + suffix < K:
                    continue
            self.add_match(seed, low, high, cost, found_seed, found_row, found_end, found_cost)

    @inline(.always)
    def agreeing_prefix(self, first: UInt64, second: UInt64, bases: Int) -> Int:
        """Leading bases two codes of `bases` bases share."""
        var bits = differing(first, second)
        if bits == 0:
            return bases
        return (Int(count_leading_zeros(bits)) - (64 - 2 * bases)) // 2

    @inline(.always)
    def agreeing_suffix(self, first: UInt64, second: UInt64, bases: Int) -> Int:
        """Trailing bases two codes of `bases` bases share."""
        var bits = differing(first, second)
        if bits == 0:
            return bases
        return Int(count_trailing_zeros(bits)) // 2

    @inline(.always)
    def reachable_seeds(self, lowest_end: Int, highest_end: Int) -> Tuple[Int, Int]:
        """The first and last inexact seed a match ending in rows `[lowest_end, highest_end]` can be
        of and still pass `add_match`.

        Seed `s` at column `16s` charges the seeds after it `2(S - s - 1)`, so its match ending at row
        `e` passes when `16s + 16 - e - 2(S - s - 1) <= C - R` and `e - 16s - 16 - 2(S - s - 1) <= R - C`,
        `S` the seeds, `C` and `R` the columns and rows: `18s <= C - R + e + 2S - 18` and
        `14s >= e + C - R - 2S - 14`, a range of seeds, here taken at its widest over the ends. With
        seeds uncounted the potential is lower and the range narrower, so this one holds it.
        """
        comptime assert INEXACT_LENGTH == 16, "the bounds below are worked out for 16-base seeds"
        var shift = self.columns - self.rows
        var first = ceildiv(lowest_end + shift - 2 * self.seeds - 14, 14)
        var last = (highest_end + shift + 2 * self.seeds - 18) // 18
        return (max(first, 0), min(last, self.seeds - 1))

    @inline(.always)
    def add_match(
        self,
        seed: Int,
        start_row: Int,
        end_row: Int,
        cost: Int,
        mut found_seed: List[Int32],
        mut found_row: List[Int32],
        mut found_end: List[Int32],
        mut found_cost: List[Int32],
    ):
        """Keeps an inexact match whose end can still reach the end: in transformed coordinates, its
        end lies at or below and left of the end's."""
        var column = seed * self.length
        var end_potential = self.potential(column) - self.cost
        var end_x = column + self.length - end_row - end_potential
        var end_y = end_row - column - self.length - end_potential
        if end_x <= self.columns - self.rows and end_y <= self.rows - self.columns:
            found_seed.append(Int32(seed))
            found_row.append(Int32(start_row))
            found_end.append(Int32(end_row))
            found_cost.append(Int32(cost))

    # Out of line: inlined into `build`, it slows the pruning loop there by a few percent.
    @inline(.never)
    def worth_keeping(
        self,
        first: ImmPointer[UInt8, _],
        second: ImmPointer[UInt8, _],
        seed: Int,
        start_row: Int,
        end_row: Int,
        match_cost: Int,
        leftmost: List[Int32],
        mut fronts: List[Int],
    ) -> Bool:
        """A*PA's local pruning: whether a match can lower the cost of some path, so dropping it would not.

        A match only helps a path that goes on to cross the next seeds for fewer edits than the
        seeds it crosses would cost unmatched. A diagonal transition from the match's end, starting
        at the match's own cost, over the next `LOOKAHEAD_SEEDS` seeds counting its own, looks for such a
        path: the match stays if a front reaches past the last of them, or slides into a match
        already kept on its diagonal, which then continues it.
        A front whose edits already equal the seeds crossed can no longer gain, and is dropped. Every
        match an optimal path's chain relies on passes, so the heuristic stays a lower bound.
        """
        var start_column = seed * self.length
        var start_potential = self.potential(start_column)
        var last = min(seed + LOOKAHEAD_SEEDS - 1, self.seeds - 1)
        var end_column = (last + 1) * self.length
        var reach = start_potential - self.potential(end_column)
        var origin = (start_column + self.length) - end_row
        # Front `d` is the diagonal `origin + d - reach`, at `fronts[d]`, its furthest column.
        var low = reach
        var high = reach + 1
        var front = fronts.unsafe_ptr()
        # A run past the last seed answers as one reaching it: in a repeat a run to its end took each match
        # thousands of bases.
        var horizon = min(end_column, self.columns)
        var reached = extend(first, second, start_column + self.length, end_row, horizon, self.rows)
        front[unsafe_offset=reach] = reached
        if reached >= end_column:
            return True
        var nearest = leftmost.unsafe_ptr()
        var kept = Int(nearest[unsafe_offset=origin + self.rows])
        if kept <= reached:
            return True
        # Past the live fronts, one on either side reads as unreachable, so every front takes the best
        # of its three sources with no test of which exist.
        comptime UNREACHED = -(1 << 40)
        for cost in range(match_cost + 1, reach):
            # One more edit: from the same diagonal, or from either neighbour, in place, the front
            # below's value before this edit carried along.
            front[unsafe_offset=low - 1] = UNREACHED
            front[unsafe_offset=high] = UNREACHED
            front[unsafe_offset=high + 1] = UNREACHED
            var below = UNREACHED
            for d in range(low - 1, high + 1):
                var current = front[unsafe_offset=d]
                front[unsafe_offset=d] = min(max(max(current, below) + 1, front[unsafe_offset=d + 1]), self.columns)
                below = current
            low -= 1
            high += 1
            # A front whose edits match the seeds it has crossed can no longer gain.
            while low < high and cost + self.potential(front[unsafe_offset=low]) >= start_potential:
                low += 1
            while high > low and cost + self.potential(front[unsafe_offset=high - 1]) >= start_potential:
                high -= 1
            if low == high:
                return False
            for d in range(low, high):
                var diagonal = origin + d - reach
                var before = front[unsafe_offset=d]
                var row = before - diagonal
                if row < 0 or row > self.rows:
                    continue
                var after = extend(first, second, before, row, horizon, self.rows)
                front[unsafe_offset=d] = after
                if after >= end_column:
                    return True
                # A point inside the matrix, so its diagonal indexes `leftmost`.
                var kept_column = Int(nearest[unsafe_offset=diagonal + self.rows])
                if before <= kept_column and kept_column <= after:
                    return True
        return False

    def add_layer(mut self):
        """One more layer, its slots empty and nothing spilled."""
        for _ in range(LAYER_SLOTS):
            self.slot_x.append(0)
            self.slot_y.append(0)
        self.counts.append(0)
        self.spill_x.append(List[Int32]())
        self.spill_y.append(List[Int32]())

    def add_point(mut self, layer: Int, x: Int, y: Int):
        """Adds the transformed start `(x, y)` to `layer`: into a slot while one is free, else onto the
        layer's staircase, unless a start there lies at or above and right of it, and in place of those it
        does."""
        var count = Int(self.counts[layer])
        if count < LAYER_SLOTS:
            self.slot_x[layer * LAYER_SLOTS + count] = Int32(x)
            self.slot_y[layer * LAYER_SLOTS + count] = Int32(y)
            self.counts[layer] = Int32(count + 1)
            return
        ref xs = self.spill_x[layer]
        ref ys = self.spill_y[layer]
        var place = first_at_least(xs, x)
        if place < len(xs) and Int(ys[place]) >= y:
            return
        # The starts this one covers: those left of it and no higher, a run just before its place as `y`
        # falls along the staircase, and one at its own `x` below it.
        var end = place + 1 if place < len(xs) and Int(xs[place]) == x else place
        var begin = place
        while begin > 0 and Int(ys[begin - 1]) <= y:
            begin -= 1
        if begin == end:
            xs.insert(begin, Int32(x))
            ys.insert(begin, Int32(y))
            return
        xs[begin] = Int32(x)
        ys[begin] = Int32(y)
        var gone = end - begin - 1
        if gone > 0:
            for index in range(end, len(xs)):
                xs[index - gone] = xs[index]
                ys[index - gone] = ys[index]
            xs.resize(len(xs) - gone, 0)
            ys.resize(len(ys) - gone, 0)

    def add_sentinel(mut self):
        """Layer zero: a point dominating everything, as no match is chained."""
        self.add_layer()
        self.add_point(0, Int(Int32.MAX), Int(Int32.MAX))

    @inline(.always)
    def contains(self, layer: Int, x: Int, y: Int) -> Bool:
        """Whether a start in `layer` lies at or above and right of `(x, y)`."""
        var count = Int(self.counts.unsafe_ptr()[unsafe_offset=layer])
        var xs = self.slot_x.unsafe_ptr().unsafe_offset(layer * LAYER_SLOTS)
        var ys = self.slot_y.unsafe_ptr().unsafe_offset(layer * LAYER_SLOTS)
        for index in range(count):
            if x <= Int(xs[unsafe_offset=index]) and y <= Int(ys[unsafe_offset=index]):
                return True
        if count == LAYER_SLOTS:
            ref step_x = self.spill_x.unsafe_ptr()[unsafe_offset=layer]
            ref step_y = self.spill_y.unsafe_ptr()[unsafe_offset=layer]
            # Of the starts right of `x`, the leftmost lies highest.
            var place = first_at_least(step_x, x)
            if place < len(step_x) and y <= Int(step_y.unsafe_ptr()[unsafe_offset=place]):
                return True
        return False

    def score(mut self, x: Int, y: Int) -> Int:
        """The most matches a chain from transformed point `(x, y)` takes.

        Layers nest, every start of layer `v + 1` lying below some start of layer `v`, so the
        layers containing a point are a prefix: a search galloping out from the last answer, as
        neighbouring queries land on neighbouring layers, finds its end.
        """
        var last = len(self.counts) - 1
        var guess = min(self.hint, last)
        var low: Int
        var high: Int
        if self.contains(guess, x, y):
            # The end lies at or above the guess.
            low = guess
            var step = 1
            high = min(guess + step, last)
            while high > low and self.contains(high, x, y):
                low = high
                step *= 2
                high = min(low + step, last)
            if high == low:
                self.hint = low
                return low
            high -= 1
        else:
            # The end lies below the guess; layer zero always contains.
            high = guess - 1
            var step = 1
            low = max(guess - step, 0)
            while low > 0 and not self.contains(low, x, y):
                high = low - 1
                step *= 2
                low = max(high - step, 0)
        while low < high:
            var middle = (low + high + 1) // 2
            if self.contains(middle, x, y):
                low = middle
            else:
                high = middle - 1
        self.hint = low
        return low

    def chains_well(mut self) -> Bool:
        """Whether at least one seed in `CHAINED_SHARE` is chained from the origin: then the bound
        there lies close to the distance, rather than at little more than an edit a seed."""
        if self.seeds == 0:
            return False
        var whole = self.potential(0)
        return (whole - self.h(0, 0)) * CHAINED_SHARE >= whole

    @inline(.always)
    def potential(self, column: Int) -> Int:
        """What the seeds starting at or after `column` cost a path matching none of them."""
        # Each length its own constant divisor: the band asks for this at every row it prunes.
        if likely(len(self.remaining) == 0):
            if self.cost == 1:
                return self.seeds - min(self.seeds, ceildiv(column, SEED_LENGTH))
            return 2 * (self.seeds - min(self.seeds, ceildiv(column, INEXACT_LENGTH)))
        return self.uncounted_potential(column)

    @inline(.always)
    def uncounted_potential(self, column: Int) -> Int:
        """`potential` when some seeds go uncounted (see `remaining`), its divisors constant too: a
        64-bit division by a variable takes tens of cycles on x86."""
        var counted = self.remaining.unsafe_ptr()
        if self.cost == 1:
            return Int(counted[unsafe_offset=min(self.seeds, ceildiv(column, SEED_LENGTH))])
        return 2 * Int(counted[unsafe_offset=min(self.seeds, ceildiv(column, INEXACT_LENGTH))])

    @inline(.always)
    def climb(self) -> Int:
        """The most the heuristic changes a row down or up a column: one with no seeds or exact ones,
        `cost` with inexact. A best chain from a row still starts, but for its first match, from the
        row above or below, as every match moves the transformed point at least one up and one
        right; the first match scores at most `cost`.

        Down a column a complete set of matches would hold the drop to one, the heuristic being
        consistent, but an exact match's one-edit neighbours are left out (see `try_windows`), and
        without them only this bound holds."""
        return self.cost

    def h(mut self, column: Int, row: Int) -> Int:
        """A lower bound on the cost from `(column, row)` to the end."""
        var gap = abs((self.columns - column) - (self.rows - row))
        if self.seeds == 0:
            return gap
        var potential = self.potential(column)
        var chained = self.score(column - row - potential, row - column - potential)
        if chained == 0:
            return max(gap, potential)
        return potential - chained

    @inline(.always)
    def bound[seeded: Bool](mut self, column: Int, row: Int) -> Int:
        """`h` in a band known to have seeds, or not: without, the gap alone, so the band's loops carry
        none of the seeds' code, which on the Skylake-X slowed pairs without seeds by 3 to 7%."""
        comptime if seeded:
            return self.h(column, row)
        return abs((self.columns - column) - (self.rows - row))
