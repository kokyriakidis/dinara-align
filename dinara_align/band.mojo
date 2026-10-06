# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Ported from `pa-bitpacking` in A*PA (https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner,
# commit bf2e14e), by Ragnar Groot Koerkamp and Pesho Ivanov, itself translated from Edlib.
"""
Band doubling, as in A*PA2-simple: guess a bound, compute only the cells a path within it could cross,
prune rows whose score already rules them out, and repeat with a larger guess until the distance fits
under it (see `pruned_distance`). Unlike A*PA2, the next guess is aimed from where the failed round died
rather than doubled (see `band_doubling`).
"""

from std.math import ceildiv

from .bit_parallel import (
    BAND_COLUMNS,
    Edge,
    Frontier,
    LANES,
    NARROW_COLUMNS,
    Profile,
    Sweep,
    tile_bounds,
    Trail,
    WORD_BITS,
    word_value,
)
from .diagonal import (
    diagonal_transition,
    DiagonalFronts,
    Probe,
    PROBE_MARGIN,
    PROJECTION_ONLY,
    SHORT_COLUMNS,
    trusted_projection,
)
from .seeds import (
    INEXACT_CHAINED,
    INEXACT_CHAINED_SPREAD,
    INEXACT_COLUMNS,
    INEXACT_DIVERGENCE,
    SEED_COLUMNS,
    SEED_DIVERGENCE,
    SEED_EDITS,
    SEED_SLACK,
    SEEDED_GROWTH,
    SEEDED_TRUST_SHARE,
    SeedHeuristic,
    SHORT_SEEDS,
)


comptime NARROW_BAND = 1024
"""The bound below which band tiles are `NARROW_COLUMNS` wide rather than `BAND_COLUMNS`."""


@fieldwise_init
struct Round(ImplicitlyCopyable, TrivialRegisterPassable):
    """What one pruned round learned: a distance, or how far across the matrix any path within the bound got."""

    var distance: Int
    """An upper bound on the distance, exact when within the round's bound; -1 when no path fit it."""
    var reached: Int
    """Matrix columns crossed before every row was pruned, or all of them when a distance was found."""
    var bound: Int
    """The bound the round ended on, which its checkpoints may have lowered (see `Band.check`)."""
    var estimate: Int
    """The distance a checkpoint projected when it gave the round up, or -1."""


struct Band(Movable):
    """One round of band doubling, advanced a tile at a time.

    `pruned_distance` documents the round; this holds its state between tiles: the frontier, the
    band's top and bottom words, the score at the top, and the deepest kept row and its score.
    """

    var frontier: Frontier
    var sweep: Sweep
    var bounds: List[Int]
    var columns: Int
    var rows: Int
    var words: Int
    var difference: Int
    var extra: Int
    var threshold: Int
    var top: Int
    var end_word: Int
    var anchor: Int
    """The score at the top of word `top`, on the left edge of the coming tile."""
    var deepest: Int
    """The deepest row a kept cell, scoring `floor` at the least, sits on, after the last tile."""
    var floor: Int
    var outcome: Round
    """Set once the round has pruned every row, the reason it stopped."""
    var adapt: Bool
    """Whether the round re-aims its bound at checkpoints (see `check`)."""
    var edge: Edge
    """The scores down the last tile's right edge."""
    var checkpoint: Int
    """Checkpoints passed."""

    def __init__(
        out self,
        mut profile: Profile,
        threshold: Int,
        adapt: Bool,
    ):
        self.columns = profile.columns
        self.rows = profile.rows
        self.words = profile.words
        self.difference = self.rows - self.columns
        self.extra = (threshold - abs(self.difference)) // 2
        self.threshold = threshold
        # One tile's width of horizontal edge, since each tile starts again from `+1` above (see
        # `Sweep.shifted`); the last tile may absorb a sliver of up to `2 * LANES` more columns.
        self.frontier = Frontier(BAND_COLUMNS + 2 * LANES, self.words)
        profile.build_planes()
        self.sweep = self.frontier.sweep(profile)
        self.bounds = tile_bounds(self.columns, NARROW_COLUMNS if threshold < NARROW_BAND else BAND_COLUMNS)
        self.top = 0
        self.end_word = 0
        self.anchor = 0
        self.deepest = 0
        self.floor = 0
        self.outcome = Round(-1, -1, threshold, -1)
        self.adapt = adapt
        self.edge = Edge()
        self.checkpoint = 0

    def tiles(self) -> Int:
        return len(self.bounds) - 1

    def first_column(self, tile: Int) -> Int:
        return self.bounds[tile]

    def end_column(self, tile: Int) -> Int:
        return self.bounds[tile + 1]

    def prepare[
        record: Bool, seeded: Bool
    ](mut self, tile: Int, mut trail: Trail, mut heuristic: SeedHeuristic) -> Bool:
        """Sets the tile's words and readies its edge; false, with `outcome` set, when no word is left."""
        var first_column = self.bounds[tile]
        var end_column = self.bounds[tile + 1]
        var width = end_column - first_column
        # Ukkonen's band bounds the bottom, and so does the deepest kept row: going `k` rows past it
        # costs at least `k` over the score there, at least `floor`, and below the diagonal that leads
        # to the end each such row also adds one to the gap still to close.
        var band_row = end_column + max(0, self.difference) + self.extra
        var slack = self.threshold - self.floor
        var end_diagonal = self.rows - (self.columns - end_column)
        var reach_row = min(
            band_row,
            self.deepest + width + slack,
            max(end_diagonal, (slack + self.deepest + width + end_diagonal) // 2),
            self.rows,
        )
        comptime if seeded:
            # The seeds bound the bottom tighter, as A*PA2 bounds a block's: a row `k` past the
            # diagonal from the deepest kept cell costs at least `floor + k` to reach, so it is in
            # reach only while that plus the heuristic there fits the bound. Going up, the first drops
            # by one a row and the heuristic by at most its `climb`, so a row `x` over rules out the
            # next `ceil(x / (1 + climb))` above it unread.
            var diagonal_row = self.deepest + width
            var step = 1 + heuristic.climb()
            while reach_row > diagonal_row:
                var over = (
                    self.floor
                    + (reach_row - diagonal_row)
                    + heuristic.bound[seeded](end_column, reach_row)
                    - self.threshold
                )
                if over <= 0:
                    break
                reach_row = max(reach_row - ceildiv(over, step), diagonal_row)
        var reach = ceildiv(reach_row, WORD_BITS)
        # Exactly the words the band reaches: `words` has a kernel for every count. The bottom never
        # rises, which would leave words behind whose differences the next tile still reads.
        self.end_word = max(self.end_word, min(reach, self.words))
        if self.end_word <= self.top:
            self.outcome = Round(-1, first_column, self.threshold, -1)
            return False
        comptime if record:
            trail.record(first_column, end_column, self.top, self.end_word, self.anchor, self.frontier)
        self.frontier.restart_horizontal(width)
        return True

    def tile_sweep(self, tile: Int) -> Sweep:
        """The sweep for one tile, its horizontal edge starting at the tile's first column."""
        return self.sweep.shifted(self.bounds[tile])

    def finish[seeded: Bool](mut self, tile: Int, mut heuristic: SeedHeuristic) -> Bool:
        """Prunes after a swept tile; false, with `outcome` set, when every row went."""
        var end_column = self.bounds[tile + 1]
        self.anchor += end_column - self.bounds[tile]
        # Read the right edge back as scores and keep, as A*PA2 does, the rows from the first to the
        # last whose score plus heuristic fits the bound. Scores change by at most one per row and the
        # heuristic by at most its `climb` either way, so a row `x` over the bound rules out the next
        # `ceil(x / (1 + climb))` rows on either side without reading them: half of `x` for exact
        # seeds or none, a third for inexact.
        self.edge.capture(self.top, self.end_word, self.anchor, self.frontier, self.rows)
        var first_kept = self.edge.low_row
        var last_kept = self.edge.high_row
        var step = 1 + heuristic.climb()
        while first_kept <= last_kept:
            var over = self.edge.score(first_kept) + heuristic.bound[seeded](end_column, first_kept) - self.threshold
            if over <= 0:
                break
            first_kept += ceildiv(over, step)
        while last_kept >= first_kept:
            var over = self.edge.score(last_kept) + heuristic.bound[seeded](end_column, last_kept) - self.threshold
            if over <= 0:
                break
            last_kept -= ceildiv(over, step)
        if first_kept > last_kept:
            self.outcome = Round(-1, end_column, self.threshold, -1)
            return False
        if self.adapt and not self.check[seeded](end_column, first_kept, last_kept, heuristic):
            return False

        # Move the top down to the word whose top row is at or above the first kept row, carrying the
        # anchor past the words it leaves.
        var new_top = first_kept // WORD_BITS
        for word in range(self.top, new_top):
            self.anchor += word_value(self.frontier.vertical_plus[word], self.frontier.vertical_minus[word])
        self.top = new_top
        # The bottom-most kept cell bounds every path below it: from there each row past the diagonal
        # costs one more, so its exact score is the floor the next tile's reach is measured from.
        self.deepest = last_kept
        self.floor = self.edge.score(last_kept)
        return True

    def check[
        seeded: Bool
    ](mut self, end_column: Int, first_kept: Int, last_kept: Int, mut heuristic: SeedHeuristic) -> Bool:
        """Re-aims the bound from the band's own climb at a checkpoint; false, with `outcome` set, to
        give the round up.

        At an eighth, a quarter and half of the columns, the least score plus heuristic down the kept
        rows, the gap to the end or the seeds still ahead, has climbed from the heuristic at the origin
        about in proportion to the columns crossed, so scaled to the whole width it projects the
        distance from hundreds of edits rather than the diagonal transition's handful: within about a
        tenth at an eighth, closer further on. The bound drops to the projection plus `CHECK_MARGIN`
        tenths for each checkpoint still ahead, and a round whose projection less half as much passes
        its bound gives up, its rest likely wasted. Lowering the bound keeps the round exact: on an
        optimal path a cell's score plus its gap to the end is at most the distance, so a distance
        within the final bound passed every column unpruned.
        """
        var passed = self.checkpoint
        while passed < CHECKPOINTS and end_column >= self.columns >> (CHECKPOINTS - passed):
            passed += 1
        if passed == self.checkpoint or end_column >= self.columns:
            return True
        self.checkpoint = passed
        # Sampled every `CHECK_ROWS` rows: the sum changes by at most two a row, a few edits at most.
        var least = Int.MAX
        var row = first_kept
        while True:
            least = min(least, self.edge.score(row) + heuristic.bound[seeded](end_column, row))
            if row == last_kept:
                break
            row = min(row + CHECK_ROWS, last_kept)
        var gap = abs(self.difference)
        var origin = heuristic.bound[seeded](0, 0)
        var estimate = origin + max(least - origin, 0) * self.columns // end_column
        var margin = estimate * (CHECKPOINTS + 1 - passed) * CHECK_MARGIN // 10
        if estimate - margin // 2 > self.threshold:
            self.outcome = Round(-1, end_column, self.threshold, estimate)
            return False
        var aim = estimate + margin + PROBE_MARGIN
        if aim < self.threshold:
            self.threshold = aim
            self.extra = (aim - gap) // 2
        return True

    def result(mut self) -> Round:
        """What the round found once every tile is swept, with the last edge captured."""
        self.edge.capture(self.top, self.end_word, self.anchor, self.frontier, self.rows)
        if self.end_word < self.words:
            return Round(-1, self.columns, self.threshold, -1)
        return Round(self.edge.score(self.rows), self.columns, self.threshold, -1)


def pruned_distance[
    record: Bool, seeded: Bool
](mut profile: Profile, threshold: Int, mut trail: Trail, mut heuristic: SeedHeuristic, adapt: Bool,) -> Round:
    """One round of band doubling with A*PA2-simple's pruning.

    Only cells a path of cost at most `threshold` could cross are computed. Ukkonen's band bounds
    them, gap from the start plus gap to the end within the bound, and pruning narrows it further:
    after each tile the absolute scores down its right edge are read back, and a word is dropped for
    good once even its lowest score plus its gap to the end exceeds `threshold`, since no path within
    the bound can pass through it, or below it from above. The bottom grows only to the deepest row
    the lowest kept score could still reach within the bound.

    Every cell outside what is computed reads as the cost of a real path (see `Frontier`), so no
    computed score is ever below the true one, and the scores along an optimal path within the
    bound are exact by induction: its predecessors are exact and kept, so its own rows are never the
    ones dropped. Hence a distance at most `threshold` is the distance.

    The band's top and bottom only move down, so no word reads differences left from an earlier
    column, and each word's differences end as the right edge of the last tile that swept it. With
    the top word reading `+1` from above in every column, the score at the corner reads back as for
    the whole matrix: along the top row, then down the right edge. Each tile sweeps exactly the
    words the band reaches.

    With `record`, every tile's left edge goes into `trail` before the tile is swept, for the
    traceback; without, the trail is untouched. Each tile matches on as many planes as its own symbols
    need (see `Profile.symbols`). `seeded` says whether `heuristic` has seeds (see `SeedHeuristic.bound`).

    With `adapt`, the round re-aims its bound as it goes (see `Band.check`), and a distance is exact
    only within the bound it ends on.
    """
    comptime if record:
        trail.clear()
    var band = Band(profile, threshold, adapt)
    for tile in range(band.tiles()):
        if not band.prepare[record, seeded](tile, trail, heuristic):
            return band.outcome
        band.tile_sweep(tile).words(
            profile.symbols(band.first_column(tile), band.end_column(tile)),
            band.top,
            band.end_word,
            band.first_column(tile),
            band.end_column(tile),
        )
        if not band.finish[seeded](tile, heuristic):
            return band.outcome
    return band.result()


comptime SHORT_AIM = 17
"""A short pair's first bound, in tenths of the projection."""


comptime SHORT_REACH = 128
"""The most a short pair's first bound reaches past the projection, so a projection already too
high, as on very divergent pairs, does not widen the band by most of the matrix."""


comptime CHECKPOINTS = 3
"""A first round without seeds re-aims its bound at `columns >> k` for `k = CHECKPOINTS ..= 1`: an
eighth, a quarter and half of the way across (see `Band.check`)."""


comptime CHECK_MARGIN = 1
"""Tenths of its projection a checkpoint's bound allows, for each checkpoint still to come: its
projection strayed by up to about an eighth at the first, a twentieth at the last."""


comptime CHECK_ROWS = 8
"""Rows between the scores a checkpoint samples down the band."""


comptime DOUBLING_START = 256
"""How far past the distance already ruled out an untrusted band's first bound reaches: A*PA2's own
first step past the heuristic at the origin."""


def band_doubling[
    record: Bool
](
    mut profile: Profile,
    give_up_wide: Bool,
    probe: Probe,
    mut trail: Trail,
    mut heuristic: SeedHeuristic,
    trusted: Bool,
) -> Int:
    """Band doubling: the distance, or -1 once the band would cover the matrix and `give_up_wide`
    hands it back for a sweep of the whole.

    An untrusted projection (see `trusted_projection`) neither aims the first bound nor sets the
    next: the first starts `DOUBLING_START` past the distance already ruled out, by the floor or the
    heuristic at the origin, and each next one at most doubles, as A*PA2's band doubling grows,
    whatever a failed round projected.

    Each round sweeps the band for one bound. A failed round leaves either a real alignment's cost,
    which caps the next bound, or the column at which it pruned every row, from which the distance is
    estimated as the bound scaled to the whole width; the next bound aims just past the estimate
    rather than doubling. Once the band would
    cover the matrix, `give_up_wide` hands back to the caller, else the bound is lifted past any path.

    The first bound aims just past where the diagonal transition's `probe` projected the distance,
    and above every score it ruled out.
    """
    var columns = profile.columns
    var rows = profile.rows
    var gap = abs(rows - columns)
    var aimed = probe.estimate + probe.estimate // 8
    if columns <= SHORT_COLUMNS:
        aimed = min(probe.estimate * SHORT_AIM // 10, probe.estimate + SHORT_REACH)
    # The heuristic at the origin is a lower bound on the distance. With seeds it lands within about a
    # sixth of it on a close pair, where the projection from a few edits strays further.
    var origin = heuristic.h(0, 0)
    var threshold = max(gap, probe.floor + 1, aimed + PROBE_MARGIN)
    if not trusted:
        # As A*PA2 starts: a step past what is known, the heuristic at the origin or the floor.
        threshold = min(threshold, max(gap, probe.floor + 1, origin) + DOUBLING_START)
    if heuristic.seeds > 0:
        if heuristic.chains_well():
            # Matches survive: the origin's bound lies within a few percent of the distance on a
            # close pair, closer than any projection, so start just past it. On a more divergent
            # pair the round dies early, a share of the way across proportional to how far the
            # bound sits above the origin's, and its death estimates the rest (see below).
            threshold = origin + SEED_SLACK
        # Otherwise few matches survive and the bound is one edit a seed, well short of a divergent
        # pair's distance: the projection leads, the bound only floors it.
        threshold = max(threshold, gap, probe.floor + 1, origin + SEED_SLACK)
    var best = Int.MAX
    # Only the first round re-aims at checkpoints and may give itself up there: on reads whose errors
    # gather at an end, a later round given up on its own climb would jump past a bound that was enough.
    var first_round = True
    while True:
        if 2 * (threshold + BAND_COLUMNS) >= rows:
            if give_up_wide:
                return -1
            threshold = max(threshold, columns + rows)
        # A round with seeds and one without each take their own copy (see `SeedHeuristic.bound`).
        var attempt: Round
        if heuristic.seeds > 0:
            attempt = pruned_distance[record, True](profile, threshold, trail, heuristic, first_round)
        else:
            attempt = pruned_distance[record, False](profile, threshold, trail, heuristic, first_round)
        first_round = False

        var found = attempt.distance
        # A checkpoint may have lowered the bound, and only a distance within the lowered one is exact.
        var bound = attempt.bound
        if found >= 0 and found <= bound:
            return found

        # Grown from the bound the round ended on, which a checkpoint may have lowered: from the one it
        # started on, a round lowered and then failed would jump back to a bound it already knew was loose.
        var next_bound = 2 * bound
        if found >= 0:
            # Still the cost of a real alignment, so it caps every later bound.
            best = min(best, found)
        else:
            # The round pruned every row `reached` columns into its sweep, where the best alignment's
            # score plus heuristic had climbed from the heuristic at the origin past the bound;
            # that climb scaled to the whole width estimates the distance.
            var estimate = Int.MAX
            if attempt.estimate >= 0:
                # A checkpoint gave the round up, with its own projection.
                estimate = attempt.estimate
            elif attempt.reached > 0 and (attempt.reached < columns or heuristic.seeds > 0):
                # A seeded round that crossed every column and still missed the end climbed past the
                # bound in its last tile: the bound itself is the projection, and the retry aims its
                # margin past it rather than doubling a bound that may be nearly enough.
                estimate = min(estimate, origin + (bound - origin) * columns // attempt.reached)
            if estimate != Int.MAX:
                next_bound = max(bound + bound // 4, estimate + estimate // 8 + PROBE_MARGIN)
                if heuristic.seeds > 0:
                    # The origin's bound is certain; only the climb above it is estimated, and the
                    # retry aims a quarter of that climb past it: half overshot the final bound by a
                    # tenth or more on the long reads, widening the band that succeeds.
                    next_bound = max(bound + SEED_SLACK, estimate + (estimate - origin) // 4 + PROBE_MARGIN)
        if found < 0 and heuristic.seeds > 0 and attempt.reached * SEEDED_TRUST_SHARE < columns:
            # A seeded round that died within its first columns projects from those alone, which on
            # real reads hold their errors gathered at the start: its margin over the origin's bound
            # grows `SEEDED_GROWTH` times instead, as A*PA2's grows, until a round gets far enough in
            # for its death to say where the distance lies.
            next_bound = origin + SEEDED_GROWTH * max(bound - origin, SEED_SLACK)
        if not trusted:
            # A round's estimate comes from where it stopped, which errors gathered at an end can
            # set as far off as the first projection: at most doubling instead.
            next_bound = min(next_bound, 2 * bound)
        threshold = min(next_bound, best)


def band_start(profile: Profile, search: Probe, mut probe: Probe, mut trusted: Bool) -> SeedHeuristic:
    """What a band starts from once diagonal transition has left the distance to it: the seed
    heuristic its gate chooses, returned, with `probe` set to the projection to aim at and `trusted`
    to whether it agrees with `search`'s (see `trusted_projection`).

    The band's bounds and the seeds' gate are tuned on one front's projection from its first
    `PROBE_START` edits, so that is the projection handed on, unless `search` went far further.
    """
    var fronts = DiagonalFronts()
    var projected = diagonal_transition(profile, PROJECTION_ONLY, fronts)
    trusted = trusted_projection(search, projected)
    probe = Probe(-1, projected.estimate, max(search.floor, projected.floor))
    var seeded = profile.columns >= SEED_COLUMNS or (
        SHORT_SEEDS and projected.estimate >= SEED_EDITS and projected.estimate * SEED_DIVERGENCE <= profile.columns
    )
    if not seeded:
        return SeedHeuristic(profile.columns, profile.rows)
    # Long pairs only: inexact seeds when the projection already says they pay, else exact ones
    # rebuilt inexact if they chain poorly (see `SeedHeuristic`).
    var long = profile.columns >= INEXACT_COLUMNS
    var inexact = long and projected.estimate * INEXACT_DIVERGENCE >= profile.columns
    var cutoff = INEXACT_CHAINED_SPREAD if trusted else INEXACT_CHAINED
    return SeedHeuristic(profile, inexact, choose=long, cutoff=min(cutoff, INEXACT_CHAINED))
