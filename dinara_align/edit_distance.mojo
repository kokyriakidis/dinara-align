# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Ported from `pa-bitpacking` in A*PA (https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner,
# commit bf2e14e), by Ragnar Groot Koerkamp and Pesho Ivanov, itself translated from Edlib.
"""
The global unit-cost edit distance between two sequences, and an optimal alignment, by bit-parallel
sweep: a substitution, an insertion and a deletion each cost one, `Costs.edit()`.

A pair first goes through diagonal transition (see `diagonal`), which settles near-identical pairs and
projects any other's distance. Band doubling (see `band`) then sweeps only the cells a path within a
bound could cross, its rows pruned with the seed heuristic on long pairs (see `seeds`), with the
bit-parallel kernels of `bit_parallel`; a pair too divergent for a band is swept whole. The alignment
is traced back by a fixed rule for ties, whichever search found the distance (see `edit_cigar`).

Each pair runs on one thread; `distances` and `alignments` in `api` spread a batch over threads a pair
at a time.
"""

from .band import band_doubling, band_start, retrace
from .bit_parallel import FIRST_SENTINEL, full_distance, LEFT, Profile, SECOND_SENTINEL, Trail, UP
from .diagonal import (
    SHORT_BAND_COLUMNS,
    full_matrix_steps,
    diagonal_transition,
    DiagonalFronts,
    FrontPair,
    reversed_codes_into,
    ALIGNMENT_COST,
    DISTANCE_COST,
    EDITS_PER_STEP,
    StepCost,
    grow_to,
    Probe,
    PROJECTION_ONLY,
    two_ended,
    two_ended_distance,
    two_ended_setup,
)
from .errors import AlignmentError
from .modes import Ties
from .probes import (
    BAND_NS,
    CIGAR_NS,
    DIAGONAL_NS,
    DIAGONAL_SCORE,
    DIAGONAL_STEPS,
    DISTANCE,
    ENCODE_NS,
    ENGINE,
    ENGINE_BAND,
    ENGINE_DIAGONAL,
    ENGINE_EMPTY,
    ENGINE_SHORT_FULL,
    ENGINE_TWO_ENDED,
    ENGINE_WIDE_FULL,
    ESTIMATE,
    FLOOR,
    FORCE_BAND,
    FORCE_DIAGONAL,
    FORCE_FULL,
    FORCED_STEP_TENTHS,
    front_steps,
    FULL_NS,
    FULL_WORD_COLUMNS,
    GROW_NS,
    GROW_STEPS,
    PROBES,
    SEEDS_NS,
    Tally,
    TRACEBACK_NS,
    TRUSTED,
    TWO_ENDED_NS,
    TWO_ENDED_SCORE,
    TWO_ENDED_STEPS,
)
from .seeds import SeedHeuristic
from .traceback import cigar_string, diagonal_cigar, EditPath


struct EditSpace(Movable):
    """A batch worker's memory for its pairs' edit distances: the profile, the reversed codes and the
    two-ended search's rings, kept from pair to pair so a warm worker allocates nothing, where each
    pair's own were a third of a short read's time and, on many threads, contended in the allocator;
    and for an alignment, the fronts it traces back."""

    var profile: Profile
    var first_back: List[UInt8]
    var second_back: List[UInt8]
    var ahead: FrontPair
    var behind: FrontPair
    var traced: DiagonalFronts
    """The fronts an alignment settled by diagonal transition traces back (see `settle`)."""
    var traced_ahead: FrontPair
    """An alignment's two-ended search from the start, recording, kept from pair to pair (see `settle`)."""
    var traced_behind: FrontPair
    """Its search from the end."""
    var tally: Tally
    """The last pair's probes, in a build with them (see `probes`)."""

    def __init__(out self):
        """Nothing profiled yet: the first pair does."""
        self.profile = Profile()
        self.first_back = List[UInt8]()
        self.second_back = List[UInt8]()
        self.ahead = FrontPair()
        self.behind = FrontPair()
        self.traced = DiagonalFronts(reserve=False)
        self.traced_ahead = FrontPair(record=True)
        self.traced_behind = FrontPair(record=True)
        self.tally = Tally()


def edit_distance(first: String, second: String) raises AlignmentError -> Int:
    """The global edit distance between two sequences (see the overload taking an `EditSpace`)."""
    var space = EditSpace()
    return edit_distance(first, second, space)


def edit_distance(first: String, second: String, mut space: EditSpace) raises AlignmentError -> Int:
    """The global edit distance between two sequences, by bit-parallel sweep, through `space`'s memory.

    Built for DNA over `ACGT`. Up to four other bytes, `N` among them, are symbols of their own,
    each matching only itself; a pair holding them sweeps a third bit plane, so runs a little slower
    than bases alone would.

    Band doubling, as in A*PA2-simple: guess a bound, sweep only the band of cells a path within it
    could cross, and raise the guess until the answer fits under it, which proves it optimal. Close
    sequences therefore cost far less than the whole matrix. Once the band would cover most of the
    matrix, the whole matrix is swept instead; a short pair's band covers it from the start, so the
    diagonal transition gives way straight to the sweep once it would cost more.
    """
    # A forced policy (see `probes.FORCE`) takes the pair before the library's own choices.
    comptime if FORCE_FULL or FORCE_BAND:
        return forced_distance(first, second, space)
    # The probes (see `probes`) are blocks a build without them compiles to nothing, the code around them
    # main's line for line: moving code into helpers for them changed what the compiler inlined.
    var started = space.tally.start()
    space.profile.reset(first, second)
    ref profile = space.profile
    space.tally.since(ENCODE_NS, started)
    if profile.columns == 0 or profile.rows == 0:
        space.tally.set(ENGINE, ENGINE_EMPTY)
        return profile.columns + profile.rows
    # A short pair's band would sweep its whole matrix (see `SHORT_BAND_COLUMNS`), so past what that
    # sweep costs the diagonal transition gives way to the sweep itself.
    var short = profile.columns <= SHORT_BAND_COLUMNS
    started = space.tally.start()
    var search = two_ended_distance(
        profile,
        StepCost(FORCED_STEP_TENTHS, EDITS_PER_STEP) if FORCE_DIAGONAL else DISTANCE_COST,
        full_matrix_steps(profile.columns, profile.rows) if short and not FORCE_DIAGONAL else Int.MAX,
        space.first_back,
        space.second_back,
        space.ahead,
        space.behind,
    )
    comptime if PROBES:
        space.tally.since(TWO_ENDED_NS, started)
        space.tally.set(TWO_ENDED_SCORE, space.ahead.score + space.behind.score)
        space.tally.set(
            TWO_ENDED_STEPS,
            front_steps(space.ahead.score, profile.columns, profile.rows)
            + front_steps(space.behind.score, profile.columns, profile.rows),
        )
        space.tally.set(ESTIMATE, search.estimate)
        space.tally.set(FLOOR, search.floor)
    if search.distance >= 0:
        space.tally.set(ENGINE, ENGINE_TWO_ENDED)
        return search.distance
    if short:
        comptime if PROBES:
            return swept_whole(profile, space.tally, ENGINE_SHORT_FULL)
        else:
            return full_distance(profile)
    var probe = search
    var trusted = True
    started = space.tally.start()
    var heuristic = band_start(profile, search, probe, trusted)
    comptime if PROBES:
        space.tally.since(SEEDS_NS, started)
        space.tally.set(TRUSTED, Int(trusted))
        space.tally.set(ESTIMATE, probe.estimate)
        heuristic.tally_into(space.tally)
    var trail = Trail()
    started = space.tally.start()
    var distance = band_doubling[False](profile, True, probe, trail, heuristic, trusted, space.tally)
    space.tally.since(BAND_NS, started)
    if distance >= 0:
        space.tally.set(ENGINE, ENGINE_BAND)
        return distance
    comptime if PROBES:
        return swept_whole(profile, space.tally, ENGINE_WIDE_FULL)
    else:
        return full_distance(profile)


def swept_whole(mut profile: Profile, mut tally: Tally, engine: Int) -> Int:
    """`full_distance`, its time and cells tallied and `engine` named as what settled the pair, in a build with
    probes or a forced policy."""
    var started = tally.start()
    var distance = full_distance(profile)
    tally.since(FULL_NS, started)
    tally.set(ENGINE, engine)
    tally.set(FULL_WORD_COLUMNS, profile.columns * profile.words)
    return distance


def projection(profile: Profile) -> Probe:
    """One front's projection of the distance from its first edits, as `band_start` takes it, settling
    nothing: what a band forced to take a pair (see `probes.FORCE`) starts from."""
    var fronts = DiagonalFronts()
    var projected = diagonal_transition(profile, PROJECTION_ONLY, fronts)
    return Probe(-1, projected.estimate, projected.floor)


def forced_distance(first: String, second: String, mut space: EditSpace) raises AlignmentError -> Int:
    """The distance under a forced policy (see `probes.FORCE`): the whole matrix swept, or band doubling from
    one front's projection, the diagonal transition settling nothing."""
    space.profile.reset(first, second)
    ref profile = space.profile
    if profile.columns == 0 or profile.rows == 0:
        return profile.columns + profile.rows
    comptime if FORCE_FULL:
        return swept_whole(profile, space.tally, ENGINE_WIDE_FULL)
    if profile.columns <= SHORT_BAND_COLUMNS:
        # The band of a short pair covers the matrix, so band doubling would hand it straight back.
        return swept_whole(profile, space.tally, ENGINE_SHORT_FULL)
    var search = projection(profile)
    var probe = search
    var trusted = True
    var heuristic = band_start(profile, search, probe, trusted)
    heuristic.tally_into(space.tally)
    var trail = Trail()
    var distance = band_doubling[False](profile, True, probe, trail, heuristic, trusted, space.tally)
    if distance >= 0:
        space.tally.set(ENGINE, ENGINE_BAND)
        return distance
    return swept_whole(profile, space.tally, ENGINE_WIDE_FULL)


@fieldwise_init
struct EditCigar(Copyable, Movable, Writable):
    """The global edit distance between two sequences and an optimal alignment as a CIGAR string, the
    first sequence the reference: `=` a match, `X` a substitution (or `M` for either), `D` a base of the
    first sequence alone, `I` a base of the second alone, each run as its length then its letter."""

    var distance: Int
    var cigar: String

    def take_cigar(deinit self) -> String:
        """The CIGAR, moved out rather than copied, the rest given up."""
        return self.cigar^


comptime EDGE_WORD_BYTES = 16
"""Bytes a kept edge word takes in a band's trail: its `+` and `-` masks."""


def edit_cigar(
    first: String, second: String, eqx: Bool = True, ties: Ties = Ties.LEFT, max_memory: Int = Int.MAX
) raises AlignmentError -> EditCigar:
    """The global edit distance between two sequences and an optimal alignment as a CIGAR string (see the
    overload taking an `EditSpace`)."""
    var space = EditSpace()
    return edit_cigar(first, second, space, eqx, ties, max_memory)


def edit_cigar(
    first: String,
    second: String,
    mut space: EditSpace,
    eqx: Bool = True,
    ties: Ties = Ties.LEFT,
    max_memory: Int = Int.MAX,
) raises AlignmentError -> EditCigar:
    """The global edit distance between two sequences, and an optimal alignment as a CIGAR string, `=`
    and `X` for matches and substitutions, or with `eqx` false `M` for both (see `EditCigar`).

    Of several optimal alignments the CIGAR is always the one `ties` names, as for the gap-affine
    alignment (see `Ties`): by default every edit as far left as it goes, indels placed as minimap2
    places them, or with `Ties.RIGHT` as far right, WFA2-lib's edit CIGAR byte for byte; whichever of
    the searches below found the distance. A distance diagonal transition found has its CIGAR written
    straight from the fronts, a run at a time (see `diagonal_cigar`).

    A band keeps its tiles' edges for the traceback within `max_memory` bytes, sweeping stretches again
    past that (see `band.retrace`): the same CIGAR, for one more sweep of the band, or a few under a cap
    far below the band's edges. The profile is `space`'s, kept from pair to pair as `edit_distance`
    keeps it."""
    # The search's rule places every edit as far right as it goes; over the reversed pair that is as far
    # left over the pair itself.
    var reverse = ties == Ties.LEFT
    var started = space.tally.start()
    space.profile.reset(first, second, reverse)
    ref profile = space.profile
    space.tally.since(ENCODE_NS, started)
    var settled = settle(
        profile,
        space.traced,
        space.first_back,
        space.second_back,
        space.traced_ahead,
        space.traced_behind,
        max_memory // EDGE_WORD_BYTES,
        space.tally,
    )
    comptime if PROBES:
        started = space.tally.start()
    if not settled.banded:
        comptime if PROBES:
            var cigar = diagonal_cigar(profile, space.traced, settled.distance, reverse, eqx)
            space.tally.since(CIGAR_NS, started)
            return EditCigar(settled.distance, cigar^)
        else:
            return EditCigar(settled.distance, diagonal_cigar(profile, space.traced, settled.distance, reverse, eqx))
    var path = oriented(settled^, profile.columns, profile.rows, reverse)
    comptime if PROBES:
        var cigar = cigar_string(first, second, path, eqx)
        space.tally.since(CIGAR_NS, started)
        return EditCigar(path.distance, cigar^)
    else:
        return EditCigar(path.distance, cigar_string(first, second, path, eqx))


def oriented(var settled: Settled, columns: Int, rows: Int, reverse: Bool) -> EditPath:
    """The moves traced, right to left from the corner, as a path over the pair: over the reversed pair
    they run left to right over the pair itself."""
    var distance = settled.distance
    if reverse:
        return EditPath(List[UInt8](), settled^.take_moves(), 0, 0, distance)
    return EditPath(settled^.take_moves(), List[UInt8](), columns, rows, distance)


struct Settled(Movable):
    """A pair's distance and how it was found: by diagonal transition, whose fronts `settle` kept hold
    the path still to trace, or by band doubling, its path already traced into `moves`, right to left."""

    var distance: Int
    var moves: List[UInt8]
    var banded: Bool

    def __init__(out self, distance: Int, var moves: List[UInt8], banded: Bool):
        """A settled pair from its parts, taking ownership of `moves`."""
        self.distance = distance
        self.moves = moves^
        self.banded = banded

    def take_moves(deinit self) -> List[UInt8]:
        """The traced moves, the rest given up."""
        return self.moves^


def settle(
    mut profile: Profile,
    mut traced: DiagonalFronts,
    mut first_back: List[UInt8],
    mut second_back: List[UInt8],
    mut ahead: FrontPair,
    mut behind: FrontPair,
    budget: Int,
    mut tally: Tally,
) raises AlignmentError -> Settled:
    """The global edit distance between a profile's two sequences, and what an optimal alignment by
    WFA2-lib's rule for ties is traced from: from the corner back, an edit into each cell whenever one
    is optimal, a substitution before a base of the first sequence alone before one of the second.

    A near-identical pair is settled by diagonal transition from the start, its fronts kept whole to be
    traced back from the corner, which is that rule. Otherwise diagonal transition from both ends finds
    the distance while it is cheaper than a band, and the forward fronts then grow on to it, pruned to
    the diagonals an optimal path passes (see `grow_to`), to be traced back the same way. A pair too
    far apart for either goes to band doubling, recording each tile's left edge in the round that
    succeeds, or within `budget` words only some, and each tile is traced by the same rule (see
    `band.retrace`).

    The fronts to trace back are left in `traced`, and the two-ended search runs over `first_back`,
    `second_back`, `ahead` and `behind`, all of whose memory a reused space keeps from pair to pair; the
    pair's probes go in `tally` (see `probes`).
    """
    var columns = profile.columns
    var rows = profile.rows
    if columns == 0 or rows == 0:
        var moves = List[UInt8](capacity=columns + rows)
        for _ in range(rows):
            moves.append(UP)
        for _ in range(columns):
            moves.append(LEFT)
        tally.set(ENGINE, ENGINE_EMPTY)
        return Settled(columns + rows, moves^, True)
    comptime if FORCE_FULL:
        return banded_alignment(profile, Probe(-1, columns + rows, 0), budget, tally, True)
    comptime if FORCE_BAND:
        return banded_alignment(profile, projection(profile), budget, tally, False)
    comptime cost = StepCost(FORCED_STEP_TENTHS, EDITS_PER_STEP) if FORCE_DIAGONAL else ALIGNMENT_COST
    # A near-identical pair: one front, kept whole, settles it before both ends' setup would pay.
    var started = tally.start()
    traced.clear()
    var close = diagonal_transition(profile, cost, traced, switch_setup=two_ended_setup(columns, rows))
    comptime if PROBES:
        tally.since(DIAGONAL_NS, started)
        var steps = 0
        for score in range(len(traced.lows)):
            steps += traced.highs[score] - traced.lows[score] + 1
        tally.set(DIAGONAL_STEPS, steps)
        tally.set(DIAGONAL_SCORE, len(traced.lows) - 1)
        tally.set(ESTIMATE, close.estimate)
    if close.distance >= 0:
        tally.set(ENGINE, ENGINE_DIAGONAL)
        tally.set(DISTANCE, close.distance)
        return Settled(close.distance, List[UInt8](), False)
    started = tally.start()
    reversed_codes_into(first_back, profile.column_codes, columns, FIRST_SENTINEL)
    reversed_codes_into(second_back, profile.row_codes, rows, SECOND_SENTINEL)
    ahead.restart(record=True)
    behind.restart(record=True)
    var meeting = two_ended(profile, first_back, second_back, cost, ahead, behind)
    comptime if PROBES:
        tally.since(TWO_ENDED_NS, started)
        tally.set(TWO_ENDED_SCORE, ahead.score + behind.score)
        tally.set(TWO_ENDED_STEPS, front_steps(ahead.score, columns, rows) + front_steps(behind.score, columns, rows))
        tally.set(ESTIMATE, meeting.probe.estimate)
        tally.set(FLOOR, meeting.probe.floor)
    if meeting.probe.distance >= 0:
        var distance = meeting.probe.distance
        started = tally.start()
        grow_to(profile, ahead.history, behind.history, distance)
        # The grown fronts become the ones traced, and the search keeps the memory they had.
        swap(traced, ahead.history)
        comptime if PROBES:
            tally.since(GROW_NS, started)
            var grown = 0
            for score in range(len(traced.lows)):
                grown += traced.highs[score] - traced.lows[score] + 1
            tally.set(GROW_STEPS, grown)
            tally.set(ENGINE, ENGINE_TWO_ENDED)
            tally.set(DISTANCE, distance)
        return Settled(distance, List[UInt8](), False)
    var probe = meeting.probe
    var trusted = True
    started = tally.start()
    var heuristic = band_start(profile, meeting.probe, probe, trusted)
    comptime if PROBES:
        tally.since(SEEDS_NS, started)
        tally.set(TRUSTED, Int(trusted))
        tally.set(ESTIMATE, probe.estimate)
        heuristic.tally_into(tally)
    var trail = Trail(columns, budget)
    started = tally.start()
    var distance = band_doubling[True](profile, False, probe, trail, heuristic, trusted, tally)
    comptime if PROBES:
        tally.since(BAND_NS, started)
        started = tally.start()
    var moves = List[UInt8](capacity=columns + rows)
    retrace(profile, trail, heuristic, columns, rows, distance, moves, tally)
    comptime if PROBES:
        tally.since(TRACEBACK_NS, started)
        tally.set(ENGINE, ENGINE_BAND)
    return Settled(distance, moves^, True)


def banded_alignment(mut profile: Profile, search: Probe, budget: Int, mut tally: Tally, whole: Bool) -> Settled:
    """Under a forced policy (see `probes.FORCE`), the distance by band doubling from `search`, its path traced
    back over the trail, its edges kept within `budget` words; with `whole`, as one band over the whole matrix
    and the gap heuristic alone."""
    var columns = profile.columns
    var rows = profile.rows
    var probe = search
    var trusted = True
    var started = tally.start()
    var heuristic: SeedHeuristic
    if whole:
        heuristic = SeedHeuristic(columns, rows)
    else:
        heuristic = band_start(profile, search, probe, trusted)
    tally.since(SEEDS_NS, started)
    tally.set(TRUSTED, Int(trusted))
    tally.set(ESTIMATE, probe.estimate)
    heuristic.tally_into(tally)
    var trail = Trail(columns, budget)
    started = tally.start()
    var distance = band_doubling[True](profile, False, probe, trail, heuristic, trusted, tally)
    tally.since(BAND_NS, started)
    started = tally.start()
    var moves = List[UInt8](capacity=columns + rows)
    retrace(profile, trail, heuristic, columns, rows, distance, moves, tally)
    tally.since(TRACEBACK_NS, started)
    tally.set(ENGINE, ENGINE_WIDE_FULL if whole else ENGINE_BAND)
    return Settled(distance, moves^, True)
