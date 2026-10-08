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

from .band import band_doubling, band_start
from .bit_parallel import FIRST_SENTINEL, full_distance, LEFT, Profile, SECOND_SENTINEL, Trail, UP
from .diagonal import (
    SHORT_BAND_COLUMNS,
    full_matrix_steps,
    diagonal_transition,
    DiagonalFronts,
    FrontPair,
    reversed_codes,
    STEP_TENTHS_ALIGNMENT,
    STEP_TENTHS_DISTANCE,
    grow_to,
    two_ended,
    two_ended_distance,
    two_ended_setup,
)
from .errors import AlignmentError
from .modes import Ties
from .traceback import cigar_string, diagonal_cigar, EditPath, trace_back


def edit_distance(first: String, second: String) raises AlignmentError -> Int:
    """The global edit distance between two sequences, by bit-parallel sweep.

    Built for DNA over `ACGT`. Up to four other bytes, `N` among them, are symbols of their own,
    each matching only itself; a pair holding them sweeps a third bit plane, so runs a little slower
    than bases alone would.

    Band doubling, as in A*PA2-simple: guess a bound, sweep only the band of cells a path within it
    could cross, and raise the guess until the answer fits under it, which proves it optimal. Close
    sequences therefore cost far less than the whole matrix. Once the band would cover most of the
    matrix, the whole matrix is swept instead; a short pair's band covers it from the start, so the
    diagonal transition gives way straight to the sweep once it would cost more.
    """
    var profile = Profile(first, second)
    if profile.columns == 0 or profile.rows == 0:
        return profile.columns + profile.rows
    # A short pair's band would sweep its whole matrix (see `SHORT_BAND_COLUMNS`), so past what that
    # sweep costs the diagonal transition gives way to the sweep itself.
    var short = profile.columns <= SHORT_BAND_COLUMNS
    var search = two_ended_distance(
        profile, STEP_TENTHS_DISTANCE, full_matrix_steps(profile.columns, profile.rows) if short else Int.MAX
    )
    if search.distance >= 0:
        return search.distance
    if short:
        return full_distance(profile)
    var probe = search
    var trusted = True
    var heuristic = band_start(profile, search, probe, trusted)
    var trail = Trail()
    var distance = band_doubling[False](profile, True, probe, trail, heuristic, trusted)
    if distance >= 0:
        return distance
    return full_distance(profile)


@fieldwise_init
struct EditCigar(Copyable, Movable, Writable):
    """The global edit distance between two sequences and an optimal alignment as a CIGAR string, the
    first sequence the reference: `=` a match, `X` a substitution (or `M` for either), `D` a base of the
    first sequence alone, `I` a base of the second alone, each run as its length then its letter."""

    var distance: Int
    var cigar: String


def edit_cigar(
    first: String, second: String, eqx: Bool = True, ties: Ties = Ties.LEFT
) raises AlignmentError -> EditCigar:
    """The global edit distance between two sequences, and an optimal alignment as a CIGAR string, `=`
    and `X` for matches and substitutions, or with `eqx` false `M` for both (see `EditCigar`).

    Of several optimal alignments the CIGAR is always the one `ties` names, as for the gap-affine
    alignment (see `Ties`): by default every edit as far left as it goes, indels placed as minimap2
    places them, or with `Ties.RIGHT` as far right, WFA2-lib's edit CIGAR byte for byte; whichever of
    the searches below found the distance. A distance diagonal transition found has its CIGAR written
    straight from the fronts, a run at a time (see `diagonal_cigar`)."""
    var reverse = ties == Ties.LEFT
    var profile = Profile(first, second, reverse)
    var settled = settle(profile)
    if not settled.banded:
        return EditCigar(settled.distance, diagonal_cigar(profile, settled.fronts, settled.distance, reverse, eqx))
    var path = oriented(settled^, profile.columns, profile.rows, reverse)
    return EditCigar(path.distance, cigar_string(first, second, path, eqx))


def oriented(var settled: Settled, columns: Int, rows: Int, reverse: Bool) -> EditPath:
    """The moves traced, right to left from the corner, as a path over the pair: over the reversed pair
    they run left to right over the pair itself."""
    var distance = settled.distance
    if reverse:
        return EditPath(List[UInt8](), settled^.take_moves(), 0, 0, distance)
    return EditPath(settled^.take_moves(), List[UInt8](), columns, rows, distance)


struct Settled(Movable):
    """A pair's distance and how it was found: by diagonal transition, whose kept `fronts` hold the path
    still to trace, or by band doubling, its path already traced into `moves`, right to left."""

    var distance: Int
    var fronts: DiagonalFronts
    var moves: List[UInt8]
    var banded: Bool

    def __init__(out self, distance: Int, var fronts: DiagonalFronts, var moves: List[UInt8], banded: Bool):
        self.distance = distance
        self.fronts = fronts^
        self.moves = moves^
        self.banded = banded

    def take_moves(deinit self) -> List[UInt8]:
        return self.moves^


def settle(mut profile: Profile) raises AlignmentError -> Settled:
    """The global edit distance between a profile's two sequences, and what an optimal alignment by
    WFA2-lib's rule for ties is traced from: from the corner back, an edit into each cell whenever one
    is optimal, a substitution before a base of the first sequence alone before one of the second.

    A near-identical pair is settled by diagonal transition from the start, its fronts kept whole to be
    traced back from the corner, which is that rule. Otherwise diagonal transition from both ends finds
    the distance while it is cheaper than a band, and the forward fronts then grow on to it, pruned to
    the diagonals an optimal path passes (see `grow_to`), to be traced back the same way. A pair too
    far apart for either goes to band doubling, recording each tile's left edge in the round that
    succeeds, and each tile is traced by the same rule (see `trace_back`).
    """
    var columns = profile.columns
    var rows = profile.rows
    var moves = List[UInt8](capacity=columns + rows)
    if columns == 0 or rows == 0:
        for _ in range(rows):
            moves.append(UP)
        for _ in range(columns):
            moves.append(LEFT)
        return Settled(columns + rows, DiagonalFronts(reserve=False), moves^, True)
    # A near-identical pair: one front, kept whole, settles it before both ends' setup would pay.
    var near = DiagonalFronts()
    var close = diagonal_transition(profile, STEP_TENTHS_ALIGNMENT, near, switch_setup=two_ended_setup(columns, rows))
    if close.distance >= 0:
        return Settled(close.distance, near^, moves^, False)
    var first_back = reversed_codes(profile.column_codes, columns, FIRST_SENTINEL)
    var second_back = reversed_codes(profile.row_codes, rows, SECOND_SENTINEL)
    var ahead = FrontPair(record=True)
    var behind = FrontPair(record=True)
    var meeting = two_ended(profile, first_back, second_back, STEP_TENTHS_ALIGNMENT, ahead, behind)
    if meeting.probe.distance >= 0:
        var distance = meeting.probe.distance
        grow_to(profile, ahead.history, behind.history, distance)
        return Settled(distance, ahead^.take_history(), moves^, False)
    var probe = meeting.probe
    var trusted = True
    var heuristic = band_start(profile, meeting.probe, probe, trusted)
    var trail = Trail(columns)
    var distance = band_doubling[True](profile, False, probe, trail, heuristic, trusted)
    trace_back(profile, trail, columns, rows, distance, moves)
    return Settled(distance, DiagonalFronts(reserve=False), moves^, True)
