# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Ported from `pa-bitpacking` in A*PA (https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner,
# commit bf2e14e), by Ragnar Groot Koerkamp and Pesho Ivanov, itself translated from Edlib.
"""
The global unit-cost edit distance between two sequences, and an optimal alignment, by bit-parallel
sweep: a substitution, an insertion and a deletion each cost one, the distance `levenshtein_alignment`
returns.

A pair first goes through diagonal transition (see `diagonal`), which settles near-identical pairs and
projects any other's distance. Band doubling (see `band`) then sweeps only the cells a path within a
bound could cross, its rows pruned with the seed heuristic on long pairs (see `seeds`), with the
bit-parallel kernels of `bit_parallel`; a pair too divergent for a band is swept whole. The alignment
is traced back by a fixed rule for ties, whichever search found the distance (see `edit_path`).

Each pair runs on one thread; `edit_distances` and `edit_alignments` in `api` spread a batch over
threads a pair at a time.
"""

from .alignment import AlignmentResult
from .band import band_doubling, band_start
from .bit_parallel import FIRST_SENTINEL, full_distance, LEFT, Profile, SECOND_SENTINEL, Trail, UP
from .diagonal import (
    diagonal_transition,
    DiagonalFronts,
    FrontPair,
    reversed_codes,
    STEP_TENTHS_ALIGNMENT,
    STEP_TENTHS_DISTANCE,
    grow_to,
    trace_diagonals,
    two_ended,
    two_ended_distance,
    TWO_ENDED_SETUP,
)
from .errors import AlignmentError
from .gap_affine import Ties
from .traceback import cigar_string, EditPath, gapped_rows, trace_back


def edit_distance(first: String, second: String) raises AlignmentError -> Int:
    """The global edit distance between two sequences, by bit-parallel sweep.

    Built for DNA over `ACGT`. Up to four other bytes, `N` among them, are symbols of their own,
    each matching only itself; a pair holding them sweeps a third bit plane, so runs a little slower
    than bases alone would.

    Band doubling, as in A*PA2-simple: guess a bound, sweep only the band of cells a path within it
    could cross, and raise the guess until the answer fits under it, which proves it optimal. Close
    sequences therefore cost far less than the whole matrix. Once the band would cover most of the
    matrix, the whole matrix is swept instead.
    """
    var profile = Profile(first, second)
    if profile.columns == 0 or profile.rows == 0:
        return profile.columns + profile.rows
    var search = two_ended_distance(profile, STEP_TENTHS_DISTANCE)
    if search.distance >= 0:
        return search.distance
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


def edit_alignment(first: String, second: String, ties: Ties = Ties.LEFT) raises AlignmentError -> AlignmentResult:
    """The global edit distance between two sequences, and an optimal alignment as two gapped rows;
    symbols past `ACGT` as `edit_distance` takes them. The score is the distance, as
    `levenshtein_alignment` reports it. Of several optimal alignments, the one `ties` names (see
    `edit_cigar`)."""
    var path = edit_path(first, second, ties)
    return gapped_rows(first, second, path.prefix, path.middle_column, path.middle_row, path.suffix, path.distance)


def edit_cigar(
    first: String, second: String, extended: Bool = True, ties: Ties = Ties.LEFT
) raises AlignmentError -> EditCigar:
    """The global edit distance between two sequences, and an optimal alignment as a CIGAR string, `=`
    and `X` for matches and substitutions, or with `extended` false `M` for both (see `EditCigar`);
    the same alignment `edit_alignment` writes out, without the gapped rows.

    Of several optimal alignments the CIGAR is always the one `ties` names, as for the gap-affine
    alignment (see `Ties`): by default every edit as far left as it goes, indels placed as minimap2
    places them, or with `Ties.RIGHT` as far right, WFA2-lib's edit CIGAR byte for byte; whichever of
    the searches below found the distance."""
    var path = edit_path(first, second, ties)
    return EditCigar(path.distance, cigar_string(first, second, path, extended))


def edit_path(first: String, second: String, ties: Ties = Ties.LEFT) raises AlignmentError -> EditPath:
    """The global edit distance between two sequences and the optimal alignment `ties` names.

    The right rule is WFA2-lib's backtrace from the corner (see `right_moves`); the left rule is that
    over both sequences reversed, its moves read the other way."""
    var moves = List[UInt8](capacity=first.byte_length() + second.byte_length())
    if ties == Ties.RIGHT:
        var distance = right_moves(first, second, False, moves)
        return EditPath(moves^, List[UInt8](), first.byte_length(), second.byte_length(), distance)
    var distance = right_moves(first, second, True, moves)
    # The reversed pair's moves, right to left over it, run left to right over the pair.
    return EditPath(List[UInt8](), moves^, 0, 0, distance)


def right_moves(first: String, second: String, reverse: Bool, mut moves: List[UInt8]) raises AlignmentError -> Int:
    """The global edit distance between two sequences, both back to front with `reverse`, appending an
    optimal alignment's moves right to left from the corner, by WFA2-lib's rule for ties: from the corner back, an edit into each cell whenever
    one is optimal, a substitution before a base of the first sequence alone before one of the second.

    A near-identical pair is settled by diagonal transition from the start, kept whole, and traced back
    from the corner, which is that rule. Otherwise diagonal transition from both ends finds the distance
    while it is cheaper than a band, and the forward fronts then grow on to it, pruned to the diagonals
    an optimal path passes (see `grow_to`), to be traced back the same way. A pair too far apart for
    either goes to band doubling, recording each tile's left edge in the round that succeeds, and each
    tile is swept again and traced cell by cell by the same rule (see `trace_back`).
    """
    var profile = Profile(first, second, reverse)
    var columns = profile.columns
    var rows = profile.rows
    if columns == 0 or rows == 0:
        for _ in range(rows):
            moves.append(UP)
        for _ in range(columns):
            moves.append(LEFT)
        return columns + rows
    # A near-identical pair: one front, kept whole, settles it before both ends' setup would pay.
    var near = DiagonalFronts()
    var close = diagonal_transition(profile, STEP_TENTHS_ALIGNMENT, near, switch_setup=TWO_ENDED_SETUP)
    if close.distance >= 0:
        trace_diagonals(profile, near, close.distance, moves)
        return close.distance
    var first_back = reversed_codes(profile.column_codes, columns, FIRST_SENTINEL)
    var second_back = reversed_codes(profile.row_codes, rows, SECOND_SENTINEL)
    var ahead = FrontPair(record=True)
    var behind = FrontPair(record=True)
    var meeting = two_ended(profile, first_back, second_back, STEP_TENTHS_ALIGNMENT, ahead, behind)
    if meeting.probe.distance >= 0:
        var distance = meeting.probe.distance
        grow_to(profile, ahead.history, behind.history, distance)
        trace_diagonals(profile, ahead.history, distance, moves)
        return distance
    var probe = meeting.probe
    var trusted = True
    var heuristic = band_start(profile, meeting.probe, probe, trusted)
    var trail = Trail(columns)
    var distance = band_doubling[True](profile, False, probe, trail, heuristic, trusted)
    trace_back(profile, trail, columns, rows, distance, moves)
    return distance
