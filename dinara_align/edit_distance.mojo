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
is traced back from the band's tile edges (see `traceback`).

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
    trace_diagonals,
    trace_from,
    two_ended,
    two_ended_distance,
    TWO_ENDED_SETUP,
)
from .errors import AlignmentError
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


def edit_alignment(first: String, second: String) raises AlignmentError -> AlignmentResult:
    """The global edit distance between two sequences, and an optimal alignment as two gapped rows;
    symbols past `ACGT` as `edit_distance` takes them. The score is the distance, as
    `levenshtein_alignment` reports it."""
    var path = edit_path(first, second)
    return gapped_rows(first, second, path.prefix, path.middle_column, path.middle_row, path.suffix, path.distance)


def edit_cigar(first: String, second: String, extended: Bool = True) raises AlignmentError -> EditCigar:
    """The global edit distance between two sequences, and an optimal alignment as a CIGAR string, `=`
    and `X` for matches and substitutions, or with `extended` false `M` for both (see `EditCigar`);
    the same alignment `edit_alignment` writes out, without the gapped rows."""
    var path = edit_path(first, second)
    return EditCigar(path.distance, cigar_string(first, second, path, extended))


def edit_path(first: String, second: String) raises AlignmentError -> EditPath:
    """The global edit distance between two sequences and an optimal alignment's moves.

    The distance comes from `edit_distance`'s band doubling, recording each tile's left edge in the
    round that succeeds; the alignment is then traced back tile by tile from those edges (see
    `trace_back`). Where the two-ended diagonal transition settles the distance, its fronts give the
    path, traced to the start and on to the end from where they met.
    """
    var profile = Profile(first, second)
    var columns = profile.columns
    var rows = profile.rows
    # Moves right to left, from the corner, or from where the fronts met, back to the origin.
    var forward_moves = List[UInt8](capacity=columns + rows)
    # Moves left to right, from where the fronts met to the corner; none without meeting.
    var backward_moves = List[UInt8]()
    var distance: Int
    if columns == 0 or rows == 0:
        distance = columns + rows
        for _ in range(rows):
            forward_moves.append(UP)
        for _ in range(columns):
            forward_moves.append(LEFT)
    else:
        # A near-identical pair: one front, kept whole, settles it before both ends' setup would pay.
        var near = DiagonalFronts()
        var close = diagonal_transition(profile, STEP_TENTHS_ALIGNMENT, near, switch_setup=TWO_ENDED_SETUP)
        if close.distance >= 0:
            trace_diagonals(profile, near, close.distance, forward_moves)
            return EditPath(forward_moves^, backward_moves^, columns, rows, close.distance)
        # Diagonal transition from both ends, keeping every front, while it is cheaper than a band:
        # where the fronts meet, the path is traced back to the start through the forward fronts
        # and on to the end through the backward ones.
        var first_back = reversed_codes(profile.column_codes, columns, FIRST_SENTINEL)
        var second_back = reversed_codes(profile.row_codes, rows, SECOND_SENTINEL)
        var ahead = FrontPair(record=True)
        var behind = FrontPair(record=True)
        var meeting = two_ended(profile, first_back, second_back, STEP_TENTHS_ALIGNMENT, ahead, behind)
        if meeting.probe.distance >= 0:
            var middle_column = meeting.column
            var middle_row = meeting.column - meeting.diagonal
            trace_from(
                ahead.history,
                columns,
                rows,
                meeting.forward_score,
                meeting.diagonal,
                middle_column,
                forward_moves,
            )
            # The backward fronts' trace runs from the meeting cell, mirrored, to the end, and its
            # moves right to left over the reversed sequences are the suffix left to right.
            backward_moves = List[UInt8](capacity=columns + rows - middle_column - middle_row)
            trace_from(
                behind.history,
                columns,
                rows,
                meeting.backward_score,
                (columns - rows) - meeting.diagonal,
                columns - middle_column,
                backward_moves,
            )
            return EditPath(forward_moves^, backward_moves^, middle_column, middle_row, meeting.probe.distance)
        var probe = meeting.probe
        var trusted = True
        var heuristic = band_start(profile, meeting.probe, probe, trusted)
        var trail = Trail(columns)
        distance = band_doubling[True](profile, False, probe, trail, heuristic, trusted)
        trace_back(profile, trail, columns, rows, distance, forward_moves)
    return EditPath(forward_moves^, backward_moves^, columns, rows, distance)
