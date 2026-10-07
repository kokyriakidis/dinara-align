# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Local alignment under `Costs` with a match reward, Smith-Waterman: the best-scoring alignment of any
part of the reference against any part of the query, as abPOA's local mode finds it.

It is ends-free alignment with all four ends free, and a reward for every match: with costs alone the
empty alignment would always win. A wavefront grows by cost, and with a reward and both ends floating
the cost no longer bounds the score, which is why WFA2-lib offers no local mode and KSW2 a score
alone. So the end is found by Smith-Waterman's own sweep, every cell's best score floored at zero, by
anti-diagonal sixteen cells at a time, as abPOA and `vector_score` vectorize it (see `best_end`), and the alignment by an extension back from that end (see
`gap_affine.extension_of`): the best alignment ending there is the local one, exact, its CIGAR chosen
by the same rule for ties as every other mode's.
"""

from .errors import AlignmentError
from .gap_affine import extension_of, extension_penalties
from .modes import Alignment, Anchor, Band, Costs, Ties


comptime WIDTH = 16
"""Cells a step of the sweep computes at once, along one anti-diagonal."""


def best_end[
    pieces: Int
](reference: Span[UInt8, _], query: Span[UInt8, _], costs: Costs, match_score: Int) -> Tuple[Int, Int, Int]:
    """The best local score, and the reference's and the query's letters up to where an alignment
    earning it ends: of several such ends the furthest along the reference, then the query.

    Gotoh's recurrence over scores, every cell's best floored at zero, where a local alignment may
    start, swept by anti-diagonal as `vector_score` sweeps it: every cell of `d = i + j` reads only
    diagonals `d - 1` and `d - 2`, so a diagonal's cells fill vector lanes, each lane keeping the best
    cell it has seen by the rule above, and the lanes compared once at the end."""
    comptime Lanes = SIMD[DType.int32, WIDTH]
    comptime LOW = Int32.MIN // 4
    var rows = len(reference)
    var columns = len(query)
    if rows == 0 or columns == 0:
        return (0, 0, 0)
    comptime two = pieces == 2
    # Lane `i` reads the reference's letter `i - 1`, stored one place on, and the query back to front,
    # so a diagonal's letters load contiguously; both padded past their ends by bytes no text holds.
    var letters = List[UInt8](length=rows + 1 + WIDTH, fill=0xFE)
    for index in range(rows):
        letters[index + 1] = reference[index]
    var reversed = List[UInt8](length=columns + WIDTH, fill=0xFF)
    for index in range(columns):
        reversed[index] = query[columns - 1 - index]
    var size = rows + 1 + WIDTH
    var two_back = List[Int32](length=size, fill=0)
    var one_back = List[Int32](length=size, fill=0)
    var current = List[Int32](length=size, fill=0)
    # A gap of reference letters alone grows down the lanes, one of query letters across diagonals.
    var deletes_back = List[Int32](length=size, fill=LOW)
    var deletes = List[Int32](length=size, fill=LOW)
    var inserts_back = List[Int32](length=size, fill=LOW)
    var inserts = List[Int32](length=size, fill=LOW)
    var deletes2_back = List[Int32](length=size, fill=LOW)
    var deletes2 = List[Int32](length=size, fill=LOW)
    var inserts2_back = List[Int32](length=size, fill=LOW)
    var inserts2 = List[Int32](length=size, fill=LOW)
    # A gap's first letter pays the opening and its extension, each further letter the extension.
    var first = Lanes(Int32(costs.opening + costs.extension))
    var further = Lanes(Int32(costs.extension))
    var first2 = Lanes(Int32(costs.opening2 + costs.extension2) if two else 0)
    var further2 = Lanes(Int32(costs.extension2) if two else 0)
    var matched = Lanes(Int32(match_score))
    var mismatched = Lanes(-Int32(costs.mismatch))
    var zero = Lanes(0)
    var lane_rows = Lanes()
    comptime for lane in range(WIDTH):
        lane_rows[lane] = Int32(lane)
    var best = Lanes(0)
    var best_row = Lanes(0)
    var best_diagonal = Lanes(0)
    for diagonal in range(2, rows + columns + 1):
        var low = max(1, diagonal - columns)
        var high = min(rows, diagonal - 1)
        var lag = columns - diagonal
        var row = low
        while row <= high:
            var above = one_back.unsafe_ptr().unsafe_offset(row - 1).unsafe_load[width=WIDTH]()
            var left = one_back.unsafe_ptr().unsafe_offset(row).unsafe_load[width=WIDTH]()
            var above_left = two_back.unsafe_ptr().unsafe_offset(row - 1).unsafe_load[width=WIDTH]()
            var mine = letters.unsafe_ptr().unsafe_offset(row).unsafe_load[width=WIDTH]()
            var theirs = reversed.unsafe_ptr().unsafe_offset(lag + row).unsafe_load[width=WIDTH]()
            var deletion = max(
                above - first, deletes_back.unsafe_ptr().unsafe_offset(row - 1).unsafe_load[width=WIDTH]() - further
            )
            var insertion = max(
                left - first, inserts_back.unsafe_ptr().unsafe_offset(row).unsafe_load[width=WIDTH]() - further
            )
            var score = max(above_left + mine.eq(theirs).select(matched, mismatched), max(deletion, insertion))
            deletes.unsafe_ptr().unsafe_offset(row).unsafe_store(deletion)
            inserts.unsafe_ptr().unsafe_offset(row).unsafe_store(insertion)
            comptime if two:
                var deletion2 = max(
                    above - first2,
                    deletes2_back.unsafe_ptr().unsafe_offset(row - 1).unsafe_load[width=WIDTH]() - further2,
                )
                var insertion2 = max(
                    left - first2, inserts2_back.unsafe_ptr().unsafe_offset(row).unsafe_load[width=WIDTH]() - further2
                )
                score = max(score, max(deletion2, insertion2))
                deletes2.unsafe_ptr().unsafe_offset(row).unsafe_store(deletion2)
                inserts2.unsafe_ptr().unsafe_offset(row).unsafe_store(insertion2)
            score = max(score, zero)
            current.unsafe_ptr().unsafe_offset(row).unsafe_store(score)
            # The furthest along the reference, then the query, of a lane's best cells: diagonals come in
            # order, so on the same row a later cell is further along the query. Lanes past the
            # diagonal's last row hold no cell.
            var at_row = lane_rows + Int32(row)
            var better = (score.gt(best) | (score.eq(best) & at_row.ge(best_row))) & at_row.le(Lanes(Int32(high)))
            best = better.select(score, best)
            best_row = better.select(at_row, best_row)
            best_diagonal = better.select(Lanes(Int32(diagonal)), best_diagonal)
            row += WIDTH
        # The border cells of this diagonal, written after the lanes that may have run over them: the
        # first row and column score zero, and no gap runs along them.
        if diagonal <= columns:
            current[0] = 0
            deletes[0] = LOW
            deletes2[0] = LOW
        if diagonal <= rows:
            current[diagonal] = 0
            inserts[diagonal] = LOW
            inserts2[diagonal] = LOW
        swap(two_back, one_back)
        swap(one_back, current)
        swap(deletes_back, deletes)
        swap(inserts_back, inserts)
        comptime if two:
            swap(deletes2_back, deletes2)
            swap(inserts2_back, inserts2)
    var top = Int32(0)
    var top_row = 0
    var top_column = 0
    comptime for lane in range(WIDTH):
        var value = best[lane]
        var at_row = Int(best_row[lane])
        var at_column = Int(best_diagonal[lane]) - at_row
        if value > top or (
            value == top and value > 0 and (at_row > top_row or (at_row == top_row and at_column > top_column))
        ):
            top = value
            top_row = at_row
            top_column = at_column
    return (Int(top), top_row, top_column)


def reversed_text(text: String) -> String:
    var bytes = text.as_bytes()
    var out = List[UInt8](capacity=len(bytes))
    for index in range(len(bytes) - 1, -1, -1):
        out.append(bytes[index])
    return String(unsafe_from_utf8=out^)


def reversed_cigar(cigar: String) -> String:
    """A CIGAR's runs in the other order: the alignment of both sequences reversed."""
    var bytes = cigar.as_bytes()
    var runs = List[String]()
    var start = 0
    for index in range(len(bytes)):
        if bytes[index] < UInt8(ord("0")) or bytes[index] > UInt8(ord("9")):
            runs.append(String(StringSlice(unsafe_from_utf8=bytes[start : index + 1])))
            start = index + 1
    var out = String()
    for index in range(len(runs) - 1, -1, -1):
        out += runs[index]
    return out


def local_alignment(
    reference: String, query: String, costs: Costs, match_score: Int, ties: Ties, extended: Bool
) raises AlignmentError -> Alignment:
    """The best local alignment, a match earning `match_score` (see `Mode.local`).

    `Ties.RIGHT` ends it as late as an equally good alignment allows (see `best_end`), starts it as
    late too, and spells its CIGAR by WFA2-lib's rule; `Ties.LEFT` is that over both sequences
    reversed, read backwards, as for every other mode: everything as early as it goes."""
    var columns = reference.byte_length()
    var rows = query.byte_length()
    if ties == Ties.LEFT:
        var mirrored = local_alignment(
            reversed_text(reference), reversed_text(query), costs, match_score, Ties.RIGHT, extended
        )
        if mirrored.score == 0:
            return mirrored^
        return Alignment(
            mirrored.cost,
            mirrored.score,
            reversed_cigar(mirrored.cigar),
            columns - mirrored.reference_end,
            columns - mirrored.reference_start,
            rows - mirrored.query_end,
            rows - mirrored.query_start,
        )
    var found = best_end[2](
        reference.as_bytes(), query.as_bytes(), costs, match_score
    ) if costs.pieces() == 2 else best_end[1](reference.as_bytes(), query.as_bytes(), costs, match_score)
    if found[0] == 0:
        return Alignment(0, 0, String(), 0, 0, 0, 0)
    var end_column = found[1]
    var end_row = found[2]
    var two = costs.pieces() == 2
    var penalties = extension_penalties(
        match_score,
        costs.mismatch,
        costs.opening,
        costs.extension,
        costs.opening2 if two else 0,
        costs.extension2 if two else 0,
    )
    var head = String(StringSlice(unsafe_from_utf8=reference.as_bytes()[:end_column]))
    var lead = String(StringSlice(unsafe_from_utf8=query.as_bytes()[:end_row]))
    # The best alignment ending at that cell, and starting wherever pays: an extension back from it.
    var back = extension_of[2](
        head, lead, penalties, extended, Anchor.END, Band(), Ties.RIGHT
    ) if two else extension_of[1](head, lead, penalties, extended, Anchor.END, Band(), Ties.RIGHT)
    return Alignment(
        match_score * back.matches - back.score,
        back.score,
        back.cigar,
        end_column - back.first_length,
        end_column,
        end_row - back.second_length,
        end_row,
    )
