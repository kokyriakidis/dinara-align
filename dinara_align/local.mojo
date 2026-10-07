# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Local and overlap alignment under `Costs` with a match reward. Local, Smith-Waterman: the best-scoring alignment of any
part of the reference against any part of the query, as abPOA's local mode finds it.

It is ends-free alignment with all four ends free, and a reward for every match: with costs alone the
empty alignment would always win. A wavefront grows by cost, and with a reward and both ends floating
the cost no longer bounds the score, which is why WFA2-lib offers no local mode and KSW2 a score
alone. So the end is found by Smith-Waterman's own sweep, every cell's best score floored at zero, by
anti-diagonal in 16-bit lanes along the shorter sequence while the scores fit (see `best_end`), and
the alignment by an extension back from that end, which stops once it earns the sweep's score (see
`gap_affine.extension_of`): the best alignment ending there is the local one, exact, its CIGAR chosen
by the same rule for ties as every other mode's.
"""

from .errors import AlignmentError
from .gap_affine import AffineExtension, EndsFree, cigar_within, extension_of, extension_penalties, traced_extension
from .modes import Alignment, Anchor, Band, Costs, Ties


comptime ANYWHERE = 0
"""A sweep for a local alignment: it starts at any cell, every cell's best floored at zero, and ends at
any cell."""
comptime FROM_EDGE = 1
"""A sweep for an overlap: it starts anywhere on the first row or column, for nothing, and ends on
the last row or column."""
comptime FROM_ORIGIN = 2
"""A sweep back from an overlap's end: it starts at the origin, a gap along either edge paid, and ends
on the last row or column."""


def best_end[
    pieces: Int, dtype: DType, width: Int, transposed: Bool, kind: Int = ANYWHERE
](reference: Span[UInt8, _], query: Span[UInt8, _], costs: Costs, match_score: Int) -> Tuple[Int, Int, Int]:
    """The best score of an alignment starting and ending where `kind` allows, and the reference's and
    the query's letters up to where it ends: of several such ends the furthest along both together,
    then along the reference.

    Gotoh's recurrence over scores, a local alignment's every cell floored at zero, swept by
    anti-diagonal as `vector_score` sweeps it: every cell of `d = i + j` reads only diagonals `d - 1`
    and `d - 2`, so a diagonal's cells fill `width` lanes of `dtype`, as narrow as the scores allow
    (see `local_alignment`). For a local alignment each lane keeps its best of the diagonal and the
    step it came at, and the diagonal's best is placed only when it meets the best so far; an overlap
    ends on the last row or column, at most two cells a diagonal, read once it is done."""
    comptime Lanes = SIMD[dtype, width]
    comptime Value = Scalar[dtype]
    comptime LOW = Value.MIN // 4
    # The lanes run along the shorter sequence, so the sweep's rows stay in the core's own cache:
    # along the reference, or with `transposed` along the query.
    var down_letters = query if transposed else reference
    var across_letters = reference if transposed else query
    var rows = len(down_letters)
    var columns = len(across_letters)
    if rows == 0 or columns == 0:
        return (0, 0, 0)
    comptime two = pieces == 2
    comptime local = kind == ANYWHERE
    # Lane `i` reads the reference's letter `i - 1`, stored one place on, and the query back to front,
    # so a diagonal's letters load contiguously; both padded past their ends by bytes no text holds.
    var letters = List[UInt8](length=rows + 1 + width, fill=0xFE)
    for index in range(rows):
        letters[index + 1] = down_letters[index]
    var reversed = List[UInt8](length=columns + width, fill=0xFF)
    for index in range(columns):
        reversed[index] = across_letters[columns - 1 - index]
    var size = rows + 1 + width
    var two_back = List[Value](length=size, fill=0)
    var one_back = List[Value](length=size, fill=0)
    var current = List[Value](length=size, fill=0)
    # A gap of reference letters alone grows down the lanes, one of query letters across diagonals.
    var deletes_back = List[Value](length=size, fill=LOW)
    var deletes = List[Value](length=size, fill=LOW)
    var inserts_back = List[Value](length=size, fill=LOW)
    var inserts = List[Value](length=size, fill=LOW)
    var deletes2_back = List[Value](length=size if two else 0, fill=LOW)
    var deletes2 = List[Value](length=size if two else 0, fill=LOW)
    var inserts2_back = List[Value](length=size if two else 0, fill=LOW)
    var inserts2 = List[Value](length=size if two else 0, fill=LOW)
    # A gap's first letter pays the opening and its extension, each further letter the extension.
    var first = Lanes(Value(costs.opening + costs.extension))
    var further = Lanes(Value(costs.extension))
    var first2 = Lanes(Value(costs.opening2 + costs.extension2) if two else 0)
    var further2 = Lanes(Value(costs.extension2) if two else 0)
    var matched = Lanes(Value(match_score))
    var mismatched = Lanes(-Value(costs.mismatch))
    var zero = Lanes(0)
    var lane_index = Lanes()
    comptime for lane in range(width):
        lane_index[lane] = Value(lane)

    @inline(.always)
    def gap(letters: Int, opening: Int, extension: Int) -> Value:
        """A gap of `letters` along an edge, as a sweep from the origin pays it."""
        return -Value(opening + extension * letters)

    var best = 0 if local else Int.MIN
    var best_row = 0
    var best_diagonal = 0

    @inline(.always)
    def ending(
        values: List[Value], diagonal: Int, mut best: Int, mut best_row: Int, mut best_diagonal: Int
    ) {imm rows, imm columns}:
        """An overlap's ends on `diagonal`: its cells on the last row and the last column, the later
        diagonal, then the further along the reference, on a tie."""
        for edge in range(2):
            var row = rows if edge == 0 else diagonal - columns
            if row < 0 or row > rows or diagonal - row < 0 or diagonal - row > columns:
                continue
            var value = Int(values[row])
            var further = (row < best_row) if transposed else (row > best_row)
            if value > best or (value == best and (diagonal > best_diagonal or further)):
                best = value
                best_row = row
                best_diagonal = diagonal

    comptime if kind == FROM_ORIGIN:
        # Diagonal one, beside the origin: a letter of either sequence against a gap.
        one_back[0] = max(
            gap(1, costs.opening, costs.extension), gap(1, costs.opening2, costs.extension2) if two else LOW
        )
        one_back[1] = one_back[0]
        inserts_back[0] = gap(1, costs.opening, costs.extension)
        deletes_back[1] = gap(1, costs.opening, costs.extension)
        comptime if two:
            inserts2_back[0] = gap(1, costs.opening2, costs.extension2)
            deletes2_back[1] = gap(1, costs.opening2, costs.extension2)
    comptime if not local:
        ending(one_back, 1, best, best_row, best_diagonal)
    for diagonal in range(2, rows + columns + 1):
        var low = max(1, diagonal - columns)
        var high = min(rows, diagonal - 1)
        var lag = columns - diagonal
        var row = low
        # Each lane's best on this diagonal, and the step it came at, the later on a tie.
        var top = zero
        var top_step = zero
        var step = Value(0)
        while row <= high:
            var above = one_back.unsafe_ptr().unsafe_offset(row - 1).unsafe_load[width=width]()
            var left = one_back.unsafe_ptr().unsafe_offset(row).unsafe_load[width=width]()
            var above_left = two_back.unsafe_ptr().unsafe_offset(row - 1).unsafe_load[width=width]()
            var mine = letters.unsafe_ptr().unsafe_offset(row).unsafe_load[width=width]()
            var theirs = reversed.unsafe_ptr().unsafe_offset(lag + row).unsafe_load[width=width]()
            var deletion = max(
                above - first, deletes_back.unsafe_ptr().unsafe_offset(row - 1).unsafe_load[width=width]() - further
            )
            var insertion = max(
                left - first, inserts_back.unsafe_ptr().unsafe_offset(row).unsafe_load[width=width]() - further
            )
            var score = max(above_left + mine.eq(theirs).select(matched, mismatched), max(deletion, insertion))
            deletes.unsafe_ptr().unsafe_offset(row).unsafe_store(deletion)
            inserts.unsafe_ptr().unsafe_offset(row).unsafe_store(insertion)
            comptime if two:
                var deletion2 = max(
                    above - first2,
                    deletes2_back.unsafe_ptr().unsafe_offset(row - 1).unsafe_load[width=width]() - further2,
                )
                var insertion2 = max(
                    left - first2, inserts2_back.unsafe_ptr().unsafe_offset(row).unsafe_load[width=width]() - further2
                )
                score = max(score, max(deletion2, insertion2))
                deletes2.unsafe_ptr().unsafe_offset(row).unsafe_store(deletion2)
                inserts2.unsafe_ptr().unsafe_offset(row).unsafe_store(insertion2)
            comptime if local:
                score = max(score, zero)
            current.unsafe_ptr().unsafe_offset(row).unsafe_store(score)
            comptime if local:
                # Lanes past the diagonal's last row hold no cell.
                var counted = score
                if row + width - 1 > high:
                    counted = lane_index.lt(Lanes(Value(high - row + 1))).select(score, zero)
                var later = counted.gt(top) if transposed else counted.ge(top)
                top = later.select(counted, top)
                top_step = later.select(Lanes(step), top_step)
                step += 1
            row += width
        comptime if local:
            # The diagonal's best, placed at its furthest row, when it meets the best so far.
            var most = Int(top.reduce_max())
            if most > 0 and most >= best:
                var furthest = 0 if not transposed else Int.MAX
                comptime for lane in range(width):
                    if Int(top[lane]) == most:
                        var at = low + Int(top_step[lane]) * width + lane
                        furthest = min(furthest, at) if transposed else max(furthest, at)
                best = most
                best_row = furthest
                best_diagonal = diagonal
        # The border cells of this diagonal, written after the lanes that may have run over them: the
        # first row and column score zero, and no gap runs along them, or from the origin a gap does.
        if diagonal <= columns:
            comptime if kind == FROM_ORIGIN:
                var along = gap(diagonal, costs.opening, costs.extension)
                current[0] = along
                inserts[0] = along
                comptime if two:
                    var along2 = gap(diagonal, costs.opening2, costs.extension2)
                    current[0] = max(along, along2)
                    inserts2[0] = along2
            else:
                current[0] = 0
            deletes[0] = LOW
            comptime if two:
                deletes2[0] = LOW
        if diagonal <= rows:
            comptime if kind == FROM_ORIGIN:
                var down = gap(diagonal, costs.opening, costs.extension)
                current[diagonal] = down
                deletes[diagonal] = down
                comptime if two:
                    var down2 = gap(diagonal, costs.opening2, costs.extension2)
                    current[diagonal] = max(down, down2)
                    deletes2[diagonal] = down2
            else:
                current[diagonal] = 0
            inserts[diagonal] = LOW
            comptime if two:
                inserts2[diagonal] = LOW
        comptime if not local:
            ending(current, diagonal, best, best_row, best_diagonal)
        swap(two_back, one_back)
        swap(one_back, current)
        swap(deletes_back, deletes)
        swap(inserts_back, inserts)
        comptime if two:
            swap(deletes2_back, deletes2)
            swap(inserts2_back, inserts2)
    # `best_row` counts the lanes' sequence, the rest of the diagonal the other's.
    if transposed:
        return (best, best_diagonal - best_row, best_row)
    return (best, best_row, best_diagonal - best_row)


def narrow_enough[kind: Int](costs: Costs, match_score: Int, rows: Int, columns: Int) -> Bool:
    """Whether every score of the sweep fits 16 bits. None passes the reward of the shorter sequence
    matched throughout. A local alignment's none falls further below zero than the dearest single
    move, as each cell takes the best of its moves from cells of zero or more; an overlap's none falls
    below every letter of both paying the dearest move. Nor does a gap's sentinel, a quarter of the
    way down, fall further than one extension below it."""
    var dearest = max(costs.mismatch, costs.opening + costs.extension)
    if costs.pieces() == 2:
        dearest = max(dearest, costs.opening2 + costs.extension2)
    var fits = match_score * (min(rows, columns) + 1) < 32000 and dearest < 4000
    comptime if kind != ANYWHERE:
        fits = fits and dearest * (rows + columns + 1) < 8000
    return fits


def swept[
    kind: Int
](reference: Span[UInt8, _], query: Span[UInt8, _], costs: Costs, match_score: Int) -> Tuple[Int, Int, Int]:
    """`best_end` with its lanes along the shorter sequence, 16 bits to a lane, thirty-two to an AVX-512
    register, while the scores fit, else 32."""
    var transposed = len(query) < len(reference)
    if narrow_enough[kind](costs, match_score, len(reference), len(query)):
        if costs.pieces() == 2:
            if transposed:
                return best_end[2, DType.int16, 32, True, kind](reference, query, costs, match_score)
            return best_end[2, DType.int16, 32, False, kind](reference, query, costs, match_score)
        if transposed:
            return best_end[1, DType.int16, 32, True, kind](reference, query, costs, match_score)
        return best_end[1, DType.int16, 32, False, kind](reference, query, costs, match_score)
    if costs.pieces() == 2:
        if transposed:
            return best_end[2, DType.int32, 16, True, kind](reference, query, costs, match_score)
        return best_end[2, DType.int32, 16, False, kind](reference, query, costs, match_score)
    if transposed:
        return best_end[1, DType.int32, 16, True, kind](reference, query, costs, match_score)
    return best_end[1, DType.int32, 16, False, kind](reference, query, costs, match_score)


def end_of(reference: Span[UInt8, _], query: Span[UInt8, _], costs: Costs, match_score: Int) -> Tuple[Int, Int, Int]:
    """A local alignment's best score and where it ends (see `best_end`)."""
    return swept[ANYWHERE](reference, query, costs, match_score)


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

    `Ties.LEFT` ends it as late as an equally good alignment allows (see `best_end`), starts it as
    late too, the shortest, and spells its CIGAR by the left rule, the extension back from its end
    traced as it searched (see `gap_affine.traced_extension`); `Ties.RIGHT` is that over both
    sequences reversed, read backwards, as for every other mode: everything as early as it goes, gaps
    right."""
    var columns = reference.byte_length()
    var rows = query.byte_length()
    if ties == Ties.RIGHT:
        var mirrored = local_alignment(
            reversed_text(reference), reversed_text(query), costs, match_score, Ties.LEFT, extended
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
    var found = end_of(reference.as_bytes(), query.as_bytes(), costs, match_score)
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
    # The best alignment ending at that cell, and starting wherever pays: an extension back from it,
    # which stops on earning the sweep's score, the best any alignment ending there earns.
    var traced = traced_extension[2](head, lead, penalties, extended, found[0]) if two else traced_extension[1](
        head, lead, penalties, extended, found[0]
    )
    var back: AffineExtension
    if traced:
        back = traced.take()
        back.cigar = reversed_cigar(back.cigar)
    else:
        # Too many fronts to keep: the extension's own search and split, by the left rule too, over
        # both sequences reversed from the end.
        var mirrored = extension_of[2](
            reversed_text(head), reversed_text(lead), penalties, extended, Anchor.START, Band(), Ties.RIGHT, found[0]
        ) if two else extension_of[1](
            reversed_text(head), reversed_text(lead), penalties, extended, Anchor.START, Band(), Ties.RIGHT, found[0]
        )
        mirrored.cigar = reversed_cigar(mirrored.cigar)
        back = mirrored^
    return Alignment(
        match_score * back.matches - back.score,
        back.score,
        back.cigar,
        end_column - back.first_length,
        end_column,
        end_row - back.second_length,
        end_row,
    )


def overlap_alignment(
    reference: String, query: String, costs: Costs, match_score: Int, ties: Ties, extended: Bool
) raises AlignmentError -> Alignment:
    """The best overlap, a match earning `match_score` (see `Mode.overlap`): every end gap free, so the
    alignment starts on the first row or column and ends on the last.

    One sweep from the edges finds the best score and its end, the furthest along both together, then
    along the reference; one back from that end over both sequences reversed, paying gaps from it,
    finds the start, chosen the same way there, the earliest; between them the alignment is a global
    one with a reward, whose letters are fixed, so the reward folds into the costs (see `gap_affine`)
    and the wavefront aligns it, its CIGAR the one `ties` names."""
    var forward = swept[FROM_EDGE](reference.as_bytes(), query.as_bytes(), costs, match_score)
    var end_column = forward[1]
    var end_row = forward[2]
    if forward[0] <= 0 or end_column == 0 or end_row == 0:
        return Alignment(0, 0, String(), end_column, end_column, end_row, end_row)
    var head = reversed_text(String(StringSlice(unsafe_from_utf8=reference.as_bytes()[:end_column])))
    var lead = reversed_text(String(StringSlice(unsafe_from_utf8=query.as_bytes()[:end_row])))
    var back = swept[FROM_ORIGIN](head.as_bytes(), lead.as_bytes(), costs, match_score)
    var start_column = end_column - back[1]
    var start_row = end_row - back[2]
    var two = costs.pieces() == 2
    var penalties = extension_penalties(
        match_score,
        costs.mismatch,
        costs.opening,
        costs.extension,
        costs.opening2 if two else 0,
        costs.extension2 if two else 0,
    )
    var part = String(StringSlice(unsafe_from_utf8=reference.as_bytes()[start_column:end_column]))
    var piece = String(StringSlice(unsafe_from_utf8=query.as_bytes()[start_row:end_row]))
    var found = cigar_within[2](
        part, piece, penalties, extended, Int.MAX, EndsFree(), Band(), ties
    ) if two else cigar_within[1](part, piece, penalties, extended, Int.MAX, EndsFree(), Band(), ties)
    var cigar = found.take().cigar
    var matches = matches_of(part, piece, cigar)
    return Alignment(
        match_score * matches - forward[0], forward[0], cigar^, start_column, end_column, start_row, end_row
    )


def matches_of(first: String, second: String, cigar: String) -> Int:
    """The equal pairs a CIGAR of `first` against `second` aligns, `M` runs compared letter by letter."""
    var a = first.as_bytes()
    var b = second.as_bytes()
    var column = 0
    var row = 0
    var length = 0
    var total = 0
    for byte in cigar.as_bytes():
        if byte >= UInt8(ord("0")) and byte <= UInt8(ord("9")):
            length = length * 10 + Int(byte - UInt8(ord("0")))
            continue
        if byte == UInt8(ord("D")):
            column += length
        elif byte == UInt8(ord("I")):
            row += length
        else:
            for _ in range(length):
                if a[column] == b[row]:
                    total += 1
                column += 1
                row += 1
        length = 0
    return total
