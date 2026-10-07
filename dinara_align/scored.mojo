"""
The modes that maximize a score under `Costs`, a match earning a reward: a local alignment,
Smith-Waterman, and free ends with a reward, semi-global as parasail and hyalite count it, an overlap
among them.

With a reward and an end floating, the cost no longer bounds the score, so a wavefront, which grows
by cost, cannot search them, which is why WFA2-lib offers no local mode and KSW2 a local score alone.
Each takes Smith-Waterman's own sweep instead, by anti-diagonal in 16-bit lanes along the shorter
sequence while the scores fit (see `best_end`), which finds the best score and where it ends, and
keeps a few rows. The alignment then comes from the wavefront, exact and in bounded memory, its CIGAR
by the same rule for ties as every other mode's: a local one by the extension back from its end,
which stops once it earns the sweep's score (see `local_alignment`), and one with free ends by a
sweep back from its end for its start, the span the one `Ties` names, the letters between them a global
alignment whose reward folds into the costs (see `rewarded_alignment`).
"""

from .cigar import cigar_cost, reversed_cigar, reversed_text
from .errors import AlignmentError
from .gap_affine import (
    AffineExtension,
    EndsFree,
    cigar_within,
    extension_of,
    extension_penalties,
    outside,
    traced_extension,
)
from .modes import Alignment, Anchor, Band, Costs, Ties


comptime ANYWHERE = 0
"""A sweep for a local alignment: it starts at any cell, every cell's best floored at zero, and ends at
any cell."""
comptime FROM_EDGE = 1
"""A sweep for ends-free alignment with a reward: it starts on the first row or column, for nothing
within the letters free there and past them paying a gap, and ends on the last row or column within
the letters free there. With none free at the start it starts at the origin."""


def best_end[
    pieces: Int, dtype: DType, width: Int, transposed: Bool, kind: Int = ANYWHERE
](
    reference: Span[UInt8, _],
    query: Span[UInt8, _],
    costs: Costs,
    match_score: Int,
    ends: EndsFree = EndsFree(),
    highest: Bool = True,
) -> Tuple[Int, Int, Int]:
    """The best score of an alignment starting and ending where `kind` allows, and where it ends (see
    `swept_cells`)."""
    var unused = List[Int32]()
    return swept_cells[pieces, dtype, width, transposed, kind, False](
        reference, query, costs, match_score, ends, highest, unused
    )


def swept_cells[
    pieces: Int, dtype: DType, width: Int, transposed: Bool, kind: Int, columns_kept: Bool
](
    reference: Span[UInt8, _],
    query: Span[UInt8, _],
    costs: Costs,
    match_score: Int,
    ends: EndsFree,
    highest: Bool,
    mut kept: List[Int32],
) -> Tuple[Int, Int, Int]:
    """The best score of an alignment starting and ending where `kind` allows, and the reference's and
    the query's letters up to where it ends: of several such ends, for a local alignment the furthest
    along both together, then along the reference; from the edges the one on the highest diagonal, the
    furthest along the reference less the query, or with `highest` false the lowest. With
    `columns_kept`, for a local alignment, `kept` holds after it the best score of every column: of
    the cells with as many reference letters, `len(reference) + 1` of them.

    Gotoh's recurrence over scores, a local alignment's every cell floored at zero, swept by
    anti-diagonal as `vector_score` sweeps it: every cell of `d = i + j` reads only diagonals `d - 1`
    and `d - 2`, so a diagonal's cells fill `width` lanes of `dtype`, as narrow as the scores allow
    (see `local_alignment`). For a local alignment each lane keeps its best of the diagonal and the
    step it came at, and the diagonal's best is placed only when it meets the best so far; one from
    the edges, `ends` saying which letters at each end are free, ends on the last row or column, at
    most two cells a diagonal, read once it is done."""
    comptime Lanes = SIMD[dtype, width]
    comptime Value = Scalar[dtype]
    comptime LOW = Value.MIN // 4
    # The lanes run along the shorter sequence, so the sweep's rows stay in the core's own cache:
    # along the reference, or with `transposed` along the query.
    var down_letters = query if transposed else reference
    var across_letters = reference if transposed else query
    var rows = len(down_letters)
    var columns = len(across_letters)
    comptime two = pieces == 2
    comptime local = kind == ANYWHERE
    # The free letters at either end of the lanes' sequence and the other's; a local alignment's start
    # is free everywhere.
    var down_start = Int.MAX if local else (ends.second_begin if transposed else ends.first_begin)
    var down_end = ends.second_end if transposed else ends.first_end
    var across_start = Int.MAX if local else (ends.first_begin if transposed else ends.second_begin)
    var across_end = ends.first_end if transposed else ends.second_end
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
    # Each column's best as the lanes find it: by the lanes' own row along the reference, or along the
    # query back from the reference's end, so a lane group's cells are always a run of the store.
    var column_best = List[Value](length=(len(reference) + 1 + width) if columns_kept else 0, fill=0)

    @inline(.always)
    def edge(letters: Int, free: Int) {imm costs} -> Value:
        """A cell on the first row or column, `letters` in: nothing within the `free` letters, past them
        a gap from where they end, at the cheaper piece."""
        if letters <= free:
            return 0
        var paid = -(costs.opening + costs.extension * (letters - free))
        comptime if two:
            paid = max(paid, -(costs.opening2 + costs.extension2 * (letters - free)))
        return Value(paid)

    if rows == 0 or columns == 0:
        comptime if local:
            comptime if columns_kept:
                kept = List[Int32](length=len(reference) + 1, fill=0)
            return (0, 0, 0)
        # One sequence empty: the alignment lies along the other's edge, ending within its free letters,
        # the furthest such end on a tie.
        var length = rows + columns
        var start = down_start if rows > 0 else across_start
        var finish = down_end if rows > 0 else across_end
        # Letters along the reference raise the diagonal, along the query lower it.
        var along_reference = (rows > 0) != transposed
        var top = Int.MIN
        var reach = 0
        for letters in range(length + 1):
            if length - letters > finish:
                continue
            var value = Int(edge(letters, start))
            if value > top or (value == top and (along_reference == highest)):
                top = value
                reach = letters
        var down = reach if rows > 0 else 0
        var across = reach if columns > 0 else 0
        if transposed:
            return (top, across, down)
        return (top, down, across)

    var best = 0 if local else Int.MIN
    var best_row = 0
    var best_diagonal = 0

    @inline(.always)
    def ending(
        values: List[Value], diagonal: Int, mut best: Int, mut best_row: Int, mut best_diagonal: Int
    ) {imm rows, imm columns, imm across_end, imm down_end, imm highest}:
        """The ends on anti-diagonal `diagonal`: its cells on the last row and the last column within the
        letters free there, the one on the highest diagonal, the reference's letters less the query's, or
        with `highest` false the lowest, on a tie."""
        for side in range(2):
            var row = rows if side == 0 else diagonal - columns
            if row < 0 or row > rows or diagonal - row < 0 or diagonal - row > columns:
                continue
            if side == 0 and columns - (diagonal - row) > across_end:
                continue
            if side == 1 and rows - row > down_end:
                continue
            var value = Int(values[row])
            # The reference's letters less the query's, from the lanes' row and the rest of the anti-diagonal.
            var lean = (diagonal - 2 * row) if transposed else (2 * row - diagonal)
            var best_lean = (best_diagonal - 2 * best_row) if transposed else (2 * best_row - best_diagonal)
            var beyond = (lean > best_lean) if highest else (lean < best_lean)
            if value > best or (value == best and beyond):
                best = value
                best_row = row
                best_diagonal = diagonal

    # Diagonal one, beside the origin: a letter of either sequence against nothing.
    one_back[0] = edge(1, across_start)
    one_back[1] = edge(1, down_start)
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
                comptime if columns_kept:
                    var at = column_best.unsafe_ptr().unsafe_offset((columns - diagonal + row) if transposed else row)
                    at.unsafe_store(max(at.unsafe_load[width=width](), counted))
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
        # The border cells of this diagonal, written after the lanes that may have run over them: free
        # within the letters free there, past them a gap. No interior cell reads a border's gap layers.
        if diagonal <= columns:
            current[0] = edge(diagonal, across_start)
            deletes[0] = LOW
            comptime if two:
                deletes2[0] = LOW
        if diagonal <= rows:
            current[diagonal] = edge(diagonal, down_start)
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
    comptime if columns_kept:
        var length = len(reference)
        kept = List[Int32](length=length + 1, fill=0)
        for letters in range(length + 1):
            kept[letters] = Int32(column_best[(columns - letters) if transposed else letters])
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
](
    reference: Span[UInt8, _],
    query: Span[UInt8, _],
    costs: Costs,
    match_score: Int,
    ends: EndsFree = EndsFree(),
    highest: Bool = True,
) -> Tuple[Int, Int, Int]:
    """`best_end` with its lanes along the shorter sequence, 16 bits to a lane, thirty-two to an AVX-512
    register, while the scores fit, else 32."""
    var transposed = len(query) < len(reference)
    if narrow_enough[kind](costs, match_score, len(reference), len(query)):
        if costs.pieces() == 2:
            if transposed:
                return best_end[2, DType.int16, 32, True, kind](reference, query, costs, match_score, ends, highest)
            return best_end[2, DType.int16, 32, False, kind](reference, query, costs, match_score, ends, highest)
        if transposed:
            return best_end[1, DType.int16, 32, True, kind](reference, query, costs, match_score, ends, highest)
        return best_end[1, DType.int16, 32, False, kind](reference, query, costs, match_score, ends, highest)
    if costs.pieces() == 2:
        if transposed:
            return best_end[2, DType.int32, 16, True, kind](reference, query, costs, match_score, ends, highest)
        return best_end[2, DType.int32, 16, False, kind](reference, query, costs, match_score, ends, highest)
    if transposed:
        return best_end[1, DType.int32, 16, True, kind](reference, query, costs, match_score, ends, highest)
    return best_end[1, DType.int32, 16, False, kind](reference, query, costs, match_score, ends, highest)


@fieldwise_init
struct LocalScores(ImplicitlyCopyable, Writable):
    """A local alignment's best score and where it ends, and the best score of an alignment ending
    elsewhere, as SSW reports them for a mapping quality: `second_score` the best of any cell whose
    reference letters lie more than the window from `reference_end`, at `second_reference_end`, the
    first such column; zero, at zero, when there is none."""

    var score: Int
    var reference_end: Int
    var query_end: Int
    var second_score: Int
    var second_reference_end: Int


def local_scores(
    reference: Span[UInt8, _], query: Span[UInt8, _], costs: Costs, match_score: Int, window: Int
) -> LocalScores:
    """The best local score and where it ends, by the sweep `end_of` runs, and the best of every
    column more than `window` reference letters from that end (see `LocalScores`), from each column's
    best the same sweep keeps."""
    var kept = List[Int32]()
    var found: Tuple[Int, Int, Int]
    var flipped = len(query) < len(reference)
    var two = costs.pieces() == 2
    var none = EndsFree()
    # As `swept` dispatches: lanes along the shorter sequence, 16 bits while the scores fit.
    if narrow_enough[ANYWHERE](costs, match_score, len(reference), len(query)):
        if two:
            found = swept_cells[2, DType.int16, 32, True, ANYWHERE, True](
                reference, query, costs, match_score, none, True, kept
            ) if flipped else swept_cells[2, DType.int16, 32, False, ANYWHERE, True](
                reference, query, costs, match_score, none, True, kept
            )
        else:
            found = swept_cells[1, DType.int16, 32, True, ANYWHERE, True](
                reference, query, costs, match_score, none, True, kept
            ) if flipped else swept_cells[1, DType.int16, 32, False, ANYWHERE, True](
                reference, query, costs, match_score, none, True, kept
            )
    elif two:
        found = swept_cells[2, DType.int32, 16, True, ANYWHERE, True](
            reference, query, costs, match_score, none, True, kept
        ) if flipped else swept_cells[2, DType.int32, 16, False, ANYWHERE, True](
            reference, query, costs, match_score, none, True, kept
        )
    else:
        found = swept_cells[1, DType.int32, 16, True, ANYWHERE, True](
            reference, query, costs, match_score, none, True, kept
        ) if flipped else swept_cells[1, DType.int32, 16, False, ANYWHERE, True](
            reference, query, costs, match_score, none, True, kept
        )
    var second = 0
    var second_end = 0
    for letters in range(len(kept)):
        if abs(letters - found[1]) <= window:
            continue
        if Int(kept[letters]) > second:
            second = Int(kept[letters])
            second_end = letters
    return LocalScores(found[0], found[1], found[2], second, second_end)


def end_of(reference: Span[UInt8, _], query: Span[UInt8, _], costs: Costs, match_score: Int) -> Tuple[Int, Int, Int]:
    """A local alignment's best score and where it ends (see `best_end`)."""
    return swept[ANYWHERE](reference, query, costs, match_score)


def latest_local(
    reference: String, query: String, costs: Costs, match_score: Int, extended: Bool, limit: Int
) raises AlignmentError -> Alignment:
    """The best local alignment ending as late as an equally good one allows (see `best_end`) and
    starting as late too, the shortest, its CIGAR by the left rule: the extension back from its end,
    traced as it searched (see `gap_affine.traced_extension`)."""
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
    var traced = traced_extension[2](head, lead, penalties, extended, found[0], limit) if two else traced_extension[1](
        head, lead, penalties, extended, found[0], limit
    )
    var back: AffineExtension
    if traced:
        back = traced.take()
        back.cigar = reversed_cigar(back.cigar)
    else:
        # Too many fronts to keep: the extension's own search and split, by the left rule too, over
        # both sequences reversed from the end.
        var mirrored = extension_of[2](
            reversed_text(head),
            reversed_text(lead),
            penalties,
            extended,
            Anchor.START,
            Band(),
            Ties.RIGHT,
            found[0],
            limit,
        ) if two else extension_of[1](
            reversed_text(head),
            reversed_text(lead),
            penalties,
            extended,
            Anchor.START,
            Band(),
            Ties.RIGHT,
            found[0],
            limit,
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


def local_alignment(
    reference: String, query: String, costs: Costs, match_score: Int, ties: Ties, extended: Bool, limit: Int
) raises AlignmentError -> Alignment:
    """The best local alignment, a match earning `match_score` (see `Mode.local`).

    `Ties.LEFT` ends it as late as an equally good alignment allows and starts it as late too, the
    shortest, its CIGAR by the left rule (see `latest_local`); `Ties.RIGHT` is that over both sequences
    reversed, read backwards, as for every other mode: it starts and ends as early as it may, gaps
    right."""
    if ties == Ties.LEFT:
        return latest_local(reference, query, costs, match_score, extended, limit)
    var columns = reference.byte_length()
    var rows = query.byte_length()
    var mirrored = latest_local(reversed_text(reference), reversed_text(query), costs, match_score, extended, limit)
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


def rewarded_span(
    reference: String, query: String, costs: Costs, match_score: Int, ends: EndsFree
) -> Tuple[Int, Int, Int, Int, Int]:
    """The best score with `ends`' letters free and a match earning `match_score`, and the span the rule
    of `Ties.LEFT` names (see `gap_affine.free_ends_alignment`): its end on the highest diagonal an
    optimum reaches, by a sweep from the edges, and its start on the highest of those ending there, by
    a sweep back from that end over both sequences reversed, to the letters free at the start. The
    score, then the start's and the end's columns and rows."""
    var forward = swept[FROM_EDGE](reference.as_bytes(), query.as_bytes(), costs, match_score, ends)
    var end_column = forward[1]
    var end_row = forward[2]
    var head = reversed_text(String(StringSlice(unsafe_from_utf8=reference.as_bytes()[:end_column])))
    var lead = reversed_text(String(StringSlice(unsafe_from_utf8=query.as_bytes()[:end_row])))
    # From the end, the highest diagonal is the lowest of the reversed sequences'.
    var back = swept[FROM_EDGE](
        head.as_bytes(), lead.as_bytes(), costs, match_score, EndsFree(0, ends.first_begin, 0, ends.second_begin), False
    )
    return (forward[0], end_column - back[1], end_row - back[2], end_column, end_row)


def rewarded_alignment(
    reference: String,
    query: String,
    costs: Costs,
    match_score: Int,
    ends: EndsFree,
    ties: Ties,
    extended: Bool,
    limit: Int,
) raises AlignmentError -> Alignment:
    """The best alignment with the letters `ends` allows free at either end, a match earning
    `match_score` (see `Mode.ends_free`): semi-global with a reward, and with every end free an
    overlap.

    Its span is the one the rule of `ties` names, as for free ends without a reward (see
    `gap_affine.free_ends_alignment`): for `Ties.LEFT` found by sweeps (see `rewarded_span`), for
    `Ties.RIGHT` by the same over both sequences reversed. Between its ends the alignment is a global one
    with a reward, whose letters are fixed, so the reward folds into the costs (see `gap_affine`) and
    the wavefront aligns it, its CIGAR the one `ties` picks. A gap past the free letters at either end
    lies inside the span, as it is paid."""
    var columns = reference.byte_length()
    var rows = query.byte_length()
    var span: Tuple[Int, Int, Int, Int, Int]
    if ties == Ties.RIGHT:
        var mirrored = rewarded_span(
            reversed_text(reference),
            reversed_text(query),
            costs,
            match_score,
            EndsFree(ends.first_end, ends.first_begin, ends.second_end, ends.second_begin),
        )
        span = (mirrored[0], columns - mirrored[3], rows - mirrored[4], columns - mirrored[1], rows - mirrored[2])
    else:
        span = rewarded_span(reference, query, costs, match_score, ends)
    var start_column = span[1]
    var start_row = span[2]
    var end_column = span[3]
    var end_row = span[4]
    if start_column == end_column and start_row == end_row:
        return Alignment(0, span[0], String(), start_column, end_column, start_row, end_row)
    var part = String(StringSlice(unsafe_from_utf8=reference.as_bytes()[start_column:end_column]))
    var piece = String(StringSlice(unsafe_from_utf8=query.as_bytes()[start_row:end_row]))
    var found = global_rewarded(part, piece, costs, match_score, Band(), ties, extended, limit)
    return Alignment(found[0], span[0], found[1], start_column, end_column, start_row, end_row)


def global_rewarded(
    reference: String,
    query: String,
    costs: Costs,
    match_score: Int,
    band: Band,
    ties: Ties,
    extended: Bool,
    limit: Int,
) raises AlignmentError -> Tuple[Int, String]:
    """A global alignment with a reward, its letters fixed, so the reward folds into the costs and the
    wavefront finds it, inside `band`, split past `limit` kept diagonals: its cost and CIGAR, or raises
    when none fits the band."""
    var two = costs.pieces() == 2
    var penalties = extension_penalties(
        match_score,
        costs.mismatch,
        costs.opening,
        costs.extension,
        costs.opening2 if two else 0,
        costs.extension2 if two else 0,
    )
    var found = cigar_within[2](
        reference, query, penalties, extended, Int.MAX, band, ties, limit
    ) if two else cigar_within[1](reference, query, penalties, extended, Int.MAX, band, ties, limit)
    if not found:
        raise outside(band)
    var cigar = found.take().cigar
    var cost = cigar_cost(
        reference, query, cigar, costs.mismatch, costs.opening, costs.extension, costs.opening2, costs.extension2
    )
    return (cost, cigar^)
