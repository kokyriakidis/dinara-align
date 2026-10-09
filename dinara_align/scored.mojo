# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
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


from .common import SubstitutionDType
from .cigar import reversed_list, cigar_counts, cigar_matches, cigar_runs, reversed_cigar, reversed_text, text_of
from .errors import AlignmentError
from .gap_affine import kept_bytes, rings_within
from .gap_affine import (
    AffineExtension,
    EndsFree,
    cigar_within,
    extension_of,
    rewarded_penalties,
    outside,
    traced_extension,
)
from .modes import Alignment, Anchor, Band, Costs, Ties
from .anti_diagonals import AntiDiagonals, GapLanes, lane_bits
from .substitutions import SubstitutionLookup


comptime ANYWHERE = 0
"""A sweep for a local alignment: it starts at any cell, every cell's best floored at zero, and ends at
any cell."""
comptime FROM_EDGE = 1
"""A sweep for ends-free alignment with a reward: it starts on the first row or column, for nothing
within the letters free there and past them paying a gap, and ends on the last row or column within
the letters free there. With none free at the start it starts at the origin."""
comptime FROM_ORIGIN = 2
"""A sweep for an extension: it starts at the origin and ends at any cell, aligning nothing scoring
zero."""


def best_end[
    pieces: Int, dtype: DType, width: Int, transposed: Bool, kind: Int = ANYWHERE
](
    reference: Span[UInt8, _],
    query: Span[UInt8, _],
    costs: Costs,
    match_score: Int,
    ends: EndsFree = EndsFree(),
    highest: Bool = True,
) -> Tuple[Int, Int, Int, Bool]:
    """The best score of an alignment starting and ending where `kind` allows, where it ends, and whether
    a Z-drop gave the sweep up (see `swept_cells`)."""
    var unused = List[Int]()
    var lookup = SubstitutionLookup.uniform_of(match_score, -costs.mismatch)
    return swept_cells[pieces, dtype, width, transposed, kind, False, True](
        reference, query, costs, ends, highest, lookup, -1, unused
    )


def swept_cells[
    pieces: Int, dtype: DType, width: Int, transposed: Bool, kind: Int, columns_kept: Bool, compared: Bool
](
    reference: Span[UInt8, _],
    query: Span[UInt8, _],
    costs: Costs,
    ends: EndsFree,
    highest: Bool,
    lookup: SubstitutionLookup,
    zdrop: Int,
    mut kept: List[Int],
) -> Tuple[Int, Int, Int, Bool]:
    """The best score of an alignment starting and ending where `kind` allows, and the reference's and
    the query's letters up to where it ends: of several such ends, for a local alignment the furthest
    along both together, then along the reference; from the edges the one on the highest diagonal, the
    furthest along the reference less the query, or with `highest` false the lowest; an extension's
    as a local alignment's. With `columns_kept`, for a local alignment, `kept` holds after it the best
    score of every column: of the cells with as many reference letters, `len(reference) + 1` of them.
    `lookup` scores each pair, the reference's letter its row: a `Scoring`'s substitutions, the sequences
    codes into them, or with `compared` a match's reward and `costs.mismatch`. An extension with a `zdrop` of zero or more
    stops as KSW2's does: once a diagonal's best lies more than `zdrop`, plus an extension a diagonal
    between them, below the best so far, the best so far stands.

    Gotoh's recurrence over scores, a local alignment's every cell floored at zero, swept by
    anti-diagonal (see `anti_diagonals`) in `width` lanes of `dtype`, as narrow as the scores allow (see
    `local_alignment`). For a local alignment each lane keeps its best of the diagonal and the step it came
    at, and the diagonal's best is placed only when it meets the best so far; one from the edges, `ends`
    saying which letters at each end are free, ends on the last row or column, at most two cells a
    diagonal, read once it is done."""
    comptime Lanes = SIMD[dtype, width]
    comptime Value = Scalar[dtype]
    comptime LOW = Value.MIN // 4
    # The lanes run along the shorter sequence, so the sweep's rows stay in the core's own cache:
    # along the reference, or with `transposed` along the query.
    var down_letters = query if transposed else reference
    var across_letters = reference if transposed else query
    var rows = len(down_letters)
    var columns = len(across_letters)
    comptime local = kind == ANYWHERE
    # A local alignment and an extension end at any cell, the best so far kept as the sweep goes.
    comptime anywhere_end = kind != FROM_EDGE
    # The free letters at either end of the lanes' sequence and the other's; a local alignment's start
    # is free everywhere.
    var down_start = Int.MAX if local else (ends.second_begin if transposed else ends.first_begin)
    var down_end = ends.second_end if transposed else ends.first_end
    var across_start = Int.MAX if local else (ends.first_begin if transposed else ends.second_begin)
    var across_end = ends.first_end if transposed else ends.second_end
    var sweep = AntiDiagonals[dtype, width, pieces](down_letters, across_letters, LOW)
    var cells = sweep.cells()
    # A gap's first letter pays the opening and its extension, each further letter the extension: a
    # gap down the lanes is of the lanes' letters, the reference's a deletion's, across them the other's.
    var down_opening = costs.opening if transposed else costs.deletion_opening
    var down_extension = costs.extension if transposed else costs.deletion_extension
    var down_opening2 = costs.opening2 if transposed else costs.deletion_opening2
    var down_extension2 = costs.extension2 if transposed else costs.deletion_extension2
    var across_opening = costs.deletion_opening if transposed else costs.opening
    var across_extension = costs.deletion_extension if transposed else costs.extension
    var across_opening2 = costs.deletion_opening2 if transposed else costs.opening2
    var across_extension2 = costs.deletion_extension2 if transposed else costs.extension2
    comptime two = pieces == 2
    var gaps = GapLanes[dtype, width](
        Lanes(Value(-down_opening - down_extension)),
        Lanes(Value(-down_extension)),
        Lanes(Value(-across_opening - across_extension)),
        Lanes(Value(-across_extension)),
        Lanes(Value(-down_opening2 - down_extension2) if two else 0),
        Lanes(Value(-down_extension2) if two else 0),
        Lanes(Value(-across_opening2 - across_extension2) if two else 0),
        Lanes(Value(-across_extension2) if two else 0),
    )
    var substitute = lookup.lanes[width, dtype, compared]()
    var zero = Lanes(0)
    # A Z-drop's slack a diagonal between two bests: the cheapest a gap grows a letter.
    var slack = costs.cheapest_extension()
    # What a lane past the diagonal's last row counts as: nothing a best could be.
    var nothing = zero if local else Lanes(LOW)
    var lane_index = Lanes()
    comptime for lane in range(width):
        lane_index[lane] = Value(lane)
    # Each column's best as the lanes find it: by the lanes' own row along the reference, or along the
    # query back from the reference's end, so a lane group's cells are always a run of the store.
    var column_best = List[Value](length=(len(reference) + 1 + width) if columns_kept else 0, fill=0)

    @inline(.always)
    def edge(letters: Int, free: Int, down: Bool) {imm costs} -> Value:
        """A cell on the first column, `letters` of the lanes' sequence in with `down`, or on the first row,
        of the other's: nothing within the `free` letters, past them a gap from where they end, at the
        cheaper piece, a deletion's costs for the reference's letters."""
        if letters <= free:
            return 0
        return Value(-costs.gap(letters - free, down != transposed))

    if rows == 0 or columns == 0:
        comptime if anywhere_end:
            comptime if columns_kept:
                kept = List[Int](length=len(reference) + 1, fill=0)
            return (0, 0, 0, False)
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
            var value = Int(edge(letters, start, rows > 0))
            if value > top or (value == top and (along_reference == highest)):
                top = value
                reach = letters
        var down = reach if rows > 0 else 0
        var across = reach if columns > 0 else 0
        if transposed:
            return (top, across, down, False)
        return (top, down, across, False)

    var best = 0 if anywhere_end else Int.MIN
    var best_row = 0
    var best_diagonal = 0
    # Whether a Z-drop gave the sweep up before it covered the matrix.
    var dropped = False

    @inline(.always)
    def ending(
        diagonal: Int, last_row: Value, last_column: Value, mut best: Int, mut best_row: Int, mut best_diagonal: Int
    ) {imm rows, imm columns, imm across_end, imm down_end, imm highest}:
        """The ends on anti-diagonal `diagonal`, its cells on the last row and the last column given: those
        within the letters free there, the one on the highest diagonal, the reference's letters less the
        query's, or with `highest` false the lowest, on a tie."""
        for side in range(2):
            var row = rows if side == 0 else diagonal - columns
            if row < 0 or row > rows or diagonal - row < 0 or diagonal - row > columns:
                continue
            if side == 0 and columns - (diagonal - row) > across_end:
                continue
            if side == 1 and rows - row > down_end:
                continue
            var value = Int(last_row if side == 0 else last_column)
            # The reference's letters less the query's, from the lanes' row and the rest of the anti-diagonal.
            var lean = (diagonal - 2 * row) if transposed else (2 * row - diagonal)
            var best_lean = (best_diagonal - 2 * best_row) if transposed else (2 * best_row - best_diagonal)
            var beyond = (lean > best_lean) if highest else (lean < best_lean)
            if value > best or (value == best and beyond):
                best = value
                best_row = row
                best_diagonal = diagonal

    # Diagonal one, beside the origin: a letter of either sequence against nothing.
    cells.one_back[unsafe_offset=0] = edge(1, across_start, False)
    cells.one_back[unsafe_offset=1] = edge(1, down_start, True)
    comptime if not anywhere_end:
        ending(1, cells.filled(rows), cells.filled(max(1 - columns, 0)), best, best_row, best_diagonal)
    for diagonal in range(2, rows + columns + 1):
        var low = max(1, diagonal - columns)
        var high = min(rows, diagonal - 1)
        var lag = columns - diagonal
        var row = low
        # Each lane's best on this diagonal, and the step it came at, the later on a tie.
        var top = nothing
        var top_step = zero
        var step = Value(0)
        while row <= high:
            var score = cells.step[local, transposed](row, lag, substitute, gaps)
            comptime if anywhere_end:
                # Lanes past the diagonal's last row hold no cell.
                var counted = score
                if row + width - 1 > high:
                    counted = lane_index.lt(Lanes(Value(high - row + 1))).select(score, nothing)
                comptime if columns_kept:
                    var at = column_best.unsafe_ptr().unsafe_offset((columns - diagonal + row) if transposed else row)
                    at.unsafe_store(max(at.unsafe_load[width=width](), counted))
                var later = counted.gt(top) if transposed else counted.ge(top)
                top = later.select(counted, top)
                top_step = later.select(Lanes(step), top_step)
                step += 1
            row += width
        comptime if anywhere_end:
            # The diagonal's best, placed at its furthest row, when it meets the best so far.
            var most = Int(top.reduce_max())
            var placed = most > 0 and most >= best
            if placed or (zdrop >= 0 and most < best):
                var furthest = 0 if not transposed else Int.MAX
                comptime for lane in range(width):
                    if Int(top[lane]) == most:
                        var at = low + Int(top_step[lane]) * width + lane
                        furthest = min(furthest, at) if transposed else max(furthest, at)
                if placed:
                    best = most
                    best_row = furthest
                    best_diagonal = diagonal
                else:
                    # A Z-drop: how far the diagonal's best has fallen, a gap's slack between them.
                    var lean = (diagonal - 2 * furthest) if transposed else (2 * furthest - diagonal)
                    var best_lean = (best_diagonal - 2 * best_row) if transposed else (2 * best_row - best_diagonal)
                    if best - most > zdrop + slack * abs(lean - best_lean):
                        dropped = True
                        break
        # The border cells of this diagonal, written after the lanes that may have run over them: free
        # within the letters free there, past them a gap. No interior cell reads a border's gap layers.
        if diagonal <= columns:
            cells.current[unsafe_offset=0] = edge(diagonal, across_start, False)
            cells.deletes[unsafe_offset=0] = LOW
            comptime if two:
                cells.deletes2[unsafe_offset=0] = LOW
        if diagonal <= rows:
            cells.current[unsafe_offset=diagonal] = edge(diagonal, down_start, True)
            cells.inserts[unsafe_offset=diagonal] = LOW
            comptime if two:
                cells.inserts2[unsafe_offset=diagonal] = LOW
        cells.advance()
        comptime if not anywhere_end:
            ending(
                diagonal,
                cells.filled(rows),
                cells.filled(max(diagonal - columns, 0)),
                best,
                best_row,
                best_diagonal,
            )
    comptime if columns_kept:
        var length = len(reference)
        kept = List[Int](length=length + 1, fill=0)
        for letters in range(length + 1):
            kept[letters] = Int(column_best[(columns - letters) if transposed else letters])
    # `best_row` counts the lanes' sequence, the rest of the diagonal the other's.
    if transposed:
        return (best, best_diagonal - best_row, best_row, dropped)
    return (best, best_row, best_diagonal - best_row, dropped)


def sweep_bits[kind: Int](costs: Costs, match_score: Int, rows: Int, columns: Int) -> Int:
    """The narrowest lanes every score of the sweep under `costs` fits (see `anti_diagonals.lane_bits`)."""
    return lane_bits[kind == ANYWHERE](match_score, costs.dearest_step(), rows, columns)


def swept[
    kind: Int
](
    reference: Span[UInt8, _],
    query: Span[UInt8, _],
    costs: Costs,
    match_score: Int,
    ends: EndsFree = EndsFree(),
    highest: Bool = True,
) -> Tuple[Int, Int, Int, Bool]:
    """`best_end` with its lanes along the shorter sequence, as narrow as the scores allow (see `sweep`)."""
    var unused = List[Int]()
    var no_table = List[Scalar[SubstitutionDType]]()
    return sweep[kind](
        reference,
        query,
        costs,
        match_score,
        ends,
        highest,
        sweep_bits[kind](costs, match_score, len(reference), len(query)),
        no_table,
        0,
        -1,
        unused,
    )


def sweep[
    kind: Int, columns_kept: Bool = False
](
    reference: Span[UInt8, _],
    query: Span[UInt8, _],
    costs: Costs,
    match_score: Int,
    ends: EndsFree,
    highest: Bool,
    bits: Int,
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet: Int,
    zdrop: Int,
    mut kept: List[Int],
) -> Tuple[Int, Int, Int, Bool]:
    """`swept_cells` with its lanes along the shorter sequence, `bits` to a lane, 16, 32 or 64 (see
    `anti_diagonals.lane_bits`), a register's 64 bytes of them, each pair's score read from `substitutions`
    over an alphabet of `alphabet` codes where it is one, else a match earning `match_score`: the one choice
    of type, pieces, orientation and table every sweep's caller makes."""
    if bits == 16:
        return sweep_in[DType.int16, 32, kind, columns_kept](
            reference, query, costs, match_score, ends, highest, substitutions, alphabet, zdrop, kept
        )
    if bits == 32:
        return sweep_in[DType.int32, 16, kind, columns_kept](
            reference, query, costs, match_score, ends, highest, substitutions, alphabet, zdrop, kept
        )
    return sweep_in[DType.int64, 8, kind, columns_kept](
        reference, query, costs, match_score, ends, highest, substitutions, alphabet, zdrop, kept
    )


def sweep_in[
    dtype: DType, width: Int, kind: Int, columns_kept: Bool
](
    reference: Span[UInt8, _],
    query: Span[UInt8, _],
    costs: Costs,
    match_score: Int,
    ends: EndsFree,
    highest: Bool,
    substitutions: ImmSpan[Scalar[SubstitutionDType], _],
    alphabet: Int,
    zdrop: Int,
    mut kept: List[Int],
) -> Tuple[Int, Int, Int, Bool]:
    """`sweep` in lanes of `dtype`, `width` of them."""
    var transposed = len(query) < len(reference)
    if alphabet > 0:
        var table = SubstitutionLookup(substitutions, alphabet)
        if transposed:
            return swept_cells[1, dtype, width, True, kind, columns_kept, False](
                reference, query, costs, ends, highest, table, zdrop, kept
            )
        return swept_cells[1, dtype, width, False, kind, columns_kept, False](
            reference, query, costs, ends, highest, table, zdrop, kept
        )
    var lookup = SubstitutionLookup.uniform_of(match_score, -costs.mismatch)
    if costs.pieces() == 2:
        if transposed:
            return swept_cells[2, dtype, width, True, kind, columns_kept, True](
                reference, query, costs, ends, highest, lookup, zdrop, kept
            )
        return swept_cells[2, dtype, width, False, kind, columns_kept, True](
            reference, query, costs, ends, highest, lookup, zdrop, kept
        )
    if transposed:
        return swept_cells[1, dtype, width, True, kind, columns_kept, True](
            reference, query, costs, ends, highest, lookup, zdrop, kept
        )
    return swept_cells[1, dtype, width, False, kind, columns_kept, True](
        reference, query, costs, ends, highest, lookup, zdrop, kept
    )


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
    var kept = List[Int]()
    var no_table = List[Scalar[SubstitutionDType]]()
    var found = sweep[ANYWHERE, True](
        reference,
        query,
        costs,
        match_score,
        EndsFree(),
        True,
        sweep_bits[ANYWHERE](costs, match_score, len(reference), len(query)),
        no_table,
        0,
        -1,
        kept,
    )
    var second = 0
    var second_end = 0
    for letters in range(len(kept)):
        if abs(letters - found[1]) <= window:
            continue
        if kept[letters] > second:
            second = kept[letters]
            second_end = letters
    return LocalScores(found[0], found[1], found[2], second, second_end)


def end_of(
    reference: Span[UInt8, _], query: Span[UInt8, _], costs: Costs, match_score: Int
) -> Tuple[Int, Int, Int, Bool]:
    """A local alignment's best score and where it ends (see `best_end`)."""
    return swept[ANYWHERE](reference, query, costs, match_score)


def latest_local(
    reference: String, query: String, costs: Costs, match_score: Int, eqx: Bool, limit: Int
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
    var penalties = rewarded_penalties(match_score, costs)
    rings_within(penalties, end_column, end_row, kept_bytes(limit))
    var head = text_of(reference.as_bytes()[:end_column])
    var lead = text_of(query.as_bytes()[:end_row])
    # The best alignment ending at that cell, and starting wherever pays: an extension back from it,
    # which stops on earning the sweep's score, the best any alignment ending there earns.
    var traced = traced_extension[2](head, lead, penalties, eqx, found[0], limit) if two else traced_extension[1](
        head, lead, penalties, eqx, found[0], limit
    )
    var back: AffineExtension
    if traced:
        back = traced.take()
        back.cigar = reversed_cigar(back.cigar)
    else:
        # Too many fronts to keep: the extension's own search and split, by the left rule too, over
        # both sequences reversed from the end.
        var mirrored = extension_of[2](
            reversed_text(head.as_bytes()),
            reversed_text(lead.as_bytes()),
            penalties,
            eqx,
            Anchor.START,
            Band(),
            Ties.RIGHT,
            found[0],
            limit,
        ) if two else extension_of[1](
            reversed_text(head.as_bytes()),
            reversed_text(lead.as_bytes()),
            penalties,
            eqx,
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
    reference: String, query: String, costs: Costs, match_score: Int, ties: Ties, eqx: Bool, limit: Int
) raises AlignmentError -> Alignment:
    """The best local alignment, a match earning `match_score` (see `Mode.local`).

    `Ties.LEFT` ends it as late as an equally good alignment allows and starts it as late too, the
    shortest, its CIGAR by the left rule (see `latest_local`); `Ties.RIGHT` is that over both sequences
    reversed, read backwards, as for every other mode: it starts and ends as early as it may, gaps
    right."""
    if ties == Ties.LEFT:
        return latest_local(reference, query, costs, match_score, eqx, limit)
    var columns = reference.byte_length()
    var rows = query.byte_length()
    var mirrored = latest_local(
        reversed_text(reference.as_bytes()), reversed_text(query.as_bytes()), costs, match_score, eqx, limit
    )
    if mirrored.score == 0:
        return mirrored^
    return mirrored.mirrored(columns, rows)


def started_span[
    B: def(List[UInt8], List[UInt8], EndsFree) -> Tuple[Int, Int, Int, Bool]
](
    first: ImmSpan[UInt8, _], second: ImmSpan[UInt8, _], ends: EndsFree, forward: Tuple[Int, Int, Int, Bool], back: B
) -> Tuple[Int, Int, Int, Int, Int]:
    """The span the rule of `Ties.LEFT` names (see `gap_affine.free_ends_alignment`) from a sweep from the edges,
    `forward`, which found the best score and the end on the highest diagonal an optimum reaches: the start on
    the highest of those ending there, by `back`, a sweep over both sequences reversed from that end to the
    letters free at the start, whose lowest diagonal it takes. The score, then the start's and the end's
    columns and rows."""
    var end_column = forward[1]
    var end_row = forward[2]
    var found = back(
        reversed_list(first[:end_column]),
        reversed_list(second[:end_row]),
        EndsFree(0, ends.first_begin, 0, ends.second_begin),
    )
    return (forward[0], end_column - found[1], end_row - found[2], end_column, end_row)


def rewarded_span(
    reference: String, query: String, costs: Costs, match_score: Int, ends: EndsFree
) -> Tuple[Int, Int, Int, Int, Int]:
    """The best score with `ends`' letters free and a match earning `match_score`, and the span the rule
    of `Ties.LEFT` names (see `gap_affine.free_ends_alignment`): its end on the highest diagonal an
    optimum reaches, by a sweep from the edges, and its start on the highest of those ending there, by
    a sweep back from that end over both sequences reversed, to the letters free at the start. The
    score, then the start's and the end's columns and rows."""
    var forward = swept[FROM_EDGE](reference.as_bytes(), query.as_bytes(), costs, match_score, ends)

    def back(
        head: List[UInt8], lead: List[UInt8], starts: EndsFree
    ) {imm costs, imm match_score} -> Tuple[Int, Int, Int, Bool]:
        """The sweep back, the lowest diagonal of the reversed sequences."""
        return swept[FROM_EDGE](Span(head), Span(lead), costs, match_score, starts, False)

    return started_span(reference.as_bytes(), query.as_bytes(), ends, forward, back)


def rewarded_alignment(
    reference: String,
    query: String,
    costs: Costs,
    match_score: Int,
    ends: EndsFree,
    ties: Ties,
    eqx: Bool,
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
            reversed_text(reference.as_bytes()),
            reversed_text(query.as_bytes()),
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
    var part = text_of(reference.as_bytes()[start_column:end_column])
    var piece = text_of(query.as_bytes()[start_row:end_row])
    var found = global_rewarded(part, piece, costs, match_score, Band(), ties, eqx, limit)
    return Alignment(found[0], span[0], found[1], start_column, end_column, start_row, end_row)


def global_rewarded(
    reference: String,
    query: String,
    costs: Costs,
    match_score: Int,
    band: Band,
    ties: Ties,
    eqx: Bool,
    limit: Int,
) raises AlignmentError -> Tuple[Int, String]:
    """A global alignment with a reward, its letters fixed, so the reward folds into the costs and the
    wavefront finds it, inside `band`, split past `limit` kept diagonals: its cost and CIGAR, or raises
    when none fits the band."""
    var two = costs.pieces() == 2
    var penalties = rewarded_penalties(match_score, costs)
    rings_within(penalties, reference.byte_length(), query.byte_length(), kept_bytes(limit))
    var found = cigar_within[2](reference, query, penalties, eqx, Int.MAX, band, ties, limit) if two else cigar_within[
        1
    ](reference, query, penalties, eqx, Int.MAX, band, ties, limit)
    if not found:
        raise outside(band)
    var cigar = found.take().cigar
    var cost = costs_of_cigar(reference, query, cigar, costs)
    return (cost, cigar^)


def costs_of_cigar(reference: String, query: String, cigar: String, costs: Costs) -> Int:
    """What a CIGAR of `reference` against `query` costs: its substitutions, `M` runs compared letter
    by letter, and each gap run at its direction's cheaper piece."""
    var runs = cigar_runs(cigar)
    var gaps = 0
    for index in range(len(runs[0])):
        var letter = runs[0][index]
        if letter == UInt8(ord("D")) or letter == UInt8(ord("I")):
            gaps += costs.gap(runs[1][index], letter == UInt8(ord("D")))
    return cigar_counts(reference.as_bytes(), query.as_bytes(), cigar).mismatches * costs.mismatch + gaps
