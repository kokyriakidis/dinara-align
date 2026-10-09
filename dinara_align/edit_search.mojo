# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Ported from `pa-bitpacking` in A*PA (https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner,
# commit bf2e14e), by Ragnar Groot Koerkamp and Pesho Ivanov, itself translated from Edlib.
"""
A pattern found inside a text, Edlib's infix mode, or at its start, its prefix mode: the least edit
distance from the pattern to any substring of the text, and an alignment there.
"""

from std.math import ceildiv

from .band import check_margin, checkpoints_passed
from .cigar import reversed_text, text_of
from .bit_parallel import ALL_ONES, BAND_COLUMNS, Frontier, Profile, WORD_BITS, word_value
from .diagonal import PROBE_MARGIN
from .errors import AlignmentError
from .modes import Ties


@fieldwise_init
struct EditHit(ImplicitlyCopyable, Writable):
    """Where a pattern best matches a text: its edit distance to `text[start:end]`."""

    var distance: Int
    var start: Int
    var end: Int


comptime SEARCH_START = 64
"""The first bound a search's band tries, doubling until the best score along the pattern's last row
falls within it."""


def last_row_scores[free_start: Bool](mut profile: Profile, latest: Bool = False) -> Tuple[Int, Int]:
    """The least score along the pattern's last row and the first column it falls in, or with `latest`
    the last, the pattern down the rows, the text across the columns; with `free_start`, the top row is
    free, a match starting anywhere in the text, else it is the global border.

    As Edlib searches: a bound guessed and doubled, each try sweeping only the band of rows some score
    within it can still reach (see `banded_last_row`), until the least score falls within the bound,
    or the bound covers every row and the band the whole matrix.
    """
    if profile.rows == 0:
        # The last row is the top: free, every column scoring nothing, or the border, the first alone.
        return (0, profile.columns if free_start and latest else 0)
    if profile.columns == 0:
        # One column, the whole pattern against nothing; the planes are never built for no text.
        return (profile.rows, 0)
    profile.build_planes()
    var bound = SEARCH_START
    while True:
        var found = banded_last_row[free_start](profile, bound, latest)
        if found[0] <= bound or bound >= profile.rows:
            return (found[0], found[1])
        # A try a checkpoint gave up on aims past the distance its climb projected, as the global
        # band's retry does (see `band.band_doubling`); any other doubles.
        var estimate = found[2]
        var next_bound = 2 * bound
        if estimate >= 0:
            next_bound = max(bound + bound // 4, estimate + estimate // 8 + PROBE_MARGIN)
        bound = next_bound


def banded_last_row[free_start: Bool](mut profile: Profile, bound: Int, latest: Bool = False) -> Tuple[Int, Int, Int]:
    """`last_row_scores` swept only where a score within `bound` can still lie, Ukkonen's cutoff: its
    least score and first column, or with `latest` its last, when that score is within the bound, else
    some score above it; and -1, or the distance a checkpoint projected when it gave the try up.

    Without `free_start`, at an eighth, a quarter and half of the pattern's length across the text, the
    least score down the band has climbed about in proportion to the columns crossed, so scaled to the
    pattern's whole length it projects the distance, as the global band's checkpoints do (see
    `band.Band.check`); a try whose projection less its margin passes the bound gives up there, where
    without it the try would sweep on until its band emptied, two thirds of the way across a 1 kbp
    read at 10% under a bound of 64. Giving up only ever loses a try that would have failed.

    A tile at a time, the band runs down to the last row scoring within the bound at the tile's left
    edge, plus the tile's width, since that row moves at most one down a column. Every word but the
    last goes through the sweep's kernels; the last, once in the band, takes Hyyrö's block step a
    column at a time, its horizontal masks read at the pattern's last row, which seldom ends a word.
    A word entering the band, or entering it again, starts from `+1` all down its left edge, the cost
    of a real path, so every score is at least the true one and those within the bound are exact.
    Without `free_start` the sweep also stops at the first column no end past it could match the best
    end found so far in, and once no row of a tile's right edge is within the bound.
    """
    var columns = profile.columns
    var rows = profile.rows
    var words = profile.words
    var last = words - 1
    var frontier = Frontier(columns, words)
    comptime if free_start:
        for column in range(columns):
            frontier.horizontal_plus[column] = 0
    var sweep = frontier.sweep(profile)
    var bit = UInt64((rows - 1) % WORD_BITS)
    var row_low = sweep.row_low[unsafe_offset=last]
    var row_high = sweep.row_high[unsafe_offset=last]
    var row_extra = sweep.row_extra[unsafe_offset=last] if profile.extended else UInt64(0)
    # The pattern's last row scores `rows` at the first column, before any of the text.
    var score = rows
    var best = score
    var best_column = 0
    # The last row reachable within the bound at the current left edge, and the words swept so far.
    var reach = min(bound, rows)
    var swept = 0
    # Without a free start the band's top moves down too: the first word it keeps, and the score at that
    # word's top on the current left edge, which the words left above it carry.
    var top = 0
    var anchor = 0
    var checkpoint = 0
    var first_column = 0
    while first_column < columns:
        comptime if not free_start:
            # From the global border, the last row at column `c` costs at least `c - rows`: the length
            # the text has run past the pattern. So no column past `rows + best` matches the best so
            # far, a later equal one included, and the sweep ends there.
            if first_column >= rows + best:
                break
        var end_column = min(first_column + BAND_COLUMNS, columns)
        var end_word = min(ceildiv(min(reach + (end_column - first_column), rows), WORD_BITS), words)
        end_word = max(end_word, 1)
        # Words left behind and now back in the band restart from `+1`, as if never swept.
        for word in range(swept, end_word):
            frontier.vertical_plus[word] = ALL_ONES
            frontier.vertical_minus[word] = 0
        swept = end_word
        comptime if not free_start:
            # A row more than `bound` above a column's diagonal scores more than `bound` there and on
            # every column after, so the words wholly above row `first_column - bound` leave the band,
            # their scores carried into the anchor. The top word then reads `+1` from above, the cost of
            # a real path, as the global band's does: every score stays at least the true one, and one
            # within the bound, whose optimal path stays inside the band, exact.
            var new_top = min(max(top, (first_column - bound) // WORD_BITS), end_word - 1, last)
            anchor += frontier.climb(top, new_top, rows)
            top = new_top
        if end_word == words:
            # The last row's score at the left edge, read down it before the tile is swept, so it
            # agrees with the words' state however long the row was out of the band: the top border,
            # or the anchor at the band's top, every word between, and the last word to the pattern's
            # last row.
            score = (anchor if top > 0 else first_column) if not free_start else 0
            score += frontier.climb(top, words, rows)
        var fast_end = min(end_word, last)
        if fast_end > top:
            sweep.words(profile.symbols(first_column, end_column), top, fast_end, first_column, end_column)
        if end_word == words:
            var vertical_plus = frontier.vertical_plus[last]
            var vertical_minus = frontier.vertical_minus[last]
            for column in range(first_column, end_column):
                var matches = (sweep.column_low[unsafe_offset=column] ^ row_low) & (
                    sweep.column_high[unsafe_offset=column] ^ row_high
                )
                if profile.extended:
                    matches &= sweep.column_extra[unsafe_offset=column] ^ row_extra
                var incoming_plus = frontier.horizontal_plus[column]
                var incoming_minus = frontier.horizontal_minus[column]
                var crossing = matches | vertical_minus
                # Myers' step takes no `-1` entering from above; Hyyrö folds one in as a match at the top.
                if incoming_minus != 0:
                    matches |= 1
                var horizontal = (((matches & vertical_plus) + vertical_plus) ^ vertical_plus) | matches
                var plus = vertical_minus | ~(horizontal | vertical_plus)
                var minus = vertical_plus & horizontal
                score += Int((plus >> bit) & 1) - Int((minus >> bit) & 1)
                if score < best or (latest and score == best):
                    best = score
                    best_column = column + 1
                plus = (plus << 1) | incoming_plus
                minus = (minus << 1) | incoming_minus
                vertical_plus = minus | ~(crossing | plus)
                vertical_minus = plus & crossing
            frontier.vertical_plus[last] = vertical_plus
            frontier.vertical_minus[last] = vertical_minus
        # Down the right edge, each word's top and bottom scores bound its least, the rows a step
        # apart: at least `(top + bottom - 64) / 2`. The band reaches the last word that may hold
        # a score within the bound.
        # The anchor before any word leaves is the top border's score at the left edge, the column.
        if top == 0:
            anchor = first_column
        anchor += end_column - first_column
        var running = anchor if not free_start else 0
        var least = running
        reach = 0
        for word in range(top, end_word):
            var top = running
            running += word_value(frontier.vertical_plus[word], frontier.vertical_minus[word])
            least = min(least, max((top + running - WORD_BITS) // 2, 0))
            if (top + running - WORD_BITS) // 2 <= bound:
                reach = min((word + 1) * WORD_BITS, rows)
        comptime if not free_start:
            # Nothing on the right edge within the bound: every later cell is reached across it, so none
            # is within the bound either, and the try has failed. A free start begins anew anywhere.
            if reach == 0:
                break
            var passed = checkpoints_passed(checkpoint, end_column, rows)
            if passed != checkpoint and end_column < rows:
                checkpoint = passed
                var estimate = least * rows // end_column
                if estimate - check_margin(estimate, passed) // 2 > bound:
                    return (bound + 1, 0, estimate)
        first_column = end_column
    return (best, best_column, -1)


def edit_search(
    pattern: String, text: String, prefix: Bool = False, ties: Ties = Ties.LEFT
) raises AlignmentError -> EditHit:
    """Where `pattern` best matches inside `text` at unit costs, Edlib's infix mode (HW), or with
    `prefix` where it best matches a prefix of the text, its prefix mode (SHW): the least edit distance
    from the pattern to any `text[start:end]`, `start` zero with `prefix`.

    Of equally good matches, the span the wavefront's rule picks (see `gap_affine.free_ends_alignment`):
    with `Ties.LEFT` the last end, then the last start for it; with `Ties.RIGHT` the first start, then
    the first end for it. One sweep across the text finds the distance and the end, and one over both
    reversed, the text cut at that end, the start; for `Ties.RIGHT` the other way round. Symbols as
    `edit_distance` takes them. Every column of the text is swept, down as many of the pattern's
    `len(pattern) / 64` words as the bound reaches, so this suits a read against a window of reference
    rather than a genome.
    """
    var length = text.byte_length()
    if ties == Ties.RIGHT and not prefix:
        # The first start: the last end of both reversed, then the first end from that start.
        var backward = Profile(reversed_text(text.as_bytes()[0:length]), reversed_text(pattern.as_bytes()))
        var found = last_row_scores[True](backward, True)
        var start = length - found[1]
        var tail = text_of(text.as_bytes()[start:])
        var forward = Profile(tail, pattern)
        var end_found = last_row_scores[False](forward)
        return EditHit(found[0], start, start + end_found[1])
    var forward = Profile(text, pattern)
    var latest = ties == Ties.LEFT
    if prefix:
        var found = last_row_scores[False](forward, latest)
        return EditHit(found[0], 0, found[1])
    var found = last_row_scores[True](forward, True)
    var end = found[1]
    var backward = Profile(reversed_text(text.as_bytes()[0:end]), reversed_text(pattern.as_bytes()))
    var start_found = last_row_scores[False](backward)
    return EditHit(found[0], end - start_found[1], end)
