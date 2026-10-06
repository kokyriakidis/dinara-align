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

from .alignment import AlignmentResult
from .bit_parallel import ALL_ONES, BAND_COLUMNS, Frontier, Profile, WORD_BITS, word_value
from .edit_distance import edit_alignment
from .errors import AlignmentError


@fieldwise_init
struct EditHit(ImplicitlyCopyable, Writable):
    """Where a pattern best matches a text: its edit distance to `text[start:end]`."""

    var distance: Int
    var start: Int
    var end: Int


def reversed_text(text: String, end: Int) -> String:
    """The first `end` bytes of `text`, back to front."""
    var bytes = text.as_bytes()[0:end]
    var out = List[UInt8](capacity=len(bytes))
    for index in range(len(bytes) - 1, -1, -1):
        out.append(bytes[index])
    return String(unsafe_from_utf8=out)


comptime SEARCH_START = 64
"""The first bound a search's band tries, doubling until the best score along the pattern's last row
falls within it."""


def last_row_scores[free_start: Bool](mut profile: Profile) -> Tuple[Int, Int]:
    """The least score along the pattern's last row and the first column it falls in, the pattern down
    the rows, the text across the columns; with `free_start`, the top row is free, a match starting
    anywhere in the text, else it is the global border.

    As Edlib searches: a bound guessed and doubled, each try sweeping only the band of rows some score
    within it can still reach (see `banded_last_row`), until the least score falls within the bound,
    or the bound covers every row and the band the whole matrix.
    """
    if profile.rows == 0:
        return (0, 0)
    profile.build_planes()
    var bound = SEARCH_START
    while True:
        var found = banded_last_row[free_start](profile, bound)
        if found[0] <= bound or bound >= profile.rows:
            return found
        bound *= 2


def banded_last_row[free_start: Bool](mut profile: Profile, bound: Int) -> Tuple[Int, Int]:
    """`last_row_scores` swept only where a score within `bound` can still lie, Ukkonen's cutoff: its
    least score and first column when that score is within the bound, else some score above it.

    A tile at a time, the band runs down to the last row scoring within the bound at the tile's left
    edge, plus the tile's width, since that row moves at most one down a column. Every word but the
    last goes through the sweep's kernels; the last, once in the band, takes Hyyrö's block step a
    column at a time, its horizontal masks read at the pattern's last row, which seldom ends a word.
    A word entering the band, or entering it again, starts from `+1` all down its left edge, the cost
    of a real path, so every score is at least the true one and those within the bound are exact.
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
    var first_column = 0
    while first_column < columns:
        var end_column = min(first_column + BAND_COLUMNS, columns)
        var end_word = min(ceildiv(min(reach + (end_column - first_column), rows), WORD_BITS), words)
        end_word = max(end_word, 1)
        # Words left behind and now back in the band restart from `+1`, as if never swept.
        for word in range(swept, end_word):
            frontier.vertical_plus[word] = ALL_ONES
            frontier.vertical_minus[word] = 0
        swept = end_word
        if end_word == words:
            # The last row's score at the left edge, read down it before the tile is swept, so it
            # agrees with the words' state however long the row was out of the band: the top border,
            # every word above, and the last word to the pattern's last row.
            score = first_column if not free_start else 0
            for word in range(last):
                score += word_value(frontier.vertical_plus[word], frontier.vertical_minus[word])
            var through = ALL_ONES if bit == UInt64(WORD_BITS - 1) else (UInt64(1) << (bit + 1)) - 1
            score += word_value(frontier.vertical_plus[last] & through, frontier.vertical_minus[last] & through)
        var fast_end = min(end_word, last)
        if fast_end > 0:
            sweep.words(profile.symbols(first_column, end_column), 0, fast_end, first_column, end_column)
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
                if incoming_minus != 0:
                    matches |= 1
                var horizontal = (((matches & vertical_plus) + vertical_plus) ^ vertical_plus) | matches
                var plus = vertical_minus | ~(horizontal | vertical_plus)
                var minus = vertical_plus & horizontal
                score += Int((plus >> bit) & 1) - Int((minus >> bit) & 1)
                if score < best:
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
        var running = end_column if not free_start else 0
        reach = 0
        for word in range(end_word):
            var top = running
            running += word_value(frontier.vertical_plus[word], frontier.vertical_minus[word])
            if (top + running - WORD_BITS) // 2 <= bound:
                reach = min((word + 1) * WORD_BITS, rows)
        first_column = end_column
    return (best, best_column)


def edit_search(pattern: String, text: String, prefix: Bool = False) raises AlignmentError -> EditHit:
    """Where `pattern` best matches inside `text` at unit costs, Edlib's infix mode (HW), or with
    `prefix` where it best matches a prefix of the text, its prefix mode (SHW): the least edit distance
    from the pattern to any `text[start:end]`, `start` zero with `prefix`.

    One sweep over the whole matrix finds the distance and the first end reaching it; the same on both
    reversed, the text cut at that end, finds the latest start reaching it. Symbols as `edit_distance`
    takes them. The whole matrix is swept, `len(text)` columns of `len(pattern) / 64` words, so this
    suits a read against a window of reference rather than a genome.
    """
    var forward = Profile(text, pattern)
    var found: Tuple[Int, Int]
    if prefix:
        found = last_row_scores[False](forward)
        return EditHit(found[0], 0, found[1])
    found = last_row_scores[True](forward)
    var end = found[1]
    var backward = Profile(reversed_text(text, end), reversed_text(pattern, pattern.byte_length()))
    var start_found = last_row_scores[False](backward)
    return EditHit(found[0], end - start_found[1], end)


def edit_search_alignment(
    pattern: String, text: String, prefix: Bool = False
) raises AlignmentError -> Tuple[EditHit, AlignmentResult]:
    """`edit_search`, and an optimal alignment of the text's matched part, `text[start:end]`, first,
    against the whole pattern, second, as `edit_alignment` gives it."""
    var hit = edit_search(pattern, text, prefix)
    var part = String(StringSlice(unsafe_from_utf8=text.as_bytes()[hit.start : hit.end]))
    var aligned = edit_alignment(part, pattern)
    return (hit, aligned^)
