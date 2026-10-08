# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
One query against many references, a database search: every reference's best score or least cost,
the hits ranked, the best kept, and alignments for those alone (see `search`).

A local alignment scores many references at once, one to a SIMD lane, as SWIPE (Rognes, 2011) and
hyalite's database mode do: the query runs down the rows of every lane's matrix together, each
reference across its own lane's columns, so a block of references of about one length costs about
what one of them costs alone, every lane busy on every cell. Each pair's own sweep keeps its lanes
along one pair's anti-diagonals, which a short pair fills only in part; across references there is
no such waste. Other modes take each pair's own search, the pairs spread over every thread.
"""

from std.atomic import Atomic

from max.algorithm import parallelize

from .common import hardware_threads
from .errors import AlignmentError, ErrorKind
from .modes import Alignment, Costs, Mode, Ties

comptime BLOCK_LOW = Int16.MIN // 4
"""A cell no path reaches, far enough below zero that a few more costs stay negative."""


@fieldwise_init
struct Hit(Copyable, Movable, Writable):
    """One reference's result in a search: its place in the list searched, its best score (minus its
    least cost for a mode with no reward), and its alignment when asked for."""

    var index: Int
    var score: Int
    var alignment: Optional[Alignment]


def lanes_fit(costs: Costs, match_score: Int, query_length: Int, longest: Int) -> Bool:
    """Whether every score of a block fits 16 bits: none passes the reward of the query matched
    throughout, and none falls further below zero than the dearest single move."""
    var dearest = max(
        costs.mismatch, max(costs.opening + costs.extension, costs.deletion_opening + costs.deletion_extension)
    )
    if costs.pieces() == 2:
        dearest = max(
            dearest, max(costs.opening2 + costs.extension2, costs.deletion_opening2 + costs.deletion_extension2)
        )
    return match_score * (min(query_length, longest) + 1) < 32000 and dearest < 4000


def block_scores[
    dtype: DType, width: Int, pieces: Int
](query: Span[UInt8, _], references: List[String], members: List[Int], costs: Costs, match_score: Int) -> List[Int]:
    """The best local score of the query against each of up to `width` references, `members`, one to a
    lane: Gotoh's recurrence floored at zero, row by row of the query, a lane's columns its reference's
    letters. Cells past a reference's end compare its padding, which matches nothing, so they score no
    more than the cells before them, and the best of a lane is its reference's."""
    comptime Lanes = SIMD[dtype, width]
    comptime Value = Scalar[dtype]
    comptime LOW = Value.MIN // 4
    comptime two = pieces == 2
    var longest = 0
    for member in members:
        longest = max(longest, references[member].byte_length())
    # The references side by side, a column's letters of every lane together; padding past each end.
    var letters = List[UInt8](length=longest * width, fill=0xFE)
    for lane in range(len(members)):
        var bytes = references[members[lane]].as_bytes()
        for column in range(len(bytes)):
            letters[column * width + lane] = bytes[column]
    var rows = len(query)
    # The previous row's cells and its gaps down the columns, each column its lanes; the gap across a
    # row runs as the row goes.
    var above = List[Value](length=(longest + 1) * width, fill=0)
    var gaps_down = List[Value](length=(longest + 1) * width, fill=LOW)
    var gaps_down2 = List[Value](length=(longest + 1) * width if two else 0, fill=LOW)
    # Down a column is a query letter alone, an insertion; across a row a reference letter alone.
    var first_down = Lanes(Value(costs.opening + costs.extension))
    var further_down = Lanes(Value(costs.extension))
    var first_across = Lanes(Value(costs.deletion_opening + costs.deletion_extension))
    var further_across = Lanes(Value(costs.deletion_extension))
    var first_down2 = Lanes(Value(costs.opening2 + costs.extension2) if two else 0)
    var further_down2 = Lanes(Value(costs.extension2) if two else 0)
    var first_across2 = Lanes(Value(costs.deletion_opening2 + costs.deletion_extension2) if two else 0)
    var further_across2 = Lanes(Value(costs.deletion_extension2) if two else 0)
    var matched = Lanes(Value(match_score))
    var mismatched = Lanes(-Value(costs.mismatch))
    var zero = Lanes(0)
    var best = zero
    for row in range(rows):
        var letter = Lanes(Value(Int(query[row])))
        var diagonal = zero
        var left = zero
        var across = Lanes(LOW)
        var across2 = Lanes(LOW)
        for column in range(1, longest + 1):
            var at = column * width
            var up = above.unsafe_ptr().unsafe_offset(at).unsafe_load[width=width]()
            var down = max(
                up - first_down, gaps_down.unsafe_ptr().unsafe_offset(at).unsafe_load[width=width]() - further_down
            )
            gaps_down.unsafe_ptr().unsafe_offset(at).unsafe_store(down)
            across = max(left - first_across, across - further_across)
            var theirs = letters.unsafe_ptr().unsafe_offset((column - 1) * width).unsafe_load[width=width]()
            var score = max(
                diagonal + theirs.cast[dtype]().eq(letter).select(matched, mismatched), max(max(down, across), zero)
            )
            comptime if two:
                var down2 = max(
                    up - first_down2,
                    gaps_down2.unsafe_ptr().unsafe_offset(at).unsafe_load[width=width]() - further_down2,
                )
                gaps_down2.unsafe_ptr().unsafe_offset(at).unsafe_store(down2)
                across2 = max(left - first_across2, across2 - further_across2)
                score = max(score, max(down2, across2))
            diagonal = up
            left = score
            above.unsafe_ptr().unsafe_offset(at).unsafe_store(score)
            best = max(best, score)
    var out = List[Int](capacity=len(members))
    for lane in range(len(members)):
        out.append(Int(best[lane]))
    return out^


def local_scores_by_block(
    query: String, references: List[String], costs: Costs, match_score: Int, threads: Int
) -> List[Int]:
    """Every reference's best local score against the query, the references in blocks of one lane
    width, each block of references of about one length so little of it is padding, the blocks over
    `threads` threads."""
    var count = len(references)
    var scores = List[Int](length=count, fill=0)
    if count == 0:
        return scores^
    # The references by length, so a block's lanes end about together.
    comptime INDEX_BITS = 25
    var keys = List[Int](capacity=count)
    var longest = 0
    for index in range(count):
        var length = references[index].byte_length()
        longest = max(longest, length)
        keys.append((length << INDEX_BITS) | index)
    sort(keys)
    var narrow = lanes_fit(costs, match_score, query.byte_length(), longest)
    var width = 32 if narrow else 16
    var blocks = (count + width - 1) // width
    var out = scores.unsafe_ptr()
    var taken = Atomic[Int64](0)
    var two = costs.pieces() == 2

    def work(slot: Int) {mut taken, imm}:
        """Scores blocks of references one after another, each claimed from the shared counter, until none
        is left."""
        while True:
            var block = Int(taken.fetch_add(1))
            if block >= blocks:
                return
            # Each key's low bits are the reference's index, its length sorted above them.
            var members = List[Int]()
            for position in range(block * width, min((block + 1) * width, count)):
                members.append(keys[position] & ((1 << INDEX_BITS) - 1))
            var found: List[Int]
            if narrow:
                found = block_scores[DType.int16, 32, 2](
                    query.as_bytes(), references, members, costs, match_score
                ) if two else block_scores[DType.int16, 32, 1](
                    query.as_bytes(), references, members, costs, match_score
                )
            else:
                found = block_scores[DType.int32, 16, 2](
                    query.as_bytes(), references, members, costs, match_score
                ) if two else block_scores[DType.int32, 16, 1](
                    query.as_bytes(), references, members, costs, match_score
                )
            for lane in range(len(members)):
                out[unsafe_offset=members[lane]] = found[lane]

    var workers = max(min(threads, blocks), 1)
    parallelize(work, workers, workers)
    return scores^
