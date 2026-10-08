# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
One query against many references, a database search: every reference's best score or least cost,
the hits ranked, the best kept, and alignments for those alone (see `search`).

A local alignment scores many references at once, one to a SIMD lane, as SWIPE (Rognes, 2011) and
hyalite's database mode do, through the lanes a batch's local scores take (see
`lanes.lane_local_scores`): each lane a reference against the query, a group of references of about one
length costing about what one of them costs alone, every lane busy on every cell. Each pair's own sweep
keeps its lanes along one pair's anti-diagonals, which a short pair fills only in part; across references
there is no such waste. Other modes take each pair's own search, the pairs spread over the threads asked
for.
"""

from .lanes import LocalCosts, StringTexts, Texts, lane_local_scores
from .modes import Alignment, Costs


@fieldwise_init
struct Hit(Copyable, Movable, Writable):
    """One reference's result in a search: its place in the list searched, its best score (minus its
    least cost for a mode with no reward), and its alignment when asked for."""

    var index: Int
    var score: Int
    var alignment: Optional[Alignment]


@fieldwise_init
struct RepeatedText(Texts, TrivialRegisterPassable):
    """One sequence for every pair: a search's query, against each reference in turn."""

    var start: ImmPointer[UInt8, ImmUntrackedOrigin]
    var count: Int

    @always_inline
    def length(self, index: Int) -> Int:
        return self.count

    @always_inline
    def letters(self, index: Int) -> ImmPointer[UInt8, ImmUntrackedOrigin]:
        return self.start


def local_scores_by_lane(
    query: String, references: List[String], costs: Costs, match_score: Int, threads: Int
) -> List[Optional[Int]]:
    """Every reference's best local score against the query, a reference a lane, as `lanes.lane_local_scores`
    scores a batch's pairs, the references by length so a group's lanes end about together; None for one the
    lanes leave, which the caller scores alone. Past a reference's or the query's end each lane reads a byte
    no UTF-8 text holds, which matches nothing."""
    var count = len(references)
    var found = List[Int32](length=count, fill=0)
    var settled = List[Bool](length=count, fill=False)
    _ = lane_local_scores(
        count,
        StringTexts.of(references),
        RepeatedText(query.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](), query.byte_length()),
        LocalCosts.of(costs, match_score),
        (UInt8(0xFE), UInt8(0xFF)),
        max(threads, 1),
        found.unsafe_ptr(),
        settled.unsafe_ptr(),
    )
    var scores = List[Optional[Int]](capacity=count)
    for index in range(count):
        scores.append(Optional[Int](Int(found[index])) if settled[index] else None)
    return scores^
