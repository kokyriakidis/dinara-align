# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
A batch of global distances many pairs at once, a pair a 16-bit lane of one SIMD register, each over a
band of diagonals its own cost proves wide enough: the inter-sequence vectorization SeqAn's and
parasail's batch modes sweep whole matrices with, here over the band alone.

The pairs of a group share a register, one lane each, the group's references down the rows and its
queries across, their letters laid side by side, position by position, so one load holds a position of
every pair. Gotoh's recurrence then runs a row at a time over the band of diagonals, column minus row,
from `low` to `high`, every lane at once; a lane's cost is its own corner, read on its own last row.
Cells past a lane's sequences hold whatever the padding gives them, and no cell past a pair's sequences
feeds its corner. The pairs are dealt by their end diagonal, columns less rows, so a group's lanes share
a band.

A band is proven by the cost found in it. A path that visits diagonal `k` past the end diagonal and the
main one must run a gap out to `k` and another back, each at least an opening and its letters'
extensions: whatever else it does, it costs at least that (`off_band`). So a cost no dearer than any
path through the diagonals just outside the band is the least; a pair whose cost passes the cap, with
every path off the band past it too, is past it. A pair the first band does not prove is scored again
over every diagonal a cheaper path could visit, which proves itself, its cost the less of the two.

A pair goes into a lane only if no cell of its matrix can reach past what 16 bits hold (see `fits`);
the rest are left to the caller.
"""

from std.bit import count_leading_zeros
from std.math import ceildiv
from std.utils import IndexList
from std.sys import simd_width_of
from std.atomic import Atomic
from max.algorithm import parallelize

from .common import next_share
from .modes import Band, Costs, Mode


trait Texts(ImplicitlyCopyable):
    """A batch's sequences, each its bytes and their count: `List[String]`'s, or a C caller's."""

    def length(self, index: Int) -> Int:
        """Sequence `index`'s bytes."""
        ...

    def letters(self, index: Int) -> ImmPointer[UInt8, ImmUntrackedOrigin]:
        """Sequence `index`'s first byte."""
        ...


@fieldwise_init
struct StringTexts(Texts, TrivialRegisterPassable):
    """A `List[String]`'s sequences, through its storage."""

    var items: ImmPointer[String, ImmUntrackedOrigin]

    @always_inline
    def length(self, index: Int) -> Int:
        return self.items[unsafe_offset=index].byte_length()

    @always_inline
    def letters(self, index: Int) -> ImmPointer[UInt8, ImmUntrackedOrigin]:
        return self.items[unsafe_offset=index].unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]()


comptime WIDTH = simd_width_of[DType.int16]()
"""Pairs a group: the lanes of one native register of 16-bit integers, 32 under AVX-512."""
comptime Lanes = SIMD[DType.int16, WIDTH]
comptime FAR = Int16(16384)
"""A cell no path inside the band reaches; any cost a lane may hold stays under it (see `fits`)."""
comptime HELD = 16000
"""The dearest path a pair may have for a lane: FAR plus the most any chain of the recurrence adds to
it stays inside 16 bits."""


@fieldwise_init
struct LaneCosts(ImplicitlyCopyable, TrivialRegisterPassable):
    """One gap piece's costs either way, and a mismatch's: a deletion is a run of reference letters
    alone, down the rows; an insertion a run of query letters alone, across."""

    var mismatch: Int
    var deletion_opening: Int
    var deletion_extension: Int
    var insertion_opening: Int
    var insertion_extension: Int

    @staticmethod
    def of(costs: Costs, mode: Mode) -> Optional[Self]:
        """The lanes' costs for `costs`, or None where they do not serve: free ends, a match's reward, or two
        pieces."""
        if not mode.is_global() or mode.is_scored() or costs.pieces() != 1:
            return None
        return Self(costs.mismatch, costs.deletion_opening, costs.deletion_extension, costs.opening, costs.extension)

    @always_inline
    def deleted(self, letters: Int) -> Int:
        """A deletion of `letters` reference letters, nothing for none."""
        return 0 if letters <= 0 else self.deletion_opening + letters * self.deletion_extension

    @always_inline
    def inserted(self, letters: Int) -> Int:
        """An insertion of `letters` query letters, nothing for none."""
        return 0 if letters <= 0 else self.insertion_opening + letters * self.insertion_extension

    def fits(self, rows: Int, columns: Int) -> Bool:
        """Whether a pair of `rows` reference letters and `columns` query letters stays inside 16 bits:
        its dearest path, all gaps, plus an opening and a mismatch more, under `HELD`."""
        var dearest = self.deleted(rows) + self.inserted(columns)
        return dearest + max(self.deletion_opening, self.insertion_opening) + self.mismatch < HELD

    def off_band(self, diagonal: Int, end: Int) -> Int:
        """The least any path visiting `diagonal` costs, its end on diagonal `end`: a gap out to the
        diagonal from the main one's side and one back to the end's; nothing between the two."""
        if diagonal > max(0, end):
            return self.inserted(diagonal) + self.deleted(diagonal - end)
        if diagonal < min(0, end):
            return self.deleted(-diagonal) + self.inserted(end - diagonal)
        return 0


struct LaneSpace(Movable):
    """A worker's memory for its groups, kept from group to group: the letters side by side, a row of
    each layer, the group's pairs, and the pairs its first bands did not prove."""

    var row_letters: List[UInt8]
    var column_letters: List[UInt8]
    var staging: List[UInt8]
    """Each lane's text whole, a lane after another, on its way to lying side by side."""
    var scores: List[Lanes]
    var deletions: List[Lanes]
    var members: List[Int]
    var retries: List[List[Int]]
    """The pairs its first bands did not prove, filed by the width of the band that will (see
    `settle`): a pair, its first cost, and that band's lowest and highest diagonal, four numbers a pair."""

    def __init__(out self):
        self.row_letters = List[UInt8]()
        self.column_letters = List[UInt8]()
        self.staging = List[UInt8]()
        self.scores = List[Lanes]()
        self.deletions = List[Lanes]()
        self.members = List[Int]()
        self.retries = List[List[Int]]()
        for _ in range(BUCKETS):
            self.retries.append(List[Int]())


comptime BLOCK = 8
"""Positions laid side by side at once: a lane's eight letters are one 64-bit load."""


def transposed_order() -> IndexList[WIDTH * BLOCK]:
    """Where each byte of a block laid position by position comes from in one laid lane by lane."""
    var mask = IndexList[WIDTH * BLOCK]()
    for position in range(BLOCK):
        for lane in range(WIDTH):
            mask[position * WIDTH + lane] = lane * BLOCK + position
    return mask


def side_by_side[
    T: Texts
](texts: T, members: List[Int], length: Int, mut staging: List[UInt8], target: MutPointer[UInt8, _]):
    """The members' texts, position `p` of every lane at `target[p WIDTH:(p + 1) WIDTH]`: each text copied
    whole into `staging`, a lane's stretch padded to whole blocks, then a block of `BLOCK` positions at
    a time, every lane's letters of the block one load, the block turned in registers."""
    comptime order = transposed_order()
    var stride = ceildiv(max(length, 1), BLOCK) * BLOCK
    staging.resize(unsafe_uninit_length=WIDTH * stride)
    var lanes = staging.unsafe_ptr()
    for lane in range(len(members)):
        var index = members[lane]
        var count = texts.length(index)
        Span(unsafe_ptr=lanes.unsafe_offset(lane * stride), length=count).copy_from(
            Span(unsafe_ptr=texts.letters(index), length=count)
        )
    var block = Array[UInt8, WIDTH * BLOCK](fill=0)
    var block_ptr = block.unsafe_ptr()
    for start in range(0, stride, BLOCK):
        comptime for lane in range(WIDTH):
            block_ptr.unsafe_offset(lane * BLOCK).unsafe_bitcast[UInt64]().unsafe_store(
                lanes.unsafe_offset(lane * stride + start).unsafe_bitcast[UInt64]().unsafe_load()
            )
        target.unsafe_offset(start * WIDTH).unsafe_store(block_ptr.unsafe_load[width=WIDTH * BLOCK]().shuffle[order]())


def band_costs[
    T: Texts
](references: T, queries: T, mut space: LaneSpace, low: Int, high: Int, costs: LaneCosts,) -> Lanes:
    """The least cost of each of `space.members`'s pairs over the paths on diagonals `low ..= high`,
    `FAR` or more for a pair none of whose paths stays on them; a lane past the members holds nothing."""
    var count = len(space.members)
    var rows = 0
    var columns = 0
    var row_ends = SIMD[DType.int16, WIDTH](-1)
    var column_ends = SIMD[DType.int16, WIDTH](0)
    for lane in range(count):
        var index = space.members[lane]
        var reference = references.length(index)
        var query = queries.length(index)
        rows = max(rows, reference)
        columns = max(columns, query)
        row_ends[lane] = Int16(reference)
        column_ends[lane] = Int16(query)
    # Each position's letters side by side. Past a lane's own letters the rows and columns hold whatever
    # was there before: no cell past a pair's sequences feeds its corner.
    space.row_letters.resize(unsafe_uninit_length=ceildiv(max(rows, 1), BLOCK) * BLOCK * WIDTH)
    space.column_letters.resize(unsafe_uninit_length=ceildiv(max(columns, 1), BLOCK) * BLOCK * WIDTH)
    side_by_side(references, space.members, rows, space.staging, space.row_letters.unsafe_ptr())
    side_by_side(queries, space.members, columns, space.staging, space.column_letters.unsafe_ptr())

    space.scores.resize(columns + 2, Lanes(FAR))
    space.deletions.resize(columns + 2, Lanes(FAR))
    var scores = space.scores.unsafe_ptr()
    var deletions = space.deletions.unsafe_ptr()
    for column in range(columns + 2):
        scores[unsafe_offset=column] = Lanes(FAR)
        deletions[unsafe_offset=column] = Lanes(FAR)
    # Row zero inside the band: an insertion of every column before.
    for column in range(max(low, 0), min(high, columns) + 1):
        scores[unsafe_offset=column] = Lanes(Int16(costs.inserted(column)))
    var found = Lanes(FAR)
    for lane in range(count):
        var column = Int(column_ends[lane])
        if row_ends[lane] == 0 and column >= low and column <= high:
            found[lane] = scores[unsafe_offset=column][lane]

    var mismatch = Lanes(Int16(costs.mismatch))
    var delete_open = Lanes(Int16(costs.deletion_opening + costs.deletion_extension))
    var delete_extend = Lanes(Int16(costs.deletion_extension))
    var insert_open = Lanes(Int16(costs.insertion_opening + costs.insertion_extension))
    var insert_extend = Lanes(Int16(costs.insertion_extension))
    var row_source = space.row_letters.unsafe_ptr()
    var column_source = space.column_letters.unsafe_ptr()
    for row in range(1, rows + 1):
        var first = max(row + low, 0)
        var last = min(row + high, columns)
        if first > last:
            continue
        var letter = row_source.unsafe_offset((row - 1) * WIDTH).unsafe_load[width=WIDTH]()
        var left = Lanes(FAR)
        var diagonal: Lanes
        if first == 0:
            # The left edge: a deletion of every row so far.
            diagonal = scores[unsafe_offset=0]
            left = Lanes(Int16(costs.deleted(row)))
            scores[unsafe_offset=0] = left
            first = 1
        else:
            # The column before the band's first: the last row's value there, which this row's band has
            # left, so it reads as reached by nothing from here on.
            diagonal = scores[unsafe_offset=first - 1]
            scores[unsafe_offset=first - 1] = Lanes(FAR)
            deletions[unsafe_offset=first - 1] = Lanes(FAR)
        var insertion = Lanes(FAR)
        for column in range(first, last + 1):
            var other = column_source.unsafe_offset((column - 1) * WIDTH).unsafe_load[width=WIDTH]()
            var above = scores[unsafe_offset=column]
            var deletion = min(above + delete_open, deletions[unsafe_offset=column] + delete_extend)
            insertion = min(left + insert_open, insertion + insert_extend)
            var score = min(min(diagonal + letter.eq(other).select(Lanes(0), mismatch), deletion), insertion)
            deletions[unsafe_offset=column] = deletion
            scores[unsafe_offset=column] = score
            diagonal = above
            left = score
        # A lane whose reference ends on this row reads its corner, if the band holds it.
        if row_ends.eq(Int16(row)).reduce_or():
            for lane in range(count):
                if Int(row_ends[lane]) == row:
                    var column = Int(column_ends[lane])
                    if column >= max(row + low, 0) and column <= last:
                        found[lane] = scores[unsafe_offset=column][lane]
    return found


def lane_distances[
    T: Texts
](
    pairs: Int,
    references: T,
    queries: T,
    costs: LaneCosts,
    reference_band: Band,
    max_cost: Int,
    workers: Int,
    costs_out: MutPointer[Optional[Int], _],
    settled: MutPointer[Bool, _],
) -> Int:
    """Every pair 16 bits hold and `settled` does not already mark, its global distance within
    `reference_band` into `costs_out`, None past `max_cost` or with no path inside the band, and `settled`
    set; the others left for the caller. The pairs settled here.

    The band is the library's, its diagonals the reference's position less the query's; the lanes run
    the reference down the rows, their diagonals the query's position less the reference's, so they
    take it mirrored."""
    var band = Band(-reference_band.high, -reference_band.low)
    if pairs == 0:
        return 0
    var stretches = max(min(workers, pairs // WIDTH), 1)
    # Each pair's end diagonal, columns less rows; `HELD` and more for a pair the lanes cannot hold.
    var ends = List[Int](capacity=pairs)
    ends.resize(unsafe_uninit_length=pairs)
    var end_ptr = ends.unsafe_ptr()
    var extremes = List[Int](length=2 * stretches, fill=0)
    var extreme_ptr = extremes.unsafe_ptr()

    def measure(
        stretch: Int,
    ) {imm references, imm queries, imm costs, imm pairs, imm stretches, imm end_ptr, imm extreme_ptr, imm settled}:
        """Stretch `stretch`'s end diagonals, and its lowest and highest."""
        var lowest = 0
        var highest = 0
        for index in range(pairs * stretch // stretches, pairs * (stretch + 1) // stretches):
            var rows = references.length(index)
            var columns = queries.length(index)
            if not settled[unsafe_offset=index] and costs.fits(rows, columns):
                var end = columns - rows
                end_ptr[unsafe_offset=index] = end
                lowest = min(lowest, end)
                highest = max(highest, end)
            else:
                end_ptr[unsafe_offset=index] = HELD
        extreme_ptr[unsafe_offset=2 * stretch] = lowest
        extreme_ptr[unsafe_offset=2 * stretch + 1] = highest

    parallelize(measure, stretches, stretches)
    var lowest = 0
    var highest = 0
    for stretch in range(stretches):
        lowest = min(lowest, extremes[2 * stretch])
        highest = max(highest, extremes[2 * stretch + 1])
    # The pairs dealt by end diagonal, each stretch counting its own and placing them, the batch's
    # order within a diagonal.
    var span = highest - lowest + 1
    var starts = List[Int](length=stretches * span, fill=0)
    var start_ptr = starts.unsafe_ptr()

    def count(stretch: Int) {imm pairs, imm stretches, imm end_ptr, imm start_ptr, imm span, imm lowest}:
        """Stretch `stretch`'s count of each end diagonal."""
        var counts = start_ptr.unsafe_offset(stretch * span)
        for index in range(pairs * stretch // stretches, pairs * (stretch + 1) // stretches):
            var end = end_ptr[unsafe_offset=index]
            if end != HELD:
                counts[unsafe_offset=end - lowest] += 1

    parallelize(count, stretches, stretches)
    var placed = 0
    for diagonal in range(span):
        for stretch in range(stretches):
            var counted = starts[stretch * span + diagonal]
            starts[stretch * span + diagonal] = placed
            placed += counted
    var order = List[Int](capacity=max(placed, 1))
    order.resize(unsafe_uninit_length=placed)
    var order_ptr = order.unsafe_ptr()

    def place(stretch: Int) {imm pairs, imm stretches, imm end_ptr, imm start_ptr, imm span, imm lowest, imm order_ptr}:
        """Stretch `stretch`'s pairs at their places."""
        var next = start_ptr.unsafe_offset(stretch * span)
        for index in range(pairs * stretch // stretches, pairs * (stretch + 1) // stretches):
            var end = end_ptr[unsafe_offset=index]
            if end != HELD:
                order_ptr[unsafe_offset=next[unsafe_offset=end - lowest]] = index
                next[unsafe_offset=end - lowest] += 1

    parallelize(place, stretches, stretches)

    # First pass: each group over the band between its end diagonals and the main one and one diagonal
    # more either side, within `band`.
    var spaces = List[LaneSpace](capacity=workers)
    for _ in range(workers):
        spaces.append(LaneSpace())
    var space_ptr = spaces.unsafe_ptr()
    var groups = ceildiv(placed, WIDTH)
    var taken = Atomic[Int64](0)
    var first_workers = min(workers, groups)

    def first_pass(
        worker: Int,
    ) {
        mut taken,
        imm references,
        imm queries,
        imm costs,
        imm band,
        imm max_cost,
        imm groups,
        imm first_workers,
        imm placed,
        imm order_ptr,
        imm end_ptr,
        imm costs_out,
        imm settled,
        imm space_ptr,
    }:
        """Takes groups until none is left, settling every pair its band proves."""
        ref space = space_ptr[unsafe_offset=worker]
        var last_share = 0
        while True:
            var share = next_share(taken, groups, first_workers, last_share)
            if share[0] >= groups:
                return
            for group in range(share[0], share[1]):
                space.members.clear()
                var low = 0
                var high = 0
                for slot in range(group * WIDTH, min(placed, (group + 1) * WIDTH)):
                    var index = order_ptr[unsafe_offset=slot]
                    space.members.append(index)
                    var end = end_ptr[unsafe_offset=index]
                    low = min(low, min(0, end) - 1)
                    high = max(high, max(0, end) + 1)
                low = max(low, band.low)
                high = min(high, band.high)
                var found = band_costs(references, queries, space, low, high, costs)
                for lane in range(len(space.members)):
                    var index = space.members[lane]
                    settle(
                        index,
                        Int(found[lane]),
                        references.length(index),
                        queries.length(index),
                        low,
                        high,
                        band,
                        costs,
                        max_cost,
                        costs_out,
                        settled,
                        space.retries,
                    )

    parallelize(first_pass, first_workers, first_workers)

    # Second pass: every pair the first band did not prove, over every diagonal a cheaper path than its
    # first cost could visit, which proves itself; gathered from the workers' buckets by the width of
    # that band, so a group's lanes need bands alike.
    var retries = List[Int]()
    for bucket in range(BUCKETS):
        for worker in range(workers):
            retries.extend(Span(spaces[worker].retries[bucket]))
    var unproven = len(retries) // 4
    if unproven == 0:
        return placed
    var retry_ptr = retries.unsafe_ptr()
    var retry_groups = ceildiv(unproven, WIDTH)
    var retaken = Atomic[Int64](0)
    var second_workers = min(workers, retry_groups)

    def second_pass(
        worker: Int,
    ) {
        mut retaken,
        imm references,
        imm queries,
        imm costs,
        imm max_cost,
        imm retry_groups,
        imm second_workers,
        imm unproven,
        imm retry_ptr,
        imm costs_out,
        imm settled,
        imm space_ptr,
    }:
        """Takes groups of unproven pairs until none is left, settling each."""
        ref space = space_ptr[unsafe_offset=worker]
        var last_share = 0
        while True:
            var share = next_share(retaken, retry_groups, second_workers, last_share)
            if share[0] >= retry_groups:
                return
            for group in range(share[0], share[1]):
                space.members.clear()
                var low = 0
                var high = 0
                for slot in range(group * WIDTH, min(unproven, (group + 1) * WIDTH)):
                    space.members.append(retry_ptr[unsafe_offset=4 * slot])
                    low = min(low, retry_ptr[unsafe_offset=4 * slot + 2])
                    high = max(high, retry_ptr[unsafe_offset=4 * slot + 3])
                var found = band_costs(references, queries, space, low, high, costs)
                for lane in range(len(space.members)):
                    var slot = group * WIDTH + lane
                    # A path as cheap as the first band's cost the first band found already.
                    var cost = min(Int(found[lane]), retry_ptr[unsafe_offset=4 * slot + 1])
                    costs_out[unsafe_offset=space.members[lane]] = Optional[Int](cost) if cost <= max_cost else None
                    settled[unsafe_offset=space.members[lane]] = True

    parallelize(second_pass, second_workers, second_workers)
    return placed


comptime BUCKETS = 24
"""The widths of second-pass bands a worker files apart, by their bit length: up to 2^23 diagonals."""


@always_inline
def settle(
    index: Int,
    found: Int,
    rows: Int,
    columns: Int,
    low: Int,
    high: Int,
    band: Band,
    costs: LaneCosts,
    max_cost: Int,
    costs_out: MutPointer[Optional[Int], _],
    settled: MutPointer[Bool, _],
    mut retries: List[List[Int]],
):
    """Settles pair `index`, of `rows` reference letters and `columns` query letters, with the cost
    `found` on diagonals `low ..= high` if no path off them could beat it, or one past the cap could not
    come under it; else files it, its cost and the band a cheaper path needs, by that band's width. A
    path cannot leave the matrix, nor the band."""
    var end = columns - rows
    if end < band.low or end > band.high or 0 < band.low or 0 > band.high:
        # Its start or its end lies outside the band: no alignment inside it.
        costs_out[unsafe_offset=index] = None
        settled[unsafe_offset=index] = True
        return
    var top = min(band.high, columns)
    var bottom = max(band.low, -rows)
    var beyond_high = costs.off_band(high + 1, end) if high < top else Int.MAX
    var beyond_low = costs.off_band(low - 1, end) if low > bottom else Int.MAX
    var off = min(beyond_high, beyond_low)
    if found <= off or off > max_cost:
        # The least cost, or, its band's past the cap as is every path off it, past the cap.
        costs_out[unsafe_offset=index] = Optional[Int](found) if found <= max_cost else None
        settled[unsafe_offset=index] = True
        return
    # Every diagonal a path cheaper than the target could visit, within the band and the matrix: a path
    # dearer than the cap need not be found, only shown past it.
    var target = found if found <= max_cost else max_cost + 1
    var wide_high = high
    while wide_high < top and costs.off_band(wide_high + 1, end) < target:
        wide_high += 1
    var wide_low = low
    while wide_low > bottom and costs.off_band(wide_low - 1, end) < target:
        wide_low -= 1
    var width_bits = 64 - Int(count_leading_zeros(UInt64(wide_high - wide_low)))
    ref bucket = retries[min(width_bits, BUCKETS - 1)]
    bucket.append(index)
    bucket.append(found)
    bucket.append(wide_low)
    bucket.append(wide_high)
