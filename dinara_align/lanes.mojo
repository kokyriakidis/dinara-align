# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
A batch of global distances many pairs at once, a pair a lane of one SIMD register, each over a band
of diagonals its own cost proves wide enough: the inter-sequence vectorization SeqAn's and parasail's
batch modes sweep whole matrices with, here over the band alone.

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

The lanes come in two sizes, as SSW's and parasail's do. First every pair goes into a byte, 64 to an
AVX-512 register, whose additions saturate: a cost under 255 is exact, 255 is 255 or more, and stands
for a cell no path inside the band reaches, so a pair whose cost passes it is left. Then every pair left
goes into 16 bits, whose dearest path must stay under what they hold (see `fits`); the rest are left to
the caller.
"""

from std.bit import count_leading_zeros
from std.math import ceildiv
from std.memory import bitcast
from std.utils import IndexList
from std.sys import llvm_intrinsic, simd_width_of
from std.atomic import Atomic
from max.algorithm import parallelize

from .common import next_share
from .gap_affine import ALIGNED, FIRST_GAP, SECOND_GAP
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


@always_inline
def lanes_of[value: DType]() -> Int:
    """Pairs a group: the lanes of one native register of `value`s, 64 bytes or 32 16-bit integers
    under AVX-512."""
    return simd_width_of[value]()


@always_inline
def far_of[value: DType]() -> Int:
    """A cell no path inside the band reaches: a byte's 255, which its additions saturate at; 16 bits'
    16384, which any cost a lane may hold stays under (see `fits`)."""
    return 255 if value == DType.uint8 else 16384


@always_inline
def added[value: DType, width: Int](left: SIMD[value, width], right: SIMD[value, width]) -> SIMD[value, width]:
    """A sum of costs: a byte's saturating at 255, 16 bits' plain."""
    comptime if value == DType.uint8:
        return llvm_intrinsic["llvm.uadd.sat", SIMD[value, width]](left, right)
    else:
        return left + right


comptime HELD = 16000
"""The dearest path a pair may have for a lane: FAR plus the most any chain of the recurrence adds to
it stays inside 16 bits."""


@fieldwise_init
struct LaneCosts(ImplicitlyCopyable, TrivialRegisterPassable):
    """One or two gap pieces' costs either way, and a mismatch's: a deletion is a run of reference letters
    alone, down the rows; an insertion a run of query letters alone, across. A gap costs the cheaper of
    its direction's pieces; with one piece the second's fields are unused."""

    var mismatch: Int
    var deletion_opening: Int
    var deletion_extension: Int
    var insertion_opening: Int
    var insertion_extension: Int
    var pieces: Int
    var deletion_opening2: Int
    var deletion_extension2: Int
    var insertion_opening2: Int
    var insertion_extension2: Int

    @staticmethod
    def one_piece(
        mismatch: Int, deletion_opening: Int, deletion_extension: Int, insertion_opening: Int, insertion_extension: Int
    ) -> Self:
        """Costs of one gap piece."""
        return Self(
            mismatch, deletion_opening, deletion_extension, insertion_opening, insertion_extension, 1, 0, 0, 0, 0
        )

    @staticmethod
    def of(costs: Costs, mode: Mode) -> Optional[Self]:
        """The lanes' costs for `costs`, or None where they do not serve: free ends or a match's reward."""
        if not mode.is_global() or mode.is_scored():
            return None
        return Self(
            costs.mismatch,
            costs.deletion_opening,
            costs.deletion_extension,
            costs.opening,
            costs.extension,
            costs.pieces(),
            costs.deletion_opening2,
            costs.deletion_extension2,
            costs.opening2,
            costs.extension2,
        )

    @always_inline
    def deleted(self, letters: Int) -> Int:
        """A deletion of `letters` reference letters at its cheaper piece, nothing for none."""
        if letters <= 0:
            return 0
        var cost = self.deletion_opening + letters * self.deletion_extension
        if self.pieces == 2:
            cost = min(cost, self.deletion_opening2 + letters * self.deletion_extension2)
        return cost

    @always_inline
    def inserted(self, letters: Int) -> Int:
        """An insertion of `letters` query letters at its cheaper piece, nothing for none."""
        if letters <= 0:
            return 0
        var cost = self.insertion_opening + letters * self.insertion_extension
        if self.pieces == 2:
            cost = min(cost, self.insertion_opening2 + letters * self.insertion_extension2)
        return cost

    def dearest_step(self) -> Int:
        """The most any one step of the recurrence adds: a mismatch, or a gap's opening letter or a further one."""
        var dearest = max(self.mismatch, self.deletion_opening + self.deletion_extension)
        dearest = max(dearest, self.insertion_opening + self.insertion_extension)
        if self.pieces == 2:
            dearest = max(dearest, self.deletion_opening2 + self.deletion_extension2)
            dearest = max(dearest, self.insertion_opening2 + self.insertion_extension2)
        return dearest

    def fits_bytes(self) -> Bool:
        """Whether every step costs under a byte's 255, which a pair's lanes need to count anything."""
        return self.dearest_step() < 255

    def fits(self, rows: Int, columns: Int) -> Bool:
        """Whether a pair of `rows` reference letters and `columns` query letters stays inside 16 bits: a
        run of every row and every column at each piece, and a step more, under `HELD`; a layer starting
        from FAR grows no faster."""
        var dearest = self.deletion_opening + rows * self.deletion_extension
        dearest += self.insertion_opening + columns * self.insertion_extension
        if self.pieces == 2:
            dearest = max(
                dearest,
                self.deletion_opening2
                + rows * self.deletion_extension2
                + self.insertion_opening2
                + columns * self.insertion_extension2,
            )
        return dearest + self.dearest_step() < HELD

    def off_band(self, diagonal: Int, end: Int) -> Int:
        """The least any path visiting `diagonal` costs, its end on diagonal `end`: a gap out to the
        diagonal from the main one's side and one back to the end's; nothing between the two."""
        if diagonal > max(0, end):
            return self.inserted(diagonal) + self.deleted(diagonal - end)
        if diagonal < min(0, end):
            return self.deleted(-diagonal) + self.inserted(end - diagonal)
        return 0


struct LaneSpace[value: DType](Movable):
    """A worker's memory for its groups, kept from group to group: the letters side by side, a row of
    each layer, the group's pairs, and the pairs its first bands did not prove."""

    var row_letters: List[UInt8]
    var column_letters: List[UInt8]
    var staging: List[UInt8]
    """Each lane's text whole, a lane after another, on its way to lying side by side."""
    var scores: List[SIMD[Self.value, lanes_of[Self.value]()]]
    var deletions: List[SIMD[Self.value, lanes_of[Self.value]()]]
    var deletions2: List[SIMD[Self.value, lanes_of[Self.value]()]]
    """The second gap piece's deletions, where the costs have one."""
    var flags: List[SIMD[DType.uint8, lanes_of[Self.value]()]]
    """For an alignment, each cell's flag on the band, a row after another (see `traced`)."""
    var members: List[Int]
    var retries: List[List[Int]]
    """The pairs its first bands did not prove, filed by the width of the band that will (see
    `proof`): a pair, its first cost, and that band's lowest and highest diagonal, four numbers a pair."""

    def __init__(out self):
        self.row_letters = List[UInt8]()
        self.column_letters = List[UInt8]()
        self.staging = List[UInt8]()
        self.scores = List[SIMD[Self.value, lanes_of[Self.value]()]]()
        self.deletions = List[SIMD[Self.value, lanes_of[Self.value]()]]()
        self.deletions2 = List[SIMD[Self.value, lanes_of[Self.value]()]]()
        self.flags = List[SIMD[DType.uint8, lanes_of[Self.value]()]]()
        self.members = List[Int]()
        self.retries = List[List[Int]]()
        for _ in range(BUCKETS):
            self.retries.append(List[Int]())


comptime BLOCK = 8
"""Positions laid side by side at once: a lane's eight letters are one 64-bit load."""


def byte_order() -> IndexList[64]:
    """An 8 by 8 transpose of bytes inside 64: eight lanes' eight letters in, eight positions' out."""
    var mask = IndexList[64]()
    for position in range(BLOCK):
        for lane in range(8):
            mask[position * 8 + lane] = lane * BLOCK + position
    return mask


def word_order[groups: Int]() -> IndexList[groups * BLOCK]:
    """Words of `groups` groups of eight lanes, a group's eight positions in turn, to each position's
    words, a group's after another."""
    var mask = IndexList[groups * BLOCK]()
    for position in range(BLOCK):
        for group in range(groups):
            mask[position * groups + group] = group * BLOCK + position
    return mask


def side_by_side[
    T: Texts, width: Int
](
    texts: T,
    members: List[Int],
    length: Int,
    mut staging: List[UInt8],
    target: MutPointer[UInt8, _],
    reverse: Bool = False,
):
    """The members' texts, each back to front with `reverse`, position `p` of every lane at
    `target[p width:(p + 1) width]`: each text copied whole into `staging`, a lane's stretch padded to
    whole blocks, then a block of `BLOCK` positions at a time, every lane's letters of the block one load,
    the block turned in registers in two steps: eight lanes' bytes inside 64 at a time, then the 64-bit
    words of all of them."""
    comptime groups = width // 8
    comptime bytes = byte_order()
    comptime words = word_order[groups]()
    var stride = ceildiv(max(length, 1), BLOCK) * BLOCK
    staging.resize(unsafe_uninit_length=width * stride)
    var lanes = staging.unsafe_ptr()
    for lane in range(len(members)):
        var index = members[lane]
        var count = texts.length(index)
        if reverse:
            var source = texts.letters(index)
            var destination = lanes.unsafe_offset(lane * stride)
            comptime CHUNK = 16
            var position = 0
            while position + CHUNK <= count:
                destination.unsafe_offset(position).unsafe_store(
                    source.unsafe_offset(count - position - CHUNK).unsafe_load[width=CHUNK]().reversed()
                )
                position += CHUNK
            while position < count:
                destination[unsafe_offset=position] = source[unsafe_offset=count - 1 - position]
                position += 1
        else:
            Span(unsafe_ptr=lanes.unsafe_offset(lane * stride), length=count).copy_from(
                Span(unsafe_ptr=texts.letters(index), length=count)
            )
    var block = Array[UInt8, width * BLOCK](fill=0)
    var block_ptr = block.unsafe_ptr()
    for start in range(0, stride, BLOCK):
        comptime for lane in range(width):
            block_ptr.unsafe_offset(lane * BLOCK).unsafe_bitcast[UInt64]().unsafe_store(
                lanes.unsafe_offset(lane * stride + start).unsafe_bitcast[UInt64]().unsafe_load()
            )
        var turned = SIMD[DType.uint64, groups * BLOCK]()
        comptime for group in range(groups):
            var inner = bitcast[DType.uint64, BLOCK](
                block_ptr.unsafe_offset(group * 64).unsafe_load[width=64]().shuffle[bytes]()
            )
            comptime for position in range(BLOCK):
                turned[group * BLOCK + position] = inner[position]
        target.unsafe_offset(start * width).unsafe_bitcast[UInt64]().unsafe_store(turned.shuffle[words]())


comptime SECOND_PIECE_DELETION = 3
"""The layer of the second gap piece's deletions, `gap_affine.gap_layer(1, True)`."""
comptime SECOND_PIECE_INSERTION = 4
"""The layer of the second gap piece's insertions, `gap_affine.gap_layer(1, False)`."""
comptime SOURCE_MASK = UInt8(7)
"""A flag's low three bits: the source its alignment layer takes."""


@always_inline
def extended_bit(layer: Int) -> UInt8:
    """A flag bit: gap layer `layer` extends, not opens, into the cell."""
    return UInt8(1) << UInt8(2 + layer)


def band_costs[
    T: Texts, value: DType, record: Bool = False
](
    references: T,
    queries: T,
    mut space: LaneSpace[value],
    low: Int,
    high: Int,
    costs: LaneCosts,
    reverse: Bool = False,
) -> SIMD[value, lanes_of[value]()]:
    """`band_costs` at the costs' own number of gap pieces."""
    if costs.pieces == 2:
        return pieced_band_costs[T, value, 2, record](references, queries, space, low, high, costs, reverse)
    return pieced_band_costs[T, value, 1, record](references, queries, space, low, high, costs, reverse)


def pieced_band_costs[
    T: Texts, value: DType, pieces: Int, record: Bool
](
    references: T,
    queries: T,
    mut space: LaneSpace[value],
    low: Int,
    high: Int,
    costs: LaneCosts,
    reverse: Bool,
) -> SIMD[value, lanes_of[value]()]:
    """The least cost of each of `space.members`'s pairs over the paths on diagonals `low ..= high`,
    `far_of[value]()` or more for a pair none of whose paths stays on them, or, in a byte, whose cost
    reaches 255; a lane past the members holds nothing. Both sequences run back to front with `reverse`.

    With `record`, each cell of the band keeps a flag in `space.flags`, for `traced`: the source its
    alignment layer takes, `ALIGNED` for the diagonal, and a bit a gap layer, `extended_bit`, where the
    layer extends. The source is the one the wavefront's backtrace takes (see `Ties`): of the moves into
    the cell as cheap as it, a substitution, then a letter of the reference alone, then one of the query,
    the second gap piece before the first, and only then a match; a gap layer's extension before its
    opening."""
    comptime WIDTH = lanes_of[value]()
    comptime Lanes = SIMD[value, WIDTH]
    comptime FAR = far_of[value]()

    @always_inline
    def held(cost: Int) -> Scalar[value]:
        """A cost as a lane holds it, 255 or more a byte's 255."""
        return Scalar[value](min(cost, FAR))

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
    side_by_side[T, WIDTH](references, space.members, rows, space.staging, space.row_letters.unsafe_ptr(), reverse)
    side_by_side[T, WIDTH](queries, space.members, columns, space.staging, space.column_letters.unsafe_ptr(), reverse)
    comptime Flags = SIMD[DType.uint8, WIDTH]
    var span = high - low + 1
    comptime if record:
        space.flags.resize(unsafe_uninit_length=rows * span)
    var flags = space.flags.unsafe_ptr()

    space.scores.resize(columns + 2, Lanes(FAR))
    space.deletions.resize(columns + 2, Lanes(FAR))
    var scores = space.scores.unsafe_ptr()
    var deletions = space.deletions.unsafe_ptr()
    for column in range(columns + 2):
        scores[unsafe_offset=column] = Lanes(FAR)
        deletions[unsafe_offset=column] = Lanes(FAR)
    comptime if pieces == 2:
        space.deletions2.resize(columns + 2, Lanes(FAR))
        for column in range(columns + 2):
            space.deletions2[column] = Lanes(FAR)
    var deletions2 = space.deletions2.unsafe_ptr()
    # Row zero inside the band: an insertion of every column before.
    for column in range(max(low, 0), min(high, columns) + 1):
        scores[unsafe_offset=column] = Lanes(held(costs.inserted(column)))
    var found = Lanes(FAR)
    for lane in range(count):
        var column = Int(column_ends[lane])
        if row_ends[lane] == 0 and column >= low and column <= high:
            found[lane] = scores[unsafe_offset=column][lane]

    var mismatch = Lanes(held(costs.mismatch))
    var delete_open = Lanes(held(costs.deletion_opening + costs.deletion_extension))
    var delete_extend = Lanes(held(costs.deletion_extension))
    var insert_open = Lanes(held(costs.insertion_opening + costs.insertion_extension))
    var insert_extend = Lanes(held(costs.insertion_extension))
    var delete_open2 = Lanes(held(costs.deletion_opening2 + costs.deletion_extension2))
    var delete_extend2 = Lanes(held(costs.deletion_extension2))
    var insert_open2 = Lanes(held(costs.insertion_opening2 + costs.insertion_extension2))
    var insert_extend2 = Lanes(held(costs.insertion_extension2))
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
            left = Lanes(held(costs.deleted(row)))
            scores[unsafe_offset=0] = left
            first = 1
        else:
            # The column before the band's first: the last row's value there, which this row's band has
            # left, so it reads as reached by nothing from here on.
            diagonal = scores[unsafe_offset=first - 1]
            scores[unsafe_offset=first - 1] = Lanes(FAR)
            deletions[unsafe_offset=first - 1] = Lanes(FAR)
            comptime if pieces == 2:
                deletions2[unsafe_offset=first - 1] = Lanes(FAR)
        var insertion = Lanes(FAR)
        var insertion2 = Lanes(FAR)
        # This row's flags, indexed by column.
        var row_flags = flags.unsafe_offset((row - 1) * span - row - low)
        for column in range(first, last + 1):
            var other = column_source.unsafe_offset((column - 1) * WIDTH).unsafe_load[width=WIDTH]()
            var above = scores[unsafe_offset=column]
            var equal = letter.eq(other)
            var substituted = added(diagonal, equal.select(Lanes(0), mismatch))
            var opened_deletion = added(above, delete_open)
            var extended_deletion = added(deletions[unsafe_offset=column], delete_extend)
            var deletion = min(opened_deletion, extended_deletion)
            var opened_insertion = added(left, insert_open)
            var extended_insertion = added(insertion, insert_extend)
            insertion = min(opened_insertion, extended_insertion)
            var score = min(min(substituted, deletion), insertion)
            var deletion2 = Lanes(FAR)
            var opened_deletion2 = Lanes(FAR)
            var extended_deletion2 = Lanes(FAR)
            var opened_insertion2 = Lanes(FAR)
            var extended_insertion2 = Lanes(FAR)
            comptime if pieces == 2:
                opened_deletion2 = added(above, delete_open2)
                extended_deletion2 = added(deletions2[unsafe_offset=column], delete_extend2)
                deletion2 = min(opened_deletion2, extended_deletion2)
                opened_insertion2 = added(left, insert_open2)
                extended_insertion2 = added(insertion2, insert_extend2)
                insertion2 = min(opened_insertion2, extended_insertion2)
                score = min(score, min(deletion2, insertion2))
                deletions2[unsafe_offset=column] = deletion2
            comptime if record:
                # The lowest priority first, each source as cheap as the cell taking over from the last.
                var flag = insertion.eq(score).select(Flags(UInt8(SECOND_GAP)), Flags(UInt8(ALIGNED)))
                comptime if pieces == 2:
                    flag = insertion2.eq(score).select(Flags(UInt8(SECOND_PIECE_INSERTION)), flag)
                flag = deletion.eq(score).select(Flags(UInt8(FIRST_GAP)), flag)
                comptime if pieces == 2:
                    flag = deletion2.eq(score).select(Flags(UInt8(SECOND_PIECE_DELETION)), flag)
                flag = (~equal & substituted.eq(score)).select(Flags(UInt8(ALIGNED)), flag)
                flag |= extended_deletion.le(opened_deletion).select(Flags(extended_bit(FIRST_GAP)), Flags(0))
                flag |= extended_insertion.le(opened_insertion).select(Flags(extended_bit(SECOND_GAP)), Flags(0))
                comptime if pieces == 2:
                    flag |= extended_deletion2.le(opened_deletion2).select(
                        Flags(extended_bit(SECOND_PIECE_DELETION)), Flags(0)
                    )
                    flag |= extended_insertion2.le(opened_insertion2).select(
                        Flags(extended_bit(SECOND_PIECE_INSERTION)), Flags(0)
                    )
                row_flags[unsafe_offset=column] = flag
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
    """Every pair `settled` does not already mark that the lanes hold, its global distance within
    `reference_band` into `costs_out`, None past `max_cost` or with no path inside the band, and `settled`
    set; the others left for the caller. The pairs settled here: first in bytes, where every step's
    cost fits one, then in 16 bits."""
    var before = 0
    for index in range(pairs):
        before += Int(settled[unsafe_offset=index])
    if costs.fits_bytes():
        lane_stage[T, DType.uint8](
            pairs, references, queries, costs, reference_band, max_cost, workers, costs_out, settled
        )
    lane_stage[T, DType.int16](pairs, references, queries, costs, reference_band, max_cost, workers, costs_out, settled)
    # A stage leaves what its lanes cannot hold, so the pairs settled are counted, not assumed.
    var after = 0
    for index in range(pairs):
        after += Int(settled[unsafe_offset=index])
    return after - before


def dealt[
    T: Texts, value: DType
](
    pairs: Int,
    references: T,
    queries: T,
    costs: LaneCosts,
    workers: Int,
    settled: MutPointer[Bool, _],
    covering: Band,
    whole: Bool,
) -> List[Int]:
    """The pairs `settled` does not mark that lanes of `value` hold, dealt by their end diagonal, columns
    less rows, the batch's order within one; with `whole`, only those both of whose sequences hold a
    letter and whose every diagonal `covering` holds."""
    comptime WIDTH = lanes_of[value]()
    var stretches = max(min(workers, pairs // WIDTH), 1)
    # Each pair's end diagonal; `HELD` and more for a pair the lanes cannot hold.
    var ends = List[Int](capacity=pairs)
    ends.resize(unsafe_uninit_length=pairs)
    var end_ptr = ends.unsafe_ptr()
    var extremes = List[Int](length=2 * stretches, fill=0)
    var extreme_ptr = extremes.unsafe_ptr()

    def measure(
        stretch: Int,
    ) {
        imm references,
        imm queries,
        imm costs,
        imm pairs,
        imm stretches,
        imm end_ptr,
        imm extreme_ptr,
        imm settled,
        imm covering,
        imm whole,
    }:
        """Stretch `stretch`'s end diagonals, and its lowest and highest."""
        var lowest = 0
        var highest = 0
        for index in range(pairs * stretch // stretches, pairs * (stretch + 1) // stretches):
            var rows = references.length(index)
            var columns = queries.length(index)
            # A byte holds any pair, its saturation telling; 16 bits only a pair whose costs fit.
            var held = value == DType.uint8 or costs.fits(rows, columns)
            if whole:
                held = held and rows > 0 and columns > 0 and covering.covers(columns, rows)
            if not settled[unsafe_offset=index] and held:
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
    # Each stretch counts its own pairs of each end diagonal and places them.
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
    return order^


def lane_stage[
    T: Texts, value: DType
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
):
    """`lane_distances` in lanes of `value`: every pair `settled` does not mark and the lanes hold settled,
    a pair whose byte saturates left.

    The band is the library's, its diagonals the reference's position less the query's; the lanes run
    the reference down the rows, their diagonals the query's position less the reference's, so they
    take it mirrored."""
    comptime WIDTH = lanes_of[value]()
    var band = Band(-reference_band.high, -reference_band.low)
    if pairs == 0:
        return
    var order = dealt[T, value](pairs, references, queries, costs, workers, settled, Band(), False)
    var placed = len(order)
    var order_ptr = order.unsafe_ptr()

    # First pass: each group over the band between its end diagonals and the main one and one diagonal
    # more either side, within `band`.
    var spaces = List[LaneSpace[value]](capacity=workers)
    for _ in range(workers):
        spaces.append(LaneSpace[value]())
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
                    var end = queries.length(index) - references.length(index)
                    low = min(low, min(0, end) - 1)
                    high = max(high, max(0, end) + 1)
                low = max(low, band.low)
                high = min(high, band.high)
                var found = band_costs(references, queries, space, low, high, costs)
                for lane in range(len(space.members)):
                    var index = space.members[lane]
                    var cost = Int(found[lane])
                    var verdict = proof[value, False](
                        index,
                        cost,
                        references.length(index),
                        queries.length(index),
                        low,
                        high,
                        band,
                        costs,
                        max_cost,
                        space.retries,
                    )
                    if verdict == PROVEN or verdict == REFUSED:
                        var kept = verdict == PROVEN and cost <= max_cost
                        costs_out[unsafe_offset=index] = Optional[Int](cost) if kept else None
                        settled[unsafe_offset=index] = True

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
        return
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
                    # A path as cheap as the first band's cost the first band found already; one cheaper is
                    # under 255, so no byte saturates on it.
                    var cost = min(Int(found[lane]), retry_ptr[unsafe_offset=4 * slot + 1])
                    costs_out[unsafe_offset=space.members[lane]] = Optional[Int](cost) if cost <= max_cost else None
                    settled[unsafe_offset=space.members[lane]] = True

    parallelize(second_pass, second_workers, second_workers)


comptime BUCKETS = 24
"""The widths of second-pass bands a worker files apart, by their bit length: up to 2^23 diagonals."""


comptime PROVEN = 0
"""`proof`'s verdict: the cost found is the pair's, or, past the cap, shows the pair past it."""
comptime REFUSED = 1
"""`proof`'s verdict: no alignment of the pair lies inside the band."""
comptime UNHELD = 2
"""`proof`'s verdict: a byte saturated, the pair left for 16 bits."""
comptime FILED = 3
"""`proof`'s verdict: filed for a wider band."""


@always_inline
def proof[
    value: DType, strict: Bool
](
    index: Int,
    found: Int,
    rows: Int,
    columns: Int,
    low: Int,
    high: Int,
    band: Band,
    costs: LaneCosts,
    max_cost: Int,
    mut retries: List[List[Int]],
) -> Int:
    """Whether the cost `found` on diagonals `low ..= high` settles pair `index`, of `rows` reference letters
    and `columns` query letters: proven when no path off them could beat it, or one past the cap could
    not come under it; else the pair is filed, its cost and the band a cheaper path needs, by that band's
    width. With `strict`, as an alignment needs, no path off them may even match it: the tie rule's
    path, an optimal one, then lies inside. A path cannot leave the matrix, nor the band. A byte's 255
    may be any cost from there up, or none inside the band: such a pair is left for 16 bits."""
    comptime if value == DType.uint8:
        if found >= far_of[value]():
            return UNHELD
    var end = columns - rows
    if end < band.low or end > band.high or 0 < band.low or 0 > band.high:
        # Its start or its end lies outside the band: no alignment inside it.
        return REFUSED
    var top = min(band.high, columns)
    var bottom = max(band.low, -rows)
    var beyond_high = costs.off_band(high + 1, end) if high < top else Int.MAX
    var beyond_low = costs.off_band(low - 1, end) if low > bottom else Int.MAX
    var off = min(beyond_high, beyond_low)
    var beaten = found < off if strict else found <= off
    if beaten or off > max_cost:
        # The least cost, or, its band's past the cap as is every path off it, past the cap.
        return PROVEN
    # Every diagonal a path cheaper than the target could visit, within the band and the matrix: a path
    # dearer than the cap need not be found, only shown past it.
    var target = min(found, max_cost) + 1 if strict else (found if found <= max_cost else max_cost + 1)
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
    return FILED


def reverse_bytes(bytes: MutPointer[UInt8, _], count: Int):
    """Turns `count` bytes back to front in place, sixteen from either end at a time."""
    comptime CHUNK = 16
    var front = 0
    var back = count
    while back - front >= 2 * CHUNK:
        var head = bytes.unsafe_offset(front).unsafe_load[width=CHUNK]()
        var tail = bytes.unsafe_offset(back - CHUNK).unsafe_load[width=CHUNK]()
        bytes.unsafe_offset(front).unsafe_store(tail.reversed())
        bytes.unsafe_offset(back - CHUNK).unsafe_store(head.reversed())
        front += CHUNK
        back -= CHUNK
    while back - front >= 2:
        back -= 1
        var kept = bytes[unsafe_offset=front]
        bytes[unsafe_offset=front] = bytes[unsafe_offset=back]
        bytes[unsafe_offset=back] = kept
        front += 1


def traced[
    value: DType
](space: LaneSpace[value], lane: Int, rows: Int, columns: Int, low: Int, high: Int, mut moves: List[UInt8]):
    """Appends, from the corner back, lane `lane`'s path through the flags `band_costs` kept over diagonals
    `low ..= high`, a pair of `rows` reference letters and `columns` query letters: at each cell the source
    its flag names, a gap layer's run until it opens. A path on the first row or column has one way left,
    a gap along it."""
    comptime WIDTH = lanes_of[value]()
    var span = high - low + 1
    var flags = space.flags.unsafe_ptr().unsafe_bitcast[UInt8]().unsafe_offset(lane)
    # The moves go straight into the list's room, at most a letter of either sequence each.
    var start = len(moves)
    moves.reserve(start + rows + columns)
    var out = moves.unsafe_ptr().unsafe_offset(start)
    var written = 0
    var row = rows
    var column = columns
    # The cell's flag, `WIDTH` bytes apart: a row up is a span less a diagonal, a column left a diagonal.
    var cell = ((row - 1) * span + column - row - low) * WIDTH
    var layer = ALIGNED
    while row > 0 and column > 0:
        var flag = flags[unsafe_offset=cell]
        if layer == ALIGNED:
            layer = Int(flag & SOURCE_MASK)
            if layer == ALIGNED:
                out[unsafe_offset=written] = UInt8(ALIGNED)
                written += 1
                row -= 1
                column -= 1
                cell -= span * WIDTH
            continue
        var extended = flag & extended_bit(layer) != 0
        if layer % 2 == 1:
            # A deletion: a reference letter alone, down a row.
            out[unsafe_offset=written] = UInt8(FIRST_GAP)
            row -= 1
            cell -= (span - 1) * WIDTH
        else:
            out[unsafe_offset=written] = UInt8(SECOND_GAP)
            column -= 1
            cell -= WIDTH
        written += 1
        if not extended:
            layer = ALIGNED
    for _ in range(row):
        out[unsafe_offset=written] = UInt8(FIRST_GAP)
        written += 1
    for _ in range(column):
        out[unsafe_offset=written] = UInt8(SECOND_GAP)
        written += 1
    moves.resize(unsafe_uninit_length=start + written)


def lane_alignments[
    T: Texts
](
    pairs: Int,
    references: T,
    queries: T,
    costs: LaneCosts,
    reference_band: Band,
    max_cost: Int,
    left: Bool,
    workers: Int,
    budget: Int,
    costs_out: MutPointer[Optional[Int], _],
    moves_out: MutPointer[List[UInt8], _],
    settled: MutPointer[Bool, _],
) -> Int:
    """`lane_distances` for alignments: every pair `settled` does not mark that the lanes hold, its global
    cost within `reference_band` into `costs_out` and its path's moves, right to left as `gap_affine.solve`
    appends them, into `moves_out`, None past `max_cost` or with no path inside the band, and `settled`
    set. The path is the one `Ties` picks, `Ties.LEFT` with `left`: the rule run over both sequences
    reversed, which takes only pairs the band leaves whole. A group whose flags would pass `budget` bytes
    is left to the caller, as is a pair with an empty side."""
    var before = 0
    for index in range(pairs):
        before += Int(settled[unsafe_offset=index])
    if costs.fits_bytes():
        lane_alignment_stage[T, DType.uint8](
            pairs,
            references,
            queries,
            costs,
            reference_band,
            max_cost,
            left,
            workers,
            budget,
            costs_out,
            moves_out,
            settled,
        )
    lane_alignment_stage[T, DType.int16](
        pairs,
        references,
        queries,
        costs,
        reference_band,
        max_cost,
        left,
        workers,
        budget,
        costs_out,
        moves_out,
        settled,
    )
    var after = 0
    for index in range(pairs):
        after += Int(settled[unsafe_offset=index])
    return after - before


@always_inline
def settle_traced[
    T: Texts, value: DType
](
    space: LaneSpace[value],
    lane: Int,
    index: Int,
    cost: Int,
    low: Int,
    high: Int,
    references: T,
    queries: T,
    max_cost: Int,
    left: Bool,
    costs_out: MutPointer[Optional[Int], _],
    moves_out: MutPointer[List[UInt8], _],
    settled: MutPointer[Bool, _],
):
    """Settles pair `index`, in lane `lane` of the band `low ..= high` just swept, with its least cost
    `cost`, and its path traced (see `traced`) if under `max_cost`."""
    if cost <= max_cost:
        var rows = references.length(index)
        var columns = queries.length(index)
        var moves = List[UInt8](capacity=rows + columns)
        traced[value](space, lane, rows, columns, low, high, moves)
        if left:
            # Traced over both sequences reversed, from the origin on: turned right to left.
            reverse_bytes(moves.unsafe_ptr(), len(moves))
        moves_out[unsafe_offset=index] = moves^
        costs_out[unsafe_offset=index] = Optional[Int](cost)
    else:
        costs_out[unsafe_offset=index] = None
    settled[unsafe_offset=index] = True


def lane_alignment_stage[
    T: Texts, value: DType
](
    pairs: Int,
    references: T,
    queries: T,
    costs: LaneCosts,
    reference_band: Band,
    max_cost: Int,
    left: Bool,
    workers: Int,
    budget: Int,
    costs_out: MutPointer[Optional[Int], _],
    moves_out: MutPointer[List[UInt8], _],
    settled: MutPointer[Bool, _],
):
    """`lane_alignments` in lanes of `value`, as `lane_stage` runs `lane_distances`'s, each band's flags
    kept and every pair it proves traced before the next group takes them; a band proves a pair only when
    no path off it costs as little (see `proof`)."""
    comptime WIDTH = lanes_of[value]()
    if pairs == 0:
        return
    # The lanes' diagonals are the query's position less the reference's; reversed, a band leaving the
    # pair whole leaves it whole again.
    var lane_band = Band(-reference_band.high, -reference_band.low)
    var band = Band() if left else lane_band
    var order = dealt[T, value](
        pairs, references, queries, costs, workers, settled, lane_band if left else Band(), True
    )
    var placed = len(order)
    var order_ptr = order.unsafe_ptr()
    var spaces = List[LaneSpace[value]](capacity=workers)
    for _ in range(workers):
        spaces.append(LaneSpace[value]())
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
        imm left,
        imm budget,
        imm groups,
        imm first_workers,
        imm placed,
        imm order_ptr,
        imm costs_out,
        imm settled,
        imm space_ptr,
        imm moves_out,
    }:
        """Takes groups until none is left, settling and tracing every pair its band proves."""
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
                var rows = 0
                for slot in range(group * WIDTH, min(placed, (group + 1) * WIDTH)):
                    var index = order_ptr[unsafe_offset=slot]
                    space.members.append(index)
                    var end = queries.length(index) - references.length(index)
                    low = min(low, min(0, end) - 1)
                    high = max(high, max(0, end) + 1)
                    rows = max(rows, references.length(index))
                low = max(low, band.low)
                high = min(high, band.high)
                if rows * (high - low + 1) * WIDTH > budget:
                    continue
                var found = band_costs[T, value, True](references, queries, space, low, high, costs, left)
                for lane in range(len(space.members)):
                    var index = space.members[lane]
                    var cost = Int(found[lane])
                    var verdict = proof[value, True](
                        index,
                        cost,
                        references.length(index),
                        queries.length(index),
                        low,
                        high,
                        band,
                        costs,
                        max_cost,
                        space.retries,
                    )
                    if verdict == PROVEN:
                        settle_traced[T, value](
                            space,
                            lane,
                            index,
                            cost,
                            low,
                            high,
                            references,
                            queries,
                            max_cost,
                            left,
                            costs_out,
                            moves_out,
                            settled,
                        )
                    elif verdict == REFUSED:
                        costs_out[unsafe_offset=index] = None
                        settled[unsafe_offset=index] = True

    parallelize(first_pass, first_workers, first_workers)

    # Second pass: every pair the first band did not prove, over every diagonal a path as cheap as its
    # first cost could visit, which proves itself and holds every optimal path.
    var retries = List[Int]()
    for bucket in range(BUCKETS):
        for worker in range(workers):
            retries.extend(Span(spaces[worker].retries[bucket]))
    var unproven = len(retries) // 4
    if unproven == 0:
        return
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
        imm left,
        imm budget,
        imm costs_out,
        imm settled,
        imm retry_groups,
        imm second_workers,
        imm unproven,
        imm retry_ptr,
        imm space_ptr,
        imm moves_out,
    }:
        """Takes groups of unproven pairs until none is left, settling and tracing each."""
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
                var rows = 0
                for slot in range(group * WIDTH, min(unproven, (group + 1) * WIDTH)):
                    var index = retry_ptr[unsafe_offset=4 * slot]
                    space.members.append(index)
                    low = min(low, retry_ptr[unsafe_offset=4 * slot + 2])
                    high = max(high, retry_ptr[unsafe_offset=4 * slot + 3])
                    rows = max(rows, references.length(index))
                if rows * (high - low + 1) * WIDTH > budget:
                    continue
                var found = band_costs[T, value, True](references, queries, space, low, high, costs, left)
                for lane in range(len(space.members)):
                    # The band holds the first one, so its cost is no dearer, and under a byte's 255.
                    settle_traced[T, value](
                        space,
                        lane,
                        space.members[lane],
                        Int(found[lane]),
                        low,
                        high,
                        references,
                        queries,
                        max_cost,
                        left,
                        costs_out,
                        moves_out,
                        settled,
                    )

    parallelize(second_pass, second_workers, second_workers)
