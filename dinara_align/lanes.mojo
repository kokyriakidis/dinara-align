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

Free ends take the lanes too, for a distance (see `LaneEnds`): each lane's first row and column cost
nothing over its free letters and a gap past them, and its cost is the least over every cell it may
end on, its last column on the rows its free reference letters reach and its last row over its free
query letters. A path off the band then runs out from the nearest start's diagonal and back to the
nearest end's. A read placed in a window may start and end on any row, and its band is the whole
matrix, swept 64 pairs at a time.

A local score takes the lanes in 16 bits, the whole matrix swept (see `lane_local_scores`): each lane
padded past its letters with bytes that match nothing, so no cell there beats its best. A table of more
than one mismatch score that fits a register, DNA's, is read by byte shuffle over the alphabet's codes
(see `CodeTexts`), each lane's pair at once.

An alignment takes the bytes alone (see `lane_alignments`). Its band must hold every optimal path, so a
path off it must cost more, not merely as much; each cell then keeps a flag naming the source the
wavefront's backtrace would take there, and a lane's path is traced from its corner through them, the
CIGAR `Ties` picks. The flags of a group's whole band would pass the caches on long pairs and fault in
fresh pages for every group, so the band's sweep keeps a row of every layer every few rows, and the
traceback sweeps each stretch between two again, from the last, keeping its flags alone.
"""

from std.bit import count_leading_zeros
from std.math import ceildiv
from std.memory import bitcast
from std.utils import IndexList
from std.sys import llvm_intrinsic, simd_width_of, size_of
from std.atomic import Atomic

from .ablation import ABLATE_CERTIFIED
from .cigar import reverse_bytes, reversed_text
from .substitutions import SHUFFLED_ENTRIES, looked_up
from .common import next_share, spread
from .gap_affine import ALIGNED, ENTRY_MASK, EndsFree, Penalties, FIRST_GAP, SECOND_GAP, gap_layer, layer_bit
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

    @staticmethod
    def of(items: List[String]) -> Self:
        """`items`' sequences, which must outlive the texts: a list read only through them is destroyed
        after its last use by name, so its owner keeps it to the end of the texts' use (`_ = items^`)."""
        return Self(items.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]())

    @always_inline
    def length(self, index: Int) -> Int:
        return self.items[unsafe_offset=index].byte_length()

    @always_inline
    def letters(self, index: Int) -> ImmPointer[UInt8, ImmUntrackedOrigin]:
        return self.items[unsafe_offset=index].unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]()


@always_inline
def lanes_of[value: DType]() -> Int:
    """Pairs a group: 64 bytes' worth of `value`s, one native register under AVX-512 and several where the
    CPU's are narrower, each step four registers' work with none waiting on another (see `local_width`): on
    the M2's 128 bits, 1 kbp pairs at 10% took 73 us a pair under affine costs where one register's 16 took
    132, and 32 at unit costs where they took 60."""
    return max(simd_width_of[value](), 64 // size_of[value]())


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


comptime TABLE_ENTRIES = SHUFFLED_ENTRIES
"""Entries a table holds at most for the lanes: one 16-byte register a byte shuffle reads, a four-letter
alphabet's, DNA's."""


@fieldwise_init
struct CodeTexts(Texts, TrivialRegisterPassable):
    """Sequences as an alphabet's codes, one after another in one buffer, sequence `i` from `starts[i]` to
    `starts[i + 1]`."""

    var codes: ImmPointer[UInt8, ImmUntrackedOrigin]
    var starts: ImmPointer[Int, ImmUntrackedOrigin]

    @always_inline
    def length(self, index: Int) -> Int:
        return self.starts[unsafe_offset=index + 1] - self.starts[unsafe_offset=index]

    @always_inline
    def letters(self, index: Int) -> ImmPointer[UInt8, ImmUntrackedOrigin]:
        return self.codes.unsafe_offset(self.starts[unsafe_offset=index])


@fieldwise_init
struct LaneEnds(ImplicitlyCopyable, TrivialRegisterPassable):
    """A pair's letters an alignment may leave unaligned for nothing, in the lanes' terms: reference
    rows and query columns at the start and at the end. All none is a global alignment."""

    var start_rows: Int
    var end_rows: Int
    var start_columns: Int
    var end_columns: Int

    @staticmethod
    def of(mode: Mode, rows: Int, columns: Int) -> Self:
        """`mode`'s free letters for a pair of `rows` reference letters and `columns` query letters (see
        `EndsFree.of`)."""
        var ends = EndsFree.of(mode, rows, columns)
        return Self(ends.first_begin, ends.first_end, ends.second_begin, ends.second_end)

    @always_inline
    def start_low(self) -> Int:
        """The lowest diagonal a path may start on: past the free reference letters, down the first column."""
        return -self.start_rows

    @always_inline
    def start_high(self) -> Int:
        """The highest: past the free query letters, along the first row."""
        return self.start_columns

    @always_inline
    def end_low(self, end: Int) -> Int:
        """The lowest diagonal a path may end on, the corner on diagonal `end`: short of the free query letters."""
        return end - self.end_columns

    @always_inline
    def end_high(self, end: Int) -> Int:
        """The highest: short of the free reference letters."""
        return end + self.end_rows


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
    var letters: Int
    """With a table, the alphabet's size, the sequences then its codes; zero to compare letters, a match
    free and a mismatch `mismatch`."""
    var table: SIMD[DType.uint8, TABLE_ENTRIES]
    """Pair `(a, b)`'s cost at `a letters + b`, `mismatch` the dearest."""

    @staticmethod
    def one_piece(
        mismatch: Int, deletion_opening: Int, deletion_extension: Int, insertion_opening: Int, insertion_extension: Int
    ) -> Self:
        """Costs of one gap piece."""
        return Self(
            mismatch,
            deletion_opening,
            deletion_extension,
            insertion_opening,
            insertion_extension,
            1,
            0,
            0,
            0,
            0,
            0,
            SIMD[DType.uint8, TABLE_ENTRIES](0),
        )

    def tabled(self, letters: Int, table: SIMD[DType.uint8, TABLE_ENTRIES]) -> Self:
        """These gap costs with each pair's cost read from `table`, over an alphabet of `letters` codes."""
        var dearest = 0
        for cell in range(letters * letters):
            dearest = max(dearest, Int(table[cell]))
        var costs = self
        costs.mismatch = dearest
        costs.letters = letters
        costs.table = table
        return costs

    @staticmethod
    def of(costs: Costs, mode: Mode, free_ends: Bool = False) -> Optional[Self]:
        """The lanes' costs for `costs`, or None where they do not serve: a match's reward, or free ends
        unless `free_ends`."""
        if mode.is_scored() or (not free_ends and not mode.is_global()):
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
            0,
            SIMD[DType.uint8, TABLE_ENTRIES](0),
        )

    @staticmethod
    def of_penalties(penalties: Penalties) -> Self:
        """The lanes' costs for a wavefront's one-piece costs, in its units."""
        return Self.one_piece(
            penalties.mismatch, penalties.opening, penalties.extension, penalties.opening, penalties.extension
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

    @always_inline
    def gap_lanes[
        value: DType, width: Int
    ](self, letters: SIMD[DType.int16, width], deletion: Bool) -> SIMD[value, width]:
        """`deleted` or `inserted` of each lane's `letters` at once, as lanes of `value` hold them (see `held`)."""
        var count = letters.cast[DType.int32]()
        var opening = self.deletion_opening if deletion else self.insertion_opening
        var extension = self.deletion_extension if deletion else self.insertion_extension
        var cost = SIMD[DType.int32, width](Int32(opening)) + count * Int32(extension)
        if self.pieces == 2:
            var opening2 = self.deletion_opening2 if deletion else self.insertion_opening2
            var extension2 = self.deletion_extension2 if deletion else self.insertion_extension2
            cost = min(cost, SIMD[DType.int32, width](Int32(opening2)) + count * Int32(extension2))
        cost = count.gt(0).select(cost, SIMD[DType.int32, width](0))
        return min(cost, SIMD[DType.int32, width](Int32(far_of[value]()))).cast[value]()

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

    def off_band(self, diagonal: Int, end: Int, ends: LaneEnds = LaneEnds(0, 0, 0, 0)) -> Int:
        """The least any path visiting `diagonal` costs, its end on diagonal `end` and its free letters
        `ends`: a gap out to the diagonal from the nearest start's side and one back to the nearest end's;
        nothing between the two."""
        var start_low = ends.start_low()
        var start_high = ends.start_high()
        var end_low = ends.end_low(end)
        var end_high = ends.end_high(end)
        if diagonal > max(start_high, end_high):
            return self.inserted(diagonal - start_high) + self.deleted(diagonal - end_high)
        if diagonal < min(start_low, end_low):
            return self.deleted(start_low - diagonal) + self.inserted(end_low - diagonal)
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
    """For an alignment, each cell's flag on the band over a stretch of rows, a row after another (see
    `band_flags`)."""
    var saved: List[SIMD[Self.value, lanes_of[Self.value]()]]
    """Every layer's row on each multiple of `every` rows, for `band_flags` to sweep again from."""
    var every: Int
    var rows: Int
    """The group's most rows and columns, and each lane's own."""
    var columns: Int
    var row_ends: SIMD[DType.int16, lanes_of[Self.value]()]
    var column_ends: SIMD[DType.int16, lanes_of[Self.value]()]
    var stop: Int
    """The cost past which a lane is done with: once every lane's row has reached it, or its last row, the
    sweep stops (see `swept`); a byte's 255, or a cap's next cost."""
    var stopped: Bool
    """Whether the last sweep stopped early so."""
    var free: Bool
    """Whether any lane has free letters at an end, which its own allowances below give."""
    var last_columns: List[Int]
    """The lanes' last columns, each once, and the lanes ending on each, for `ends_reached`."""
    var last_lanes: List[SIMD[DType.bool, lanes_of[Self.value]()]]
    var start_rows: SIMD[DType.int16, lanes_of[Self.value]()]
    var end_rows: SIMD[DType.int16, lanes_of[Self.value]()]
    var start_columns: SIMD[DType.int16, lanes_of[Self.value]()]
    var end_columns: SIMD[DType.int16, lanes_of[Self.value]()]
    var found_diagonals: SIMD[DType.int16, lanes_of[Self.value]()]
    """With free ends, the diagonal of the cell each lane's least cost so far ends on (see `ends_reached`)."""
    var lowest_end: Bool
    """Of the cells a lane may end on at its least cost, whether `found_diagonals` keeps the one on the
    lowest diagonal, or the highest."""
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
        self.saved = List[SIMD[Self.value, lanes_of[Self.value]()]]()
        self.every = 0
        self.rows = 0
        self.columns = 0
        self.row_ends = SIMD[DType.int16, lanes_of[Self.value]()](-1)
        self.column_ends = SIMD[DType.int16, lanes_of[Self.value]()](0)
        self.stop = Int.MAX
        self.stopped = False
        self.free = False
        self.last_columns = List[Int]()
        self.last_lanes = List[SIMD[DType.bool, lanes_of[Self.value]()]]()
        self.start_rows = SIMD[DType.int16, lanes_of[Self.value]()](0)
        self.end_rows = SIMD[DType.int16, lanes_of[Self.value]()](0)
        self.start_columns = SIMD[DType.int16, lanes_of[Self.value]()](0)
        self.end_columns = SIMD[DType.int16, lanes_of[Self.value]()](0)
        self.found_diagonals = SIMD[DType.int16, lanes_of[Self.value]()](0)
        self.lowest_end = True
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
    pad: Int = -1,
):
    """The members' texts, each back to front with `reverse`, position `p` of every lane at
    `target[p width:(p + 1) width]`: each text copied whole into `staging`, a lane's stretch padded to
    whole blocks, then a block of `BLOCK` positions at a time, every lane's letters of the block one load,
    the block turned in registers in two steps: eight lanes' bytes inside 64 at a time, then the 64-bit
    words of all of them. Past a lane's own letters each position holds `pad`, if one is given, else
    whatever was there."""
    comptime groups = width // 8
    comptime bytes = byte_order()
    comptime words = word_order[groups]()
    var stride = ceildiv(max(length, 1), BLOCK) * BLOCK
    staging.resize(unsafe_uninit_length=width * stride)
    var lanes = staging.unsafe_ptr()
    for lane in range(len(members)):
        var index = members[lane]
        var count = texts.length(index)
        if pad >= 0:
            for position in range(count, stride):
                lanes[unsafe_offset=lane * stride + position] = UInt8(pad)
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
    comptime if width < 8:
        # Too few lanes for a block's 64-bit words, which 32-bit lanes on NEON's 128 bits leave: a byte at a
        # time.
        for start in range(stride):
            comptime for lane in range(width):
                target[unsafe_offset=start * width + lane] = lanes[unsafe_offset=lane * stride + start]
    else:
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


comptime SECOND_PIECE_DELETION = gap_layer(1, True)
"""The layer of the second gap piece's deletions."""
comptime SECOND_PIECE_INSERTION = gap_layer(1, False)
"""The layer of the second gap piece's insertions."""
comptime SOURCE_MASK = ENTRY_MASK
"""A flag's low three bits: the source its alignment layer takes, as the wavefront's flags keep it."""


@always_inline
def extended_bit(layer: Int) -> UInt8:
    """A flag bit: gap layer `layer` extends, not opens, into the cell; the wavefront's same bit says the
    opposite (see `gap_affine.opened_bit`)."""
    return layer_bit(layer)


def band_costs[
    T: Texts, value: DType
](
    references: T,
    queries: T,
    mut space: LaneSpace[value],
    low: Int,
    high: Int,
    costs: LaneCosts,
    reverse: Bool = False,
    every: Int = 0,
    mode: Mode = Mode.GLOBAL,
) -> SIMD[value, lanes_of[value]()]:
    """The least cost of each of `space.members`'s pairs over the paths on diagonals `low ..= high`, the
    letters `mode` leaves free at either end unaligned for nothing, `far_of[value]()` or more for a pair
    none of whose paths stays on them, or, in a byte, whose cost reaches 255; a lane past the members
    holds nothing. References run down the rows and queries across, both back to front with `reverse`.
    With `every` rows, the rows on a multiple of it are kept in `space.saved`, for `band_flags` to sweep
    a stretch again from."""
    var count = len(space.members)
    var rows = 0
    var columns = 0
    space.row_ends = SIMD[DType.int16, lanes_of[value]()](-1)
    space.column_ends = SIMD[DType.int16, lanes_of[value]()](0)
    space.free = False
    for lane in range(count):
        var index = space.members[lane]
        var reference = references.length(index)
        var query = queries.length(index)
        rows = max(rows, reference)
        columns = max(columns, query)
        space.row_ends[lane] = Int16(reference)
        space.column_ends[lane] = Int16(query)
        var ends = LaneEnds.of(mode, reference, query)
        space.start_rows[lane] = Int16(ends.start_rows)
        space.end_rows[lane] = Int16(ends.end_rows)
        space.start_columns[lane] = Int16(ends.start_columns)
        space.end_columns[lane] = Int16(ends.end_columns)
        space.free = space.free or ends.start_rows + ends.end_rows + ends.start_columns + ends.end_columns > 0
    if space.free:
        space.last_columns.clear()
        space.last_lanes.clear()
        for lane in range(count):
            var column = Int(space.column_ends[lane])
            if column not in space.last_columns:
                space.last_columns.append(column)
                space.last_lanes.append(space.column_ends.eq(Int16(column)))
    space.rows = rows
    space.columns = columns
    space.every = every
    # Each position's letters side by side. Past a lane's own letters the rows and columns hold whatever
    # was there before: no cell past a pair's sequences feeds its corner.
    comptime WIDTH = lanes_of[value]()
    space.row_letters.resize(unsafe_uninit_length=ceildiv(max(rows, 1), BLOCK) * BLOCK * WIDTH)
    space.column_letters.resize(unsafe_uninit_length=ceildiv(max(columns, 1), BLOCK) * BLOCK * WIDTH)
    side_by_side[T, WIDTH](references, space.members, rows, space.staging, space.row_letters.unsafe_ptr(), reverse)
    side_by_side[T, WIDTH](queries, space.members, columns, space.staging, space.column_letters.unsafe_ptr(), reverse)
    if costs.letters > 0:
        if costs.pieces == 2:
            return pieced_band_costs[value, 2, True](space, low, high, costs, every)
        return pieced_band_costs[value, 1, True](space, low, high, costs, every)
    if costs.pieces == 2:
        return pieced_band_costs[value, 2](space, low, high, costs, every)
    return pieced_band_costs[value, 1](space, low, high, costs, every)


def pieced_band_costs[
    value: DType, pieces: Int, tabled: Bool = False
](mut space: LaneSpace[value], low: Int, high: Int, costs: LaneCosts, every: Int) -> SIMD[value, lanes_of[value]()]:
    """`band_costs` at `pieces` gap pieces, each pair's cost from `costs.table` with `tabled`, the letters
    already side by side."""
    comptime WIDTH = lanes_of[value]()
    comptime Lanes = SIMD[value, WIDTH]
    comptime FAR = far_of[value]()
    var columns = space.columns
    space.stopped = False
    space.scores.resize(columns + 2, Lanes(FAR))
    space.deletions.resize(columns + 2, Lanes(FAR))
    comptime if pieces == 2:
        space.deletions2.resize(columns + 2, Lanes(FAR))
    first_row[value, pieces](space, low, high, costs)
    var found = Lanes(FAR)
    if space.free:
        space.found_diagonals = SIMD[DType.int16, lanes_of[value]()](0)
        ends_reached[value](space, low, high, 0, found)
        swept[value, pieces, False, True, tabled](space, low, high, costs, 0, space.rows, found)
        return found
    for lane in range(len(space.members)):
        var column = Int(space.column_ends[lane])
        if space.row_ends[lane] == 0 and column >= low and column <= high:
            found[lane] = space.scores[column][lane]
    if every <= 0:
        swept[value, pieces, False, False, tabled](space, low, high, costs, 0, space.rows, found)
        return found
    space.saved.clear()
    saved_row[value, pieces](space, low, high, 0)
    var row = 0
    while row < space.rows:
        var next = min(row + every, space.rows)
        swept[value, pieces, False, False, tabled](space, low, high, costs, row, next, found)
        if space.stopped:
            # Every lane done: its walk, if any, stays in the rows kept so far.
            break
        if next % every == 0 and next < space.rows:
            saved_row[value, pieces](space, low, high, next)
        row = next
    return found


@always_inline
def held[value: DType](cost: Int) -> Scalar[value]:
    """A cost as a lane holds it, 255 or more a byte's 255."""
    return Scalar[value](min(cost, far_of[value]()))


def first_row[value: DType, pieces: Int](mut space: LaneSpace[value], low: Int, high: Int, costs: LaneCosts):
    """Every layer's row as row zero leaves it inside the band, an insertion of every column before, and
    reached by nothing elsewhere."""
    comptime Lanes = SIMD[value, lanes_of[value]()]
    comptime FAR = far_of[value]()
    var scores = space.scores.unsafe_ptr()
    var deletions = space.deletions.unsafe_ptr()
    var deletions2 = space.deletions2.unsafe_ptr()
    for column in range(space.columns + 2):
        scores[unsafe_offset=column] = Lanes(FAR)
        deletions[unsafe_offset=column] = Lanes(FAR)
        comptime if pieces == 2:
            deletions2[unsafe_offset=column] = Lanes(FAR)
    if space.free:
        # Each lane's free query letters cost nothing; past them, a gap from the last.
        for column in range(max(low, 0), min(high, space.columns) + 1):
            scores[unsafe_offset=column] = costs.gap_lanes[value](Int16(column) - space.start_columns, False)
        return
    for column in range(max(low, 0), min(high, space.columns) + 1):
        scores[unsafe_offset=column] = Lanes(held[value](costs.inserted(column)))


@always_inline
def saved_span(row: Int, low: Int, high: Int, columns: Int) -> Tuple[Int, Int]:
    """The columns of row `row` a sweep from it reads again: its band and the column before."""
    return (max(max(row + low, 0) - 1, 0), min(row + high, columns))


def saved_row[value: DType, pieces: Int](mut space: LaneSpace[value], low: Int, high: Int, row: Int):
    """Keeps row `row` of every layer in `space.saved`, as a sweep from it reads it, `high - low + 2`
    entries a layer."""
    comptime Lanes = SIMD[value, lanes_of[value]()]
    var stride = high - low + 2
    var bounds = saved_span(row, low, high, space.columns)
    var start = len(space.saved)
    space.saved.resize(start + (1 + pieces) * stride, Lanes(0))
    var target = space.saved.unsafe_ptr().unsafe_offset(start)
    for column in range(bounds[0], bounds[1] + 1):
        target[unsafe_offset=column - bounds[0]] = space.scores[column]
        target[unsafe_offset=stride + column - bounds[0]] = space.deletions[column]
        comptime if pieces == 2:
            target[unsafe_offset=2 * stride + column - bounds[0]] = space.deletions2[column]


def restored_row[
    value: DType, pieces: Int
](mut space: LaneSpace[value], low: Int, high: Int, row: Int, slot: Int, rows: Int):
    """Every layer's row as the sweep left it after row `row`, from `space.saved`'s slot `slot`: what a
    sweep of the next `rows` rows from it reads, the band and the column before, and the columns past the
    band they reach, one a row, reached by nothing, as no later row had written there yet."""
    comptime Lanes = SIMD[value, lanes_of[value]()]
    comptime FAR = far_of[value]()
    var stride = high - low + 2
    var bounds = saved_span(row, low, high, space.columns)
    var source = space.saved.unsafe_ptr().unsafe_offset(slot * (1 + pieces) * stride)
    for column in range(bounds[0], bounds[1] + 1):
        space.scores[column] = source[unsafe_offset=column - bounds[0]]
        space.deletions[column] = source[unsafe_offset=stride + column - bounds[0]]
        comptime if pieces == 2:
            space.deletions2[column] = source[unsafe_offset=2 * stride + column - bounds[0]]
    for column in range(bounds[1] + 1, min(bounds[1] + rows + 1, space.columns + 2)):
        space.scores[column] = Lanes(FAR)
        space.deletions[column] = Lanes(FAR)
        comptime if pieces == 2:
            space.deletions2[column] = Lanes(FAR)


def band_flags[value: DType](mut space: LaneSpace[value], low: Int, high: Int, costs: LaneCosts, start: Int, end: Int):
    """Sweeps rows `start + 1 ..= end` again, from the row `band_costs` kept at `start`, keeping each cell's
    flag in `space.flags`, a row after another from `start + 1` (see `swept`)."""
    var found = SIMD[value, lanes_of[value]()](0)
    if costs.pieces == 2:
        if start == 0:
            first_row[value, 2](space, low, high, costs)
        else:
            restored_row[value, 2](space, low, high, start, start // space.every, end - start)
        if costs.letters > 0:
            swept[value, 2, True, False, True](space, low, high, costs, start, end, found)
        else:
            swept[value, 2, True](space, low, high, costs, start, end, found)
    else:
        if start == 0:
            first_row[value, 1](space, low, high, costs)
        else:
            restored_row[value, 1](space, low, high, start, start // space.every, end - start)
        if costs.letters > 0:
            swept[value, 1, True, False, True](space, low, high, costs, start, end, found)
        else:
            swept[value, 1, True](space, low, high, costs, start, end, found)


def ends_reached[
    value: DType
](mut space: LaneSpace[value], low: Int, high: Int, row: Int, mut found: SIMD[value, lanes_of[value]()]):
    """Lowers each lane's `found` to its cells on row `row` an alignment may end on: its last column, on
    a row its free reference letters reach, and on its last row, every column its free query letters do;
    inside the band. The last columns a column at a time, every lane ending on it at once. Each lane's
    `space.found_diagonals` follows the cell, of those at its least cost the one on the lowest diagonal,
    or with `space.lowest_end` off the highest."""
    comptime WIDTH = lanes_of[value]()
    comptime Lanes = SIMD[value, WIDTH]
    comptime Diagonals = SIMD[DType.int16, WIDTH]
    var last_rows = space.row_ends
    var reached = (last_rows - space.end_rows).le(Int16(row)) & last_rows.ge(Int16(row))
    if not reached.reduce_or():
        return
    var first = max(row + low, 0)
    var last = min(row + high, space.columns)
    var lowest = space.lowest_end
    for slot in range(len(space.last_columns)):
        var column = space.last_columns[slot]
        if column >= first and column <= last:
            var cell = space.scores[column]
            var diagonal = Diagonals(Int16(column - row))
            var nearer = diagonal.lt(space.found_diagonals) if lowest else diagonal.gt(space.found_diagonals)
            var taken = space.last_lanes[slot] & reached & (cell.lt(found) | (cell.eq(found) & nearer))
            found = taken.select(cell, found)
            space.found_diagonals = taken.select(diagonal, space.found_diagonals)
    # A lane's free query letters on its last row, once.
    var finishing = last_rows.eq(Int16(row)) & space.end_columns.gt(0)
    if finishing.reduce_or():
        for lane in range(len(space.members)):
            if finishing[lane]:
                var column = Int(space.column_ends[lane])
                var best = found[lane]
                var best_diagonal = Int(space.found_diagonals[lane])
                for end in range(max(column - Int(space.end_columns[lane]), first), min(column, last) + 1):
                    var cell = space.scores[end][lane]
                    var diagonal = end - row
                    var nearer = diagonal < best_diagonal if lowest else diagonal > best_diagonal
                    if cell < best or (cell == best and nearer):
                        best = cell
                        best_diagonal = diagonal
                found[lane] = best
                space.found_diagonals[lane] = Int16(best_diagonal)


def swept[
    value: DType, pieces: Int, record: Bool, free: Bool = False, tabled: Bool = False
](
    mut space: LaneSpace[value],
    low: Int,
    high: Int,
    costs: LaneCosts,
    start: Int,
    end: Int,
    mut found: SIMD[value, lanes_of[value]()],
):
    """Gotoh's recurrence over rows `start + 1 ..= end` of the band, every lane at once, from the rows
    `space` holds for row `start`; a lane whose reference ends on one of them reads its corner into
    `found`, or with `free`, every cell it may end on (see `ends_reached`).

    With `record`, each cell keeps a flag in `space.flags`, for `traced`: the source its alignment layer
    takes, `ALIGNED` for the diagonal, and a bit a gap layer, `extended_bit`, where the layer extends. The
    source is the one the wavefront's backtrace takes (see `Ties`): of the moves into the cell as cheap as
    it, a substitution, then a letter of the reference alone, then one of the query, the second gap piece
    before the first, and only then a match; a gap layer's extension before its opening."""
    comptime WIDTH = lanes_of[value]()
    comptime Lanes = SIMD[value, WIDTH]
    comptime FAR = far_of[value]()
    comptime Flags = SIMD[DType.uint8, WIDTH]
    var columns = space.columns
    var span = high - low + 1
    comptime if record:
        space.flags.resize(unsafe_uninit_length=(end - start) * span)
    var flags = space.flags.unsafe_ptr()
    var scores = space.scores.unsafe_ptr()
    var deletions = space.deletions.unsafe_ptr()
    var deletions2 = space.deletions2.unsafe_ptr()
    var row_ends = space.row_ends
    var column_ends = space.column_ends
    var mismatch = Lanes(held[value](costs.mismatch))
    var delete_open = Lanes(held[value](costs.deletion_opening + costs.deletion_extension))
    var delete_extend = Lanes(held[value](costs.deletion_extension))
    var insert_open = Lanes(held[value](costs.insertion_opening + costs.insertion_extension))
    var insert_extend = Lanes(held[value](costs.insertion_extension))
    var delete_open2 = Lanes(held[value](costs.deletion_opening2 + costs.deletion_extension2))
    var delete_extend2 = Lanes(held[value](costs.deletion_extension2))
    var insert_open2 = Lanes(held[value](costs.insertion_opening2 + costs.insertion_extension2))
    var insert_extend2 = Lanes(held[value](costs.insertion_extension2))
    var row_source = space.row_letters.unsafe_ptr()
    var column_source = space.column_letters.unsafe_ptr()
    var stopping = space.stop <= FAR
    var stop = Lanes(Scalar[value](min(space.stop, FAR)))
    for row in range(start + 1, end + 1):
        var first = max(row + low, 0)
        var last = min(row + high, columns)
        if first > last:
            continue
        var letter = row_source.unsafe_offset((row - 1) * WIDTH).unsafe_load[width=WIDTH]()
        # With a table, the row's letters as the first half of each pair's entry.
        var row_entries = letter * UInt8(costs.letters)
        var left = Lanes(FAR)
        var diagonal: Lanes
        if first == 0:
            # The left edge: a deletion of every row so far, past the lane's free reference letters.
            diagonal = scores[unsafe_offset=0]
            comptime if free:
                left = costs.gap_lanes[value](Int16(row) - space.start_rows, True)
            else:
                left = Lanes(held[value](costs.deleted(row)))
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
        # The row's least cost so far, the left edge's among them.
        var least = left
        # This row's flags, indexed by column.
        var row_flags = flags.unsafe_offset((row - start - 1) * span - row - low)
        for column in range(first, last + 1):
            var other = column_source.unsafe_offset((column - 1) * WIDTH).unsafe_load[width=WIDTH]()
            var above = scores[unsafe_offset=column]
            var equal = letter.eq(other)
            var substituted: Lanes
            comptime if tabled:
                substituted = added(diagonal, looked_up[value, WIDTH, False](costs.table, row_entries + other))
            else:
                substituted = added(diagonal, equal.select(Lanes(0), mismatch))
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
            comptime if not record:
                least = min(least, score)
            diagonal = above
            left = score
        # A lane whose reference ends on this row reads its corner, if the band holds it; with free ends,
        # every cell of the row it may end on.
        comptime if free:
            ends_reached[value](space, low, high, row, found)
        elif not record:
            if row_ends.eq(Int16(row)).reduce_or():
                for lane in range(len(space.members)):
                    if Int(row_ends[lane]) == row:
                        var column = Int(column_ends[lane])
                        if column >= max(row + low, 0) and column <= last:
                            found[lane] = scores[unsafe_offset=column][lane]
        # No row costs less than the one before, every cell coming from it or from its left at no saving:
        # once every lane's row has reached `stop`, or the lane has passed its last row, no lane can end
        # under it, and the sweep has nothing left to find.
        comptime if not record:
            if stopping and (least.ge(stop) | row_ends.le(Int16(row))).reduce_and():
                space.stopped = True
                return


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
    mode: Mode = Mode.GLOBAL,
) -> Int:
    """Every pair `settled` does not already mark that the lanes hold, its distance within `reference_band`
    with the letters `mode` leaves free at either end into `costs_out`, None past `max_cost` or with no
    path inside the band, and `settled` set; the others left for the caller, a pair with an empty side
    among them when an end is free. The pairs settled here: first in bytes, where every step's cost fits
    one, then in 16 bits."""
    var before = 0
    for index in range(pairs):
        before += Int(settled[unsafe_offset=index])
    var unlocated = List[Int](length=1, fill=0)
    if costs.fits_bytes():
        lane_stage[T, DType.uint8](
            pairs,
            references,
            queries,
            costs,
            reference_band,
            max_cost,
            workers,
            costs_out,
            settled,
            mode,
            unlocated.unsafe_ptr(),
        )
    lane_stage[T, DType.int16](
        pairs,
        references,
        queries,
        costs,
        reference_band,
        max_cost,
        workers,
        costs_out,
        settled,
        mode,
        unlocated.unsafe_ptr(),
    )
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
            # A byte holds any pair its lanes' 16-bit coordinates reach, its saturation telling; 16 bits only
            # a pair whose costs fit, which `fits` keeps far shorter.
            var held = (value == DType.uint8 and max(rows, columns) <= Int(Int16.MAX)) or costs.fits(rows, columns)
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

    spread(measure, stretches, stretches)
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

    spread(count, stretches, stretches)
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

    spread(place, stretches, stretches)
    return order^


def widest_last[value: DType](spaces: List[LaneSpace[value]], workers: Int) -> List[Int]:
    """Every worker's unproven pairs, four numbers a pair (see `proof`), the narrowest band first, so a
    second-pass group's lanes need bands alike and its band, the widest of theirs, wastes little: by the
    bands' bit lengths alone, a group swept up to twice the band most of its lanes needed."""
    var gathered = List[Int]()
    for bucket in range(BUCKETS):
        for worker in range(workers):
            gathered.extend(Span(spaces[worker].retries[bucket]))
    var keys = List[Int](capacity=len(gathered) // 4)
    for slot in range(len(gathered) // 4):
        keys.append(((gathered[4 * slot + 3] - gathered[4 * slot + 2]) << 32) | slot)
    sort(keys)
    var retries = List[Int](capacity=len(gathered))
    for key in keys:
        var slot = key & ((1 << 32) - 1)
        for field in range(4):
            retries.append(gathered[4 * slot + field])
    return retries^


def lane_stage[
    T: Texts, value: DType, located: Bool = False
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
    mode: Mode,
    ends_out: MutPointer[Int, _],
    lowest_end: Bool = True,
):
    """`lane_distances` in lanes of `value`: every pair `settled` does not mark and the lanes hold settled,
    a pair whose byte saturates left. With `located`, as `lane_ends` needs, each pair's band proves its
    cost only when no path off it costs as little, so every optimal path lies inside, and the diagonal of
    the cell its alignment ends on goes to `ends_out`: of those at its least cost, the lowest, or with
    `lowest_end` off the highest.

    The band is the library's, its diagonals the reference's position less the query's; the lanes run
    the reference down the rows, their diagonals the query's position less the reference's, so they
    take it mirrored."""
    comptime WIDTH = lanes_of[value]()
    var band = Band(-reference_band.high, -reference_band.low)
    if pairs == 0:
        return
    var free = not mode.is_global()
    var order = dealt[T, value](pairs, references, queries, costs, workers, settled, Band(), free)
    var placed = len(order)
    var order_ptr = order.unsafe_ptr()

    # First pass: each group over the band between its starts' and its ends' diagonals and one diagonal
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
        imm mode,
        imm groups,
        imm first_workers,
        imm placed,
        imm order_ptr,
        imm costs_out,
        imm settled,
        imm space_ptr,
        imm ends_out,
        imm lowest_end,
    }:
        """Takes groups until none is left, settling every pair its band proves."""
        ref space = space_ptr[unsafe_offset=worker]
        space.lowest_end = lowest_end
        # A byte is done at its 255; 16 bits at the cap's next cost, if there is a cap they hold.
        space.stop = far_of[value]() if value == DType.uint8 or max_cost >= far_of[value]() - 1 else max_cost + 1
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
                    var rows = references.length(index)
                    var columns = queries.length(index)
                    var end = columns - rows
                    var ends = LaneEnds.of(mode, rows, columns)
                    low = min(low, min(ends.start_low(), ends.end_low(end)) - 1)
                    high = max(high, max(ends.start_high(), ends.end_high(end)) + 1)
                comptime if ABLATE_CERTIFIED:
                    # The ablation: every lane's whole matrix, which proves itself.
                    for lane in range(len(space.members)):
                        low = min(low, -references.length(space.members[lane]))
                        high = max(high, queries.length(space.members[lane]))
                low = max(low, band.low)
                high = min(high, band.high)
                # A band holding none of the group's diagonals leaves its pairs to their own searches.
                if low > high:
                    continue
                var found = band_costs(references, queries, space, low, high, costs, False, 0, mode)
                for lane in range(len(space.members)):
                    var index = space.members[lane]
                    var cost = Int(found[lane])
                    var rows = references.length(index)
                    var columns = queries.length(index)
                    var verdict = proof[value, located](
                        index,
                        cost,
                        rows,
                        columns,
                        low,
                        high,
                        band,
                        costs,
                        max_cost,
                        space.retries,
                        LaneEnds.of(mode, rows, columns),
                    )
                    if verdict == PROVEN or verdict == REFUSED:
                        var kept = verdict == PROVEN and cost <= max_cost
                        costs_out[unsafe_offset=index] = Optional[Int](cost) if kept else None
                        settled[unsafe_offset=index] = True
                        comptime if located:
                            ends_out[unsafe_offset=index] = Int(space.found_diagonals[lane])

    spread(first_pass, first_workers, first_workers)

    # Second pass: every pair the first band did not prove, over every diagonal a cheaper path than its
    # first cost could visit, which proves itself; gathered from the workers' buckets by the width of
    # that band, so a group's lanes need bands alike.
    var retries = widest_last(spaces, workers)
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
        imm band,
        imm references,
        imm queries,
        imm costs,
        imm max_cost,
        imm mode,
        imm retry_groups,
        imm second_workers,
        imm unproven,
        imm retry_ptr,
        imm costs_out,
        imm settled,
        imm space_ptr,
        imm ends_out,
        imm lowest_end,
    }:
        """Takes groups of unproven pairs until none is left, settling each."""
        ref space = space_ptr[unsafe_offset=worker]
        space.lowest_end = lowest_end
        # A byte is done at its 255; 16 bits at the cap's next cost, if there is a cap they hold.
        space.stop = far_of[value]() if value == DType.uint8 or max_cost >= far_of[value]() - 1 else max_cost + 1
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
                # Within the band as the first pass is: the diagonal the group started from may lie outside it.
                low = max(low, band.low)
                high = min(high, band.high)
                var found = band_costs(references, queries, space, low, high, costs, False, 0, mode)
                for lane in range(len(space.members)):
                    var slot = group * WIDTH + lane
                    # A path as cheap as the first band's cost the first band found already; one cheaper is
                    # under 255, so no byte saturates on it. Located, the band holds the first one, every
                    # optimal path with it.
                    var cost = min(Int(found[lane]), retry_ptr[unsafe_offset=4 * slot + 1])
                    costs_out[unsafe_offset=space.members[lane]] = Optional[Int](cost) if cost <= max_cost else None
                    settled[unsafe_offset=space.members[lane]] = True
                    comptime if located:
                        ends_out[unsafe_offset=space.members[lane]] = Int(space.found_diagonals[lane])

    spread(second_pass, second_workers, second_workers)


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
    ends: LaneEnds = LaneEnds(0, 0, 0, 0),
) -> Int:
    """Whether the cost `found` on diagonals `low ..= high` settles pair `index`, of `rows` reference letters
    and `columns` query letters and the free letters `ends`: proven when no path off them could beat it,
    or one past the cap could not come under it; else the pair is filed, its cost and the band a cheaper
    path needs, by that band's width. With `strict`, as an alignment needs, no path off them may even
    match it: the tie rule's path, an optimal one, then lies inside. A path cannot leave the matrix, nor
    the band. A cost at FAR may be any from there up, or none inside the band: such a pair is left, by a
    byte's lanes for 16 bits, by 16 bits' for the pair's own search. 16 bits reach FAR only where a band
    rules out the path `LaneCosts.fits` bounds a pair's cost by."""
    if found >= far_of[value]():
        return UNHELD
    var end = columns - rows
    if max(ends.start_low(), band.low) > min(ends.start_high(), band.high) or max(ends.end_low(end), band.low) > min(
        ends.end_high(end), band.high
    ):
        # Every start or every end lies outside the band: no alignment inside it.
        return REFUSED
    var top = min(band.high, columns)
    var bottom = max(band.low, -rows)
    var beyond_high = costs.off_band(high + 1, end, ends) if high < top else Int.MAX
    var beyond_low = costs.off_band(low - 1, end, ends) if low > bottom else Int.MAX
    var off = min(beyond_high, beyond_low)
    var beaten = found < off if strict else found <= off
    if beaten or off > max_cost:
        # The least cost, or, its band's past the cap as is every path off it, past the cap.
        return PROVEN
    # Every diagonal a path cheaper than the target could visit, within the band and the matrix: a path
    # dearer than the cap need not be found, only shown past it.
    var target = min(found, max_cost) + 1 if strict else (found if found <= max_cost else max_cost + 1)
    var wide_high = high
    while wide_high < top and costs.off_band(wide_high + 1, end, ends) < target:
        wide_high += 1
    var wide_low = low
    while wide_low > bottom and costs.off_band(wide_low - 1, end, ends) < target:
        wide_low -= 1
    var width_bits = 64 - Int(count_leading_zeros(UInt64(wide_high - wide_low)))
    ref bucket = retries[min(width_bits, BUCKETS - 1)]
    bucket.append(index)
    bucket.append(found)
    bucket.append(wide_low)
    bucket.append(wide_high)
    return FILED


@fieldwise_init
struct Walk(Movable):
    """One lane's traceback as it goes up the band a stretch of rows at a time (see `walked`): the cell it
    stands at, the layer it is in, and the moves so far, right to left from its corner."""

    var lane: Int
    var index: Int
    """The pair's place in the batch."""
    var cost: Int
    var row: Int
    var column: Int
    var layer: Int
    var moves: List[UInt8]


def walked[value: DType](space: LaneSpace[value], mut walk: Walk, start: Int, low: Int, high: Int):
    """Takes `walk` up through the flags `band_flags` kept over rows `start + 1` on, until it stands on row
    `start`, or ends on the first row or column: at each cell the source its flag names, a gap layer's run
    until it opens. A path on the first row or column has one way left, a gap along it."""
    comptime WIDTH = lanes_of[value]()
    var span = high - low + 1
    var flags = space.flags.unsafe_ptr().unsafe_bitcast[UInt8]().unsafe_offset(walk.lane)
    var row = walk.row
    var column = walk.column
    var layer = walk.layer
    var out = walk.moves.unsafe_ptr()
    var written = len(walk.moves)
    # The cell's flag, `WIDTH` bytes apart: a row up is a span less a diagonal, a column left a diagonal.
    var cell = ((row - start - 1) * span + column - row - low) * WIDTH
    while row > start and column > 0:
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
    if row == 0 or column == 0:
        for _ in range(row):
            out[unsafe_offset=written] = UInt8(FIRST_GAP)
            written += 1
        for _ in range(column):
            out[unsafe_offset=written] = UInt8(SECOND_GAP)
            written += 1
        row = 0
        column = 0
    walk.moves.resize(unsafe_uninit_length=written)
    walk.row = row
    walk.column = column
    walk.layer = layer


@always_inline
def stretch_rows[value: DType](rows: Int, pieces: Int) -> Int:
    """The rows between the ones a traceback keeps (see `band_costs`' `every`): the square root of a kept
    row's bytes a cell over its flag's, which keeps the least memory, the kept rows and one stretch's
    flags weighing the same."""
    var kept = (1 + pieces) * size_of[value]()
    var every = 1
    while every * every < kept * rows:
        every += 1
    return every


def walk_of[T: Texts](lane: Int, index: Int, cost: Int, references: T, queries: T) -> Walk:
    """A walk for pair `index` in lane `lane`, at its corner, with room for its moves."""
    var rows = references.length(index)
    var columns = queries.length(index)
    return Walk(lane, index, cost, rows, columns, ALIGNED, List[UInt8](capacity=rows + columns))


def traced_bytes[value: DType](rows: Int, span: Int, every: Int, pieces: Int) -> Int:
    """The memory a group's traceback takes over `rows` rows of a band `span` diagonals wide, a row kept
    every `every`: those rows of every layer, and one stretch's flags."""
    comptime WIDTH = lanes_of[value]()
    var kept = (rows // every + 1) * (span + 1) * (1 + pieces) * size_of[value]() * WIDTH
    return kept + every * span * WIDTH


def traced_group[
    T: Texts, value: DType
](
    mut space: LaneSpace[value],
    mut walks: List[Walk],
    low: Int,
    high: Int,
    costs: LaneCosts,
    references: T,
    queries: T,
    left: Bool,
    costs_out: MutPointer[Optional[Int], _],
    moves_out: MutPointer[List[UInt8], _],
    settled: MutPointer[Bool, _],
):
    """Traces each of `walks`, the lanes of the group `band_costs` just swept with its kept rows, a stretch
    of rows at a time from the last: each stretch's flags swept again from the row kept above it, and
    every walk standing in it taken up through it. Each pair is then settled with its cost and its path's
    moves, right to left; with `left`, traced over both sequences reversed, turned around."""
    var every = space.every
    var start = ((space.rows - 1) // every) * every
    while start >= 0:
        var needed = False
        for slot in range(len(walks)):
            needed = needed or walks[slot].row > start
        if needed:
            band_flags[value](space, low, high, costs, start, min(start + every, space.rows))
            for slot in range(len(walks)):
                if walks[slot].row > start:
                    walked[value](space, walks[slot], start, low, high)
        start -= every
    for slot in range(len(walks)):
        ref walk = walks[slot]
        if left:
            # Traced over both sequences reversed, from the origin on: turned right to left.
            reverse_bytes(walk.moves.unsafe_ptr(), len(walk.moves))
        var moves = List[UInt8]()
        swap(moves, walk.moves)
        moves_out[unsafe_offset=walk.index] = moves^
        costs_out[unsafe_offset=walk.index] = Optional[Int](walk.cost)
        settled[unsafe_offset=walk.index] = True


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
    """`lane_distances` for alignments: every pair `settled` does not mark that a byte's lanes hold, its
    global cost within `reference_band` into `costs_out` and its path's moves, right to left as
    `gap_affine.solve` appends them, into `moves_out`, None past `max_cost` or with no path inside the
    band, and `settled` set. The path is the one `Ties` picks, `Ties.LEFT` with `left`: the rule run over
    both sequences reversed, which takes only pairs the band leaves whole. A group whose traceback would
    pass `budget` bytes is left to the caller, as is a pair with an empty side.

    Bytes alone, 64 pairs a register under AVX-512: a pair whose cost a byte cannot hold is left too. In
    16 bits, half as many lanes sweep the wide band such a pair's proof needs, which on the Skylake-X
    took as long as the searches on 1 kbp reads at 10%, and twice as long on NEON's eight lanes."""
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
    var after = 0
    for index in range(pairs):
        after += Int(settled[unsafe_offset=index])
    return after - before


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
        # A byte is done at its 255; 16 bits at the cap's next cost, if there is a cap they hold.
        space.stop = far_of[value]() if value == DType.uint8 or max_cost >= far_of[value]() - 1 else max_cost + 1
        var walks = List[Walk]()
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
                # A band holding none of the group's diagonals leaves its pairs to their own searches, which
                # find no alignment inside it.
                if low > high:
                    continue
                var every = stretch_rows[value](rows, costs.pieces)
                if traced_bytes[value](rows, high - low + 1, every, costs.pieces) > budget:
                    continue
                var found = band_costs[T, value](references, queries, space, low, high, costs, left, every)
                walks.clear()
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
                    if verdict == PROVEN and cost <= max_cost:
                        walks.append(walk_of(lane, index, cost, references, queries))
                    elif verdict == PROVEN or verdict == REFUSED:
                        costs_out[unsafe_offset=index] = None
                        settled[unsafe_offset=index] = True
                if len(walks) > 0:
                    traced_group[T, value](
                        space, walks, low, high, costs, references, queries, left, costs_out, moves_out, settled
                    )

    spread(first_pass, first_workers, first_workers)

    # Second pass: every pair the first band did not prove, over every diagonal a path as cheap as its
    # first cost could visit, which proves itself and holds every optimal path.
    var retries = widest_last(spaces, workers)
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
        # A byte is done at its 255; 16 bits at the cap's next cost, if there is a cap they hold.
        space.stop = far_of[value]() if value == DType.uint8 or max_cost >= far_of[value]() - 1 else max_cost + 1
        var walks = List[Walk]()
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
                var every = stretch_rows[value](rows, costs.pieces)
                if traced_bytes[value](rows, high - low + 1, every, costs.pieces) > budget:
                    continue
                var found = band_costs[T, value](references, queries, space, low, high, costs, left, every)
                walks.clear()
                for lane in range(len(space.members)):
                    # The band holds the first one, so its cost is no dearer, and under a byte's 255.
                    var index = space.members[lane]
                    var cost = Int(found[lane])
                    if cost <= max_cost:
                        walks.append(walk_of(lane, index, cost, references, queries))
                    else:
                        costs_out[unsafe_offset=index] = None
                        settled[unsafe_offset=index] = True
                if len(walks) > 0:
                    traced_group[T, value](
                        space, walks, low, high, costs, references, queries, left, costs_out, moves_out, settled
                    )

    spread(second_pass, second_workers, second_workers)


def lane_ends[
    T: Texts
](
    pairs: Int,
    references: T,
    queries: T,
    costs: LaneCosts,
    max_cost: Int,
    workers: Int,
    costs_out: MutPointer[Optional[Int], _],
    ends_out: MutPointer[Int, _],
    settled: MutPointer[Bool, _],
    mode: Mode,
    lowest_end: Bool,
):
    """`lane_distances` with `mode`'s free ends and no band, each pair's cost proven only where every optimal
    path lies in its band, and the diagonal of the cell its alignment ends on into `ends_out`: of those at
    its least cost, the lowest, or with `lowest_end` off the highest. A pair with an empty side is left."""
    if costs.fits_bytes():
        lane_stage[T, DType.uint8, True](
            pairs,
            references,
            queries,
            costs,
            Band(),
            max_cost,
            workers,
            costs_out,
            settled,
            mode,
            ends_out,
            lowest_end,
        )
    lane_stage[T, DType.int16, True](
        pairs, references, queries, costs, Band(), max_cost, workers, costs_out, settled, mode, ends_out, lowest_end
    )


@fieldwise_init
struct FramedSpans(Movable):
    """What `framed_spans` finds in the rule's frame: each pair's cost, its span, rows then columns, and
    whether the lanes found the span whole."""

    var found: List[Optional[Int]]
    var spans: List[Int]
    var spanning: List[Bool]


def framed_spans[
    F: Texts
](
    pairs: Int,
    firsts: F,
    seconds: F,
    costs: LaneCosts,
    frame_mode: Mode,
    max_cost: Int,
    workers: Int,
    costs_out: MutPointer[Optional[Int], _],
    settled: MutPointer[Bool, _],
) -> FramedSpans:
    """Each pair's span by the rule in its frame, `firsts` and `seconds` as the frame reads them (see
    `lane_free_alignments`): the end by `lane_ends` from the starts, then the start back from it. A pair
    past `max_cost` is settled here, None."""
    # The end: its cost and its diagonal, the lowest at that cost.
    var found = List[Optional[Int]](length=pairs, fill=None)
    var ends = List[Int](length=pairs, fill=0)
    var reached = List[Bool](length=pairs, fill=False)
    lane_ends(
        pairs,
        firsts,
        seconds,
        costs,
        max_cost,
        workers,
        found.unsafe_ptr(),
        ends.unsafe_ptr(),
        reached.unsafe_ptr(),
        frame_mode,
        True,
    )
    # Each pair's span in the frame: rows `0`, then `1`, columns `2` and `3`.
    var spans = List[Int](length=4 * pairs, fill=0)
    var spanning = List[Bool](length=pairs, fill=False)
    var starts_free = frame_mode.reference_start > 0 or frame_mode.query_start > 0
    var back_references = List[String](length=pairs, fill=String())
    var back_queries = List[String](length=pairs, fill=String())
    var backed = List[Bool](length=pairs, fill=True)
    for index in range(pairs):
        if not reached[index]:
            continue
        if not found[index]:
            costs_out[unsafe_offset=index] = None
            settled[unsafe_offset=index] = True
            continue
        var rows = firsts.length(index)
        var columns = seconds.length(index)
        var diagonal = ends[index]
        var end_row = rows if diagonal <= columns - rows else columns - diagonal
        var end_column = rows + diagonal if diagonal <= columns - rows else columns
        if end_row == 0 or end_column == 0:
            continue
        spanning[index] = True
        spans[4 * index + 1] = end_row
        spans[4 * index + 3] = end_column
        if starts_free:
            back_references[index] = reversed_text(Span(unsafe_ptr=firsts.letters(index), length=end_row))
            back_queries[index] = reversed_text(Span(unsafe_ptr=seconds.letters(index), length=end_column))
            backed[index] = False

    # The start, back from the end over its reversed prefixes, the frame's free starts their free ends: of
    # the cells an optimal alignment leaves from, the one on the lowest diagonal in the frame, the highest
    # here.
    if starts_free:
        var back_mode = Mode(
            Mode.ENDS, 0, frame_mode.reference_start, 0, frame_mode.query_start, 0, frame_mode.anchor, -1, -1
        )
        var back_found = List[Optional[Int]](length=pairs, fill=None)
        var back_ends = List[Int](length=pairs, fill=0)
        lane_ends(
            pairs,
            StringTexts.of(back_references),
            StringTexts.of(back_queries),
            costs,
            max_cost,
            workers,
            back_found.unsafe_ptr(),
            back_ends.unsafe_ptr(),
            backed.unsafe_ptr(),
            back_mode,
            False,
        )
        _ = back_references^
        _ = back_queries^
        for index in range(pairs):
            if not spanning[index]:
                continue
            spanning[index] = False
            if not backed[index] or not back_found[index] or back_found[index].value() != found[index].value():
                continue
            var rows = spans[4 * index + 1]
            var columns = spans[4 * index + 3]
            var diagonal = back_ends[index]
            var back_row = rows if diagonal <= columns - rows else columns - diagonal
            var back_column = rows + diagonal if diagonal <= columns - rows else columns
            spans[4 * index] = rows - back_row
            spans[4 * index + 2] = columns - back_column
            spanning[index] = True
    return FramedSpans(found^, spans^, spanning^)


def lane_free_alignments[
    T: Texts
](
    pairs: Int,
    references: T,
    queries: T,
    costs: LaneCosts,
    mode: Mode,
    max_cost: Int,
    right: Bool,
    workers: Int,
    budget: Int,
    costs_out: MutPointer[Optional[Int], _],
    moves_out: MutPointer[List[UInt8], _],
    spans_out: MutPointer[Int, _],
    settled: MutPointer[Bool, _],
):
    """Every pair's alignment with `mode`'s free ends, no band, in lanes, as `gap_affine.free_ends_alignment`
    finds it one pair at a time: its cost into `costs_out`, None past `max_cost`; the span it aligns into
    `spans_out`, four numbers a pair, the reference's first letter and the one past its last, then the
    query's; and the moves of the global alignment of that span, right to left, into `moves_out`; and
    `settled` set. A pair the lanes do not settle at every step is left to the caller.

    The span is the rule's (see `Ties`), run over both sequences reversed with `right`: the end on the
    lowest of the lanes' diagonals, query less reference, an optimal alignment reaches, found by `lane_ends`
    from the starts; then the start on the lowest such diagonal too, of those an optimal alignment ending
    there leaves, by `lane_ends` again, back from that end over its reversed prefixes; then the span's own
    global alignment by `lane_alignments`, the tie rule the other way round with `right`, at the same
    cost. A pair whose end lies on the first row or column, a path running along it, is left."""
    if pairs == 0:
        return
    var framed: FramedSpans
    if right:
        # The rule's frame: both sequences reversed, the free letters swapped end for end.
        var frame_mode = Mode(
            Mode.ENDS,
            mode.reference_end,
            mode.reference_start,
            mode.query_end,
            mode.query_start,
            0,
            mode.anchor,
            -1,
            -1,
        )
        var frame_references = List[String](length=pairs, fill=String())
        var frame_queries = List[String](length=pairs, fill=String())
        var reference_ptr = frame_references.unsafe_ptr()
        var query_ptr = frame_queries.unsafe_ptr()

        def turn(stretch: Int) {imm references, imm queries, imm reference_ptr, imm query_ptr, imm pairs, imm workers}:
            """Stretch `stretch`'s pairs back to front."""
            for index in range(pairs * stretch // workers, pairs * (stretch + 1) // workers):
                reference_ptr[unsafe_offset=index] = reversed_text(
                    Span(unsafe_ptr=references.letters(index), length=references.length(index))
                )
                query_ptr[unsafe_offset=index] = reversed_text(
                    Span(unsafe_ptr=queries.letters(index), length=queries.length(index))
                )

        spread(turn, workers, workers)
        framed = framed_spans(
            pairs,
            StringTexts.of(frame_references),
            StringTexts.of(frame_queries),
            costs,
            frame_mode,
            max_cost,
            workers,
            costs_out,
            settled,
        )
        _ = frame_references^
        _ = frame_queries^
    else:
        framed = framed_spans(pairs, references, queries, costs, mode, max_cost, workers, costs_out, settled)

    # The span's global alignment, at the same cost, the tie rule its own way.
    ref spans = framed.spans
    var span_references = List[String](length=pairs, fill=String())
    var span_queries = List[String](length=pairs, fill=String())
    var span_settled = List[Bool](length=pairs, fill=True)
    for index in range(pairs):
        if not framed.spanning[index]:
            continue
        var rows = references.length(index)
        var columns = queries.length(index)
        # In the batch's own order, unreversed.
        var row_start = rows - spans[4 * index + 1] if right else spans[4 * index]
        var row_end = rows - spans[4 * index] if right else spans[4 * index + 1]
        var column_start = columns - spans[4 * index + 3] if right else spans[4 * index + 2]
        var column_end = columns - spans[4 * index + 2] if right else spans[4 * index + 3]
        spans[4 * index] = row_start
        spans[4 * index + 1] = row_end
        spans[4 * index + 2] = column_start
        spans[4 * index + 3] = column_end
        span_references[index] = String(
            unsafe_from_utf8=Span(
                unsafe_ptr=references.letters(index).unsafe_offset(row_start), length=row_end - row_start
            )
        )
        span_queries[index] = String(
            unsafe_from_utf8=Span(
                unsafe_ptr=queries.letters(index).unsafe_offset(column_start), length=column_end - column_start
            )
        )
        span_settled[index] = False
    var span_costs = List[Optional[Int]](length=pairs, fill=None)
    _ = lane_alignments(
        pairs,
        StringTexts.of(span_references),
        StringTexts.of(span_queries),
        costs,
        Band(),
        max_cost,
        not right,
        workers,
        budget,
        span_costs.unsafe_ptr(),
        moves_out,
        span_settled.unsafe_ptr(),
    )
    _ = span_references^
    _ = span_queries^
    for index in range(pairs):
        if not framed.spanning[index] or not span_settled[index] or not span_costs[index]:
            continue
        if span_costs[index].value() != framed.found[index].value():
            continue
        costs_out[unsafe_offset=index] = span_costs[index]
        comptime for field in range(4):
            spans_out[unsafe_offset=4 * index + field] = spans[4 * index + field]
        settled[unsafe_offset=index] = True


@fieldwise_init
struct LocalCosts(ImplicitlyCopyable, TrivialRegisterPassable):
    """A local alignment's scores for the lanes: a match's and a mismatch's, and each gap piece's first
    letter's and every further letter's, a deletion's, a reference letter alone, apart from an
    insertion's; every one below zero but a match's. A table's lanes read their pairs from it, `mismatch`
    then its least entry."""

    var hit: Int
    var mismatch: Int
    var deletion_open: Int
    var deletion_extend: Int
    var insertion_open: Int
    var insertion_extend: Int
    var pieces: Int
    var deletion_open2: Int
    var deletion_extend2: Int
    var insertion_open2: Int
    var insertion_extend2: Int

    @staticmethod
    def symmetric(hit: Int, mismatch: Int, open: Int, extend: Int) -> Self:
        """One gap piece, a gap of `k` letters `open + (k - 1) extend` either way."""
        return Self(hit, mismatch, open, extend, open, extend, 1, 0, 0, 0, 0)

    @staticmethod
    def of(costs: Costs, match_score: Int) -> Self:
        """`costs` as scores, every match earning `match_score`."""
        return Self(
            match_score,
            -costs.mismatch,
            -(costs.deletion_opening + costs.deletion_extension),
            -costs.deletion_extension,
            -(costs.opening + costs.extension),
            -costs.extension,
            costs.pieces(),
            -(costs.deletion_opening2 + costs.deletion_extension2),
            -costs.deletion_extension2,
            -(costs.opening2 + costs.extension2),
            -costs.extension2,
        )

    def dearest(self) -> Int:
        """The dearest single move, a mismatch or a gap's first letter."""
        var dearest = max(-self.mismatch, max(-self.deletion_open, -self.insertion_open))
        if self.pieces == 2:
            dearest = max(dearest, max(-self.deletion_open2, -self.insertion_open2))
        return dearest

    def pads_lose(self) -> Bool:
        """Whether no move but a match earns, as padding needs: a cell past a lane's letters then scores no
        more than the cell it came from, so never more than the lane's best."""
        var losing = self.mismatch <= 0 and max(self.deletion_open, self.insertion_open) <= 0
        losing = losing and max(self.deletion_extend, self.insertion_extend) <= 0
        if self.pieces == 2:
            losing = losing and max(self.deletion_open2, self.insertion_open2) <= 0
            losing = losing and max(self.deletion_extend2, self.insertion_extend2) <= 0
        return losing


def lane_local_scores[
    R: Texts, Q: Texts
](
    pairs: Int,
    references: R,
    queries: Q,
    costs: LocalCosts,
    pads: Tuple[UInt8, UInt8],
    workers: Int,
    scores_out: MutPointer[Int32, _],
    settled: MutPointer[Bool, _],
    letters: Int = 0,
    table: SIMD[DType.uint8, TABLE_ENTRIES] = SIMD[DType.uint8, TABLE_ENTRIES](0),
) -> Int:
    """Every pair `settled` does not mark, its best local alignment's score into `scores_out` and `settled`
    set, many pairs at once, a pair a lane, both sequences swept whole, as Smith and Waterman's recurrence
    with Gotoh's gaps runs, as SWIPE scores a database (Rognes, 2011): first in 16 bits, every pair whose
    scores they hold, then the rest in 32. The pairs settled are counted.

    A lane past its own letters reads `pads`, a byte for the references and one for the queries that
    match no letter and not each other, so a cell past them scores no more than the cell it came from and
    never beats the lane's best; which needs no move but a match to earn. Otherwise every pair is left.

    With `letters`, the sequences are an alphabet's codes and each pair scores `table`'s signed byte at
    `a letters + b`, `costs.mismatch` its least; a pad is then any code past the alphabet, and scores
    that."""
    if not costs.pads_lose() or pairs == 0:
        return 0
    var narrow = local_stage[R, Q, DType.int16](
        pairs, references, queries, costs, pads, workers, scores_out, settled, letters, table
    )
    var wide = local_stage[R, Q, DType.int32](
        pairs, references, queries, costs, pads, workers, scores_out, settled, letters, table
    )
    return narrow + wide


def local_width[value: DType]() -> Int:
    """A local score's lanes a group: 32 in 16 bits and 16 in 32, several registers' worth where the CPU's
    are narrower: on the M2's 128 bits a step of 32 lanes, four registers' work with none waiting on
    another, scored 20,000 references against a read in 83 ms where one register's 8 took 160, and 20,000
    pairs under a table in 108 where they took 230."""
    return max(lanes_of[value](), 32 if value == DType.int16 else 16)


def local_stage[
    R: Texts, Q: Texts, value: DType
](
    pairs: Int,
    references: R,
    queries: Q,
    costs: LocalCosts,
    pads: Tuple[UInt8, UInt8],
    workers: Int,
    scores_out: MutPointer[Int32, _],
    settled: MutPointer[Bool, _],
    letters: Int,
    table: SIMD[DType.uint8, TABLE_ENTRIES],
) -> Int:
    """`lane_local_scores` in lanes of `value`, for every pair `settled` does not mark whose scores they
    hold: a match's reward over the shorter sequence and the dearest move each well inside them."""
    comptime WIDTH = local_width[value]()
    comptime LIMIT = 16000 if value == DType.int16 else (1 << 29)
    if costs.dearest() >= LIMIT // 4:
        return 0
    # The pairs a lane holds, by their references' lengths, so a group's lanes pad few rows.
    var keys = List[Int](capacity=pairs)
    for index in range(pairs):
        var rows = references.length(index)
        var columns = queries.length(index)
        if not settled[unsafe_offset=index] and max(costs.hit, 0) * min(rows, columns) < LIMIT:
            keys.append((min(rows, (1 << 30) - 1) << 32) | index)
    sort(keys)
    var held_pairs = len(keys)
    if held_pairs == 0:
        return 0
    var order = List[Int](capacity=held_pairs)
    for key in keys:
        order.append(key & ((1 << 32) - 1))
    var order_ptr = order.unsafe_ptr()
    var spaces = List[LaneSpace[value]](capacity=workers)
    for _ in range(workers):
        spaces.append(LaneSpace[value]())
    var space_ptr = spaces.unsafe_ptr()
    var groups = ceildiv(held_pairs, WIDTH)
    var taken = Atomic[Int64](0)
    var group_workers = min(workers, groups)

    def sweep(
        worker: Int,
    ) {
        mut taken,
        imm references,
        imm queries,
        imm costs,
        imm pads,
        imm letters,
        imm table,
        imm groups,
        imm group_workers,
        imm held_pairs,
        imm order_ptr,
        imm space_ptr,
        imm scores_out,
        imm settled,
    }:
        """Takes groups until none is left, scoring each pair."""
        ref space = space_ptr[unsafe_offset=worker]
        var last_share = 0
        while True:
            var share = next_share(taken, groups, group_workers, last_share)
            if share[0] >= groups:
                return
            for group in range(share[0], share[1]):
                space.members.clear()
                for slot in range(group * WIDTH, min(held_pairs, (group + 1) * WIDTH)):
                    space.members.append(order_ptr[unsafe_offset=slot])
                var best: SIMD[value, WIDTH]
                if letters > 0:
                    if costs.pieces == 2:
                        best = local_band_scores[R, Q, value, 2, True](
                            references, queries, space, costs, pads, letters, table
                        )
                    else:
                        best = local_band_scores[R, Q, value, 1, True](
                            references, queries, space, costs, pads, letters, table
                        )
                elif costs.pieces == 2:
                    best = local_band_scores[R, Q, value, 2](references, queries, space, costs, pads)
                else:
                    best = local_band_scores[R, Q, value, 1](references, queries, space, costs, pads)
                for lane in range(len(space.members)):
                    scores_out[unsafe_offset=space.members[lane]] = Int32(best[lane])
                    settled[unsafe_offset=space.members[lane]] = True

    spread(sweep, group_workers, group_workers)
    return held_pairs


def local_band_scores[
    R: Texts, Q: Texts, value: DType, pieces: Int, tabled: Bool = False
](
    references: R,
    queries: Q,
    mut space: LaneSpace[value],
    costs: LocalCosts,
    pads: Tuple[UInt8, UInt8],
    letters: Int = 0,
    table: SIMD[DType.uint8, TABLE_ENTRIES] = SIMD[DType.uint8, TABLE_ENTRIES](0),
) -> SIMD[value, local_width[value]()]:
    """The best local score of each of `space.members`'s pairs, the whole matrix swept a row at a time
    for every lane at once (see `lane_local_scores`): a row a reference letter, so a step down a row is a
    deletion and one across an insertion. Its rows of each layer are its own, `local_width` lanes wide."""
    comptime WIDTH = local_width[value]()
    comptime Lanes = SIMD[value, WIDTH]
    comptime two = pieces == 2
    var rows = 0
    var columns = 0
    for lane in range(len(space.members)):
        rows = max(rows, references.length(space.members[lane]))
        columns = max(columns, queries.length(space.members[lane]))
    space.row_letters.resize(unsafe_uninit_length=ceildiv(max(rows, 1), BLOCK) * BLOCK * WIDTH)
    space.column_letters.resize(unsafe_uninit_length=ceildiv(max(columns, 1), BLOCK) * BLOCK * WIDTH)
    side_by_side[R, WIDTH](
        references, space.members, rows, space.staging, space.row_letters.unsafe_ptr(), False, Int(pads[0])
    )
    side_by_side[Q, WIDTH](
        queries, space.members, columns, space.staging, space.column_letters.unsafe_ptr(), False, Int(pads[1])
    )
    # A gap layer no path has entered: below any score a gap opened from a cell, which is at least the
    # dearest gap's first letter.
    comptime NONE = -16000 if value == DType.int16 else -(1 << 29)
    var none = Lanes(Scalar[value](NONE))
    var score_row = List[Lanes](length=columns + 1, fill=Lanes(0))
    var deletion_row = List[Lanes](length=columns + 1, fill=none)
    var deletion2_row = List[Lanes](length=columns + 1 if two else 0, fill=none)
    var scores = score_row.unsafe_ptr()
    var deletions = deletion_row.unsafe_ptr()
    var deletions2 = deletion2_row.unsafe_ptr()
    var matched = Lanes(Scalar[value](costs.hit))
    var mismatched = Lanes(Scalar[value](costs.mismatch))
    var delete_open = Lanes(Scalar[value](costs.deletion_open))
    var delete_extend = Lanes(Scalar[value](costs.deletion_extend))
    var insert_open = Lanes(Scalar[value](costs.insertion_open))
    var insert_extend = Lanes(Scalar[value](costs.insertion_extend))
    var delete_open2 = Lanes(Scalar[value](costs.deletion_open2))
    var delete_extend2 = Lanes(Scalar[value](costs.deletion_extend2))
    var insert_open2 = Lanes(Scalar[value](costs.insertion_open2))
    var insert_extend2 = Lanes(Scalar[value](costs.insertion_extend2))
    var zero = Lanes(0)
    var best = Lanes(0)
    var row_source = space.row_letters.unsafe_ptr()
    var column_source = space.column_letters.unsafe_ptr()
    var alphabet = SIMD[DType.uint8, WIDTH](UInt8(letters))
    for row in range(1, rows + 1):
        var letter = row_source.unsafe_offset((row - 1) * WIDTH).unsafe_load[width=WIDTH]()
        var row_entries = letter * UInt8(letters)
        var row_padded = letter.ge(alphabet)
        var diagonal = Lanes(0)
        var left = Lanes(0)
        var insertion = none
        var insertion2 = none
        for column in range(1, columns + 1):
            var other = column_source.unsafe_offset((column - 1) * WIDTH).unsafe_load[width=WIDTH]()
            var above = scores[unsafe_offset=column]
            var deletion = max(above + delete_open, deletions[unsafe_offset=column] + delete_extend)
            insertion = max(left + insert_open, insertion + insert_extend)
            var paired: Lanes
            comptime if tabled:
                paired = (row_padded | other.ge(alphabet)).select(
                    mismatched, looked_up[value, WIDTH](table, row_entries + other)
                )
            else:
                paired = letter.eq(other).select(matched, mismatched)
            var score = max(max(diagonal + paired, deletion), max(insertion, zero))
            comptime if two:
                var deletion2 = max(above + delete_open2, deletions2[unsafe_offset=column] + delete_extend2)
                insertion2 = max(left + insert_open2, insertion2 + insert_extend2)
                deletions2[unsafe_offset=column] = deletion2
                score = max(score, max(deletion2, insertion2))
            deletions[unsafe_offset=column] = deletion
            scores[unsafe_offset=column] = score
            best = max(best, score)
            diagonal = above
            left = score
    return best
