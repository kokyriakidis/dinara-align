# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Ported from `pa-bitpacking` in A*PA (https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner,
# commit bf2e14e), by Ragnar Groot Koerkamp and Pesho Ivanov, itself translated from Edlib.
"""
Bit-parallel global edit distance over DNA, after Myers (1999) as A*PA2 implements it.

The second sequence is packed 64 bases to a word, so one machine word holds a 64-row slice of a
column, and a column of the matrix advances with twenty bitwise operations per word. Each word
carries the vertical differences down its slice as a `+` mask and a `-` mask, and a single
horizontal difference crosses from one word to the word below.

Eight words go through one SIMD vector at a time, staggered by a column each: the lane holding the
lower word works one column behind the lane above it, so the difference the upper lane emits is
exactly the one the lower lane needs next, handed over by rotating the lanes. The stagger fills and
drains with the lanes outside the tile masked off, rather than a word at a time. A band only two or
three words tall runs the same stagger in scalar registers instead, where an operation takes one
cycle rather than a vector's two, since a band that narrow is bound by the chain, not the work.

Myers' step is regrouped so the next column waits on fewer operations than in A*PA2 (see `advance`).

The matrix is rarely swept whole. As in A*PA2-simple, band doubling guesses a bound, computes only
the cells a path within it could cross, prunes rows whose score already rules them out, and repeats
with a larger guess until the distance fits under it (see `pruned_distance`). Unlike A*PA2, the next
guess is aimed from where the failed round died rather than doubled (see `band_doubling`).

Before any band, diagonal transition, as WFA runs it, searches from the start for as long as it
stays cheaper than a band would be (see `diagonal_transition`). Near-identical pairs, which the post
on A*PA2 names as its weak spot, finish there, traceback included; any other pair leaves with a
projection of its distance, which becomes the band's first bound.

A pair too divergent for a band is swept whole, cut into tiles of `LANES` words by `TILE_COLUMNS`
columns. A tile needs only the tile above it and the tile to its left, so every tile on one
anti-diagonal of tiles runs on its own thread; they touch disjoint columns of the horizontal edge
and disjoint words of the vertical one. A*PA sweeps on one thread.

Unit costs only: a substitution, an insertion and a deletion each cost one, so this is the
distance `levenshtein_alignment` returns, without the alignment.
"""

from std.atomic import Atomic
from std.bit import byte_swap, count_leading_zeros, count_trailing_zeros, pop_count
from std.math import ceildiv, sqrt
from std.sys import inlined_assembly
from std.time import perf_counter_ns

from max.algorithm import parallelize

from .alignment import AlignmentResult
from .common import hardware_threads
from .errors import AlignmentError, ErrorKind

comptime WORD_BITS = 64
"""Rows one word of the second sequence covers."""

comptime LANES = 8
"""
Words advanced together in one SIMD vector, and so the rows of a block: `LANES * WORD_BITS`.

Eight measured fastest on an M2, where four left the pipeline idle and sixteen ran out of registers.
That matches A*PA2's own choice of two four-lane vectors side by side.
"""

comptime NARROW_LANES = 4
"""
The smaller group a sweep falls back to below `LANES` words, A*PA2's own width, so a narrow band rounds
its bottom up to four words rather than eight.
"""

comptime TILE_COLUMNS = 2048
"""
Columns of one parallel tile. Each tile repeats the stagger's two triangles, about `LANES * LANES`
single-word steps, so a tile must be wide enough for those to vanish beside its `LANES * width`
vector steps, and narrow enough that a long pair still yields more tiles per diagonal than threads.
"""

comptime PARALLEL_CELLS = 64_000_000
"""Cells below which one thread finishes before a fork would pay for itself."""

comptime BAND_COLUMNS = 256
"""
Columns of one band tile, A*PA2's block width. A tile computes every word its columns reach, so a
diagonal band costs about this many extra rows per tile, against the stagger's triangles repaid on
every narrower one.
"""

comptime DIAGONAL = UInt8(0)
"""A traceback move that consumes a base of each sequence, a match or a substitution."""
comptime LEFT = UInt8(1)
"""A move that consumes a base of the first sequence against a gap."""
comptime UP = UInt8(2)
"""A move that consumes a base of the second sequence against a gap."""

comptime WAVEFRONT_FLOOR = 32
"""The fewest edits a tile's wavefront search may spend before the tile is recomputed instead."""

comptime TRACE_PADDING = 2
"""Unreached diagonals stored either side of a traceback front, so the next reads its neighbours unchecked."""

comptime FRONT_DROP = 20
"""
How far behind the furthest front, in column plus row, a front may fall before the search drops it;
about ten diagonal steps, A*PA2's `fr_drop`.
"""

comptime NARROW_COLUMNS = 64
"""
Columns of one band tile when the band itself is narrow. A tile computes every word its columns reach,
so a band slopes down by its own width across each tile; on a band only a few hundred rows tall that
slope is most of the work, and narrower tiles cut it at the cost of more triangles.
"""

comptime NARROW_BAND = 1024
"""The bound below which band tiles are `NARROW_COLUMNS` wide rather than `BAND_COLUMNS`."""

comptime MEET_COLUMNS = 2048
"""Columns below which one direction finishes before a second thread would pay for itself."""

comptime CODE_PADDING = 16
"""Sentinel bytes after each sequence's codes, enough for a sixteen-byte comparison at the last base."""

comptime FIRST_SENTINEL = UInt8(0xFE)
"""Past the first sequence's last base: no code, and unequal to `SECOND_SENTINEL`."""

comptime SECOND_SENTINEL = UInt8(0xFF)
"""Past the second sequence's last base."""

comptime COLUMN_PADDING = LANES
"""Columns of zero bases before and after the profile's planes, for lanes standing outside a tile."""

comptime Words = SIMD[DType.uint64, LANES]
comptime ALL_ONES = ~UInt64(0)


@always_inline
def opaque[width: Int](value: SIMD[DType.uint64, width]) -> SIMD[DType.uint64, width]:
    """`value`, unchanged, hidden from the optimizer so it cannot fold it back into a longer chain.

    Only up to four lanes: eight are bound by how many vector operations issue, not by the chain,
    and there the folded form is the cheaper one.
    """
    comptime if width == 1:
        return inlined_assembly["", SIMD[DType.uint64, width], constraints="=r,0", has_side_effect=False](value)
    elif width == 2:
        return inlined_assembly["", SIMD[DType.uint64, width], constraints="=w,0", has_side_effect=False](value)
    elif width > 4:
        return value
    else:
        comptime half = width // 2
        var low = opaque[half](value.slice[half, offset=0]())
        var high = opaque[half](value.slice[half, offset=half]())
        return rebind[SIMD[DType.uint64, width]](low.join(high))


@always_inline
def advance[
    width: Int
](
    mut horizontal_plus: SIMD[DType.uint64, width],
    mut horizontal_minus: SIMD[DType.uint64, width],
    mut vertical_plus: SIMD[DType.uint64, width],
    mut vertical_minus: SIMD[DType.uint64, width],
    matches: SIMD[DType.uint64, width],
):
    """One column of one word: Myers' step, as Edlib and A*PA write it, regrouped for latency.

    The horizontal difference enters at the top of the word in bit zero and the one leaving at the
    bottom replaces it; the vertical masks are updated in place. The addition is what lets a
    vertical run of matches propagate down the whole word at once.

    A narrow band is bound by how long one column waits on the last, not by how many operations it
    takes, so the step is regrouped to shorten that wait, and `opaque` keeps the compiler from
    folding the regrouping back: a single word went from fourteen cycles a column to eight.
    """
    var crossing = matches | vertical_minus
    # Myers assumes the incoming horizontal difference is never -1; A*PA folds it into the matches.
    var equal = matches | horizontal_minus
    # Myers' horizontal mask is `((sum ^ vp) | equal)`; both its uses below are rewritten on `sum`
    # directly, `horizontal | vp` as `sum | equal | vp` and `vp & horizontal` as
    # `(vp & ~sum) | (vp & equal)`, which takes two operations off the path one column waits on.
    var sum = (equal & vertical_plus) + vertical_plus
    var plus = vertical_minus | ~(sum | (equal | vertical_plus))
    var minus = (vertical_plus & ~sum) | (vertical_plus & equal)
    var plus_out = plus >> (WORD_BITS - 1)
    var minus_out = minus >> (WORD_BITS - 1)
    var plus_shifted = plus << 1
    # `~(crossing | plus)` with the incoming `+` folded into `crossing`, which is known early.
    var open = opaque(~(crossing | horizontal_plus))
    minus = opaque((minus << 1) | horizontal_minus)
    plus = plus_shifted | horizontal_plus
    horizontal_plus = plus_out
    horizontal_minus = minus_out
    vertical_plus = minus | (open & ~plus_shifted)
    vertical_minus = plus & crossing


def base_code(byte: UInt8) raises AlignmentError -> Int:
    """`A`, `C`, `G` and `T` as zero to three, the two bits the profile is built from."""
    if byte == UInt8(ord("A")):
        return 0
    if byte == UInt8(ord("C")):
        return 1
    if byte == UInt8(ord("G")):
        return 2
    if byte == UInt8(ord("T")):
        return 3
    raise AlignmentError(ErrorKind.UNKNOWN_SYMBOL, String("bit-parallel edit distance takes ACGT, not byte ", byte))


struct Sweep(ImplicitlyCopyable, TrivialRegisterPassable):
    """Pointers to the profile and the frontier, which the caller's lists own and outlive.

    The profile holds both sequences as bit planes, so a word of matches is two XORs and an AND. A
    base of the first sequence becomes two whole-word masks, all ones where its code has that bit;
    the second is packed 64 bases to a word with each plane stored negated, so a base equals a
    column's base exactly where both planes XOR to ones.

    The frontier holds the differences along the two edges still being computed: one horizontal
    difference per column, in bit zero, along the bottom of the rows done so far, and one word of
    vertical differences per row word, down the right edge of the columns done so far.

    Plain pointers rather than the lists themselves, so tiles on different threads can each write
    their own disjoint slice of the same frontier.
    """

    var column_low: MutPointer[UInt64, MutUntrackedOrigin]
    var column_high: MutPointer[UInt64, MutUntrackedOrigin]
    var row_low: MutPointer[UInt64, MutUntrackedOrigin]
    var row_high: MutPointer[UInt64, MutUntrackedOrigin]
    var horizontal_plus: MutPointer[UInt64, MutUntrackedOrigin]
    var horizontal_minus: MutPointer[UInt64, MutUntrackedOrigin]
    var vertical_plus: MutPointer[UInt64, MutUntrackedOrigin]
    var vertical_minus: MutPointer[UInt64, MutUntrackedOrigin]

    def __init__(
        out self,
        mut column_low: List[UInt64],
        mut column_high: List[UInt64],
        mut row_low: List[UInt64],
        mut row_high: List[UInt64],
        mut horizontal_plus: List[UInt64],
        mut horizontal_minus: List[UInt64],
        mut vertical_plus: List[UInt64],
        mut vertical_minus: List[UInt64],
    ):
        # Past the padding, so column `c` of the matrix is element `c` here.
        self.column_low = column_low.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]().unsafe_offset(COLUMN_PADDING)
        self.column_high = (
            column_high.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]().unsafe_offset(COLUMN_PADDING)
        )
        self.row_low = row_low.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        self.row_high = row_high.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        self.horizontal_plus = horizontal_plus.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        self.horizontal_minus = horizontal_minus.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        self.vertical_plus = vertical_plus.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()
        self.vertical_minus = vertical_minus.unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]()

    def shifted(self, first_column: Int) -> Self:
        """The same sweep reading the horizontal edge from a buffer that starts at `first_column`.

        A band tile's top word always reads `+1` from above, so the horizontal edge only matters within
        one tile; a buffer one tile wide, reset per tile, replaces an array as long as the sequence.
        """
        var moved = self
        moved.horizontal_plus = self.horizontal_plus.unsafe_offset(-first_column)
        moved.horizontal_minus = self.horizontal_minus.unsafe_offset(-first_column)
        return moved

    @always_inline
    def run(self, word: Int, first_column: Int, end_column: Int):
        """One word through `[first_column, end_column)`, its vertical differences held in registers.

        Only the horizontal edge goes through memory, and each column's is a different address, so
        consecutive columns do not wait on a store being read back.
        """
        var vp = self.vertical_plus[unsafe_offset=word]
        var vm = self.vertical_minus[unsafe_offset=word]
        var row_low = self.row_low[unsafe_offset=word]
        var row_high = self.row_high[unsafe_offset=word]
        for column in range(first_column, end_column):
            var hp = self.horizontal_plus[unsafe_offset=column]
            var hm = self.horizontal_minus[unsafe_offset=column]
            var matches = (self.column_low[unsafe_offset=column] ^ row_low) & (
                self.column_high[unsafe_offset=column] ^ row_high
            )
            advance[1](hp, hm, vp, vm, matches)
            self.horizontal_plus[unsafe_offset=column] = hp
            self.horizontal_minus[unsafe_offset=column] = hm
        self.vertical_plus[unsafe_offset=word] = vp
        self.vertical_minus[unsafe_offset=word] = vm

    def words(self, first_word: Int, end_word: Int, first_column: Int, end_column: Int):
        """Words `[first_word, end_word)` through columns `[first_column, end_column)`.

        Full groups of `LANES` take the staggered vector sweep when the span is wide enough for it.
        Of what is left, four words take a narrow vector, five to seven a narrow vector with the rest
        in scalar registers beside it, two or three scalar registers alone, and one word goes on its own.
        """
        var word = first_word
        if end_column - first_column >= 2 * LANES:
            while word + LANES <= end_word:
                self.block[LANES](word, first_column, end_column)
                word += LANES
            # Five to seven words left: a narrow vector with the rest in scalar registers in its shadow.
            var left = end_word - word
            if left == 7:
                self.hybrid_block[NARROW_LANES, 3](word, first_column, end_column)
                word += 7
            elif left == 6:
                self.hybrid_block[NARROW_LANES, 2](word, first_column, end_column)
                word += 6
            elif left == 5:
                self.hybrid_block[NARROW_LANES, 1](word, first_column, end_column)
                word += 5
            elif left == 4:
                self.block[NARROW_LANES](word, first_column, end_column)
                word += NARROW_LANES
            # Two or three words left run side by side in scalar registers, which beats a vector
            # this narrow: its chain waits two cycles an operation, a scalar one.
            if end_word - word == 3:
                self.scalar_block[3](word, first_column, end_column)
                word += 3
            elif end_word - word == 2:
                self.scalar_block[2](word, first_column, end_column)
                word += 2
        while word < end_word:
            self.run(word, first_column, end_column)
            word += 1

    def block[lanes: Int](self, first_word: Int, first_column: Int, end_column: Int):
        """`lanes` words through `[first_column, end_column)` in one vector, staggered (see `VectorGroup`)."""
        var group = VectorGroup[lanes](self, first_word, first_column, end_column)
        stagger(group)
        group.finish()

    def scalar_block[lanes: Int](self, first_word: Int, first_column: Int, end_column: Int):
        """`lanes` words through `[first_column, end_column)` in scalar registers, staggered (see `ScalarGroup`)."""
        var group = ScalarGroup[lanes](self, first_word, first_column, end_column)
        stagger(group)
        group.finish()

    def hybrid_block[lanes: Int, below: Int](self, first_word: Int, first_column: Int, end_column: Int):
        """`lanes` words in a vector and the `below` words under them in scalar registers, in one loop.

        A narrow vector's chain leaves the integer units idle, so the scalar words run in its
        shadow, a few columns behind so the differences the vector sends down are already stored.
        """
        var top = VectorGroup[lanes](self, first_word, first_column, end_column)
        var bottom = ScalarGroup[below](self, first_word + lanes, first_column, end_column)
        stagger_pair(top, bottom, below + 2)
        top.finish()
        bottom.finish()


trait Staggered(Movable):
    """A group of words swept through one tile, each a column behind the word above it.

    At step `offset`, the group's `k`-th word from the bottom works column `offset + 1 + k`, so a
    word takes the difference the word above it sent down on the step before. The group takes the
    difference entering its top word from the horizontal edge, and leaves its bottom word's there.
    """

    def width(self) -> Int:
        """How many words, and so how many steps the stagger takes to fill."""
        ...

    def first_column(self) -> Int:
        ...

    def end_column(self) -> Int:
        ...

    def step[masked: Bool](mut self, offset: Int):
        """One step. `masked` when some word stands outside the tile: such a word keeps its state,
        and the edge is read and written only for columns inside it."""
        ...

    def finish(self):
        """Writes the words' vertical differences back to the frontier."""
        ...


struct VectorGroup[lanes: Int](Staggered, TrivialRegisterPassable):
    """`lanes` words in one SIMD vector: lane `k` holds the `k`-th word from the bottom.

    Each lane takes the difference the lane above sent by rotating the lanes, and the top lane
    takes the one entering from the edge. A lane outside the tile reads bases from the profile's
    padding and keeps its state.
    """

    var sweep: Sweep
    var first_word: Int
    var first: Int
    var end: Int
    var row_low: SIMD[DType.uint64, Self.lanes]
    var row_high: SIMD[DType.uint64, Self.lanes]
    var vertical_plus: SIMD[DType.uint64, Self.lanes]
    var vertical_minus: SIMD[DType.uint64, Self.lanes]
    var horizontal_plus: SIMD[DType.uint64, Self.lanes]
    var horizontal_minus: SIMD[DType.uint64, Self.lanes]
    var lane_columns: SIMD[DType.int, Self.lanes]

    @always_inline
    def __init__(out self, sweep: Sweep, first_word: Int, first_column: Int, end_column: Int):
        self.sweep = sweep
        self.first_word = first_word
        self.first = first_column
        self.end = end_column
        self.row_low = SIMD[DType.uint64, Self.lanes]()
        self.row_high = SIMD[DType.uint64, Self.lanes]()
        self.vertical_plus = SIMD[DType.uint64, Self.lanes]()
        self.vertical_minus = SIMD[DType.uint64, Self.lanes]()
        self.horizontal_plus = SIMD[DType.uint64, Self.lanes]()
        self.horizontal_minus = SIMD[DType.uint64, Self.lanes]()
        self.lane_columns = SIMD[DType.int, Self.lanes]()
        comptime for lane in range(Self.lanes):
            var word = first_word + Self.lanes - 1 - lane
            self.row_low[lane] = sweep.row_low[unsafe_offset=word]
            self.row_high[lane] = sweep.row_high[unsafe_offset=word]
            self.vertical_plus[lane] = sweep.vertical_plus[unsafe_offset=word]
            self.vertical_minus[lane] = sweep.vertical_minus[unsafe_offset=word]
            self.lane_columns[lane] = lane + 1

    @always_inline
    def width(self) -> Int:
        return Self.lanes

    @always_inline
    def first_column(self) -> Int:
        return self.first

    @always_inline
    def end_column(self) -> Int:
        return self.end

    @always_inline
    def step[masked: Bool](mut self, offset: Int):
        self.horizontal_plus = self.horizontal_plus.rotate_left[1]()
        self.horizontal_minus = self.horizontal_minus.rotate_left[1]()
        var top = offset + Self.lanes
        if not masked or (top >= self.first and top < self.end):
            self.horizontal_plus[Self.lanes - 1] = self.sweep.horizontal_plus[unsafe_offset=top]
            self.horizontal_minus[Self.lanes - 1] = self.sweep.horizontal_minus[unsafe_offset=top]
        # The bases lane `k` needs sit side by side, from column `offset + 1`.
        var low = self.sweep.column_low.unsafe_offset(offset + 1).unsafe_load[width=Self.lanes]()
        var high = self.sweep.column_high.unsafe_offset(offset + 1).unsafe_load[width=Self.lanes]()
        var matches = (low ^ self.row_low) & (high ^ self.row_high)
        comptime if masked:
            var kept_plus = self.vertical_plus
            var kept_minus = self.vertical_minus
            advance[Self.lanes](
                self.horizontal_plus, self.horizontal_minus, self.vertical_plus, self.vertical_minus, matches
            )
            var columns = self.lane_columns + offset
            var inside = columns.ge(self.first) & columns.lt(self.end)
            self.vertical_plus = inside.select(self.vertical_plus, kept_plus)
            self.vertical_minus = inside.select(self.vertical_minus, kept_minus)
        else:
            advance[Self.lanes](
                self.horizontal_plus, self.horizontal_minus, self.vertical_plus, self.vertical_minus, matches
            )
        # The bottom word's difference is final, and leaves the group at its column.
        var bottom = offset + 1
        if not masked or (bottom >= self.first and bottom < self.end):
            self.sweep.horizontal_plus[unsafe_offset=bottom] = self.horizontal_plus[0]
            self.sweep.horizontal_minus[unsafe_offset=bottom] = self.horizontal_minus[0]

    @always_inline
    def finish(self):
        comptime for lane in range(Self.lanes):
            var word = self.first_word + Self.lanes - 1 - lane
            self.sweep.vertical_plus[unsafe_offset=word] = self.vertical_plus[lane]
            self.sweep.vertical_minus[unsafe_offset=word] = self.vertical_minus[lane]


struct ScalarGroup[lanes: Int](Staggered):
    """`lanes` words in general-purpose registers, one scalar chain each, staggered as `VectorGroup`.

    Each word's difference passes down through a register rather than a rotation, and a scalar
    operation takes one cycle where a vector one takes two, so a band a few words tall, which is
    all latency, runs faster here. Word `j` is the `j`-th from the top.
    """

    var sweep: Sweep
    var first_word: Int
    var first: Int
    var end: Int
    var row_low: Array[UInt64, Self.lanes]
    var row_high: Array[UInt64, Self.lanes]
    var vertical_plus: Array[UInt64, Self.lanes]
    var vertical_minus: Array[UInt64, Self.lanes]
    var horizontal_plus: Array[UInt64, Self.lanes]
    var horizontal_minus: Array[UInt64, Self.lanes]

    @always_inline
    def __init__(out self, sweep: Sweep, first_word: Int, first_column: Int, end_column: Int):
        self.sweep = sweep
        self.first_word = first_word
        self.first = first_column
        self.end = end_column
        self.row_low = Array[UInt64, Self.lanes](fill=0)
        self.row_high = Array[UInt64, Self.lanes](fill=0)
        self.vertical_plus = Array[UInt64, Self.lanes](fill=0)
        self.vertical_minus = Array[UInt64, Self.lanes](fill=0)
        self.horizontal_plus = Array[UInt64, Self.lanes](fill=0)
        self.horizontal_minus = Array[UInt64, Self.lanes](fill=0)
        comptime for j in range(Self.lanes):
            self.row_low[j] = sweep.row_low[unsafe_offset=first_word + j]
            self.row_high[j] = sweep.row_high[unsafe_offset=first_word + j]
            self.vertical_plus[j] = sweep.vertical_plus[unsafe_offset=first_word + j]
            self.vertical_minus[j] = sweep.vertical_minus[unsafe_offset=first_word + j]

    @always_inline
    def width(self) -> Int:
        return Self.lanes

    @always_inline
    def first_column(self) -> Int:
        return self.first

    @always_inline
    def end_column(self) -> Int:
        return self.end

    @always_inline
    def step[masked: Bool](mut self, offset: Int):
        # Bottom word first, so each word reads what the one above sent before it is replaced.
        comptime for i in range(Self.lanes):
            comptime j = Self.lanes - 1 - i
            var column = offset + Self.lanes - j
            var inside = column >= self.first and column < self.end
            var hp: UInt64
            var hm: UInt64
            comptime if j == 0:
                hp = 0
                hm = 0
                if not masked or inside:
                    hp = self.sweep.horizontal_plus[unsafe_offset=column]
                    hm = self.sweep.horizontal_minus[unsafe_offset=column]
            else:
                hp = self.horizontal_plus[j - 1]
                hm = self.horizontal_minus[j - 1]
            var vp = self.vertical_plus[j]
            var vm = self.vertical_minus[j]
            var matches = (self.sweep.column_low[unsafe_offset=column] ^ self.row_low[j]) & (
                self.sweep.column_high[unsafe_offset=column] ^ self.row_high[j]
            )
            advance[1](hp, hm, vp, vm, matches)
            comptime if masked:
                self.vertical_plus[j] = vp if inside else self.vertical_plus[j]
                self.vertical_minus[j] = vm if inside else self.vertical_minus[j]
            else:
                self.vertical_plus[j] = vp
                self.vertical_minus[j] = vm
            comptime if j == Self.lanes - 1:
                if not masked or inside:
                    self.sweep.horizontal_plus[unsafe_offset=column] = hp
                    self.sweep.horizontal_minus[unsafe_offset=column] = hm
            else:
                self.horizontal_plus[j] = hp
                self.horizontal_minus[j] = hm

    @always_inline
    def finish(self):
        comptime for j in range(Self.lanes):
            self.sweep.vertical_plus[unsafe_offset=self.first_word + j] = self.vertical_plus[j]
            self.sweep.vertical_minus[unsafe_offset=self.first_word + j] = self.vertical_minus[j]


@always_inline
def stagger[G: Staggered](mut group: G):
    """One group through its tile: the stagger fills with the lower words masked, runs, then drains.

    Steps `first - width ..< end - 1` take every word through every column; only between
    `first - 1` and `end - width` is every word inside the tile.
    """
    var first = group.first_column()
    var end = group.end_column()
    var width = group.width()
    for offset in range(first - width, first - 1):
        group.step[True](offset)
    for offset in range(first - 1, end - width):
        group.step[False](offset)
    for offset in range(max(end - width, first - 1), end - 1):
        group.step[True](offset)


@always_inline
def masked_step[G: Staggered](mut group: G, offset: Int):
    """A masked step, or none when the group has not started or has already finished."""
    if offset >= group.first_column() - group.width() and offset < group.end_column() - 1:
        group.step[True](offset)


@always_inline
def stagger_pair[A: Staggered, B: Staggered](mut first: A, mut second: B, lag: Int):
    """Two groups through the same tile in one loop, `second` running `lag` steps behind `first`.

    Two chains in flight hide each other's latency. The groups may be independent, or `second`
    may be the words just under `first`, as long as `lag` lets every difference `first` sends down
    be stored before `second` reads it: at least `second`'s width.
    """
    var first_column = first.first_column()
    var end_column = first.end_column()
    var start = min(first_column - first.width(), first_column - second.width() + lag)
    var stop = max(end_column - 1, end_column - 1 + lag)
    var steady = first_column - 1 + lag
    var steady_end = min(end_column - first.width(), end_column - second.width() + lag)
    if steady_end < steady:
        steady_end = steady
    for offset in range(start, steady):
        masked_step(first, offset)
        masked_step(second, offset - lag)
    for offset in range(steady, steady_end):
        first.step[False](offset)
        second.step[False](offset - lag)
    for offset in range(steady_end, stop):
        masked_step(first, offset)
        masked_step(second, offset - lag)


def tile_bounds(columns: Int, width: Int) -> List[Int]:
    """Where each column tile starts, plus `columns`; the last tile absorbs a sliver too thin to stagger."""
    var bounds: List[Int] = [0]
    while bounds[len(bounds) - 1] < columns:
        var end = min(bounds[len(bounds) - 1] + width, columns)
        if columns - end < 2 * LANES:
            end = columns
        bounds.append(end)
    return bounds^


@always_inline
def word_value(plus: UInt64, minus: UInt64) -> Int:
    """How much the score grows down a word's 64 rows, from its `+` and `-` masks."""
    return Int(pop_count(plus)) - Int(pop_count(minus))


def all_bases(text: String) -> Bool:
    """Whether every byte is `A`, `C`, `G` or `T`, sixteen at a time."""
    comptime CHUNK = 16
    var bytes = text.unsafe_ptr()
    var length = text.byte_length()
    var other = SIMD[DType.bool, CHUNK](fill=False)
    var index = 0
    while index + CHUNK <= length:
        var chunk = bytes.unsafe_offset(index).unsafe_load[width=CHUNK]()
        other |= ~(
            chunk.eq(UInt8(ord("A")))
            | chunk.eq(UInt8(ord("C")))
            | chunk.eq(UInt8(ord("G")))
            | chunk.eq(UInt8(ord("T")))
        )
        index += CHUNK
    var found = other.reduce_or()
    while index < length:
        var byte = bytes[unsafe_offset=index]
        found |= not (
            byte == UInt8(ord("A")) or byte == UInt8(ord("C")) or byte == UInt8(ord("G")) or byte == UInt8(ord("T"))
        )
        index += 1
    return not found


struct Profile(Movable):
    """Both sequences as codes, and once a band needs them, as the bit planes `Sweep` reads; see
    `Sweep` for the encoding."""

    var columns: Int
    var rows: Int
    var words: Int
    var column_low: List[UInt64]
    var column_high: List[UInt64]
    var row_low: List[UInt64]
    var row_high: List[UInt64]
    var column_codes: List[UInt8]
    """The first sequence as codes, which the traceback compares base by base."""
    var row_codes: List[UInt8]
    """The second sequence as codes."""

    def __init__(out self, first: String, second: String, reverse: Bool = False) raises AlignmentError:
        """Both sequences as bit planes, or, with `reverse`, both read back to front.

        Reversed, a sweep from the start of this profile is a sweep from the end of the original,
        with no reversed copy of either string.
        """
        self.columns = first.byte_length()
        self.rows = second.byte_length()
        self.words = ceildiv(self.rows, WORD_BITS)
        if not all_bases(first) or not all_bases(second):
            for byte in first.as_bytes():
                _ = base_code(byte)
            for byte in second.as_bytes():
                _ = base_code(byte)

        # Every byte is now `A`, `C`, `G` or `T`, whose ASCII bits give the code directly: bit 2 is
        # set for `G` and `T`, the code's high bit, and bit 1 differs from bit 2 for `C` and `T`, its
        # low bit. Only the codes are built here; the planes wait for a band (see `build_planes`).
        self.column_low = List[UInt64]()
        self.column_high = List[UInt64]()
        self.row_low = List[UInt64]()
        self.row_high = List[UInt64]()
        self.column_codes = List[UInt8](capacity=self.columns + CODE_PADDING)
        self.column_codes.resize(unsafe_uninit_length=self.columns)
        var first_bytes = first.unsafe_ptr()
        var column_codes = self.column_codes.unsafe_ptr()
        comptime CHUNK = 16
        var column = 0
        while column + CHUNK <= self.columns:
            var bytes: SIMD[DType.uint8, CHUNK]
            if reverse:
                bytes = first_bytes.unsafe_offset(self.columns - CHUNK - column).unsafe_load[width=CHUNK]().reversed()
            else:
                bytes = first_bytes.unsafe_offset(column).unsafe_load[width=CHUNK]()
            column_codes.unsafe_offset(column).unsafe_store((((bytes >> 1) ^ (bytes >> 2)) & 1) | ((bytes >> 1) & 2))
            column += CHUNK
        while column < self.columns:
            var byte = first_bytes[unsafe_offset=self.columns - 1 - column if reverse else column]
            column_codes[unsafe_offset=column] = (((byte >> 1) ^ (byte >> 2)) & 1) | ((byte >> 1) & 2)
            column += 1

        comptime ONES = UInt64(0x0101010101010101)
        var second_bytes = second.unsafe_ptr()
        self.row_codes = List[UInt8](capacity=self.rows + CODE_PADDING)
        self.row_codes.resize(unsafe_uninit_length=self.rows)
        var row_codes = self.row_codes.unsafe_ptr()
        var row = 0
        while row + 8 <= self.rows:
            var eight: UInt64
            if reverse:
                eight = byte_swap(
                    second_bytes.unsafe_offset(self.rows - 8 - row).unsafe_bitcast[UInt64]().unsafe_load()
                )
            else:
                eight = second_bytes.unsafe_offset(row).unsafe_bitcast[UInt64]().unsafe_load()
            var low_bits = ((eight >> 1) ^ (eight >> 2)) & ONES
            var high_bits = (eight >> 2) & ONES
            row_codes.unsafe_offset(row).unsafe_bitcast[UInt64]().unsafe_store(low_bits | (high_bits << 1))
            row += 8
        while row < self.rows:
            var byte = second_bytes[unsafe_offset=self.rows - 1 - row if reverse else row]
            row_codes[unsafe_offset=row] = (((byte >> 1) ^ (byte >> 2)) & 1) | ((byte >> 1) & 2)
            row += 1
        # Past the last base of each, sentinels that match nothing, the two of them distinct, so a
        # match extension stops at the matrix's edge without checking it (see `slide_forward`).
        for _ in range(CODE_PADDING):
            self.column_codes.append(FIRST_SENTINEL)
            self.row_codes.append(SECOND_SENTINEL)

    def build_planes(mut self):
        """The bit planes a band sweeps, from the codes, once; a pair the diagonal transition settles
        never pays for them.

        The column planes run `COLUMN_PADDING` columns past either end, zero, for the lanes a
        staggered block holds outside its tile (see `Sweep.block`). The row planes take eight bases
        at a time, each byte's bit packed into the plane by a multiply whose partial products never
        overlap, stored negated so a row matches a column where both planes XOR to ones.
        """
        if len(self.column_low) > 0 or self.columns == 0:
            return
        var padded = self.columns + 2 * COLUMN_PADDING
        self.column_low = List[UInt64](length=padded, fill=0)
        self.column_high = List[UInt64](length=padded, fill=0)
        var low = self.column_low.unsafe_ptr().unsafe_offset(COLUMN_PADDING)
        var high = self.column_high.unsafe_ptr().unsafe_offset(COLUMN_PADDING)
        var column_codes = self.column_codes.unsafe_ptr()
        comptime CHUNK = 16
        var column = 0
        while column + CHUNK <= self.columns:
            var codes = column_codes.unsafe_offset(column).unsafe_load[width=CHUNK]()
            low.unsafe_offset(column).unsafe_store(UInt64(0) - (codes & 1).cast[DType.uint64]())
            high.unsafe_offset(column).unsafe_store(UInt64(0) - (codes >> 1).cast[DType.uint64]())
            column += CHUNK
        while column < self.columns:
            var code = column_codes[unsafe_offset=column]
            low[unsafe_offset=column] = UInt64(0) - (code & 1).cast[DType.uint64]()
            high[unsafe_offset=column] = UInt64(0) - (code >> 1).cast[DType.uint64]()
            column += 1

        comptime ONES = UInt64(0x0101010101010101)
        comptime GATHER = UInt64(0x0102040810204080)
        self.row_low = List[UInt64](length=self.words, fill=0)
        self.row_high = List[UInt64](length=self.words, fill=0)
        var row_codes = self.row_codes.unsafe_ptr()
        var row = 0
        while row + 8 <= self.rows:
            var eight = row_codes.unsafe_offset(row).unsafe_bitcast[UInt64]().unsafe_load()
            var shift = UInt64(row % WORD_BITS)
            self.row_low[row // WORD_BITS] |= ((((eight & ONES) * GATHER) >> 56) ^ 0xFF) << shift
            self.row_high[row // WORD_BITS] |= (((((eight >> 1) & ONES) * GATHER) >> 56) ^ 0xFF) << shift
            row += 8
        while row < self.rows:
            var code = row_codes[unsafe_offset=row]
            var shift = UInt64(row % WORD_BITS)
            self.row_low[row // WORD_BITS] |= ((code & 1) ^ 1).cast[DType.uint64]() << shift
            self.row_high[row // WORD_BITS] |= ((code >> 1) ^ 1).cast[DType.uint64]() << shift
            row += 1

    @always_inline
    def matches(self, column: Int, word: Int) -> UInt64:
        """Which rows of a word hold the same base as a column."""
        var padded = column + COLUMN_PADDING
        return (self.column_low[padded] ^ self.row_low[word]) & (self.column_high[padded] ^ self.row_high[word])


struct Frontier(Movable):
    """The two edges still being computed, starting from the global borders.

    Each step along the top row or down the left column costs one, so every horizontal difference
    starts at `+1` and every vertical word at all `+1`. A banded sweep relies on both: a word that
    first enters the band still reads `+1` down its left edge, and the top word of a tile still reads
    `+1` from above, and each is the cost of a real path, so never below the true score.
    """

    var horizontal_plus: List[UInt64]
    var horizontal_minus: List[UInt64]
    var vertical_plus: List[UInt64]
    var vertical_minus: List[UInt64]

    def __init__(out self, columns: Int, words: Int):
        self.horizontal_plus = List[UInt64](length=columns, fill=1)
        self.horizontal_minus = List[UInt64](length=columns, fill=0)
        self.vertical_plus = List[UInt64](length=words, fill=ALL_ONES)
        self.vertical_minus = List[UInt64](length=words, fill=0)

    def restart_horizontal(mut self, count: Int):
        """The first `count` horizontal differences back to the border's `+1`, a vector at a time."""
        var plus = self.horizontal_plus.unsafe_ptr()
        var minus = self.horizontal_minus.unsafe_ptr()
        var index = 0
        while index + LANES <= count:
            plus.unsafe_offset(index).unsafe_store(Words(1))
            minus.unsafe_offset(index).unsafe_store(Words(0))
            index += LANES
        while index < count:
            plus[unsafe_offset=index] = 1
            minus[unsafe_offset=index] = 0
            index += 1

    def sweep(mut self, mut profile: Profile) -> Sweep:
        """Pointers into both, valid while both outlive the sweep."""
        return Sweep(
            profile.column_low,
            profile.column_high,
            profile.row_low,
            profile.row_high,
            self.horizontal_plus,
            self.horizontal_minus,
            self.vertical_plus,
            self.vertical_minus,
        )

    def down_from(self, word: Int, rows: Int) -> Int:
        """How much the score grows from the top of `word` to the bottom row, down the right edge.

        Rows past the end of the last word only padded it and are masked out.
        """
        var words = len(self.vertical_plus)
        var total = 0
        for index in range(word, words):
            var plus = self.vertical_plus[index]
            var minus = self.vertical_minus[index]
            if index == words - 1 and rows % WORD_BITS != 0:
                var kept = (UInt64(1) << UInt64(rows % WORD_BITS)) - 1
                plus &= kept
                minus &= kept
            total += word_value(plus, minus)
        return total


def full_distance(mut profile: Profile, workers: Int) -> Int:
    """The whole matrix, tiled across `workers` threads once it is large enough to pay for them."""
    profile.build_planes()
    var frontier = Frontier(profile.columns, profile.words)
    var sweep = frontier.sweep(profile)
    var words = profile.words
    if workers == 1 or profile.columns * profile.rows < PARALLEL_CELLS:
        sweep.words(0, words, 0, profile.columns)
    else:
        # Tiles of `LANES` words by a column tile, swept one anti-diagonal of tiles at a time; the
        # words left over below the last full block form one more row of tiles.
        var bounds = tile_bounds(profile.columns, TILE_COLUMNS)
        var column_tiles = len(bounds) - 1
        var row_tiles = ceildiv(words, LANES)
        for diagonal in range(row_tiles + column_tiles - 1):
            var first_row_tile = max(0, diagonal - column_tiles + 1)
            var last_row_tile = min(diagonal, row_tiles - 1)

            def tile(slot: Int) {imm}:
                var row_tile = first_row_tile + slot
                var column_tile = diagonal - row_tile
                sweep.words(
                    row_tile * LANES,
                    min(row_tile * LANES + LANES, words),
                    bounds[column_tile],
                    bounds[column_tile + 1],
                )

            parallelize(tile, last_row_tile - first_row_tile + 1, workers)

    # The top-right corner is `columns`; walking down the right edge adds each vertical difference.
    return profile.columns + frontier.down_from(0, profile.rows)


struct Trail(Movable):
    """The left edge of every tile a round swept, which is all the traceback needs to retrace it.

    Per tile: its columns, the band's words, the score at the band's top on the left edge, and the
    vertical differences down that edge. A few words per tile rather than every column, so recording
    costs a copy per tile, not the matrix.
    """

    var first_columns: List[Int]
    var end_columns: List[Int]
    var tops: List[Int]
    var ends: List[Int]
    var anchors: List[Int]
    """The score at the left edge's top, row `tops[tile] * WORD_BITS`."""
    var offsets: List[Int]
    """Where each tile's edge words start in `edge_plus` and `edge_minus`."""
    var edge_plus: List[UInt64]
    var edge_minus: List[UInt64]

    def __init__(out self, columns: Int = 0):
        """Room for a round over `columns` columns in its narrowest tiles, a vector's words each."""
        var tiles = columns // NARROW_COLUMNS + 2
        self.first_columns = List[Int](capacity=tiles)
        self.end_columns = List[Int](capacity=tiles)
        self.tops = List[Int](capacity=tiles)
        self.ends = List[Int](capacity=tiles)
        self.anchors = List[Int](capacity=tiles)
        self.offsets = List[Int](capacity=tiles)
        self.edge_plus = List[UInt64](capacity=tiles * LANES)
        self.edge_minus = List[UInt64](capacity=tiles * LANES)

    def clear(mut self):
        self.first_columns.clear()
        self.end_columns.clear()
        self.tops.clear()
        self.ends.clear()
        self.anchors.clear()
        self.offsets.clear()
        self.edge_plus.clear()
        self.edge_minus.clear()

    def record(mut self, first_column: Int, end_column: Int, top: Int, end: Int, anchor: Int, frontier: Frontier):
        """A tile's left edge, read from the frontier just before the tile is swept."""
        self.first_columns.append(first_column)
        self.end_columns.append(end_column)
        self.tops.append(top)
        self.ends.append(end)
        self.anchors.append(anchor)
        self.offsets.append(len(self.edge_plus))
        for word in range(top, end):
            self.edge_plus.append(frontier.vertical_plus[word])
            self.edge_minus.append(frontier.vertical_minus[word])


@fieldwise_init
struct Round(ImplicitlyCopyable, TrivialRegisterPassable):
    """What one pruned round learned: a distance, or how far across the matrix any path within the bound got."""

    var distance: Int
    """An upper bound on the distance, exact when within the round's bound; -1 when no path fit it."""
    var reached: Int
    """Matrix columns crossed before every row was pruned, or all of them when a distance was found."""
    var bound: Int
    """The bound the round ended on, which its checkpoints may have lowered (see `HalfBand.check`)."""
    var estimate: Int
    """The distance a checkpoint projected when it gave the round up, or -1."""


# region Seed heuristic

comptime SEED_LENGTH = 12
"""Bases per seed: A*PA2-full's `k`, long enough that a random match is rare, short enough that most
seeds of a moderately divergent pair still match exactly."""

comptime INEXACT_LENGTH = 16
"""Bases per seed when seeds may match with one edit: A*PA's inexact seeds, long enough that a random
one-edit match is still rare, each seed then charging two edits to a path that matches it nowhere."""

comptime INEXACT_COLUMNS = 65_536
"""Columns from which the seeds may match with one edit. Their setup grows with the length and the band
they save with its square: on real reads of 16 to 70 kbp exact seeds win however divergent the pair,
and from about 100 kbp, and on uniform errors of one in seven, inexact ones do."""

comptime INEXACT_DIVERGENCE = 10
"""Columns per projected edit at or below which the seeds may match with one edit: past about one edit
in ten bases most exact seeds are broken, while most still match within one edit."""

comptime INEXACT_CHAINED = 40
"""Percent of the exact seeds chained from the origin below which the seeds are rebuilt to match within
one edit: about one edit in fifteen bases on real reads, whose errors gather, and one in fifteen to
twenty on uniform ones. Above it the exact seeds' bound already lies within a tenth or so of the
distance and inexact seeds narrow the band by less than their setup costs."""

comptime HALF_BITS = INEXACT_LENGTH
"""Bits in half an inexact seed's two-bit code: a one-edit match matches one half exactly, so each half
indexes a table of this many bits."""

comptime SEED_EDITS = 1500
"""Projected edits from which a band prunes with the seed heuristic. Its setup costs about ten
nanoseconds a column; what it saves grows with the distance, the band otherwise sweeping rows in
proportion to it, and passes the setup about here."""

comptime SEED_DIVERGENCE = 7
"""Columns per projected edit below which the seeds are left out: past about one edit in seven bases,
almost no seed survives local pruning, the heuristic is little more than an edit a seed, and its
setup is not repaid."""

comptime SEED_COLUMNS = 16_384
"""Columns from which a band prunes with the seeds whatever the projection: its setup is a small share
of any band this long, and on real reads, whose errors gather at the ends, the projection that gates
shorter pairs can put a divergence of one edit in ten at one in two."""

comptime SEEDED_GROWTH = 4
"""How many times a seeded band's margin over the heuristic at the origin grows after a round that
died early."""

comptime SEEDED_TRUST_SHARE = 8
"""A seeded round must cross one part in this many of the columns before its death projects the next
bound."""

comptime CHAINED_SHARE = 32
"""With seeds, the band's first bound starts at the origin's when at least one seed in this many is
chained there; with fewer, the bound is little more than an edit a seed, and the projection leads."""

comptime SEED_SLACK = 64
"""With seeds, the band's first bound reaches at least this far past the heuristic at the origin."""

comptime LOOKAHEAD_SEEDS = 14
"""Seeds a match's local pruning looks ahead, its own included: A*PA2-full's `p`."""

comptime LAYER_SLOTS = 8
"""Starts a layer holds in place before the rest spill into a chain of its own."""


struct SeedHeuristic(Movable):
    """A*PA2-full's gap-chaining seed heuristic, without match pruning.

    The first sequence is cut into disjoint seeds of `SEED_LENGTH` bases, and every exact occurrence
    of a seed in the second sequence is a match. A path that misses a seed must spend an edit on it,
    so the seeds still ahead of a cell, its potential `P`, bound the cost to the end, less one for
    each match a path can still chain. Chaining charges the gap between matches: in the coordinates
    `T(i, j) = (i - j - P(i), j - i - P(i))` one match can follow another exactly when the second's
    start lies above and right of the first's end, so the longest chain from a cell is a dominance
    query, answered from layers of match starts, layer `v` holding the starts from which `v` matches
    chain (see `score`).

    With inexact seeds, as A*PA's `r = 2`, seeds are `INEXACT_LENGTH` bases, a match may cost one
    edit, and a seed matched nowhere costs a path two: the potential counts two a seed, a match
    scores two less its cost, and the transform is the same in those units. A match scoring two
    starts in two layers, so the layers still nest.

    The heuristic never overestimates the cost to the end. Down a column it drops by at most one a
    row, as it is consistent, and up a column by at most `cost` (see `climb`), so the band's
    pruning stays exact and its jumps stay valid. With no seeds it is the plain gap to the end's
    diagonal, which is what every band used before.
    """

    var columns: Int
    var rows: Int
    var length: Int
    """Bases per seed."""
    var cost: Int
    """What a path pays crossing a seed it matches nowhere: one for exact seeds, two for inexact."""
    var seeds: Int
    var slot_x: List[Int32]
    """Per layer, `LAYER_SLOTS` slots for the transformed starts of the matches from which that many
    matches chain; a layer rarely holds more than three, and any past the slots spill over."""
    var slot_y: List[Int32]
    var counts: List[Int32]
    var spill_head: List[Int32]
    """Per layer, the last of its spilled starts, -1 for none; each spilled start links to the one
    spilled before it in `spill_next`."""
    var spill_next: List[Int32]
    var spill_x: List[Int32]
    var spill_y: List[Int32]
    var hint: Int
    """The layer the last query ended on: neighbouring queries land near it."""

    def __init__(out self, columns: Int, rows: Int):
        """No seeds: the gap heuristic."""
        self.columns = columns
        self.rows = rows
        self.length = SEED_LENGTH
        self.cost = 1
        self.seeds = 0
        self.slot_x = List[Int32]()
        self.slot_y = List[Int32]()
        self.counts = List[Int32]()
        self.spill_head = List[Int32]()
        self.spill_next = List[Int32]()
        self.spill_x = List[Int32]()
        self.spill_y = List[Int32]()
        self.hint = 0

    def __init__(
        out self, profile: Profile, inexact: Bool = False, choose: Bool = False, lookahead: Int = LOOKAHEAD_SEEDS
    ):
        """Seeds of the profile's first sequence, matched exactly in its second, or with `inexact`
        within one edit. With `choose`, exact seeds, rebuilt inexact when fewer than
        `INEXACT_CHAINED` percent of them chain from the origin.

        With `lookahead`, a match is kept only if a path from its start crosses the next `lookahead`
        seeds for less than they would cost unmatched (see `worth_keeping`).
        """
        self.columns = profile.columns
        self.rows = profile.rows
        self.length = SEED_LENGTH
        self.cost = 1
        self.seeds = 0
        self.slot_x = List[Int32]()
        self.slot_y = List[Int32]()
        self.counts = List[Int32]()
        self.spill_head = List[Int32]()
        self.spill_next = List[Int32]()
        self.spill_x = List[Int32]()
        self.spill_y = List[Int32]()
        self.hint = 0
        self.build(profile, inexact, lookahead)
        if choose and not inexact and self.seeds > 0:
            var chained = self.seeds - self.h(0, 0)
            if chained * 100 < INEXACT_CHAINED * self.seeds:
                self.build(profile, True, lookahead)

    def build(mut self, profile: Profile, inexact: Bool, lookahead: Int):
        """The seeds, their matches, and the layers, from scratch."""
        self.length = INEXACT_LENGTH if inexact else SEED_LENGTH
        self.cost = 2 if inexact else 1
        self.seeds = profile.columns // self.length
        self.slot_x.clear()
        self.slot_y.clear()
        self.counts.clear()
        self.spill_head.clear()
        self.spill_next.clear()
        self.spill_x.clear()
        self.spill_y.clear()
        self.slot_x.reserve(LAYER_SLOTS * (self.seeds + 1))
        self.slot_y.reserve(LAYER_SLOTS * (self.seeds + 1))
        self.counts.reserve(self.seeds + 1)
        self.spill_head.reserve(self.seeds + 1)
        self.hint = 0
        self.add_sentinel()
        if self.seeds == 0 or profile.rows < self.length + 1:
            return
        var first = profile.column_codes.unsafe_ptr()
        var second = profile.row_codes.unsafe_ptr()
        # Each match a seed and its start row, and for inexact seeds its end row and cost; an exact
        # match ends a seed's length further down, at no cost.
        var found_seed = List[Int32]()
        var found_row = List[Int32]()
        var found_end = List[Int32]()
        var found_cost = List[Int32]()
        if inexact:
            self.inexact_matches(first, second, found_seed, found_row, found_end, found_cost)
        else:
            self.exact_matches(first, second, found_seed, found_row)

        # Bucketed by seed: a start can dominate another match's end only from a later seed, so taking
        # the seeds last first is an order the layers can be built in, no sort.
        var firsts = List[Int32](length=self.seeds + 1, fill=0)
        for index in range(len(found_seed)):
            firsts[Int(found_seed[index]) + 1] += 1
        for seed in range(self.seeds):
            firsts[seed + 1] += firsts[seed]
        var rows_by_seed = List[Int32](length=len(found_seed), fill=0)
        var ends_by_seed = List[Int32](length=len(found_end), fill=0)
        var costs_by_seed = List[Int32](length=len(found_cost), fill=0)
        var filled = firsts.copy()
        for index in range(len(found_seed)):
            var seed = Int(found_seed[index])
            var slot = Int(filled[seed])
            rows_by_seed[slot] = found_row[index]
            if inexact:
                ends_by_seed[slot] = found_end[index]
                costs_by_seed[slot] = found_cost[index]
            filled[seed] += 1

        # A match starts `cost` less its own cost layers above the best its end can chain on to. Right
        # to left on every diagonal, as the seeds go last first, so local pruning sees the matches
        # kept after it.
        var leftmost = List[Int32](length=self.columns + self.rows + 1, fill=Int32.MAX)
        var fronts = List[Int](length=4 * self.cost * LOOKAHEAD_SEEDS + 3, fill=0)
        var spare = List[Int](length=4 * self.cost * LOOKAHEAD_SEEDS + 3, fill=0)
        for seed in range(self.seeds - 1, -1, -1):
            var column = seed * self.length
            var potential = self.potential(column)
            var end_potential = potential - self.cost
            for slot in range(Int(firsts[seed]), Int(firsts[seed + 1])):
                var start_row = Int(rows_by_seed[slot])
                var end_row = Int(ends_by_seed[slot]) if inexact else start_row + SEED_LENGTH
                var match_cost = Int(costs_by_seed[slot]) if inexact else 0
                if lookahead > 0 and not self.worth_keeping(
                    first, second, seed, start_row, end_row, match_cost, lookahead, leftmost, fronts, spare
                ):
                    continue
                leftmost[column - start_row + self.rows] = Int32(column)
                var x = column - start_row - potential
                var y = start_row - column - potential
                var end_x = column + self.length - end_row - end_potential
                var end_y = end_row - column - self.length - end_potential
                var score = self.cost - match_cost
                var layer = self.score(end_x, end_y) + score
                while layer >= len(self.counts):
                    self.add_layer()
                for below in range(score):
                    self.add_point(layer - below, x, y)

    def exact_matches(
        self,
        first: ImmPointer[UInt8, _],
        second: ImmPointer[UInt8, _],
        mut found_seed: List[Int32],
        mut found_row: List[Int32],
    ):
        """Every exact occurrence of a seed whose chain can still reach the end.

        Every seed's two-bit code is hashed by open addressing on the multiply's top bits; a slot
        holds its code and the first seed with it in one word, and seeds sharing a code chain on.
        """
        comptime MASK = (1 << (2 * SEED_LENGTH)) - 1
        comptime EMPTY = Int64(-1)
        var bits = 1
        while (1 << bits) < 2 * self.seeds:
            bits += 1
        var size = 1 << bits
        var table = List[Int64](length=size, fill=EMPTY)
        var next = List[Int32](length=self.seeds, fill=-1)
        var slots = table.unsafe_ptr()
        for seed in range(self.seeds):
            var code = 0
            for offset in range(SEED_LENGTH):
                code = (code << 2) | Int(first[unsafe_offset=seed * SEED_LENGTH + offset])
            var slot = Int((UInt64(code) * 0x9E3779B97F4A7C15) >> UInt64(64 - bits))
            while slots[unsafe_offset=slot] != EMPTY and Int(slots[unsafe_offset=slot] >> 32) != code:
                slot = (slot + 1) & (size - 1)
            var held = slots[unsafe_offset=slot]
            next[seed] = Int32(held & 0xFFFFFFFF) if held != EMPTY else -1
            slots[unsafe_offset=slot] = (Int64(code) << 32) | Int64(seed)

        # Every window of the second sequence looked up. A match's end, one seed on, is `(x + 1, y + 1)`
        # in transformed coordinates, and must lie at or below and left of the end's.
        var target_x = self.columns - self.rows
        var target_y = self.rows - self.columns
        var code = 0
        for row in range(self.rows):
            code = ((code << 2) | Int(second[unsafe_offset=row])) & MASK
            if row + 1 < SEED_LENGTH:
                continue
            var slot = Int((UInt64(code) * 0x9E3779B97F4A7C15) >> UInt64(64 - bits))
            var held = slots[unsafe_offset=slot]
            while held != EMPTY and Int(held >> 32) != code:
                slot = (slot + 1) & (size - 1)
                held = slots[unsafe_offset=slot]
            if held == EMPTY:
                continue
            var start_row = row + 1 - SEED_LENGTH
            var seed = Int(held & 0xFFFFFFFF)
            while seed >= 0:
                var column = seed * SEED_LENGTH
                var potential = self.seeds - seed
                if column - start_row - potential + 1 <= target_x and start_row - column - potential + 1 <= target_y:
                    found_seed.append(Int32(seed))
                    found_row.append(Int32(start_row))
                seed = Int(next[seed])

    def inexact_matches(
        self,
        first: ImmPointer[UInt8, _],
        second: ImmPointer[UInt8, _],
        mut found_seed: List[Int32],
        mut found_row: List[Int32],
        mut found_end: List[Int32],
        mut found_cost: List[Int32],
    ):
        """Every occurrence of a seed within one edit whose chain can still reach the end.

        One edit leaves one half of the seed matching exactly, at the window's start or its end. Each
        half's code indexes a table of the seeds with it, and every window of the second sequence is
        looked up as either half: a left half found at a row tries the windows starting there, one
        base shorter, as long, or one longer than the seed; a right half found tries those ending
        where it ends. A window both halves find is taken from the left half alone.
        """
        comptime HALF = INEXACT_LENGTH // 2
        comptime HALF_MASK = (1 << HALF_BITS) - 1
        comptime BUCKETS = 1 << HALF_BITS
        # The quarters below are a code's top, second, third and bottom bytes.
        comptime assert INEXACT_LENGTH == 16, "an inexact seed's quarters are its code's bytes"
        var codes = List[UInt64](length=self.seeds, fill=0)
        # Each half's seeds bucketed by its code, contiguous, the codes beside them.
        var left_start = List[Int32](length=BUCKETS + 1, fill=0)
        var right_start = List[Int32](length=BUCKETS + 1, fill=0)
        for seed in range(self.seeds):
            var code = UInt64(0)
            for offset in range(INEXACT_LENGTH):
                code = (code << 2) | UInt64(first[unsafe_offset=seed * INEXACT_LENGTH + offset])
            codes[seed] = code
            left_start[Int(code >> UInt64(HALF_BITS)) + 1] += 1
            right_start[(Int(code) & HALF_MASK) + 1] += 1
        for bucket in range(BUCKETS):
            left_start[bucket + 1] += left_start[bucket]
            right_start[bucket + 1] += right_start[bucket]
        var left_seeds = List[Int32](length=self.seeds, fill=0)
        var left_codes = List[UInt64](length=self.seeds, fill=0)
        var right_seeds = List[Int32](length=self.seeds, fill=0)
        var right_codes = List[UInt64](length=self.seeds, fill=0)
        var left_fill = left_start.copy()
        var right_fill = right_start.copy()
        for seed in range(self.seeds):
            var code = codes[seed]
            var left = Int(left_fill[Int(code >> UInt64(HALF_BITS))])
            left_seeds[left] = Int32(seed)
            left_codes[left] = code
            left_fill[Int(code >> UInt64(HALF_BITS))] += 1
            var right = Int(right_fill[Int(code) & HALF_MASK])
            right_seeds[right] = Int32(seed)
            right_codes[right] = code
            right_fill[Int(code) & HALF_MASK] += 1

        # Every window of `INEXACT_LENGTH` bases of the second sequence as one code, the first base in
        # the high bits, and its first four bases alone; past the end, the bases read as zero, and no
        # window reaching there is tried.
        var windows = List[UInt64](length=self.rows + 1, fill=0)
        var quarters = List[UInt8](length=self.rows + 1, fill=0)
        var rolling = UInt64(0)
        for row in range(self.rows + INEXACT_LENGTH - 1, -1, -1):
            var base = UInt64(second[unsafe_offset=row]) if row < self.rows else UInt64(0)
            rolling = (rolling >> 2) | (base << UInt64(2 * INEXACT_LENGTH - 2))
            if row <= self.rows:
                windows[row] = rolling
                quarters[row] = UInt8(rolling >> UInt64(2 * INEXACT_LENGTH - 8))
        var window = windows.unsafe_ptr()
        var quarter = quarters.unsafe_ptr()
        var last = self.rows
        for row in range(self.rows - HALF + 1):
            # A left half found here leaves the right half within one edit of the bases after it: the
            # seed's third quarter matches at `row + 8`, or its last quarter at `row + 11`, `row + 12`
            # or `row + 13`, as the edit falls after or before the last quarter. A right half ending
            # at `end` likewise puts the second quarter at `end - 12` or the first at `end - 15`,
            # `end - 16` or `end - 17`. Most chance lookups miss all four.
            var half = Int(window[unsafe_offset=row] >> UInt64(HALF_BITS))
            var third = quarter[unsafe_offset=min(row + 8, last)]
            var at_11 = quarter[unsafe_offset=min(row + 11, last)]
            var at_12 = quarter[unsafe_offset=min(row + 12, last)]
            var at_13 = quarter[unsafe_offset=min(row + 13, last)]
            for slot in range(Int(left_start[half]), Int(left_start[half + 1])):
                var code = left_codes[slot]
                var fourth = UInt8(code & 0xFF)
                if UInt8((code >> 8) & 0xFF) != third and fourth != at_11 and fourth != at_12 and fourth != at_13:
                    continue
                self.try_windows(
                    code, Int(left_seeds[slot]), row, -1, window, found_seed, found_row, found_end, found_cost
                )
            var end = row + HALF
            var second_at = quarter[unsafe_offset=max(end - 12, 0)]
            var at_15 = quarter[unsafe_offset=max(end - 15, 0)]
            var at_16 = quarter[unsafe_offset=max(end - 16, 0)]
            var at_17 = quarter[unsafe_offset=max(end - 17, 0)]
            for slot in range(Int(right_start[half]), Int(right_start[half + 1])):
                var code = right_codes[slot]
                var opening = UInt8(code >> 24)
                if (
                    UInt8((code >> 16) & 0xFF) != second_at
                    and opening != at_15
                    and opening != at_16
                    and opening != at_17
                ):
                    continue
                self.try_windows(
                    code, Int(right_seeds[slot]), -1, end, window, found_seed, found_row, found_end, found_cost
                )

    @always_inline
    def try_windows(
        self,
        code: UInt64,
        seed: Int,
        start: Int,
        end: Int,
        window: ImmPointer[UInt64, _],
        mut found_seed: List[Int32],
        mut found_row: List[Int32],
        mut found_end: List[Int32],
        mut found_cost: List[Int32],
    ):
        """The windows within one edit of a seed that start at `start`, or end at `end` when `start`
        is -1: one base shorter than the seed, as long, and one longer.

        Two codes of equal length agree base by base where their XOR has neither bit of a pair set;
        the leading agreeing bases are the common prefix, the trailing ones the common suffix. A
        window one base shorter is the seed less one base when the two cover it between them, and
        one longer is the seed plus one base when they cover the seed.

        An exact match's windows one base shorter or longer, sharing its start or its end, are left
        out: a path crossing the seed along one costs at least one edit there, and taking the exact
        match in its place moves the chain's gaps by one edit at most, so the heuristic stays a
        lower bound without them.
        """
        comptime K = INEXACT_LENGTH
        comptime PAIRS = UInt64(0x5555555555555555)
        var exact = False
        for step in range(3):
            # As long as the seed first, so an exact match is known before its neighbours.
            var extra = 0 if step == 0 else (-1 if step == 1 else 1)
            if exact and extra != 0:
                return
            var size = K + extra
            var low = start if start >= 0 else end - size
            var high = low + size
            if low < 0 or high > self.rows:
                continue
            # A window whose left half matches too is the left half's lookup's to take.
            var taken = start < 0 and (window[unsafe_offset=low] >> UInt64(HALF_BITS)) == (code >> UInt64(HALF_BITS))
            var cost = 1
            if extra == 0:
                var mismatched = window[unsafe_offset=low] ^ code
                var bases = Int(pop_count((mismatched | (mismatched >> 1)) & PAIRS))
                exact = bases == 0
                if bases > 1 or taken:
                    continue
                cost = bases
            elif taken:
                continue
            elif extra == -1:
                # The window's bases, `K - 1` of them, against the seed's first and last `K - 1`.
                var shorter = window[unsafe_offset=low] >> 2
                var prefix = self.agreeing_prefix(shorter, code >> 2, K - 1)
                var suffix = self.agreeing_suffix(shorter, code & ((UInt64(1) << UInt64(2 * K - 2)) - 1), K - 1)
                if prefix + suffix < K - 1:
                    continue
            else:
                # The window's first and last `K` bases against the seed.
                var head = window[unsafe_offset=low]
                var tail = window[unsafe_offset=low + 1]
                var prefix = self.agreeing_prefix(head, code, K)
                var suffix = self.agreeing_suffix(tail, code, K)
                if prefix + suffix < K:
                    continue
            self.add_match(seed, low, high, cost, found_seed, found_row, found_end, found_cost)

    @always_inline
    def agreeing_prefix(self, first: UInt64, second: UInt64, bases: Int) -> Int:
        """Leading bases two codes of `bases` bases share."""
        comptime PAIRS = UInt64(0x5555555555555555)
        var mismatched = first ^ second
        var differing = (mismatched | (mismatched >> 1)) & PAIRS
        if differing == 0:
            return bases
        return (Int(count_leading_zeros(differing)) - (64 - 2 * bases)) // 2

    @always_inline
    def agreeing_suffix(self, first: UInt64, second: UInt64, bases: Int) -> Int:
        """Trailing bases two codes of `bases` bases share."""
        comptime PAIRS = UInt64(0x5555555555555555)
        var mismatched = first ^ second
        var differing = (mismatched | (mismatched >> 1)) & PAIRS
        if differing == 0:
            return bases
        return Int(count_trailing_zeros(differing)) // 2

    @always_inline
    def add_match(
        self,
        seed: Int,
        start_row: Int,
        end_row: Int,
        cost: Int,
        mut found_seed: List[Int32],
        mut found_row: List[Int32],
        mut found_end: List[Int32],
        mut found_cost: List[Int32],
    ):
        """Keeps an inexact match whose end can still reach the end: in transformed coordinates, its
        end lies at or below and left of the end's."""
        var column = seed * self.length
        var end_potential = self.potential(column) - self.cost
        var end_x = column + self.length - end_row - end_potential
        var end_y = end_row - column - self.length - end_potential
        if end_x <= self.columns - self.rows and end_y <= self.rows - self.columns:
            found_seed.append(Int32(seed))
            found_row.append(Int32(start_row))
            found_end.append(Int32(end_row))
            found_cost.append(Int32(cost))

    def worth_keeping(
        self,
        first: ImmPointer[UInt8, _],
        second: ImmPointer[UInt8, _],
        seed: Int,
        start_row: Int,
        end_row: Int,
        match_cost: Int,
        lookahead: Int,
        leftmost: List[Int32],
        mut fronts: List[Int],
        mut spare: List[Int],
    ) -> Bool:
        """A*PA's local pruning: whether a match can lower the cost of some path, so dropping it would not.

        A match only helps a path that goes on to cross the next seeds for fewer edits than the
        seeds it crosses would cost unmatched. A diagonal transition from the match's end, starting
        at the match's own cost, over the next `lookahead` seeds counting its own, looks for such a
        path: the match stays if a front reaches past the last of them, or slides into a match
        already kept on its diagonal, which then continues it.
        A front whose edits already equal the seeds crossed can no longer gain, and is dropped. Every
        match an optimal path's chain relies on passes, so the heuristic stays a lower bound.
        """
        var start_column = seed * self.length
        var start_potential = self.potential(start_column)
        var last = min(seed + lookahead - 1, self.seeds - 1)
        var end_column = (last + 1) * self.length
        var reach = start_potential - self.potential(end_column)
        var origin = (start_column + self.length) - end_row
        # Front `d` is the diagonal `origin + d - reach`, at `fronts[d]`, its furthest column.
        var low = reach
        var high = reach + 1
        fronts[reach] = start_column + self.length
        fronts[reach] = extend(first, second, fronts[reach], end_row, self.columns, self.rows)
        if fronts[reach] >= end_column:
            return True
        var kept = Int(leftmost[origin + self.rows])
        if kept <= fronts[reach]:
            return True
        for cost in range(match_cost + 1, reach):
            # One more edit: from the same diagonal, or from either neighbour.
            for d in range(low - 1, high + 1):
                var best = -1
                if d >= low and d < high:
                    best = fronts[d] + 1
                if d + 1 >= low and d + 1 < high:
                    best = max(best, fronts[d + 1])
                if d - 1 >= low and d - 1 < high:
                    best = max(best, fronts[d - 1] + 1)
                spare[d] = min(best, self.columns)
            for d in range(low - 1, high + 1):
                fronts[d] = spare[d]
            low -= 1
            high += 1
            # A front whose edits match the seeds it has crossed can no longer gain.
            while low < high and cost + self.potential(fronts[low]) >= start_potential:
                low += 1
            while high > low and cost + self.potential(fronts[high - 1]) >= start_potential:
                high -= 1
            if low == high:
                return False
            for d in range(low, high):
                var diagonal = origin + d - reach
                var before = fronts[d]
                var row = before - diagonal
                if row < 0 or row > self.rows:
                    continue
                fronts[d] = extend(first, second, before, row, self.columns, self.rows)
                if fronts[d] >= end_column:
                    return True
                var next = Int(
                    leftmost[diagonal + self.rows]
                ) if diagonal + self.rows >= 0 and diagonal + self.rows < len(leftmost) else Int(Int32.MAX)
                if before <= next and next <= fronts[d]:
                    return True
        return False

    def add_layer(mut self):
        for _ in range(LAYER_SLOTS):
            self.slot_x.append(0)
            self.slot_y.append(0)
        self.counts.append(0)
        self.spill_head.append(-1)

    def add_point(mut self, layer: Int, x: Int, y: Int):
        var count = Int(self.counts[layer])
        if count < LAYER_SLOTS:
            self.slot_x[layer * LAYER_SLOTS + count] = Int32(x)
            self.slot_y[layer * LAYER_SLOTS + count] = Int32(y)
            self.counts[layer] = Int32(count + 1)
        else:
            self.spill_next.append(self.spill_head[layer])
            self.spill_head[layer] = Int32(len(self.spill_x))
            self.spill_x.append(Int32(x))
            self.spill_y.append(Int32(y))

    def add_sentinel(mut self):
        """Layer zero: a point dominating everything, as no match is chained."""
        self.add_layer()
        self.add_point(0, Int(Int32.MAX), Int(Int32.MAX))

    @always_inline
    def contains(self, layer: Int, x: Int, y: Int) -> Bool:
        """Whether a start in `layer` lies at or above and right of `(x, y)`."""
        var count = Int(self.counts[layer])
        var xs = self.slot_x.unsafe_ptr().unsafe_offset(layer * LAYER_SLOTS)
        var ys = self.slot_y.unsafe_ptr().unsafe_offset(layer * LAYER_SLOTS)
        for index in range(count):
            if x <= Int(xs[unsafe_offset=index]) and y <= Int(ys[unsafe_offset=index]):
                return True
        if count == LAYER_SLOTS:
            var index = Int(self.spill_head[layer])
            while index >= 0:
                if x <= Int(self.spill_x[index]) and y <= Int(self.spill_y[index]):
                    return True
                index = Int(self.spill_next[index])
        return False

    def score(mut self, x: Int, y: Int) -> Int:
        """The most matches a chain from transformed point `(x, y)` takes.

        Layers nest, every start of layer `v + 1` lying below some start of layer `v`, so the
        layers containing a point are a prefix: a search galloping out from the last answer, as
        neighbouring queries land on neighbouring layers, finds its end.
        """
        var last = len(self.counts) - 1
        var guess = min(self.hint, last)
        var low: Int
        var high: Int
        if self.contains(guess, x, y):
            # The end lies at or above the guess.
            low = guess
            var step = 1
            high = min(guess + step, last)
            while high > low and self.contains(high, x, y):
                low = high
                step *= 2
                high = min(low + step, last)
            if high == low:
                self.hint = low
                return low
            high -= 1
        else:
            # The end lies below the guess; layer zero always contains.
            high = guess - 1
            var step = 1
            low = max(guess - step, 0)
            while low > 0 and not self.contains(low, x, y):
                high = low - 1
                step *= 2
                low = max(high - step, 0)
        while low < high:
            var middle = (low + high + 1) // 2
            if self.contains(middle, x, y):
                low = middle
            else:
                high = middle - 1
        self.hint = low
        return low

    def chains_well(mut self) -> Bool:
        """Whether at least one seed in `CHAINED_SHARE` is chained from the origin: then the bound
        there lies close to the distance, rather than at little more than an edit a seed."""
        if self.seeds == 0:
            return False
        var whole = self.potential(0)
        return (whole - self.h(0, 0)) * CHAINED_SHARE >= whole

    @always_inline
    def potential(self, column: Int) -> Int:
        """What the seeds starting at or after `column` cost a path matching none of them."""
        # Each length its own constant divisor: the band asks for this at every row it prunes.
        if self.cost == 1:
            return self.seeds - min(self.seeds, ceildiv(column, SEED_LENGTH))
        return 2 * (self.seeds - min(self.seeds, ceildiv(column, INEXACT_LENGTH)))

    @always_inline
    def climb(self) -> Int:
        """The most the heuristic drops a row going up a column: one with no seeds or exact ones,
        `cost` with inexact. A best chain from a row still starts, but for its first match, from
        the row below: every match moves the transformed point at least one up and one right."""
        return self.cost

    def h(mut self, column: Int, row: Int) -> Int:
        """A lower bound on the cost from `(column, row)` to the end."""
        var gap = abs((self.columns - column) - (self.rows - row))
        if self.seeds == 0:
            return gap
        var potential = self.potential(column)
        var chained = self.score(column - row - potential, row - column - potential)
        if chained == 0:
            return max(gap, potential)
        return potential - chained


# endregion Seed heuristic


struct HalfBand(Movable):
    """One round of band doubling, advanced a tile at a time, so two rounds can share a loop.

    `pruned_distance` documents the round; this holds its state between tiles: the frontier, the
    band's top and bottom words, the score at the top, and the deepest kept row and its score.
    """

    var frontier: Frontier
    var sweep: Sweep
    var bounds: List[Int]
    var columns: Int
    var rows: Int
    var words: Int
    var difference: Int
    var extra: Int
    var threshold: Int
    var stop_column: Int
    var top: Int
    var end_word: Int
    var anchor: Int
    """The score at the top of word `top`, on the left edge of the coming tile."""
    var deepest: Int
    """The deepest row a kept cell, scoring `floor` at the least, sits on, after the last tile."""
    var floor: Int
    var outcome: Round
    """Set once the round has pruned every row, the reason it stopped."""
    var adapt: Bool
    """Whether the round re-aims its bound at checkpoints (see `check`)."""
    var lower: Bool
    """Whether a checkpoint may give the round up as well as lower its bound."""
    var lowest: Int
    """The least a checkpoint may lower the bound to."""
    var checkpoint: Int
    """Checkpoints passed."""

    def __init__(
        out self,
        mut profile: Profile,
        threshold: Int,
        stop_column: Int,
        adapt: Bool = False,
        lower: Bool = True,
        lowest: Int = 0,
    ):
        self.columns = profile.columns
        self.rows = profile.rows
        self.words = profile.words
        self.difference = self.rows - self.columns
        self.extra = (threshold - abs(self.difference)) // 2
        self.threshold = threshold
        self.stop_column = stop_column
        # One tile's width of horizontal edge, since each tile starts again from `+1` above (see
        # `Sweep.shifted`); the last tile may absorb a sliver of up to `2 * LANES` more columns.
        self.frontier = Frontier(BAND_COLUMNS + 2 * LANES, self.words)
        profile.build_planes()
        self.sweep = self.frontier.sweep(profile)
        self.bounds = tile_bounds(stop_column, NARROW_COLUMNS if threshold < NARROW_BAND else BAND_COLUMNS)
        self.top = 0
        self.end_word = 0
        self.anchor = 0
        self.deepest = 0
        self.floor = 0
        self.outcome = Round(-1, -1, threshold, -1)
        self.adapt = adapt and stop_column == self.columns
        self.lower = lower
        self.lowest = lowest
        self.checkpoint = 0

    def tiles(self) -> Int:
        return len(self.bounds) - 1

    def first_column(self, tile: Int) -> Int:
        return self.bounds[tile]

    def end_column(self, tile: Int) -> Int:
        return self.bounds[tile + 1]

    def prepare[record: Bool](mut self, tile: Int, mut trail: Trail, mut heuristic: SeedHeuristic) -> Bool:
        """Sets the tile's words and readies its edge; false, with `outcome` set, when no word is left."""
        var first_column = self.bounds[tile]
        var end_column = self.bounds[tile + 1]
        var width = end_column - first_column
        # Ukkonen's band bounds the bottom, and so does the deepest kept row: going `k` rows past it
        # costs at least `k` over the score there, at least `floor`, and below the diagonal that leads
        # to the end each such row also adds one to the gap still to close.
        var band_row = end_column + max(0, self.difference) + self.extra
        var slack = self.threshold - self.floor
        var end_diagonal = self.rows - (self.columns - end_column)
        var reach_row = min(
            band_row,
            self.deepest + width + slack,
            max(end_diagonal, (slack + self.deepest + width + end_diagonal) // 2),
            self.rows,
        )
        if heuristic.seeds > 0:
            # The seeds bound the bottom tighter, as A*PA2 bounds a block's: a row `k` past the
            # diagonal from the deepest kept cell costs at least `floor + k` to reach, so it is in
            # reach only while that plus the heuristic there fits the bound. Going up, the first drops
            # by one a row and the heuristic by at most its `climb`, so a row `x` over rules out the
            # next `ceil(x / (1 + climb))` above it unread.
            var diagonal_row = self.deepest + width
            var step = 1 + heuristic.climb()
            while reach_row > diagonal_row:
                var over = self.floor + (reach_row - diagonal_row) + heuristic.h(end_column, reach_row) - self.threshold
                if over <= 0:
                    break
                reach_row = max(reach_row - (over + step - 1) // step, diagonal_row)
        var reach = (reach_row - 1) // WORD_BITS + 1
        # Exactly the words the band reaches: `words` has a kernel for every count. The bottom never
        # rises, which would leave words behind whose differences the next tile still reads.
        self.end_word = max(self.end_word, min(reach, self.words))
        if self.end_word <= self.top:
            self.outcome = Round(-1, first_column, self.threshold, -1)
            return False
        comptime if record:
            trail.record(first_column, end_column, self.top, self.end_word, self.anchor, self.frontier)
        self.frontier.restart_horizontal(width)
        return True

    def tile_sweep(self, tile: Int) -> Sweep:
        """The sweep for one tile, its horizontal edge starting at the tile's first column."""
        return self.sweep.shifted(self.bounds[tile])

    def finish(mut self, tile: Int, mut edge: Edge, mut heuristic: SeedHeuristic) -> Bool:
        """Prunes after a swept tile; false, with `outcome` set, when every row went."""
        var end_column = self.bounds[tile + 1]
        self.anchor += end_column - self.bounds[tile]
        # Read the right edge back as scores and keep, as A*PA2 does, the rows from the first to the
        # last whose score plus heuristic fits the bound. Scores change by at most one per row and the
        # heuristic by at most one going down, so their sum by at most two, and a row `x` over the
        # bound rules out the next `ceil(x / 2)` rows below without reading them; going up, the
        # heuristic may drop by its `climb`, and a row rules out `ceil(x / (1 + climb))` above.
        edge.capture(self.top, self.end_word, self.anchor, self.frontier, self.rows)
        var first_kept = edge.low_row
        var last_kept = edge.high_row
        while first_kept <= last_kept:
            var over = edge.score(first_kept) + heuristic.h(end_column, first_kept) - self.threshold
            if over <= 0:
                break
            first_kept += (over + 1) // 2
        var step = 1 + heuristic.climb()
        while last_kept >= first_kept:
            var over = edge.score(last_kept) + heuristic.h(end_column, last_kept) - self.threshold
            if over <= 0:
                break
            last_kept -= (over + step - 1) // step
        if first_kept > last_kept:
            self.outcome = Round(-1, end_column, self.threshold, -1)
            return False
        if self.adapt and not self.check(end_column, first_kept, last_kept, edge, heuristic):
            return False

        # Move the top down to the word whose top row is at or above the first kept row, carrying the
        # anchor past the words it leaves.
        var new_top = first_kept // WORD_BITS
        for word in range(self.top, new_top):
            self.anchor += word_value(self.frontier.vertical_plus[word], self.frontier.vertical_minus[word])
        self.top = new_top
        # The bottom-most kept cell bounds every path below it: from there each row past the diagonal
        # costs one more, so its exact score is the floor the next tile's reach is measured from.
        self.deepest = last_kept
        self.floor = edge.score(last_kept)
        return True

    def check(
        mut self, end_column: Int, first_kept: Int, last_kept: Int, edge: Edge, mut heuristic: SeedHeuristic
    ) -> Bool:
        """Re-aims the bound from the band's own climb at a checkpoint; false, with `outcome` set, to
        give the round up.

        At an eighth, a quarter and half of the columns, the least score plus heuristic down the kept
        rows, the gap to the end or the seeds still ahead, has climbed from the heuristic at the origin
        about in proportion to the columns crossed,
        so scaled to the whole width it projects the distance from hundreds of edits rather than the
        diagonal transition's handful: within about a tenth at an eighth, closer further on. The
        bound drops to the projection plus `CHECK_MARGIN` tenths for each checkpoint still ahead,
        and a round whose projection less half as much passes its bound gives up, its rest likely
        wasted. Lowering
        the bound keeps the round exact: on an optimal path a cell's score plus its gap to the end is at
        most the distance, so a distance within the final bound passed every column unpruned.
        """
        var passed = self.checkpoint
        while passed < CHECKPOINTS and end_column >= self.columns >> (CHECKPOINTS - passed):
            passed += 1
        if passed == self.checkpoint or end_column >= self.columns:
            return True
        self.checkpoint = passed
        # Sampled every `CHECK_ROWS` rows: the sum changes by at most two a row, a few edits at most.
        var least = Int.MAX
        var row = first_kept
        while True:
            least = min(least, edge.score(row) + heuristic.h(end_column, row))
            if row == last_kept:
                break
            row = min(row + CHECK_ROWS, last_kept)
        var gap = abs(self.difference)
        var origin = heuristic.h(0, 0)
        var estimate = origin + max(least - origin, 0) * self.columns // end_column
        var margin = estimate * (CHECKPOINTS + 1 - passed) * CHECK_MARGIN // 10
        if self.lower and estimate - margin // 2 > self.threshold:
            self.outcome = Round(-1, end_column, self.threshold, estimate)
            return False
        var aim = max(estimate + margin + PROBE_MARGIN, self.lowest)
        if aim < self.threshold:
            self.threshold = aim
            self.extra = (aim - gap) // 2
        return True

    def result(mut self, mut edge: Edge) -> Round:
        """What the round found once every tile is swept, with the last edge captured."""
        edge.capture(self.top, self.end_word, self.anchor, self.frontier, self.rows)
        if self.stop_column < self.columns:
            return Round(-1, self.stop_column, self.threshold, -1)
        if self.end_word < self.words:
            return Round(-1, self.columns, self.threshold, -1)
        return Round(edge.score(self.rows), self.columns, self.threshold, -1)


def pruned_distance[
    record: Bool
](
    mut profile: Profile,
    threshold: Int,
    stop_column: Int,
    mut trail: Trail,
    mut edge: Edge,
    mut heuristic: SeedHeuristic,
    adapt: Bool = False,
    lower: Bool = True,
    lowest: Int = 0,
) -> Round:
    """One round of band doubling with A*PA2-simple's pruning, on one thread.

    Only cells a path of cost at most `threshold` could cross are computed. Ukkonen's band bounds
    them, gap from the start plus gap to the end within the bound, and pruning narrows it further:
    after each tile the absolute scores down its right edge are read back, and a word is dropped for
    good once even its lowest score plus its gap to the end exceeds `threshold`, since no path within
    the bound can pass through it, or below it from above. The bottom grows only to the deepest row
    the lowest kept score could still reach within the bound.

    Every cell outside what is computed reads as the cost of a real path (see `Frontier`), so no
    computed score is ever below the true one, and the scores along an optimal path within the
    bound are exact by induction: its predecessors are exact and kept, so its own rows are never the
    ones dropped. Hence a distance at most `threshold` is the distance.

    The band's top and bottom only move down, so no word reads differences left from an earlier
    column, and each word's differences end as the right edge of the last tile that swept it. With
    the top word reading `+1` from above in every column, the score at the corner reads back as for
    the whole matrix: along the top row, then down the right edge. Tiles cover whole groups of
    `LANES` words from the band's top, growing the bottom to fill the last group, so every word
    takes the vector path.

    With `record`, every tile's left edge goes into `trail` before the tile is swept, for the
    traceback; without, the trail is untouched.

    The round may stop short at `stop_column`, still pruning against the whole problem's end, and
    leaves the scores down that column in `edge`; meeting in the middle joins two such half rounds.
    The distance is reported only when the round runs to the last column.

    With `adapt`, a round over the whole width without seeds re-aims its bound as it goes (see
    `HalfBand.check`), and a distance is exact only within the bound it ends on.
    """
    comptime if record:
        trail.clear()
    var band = HalfBand(profile, threshold, stop_column, adapt, lower, lowest)
    for tile in range(band.tiles()):
        if not band.prepare[record](tile, trail, heuristic):
            return band.outcome
        band.tile_sweep(tile).words(band.top, band.end_word, band.first_column(tile), band.end_column(tile))
        if not band.finish(tile, edge, heuristic):
            return band.outcome
    return band.result(edge)


# region Diagonal transition

comptime UNREACHED_OFFSET = Int32(-(1 << 28))
"""A diagonal no path of the score reaches: far enough below zero that one more column stays negative."""

comptime FRONT_PADDING = 2
"""Unreached diagonals stored either side of a front, so the next front reads its neighbours unchecked."""

comptime FRONT_LANES = 8
"""Diagonals a front step computes at once."""

comptime PROBE_STRIDE = 4
"""Scores between the diagonal transition's checks of its projection."""

comptime PROBE_START = 8
"""The score from which the diagonal transition judges whether to go on, so the projection has a few
edits to go on."""

comptime STEP_TENTHS_DISTANCE = 60
"""Diagonal steps, in tenths, a distance's diagonal transition may still take per column before a band
would be cheaper, at no distance: a step costs about 1.3 ns, a band column about 6 ns."""

comptime STEP_TENTHS_ALIGNMENT = 50
"""`STEP_TENTHS_DISTANCE` when an alignment is wanted: the diagonal transition then keeps every front
and the band records and retraces, which roughly cancel, so the budget sits a little under the
distance's."""

comptime TWO_ENDED_PERCENT = 57
"""The two-ended search's cost per square edit of distance, in percent of one front's step."""

comptime LOCKSTEP_SCORE = 128
"""The score each direction of the two-ended search reaches before a second thread may take one of
them: a front of this many diagonals takes a few times a spin barrier's round trip."""

comptime LOCKSTEP_SETUP = 30_000
"""What starting the second thread costs, in one front's steps of about a nanosecond: a fork and join
of two tasks, with room for one waking from sleep and for a projection that ran high, as one from a
couple of hundred edits can by a fifth."""

comptime CHECK_IN_NANOSECONDS = 100_000
"""How long the first of the lockstep's two tasks waits for the second before going on alone."""

comptime ALONE = Int64(1 << 20)
"""Marks the lockstep claimed by one task alone; past any count of arrivals."""

comptime BARRIER_STEPS = 150
"""What the lockstep's spin barrier costs a score, in one front's steps."""

comptime TWO_ENDED_SETUP = 600
"""What setting up the two-ended search's second front and histories costs, in one front's steps: an
alignment's search runs one front while what it has left costs less than starting from both ends."""

comptime NO_DIAGONAL = Int.MIN
"""No diagonal: the fronts have not met."""

comptime GIVE_UP_SHARE = 24
"""A diagonal transition gives up on its projection only once it has spent at least one part in this
many of what the cheapest band could cost, the band whose bound is the score already searched, since
the distance is at least that.

The projection from a front's first edits assumes them spread along the pair. Real reads and genomes
gather theirs at the ends, a primer or a poly-A tail apart, and a projection from those alone can run
hundreds of times over a distance of a few dozen; giving up on it hands a pair the search would finish
in microseconds to a band over the whole matrix. Spending this share first bounds what a search that
was going to give up anyway wastes to it, of a band's cost."""

comptime PROBE_SPREAD = 1.5
"""Spreads of a projection's noise the search allows for before it gives up (see `noisy_budget`)."""

comptime PROJECTION_ONLY = -(1 << 40)
"""A budget no search fits: the diagonal transition stops at its first check, with its projection."""

comptime EDITS_PER_STEP = 83
"""A band's column grows about 1.3 ns, a step, per this many edits of projected distance, its band
growing taller with the distance."""

comptime UNSEEDED_EDITS_PER_STEP = 83
"""`EDITS_PER_STEP` for a distance whose band runs without seeds over more than `SHORT_COLUMNS`
columns. Measured since the band re-aims its bound at checkpoints, it costs what a seeded one does
per edit, about a step a column per 80 to 100 edits."""

comptime PROBE_MARGIN = 16
"""How far past the projected distance the band's first bound reaches, on top of an eighth of it."""

comptime SHORT_COLUMNS = 4096
"""Pairs up to this many columns aim the band's first bound further past the projection.

On a short pair the projection from a handful of edits ran as much as 1.7 times under the distance,
and a first round that fails there costs nearly a whole round more, where aiming wide costs a band
some rows taller; a long pair's projection landed within about a sixth either way."""

comptime SHORT_AIM = 17
"""A short pair's first bound, in tenths of the projection."""

comptime SHORT_REACH = 128
"""The most a short pair's first bound reaches past the projection, so a projection already too
high, as on very divergent pairs, does not widen the band by most of the matrix."""

comptime CHECKPOINTS = 3
"""A first round without seeds re-aims its bound at `columns >> k` for `k = CHECKPOINTS ..= 1`: an
eighth, a quarter and half of the way across (see `HalfBand.check`)."""

comptime CHECK_MARGIN = 1
"""Tenths of its projection a checkpoint's bound allows, for each checkpoint still to come: its
projection strayed by up to about an eighth at the first, a twentieth at the last."""

comptime CHECK_ROWS = 8
"""Rows between the scores a checkpoint samples down the band."""

comptime PROBE_CEILING = 2048
"""The highest score the diagonal transition reaches, which bounds its memory to a few megabytes."""


struct DiagonalFronts(Movable):
    """Every score's wavefront: for score `s`, the furthest column each diagonal reaches.

    Diagonal `k` holds the cells whose column minus row is `k`. Score `s` keeps diagonals
    `lows[s] ..= highs[s]`, with `FRONT_PADDING` unreached diagonals either side, from position
    `starts[s]` of `offsets`; all of them stay, since the traceback reads them back.
    """

    var offsets: List[Int32]
    var starts: List[Int]
    var lows: List[Int]
    var highs: List[Int]

    def __init__(out self, reserve: Bool = True):
        """Room for a typical search up front, unless `reserve` is false and none is kept."""
        var offsets = 4096 if reserve else 0
        var scores = 64 if reserve else 0
        self.offsets = List[Int32](capacity=offsets)
        self.starts = List[Int](capacity=scores)
        self.lows = List[Int](capacity=scores)
        self.highs = List[Int](capacity=scores)

    @always_inline
    def at(self, score: Int, diagonal: Int) -> Int:
        """The furthest column of `diagonal` at `score`, or -1 when no path of that score reaches it."""
        if diagonal < self.lows[score] or diagonal > self.highs[score]:
            return -1
        var column = Int(self.offsets[self.starts[score] + FRONT_PADDING + diagonal - self.lows[score]])
        return column if column >= 0 else -1


@fieldwise_init
struct Probe(ImplicitlyCopyable, TrivialRegisterPassable):
    """What the diagonal transition learned: the distance, or a projection of it and a floor under it."""

    var distance: Int
    """The edit distance, or -1 when the search stopped first."""
    var estimate: Int
    """Where the distance was heading when the search stopped: its score scaled by the progress made."""
    var floor: Int
    """Every score up to this was searched, so the distance exceeds it."""


@always_inline
def slide_forward(first: ImmPointer[UInt8, _], second: ImmPointer[UInt8, _], column: Int, row: Int) -> Int:
    """How far matches carry `(column, row)` along its diagonal, sixteen bases at a time.

    The sequences' sentinels differ from every base and from each other, so the run stops at the
    matrix's edge by itself (see `CODE_PADDING`).
    """
    var at = column
    var down = row
    while True:
        var low = (
            first.unsafe_offset(at).unsafe_bitcast[UInt64]().unsafe_load()
            ^ second.unsafe_offset(down).unsafe_bitcast[UInt64]().unsafe_load()
        )
        if low != 0:
            return at + Int(count_trailing_zeros(low)) // 8
        var high = (
            first.unsafe_offset(at + 8).unsafe_bitcast[UInt64]().unsafe_load()
            ^ second.unsafe_offset(down + 8).unsafe_bitcast[UInt64]().unsafe_load()
        )
        if high != 0:
            return at + 8 + Int(count_trailing_zeros(high)) // 8
        at += 16
        down += 16


@always_inline
def extend(
    first: ImmPointer[UInt8, _], second: ImmPointer[UInt8, _], column: Int, row: Int, columns: Int, rows: Int
) -> Int:
    """How far matches carry `(column, row)` along its diagonal, eight bases at a time while both have eight."""
    var at = column
    var down = row
    while at + 8 <= columns and down + 8 <= rows:
        var mismatches = (
            first.unsafe_offset(at).unsafe_bitcast[UInt64]().unsafe_load()
            ^ second.unsafe_offset(down).unsafe_bitcast[UInt64]().unsafe_load()
        )
        if mismatches != 0:
            return at + Int(count_trailing_zeros(mismatches)) // 8
        at += 8
        down += 8
    while at < columns and down < rows and first[unsafe_offset=at] == second[unsafe_offset=down]:
        at += 1
        down += 1
    return at


@always_inline
def best_source(fronts: DiagonalFronts, score: Int, diagonal: Int, columns: Int, rows: Int) -> Tuple[Int, UInt8]:
    """The furthest column a path of `score` reaches on `diagonal` before sliding over matches, and its last move.

    One more edit after the furthest point of score `score - 1` on this diagonal or a neighbour,
    each only where that edit stays inside the matrix; -1 when none does.
    """
    var best = -1
    var move = DIAGONAL
    var same = fronts.at(score - 1, diagonal)
    if same >= 0 and same < columns and same - diagonal < rows:
        best = same + 1
    # A base of the first sequence against a gap, from the diagonal below.
    var left = fronts.at(score - 1, diagonal - 1)
    if left >= 0 and left < columns and left + 1 > best:
        best = left + 1
        move = LEFT
    # A base of the second sequence against a gap, from the diagonal above.
    var up = fronts.at(score - 1, diagonal + 1)
    if up >= 0 and up - diagonal - 1 < rows and up > best:
        best = up
        move = UP
    return (best, move)


@always_inline
def noisy_budget(budget: Int, edits: Int) -> Int:
    """`budget` widened for a projection from `edits` edits, about one in `sqrt(edits)` off.

    A projection from a few edits can land far above the distance, and giving up on it then hands a
    pair the band would align more slowly; so an alignment's search goes on until the projected work
    passes the budget by more than the projection's noise, `PROBE_SPREAD` of its spreads. A distance
    has a cheaper band to fall back on, and gives up on the plain budget.
    """
    if budget < 0:
        return budget
    return Int(Float64(budget) * (1.0 + PROBE_SPREAD / sqrt(Float64(max(edits, 1)))))


@always_inline
def step_budget(columns: Int, step_tenths: Int, estimate: Int, edits_per_step: Int = EDITS_PER_STEP) -> Int:
    """How many diagonal steps cost what a band over `columns` columns would, at a projected distance.

    A band's column costs a fixed part plus a part growing with the distance, its band taller; in
    steps that is `step_tenths / 10 + estimate / edits_per_step` a column.
    """
    return columns * step_tenths // 10 + columns * estimate // edits_per_step


@always_inline
def slide_from(first: ImmPointer[UInt8, _], second: ImmPointer[UInt8, _], start: Int, diagonal: Int) -> Int:
    """How far matches carry column `start` of `diagonal`, eight bases at a time, to the sentinels."""
    var column = start
    var lag = second.unsafe_offset(-diagonal)
    var mismatches = (
        first.unsafe_offset(column).unsafe_bitcast[UInt64]().unsafe_load()
        ^ lag.unsafe_offset(column).unsafe_bitcast[UInt64]().unsafe_load()
    )
    while mismatches == 0:
        column += 8
        mismatches = (
            first.unsafe_offset(column).unsafe_bitcast[UInt64]().unsafe_load()
            ^ lag.unsafe_offset(column).unsafe_bitcast[UInt64]().unsafe_load()
        )
    return column + (Int(count_trailing_zeros(mismatches)) >> 3)


@no_inline
def step_front[
    measure: Bool
](
    previous: MutPointer[Int32, MutUntrackedOrigin],
    current: MutPointer[Int32, MutUntrackedOrigin],
    low: Int,
    high: Int,
    columns: Int,
    rows: Int,
    first: ImmPointer[UInt8, _],
    second: ImmPointer[UInt8, _],
) -> Int:
    """One score's front from the last: every diagonal `low ..= high` of `current`, both indexed by diagonal.

    A function of its own, so its handful of values stay in registers. With `measure`, returns
    the furthest anti-diagonal, column plus row, any diagonal reached; otherwise zero.
    """
    comptime Lanes = SIMD[DType.int32, FRONT_LANES]
    # Every diagonal's start, eight at a time: one more edit after the previous front on the same
    # diagonal or either neighbour, each only where that edit stays inside the matrix. The previous
    # front's padding reads as unreached, far below any column.
    var lane_diagonals = Lanes()
    comptime for lane in range(FRONT_LANES):
        lane_diagonals[lane] = Int32(lane)
    var column_limit = Lanes(Int32(columns))
    var row_limit = Lanes(Int32(rows))
    var unreached = Lanes(UNREACHED_OFFSET)
    var diagonal = low
    while diagonal <= high:
        var diagonals = lane_diagonals + Int32(diagonal)
        var same = previous.unsafe_offset(diagonal).unsafe_load[width=FRONT_LANES]()
        var below = previous.unsafe_offset(diagonal - 1).unsafe_load[width=FRONT_LANES]()
        var above = previous.unsafe_offset(diagonal + 1).unsafe_load[width=FRONT_LANES]()
        var substituted = (same.lt(column_limit) & (same - diagonals).lt(row_limit)).select(same + 1, unreached)
        # A base of the first sequence against a gap, from the diagonal below.
        var deleted = below.lt(column_limit).select(below + 1, unreached)
        # A base of the second sequence against a gap, from the diagonal above.
        var inserted = (above - diagonals).le(row_limit).select(above, unreached)
        current.unsafe_offset(diagonal).unsafe_store(max(substituted, max(deleted, inserted)))
        diagonal += FRONT_LANES
    for index in range(FRONT_PADDING):
        current[unsafe_offset=high + 1 + index] = UNREACHED_OFFSET

    # Then each slides over its matches, eight bases at a time; the sentinels stop it at the edge.
    # Two diagonals a turn, so their loads overlap and they share the loop's bookkeeping.
    var furthest = 0
    diagonal = low
    while diagonal + 1 <= high:
        var left = Int(current[unsafe_offset=diagonal])
        var right = Int(current[unsafe_offset=diagonal + 1])
        if left >= 0:
            left = slide_from(first, second, left, diagonal)
            current[unsafe_offset=diagonal] = Int32(left)
            comptime if measure:
                furthest = max(furthest, 2 * left - diagonal)
        if right >= 0:
            right = slide_from(first, second, right, diagonal + 1)
            current[unsafe_offset=diagonal + 1] = Int32(right)
            comptime if measure:
                furthest = max(furthest, 2 * right - diagonal - 1)
        diagonal += 2
    if diagonal <= high:
        var column = Int(current[unsafe_offset=diagonal])
        if column >= 0:
            column = slide_from(first, second, column, diagonal)
            current[unsafe_offset=diagonal] = Int32(column)
            comptime if measure:
                furthest = max(furthest, 2 * column - diagonal)
    return furthest


def diagonal_transition(
    profile: Profile,
    step_tenths: Int,
    mut fronts: DiagonalFronts,
    ceiling: Int = PROBE_CEILING,
    switch_setup: Int = -1,
) -> Probe:
    """The edit distance by diagonal transition, as WFA computes it, while it stays cheaper than a band.

    Score `s` reaches, on every diagonal, the furthest cell some path of `s` edits reaches; matches
    are free, so each front slides as far as they carry it. The distance is the first score whose
    front reaches the corner, after about `d²` diagonals and the matches along the way, where the
    band's sweep pays for every column whatever the distance. So near-identical pairs, the case the
    post on A*PA2 names as its weak spot, finish here.

    From `PROBE_START` on, the furthest anti-diagonal reached projects the distance, and once the
    search still to do, the projection squared less the score squared, passes `budget`, it stops
    and hands the projection on as the band's first bound.
    """
    var columns = profile.columns
    var rows = profile.rows
    var first = profile.column_codes.unsafe_ptr()
    var second = profile.row_codes.unsafe_ptr()
    var target = columns - rows
    fronts.offsets.clear()
    fronts.starts.clear()
    fronts.lows.clear()
    fronts.highs.clear()

    var start = slide_forward(first, second, 0, 0)
    fronts.starts.append(0)
    fronts.lows.append(0)
    fronts.highs.append(0)
    for _ in range(FRONT_PADDING):
        fronts.offsets.append(UNREACHED_OFFSET)
    fronts.offsets.append(Int32(start))
    for _ in range(FRONT_PADDING):
        fronts.offsets.append(UNREACHED_OFFSET)
    if target == 0 and start == columns:
        return Probe(0, 0, 0)
    var score = 0
    while score < ceiling:
        score += 1
        var previous_low = fronts.lows[score - 1]
        var previous_start = fronts.starts[score - 1]
        var low = max(-score, -rows)
        var high = min(score, columns)
        var count = high - low + 1
        var row_start = len(fronts.offsets)
        fronts.starts.append(row_start)
        fronts.lows.append(low)
        fronts.highs.append(high)
        # Room for the front, its padding, and a last vector's spill past the end.
        fronts.offsets.resize(unsafe_uninit_length=row_start + count + 2 * FRONT_PADDING + FRONT_LANES)
        var offsets = fronts.offsets.unsafe_ptr()
        for index in range(FRONT_PADDING):
            offsets[unsafe_offset=row_start + index] = UNREACHED_OFFSET
        # The previous front and the new one, both indexed by diagonal.
        var previous = offsets.unsafe_origin_cast[MutUntrackedOrigin]().unsafe_offset(
            previous_start + FRONT_PADDING - previous_low
        )
        var current = offsets.unsafe_origin_cast[MutUntrackedOrigin]().unsafe_offset(row_start + FRONT_PADDING - low)

        # The projection is checked every `PROBE_STRIDE` scores, so most fronts skip measuring it.
        var checking = score >= PROBE_START and score % PROBE_STRIDE == 0
        var furthest: Int
        if checking:
            furthest = step_front[True](previous, current, low, high, columns, rows, first, second)
        else:
            furthest = step_front[False](previous, current, low, high, columns, rows, first, second)
        if target >= low and target <= high and Int(current[unsafe_offset=target]) == columns:
            fronts.offsets.resize(unsafe_uninit_length=row_start + count + 2 * FRONT_PADDING)
            return Probe(score, score, score - 1)
        fronts.offsets.resize(unsafe_uninit_length=row_start + count + 2 * FRONT_PADDING)
        if checking:
            var estimate = score * (columns + rows) // max(furthest, 1)
            # With `switch_setup`, give way to the two-ended search once what this front has left costs
            # more than that search would from scratch, setup included; the work done is spent.
            if switch_setup >= 0 and estimate * estimate - score * score > (
                estimate * estimate * TWO_ENDED_PERCENT // 100 + switch_setup
            ):
                return Probe(-1, max(estimate, score + 1), score)
            # What is left to search, about `estimate² - score²` diagonals, against what a band
            # would cost; the work already done is spent either way. Only once the search has spent
            # its share of the cheapest band there could be (see `GIVE_UP_SHARE`).
            if estimate * estimate - score * score > noisy_budget(
                step_budget(columns, step_tenths, estimate), score
            ) and score * score * GIVE_UP_SHARE >= step_budget(columns, step_tenths, score):
                return Probe(-1, max(estimate, score + 1), score)
    return Probe(-1, ceiling + 1, ceiling)


def reversed_codes(codes: List[UInt8], count: Int, sentinel: UInt8) -> List[UInt8]:
    """The first `count` codes back to front, with `CODE_PADDING` sentinels after them."""
    var flipped = List[UInt8](capacity=count + CODE_PADDING)
    flipped.resize(unsafe_uninit_length=count)
    var source = codes.unsafe_ptr()
    var target = flipped.unsafe_ptr()
    comptime CHUNK = 16
    var index = 0
    while index + CHUNK <= count:
        target.unsafe_offset(index).unsafe_store(
            source.unsafe_offset(count - CHUNK - index).unsafe_load[width=CHUNK]().reversed()
        )
        index += CHUNK
    while index < count:
        target[unsafe_offset=index] = source[unsafe_offset=count - 1 - index]
        index += 1
    for _ in range(CODE_PADDING):
        flipped.append(sentinel)
    return flipped^


comptime FRONT_RING = 3
"""Fronts a direction keeps in turn: the latest, the one before, which an overlap check reads, and the
next, which a second thread may already be growing while the first checks the other two."""


struct FrontPair(Movable):
    """One direction's latest fronts, in a ring of `FRONT_RING`, each indexed by diagonal."""

    var buffers: List[Int32]
    var width: Int
    var low: Int
    var high: Int
    var previous_low: Int
    """The previous front's diagonals, which the ring still holds."""
    var previous_high: Int
    var score: Int
    var slot: Int
    """Which of the ring holds the latest front."""
    var furthest: Int
    var history: DiagonalFronts
    """Every front so far, when `record`, laid out as `diagonal_transition` keeps them, for a traceback."""
    var record: Bool

    def __init__(out self, limit: Int, record: Bool = False):
        self.width = 2 * (limit + FRONT_PADDING) + FRONT_LANES + 1
        # Every front writes its own diagonals and the padding either side before the next reads it,
        # so only the first front's surroundings need setting.
        self.buffers = List[Int32](capacity=FRONT_RING * self.width)
        self.buffers.resize(unsafe_uninit_length=FRONT_RING * self.width)
        self.low = 0
        self.high = 0
        self.previous_low = 0
        self.previous_high = -1
        self.score = 0
        self.slot = 0
        self.furthest = 0
        self.history = DiagonalFronts(reserve=record)
        self.record = record
        var first = self.front(0)
        for diagonal in range(-FRONT_PADDING, FRONT_PADDING + 1):
            first[unsafe_offset=diagonal] = UNREACHED_OFFSET

    def keep(mut self):
        """Copies the latest front, with its padding, onto the history, when recording."""
        if not self.record:
            return
        var source = self.front(self.slot)
        var start = len(self.history.offsets)
        self.history.starts.append(start)
        self.history.lows.append(self.low)
        self.history.highs.append(self.high)
        var count = self.high - self.low + 1 + 2 * FRONT_PADDING
        self.history.offsets.resize(unsafe_uninit_length=start + count)
        var target = self.history.offsets.unsafe_ptr().unsafe_offset(start)
        for index in range(count):
            target[unsafe_offset=index] = source[unsafe_offset=self.low - FRONT_PADDING + index]

    @always_inline
    def previous(self) -> Int:
        """The ring slot of the front before the latest."""
        return self.slot - 1 if self.slot > 0 else FRONT_RING - 1

    @always_inline
    def state(self) -> FrontState:
        """What another thread needs to read this direction's two latest fronts."""
        return FrontState(
            self.slot,
            self.low,
            self.high,
            self.previous(),
            self.previous_low,
            self.previous_high,
            self.score,
            self.furthest,
        )

    @always_inline
    def advance(
        mut self, codes: ImmPointer[UInt8, _], others: ImmPointer[UInt8, _], measure: Bool, columns: Int, rows: Int
    ):
        """One more score: the next front in the ring, grown from the latest."""
        var previous = self.front(self.slot)
        var next = self.slot + 1 if self.slot + 1 < FRONT_RING else 0
        var current = self.front(next)
        self.score += 1
        var low = max(-self.score, -rows)
        var high = min(self.score, columns)
        for index in range(1, FRONT_PADDING + 1):
            current[unsafe_offset=low - index] = UNREACHED_OFFSET
        if measure:
            self.furthest = step_front[True](previous, current, low, high, columns, rows, codes, others)
        else:
            _ = step_front[False](previous, current, low, high, columns, rows, codes, others)
        self.previous_low = self.low
        self.previous_high = self.high
        self.low = low
        self.high = high
        self.slot = next
        self.keep()

    @always_inline
    def front(self, which: Int) -> MutPointer[Int32, MutUntrackedOrigin]:
        """Front `which` of the ring, as a pointer indexed by diagonal.

        Writable from a shared reference: in the lockstep each thread writes only its own fronts, and
        reads the other's only where their owner no longer writes (see `lockstep`).
        """
        return (
            self.buffers.unsafe_ptr()
            .unsafe_mut_cast[True]()
            .unsafe_origin_cast[MutUntrackedOrigin]()
            .unsafe_offset(which * self.width + self.width // 2)
        )


@fieldwise_init
struct FrontState(ImplicitlyCopyable, TrivialRegisterPassable):
    """A direction's two latest fronts as one moment saw them: their ring slots and diagonals, its
    score and its progress, so a thread can check them while their owner grows the next."""

    var slot: Int
    var low: Int
    var high: Int
    var previous_slot: Int
    var previous_low: Int
    var previous_high: Int
    var score: Int
    var furthest: Int


@always_inline
def overlap(
    forward: MutPointer[Int32, MutUntrackedOrigin],
    forward_low: Int,
    forward_high: Int,
    backward: MutPointer[Int32, MutUntrackedOrigin],
    back_low: Int,
    back_high: Int,
    target: Int,
    columns: Int,
) -> Int:
    """A diagonal whose forward front reaches past where the backward front comes back to, or
    `NO_DIAGONAL`; the backward front's diagonal `target - k` mirrors the forward's `k`."""
    var low = max(forward_low, target - back_high)
    var high = min(forward_high, target - back_low)
    # An unreached diagonal holds a value far below zero, so the sum alone rules it out.
    var needed = SIMD[DType.int32, FRONT_LANES](Int32(columns))
    var diagonal = low
    while diagonal + FRONT_LANES - 1 <= high:
        var reached = forward.unsafe_offset(diagonal).unsafe_load[width=FRONT_LANES]()
        var back = (
            backward.unsafe_offset(target - diagonal - FRONT_LANES + 1).unsafe_load[width=FRONT_LANES]().reversed()
        )
        if (reached + back).ge(needed).reduce_or():
            break
        diagonal += FRONT_LANES
    while diagonal <= high:
        if Int(forward[unsafe_offset=diagonal]) + Int(backward[unsafe_offset=target - diagonal]) >= columns:
            return diagonal
        diagonal += 1
    return NO_DIAGONAL


@fieldwise_init
struct Meeting(ImplicitlyCopyable, TrivialRegisterPassable):
    """Where the two-ended search's fronts met: the forward front's furthest cell on a diagonal both
    reached, and each side's score there; or, when it gave up, `probe` alone."""

    var probe: Probe
    var diagonal: Int
    var column: Int
    var forward_score: Int
    var backward_score: Int


def two_ended_distance(profile: Profile, step_tenths: Int, workers: Int = 1) -> Probe:
    """The edit distance by `two_ended`, keeping no history."""
    var first_back = reversed_codes(profile.column_codes, profile.columns, FIRST_SENTINEL)
    var second_back = reversed_codes(profile.row_codes, profile.rows, SECOND_SENTINEL)
    var ahead = FrontPair(PROBE_CEILING // 2 + 1)
    var behind = FrontPair(PROBE_CEILING // 2 + 1)
    return two_ended(
        profile, first_back, second_back, step_tenths, ahead, behind, UNSEEDED_EDITS_PER_STEP, workers
    ).probe


def two_ended(
    profile: Profile,
    first_back: List[UInt8],
    second_back: List[UInt8],
    step_tenths: Int,
    mut ahead: FrontPair,
    mut behind: FrontPair,
    unseeded_edits_per_step: Int = EDITS_PER_STEP,
    workers: Int = 1,
) -> Meeting:
    """The edit distance by diagonal transition from both ends at once, as BiWFA scores, while cheap.

    One front grows from the start and one from the end, over the reversed sequences, a score at a
    time each in turn, and the distance is the first total score at which they overlap: on some
    diagonal the forward front reaches at least as far as the backward one comes back. An overlap
    joins a real path from the start to one to the end, so the total is at least the distance; and
    every cost level an optimal path passes splits it in two whose halves the fronts have reached
    by the time their scores add up to the distance, edit costs never falling along a diagonal. So
    the first overlap is exact, after about half the diagonals of one front grown alone, keeping
    only the last two fronts each way.

    It stops as `diagonal_transition` does, once the search still to do passes the budget (see
    `step_budget`), projecting the distance from both fronts' progress; against a band the seeds
    would not narrow, at `unseeded_edits_per_step`. With more than one worker, once both fronts
    reach `LOCKSTEP_SCORE` the two directions go on in parallel (see `lockstep`); `workers` of zero
    means every thread this process may use, asked for only then.
    """
    var columns = profile.columns
    var rows = profile.rows
    var target = columns - rows
    var first = profile.column_codes.unsafe_ptr()
    var second = profile.row_codes.unsafe_ptr()
    var first_reversed = first_back.unsafe_ptr()
    var second_reversed = second_back.unsafe_ptr()
    var limit = PROBE_CEILING // 2 + 1
    ahead.front(0)[unsafe_offset=0] = Int32(slide_forward(first, second, 0, 0))
    behind.front(0)[unsafe_offset=0] = Int32(slide_forward(first_reversed, second_reversed, 0, 0))
    ahead.keep()
    behind.keep()

    @always_inline
    def overlapping(
        ahead: FrontPair, back: FrontState, forward_earlier: Bool, backward_earlier: Bool, backward_buffers: FrontPair
    ) {imm target, imm columns} -> Int:
        """`overlap` of the forward front, or the one before it, and the backward front `back` saw,
        or the one before it."""
        return overlap(
            ahead.front(ahead.previous() if forward_earlier else ahead.slot),
            ahead.previous_low if forward_earlier else ahead.low,
            ahead.previous_high if forward_earlier else ahead.high,
            backward_buffers.front(back.previous_slot if backward_earlier else back.slot),
            back.previous_low if backward_earlier else back.low,
            back.previous_high if backward_earlier else back.high,
            target,
            columns,
        )

    @always_inline
    def met(
        ahead: FrontPair, back: FrontState, diagonal: Int, forward_earlier: Bool, backward_earlier: Bool
    ) -> Meeting:
        var column = Int(ahead.front(ahead.previous() if forward_earlier else ahead.slot)[unsafe_offset=diagonal])
        var forward_score = ahead.score - 1 if forward_earlier else ahead.score
        var backward_score = back.score - 1 if backward_earlier else back.score
        var total = forward_score + backward_score
        return Meeting(Probe(total, total, total - 1), diagonal, column, forward_score, backward_score)

    var meeting = overlapping(ahead, behind.state(), False, False, behind)
    if meeting != NO_DIAGONAL:
        return met(ahead, behind.state(), meeting, False, False)
    var total = 0
    while ahead.score < limit and behind.score < limit:
        total += 1
        # Each side measures its progress on its step just before a check, every `PROBE_STRIDE`, from
        # when each has taken `PROBE_START` edits, so the projection has as many to go on as one front's.
        var checking = total >= 2 * PROBE_START and total % PROBE_STRIDE == 0
        var measuring = total + 1 >= 2 * PROBE_START and (total + 1) % PROBE_STRIDE == 0
        if total % 2 == 1:
            ahead.advance(first, second, checking or measuring, columns, rows)
        else:
            behind.advance(first_reversed, second_reversed, checking or measuring, columns, rows)
        # Checked after the backward front's steps alone: a first overlap one score sooner is found
        # then, against the backward front's previous score, which the ring still holds.
        if total % 2 == 0:
            var back = behind.state()
            meeting = overlapping(ahead, back, False, False, behind)
            if meeting != NO_DIAGONAL:
                var sooner = overlapping(ahead, back, False, True, behind)
                if sooner != NO_DIAGONAL:
                    return met(ahead, back, sooner, False, True)
                return met(ahead, back, meeting, False, False)
        if checking:
            var estimate = two_ended_gives_up(
                total, ahead.furthest + behind.furthest, columns, rows, step_tenths, unseeded_edits_per_step
            )
            if estimate >= 0:
                return Meeting(Probe(-1, estimate, total), 0, 0, 0, 0)
            # Both fronts stand at the same score after a check; a second thread takes one of them
            # once half the steps projected to be left outweigh its start and a barrier a score.
            if workers != 1 and total >= 2 * LOCKSTEP_SCORE:
                var projected = total * (columns + rows) // max(ahead.furthest + behind.furthest, 1)
                var saved = (projected * projected - total * total) * TWO_ENDED_PERCENT // 200
                # Asking for the thread count is a system call, so only a search this far along asks.
                var available = workers if workers > 0 else max(hardware_threads(), 1)
                if available > 1 and saved > LOCKSTEP_SETUP + (projected - total) // 2 * BARRIER_STEPS:
                    return lockstep(
                        profile, first_back, second_back, ahead, behind, limit, step_tenths, unseeded_edits_per_step
                    )
    return Meeting(Probe(-1, 2 * limit + 1, total), 0, 0, 0, 0)


comptime DOUBLING_START = 256
"""How far past the distance already ruled out an untrusted band's first bound reaches: A*PA2's own
first step past the heuristic at the origin."""

comptime AGREEMENT = 3
"""Two projections agree when neither exceeds the other by more than one part in this many."""


def trusted_projection(search: Probe, projected: Probe) -> Bool:
    """Whether one front's projection from its first edits agrees with the search's from many more.

    On pairs whose errors spread along them, as mutated sequences' do, the two land within a few
    percent of each other and of the distance, and the band aims straight at them. Real reads gather
    errors at their ends, a primer or a tail, and there the two disagree several times over; neither
    is then worth aiming at.
    """
    if search.floor < 2 * PROBE_START or search.estimate <= search.floor:
        return True
    var first = projected.estimate
    var later = search.estimate
    return first * AGREEMENT <= later * (AGREEMENT + 1) and later * AGREEMENT <= first * (AGREEMENT + 1)


def two_ended_gives_up(
    total: Int, reached: Int, columns: Int, rows: Int, step_tenths: Int, unseeded_edits_per_step: Int
) -> Int:
    """The two-ended search's projected distance once what is left of it passes the budget, or -1."""
    var estimate = total * (columns + rows) // max(reached, 1)
    # What is left, about half of `estimate² - total²` diagonals, against what a band costs.
    # Both fronts' steps, about 0.57 ns per square edit with the overlap check, in one-front steps.
    var seeded = estimate >= SEED_EDITS and estimate * SEED_DIVERGENCE <= columns
    var per_step = EDITS_PER_STEP if seeded or columns <= SHORT_COLUMNS else unseeded_edits_per_step
    if (estimate * estimate - total * total) * TWO_ENDED_PERCENT // 100 > step_budget(
        columns, step_tenths, estimate, per_step
    ) and total * total * TWO_ENDED_PERCENT // 100 * GIVE_UP_SHARE >= step_budget(
        columns, step_tenths, total, per_step
    ):
        return max(estimate, total + 1)
    return -1


def lockstep(
    profile: Profile,
    first_back: List[UInt8],
    second_back: List[UInt8],
    mut ahead: FrontPair,
    mut behind: FrontPair,
    limit: Int,
    step_tenths: Int,
    unseeded_edits_per_step: Int,
) -> Meeting:
    """`two_ended` from equal scores on, with one thread a direction, a score at a time each.

    After both have grown the fronts of score `s`, the first thread checks the totals `2s - 1`, one
    front of each pair against the other's previous, then `2s`, in that order, so the first overlap
    is still the least total; the second thread meanwhile grows its next front into the third slot
    of its ring, having published where its two latest lie (see `FrontState`). One spin barrier a
    score: a few dozen nanoseconds, against the hundreds a front of a few hundred diagonals takes.

    The runtime may run the two tasks one after the other, and a task spinning for a partner that
    starts only once it returns would never return. So each checks in first: the first to arrive
    waits `CHECK_IN_NANOSECONDS` for the other, and if none comes, claims the search and finishes it
    alone, both directions in turn, while the late task finds it claimed and returns.
    """
    var columns = profile.columns
    var rows = profile.rows
    var first = profile.column_codes.unsafe_ptr()
    var second = profile.row_codes.unsafe_ptr()
    var first_reversed = first_back.unsafe_ptr()
    var second_reversed = second_back.unsafe_ptr()
    var target = columns - rows
    var arrived = Atomic[Int64](0)
    var generation = Atomic[Int64](0)
    var done = Atomic[Int64](0)
    var present = Atomic[Int64](0)
    var start = ahead.score
    # What a search that ran out of scores reports; any other end overwrites it.
    var outcome = Meeting(Probe(-1, 2 * limit + 1, 2 * (limit - 1)), 0, 0, 0, 0)
    # The backward direction's latest state, one record per score parity, so publishing the next
    # never overwrites the one being checked.
    var published = List[FrontState](length=2, fill=behind.state())

    @always_inline
    def wait(mut arrived: Atomic[Int64], mut generation: Atomic[Int64]):
        var seen = generation.load()
        if arrived.fetch_add(1) == 1:
            arrived.store(0)
            _ = generation.fetch_add(1)
        else:
            while generation.load() == seen:
                pass

    def side(
        index: Int,
    ) {
        mut ahead,
        mut behind,
        mut arrived,
        mut generation,
        mut done,
        mut present,
        mut outcome,
        mut published,
        imm start,
        imm limit,
        imm columns,
        imm rows,
        imm target,
        imm step_tenths,
        imm unseeded_edits_per_step,
        imm first,
        imm second,
        imm first_reversed,
        imm second_reversed,
    }:
        @always_inline
        def settled(
            score: Int, back: FrontState, checking: Bool
        ) {
            mut outcome,
            imm ahead,
            imm behind,
            imm target,
            imm columns,
            imm rows,
            imm step_tenths,
            imm unseeded_edits_per_step,
        } -> Bool:
            """Checks the totals `2 score - 1` then `2 score`, and the budget; true once the search ends."""
            for earlier in range(3):
                # Either front against the other's previous first, then both latest.
                var forward_earlier = earlier == 1
                var backward_earlier = earlier == 0
                var diagonal = overlap(
                    ahead.front(ahead.previous() if forward_earlier else ahead.slot),
                    ahead.previous_low if forward_earlier else ahead.low,
                    ahead.previous_high if forward_earlier else ahead.high,
                    behind.front(back.previous_slot if backward_earlier else back.slot),
                    back.previous_low if backward_earlier else back.low,
                    back.previous_high if backward_earlier else back.high,
                    target,
                    columns,
                )
                if diagonal != NO_DIAGONAL:
                    var column = Int(
                        ahead.front(ahead.previous() if forward_earlier else ahead.slot)[unsafe_offset=diagonal]
                    )
                    var forward_score = ahead.score - 1 if forward_earlier else ahead.score
                    var backward_score = back.score - 1 if backward_earlier else back.score
                    var total = forward_score + backward_score
                    outcome = Meeting(Probe(total, total, total - 1), diagonal, column, forward_score, backward_score)
                    return True
            if checking:
                var estimate = two_ended_gives_up(
                    2 * score, ahead.furthest + back.furthest, columns, rows, step_tenths, unseeded_edits_per_step
                )
                if estimate >= 0:
                    outcome = Meeting(Probe(-1, estimate, 2 * score), 0, 0, 0, 0)
                    return True
            return False

        var arrival = present.fetch_add(1)
        if arrival >= ALONE:
            return
        if arrival == 0:
            var waited_from = perf_counter_ns()
            while present.load() == 1:
                if Int(perf_counter_ns() - waited_from) > CHECK_IN_NANOSECONDS:
                    var expected = Int64(1)
                    if present.compare_exchange(expected, ALONE):
                        break
            if present.load() >= ALONE:
                # No partner came: both directions on this thread, in turn.
                var score = start
                while score + 1 < limit:
                    score += 1
                    var checking = (2 * score) % PROBE_STRIDE == 0
                    ahead.advance(first, second, checking, columns, rows)
                    behind.advance(first_reversed, second_reversed, checking, columns, rows)
                    if settled(score, behind.state(), checking):
                        return
                return
        var score = start
        while score + 1 < limit:
            score += 1
            var checking = (2 * score) % PROBE_STRIDE == 0
            if index == 1:
                behind.advance(first_reversed, second_reversed, checking, columns, rows)
                published[score % 2] = behind.state()
                wait(arrived, generation)
                # The first thread decides a score's end only after that score's barrier, so an end
                # at this very score may already show; this thread then grows one more front and
                # meets it at the next barrier, where it waits, and an earlier end lets it go.
                var ended = Int(done.load())
                if ended != 0 and ended < score:
                    return
                continue
            ahead.advance(first, second, checking, columns, rows)
            wait(arrived, generation)
            if settled(score, published[score % 2], checking):
                done.store(Int64(score))
                # The other thread is growing its next front, if there is a next, and waits there.
                if score + 1 < limit:
                    wait(arrived, generation)
                return

    parallelize(side, 2, 2)
    return outcome


def trace_from(
    fronts: DiagonalFronts,
    columns: Int,
    rows: Int,
    score: Int,
    start_diagonal: Int,
    start_column: Int,
    mut moves: List[UInt8],
):
    """An optimal path from the origin to any cell a front of `score` reached, as moves right to left.

    Edit costs never fall along a diagonal, so a cell's cost is the least score whose front on its
    diagonal reaches it, and every cell a front of score `s` reaches costs at most `s`. At a cell of
    cost `s`, an edit can enter it from a cell of cost `s - 1` only up to the furthest of the
    previous front's three starts there; past that, the cells of cost `s` are reached by matches
    alone. So each step slides back to `min(cell, furthest start)` and takes an edit whose source
    reaches that far, landing on a cell of cost at most `s - 1`; the cell need not be a front's
    furthest, as where two fronts met.
    """
    var diagonal = start_diagonal
    var column = start_column
    var level = score
    while True:
        while level > 0 and fronts.at(level - 1, diagonal) >= column:
            level -= 1
        if level == 0:
            break
        var same = fronts.at(level - 1, diagonal)
        var below = fronts.at(level - 1, diagonal - 1)
        var above = fronts.at(level - 1, diagonal + 1)
        # One more edit after the previous front, only where it stays inside the matrix.
        var substituted = same + 1 if same >= 0 and same < columns and same - diagonal < rows else -1
        var deleted = below + 1 if below >= 0 and below < columns else -1
        var inserted = above if above >= 0 and above - diagonal - 1 < rows else -1
        var entry = min(column, max(substituted, max(deleted, inserted)))
        for _ in range(column - entry):
            moves.append(DIAGONAL)
        if substituted >= entry and entry >= 1 and entry - diagonal >= 1:
            moves.append(DIAGONAL)
            column = entry - 1
        elif deleted >= entry and entry >= 1:
            moves.append(LEFT)
            column = entry - 1
            diagonal -= 1
        else:
            moves.append(UP)
            column = entry
            diagonal += 1
        level -= 1
    for _ in range(column):
        moves.append(DIAGONAL)


def trace_diagonals(profile: Profile, fronts: DiagonalFronts, distance: Int, mut moves: List[UInt8]):
    """The optimal path `diagonal_transition` found, as moves right to left, like `trace_back`.

    From the corner, each score's front is undone: the matches it slid over, then the edit that
    reached its start from the furthest front of one score less, which `best_source` finds again.
    """
    var columns = profile.columns
    var rows = profile.rows
    var diagonal = columns - rows
    var column = columns
    var score = distance
    while score > 0:
        var source = best_source(fronts, score, diagonal, columns, rows)
        var best = source[0]
        for _ in range(column - best):
            moves.append(DIAGONAL)
        var move = source[1]
        moves.append(move)
        if move == DIAGONAL:
            column = best - 1
        elif move == LEFT:
            column = best - 1
            diagonal -= 1
        else:
            column = best
            diagonal += 1
        score -= 1
    for _ in range(column):
        moves.append(DIAGONAL)


# endregion Diagonal transition


@fieldwise_init
struct Outcome(ImplicitlyCopyable, TrivialRegisterPassable):
    """Where band doubling ended: the distance, and the cell at which the two halves met."""

    var distance: Int
    """The edit distance, or -1 when the band grew to the whole matrix and the caller should sweep it."""
    var middle: Int
    """The column the forward half stopped at; every column when it ran alone."""
    var meeting_row: Int
    """The row of `middle` an optimal path crosses."""


def band_doubling[
    record: Bool
](
    mut forward: Profile,
    mut backward: Profile,
    meet: Bool,
    give_up_wide: Bool,
    probe: Probe,
    mut forward_trail: Trail,
    mut backward_trail: Trail,
    mut forward_edge: Edge,
    mut backward_edge: Edge,
    mut forward_heuristic: SeedHeuristic,
    mut backward_heuristic: SeedHeuristic,
    trusted: Bool = True,
) -> Outcome:
    """Band doubling, from the start alone or from both ends at once to meet in the middle.

    An untrusted projection (see `trusted_projection`) neither aims the first bound nor sets the
    next: the first starts `DOUBLING_START` past the distance already ruled out, by the floor or the
    heuristic at the origin, and each next one at
    most doubles, as A*PA2's band doubling grows, whatever a failed round projected.

    Each round sweeps the band for one bound. Meeting in the middle, one thread sweeps forward to
    the middle column while another sweeps the reversed sequences forward to the same column, which
    is the original's suffix swept from its end; the distance is the least, over the middle's rows,
    of the two scores there. Every computed score is the cost of a real path and the optimal path's
    are exact (see `pruned_distance`), so that least is exact whenever it is within the bound.

    A failed round leaves either a real alignment's cost, which caps the next bound, or the column
    at which a half pruned every row, from which the distance is estimated as the bound scaled to the
    whole width; the next bound aims just past the estimate rather than doubling. Once the band would
    cover the matrix, `give_up_wide` hands back to the caller, else the bound is lifted past any path.

    The first bound aims just past where the diagonal transition's `probe` projected the distance,
    and above every score it ruled out.
    """
    var columns = forward.columns
    var rows = forward.rows
    var middle = columns // 2 if meet else columns
    var gap = abs(rows - columns)
    var aimed = probe.estimate + probe.estimate // 8
    if columns <= SHORT_COLUMNS:
        aimed = min(probe.estimate * SHORT_AIM // 10, probe.estimate + SHORT_REACH)
    # The heuristic at the origin is a lower bound on the distance. With seeds it lands within about a
    # sixth of it on a close pair, where the projection from a few edits strays further.
    var origin = forward_heuristic.h(0, 0)
    var threshold = max(gap, probe.floor + 1, aimed + PROBE_MARGIN)
    if not trusted:
        # As A*PA2 starts: a step past what is known, the heuristic at the origin or the floor.
        threshold = min(threshold, max(gap, probe.floor + 1, origin) + DOUBLING_START)
    if forward_heuristic.seeds > 0:
        if forward_heuristic.chains_well():
            # Matches survive: the origin's bound lies within a few percent of the distance on a
            # close pair, closer than any projection, so start just past it. On a more divergent
            # pair the round dies early, a share of the way across proportional to how far the
            # bound sits above the origin's, and its death estimates the rest (see below).
            threshold = origin + SEED_SLACK
        # Otherwise few matches survive and the bound is one edit a seed, well short of a divergent
        # pair's distance: the projection leads, the bound only floors it.
        threshold = max(threshold, gap, probe.floor + 1, origin + SEED_SLACK)
    var best = Int.MAX
    # Every round re-aims at checkpoints, a later one never below a quarter of the last bound's climb
    # above the origin's past it, so a lowered round that falls short still leaves the next one higher. Only the first may give itself
    # up there: on reads whose errors gather at an end, a later round given up on its own climb would
    # jump past a bound that was enough.
    var first_round = True
    var last_bound = 0
    while True:
        if 2 * (threshold + BAND_COLUMNS) >= rows:
            if give_up_wide:
                return Outcome(-1, middle, rows)
            threshold = max(threshold, columns + rows)
        var forward_round = Round(-1, 0, threshold, -1)
        var backward_round = Round(-1, 0, threshold, -1)
        if meet:

            def half(
                index: Int,
            ) {
                mut forward,
                mut backward,
                mut forward_trail,
                mut backward_trail,
                mut forward_edge,
                mut backward_edge,
                mut forward_round,
                mut backward_round,
                mut forward_heuristic,
                mut backward_heuristic,
                imm threshold,
                imm middle,
                imm columns,
            }:
                if index == 0:
                    forward_round = pruned_distance[record](
                        forward, threshold, middle, forward_trail, forward_edge, forward_heuristic
                    )
                else:
                    backward_round = pruned_distance[record](
                        backward, threshold, columns - middle, backward_trail, backward_edge, backward_heuristic
                    )

            parallelize(half, 2, 2)
        else:
            forward_round = pruned_distance[record](
                forward,
                threshold,
                columns,
                forward_trail,
                forward_edge,
                forward_heuristic,
                first_round,
                first_round,
                0 if first_round else last_bound + max((last_bound - origin) // 4, SEED_SLACK),
            )
        first_round = False

        var found = -1
        var meeting_row = rows
        if not meet:
            found = forward_round.distance
        elif forward_round.reached == middle and backward_round.reached == columns - middle:
            # Join at the middle: the forward score at a row plus the backward score at its mirror.
            var low = max(forward_edge.low_row, rows - backward_edge.high_row)
            var high = min(forward_edge.high_row, rows - backward_edge.low_row)
            for row in range(low, high + 1):
                var total = forward_edge.score(row) + backward_edge.score(rows - row)
                if found < 0 or total < found:
                    found = total
                    meeting_row = row
        # A checkpoint may have lowered the bound, and only a distance within the lowered one is exact.
        var bound = min(forward_round.bound, backward_round.bound) if meet else forward_round.bound
        if found >= 0 and found <= bound:
            return Outcome(found, middle, meeting_row)

        # Grown from the bound the round ended on, which a checkpoint may have lowered: from the one it
        # started on, a round lowered and then failed would jump back to a bound it already knew was loose.
        var next = 2 * bound
        if found >= 0:
            # Still the cost of a real alignment, so it caps every later bound.
            best = min(best, found)
        else:
            # A half pruned every row `reached` columns into its sweep, where the best alignment's
            # score plus heuristic had climbed from the heuristic at the origin past the bound;
            # that climb scaled to the whole width estimates the distance.
            var estimate = Int.MAX
            if forward_round.estimate >= 0:
                # A checkpoint gave the round up, with its own projection.
                estimate = forward_round.estimate
            elif forward_round.reached > 0 and (
                forward_round.reached < middle
                or (not meet and forward_round.reached == middle and forward_heuristic.seeds > 0)
            ):
                # A seeded round that crossed every column and still missed the end climbed past the
                # bound in its last tile: the bound itself is the projection, and the retry aims its
                # margin past it rather than doubling a bound that may be nearly enough.
                estimate = min(estimate, origin + (bound - origin) * columns // forward_round.reached)
            if meet and backward_round.reached < columns - middle and backward_round.reached > 0:
                estimate = min(estimate, origin + (threshold - origin) * columns // backward_round.reached)
            if estimate != Int.MAX:
                next = max(bound + bound // 4, estimate + estimate // 8 + PROBE_MARGIN)
                if forward_heuristic.seeds > 0:
                    # The origin's bound is certain; only the climb above it is estimated.
                    next = max(bound + SEED_SLACK, estimate + (estimate - origin) // 2 + PROBE_MARGIN)
        last_bound = bound
        if found < 0 and forward_heuristic.seeds > 0 and forward_round.reached * SEEDED_TRUST_SHARE < middle:
            # A seeded round that died within its first columns projects from those alone, which on
            # real reads hold their errors gathered at the start: its margin over the origin's bound
            # grows `SEEDED_GROWTH` times instead, as A*PA2's grows, until a round gets far enough in
            # for its death to say where the distance lies.
            next = origin + SEEDED_GROWTH * max(bound - origin, SEED_SLACK)
        if not trusted:
            # A round's estimate comes from where it stopped, which errors gathered at an end can
            # set as far off as the first projection: at most doubling instead.
            next = min(next, 2 * bound)
        threshold = min(next, best)


def edit_distance(first: String, second: String, threads: Optional[Int] = None) raises AlignmentError -> Int:
    """The global edit distance between two DNA sequences over `ACGT`, by bit-parallel sweep.

    Band doubling, as in A*PA2-simple: guess a bound, sweep only the band of cells a path within it
    could cross, and raise the guess until the answer fits under it, which proves it optimal. Close
    sequences therefore cost far less than the whole matrix. With more than one thread, a long pair
    is swept from both ends at once and joined in the middle (see `band_doubling`), which A*PA2 does
    not do. Once the band would cover most of the matrix, the whole matrix is swept instead, tiled
    across `threads` threads, every thread this process may use by default.
    """
    var forward = Profile(first, second)
    if forward.columns == 0 or forward.rows == 0:
        return forward.columns + forward.rows
    var workers = 1
    if forward.columns >= MEET_COLUMNS or forward.columns * forward.rows >= PARALLEL_CELLS:
        # A system call, so only a pair long enough to use more than one thread asks.
        workers = max(threads.or_else(hardware_threads()), 1)
    var probe = two_ended_distance(forward, STEP_TENTHS_DISTANCE, workers)
    if probe.distance >= 0:
        return probe.distance
    # The band's bounds and the seeds' gate are tuned on one front's projection from its first
    # `PROBE_START` edits, so that is the projection handed on, unless the search went far further.
    var fronts = DiagonalFronts()
    var projected = diagonal_transition(forward, PROJECTION_ONLY, fronts)
    var trusted = trusted_projection(probe, projected)
    probe = Probe(-1, projected.estimate, max(probe.floor, projected.floor))
    # The seeds' gate stays on the first projection, which it is tuned on.
    var seeded = forward.columns >= SEED_COLUMNS or (
        projected.estimate >= SEED_EDITS and projected.estimate * SEED_DIVERGENCE <= forward.columns
    )
    # Long pairs only: inexact seeds when the projection already says they pay, else exact ones
    # rebuilt inexact if they chain poorly (see `SeedHeuristic`).
    var long = forward.columns >= INEXACT_COLUMNS
    var inexact = long and projected.estimate * INEXACT_DIVERGENCE >= forward.columns
    var forward_heuristic = SeedHeuristic(forward, inexact, choose=long) if seeded else SeedHeuristic(
        forward.columns, forward.rows
    )
    inexact = forward_heuristic.cost > 1
    # A band narrowed by many chained seeds gains less from a second thread than that thread's half,
    # with its own heuristic to build and a weaker start, costs.
    var meet = workers > 1 and forward.columns >= MEET_COLUMNS and not forward_heuristic.chains_well()
    var backward = Profile(first, second, reverse=True) if meet else Profile(String(), String())
    var forward_trail = Trail()
    var backward_trail = Trail()
    var forward_edge = Edge()
    var backward_edge = Edge()
    var backward_heuristic = SeedHeuristic(backward, inexact) if seeded and meet else SeedHeuristic(
        backward.columns, backward.rows
    )
    var outcome = band_doubling[False](
        forward,
        backward,
        meet,
        True,
        probe,
        forward_trail,
        backward_trail,
        forward_edge,
        backward_edge,
        forward_heuristic,
        backward_heuristic,
        trusted,
    )
    if outcome.distance >= 0:
        return outcome.distance
    return full_distance(forward, workers)


# region Traceback


struct Edge(Movable):
    """A column's scores, read back from its differences, so a row's score costs one lookup and a bit count.

    `bases[k]` is the score at the top of the edge's `k`-th word; a row inside a word adds the
    differences of the word's rows above it.
    """

    var top: Int
    var low_row: Int
    """The first row with a known score: the band's top on this column."""
    var high_row: Int
    """The last row with a known score."""
    var bases: List[Int]
    var plus: List[UInt64]
    var minus: List[UInt64]

    def __init__(out self, words: Int = 64):
        """Room for `words` words of edge, which the edge grows past only when it must."""
        self.top = 0
        self.low_row = 0
        self.high_row = 0
        self.bases = List[Int](capacity=words + 1)
        self.plus = List[UInt64](capacity=words)
        self.minus = List[UInt64](capacity=words)

    def capture(mut self, top: Int, end: Int, anchor: Int, frontier: Frontier, rows: Int):
        """The right edge a round finished on: words `top` to `end` of the frontier, scoring `anchor` at the top."""
        self.top = top
        self.low_row = top * WORD_BITS
        self.high_row = min(end * WORD_BITS, rows)
        self.bases.clear()
        self.plus.clear()
        self.minus.clear()
        var running = anchor
        self.bases.append(running)
        for word in range(top, end):
            var plus = frontier.vertical_plus[word]
            var minus = frontier.vertical_minus[word]
            self.plus.append(plus)
            self.minus.append(minus)
            running += word_value(plus, minus)
            self.bases.append(running)

    def load(mut self, trail: Trail, tile: Int, rows: Int):
        """Reads one tile's left edge from the trail, reusing this edge's buffers.

        A segment may end anywhere from the band's top down to the deepest row the tile to the left
        computed; on the first column, the global border, every row is known.
        """
        self.top = trail.tops[tile]
        var count = trail.ends[tile] - self.top
        var offset = trail.offsets[tile]
        self.low_row = self.top * WORD_BITS
        self.high_row = rows if tile == 0 else min(trail.ends[tile - 1] * WORD_BITS, rows)
        self.bases.clear()
        self.plus.clear()
        self.minus.clear()
        var running = trail.anchors[tile]
        self.bases.append(running)
        for word in range(count):
            var plus = trail.edge_plus[offset + word]
            var minus = trail.edge_minus[offset + word]
            self.plus.append(plus)
            self.minus.append(minus)
            running += word_value(plus, minus)
            self.bases.append(running)

    @always_inline
    def score(self, row: Int) -> Int:
        """The score at `row`, which must lie within `low_row ..= high_row`."""
        if row == self.low_row:
            return self.bases[0]
        var word = (row - 1) // WORD_BITS - self.top
        var bits = row - (self.top + word) * WORD_BITS
        var kept = ALL_ONES if bits == WORD_BITS else (UInt64(1) << UInt64(bits)) - 1
        return self.bases[word] + word_value(self.plus[word] & kept, self.minus[word] & kept)


struct Wavefronts(Movable):
    """The wavefront search's buffers, kept across tiles so a tile allocates nothing.

    Cost `s` keeps its diagonals from offset `lows[s]` from the start's, padding included, at position
    `starts[s]` of the flat buffers on.
    """

    var landed: List[Int]
    var reached: List[Int]
    var how: List[UInt8]
    var starts: List[Int]
    var lows: List[Int]

    def __init__(out self):
        # Room for a typical tile up front: growing from empty cost more than the search itself.
        comptime DIAGONALS = 1024
        comptime COSTS = 64
        self.landed = List[Int](capacity=DIAGONALS)
        self.reached = List[Int](capacity=DIAGONALS)
        self.how = List[UInt8](capacity=DIAGONALS)
        self.starts = List[Int](capacity=COSTS)
        self.lows = List[Int](capacity=COSTS)


@always_inline
def slide(
    first: ImmPointer[UInt8, _], second: ImmPointer[UInt8, _], column: Int, diagonal: Int, first_column: Int
) -> Int:
    """How far back along a diagonal matches carry `column`, eight bases at a time while both have eight.

    The eight bases ending just before the position are loaded from each sequence as one word; their
    XOR is zero in every matching byte, so the matching run is the count of zero bytes from the top.
    """
    var at = column
    var row = at + diagonal
    while at - first_column >= 8 and row >= 8:
        var mismatches = (
            first.unsafe_offset(at - 8).unsafe_bitcast[UInt64]().unsafe_load()
            ^ second.unsafe_offset(row - 8).unsafe_bitcast[UInt64]().unsafe_load()
        )
        if mismatches != 0:
            var run = Int(count_leading_zeros(mismatches)) // 8
            return at - run
        at -= 8
        row -= 8
    while at > first_column and row > 0 and first[unsafe_offset=at - 1] == second[unsafe_offset=row - 1]:
        at -= 1
        row -= 1
    return at


def wavefront_segment(
    profile: Profile,
    edge: Edge,
    first_column: Int,
    end_column: Int,
    end_row: Int,
    score: Int,
    limit: Int,
    mut fronts: Wavefronts,
    mut moves: List[UInt8],
) -> Int:
    """Traces from `(end_column, end_row)`, scoring `score`, back to the tile's left edge by wavefront.

    Diagonal transition run backwards, as A*PA2's trace does: for each cost `s`, the furthest a path
    of `s` edits reaches back along every live diagonal, sliding over matches for free. A path is
    accepted only when it reaches the left edge at a row whose score plus `s` is `score`, so it is
    optimal and joins the trace exactly.

    As in A*PA2, the live diagonals shrink from both ends while the outermost front has fallen more
    than `FRONT_DROP` behind the furthest or already stands on the left edge, and the search gives up
    when half its budget has not carried any front halfway across. Each only ever turns a success
    into a miss, and a miss falls back to recomputing the tile. Appends the segment's moves right to
    left and returns its left-edge row, or -1 on a miss.
    """
    var first = profile.column_codes.unsafe_ptr()
    var second = profile.row_codes.unsafe_ptr()
    var home = end_row - end_column
    # Unreached reads as a column far left of any tile, so every move from it fails its own test
    # unchecked; each stored row carries `TRACE_PADDING` of it either side, so the next row reads its
    # neighbours without range checks.
    comptime UNREACHED = -(1 << 40)
    comptime NONE = Int.MAX

    fronts.landed.clear()
    fronts.reached.clear()
    fronts.how.clear()
    fronts.starts.clear()
    fronts.lows.clear()

    @always_inline
    def ends_here(column: Int, offset: Int, cost: Int) {imm edge, imm first_column, imm score, imm home} -> Bool:
        if column != first_column:
            return False
        var row = column + home + offset
        return row >= edge.low_row and row <= edge.high_row and edge.score(row) + cost == score

    @always_inline
    def row(mut fronts: Wavefronts, low: Int, high: Int) -> Int:
        """Lays out the next cost's row over diagonals `low ..= high`, its padding unreached, and
        returns where the row's diagonal zero would sit, so diagonal `k` is that plus `k`."""
        var start = len(fronts.reached)
        var width = high - low + 1 + 2 * TRACE_PADDING
        fronts.starts.append(start)
        fronts.lows.append(low - TRACE_PADDING)
        fronts.reached.resize(unsafe_uninit_length=start + width)
        fronts.landed.resize(unsafe_uninit_length=start + width)
        fronts.how.resize(unsafe_uninit_length=start + width)
        var zero = start + TRACE_PADDING - low
        var reached = fronts.reached.unsafe_ptr()
        for index in range(1, TRACE_PADDING + 1):
            reached[unsafe_offset=zero + low - index] = UNREACHED
            reached[unsafe_offset=zero + high + index] = UNREACHED
        return zero

    var finish_cost = -1
    var finish_offset = 0
    var zero = row(fronts, 0, 0)
    fronts.landed[zero] = end_column
    fronts.reached[zero] = slide(first, second, end_column, home, first_column)
    fronts.how[zero] = DIAGONAL
    if ends_here(fronts.reached[zero], 0, 0):
        finish_cost = 0
    var low = 0
    var high = 0
    var cost = 1
    var halfway = (first_column + end_column) // 2
    while finish_cost < 0 and cost <= limit:
        low -= 1
        high += 1
        var current_zero = row(fronts, low, high)
        # The previous row is laid out from its own first stored diagonal, unreached outside its live ones.
        var previous = fronts.reached.unsafe_ptr().unsafe_offset(zero)
        var current = fronts.reached.unsafe_ptr().unsafe_offset(current_zero)
        var landed = fronts.landed.unsafe_ptr().unsafe_offset(current_zero)
        var how = fronts.how.unsafe_ptr().unsafe_offset(current_zero)
        zero = current_zero
        var furthest = NONE
        var nearest_column = NONE
        for offset in range(low, high + 1):
            var same = previous[unsafe_offset=offset]
            var below = previous[unsafe_offset=offset - 1]
            var above = previous[unsafe_offset=offset + 1]
            var best = NONE
            var move = DIAGONAL
            var diagonal = home + offset
            # Substitution, along the same diagonal.
            if same > first_column and same + diagonal > 0:
                best = same - 1
            # A base of the first sequence against a gap, from the diagonal below.
            if below > first_column and below - 1 < best:
                best = below - 1
                move = LEFT
            # A base of the second sequence against a gap, from the diagonal above.
            if above + diagonal + 1 > 0 and above < best:
                best = above
                move = UP
            how[unsafe_offset=offset] = move
            if best == NONE:
                landed[unsafe_offset=offset] = UNREACHED
                current[unsafe_offset=offset] = UNREACHED
                continue
            landed[unsafe_offset=offset] = best
            var slid = slide(first, second, best, diagonal, first_column)
            current[unsafe_offset=offset] = slid
            furthest = min(furthest, 2 * slid - offset)
            nearest_column = min(nearest_column, slid)
            if finish_cost < 0 and ends_here(slid, offset, cost):
                finish_cost = cost
                finish_offset = offset
        if finish_cost >= 0:
            break
        if furthest == NONE:
            return -1
        if 2 * cost >= limit and nearest_column > halfway:
            return -1

        @always_inline
        def dropped(offset: Int) {imm current, imm furthest, imm first_column} -> Bool:
            var at = current[unsafe_offset=offset]
            return at <= first_column or 2 * at - offset > furthest + FRONT_DROP

        # Shrink the live diagonals from both ends, marking the dropped ones unreached in place.
        var new_low = low
        var new_high = high
        while new_low < new_high and dropped(new_low):
            current[unsafe_offset=new_low] = UNREACHED
            new_low += 1
        while new_high > new_low and dropped(new_high):
            current[unsafe_offset=new_high] = UNREACHED
            new_high -= 1
        if dropped(new_low):
            return -1
        low = new_low
        high = new_high
        cost += 1
    if finish_cost < 0:
        return -1

    # Walk the costs back down, left to right: once to size the segment, an edit per cost plus the
    # matches each slid over, then again writing its moves from the block's far end, since they go
    # after the moves already traced right to left.
    var length = finish_cost
    var offset = finish_offset
    var at_cost = finish_cost
    while True:
        var index = fronts.starts[at_cost] + offset - fronts.lows[at_cost]
        length += fronts.landed[index] - fronts.reached[index]
        if at_cost == 0:
            break
        var move = fronts.how[index]
        if move == LEFT:
            offset -= 1
        elif move == UP:
            offset += 1
        at_cost -= 1
    var block = len(moves)
    moves.resize(unsafe_uninit_length=block + length)
    var out = moves.unsafe_ptr().unsafe_offset(block + length)
    offset = finish_offset
    at_cost = finish_cost
    while True:
        var index = fronts.starts[at_cost] + offset - fronts.lows[at_cost]
        var run = fronts.landed[index] - fronts.reached[index]
        out = out.unsafe_offset(-run)
        for step in range(run):
            out[unsafe_offset=step] = DIAGONAL
        if at_cost == 0:
            break
        var move = fronts.how[index]
        out = out.unsafe_offset(-1)
        out[] = move
        if move == LEFT:
            offset -= 1
        elif move == UP:
            offset += 1
        at_cost -= 1
    return first_column + home + finish_offset


comptime RECOMPUTE_WORDS = 4
"""Words a tile's recompute first takes above the point it traces from, doubling while that falls short."""


def recomputed_segment(
    profile: Profile, trail: Trail, tile: Int, end_column: Int, end_row: Int, score: Int, mut moves: List[UInt8]
) -> Int:
    """Traces a tile back cell by cell, after sweeping it again and keeping every column's differences.

    The fallback for a tile with more edits than the wavefront search is allowed. As A*PA2 does, it
    first recomputes only a window of rows above the point it traces from, `RECOMPUTE_WORDS` words
    and doubling, rather than the band's whole height (see `window_segment`). Appends the moves
    right to left and returns the left-edge row.
    """
    var top = trail.tops[tile]
    var end_word = min(max(ceildiv(end_row, WORD_BITS), top + 1), trail.ends[tile])
    var window = RECOMPUTE_WORDS
    while True:
        var first_word = max(top, end_word - window)
        var left = window_segment(profile, trail, tile, end_column, end_row, score, first_word, end_word, moves)
        if left >= 0 or first_word == top:
            return left
        window *= 2


def window_segment(
    profile: Profile,
    trail: Trail,
    tile: Int,
    end_column: Int,
    end_row: Int,
    score: Int,
    first_word: Int,
    end_word: Int,
    mut moves: List[UInt8],
) -> Int:
    """`recomputed_segment` over words `[first_word, end_word)` alone, or -1 when they do not hold the path.

    The window's left edge comes from the recorded one, and its top reads `+1` from above, as a
    band's top does, so every score is a real path's and never below the true one. When the
    recomputed score at the traced cell is its exact `score`, a path back to the left edge along
    those scores is optimal; the window falls short when the scores differ, or when the path would
    leave through its top. Below the band's own top, which is the full recompute, neither can happen.
    """
    var first_column = trail.first_columns[tile]
    var top = first_word
    var count = end_word - first_word
    var width = end_column - first_column
    var offset = trail.offsets[tile] + first_word - trail.tops[tile]
    var plus = List[UInt64](length=(width + 1) * count, fill=0)
    var minus = List[UInt64](length=(width + 1) * count, fill=0)
    # The left edge's score at the window's top, carried down from the band's.
    var anchor = trail.anchors[tile]
    for word in range(trail.offsets[tile], offset):
        anchor += word_value(trail.edge_plus[word], trail.edge_minus[word])
    for word in range(count):
        plus[word] = trail.edge_plus[offset + word]
        minus[word] = trail.edge_minus[offset + word]
    for step in range(1, width + 1):
        var horizontal_plus = UInt64(1)
        var horizontal_minus = UInt64(0)
        for word in range(count):
            var vertical_plus = plus[(step - 1) * count + word]
            var vertical_minus = minus[(step - 1) * count + word]
            advance[1](
                horizontal_plus,
                horizontal_minus,
                vertical_plus,
                vertical_minus,
                profile.matches(first_column + step - 1, top + word),
            )
            plus[step * count + word] = vertical_plus
            minus[step * count + word] = vertical_minus

    # Scores at each word's top on every column: the window's top scores the anchor plus one per column.
    var bases = List[Int](length=(width + 1) * (count + 1), fill=0)
    for step in range(width + 1):
        var running = anchor + step
        bases[step * (count + 1)] = running
        for word in range(count):
            running += word_value(plus[step * count + word], minus[step * count + word])
            bases[step * (count + 1) + word + 1] = running

    @always_inline
    def score_at(step: Int, row: Int) {imm bases, imm plus, imm minus, imm count, imm top} -> Int:
        var word = (row - 1) // WORD_BITS - top if row > top * WORD_BITS else 0
        if row == top * WORD_BITS:
            return bases[step * (count + 1)]
        var bits = row - (top + word) * WORD_BITS
        var kept = ALL_ONES if bits == WORD_BITS else (UInt64(1) << UInt64(bits)) - 1
        return bases[step * (count + 1) + word] + word_value(
            plus[step * count + word] & kept, minus[step * count + word] & kept
        )

    var whole = first_word == trail.tops[tile]
    if not whole and score_at(width, end_row) != score:
        return -1
    var lowest = top * WORD_BITS
    var start = len(moves)
    var step = width
    var row = end_row
    var current = score
    while step > 0:
        var column = first_column + step - 1
        if row > lowest:
            var change = 0 if profile.column_codes[column] == profile.row_codes[row - 1] else 1
            if score_at(step - 1, row - 1) + change == current:
                moves.append(DIAGONAL)
                current -= change
                step -= 1
                row -= 1
                continue
        if score_at(step - 1, row) + 1 == current:
            moves.append(LEFT)
            current -= 1
            step -= 1
            continue
        if not whole and row == lowest:
            # The path leaves through the window's top: a taller window must hold it.
            moves.resize(start, 0)
            return -1
        moves.append(UP)
        current -= 1
        row -= 1
    return row


def trace_back(profile: Profile, trail: Trail, start_column: Int, start_row: Int, score: Int, mut moves: List[UInt8]):
    """An optimal path from `(start_column, start_row)`, scoring `score`, back to the origin, as moves right to left.

    One tile at a time from the last: each segment ends at a left-edge row whose recorded score is
    exact, because it plus the segment's cost is the exact score it started from, so the next tile
    starts from a known score. The wavefront search allows a few times the tile's share of the
    score in edits before handing the tile to the recompute.
    """
    var edge = Edge()
    var fronts = Wavefronts()
    var column = start_column
    var row = start_row
    var current = score
    for tile in range(len(trail.first_columns) - 1, -1, -1):
        var first_column = trail.first_columns[tile]
        edge.load(trail, tile, profile.rows)
        var share = ceildiv(score * (column - first_column), max(start_column, 1))
        var limit = max(WAVEFRONT_FLOOR, 3 * share)
        var left = wavefront_segment(profile, edge, first_column, column, row, current, limit, fronts, moves)
        if left < 0:
            left = recomputed_segment(profile, trail, tile, column, row, current, moves)
        # Exact, since it plus the segment's cost is the exact score the segment started from.
        current = edge.score(left)
        column = first_column
        row = left
    # The first column is the border: the rest of the way up is gaps against the second sequence.
    for _ in range(row):
        moves.append(UP)


def edit_alignment(
    first: String, second: String, threads: Optional[Int] = None
) raises AlignmentError -> AlignmentResult:
    """The global edit distance between two DNA sequences over `ACGT`, and an optimal alignment.

    The distance comes from `edit_distance`'s band doubling, recording each tile's left edge in the
    round that succeeds; the alignment is then traced back tile by tile from those edges (see
    `trace_back`). Meeting in the middle, both halves are traced at once from the cell where they
    met. The score is the distance, as `levenshtein_alignment` reports it.
    """
    var forward = Profile(first, second)
    var columns = forward.columns
    var rows = forward.rows
    # Moves right to left, from the corner, or from where the halves met, back to the origin.
    var forward_moves = List[UInt8](capacity=columns + rows)
    # Moves left to right, from where the halves met to the corner; none without meeting.
    var backward_moves = List[UInt8]()
    var distance: Int
    var middle = columns
    var meeting_row = rows
    if columns == 0 or rows == 0:
        distance = columns + rows
        for _ in range(rows):
            forward_moves.append(UP)
        for _ in range(columns):
            forward_moves.append(LEFT)
    else:
        # A near-identical pair: one front, kept whole, settles it before both ends' setup would pay.
        var near = DiagonalFronts()
        var close = diagonal_transition(forward, STEP_TENTHS_ALIGNMENT, near, switch_setup=TWO_ENDED_SETUP)
        if close.distance >= 0:
            trace_diagonals(forward, near, close.distance, forward_moves)
            return gapped_rows(first, second, forward_moves, columns, rows, backward_moves, close.distance)
        # Diagonal transition from both ends, keeping every front, while it is cheaper than a band:
        # where the fronts meet, the path is traced back to the start through the forward fronts
        # and on to the end through the backward ones.
        var first_back = reversed_codes(forward.column_codes, columns, FIRST_SENTINEL)
        var second_back = reversed_codes(forward.row_codes, rows, SECOND_SENTINEL)
        var ahead = FrontPair(PROBE_CEILING // 2 + 1, record=True)
        var behind = FrontPair(PROBE_CEILING // 2 + 1, record=True)
        var workers = max(threads.or_else(0), 0) if columns >= MEET_COLUMNS else 1
        var meeting = two_ended(
            forward, first_back, second_back, STEP_TENTHS_ALIGNMENT, ahead, behind, EDITS_PER_STEP, workers
        )
        if meeting.probe.distance >= 0:
            var middle_column = meeting.column
            var middle_row = meeting.column - meeting.diagonal
            trace_from(
                ahead.history,
                columns,
                rows,
                meeting.forward_score,
                meeting.diagonal,
                middle_column,
                forward_moves,
            )
            # The backward fronts' trace runs from the meeting cell, mirrored, to the end, and its
            # moves right to left over the reversed sequences are the suffix left to right.
            backward_moves = List[UInt8](capacity=columns + rows - middle_column - middle_row)
            trace_from(
                behind.history,
                columns,
                rows,
                meeting.backward_score,
                (columns - rows) - meeting.diagonal,
                columns - middle_column,
                backward_moves,
            )
            return gapped_rows(
                first, second, forward_moves, middle_column, middle_row, backward_moves, meeting.probe.distance
            )
        # The band's bounds and the seeds' gate are tuned on one front's projection from its first
        # `PROBE_START` edits, so that is the projection handed on, unless the search went far further.
        var fronts = DiagonalFronts()
        var projected = diagonal_transition(forward, PROJECTION_ONLY, fronts)
        var trusted = trusted_projection(meeting.probe, projected)
        var probe = Probe(-1, projected.estimate, max(meeting.probe.floor, projected.floor))
        # Asking for the thread count is a system call, so only a pair long enough to split asks.
        # The seeds' gate stays on the first projection, which it is tuned on.
        var seeded = columns >= SEED_COLUMNS or (
            projected.estimate >= SEED_EDITS and projected.estimate * SEED_DIVERGENCE <= columns
        )
        # Long pairs only: inexact seeds when the projection already says they pay, else exact ones
        # rebuilt inexact if they chain poorly (see `SeedHeuristic`).
        var long = columns >= INEXACT_COLUMNS
        var inexact = long and projected.estimate * INEXACT_DIVERGENCE >= columns
        var forward_heuristic = SeedHeuristic(forward, inexact, choose=long) if seeded else SeedHeuristic(columns, rows)
        inexact = forward_heuristic.cost > 1
        # A band narrowed by many chained seeds gains less from a second thread than its half costs.
        # Asking for the thread count is a system call, so only a pair long enough to split asks.
        var meet = (
            not forward_heuristic.chains_well()
            and columns >= MEET_COLUMNS
            and max(threads.or_else(hardware_threads()), 1) > 1
        )
        var backward = Profile(first, second, reverse=True) if meet else Profile(String(), String())
        var forward_trail = Trail(columns)
        var backward_trail = Trail(columns)
        var forward_edge = Edge()
        var backward_edge = Edge()
        var backward_heuristic = SeedHeuristic(backward, inexact) if seeded and meet else SeedHeuristic(
            backward.columns, backward.rows
        )
        var outcome = band_doubling[True](
            forward,
            backward,
            meet,
            False,
            probe,
            forward_trail,
            backward_trail,
            forward_edge,
            backward_edge,
            forward_heuristic,
            backward_heuristic,
            trusted,
        )
        distance = outcome.distance
        if meet:
            middle = outcome.middle
            meeting_row = outcome.meeting_row
            var forward_score = forward_edge.score(meeting_row)
            var backward_score = backward_edge.score(rows - meeting_row)

            def trace(
                index: Int,
            ) {
                imm forward,
                imm backward,
                imm forward_trail,
                imm backward_trail,
                mut forward_moves,
                mut backward_moves,
                imm middle,
                imm meeting_row,
                imm forward_score,
                imm backward_score,
                imm columns,
                imm rows,
            }:
                if index == 0:
                    trace_back(forward, forward_trail, middle, meeting_row, forward_score, forward_moves)
                else:
                    trace_back(
                        backward, backward_trail, columns - middle, rows - meeting_row, backward_score, backward_moves
                    )

            parallelize(trace, 2, 2)
        else:
            trace_back(forward, forward_trail, columns, rows, distance, forward_moves)
        # The forward half runs right to left; the backward half, traced on the reversed sequences
        # right to left, is already the original suffix left to right, move for move.
    return gapped_rows(first, second, forward_moves, middle, meeting_row, backward_moves, distance)


@always_inline
def copy_bytes(destination: MutPointer[UInt8, _], source: ImmPointer[UInt8, _], count: Int):
    """`count` bytes, sixteen at a time while that many are left."""
    var index = 0
    while index + 16 <= count:
        destination.unsafe_offset(index).unsafe_store(source.unsafe_offset(index).unsafe_load[width=16]())
        index += 16
    while index < count:
        destination[unsafe_offset=index] = source[unsafe_offset=index]
        index += 1


@always_inline
def diagonal_run(moves: ImmPointer[UInt8, _], start: Int, end: Int) -> Int:
    """How many moves from `start` on, before `end`, are `DIAGONAL`, eight at a time."""
    var index = start
    while index + 8 <= end:
        var eight = moves.unsafe_offset(index).unsafe_bitcast[UInt64]().unsafe_load()
        if eight != 0:
            return index - start + Int(count_trailing_zeros(eight)) // 8
        index += 8
    while index < end and moves[unsafe_offset=index] == DIAGONAL:
        index += 1
    return index - start


def gapped_rows(
    first: String,
    second: String,
    prefix: List[UInt8],
    middle: Int,
    meeting_row: Int,
    suffix: List[UInt8],
    distance: Int,
) -> AlignmentResult:
    """An alignment written out as the two gapped rows.

    `prefix` holds the moves from the origin to `(middle, meeting_row)` right to left, as the
    traceback appends them, and `suffix` the moves from there to the corner left to right, so
    neither is reversed first. A run of diagonal moves, most of any alignment, copies both
    sequences' bases a block at a time; a gap move writes one base against `-`.
    """
    comptime GAP = UInt8(ord("-"))
    var length = len(prefix) + len(suffix)
    var top_row = List[UInt8](capacity=length)
    var bottom_row = List[UInt8](capacity=length)
    top_row.resize(unsafe_uninit_length=length)
    bottom_row.resize(unsafe_uninit_length=length)
    var first_bytes = first.unsafe_ptr()
    var second_bytes = second.unsafe_ptr()
    var top = top_row.unsafe_ptr()
    var bottom = bottom_row.unsafe_ptr()

    # The prefix, from its last move back to its first, fills the rows from `len(prefix)` down.
    var moves = prefix.unsafe_ptr()
    var count = len(prefix)
    var column = middle
    var row = meeting_row
    var at = count
    var index = 0
    while index < count:
        var run = diagonal_run(moves, index, count)
        if run > 0:
            at -= run
            column -= run
            row -= run
            copy_bytes(top.unsafe_offset(at), first_bytes.unsafe_offset(column), run)
            copy_bytes(bottom.unsafe_offset(at), second_bytes.unsafe_offset(row), run)
            index += run
            continue
        var move = moves[unsafe_offset=index]
        at -= 1
        if move == LEFT:
            column -= 1
            top[unsafe_offset=at] = first_bytes[unsafe_offset=column]
            bottom[unsafe_offset=at] = GAP
        else:
            row -= 1
            top[unsafe_offset=at] = GAP
            bottom[unsafe_offset=at] = second_bytes[unsafe_offset=row]
        index += 1

    # The suffix, first move first, fills the rest.
    var later = suffix.unsafe_ptr()
    count = len(suffix)
    column = middle
    row = meeting_row
    at = len(prefix)
    index = 0
    while index < count:
        var run = diagonal_run(later, index, count)
        if run > 0:
            copy_bytes(top.unsafe_offset(at), first_bytes.unsafe_offset(column), run)
            copy_bytes(bottom.unsafe_offset(at), second_bytes.unsafe_offset(row), run)
            at += run
            column += run
            row += run
            index += run
            continue
        if later[unsafe_offset=index] == LEFT:
            top[unsafe_offset=at] = first_bytes[unsafe_offset=column]
            bottom[unsafe_offset=at] = GAP
            column += 1
        else:
            top[unsafe_offset=at] = GAP
            bottom[unsafe_offset=at] = second_bytes[unsafe_offset=row]
            row += 1
        at += 1
        index += 1
    return AlignmentResult(Int32(distance), String(unsafe_from_utf8=top_row), String(unsafe_from_utf8=bottom_row))


# endregion Traceback
