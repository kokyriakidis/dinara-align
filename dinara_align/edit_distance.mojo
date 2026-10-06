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

A pair too divergent for a band is swept whole. Each pair runs on one thread; `edit_distances` and
`edit_alignments` in `api` spread a batch over threads a pair at a time.

Unit costs only: a substitution, an insertion and a deletion each cost one, so this is the
distance `levenshtein_alignment` returns, without the alignment.
"""

from std.bit import count_leading_zeros, count_trailing_zeros, pop_count
from std.math import ceildiv, clamp, sqrt
from std.sys import inlined_assembly, simd_width_of
from std.sys.info import CompilationTarget
from std.sys.intrinsics import PrefetchOptions, prefetch

from .alignment import AlignmentResult
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


@inline(.always)
def opaque[width: Int](value: SIMD[DType.uint64, width]) -> SIMD[DType.uint64, width]:
    """`value`, unchanged, hidden from the optimizer so it cannot fold it back into a longer chain.

    Only up to four lanes: eight are bound by how many vector operations issue, not by the chain,
    and there the folded form is the cheaper one. Two lanes sit in a vector register, which the
    constraint names per target: `w` for NEON, `x` for SSE; elsewhere the value goes through as is.
    """
    comptime if width == 1:
        return inlined_assembly["", SIMD[DType.uint64, width], constraints="=r,0", has_side_effect=False](value)
    elif width == 2:
        comptime if CompilationTarget.has_neon():
            return inlined_assembly["", SIMD[DType.uint64, width], constraints="=w,0", has_side_effect=False](value)
        elif CompilationTarget.is_x86():
            return inlined_assembly["", SIMD[DType.uint64, width], constraints="=x,0", has_side_effect=False](value)
        else:
            return value
    elif width > 4:
        return value
    else:
        comptime half = width // 2
        var low = opaque[half](value.slice[half, offset=0]())
        var high = opaque[half](value.slice[half, offset=half]())
        return rebind[SIMD[DType.uint64, width]](low.join(high))


@inline(.always)
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
    # Myers' horizontal mask is `((carried ^ vp) | equal)`; both its uses below are rewritten on
    # `carried` directly, `horizontal | vp` as `carried | equal | vp` and `vp & horizontal` as
    # `(vp & ~carried) | (vp & equal)`, which takes two operations off the path one column waits on.
    var carried = (equal & vertical_plus) + vertical_plus
    var plus = vertical_minus | ~(carried | (equal | vertical_plus))
    var minus = (vertical_plus & ~carried) | (vertical_plus & equal)
    var plus_out = plus >> (WORD_BITS - 1)
    var minus_out = minus >> (WORD_BITS - 1)
    var plus_shifted = plus << 1
    # `~(crossing | plus)` with the incoming `+` folded into `crossing`, which is known early.
    var unblocked = opaque(~(crossing | horizontal_plus))
    minus = opaque((minus << 1) | horizontal_minus)
    plus = plus_shifted | horizontal_plus
    horizontal_plus = plus_out
    horizontal_minus = minus_out
    vertical_plus = minus | (unblocked & ~plus_shifted)
    vertical_minus = plus & crossing


struct Sweep(ImplicitlyCopyable, TrivialRegisterPassable):
    """Pointers to the profile and the frontier, which the caller's lists own and outlive.

    The profile holds both sequences as bit planes, so a word of matches is two XORs and an AND. A
    base of the first sequence becomes two whole-word masks, all ones where its code has that bit;
    the second is packed 64 bases to a word with each plane stored negated, so a base equals a
    column's base exactly where both planes XOR to ones.

    The frontier holds the differences along the two edges still being computed: one horizontal
    difference per column, in bit zero, along the bottom of the rows done so far, and one word of
    vertical differences per row word, down the right edge of the columns done so far.

    Plain pointers rather than the lists themselves, so a sweep can write the frontier it was handed
    while it reads the profile.
    """

    var column_low: ImmPointer[UInt64, ImmUntrackedOrigin]
    var column_high: ImmPointer[UInt64, ImmUntrackedOrigin]
    var column_extra: ImmPointer[UInt64, ImmUntrackedOrigin]
    var row_low: ImmPointer[UInt64, ImmUntrackedOrigin]
    var row_high: ImmPointer[UInt64, ImmUntrackedOrigin]
    var row_extra: ImmPointer[UInt64, ImmUntrackedOrigin]
    """The third plane, read only by the kernels a profile with symbols past `ACGT` takes."""
    var horizontal_plus: MutPointer[UInt64, MutUntrackedOrigin]
    var horizontal_minus: MutPointer[UInt64, MutUntrackedOrigin]
    var vertical_plus: MutPointer[UInt64, MutUntrackedOrigin]
    var vertical_minus: MutPointer[UInt64, MutUntrackedOrigin]

    def __init__(
        out self,
        column_low: List[UInt64],
        column_high: List[UInt64],
        column_extra: List[UInt64],
        row_low: List[UInt64],
        row_high: List[UInt64],
        row_extra: List[UInt64],
        mut horizontal_plus: List[UInt64],
        mut horizontal_minus: List[UInt64],
        mut vertical_plus: List[UInt64],
        mut vertical_minus: List[UInt64],
    ):
        # Past the padding, so column `c` of the matrix is element `c` here.
        self.column_low = column_low.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]().unsafe_offset(COLUMN_PADDING)
        self.column_high = (
            column_high.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]().unsafe_offset(COLUMN_PADDING)
        )
        self.column_extra = (
            column_extra.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]().unsafe_offset(COLUMN_PADDING)
        )
        self.row_low = row_low.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]()
        self.row_high = row_high.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]()
        self.row_extra = row_extra.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]()
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

    @inline(.always)
    def run[extended: Bool](self, word: Int, first_column: Int, end_column: Int):
        """One word through `[first_column, end_column)`, its vertical differences held in registers.

        Only the horizontal edge goes through memory, and each column's is a different address, so
        consecutive columns do not wait on a store being read back.
        """
        var vp = self.vertical_plus[unsafe_offset=word]
        var vm = self.vertical_minus[unsafe_offset=word]
        var row_low = self.row_low[unsafe_offset=word]
        var row_high = self.row_high[unsafe_offset=word]
        var row_extra = self.row_extra[unsafe_offset=word] if extended else UInt64(0)
        for column in range(first_column, end_column):
            var hp = self.horizontal_plus[unsafe_offset=column]
            var hm = self.horizontal_minus[unsafe_offset=column]
            var matches = (self.column_low[unsafe_offset=column] ^ row_low) & (
                self.column_high[unsafe_offset=column] ^ row_high
            )
            comptime if extended:
                matches &= self.column_extra[unsafe_offset=column] ^ row_extra
            advance[1](hp, hm, vp, vm, matches)
            self.horizontal_plus[unsafe_offset=column] = hp
            self.horizontal_minus[unsafe_offset=column] = hm
        self.vertical_plus[unsafe_offset=word] = vp
        self.vertical_minus[unsafe_offset=word] = vm

    def words(self, extended: Bool, first_word: Int, end_word: Int, first_column: Int, end_column: Int):
        """Words `[first_word, end_word)` through columns `[first_column, end_column)`, matching on the
        third plane when `extended`, as a profile with symbols past `ACGT` needs (see `words_matching`)."""
        if extended:
            self.words_matching[True](first_word, end_word, first_column, end_column)
        else:
            self.words_matching[False](first_word, end_word, first_column, end_column)

    # Out of line, each kernel its own function: inlined side by side, or chosen per band round, the
    # two slowed the bases' one by 1 to 4%.
    @inline(.never)
    def words_matching[extended: Bool](self, first_word: Int, end_word: Int, first_column: Int, end_column: Int):
        """`words`, its match test fixed by `extended`.

        Full groups of `LANES` take the staggered vector sweep when the span is wide enough for it.
        Of what is left, four words take a narrow vector, five to seven a narrow vector with the rest
        in scalar registers beside it, two or three scalar registers alone, and one word goes on its own.
        """
        var word = first_word
        if end_column - first_column >= 2 * LANES:
            while word + LANES <= end_word:
                self.block[LANES, extended](word, first_column, end_column)
                word += LANES
            # Five to seven words left: a narrow vector with the rest in scalar registers in its shadow.
            var left = end_word - word
            if left == 7:
                self.hybrid_block[NARROW_LANES, 3, extended](word, first_column, end_column)
                word += 7
            elif left == 6:
                self.hybrid_block[NARROW_LANES, 2, extended](word, first_column, end_column)
                word += 6
            elif left == 5:
                self.hybrid_block[NARROW_LANES, 1, extended](word, first_column, end_column)
                word += 5
            elif left == 4:
                self.block[NARROW_LANES, extended](word, first_column, end_column)
                word += NARROW_LANES
            # Two or three words left run side by side in scalar registers, which beats a vector
            # this narrow: its chain waits two cycles an operation, a scalar one.
            if end_word - word == 3:
                self.scalar_block[3, extended](word, first_column, end_column)
                word += 3
            elif end_word - word == 2:
                self.scalar_block[2, extended](word, first_column, end_column)
                word += 2
        while word < end_word:
            self.run[extended](word, first_column, end_column)
            word += 1

    def block[lanes: Int, extended: Bool](self, first_word: Int, first_column: Int, end_column: Int):
        """`lanes` words through `[first_column, end_column)` in one vector, staggered (see `VectorGroup`)."""
        var group = VectorGroup[lanes, extended](self, first_word, first_column, end_column)
        stagger(group)
        group.finish()

    def scalar_block[lanes: Int, extended: Bool](self, first_word: Int, first_column: Int, end_column: Int):
        """`lanes` words through `[first_column, end_column)` in scalar registers, staggered (see `ScalarGroup`)."""
        var group = ScalarGroup[lanes, extended](self, first_word, first_column, end_column)
        stagger(group)
        group.finish()

    def hybrid_block[lanes: Int, below: Int, extended: Bool](self, first_word: Int, first_column: Int, end_column: Int):
        """`lanes` words in a vector and the `below` words under them in scalar registers, in one loop.

        A narrow vector's chain leaves the integer units idle, so the scalar words run in its
        shadow, a few columns behind so the differences the vector sends down are already stored.
        """
        var top = VectorGroup[lanes, extended](self, first_word, first_column, end_column)
        var bottom = ScalarGroup[below, extended](self, first_word + lanes, first_column, end_column)
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


struct VectorGroup[lanes: Int, extended: Bool](Staggered, TrivialRegisterPassable):
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
    var row_extra: SIMD[DType.uint64, Self.lanes]
    var vertical_plus: SIMD[DType.uint64, Self.lanes]
    var vertical_minus: SIMD[DType.uint64, Self.lanes]
    var horizontal_plus: SIMD[DType.uint64, Self.lanes]
    var horizontal_minus: SIMD[DType.uint64, Self.lanes]
    var lane_columns: SIMD[DType.int, Self.lanes]

    @inline(.always)
    def __init__(out self, sweep: Sweep, first_word: Int, first_column: Int, end_column: Int):
        self.sweep = sweep
        self.first_word = first_word
        self.first = first_column
        self.end = end_column
        self.row_low = SIMD[DType.uint64, Self.lanes]()
        self.row_high = SIMD[DType.uint64, Self.lanes]()
        self.row_extra = SIMD[DType.uint64, Self.lanes]()
        self.vertical_plus = SIMD[DType.uint64, Self.lanes]()
        self.vertical_minus = SIMD[DType.uint64, Self.lanes]()
        self.horizontal_plus = SIMD[DType.uint64, Self.lanes]()
        self.horizontal_minus = SIMD[DType.uint64, Self.lanes]()
        self.lane_columns = SIMD[DType.int, Self.lanes]()
        comptime for lane in range(Self.lanes):
            var word = first_word + Self.lanes - 1 - lane
            self.row_low[lane] = sweep.row_low[unsafe_offset=word]
            self.row_high[lane] = sweep.row_high[unsafe_offset=word]
            comptime if Self.extended:
                self.row_extra[lane] = sweep.row_extra[unsafe_offset=word]
            self.vertical_plus[lane] = sweep.vertical_plus[unsafe_offset=word]
            self.vertical_minus[lane] = sweep.vertical_minus[unsafe_offset=word]
            self.lane_columns[lane] = lane + 1

    @inline(.always)
    def width(self) -> Int:
        return Self.lanes

    @inline(.always)
    def first_column(self) -> Int:
        return self.first

    @inline(.always)
    def end_column(self) -> Int:
        return self.end

    @inline(.always)
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
        comptime if Self.extended:
            matches &= (
                self.sweep.column_extra.unsafe_offset(offset + 1).unsafe_load[width=Self.lanes]() ^ self.row_extra
            )
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

    @inline(.always)
    def finish(self):
        comptime for lane in range(Self.lanes):
            var word = self.first_word + Self.lanes - 1 - lane
            self.sweep.vertical_plus[unsafe_offset=word] = self.vertical_plus[lane]
            self.sweep.vertical_minus[unsafe_offset=word] = self.vertical_minus[lane]


struct ScalarGroup[lanes: Int, extended: Bool](Staggered):
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
    var row_extra: Array[UInt64, Self.lanes]
    var vertical_plus: Array[UInt64, Self.lanes]
    var vertical_minus: Array[UInt64, Self.lanes]
    var horizontal_plus: Array[UInt64, Self.lanes]
    var horizontal_minus: Array[UInt64, Self.lanes]

    @inline(.always)
    def __init__(out self, sweep: Sweep, first_word: Int, first_column: Int, end_column: Int):
        self.sweep = sweep
        self.first_word = first_word
        self.first = first_column
        self.end = end_column
        self.row_low = Array[UInt64, Self.lanes](fill=0)
        self.row_high = Array[UInt64, Self.lanes](fill=0)
        self.row_extra = Array[UInt64, Self.lanes](fill=0)
        self.vertical_plus = Array[UInt64, Self.lanes](fill=0)
        self.vertical_minus = Array[UInt64, Self.lanes](fill=0)
        self.horizontal_plus = Array[UInt64, Self.lanes](fill=0)
        self.horizontal_minus = Array[UInt64, Self.lanes](fill=0)
        comptime for j in range(Self.lanes):
            self.row_low[j] = sweep.row_low[unsafe_offset=first_word + j]
            self.row_high[j] = sweep.row_high[unsafe_offset=first_word + j]
            comptime if Self.extended:
                self.row_extra[j] = sweep.row_extra[unsafe_offset=first_word + j]
            self.vertical_plus[j] = sweep.vertical_plus[unsafe_offset=first_word + j]
            self.vertical_minus[j] = sweep.vertical_minus[unsafe_offset=first_word + j]

    @inline(.always)
    def width(self) -> Int:
        return Self.lanes

    @inline(.always)
    def first_column(self) -> Int:
        return self.first

    @inline(.always)
    def end_column(self) -> Int:
        return self.end

    @inline(.always)
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
            comptime if Self.extended:
                matches &= self.sweep.column_extra[unsafe_offset=column] ^ self.row_extra[j]
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

    @inline(.always)
    def finish(self):
        comptime for j in range(Self.lanes):
            self.sweep.vertical_plus[unsafe_offset=self.first_word + j] = self.vertical_plus[j]
            self.sweep.vertical_minus[unsafe_offset=self.first_word + j] = self.vertical_minus[j]


@inline(.always)
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


@inline(.always)
def masked_step[G: Staggered](mut group: G, offset: Int):
    """A masked step, or none when the group has not started or has already finished."""
    if offset >= group.first_column() - group.width() and offset < group.end_column() - 1:
        group.step[True](offset)


@inline(.always)
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


@inline(.always)
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


def symbol_codes(
    first: String, second: String, mut column_codes: List[UInt8], mut row_codes: List[UInt8]
) raises AlignmentError:
    """Both sequences' codes when some byte is not a base: `ACGT` zero to three, every other byte the
    next free code, in the order the bytes first appear; past four such symbols, three bits run out."""
    comptime UNSEEN = Int16(-1)
    var table = List[Int16](length=256, fill=UNSEEN)
    table[ord("A")] = 0
    table[ord("C")] = 1
    table[ord("G")] = 2
    table[ord("T")] = 3
    var next_code = 4

    @inline(.always)
    def code(byte: UInt8, mut table: List[Int16], mut next_code: Int) raises AlignmentError -> UInt8:
        if table[Int(byte)] == UNSEEN:
            if next_code == 8:
                raise AlignmentError(
                    ErrorKind.UNKNOWN_SYMBOL, "bit-parallel edit distance takes ACGT and at most four other symbols"
                )
            table[Int(byte)] = Int16(next_code)
            next_code += 1
        return UInt8(table[Int(byte)])

    for byte in first.as_bytes():
        column_codes.append(code(byte, table, next_code))
    for byte in second.as_bytes():
        row_codes.append(code(byte, table, next_code))


struct Profile(Movable):
    """Both sequences as codes, and once a band needs them, as the bit planes `Sweep` reads; see
    `Sweep` for the encoding."""

    var columns: Int
    var rows: Int
    var words: Int
    var column_low: List[UInt64]
    var column_high: List[UInt64]
    var column_extra: List[UInt64]
    var row_low: List[UInt64]
    var row_high: List[UInt64]
    var row_extra: List[UInt64]
    var extended: Bool
    """Whether either sequence holds a symbol past `ACGT`, coded from four up, so a match also needs
    the codes' third bit, `column_extra` and `row_extra`; the seeds, packed two bits a base, stay out."""
    var column_codes: List[UInt8]
    """The first sequence as codes, which the traceback compares base by base."""
    var row_codes: List[UInt8]
    """The second sequence as codes."""

    def __init__(out self, first: String, second: String) raises AlignmentError:
        """Both sequences as codes: `A`, `C`, `G` and `T` zero to three, and up to four other bytes the
        codes four to seven, in the order they first appear, each matching only itself."""
        self.columns = first.byte_length()
        self.rows = second.byte_length()
        self.words = ceildiv(self.rows, WORD_BITS)
        self.column_low = List[UInt64]()
        self.column_high = List[UInt64]()
        self.column_extra = List[UInt64]()
        self.row_low = List[UInt64]()
        self.row_high = List[UInt64]()
        self.row_extra = List[UInt64]()
        self.extended = not all_bases(first) or not all_bases(second)
        if self.extended:
            self.column_codes = List[UInt8](capacity=self.columns + CODE_PADDING)
            self.row_codes = List[UInt8](capacity=self.rows + CODE_PADDING)
            symbol_codes(first, second, self.column_codes, self.row_codes)
            for _ in range(CODE_PADDING):
                self.column_codes.append(FIRST_SENTINEL)
                self.row_codes.append(SECOND_SENTINEL)
            return

        # Every byte is `A`, `C`, `G` or `T`, whose ASCII bits give the code directly: bit 2 is set
        # for `G` and `T`, the code's high bit, and bit 1 differs from bit 2 for `C` and `T`, its low
        # bit. Only the codes are built here; the planes wait for a band (see `build_planes`).
        self.column_codes = List[UInt8](capacity=self.columns + CODE_PADDING)
        self.column_codes.resize(unsafe_uninit_length=self.columns)
        var first_bytes = first.unsafe_ptr()
        var column_codes = self.column_codes.unsafe_ptr()
        comptime CHUNK = 16
        var column = 0
        while column + CHUNK <= self.columns:
            var bytes = first_bytes.unsafe_offset(column).unsafe_load[width=CHUNK]()
            column_codes.unsafe_offset(column).unsafe_store((((bytes >> 1) ^ (bytes >> 2)) & 1) | ((bytes >> 1) & 2))
            column += CHUNK
        while column < self.columns:
            var byte = first_bytes[unsafe_offset=column]
            column_codes[unsafe_offset=column] = (((byte >> 1) ^ (byte >> 2)) & 1) | ((byte >> 1) & 2)
            column += 1

        comptime ONES = UInt64(0x0101010101010101)
        var second_bytes = second.unsafe_ptr()
        self.row_codes = List[UInt8](capacity=self.rows + CODE_PADDING)
        self.row_codes.resize(unsafe_uninit_length=self.rows)
        var row_codes = self.row_codes.unsafe_ptr()
        var row = 0
        while row + 8 <= self.rows:
            var eight = second_bytes.unsafe_offset(row).unsafe_bitcast[UInt64]().unsafe_load()
            var low_bits = ((eight >> 1) ^ (eight >> 2)) & ONES
            var high_bits = (eight >> 2) & ONES
            row_codes.unsafe_offset(row).unsafe_bitcast[UInt64]().unsafe_store(low_bits | (high_bits << 1))
            row += 8
        while row < self.rows:
            var byte = second_bytes[unsafe_offset=row]
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
            high.unsafe_offset(column).unsafe_store(UInt64(0) - ((codes >> 1) & 1).cast[DType.uint64]())
            column += CHUNK
        while column < self.columns:
            var code = column_codes[unsafe_offset=column]
            low[unsafe_offset=column] = UInt64(0) - (code & 1).cast[DType.uint64]()
            high[unsafe_offset=column] = UInt64(0) - ((code >> 1) & 1).cast[DType.uint64]()
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
            self.row_high[row // WORD_BITS] |= (((code >> 1) & 1) ^ 1).cast[DType.uint64]() << shift
            row += 1
        if not self.extended:
            return
        # The codes' third bit, as the first two: whole-word masks per column, and per row stored
        # negated. The eight-at-a-time row packing above reads each byte's two low bits alone, so
        # codes from four up leave the first two planes exact.
        self.column_extra = List[UInt64](length=padded, fill=0)
        self.row_extra = List[UInt64](length=self.words, fill=0)
        var extra = self.column_extra.unsafe_ptr().unsafe_offset(COLUMN_PADDING)
        for index in range(self.columns):
            extra[unsafe_offset=index] = UInt64(0) - ((column_codes[unsafe_offset=index] >> 2) & 1).cast[DType.uint64]()
        for index in range(self.rows):
            var bit = (((row_codes[unsafe_offset=index] >> 2) & 1) ^ 1).cast[DType.uint64]()
            self.row_extra[index // WORD_BITS] |= bit << UInt64(index % WORD_BITS)


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
            profile.column_extra,
            profile.row_low,
            profile.row_high,
            profile.row_extra,
            self.horizontal_plus,
            self.horizontal_minus,
            self.vertical_plus,
            self.vertical_minus,
        )

    def down_right_edge(self, rows: Int) -> Int:
        """How much the score grows from the top row to the bottom row, down the right edge.

        Rows past the end of the last word only padded it and are masked out.
        """
        var words = len(self.vertical_plus)
        var total = 0
        for index in range(words):
            var plus = self.vertical_plus[index]
            var minus = self.vertical_minus[index]
            if index == words - 1 and rows % WORD_BITS != 0:
                var kept = (UInt64(1) << UInt64(rows % WORD_BITS)) - 1
                plus &= kept
                minus &= kept
            total += word_value(plus, minus)
        return total


def full_distance(mut profile: Profile) -> Int:
    """The whole matrix."""
    profile.build_planes()
    var frontier = Frontier(profile.columns, profile.words)
    var sweep = frontier.sweep(profile)
    sweep.words(profile.extended, 0, profile.words, 0, profile.columns)

    # The top-right corner is `columns`; walking down the right edge adds each vertical difference.
    return profile.columns + frontier.down_right_edge(profile.rows)


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
    """The bound the round ended on, which its checkpoints may have lowered (see `Band.check`)."""
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

comptime INEXACT_CHAINED = 40 if simd_width_of[DType.uint64]() <= 2 else (
    30 if simd_width_of[DType.uint64]() <= 4 else 20
)
"""Percent of the exact seeds chained from the origin below which the seeds are rebuilt to match within
one edit. Where that pays depends on the machine: the band is vector work and the seeds' setup scalar,
memory-bound lookups, so the wider the vectors beside the scalar core, the less the band inexact seeds
narrow is worth beside the setup they cost. On an M2, two 64-bit lanes to a NEON register, 40 percent,
about one edit in fifteen bases on real reads; on a Skylake-X at a fixed 3.3 GHz, eight to an AVX-512
register, its seeds' setup twice as slow as the M2's and its band no slower, 20 percent, short of which
the long reads ran faster on exact seeds. Four lanes, AVX2, are between the two and unmeasured."""

comptime INEXACT_CHAINED_SPREAD = 20
"""`INEXACT_CHAINED` for a pair whose two projections agree (see `trusted_projection`), its errors
spread along it as mutated sequences' are, not gathered at an end as real reads' are. Spread errors
break fewer exact seeds' chains than gathered ones for the same narrowing of the band, so exact seeds
still pay down to this share: on the M2, 100 kbp pairs at 6 to 8% aligned in 3.7 to 4.3 ms on exact
seeds against 5.7 ms rebuilt, where the real reads still wanted `INEXACT_CHAINED`'s 40 percent."""

comptime HALF_BITS = INEXACT_LENGTH
"""Bits in half an inexact seed's two-bit code: a one-edit match matches one half exactly, so each half
indexes a table of this many bits."""

comptime ENTRY_BITS = 32
"""Bits below an inexact seed's index in its half tables' entries, its two-bit code there: one load
gives a lookup both."""

comptime SCAN_BATCH = 16
"""Rows whose half tables' buckets are looked up together before any is searched (see
`inexact_scan`): enough misses in flight at once to hide most of their wait."""

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


@inline(.always)
def differing(first: UInt64, second: UInt64) -> UInt64:
    """The bases two two-bit codes disagree on, each as the low bit of its pair."""
    var mismatched = first ^ second
    return (mismatched | (mismatched >> 1)) & UInt64(0x5555555555555555)


@inline(.always)
def covers(prefix: UInt64, suffix: UInt64) -> Bool:
    """Whether the leading bases one comparison agrees on and the trailing bases another agrees on,
    both of the same bases, cover all of them between them.

    The leading run reaches past the trailing run's start when every base the first comparison
    differs on lies below the lowest the second differs on, an unsigned comparison of the first's
    differing bits against the second's lowest, with no count of leading or trailing zeros.
    """
    return suffix == 0 or prefix < (suffix & (0 - suffix))


@inline(.always)
def ends_within_one_edit(shorter: UInt64, level: UInt64, longer: UInt64, code: UInt64) -> Bool:
    """`within_one_edit` for the windows ending a seed's length on from `level`'s start, which start
    at `shorter`, `level` and `longer`, less any whose left half matches the seed's too: that one is
    the left half's lookup's to take."""
    var left = code >> UInt64(HALF_BITS)
    return within_one_edit(
        shorter >> 2,
        level,
        longer,
        level,
        code,
        left != shorter >> UInt64(HALF_BITS),
        left != level >> UInt64(HALF_BITS),
        left != longer >> UInt64(HALF_BITS),
    )


@inline(.always)
def first_reachable(entries: ImmPointer[UInt64, _], start: Int, end: Int, seed: Int) -> Int:
    """The first of a bucket's entries, `[start, end)` in seed order, of `seed` or a later seed: a
    step at a time in a small bucket, by bisection in a large one, a repeat's."""
    comptime LINEAR = 8
    var low = start
    var high = end
    while high - low > LINEAR:
        var middle = (low + high) // 2
        if Int(entries[unsafe_offset=middle] >> UInt64(ENTRY_BITS)) < seed:
            low = middle + 1
        else:
            high = middle
    while low < high and Int(entries[unsafe_offset=low] >> UInt64(ENTRY_BITS)) < seed:
        low += 1
    return low


@inline(.always)
def within_one_edit(
    shorter: UInt64,
    level: UInt64,
    head: UInt64,
    tail: UInt64,
    code: UInt64,
    shorter_open: Bool = True,
    level_open: Bool = True,
    longer_open: Bool = True,
) -> Bool:
    """Whether a seed's code is within one edit of any of three windows `try_windows` would keep:
    `shorter`, one base shorter than the seed, `level`, as long, or the window one longer, whose
    first and last `INEXACT_LENGTH` bases are `head` and `tail`, each only while open.

    The three tests folded together, with no branch to mispredict: an equal-length window within
    one substitution differs on at most one base, a shorter one is the seed less a base when the
    seed's first and last bases but one cover it between them, and a longer one the seed plus a base
    when its first and last `INEXACT_LENGTH` bases cover the seed. An exact match leaves out the
    shorter and longer windows, as `try_windows` does.
    """
    comptime K = INEXACT_LENGTH
    comptime SHORTER = (UInt64(1) << UInt64(2 * K - 2)) - 1
    var whole = differing(level, code)
    var substituted = (whole & (whole - 1)) == 0
    var deleted = covers(differing(shorter, code >> 2), differing(shorter, code & SHORTER))
    var inserted = covers(differing(head, code), differing(tail, code))
    return (substituted & level_open) | (((deleted & shorter_open) | (inserted & longer_open)) & (whole != 0))


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

    The heuristic never overestimates the cost to the end, and changes by at most `cost` a row down a
    column or up it (see `climb`), so the band's pruning stays exact and its jumps stay valid. With no
    seeds it is the plain gap to the end's diagonal.
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
        out self,
        profile: Profile,
        inexact: Bool = False,
        choose: Bool = False,
        cutoff: Int = INEXACT_CHAINED,
    ):
        """Seeds of the profile's first sequence, matched exactly in its second, or with `inexact`
        within one edit. With `choose`, exact seeds, rebuilt inexact when fewer than `cutoff` percent of
        them chain from the origin.

        A match is kept only if a path from its start crosses the next `LOOKAHEAD_SEEDS` seeds for
        less than they would cost unmatched (see `worth_keeping`).
        """
        self = Self(profile.columns, profile.rows)
        self.build(profile, inexact)
        if choose and not inexact and self.seeds > 0:
            var chained = self.seeds - self.h(0, 0)
            if chained * 100 < cutoff * self.seeds:
                self.build(profile, True)

    def build(mut self, profile: Profile, inexact: Bool):
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
        for seed in range(self.seeds - 1, -1, -1):
            var column = seed * self.length
            var potential = self.potential(column)
            var end_potential = potential - self.cost
            for slot in range(Int(firsts[seed]), Int(firsts[seed + 1])):
                var start_row = Int(rows_by_seed[slot])
                var end_row = Int(ends_by_seed[slot]) if inexact else start_row + SEED_LENGTH
                var match_cost = Int(costs_by_seed[slot]) if inexact else 0
                if not self.worth_keeping(first, second, seed, start_row, end_row, match_cost, leftmost, fronts):
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

    # Out of line: inlined into `build`, it slows the pruning loop there by a few percent.
    @inline(.never)
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
        var chained = List[Int32](length=self.seeds, fill=-1)
        var slots = table.unsafe_ptr()
        var chain = chained.unsafe_ptr()
        for seed in range(self.seeds):
            var code = 0
            for offset in range(SEED_LENGTH):
                code = (code << 2) | Int(first[unsafe_offset=seed * SEED_LENGTH + offset])
            var slot = Int((UInt64(code) * 0x9E3779B97F4A7C15) >> UInt64(64 - bits))
            while slots[unsafe_offset=slot] != EMPTY and Int(slots[unsafe_offset=slot] >> 32) != code:
                slot = (slot + 1) & (size - 1)
            var held = slots[unsafe_offset=slot]
            chain[unsafe_offset=seed] = Int32(held & 0xFFFFFFFF) if held != EMPTY else -1
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
                seed = Int(chain[unsafe_offset=seed])

    # Out of line: inlined into `build`, it slows the pruning loop there by a few percent.
    @inline(.never)
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

        Within a bucket the seeds go in order, and only a range of them can be found at any row with
        their chain reaching the end (see `reachable_seeds`), so each lookup tries that range alone.
        """
        comptime HALF = INEXACT_LENGTH // 2
        comptime HALF_MASK = (1 << HALF_BITS) - 1
        comptime BUCKETS = 1 << HALF_BITS
        var codes = List[UInt64](length=self.seeds, fill=0)
        # Each half's seeds bucketed by its code, contiguous, in order, each with its code below it in
        # one word (see `ENTRY_BITS`).
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
        var left_entries = List[UInt64](length=self.seeds, fill=0)
        var right_entries = List[UInt64](length=self.seeds, fill=0)
        var left_fill = left_start.copy()
        var right_fill = right_start.copy()
        for seed in range(self.seeds):
            var code = codes[seed]
            var entry = (UInt64(seed) << UInt64(ENTRY_BITS)) | code
            left_entries[Int(left_fill[Int(code >> UInt64(HALF_BITS))])] = entry
            left_fill[Int(code >> UInt64(HALF_BITS))] += 1
            right_entries[Int(right_fill[Int(code) & HALF_MASK])] = entry
            right_fill[Int(code) & HALF_MASK] += 1

        # Every window of `INEXACT_LENGTH` bases of the second sequence as one code, the first base in
        # the high bits; past the end, the bases read as zero, and no window reaching there is tried.
        var windows = List[UInt64](length=self.rows + 1, fill=0)
        var rolling = UInt64(0)
        for row in range(self.rows + INEXACT_LENGTH - 1, -1, -1):
            var base = UInt64(second[unsafe_offset=row]) if row < self.rows else UInt64(0)
            rolling = (rolling >> 2) | (base << UInt64(2 * INEXACT_LENGTH - 2))
            if row <= self.rows:
                windows[row] = rolling
        self.inexact_scan(
            left_start.unsafe_ptr(),
            right_start.unsafe_ptr(),
            left_entries.unsafe_ptr(),
            right_entries.unsafe_ptr(),
            windows.unsafe_ptr(),
            found_seed,
            found_row,
            found_end,
            found_cost,
        )

    def inexact_scan(
        self,
        left_start: ImmPointer[Int32, _],
        right_start: ImmPointer[Int32, _],
        left_entries: ImmPointer[UInt64, _],
        right_entries: ImmPointer[UInt64, _],
        window: ImmPointer[UInt64, _],
        mut found_seed: List[Int32],
        mut found_row: List[Int32],
        mut found_end: List[Int32],
        mut found_cost: List[Int32],
    ):
        """The inexact matches every row of the second sequence finds, from the half tables (each
        half's bucket bounds and entries) and the windows.

        Every seed a half finds is tested against the three windows it could match without a branch
        (see `within_one_edit`), most failing all three, and only one that passes is tried in full.

        The tables are larger than a core's own cache on some machines, and every lookup goes from
        a bucket's bounds to its entries, each load a miss waiting on the last. So the rows go in
        batches of `SCAN_BATCH`: first every row's bounds, the loads independent and overlapping,
        their entries fetched ahead, and then the tests, the entries on their way or arrived.
        """
        comptime HALF = INEXACT_LENGTH // 2
        comptime K = INEXACT_LENGTH
        comptime CODE = (UInt64(1) << UInt64(ENTRY_BITS)) - 1
        comptime AHEAD = PrefetchOptions().for_read().high_locality()
        var batch = List[Int32](length=4 * SCAN_BATCH, fill=0)
        var bounds = batch.unsafe_ptr()
        # A left half needs half a seed after it.
        var end_row = self.rows - HALF + 1
        for batch_row in range(0, end_row, SCAN_BATCH):
            var batch_end = min(batch_row + SCAN_BATCH, end_row)
            for row in range(batch_row, batch_end):
                var half = Int(window[unsafe_offset=row] >> UInt64(HALF_BITS))
                var at = 4 * (row - batch_row)
                bounds[unsafe_offset=at] = left_start[unsafe_offset=half]
                bounds[unsafe_offset=at + 1] = left_start[unsafe_offset=half + 1]
                bounds[unsafe_offset=at + 2] = right_start[unsafe_offset=half]
                bounds[unsafe_offset=at + 3] = right_start[unsafe_offset=half + 1]
                prefetch[AHEAD](left_entries.unsafe_offset(Int(bounds[unsafe_offset=at])))
                prefetch[AHEAD](right_entries.unsafe_offset(Int(bounds[unsafe_offset=at + 2])))
            for row in range(batch_row, batch_end):
                var at = 4 * (row - batch_row)
                # A left half found here: the windows start at `row`, one base shorter than the seed,
                # as long, and one longer, whose last `K` bases start a row on.
                var head = window[unsafe_offset=row]
                var tail = window[unsafe_offset=row + 1]
                var starting = self.reachable_seeds(row + K - 1, row + K + 1)
                var after = Int(bounds[unsafe_offset=at + 1])
                for slot in range(
                    first_reachable(left_entries, Int(bounds[unsafe_offset=at]), after, starting[0]), after
                ):
                    var seed = Int(left_entries[unsafe_offset=slot] >> UInt64(ENTRY_BITS))
                    if seed > starting[1]:
                        break
                    var code = left_entries[unsafe_offset=slot] & CODE
                    if within_one_edit(head >> 2, head, head, tail, code):
                        self.try_windows(code, seed, row, -1, window, found_seed, found_row, found_end, found_cost)
                # A right half ending at `end`: the windows end there, so start a base later, at the
                # same row, or a base earlier than the seed's length back.
                var end = row + HALF
                var shorter = window[unsafe_offset=max(end - K + 1, 0)]
                var level = window[unsafe_offset=max(end - K, 0)]
                var longer = window[unsafe_offset=max(end - K - 1, 0)]
                var ending = self.reachable_seeds(end, end)
                after = Int(bounds[unsafe_offset=at + 3])
                for slot in range(
                    first_reachable(right_entries, Int(bounds[unsafe_offset=at + 2]), after, ending[0]), after
                ):
                    var seed = Int(right_entries[unsafe_offset=slot] >> UInt64(ENTRY_BITS))
                    if seed > ending[1]:
                        break
                    var code = right_entries[unsafe_offset=slot] & CODE
                    if ends_within_one_edit(shorter, level, longer, code):
                        self.try_windows(code, seed, -1, end, window, found_seed, found_row, found_end, found_cost)

    @inline(.always)
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
                var bases = Int(pop_count(differing(window[unsafe_offset=low], code)))
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

    @inline(.always)
    def agreeing_prefix(self, first: UInt64, second: UInt64, bases: Int) -> Int:
        """Leading bases two codes of `bases` bases share."""
        var bits = differing(first, second)
        if bits == 0:
            return bases
        return (Int(count_leading_zeros(bits)) - (64 - 2 * bases)) // 2

    @inline(.always)
    def agreeing_suffix(self, first: UInt64, second: UInt64, bases: Int) -> Int:
        """Trailing bases two codes of `bases` bases share."""
        var bits = differing(first, second)
        if bits == 0:
            return bases
        return Int(count_trailing_zeros(bits)) // 2

    @inline(.always)
    def reachable_seeds(self, lowest_end: Int, highest_end: Int) -> Tuple[Int, Int]:
        """The first and last inexact seed a match ending in rows `[lowest_end, highest_end]` can be
        of and still pass `add_match`.

        Seed `s` at column `16s` charges the seeds after it `2(S - s - 1)`, so its match ending at row
        `e` passes when `16s + 16 - e - 2(S - s - 1) <= C - R` and `e - 16s - 16 - 2(S - s - 1) <= R - C`,
        `S` the seeds, `C` and `R` the columns and rows: `18s <= C - R + e + 2S - 18` and
        `14s >= e + C - R - 2S - 14`, a range of seeds, here taken at its widest over the ends.
        """
        comptime assert INEXACT_LENGTH == 16, "the bounds below are worked out for 16-base seeds"
        var shift = self.columns - self.rows
        var first = ceildiv(lowest_end + shift - 2 * self.seeds - 14, 14)
        var last = (highest_end + shift + 2 * self.seeds - 18) // 18
        return (max(first, 0), min(last, self.seeds - 1))

    @inline(.always)
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

    # Out of line: inlined into `build`, it slows the pruning loop there by a few percent.
    @inline(.never)
    def worth_keeping(
        self,
        first: ImmPointer[UInt8, _],
        second: ImmPointer[UInt8, _],
        seed: Int,
        start_row: Int,
        end_row: Int,
        match_cost: Int,
        leftmost: List[Int32],
        mut fronts: List[Int],
    ) -> Bool:
        """A*PA's local pruning: whether a match can lower the cost of some path, so dropping it would not.

        A match only helps a path that goes on to cross the next seeds for fewer edits than the
        seeds it crosses would cost unmatched. A diagonal transition from the match's end, starting
        at the match's own cost, over the next `LOOKAHEAD_SEEDS` seeds counting its own, looks for such a
        path: the match stays if a front reaches past the last of them, or slides into a match
        already kept on its diagonal, which then continues it.
        A front whose edits already equal the seeds crossed can no longer gain, and is dropped. Every
        match an optimal path's chain relies on passes, so the heuristic stays a lower bound.
        """
        var start_column = seed * self.length
        var start_potential = self.potential(start_column)
        var last = min(seed + LOOKAHEAD_SEEDS - 1, self.seeds - 1)
        var end_column = (last + 1) * self.length
        var reach = start_potential - self.potential(end_column)
        var origin = (start_column + self.length) - end_row
        # Front `d` is the diagonal `origin + d - reach`, at `fronts[d]`, its furthest column.
        var low = reach
        var high = reach + 1
        var front = fronts.unsafe_ptr()
        var reached = extend(first, second, start_column + self.length, end_row, self.columns, self.rows)
        front[unsafe_offset=reach] = reached
        if reached >= end_column:
            return True
        var nearest = leftmost.unsafe_ptr()
        var kept = Int(nearest[unsafe_offset=origin + self.rows])
        if kept <= reached:
            return True
        # Past the live fronts, one on either side reads as unreachable, so every front takes the best
        # of its three sources with no test of which exist.
        comptime UNREACHED = -(1 << 40)
        for cost in range(match_cost + 1, reach):
            # One more edit: from the same diagonal, or from either neighbour, in place, the front
            # below's value before this edit carried along.
            front[unsafe_offset=low - 1] = UNREACHED
            front[unsafe_offset=high] = UNREACHED
            front[unsafe_offset=high + 1] = UNREACHED
            var below = UNREACHED
            for d in range(low - 1, high + 1):
                var current = front[unsafe_offset=d]
                front[unsafe_offset=d] = min(max(max(current, below) + 1, front[unsafe_offset=d + 1]), self.columns)
                below = current
            low -= 1
            high += 1
            # A front whose edits match the seeds it has crossed can no longer gain.
            while low < high and cost + self.potential(front[unsafe_offset=low]) >= start_potential:
                low += 1
            while high > low and cost + self.potential(front[unsafe_offset=high - 1]) >= start_potential:
                high -= 1
            if low == high:
                return False
            for d in range(low, high):
                var diagonal = origin + d - reach
                var before = front[unsafe_offset=d]
                var row = before - diagonal
                if row < 0 or row > self.rows:
                    continue
                var after = extend(first, second, before, row, self.columns, self.rows)
                front[unsafe_offset=d] = after
                if after >= end_column:
                    return True
                # A point inside the matrix, so its diagonal indexes `leftmost`.
                var kept_column = Int(nearest[unsafe_offset=diagonal + self.rows])
                if before <= kept_column and kept_column <= after:
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

    @inline(.always)
    def contains(self, layer: Int, x: Int, y: Int) -> Bool:
        """Whether a start in `layer` lies at or above and right of `(x, y)`."""
        var count = Int(self.counts.unsafe_ptr()[unsafe_offset=layer])
        var xs = self.slot_x.unsafe_ptr().unsafe_offset(layer * LAYER_SLOTS)
        var ys = self.slot_y.unsafe_ptr().unsafe_offset(layer * LAYER_SLOTS)
        for index in range(count):
            if x <= Int(xs[unsafe_offset=index]) and y <= Int(ys[unsafe_offset=index]):
                return True
        if count == LAYER_SLOTS:
            var spill_x = self.spill_x.unsafe_ptr()
            var spill_y = self.spill_y.unsafe_ptr()
            var spill_next = self.spill_next.unsafe_ptr()
            var index = Int(self.spill_head.unsafe_ptr()[unsafe_offset=layer])
            while index >= 0:
                if x <= Int(spill_x[unsafe_offset=index]) and y <= Int(spill_y[unsafe_offset=index]):
                    return True
                index = Int(spill_next[unsafe_offset=index])
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

    @inline(.always)
    def potential(self, column: Int) -> Int:
        """What the seeds starting at or after `column` cost a path matching none of them."""
        # Each length its own constant divisor: the band asks for this at every row it prunes.
        if self.cost == 1:
            return self.seeds - min(self.seeds, ceildiv(column, SEED_LENGTH))
        return 2 * (self.seeds - min(self.seeds, ceildiv(column, INEXACT_LENGTH)))

    @inline(.always)
    def climb(self) -> Int:
        """The most the heuristic changes a row down or up a column: one with no seeds or exact ones,
        `cost` with inexact. A best chain from a row still starts, but for its first match, from the
        row above or below, as every match moves the transformed point at least one up and one
        right; the first match scores at most `cost`.

        Down a column a complete set of matches would hold the drop to one, the heuristic being
        consistent, but an exact match's one-edit neighbours are left out (see `try_windows`), and
        without them only this bound holds."""
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


struct Band(Movable):
    """One round of band doubling, advanced a tile at a time.

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
    var edge: Edge
    """The scores down the last tile's right edge."""
    var checkpoint: Int
    """Checkpoints passed."""

    def __init__(
        out self,
        mut profile: Profile,
        threshold: Int,
        adapt: Bool,
    ):
        self.columns = profile.columns
        self.rows = profile.rows
        self.words = profile.words
        self.difference = self.rows - self.columns
        self.extra = (threshold - abs(self.difference)) // 2
        self.threshold = threshold
        # One tile's width of horizontal edge, since each tile starts again from `+1` above (see
        # `Sweep.shifted`); the last tile may absorb a sliver of up to `2 * LANES` more columns.
        self.frontier = Frontier(BAND_COLUMNS + 2 * LANES, self.words)
        profile.build_planes()
        self.sweep = self.frontier.sweep(profile)
        self.bounds = tile_bounds(self.columns, NARROW_COLUMNS if threshold < NARROW_BAND else BAND_COLUMNS)
        self.top = 0
        self.end_word = 0
        self.anchor = 0
        self.deepest = 0
        self.floor = 0
        self.outcome = Round(-1, -1, threshold, -1)
        self.adapt = adapt
        self.edge = Edge()
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
                reach_row = max(reach_row - ceildiv(over, step), diagonal_row)
        var reach = ceildiv(reach_row, WORD_BITS)
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

    def finish(mut self, tile: Int, mut heuristic: SeedHeuristic) -> Bool:
        """Prunes after a swept tile; false, with `outcome` set, when every row went."""
        var end_column = self.bounds[tile + 1]
        self.anchor += end_column - self.bounds[tile]
        # Read the right edge back as scores and keep, as A*PA2 does, the rows from the first to the
        # last whose score plus heuristic fits the bound. Scores change by at most one per row and the
        # heuristic by at most its `climb` either way, so a row `x` over the bound rules out the next
        # `ceil(x / (1 + climb))` rows on either side without reading them: half of `x` for exact
        # seeds or none, a third for inexact.
        self.edge.capture(self.top, self.end_word, self.anchor, self.frontier, self.rows)
        var first_kept = self.edge.low_row
        var last_kept = self.edge.high_row
        var step = 1 + heuristic.climb()
        while first_kept <= last_kept:
            var over = self.edge.score(first_kept) + heuristic.h(end_column, first_kept) - self.threshold
            if over <= 0:
                break
            first_kept += ceildiv(over, step)
        while last_kept >= first_kept:
            var over = self.edge.score(last_kept) + heuristic.h(end_column, last_kept) - self.threshold
            if over <= 0:
                break
            last_kept -= ceildiv(over, step)
        if first_kept > last_kept:
            self.outcome = Round(-1, end_column, self.threshold, -1)
            return False
        if self.adapt and not self.check(end_column, first_kept, last_kept, heuristic):
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
        self.floor = self.edge.score(last_kept)
        return True

    def check(mut self, end_column: Int, first_kept: Int, last_kept: Int, mut heuristic: SeedHeuristic) -> Bool:
        """Re-aims the bound from the band's own climb at a checkpoint; false, with `outcome` set, to
        give the round up.

        At an eighth, a quarter and half of the columns, the least score plus heuristic down the kept
        rows, the gap to the end or the seeds still ahead, has climbed from the heuristic at the origin
        about in proportion to the columns crossed, so scaled to the whole width it projects the
        distance from hundreds of edits rather than the diagonal transition's handful: within about a
        tenth at an eighth, closer further on. The bound drops to the projection plus `CHECK_MARGIN`
        tenths for each checkpoint still ahead, and a round whose projection less half as much passes
        its bound gives up, its rest likely wasted. Lowering the bound keeps the round exact: on an
        optimal path a cell's score plus its gap to the end is at most the distance, so a distance
        within the final bound passed every column unpruned.
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
            least = min(least, self.edge.score(row) + heuristic.h(end_column, row))
            if row == last_kept:
                break
            row = min(row + CHECK_ROWS, last_kept)
        var gap = abs(self.difference)
        var origin = heuristic.h(0, 0)
        var estimate = origin + max(least - origin, 0) * self.columns // end_column
        var margin = estimate * (CHECKPOINTS + 1 - passed) * CHECK_MARGIN // 10
        if estimate - margin // 2 > self.threshold:
            self.outcome = Round(-1, end_column, self.threshold, estimate)
            return False
        var aim = estimate + margin + PROBE_MARGIN
        if aim < self.threshold:
            self.threshold = aim
            self.extra = (aim - gap) // 2
        return True

    def result(mut self) -> Round:
        """What the round found once every tile is swept, with the last edge captured."""
        self.edge.capture(self.top, self.end_word, self.anchor, self.frontier, self.rows)
        if self.end_word < self.words:
            return Round(-1, self.columns, self.threshold, -1)
        return Round(self.edge.score(self.rows), self.columns, self.threshold, -1)


def pruned_distance[
    record: Bool
](mut profile: Profile, threshold: Int, mut trail: Trail, mut heuristic: SeedHeuristic, adapt: Bool,) -> Round:
    """One round of band doubling with A*PA2-simple's pruning.

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
    the whole matrix: along the top row, then down the right edge. Each tile sweeps exactly the
    words the band reaches.

    With `record`, every tile's left edge goes into `trail` before the tile is swept, for the
    traceback; without, the trail is untouched. With `extended`, matches take the third plane (see
    `Profile.extended`).

    With `adapt`, the round re-aims its bound as it goes (see `Band.check`), and a distance is exact
    only within the bound it ends on.
    """
    comptime if record:
        trail.clear()
    var band = Band(profile, threshold, adapt)
    for tile in range(band.tiles()):
        if not band.prepare[record](tile, trail, heuristic):
            return band.outcome
        band.tile_sweep(tile).words(
            profile.extended, band.top, band.end_word, band.first_column(tile), band.end_column(tile)
        )
        if not band.finish(tile, heuristic):
            return band.outcome
    return band.result()


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

comptime PROBE_MARGIN = 16
"""How far past the projected distance the band's first bound reaches, on top of an eighth of it."""

comptime SHORT_COLUMNS = 4096
"""Pairs up to this many columns aim the band's first bound further past the projection.

On a short pair the projection from a handful of edits ran as much as 1.7 times under the distance,
and a first round that fails there costs nearly a whole round more, where aiming wide costs a band
some rows taller; a long pair's projection landed within about a sixth either way."""

comptime SHORT_BAND_COLUMNS = 2048
"""Pairs up to this many columns whose band costs more than its columns: aimed wide (see `SHORT_AIM`),
its first bound lifts it over the whole matrix, which it sweeps and then traces back through
full-height tiles. On ONT reads of about a kilobase that came to some 40 µs a read, several times what
the columns count, so a diagonal transition the per-column budget gave up on still beat it."""

comptime SHORT_BAND_SETUP = 30_000
"""What a band costs a pair of up to `SHORT_BAND_COLUMNS` columns beyond its columns, in diagonal
steps: ONT reads of about a kilobase aligned fastest from 20,000 to 40,000, a quarter faster than with
none, and longer pairs, whose bands stay narrow, are left out."""

comptime SHORT_AIM = 17
"""A short pair's first bound, in tenths of the projection."""

comptime SHORT_REACH = 128
"""The most a short pair's first bound reaches past the projection, so a projection already too
high, as on very divergent pairs, does not widen the band by most of the matrix."""

comptime CHECKPOINTS = 3
"""A first round without seeds re-aims its bound at `columns >> k` for `k = CHECKPOINTS ..= 1`: an
eighth, a quarter and half of the way across (see `Band.check`)."""

comptime CHECK_MARGIN = 1
"""Tenths of its projection a checkpoint's bound allows, for each checkpoint still to come: its
projection strayed by up to about an eighth at the first, a twentieth at the last."""

comptime CHECK_ROWS = 8
"""Rows between the scores a checkpoint samples down the band."""

comptime PROBE_CEILING = 2048
"""The highest score the diagonal transition reaches, which bounds its memory to a few megabytes."""

comptime TWO_ENDED_LIMIT = PROBE_CEILING // 2 + 1
"""Scores each direction of the two-ended search keeps room for, the two together the ceiling."""

comptime FRONT_CENTER = TWO_ENDED_LIMIT + FRONT_PADDING
"""Where diagonal zero sits in a two-ended front's slot: room on the left for its lowest diagonal and
padding, on the right also for a step's last vector, which reads `FRONT_LANES` past the front."""

comptime FRONT_WIDTH = 2 * FRONT_CENTER + FRONT_LANES + 1
"""A two-ended front's slot."""


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

    @inline(.always)
    def at(self, score: Int, diagonal: Int) -> Int:
        """The furthest column of `diagonal` at `score`, or -1 when no path of that score reaches it."""
        var low = Int(self.lows.unsafe_ptr()[unsafe_offset=score])
        if diagonal < low or diagonal > Int(self.highs.unsafe_ptr()[unsafe_offset=score]):
            return -1
        var start = Int(self.starts.unsafe_ptr()[unsafe_offset=score])
        var column = Int(self.offsets.unsafe_ptr()[unsafe_offset=start + FRONT_PADDING + diagonal - low])
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


@inline(.always)
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


@inline(.always)
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


@inline(.always)
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


@inline(.always)
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


@inline(.always)
def step_budget(columns: Int, step_tenths: Int, estimate: Int) -> Int:
    """How many diagonal steps cost what a band over `columns` columns would, at a projected distance.

    A band's column costs a fixed part plus a part growing with the distance, its band taller; in
    steps that is `step_tenths / 10 + estimate / EDITS_PER_STEP` a column. A short pair's band costs
    `SHORT_BAND_SETUP` more (see `SHORT_BAND_COLUMNS`).
    """
    var setup = SHORT_BAND_SETUP if columns <= SHORT_BAND_COLUMNS else 0
    return setup + columns * step_tenths // 10 + columns * estimate // EDITS_PER_STEP


@inline(.always)
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


@inline(.never)
def step_front[
    measure: Bool
](
    previous: ImmPointer[Int32, _],
    current: MutPointer[Int32, _],
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
    while score < PROBE_CEILING:
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
    return Probe(-1, PROBE_CEILING + 1, PROBE_CEILING)


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


comptime FRONT_RING = 2
"""Fronts a direction keeps in turn: the latest, and the one before, which an overlap check reads; the
next grows over the one before."""


struct FrontPair(Movable):
    """One direction's latest fronts, in a ring of `FRONT_RING` slots of `FRONT_WIDTH`, each indexed
    by diagonal, for up to `TWO_ENDED_LIMIT` scores."""

    var buffers: List[Int32]
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

    def __init__(out self, record: Bool = False):
        # Every front writes its own diagonals and the padding either side before the next reads it,
        # so only the first front's surroundings need setting.
        self.buffers = List[Int32](capacity=FRONT_RING * FRONT_WIDTH)
        self.buffers.resize(unsafe_uninit_length=FRONT_RING * FRONT_WIDTH)
        self.low = 0
        self.high = 0
        self.previous_low = 0
        self.previous_high = -1
        self.score = 0
        self.slot = 0
        self.furthest = 0
        self.history = DiagonalFronts(reserve=record)
        self.record = record
        var first = self.front_mut(0)
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

    @inline(.always)
    def previous(self) -> Int:
        """The ring slot of the front before the latest."""
        return self.slot - 1 if self.slot > 0 else FRONT_RING - 1

    @inline(.always)
    def advance(
        mut self, codes: ImmPointer[UInt8, _], others: ImmPointer[UInt8, _], measure: Bool, columns: Int, rows: Int
    ):
        """One more score: the next front in the ring, grown from the latest."""
        var previous = self.front(self.slot)
        var next_slot = self.slot + 1 if self.slot + 1 < FRONT_RING else 0
        var current = self.front_mut(next_slot)
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
        self.slot = next_slot
        self.keep()

    @inline(.always)
    def front(self, which: Int) -> ImmPointer[Int32, ImmUntrackedOrigin]:
        """Front `which` of the ring, to read, indexed by diagonal."""
        return (
            self.buffers.unsafe_ptr()
            .unsafe_origin_cast[ImmUntrackedOrigin]()
            .unsafe_offset(which * FRONT_WIDTH + FRONT_CENTER)
        )

    @inline(.always)
    def front_mut(mut self, which: Int) -> MutPointer[Int32, MutUntrackedOrigin]:
        """Front `which` of the ring, to write, indexed by diagonal. Untracked, as a front is written
        while `front` reads the one before it from the same buffer."""
        return (
            self.buffers.unsafe_ptr()
            .unsafe_origin_cast[MutUntrackedOrigin]()
            .unsafe_offset(which * FRONT_WIDTH + FRONT_CENTER)
        )


@inline(.always)
def overlap(
    forward: ImmPointer[Int32, _],
    forward_low: Int,
    forward_high: Int,
    backward: ImmPointer[Int32, _],
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


def two_ended_distance(profile: Profile, step_tenths: Int) -> Probe:
    """The edit distance by `two_ended`, keeping no history."""
    var first_back = reversed_codes(profile.column_codes, profile.columns, FIRST_SENTINEL)
    var second_back = reversed_codes(profile.row_codes, profile.rows, SECOND_SENTINEL)
    var ahead = FrontPair()
    var behind = FrontPair()
    return two_ended(profile, first_back, second_back, step_tenths, ahead, behind).probe


def two_ended(
    profile: Profile,
    first_back: List[UInt8],
    second_back: List[UInt8],
    step_tenths: Int,
    mut ahead: FrontPair,
    mut behind: FrontPair,
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
    `step_budget`), projecting the distance from both fronts' progress.
    """
    var columns = profile.columns
    var rows = profile.rows
    var target = columns - rows
    var first = profile.column_codes.unsafe_ptr()
    var second = profile.row_codes.unsafe_ptr()
    var first_reversed = first_back.unsafe_ptr()
    var second_reversed = second_back.unsafe_ptr()
    var limit = TWO_ENDED_LIMIT
    ahead.front_mut(0)[unsafe_offset=0] = Int32(slide_forward(first, second, 0, 0))
    behind.front_mut(0)[unsafe_offset=0] = Int32(slide_forward(first_reversed, second_reversed, 0, 0))
    ahead.keep()
    behind.keep()

    @inline(.always)
    def overlapping(
        ahead: FrontPair, behind: FrontPair, forward_earlier: Bool, backward_earlier: Bool
    ) {imm target, imm columns} -> Int:
        """`overlap` of the forward front, or the one before it, and the backward front, or the one
        before it."""
        return overlap(
            ahead.front(ahead.previous() if forward_earlier else ahead.slot),
            ahead.previous_low if forward_earlier else ahead.low,
            ahead.previous_high if forward_earlier else ahead.high,
            behind.front(behind.previous() if backward_earlier else behind.slot),
            behind.previous_low if backward_earlier else behind.low,
            behind.previous_high if backward_earlier else behind.high,
            target,
            columns,
        )

    @inline(.always)
    def met(
        ahead: FrontPair, behind: FrontPair, diagonal: Int, forward_earlier: Bool, backward_earlier: Bool
    ) -> Meeting:
        var column = Int(ahead.front(ahead.previous() if forward_earlier else ahead.slot)[unsafe_offset=diagonal])
        var forward_score = ahead.score - 1 if forward_earlier else ahead.score
        var backward_score = behind.score - 1 if backward_earlier else behind.score
        var total = forward_score + backward_score
        return Meeting(Probe(total, total, total - 1), diagonal, column, forward_score, backward_score)

    var meeting = overlapping(ahead, behind, False, False)
    if meeting != NO_DIAGONAL:
        return met(ahead, behind, meeting, False, False)
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
            meeting = overlapping(ahead, behind, False, False)
            if meeting != NO_DIAGONAL:
                var sooner = overlapping(ahead, behind, False, True)
                if sooner != NO_DIAGONAL:
                    return met(ahead, behind, sooner, False, True)
                return met(ahead, behind, meeting, False, False)
        if checking:
            var estimate = two_ended_gives_up(total, ahead.furthest + behind.furthest, columns, rows, step_tenths)
            if estimate >= 0:
                return Meeting(Probe(-1, estimate, total), 0, 0, 0, 0)
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


def two_ended_gives_up(total: Int, reached: Int, columns: Int, rows: Int, step_tenths: Int) -> Int:
    """The two-ended search's projected distance once what is left of it passes the budget, or -1."""
    var estimate = total * (columns + rows) // max(reached, 1)
    # What is left, about half of `estimate² - total²` diagonals, against what a band costs.
    # Both fronts' steps, about 0.57 ns per square edit with the overlap check, in one-front steps.
    if (estimate * estimate - total * total) * TWO_ENDED_PERCENT // 100 > step_budget(
        columns, step_tenths, estimate
    ) and total * total * TWO_ENDED_PERCENT // 100 * GIVE_UP_SHARE >= step_budget(columns, step_tenths, total):
        return max(estimate, total + 1)
    return -1


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


def band_doubling[
    record: Bool
](
    mut profile: Profile,
    give_up_wide: Bool,
    probe: Probe,
    mut trail: Trail,
    mut heuristic: SeedHeuristic,
    trusted: Bool,
) -> Int:
    """Band doubling: the distance, or -1 once the band would cover the matrix and `give_up_wide`
    hands it back for a sweep of the whole.

    An untrusted projection (see `trusted_projection`) neither aims the first bound nor sets the
    next: the first starts `DOUBLING_START` past the distance already ruled out, by the floor or the
    heuristic at the origin, and each next one at most doubles, as A*PA2's band doubling grows,
    whatever a failed round projected.

    Each round sweeps the band for one bound. A failed round leaves either a real alignment's cost,
    which caps the next bound, or the column at which it pruned every row, from which the distance is
    estimated as the bound scaled to the whole width; the next bound aims just past the estimate
    rather than doubling. Once the band would
    cover the matrix, `give_up_wide` hands back to the caller, else the bound is lifted past any path.

    The first bound aims just past where the diagonal transition's `probe` projected the distance,
    and above every score it ruled out.
    """
    var columns = profile.columns
    var rows = profile.rows
    var gap = abs(rows - columns)
    var aimed = probe.estimate + probe.estimate // 8
    if columns <= SHORT_COLUMNS:
        aimed = min(probe.estimate * SHORT_AIM // 10, probe.estimate + SHORT_REACH)
    # The heuristic at the origin is a lower bound on the distance. With seeds it lands within about a
    # sixth of it on a close pair, where the projection from a few edits strays further.
    var origin = heuristic.h(0, 0)
    var threshold = max(gap, probe.floor + 1, aimed + PROBE_MARGIN)
    if not trusted:
        # As A*PA2 starts: a step past what is known, the heuristic at the origin or the floor.
        threshold = min(threshold, max(gap, probe.floor + 1, origin) + DOUBLING_START)
    if heuristic.seeds > 0:
        if heuristic.chains_well():
            # Matches survive: the origin's bound lies within a few percent of the distance on a
            # close pair, closer than any projection, so start just past it. On a more divergent
            # pair the round dies early, a share of the way across proportional to how far the
            # bound sits above the origin's, and its death estimates the rest (see below).
            threshold = origin + SEED_SLACK
        # Otherwise few matches survive and the bound is one edit a seed, well short of a divergent
        # pair's distance: the projection leads, the bound only floors it.
        threshold = max(threshold, gap, probe.floor + 1, origin + SEED_SLACK)
    var best = Int.MAX
    # Only the first round re-aims at checkpoints and may give itself up there: on reads whose errors
    # gather at an end, a later round given up on its own climb would jump past a bound that was enough.
    var first_round = True
    while True:
        if 2 * (threshold + BAND_COLUMNS) >= rows:
            if give_up_wide:
                return -1
            threshold = max(threshold, columns + rows)
        # Symbols past `ACGT` take their own copy of the round, so the bases' copy matches on two planes.
        var attempt = pruned_distance[record](profile, threshold, trail, heuristic, first_round)
        first_round = False

        var found = attempt.distance
        # A checkpoint may have lowered the bound, and only a distance within the lowered one is exact.
        var bound = attempt.bound
        if found >= 0 and found <= bound:
            return found

        # Grown from the bound the round ended on, which a checkpoint may have lowered: from the one it
        # started on, a round lowered and then failed would jump back to a bound it already knew was loose.
        var next_bound = 2 * bound
        if found >= 0:
            # Still the cost of a real alignment, so it caps every later bound.
            best = min(best, found)
        else:
            # The round pruned every row `reached` columns into its sweep, where the best alignment's
            # score plus heuristic had climbed from the heuristic at the origin past the bound;
            # that climb scaled to the whole width estimates the distance.
            var estimate = Int.MAX
            if attempt.estimate >= 0:
                # A checkpoint gave the round up, with its own projection.
                estimate = attempt.estimate
            elif attempt.reached > 0 and (attempt.reached < columns or heuristic.seeds > 0):
                # A seeded round that crossed every column and still missed the end climbed past the
                # bound in its last tile: the bound itself is the projection, and the retry aims its
                # margin past it rather than doubling a bound that may be nearly enough.
                estimate = min(estimate, origin + (bound - origin) * columns // attempt.reached)
            if estimate != Int.MAX:
                next_bound = max(bound + bound // 4, estimate + estimate // 8 + PROBE_MARGIN)
                if heuristic.seeds > 0:
                    # The origin's bound is certain; only the climb above it is estimated, and the
                    # retry aims a quarter of that climb past it: half overshot the final bound by a
                    # tenth or more on the long reads, widening the band that succeeds.
                    next_bound = max(bound + SEED_SLACK, estimate + (estimate - origin) // 4 + PROBE_MARGIN)
        if found < 0 and heuristic.seeds > 0 and attempt.reached * SEEDED_TRUST_SHARE < columns:
            # A seeded round that died within its first columns projects from those alone, which on
            # real reads hold their errors gathered at the start: its margin over the origin's bound
            # grows `SEEDED_GROWTH` times instead, as A*PA2's grows, until a round gets far enough in
            # for its death to say where the distance lies.
            next_bound = origin + SEEDED_GROWTH * max(bound - origin, SEED_SLACK)
        if not trusted:
            # A round's estimate comes from where it stopped, which errors gathered at an end can
            # set as far off as the first projection: at most doubling instead.
            next_bound = min(next_bound, 2 * bound)
        threshold = min(next_bound, best)


def edit_distance(first: String, second: String) raises AlignmentError -> Int:
    """The global edit distance between two sequences, by bit-parallel sweep.

    Built for DNA over `ACGT`. Up to four other bytes, `N` among them, are symbols of their own,
    each matching only itself; a pair holding them is aligned without the seed heuristic, so a long
    divergent one runs slower than bases alone would.

    Band doubling, as in A*PA2-simple: guess a bound, sweep only the band of cells a path within it
    could cross, and raise the guess until the answer fits under it, which proves it optimal. Close
    sequences therefore cost far less than the whole matrix. Once the band would cover most of the
    matrix, the whole matrix is swept instead.
    """
    var profile = Profile(first, second)
    if profile.columns == 0 or profile.rows == 0:
        return profile.columns + profile.rows
    var search = two_ended_distance(profile, STEP_TENTHS_DISTANCE)
    if search.distance >= 0:
        return search.distance
    var probe = search
    var trusted = True
    var heuristic = band_start(profile, search, probe, trusted)
    var trail = Trail()
    var distance = band_doubling[False](profile, True, probe, trail, heuristic, trusted)
    if distance >= 0:
        return distance
    return full_distance(profile)


def band_start(profile: Profile, search: Probe, mut probe: Probe, mut trusted: Bool) -> SeedHeuristic:
    """What a band starts from once diagonal transition has left the distance to it: the seed
    heuristic its gate chooses, returned, with `probe` set to the projection to aim at and `trusted`
    to whether it agrees with `search`'s (see `trusted_projection`).

    The band's bounds and the seeds' gate are tuned on one front's projection from its first
    `PROBE_START` edits, so that is the projection handed on, unless `search` went far further.
    """
    var fronts = DiagonalFronts()
    var projected = diagonal_transition(profile, PROJECTION_ONLY, fronts)
    trusted = trusted_projection(search, projected)
    probe = Probe(-1, projected.estimate, max(search.floor, projected.floor))
    var seeded = profile.columns >= SEED_COLUMNS or (
        projected.estimate >= SEED_EDITS and projected.estimate * SEED_DIVERGENCE <= profile.columns
    )
    # Seeds are packed two bits a base, which symbols past `ACGT` do not fit; such a pair sweeps its
    # band on the gap heuristic, exact still, only less narrowed.
    if not seeded or profile.extended:
        return SeedHeuristic(profile.columns, profile.rows)
    # Long pairs only: inexact seeds when the projection already says they pay, else exact ones
    # rebuilt inexact if they chain poorly (see `SeedHeuristic`).
    var long = profile.columns >= INEXACT_COLUMNS
    var inexact = long and projected.estimate * INEXACT_DIVERGENCE >= profile.columns
    var cutoff = INEXACT_CHAINED_SPREAD if trusted else INEXACT_CHAINED
    return SeedHeuristic(profile, inexact, choose=long, cutoff=min(cutoff, INEXACT_CHAINED))


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
        self.fill(
            frontier.vertical_plus.unsafe_ptr().unsafe_offset(top),
            frontier.vertical_minus.unsafe_ptr().unsafe_offset(top),
            end - top,
            anchor,
        )

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
        self.fill(
            trail.edge_plus.unsafe_ptr().unsafe_offset(offset),
            trail.edge_minus.unsafe_ptr().unsafe_offset(offset),
            count,
            trail.anchors[tile],
        )

    @inline(.always)
    def fill(mut self, plus: ImmPointer[UInt64, _], minus: ImmPointer[UInt64, _], count: Int, anchor: Int):
        """`count` words of differences, and the score at each word's top from `anchor` on, the buffers
        sized once and written through pointers rather than grown an element at a time."""
        self.plus.resize(unsafe_uninit_length=count)
        self.minus.resize(unsafe_uninit_length=count)
        self.bases.resize(unsafe_uninit_length=count + 1)
        var plus_out = self.plus.unsafe_ptr()
        var minus_out = self.minus.unsafe_ptr()
        var bases = self.bases.unsafe_ptr()
        var running = anchor
        bases[unsafe_offset=0] = running
        for word in range(count):
            var up = plus[unsafe_offset=word]
            var down = minus[unsafe_offset=word]
            plus_out[unsafe_offset=word] = up
            minus_out[unsafe_offset=word] = down
            running += word_value(up, down)
            bases[unsafe_offset=word + 1] = running

    @inline(.always)
    def score(self, row: Int) -> Int:
        """The score at `row`, which must lie within `low_row ..= high_row`."""
        var bases = self.bases.unsafe_ptr()
        if row == self.low_row:
            return bases[unsafe_offset=0]
        var word = (row - 1) // WORD_BITS - self.top
        var bits = row - (self.top + word) * WORD_BITS
        var kept = ALL_ONES if bits == WORD_BITS else (UInt64(1) << UInt64(bits)) - 1
        return bases[unsafe_offset=word] + word_value(
            self.plus.unsafe_ptr()[unsafe_offset=word] & kept, self.minus.unsafe_ptr()[unsafe_offset=word] & kept
        )


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


@inline(.always)
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

    @inline(.always)
    def ends_here(column: Int, offset: Int, cost: Int) {imm edge, imm first_column, imm score, imm home} -> Bool:
        if column != first_column:
            return False
        var row = column + home + offset
        return row >= edge.low_row and row <= edge.high_row and edge.score(row) + cost == score

    @inline(.always)
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

        @inline(.always)
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


struct Recompute(Movable):
    """A tile recompute's buffers, kept across tiles so a recompute neither allocates nor zeroes: every
    slot it reads it has written first."""

    var plus: List[UInt64]
    var minus: List[UInt64]
    var bases: List[Int]

    def __init__(out self):
        self.plus = List[UInt64]()
        self.minus = List[UInt64]()
        self.bases = List[Int]()


def recomputed_segment(
    profile: Profile,
    trail: Trail,
    tile: Int,
    end_column: Int,
    end_row: Int,
    score: Int,
    mut buffers: Recompute,
    mut moves: List[UInt8],
) -> Int:
    """Traces a tile back cell by cell, after sweeping it again and keeping every column's differences.

    The fallback for a tile with more edits than the wavefront search is allowed. As A*PA2 does, it
    first recomputes only a window of rows above the point it traces from, `RECOMPUTE_WORDS` words
    and doubling, rather than the band's whole height (see `window_segment`). Appends the moves
    right to left and returns the left-edge row.
    """
    var top = trail.tops[tile]
    # A tile is recorded only with a word in it, so its end is past its top.
    var end_word = clamp(ceildiv(end_row, WORD_BITS), top + 1, trail.ends[tile])
    var window = RECOMPUTE_WORDS
    while True:
        var first_word = max(top, end_word - window)
        var left: Int
        if profile.extended:
            left = window_segment[True](
                profile, trail, tile, end_column, end_row, score, first_word, end_word, buffers, moves
            )
        else:
            left = window_segment[False](
                profile, trail, tile, end_column, end_row, score, first_word, end_word, buffers, moves
            )
        if left >= 0 or first_word == top:
            return left
        window *= 2


def window_segment[
    extended: Bool
](
    profile: Profile,
    trail: Trail,
    tile: Int,
    end_column: Int,
    end_row: Int,
    score: Int,
    first_word: Int,
    end_word: Int,
    mut buffers: Recompute,
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
    buffers.plus.resize(unsafe_uninit_length=(width + 1) * count)
    buffers.minus.resize(unsafe_uninit_length=(width + 1) * count)
    buffers.bases.resize(unsafe_uninit_length=(width + 1) * (count + 1))
    var plus = buffers.plus.unsafe_ptr()
    var minus = buffers.minus.unsafe_ptr()
    var bases = buffers.bases.unsafe_ptr()
    # The left edge's score at the window's top, carried down from the band's.
    var anchor = trail.anchors[tile]
    for word in range(trail.offsets[tile], offset):
        anchor += word_value(trail.edge_plus[word], trail.edge_minus[word])
    for word in range(count):
        plus[unsafe_offset=word] = trail.edge_plus[offset + word]
        minus[unsafe_offset=word] = trail.edge_minus[offset + word]
    var column_low = profile.column_low.unsafe_ptr().unsafe_offset(COLUMN_PADDING + first_column - 1)
    var column_high = profile.column_high.unsafe_ptr().unsafe_offset(COLUMN_PADDING + first_column - 1)
    var row_low = profile.row_low.unsafe_ptr().unsafe_offset(top)
    var row_high = profile.row_high.unsafe_ptr().unsafe_offset(top)
    # The third plane exists only for symbols past `ACGT` (see `Profile.extended`).
    var column_extra = profile.column_extra.unsafe_ptr().unsafe_offset(COLUMN_PADDING + first_column - 1)
    var row_extra = profile.row_extra.unsafe_ptr().unsafe_offset(top)
    for step in range(1, width + 1):
        var horizontal_plus = UInt64(1)
        var horizontal_minus = UInt64(0)
        var low = column_low[unsafe_offset=step]
        var high = column_high[unsafe_offset=step]
        var extra = column_extra[unsafe_offset=step] if extended else UInt64(0)
        for word in range(count):
            var vertical_plus = plus[unsafe_offset=(step - 1) * count + word]
            var vertical_minus = minus[unsafe_offset=(step - 1) * count + word]
            var matches = (low ^ row_low[unsafe_offset=word]) & (high ^ row_high[unsafe_offset=word])
            comptime if extended:
                matches &= extra ^ row_extra[unsafe_offset=word]
            advance[1](horizontal_plus, horizontal_minus, vertical_plus, vertical_minus, matches)
            plus[unsafe_offset=step * count + word] = vertical_plus
            minus[unsafe_offset=step * count + word] = vertical_minus

    # Scores at each word's top on every column: the window's top scores the anchor plus one per column.
    for step in range(width + 1):
        var running = anchor + step
        bases[unsafe_offset=step * (count + 1)] = running
        for word in range(count):
            running += word_value(plus[unsafe_offset=step * count + word], minus[unsafe_offset=step * count + word])
            bases[unsafe_offset=step * (count + 1) + word + 1] = running

    @inline(.always)
    def score_at(step: Int, row: Int) {imm bases, imm plus, imm minus, imm count, imm top} -> Int:
        var word = (row - 1) // WORD_BITS - top if row > top * WORD_BITS else 0
        if row == top * WORD_BITS:
            return bases[unsafe_offset=step * (count + 1)]
        var bits = row - (top + word) * WORD_BITS
        var kept = ALL_ONES if bits == WORD_BITS else (UInt64(1) << UInt64(bits)) - 1
        return bases[unsafe_offset=step * (count + 1) + word] + word_value(
            plus[unsafe_offset=step * count + word] & kept, minus[unsafe_offset=step * count + word] & kept
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
    var buffers = Recompute()
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
            left = recomputed_segment(profile, trail, tile, column, row, current, buffers, moves)
        # Exact, since it plus the segment's cost is the exact score the segment started from.
        current = edge.score(left)
        column = first_column
        row = left
    # The first column is the border: the rest of the way up is gaps against the second sequence.
    for _ in range(row):
        moves.append(UP)


def edit_alignment(first: String, second: String) raises AlignmentError -> AlignmentResult:
    """The global edit distance between two sequences, and an optimal alignment; symbols past `ACGT`
    as `edit_distance` takes them.

    The distance comes from `edit_distance`'s band doubling, recording each tile's left edge in the
    round that succeeds; the alignment is then traced back tile by tile from those edges (see
    `trace_back`). Where the two-ended diagonal transition settles the distance, its fronts give the
    path, traced to the start and on to the end from where they met. The score is the distance, as
    `levenshtein_alignment` reports it.
    """
    var profile = Profile(first, second)
    var columns = profile.columns
    var rows = profile.rows
    # Moves right to left, from the corner, or from where the fronts met, back to the origin.
    var forward_moves = List[UInt8](capacity=columns + rows)
    # Moves left to right, from where the fronts met to the corner; none without meeting.
    var backward_moves = List[UInt8]()
    var distance: Int
    if columns == 0 or rows == 0:
        distance = columns + rows
        for _ in range(rows):
            forward_moves.append(UP)
        for _ in range(columns):
            forward_moves.append(LEFT)
    else:
        # A near-identical pair: one front, kept whole, settles it before both ends' setup would pay.
        var near = DiagonalFronts()
        var close = diagonal_transition(profile, STEP_TENTHS_ALIGNMENT, near, switch_setup=TWO_ENDED_SETUP)
        if close.distance >= 0:
            trace_diagonals(profile, near, close.distance, forward_moves)
            return gapped_rows(first, second, forward_moves, columns, rows, backward_moves, close.distance)
        # Diagonal transition from both ends, keeping every front, while it is cheaper than a band:
        # where the fronts meet, the path is traced back to the start through the forward fronts
        # and on to the end through the backward ones.
        var first_back = reversed_codes(profile.column_codes, columns, FIRST_SENTINEL)
        var second_back = reversed_codes(profile.row_codes, rows, SECOND_SENTINEL)
        var ahead = FrontPair(record=True)
        var behind = FrontPair(record=True)
        var meeting = two_ended(profile, first_back, second_back, STEP_TENTHS_ALIGNMENT, ahead, behind)
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
        var probe = meeting.probe
        var trusted = True
        var heuristic = band_start(profile, meeting.probe, probe, trusted)
        var trail = Trail(columns)
        distance = band_doubling[True](profile, False, probe, trail, heuristic, trusted)
        trace_back(profile, trail, columns, rows, distance, forward_moves)
    return gapped_rows(first, second, forward_moves, columns, rows, backward_moves, distance)


@inline(.always)
def copy_bytes(destination: MutPointer[UInt8, _], source: ImmPointer[UInt8, _], count: Int):
    """`count` bytes, sixteen at a time while that many are left."""
    var index = 0
    while index + 16 <= count:
        destination.unsafe_offset(index).unsafe_store(source.unsafe_offset(index).unsafe_load[width=16]())
        index += 16
    while index < count:
        destination[unsafe_offset=index] = source[unsafe_offset=index]
        index += 1


@inline(.always)
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


# region Semi-global

comptime HIGH_BIT = UInt64(1) << UInt64(WORD_BITS - 1)
"""A word's last row."""


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
            sweep.words(profile.extended, 0, fast_end, first_column, end_column)
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


# endregion Semi-global
