# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Ported from `pa-bitpacking` in A*PA (https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner,
# commit bf2e14e), by Ragnar Groot Koerkamp and Pesho Ivanov, itself translated from Edlib.
"""
Bit-parallel unit-cost edit distance over DNA, after Myers (1999) as A*PA2 implements it: the sequences'
encoding, the sweep's kernels, and the whole matrix.

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
"""

from std.bit import pop_count
from std.math import ceildiv
from std.sys import inlined_assembly, simd_width_of
from std.sys.info import CompilationTarget

from .errors import AlignmentError, ErrorKind


comptime WORD_BITS = 64
"""Rows one word of the second sequence covers."""


comptime LANES = 8
"""
Words advanced together in one SIMD vector, and so the rows of a block: `LANES * WORD_BITS`.

Eight measured fastest on an M2, where four left the pipeline idle and sixteen ran out of registers.
That matches A*PA2's own choice of two four-lane vectors side by side.
"""


comptime PAIRED_GROUPS = simd_width_of[DType.uint64]() >= 8
"""
Whether a sweep runs its groups two at a time, one under the other in one loop (see `pair_block`).

On AVX-512 the pair took 5 to 7% off the long reads and left the rest within noise; sixteen lanes in one
vector did about as well but cost up to 2% on short divergent pairs. On an M2 both were slower.
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


@inline(.always)
def append_diagonals(mut moves: List[UInt8], count: Int):
    """`count` diagonal moves, sixteen a store rather than an append each: a run of matches is most of
    any traceback's moves."""
    var start = len(moves)
    moves.resize(unsafe_uninit_length=start + count)
    var out = moves.unsafe_ptr().unsafe_offset(start)
    var index = 0
    while index + 16 <= count:
        out.unsafe_offset(index).unsafe_store(SIMD[DType.uint8, 16](DIAGONAL))
        index += 16
    while index < count:
        out[unsafe_offset=index] = DIAGONAL
        index += 1


comptime NARROW_COLUMNS = 64
"""
Columns of one band tile when the band itself is narrow. A tile computes every word its columns reach,
so a band slopes down by its own width across each tile; on a band only a few hundred rows tall that
slope is most of the work, and narrower tiles cut it at the cost of more triangles.
"""


comptime CODE_PADDING = 16
"""Sentinel bytes after each sequence's codes, enough for a sixteen-byte comparison at the last base."""


comptime FIRST_SENTINEL = UInt8(0xFE)
"""Past the first sequence's last base: no code, and unequal to `SECOND_SENTINEL`."""


comptime SECOND_SENTINEL = UInt8(0xFF)
"""Past the second sequence's last base."""


comptime COLUMN_PADDING = LANES
"""Columns of zero bases before and after the profile's planes, for lanes standing outside a tile."""


comptime Words = SIMD[DType.uint64, LANES]
"""One vector of `LANES` words."""


comptime ALL_ONES = ~UInt64(0)
"""A word with every bit set."""


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


comptime BASES_ONLY = 0
"""A tile's match test when neither its columns nor any row holds a symbol past `ACGT`: two planes."""


comptime ROW_SYMBOLS = 1
"""A tile's match test when rows hold symbols past `ACGT` but its columns none: those rows match no
column, and the third row plane, stored negated, is all ones but at them, so an AND with it masks
them out, the full test less a load and an XOR."""


comptime ALL_SYMBOLS = 2
"""A tile's match test when its columns hold symbols past `ACGT`: the third plane in full."""


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
        """Pointers into the profile's planes and the frontier's edges, which must outlive the sweep."""
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
    def run[symbols: Int](self, word: Int, first_column: Int, end_column: Int):
        """One word through `[first_column, end_column)`, its vertical differences held in registers.

        Only the horizontal edge goes through memory, and each column's is a different address, so
        consecutive columns do not wait on a store being read back.
        """
        var vp = self.vertical_plus[unsafe_offset=word]
        var vm = self.vertical_minus[unsafe_offset=word]
        var row_low = self.row_low[unsafe_offset=word]
        var row_high = self.row_high[unsafe_offset=word]
        var row_extra = self.row_extra[unsafe_offset=word] if symbols != BASES_ONLY else UInt64(0)
        for column in range(first_column, end_column):
            var hp = self.horizontal_plus[unsafe_offset=column]
            var hm = self.horizontal_minus[unsafe_offset=column]
            var matches = (self.column_low[unsafe_offset=column] ^ row_low) & (
                self.column_high[unsafe_offset=column] ^ row_high
            )
            comptime if symbols == ALL_SYMBOLS:
                matches &= self.column_extra[unsafe_offset=column] ^ row_extra
            elif symbols == ROW_SYMBOLS:
                matches &= row_extra
            advance[1](hp, hm, vp, vm, matches)
            self.horizontal_plus[unsafe_offset=column] = hp
            self.horizontal_minus[unsafe_offset=column] = hm
        self.vertical_plus[unsafe_offset=word] = vp
        self.vertical_minus[unsafe_offset=word] = vm

    def words(self, symbols: Int, first_word: Int, end_word: Int, first_column: Int, end_column: Int):
        """Words `[first_word, end_word)` through columns `[first_column, end_column)`, matching as
        `symbols` says, `BASES_ONLY`, `ROW_SYMBOLS` or `ALL_SYMBOLS` (see `Profile.symbols`)."""
        if symbols == ALL_SYMBOLS:
            self.words_matching[ALL_SYMBOLS](first_word, end_word, first_column, end_column)
        elif symbols == ROW_SYMBOLS:
            self.words_matching[ROW_SYMBOLS](first_word, end_word, first_column, end_column)
        else:
            self.words_matching[BASES_ONLY](first_word, end_word, first_column, end_column)

    # Out of line, each kernel its own function: inlined side by side, or chosen per band round, the
    # two slowed the bases' one by 1 to 4%.
    @inline(.never)
    def words_matching[symbols: Int](self, first_word: Int, end_word: Int, first_column: Int, end_column: Int):
        """`words`, its match test fixed by `symbols`.

        Full groups of `LANES` take the staggered vector sweep when the span is wide enough for it, two
        at a time where `PAIRED_GROUPS`.
        Of what is left, four words take a narrow vector, five to seven a narrow vector with the rest
        in scalar registers beside it, two or three scalar registers alone, and one word goes on its own.
        """
        var word = first_word
        if end_column - first_column >= 2 * LANES:
            comptime if PAIRED_GROUPS:
                while word + 2 * LANES <= end_word:
                    self.pair_block[symbols](word, first_column, end_column)
                    word += 2 * LANES
            while word + LANES <= end_word:
                self.block[LANES, symbols](word, first_column, end_column)
                word += LANES
            # Five to seven words left: a narrow vector with the rest in scalar registers in its shadow.
            var left = end_word - word
            if left == 7:
                self.hybrid_block[NARROW_LANES, 3, symbols](word, first_column, end_column)
                word += 7
            elif left == 6:
                self.hybrid_block[NARROW_LANES, 2, symbols](word, first_column, end_column)
                word += 6
            elif left == 5:
                self.hybrid_block[NARROW_LANES, 1, symbols](word, first_column, end_column)
                word += 5
            elif left == 4:
                self.block[NARROW_LANES, symbols](word, first_column, end_column)
                word += NARROW_LANES
            # Two or three words left run side by side in scalar registers, which beats a vector
            # this narrow: its chain waits two cycles an operation, a scalar one.
            if end_word - word == 3:
                self.scalar_block[3, symbols](word, first_column, end_column)
                word += 3
            elif end_word - word == 2:
                self.scalar_block[2, symbols](word, first_column, end_column)
                word += 2
        while word < end_word:
            self.run[symbols](word, first_column, end_column)
            word += 1

    def block[lanes: Int, symbols: Int](self, first_word: Int, first_column: Int, end_column: Int):
        """`lanes` words through `[first_column, end_column)` in one vector, staggered (see `VectorGroup`)."""
        var group = VectorGroup[lanes, symbols](self, first_word, first_column, end_column)
        stagger(group)
        group.finish()

    def pair_block[symbols: Int](self, first_word: Int, first_column: Int, end_column: Int):
        """Two groups of `LANES` words, one under the other, in one loop, the lower a little behind."""
        var top = VectorGroup[LANES, symbols](self, first_word, first_column, end_column)
        var bottom = VectorGroup[LANES, symbols](self, first_word + LANES, first_column, end_column)
        stagger_pair(top, bottom, LANES + 2)
        top.finish()
        bottom.finish()

    def scalar_block[lanes: Int, symbols: Int](self, first_word: Int, first_column: Int, end_column: Int):
        """`lanes` words through `[first_column, end_column)` in scalar registers, staggered (see `ScalarGroup`)."""
        var group = ScalarGroup[lanes, symbols](self, first_word, first_column, end_column)
        stagger(group)
        group.finish()

    def hybrid_block[lanes: Int, below: Int, symbols: Int](self, first_word: Int, first_column: Int, end_column: Int):
        """`lanes` words in a vector and the `below` words under them in scalar registers, in one loop.

        A narrow vector's chain leaves the integer units idle, so the scalar words run in its
        shadow, a few columns behind so the differences the vector sends down are already stored.
        """
        var top = VectorGroup[lanes, symbols](self, first_word, first_column, end_column)
        var bottom = ScalarGroup[below, symbols](self, first_word + lanes, first_column, end_column)
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
        """The tile's first column."""
        ...

    def end_column(self) -> Int:
        """The column just past the tile's last."""
        ...

    def step[masked: Bool](mut self, offset: Int):
        """One step. `masked` when some word stands outside the tile: such a word keeps its state,
        and the edge is read and written only for columns inside it."""
        ...

    def finish(self):
        """Writes the words' vertical differences back to the frontier."""
        ...


struct VectorGroup[lanes: Int, symbols: Int](Staggered, TrivialRegisterPassable):
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
        """Words `first_word` on, loaded bottom lane first, for columns `[first_column, end_column)`."""
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
            comptime if Self.symbols != BASES_ONLY:
                self.row_extra[lane] = sweep.row_extra[unsafe_offset=word]
            self.vertical_plus[lane] = sweep.vertical_plus[unsafe_offset=word]
            self.vertical_minus[lane] = sweep.vertical_minus[unsafe_offset=word]
            self.lane_columns[lane] = lane + 1

    @inline(.always)
    def width(self) -> Int:
        """The group's `lanes` words."""
        return Self.lanes

    @inline(.always)
    def first_column(self) -> Int:
        """The tile's first column."""
        return self.first

    @inline(.always)
    def end_column(self) -> Int:
        """The column just past the tile's last."""
        return self.end

    @inline(.always)
    def step[masked: Bool](mut self, offset: Int):
        """One step: the lanes pass their differences down a lane, the top lane takes column
        `offset + lanes` from the edge, and the bottom lane leaves column `offset + 1`'s there."""
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
        comptime if Self.symbols == ALL_SYMBOLS:
            matches &= (
                self.sweep.column_extra.unsafe_offset(offset + 1).unsafe_load[width=Self.lanes]() ^ self.row_extra
            )
        elif Self.symbols == ROW_SYMBOLS:
            matches &= self.row_extra
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
        """Writes each lane's vertical differences back to its word in the frontier."""
        comptime for lane in range(Self.lanes):
            var word = self.first_word + Self.lanes - 1 - lane
            self.sweep.vertical_plus[unsafe_offset=word] = self.vertical_plus[lane]
            self.sweep.vertical_minus[unsafe_offset=word] = self.vertical_minus[lane]


struct ScalarGroup[lanes: Int, symbols: Int](Staggered):
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
        """Words `first_word` on, top word first, for columns `[first_column, end_column)`."""
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
            comptime if Self.symbols != BASES_ONLY:
                self.row_extra[j] = sweep.row_extra[unsafe_offset=first_word + j]
            self.vertical_plus[j] = sweep.vertical_plus[unsafe_offset=first_word + j]
            self.vertical_minus[j] = sweep.vertical_minus[unsafe_offset=first_word + j]

    @inline(.always)
    def width(self) -> Int:
        """The group's `lanes` words."""
        return Self.lanes

    @inline(.always)
    def first_column(self) -> Int:
        """The tile's first column."""
        return self.first

    @inline(.always)
    def end_column(self) -> Int:
        """The column just past the tile's last."""
        return self.end

    @inline(.always)
    def step[masked: Bool](mut self, offset: Int):
        """One step: word `j` from the top works column `offset + lanes - j`, the top word taking its
        difference from the edge and the bottom word leaving its own there."""
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
            comptime if Self.symbols == ALL_SYMBOLS:
                matches &= self.sweep.column_extra[unsafe_offset=column] ^ self.row_extra[j]
            elif Self.symbols == ROW_SYMBOLS:
                matches &= self.row_extra[j]
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
        """Writes each word's vertical differences back to the frontier."""
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


@always_inline
def base_codes[width: Int](bytes: SIMD[DType.uint8, width]) -> SIMD[DType.uint8, width]:
    """Each byte's code, right for `A`, `C`, `G` and `T` alone (see `encoded_bases`)."""
    return (((bytes >> 1) ^ (bytes >> 2)) & 1) | ((bytes >> 1) & 2)


@always_inline
def not_bases[width: Int](bytes: SIMD[DType.uint8, width]) -> SIMD[DType.bool, width]:
    """Which bytes are not `A`, `C`, `G` or `T`."""
    return ~(
        bytes.eq(UInt8(ord("A"))) | bytes.eq(UInt8(ord("C"))) | bytes.eq(UInt8(ord("G"))) | bytes.eq(UInt8(ord("T")))
    )


def encoded_bases(source: ImmPointer[UInt8, _], count: Int, target: MutPointer[UInt8, _]) -> Bool:
    """Writes each of `count` bytes' base code to `target`, sixteen at a time, and whether every byte
    was `A`, `C`, `G` or `T`: a byte of these gives its code from its ASCII bits, bit 2 set for `G` and
    `T`, the code's high bit, and bit 1 differing from bit 2 for `C` and `T`, its low bit. Any other
    byte's code is wrong until recoded (see `symbol_codes`). One pass for both, where a short read's
    profile once read each sequence twice."""
    comptime CHUNK = 16
    var others = SIMD[DType.bool, CHUNK](fill=False)
    var index = 0
    while index + CHUNK <= count:
        var bytes = source.unsafe_offset(index).unsafe_load[width=CHUNK]()
        target.unsafe_offset(index).unsafe_store(base_codes[CHUNK](bytes))
        others |= not_bases[CHUNK](bytes)
        index += CHUNK
    var found = others.reduce_or()
    while index < count:
        var byte = source[unsafe_offset=index]
        target[unsafe_offset=index] = base_codes[1](byte)
        found |= not_bases[1](byte)[0]
        index += 1
    return not found


def symbol_codes(
    first: String, second: String, mut column_codes: List[UInt8], mut row_codes: List[UInt8]
) raises AlignmentError:
    """Recodes the bytes that are not bases, which the bases' codes read as one: each the next free
    code from four, in the order the bytes first appear; past four such symbols, three bits run out.

    Only the sixteen-byte chunks holding such a byte are gone through byte by byte: coded one byte at
    a time, a 10 kbp read with one `N` took five times as long to code as without.
    """
    comptime UNSEEN = Int16(-1)
    var table = List[Int16](length=256, fill=UNSEEN)
    table[ord("A")] = 0
    table[ord("C")] = 1
    table[ord("G")] = 2
    table[ord("T")] = 3
    var next_code = 4
    recode_symbols(first, column_codes, table, next_code)
    recode_symbols(second, row_codes, table, next_code)


def recode_symbols(
    text: String, mut codes: List[UInt8], mut table: List[Int16], mut next_code: Int
) raises AlignmentError:
    """`symbol_codes` for one sequence, `table` and `next_code` carried from the one before."""
    comptime CHUNK = 16
    var bytes = text.unsafe_ptr()
    var out = codes.unsafe_ptr()
    var length = text.byte_length()
    var index = 0
    while index < length:
        var size = min(CHUNK, length - index)
        if size == CHUNK:
            var chunk = bytes.unsafe_offset(index).unsafe_load[width=CHUNK]()
            var bases = (
                chunk.eq(UInt8(ord("A")))
                | chunk.eq(UInt8(ord("C")))
                | chunk.eq(UInt8(ord("G")))
                | chunk.eq(UInt8(ord("T")))
            )
            if bases.reduce_and():
                index += CHUNK
                continue
        for at in range(index, index + size):
            var byte = Int(bytes[unsafe_offset=at])
            if table[byte] >= 0 and table[byte] < 4:
                continue
            if table[byte] < 0:
                if next_code == 8:
                    raise AlignmentError(
                        ErrorKind.UNKNOWN_SYMBOL, "bit-parallel edit distance takes ACGT and at most four other symbols"
                    )
                table[byte] = Int16(next_code)
                next_code += 1
            out[unsafe_offset=at] = UInt8(table[byte])
        index += size


def folded(codes: List[UInt8]) -> List[UInt8]:
    """`codes` with each symbol past `ACGT` read as the base its two low bits name, for the seeds, which
    pack two bits a base; the sentinels past the end stay.

    Folding never parts two equal symbols, so it never raises an edit distance: a match found in the
    folded second sequence costs no more than the real one, and every real one is found. It only adds
    matches, which weaken the bound and never break it. The first sequence's seeds that hold such a
    symbol go uncounted instead (see `SeedHeuristic.remaining`).
    """
    var out = List[UInt8](capacity=len(codes))
    out.resize(unsafe_uninit_length=len(codes))
    var source = codes.unsafe_ptr()
    var target = out.unsafe_ptr()
    comptime CHUNK = 16
    var index = 0
    while index + CHUNK <= len(codes):
        var chunk = source.unsafe_offset(index).unsafe_load[width=CHUNK]()
        target.unsafe_offset(index).unsafe_store(chunk.lt(8).select(chunk & 3, chunk))
        index += CHUNK
    while index < len(codes):
        var code = source[unsafe_offset=index]
        target[unsafe_offset=index] = code & 3 if code < 8 else code
        index += 1
    return out^


def reverse_in_place(codes: MutPointer[UInt8, _], count: Int):
    """`count` codes back to front, sixteen from each end at a time."""
    comptime CHUNK = 16
    var low = 0
    var high = count
    while high - low >= 2 * CHUNK:
        var front = codes.unsafe_offset(low).unsafe_load[width=CHUNK]()
        var back = codes.unsafe_offset(high - CHUNK).unsafe_load[width=CHUNK]()
        codes.unsafe_offset(low).unsafe_store(back.reversed())
        codes.unsafe_offset(high - CHUNK).unsafe_store(front.reversed())
        low += CHUNK
        high -= CHUNK
    high -= 1
    while low < high:
        var swapped = codes[unsafe_offset=low]
        codes[unsafe_offset=low] = codes[unsafe_offset=high]
        codes[unsafe_offset=high] = swapped
        low += 1
        high -= 1


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
    the codes' third bit, `column_extra` and `row_extra`; the seeds, packed two bits a base, read the
    codes folded (see `folded`)."""
    var row_symbols: Bool
    """Whether the second sequence holds a symbol past `ACGT`, so that no tile matches on two planes alone."""
    var column_symbols: List[Int32]
    """Per column, the columns before it holding a symbol past `ACGT`, one past the last, built with the
    third plane: a tile whose columns hold none needs no more than `ROW_SYMBOLS` (see `symbols`)."""
    var column_codes: List[UInt8]
    """The first sequence as codes, which the traceback compares base by base."""
    var row_codes: List[UInt8]
    """The second sequence as codes."""

    def __init__(out self, first: String, second: String, reverse: Bool = False) raises AlignmentError:
        """Both sequences as codes: `A`, `C`, `G` and `T` zero to three, and up to four other bytes the
        codes four to seven, in the order they first appear, each matching only itself; with
        `reverse`, both back to front, the profile of the reversed pair."""
        self.columns = 0
        self.rows = 0
        self.words = 0
        self.column_low = List[UInt64]()
        self.column_high = List[UInt64]()
        self.column_extra = List[UInt64]()
        self.row_low = List[UInt64]()
        self.row_high = List[UInt64]()
        self.row_extra = List[UInt64]()
        self.row_symbols = False
        self.column_symbols = List[Int32]()
        self.extended = False
        self.column_codes = List[UInt8]()
        self.row_codes = List[UInt8]()
        self.reset(first, second, reverse)

    def reset(mut self, first: String, second: String, reverse: Bool = False) raises AlignmentError:
        """As new for this pair, the memory the last pair's codes took kept: a batch's worker profiles
        each of its pairs in the same lists (see `edit_distance.EditSpace`)."""
        self.columns = first.byte_length()
        self.rows = second.byte_length()
        self.words = ceildiv(self.rows, WORD_BITS)
        # The planes wait for a band, which builds them anew (see `build_planes`).
        self.column_low.clear()
        self.column_high.clear()
        self.column_extra.clear()
        self.row_low.clear()
        self.row_high.clear()
        self.row_extra.clear()
        self.row_symbols = False
        self.column_symbols.clear()
        # Only the codes are built here; the planes wait for a band (see `build_planes`).
        self.column_codes.resize(unsafe_uninit_length=self.columns)
        self.row_codes.resize(unsafe_uninit_length=self.rows)
        var bases = encoded_bases(first.unsafe_ptr(), self.columns, self.column_codes.unsafe_ptr())
        bases = encoded_bases(second.unsafe_ptr(), self.rows, self.row_codes.unsafe_ptr()) and bases
        self.extended = not bases
        if self.extended:
            symbol_codes(first, second, self.column_codes, self.row_codes)
        if reverse:
            reverse_in_place(self.column_codes.unsafe_ptr(), self.columns)
            reverse_in_place(self.row_codes.unsafe_ptr(), self.rows)
        # Past the last base of each, sentinels that match nothing, the two of them distinct, so a
        # match extension stops at the matrix's edge without checking it (see `slide_forward`).
        self.column_codes.resize(unsafe_uninit_length=self.columns + CODE_PADDING)
        self.row_codes.resize(unsafe_uninit_length=self.rows + CODE_PADDING)
        self.column_codes.unsafe_ptr().unsafe_offset(self.columns).unsafe_store(
            SIMD[DType.uint8, CODE_PADDING](FIRST_SENTINEL)
        )
        self.row_codes.unsafe_ptr().unsafe_offset(self.rows).unsafe_store(
            SIMD[DType.uint8, CODE_PADDING](SECOND_SENTINEL)
        )

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
        self.column_symbols = List[Int32](length=self.columns + 1, fill=0)
        var extra = self.column_extra.unsafe_ptr().unsafe_offset(COLUMN_PADDING)
        var before = self.column_symbols.unsafe_ptr()
        column = 0
        while column + CHUNK <= self.columns:
            var symbols = (column_codes.unsafe_offset(column).unsafe_load[width=CHUNK]() >> 2) & 1
            extra.unsafe_offset(column).unsafe_store(UInt64(0) - symbols.cast[DType.uint64]())
            var count = before[unsafe_offset=column]
            if symbols.reduce_or() == 0:
                for at in range(column + 1, column + CHUNK + 1):
                    before[unsafe_offset=at] = count
            else:
                for at in range(CHUNK):
                    count += Int32(symbols[at])
                    before[unsafe_offset=column + at + 1] = count
            column += CHUNK
        while column < self.columns:
            var symbol = (column_codes[unsafe_offset=column] >> 2) & 1
            extra[unsafe_offset=column] = UInt64(0) - symbol.cast[DType.uint64]()
            before[unsafe_offset=column + 1] = before[unsafe_offset=column] + Int32(symbol)
            column += 1
        var row_extra = self.row_extra.unsafe_ptr()
        var seen = UInt64(0)
        row = 0
        while row + 8 <= self.rows:
            var symbols = (row_codes.unsafe_offset(row).unsafe_bitcast[UInt64]().unsafe_load() >> 2) & ONES
            seen |= symbols
            row_extra[unsafe_offset=row // WORD_BITS] |= (((symbols * GATHER) >> 56) ^ 0xFF) << UInt64(row % WORD_BITS)
            row += 8
        while row < self.rows:
            var symbol = ((row_codes[unsafe_offset=row] >> 2) & 1).cast[DType.uint64]()
            seen |= symbol
            row_extra[unsafe_offset=row // WORD_BITS] |= (symbol ^ 1) << UInt64(row % WORD_BITS)
            row += 1
        self.row_symbols = seen != 0

    def symbols(self, first_column: Int, end_column: Int) -> Int:
        """The match test a sweep of columns `[first_column, end_column)` takes, once the planes are
        built: the third plane in full only where those columns hold a symbol past `ACGT`, a mask of
        the rows holding one elsewhere, and two planes alone when neither does.

        A pair with one `N` used to sweep every tile on the third plane.
        """
        if not self.extended or first_column == end_column:
            return BASES_ONLY
        if self.column_symbols[end_column] != self.column_symbols[first_column]:
            return ALL_SYMBOLS
        return ROW_SYMBOLS if self.row_symbols else BASES_ONLY


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
        """Both edges at the global borders: `columns` horizontal differences and `words` vertical words."""
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
    sweep.words(profile.symbols(0, profile.columns), 0, profile.words, 0, profile.columns)

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
        """Empties the trail for a new round, keeping its buffers."""
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
        self.edge_plus.extend(Span(frontier.vertical_plus)[top:end])
        self.edge_minus.extend(Span(frontier.vertical_minus)[top:end])


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
        """The right edge of the last tile swept: words `top` to `end` of the frontier, scoring `anchor` at
        the top."""
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
        # A row on a word's boundary is read as the bottom of the word above, all 64 of its rows
        # counted, so the last row needs no base past the last word.
        var word = (row - 1) // WORD_BITS - self.top
        var bits = row - (self.top + word) * WORD_BITS
        var kept = ALL_ONES if bits == WORD_BITS else (UInt64(1) << UInt64(bits)) - 1
        return bases[unsafe_offset=word] + word_value(
            self.plus.unsafe_ptr()[unsafe_offset=word] & kept, self.minus.unsafe_ptr()[unsafe_offset=word] & kept
        )
