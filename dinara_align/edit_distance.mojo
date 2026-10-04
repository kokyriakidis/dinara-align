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
guess is aimed from where the failed round died rather than doubled (see `edit_distance`).

A pair too divergent for a band is swept whole, cut into tiles of `LANES` words by `TILE_COLUMNS`
columns. A tile needs only the tile above it and the tile to its left, so every tile on one
anti-diagonal of tiles runs on its own thread; they touch disjoint columns of the horizontal edge
and disjoint words of the vertical one. A*PA sweeps on one thread.

Unit costs only: a substitution, an insertion and a deletion each cost one, so this is the
distance `levenshtein_alignment` returns, without the alignment.
"""

from std.bit import byte_swap, count_leading_zeros, pop_count
from std.math import ceildiv
from std.sys import inlined_assembly

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

comptime BAND_START = 64
"""
The least the first guess reaches past the length difference. A round costs at least a tile's height
of rows across the matrix whatever the guess, so starting from the bare difference only adds rounds.
"""

comptime FIRST_GUESS_DIVISOR = 64
"""
The first guess allows one edit per this many columns, about 1.6%. Close pairs then finish in one
round, and a more distant pair fails it late enough for its estimate (see `edit_distance`) to land
near the distance; starting from `BAND_START` alone failed early, on a noisier estimate.
"""

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

        Full groups of `LANES` take the staggered vector sweep when the span is wide enough for it,
        then one group of `NARROW_LANES` if that many are left, then two or three words staggered in
        scalar registers; a single word left goes on its own.
        """
        var word = first_word
        if end_column - first_column >= 2 * LANES:
            while word + LANES <= end_word:
                self.block[LANES](word, first_column, end_column)
                word += LANES
            if word + NARROW_LANES <= end_word:
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
        """`lanes` words through `[first_column, end_column)`, staggered a column apart.

        Lane `k` holds word `first_word + lanes - 1 - k` and works column `offset + 1 + k`, so the
        lowest word sits in lane zero and works furthest behind. The stagger fills over the first
        `lanes - 1` steps, while the lower lanes still stand left of the tile, and drains over the
        last ones; a lane outside the tile keeps its state, and reads bases from the profile's
        padding. Needs at least `2 * lanes` columns, which `words` checks.
        """
        var row_low = SIMD[DType.uint64, lanes]()
        var row_high = SIMD[DType.uint64, lanes]()
        var vertical_plus = SIMD[DType.uint64, lanes]()
        var vertical_minus = SIMD[DType.uint64, lanes]()
        var horizontal_plus = SIMD[DType.uint64, lanes]()
        var horizontal_minus = SIMD[DType.uint64, lanes]()
        var lane_columns = SIMD[DType.int, lanes]()
        comptime for lane in range(lanes):
            var word = first_word + lanes - 1 - lane
            row_low[lane] = self.row_low[unsafe_offset=word]
            row_high[lane] = self.row_high[unsafe_offset=word]
            vertical_plus[lane] = self.vertical_plus[unsafe_offset=word]
            vertical_minus[lane] = self.vertical_minus[unsafe_offset=word]
            lane_columns[lane] = lane + 1

        @always_inline
        def stagger[
            masked: Bool, enters: Bool, leaves: Bool
        ](offset: Int) {
            imm self,
            imm row_low,
            imm row_high,
            mut vertical_plus,
            mut vertical_minus,
            mut horizontal_plus,
            mut horizontal_minus,
            imm lane_columns,
            imm first_column,
            imm end_column,
        }:
            # Each lane takes the difference the lane above sent last step; the top lane takes the
            # one entering the block from above, while its column is still inside the tile.
            horizontal_plus = horizontal_plus.rotate_left[1]()
            horizontal_minus = horizontal_minus.rotate_left[1]()
            comptime if enters:
                horizontal_plus[lanes - 1] = self.horizontal_plus[unsafe_offset=offset + lanes]
                horizontal_minus[lanes - 1] = self.horizontal_minus[unsafe_offset=offset + lanes]
            # The bases lane `k` needs sit side by side, from column `offset + 1`.
            var low = self.column_low.unsafe_offset(offset + 1).unsafe_load[width=lanes]()
            var high = self.column_high.unsafe_offset(offset + 1).unsafe_load[width=lanes]()
            var matches = (low ^ row_low) & (high ^ row_high)
            comptime if masked:
                var kept_plus = vertical_plus
                var kept_minus = vertical_minus
                advance[lanes](horizontal_plus, horizontal_minus, vertical_plus, vertical_minus, matches)
                var columns = lane_columns + offset
                var inside = columns.ge(first_column) & columns.lt(end_column)
                vertical_plus = inside.select(vertical_plus, kept_plus)
                vertical_minus = inside.select(vertical_minus, kept_minus)
            else:
                advance[lanes](horizontal_plus, horizontal_minus, vertical_plus, vertical_minus, matches)
            # The lowest word's difference is final, and leaves the block at its column.
            comptime if leaves:
                self.horizontal_plus[unsafe_offset=offset + 1] = horizontal_plus[0]
                self.horizontal_minus[unsafe_offset=offset + 1] = horizontal_minus[0]

        for offset in range(first_column - lanes, first_column - 1):
            stagger[True, True, False](offset)
        for offset in range(first_column - 1, end_column - lanes):
            stagger[False, True, True](offset)
        for offset in range(end_column - lanes, end_column - 1):
            stagger[True, False, True](offset)

        comptime for lane in range(lanes):
            var word = first_word + lanes - 1 - lane
            self.vertical_plus[unsafe_offset=word] = vertical_plus[lane]
            self.vertical_minus[unsafe_offset=word] = vertical_minus[lane]

    def scalar_block[lanes: Int](self, first_word: Int, first_column: Int, end_column: Int):
        """`block` in general-purpose registers, one scalar chain per word, staggered the same way.

        Each word's horizontal difference passes to the word below through a register rather than a
        vector rotation, and scalar operations take one cycle where vector ones take two, so a band
        a few words tall, which is all latency, runs faster here than across vector lanes.
        """
        for j in range(lanes):
            self.run(first_word + j, first_column, first_column + lanes - j)

        # Word `j` from the top: its row planes, its vertical differences, and the horizontal
        # difference it last sent down, which the word below takes next.
        var row_low = Array[UInt64, lanes](fill=0)
        var row_high = Array[UInt64, lanes](fill=0)
        var vertical_plus = Array[UInt64, lanes](fill=0)
        var vertical_minus = Array[UInt64, lanes](fill=0)
        var horizontal_plus = Array[UInt64, lanes](fill=0)
        var horizontal_minus = Array[UInt64, lanes](fill=0)
        comptime for j in range(lanes):
            row_low[j] = self.row_low[unsafe_offset=first_word + j]
            row_high[j] = self.row_high[unsafe_offset=first_word + j]
            vertical_plus[j] = self.vertical_plus[unsafe_offset=first_word + j]
            vertical_minus[j] = self.vertical_minus[unsafe_offset=first_word + j]
            horizontal_plus[j] = self.horizontal_plus[unsafe_offset=first_column + lanes - 1 - j]
            horizontal_minus[j] = self.horizontal_minus[unsafe_offset=first_column + lanes - 1 - j]

        for offset in range(first_column, end_column - lanes):
            # Bottom word first, so each word reads what the one above sent before it is replaced.
            comptime for i in range(lanes):
                comptime j = lanes - 1 - i
                var column = offset + lanes - j
                var hp: UInt64
                var hm: UInt64
                comptime if j == 0:
                    hp = self.horizontal_plus[unsafe_offset=column]
                    hm = self.horizontal_minus[unsafe_offset=column]
                else:
                    hp = horizontal_plus[j - 1]
                    hm = horizontal_minus[j - 1]
                var vp = vertical_plus[j]
                var vm = vertical_minus[j]
                var matches = (self.column_low[unsafe_offset=column] ^ row_low[j]) & (
                    self.column_high[unsafe_offset=column] ^ row_high[j]
                )
                advance[1](hp, hm, vp, vm, matches)
                vertical_plus[j] = vp
                vertical_minus[j] = vm
                comptime if j == lanes - 1:
                    self.horizontal_plus[unsafe_offset=column] = hp
                    self.horizontal_minus[unsafe_offset=column] = hm
                else:
                    horizontal_plus[j] = hp
                    horizontal_minus[j] = hm

        comptime for j in range(lanes):
            comptime if j != lanes - 1:
                self.horizontal_plus[unsafe_offset=end_column - 1 - j] = horizontal_plus[j]
                self.horizontal_minus[unsafe_offset=end_column - 1 - j] = horizontal_minus[j]
            self.vertical_plus[unsafe_offset=first_word + j] = vertical_plus[j]
            self.vertical_minus[unsafe_offset=first_word + j] = vertical_minus[j]

        for j in range(lanes):
            self.run(first_word + j, end_column - j, end_column)


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
    """Both sequences as the bit planes `Sweep` reads; see `Sweep` for the encoding."""

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
        # set for `G` and `T`, the code's high bit, and bit 1 differs from bit 2 for `C` and `T`, its low bit.
        # The planes run `COLUMN_PADDING` columns past either end, zero, for the lanes a staggered
        # block holds outside its tile (see `Sweep.block`).
        var padded = self.columns + 2 * COLUMN_PADDING
        self.column_low = List[UInt64](capacity=padded)
        self.column_high = List[UInt64](capacity=padded)
        self.column_codes = List[UInt8](capacity=self.columns)
        self.column_low.resize(unsafe_uninit_length=padded)
        self.column_high.resize(unsafe_uninit_length=padded)
        self.column_codes.resize(unsafe_uninit_length=self.columns)
        for index in range(COLUMN_PADDING):
            self.column_low[index] = 0
            self.column_high[index] = 0
            self.column_low[padded - 1 - index] = 0
            self.column_high[padded - 1 - index] = 0
        var first_bytes = first.unsafe_ptr()
        var low = self.column_low.unsafe_ptr().unsafe_offset(COLUMN_PADDING)
        var high = self.column_high.unsafe_ptr().unsafe_offset(COLUMN_PADDING)
        var column_codes = self.column_codes.unsafe_ptr()
        comptime CHUNK = 16
        var column = 0
        while column + CHUNK <= self.columns:
            var bytes: SIMD[DType.uint8, CHUNK]
            if reverse:
                bytes = first_bytes.unsafe_offset(self.columns - CHUNK - column).unsafe_load[width=CHUNK]().reversed()
            else:
                bytes = first_bytes.unsafe_offset(column).unsafe_load[width=CHUNK]()
            var low_bits = ((bytes >> 1) ^ (bytes >> 2)) & 1
            var high_bits = (bytes >> 2) & 1
            column_codes.unsafe_offset(column).unsafe_store(low_bits | (high_bits << 1))
            low.unsafe_offset(column).unsafe_store(UInt64(0) - low_bits.cast[DType.uint64]())
            high.unsafe_offset(column).unsafe_store(UInt64(0) - high_bits.cast[DType.uint64]())
            column += CHUNK
        while column < self.columns:
            var byte = first_bytes[unsafe_offset=self.columns - 1 - column if reverse else column]
            var low_bit = ((byte >> 1) ^ (byte >> 2)) & 1
            var high_bit = (byte >> 2) & 1
            column_codes[unsafe_offset=column] = low_bit | (high_bit << 1)
            low[unsafe_offset=column] = UInt64(0) - low_bit.cast[DType.uint64]()
            high[unsafe_offset=column] = UInt64(0) - high_bit.cast[DType.uint64]()
            column += 1

        # Eight bases at a time as one word, each byte's bit packed into a byte of the plane by a
        # multiply whose partial products never overlap.
        comptime ONES = UInt64(0x0101010101010101)
        comptime GATHER = UInt64(0x0102040810204080)
        var second_bytes = second.unsafe_ptr()
        self.row_low = List[UInt64](length=self.words, fill=0)
        self.row_high = List[UInt64](length=self.words, fill=0)
        self.row_codes = List[UInt8](capacity=self.rows)
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
            var shift = UInt64(row % WORD_BITS)
            # Stored negated, so a row matches a column where both planes XOR to ones.
            self.row_low[row // WORD_BITS] |= ((((low_bits * GATHER) >> 56) ^ 0xFF)) << shift
            self.row_high[row // WORD_BITS] |= ((((high_bits * GATHER) >> 56) ^ 0xFF)) << shift
            row += 8
        while row < self.rows:
            var byte = second_bytes[unsafe_offset=self.rows - 1 - row if reverse else row]
            var low_bit = ((byte >> 1) ^ (byte >> 2)) & 1
            var high_bit = (byte >> 2) & 1
            row_codes[unsafe_offset=row] = low_bit | (high_bit << 1)
            var shift = UInt64(row % WORD_BITS)
            self.row_low[row // WORD_BITS] |= (low_bit ^ 1).cast[DType.uint64]() << shift
            self.row_high[row // WORD_BITS] |= (high_bit ^ 1).cast[DType.uint64]() << shift
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


def pruned_distance[
    record: Bool
](mut profile: Profile, threshold: Int, stop_column: Int, mut trail: Trail, mut edge: Edge) -> Round:
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
    """
    comptime if record:
        trail.clear()
    var columns = profile.columns
    var rows = profile.rows
    var words = profile.words
    var difference = rows - columns
    var extra = (threshold - abs(difference)) // 2
    # One tile's width of horizontal edge, since each tile starts again from `+1` above (see `Sweep.shifted`);
    # the last tile may absorb a sliver of up to `2 * LANES` more columns.
    var frontier = Frontier(BAND_COLUMNS + 2 * LANES, words)
    var sweep = frontier.sweep(profile)
    var bounds = tile_bounds(stop_column, NARROW_COLUMNS if threshold < NARROW_BAND else BAND_COLUMNS)

    var top = 0
    var end_word = 0
    # The score at the top of word `top`, on the left edge of the coming tile.
    var anchor = 0
    # The deepest row a kept cell, scoring `floor` at the least, sits on, after the last tile.
    var deepest = 0
    var floor = 0
    for tile in range(len(bounds) - 1):
        var first_column = bounds[tile]
        var end_column = bounds[tile + 1]
        var width = end_column - first_column
        # Ukkonen's band bounds the bottom, and so does the deepest kept row: going `k` rows past it
        # costs at least `k` over the score there, at least `floor`, and below the diagonal that leads
        # to the end each such row also adds one to the gap still to close.
        var band_row = end_column + max(0, difference) + extra
        var slack = threshold - floor
        var end_diagonal = rows - (columns - end_column)
        var reach_row = min(
            band_row, deepest + width + slack, max(end_diagonal, (slack + deepest + width + end_diagonal) // 2), rows
        )
        var reach = (reach_row - 1) // WORD_BITS + 1
        # Fewer words than a narrow vector go through scalar registers as they are; more grow to whole
        # groups of `NARROW_LANES` from the top, which keeps every word on the vector path, since
        # `words` takes groups of `LANES` and then at most one of `NARROW_LANES`. The bottom grows
        # to fill the last group rather than the top rising, which would revisit words already left.
        var span = reach - top
        if span > NARROW_LANES:
            span = ceildiv(span, NARROW_LANES) * NARROW_LANES
        end_word = max(end_word, min(top + span, words))
        if end_word <= top:
            return Round(-1, first_column)
        comptime if record:
            trail.record(first_column, end_column, top, end_word, anchor, frontier)
        frontier.restart_horizontal(width)
        sweep.shifted(first_column).words(top, end_word, first_column, end_column)
        anchor += width

        # Read the right edge back as scores and keep, as A*PA2 does, the rows from the first to the
        # last whose score plus gap to the end fits the bound. Both scores and gaps change by at most
        # one per row, so their sum by at most two, and a row `x` over the bound rules out the next
        # `ceil(x / 2)` rows without reading them.
        edge.capture(top, end_word, anchor, frontier, rows)
        var to_end = rows - (columns - end_column)
        var first_kept = edge.low_row
        var last_kept = edge.high_row
        while first_kept <= last_kept:
            var over = edge.score(first_kept) + abs(first_kept - to_end) - threshold
            if over <= 0:
                break
            first_kept += (over + 1) // 2
        while last_kept >= first_kept:
            var over = edge.score(last_kept) + abs(last_kept - to_end) - threshold
            if over <= 0:
                break
            last_kept -= (over + 1) // 2
        if first_kept > last_kept:
            return Round(-1, end_column)

        # Move the top down to the word whose top row is at or above the first kept row, carrying the
        # anchor past the words it leaves.
        var new_top = first_kept // WORD_BITS
        for word in range(top, new_top):
            anchor += word_value(frontier.vertical_plus[word], frontier.vertical_minus[word])
        top = new_top
        # The bottom-most kept cell bounds every path below it: from there each row past the diagonal
        # costs one more, so its exact score is the floor the next tile's reach is measured from.
        deepest = last_kept
        floor = edge.score(last_kept)

    edge.capture(top, end_word, anchor, frontier, rows)
    if stop_column < columns:
        return Round(-1, stop_column)
    if end_word < words:
        return Round(-1, columns)
    return Round(edge.score(rows), columns)


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
    mut forward_trail: Trail,
    mut backward_trail: Trail,
    mut forward_edge: Edge,
    mut backward_edge: Edge,
) -> Outcome:
    """Band doubling, from the start alone or from both ends at once to meet in the middle.

    Each round sweeps the band for one bound. Meeting in the middle, one thread sweeps forward to
    the middle column while another sweeps the reversed sequences forward to the same column, which
    is the original's suffix swept from its end; the distance is the least, over the middle's rows,
    of the two scores there. Every computed score is the cost of a real path and the optimal path's
    are exact (see `pruned_distance`), so that least is exact whenever it is within the bound.

    A failed round leaves either a real alignment's cost, which caps the next bound, or the column
    at which a half pruned every row, from which the distance is estimated as the bound scaled to the
    whole width; the next bound aims just past the estimate rather than doubling. Once the band would
    cover the matrix, `give_up_wide` hands back to the caller, else the bound is lifted past any path.
    """
    var columns = forward.columns
    var rows = forward.rows
    var middle = columns // 2 if meet else columns
    var gap = abs(rows - columns)
    var threshold = gap + max(BAND_START, columns // FIRST_GUESS_DIVISOR)
    var best = Int.MAX
    while True:
        if 2 * (threshold + BAND_COLUMNS) >= rows:
            if give_up_wide:
                return Outcome(-1, middle, rows)
            threshold = max(threshold, columns + rows)
        var forward_round = Round(-1, 0)
        var backward_round = Round(-1, 0)
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
                imm threshold,
                imm middle,
                imm columns,
            }:
                if index == 0:
                    forward_round = pruned_distance[record](forward, threshold, middle, forward_trail, forward_edge)
                else:
                    backward_round = pruned_distance[record](
                        backward, threshold, columns - middle, backward_trail, backward_edge
                    )

            parallelize(half, 2, 2)
        else:
            forward_round = pruned_distance[record](forward, threshold, columns, forward_trail, forward_edge)

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
        if found >= 0 and found <= threshold:
            return Outcome(found, middle, meeting_row)

        var next = 2 * threshold
        if found >= 0:
            # Still the cost of a real alignment, so it caps every later bound.
            best = min(best, found)
        else:
            # A half pruned every row `reached` columns into its sweep, where the best alignment's
            # score had passed the bound; scaled to the whole width, that estimates the distance.
            var estimate = Int.MAX
            if forward_round.reached < middle and forward_round.reached > 0:
                estimate = min(estimate, threshold * columns // forward_round.reached)
            if meet and backward_round.reached < columns - middle and backward_round.reached > 0:
                estimate = min(estimate, threshold * columns // backward_round.reached)
            if estimate != Int.MAX:
                next = max(threshold + threshold // 4, estimate + estimate // 8 + BAND_START)
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
    var meet = workers > 1 and forward.columns >= MEET_COLUMNS
    var backward = Profile(first, second, reverse=True) if meet else Profile(String(), String())
    var forward_trail = Trail()
    var backward_trail = Trail()
    var forward_edge = Edge()
    var backward_edge = Edge()
    var outcome = band_doubling[False](
        forward, backward, meet, True, forward_trail, backward_trail, forward_edge, backward_edge
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

    Cost `s` keeps its live diagonals, offsets `lows[s] ..= lows[s] + widths[s] - 1` from the start's,
    from position `starts[s]` of the flat buffers.
    """

    var landed: List[Int]
    var reached: List[Int]
    var how: List[UInt8]
    var starts: List[Int]
    var lows: List[Int]
    var widths: List[Int]
    var segment: List[UInt8]

    def __init__(out self):
        # Room for a typical tile up front: growing from empty cost more than the search itself.
        comptime DIAGONALS = 1024
        comptime COSTS = 64
        self.landed = List[Int](capacity=DIAGONALS)
        self.reached = List[Int](capacity=DIAGONALS)
        self.how = List[UInt8](capacity=DIAGONALS)
        self.starts = List[Int](capacity=COSTS)
        self.lows = List[Int](capacity=COSTS)
        self.widths = List[Int](capacity=COSTS)
        self.segment = List[UInt8](capacity=2 * TILE_COLUMNS)


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
    comptime UNREACHED = Int.MAX

    fronts.landed.clear()
    fronts.reached.clear()
    fronts.how.clear()
    fronts.starts.clear()
    fronts.lows.clear()
    fronts.widths.clear()

    @always_inline
    def ends_here(column: Int, offset: Int, cost: Int) {imm edge, imm first_column, imm score, imm home} -> Bool:
        if column != first_column:
            return False
        var row = column + home + offset
        return row >= edge.low_row and row <= edge.high_row and edge.score(row) + cost == score

    var finish_cost = -1
    var finish_offset = 0
    fronts.starts.append(0)
    fronts.lows.append(0)
    fronts.widths.append(1)
    fronts.landed.append(end_column)
    fronts.reached.append(slide(first, second, end_column, home, first_column))
    fronts.how.append(DIAGONAL)
    if ends_here(fronts.reached[0], 0, 0):
        finish_cost = 0
    var low = 0
    var high = 0
    var cost = 1
    var halfway = (first_column + end_column) // 2
    while finish_cost < 0 and cost <= limit:
        var previous_start = fronts.starts[cost - 1]
        # The previous row is laid out from its own first diagonal; only `low ..= high` of it is live.
        var previous_base = fronts.lows[cost - 1]
        var previous_low = low
        var previous_high = high
        low -= 1
        high += 1
        var row_start = len(fronts.reached)
        fronts.starts.append(row_start)
        fronts.lows.append(low)
        fronts.widths.append(high - low + 1)
        var furthest = UNREACHED
        var nearest_column = UNREACHED
        for offset in range(low, high + 1):
            var best = UNREACHED
            var move = DIAGONAL
            var diagonal = home + offset
            # Substitution, along the same diagonal.
            if offset >= previous_low and offset <= previous_high:
                var source = fronts.reached[previous_start + offset - previous_base]
                if source != UNREACHED and source > first_column and source + diagonal > 0:
                    best = source - 1
            # A base of the first sequence against a gap, from the diagonal below.
            if offset - 1 >= previous_low and offset - 1 <= previous_high:
                var source = fronts.reached[previous_start + offset - 1 - previous_base]
                if source != UNREACHED and source > first_column and source - 1 < best:
                    best = source - 1
                    move = LEFT
            # A base of the second sequence against a gap, from the diagonal above.
            if offset + 1 >= previous_low and offset + 1 <= previous_high:
                var source = fronts.reached[previous_start + offset + 1 - previous_base]
                if source != UNREACHED and source + diagonal + 1 > 0 and source < best:
                    best = source
                    move = UP
            fronts.landed.append(best)
            fronts.how.append(move)
            if best == UNREACHED:
                fronts.reached.append(UNREACHED)
                continue
            var slid = slide(first, second, best, diagonal, first_column)
            fronts.reached.append(slid)
            furthest = min(furthest, 2 * slid - offset)
            nearest_column = min(nearest_column, slid)
            if finish_cost < 0 and ends_here(slid, offset, cost):
                finish_cost = cost
                finish_offset = offset
        if finish_cost >= 0:
            break
        if furthest == UNREACHED:
            return -1
        if 2 * cost >= limit and nearest_column > halfway:
            return -1

        @always_inline
        def dropped(offset: Int) {imm fronts, imm row_start, imm low, imm furthest, imm first_column} -> Bool:
            var at = fronts.reached[row_start + offset - low]
            return at == UNREACHED or at <= first_column or 2 * at - offset > furthest + FRONT_DROP

        # Shrink the live diagonals from both ends, keeping the stored row's layout intact.
        var new_low = low
        var new_high = high
        while new_low < new_high and dropped(new_low):
            new_low += 1
        while new_high > new_low and dropped(new_high):
            new_high -= 1
        if dropped(new_low):
            return -1
        low = new_low
        high = new_high
        cost += 1
    if finish_cost < 0:
        return -1

    # Walk the costs back down, collecting moves left to right, then hand them over right to left.
    fronts.segment.clear()
    var offset = finish_offset
    var at_cost = finish_cost
    while True:
        var index = fronts.starts[at_cost] + offset - fronts.lows[at_cost]
        for _ in range(fronts.landed[index] - fronts.reached[index]):
            fronts.segment.append(DIAGONAL)
        if at_cost == 0:
            break
        var move = fronts.how[index]
        fronts.segment.append(move)
        if move == LEFT:
            offset -= 1
        elif move == UP:
            offset += 1
        at_cost -= 1
    for index in range(len(fronts.segment) - 1, -1, -1):
        moves.append(fronts.segment[index])
    return first_column + home + finish_offset


def recomputed_segment(
    profile: Profile, trail: Trail, tile: Int, end_column: Int, end_row: Int, score: Int, mut moves: List[UInt8]
) -> Int:
    """Traces a tile back cell by cell, after sweeping it again and keeping every column's differences.

    The fallback for a tile with more edits than the wavefront search is allowed: the tile is
    recomputed one column at a time from its recorded left edge, its scores are read back from the
    differences, and each step takes a neighbour whose score plus the step's cost is the current
    score. Appends the moves right to left and returns the left-edge row.
    """
    var first_column = trail.first_columns[tile]
    var top = trail.tops[tile]
    var count = trail.ends[tile] - top
    var width = end_column - first_column
    var offset = trail.offsets[tile]
    var plus = List[UInt64](length=(width + 1) * count, fill=0)
    var minus = List[UInt64](length=(width + 1) * count, fill=0)
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

    # Scores at each word's top on every column: the band's top scores the anchor plus one per column.
    var bases = List[Int](length=(width + 1) * (count + 1), fill=0)
    for step in range(width + 1):
        var running = trail.anchors[tile] + step
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

    var lowest = top * WORD_BITS
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
    # Moves left to right, from the origin to the corner.
    var path = List[UInt8](capacity=columns + rows)
    var distance: Int
    if columns == 0 or rows == 0:
        distance = columns + rows
        for _ in range(columns):
            path.append(LEFT)
        for _ in range(rows):
            path.append(UP)
    else:
        # Asking for the thread count is a system call, so only a pair long enough to split asks.
        var meet = columns >= MEET_COLUMNS and max(threads.or_else(hardware_threads()), 1) > 1
        var backward = Profile(first, second, reverse=True) if meet else Profile(String(), String())
        var forward_trail = Trail(columns)
        var backward_trail = Trail(columns)
        var forward_edge = Edge()
        var backward_edge = Edge()
        var outcome = band_doubling[True](
            forward, backward, meet, False, forward_trail, backward_trail, forward_edge, backward_edge
        )
        distance = outcome.distance
        var forward_moves = List[UInt8](capacity=columns + rows)
        var backward_moves = List[UInt8]()
        if meet:
            var middle = outcome.middle
            var meeting_row = outcome.meeting_row
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
        for index in range(len(forward_moves) - 1, -1, -1):
            path.append(forward_moves[index])
        for move in backward_moves:
            path.append(move)

    # Branch-free: a move other than `UP` takes the next base of the first sequence, one other than
    # `LEFT` the next base of the second, and the other row of the pair gets a gap.
    var length = len(path)
    var top_row = List[UInt8](capacity=length)
    var bottom_row = List[UInt8](capacity=length)
    top_row.resize(unsafe_uninit_length=length)
    bottom_row.resize(unsafe_uninit_length=length)
    var first_bytes = first.unsafe_ptr()
    var second_bytes = second.unsafe_ptr()
    var top = top_row.unsafe_ptr()
    var bottom = bottom_row.unsafe_ptr()
    var moves = path.unsafe_ptr()
    var column = 0
    var row = 0
    comptime GAP = UInt8(ord("-"))
    if columns == 0 or rows == 0:
        # All gaps against one sequence, with nothing of the other to read.
        for index in range(length):
            top[unsafe_offset=index] = first_bytes[unsafe_offset=index] if columns > 0 else GAP
            bottom[unsafe_offset=index] = second_bytes[unsafe_offset=index] if rows > 0 else GAP
        length = 0
    for index in range(length):
        var move = moves[unsafe_offset=index]
        var takes_first = move != UP
        var takes_second = move != LEFT
        # Past either end only on the move that does not read it, so the read stays in bounds.
        var first_byte = first_bytes[unsafe_offset=min(column, columns - 1)]
        var second_byte = second_bytes[unsafe_offset=min(row, rows - 1)]
        top[unsafe_offset=index] = first_byte if takes_first else GAP
        bottom[unsafe_offset=index] = second_byte if takes_second else GAP
        column += Int(takes_first)
        row += Int(takes_second)
    return AlignmentResult(Int32(distance), String(unsafe_from_utf8=top_row), String(unsafe_from_utf8=bottom_row))


# endregion Traceback
