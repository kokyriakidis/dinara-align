"""
The global affine-gap score by wavefront, as WFA computes it, for a table with one match score and
one mismatch score.

The Gotoh recurrence maximizes a score with a reward for every match, which a wavefront cannot
search: it grows by cost, and a match must cost nothing. For a global alignment the reward folds
away (Eizenga and Lindquist, 2022). A path over sequences of `n` and `m` letters with `M`
matches, `X` mismatches, `g` gap runs and `G` gapped letters spends `2M + 2X + G = n + m`, so with
match `a`, mismatch `b`, and a run of `k` costing `o + k e` under the gap penalties,

    a (n + m) - 2 score = 2 (a - b) X + 2 o g + (2 e + a) G,

a sum of non-negative penalties that are zero for a match. Minimizing it is maximizing the score:
mismatch `x = 2 (a - b)`, opening `2 o` and extension `2 e + a`, all divided by their common factor.

The wavefront then runs as in WFA (Marco-Sola et al., 2021): for each cost `s`, three fronts hold
the furthest column every diagonal reaches at exactly that cost, ending in an alignment of two
letters or a gap either way, each built from the fronts `x`, `o + e` and `e` costs back, and the
alignment front slides over matches for free. The work grows with the square of the cost, not with
the matrix, so a close pair costs a small share of a full sweep.

It runs from both ends at once, as BiWFA does (Marco-Sola et al., 2023): a search from the origin
and one from the corner over the reversed sequences, a cost at a time each, until where they
overlap proves the optimum (see `bidirectional`). Each grows to about half the cost, so the two
together step about half the diagonals one search would.

A gap may also cost the least of two affine costs, WFA's two-piece gap-affine model, where one piece
opens cheaply and extends dearly and the other the reverse, so a long gap costs about its length
rather than a short gap's extension a letter. Each piece then has its two gap fronts of its own, built
from its own opening and extension, and the alignment front takes the best of all four; the search,
the meeting inside a gap, which pays that piece's opening once, and the split are otherwise the same.
The number of pieces is a parameter of the search (see `Wavefront`), so a single piece compiles to the
three fronts alone.

For an alignment both searches keep, of every cost, the alignment front's columns and a byte of
which source won each layer (see `History`), five bytes a diagonal where WFA's high-memory mode
keeps twelve, and the path is traced from where they met back to each end (see `trace`). A pair
whose fronts would grow past `HISTORY_LIMIT` is split where an optimal path crosses, which the two
searches find keeping only their last few costs, and each piece is aligned the same way (see
`solve`), so the memory stays bounded whatever the pair.
"""

from std.bit import count_trailing_zeros
from std.math import gcd

from .errors import AlignmentError, ErrorKind
from .slides import GATHERED_SLIDES, gathered_slides, slide
from .traceback import cigar_string, EditPath

comptime Slot = MutPointer[Int32, MutUntrackedOrigin]
"""A front read or written by `step`, indexed by diagonal: a source and its destination may share a row."""

comptime UNREACHED = Int32(-(1 << 28))
"""A diagonal no path of the cost reaches: far enough below zero that a few more columns stay negative."""

comptime PADDING = 16
"""Sentinels after each sequence's codes, so a slide reads eight letters at a time past the end."""

comptime FIRST_SENTINEL = UInt8(0xFE)
"""Ends the first sequence's codes; never a letter, and never the second's sentinel."""

comptime SECOND_SENTINEL = UInt8(0xFF)
"""Ends the second sequence's codes."""

comptime LANES = 8
"""Diagonals a front step computes at once."""

comptime CHECK_STRIDE = 64
"""Costs between the search's checks of its projected work."""

comptime CHECK_START = 128
"""The cost from which the search judges its projection, so it has a few edits to go on."""

comptime HISTORY_LIMIT = 1 << 24
"""Diagonals the two searches' kept fronts may hold together, 80 MB, past which the pair is split
(see `solve`)."""

comptime HISTORY_KEPT_PER_LETTER = 8
"""Diagonals a recording search makes room for per letter of the pair before its kept fronts first grow."""

comptime ALIGNED = 0
"""The front ending in two letters aligned, a match or a substitution."""
comptime FIRST_GAP = 1
"""The front ending in a letter of the first sequence against a gap of the first piece; a move that
consumes a letter of the first sequence alone."""
comptime SECOND_GAP = 2
"""The front ending in a letter of the second sequence against a gap of the first piece; a move that
consumes a letter of the second sequence alone."""
comptime MAX_PIECES = 2
"""Gap pieces a cost may have: a gap of `k` letters costs the least of `opening + k extension` over them."""


@always_inline
def layers_of[pieces: Int]() -> Int:
    """Fronts a cost keeps: the alignment front and two gap layers a piece, one each way."""
    return 1 + 2 * pieces


@always_inline
def gap_layer(piece: Int, along_first: Bool) -> Int:
    """The layer of `piece`'s gaps of the first sequence's letters, or of the second's: `FIRST_GAP` and
    `SECOND_GAP` for the first piece, the next two for the second."""
    return 1 + 2 * piece + (0 if along_first else 1)


@always_inline
def piece_of(layer: Int) -> Int:
    return (layer - 1) // 2


@always_inline
def along_first(layer: Int) -> Bool:
    """Whether a gap layer's letters are the first sequence's, its moves `FIRST_GAP`."""
    return (layer - 1) % 2 == 0


comptime CELLS_PER_STEP = 4
"""Cells of the vectorized full sweep (see `vector_score`) one diagonal step of the three fronts
costs about as much as: about 1.5 against 0.4 ns."""

comptime MET = 0
"""`bidirectional`'s answer when the searches proved where an optimal path splits."""
comptime OVER = 1
"""Its answer when they proved every path dearer than the ceiling."""
comptime HALTED = 2
"""Its answer when they stopped first: a full sweep would be cheaper, or the kept fronts too large."""

comptime FREE_START = 0
"""A search's origin as a whole alignment's: any first move, each at its own cost. An origin of a gap
layer `g` instead lies inside a gap of that layer, which a first move along it continues for an
extension alone: the piece after a split inside that gap."""
comptime OPENING = 8
"""An origin of `OPENING + g` must open a gap of layer `g` with its first move: the backward search of
the piece before a split inside that gap, which must end in it."""


@fieldwise_init
struct Penalties(ImplicitlyCopyable, TrivialRegisterPassable):
    """The wavefront's costs for one scoring, and what turns a cost back into a score. A second gap
    piece, `opening2` and `extension2`, counts only where a search runs two (see `MAX_PIECES`)."""

    var mismatch: Int
    var opening: Int
    """Charged once a gap run, on top of `extension` for its first letter."""
    var extension: Int
    var scale: Int
    """The common factor the costs were divided by."""
    var reward: Int
    """The match score, which every letter of both sequences earns half of before the costs."""
    var opening2: Int
    var extension2: Int

    def score(self, cost: Int, letters: Int) -> Int:
        """The Gotoh score of an alignment over `letters` letters in all, costing `cost`."""
        return (self.reward * letters - cost * self.scale) // 2

    @always_inline
    def opening_of(self, piece: Int) -> Int:
        return self.opening if piece == 0 else self.opening2

    @always_inline
    def extension_of(self, piece: Int) -> Int:
        return self.extension if piece == 0 else self.extension2

    def window[pieces: Int](self) -> Int:
        """The dearest single move: a substitution, or a gap's first letter in any piece."""
        var dearest = max(self.mismatch, self.opening + self.extension)
        comptime if pieces == 2:
            dearest = max(dearest, self.opening2 + self.extension2)
        return dearest

    def widest_opening[pieces: Int](self) -> Int:
        """The dearest opening, what a meeting inside a gap pays once for both halves at most."""
        comptime if pieces == 2:
            return max(self.opening, self.opening2)
        return self.opening


@fieldwise_init
struct EndsFree(ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """How many letters at each end of each sequence an alignment may leave unaligned for nothing, as
    WFA2-lib's ends-free alignment counts them: all zero is a global alignment, the first sequence's
    both ends at its length a read placed anywhere inside it, and one sequence's end with the other's
    start an overlap. Letters past an allowance pay as a gap would."""

    var first_begin: Int
    var first_end: Int
    var second_begin: Int
    var second_end: Int

    def __init__(out self):
        """A global alignment: every letter aligned or paid for."""
        self.first_begin = 0
        self.first_end = 0
        self.second_begin = 0
        self.second_end = 0

    def leading(self) -> EndsFree:
        """The allowances at the start alone, for the piece of a split before its crossing."""
        return EndsFree(self.first_begin, 0, self.second_begin, 0)

    def trailing(self) -> EndsFree:
        """The allowances at the end alone, for the piece after it."""
        return EndsFree(0, self.first_end, 0, self.second_end)


def wavefront_penalties(
    substitutions: List[Scalar[DType.int8]], alphabet_size: Int, open: Int, extend: Int
) -> Optional[Penalties]:
    """The wavefront costs for a table of one match and one mismatch score, if it is one.

    `open` and `extend` are the Gotoh gap scores, a run of `k` scoring `open + (k - 1) extend`.
    None when the table is not uniform, or when the costs would not all be positive: a free
    mismatch or extension would let a front grow without its cost growing.
    """
    if alphabet_size < 2:
        return None
    var reward = Int(substitutions[0])
    var mismatch = Int(substitutions[1])
    for row in range(alphabet_size):
        for column in range(alphabet_size):
            var expected = reward if row == column else mismatch
            if Int(substitutions[row * alphabet_size + column]) != expected:
                return None
    # A run of `k` scores `open + (k - 1) extend`, so it costs `(extend - open) - k extend`.
    var gap_opening = extend - open
    var gap_extension = -extend
    var x = 2 * (reward - mismatch)
    var o = 2 * gap_opening
    var e = 2 * gap_extension + reward
    if x <= 0 or e <= 0 or o < 0:
        return None
    var scale = gcd(gcd(x, e), o)
    return Penalties(x // scale, o // scale, e // scale, scale, reward, 0, 0)


def padded(codes: Span[UInt8, _], sentinel: UInt8, reverse: Bool) -> List[UInt8]:
    """The codes, back to front with `reverse`, and `PADDING` sentinels after them."""
    var count = len(codes)
    var out = List[UInt8](length=count + PADDING, fill=sentinel)
    var source = codes.unsafe_ptr()
    var destination = out.unsafe_ptr()
    if not reverse:
        Span(unsafe_ptr=destination, length=count).copy_from(codes)
        return out^
    comptime CHUNK = 16
    var index = 0
    while index + CHUNK <= count:
        var chunk = source.unsafe_offset(count - index - CHUNK).unsafe_load[width=CHUNK]()
        destination.unsafe_offset(index).unsafe_store(chunk.reversed())
        index += CHUNK
    while index < count:
        destination[unsafe_offset=index] = source[unsafe_offset=count - 1 - index]
        index += 1
    return out^


comptime ENTRY_MASK = UInt8(7)
"""A flag's low three bits: the layer the alignment front's entry came from, zero for a substitution."""


@always_inline
def opened_bit(layer: Int) -> UInt8:
    """A flag bit: gap layer `layer` came from an opening, not an extension."""
    return UInt8(1) << UInt8(2 + layer)


struct History(Movable):
    """What the traceback needs of every cost's fronts: the alignment front's column on each diagonal
    kept, and a flag of which source won each layer there (see `ENTRY_MASK` and `opened_bit`). Cost `s` holds diagonals `lows[s] ..= highs[s]`, its columns and flags from `starts[s]`;
    elsewhere it reads as unreached. The step writes both as it goes (see `begin`), so the fronts are
    never copied out of the rings, and they go into blocks that are never moved either: a list grown
    by doubling would copy everything kept so far each time, a tenth of an alignment's time.
    """

    var starts: List[Int]
    """Each cost's block, above `BLOCK_SHIFT`, and its first diagonal's place in the block below."""
    var lows: List[Int]
    var highs: List[Int]
    var columns: List[List[Int32]]
    var flags: List[List[UInt8]]
    var used: Int
    """Entries taken from the last block."""
    var kept: Int
    """Diagonals kept in all, what `HISTORY_LIMIT` bounds."""
    var next_block: Int
    """The entries the next block holds, doubling up to `HISTORY_BLOCK`."""
    var pending_start: Int
    """Where the cost `begin` made room for starts in the last block, and the diagonal it starts at."""
    var pending_low: Int

    def __init__(out self, capacity: Int = 0):
        """Kept fronts whose first block holds `capacity` diagonals."""
        self.starts = List[Int]()
        self.lows = List[Int]()
        self.highs = List[Int]()
        self.columns = List[List[Int32]]()
        self.flags = List[List[UInt8]]()
        self.used = 0
        self.kept = 0
        self.next_block = max(capacity, 1024)
        self.pending_start = 0
        self.pending_low = 0

    def begin(mut self, low: Int, high: Int) -> Tuple[Slot, MutPointer[UInt8, MutUntrackedOrigin]]:
        """Room for the next cost's columns and flags on `low ..= high` and a lane group past, each
        indexed by diagonal, which the step fills before `finish` keeps them."""
        var needed = high - low + 1 + LANES
        if len(self.columns) == 0 or self.used + needed > len(self.columns[len(self.columns) - 1]):
            # A block is only ever written, never filled first, so its pages arrive as the step reaches them.
            var size = max(self.next_block, needed)
            var block = List[Int32](capacity=size)
            block.resize(unsafe_uninit_length=size)
            var block_flags = List[UInt8](capacity=size)
            block_flags.resize(unsafe_uninit_length=size)
            self.columns.append(block^)
            self.flags.append(block_flags^)
            self.used = 0
            self.next_block = min(2 * self.next_block, HISTORY_BLOCK)
        var last = len(self.columns) - 1
        self.pending_start = self.used
        self.pending_low = low
        return (
            self.columns[last].unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]().unsafe_offset(self.used - low),
            self.flags[last].unsafe_ptr().unsafe_origin_cast[MutUntrackedOrigin]().unsafe_offset(self.used - low),
        )

    def finish(mut self, low: Int, high: Int):
        """Keeps the cost `begin` made room for on `low ..= high`, inside that room; none when `low > high`."""
        var block = (len(self.columns) - 1) << BLOCK_SHIFT
        self.lows.append(low)
        self.highs.append(high)
        if low > high:
            self.starts.append(block)
            return
        self.starts.append(block + self.pending_start + low - self.pending_low)
        self.used = self.pending_start + high - self.pending_low + 1
        self.kept += high - low + 1

    def skip(mut self):
        """The next cost reached nothing."""
        self.starts.append(0)
        self.lows.append(1)
        self.highs.append(0)

    @inline(.always)
    def column(self, cost: Int, diagonal: Int) -> Int:
        """The alignment front's furthest column on `diagonal` at `cost`, or `UNREACHED`."""
        if cost < 0 or cost >= len(self.lows) or diagonal < self.lows[cost] or diagonal > self.highs[cost]:
            return Int(UNREACHED)
        var start = self.starts[cost]
        return Int(self.columns[start >> BLOCK_SHIFT][(start & BLOCK_MASK) + diagonal - self.lows[cost]])

    @inline(.always)
    def flag(self, cost: Int, diagonal: Int) -> UInt8:
        """The flag of `diagonal` at `cost`, which the traceback reads only where a front reached."""
        var start = self.starts[cost]
        return self.flags[start >> BLOCK_SHIFT][(start & BLOCK_MASK) + diagonal - self.lows[cost]]


comptime ROW_PADDING = 16
"""Diagonals, one cache line, between the end of one ring row and the start of the next.

The rows' widths double from a power of two, so without it every row of a wide ring started at the
same place modulo 4 KB: a step's stores and its loads from the sources a few rows away then matched
in their low twelve address bits, which x86 checks first for a store a load might depend on, and
the rows shared their L1 sets. On the Skylake-X that blocked 1.5 billion loads over ten 100 kbp
pairs, as `ld_blocks_partial.address_alias` counts them."""

comptime HISTORY_BLOCK = 1 << 22
"""The most diagonals one block of kept fronts holds, 20 MB, unless one cost needs more."""
comptime BLOCK_SHIFT = 40
"""Where a kept cost's block number starts in its `starts` entry, above its place in the block."""
comptime BLOCK_MASK = (1 << BLOCK_SHIFT) - 1


struct Fronts[layers: Int](Movable):
    """The last `slots` costs' fronts, `layers` a cost, in rings indexed by cost, each over the
    diagonals reached.

    Diagonal `k` holds the cells whose column minus row is `k`, from `-rows` to `columns`, stored
    at `k + base` of each of the `layers slots` rows of one buffer, a slot's layers side by side. The rows
    cover only the diagonals the fronts have reached, about `2 s / e` of them at cost `s`, and double
    towards whichever side they outgrow, so the memory grows with the score rather than the
    sequences. A slot holds its cost's values on `lows[slot] ..= highs[slot]` and reads as unreached
    everywhere else, so a front reads its sources' neighbours unchecked: the rows start unreached,
    and a slot reused for a later cost has the old cost's diagonals outside the new range cleared.
    """

    var slots: Int
    var current: Int
    """The slot of the last cost grown."""
    var width: Int
    var stride: Int
    """Where each row starts after the last: `width` and a cache line more (see `ROW_PADDING`)."""
    var base: Int
    var buffer: List[Int32]
    var lows: List[Int]
    var highs: List[Int]
    var reach: List[Int]
    """Each slot's furthest anti-diagonal, `2 column - diagonal`, over its alignment front."""
    var least: Int
    """The diagonals a front may ever read: the matrix's and two either side, and a lane group past."""
    var most: Int

    def __init__(out self, slots: Int, columns: Int, rows: Int):
        self.slots = slots
        self.current = 0
        self.least = -rows - 2
        self.most = columns + LANES + 2
        self.width = 4 * LANES
        self.stride = self.width + ROW_PADDING
        self.base = self.width // 2
        self.buffer = List[Int32](length=Self.layers * slots * self.stride, fill=UNREACHED)
        self.lows = List[Int](length=slots, fill=1)
        self.highs = List[Int](length=slots, fill=0)
        self.reach = List[Int](length=slots, fill=Int.MIN // 2)

    @inline(.always)
    def back(self, lag: Int) -> Int:
        """The slot of the cost `lag` before the last one grown, `lag < slots`; before the first cost,
        one no cost has taken yet."""
        var slot = self.current - lag
        return slot + self.slots if slot < 0 else slot

    @inline(.always)
    def row(mut self, slot: Int, layer: Int) -> Slot:
        """A slot's front of one layer, indexed by diagonal."""
        return (
            self.buffer.unsafe_ptr()
            .unsafe_origin_cast[MutUntrackedOrigin]()
            .unsafe_offset((Self.layers * slot + layer) * self.stride + self.base)
        )

    def ready(mut self, low: Int, high: Int):
        """Makes every row cover `low ..= high`, within the diagonals a front may read."""
        if max(low, self.least) + self.base < 0 or min(high, self.most) + self.base >= self.width:
            self.grow(max(low, self.least), min(high, self.most))

    def grow(mut self, low: Int, high: Int):
        """Widens every row to cover `low ..= high`, doubling towards the side outgrown, and moves
        what they hold with it; the rest reads as unreached."""
        var size = self.width
        var first = -self.base
        var last = size - 1 - self.base
        var new_first = first
        var new_last = last
        if low < first:
            new_first = max(min(low, first - size), self.least)
        if high > last:
            new_last = min(max(high, last + size), self.most)
        var new_size = new_last - new_first + 1
        var shift = first - new_first
        var rows = Self.layers * self.slots
        var new_stride = new_size + ROW_PADDING
        var wider = List[Int32](length=rows * new_stride, fill=UNREACHED)
        var source = self.buffer.unsafe_ptr()
        var destination = wider.unsafe_ptr()
        for row in range(rows):
            Span(unsafe_ptr=destination.unsafe_offset(row * new_stride + shift), length=size).copy_from(
                Span(unsafe_ptr=source.unsafe_offset(row * self.stride), length=size)
            )
        self.buffer = wider^
        self.width = new_size
        self.stride = new_stride
        self.base = -new_first

    def claim(mut self, slot: Int, low: Int, high: Int):
        """Hands `slot` to a new cost over `low ..= high`, clearing what its old cost left outside."""
        var old_low = self.lows[slot]
        var old_high = self.highs[slot]
        if old_low <= old_high:
            # Only the old diagonals outside the new range: the new cost writes over the rest.
            var below_end = min(old_high, low - 1) if low <= high else old_high
            var above_start = max(old_low, high + 1) if low <= high else old_high + 1
            comptime for layer in range(Self.layers):
                var values = self.row(slot, layer)
                for diagonal in range(old_low, below_end + 1):
                    values[unsafe_offset=diagonal] = UNREACHED
                for diagonal in range(above_start, old_high + 1):
                    values[unsafe_offset=diagonal] = UNREACHED
        self.lows[slot] = low
        self.highs[slot] = high


@inline(.never)
def step[
    record: Bool, pieces: Int
](
    mismatched: ImmPointer[Int32, _],
    opening: ImmPointer[Int32, _],
    first_gaps: ImmPointer[Int32, _],
    second_gaps: ImmPointer[Int32, _],
    opening2: ImmPointer[Int32, _],
    first_gaps2: ImmPointer[Int32, _],
    second_gaps2: ImmPointer[Int32, _],
    aligned: Slot,
    opened_first: Slot,
    opened_second: Slot,
    opened_first2: Slot,
    opened_second2: Slot,
    kept: Slot,
    flags: MutPointer[UInt8, MutUntrackedOrigin],
    first: ImmPointer[UInt8, _],
    second: ImmPointer[UInt8, _],
    low: Int,
    high: Int,
    columns: Int,
    rows: Int,
) -> Int:
    """One cost's fronts on diagonals `low ..= high`, every pointer indexed by diagonal, and with
    `record` each diagonal's alignment column again in `kept` and its flag (see `ENTRY_MASK`).

    `mismatched` is the alignment front a mismatch back, `opening` the one an opened gap back, and
    the gap fronts are their own layers an extension back. A letter of the first sequence against
    a gap comes from the diagonal below and moves one column; one of the second, from the diagonal
    above, stays in its column and moves one row. Each only where it stays inside the matrix. With
    two `pieces` the second piece's sources and layers, the ones ending in `2`, step the same way;
    with one they are never read.

    Each lane group's alignment front then slides over its matches straight away, the eight slides
    independent of each other, and the furthest anti-diagonal, `2 column - diagonal`, comes back.
    """
    comptime Lanes = SIMD[DType.int32, LANES]
    var lane_diagonals = Lanes()
    comptime for lane in range(LANES):
        lane_diagonals[lane] = Int32(lane)
    var column_limit = Lanes(Int32(columns))
    var row_limit = Lanes(Int32(rows))
    var unreached = Lanes(UNREACHED)
    var reach = Int.MIN // 2
    # The gathered slides keep their furthest anti-diagonals a lane each, reduced once after the loop:
    # reducing every group put a chain of shuffles and maxima on each one's path.
    var reaches = Lanes(Int32.MIN // 2)
    var diagonal = low
    while diagonal <= high:
        var diagonals = lane_diagonals + Int32(diagonal)
        var same = mismatched.unsafe_offset(diagonal).unsafe_load[width=LANES]()
        var opened_below = opening.unsafe_offset(diagonal - 1).unsafe_load[width=LANES]()
        var extended_below = first_gaps.unsafe_offset(diagonal - 1).unsafe_load[width=LANES]()
        var opened_above = opening.unsafe_offset(diagonal + 1).unsafe_load[width=LANES]()
        var extended_above = second_gaps.unsafe_offset(diagonal + 1).unsafe_load[width=LANES]()
        var below = max(opened_below, extended_below)
        var above = max(opened_above, extended_above)
        var first_gap = below.lt(column_limit).select(below + 1, unreached)
        var second_gap = (above - diagonals).le(row_limit).select(above, unreached)
        var substituted = (same.lt(column_limit) & (same - diagonals).lt(row_limit)).select(same + 1, unreached)
        var gapped = max(first_gap, second_gap)
        opened_first.unsafe_offset(diagonal).unsafe_store(first_gap)
        opened_second.unsafe_offset(diagonal).unsafe_store(second_gap)
        var opened_below2 = unreached
        var extended_below2 = unreached
        var opened_above2 = unreached
        var extended_above2 = unreached
        var first_gap2 = unreached
        var second_gap2 = unreached
        comptime if pieces == 2:
            opened_below2 = opening2.unsafe_offset(diagonal - 1).unsafe_load[width=LANES]()
            extended_below2 = first_gaps2.unsafe_offset(diagonal - 1).unsafe_load[width=LANES]()
            opened_above2 = opening2.unsafe_offset(diagonal + 1).unsafe_load[width=LANES]()
            extended_above2 = second_gaps2.unsafe_offset(diagonal + 1).unsafe_load[width=LANES]()
            var below2 = max(opened_below2, extended_below2)
            var above2 = max(opened_above2, extended_above2)
            first_gap2 = below2.lt(column_limit).select(below2 + 1, unreached)
            second_gap2 = (above2 - diagonals).le(row_limit).select(above2, unreached)
            gapped = max(gapped, max(first_gap2, second_gap2))
            opened_first2.unsafe_offset(diagonal).unsafe_store(first_gap2)
            opened_second2.unsafe_offset(diagonal).unsafe_store(second_gap2)
        var front = aligned.unsafe_offset(diagonal)
        comptime if GATHERED_SLIDES:
            var slid = gathered_slides(first, second, max(substituted, gapped), diagonals)
            front.unsafe_store(slid)
            reaches = max(reaches, slid + slid - diagonals)
        else:
            front.unsafe_store(max(substituted, gapped))
            comptime for lane in range(LANES):
                var column = Int(front[unsafe_offset=lane])
                if column >= 0:
                    column = slide(first, second, column, diagonal + lane)
                    front[unsafe_offset=lane] = Int32(column)
                    reach = max(reach, 2 * column - diagonal - lane)
        comptime if record:
            kept.unsafe_offset(diagonal).unsafe_store(front.unsafe_load[width=LANES]())
        comptime if record:
            # Worked out in the fronts' own lanes and narrowed once.
            @always_inline
            def bit(won: SIMD[DType.bool, LANES], layer: Int) -> Lanes:
                return won.select(Lanes(Int32(opened_bit(layer))), Lanes(0))

            var entry = first_gap.ge(second_gap).select(Lanes(Int32(FIRST_GAP)), Lanes(Int32(SECOND_GAP)))
            var opened = bit(opened_below.ge(extended_below), FIRST_GAP) | bit(
                opened_above.ge(extended_above), SECOND_GAP
            )
            comptime if pieces == 2:
                var second_entry = first_gap2.ge(second_gap2).select(
                    Lanes(Int32(gap_layer(1, True))), Lanes(Int32(gap_layer(1, False)))
                )
                entry = max(first_gap, second_gap).ge(max(first_gap2, second_gap2)).select(entry, second_entry)
                opened |= bit(opened_below2.ge(extended_below2), gap_layer(1, True)) | bit(
                    opened_above2.ge(extended_above2), gap_layer(1, False)
                )
            entry = substituted.ge(gapped).select(Lanes(0), entry)
            flags.unsafe_offset(diagonal).unsafe_store((entry | opened).cast[DType.uint8]())
        diagonal += LANES
    return max(reach, Int(reaches.reduce_max()))


struct Wavefront[pieces: Int](Movable):
    """One direction's search over two encoded sequences: the fronts of each cost in turn, from cost
    zero at the origin, and, for the costs grown recording, what the traceback needs of each (see
    `History`). A gap costs the least over its `pieces`, one or two (see `MAX_PIECES`).

    The backward search is the same search over both sequences reversed, whose origin is the corner.
    """

    var first: List[UInt8]
    var second: List[UInt8]
    var columns: Int
    var rows: Int
    var penalties: Penalties
    var fronts: Fronts[layers_of[Self.pieces]()]
    var history: History
    var cost: Int
    """The last cost whose fronts are grown."""
    var furthest: Int
    """The furthest anti-diagonal, `2 column - diagonal`, any alignment front has reached."""
    var work: Int
    """Diagonals stepped so far."""
    var origin: Int
    """What the origin allows, `FREE_START` or one of the constants after it."""

    def __init__(
        out self,
        first: Span[UInt8, _],
        second: Span[UInt8, _],
        penalties: Penalties,
        origin: Int,
        record: Bool,
        reverse: Bool,
        free_first: Int = 0,
        free_second: Int = 0,
    ):
        """The fronts of cost zero over two encoded sequences, both back to front with `reverse`: the
        matches from the origin, unless its first move must open a gap, kept with `record`. Both
        sequences must hold a letter.

        With up to `free_first` letters of the first sequence or `free_second` of the second left
        unaligned at the origin's end for nothing, cost zero holds every diagonal such a start lies
        on, each from its cell on the first row or column, as WFA2-lib's ends-free alignment starts."""
        self.columns = len(first)
        self.rows = len(second)
        self.first = padded(first, FIRST_SENTINEL, reverse)
        self.second = padded(second, SECOND_SENTINEL, reverse)
        self.penalties = penalties
        self.origin = origin
        # Every source a cost reads lies at most this far back, and a slot is reused after as many.
        self.fronts = Fronts[layers_of[Self.pieces]()](penalties.window[Self.pieces]() + 1, self.columns, self.rows)
        # A search keeps a few diagonals a letter on close pairs, so its lists rarely grow on them.
        self.history = History(
            min(HISTORY_KEPT_PER_LETTER * (self.columns + self.rows), HISTORY_LIMIT // 2) if record else 0
        )
        self.cost = 0
        self.work = 0
        self.furthest = Int.MIN // 2
        self.fronts.ready(-1 - LANES, 1 + LANES)
        if origin > OPENING:
            # Nothing at cost zero: the opening gap enters at its own cost (see `advance`). The kept
            # fronts still hold the origin at column zero, where that gap's walk back ends.
            if record:
                self.history.begin(0, 0)[0][unsafe_offset=0] = 0
                self.history.finish(0, 0)
            return
        var low = -min(free_second, self.rows)
        var high = min(free_first, self.columns)
        self.fronts.ready(low - 1 - LANES, high + 1 + LANES)
        var front = self.fronts.row(0, ALIGNED)
        self.fronts.claim(0, low, high)
        # The flags of cost zero are never read: the traceback stops there.
        var kept = front
        if record:
            kept = self.history.begin(low, high)[0]
        var reach = Int.MIN // 2
        for diagonal in range(low, high + 1):
            var column = slide(self.first.unsafe_ptr(), self.second.unsafe_ptr(), max(diagonal, 0), diagonal)
            front[unsafe_offset=diagonal] = Int32(column)
            kept[unsafe_offset=diagonal] = Int32(column)
            reach = max(reach, 2 * column - diagonal)
        # Inside a gap, its layer holds the origin too, which extends without a second opening.
        if origin != FREE_START:
            self.fronts.row(0, origin)[unsafe_offset=0] = 0
        self.fronts.reach[0] = reach
        self.furthest = reach
        if record:
            self.history.finish(low, high)

    def advance[record: Bool](mut self):
        """Grows the next cost's fronts from the ring, and with `record` keeps what the traceback needs
        of them."""
        self.cost += 1
        var cost = self.cost
        var x = self.penalties.mismatch
        var columns = self.columns
        var rows = self.rows
        self.fronts.current = self.fronts.back(self.fronts.slots - 1)
        var slot = self.fronts.current
        # A source before the first cost reads a slot no cost has taken yet, which reads as unreached.
        var mismatch_slot = self.fronts.back(x)
        var low = Int.MAX
        var high = Int.MIN
        if cost >= x and self.fronts.lows[mismatch_slot] <= self.fronts.highs[mismatch_slot]:
            low = min(low, self.fronts.lows[mismatch_slot])
            high = max(high, self.fronts.highs[mismatch_slot])
        # Each piece's opening and extension sources, their gaps a diagonal either side.
        var opening_slots = Array[Int, MAX_PIECES](fill=0)
        var extension_slots = Array[Int, MAX_PIECES](fill=0)
        comptime for piece in range(Self.pieces):
            var o = self.penalties.opening_of(piece)
            var e = self.penalties.extension_of(piece)
            opening_slots[piece] = self.fronts.back(o + e)
            extension_slots[piece] = self.fronts.back(e)
            var source = opening_slots[piece]
            if cost >= o + e and self.fronts.lows[source] <= self.fronts.highs[source]:
                low = min(low, self.fronts.lows[source] - 1)
                high = max(high, self.fronts.highs[source] + 1)
            source = extension_slots[piece]
            if cost >= e and self.fronts.lows[source] <= self.fronts.highs[source]:
                low = min(low, self.fronts.lows[source] - 1)
                high = max(high, self.fronts.highs[source] + 1)
        # An origin that must open a gap reaches the gap's first letter at the opening's cost.
        var opened = 0
        var opened_layer = self.origin - OPENING
        if self.origin > OPENING:
            var piece = piece_of(opened_layer)
            if cost == self.penalties.opening_of(piece) + self.penalties.extension_of(piece):
                opened = 1 if along_first(opened_layer) else -1
                low = min(low, opened)
                high = max(high, opened)
        low = max(low, -rows)
        high = min(high, columns)
        if low > high:
            self.empty[record](slot)
            return
        # The step reads a lane group past `high` and a diagonal either side of the range.
        self.fronts.ready(low - 1, high + LANES + 1)
        self.fronts.claim(slot, low, high)
        var front = self.fronts.row(slot, ALIGNED)
        # Without `record` the step writes neither, and these stand for nothing.
        var kept = front
        var flags = front.unsafe_bitcast[UInt8]()
        comptime if record:
            var room = self.history.begin(low, high)
            kept = room[0]
            flags = room[1]
        # With one piece the second's sources and layers stand for the first's, and the step reads none.
        comptime last = Self.pieces - 1
        var reach = step[record, Self.pieces](
            self.fronts.row(mismatch_slot, ALIGNED),
            self.fronts.row(opening_slots[0], ALIGNED),
            self.fronts.row(extension_slots[0], FIRST_GAP),
            self.fronts.row(extension_slots[0], SECOND_GAP),
            self.fronts.row(opening_slots[last], ALIGNED),
            self.fronts.row(extension_slots[last], gap_layer(last, True)),
            self.fronts.row(extension_slots[last], gap_layer(last, False)),
            front,
            self.fronts.row(slot, FIRST_GAP),
            self.fronts.row(slot, SECOND_GAP),
            self.fronts.row(slot, gap_layer(last, True)),
            self.fronts.row(slot, gap_layer(last, False)),
            kept,
            flags,
            self.first.unsafe_ptr(),
            self.second.unsafe_ptr(),
            low,
            high,
            columns,
            rows,
        )
        # The step wrote whole lane groups; past `high` the slot must read unreached again, a lane group
        # a layer, within the room `ready` made.
        comptime for layer in range(layers_of[Self.pieces]()):
            self.fronts.row(slot, layer).unsafe_offset(high + 1).unsafe_store(SIMD[DType.int32, LANES](UNREACHED))
        # The gap an origin must open enters past the step, and slides as the step's columns did.
        if opened != 0:
            var column = 1 if opened == 1 else 0
            self.fronts.row(slot, opened_layer)[unsafe_offset=opened] = Int32(column)
            comptime if record:
                flags[unsafe_offset=opened] = UInt8(opened_layer) | opened_bit(opened_layer)
            if Int(front[unsafe_offset=opened]) < column:
                column = slide(self.first.unsafe_ptr(), self.second.unsafe_ptr(), column, opened)
                front[unsafe_offset=opened] = Int32(column)
                reach = max(reach, 2 * column - opened)
                comptime if record:
                    kept[unsafe_offset=opened] = Int32(column)
        self.fronts.reach[slot] = reach
        self.furthest = max(self.furthest, reach)
        self.work += high - low + 1

        # Diagonals no layer reached at either end are dropped, as WFA trims them, so the range
        # tracks the paths alive rather than every diagonal the gap costs allow.
        var base = self.fronts.row(slot, ALIGNED)
        var stride = self.fronts.stride

        @inline(.always)
        def dead(diagonal: Int) {imm base, imm stride} -> Bool:
            var reached = base[unsafe_offset=diagonal]
            comptime for layer in range(1, layers_of[Self.pieces]()):
                reached = max(reached, base[unsafe_offset=layer * stride + diagonal])
            return reached < 0

        var kept_low = low
        var kept_high = high
        while kept_low <= kept_high and dead(kept_low):
            kept_low += 1
        while kept_high >= kept_low and dead(kept_high):
            kept_high -= 1
        if kept_low > kept_high:
            self.fronts.claim(slot, 1, 0)
            self.fronts.reach[slot] = Int.MIN // 2
            comptime if record:
                self.history.finish(1, 0)
            return
        self.fronts.lows[slot] = kept_low
        self.fronts.highs[slot] = kept_high
        comptime if record:
            self.history.finish(kept_low, kept_high)

    def empty[record: Bool](mut self, slot: Int):
        """Leaves the cost just begun with no diagonal reached."""
        self.fronts.claim(slot, 1, 0)
        self.fronts.reach[slot] = Int.MIN // 2
        comptime if record:
            self.history.skip()


@fieldwise_init
struct Meeting(ImplicitlyCopyable, TrivialRegisterPassable):
    """A cell an optimal path passes, in one layer, and what it costs from each end: the forward front
    of `forward_cost` reaches it on `diagonal`, and the backward front of `backward_cost` comes back
    to it, in the same layer for a gap, which then pays its opening once for both halves."""

    var cost: Int
    var layer: Int
    var diagonal: Int
    var column: Int
    var forward_cost: Int
    var backward_cost: Int

    @staticmethod
    def none() -> Meeting:
        """No meeting yet, dearer than any."""
        return Meeting(Int.MAX, ALIGNED, 0, 0, 0, 0)


def meet[
    pieces: Int
](
    mut forward: Wavefront[pieces],
    forward_cost: Int,
    mut backward: Wavefront[pieces],
    backward_cost: Int,
    mut best: Meeting,
):
    """Lowers `best` to where the forward fronts of one cost and the backward fronts of another
    overlap, if that costs less.

    The backward front's diagonal `target - k` mirrors the forward's `k`, and its column `c` the
    forward's `columns - c`. Two alignment fronts overlap where their columns add up to `columns` or
    more: the forward path reaches a cell the backward one comes back past, and edit costs never fall
    along a diagonal, so the cell costs at most the sum. Two gap fronts of one layer overlap the same
    way, at a cell inside the gap both may hold, and join into one gap of one opening, its piece's.
    """
    var columns = forward.columns
    var rows = forward.rows
    var penalties = forward.penalties
    var total = forward_cost + backward_cost
    if total - penalties.widest_opening[pieces]() >= best.cost:
        return
    var ahead = forward.fronts.back(forward.cost - forward_cost)
    var behind = backward.fronts.back(backward.cost - backward_cost)
    # A front's furthest anti-diagonal bounds all its layers': the two must add up to the matrix's.
    if forward.fronts.reach[ahead] + backward.fronts.reach[behind] < columns + rows:
        return
    var target = columns - rows
    var low = max(forward.fronts.lows[ahead], target - backward.fronts.highs[behind])
    var high = min(forward.fronts.highs[ahead], target - backward.fronts.lows[behind])
    var aligned = forward.fronts.row(ahead, ALIGNED)
    var back_aligned = backward.fronts.row(behind, ALIGNED)
    var stride = forward.fronts.stride
    var back_stride = backward.fronts.stride

    @inline(.always)
    def check(
        diagonal: Int,
    ) {
        mut best,
        imm aligned,
        imm back_aligned,
        imm stride,
        imm back_stride,
        imm target,
        imm columns,
        imm rows,
        imm total,
        imm penalties,
        imm forward_cost,
        imm backward_cost,
    }:
        var mirrored = target - diagonal
        var column = Int(aligned[unsafe_offset=diagonal])
        if total < best.cost and column + Int(back_aligned[unsafe_offset=mirrored]) >= columns:
            best = Meeting(total, ALIGNED, diagonal, column, forward_cost, backward_cost)
        comptime for layer in range(1, layers_of[pieces]()):
            var joined = total - penalties.opening_of(piece_of(layer))
            if joined < best.cost:
                var gap = Int(aligned[unsafe_offset=layer * stride + diagonal])
                var back_gap = Int(back_aligned[unsafe_offset=layer * back_stride + mirrored])
                if gap + back_gap >= columns:
                    comptime if along_first(layer):
                        # A cell inside a gap of the first sequence's letters has consumed one of them either way.
                        var cell = min(gap, min(columns - 1, rows + diagonal))
                        if cell >= max(max(1, diagonal), columns - back_gap):
                            best = Meeting(joined, layer, diagonal, cell, forward_cost, backward_cost)
                    else:
                        # And inside a gap of the second's, one of its rows either way.
                        var cell = min(gap, rows + diagonal - 1)
                        if cell >= max(diagonal + 1, columns - back_gap):
                            best = Meeting(joined, layer, diagonal, cell, forward_cost, backward_cost)

    # A lane group at a time, and the cells of a group that might meet one by one. An alignment front
    # reaches as far as either gap front of its cost, so where no two alignment fronts overlap, no
    # two gap fronts do.
    comptime Lanes = SIMD[DType.int32, LANES]
    var needed = Lanes(Int32(columns))
    var diagonal = low
    while diagonal + LANES - 1 <= high:
        var mirrored = target - diagonal - LANES + 1
        var meets = (
            aligned.unsafe_offset(diagonal).unsafe_load[width=LANES]()
            + back_aligned.unsafe_offset(mirrored).unsafe_load[width=LANES]().reversed()
        ).ge(needed)
        if meets.reduce_or():
            for lane in range(LANES):
                check(diagonal + lane)
        diagonal += LANES
    while diagonal <= high:
        check(diagonal)
        diagonal += 1


def bidirectional[
    pieces: Int, record: Bool
](
    mut forward: Wavefront[pieces],
    mut backward: Wavefront[pieces],
    mut best: Meeting,
    give_up: Bool,
    limit: Int,
    ceiling: Int = Int.MAX,
) -> Int:
    """Grows the two searches a cost at a time each in turn, lowering `best` to the cheapest place an
    optimal path splits between them, until no cheaper one is left unchecked: `MET`. `OVER` once no
    path can cost `ceiling` or less; `HALTED` once a full sweep would be cheaper, with `give_up`, or,
    with `record`, once the kept fronts would pass `limit` entries, and the searches may resume.

    Each new front is checked against the other side's last `window` costs, which its ring holds. Every
    optimal path has a cell, between two moves or inside a gap, whose costs from the two ends differ by
    at most `window`, the dearest single move: the difference rises from `-cost` to `cost` in steps of
    at most twice that. So once both sides have passed half of `best`, plus the dearest gap opening the
    two halves both paid and `window`, no cheaper path is left unchecked. Each side then grows to about
    half the optimum, and the two together about half the diagonals one search grows alone. By the
    same count, once both have passed half of `ceiling` with that margin and found nothing within it,
    nothing is, and a capped search stops there, short of the optimum.
    """
    var columns = forward.columns
    var rows = forward.rows
    var letters = columns + rows
    var o = forward.penalties.widest_opening[pieces]()
    var window = forward.penalties.window[pieces]()
    if forward.cost == 0 and backward.cost == 0:
        meet(forward, 0, backward, 0, best)
    var budget = columns * rows // CELLS_PER_STEP
    var next_check = CHECK_START
    while True:
        var reached = 2 * min(forward.cost, backward.cost)
        # The cap first: an optimum just past it proves itself at the same reach.
        if best.cost > ceiling and reached >= ceiling + o + window:
            return OVER
        if best.cost != Int.MAX and reached >= best.cost - 1 + o + window:
            return MET
        if forward.cost <= backward.cost:
            forward.advance[record]()
            var ahead = forward.cost
            for behind in range(max(0, ahead - window, backward.cost - window), backward.cost + 1):
                meet(forward, ahead, backward, behind, best)
        else:
            backward.advance[record]()
            var behind = backward.cost
            for ahead in range(max(0, behind - window, forward.cost - window), forward.cost + 1):
                meet(forward, ahead, backward, behind, best)
        comptime if record:
            if forward.history.kept + backward.history.kept > limit:
                return HALTED
        var spent = forward.cost + backward.cost
        if give_up and spent >= next_check:
            next_check = spent + CHECK_STRIDE
            # The fronts widen with the cost, so the work grows with its square; projected from how
            # far across the two furthest anti-diagonals have come together.
            var projected = spent * letters // max(forward.furthest + backward.furthest, 1)
            var ratio = Float64(projected) / Float64(spent)
            var work = forward.work + backward.work
            if Float64(work) * (ratio * ratio - 1.0) > Float64(budget - work):
                return HALTED


def wavefront_score(
    first: List[UInt8], second: List[UInt8], penalties: Penalties, give_up: Bool = True
) -> Optional[Int]:
    """The optimal global score of two encoded sequences, or None once a full sweep would be cheaper.

    With `give_up` off, the search runs to the end whatever it costs.
    """
    var letters = len(first) + len(second)
    if len(first) == 0 or len(second) == 0:
        return penalties.score(gapped_cost[1](penalties, letters), letters)
    var forward = Wavefront[1](Span(first), Span(second), penalties, FREE_START, False, False)
    var backward = Wavefront[1](Span(first), Span(second), penalties, FREE_START, False, True)
    var best = Meeting.none()
    if bidirectional[1, False](forward, backward, best, give_up, Int.MAX) != MET:
        return None
    return penalties.score(best.cost, letters)


def wavefront_distance[
    pieces: Int
](
    first: Span[UInt8, _],
    second: Span[UInt8, _],
    penalties: Penalties,
    ceiling: Int,
    ends_free: EndsFree = EndsFree(),
) -> Int:
    """The least cost of a global alignment of two encoded sequences, the letters `ends_free` allows
    left unaligned for nothing, or -1 when every one costs more than `ceiling`: the two searches
    keeping only their rings, no fronts for a traceback."""
    if len(first) == 0 or len(second) == 0:
        var cost = gapped_cost[pieces](penalties, unpaid_letters(len(first), len(second), ends_free))
        return cost if cost <= ceiling else -1
    var forward = Wavefront[pieces](
        first, second, penalties, FREE_START, False, False, ends_free.first_begin, ends_free.second_begin
    )
    var backward = Wavefront[pieces](
        first, second, penalties, FREE_START, False, True, ends_free.first_end, ends_free.second_end
    )
    var best = Meeting.none()
    if bidirectional[pieces, False](forward, backward, best, False, Int.MAX, ceiling) != MET:
        return -1
    return best.cost


@inline(.always)
def unpaid_letters(columns: Int, rows: Int, ends_free: EndsFree) -> Int:
    """The letters a pair with one side empty must pay a gap for: the other side's, less what its two
    ends may leave unaligned."""
    if rows == 0:
        return max(columns - ends_free.first_begin - ends_free.first_end, 0)
    return max(rows - ends_free.second_begin - ends_free.second_end, 0)


@inline(.always)
def gapped_cost[pieces: Int](penalties: Penalties, letters: Int, continued: Int = -1) -> Int:
    """The cost of aligning `letters` letters against nothing: one gap of the cheaper piece, or none
    for no letters; or, with `continued` a piece, the extensions alone of the gap the letters continue
    if that is cheaper."""
    if letters == 0:
        return 0
    var cost = penalties.opening + penalties.extension * letters
    comptime if pieces == 2:
        cost = min(cost, penalties.opening2 + penalties.extension2 * letters)
    if continued >= 0:
        cost = min(cost, penalties.extension_of(continued) * letters)
    return cost


def trace(
    history: History,
    penalties: Penalties,
    layer: Int,
    cost: Int,
    diagonal: Int,
    column: Int,
    mut moves: List[UInt8],
):
    """An optimal path from the origin to a cell in `layer` that the front of `cost` in that layer
    reaches on `diagonal`, as moves right to left, through every cost's kept fronts.

    The cell need not be the front's furthest, as where the two searches met: every cell short of a
    front on its diagonal costs no more than the front, edit costs never falling along a diagonal. So
    the walk holds a layer, a cost and a cell its front of that cost reaches, and each step follows
    the source that won the front there, which reaches at least as far as any other. In the
    alignment layer that is a substitution from the front `x` back, or a gap layer at the same cost,
    entered after the matches back to its column when the cell lies past it. A gap cell came by a
    letter of its sequence from the diagonal beside it, opened from the alignment front `o + e` back
    or extended from its own layer `e` back, at its piece's `o` and `e`. A cell on the first row or
    column has one path left, a gap along it.
    """
    var x = penalties.mismatch
    var spent = cost
    var current = layer
    var k = diagonal
    var c = column
    while True:
        if current == ALIGNED:
            var row = c - k
            if c == 0:
                for _ in range(row):
                    moves.append(UInt8(SECOND_GAP))
                return
            if row == 0:
                for _ in range(c):
                    moves.append(UInt8(FIRST_GAP))
                return
            if spent == 0:
                # Matches back to where the diagonal starts on the first row or column, whose free
                # letters the walk then takes along that edge.
                var start = max(k, 0)
                for _ in range(c - start):
                    moves.append(UInt8(ALIGNED))
                c = start
                continue
            var source = Int(history.flag(spent, k) & ENTRY_MASK)
            # The column the winning source brought the front to, before its matches.
            var entry: Int
            if source == ALIGNED:
                entry = history.column(spent - x, k) + 1
            else:
                # A gap layer's column is where its run opened, and one letter on per extension.
                var o = penalties.opening_of(piece_of(source))
                var e = penalties.extension_of(piece_of(source))
                var opened = opened_bit(source)
                var step = 1 if along_first(source) else -1
                var gap_cost = spent
                var gap_diagonal = k
                var letters = 0
                # A run that reaches cost zero continues a gap the origin was inside, at column zero.
                while gap_cost > 0 and history.flag(gap_cost, gap_diagonal) & opened == 0:
                    gap_cost -= e
                    gap_diagonal -= step
                    letters += 1
                if along_first(source):
                    entry = (
                        letters if gap_cost == 0 else history.column(gap_cost - o - e, gap_diagonal - 1) + letters + 1
                    )
                else:
                    entry = 0 if gap_cost == 0 else history.column(gap_cost - o - e, gap_diagonal + 1)
            if entry <= c:
                for _ in range(c - entry):
                    moves.append(UInt8(ALIGNED))
                c = entry
            if source == ALIGNED:
                moves.append(UInt8(ALIGNED))
                c -= 1
                spent -= x
            else:
                current = source
        else:
            var piece = piece_of(current)
            var opened = history.flag(spent, k) & opened_bit(current) != 0
            if along_first(current):
                moves.append(UInt8(FIRST_GAP))
                c -= 1
                k -= 1
                if c == 0:
                    for _ in range(-k):
                        moves.append(UInt8(SECOND_GAP))
                    return
            else:
                moves.append(UInt8(SECOND_GAP))
                k += 1
                if c - k == 0:
                    for _ in range(c):
                        moves.append(UInt8(FIRST_GAP))
                    return
            if opened:
                spent -= penalties.opening_of(piece) + penalties.extension_of(piece)
                current = ALIGNED
            else:
                spent -= penalties.extension_of(piece)


def solve[
    pieces: Int
](
    first: Span[UInt8, _],
    second: Span[UInt8, _],
    penalties: Penalties,
    start: Int,
    finish: Int,
    limit: Int,
    mut moves: List[UInt8],
    keep: Bool = True,
    ceiling: Int = Int.MAX,
    ends_free: EndsFree = EndsFree(),
) -> Int:
    """Appends an optimal path's moves right to left and returns its cost, from an origin as `start`
    allows and to a corner the backward search's origin `finish` allows (see `FREE_START`), the
    letters `ends_free` allows left unaligned for nothing; or, when every path costs more than
    `ceiling`, appends nothing and returns -1.

    With `keep`, both searches keep every cost's fronts while they stay within `limit` entries, and
    the path is traced from where they met: back to the origin through the forward fronts, and to the
    corner through the backward ones. A pair too large is split instead where an optimal path
    crosses, which the two searches find keeping only their rings, as BiWFA does; a crossing inside a
    gap leaves the piece before it to end in that gap and the piece after it to begin there, the
    opening paid once, its piece's. A piece is about a quarter of the pair, its two searches about half of the
    diagonals the pair's search grew from its end, so one that cannot fit skips keeping at once.
    """
    var columns = len(first)
    var rows = len(second)
    if columns == 0 or rows == 0:
        # An origin inside a gap along the letters left continues it without a second opening.
        var continued = -1
        if start != FREE_START and start < OPENING and along_first(start) == (rows == 0):
            continued = piece_of(start)
        var cost = gapped_cost[pieces](penalties, unpaid_letters(columns, rows, ends_free), continued)
        if cost > ceiling:
            return -1
        for _ in range(columns):
            moves.append(UInt8(FIRST_GAP))
        for _ in range(rows):
            moves.append(UInt8(SECOND_GAP))
        return cost
    var forward = Wavefront[pieces](
        first, second, penalties, start, keep, False, ends_free.first_begin, ends_free.second_begin
    )
    var backward = Wavefront[pieces](
        first, second, penalties, finish, keep, True, ends_free.first_end, ends_free.second_end
    )
    var best = Meeting.none()
    if keep:
        var status = bidirectional[pieces, True](forward, backward, best, False, limit, ceiling)
        if status == OVER:
            return -1
        if status == MET:
            # The backward walk runs from the meeting to the corner, left to right as the forward path goes.
            var behind = List[UInt8](capacity=columns + rows)
            trace(
                backward.history,
                penalties,
                best.layer,
                best.backward_cost,
                columns - rows - best.diagonal,
                columns - best.column,
                behind,
            )
            for index in range(len(behind) - 1, -1, -1):
                moves.append(behind[index])
            trace(forward.history, penalties, best.layer, best.forward_cost, best.diagonal, best.column, moves)
            return best.cost
        # Too large to keep: the searches go on from where they stopped keeping only their rings.
        forward.history = History()
        backward.history = History()
    if bidirectional[pieces, False](forward, backward, best, False, Int.MAX, ceiling) == OVER:
        return -1
    var column = best.column
    var row = column - best.diagonal
    # A crossing at either end splits nothing: such a pair costs too little for its fronts not to fit.
    if (column == 0 and row == 0) or (column == columns and row == rows):
        return solve[pieces](first, second, penalties, start, finish, Int.MAX, moves, True, Int.MAX, ends_free)
    # A crossing inside a gap: the piece after begins in it, and the piece before must end opening it.
    var after_start = best.layer
    var before_finish = FREE_START if best.layer == ALIGNED else OPENING + best.layer
    _ = solve[pieces](
        first[column:],
        second[row:],
        penalties,
        after_start,
        finish,
        limit,
        moves,
        backward.work // 2 <= limit,
        Int.MAX,
        ends_free.trailing(),
    )
    _ = solve[pieces](
        first[:column],
        second[:row],
        penalties,
        start,
        before_finish,
        limit,
        moves,
        forward.work // 2 <= limit,
        Int.MAX,
        ends_free.leading(),
    )
    return best.cost


def wavefront_align(
    first: List[UInt8], second: List[UInt8], penalties: Penalties, alphabet: String, limit: Int = HISTORY_LIMIT
) -> Tuple[Int, String, String]:
    """The optimal global score of two encoded sequences and the gapped rows of an alignment that earns
    it, keeping at most `limit` entries of fronts at once (see `solve`)."""
    var letters = len(first) + len(second)
    var moves = List[UInt8](capacity=letters)
    var cost = solve[1](Span(first), Span(second), penalties, FREE_START, FREE_START, limit, moves)
    comptime GAP = UInt8(ord("-"))
    var symbols = alphabet.as_bytes()
    var top = List[UInt8](capacity=len(moves))
    var bottom = List[UInt8](capacity=len(moves))
    var column = 0
    var row = 0
    for index in range(len(moves) - 1, -1, -1):
        var move = Int(moves[index])
        if move == ALIGNED:
            top.append(symbols[Int(first[column])])
            bottom.append(symbols[Int(second[row])])
            column += 1
            row += 1
        elif move == FIRST_GAP:
            top.append(symbols[Int(first[column])])
            bottom.append(GAP)
            column += 1
        else:
            top.append(GAP)
            bottom.append(symbols[Int(second[row])])
            row += 1
    return (penalties.score(cost, letters), String(unsafe_from_utf8=top), String(unsafe_from_utf8=bottom))


@fieldwise_init
struct AffineCigar(Copyable, Movable, Writable):
    """The least cost of a global alignment of two sequences under gap-affine costs and an optimal
    alignment as a CIGAR string, written as `EditCigar`'s is."""

    var cost: Int
    var cigar: String


def affine_penalties(mismatch: Int, opening: Int, extension: Int) raises AlignmentError -> Penalties:
    """The wavefront's costs for gap-affine costs as WFA counts them, divided by their common factor."""
    if mismatch <= 0 or extension <= 0 or opening < 0:
        raise AlignmentError(
            ErrorKind.INVALID_SCORING,
            String("costs ", mismatch, ", ", opening, ", ", extension, ": a mismatch and an extension must cost"),
        )
    var scale = gcd(gcd(mismatch, extension), opening)
    return Penalties(mismatch // scale, opening // scale, extension // scale, scale, 0, 0, 0)


def affine2p_penalties(
    mismatch: Int, opening1: Int, extension1: Int, opening2: Int, extension2: Int
) raises AlignmentError -> Penalties:
    """The wavefront's costs for two-piece gap-affine costs as WFA counts them, divided by their common
    factor."""
    if mismatch <= 0 or extension1 <= 0 or extension2 <= 0 or opening1 < 0 or opening2 < 0:
        raise AlignmentError(
            ErrorKind.INVALID_SCORING,
            String(
                "costs ",
                mismatch,
                ", ",
                opening1,
                ", ",
                extension1,
                ", ",
                opening2,
                ", ",
                extension2,
                ": a mismatch and an extension must cost",
            ),
        )
    var scale = gcd(gcd(gcd(mismatch, extension1), gcd(opening1, extension2)), opening2)
    return Penalties(
        mismatch // scale, opening1 // scale, extension1 // scale, scale, 0, opening2 // scale, extension2 // scale
    )


def affine_distance(
    first: String,
    second: String,
    mismatch: Int,
    opening: Int,
    extension: Int,
    *,
    ends_free: EndsFree = EndsFree(),
) raises AlignmentError -> Int:
    """The least cost of a global alignment of two sequences under gap-affine costs as WFA counts them,
    a substitution `mismatch` and a gap of `k` letters `opening + k extension`, without the alignment;
    the letters `ends_free` allows at either end left unaligned for nothing.

    `affine_cigar`'s cost, by the same two searches keeping only their last few costs' fronts: about
    the same search, none of the traceback, and a few rows of memory however long the pair.
    """
    var penalties = affine_penalties(mismatch, opening, extension)
    return wavefront_distance[1](first.as_bytes(), second.as_bytes(), penalties, Int.MAX, ends_free) * penalties.scale


def affine_distance(
    first: String,
    second: String,
    mismatch: Int,
    opening: Int,
    extension: Int,
    *,
    max_cost: Int,
    ends_free: EndsFree = EndsFree(),
) raises AlignmentError -> Optional[Int]:
    """`affine_distance`, or None when it would pass `max_cost`: the searches stop once each has grown
    to about half of `max_cost` without the two meeting within it, so a pair far over costs a fraction
    of its full search."""
    return distance_within[1](first, second, affine_penalties(mismatch, opening, extension), max_cost, ends_free)


def affine_cigar(
    first: String,
    second: String,
    mismatch: Int,
    opening: Int,
    extension: Int,
    extended: Bool = True,
    *,
    ends_free: EndsFree = EndsFree(),
) raises AlignmentError -> AffineCigar:
    """The least cost of a global alignment of two sequences under gap-affine costs as WFA counts them,
    a substitution `mismatch` and a gap of `k` letters `opening + k extension`, and an optimal
    alignment as a CIGAR string, `=` and `X`, or with `extended` false `M` for both (see `EditCigar`).
    The letters `ends_free` allows at either end may go unaligned for nothing, as `D` or `I` runs.

    Every byte is a symbol matching only itself, so DNA in either case, or any other text, needs no
    alphabet. The two-ended wavefront finds it (see the module): its work grows with the square of
    the cost rather than with the matrix, and its memory stays bounded.

    At costs (2, 0, 1) the cost is the indel distance, but a substitution costs there what a deletion
    and an insertion do, and the CIGAR may write one as `X` where an indel alignment has `1D1I`.
    """
    var found = cigar_within[1](
        first, second, affine_penalties(mismatch, opening, extension), extended, Int.MAX, ends_free
    )
    return found.take()


def affine_cigar(
    first: String,
    second: String,
    mismatch: Int,
    opening: Int,
    extension: Int,
    extended: Bool = True,
    *,
    max_cost: Int,
    ends_free: EndsFree = EndsFree(),
) raises AlignmentError -> Optional[AffineCigar]:
    """`affine_cigar`, or None when the cost would pass `max_cost`, found as `affine_distance` finds
    that, with no fronts traced."""
    var penalties = affine_penalties(mismatch, opening, extension)
    if max_cost < 0:
        return None
    return cigar_within[1](first, second, penalties, extended, max_cost // penalties.scale, ends_free)


def affine2p_distance(
    first: String,
    second: String,
    mismatch: Int,
    opening1: Int,
    extension1: Int,
    opening2: Int,
    extension2: Int,
    *,
    ends_free: EndsFree = EndsFree(),
) raises AlignmentError -> Int:
    """The least cost of a global alignment of two sequences under two-piece gap-affine costs as WFA
    counts them, a substitution `mismatch` and a gap of `k` letters the least of `opening1 + k
    extension1` and `opening2 + k extension2`, without the alignment; the letters `ends_free` allows
    at either end left unaligned for nothing.

    Usually one piece opens cheaply and extends dearly and the other the reverse, so a long gap costs
    less than one affine cost would charge it. The same two searches as `affine_distance`'s, each
    cost with two more fronts.
    """
    var penalties = affine2p_penalties(mismatch, opening1, extension1, opening2, extension2)
    return wavefront_distance[2](first.as_bytes(), second.as_bytes(), penalties, Int.MAX, ends_free) * penalties.scale


def affine2p_distance(
    first: String,
    second: String,
    mismatch: Int,
    opening1: Int,
    extension1: Int,
    opening2: Int,
    extension2: Int,
    *,
    max_cost: Int,
    ends_free: EndsFree = EndsFree(),
) raises AlignmentError -> Optional[Int]:
    """`affine2p_distance`, or None when it would pass `max_cost`, as `affine_distance` caps it."""
    var penalties = affine2p_penalties(mismatch, opening1, extension1, opening2, extension2)
    return distance_within[2](first, second, penalties, max_cost, ends_free)


def affine2p_cigar(
    first: String,
    second: String,
    mismatch: Int,
    opening1: Int,
    extension1: Int,
    opening2: Int,
    extension2: Int,
    extended: Bool = True,
    *,
    ends_free: EndsFree = EndsFree(),
) raises AlignmentError -> AffineCigar:
    """`affine2p_distance` and an optimal alignment as a CIGAR string, written as `affine_cigar`'s is."""
    var penalties = affine2p_penalties(mismatch, opening1, extension1, opening2, extension2)
    var found = cigar_within[2](first, second, penalties, extended, Int.MAX, ends_free)
    return found.take()


def affine2p_cigar(
    first: String,
    second: String,
    mismatch: Int,
    opening1: Int,
    extension1: Int,
    opening2: Int,
    extension2: Int,
    extended: Bool = True,
    *,
    max_cost: Int,
    ends_free: EndsFree = EndsFree(),
) raises AlignmentError -> Optional[AffineCigar]:
    """`affine2p_cigar`, or None when the cost would pass `max_cost`, as `affine_cigar` caps it."""
    var penalties = affine2p_penalties(mismatch, opening1, extension1, opening2, extension2)
    if max_cost < 0:
        return None
    return cigar_within[2](first, second, penalties, extended, max_cost // penalties.scale, ends_free)


def distance_within[
    pieces: Int
](first: String, second: String, penalties: Penalties, max_cost: Int, ends_free: EndsFree) -> Optional[Int]:
    """The least cost in the caller's units, or None when it would pass `max_cost`."""
    if max_cost < 0:
        return None
    var cost = wavefront_distance[pieces](
        first.as_bytes(), second.as_bytes(), penalties, max_cost // penalties.scale, ends_free
    )
    if cost < 0:
        return None
    return cost * penalties.scale


def cigar_within[
    pieces: Int
](
    first: String,
    second: String,
    penalties: Penalties,
    extended: Bool,
    ceiling: Int,
    ends_free: EndsFree = EndsFree(),
) -> Optional[AffineCigar]:
    """An optimal alignment's cost and CIGAR, or None when its cost, in `penalties`' units, would pass
    `ceiling`."""
    var columns = first.byte_length()
    var rows = second.byte_length()
    # UTF-8 never holds the sentinels' bytes, so the text's own bytes serve as codes.
    var moves = List[UInt8](capacity=columns + rows)
    var cost = solve[pieces](
        first.as_bytes(),
        second.as_bytes(),
        penalties,
        FREE_START,
        FREE_START,
        HISTORY_LIMIT,
        moves,
        True,
        ceiling,
        ends_free,
    )
    if cost < 0:
        return None
    # The CIGAR's room is bounded by the edits: every gapped letter, and a substitution per mismatch cost.
    var gapped = 0
    for move in moves:
        if move != UInt8(ALIGNED):
            gapped += 1
    var path = EditPath(moves^, List[UInt8](), columns, rows, gapped + cost // penalties.mismatch)
    return AffineCigar(cost * penalties.scale, cigar_string(first, second, path, extended))
