# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
A substitution table as the anti-diagonal sweeps read it, a vector of letter pairs at a time.

A table of one match and one mismatch score scores a vector of pairs by one comparison. A table of up to
sixteen entries, a four-letter alphabet's, sits in one vector register and is read by byte shuffle; any
larger one, BLOSUM62's say, is gathered from memory, every lane's pair at once. So a table telling
transitions from transversions sweeps sixteen cells a step about as fast as a uniform one, to the same
scores a cell-by-cell sweep gives.
"""

from std.memory import bitcast

from .common import SubstitutionDType


comptime SHUFFLED_ENTRIES = 16
"""Entries a table holds at most to be looked up by byte shuffle: one 16-byte vector, so a four-letter
alphabet's, DNA's."""


@inline(.always)
@always_inline
def looked_up[
    value: DType, width: Int
](table: SIMD[DType.uint8, SHUFFLED_ENTRIES], index: SIMD[DType.uint8, width]) -> SIMD[value, width]:
    """Each lane's entry of `table`, `pshufb` or `tbl` a sixteen lanes, widened to `value`, signed for a
    signed `value`."""
    var out = SIMD[value, width]()
    comptime if width < SHUFFLED_ENTRIES:
        # Fewer lanes than a shuffle takes, as NEON's eight 16-bit ones: one shuffle, its first lanes kept.
        var part = table._dynamic_shuffle(SIMD[DType.uint8, SHUFFLED_ENTRIES](0).insert[offset=0](index)).slice[width]()
        comptime if value.is_signed():
            return bitcast[DType.int8, width](part).cast[value]()
        else:
            return part.cast[value]()
    comptime for chunk in range(width // SHUFFLED_ENTRIES):
        var part = table._dynamic_shuffle(index.slice[SHUFFLED_ENTRIES, offset=chunk * SHUFFLED_ENTRIES]())
        comptime if value.is_signed():
            out = out.insert[offset=chunk * SHUFFLED_ENTRIES](bitcast[DType.int8, SHUFFLED_ENTRIES](part).cast[value]())
        else:
            out = out.insert[offset=chunk * SHUFFLED_ENTRIES](part.cast[value]())
    return out


def table_extremes(substitutions: ImmSpan[Scalar[SubstitutionDType], _], alphabet_size: Int) -> Tuple[Int, Int]:
    """A square table's largest score and its least."""
    var best = Int.MIN
    var least = Int.MAX
    for cell in range(alphabet_size * alphabet_size):
        best = max(best, Int(substitutions[cell]))
        least = min(least, Int(substitutions[cell]))
    return (best, least)


def uniform_pair(substitutions: ImmSpan[Scalar[SubstitutionDType], _], alphabet_size: Int) -> Optional[Tuple[Int, Int]]:
    """The match and the mismatch score of a table of one each, if it is one; a one-letter alphabet's
    table never is, having no mismatch to read."""
    if alphabet_size < 2:
        return None
    var hit = Int(substitutions[0])
    var mismatch = Int(substitutions[1])
    for row in range(alphabet_size):
        for column in range(alphabet_size):
            if Int(substitutions[row * alphabet_size + column]) != (hit if row == column else mismatch):
                return None
    return (hit, mismatch)


def shuffled_table(substitutions: ImmSpan[Scalar[SubstitutionDType], _], alphabet_size: Int) -> SIMD[DType.uint8, 16]:
    """A table of at most `SHUFFLED_ENTRIES` entries as the bytes `looked_up` reads, row by row; the
    rest zero."""
    var table = SIMD[DType.uint8, SHUFFLED_ENTRIES](0)
    for cell in range(min(alphabet_size * alphabet_size, SHUFFLED_ENTRIES)):
        table[cell] = bitcast[DType.uint8](substitutions[cell])
    return table


struct SubstitutionLookup(Copyable, Movable):
    """A substitution table laid out for the sweeps' three ways of reading it, and what a sweep needs to
    know of it beside."""

    var cells: List[Int32]
    """Pair `(a, b)`'s score at `a * stride + b`, for a gather."""
    var stride: Int
    """The alphabet's size. A code past it, the padding past either sequence's end, reads as the last
    letter: no sweep keeps a lane that reads padding."""
    var uniform: Bool
    """Whether every match scores `reward` and every mismatch `mismatch`, so one comparison scores a pair."""
    var reward: Int32
    var mismatch: Int32
    var shuffled: SIMD[DType.uint8, SHUFFLED_ENTRIES]
    """The table as bytes for `looked_up`, when it has no more than `SHUFFLED_ENTRIES` entries."""
    var small: Bool
    var best: Int
    """The table's largest score: no substitution earns more."""

    def __init__(out self, substitutions: ImmSpan[Scalar[SubstitutionDType], _], alphabet_size: Int):
        """The lookup of a square table of `alphabet_size` letters, row by row."""
        self.stride = alphabet_size
        self.cells = List[Int32](capacity=alphabet_size * alphabet_size)
        self.best = Int(Int32.MIN)
        for cell in range(alphabet_size * alphabet_size):
            self.cells.append(Int32(substitutions[cell]))
            self.best = max(self.best, Int(substitutions[cell]))
        self.small = alphabet_size * alphabet_size <= SHUFFLED_ENTRIES
        self.shuffled = shuffled_table(substitutions, alphabet_size)
        # A one-letter alphabet has no mismatch; any value below its match keeps the comparison honest.
        var reward = Int(substitutions[0]) if alphabet_size > 0 else 0
        var mismatch = Int(substitutions[1]) if alphabet_size > 1 else reward - 1
        self.uniform = alphabet_size > 0
        for row in range(alphabet_size):
            for column in range(alphabet_size):
                if Int(substitutions[row * alphabet_size + column]) != (reward if row == column else mismatch):
                    self.uniform = False
        self.reward = Int32(reward)
        self.mismatch = Int32(mismatch)

    def lanes[width: Int](self) -> SubstitutionLanes[width]:
        """The table as a sweep of `width` lanes reads it, held in registers for the sweep's length; the
        lookup must outlive it."""
        return SubstitutionLanes[width](
            self.cells.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
            self.stride,
            self.uniform,
            self.reward,
            self.mismatch,
            self.small,
            self.shuffled,
        )


@fieldwise_init
struct SubstitutionLanes[width: Int](ImplicitlyCopyable, TrivialRegisterPassable):
    """A `SubstitutionLookup` as one sweep reads it, `width` pairs a call (see `SubstitutionLookup.lanes`)."""

    var cells: ImmPointer[Int32, ImmUntrackedOrigin]
    var stride: Int
    var uniform: Bool
    var reward: Int32
    var mismatch: Int32
    var small: Bool
    var shuffled: SIMD[DType.uint8, SHUFFLED_ENTRIES]

    @inline(.always)
    def __call__(
        self, mine: SIMD[DType.uint8, Self.width], theirs: SIMD[DType.uint8, Self.width]
    ) -> SIMD[DType.int32, Self.width]:
        """Each lane's pair's score: by comparison for a uniform table, by byte shuffle for a small one,
        else gathered from the table, a code past the alphabet, the sweeps' padding, read as its last
        letter."""
        comptime Scores = SIMD[DType.int32, Self.width]
        if self.uniform:
            return mine.eq(theirs).select(Scores(self.reward), Scores(self.mismatch))
        var last = SIMD[DType.uint8, Self.width](UInt8(self.stride - 1))
        var row = min(mine, last)
        var column = min(theirs, last)
        if self.small:
            return looked_up[DType.int32, Self.width](self.shuffled, row * UInt8(self.stride) + column)
        return self.cells.unsafe_gather[width=Self.width](
            row.cast[DType.int32]() * Int32(self.stride) + column.cast[DType.int32]()
        )
