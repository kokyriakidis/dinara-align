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
    value: DType, width: Int, signed: Bool = value.is_signed()
](table: SIMD[DType.uint8, SHUFFLED_ENTRIES], index: SIMD[DType.uint8, width]) -> SIMD[value, width]:
    """Each lane's entry of `table`, `pshufb` or `tbl` a sixteen lanes, widened to `value`: a score signed,
    as for a signed `value`, or with `signed` false a cost of up to 255, whatever `value` is."""
    var out = SIMD[value, width]()
    comptime if width < SHUFFLED_ENTRIES:
        # Fewer lanes than a shuffle takes, as NEON's eight 16-bit ones: one shuffle, its first lanes kept.
        var part = table._dynamic_shuffle(SIMD[DType.uint8, SHUFFLED_ENTRIES](0).insert[offset=0](index)).slice[width]()
        comptime if signed:
            return bitcast[DType.int8, width](part).cast[value]()
        else:
            return part.cast[value]()
    comptime for chunk in range(width // SHUFFLED_ENTRIES):
        var part = table._dynamic_shuffle(index.slice[SHUFFLED_ENTRIES, offset=chunk * SHUFFLED_ENTRIES]())
        comptime if signed:
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
    """The alphabet's size, which every code read lies within: the sweeps pad their sequences with code zero."""
    var uniform: Bool
    """Whether every match scores `reward` and every mismatch `mismatch`, so one comparison scores a pair."""
    var reward: Int
    var mismatch: Int
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
        self.reward = reward
        self.mismatch = mismatch

    @staticmethod
    def uniform_of(reward: Int, mismatch: Int) -> Self:
        """The lookup of every match scoring `reward` and every mismatch `mismatch`, over any letters, which a
        sweep reads by comparison alone (see `lanes`), there being no table to read."""
        var table: List[Scalar[SubstitutionDType]] = [Scalar[SubstitutionDType](0)]
        var lookup = Self(table, 1)
        lookup.reward = reward
        lookup.mismatch = mismatch
        lookup.best = reward
        return lookup^

    def lanes[
        width: Int, dtype: DType = DType.int32, compared: Bool = False
    ](self) -> SubstitutionLanes[width, dtype, compared]:
        """The table as a sweep of `width` lanes of `dtype` reads it, held in registers for the sweep's
        length; the lookup must outlive it. `compared` says the table is uniform, read by comparison alone."""
        return SubstitutionLanes[width, dtype, compared](
            self.cells.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
            self.stride,
            self.uniform,
            Scalar[dtype](self.reward),
            Scalar[dtype](self.mismatch),
            self.small,
            self.shuffled,
        )


@fieldwise_init
struct SubstitutionLanes[width: Int, dtype: DType = DType.int32, compared: Bool = False](
    ImplicitlyCopyable, TrivialRegisterPassable
):
    """A `SubstitutionLookup` as one sweep reads it, `width` pairs a call (see `SubstitutionLookup.lanes`)."""

    var cells: ImmPointer[Int32, ImmUntrackedOrigin]
    var stride: Int
    var uniform: Bool
    var reward: Scalar[Self.dtype]
    var mismatch: Scalar[Self.dtype]
    var small: Bool
    var shuffled: SIMD[DType.uint8, SHUFFLED_ENTRIES]

    @inline(.always)
    def __call__(
        self, mine: SIMD[DType.uint8, Self.width], theirs: SIMD[DType.uint8, Self.width]
    ) -> SIMD[Self.dtype, Self.width]:
        """Each lane's pair's score, `mine` the table's row: by comparison for a uniform table, by byte
        shuffle for a small one, else gathered from the table."""
        comptime Scores = SIMD[Self.dtype, Self.width]
        comptime if Self.compared:
            return mine.eq(theirs).select(Scores(self.reward), Scores(self.mismatch))
        if self.small:
            return looked_up[Self.dtype, Self.width](self.shuffled, mine * UInt8(self.stride) + theirs)
        return self.cells.unsafe_gather[width=Self.width](
            mine.cast[DType.int32]() * Int32(self.stride) + theirs.cast[DType.int32]()
        ).cast[Self.dtype]()
