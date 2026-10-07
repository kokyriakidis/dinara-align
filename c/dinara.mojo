"""
The C API of dinara-align: the unit-cost edit distance, and the least gap-affine cost as WFA counts it,
with one gap piece or two, alone, capped, banded or with an optimal alignment as a CIGAR, and the best
extension from one end, for C, C++ and any language with a C foreign-function interface. `dinara.h` declares these, with a C++ wrapper.

    pixi run build-c [target-cpu]   # build/c: libdinara, its runtime libraries and the header

Each function takes the two sequences as bytes and lengths, holds nothing between calls and so may be
called from many threads at once, and reports failure as a negative result (`DINARA_*` in the header).
"""

from std.ffi import external_call

from dinara_align import (
    AffineCigar,
    AffineExtension,
    AlignmentError,
    Anchor,
    Band,
    EndsFree,
    ErrorKind,
    affine2p_cigar,
    affine2p_distance,
    affine2p_extension,
    affine_cigar,
    affine_distance,
    affine_extension,
    edit_cigar,
    edit_distance,
)

comptime UNSUPPORTED_SYMBOLS = -1
"""More than four symbols past `ACGT` between the two sequences."""
comptime OUT_OF_MEMORY = -2
"""The CIGAR's memory could not be allocated."""
comptime INVALID_COSTS = -3
"""Gap-affine costs no alignment can be searched by: a free mismatch or extension, or a negative cost."""
comptime ABOVE_MAX = -4
"""Every alignment costs more than the `max_cost` asked for."""
comptime OUTSIDE_BAND = -5
"""No alignment stays inside the band asked for."""


def sequence(bytes: ImmPointer[UInt8, MutAnyOrigin], length: Int) -> String:
    """A sequence from C bytes; an empty one may come as a null pointer, never read."""
    if length <= 0:
        return String()
    return String(StringSlice(unsafe_from_utf8=Span(unsafe_ptr=bytes, length=length)))


@export("dinara_edit_distance")
def dinara_edit_distance(
    first: ImmPointer[UInt8, MutAnyOrigin],
    first_length: Int,
    second: ImmPointer[UInt8, MutAnyOrigin],
    second_length: Int,
) abi("C") -> Int:
    """The global edit distance between two sequences."""
    try:
        return edit_distance(sequence(first, first_length), sequence(second, second_length))
    except:
        return UNSUPPORTED_SYMBOLS


@export("dinara_edit_cigar")
def dinara_edit_cigar(
    first: ImmPointer[UInt8, MutAnyOrigin],
    first_length: Int,
    second: ImmPointer[UInt8, MutAnyOrigin],
    second_length: Int,
    extended: Int32,
    cigar: MutPointer[MutPointer[UInt8, MutAnyOrigin], MutAnyOrigin],
    cigar_length: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """The global edit distance and an optimal alignment's CIGAR, `=` and `X` when `extended` is nonzero,
    else `M`. The CIGAR, NUL-terminated, is allocated with C's `malloc`, as its length is known only once
    the alignment is, and goes to the caller, who frees it with `dinara_free`."""
    try:
        var aligned = edit_cigar(sequence(first, first_length), sequence(second, second_length), extended != 0)
        hand_over(aligned.cigar, cigar, cigar_length)
        return aligned.distance
    except:
        return UNSUPPORTED_SYMBOLS


@export("dinara_affine_distance")
def dinara_affine_distance(
    first: ImmPointer[UInt8, MutAnyOrigin],
    first_length: Int,
    second: ImmPointer[UInt8, MutAnyOrigin],
    second_length: Int,
    mismatch: Int,
    opening: Int,
    extension: Int,
    max_cost: Int,
    first_begin_free: Int,
    first_end_free: Int,
    second_begin_free: Int,
    second_end_free: Int,
    band_low: Int,
    band_high: Int,
) abi("C") -> Int:
    """The least global cost under gap-affine costs, with no alignment, or `ABOVE_MAX` when it passes a
    `max_cost` of zero or more; a negative `max_cost` caps nothing. The four `free` counts are the
    letters at each end of each sequence left unaligned for nothing, all zero for a global alignment.
    Every move stays on the diagonals `band_low ..= band_high` (see `Band`), the C integer limits for
    none; `OUTSIDE_BAND` when no alignment does, or `ABOVE_MAX` under a cap."""
    var pair = Pair(first, first_length, second, second_length, first_begin_free, first_end_free)
    return distance(
        pair, mismatch, opening, extension, -1, 0, max_cost, second_begin_free, second_end_free, band_low, band_high
    )


@export("dinara_affine_cigar")
def dinara_affine_cigar(
    first: ImmPointer[UInt8, MutAnyOrigin],
    first_length: Int,
    second: ImmPointer[UInt8, MutAnyOrigin],
    second_length: Int,
    mismatch: Int,
    opening: Int,
    extension: Int,
    max_cost: Int,
    first_begin_free: Int,
    first_end_free: Int,
    second_begin_free: Int,
    second_end_free: Int,
    band_low: Int,
    band_high: Int,
    extended: Int32,
    cigar: MutPointer[MutPointer[UInt8, MutAnyOrigin], MutAnyOrigin],
    cigar_length: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """The least global cost under gap-affine costs, a substitution `mismatch` and a gap of `k` letters
    `opening + k extension`, and an optimal alignment's CIGAR, handed over as `dinara_edit_cigar`'s; or
    `ABOVE_MAX`, and no CIGAR, when the cost passes a `max_cost` of zero or more. Every byte is a symbol
    matching only itself, save the two UTF-8 never holds, which mark the ends; the free counts and the
    band as for `dinara_affine_distance`, the free letters `D` and `I` runs."""
    var pair = Pair(first, first_length, second, second_length, first_begin_free, first_end_free)
    return aligned(
        pair,
        mismatch,
        opening,
        extension,
        -1,
        0,
        max_cost,
        second_begin_free,
        second_end_free,
        band_low,
        band_high,
        extended,
        cigar,
        cigar_length,
    )


@export("dinara_affine2p_distance")
def dinara_affine2p_distance(
    first: ImmPointer[UInt8, MutAnyOrigin],
    first_length: Int,
    second: ImmPointer[UInt8, MutAnyOrigin],
    second_length: Int,
    mismatch: Int,
    opening1: Int,
    extension1: Int,
    opening2: Int,
    extension2: Int,
    max_cost: Int,
    first_begin_free: Int,
    first_end_free: Int,
    second_begin_free: Int,
    second_end_free: Int,
    band_low: Int,
    band_high: Int,
) abi("C") -> Int:
    """`dinara_affine_distance` under two-piece gap-affine costs, a gap of `k` letters the less of
    `opening1 + k extension1` and `opening2 + k extension2`."""
    if opening2 < 0:
        return INVALID_COSTS
    var pair = Pair(first, first_length, second, second_length, first_begin_free, first_end_free)
    return distance(
        pair,
        mismatch,
        opening1,
        extension1,
        opening2,
        extension2,
        max_cost,
        second_begin_free,
        second_end_free,
        band_low,
        band_high,
    )


@export("dinara_affine2p_cigar")
def dinara_affine2p_cigar(
    first: ImmPointer[UInt8, MutAnyOrigin],
    first_length: Int,
    second: ImmPointer[UInt8, MutAnyOrigin],
    second_length: Int,
    mismatch: Int,
    opening1: Int,
    extension1: Int,
    opening2: Int,
    extension2: Int,
    max_cost: Int,
    first_begin_free: Int,
    first_end_free: Int,
    second_begin_free: Int,
    second_end_free: Int,
    band_low: Int,
    band_high: Int,
    extended: Int32,
    cigar: MutPointer[MutPointer[UInt8, MutAnyOrigin], MutAnyOrigin],
    cigar_length: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """`dinara_affine_cigar` under two-piece gap-affine costs, as `dinara_affine2p_distance` counts them."""
    if opening2 < 0:
        return INVALID_COSTS
    var pair = Pair(first, first_length, second, second_length, first_begin_free, first_end_free)
    return aligned(
        pair,
        mismatch,
        opening1,
        extension1,
        opening2,
        extension2,
        max_cost,
        second_begin_free,
        second_end_free,
        band_low,
        band_high,
        extended,
        cigar,
        cigar_length,
    )


@export("dinara_affine_extension")
def dinara_affine_extension(
    first: ImmPointer[UInt8, MutAnyOrigin],
    first_length: Int,
    second: ImmPointer[UInt8, MutAnyOrigin],
    second_length: Int,
    match_score: Int,
    mismatch: Int,
    opening: Int,
    extension: Int,
    at_end: Int32,
    band_low: Int,
    band_high: Int,
    extended: Int32,
    cigar: MutPointer[MutPointer[UInt8, MutAnyOrigin], MutAnyOrigin],
    cigar_length: MutPointer[Int, MutAnyOrigin],
    first_covered: MutPointer[Int, MutAnyOrigin],
    second_covered: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """The best score of an alignment fixed at both sequences' starts, or with `at_end` nonzero their
    ends, and free to stop anywhere (see `affine_extension`), or a negative code: a match earns
    `match_score`, the costs as for `dinara_affine_cigar`. The letters of each sequence it covers from
    that end go to `first_covered` and `second_covered`, and its CIGAR over them is handed over as
    `dinara_edit_cigar`'s. The band counts diagonals from the anchor."""
    var pair = Pair(first, first_length, second, second_length, 0, 0)
    return extended_from(
        pair,
        match_score,
        mismatch,
        opening,
        extension,
        -1,
        0,
        at_end,
        band_low,
        band_high,
        extended,
        cigar,
        cigar_length,
        first_covered,
        second_covered,
    )


@export("dinara_affine2p_extension")
def dinara_affine2p_extension(
    first: ImmPointer[UInt8, MutAnyOrigin],
    first_length: Int,
    second: ImmPointer[UInt8, MutAnyOrigin],
    second_length: Int,
    match_score: Int,
    mismatch: Int,
    opening1: Int,
    extension1: Int,
    opening2: Int,
    extension2: Int,
    at_end: Int32,
    band_low: Int,
    band_high: Int,
    extended: Int32,
    cigar: MutPointer[MutPointer[UInt8, MutAnyOrigin], MutAnyOrigin],
    cigar_length: MutPointer[Int, MutAnyOrigin],
    first_covered: MutPointer[Int, MutAnyOrigin],
    second_covered: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """`dinara_affine_extension` under two-piece gap-affine costs."""
    if opening2 < 0:
        return INVALID_COSTS
    var pair = Pair(first, first_length, second, second_length, 0, 0)
    return extended_from(
        pair,
        match_score,
        mismatch,
        opening1,
        extension1,
        opening2,
        extension2,
        at_end,
        band_low,
        band_high,
        extended,
        cigar,
        cigar_length,
        first_covered,
        second_covered,
    )


struct Pair(Movable):
    """Two sequences C handed over, as text, whether neither holds a byte the wavefront reserves, and
    the first sequence's free letters at either end."""

    var plain: Bool
    var first: String
    var second: String
    var first_begin_free: Int
    var first_end_free: Int

    def __init__(
        out self,
        first: ImmPointer[UInt8, MutAnyOrigin],
        first_length: Int,
        second: ImmPointer[UInt8, MutAnyOrigin],
        second_length: Int,
        first_begin_free: Int,
        first_end_free: Int,
    ):
        self.plain = plain_bytes(first, first_length) and plain_bytes(second, second_length)
        self.first = sequence(first, first_length) if self.plain else String()
        self.second = sequence(second, second_length) if self.plain else String()
        self.first_begin_free = first_begin_free
        self.first_end_free = first_end_free


def band_of(low: Int, high: Int) -> Band:
    """A band from C's two edges, the integer limits standing for none, kept clear of overflow."""
    comptime EDGE = 1 << 60
    return Band(max(low, -EDGE), min(high, EDGE))


def failure(error: AlignmentError) -> Int:
    """The code for what the library raised: a band no alignment fits, or costs it cannot search by."""
    return OUTSIDE_BAND if error.kind == ErrorKind.INVALID_ARGUMENT else INVALID_COSTS


def distance(
    pair: Pair,
    mismatch: Int,
    opening: Int,
    extension: Int,
    opening2: Int,
    extension2: Int,
    max_cost: Int,
    second_begin_free: Int,
    second_end_free: Int,
    band_low: Int,
    band_high: Int,
) -> Int:
    """The cost for the exports above: one gap piece, or a second with `opening2` not negative."""
    if not pair.plain:
        return UNSUPPORTED_SYMBOLS
    var ends = EndsFree(pair.first_begin_free, pair.first_end_free, second_begin_free, second_end_free)
    var band = band_of(band_low, band_high)
    var a = pair.first
    var b = pair.second
    try:
        if max_cost < 0:
            if opening2 >= 0:
                return affine2p_distance(
                    a, b, mismatch, opening, extension, opening2, extension2, ends_free=ends, band=band
                )
            return affine_distance(a, b, mismatch, opening, extension, ends_free=ends, band=band)
        var found: Optional[Int]
        if opening2 >= 0:
            found = affine2p_distance(
                a, b, mismatch, opening, extension, opening2, extension2, max_cost=max_cost, ends_free=ends, band=band
            )
        else:
            found = affine_distance(a, b, mismatch, opening, extension, max_cost=max_cost, ends_free=ends, band=band)
        return found.value() if found else ABOVE_MAX
    except error:
        return failure(error)


def aligned(
    pair: Pair,
    mismatch: Int,
    opening: Int,
    extension: Int,
    opening2: Int,
    extension2: Int,
    max_cost: Int,
    second_begin_free: Int,
    second_end_free: Int,
    band_low: Int,
    band_high: Int,
    extended: Int32,
    cigar: MutPointer[MutPointer[UInt8, MutAnyOrigin], MutAnyOrigin],
    cigar_length: MutPointer[Int, MutAnyOrigin],
) -> Int:
    """The cost and CIGAR for the exports above, as `distance` takes the costs."""
    if not pair.plain:
        return UNSUPPORTED_SYMBOLS
    var ends = EndsFree(pair.first_begin_free, pair.first_end_free, second_begin_free, second_end_free)
    var band = band_of(band_low, band_high)
    var a = pair.first
    var b = pair.second
    try:
        var found: Optional[AffineCigar]
        if max_cost < 0:
            if opening2 >= 0:
                found = affine2p_cigar(
                    a, b, mismatch, opening, extension, opening2, extension2, extended != 0, ends_free=ends, band=band
                )
            else:
                found = affine_cigar(a, b, mismatch, opening, extension, extended != 0, ends_free=ends, band=band)
        elif opening2 >= 0:
            found = affine2p_cigar(
                a,
                b,
                mismatch,
                opening,
                extension,
                opening2,
                extension2,
                extended != 0,
                max_cost=max_cost,
                ends_free=ends,
                band=band,
            )
        else:
            found = affine_cigar(
                a, b, mismatch, opening, extension, extended != 0, max_cost=max_cost, ends_free=ends, band=band
            )
        if not found:
            return ABOVE_MAX
        hand_over(found.value().cigar, cigar, cigar_length)
        return found.value().cost
    except error:
        return failure(error)


def extended_from(
    pair: Pair,
    match_score: Int,
    mismatch: Int,
    opening: Int,
    extension: Int,
    opening2: Int,
    extension2: Int,
    at_end: Int32,
    band_low: Int,
    band_high: Int,
    extended: Int32,
    cigar: MutPointer[MutPointer[UInt8, MutAnyOrigin], MutAnyOrigin],
    cigar_length: MutPointer[Int, MutAnyOrigin],
    first_covered: MutPointer[Int, MutAnyOrigin],
    second_covered: MutPointer[Int, MutAnyOrigin],
) -> Int:
    """The extension for the exports above, as `distance` takes the costs."""
    if not pair.plain:
        return UNSUPPORTED_SYMBOLS
    var band = band_of(band_low, band_high)
    var anchor = Anchor.END if at_end != 0 else Anchor.START
    try:
        var found: AffineExtension
        if opening2 >= 0:
            found = affine2p_extension(
                pair.first,
                pair.second,
                match_score,
                mismatch,
                opening,
                extension,
                opening2,
                extension2,
                extended != 0,
                anchor=anchor,
                band=band,
            )
        else:
            found = affine_extension(
                pair.first,
                pair.second,
                match_score,
                mismatch,
                opening,
                extension,
                extended != 0,
                anchor=anchor,
                band=band,
            )
        hand_over(found.cigar, cigar, cigar_length)
        first_covered[] = found.first_length
        second_covered[] = found.second_length
        return found.score
    except error:
        return failure(error)


def plain_bytes(bytes: ImmPointer[UInt8, MutAnyOrigin], length: Int) -> Bool:
    """Whether no byte is one of the two UTF-8 never holds, which the wavefront's sentinels are."""
    for index in range(length):
        if bytes[unsafe_offset=index] >= 0xFE:
            return False
    return True


def hand_over(
    text: String,
    cigar: MutPointer[MutPointer[UInt8, MutAnyOrigin], MutAnyOrigin],
    cigar_length: MutPointer[Int, MutAnyOrigin],
):
    """`text`, NUL-terminated, in memory from C's `malloc`, as its length is known only once the alignment
    is, for the caller to free with `dinara_free`."""
    var bytes = text.as_bytes()
    var copy = external_call["malloc", MutPointer[UInt8, MutAnyOrigin]](len(bytes) + 1)
    for index in range(len(bytes)):
        copy[unsafe_offset=index] = bytes[index]
    copy[unsafe_offset=len(bytes)] = 0
    cigar[] = copy
    cigar_length[] = len(bytes)


@export("dinara_free")
def dinara_free(text: MutPointer[UInt8, MutAnyOrigin]) abi("C"):
    """Frees a CIGAR a function here allocated, with the C library's `free` that matches its `malloc`."""
    external_call["free", NoneType](text)
