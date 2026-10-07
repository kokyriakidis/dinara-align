"""
The C API of dinara-align: the unit-cost edit distance, and the least gap-affine cost as WFA counts it,
with one gap piece or two, alone, capped or with an optimal alignment as a CIGAR, for C, C++ and any language with a C foreign-function
interface. `dinara.h` declares these, with a C++ wrapper.

    pixi run build-c [target-cpu]   # build/c: libdinara, its runtime libraries and the header

Each function takes the two sequences as bytes and lengths, holds nothing between calls and so may be
called from many threads at once, and reports failure as a negative result (`DINARA_*` in the header).
"""

from std.ffi import external_call

from dinara_align import (
    EndsFree,
    affine2p_cigar,
    affine2p_distance,
    affine_cigar,
    affine_distance,
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
) abi("C") -> Int:
    """The least global cost under gap-affine costs, with no alignment, or `ABOVE_MAX` when it passes a
    `max_cost` of zero or more; a negative `max_cost` caps nothing. The four `free` counts are the
    letters at each end of each sequence left unaligned for nothing, all zero for a global alignment."""
    var ends = EndsFree(first_begin_free, first_end_free, second_begin_free, second_end_free)
    if not plain_bytes(first, first_length) or not plain_bytes(second, second_length):
        return UNSUPPORTED_SYMBOLS
    try:
        var a = sequence(first, first_length)
        var b = sequence(second, second_length)
        if max_cost < 0:
            return affine_distance(a, b, mismatch, opening, extension, ends_free=ends)
        var found = affine_distance(a, b, mismatch, opening, extension, max_cost=max_cost, ends_free=ends)
        return found.value() if found else ABOVE_MAX
    except:
        return INVALID_COSTS


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
    extended: Int32,
    cigar: MutPointer[MutPointer[UInt8, MutAnyOrigin], MutAnyOrigin],
    cigar_length: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """The least global cost under gap-affine costs, a substitution `mismatch` and a gap of `k` letters
    `opening + k extension`, and an optimal alignment's CIGAR, handed over as `dinara_edit_cigar`'s; or
    `ABOVE_MAX`, and no CIGAR, when the cost passes a `max_cost` of zero or more. Every byte is a symbol
    matching only itself, save the two UTF-8 never holds, which mark the ends; the free counts as for
    `dinara_affine_distance`, their letters `D` and `I` runs."""
    var ends = EndsFree(first_begin_free, first_end_free, second_begin_free, second_end_free)
    if not plain_bytes(first, first_length) or not plain_bytes(second, second_length):
        return UNSUPPORTED_SYMBOLS
    try:
        var a = sequence(first, first_length)
        var b = sequence(second, second_length)
        if max_cost < 0:
            var aligned = affine_cigar(a, b, mismatch, opening, extension, extended != 0, ends_free=ends)
            hand_over(aligned.cigar, cigar, cigar_length)
            return aligned.cost
        var found = affine_cigar(a, b, mismatch, opening, extension, extended != 0, max_cost=max_cost, ends_free=ends)
        if not found:
            return ABOVE_MAX
        hand_over(found.value().cigar, cigar, cigar_length)
        return found.value().cost
    except:
        return INVALID_COSTS


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
) abi("C") -> Int:
    """`dinara_affine_distance` under two-piece gap-affine costs, a gap of `k` letters the less of
    `opening1 + k extension1` and `opening2 + k extension2`."""
    var ends = EndsFree(first_begin_free, first_end_free, second_begin_free, second_end_free)
    if not plain_bytes(first, first_length) or not plain_bytes(second, second_length):
        return UNSUPPORTED_SYMBOLS
    try:
        var a = sequence(first, first_length)
        var b = sequence(second, second_length)
        if max_cost < 0:
            return affine2p_distance(a, b, mismatch, opening1, extension1, opening2, extension2, ends_free=ends)
        var found = affine2p_distance(
            a, b, mismatch, opening1, extension1, opening2, extension2, max_cost=max_cost, ends_free=ends
        )
        return found.value() if found else ABOVE_MAX
    except:
        return INVALID_COSTS


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
    extended: Int32,
    cigar: MutPointer[MutPointer[UInt8, MutAnyOrigin], MutAnyOrigin],
    cigar_length: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """`dinara_affine_cigar` under two-piece gap-affine costs, as `dinara_affine2p_distance` counts them."""
    var ends = EndsFree(first_begin_free, first_end_free, second_begin_free, second_end_free)
    if not plain_bytes(first, first_length) or not plain_bytes(second, second_length):
        return UNSUPPORTED_SYMBOLS
    try:
        var a = sequence(first, first_length)
        var b = sequence(second, second_length)
        if max_cost < 0:
            var aligned = affine2p_cigar(
                a, b, mismatch, opening1, extension1, opening2, extension2, extended != 0, ends_free=ends
            )
            hand_over(aligned.cigar, cigar, cigar_length)
            return aligned.cost
        var found = affine2p_cigar(
            a, b, mismatch, opening1, extension1, opening2, extension2, extended != 0, max_cost=max_cost, ends_free=ends
        )
        if not found:
            return ABOVE_MAX
        hand_over(found.value().cigar, cigar, cigar_length)
        return found.value().cost
    except:
        return INVALID_COSTS


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
