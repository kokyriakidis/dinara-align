"""
The C API of dinara-align: the unit-cost edit distance, and the least gap-affine cost as WFA counts it,
each with an optimal alignment as a CIGAR, for C, C++ and any language with a C foreign-function
interface. `dinara.h` declares these, with a C++ wrapper.

    pixi run build-c [target-cpu]   # build/c: libdinara, its runtime libraries and the header

Each function takes the two sequences as bytes and lengths, holds nothing between calls and so may be
called from many threads at once, and reports failure as a negative result (`DINARA_*` in the header).
"""

from std.ffi import external_call

from dinara_align import affine_cigar, edit_cigar, edit_distance

comptime UNSUPPORTED_SYMBOLS = -1
"""More than four symbols past `ACGT` between the two sequences."""
comptime OUT_OF_MEMORY = -2
"""The CIGAR's memory could not be allocated."""
comptime INVALID_COSTS = -3
"""Gap-affine costs no alignment can be searched by: a free mismatch or extension, or a negative cost."""


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


@export("dinara_affine_cigar")
def dinara_affine_cigar(
    first: ImmPointer[UInt8, MutAnyOrigin],
    first_length: Int,
    second: ImmPointer[UInt8, MutAnyOrigin],
    second_length: Int,
    mismatch: Int,
    opening: Int,
    extension: Int,
    extended: Int32,
    cigar: MutPointer[MutPointer[UInt8, MutAnyOrigin], MutAnyOrigin],
    cigar_length: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """The least global cost under gap-affine costs, a substitution `mismatch` and a gap of `k` letters
    `opening + k extension`, and an optimal alignment's CIGAR, handed over as `dinara_edit_cigar`'s.
    Every byte is a symbol matching only itself, save the two UTF-8 never holds, which mark the ends."""
    # The wavefront's sentinels are the two bytes UTF-8 never uses.
    for index in range(first_length):
        if first[unsafe_offset=index] >= 0xFE:
            return UNSUPPORTED_SYMBOLS
    for index in range(second_length):
        if second[unsafe_offset=index] >= 0xFE:
            return UNSUPPORTED_SYMBOLS
    try:
        var aligned = affine_cigar(
            sequence(first, first_length), sequence(second, second_length), mismatch, opening, extension, extended != 0
        )
        hand_over(aligned.cigar, cigar, cigar_length)
        return aligned.cost
    except:
        return INVALID_COSTS


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
