"""
The C API of dinara-align: the unit-cost edit distance and an optimal alignment as a CIGAR, for C, C++
and any language with a C foreign-function interface. `dinara.h` declares these, with a C++ wrapper.

    pixi run build-c [target-cpu]   # build/c: libdinara, its runtime libraries and the header

Each function takes the two sequences as bytes and lengths, holds nothing between calls and so may be
called from many threads at once, and reports failure as a negative result (`DINARA_*` in the header).
"""

from std.ffi import external_call

from dinara_align import edit_cigar, edit_distance

comptime UNSUPPORTED_SYMBOLS = -1
"""More than four symbols past `ACGT` between the two sequences."""
comptime OUT_OF_MEMORY = -2
"""The CIGAR's memory could not be allocated."""


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
        var bytes = aligned.cigar.as_bytes()
        var text = external_call["malloc", MutPointer[UInt8, MutAnyOrigin]](len(bytes) + 1)
        for index in range(len(bytes)):
            text[unsafe_offset=index] = bytes[index]
        text[unsafe_offset=len(bytes)] = 0
        cigar[] = text
        cigar_length[] = len(bytes)
        return aligned.distance
    except:
        return UNSUPPORTED_SYMBOLS


@export("dinara_free")
def dinara_free(text: MutPointer[UInt8, MutAnyOrigin]) abi("C"):
    """Frees what `dinara_edit_cigar` allocated, with the C library's `free` that matches its `malloc`."""
    external_call["free", NoneType](text)
