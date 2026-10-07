"""
The C API of dinara-align: the least cost of aligning a query to a reference, and an optimal alignment
as a CIGAR, under any `Costs` and `Mode`, for C, C++ and any language with a C foreign-function
interface. `dinara.h` declares these, with a C++ wrapper.

    pixi run build-c [target-cpu]   # build/c: libdinara, its runtime libraries and the header

Each function takes the two sequences as bytes and lengths, and the costs, the mode and the options as
structs of 64-bit integers, any of them null for its default. It holds nothing between calls, so may
be called from many threads at once, and reports failure as a negative result (`DINARA_*` in the
header).
"""

from std.ffi import external_call

from dinara_align import Alignment, AlignmentError, Anchor, Band, Costs, ErrorKind, Mode, Ties, align, distance

comptime UNSUPPORTED_SYMBOLS = -1
"""A 0xFE or 0xFF byte, which the wavefront's sentinels are and UTF-8 never holds."""
comptime OUT_OF_MEMORY = -2
"""The CIGAR's memory could not be allocated."""
comptime INVALID_COSTS = -3
"""Costs no alignment can be searched by: a free mismatch or extension, or a negative cost."""
comptime ABOVE_MAX = -4
"""Every alignment costs more than the `max_cost` asked for."""
comptime OUTSIDE_BAND = -5
"""No alignment stays inside the band asked for."""
comptime INVALID_MODE = -6
"""A mode that cannot serve what was asked: the least cost of an extension or a local alignment, which
maximize a score, a cap on either, or a band on a local alignment."""

comptime C_ENDS_FREE = 0
comptime C_EXTENSION = 1
comptime C_LOCAL = 2

comptime CInts = ImmPointer[Int, MutAnyOrigin]
"""A C struct of `int64_t` fields, read by index; null for the default."""


def sequence(bytes: ImmPointer[UInt8, MutAnyOrigin], length: Int) -> String:
    """A sequence from C bytes; an empty one may come as a null pointer, never read."""
    if length <= 0:
        return String()
    return String(StringSlice(unsafe_from_utf8=Span(unsafe_ptr=bytes, length=length)))


def plain_bytes(bytes: ImmPointer[UInt8, MutAnyOrigin], length: Int) -> Bool:
    """Whether no byte is one of the two UTF-8 never holds, which the wavefront's sentinels are."""
    for index in range(length):
        if bytes[unsafe_offset=index] >= 0xFE:
            return False
    return True


def costs_of(fields: OptionalPointer[Int, MutAnyOrigin]) raises AlignmentError -> Costs:
    """`dinara_costs`: mismatch, opening, extension, opening2, extension2, a negative opening2 for one
    piece. Null is unit costs."""
    if not fields:
        return Costs.edit()
    var at = fields.value()
    if at[unsafe_offset=3] < 0:
        return Costs.affine(at[unsafe_offset=0], at[unsafe_offset=1], at[unsafe_offset=2])
    return Costs.two_piece(
        at[unsafe_offset=0], at[unsafe_offset=1], at[unsafe_offset=2], at[unsafe_offset=3], at[unsafe_offset=4]
    )


def mode_of(fields: OptionalPointer[Int, MutAnyOrigin]) raises AlignmentError -> Mode:
    """`dinara_mode`: kind, the four free letter counts, the match score and the anchor. Null is global."""
    if not fields:
        return Mode.GLOBAL
    var at = fields.value()
    if at[unsafe_offset=0] == C_EXTENSION:
        return Mode.extension(at[unsafe_offset=5], Anchor.END if at[unsafe_offset=6] != 0 else Anchor.START)
    if at[unsafe_offset=0] == C_LOCAL:
        return Mode.local(at[unsafe_offset=5])
    if at[unsafe_offset=0] != C_ENDS_FREE:
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "an unknown mode")
    return Mode.ends_free(
        reference_start=at[unsafe_offset=1],
        reference_end=at[unsafe_offset=2],
        query_start=at[unsafe_offset=3],
        query_end=at[unsafe_offset=4],
    )


@fieldwise_init
struct Options(ImplicitlyCopyable):
    """`dinara_options`: the band's two edges, the cost cap, `extended` and `right_ties`. Null is no band,
    no cap, `=` and `X`, and indels placed left."""

    var band: Band
    var max_cost: Int
    var extended: Bool
    var ties: Ties


def options_of(fields: OptionalPointer[Int, MutAnyOrigin]) -> Options:
    if not fields:
        return Options(Band(), -1, True, Ties.LEFT)
    var at = fields.value()
    comptime EDGE = 1 << 60
    # The C integer limits stand for no band, kept clear of overflow.
    var band = Band(max(at[unsafe_offset=0], -EDGE), min(at[unsafe_offset=1], EDGE))
    var ties = Ties.RIGHT if at[unsafe_offset=4] != 0 else Ties.LEFT
    return Options(band, at[unsafe_offset=2], at[unsafe_offset=3] != 0, ties)


def failure(error: AlignmentError) -> Int:
    """The code for what the library raised."""
    if error.kind == ErrorKind.OUTSIDE_BAND:
        return OUTSIDE_BAND
    if error.kind == ErrorKind.INVALID_SCORING:
        return INVALID_COSTS
    if error.kind == ErrorKind.INVALID_ARGUMENT:
        return INVALID_MODE
    return UNSUPPORTED_SYMBOLS


@export("dinara_distance")
def dinara_distance(
    reference: ImmPointer[UInt8, MutAnyOrigin],
    reference_length: Int,
    query: ImmPointer[UInt8, MutAnyOrigin],
    query_length: Int,
    costs: OptionalPointer[Int, MutAnyOrigin],
    mode: OptionalPointer[Int, MutAnyOrigin],
    options: OptionalPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """The least cost of aligning the query to the reference (see `distance`), or a negative code:
    `ABOVE_MAX` when it passes a `max_cost` of zero or more, `OUTSIDE_BAND` when no alignment fits the
    band, or under a cap `ABOVE_MAX` again. The options' `extended` and `right_ties` change nothing."""
    if not plain_bytes(reference, reference_length) or not plain_bytes(query, query_length):
        return UNSUPPORTED_SYMBOLS
    var asked = options_of(options)
    try:
        var first = sequence(reference, reference_length)
        var second = sequence(query, query_length)
        var cap = asked.max_cost if asked.max_cost >= 0 else Int.MAX
        var found = distance(first, second, costs_of(costs), mode_of(mode), max_cost=cap, band=asked.band)
        if not found:
            # Under a cap, a band no alignment fits also leaves nothing within it.
            return OUTSIDE_BAND if cap == Int.MAX else ABOVE_MAX
        return found.value()
    except error:
        return failure(error)


@export("dinara_align")
def dinara_align(
    reference: ImmPointer[UInt8, MutAnyOrigin],
    reference_length: Int,
    query: ImmPointer[UInt8, MutAnyOrigin],
    query_length: Int,
    costs: OptionalPointer[Int, MutAnyOrigin],
    mode: OptionalPointer[Int, MutAnyOrigin],
    options: OptionalPointer[Int, MutAnyOrigin],
    alignment: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """An optimal alignment of the query to the reference (see `align`) into `alignment`, a
    `dinara_alignment`: its cost, score, the reference's and the query's aligned spans, and the CIGAR,
    NUL-terminated, in memory from C's `malloc` as its length is known only once the alignment is, which
    the caller frees with `dinara_free`. Returns zero, or a negative code, and then no CIGAR."""
    if not plain_bytes(reference, reference_length) or not plain_bytes(query, query_length):
        return UNSUPPORTED_SYMBOLS
    var asked = options_of(options)
    try:
        var first = sequence(reference, reference_length)
        var second = sequence(query, query_length)
        var wanted_mode = mode_of(mode)
        var found: Optional[Alignment]
        if asked.max_cost < 0 or wanted_mode.kind != Mode.ENDS:
            if asked.max_cost >= 0:
                return INVALID_MODE
            found = align(
                first, second, costs_of(costs), wanted_mode, band=asked.band, ties=asked.ties, extended=asked.extended
            )
        else:
            found = align(
                first,
                second,
                costs_of(costs),
                wanted_mode,
                max_cost=asked.max_cost,
                band=asked.band,
                ties=asked.ties,
                extended=asked.extended,
            )
            if not found:
                return ABOVE_MAX
        ref result = found.value()
        var text = copied(result.cigar)
        if not text:
            return OUT_OF_MEMORY
        alignment[unsafe_offset=0] = result.cost
        alignment[unsafe_offset=1] = result.score
        alignment[unsafe_offset=2] = result.reference_start
        alignment[unsafe_offset=3] = result.reference_end
        alignment[unsafe_offset=4] = result.query_start
        alignment[unsafe_offset=5] = result.query_end
        alignment.unsafe_offset(6).unsafe_bitcast[MutPointer[UInt8, MutAnyOrigin]]()[] = text.value()
        alignment[unsafe_offset=7] = result.cigar.byte_length()
        return 0
    except error:
        return failure(error)


def copied(text: String) -> OptionalPointer[UInt8, MutAnyOrigin]:
    """`text`, NUL-terminated, in memory from C's `malloc`, for the caller to free with `dinara_free`."""
    var bytes = text.as_bytes()
    var copy = external_call["malloc", OptionalPointer[UInt8, MutAnyOrigin]](len(bytes) + 1)
    if not copy:
        return None
    var at = copy.value()
    for index in range(len(bytes)):
        at[unsafe_offset=index] = bytes[index]
    at[unsafe_offset=len(bytes)] = 0
    return copy


@export("dinara_free")
def dinara_free(text: OptionalPointer[UInt8, MutAnyOrigin]) abi("C"):
    """Frees a CIGAR a function here allocated, with the C library's `free` that matches its `malloc`;
    null is nothing to free."""
    external_call["free", NoneType](text)
