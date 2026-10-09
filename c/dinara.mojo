# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
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

from std.atomic import Atomic
from std.ffi import c_int, c_size_t, external_call
from std.sys import size_of


from dinara_align import (
    DEFAULT_MAX_MEMORY,
    Alignment,
    AlignmentError,
    Anchor,
    Band,
    Costs,
    ErrorKind,
    Hit,
    LocalScores,
    Mode,
    Ties,
    distance,
    local_scores,
    score,
    search,
)
from dinara_align.api import (
    aligned_within,
    alignment_from_lanes,
    alignments_in_lanes,
    cost_within,
    distances_in_lanes,
    each_pair,
    longest_first,
)
from dinara_align.cigar import text_of
from dinara_align.common import FIRST_SENTINEL, spread
from dinara_align.gap_affine import KEPT_BYTES, Penalties, SearchSpace, cigar_of, penalties_of
from dinara_align.lanes import LaneCosts, Texts, lane_alignments, lane_distances, lane_free_alignments

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
"""A mode that cannot serve what was asked: the least cost of an extension, a local alignment or an
overlap, which maximize a score, a cap on any of them, or a band on the last two."""
comptime INVALID_LENGTH = -7
"""A sequence of fewer than no letters."""

comptime C_ENDS_FREE = 0
"""`DINARA_ENDS_FREE`: a `dinara_mode`'s kind for free ends, global among them."""
comptime C_EXTENSION = 1
"""`DINARA_EXTENSION`: a `dinara_mode`'s kind for an extension from one end."""
comptime C_LOCAL = 2
"""`DINARA_LOCAL`: a `dinara_mode`'s kind for a local alignment."""
comptime C_OVERLAP = 3
"""`DINARA_OVERLAP`: a `dinara_mode`'s kind for an overlap, each sequence free at one end."""

comptime CInts = ImmPointer[Int, MutAnyOrigin]
"""A C struct of `int64_t` fields, read by index; null for the default."""


def sequence(bytes: ImmPointer[UInt8, MutAnyOrigin], length: Int) -> String:
    """A sequence from C bytes; an empty one may come as a null pointer, never read."""
    if length <= 0:
        return String()
    return text_of(Span(unsafe_ptr=bytes, length=length))


def plain_bytes(bytes: ImmPointer[UInt8, MutAnyOrigin], length: Int) -> Bool:
    """Whether no byte is one of the two UTF-8 never holds, which the wavefront's sentinels are."""
    for index in range(length):
        if bytes[unsafe_offset=index] >= FIRST_SENTINEL:
            return False
    return True


def pair_code(
    reference: ImmPointer[UInt8, MutAnyOrigin],
    reference_length: Int,
    query: ImmPointer[UInt8, MutAnyOrigin],
    query_length: Int,
) -> Int:
    """Zero for a pair the library takes, else why not: a negative length, or a byte it refuses (see
    `plain_bytes`)."""
    if reference_length < 0 or query_length < 0:
        return INVALID_LENGTH
    if not plain_bytes(reference, reference_length) or not plain_bytes(query, query_length):
        return UNSUPPORTED_SYMBOLS
    return 0


def costs_of(fields: OptionalPointer[Int, MutAnyOrigin]) raises AlignmentError -> Costs:
    """`dinara_costs`: mismatch, opening, extension, opening2, extension2, a negative opening2 for one
    piece, then a deletion's four, which count only with its extension above zero. Null is unit costs."""
    if not fields:
        return Costs.edit()
    var at = fields.value()
    var costs: Costs
    if at[unsafe_offset=3] < 0:
        costs = Costs.affine(at[unsafe_offset=0], at[unsafe_offset=1], at[unsafe_offset=2])
    else:
        costs = Costs.two_piece(
            at[unsafe_offset=0], at[unsafe_offset=1], at[unsafe_offset=2], at[unsafe_offset=3], at[unsafe_offset=4]
        )
    if at[unsafe_offset=6] > 0:
        costs = costs.with_deletions(at[unsafe_offset=5], at[unsafe_offset=6], at[unsafe_offset=7], at[unsafe_offset=8])
    return costs


def mode_of(fields: OptionalPointer[Int, MutAnyOrigin]) raises AlignmentError -> Mode:
    """`dinara_mode`: kind, the four free letter counts, the match score, the anchor, and the Z-drop and
    the end bonus, each none unless above zero. Null is global.
    Free ends with a match score of zero minimize the costs; above zero they maximize the score."""
    if not fields:
        return Mode.GLOBAL
    var at = fields.value()
    if at[unsafe_offset=0] == C_EXTENSION:
        var anchor = Anchor.END if at[unsafe_offset=6] != 0 else Anchor.START
        var zdrop = Optional[Int](at[unsafe_offset=7]) if at[unsafe_offset=7] > 0 else None
        var bonus = Optional[Int](at[unsafe_offset=8]) if at[unsafe_offset=8] > 0 else None
        return Mode.extension(at[unsafe_offset=5], anchor, zdrop=zdrop, end_bonus=bonus)
    if at[unsafe_offset=0] == C_LOCAL:
        return Mode.local(at[unsafe_offset=5])
    if at[unsafe_offset=0] == C_OVERLAP:
        return Mode.overlap(at[unsafe_offset=5])
    if at[unsafe_offset=0] != C_ENDS_FREE:
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "an unknown mode")
    return Mode.ends_free(
        reference_start=at[unsafe_offset=1],
        reference_end=at[unsafe_offset=2],
        query_start=at[unsafe_offset=3],
        query_end=at[unsafe_offset=4],
        match_score=at[unsafe_offset=5],
    )


@fieldwise_init
struct Options(ImplicitlyCopyable):
    """`dinara_options`: the band's two edges, the cost cap, `eqx`, `right_ties` and the memory for
    kept fronts. Null is no band, no cap, `=` and `X`, indels placed left, and the default memory."""

    var band: Band
    var max_cost: Int
    var eqx: Bool
    var ties: Ties
    var max_memory: Int


def options_of(fields: OptionalPointer[Int, MutAnyOrigin]) -> Options:
    """The `Options` a `dinara_options` asks for: a negative cap for none, and the default memory for zero or
    less."""
    if not fields:
        return Options(Band(), -1, True, Ties.LEFT, DEFAULT_MAX_MEMORY)
    var at = fields.value()
    # The C integer limits stand for no band.
    var band = Band(at[unsafe_offset=0], at[unsafe_offset=1])
    var ties = Ties.RIGHT if at[unsafe_offset=4] != 0 else Ties.LEFT
    var memory = at[unsafe_offset=5] if at[unsafe_offset=5] > 0 else DEFAULT_MAX_MEMORY
    return Options(band, at[unsafe_offset=2], at[unsafe_offset=3] != 0, ties, memory)


def failure(error: AlignmentError) -> Int:
    """The code for what the library raised."""
    if error.kind == ErrorKind.OUTSIDE_BAND:
        return OUTSIDE_BAND
    if error.kind == ErrorKind.INVALID_SCORING:
        return INVALID_COSTS
    if error.kind == ErrorKind.INVALID_ARGUMENT:
        return INVALID_MODE
    return UNSUPPORTED_SYMBOLS


def distance_code(
    reference: ImmPointer[UInt8, MutAnyOrigin],
    reference_length: Int,
    query: ImmPointer[UInt8, MutAnyOrigin],
    query_length: Int,
    costs: Costs,
    mode: Mode,
    asked: Options,
    mut space: SearchSpace,
) -> Int:
    """One pair's least cost, or its code (see `dinara_distance`), through `space`'s searches."""
    var refused = pair_code(reference, reference_length, query, query_length)
    if refused != 0:
        return refused
    try:
        var first = sequence(reference, reference_length)
        var second = sequence(query, query_length)
        var cap = asked.max_cost if asked.max_cost >= 0 else Int.MAX
        var found = cost_within(first, second, costs, mode, asked.band, cap, space)
        if not found:
            # Under a cap, a band no alignment fits also leaves nothing within it.
            return OUTSIDE_BAND if cap == Int.MAX else ABOVE_MAX
        return found.value()
    except error:
        return failure(error)


def align_into(
    reference: ImmPointer[UInt8, MutAnyOrigin],
    reference_length: Int,
    query: ImmPointer[UInt8, MutAnyOrigin],
    query_length: Int,
    costs: Costs,
    mode: Mode,
    asked: Options,
    alignment: MutPointer[Int, MutAnyOrigin],
    mut space: SearchSpace,
) -> Int:
    """One pair's optimal alignment into `alignment`, zero, or its code (see `dinara_align`), through
    `space`'s searches, which a batch's worker keeps from pair to pair."""
    var refused = pair_code(reference, reference_length, query, query_length)
    if refused != 0:
        return refused
    var capped = asked.max_cost >= 0
    if capped and mode.kind != Mode.ENDS:
        return INVALID_MODE
    try:
        var first = sequence(reference, reference_length)
        var second = sequence(query, query_length)
        # As `align` finds it, with or without the cap.
        var found = aligned_within(
            first,
            second,
            costs,
            mode,
            asked.band,
            asked.max_cost if capped else Int.MAX,
            asked.ties,
            asked.eqx,
            asked.max_memory // KEPT_BYTES,
            space,
        )
        if not found:
            return ABOVE_MAX if capped else OUTSIDE_BAND
        return written(found.value(), alignment)
    except error:
        return failure(error)


def written(result: Alignment, alignment: MutPointer[Int, MutAnyOrigin]) -> Int:
    """`result` into `alignment`, a `dinara_alignment`, its CIGAR copied for C; zero, or `OUT_OF_MEMORY`."""
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


comptime ALIGNMENT_FIELDS = 8
"""`int64_t`s a `dinara_alignment` spans, its CIGAR's pointer one of them."""


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
    band, or under a cap `ABOVE_MAX` again. The options' `eqx` and `right_ties` change nothing."""
    try:
        var space = SearchSpace()
        return distance_code(
            reference,
            reference_length,
            query,
            query_length,
            costs_of(costs),
            mode_of(mode),
            options_of(options),
            space,
        )
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
    try:
        var space = SearchSpace()
        return align_into(
            reference,
            reference_length,
            query,
            query_length,
            costs_of(costs),
            mode_of(mode),
            options_of(options),
            alignment,
            space,
        )
    except error:
        return failure(error)


@export("dinara_aligner_new")
def dinara_aligner_new() abi("C") -> OptionalPointer[SearchSpace, MutAnyOrigin]:
    """One thread's aligner, its searches' memory kept from call to call, in memory from C's `malloc`, for
    `dinara_aligner_free`; null when there is none to be had. It keeps nothing a call depends on, so it
    gives the same answers as the functions without one, and it is never shared between threads."""
    var aligner = external_call["malloc", OptionalPointer[SearchSpace, MutAnyOrigin]](size_of[SearchSpace]())
    if aligner:
        aligner.value().unsafe_write(SearchSpace())
    return aligner


@export("dinara_aligner_free")
def dinara_aligner_free(aligner: OptionalPointer[SearchSpace, MutAnyOrigin]) abi("C"):
    """Frees an aligner `dinara_aligner_new` made, and the memory it kept; null is nothing to free."""
    if aligner:
        _ = aligner.value().unsafe_take_pointee()
        external_call["free", NoneType](aligner)


@export("dinara_aligner_distance")
def dinara_aligner_distance(
    aligner: MutPointer[SearchSpace, MutAnyOrigin],
    reference: ImmPointer[UInt8, MutAnyOrigin],
    reference_length: Int,
    query: ImmPointer[UInt8, MutAnyOrigin],
    query_length: Int,
    costs: OptionalPointer[Int, MutAnyOrigin],
    mode: OptionalPointer[Int, MutAnyOrigin],
    options: OptionalPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """`dinara_distance` through `aligner`'s memory."""
    try:
        return distance_code(
            reference,
            reference_length,
            query,
            query_length,
            costs_of(costs),
            mode_of(mode),
            options_of(options),
            aligner[],
        )
    except error:
        return failure(error)


@export("dinara_aligner_align")
def dinara_aligner_align(
    aligner: MutPointer[SearchSpace, MutAnyOrigin],
    reference: ImmPointer[UInt8, MutAnyOrigin],
    reference_length: Int,
    query: ImmPointer[UInt8, MutAnyOrigin],
    query_length: Int,
    costs: OptionalPointer[Int, MutAnyOrigin],
    mode: OptionalPointer[Int, MutAnyOrigin],
    options: OptionalPointer[Int, MutAnyOrigin],
    alignment: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """`dinara_align` through `aligner`'s memory."""
    try:
        return align_into(
            reference,
            reference_length,
            query,
            query_length,
            costs_of(costs),
            mode_of(mode),
            options_of(options),
            alignment,
            aligner[],
        )
    except error:
        return failure(error)


@export("dinara_score")
def dinara_score(
    reference: ImmPointer[UInt8, MutAnyOrigin],
    reference_length: Int,
    query: ImmPointer[UInt8, MutAnyOrigin],
    query_length: Int,
    costs: OptionalPointer[Int, MutAnyOrigin],
    mode: OptionalPointer[Int, MutAnyOrigin],
    options: OptionalPointer[Int, MutAnyOrigin],
    found: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """The best score `dinara_align` would return, with no alignment traced (see `score`), into
    `found`: zero, or a negative code."""
    var refused = pair_code(reference, reference_length, query, query_length)
    if refused != 0:
        return refused
    try:
        var asked = options_of(options)
        found[] = score(
            sequence(reference, reference_length),
            sequence(query, query_length),
            costs_of(costs),
            mode_of(mode),
            band=asked.band,
        )
        return 0
    except error:
        return failure(error)


@export("dinara_local_scores")
def dinara_local_scores(
    reference: ImmPointer[UInt8, MutAnyOrigin],
    reference_length: Int,
    query: ImmPointer[UInt8, MutAnyOrigin],
    query_length: Int,
    costs: OptionalPointer[Int, MutAnyOrigin],
    mode: OptionalPointer[Int, MutAnyOrigin],
    window: Int,
    found: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """A local alignment's best score, its end and SSW's second best (see `local_scores`) into
    `found`, a `dinara_local_scores_result`: zero, or a negative code."""
    var refused = pair_code(reference, reference_length, query, query_length)
    if refused != 0:
        return refused
    try:
        var scores: LocalScores
        if window >= 0:
            scores = local_scores(
                sequence(reference, reference_length),
                sequence(query, query_length),
                costs_of(costs),
                mode_of(mode),
                window=window,
            )
        else:
            scores = local_scores(
                sequence(reference, reference_length), sequence(query, query_length), costs_of(costs), mode_of(mode)
            )
        found[unsafe_offset=0] = scores.score
        found[unsafe_offset=1] = scores.reference_end
        found[unsafe_offset=2] = scores.query_end
        found[unsafe_offset=3] = scores.second_score
        found[unsafe_offset=4] = scores.second_reference_end
        return 0
    except error:
        return failure(error)


comptime CSequences = ImmPointer[ImmPointer[UInt8, MutAnyOrigin], MutAnyOrigin]
"""A C array of sequences' first bytes."""


@fieldwise_init
struct ByteTexts(Texts, TrivialRegisterPassable):
    """A C caller's sequences for the lanes (see `lanes.Texts`): its arrays of first bytes, held as
    addresses, and of lengths."""

    var starts: ImmPointer[Int, ImmUntrackedOrigin]
    var lengths: ImmPointer[Int, ImmUntrackedOrigin]

    @staticmethod
    def of(sequences: CSequences, lengths: CInts) -> Self:
        return Self(
            sequences.unsafe_bitcast[Int]().unsafe_origin_cast[ImmUntrackedOrigin](),
            lengths.unsafe_origin_cast[ImmUntrackedOrigin](),
        )

    @always_inline
    def length(self, index: Int) -> Int:
        """Pair `index`'s length, none for a negative one, which `refused_pairs` keeps from the lanes."""
        return max(self.lengths[unsafe_offset=index], 0)

    @always_inline
    def letters(self, index: Int) -> ImmPointer[UInt8, ImmUntrackedOrigin]:
        return ImmPointer[UInt8, ImmUntrackedOrigin](unsafe_from_address=self.starts[unsafe_offset=index])


def workers_for(pairs: Int, threads: Int) -> Int:
    """Threads for a batch of `pairs`: `threads`, or the caller's own alone for zero or fewer, and no more
    than there are pairs. The caller spreads its own calls over its threads; the library starts more only
    when asked."""
    return max(min(threads, pairs), 1)


def refused_pairs(
    pairs: Int,
    references: CSequences,
    reference_lengths: CInts,
    queries: CSequences,
    query_lengths: CInts,
    workers: Int,
) -> List[Bool]:
    """Which pairs the library refuses (see `pair_code`), on `workers` threads: the lanes leave them, and
    each is answered with its own code."""
    var refused = List[Bool](length=max(pairs, 1), fill=False)
    var refused_ptr = refused.unsafe_ptr()

    def refuse(stretch: Int) {imm}:
        """Marks stretch `stretch`'s pairs holding such a byte."""
        for index in range(pairs * stretch // workers, pairs * (stretch + 1) // workers):
            if (
                pair_code(
                    references[unsafe_offset=index],
                    reference_lengths[unsafe_offset=index],
                    queries[unsafe_offset=index],
                    query_lengths[unsafe_offset=index],
                )
                != 0
            ):
                refused_ptr[unsafe_offset=index] = True

    spread(refuse, workers, workers)
    return refused^


@export("dinara_distances")
def dinara_distances(
    pairs: Int,
    references: CSequences,
    reference_lengths: CInts,
    queries: CSequences,
    query_lengths: CInts,
    costs: OptionalPointer[Int, MutAnyOrigin],
    mode: OptionalPointer[Int, MutAnyOrigin],
    options: OptionalPointer[Int, MutAnyOrigin],
    threads: Int,
    results: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """Every pair's `dinara_distance` into `results`, each its least cost or its own code, the pairs
    spread over `threads` threads; zero, or a code that fails them all, and then no result written."""
    var wanted_costs: Costs
    var wanted_mode: Mode
    try:
        wanted_costs = costs_of(costs)
        wanted_mode = mode_of(mode)
    except error:
        return failure(error)
    var asked = options_of(options)

    if pairs <= 0:
        return 0
    var workers = workers_for(pairs, threads)
    if workers == 1:
        distances_part(
            pairs, references, reference_lengths, queries, query_lengths, wanted_costs, wanted_mode, asked, results
        )
        return 0
    var job = SharedBatch(
        DISTANCES,
        dealt_chunks(pair_weights(pairs, reference_lengths, query_lengths), workers),
        references,
        reference_lengths,
        queries,
        query_lengths,
        wanted_costs,
        wanted_mode,
        asked,
        results,
        results,
    )
    shared_out(job, workers)
    return 0


def distances_part(
    pairs: Int,
    references: CSequences,
    reference_lengths: CInts,
    queries: CSequences,
    query_lengths: CInts,
    wanted_costs: Costs,
    wanted_mode: Mode,
    asked: Options,
    results: MutPointer[Int, MutAnyOrigin],
):
    """`dinara_distances` of `pairs` pairs on the caller's thread alone."""
    var workers = 1
    # As `distances` sends them (see `api.distances_in_lanes`): the lanes first, every pair whose bytes the library
    # takes; a pair it refuses, and every pair of costs it refuses, one at a time for its own code, the longest
    # first.
    var cap = asked.max_cost if asked.max_cost >= 0 else Int.MAX
    var settled = refused_pairs(pairs, references, reference_lengths, queries, query_lengths, workers)
    var before = settled.copy()
    var found = List[Optional[Int]](length=pairs, fill=None)
    var found_ptr = found.unsafe_ptr()
    var reference_texts = ByteTexts.of(references, reference_lengths)
    var query_texts = ByteTexts.of(queries, query_lengths)
    _ = distances_in_lanes(
        pairs,
        reference_texts,
        query_texts,
        wanted_costs,
        wanted_mode,
        asked.band,
        cap,
        workers,
        found_ptr,
        settled.unsafe_ptr(),
    )
    var laned = List[Bool](length=pairs, fill=False)
    for index in range(pairs):
        laned[index] = settled[index] and not before[index]
    var laned_ptr = laned.unsafe_ptr()

    def one(index: Int, mut space: SearchSpace) {imm}:
        """Pair `index`'s least cost or its code, from the lanes where they settled it."""
        if laned_ptr[unsafe_offset=index]:
            var cost = found_ptr[unsafe_offset=index]
            # Under a cap, a band no alignment fits also leaves nothing within it.
            results[unsafe_offset=index] = cost.value() if cost else (OUTSIDE_BAND if cap == Int.MAX else ABOVE_MAX)
            return
        results[unsafe_offset=index] = distance_code(
            references[unsafe_offset=index],
            reference_lengths[unsafe_offset=index],
            queries[unsafe_offset=index],
            query_lengths[unsafe_offset=index],
            wanted_costs,
            wanted_mode,
            asked,
            space,
        )

    each_pair(one, longest_first(pairs, reference_texts, query_texts, workers), workers)


@export("dinara_alignments")
def dinara_alignments(
    pairs: Int,
    references: CSequences,
    reference_lengths: CInts,
    queries: CSequences,
    query_lengths: CInts,
    costs: OptionalPointer[Int, MutAnyOrigin],
    mode: OptionalPointer[Int, MutAnyOrigin],
    options: OptionalPointer[Int, MutAnyOrigin],
    threads: Int,
    alignments: MutPointer[Int, MutAnyOrigin],
    statuses: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """Every pair's `dinara_align` into `alignments`, an array of `dinara_alignment`, and its result
    into `statuses`: zero and a CIGAR to free, or the pair's code and none. The pairs spread over
    `threads` threads; returns zero, or a code that fails them all, and then nothing written."""
    var wanted_costs: Costs
    var wanted_mode: Mode
    try:
        wanted_costs = costs_of(costs)
        wanted_mode = mode_of(mode)
    except error:
        return failure(error)
    var asked = options_of(options)

    if pairs <= 0:
        return 0
    var workers = workers_for(pairs, threads)
    if workers == 1:
        alignments_part(
            pairs,
            references,
            reference_lengths,
            queries,
            query_lengths,
            wanted_costs,
            wanted_mode,
            asked,
            alignments,
            statuses,
        )
        return 0
    var job = SharedBatch(
        ALIGNMENTS,
        dealt_chunks(pair_weights(pairs, reference_lengths, query_lengths), workers),
        references,
        reference_lengths,
        queries,
        query_lengths,
        wanted_costs,
        wanted_mode,
        asked,
        alignments,
        statuses,
    )
    shared_out(job, workers)
    return 0


def alignments_part(
    pairs: Int,
    references: CSequences,
    reference_lengths: CInts,
    queries: CSequences,
    query_lengths: CInts,
    wanted_costs: Costs,
    wanted_mode: Mode,
    asked: Options,
    alignments: MutPointer[Int, MutAnyOrigin],
    statuses: MutPointer[Int, MutAnyOrigin],
):
    """`dinara_alignments` of `pairs` pairs on the caller's thread alone."""
    var workers = 1
    # As `alignments` sends them (see `api.alignments_in_lanes`): the lanes first, every pair whose bytes the
    # library takes, each one's CIGAR spelled from their path; a pair the library refuses, every pair of costs it
    # refuses or of a mode that scores, and a pair the searches might split, one at a time for its own code, the
    # longest first.
    var capped = asked.max_cost >= 0
    var settled = refused_pairs(pairs, references, reference_lengths, queries, query_lengths, workers)
    var before = settled.copy()
    var found = List[Optional[Int]](length=pairs, fill=None)
    var paths = List[List[UInt8]](capacity=pairs)
    for _ in range(pairs):
        paths.append(List[UInt8]())
    var spans = List[Int](length=4 * pairs, fill=0)
    var found_ptr = found.unsafe_ptr()
    var path_ptr = paths.unsafe_ptr()
    var span_ptr = spans.unsafe_ptr()
    var limit = asked.max_memory // KEPT_BYTES
    var reference_texts = ByteTexts.of(references, reference_lengths)
    var query_texts = ByteTexts.of(queries, query_lengths)
    var penalties = alignments_in_lanes(
        pairs,
        reference_texts,
        query_texts,
        wanted_costs,
        wanted_mode,
        asked.band,
        asked.max_cost if capped else Int.MAX,
        asked.ties,
        workers,
        limit,
        found_ptr,
        path_ptr,
        span_ptr,
        settled.unsafe_ptr(),
    )
    var laned = List[Bool](length=pairs, fill=False)
    for index in range(pairs):
        laned[index] = settled[index] and not before[index]
    var laned_ptr = laned.unsafe_ptr()

    def one(index: Int, mut space: SearchSpace) {imm}:
        """Pair `index`'s alignment and status, a pair the lanes settled spelled from their path."""
        if laned_ptr[unsafe_offset=index]:
            if not found_ptr[unsafe_offset=index]:
                statuses[unsafe_offset=index] = ABOVE_MAX if capped else OUTSIDE_BAND
                return
            var moves = List[UInt8]()
            swap(moves, path_ptr[unsafe_offset=index])
            var spelled = alignment_from_lanes(
                Span(unsafe_ptr=references[unsafe_offset=index], length=reference_lengths[unsafe_offset=index]),
                Span(unsafe_ptr=queries[unsafe_offset=index], length=query_lengths[unsafe_offset=index]),
                found_ptr[unsafe_offset=index].value(),
                moves^,
                span_ptr.unsafe_offset(4 * index),
                penalties.value(),
                asked.eqx,
                limit,
            )
            if spelled:
                statuses[unsafe_offset=index] = written(
                    spelled.value(), alignments.unsafe_offset(index * ALIGNMENT_FIELDS)
                )
                return
        statuses[unsafe_offset=index] = align_into(
            references[unsafe_offset=index],
            reference_lengths[unsafe_offset=index],
            queries[unsafe_offset=index],
            query_lengths[unsafe_offset=index],
            wanted_costs,
            wanted_mode,
            asked,
            alignments.unsafe_offset(index * ALIGNMENT_FIELDS),
            space,
        )

    each_pair(one, longest_first(pairs, reference_texts, query_texts, workers), workers)


comptime DISTANCES = 0
"""A `SharedBatch` of `dinara_distances`."""
comptime ALIGNMENTS = 1
"""A `SharedBatch` of `dinara_alignments`."""

comptime THREAD_STACK = 8 << 20
"""The stack a thread the library starts gets: a process's main thread's, where a platform's default for
other threads, 512 KiB on macOS, is far less."""


def pair_weights(pairs: Int, reference_lengths: CInts, query_lengths: CInts) -> List[Int]:
    """Each pair's letters, by which `dealt_chunks` orders the work."""
    var weights = List[Int](capacity=pairs)
    for index in range(pairs):
        weights.append(max(reference_lengths[unsafe_offset=index], 0) + max(query_lengths[unsafe_offset=index], 0))
    return weights^


def dealt_chunks(weights: List[Int], workers: Int) -> List[Int]:
    """The items in eight chunks a worker, or one an item when there are fewer, as `api.scores` deals a
    batch: each chunk's first item and the one past its last, the heaviest chunk first, so the longest
    work starts first and the rest fills in around it."""
    var items = len(weights)
    var count = max(min(items, workers * 8), 1)
    var totals = List[Int](length=count, fill=0)
    for chunk in range(count):
        for index in range(items * chunk // count, items * (chunk + 1) // count):
            totals[chunk] += weights[index]
    var order = List[Int](capacity=count)
    for chunk in range(count):
        order.append(chunk)

    def heavier(left: Int, right: Int) {imm totals} -> Bool:
        """Whether chunk `left` goes before `right`: the heavier, then the earlier."""
        if totals[left] != totals[right]:
            return totals[left] > totals[right]
        return left < right

    sort(order, heavier)
    var bounds = List[Int](capacity=2 * count)
    for chunk in order:
        bounds.append(items * chunk // count)
        bounds.append(items * (chunk + 1) // count)
    return bounds^


struct SharedBatch(Movable):
    """A C batch dealt out over threads the library starts itself (see `shared_out`), each chunk of it a
    batch of its own on one thread.

    Not Mojo's own threads: two of a caller's threads each asking a batch for several at once could leave
    one waiting on the runtime for good, and a C caller's threads are its own to run as it likes."""

    var kind: Int
    """`DISTANCES` or `ALIGNMENTS`."""
    var bounds: List[Int]
    """The chunks, in the order they are taken (see `dealt_chunks`)."""
    var next: Atomic[Int64]
    """The next chunk to take."""
    var references: Int
    """The caller's arrays, as addresses: what a C caller hands over outlives its call."""
    var reference_lengths: Int
    var queries: Int
    var query_lengths: Int
    var costs: Costs
    var mode: Mode
    var asked: Options
    var results: Int
    """Distances' results, or alignments' records."""
    var statuses: Int
    """Alignments' statuses."""

    def __init__(
        out self,
        kind: Int,
        var bounds: List[Int],
        references: CSequences,
        reference_lengths: CInts,
        queries: CSequences,
        query_lengths: CInts,
        costs: Costs,
        mode: Mode,
        asked: Options,
        results: MutPointer[Int, MutAnyOrigin],
        statuses: MutPointer[Int, MutAnyOrigin],
    ):
        self.kind = kind
        self.bounds = bounds^
        self.next = Atomic[Int64](0)
        self.references = Int(references)
        self.reference_lengths = Int(reference_lengths)
        self.queries = Int(queries)
        self.query_lengths = Int(query_lengths)
        self.costs = costs
        self.mode = mode
        self.asked = asked
        self.results = Int(results)
        self.statuses = Int(statuses)

    def run(self, chunk: Int):
        """Chunk `chunk`'s pairs, its own batch."""
        var first = self.bounds[2 * chunk]
        var pairs = self.bounds[2 * chunk + 1] - first
        if pairs <= 0:
            return
        var references = CSequences(unsafe_from_address=self.references).unsafe_offset(first)
        var reference_lengths = CInts(unsafe_from_address=self.reference_lengths).unsafe_offset(first)
        var queries = CSequences(unsafe_from_address=self.queries).unsafe_offset(first)
        var query_lengths = CInts(unsafe_from_address=self.query_lengths).unsafe_offset(first)
        var results = MutPointer[Int, MutAnyOrigin](unsafe_from_address=self.results)
        if self.kind == DISTANCES:
            distances_part(
                pairs,
                references,
                reference_lengths,
                queries,
                query_lengths,
                self.costs,
                self.mode,
                self.asked,
                results.unsafe_offset(first),
            )
        else:
            alignments_part(
                pairs,
                references,
                reference_lengths,
                queries,
                query_lengths,
                self.costs,
                self.mode,
                self.asked,
                results.unsafe_offset(first * ALIGNMENT_FIELDS),
                MutPointer[Int, MutAnyOrigin](unsafe_from_address=self.statuses).unsafe_offset(first),
            )


comptime Raw = MutPointer[NoneType, MutAnyOrigin]


def batch_worker(argument: Raw) abi("C") -> Raw:
    """A thread's share of a `SharedBatch`: the next chunk until none is left."""
    ref job = argument.unsafe_bitcast[SharedBatch]()[]
    var chunks = len(job.bounds) // 2
    while True:
        var chunk = Int(job.next.fetch_add(1))
        if chunk >= chunks:
            return argument
        job.run(chunk)


def shared_out(mut job: SharedBatch, workers: Int):
    """`job` on `workers` threads (see `start_threads`)."""
    start_threads(batch_worker, Pointer(to=job).unsafe_origin_cast[MutAnyOrigin]().unsafe_bitcast[NoneType](), workers)


def start_threads(entry: def(Raw) abi("C") thin -> Raw, argument: Raw, workers: Int):
    """`entry(argument)` on `workers` threads: the caller's own and as many more as the library can start,
    each taking the next chunk until none is left, so a thread that fails to start leaves its share to
    the rest."""
    if workers <= 1:
        _ = entry(argument)
        return
    # Room for any platform's `pthread_attr_t`.
    var attributes = List[UInt64](length=16, fill=0)
    var attribute = attributes.unsafe_ptr().unsafe_bitcast[NoneType]().unsafe_origin_cast[MutAnyOrigin]()
    var handles = List[UInt64](length=workers, fill=0)
    var started = List[Bool](length=workers, fill=False)
    if Int(external_call["pthread_attr_init", c_int](attribute)) == 0:
        _ = external_call["pthread_attr_setstacksize", c_int](attribute, c_size_t(THREAD_STACK))
        for index in range(workers - 1):
            var created = external_call["pthread_create", c_int](
                handles.unsafe_ptr().unsafe_offset(index), attribute, entry, argument
            )
            started[index] = Int(created) == 0
        _ = external_call["pthread_attr_destroy", c_int](attribute)
    _ = entry(argument)
    for index in range(workers - 1):
        if started[index]:
            _ = external_call["pthread_join", c_int](handles[index], OptionalPointer[NoneType, MutAnyOrigin]())
    _ = attributes^


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


@export("dinara_search")
def dinara_search(
    count: Int,
    references: CSequences,
    reference_lengths: CInts,
    query: ImmPointer[UInt8, MutAnyOrigin],
    query_length: Int,
    costs: OptionalPointer[Int, MutAnyOrigin],
    mode: OptionalPointer[Int, MutAnyOrigin],
    options: OptionalPointer[Int, MutAnyOrigin],
    best: Int,
    threads: Int,
    indices: MutPointer[Int, MutAnyOrigin],
    scores: MutPointer[Int, MutAnyOrigin],
) abi("C") -> Int:
    """The query against `count` references (see `search`): the hits, best first, their places into
    `indices` and their scores into `scores`, each with room for `count`, and their number back, or a
    negative code. `best` above zero keeps that many; the options' cap drops what passes it. The
    references spread over `threads` threads the library starts itself (see `SharedSearch`)."""
    if count <= 0:
        return 0
    for index in range(count):
        var refused = pair_code(
            references[unsafe_offset=index], reference_lengths[unsafe_offset=index], query, query_length
        )
        if refused != 0:
            return refused
    var wanted_costs: Costs
    var wanted_mode: Mode
    try:
        wanted_costs = costs_of(costs)
        wanted_mode = mode_of(mode)
    except error:
        return failure(error)
    var asked = options_of(options)
    var texts = List[String](capacity=count)
    for index in range(count):
        texts.append(sequence(references[unsafe_offset=index], reference_lengths[unsafe_offset=index]))
    var workers = workers_for(count, threads)
    var weights = List[Int](capacity=count)
    for index in range(count):
        weights.append(reference_lengths[unsafe_offset=index])
    var job = SharedSearch(
        dealt_chunks(weights, workers),
        texts^,
        sequence(query, query_length),
        wanted_costs,
        wanted_mode,
        asked.max_cost,
        best,
    )
    start_threads(search_worker, Pointer(to=job).unsafe_origin_cast[MutAnyOrigin]().unsafe_bitcast[NoneType](), workers)
    # The code a serial search raises: the first chunk's, in the references' order, holding one.
    var chunks = len(job.codes)
    var first_code = 0
    var first_at = Int.MAX
    for chunk in range(chunks):
        if job.codes[chunk] != 0 and job.bounds[2 * chunk] < first_at:
            first_code = job.codes[chunk]
            first_at = job.bounds[2 * chunk]
    if first_code != 0:
        return first_code
    # The best first, ties by the references' order, as `search` ranks them.
    var places = List[Int]()
    var values = List[Int]()
    for chunk in range(chunks):
        for hit in job.hits[chunk]:
            places.append(hit.index + job.bounds[2 * chunk])
            values.append(hit.score)
    var order = List[Int](capacity=len(places))
    for rank in range(len(places)):
        order.append(rank)

    def ahead(left: Int, right: Int) {imm places, imm values} -> Bool:
        """Whether hit `left` ranks before `right`: the higher score, then the earlier reference."""
        if values[left] != values[right]:
            return values[left] > values[right]
        return places[left] < places[right]

    sort(order, ahead)
    var kept = min(len(order), best) if best > 0 else len(order)
    for rank in range(kept):
        indices[unsafe_offset=rank] = places[order[rank]]
        scores[unsafe_offset=rank] = values[order[rank]]
    return kept


struct SharedSearch(Movable):
    """`dinara_search` dealt out as `SharedBatch` deals a batch: each chunk of the references searched on
    its own, its hits kept for the merge, which ranks them as one search does; a chunk's best `best` hold
    every hit of the search's best."""

    var bounds: List[Int]
    var next: Atomic[Int64]
    var texts: List[String]
    var query: String
    var costs: Costs
    var mode: Mode
    var max_cost: Int
    var best: Int
    var hits: List[List[Hit]]
    """Each chunk's hits, its own references numbered from zero."""
    var codes: List[Int]
    """Each chunk's code, zero when it searched."""

    def __init__(
        out self,
        var bounds: List[Int],
        var texts: List[String],
        var query: String,
        costs: Costs,
        mode: Mode,
        max_cost: Int,
        best: Int,
    ):
        var chunks = len(bounds) // 2
        self.bounds = bounds^
        self.next = Atomic[Int64](0)
        self.texts = texts^
        self.query = query^
        self.costs = costs
        self.mode = mode
        self.max_cost = max_cost
        self.best = best
        self.hits = List[List[Hit]](capacity=chunks)
        for _ in range(chunks):
            self.hits.append(List[Hit]())
        self.codes = List[Int](length=chunks, fill=0)

    def run(mut self, chunk: Int):
        """Chunk `chunk`'s references searched, its hits or its code kept."""
        var part = List[String](capacity=self.bounds[2 * chunk + 1] - self.bounds[2 * chunk])
        for index in range(self.bounds[2 * chunk], self.bounds[2 * chunk + 1]):
            part.append(self.texts[index])
        try:
            self.hits[chunk] = search(
                part,
                self.query,
                self.costs,
                self.mode,
                best=Optional[Int](self.best) if self.best > 0 else None,
                max_cost=Optional[Int](self.max_cost) if self.max_cost >= 0 else None,
            )
        except error:
            self.codes[chunk] = failure(error)


def search_worker(argument: Raw) abi("C") -> Raw:
    """A thread's share of a `SharedSearch`: the next chunk until none is left."""
    ref job = argument.unsafe_bitcast[SharedSearch]()[]
    var chunks = len(job.bounds) // 2
    while True:
        var chunk = Int(job.next.fetch_add(1))
        if chunk >= chunks:
            return argument
        job.run(chunk)
