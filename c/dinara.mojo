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
from std.ffi import external_call

from max.algorithm import parallelize

from dinara_align import (
    DEFAULT_MAX_MEMORY,
    Alignment,
    AlignmentError,
    Anchor,
    Band,
    Costs,
    ErrorKind,
    LocalScores,
    Mode,
    Ties,
    distance,
    local_scores,
    score,
    search,
)
from dinara_align.api import aligned_within
from dinara_align.common import hardware_threads, next_share
from dinara_align.gap_affine import KEPT_BYTES, Penalties, SearchSpace, cigar_of, penalties_of
from dinara_align.lanes import LaneCosts, Texts, lane_alignments, lane_distances

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
    return String(StringSlice(unsafe_from_utf8=Span(unsafe_ptr=bytes, length=length)))


def plain_bytes(bytes: ImmPointer[UInt8, MutAnyOrigin], length: Int) -> Bool:
    """Whether no byte is one of the two UTF-8 never holds, which the wavefront's sentinels are."""
    for index in range(length):
        if bytes[unsafe_offset=index] >= 0xFE:
            return False
    return True


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
    comptime EDGE = 1 << 60
    # The C integer limits stand for no band, kept clear of overflow.
    var band = Band(max(at[unsafe_offset=0], -EDGE), min(at[unsafe_offset=1], EDGE))
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
) -> Int:
    """One pair's least cost, or its code (see `dinara_distance`)."""
    if not plain_bytes(reference, reference_length) or not plain_bytes(query, query_length):
        return UNSUPPORTED_SYMBOLS
    try:
        var first = sequence(reference, reference_length)
        var second = sequence(query, query_length)
        var cap = asked.max_cost if asked.max_cost >= 0 else Int.MAX
        var found = distance(first, second, costs, mode, max_cost=cap, band=asked.band)
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
    if not plain_bytes(reference, reference_length) or not plain_bytes(query, query_length):
        return UNSUPPORTED_SYMBOLS
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
        return distance_code(
            reference, reference_length, query, query_length, costs_of(costs), mode_of(mode), options_of(options)
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
    if not plain_bytes(reference, reference_length) or not plain_bytes(query, query_length):
        return UNSUPPORTED_SYMBOLS
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
    if not plain_bytes(reference, reference_length) or not plain_bytes(query, query_length):
        return UNSUPPORTED_SYMBOLS
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
        return self.lengths[unsafe_offset=index]

    @always_inline
    def letters(self, index: Int) -> ImmPointer[UInt8, ImmUntrackedOrigin]:
        return ImmPointer[UInt8, ImmUntrackedOrigin](unsafe_from_address=self.starts[unsafe_offset=index])


def workers_for(pairs: Int, threads: Int) -> Int:
    """Threads for a batch of `pairs`: `threads`, or every thread this process may use for zero or
    fewer, and no more than there are pairs."""
    var workers = threads if threads > 0 else hardware_threads()
    return max(min(workers, pairs), 1)


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
    var taken = Atomic[Int64](0)
    var workers = workers_for(pairs, threads)

    # Costs with no reward, globally or with free ends, go many pairs at once into the lanes of a register,
    # as `distances` sends them (see `lanes`): every pair whose bytes the library takes and 16 bits hold. A pair the
    # library refuses, and every pair of costs it refuses, goes one at a time for its own code.
    var settled = List[Bool](length=pairs, fill=False)
    var found = List[Optional[Int]](length=pairs, fill=None)
    var laned = List[Bool](length=pairs, fill=False)
    var settled_ptr = settled.unsafe_ptr()
    var found_ptr = found.unsafe_ptr()
    var laned_ptr = laned.unsafe_ptr()
    var lane_costs = LaneCosts.of(wanted_costs, wanted_mode, free_ends=True)
    if lane_costs:
        try:
            _ = penalties_of(wanted_costs)
        except:
            lane_costs = None
    var cap = asked.max_cost if asked.max_cost >= 0 else Int.MAX
    if lane_costs:

        def refuse(stretch: Int) {imm}:
            """Marks stretch `stretch`'s pairs holding bytes the library refuses, which the lanes leave."""
            for index in range(pairs * stretch // workers, pairs * (stretch + 1) // workers):
                if not plain_bytes(references[unsafe_offset=index], reference_lengths[unsafe_offset=index]) or not (
                    plain_bytes(queries[unsafe_offset=index], query_lengths[unsafe_offset=index])
                ):
                    settled_ptr[unsafe_offset=index] = True

        parallelize(refuse, workers, workers)
        var before = settled.copy()
        _ = lane_distances(
            pairs,
            ByteTexts.of(references, reference_lengths),
            ByteTexts.of(queries, query_lengths),
            lane_costs.value(),
            asked.band,
            cap,
            workers,
            found_ptr,
            settled_ptr,
            wanted_mode,
        )
        for index in range(pairs):
            laned[index] = settled[index] and not before[index]

    def work(slot: Int) {mut taken, imm}:
        """Takes the next pairs not yet taken and writes each one's least cost or code, until none is left."""
        var last = 0
        while True:
            var share = next_share(taken, pairs, workers, last)
            if share[0] >= pairs:
                return
            for index in range(share[0], share[1]):
                if laned_ptr[unsafe_offset=index]:
                    var cost = found_ptr[unsafe_offset=index]
                    # Under a cap, a band no alignment fits also leaves nothing within it.
                    results[unsafe_offset=index] = cost.value() if cost else (
                        OUTSIDE_BAND if cap == Int.MAX else ABOVE_MAX
                    )
                    continue
                results[unsafe_offset=index] = distance_code(
                    references[unsafe_offset=index],
                    reference_lengths[unsafe_offset=index],
                    queries[unsafe_offset=index],
                    query_lengths[unsafe_offset=index],
                    wanted_costs,
                    wanted_mode,
                    asked,
                )

    parallelize(work, workers, workers)
    return 0


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
    var taken = Atomic[Int64](0)
    var workers = workers_for(pairs, threads)

    # Global alignments go many pairs at once into the lanes of a register, as `alignments` sends them
    # (see `lanes.lane_alignments`): every pair whose bytes the library takes. A pair the library refuses,
    # and every pair of costs it refuses or of other modes, goes one at a time for its own code.
    var capped = asked.max_cost >= 0
    var settled = List[Bool](length=pairs, fill=False)
    var found = List[Optional[Int]](length=pairs, fill=None)
    var laned = List[Bool](length=pairs, fill=False)
    var paths = List[List[UInt8]](capacity=pairs)
    for _ in range(pairs):
        paths.append(List[UInt8]())
    var settled_ptr = settled.unsafe_ptr()
    var found_ptr = found.unsafe_ptr()
    var laned_ptr = laned.unsafe_ptr()
    var path_ptr = paths.unsafe_ptr()
    var lane_costs = LaneCosts.of(wanted_costs, wanted_mode)
    var penalties: Optional[Penalties] = None
    if lane_costs:
        try:
            penalties = penalties_of(wanted_costs)
        except:
            lane_costs = None
    var limit = asked.max_memory // KEPT_BYTES
    if lane_costs:

        def refuse(stretch: Int) {imm}:
            """Marks stretch `stretch`'s pairs holding bytes the library refuses, which the lanes leave."""
            for index in range(pairs * stretch // workers, pairs * (stretch + 1) // workers):
                if not plain_bytes(references[unsafe_offset=index], reference_lengths[unsafe_offset=index]) or not (
                    plain_bytes(queries[unsafe_offset=index], query_lengths[unsafe_offset=index])
                ):
                    settled_ptr[unsafe_offset=index] = True

        parallelize(refuse, workers, workers)
        var before = settled.copy()
        _ = lane_alignments(
            pairs,
            ByteTexts.of(references, reference_lengths),
            ByteTexts.of(queries, query_lengths),
            lane_costs.value(),
            asked.band,
            asked.max_cost if capped else Int.MAX,
            asked.ties == Ties.LEFT,
            workers,
            asked.max_memory,
            found_ptr,
            path_ptr,
            settled_ptr,
        )
        for index in range(pairs):
            laned[index] = settled[index] and not before[index]

    def work(slot: Int) {mut taken, imm}:
        """Takes the next pairs not yet taken and writes each one's alignment and status, until none is left;
        a pair the lanes settled has its CIGAR spelled from their path."""
        var last = 0
        var space = SearchSpace()
        while True:
            var share = next_share(taken, pairs, workers, last)
            if share[0] >= pairs:
                return
            for index in range(share[0], share[1]):
                if laned_ptr[unsafe_offset=index]:
                    if not found_ptr[unsafe_offset=index]:
                        statuses[unsafe_offset=index] = ABOVE_MAX if capped else OUTSIDE_BAND
                        continue
                    var cost = found_ptr[unsafe_offset=index].value()
                    var columns = reference_lengths[unsafe_offset=index]
                    var rows = query_lengths[unsafe_offset=index]
                    ref scaled = penalties.value()
                    # As `alignments` keeps them: only a pair the searches would never split for memory.
                    if 2 * (cost // scaled.scale + 1) * (columns + rows + 1) <= limit:
                        var first = sequence(references[unsafe_offset=index], columns)
                        var second = sequence(queries[unsafe_offset=index], rows)
                        var moves = List[UInt8]()
                        swap(moves, path_ptr[unsafe_offset=index])
                        var cigar = cigar_of(first, second, moves^, cost // scaled.scale, scaled, asked.eqx)
                        statuses[unsafe_offset=index] = written(
                            Alignment(cost, -cost, cigar^, 0, columns, 0, rows),
                            alignments.unsafe_offset(index * ALIGNMENT_FIELDS),
                        )
                        continue
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

    parallelize(work, workers, workers)
    return 0


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
    negative code. `best` above zero keeps that many; the options' cap drops what passes it."""
    if not plain_bytes(query, query_length):
        return UNSUPPORTED_SYMBOLS
    try:
        var asked = options_of(options)
        var texts = List[String](capacity=count)
        for index in range(count):
            if not plain_bytes(references[unsafe_offset=index], reference_lengths[unsafe_offset=index]):
                return UNSUPPORTED_SYMBOLS
            texts.append(sequence(references[unsafe_offset=index], reference_lengths[unsafe_offset=index]))
        var hits = search(
            texts,
            sequence(query, query_length),
            costs_of(costs),
            mode_of(mode),
            best=Optional[Int](best) if best > 0 else None,
            max_cost=Optional[Int](asked.max_cost) if asked.max_cost >= 0 else None,
            threads=Optional[Int](threads) if threads > 0 else None,
        )
        for rank in range(len(hits)):
            indices[unsafe_offset=rank] = hits[rank].index
            scores[unsafe_offset=rank] = hits[rank].score
        return len(hits)
    except error:
        return failure(error)
