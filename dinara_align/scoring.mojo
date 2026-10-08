# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# Derived from AffineGaps (https://github.com/unum-science/AffineGaps), Copyright Ash Vardanian, under the
# Apache License, Version 2.0, and changed since: see LICENSES/Apache-2.0.txt and NOTICE.
"""
Alignment under a `Scoring`, an alphabet's substitution table and affine gap scores, which a score
maximizes: Gotoh's Needleman-Wunsch and Smith-Waterman, on the host or the device, each pair served by
whichever kernel it suits so the caller never has to choose.

Two choices are made here and nowhere else, and both are cost rather than correctness:

- The traceback keeps a decision per cell while the matrix fits `stored_cells`, and recurses in
  linear space once it does not. Both return an optimal path for the same score.
- On the device, a pair short enough for one block's shared-memory carry takes the banded strip
  sweep, and a taller one tiles over global memory. The bound comes from what the card reports.

A batch on the device goes out as one launch for every pair both bounds admit; the rest are served
one by one, so a single oversized pair never sinks the batch it arrived in.
"""

from .alignment import (
    AffineGapCosts,
    AlignmentMode,
    GappedAlignment,
    DEFAULT_LEAF_CELLS,
    DEVICE_STORED_CELLS,
    Layer,
    Space,
    SweepHalf,
    band_length,
    device_align,
    device_alignments,
    device_score,
    device_scores,
    expand_path,
    score_path,
    serial_align,
    serial_hirschberg,
    serving_space,
)
from .common import (
    spread,
    Device,
    DeviceScope,
    MAX_ALPHABET_SIZE,
    OffsetDType,
    Placement,
    SubstitutionDType,
    SymbolDType,
    UNKNOWN_SYMBOL,
    translate,
    uniform_matrix,
)
from .errors import AlignmentError, ErrorKind
from .gap_affine import (
    DEFAULT_MAX_MEMORY,
    FREE_START,
    HISTORY_LIMIT,
    EndsFree,
    Penalties,
    cigar_of,
    SearchSpace,
    solve,
    wavefront_align,
    wavefront_penalties,
    wavefront_score,
)
from .lanes import (
    TABLE_ENTRIES,
    CodeTexts,
    LaneCosts,
    StringTexts,
    lane_alignments,
    lane_distances,
    lane_local_scores,
)
from .modes import Alignment, Anchor, Band, Costs, Mode
from .scored import ANYWHERE, FROM_EDGE, FROM_ORIGIN, swept_cells
from .band_groups import banded_scores
from .score_groups import grouped_scores
from .substitutions import SubstitutionLookup
from .vector_score import optimal_band, reach_back, vector_align, vector_score


from std.memory import bitcast
from std.math import gcd

comptime STORED_CELL_BYTES = 12
"""Bytes the host traceback keeps a stored cell in: three `int32` layers. The device packs a nibble per cell
and is capped again by `DEVICE_STORED_CELLS`."""


def cells_within(max_memory: Int) -> Int:
    """Cells a traceback may store within `max_memory` bytes, above which it recurses in linear space."""
    return max(max_memory, 0) // STORED_CELL_BYTES


# region Scoring

comptime DNA_ALPHABET = "ACGT"
"""The four bases. An `N` or a soft-masked lowercase base must be added to an alphabet explicitly."""

comptime DEFAULT_MATCH = 2
"""Minimap2's match score, `-A2`."""
comptime DEFAULT_MISMATCH = -4
"""Minimap2's mismatch score, `-B4`."""
comptime DEFAULT_GAP_OPENING = -4
"""Minimap2's gap opening, `-O4`: a gap of `k` letters scores `-(4 + 2k)`."""
comptime DEFAULT_GAP_EXTENSION = -2
"""Minimap2's gap extension, `-E2`."""


struct Scoring(Copyable, Movable):
    """An alphabet, the substitution table it indexes, and the affine gap model, which travel together.

    Built through `dna`, `edit_distance`, `uniform` or `tabulated` rather than field by field, so a
    table whose shape disagrees with its alphabet cannot be expressed. A gap of `k` letters scores
    `opening + k extension`, both scores zero or less, as `Costs` counts a gap's cost.
    """

    var alphabet: String
    """The letters a sequence may hold, in the order the table is indexed by."""
    var substitutions: List[Scalar[SubstitutionDType]]
    """Row-major, one row per letter of `alphabet`."""
    var gaps: AffineGapCosts
    """The gap scores as the kernels take them: a gap's first letter scores `open`, each further one `extend`."""

    def __init__(
        out self, var alphabet: String, var substitutions: List[Scalar[SubstitutionDType]], gaps: AffineGapCosts
    ):
        """Trusts its arguments; the factories are the checked way in."""
        self.alphabet = alphabet^
        self.substitutions = substitutions^
        self.gaps = gaps

    @staticmethod
    def dna() raises AlignmentError -> Self:
        """Minimap2's scoring over `ACGT`: match 2, mismatch -4, a gap of `k` letters `-(4 + 2k)`."""
        return Self.uniform(DEFAULT_MATCH, DEFAULT_MISMATCH, DEFAULT_GAP_OPENING, DEFAULT_GAP_EXTENSION)

    @staticmethod
    def edit_distance(alphabet: String = String(DNA_ALPHABET)) raises AlignmentError -> Self:
        """Unit costs, under which a global score is the negated Levenshtein distance."""
        return Self.uniform(0, -1, 0, -1, alphabet)

    @staticmethod
    def uniform(
        match_score: Int,
        mismatch_score: Int,
        opening: Int = DEFAULT_GAP_OPENING,
        extension: Int = DEFAULT_GAP_EXTENSION,
        alphabet: String = String(DNA_ALPHABET),
    ) raises AlignmentError -> Self:
        """One score for equal letters and one for unequal, over any alphabet."""
        var size = alphabet.byte_length()
        if size == 0 or size > MAX_ALPHABET_SIZE:
            raise AlignmentError(ErrorKind.ALPHABET_TOO_LARGE, String(size, " letters"))
        return Self(
            alphabet,
            uniform_matrix(size, match_score, mismatch_score),
            gap_scores(opening, extension),
        )

    @staticmethod
    def tabulated(
        alphabet: String,
        var substitutions: List[Scalar[SubstitutionDType]],
        opening: Int = DEFAULT_GAP_OPENING,
        extension: Int = DEFAULT_GAP_EXTENSION,
    ) raises AlignmentError -> Self:
        """A caller's own table, refused unless it is square in the alphabet that indexes it."""
        var size = alphabet.byte_length()
        if size == 0 or size > MAX_ALPHABET_SIZE:
            raise AlignmentError(ErrorKind.ALPHABET_TOO_LARGE, String(size, " letters"))
        if len(substitutions) != size * size:
            raise AlignmentError(ErrorKind.INVALID_SCORING, String(len(substitutions), " cells for ", size, " letters"))
        return Self(alphabet, substitutions^, gap_scores(opening, extension))

    def alphabet_size(self) -> Int:
        """Letters the table is indexed by, which is its stride."""
        return self.alphabet.byte_length()


def gap_scores(opening: Int, extension: Int) raises AlignmentError -> AffineGapCosts:
    """A gap scoring `opening + k extension` over `k` letters, as the kernels take it: its first letter
    `opening + extension`, each further one `extension`."""
    if opening > 0:
        raise AlignmentError(ErrorKind.INVALID_SCORING, "a rewarded gap opening")
    return AffineGapCosts.checked(Int32(opening + extension), Int32(extension))


# endregion Scoring

# region Batch Tape


@fieldwise_init
struct BatchTape(Movable):
    """Both sides of a batch on one tape, which is the shape a kernel can take a pointer to."""

    var sequences: List[Scalar[SymbolDType]]
    """Every sequence concatenated, first and second of each pair alternating."""
    var offsets: List[Scalar[OffsetDType]]
    """Where each sequence begins, so a block can find its own pair."""


def pack_batch(
    firsts: List[String], seconds: List[String], indices: List[Int], alphabet: String, threads: Int = 1
) raises AlignmentError -> BatchTape:
    """The named pairs concatenated onto one tape, with offsets marking where each sequence begins.

    Where each sequence goes is worked out first, and then each of `threads` writes its stretch of the
    pairs' codes straight onto the tape: one thread translating a sequence at a time into a list of
    its own took three times as long as the kernel that scored 500,000 short reads. A letter outside
    the alphabet raises as `translate` does, for the first pair in order holding one."""
    var codes_by_byte = Array[UInt8, 256](fill=UNKNOWN_SYMBOL)
    var alphabet_bytes = alphabet.as_bytes()
    for index in range(len(alphabet_bytes)):
        codes_by_byte[Int(alphabet_bytes[index])] = UInt8(index)
    var pairs = len(indices)
    var offsets = List[Scalar[OffsetDType]](capacity=2 * pairs + 1)
    offsets.append(0)
    var length = 0
    for index in indices:
        length += firsts[index].byte_length()
        offsets.append(Scalar[OffsetDType](length))
        length += seconds[index].byte_length()
        offsets.append(Scalar[OffsetDType](length))
    var sequences = List[Scalar[SymbolDType]](capacity=length)
    sequences.resize(unsafe_uninit_length=length)
    var failed = List[Bool](length=pairs, fill=False)
    var tape = sequences.unsafe_ptr()
    var places = offsets.unsafe_ptr()
    var flags = failed.unsafe_ptr()
    var stretches = max(min(threads, pairs), 1)

    def encode(
        stretch: Int,
    ) {
        imm firsts,
        imm seconds,
        imm indices,
        imm codes_by_byte,
        imm pairs,
        imm stretches,
        imm tape,
        imm places,
        imm flags,
    }:
        """The codes of stretch `stretch`'s pairs onto the tape, flagging a pair holding a letter
        outside the alphabet."""
        for slot in range(pairs * stretch // stretches, pairs * (stretch + 1) // stretches):
            var index = indices[slot]
            var known = encoded_into(
                firsts[index], codes_by_byte, tape.unsafe_offset(Int(places[unsafe_offset=2 * slot]))
            )
            known = (
                encoded_into(seconds[index], codes_by_byte, tape.unsafe_offset(Int(places[unsafe_offset=2 * slot + 1])))
                and known
            )
            if not known:
                flags[unsafe_offset=slot] = True

    spread(encode, stretches, stretches)
    for slot in range(pairs):
        if failed[slot]:
            # The pair's own translation raises the error the serial packing raised.
            _ = translate(firsts[indices[slot]], alphabet)
            _ = translate(seconds[indices[slot]], alphabet)
    return BatchTape(sequences^, offsets^)


@always_inline
def encoded_into(text: String, codes_by_byte: Array[UInt8, 256], target: MutPointer[Scalar[SymbolDType], _]) -> Bool:
    """Writes `text`'s codes from `target` on, and whether every letter had one."""
    var bytes = text.unsafe_ptr()
    var unknown = False
    for position in range(text.byte_length()):
        var code = codes_by_byte[Int(bytes[unsafe_offset=position])]
        unknown = unknown or code == UNKNOWN_SYMBOL
        target[unsafe_offset=position] = Scalar[SymbolDType](code)
    return not unknown


def paired_length(firsts: List[String], seconds: List[String]) raises AlignmentError -> Int:
    """The number of pairs, refusing two sides that do not line up."""
    if len(firsts) != len(seconds):
        raise AlignmentError(ErrorKind.LENGTH_MISMATCH, String(len(firsts), " firsts, ", len(seconds), " seconds"))
    return len(firsts)


# endregion Batch Tape

# region Linear-Space Host Traceback


def global_linear(
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    scoring: Scoring,
    vectorized: Bool = True,
) raises -> GappedAlignment:
    """Global alignment in linear space, splitting rows and joining halves Myers-Miller style.

    `vectorized` lets the sweeps run sixteen cells at a time, which computes the same rows; off, every
    sweep runs cell by cell, as a reference to check it against.
    """
    var path_columns = List[Int32](length=len(first) + 1, fill=Int32(0))
    var path_layers = List[Layer](length=len(first) + 1, fill=Layer.ALIGNING)
    serial_hirschberg(
        first,
        second,
        0,
        len(first),
        0,
        len(second),
        scoring.substitutions,
        scoring.alphabet_size(),
        scoring.gaps,
        DEFAULT_LEAF_CELLS,
        path_columns,
        path_layers,
        vectorized,
    )
    var score = score_path(
        first,
        second,
        path_columns,
        path_layers,
        scoring.substitutions,
        scoring.alphabet_size(),
        scoring.gaps,
        len(first),
    )
    var expanded = expand_path(
        first, second, path_columns, path_layers, scoring.alphabet, AlignmentMode.GLOBAL, 0, len(first)
    )
    return GappedAlignment(score, expanded[0], expanded[1])


# endregion Linear-Space Host Traceback

# region Routing


def host_score[
    mode: AlignmentMode
](first: ImmSpan[Scalar[SymbolDType], _], second: ImmSpan[Scalar[SymbolDType], _], scoring: Scoring) -> Int32:
    """The optimal score on the host: by wavefront for a global one under a table of one match and
    one mismatch score while it is cheaper (see `gap_affine`), else by full sweep, sixteen cells at a
    time under any table (see `vector_score`)."""
    var codes_first = List[UInt8](first)
    var codes_second = List[UInt8](second)
    comptime if mode == AlignmentMode.GLOBAL:
        # A table of one match and one mismatch score has a wavefront, whose work grows with the
        # score rather than the matrix; it hands back a pair a full sweep would serve sooner.
        var penalties = wavefront_penalties(
            scoring.substitutions, scoring.alphabet_size(), Int(scoring.gaps.open), Int(scoring.gaps.extend)
        )
        if penalties:
            var found = wavefront_score(codes_first, codes_second, penalties.value())
            if found:
                return Int32(found.value())
    var lookup = SubstitutionLookup(scoring.substitutions, scoring.alphabet_size())
    return vector_score[mode](codes_first, codes_second, lookup, scoring.gaps)


def align_on_host[
    mode: AlignmentMode
](
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    scoring: Scoring,
    stored_cells: Int,
) raises -> GappedAlignment:
    """One pair on the host, in memory that grows with the alignment rather than the matrix.

    A local alignment is aligned by its span (see `local_by_span`). A global one under a table of one
    match and one mismatch score, whose gap costs a wavefront can grow by, goes to the two-ended
    wavefront, traced back through its own fronts and split where they would grow too large (see
    `gap_affine.wavefront_align`). Under any other table its matrix is swept sixteen cells at a time,
    storing only the band of diagonals every optimal path stays on, which its score bounds (see
    `optimal_band`, `vector_align`); past `stored_cells` the alignment recurses in linear space.
    """
    comptime if mode == AlignmentMode.LOCAL:
        return local_by_span(first, second, scoring, stored_cells)
    return global_on_host(first, second, scoring, stored_cells)


def global_on_host(
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    scoring: Scoring,
    stored_cells: Int,
    known: Optional[Int] = None,
) raises -> GappedAlignment:
    """`align_on_host`'s global alignment; with the score `known` beforehand, as a local alignment's span
    knows it, the sweep that would find it to bound the band is skipped."""
    if len(first) == 0 or len(second) == 0:
        return serial_align[AlignmentMode.GLOBAL](
            first, second, scoring.substitutions, scoring.alphabet_size(), scoring.gaps, scoring.alphabet
        )
    var codes_first = List[UInt8](first)
    var codes_second = List[UInt8](second)
    # The wavefront's work grows with the score rather than the matrix, and it never hands a pair back.
    var penalties = wavefront_penalties(
        scoring.substitutions, scoring.alphabet_size(), Int(scoring.gaps.open), Int(scoring.gaps.extend)
    )
    if penalties:
        var traced = wavefront_align(codes_first, codes_second, penalties.value(), scoring.alphabet)
        return GappedAlignment(Int32(traced[0]), traced[1], traced[2])
    var lookup = SubstitutionLookup(scoring.substitutions, scoring.alphabet_size())
    # Not `known.or_else(...)`, which would sweep for the score whether it is known or not.
    var best: Int
    if known:
        best = known.value()
    else:
        best = Int(vector_score[AlignmentMode.GLOBAL](codes_first, codes_second, lookup, scoring.gaps))
    var band = optimal_band(len(first), len(second), lookup.best, scoring.gaps, best)
    var width = min(band[1], len(second)) - max(band[0], -len(first)) + 1
    if (len(first) + 1) * width <= stored_cells:
        return vector_align(
            codes_first,
            codes_second,
            lookup,
            scoring.gaps,
            scoring.substitutions,
            scoring.alphabet_size(),
            scoring.alphabet,
            band[0],
            band[1],
        )
    return global_linear(first, second, scoring)


def local_by_span(
    first: ImmSpan[Scalar[SymbolDType], _], second: ImmSpan[Scalar[SymbolDType], _], scoring: Scoring, stored_cells: Int
) raises -> GappedAlignment:
    """A local alignment by its span, as SSW finds one: the forward sweep's best cell ends it, in lanes
    of 16 bits while the scores fit (see `tabulated_end`), a sweep back from there to where it earns its
    score starts it (see `reach_back`), and the letters between are aligned globally, which scores the
    same. Only the forward sweep covers the matrix; the rest grows with the alignment."""
    var codes_first = List[UInt8](first)
    var codes_second = List[UInt8](second)
    var end = tabulated_end[ANYWHERE](Span(codes_first), Span(codes_second), scoring, EndsFree(), True)
    var best = end[0]
    if best <= 0:
        return GappedAlignment(0, String(), String())
    var lookup = SubstitutionLookup(scoring.substitutions, scoring.alphabet_size())
    var start = reach_back(codes_first, codes_second, end[1], end[2], lookup, scoring.gaps, Int32(best))
    var inner = global_on_host(first[start[0] : end[1]], second[start[1] : end[2]], scoring, stored_cells, best)
    # A global alignment of the span scores what the local one does, and is one.
    return inner^


def align_on_device[
    mode: AlignmentMode
](
    scope: DeviceScope,
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    scoring: Scoring,
    stored_cells: Int,
    placement: Placement,
) raises -> GappedAlignment:
    """One pair on the device, on the sweep its height and its matrix can afford.

    Both bounds are real and independent: the stored kernel indexes its carry by the first
    sequence, so a tall pair fails it even when the whole matrix would fit.
    """
    var stored = len(first) * len(second) <= min(stored_cells, DEVICE_STORED_CELLS)
    if not stored or serving_space(len(first), band_length(scope.specs)) == Space.TILED:
        return device_align[mode](
            scope,
            first,
            second,
            scoring.substitutions,
            scoring.alphabet_size(),
            scoring.gaps,
            scoring.alphabet,
            DEFAULT_LEAF_CELLS,
            placement,
        )
    # No single-pair device entry exists for the stored traceback, so this is a batch of one.
    var sequences = List[Scalar[SymbolDType]](capacity=len(first) + len(second))
    sequences.extend(first)
    sequences.extend(second)
    var offsets: List[Scalar[OffsetDType]] = [0, Scalar[OffsetDType](len(first)), Scalar[OffsetDType](len(sequences))]
    var aligned = device_alignments[mode](
        scope, sequences, offsets, scoring.substitutions, scoring.alphabet, scoring.gaps
    )
    return aligned.pop(0)


def score_on_device[
    mode: AlignmentMode
](
    scope: DeviceScope,
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    scoring: Scoring,
) raises -> Int32:
    """One pair scored on the device, banded when one block's carry can index it and tiled otherwise."""
    if serving_space(len(first), band_length(scope.specs)) == Space.TILED:
        return device_score[mode](scope, first, second, scoring.substitutions, scoring.alphabet_size(), scoring.gaps)
    var sequences = List[Scalar[SymbolDType]](capacity=len(first) + len(second))
    sequences.extend(first)
    sequences.extend(second)
    var offsets: List[Scalar[OffsetDType]] = [0, Scalar[OffsetDType](len(first)), Scalar[OffsetDType](len(sequences))]
    return device_scores[mode](scope, sequences, offsets, scoring.substitutions, scoring.alphabet_size(), scoring.gaps)[
        0
    ]


# endregion Routing

# region Every Mode


def gap_costs(scoring: Scoring) -> Costs:
    """A `Scoring`'s gaps as `Costs` price them, the sweep's own terms: a gap of `k` letters scoring
    `opening + k extension` costs minus that. Its mismatch stands for nothing: a table scores pairs."""
    return Costs(1, Int(scoring.gaps.extend - scoring.gaps.open), Int(-scoring.gaps.extend), -1, 0)


def table_fits[kind: Int](scoring: Scoring, rows: Int, columns: Int) -> Bool:
    """Whether every score of a sweep under the table fits 16 bits, as `narrow_enough` asks of costs:
    the best pair over the shorter sequence, and the dearest move over both for a sweep from the edges."""
    var most = 0
    var least = 0
    for value in scoring.substitutions:
        most = max(most, Int(value))
        least = min(least, Int(value))
    var dearest = max(-least, Int(-scoring.gaps.open))
    var fits = most * (min(rows, columns) + 1) < 32000 and dearest < 4000
    comptime if kind != ANYWHERE:
        fits = fits and dearest * (rows + columns + 1) < 8000
    return fits


def tabulated_end[
    kind: Int
](
    first: Span[UInt8, _], second: Span[UInt8, _], scoring: Scoring, ends: EndsFree, highest: Bool, zdrop: Int = -1
) -> Tuple[Int, Int, Int, Bool]:
    """The best score under `scoring` of an alignment `kind` allows, where it ends, and whether a Z-drop
    gave the sweep up, by the sweep
    `scored.swept_cells` runs, each pair's score read from the table: `first` and `second` are codes
    into the alphabet. Lanes along the shorter sequence, 16 bits while the scores fit."""
    var gaps = gap_costs(scoring)
    var size = scoring.alphabet_size()
    var transposed = len(second) < len(first)
    var unused = List[Int32]()
    if table_fits[kind](scoring, len(first), len(second)):
        var table = List[Int16](capacity=len(scoring.substitutions))
        for value in scoring.substitutions:
            table.append(Int16(value))
        if transposed:
            return swept_cells[1, DType.int16, 32, True, kind, False, True](
                first, second, gaps, 0, ends, highest, table, size, zdrop, unused
            )
        return swept_cells[1, DType.int16, 32, False, kind, False, True](
            first, second, gaps, 0, ends, highest, table, size, zdrop, unused
        )
    var table = List[Int32](capacity=len(scoring.substitutions))
    for value in scoring.substitutions:
        table.append(Int32(value))
    if transposed:
        return swept_cells[1, DType.int32, 16, True, kind, False, True](
            first, second, gaps, 0, ends, highest, table, size, zdrop, unused
        )
    return swept_cells[1, DType.int32, 16, False, kind, False, True](
        first, second, gaps, 0, ends, highest, table, size, zdrop, unused
    )


def ends_of(mode: Mode, columns: Int, rows: Int) -> EndsFree:
    """A mode's free letters, none past its sequence's length."""
    return EndsFree(
        min(mode.reference_start, columns),
        min(mode.reference_end, columns),
        min(mode.query_start, rows),
        min(mode.query_end, rows),
    )


def extension_span(
    first: List[Scalar[SymbolDType]], second: List[Scalar[SymbolDType]], scoring: Scoring, mode: Mode
) -> Tuple[Int, Int, Int, Int, Int, Bool]:
    """An extension's best stop under `scoring` and the span it covers, as `mode_span` gives them, and
    whether the Z-drop gave the sweep up: by a sweep from its anchor for the end as late as an equally
    good one allows."""
    var columns = len(first)
    var rows = len(second)
    if mode.anchor == Anchor.END:
        var back_first = List[Scalar[SymbolDType]](capacity=columns)
        for index in range(columns - 1, -1, -1):
            back_first.append(first[index])
        var back_second = List[Scalar[SymbolDType]](capacity=rows)
        for index in range(rows - 1, -1, -1):
            back_second.append(second[index])
        var found = tabulated_end[FROM_ORIGIN](back_first, back_second, scoring, EndsFree(), True, mode.zdrop)
        return (found[0], columns - found[1], rows - found[2], columns, rows, found[3])
    var found = tabulated_end[FROM_ORIGIN](first, second, scoring, EndsFree(), True, mode.zdrop)
    return (found[0], 0, 0, found[1], found[2], found[3])


def mode_span(
    first: List[Scalar[SymbolDType]], second: List[Scalar[SymbolDType]], scoring: Scoring, mode: Mode
) -> Tuple[Int, Int, Int, Int, Int]:
    """The best score under `scoring` with free ends or as an extension, and the span it covers, as
    `Costs` place it under `Ties.LEFT` (see `gap_affine.free_ends_alignment`): free ends by a sweep
    from the edges for the end on the highest diagonal, then one back from it for the start on the
    highest too; an extension by a sweep from its anchor for the end as late as an equally good one
    allows. The score, then the start's and the end's letters of each sequence."""
    var columns = len(first)
    var rows = len(second)
    if mode.kind == Mode.EXTENSION:
        var stop = extension_span(first, second, scoring, mode)
        # The end bonus prefers the best extension reaching the query's far end, unless the Z-drop gave up.
        if mode.end_bonus > 0 and not stop[5]:
            var reaching = mode_span(first, second, scoring, mode.reaching_end())
            if reaching[0] + mode.end_bonus > stop[0]:
                return reaching
        return (stop[0], stop[1], stop[2], stop[3], stop[4])
    var ends = ends_of(mode, columns, rows)
    var forward = tabulated_end[FROM_EDGE](first, second, scoring, ends, True)
    var end_column = forward[1]
    var end_row = forward[2]
    var head = List[Scalar[SymbolDType]](capacity=end_column)
    for index in range(end_column - 1, -1, -1):
        head.append(first[index])
    var lead = List[Scalar[SymbolDType]](capacity=end_row)
    for index in range(end_row - 1, -1, -1):
        lead.append(second[index])
    var back = tabulated_end[FROM_EDGE](head, lead, scoring, EndsFree(0, ends.first_begin, 0, ends.second_begin), False)
    return (forward[0], end_column - back[1], end_row - back[2], end_column, end_row)


def as_alignment(gapped: GappedAlignment, first: String, second: String, whole: Bool, eqx: Bool) -> Alignment:
    """Gotoh's gapped rows as an `Alignment`: spanning both sequences when `whole`, else a local
    alignment's, placed where its letters lie in each, the last place: any place both lie aligns the
    same pairs of letters for the same score."""
    var score = Int(gapped.score)
    var cigar = gapped.cigar(eqx)
    if whole:
        return Alignment(-score, score, cigar, 0, first.byte_length(), 0, second.byte_length())
    var part = gapped.first_gapped.replace("-", "")
    var piece = gapped.second_gapped.replace("-", "")
    var first_start = first.rfind(part) if part.byte_length() > 0 else 0
    var second_start = second.rfind(piece) if piece.byte_length() > 0 else 0
    return Alignment(
        -score,
        score,
        cigar,
        first_start,
        first_start + part.byte_length(),
        second_start,
        second_start + piece.byte_length(),
    )


def scoring_alignment(
    first: String,
    second: String,
    scoring: Scoring,
    mode: Mode,
    placement: Optional[Placement],
    stored_cells: Int,
    eqx: Bool,
) raises -> Alignment:
    """An optimal alignment under `scoring` as `mode` asks, as an `Alignment`: its CIGAR, its spans and
    its score, its cost minus the score. Global and local alignments take Gotoh's sweeps on either
    device (see `align_with`), and a local one's spans are found where its letters lie; free ends and
    extensions, on the host, find their span by sweep (see `mode_span`) and align it globally."""
    var space = SearchSpace()
    return scoring_alignment(first, second, scoring, mode, placement, stored_cells, eqx, space)


def scoring_alignment(
    first: String,
    second: String,
    scoring: Scoring,
    mode: Mode,
    placement: Optional[Placement],
    stored_cells: Int,
    eqx: Bool,
    mut space: SearchSpace,
) raises -> Alignment:
    """`scoring_alignment` through `space`'s searches, which a batch's worker keeps from pair to pair."""
    if mode.match_score > 0:
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a Scoring's table holds what a match earns")
    var on_host = not placement or placement.value().device != Device.GPU
    if on_host and mode.is_global() and first.byte_length() > 0 and second.byte_length() > 0:
        # A table of one match and one mismatch score: the wavefront's own moves spell the CIGAR, with
        # no rows between, the same alignment `align_on_host` would give (see `wavefront_align`).
        var penalties = wavefront_penalties(
            scoring.substitutions, scoring.alphabet_size(), Int(scoring.gaps.open), Int(scoring.gaps.extend)
        )
        if penalties:
            var codes_first = translate(first, scoring.alphabet)
            var codes_second = translate(second, scoring.alphabet)
            var moves = List[UInt8](capacity=len(codes_first) + len(codes_second))
            var cost = solve[1](
                space.forward,
                space.backward,
                Span(codes_first),
                Span(codes_second),
                penalties.value(),
                FREE_START,
                FREE_START,
                HISTORY_LIMIT,
                moves,
            )
            var score = penalties.value().score(cost, len(codes_first) + len(codes_second))
            var cigar = cigar_of(first, second, moves^, cost, penalties.value(), eqx)
            return Alignment(-score, score, cigar, 0, len(codes_first), 0, len(codes_second))
    if mode.kind == Mode.SMITH_WATERMAN or mode.is_global():
        var gapped: GappedAlignment
        if mode.kind == Mode.SMITH_WATERMAN:
            gapped = align_with[AlignmentMode.LOCAL](first, second, scoring, placement, stored_cells)
        else:
            gapped = align_with[AlignmentMode.GLOBAL](first, second, scoring, placement, stored_cells)
        return as_alignment(gapped, first, second, mode.is_global(), eqx)
    if placement and placement.value().device == Device.GPU:
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "on the GPU a Scoring aligns globally or locally")
    var codes_first = translate(first, scoring.alphabet)
    var codes_second = translate(second, scoring.alphabet)
    var span = mode_span(codes_first, codes_second, scoring, mode)
    var start_column = span[1]
    var start_row = span[2]
    var end_column = span[3]
    var end_row = span[4]
    var inner = align_on_host[AlignmentMode.GLOBAL](
        Span(codes_first)[start_column:end_column], Span(codes_second)[start_row:end_row], scoring, stored_cells
    )
    return Alignment(-span[0], span[0], inner.cigar(eqx), start_column, end_column, start_row, end_row)


def scoring_score(
    first: String, second: String, scoring: Scoring, mode: Mode, placement: Optional[Placement]
) raises -> Int:
    """The best score under `scoring` as `mode` asks, with no alignment traced: global and local on
    either device (see `score_with`), free ends and extensions by their sweep on the host."""
    if mode.match_score > 0:
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a Scoring's table holds what a match earns")
    if mode.kind == Mode.SMITH_WATERMAN:
        return Int(score_with[AlignmentMode.LOCAL](first, second, scoring, placement))
    if mode.is_global():
        return Int(score_with[AlignmentMode.GLOBAL](first, second, scoring, placement))
    if placement and placement.value().device == Device.GPU:
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "on the GPU a Scoring aligns globally or locally")
    return mode_span(translate(first, scoring.alphabet), translate(second, scoring.alphabet), scoring, mode)[0]


# endregion Every Mode

# region Entry Points


def score_with[
    mode: AlignmentMode
](first: String, second: String, scoring: Scoring, placement: Optional[Placement] = None) raises -> Int32:
    """The optimal score alone, in two rows of memory on either device.

    On the host, a global score under a table of one match and one mismatch score runs the
    wavefront first (see `gap_affine`), and the full sweep only when the wavefront gives up; under
    such a table the sweep runs sixteen cells at a time (see `vector_score`).
    """
    var resolved = placement.or_else(Placement.default())
    var encoded_first = translate(first, scoring.alphabet)
    var encoded_second = translate(second, scoring.alphabet)
    if resolved.device == Device.GPU:
        return score_on_device[mode](DeviceScope(resolved.gpu_id), encoded_first, encoded_second, scoring)
    return host_score[mode](encoded_first, encoded_second, scoring)


def align_with[
    mode: AlignmentMode
](
    first: String,
    second: String,
    scoring: Scoring,
    placement: Optional[Placement] = None,
    stored_cells: Int = cells_within(DEFAULT_MAX_MEMORY),
) raises -> GappedAlignment:
    """The optimal score and the two gapped strings that realize it.

    A global alignment spans both sequences; a local one returns only the best-scoring window,
    trimmed at both ends.
    """
    var resolved = placement.or_else(Placement.default())
    var encoded_first = translate(first, scoring.alphabet)
    var encoded_second = translate(second, scoring.alphabet)
    if resolved.device == Device.GPU:
        return align_on_device[mode](
            DeviceScope(resolved.gpu_id), encoded_first, encoded_second, scoring, stored_cells, resolved
        )
    return align_on_host[mode](encoded_first, encoded_second, scoring, stored_cells)


comptime CHUNKS_PER_THREAD = 8
"""Chunks a host batch is cut into per thread, so a thread that draws long pairs does not hold up
the rest while the others sit idle."""


def chunk_count(pairs: Int, threads: Int) -> Int:
    """How many contiguous chunks a host batch of `pairs` is cut into for `threads` threads."""
    return min(pairs, max(threads, 1) * CHUNKS_PER_THREAD)


def unknown_letters(firsts: List[String], seconds: List[String], alphabet: String, workers: Int) -> List[Bool]:
    """Which pairs hold a letter outside `alphabet`, checked over `workers` threads: the lanes leave those
    for their own calls to raise as they do."""
    var pairs = len(firsts)
    var held = Array[Bool, 256](fill=False)
    for letter in alphabet.as_bytes():
        held[Int(letter)] = True
    var refused = List[Bool](length=pairs, fill=False)
    var refused_ptr = refused.unsafe_ptr()

    def refuse(stretch: Int) {imm firsts, imm seconds, imm held, imm pairs, imm workers, imm refused_ptr}:
        """Marks stretch `stretch`'s pairs holding a letter outside the alphabet."""
        for index in range(pairs * stretch // workers, pairs * (stretch + 1) // workers):
            var known = True
            for letter in firsts[index].as_bytes():
                known = known and held[Int(letter)]
            for letter in seconds[index].as_bytes():
                known = known and held[Int(letter)]
            if not known:
                refused_ptr[unsafe_offset=index] = True

    spread(refuse, workers, workers)
    return refused^


def laned_scores(
    firsts: List[String],
    seconds: List[String],
    alphabet: String,
    costs: LaneCosts,
    penalties: Penalties,
    threads: Int,
    scores_out: MutPointer[Int32, _],
) -> List[Bool]:
    """The global scores of every pair the lanes take into `scores_out`, and which they were: those whose letters
    `alphabet` holds, a pair holding another left for its own call to raise as it does."""
    var pairs = len(firsts)
    var workers = max(threads, 1)
    # A pair holding a letter outside the alphabet is marked settled beforehand, so the lanes leave it.
    var refused = unknown_letters(firsts, seconds, alphabet, workers)
    var settled = refused.copy()
    var found = List[Optional[Int]](length=pairs, fill=None)
    _ = lane_distances(
        pairs,
        StringTexts(firsts.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]()),
        StringTexts(seconds.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]()),
        costs,
        Band(),
        Int.MAX,
        workers,
        found.unsafe_ptr(),
        settled.unsafe_ptr(),
    )
    for index in range(pairs):
        if refused[index]:
            settled[index] = False
        elif settled[index]:
            var letters = firsts[index].byte_length() + seconds[index].byte_length()
            scores_out[unsafe_offset=index] = Int32(penalties.score(found[index].value(), letters))
    return settled^


def laned_alignments(
    firsts: List[String],
    seconds: List[String],
    scoring: Scoring,
    eqx: Bool,
    threads: Int,
    budget: Int,
    alignments_out: MutPointer[Alignment, _],
) -> List[Bool]:
    """The global alignments of every pair the lanes take into `alignments_out`, and which they were, under
    a table of one match and one mismatch score: the alignment `scoring_alignment` gives, the wavefront's
    path by `Ties.LEFT`, traced from the flags a group's band kept within `budget` bytes (see
    `lanes.lane_alignments`). A pair holding a letter outside the alphabet, with an empty side, or one the
    searches might split for memory, which would follow the tie rule within each piece, is left for its own
    call."""
    var pairs = len(firsts)
    var workers = max(threads, 1)
    var penalties = wavefront_penalties(
        scoring.substitutions, scoring.alphabet_size(), Int(scoring.gaps.open), Int(scoring.gaps.extend)
    )
    if not penalties:
        return List[Bool](length=pairs, fill=False)
    var found = penalties.value()
    var costs = LaneCosts.one_piece(found.mismatch, found.opening, found.extension, found.opening, found.extension)
    var refused = unknown_letters(firsts, seconds, scoring.alphabet, workers)
    var settled = refused.copy()
    var laned = List[Optional[Int]](length=pairs, fill=None)
    var paths = List[List[UInt8]](capacity=pairs)
    for _ in range(pairs):
        paths.append(List[UInt8]())
    _ = lane_alignments(
        pairs,
        StringTexts(firsts.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]()),
        StringTexts(seconds.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]()),
        costs,
        Band(),
        Int.MAX,
        True,
        workers,
        budget,
        laned.unsafe_ptr(),
        paths.unsafe_ptr(),
        settled.unsafe_ptr(),
    )
    var settled_ptr = settled.unsafe_ptr()
    var laned_ptr = laned.unsafe_ptr()
    var path_ptr = paths.unsafe_ptr()
    var refused_ptr = refused.unsafe_ptr()

    def spell(
        stretch: Int,
    ) {
        imm firsts,
        imm seconds,
        imm found,
        imm eqx,
        imm pairs,
        imm workers,
        imm settled_ptr,
        imm laned_ptr,
        imm path_ptr,
        imm refused_ptr,
        imm alignments_out,
    }:
        """Spells the CIGARs of stretch `stretch`'s pairs the lanes traced, and unsettles the rest."""
        for index in range(pairs * stretch // workers, pairs * (stretch + 1) // workers):
            if refused_ptr[unsafe_offset=index] or not settled_ptr[unsafe_offset=index]:
                settled_ptr[unsafe_offset=index] = False
                continue
            var cost = laned_ptr[unsafe_offset=index].value()
            var columns = firsts[index].byte_length()
            var rows = seconds[index].byte_length()
            # The searches split a pair whose fronts pass `HISTORY_LIMIT`, which they cannot below this.
            if 2 * (cost + 1) * (columns + rows + 1) > HISTORY_LIMIT:
                settled_ptr[unsafe_offset=index] = False
                continue
            var moves = List[UInt8]()
            swap(moves, path_ptr[unsafe_offset=index])
            var score = found.score(cost, columns + rows)
            var cigar = cigar_of(firsts[index], seconds[index], moves^, cost, found, eqx)
            alignments_out[unsafe_offset=index] = Alignment(-score, score, cigar^, 0, columns, 0, rows)

    spread(spell, workers, workers)
    return settled^


def uniform_scores(scoring: Scoring) -> Optional[Tuple[Int, Int]]:
    """The match and the mismatch score of a table of one each, if it is one."""
    var size = scoring.alphabet_size()
    if size < 2:
        return None
    var hit = Int(scoring.substitutions[0])
    var mismatch = Int(scoring.substitutions[1])
    for row in range(size):
        for column in range(size):
            if Int(scoring.substitutions[row * size + column]) != (hit if row == column else mismatch):
                return None
    return (hit, mismatch)


def laned_local_scores(
    firsts: List[String], seconds: List[String], scoring: Scoring, threads: Int, scores_out: MutPointer[Int32, _]
) -> List[Bool]:
    """The local scores of every pair the lanes take into `scores_out`, and which they were, under a table
    of one match and one mismatch score (see `lanes.lane_local_scores`): those whose letters the alphabet
    holds, padded past their ends with two bytes it does not."""
    var pairs = len(firsts)
    var workers = max(threads, 1)
    var uniform = uniform_scores(scoring)
    if not uniform:
        return List[Bool](length=pairs, fill=False)
    var held = Array[Bool, 256](fill=False)
    for letter in scoring.alphabet.as_bytes():
        held[Int(letter)] = True
    var pads = List[UInt8]()
    for byte in range(256):
        if not held[byte] and len(pads) < 2:
            pads.append(UInt8(byte))
    if len(pads) < 2:
        return List[Bool](length=pairs, fill=False)
    var refused = unknown_letters(firsts, seconds, scoring.alphabet, workers)
    var settled = refused.copy()
    _ = lane_local_scores(
        pairs,
        StringTexts(firsts.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]()),
        StringTexts(seconds.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]()),
        uniform.value()[0],
        uniform.value()[1],
        Int(scoring.gaps.open),
        Int(scoring.gaps.extend),
        (pads[0], pads[1]),
        workers,
        scores_out,
        settled.unsafe_ptr(),
    )
    for index in range(pairs):
        if refused[index]:
            settled[index] = False
    return settled^


@fieldwise_init
struct CodedBatch(Movable):
    """A batch's sequences as an alphabet's codes for the lanes (see `lanes.CodeTexts`), the firsts' after
    the seconds', and the pairs holding a letter outside the alphabet, whose codes are left zero."""

    var codes: List[UInt8]
    var starts: List[Int]
    var refused: List[Bool]

    def firsts(self) -> CodeTexts:
        return CodeTexts(
            self.codes.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
            self.starts.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
        )

    def seconds(self) -> CodeTexts:
        return CodeTexts(
            self.codes.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
            self.starts.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]().unsafe_offset(len(self.refused)),
        )


def coded_batch(firsts: List[String], seconds: List[String], alphabet: String, workers: Int) -> CodedBatch:
    """Every pair's letters as `alphabet`'s codes, translated over `workers` threads."""
    var pairs = len(firsts)
    var codes_by_byte = Array[UInt8, 256](fill=UInt8(255))
    var alphabet_bytes = alphabet.as_bytes()
    for index in range(len(alphabet_bytes)):
        codes_by_byte[Int(alphabet_bytes[index])] = UInt8(index)
    # Sequence `i` of the firsts, then of the seconds, from `starts[i]`, and one entry past the last.
    var starts = List[Int](capacity=2 * pairs + 1)
    var total = 0
    for index in range(pairs):
        starts.append(total)
        total += firsts[index].byte_length()
    for index in range(pairs):
        starts.append(total)
        total += seconds[index].byte_length()
    starts.append(total)
    var codes = List[UInt8](capacity=max(total, 1))
    codes.resize(unsafe_uninit_length=total)
    var refused = List[Bool](length=pairs, fill=False)
    var code_ptr = codes.unsafe_ptr()
    var start_ptr = starts.unsafe_ptr()
    var refused_ptr = refused.unsafe_ptr()

    def translate_stretch(
        stretch: Int,
    ) {
        imm firsts, imm seconds, imm codes_by_byte, imm pairs, imm workers, imm code_ptr, imm start_ptr, imm refused_ptr
    }:
        """Translates stretch `stretch`'s pairs, marking any holding a letter outside the alphabet."""
        for index in range(pairs * stretch // workers, pairs * (stretch + 1) // workers):
            var unknown = False
            var first = firsts[index].as_bytes()
            var target = code_ptr.unsafe_offset(start_ptr[unsafe_offset=index])
            for position in range(len(first)):
                var code = codes_by_byte[Int(first[position])]
                unknown = unknown or code == 255
                target[unsafe_offset=position] = 0 if code == 255 else code
            var second = seconds[index].as_bytes()
            target = code_ptr.unsafe_offset(start_ptr[unsafe_offset=index + pairs])
            for position in range(len(second)):
                var code = codes_by_byte[Int(second[position])]
                unknown = unknown or code == 255
                target[unsafe_offset=position] = 0 if code == 255 else code
            refused_ptr[unsafe_offset=index] = unknown

    spread(translate_stretch, workers, workers)
    return CodedBatch(codes^, starts^, refused^)


def lane_table(scoring: Scoring) -> SIMD[DType.uint8, TABLE_ENTRIES]:
    """The table's signed scores as bytes, row by row, for the lanes' byte shuffle."""
    var table = SIMD[DType.uint8, TABLE_ENTRIES](0)
    var size = scoring.alphabet_size()
    for cell in range(min(size * size, TABLE_ENTRIES)):
        table[cell] = bitcast[DType.uint8](Int8(scoring.substitutions[cell]))
    return table


def tabled_scores[
    mode: AlignmentMode
](
    firsts: List[String], seconds: List[String], scoring: Scoring, threads: Int, scores_out: MutPointer[Int32, _]
) -> List[Bool]:
    """The global or local scores of every pair the lanes take into `scores_out`, and which they were,
    under a table of up to `TABLE_ENTRIES` entries with more than one mismatch score, its sequences as the
    alphabet's codes. A global score is a cost, as `wavefront_penalties` folds a uniform table's reward:
    each pair costing twice the table's best score less its own, a gap letter that best less twice its
    extension, and a gap's opening twice its extension less its opening, all over their common factor."""
    var pairs = len(firsts)
    var workers = max(threads, 1)
    var size = scoring.alphabet_size()
    var none = List[Bool](length=pairs, fill=False)
    if size < 2 or size * size > TABLE_ENTRIES:
        return none^
    var best = Int.MIN
    var least = Int.MAX
    for cell in range(size * size):
        best = max(best, Int(scoring.substitutions[cell]))
        least = min(least, Int(scoring.substitutions[cell]))
    var open = Int(scoring.gaps.open)
    var extend = Int(scoring.gaps.extend)
    var batch = coded_batch(firsts, seconds, scoring.alphabet, workers)
    var settled = batch.refused.copy()
    comptime if mode == AlignmentMode.GLOBAL:
        var opening = 2 * (extend - open)
        var extension = best - 2 * extend
        if opening < 0 or extension <= 0:
            return none^
        var scale = gcd(opening, extension)
        for cell in range(size * size):
            scale = gcd(scale, 2 * (best - Int(scoring.substitutions[cell])))
        var table = SIMD[DType.uint8, TABLE_ENTRIES](0)
        for cell in range(size * size):
            var cost = 2 * (best - Int(scoring.substitutions[cell])) // scale
            if cost > 255:
                return none^
            table[cell] = UInt8(cost)
        var penalties = Penalties(0, opening // scale, extension // scale, scale, best, 0, 0)
        var costs = LaneCosts.one_piece(
            0, opening // scale, extension // scale, opening // scale, extension // scale
        ).tabled(size, table)
        var found = List[Optional[Int]](length=pairs, fill=None)
        _ = lane_distances(
            pairs,
            batch.firsts(),
            batch.seconds(),
            costs,
            Band(),
            Int.MAX,
            workers,
            found.unsafe_ptr(),
            settled.unsafe_ptr(),
        )
        for index in range(pairs):
            if batch.refused[index]:
                settled[index] = False
            elif settled[index]:
                var letters = firsts[index].byte_length() + seconds[index].byte_length()
                scores_out[unsafe_offset=index] = Int32(penalties.score(found[index].value(), letters))
    else:
        # Codes past the alphabet pad the lanes.
        _ = lane_local_scores(
            pairs,
            batch.firsts(),
            batch.seconds(),
            best,
            least,
            open,
            extend,
            (UInt8(size), UInt8(size + 1)),
            workers,
            scores_out,
            settled.unsafe_ptr(),
            size,
            lane_table(scoring),
        )
        for index in range(pairs):
            if batch.refused[index]:
                settled[index] = False
    return settled^


def scores_with[
    mode: AlignmentMode
](firsts: List[String], seconds: List[String], scoring: Scoring, placement: Optional[Placement] = None) raises -> List[
    Int32
]:
    """Scores every pair; on the device, every pair one block can carry goes out in one launch."""
    var resolved = placement.or_else(Placement.default())
    var pairs = paired_length(firsts, seconds)
    var results = List[Int32](length=pairs, fill=0)
    if pairs == 0:
        return results^
    if resolved.device != Device.GPU:
        var out = results.unsafe_ptr()
        var failed = List[Bool](length=pairs, fill=False)
        var flags = failed.unsafe_ptr()
        var single = Placement.on_cpu(1)
        # A global score under a table of one match and one mismatch score is a cost the wavefront's
        # penalties count, the reward folded in (see `wavefront_penalties`): many pairs at once in the
        # lanes of a register (see `lanes`), every pair whose letters the alphabet holds and 16 bits its
        # cost. The rest, and every other table and mode, one at a time.
        var settled = List[Bool](length=pairs, fill=False)
        comptime if mode == AlignmentMode.GLOBAL:
            var penalties = wavefront_penalties(
                scoring.substitutions, scoring.alphabet_size(), Int(scoring.gaps.open), Int(scoring.gaps.extend)
            )
            if penalties:
                var found = penalties.value()
                var lane_costs = LaneCosts.one_piece(
                    found.mismatch, found.opening, found.extension, found.opening, found.extension
                )
                settled = laned_scores(firsts, seconds, scoring.alphabet, lane_costs, found, resolved.threads, out)
        comptime if mode == AlignmentMode.LOCAL:
            settled = laned_local_scores(firsts, seconds, scoring, resolved.threads, out)
        # A table of more than one mismatch score, small enough for a register, over its codes.
        if not uniform_scores(scoring):
            settled = tabled_scores[mode](firsts, seconds, scoring, resolved.threads, out)
        var settled_ptr = settled.unsafe_ptr()

        # The pairs are independent, so each is aligned on one thread start to finish, in
        # contiguous chunks, several a thread so one that draws long pairs does not hold up the rest.
        var chunks = chunk_count(pairs, resolved.threads)

        def score_range(slot: Int) {imm}:
            """Scores chunk `slot` of the pairs, flagging any pair that raised for the serial retry below."""
            for index in range(pairs * slot // chunks, pairs * (slot + 1) // chunks):
                if settled_ptr[unsafe_offset=index]:
                    continue
                try:
                    out[unsafe_offset=index] = score_with[mode](firsts[index], seconds[index], scoring, single)
                except:
                    flags[unsafe_offset=index] = True

        spread(score_range, chunks, max(resolved.threads, 1))
        # A pair that failed raises here, the same error a serial loop would have raised first.
        for index in range(pairs):
            if failed[index]:
                results[index] = score_with[mode](firsts[index], seconds[index], scoring, single)
        return results^

    var scope = DeviceScope(resolved.gpu_id)
    # A global score over the band its cost proves, a thread a pair, or over its whole matrix where that
    # band does not prove it (see `band_groups`).
    var remaining = List[Int](capacity=pairs)
    for index in range(pairs):
        remaining.append(index)
    comptime if mode == AlignmentMode.GLOBAL:
        var swept = banded_scores(
            scope,
            firsts,
            seconds,
            scoring.alphabet,
            scoring.substitutions,
            scoring.gaps,
            resolved.threads,
            resolved.gpu_id,
        )
        if swept:
            ref answer = swept.value()
            results = answer[0].copy()
            remaining.clear()
            for index in range(pairs):
                if not answer[1][index]:
                    remaining.append(index)
            if len(remaining) == 0:
                return results^
    # The rest whole, several pairs a warp, whatever their rows, wherever every second sequence fits a
    # shape; otherwise a warp a pair where one block's carry can index the pair, and a tiled sweep where
    # it cannot.
    var grouped = grouped_scores[mode](
        scope,
        firsts,
        seconds,
        remaining,
        scoring.alphabet,
        scoring.substitutions,
        scoring.gaps,
        resolved.threads,
        resolved.gpu_id,
    )
    if grouped:
        var scored = grouped.take()
        for slot in range(len(remaining)):
            results[remaining[slot]] = scored[slot]
        return results^
    var band = band_length(scope.specs)
    var banded = List[Int]()
    for index in remaining:
        if serving_space(firsts[index].byte_length(), band) == Space.BANDED:
            banded.append(index)
        else:
            results[index] = device_score[mode](
                scope,
                translate(firsts[index], scoring.alphabet),
                translate(seconds[index], scoring.alphabet),
                scoring.substitutions,
                scoring.alphabet_size(),
                scoring.gaps,
            )
    if len(banded) > 0:
        var tape = pack_batch(firsts, seconds, banded, scoring.alphabet, resolved.threads)
        var scored = device_scores[mode](
            scope, tape.sequences, tape.offsets, scoring.substitutions, scoring.alphabet_size(), scoring.gaps
        )
        for slot in range(len(banded)):
            results[banded[slot]] = scored[slot]
    return results^


def alignments_with[
    mode: AlignmentMode
](
    firsts: List[String],
    seconds: List[String],
    scoring: Scoring,
    placement: Optional[Placement] = None,
    stored_cells: Int = cells_within(DEFAULT_MAX_MEMORY),
) raises -> List[GappedAlignment]:
    """Aligns every pair; on the device, every pair both bounds admit goes out in one launch."""
    var resolved = placement.or_else(Placement.default())
    var pairs = paired_length(firsts, seconds)
    var results = List[GappedAlignment](capacity=pairs)
    for _ in range(pairs):
        results.append(GappedAlignment(0, String(), String()))
    if pairs == 0:
        return results^
    if resolved.device != Device.GPU:
        var out = results.unsafe_ptr()
        var failed = List[Bool](length=pairs, fill=False)
        var flags = failed.unsafe_ptr()
        var single = Placement.on_cpu(1)

        # The pairs are independent, so each is aligned on one thread start to finish, in
        # contiguous chunks, several a thread so one that draws long pairs does not hold up the rest.
        var chunks = chunk_count(pairs, resolved.threads)

        def align_range(slot: Int) {imm}:
            """Aligns chunk `slot` of the pairs, flagging any pair that raised for the serial retry below."""
            for index in range(pairs * slot // chunks, pairs * (slot + 1) // chunks):
                try:
                    out[unsafe_offset=index] = align_with[mode](
                        firsts[index], seconds[index], scoring, single, stored_cells
                    )
                except:
                    flags[unsafe_offset=index] = True

        spread(align_range, chunks, max(resolved.threads, 1))
        # A pair that failed raises here, the same error a serial loop would have raised first.
        for index in range(pairs):
            if failed[index]:
                results[index] = align_with[mode](firsts[index], seconds[index], scoring, single, stored_cells)
        return results^

    # A batch of one is a single pair however it arrived, and gets the single pair's crossover.
    var limit = stored_cells if pairs > 1 else min(stored_cells, DEVICE_STORED_CELLS)
    var scope = DeviceScope(resolved.gpu_id)
    var band = band_length(scope.specs)
    var batchable = List[Int]()
    for index in range(pairs):
        var rows = firsts[index].byte_length()
        if rows * seconds[index].byte_length() <= limit and serving_space(rows, band) == Space.BANDED:
            batchable.append(index)
        else:
            results[index] = align_on_device[mode](
                scope,
                translate(firsts[index], scoring.alphabet),
                translate(seconds[index], scoring.alphabet),
                scoring,
                stored_cells,
                resolved,
            )
    if len(batchable) > 0:
        var tape = pack_batch(firsts, seconds, batchable, scoring.alphabet, resolved.threads)
        var aligned = device_alignments[mode](
            scope, tape.sequences, tape.offsets, scoring.substitutions, scoring.alphabet, scoring.gaps
        )
        # Each alignment moved out rather than copied, the last slot first as `pop` hands them back.
        for slot in range(len(batchable) - 1, -1, -1):
            results[batchable[slot]] = aligned.pop()
    return results^


# endregion Entry Points
