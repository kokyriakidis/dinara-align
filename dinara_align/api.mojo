"""
The entry points a caller uses, each deciding which kernel serves a pair so the caller never has to.

Two choices are made here and nowhere else, and both are cost rather than correctness:

- The traceback keeps a decision per cell while the matrix fits `stored_budget`, and recurses in
  linear space once it does not. Both return an optimal path for the same score.
- On the device, a pair short enough for one block's shared-memory carry takes the banded strip
  sweep, and a taller one tiles over global memory. The bound comes from what the card reports.

A batch on the device goes out as one launch for every pair both bounds admit; the rest are served
one by one, so a single oversized pair never sinks the batch it arrived in.
"""

from .alignment import (
    AffineGapCosts,
    AlignmentMode,
    AlignmentResult,
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
    serial_local_extremum,
    serial_score,
    serving_space,
)
from .common import (
    Device,
    DeviceScope,
    FALLBACK_LETTER,
    MAX_ALPHABET_SIZE,
    OffsetDType,
    Placement,
    SubstitutionDType,
    SymbolDType,
    translate,
    uniform_matrix,
)
from .errors import AlignmentError, ErrorKind
from .gap_affine import wavefront_penalties, wavefront_score
from .vector_score import uniform_table, vector_align, vector_score

from max.algorithm import parallelize

comptime STORED_MATRIX_BUDGET = 6_000_000
"""
Cells above which the host traceback switches to the linear-space recursion. Three `int32` layers at twelve bytes a
cell keep a stored host alignment near 72 MB; the device packs a nibble per cell and is capped again by
`DEVICE_STORED_CELLS`.
"""

# region Scoring

comptime DNA_ALPHABET = "ACGT"
"""The four bases. An `N` or a soft-masked lowercase base must be added to an alphabet explicitly."""

comptime DEFAULT_MATCH = 2
comptime DEFAULT_MISMATCH = -4
comptime DEFAULT_GAP_OPENING = -6
comptime DEFAULT_GAP_EXTENSION = -2
"""
Minimap2's defaults (`-A2 -B4 -O4 -E2`). Minimap2 charges a gap of length `k` as `4 + 2k`, while this package charges
its first letter the opening and every further letter the extension, so the same model is an opening of six.
"""


struct Scoring(Copyable, Movable):
    """An alphabet, the substitution table it indexes, and the affine gap model, which travel together.

    Built through `dna`, `edit_distance`, `uniform` or `tabulated` rather than field by field, so a
    table whose shape disagrees with its alphabet cannot be expressed. A gap of length `k` costs
    `opening + (k - 1) * extension`.
    """

    var alphabet: String
    """The letters a sequence may hold, in the order the table is indexed by."""
    var substitutions: List[Scalar[SubstitutionDType]]
    """Row-major, one row per letter of `alphabet`."""
    var gaps: AffineGapCosts
    """The opening and extension penalties, both non-positive."""

    def __init__(
        out self, var alphabet: String, var substitutions: List[Scalar[SubstitutionDType]], gaps: AffineGapCosts
    ):
        """Trusts its arguments; the factories are the checked way in."""
        self.alphabet = alphabet^
        self.substitutions = substitutions^
        self.gaps = gaps

    @staticmethod
    def dna() raises AlignmentError -> Self:
        """Minimap2's scoring over `ACGT`: match 2, mismatch -4, gap opening -6 and extension -2."""
        return Self.uniform(DEFAULT_MATCH, DEFAULT_MISMATCH, DEFAULT_GAP_OPENING, DEFAULT_GAP_EXTENSION)

    @staticmethod
    def edit_distance(alphabet: String = String(DNA_ALPHABET)) raises AlignmentError -> Self:
        """Unit costs, under which a global score is the negated Levenshtein distance."""
        return Self.uniform(0, -1, -1, -1, alphabet)

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
            AffineGapCosts.checked(Int32(opening), Int32(extension)),
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
        return Self(alphabet, substitutions^, AffineGapCosts.checked(Int32(opening), Int32(extension)))

    def alphabet_size(self) -> Int:
        """Letters the table is indexed by, which is its stride."""
        return self.alphabet.byte_length()


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
    firsts: List[String], seconds: List[String], indices: List[Int], alphabet: String
) raises AlignmentError -> BatchTape:
    """The named pairs concatenated onto one tape, with offsets marking where each sequence begins."""
    var sequences = List[Scalar[SymbolDType]]()
    var offsets = List[Scalar[OffsetDType]]()
    offsets.append(0)
    for index in indices:
        sequences.extend(translate(firsts[index], alphabet))
        offsets.append(Scalar[OffsetDType](len(sequences)))
        sequences.extend(translate(seconds[index], alphabet))
        offsets.append(Scalar[OffsetDType](len(sequences)))
    return BatchTape(sequences^, offsets^)


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
) raises -> AlignmentResult:
    """Global alignment in linear space, splitting rows and joining halves Myers-Miller style.

    `vectorized` lets a uniform table's sweeps run sixteen cells at a time, which computes the same
    rows; off, every sweep runs cell by cell, as a reference to check it against.
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
        uniform_table(scoring.substitutions, scoring.alphabet_size()) if vectorized else None,
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
    return AlignmentResult(score, expanded[0], expanded[1])


def local_linear(
    first: ImmSpan[Scalar[SymbolDType], _], second: ImmSpan[Scalar[SymbolDType], _], scoring: Scoring
) raises -> AlignmentResult:
    """Local alignment in linear space, by reduction to the global problem.

    A forward local sweep finds where the best alignment ends, a backward sweep over those prefixes
    finds where it starts, and the global recursion then runs on that rectangle alone. Keeping one
    well-understood global recursion rather than four subproblem modes is what Myers and Miller
    themselves prescribe.
    """
    var alphabet_size = scoring.alphabet_size()
    var last_row, last_column, score = serial_local_extremum[SweepHalf.FORWARD](
        first, second, len(first), len(second), scoring.substitutions, alphabet_size, scoring.gaps
    )
    var path_columns = List[Int32](length=len(first) + 1, fill=Int32(0))
    var path_layers = List[Layer](length=len(first) + 1, fill=Layer.ALIGNING)

    var first_row = last_row
    if score > 0:
        var back_rows, back_columns, _ = serial_local_extremum[SweepHalf.REVERSE](
            first, second, last_row, last_column, scoring.substitutions, alphabet_size, scoring.gaps
        )
        first_row = last_row - back_rows
        var first_column = last_column - back_columns
        path_columns[last_row] = Int32(last_column)
        # The recursion runs on a copy, and only the rows of the core it solved are taken back.
        var core_columns = path_columns.copy()
        var core_layers = path_layers.copy()
        serial_hirschberg(
            first,
            second,
            first_row,
            last_row,
            first_column,
            last_column,
            scoring.substitutions,
            alphabet_size,
            scoring.gaps,
            DEFAULT_LEAF_CELLS,
            core_columns,
            core_layers,
            uniform_table(scoring.substitutions, scoring.alphabet_size()),
        )
        for index in range(first_row, last_row + 1):
            path_columns[index] = core_columns[index]
            path_layers[index] = core_layers[index]

    var expanded = expand_path(
        first, second, path_columns, path_layers, scoring.alphabet, AlignmentMode.LOCAL, first_row, last_row
    )
    return AlignmentResult(score, expanded[0], expanded[1])


# endregion Linear-Space Host Traceback

# region Routing


def align_on_host[
    mode: AlignmentMode
](
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    scoring: Scoring,
    stored_budget: Int,
) raises -> AlignmentResult:
    """One pair on the host, stored while its matrix fits the budget and linear once it does not.

    Stored under a table of one match and one mismatch score, the matrix is swept sixteen cells at a
    time, the same cells and so the same alignment (see `vector_align`).
    """
    if len(first) * len(second) > stored_budget:
        comptime if mode == AlignmentMode.LOCAL:
            return local_linear(first, second, scoring)
        else:
            return global_linear(first, second, scoring)
    var table = uniform_table(scoring.substitutions, scoring.alphabet_size())
    if table:
        var codes_first = List[UInt8](first)
        var codes_second = List[UInt8](second)
        return vector_align[mode](
            codes_first,
            codes_second,
            table.value()[0],
            table.value()[1],
            scoring.gaps,
            scoring.substitutions,
            scoring.alphabet_size(),
            scoring.alphabet,
        )
    return serial_align[mode](
        first, second, scoring.substitutions, scoring.alphabet_size(), scoring.gaps, scoring.alphabet
    )


def align_on_device[
    mode: AlignmentMode
](
    scope: DeviceScope,
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    scoring: Scoring,
    stored_budget: Int,
    placement: Placement,
) raises -> AlignmentResult:
    """One pair on the device, on the sweep its height and its matrix can afford.

    Both bounds are real and independent: the stored kernel indexes its carry by the first
    sequence, so a tall pair fails it even when the whole matrix would fit.
    """
    var stored = len(first) * len(second) <= min(stored_budget, DEVICE_STORED_CELLS)
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
    return aligned[0].copy()


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

# region Entry Points


def score[
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
    comptime if mode == AlignmentMode.GLOBAL:
        # A table of one match and one mismatch score has a wavefront, whose work grows with the
        # score rather than the matrix; it hands back a pair a full sweep would serve sooner.
        var penalties = wavefront_penalties(
            scoring.substitutions, scoring.alphabet_size(), Int(scoring.gaps.open), Int(scoring.gaps.extend)
        )
        if penalties:
            var found = wavefront_score(encoded_first, encoded_second, penalties.value())
            if found:
                return Int32(found.value())
    var table = uniform_table(scoring.substitutions, scoring.alphabet_size())
    if table:
        return vector_score[mode](encoded_first, encoded_second, table.value()[0], table.value()[1], scoring.gaps)
    return serial_score[mode](
        encoded_first, encoded_second, scoring.substitutions, scoring.alphabet_size(), scoring.gaps
    )


def align[
    mode: AlignmentMode
](
    first: String,
    second: String,
    scoring: Scoring,
    placement: Optional[Placement] = None,
    stored_budget: Int = STORED_MATRIX_BUDGET,
) raises -> AlignmentResult:
    """The optimal score and the two gapped strings that realize it.

    A global alignment spans both sequences; a local one returns only the best-scoring window,
    trimmed at both ends.
    """
    var resolved = placement.or_else(Placement.default())
    var encoded_first = translate(first, scoring.alphabet)
    var encoded_second = translate(second, scoring.alphabet)
    if resolved.device == Device.GPU:
        return align_on_device[mode](
            DeviceScope(resolved.gpu_id), encoded_first, encoded_second, scoring, stored_budget, resolved
        )
    return align_on_host[mode](encoded_first, encoded_second, scoring, stored_budget)


comptime CHUNKS_PER_THREAD = 8
"""Chunks a host batch is cut into per thread, so a thread that draws long pairs does not hold up
the rest while the others sit idle."""


def chunk_count(pairs: Int, threads: Int) -> Int:
    """How many contiguous chunks a host batch of `pairs` is cut into for `threads` threads."""
    return min(pairs, max(threads, 1) * CHUNKS_PER_THREAD)


def scores[
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

        # The pairs are independent, so each is aligned on one thread start to finish, in
        # contiguous chunks, several a thread so one that draws long pairs does not hold up the rest.
        var chunks = chunk_count(pairs, resolved.threads)

        def score_range(slot: Int) {imm}:
            for index in range(pairs * slot // chunks, pairs * (slot + 1) // chunks):
                try:
                    out[unsafe_offset=index] = score[mode](firsts[index], seconds[index], scoring, single)
                except:
                    flags[unsafe_offset=index] = True

        parallelize(score_range, chunks, max(resolved.threads, 1))
        # A pair that failed raises here, the same error a serial loop would have raised first.
        for index in range(pairs):
            if failed[index]:
                results[index] = score[mode](firsts[index], seconds[index], scoring, single)
        return results^

    var scope = DeviceScope(resolved.gpu_id)
    var band = band_length(scope.specs)
    var banded = List[Int]()
    for index in range(pairs):
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
        var tape = pack_batch(firsts, seconds, banded, scoring.alphabet)
        var scored = device_scores[mode](
            scope, tape.sequences, tape.offsets, scoring.substitutions, scoring.alphabet_size(), scoring.gaps
        )
        for slot in range(len(banded)):
            results[banded[slot]] = scored[slot]
    return results^


def alignments[
    mode: AlignmentMode
](
    firsts: List[String],
    seconds: List[String],
    scoring: Scoring,
    placement: Optional[Placement] = None,
    stored_budget: Int = STORED_MATRIX_BUDGET,
) raises -> List[AlignmentResult]:
    """Aligns every pair; on the device, every pair both bounds admit goes out in one launch."""
    var resolved = placement.or_else(Placement.default())
    var pairs = paired_length(firsts, seconds)
    var results = List[AlignmentResult](capacity=pairs)
    for _ in range(pairs):
        results.append(AlignmentResult(0, String(), String()))
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
            for index in range(pairs * slot // chunks, pairs * (slot + 1) // chunks):
                try:
                    out[unsafe_offset=index] = align[mode](
                        firsts[index], seconds[index], scoring, single, stored_budget
                    )
                except:
                    flags[unsafe_offset=index] = True

        parallelize(align_range, chunks, max(resolved.threads, 1))
        # A pair that failed raises here, the same error a serial loop would have raised first.
        for index in range(pairs):
            if failed[index]:
                results[index] = align[mode](firsts[index], seconds[index], scoring, single, stored_budget)
        return results^

    # A batch of one is a single pair however it arrived, and gets the single pair's crossover.
    var limit = stored_budget if pairs > 1 else min(stored_budget, DEVICE_STORED_CELLS)
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
                stored_budget,
                resolved,
            )
    if len(batchable) > 0:
        var tape = pack_batch(firsts, seconds, batchable, scoring.alphabet)
        var aligned = device_alignments[mode](
            scope, tape.sequences, tape.offsets, scoring.substitutions, scoring.alphabet, scoring.gaps
        )
        for slot in range(len(batchable)):
            results[batchable[slot]] = aligned[slot].copy()
    return results^


def needleman_wunsch_gotoh_score(
    first: String, second: String, scoring: Scoring, placement: Optional[Placement] = None
) raises -> Int32:
    """Global score: the path spans both sequences end to end."""
    return score[AlignmentMode.GLOBAL](first, second, scoring, placement)


def smith_waterman_gotoh_score(
    first: String, second: String, scoring: Scoring, placement: Optional[Placement] = None
) raises -> Int32:
    """Local score: the best-scoring window, never below zero."""
    return score[AlignmentMode.LOCAL](first, second, scoring, placement)


def needleman_wunsch_gotoh_alignment(
    first: String, second: String, scoring: Scoring, placement: Optional[Placement] = None
) raises -> AlignmentResult:
    """Global alignment: both sequences, gapped to one length."""
    return align[AlignmentMode.GLOBAL](first, second, scoring, placement)


def smith_waterman_gotoh_alignment(
    first: String, second: String, scoring: Scoring, placement: Optional[Placement] = None
) raises -> AlignmentResult:
    """Local alignment: the best-scoring window of each, gapped to one length."""
    return align[AlignmentMode.LOCAL](first, second, scoring, placement)


def combined_alphabet(first: String, second: String) -> String:
    """The distinct bytes of both strings, so unit-cost alignment needs no fixed alphabet."""
    var seen = List[Bool](length=256, fill=False)
    var letters = List[UInt8]()
    for text in [first, second]:
        for byte in text.as_bytes():
            if not seen[Int(byte)]:
                seen[Int(byte)] = True
                letters.append(byte)
    if len(letters) == 0:
        letters.append(FALLBACK_LETTER)
    return String(unsafe_from_utf8=letters)


def levenshtein_alignment(first: String, second: String) raises -> AlignmentResult:
    """Unit-cost edit distance with its alignment, on the host; the score is the distance, not its negation.

    Both gap penalties at minus one and substitutions at zero and minus one turn the maximizing
    Gotoh recurrence into the negated Levenshtein minimization, tie-break chain included, so this
    needs no kernel of its own. Any ASCII is accepted, since the alphabet is read off the inputs.
    """
    for text in [first, second]:
        for byte in text.as_bytes():
            if byte >= 0x80:
                raise AlignmentError(ErrorKind.NOT_ASCII, "unit-cost alignment")
    var alphabet = combined_alphabet(first, second)
    var size = alphabet.byte_length()
    var scoring = Scoring(alphabet, uniform_matrix(size, 0, -1), AffineGapCosts.checked(Int32(-1), Int32(-1)))
    var aligned = align_on_host[AlignmentMode.GLOBAL](
        translate(first, alphabet), translate(second, alphabet), scoring, STORED_MATRIX_BUDGET
    )
    aligned.score = -aligned.score
    return aligned^


# endregion Entry Points
