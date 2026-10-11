# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Alignment that maximizes a score under a `Scoring`: a substitution table over an alphabet, and affine gap
scores. Every call routes its pairs to whichever engine serves them fastest, all of them exact.

- On the host, a table of one match and one mismatch score goes to the wavefront (see `gap_affine`), or a
  batch's pairs many at once to the lanes of a register (see `lanes`). Any other table sweeps the band of
  diagonals its score bounds (see `vector_score`), and past the memory allowed takes the linear-space path
  (see `gotoh.linear_path`).
- On the device, a pair whose first sequence one block's carry can index joins a batch launch, and a taller
  one tiles over global memory (see `device_align`). A batch is cut into as many launches as the device's
  buffers need, and a pair no launch admits is served alone, so it never holds the others back.

Free ends and extensions find their span by sweep (see `scored.sweep`) and align the letters between
globally, which scores the same.
"""

from std.math import gcd

from .anti_diagonals import lane_bits
from .band_groups import banded_scores
from .cigar import reversed_list
from .common import (
    GAP_BYTE,
    MAX_ALPHABET_SIZE,
    UNKNOWN_SYMBOL,
    Device,
    DeviceScope,
    OffsetDType,
    Placement,
    SubstitutionDType,
    SymbolDType,
    code_table,
    raise_unknown,
    spread,
    translate,
    uniform_matrix,
)
from .device_align import (
    DEFAULT_LEAF_CELLS,
    DEVICE_STORED_CELLS,
    Space,
    band_length,
    device_align,
    device_alignments,
    device_score,
    device_scores,
    launch_bytes,
    serving_space,
)
from .errors import AlignmentError, ErrorKind
from .gap_affine import (
    DEFAULT_MAX_MEMORY,
    FREE_START,
    KEPT_BYTES,
    EndsFree,
    Penalties,
    SearchSpace,
    cigar_of,
    rings_fit,
    solve,
    wavefront_align,
    wavefront_penalties,
    wavefront_score,
)
from .gotoh import AffineGapCosts, AlignmentMode, GappedAlignment, Rectangle, RowPath, linear_path, serial_align
from .lanes import (
    TABLE_ENTRIES,
    CodeTexts,
    LaneCosts,
    LocalCosts,
    StringTexts,
    lane_alignments,
    lane_distances,
    lane_local_scores,
)
from .modes import Alignment, Anchor, Band, Costs, Mode
from .score_groups import grouped_scores
from .scored import ANYWHERE, FROM_EDGE, FROM_ORIGIN, costs_deficit, started_span, sweep
from .substitutions import SubstitutionLookup, shuffled_table, table_extremes, uniform_pair
from .vector_score import WIDTH as VECTOR_WIDTH, optimal_band, reach_back, vector_align, vector_bytes

# region Memory

comptime STORED_CELL_BYTES = 12
"""The unit `max_memory` is counted in for the host's tracebacks: a score in each of three layers, 32 bits
each, as the linear-space path's leaves store them. The banded traceback keeps a byte a cell (see
`vector_score.vector_bytes`)."""


def cells_bytes(cells: Int) -> Int:
    """The bytes `cells` stored cells take, saturating at `Int.MAX`."""
    if cells >= Int.MAX // STORED_CELL_BYTES:
        return Int.MAX
    return cells * STORED_CELL_BYTES


def cells_within(max_memory: Int) -> Int:
    """How many cells a traceback may store in `max_memory` bytes; past them it recurses in linear space."""
    return max(max_memory, 0) // STORED_CELL_BYTES


def fronts_within(stored_cells: Int) -> Int:
    """The wavefront's entries that fit where `stored_cells` cells would, past which it splits the pair where
    an optimal path crosses, as a `Costs` alignment does within its `max_memory`."""
    return max(stored_cells, 0) * STORED_CELL_BYTES // KEPT_BYTES


# endregion Memory

# region Scoring

comptime DNA_ALPHABET = "ACGT"
"""The four bases. An alphabet that should also take `N`, or soft-masked lowercase, must name those letters."""

comptime MINIMAP2_MATCH = 2
"""Minimap2's `-A2`."""
comptime MINIMAP2_MISMATCH = -4
"""Minimap2's `-B4`."""
comptime MINIMAP2_OPENING = -4
"""Minimap2's `-O4`; with `-E2`, a gap of `k` letters scores `-(4 + 2k)`."""
comptime MINIMAP2_EXTENSION = -2
"""Minimap2's `-E2`."""


@fieldwise_init
struct Scoring(Copyable, Movable):
    """What an alignment earns: a score for every pair of letters, and a score for every gap.

    The factories, `dna`, `edit_distance`, `uniform` and `tabulated`, check what they build: an alphabet
    the kernels can index, a square table over it, and gap scores that cost. The fieldwise constructor
    checks nothing. A gap of `k` letters scores `opening + k extension`, the way `Costs` prices one.
    """

    var alphabet: String
    """The letters a sequence may hold; letter `i` indexes row and column `i` of the table."""
    var substitutions: List[Scalar[SubstitutionDType]]
    """The table, row by row: what aligning letter `i` of the first sequence to letter `j` of the second
    scores sits at `i * alphabet_size() + j`."""
    var gaps: AffineGapCosts
    """The gap scores in the kernels' terms: `open` for a gap's first letter, `extend` for each one after."""

    @staticmethod
    def dna() raises AlignmentError -> Self:
        """Minimap2's defaults over `ACGT`: a match scores 2, a mismatch -4, a gap of `k` letters `-(4 + 2k)`."""
        return Self.uniform(MINIMAP2_MATCH, MINIMAP2_MISMATCH, MINIMAP2_OPENING, MINIMAP2_EXTENSION)

    @staticmethod
    def edit_distance(alphabet: String = DNA_ALPHABET) raises AlignmentError -> Self:
        """Every edit scoring -1 and a match 0, so a global score is minus the Levenshtein distance."""
        return Self.uniform(0, -1, 0, -1, alphabet)

    @staticmethod
    def uniform(
        match_score: Int,
        mismatch_score: Int,
        opening: Int = MINIMAP2_OPENING,
        extension: Int = MINIMAP2_EXTENSION,
        alphabet: String = DNA_ALPHABET,
    ) raises AlignmentError -> Self:
        """`match_score` for two equal letters and `mismatch_score` for two different ones, over `alphabet`."""
        var letters = letter_count(alphabet)
        return Self(alphabet, uniform_matrix(letters, match_score, mismatch_score), gap_scores(opening, extension))

    @staticmethod
    def tabulated(
        alphabet: String,
        var substitutions: List[Scalar[SubstitutionDType]],
        opening: Int = MINIMAP2_OPENING,
        extension: Int = MINIMAP2_EXTENSION,
    ) raises AlignmentError -> Self:
        """A table of the caller's, row by row, which must hold a score for every pair of `alphabet`'s letters."""
        var letters = letter_count(alphabet)
        if len(substitutions) != letters * letters:
            raise AlignmentError(
                ErrorKind.INVALID_SCORING,
                String(len(substitutions), " scores for ", letters * letters, " pairs of letters"),
            )
        return Self(alphabet, substitutions^, gap_scores(opening, extension))

    def alphabet_size(self) -> Int:
        """How many letters the table covers, the length of each of its rows."""
        return self.alphabet.byte_length()

    def penalties(self) -> Optional[Penalties]:
        """The wavefront's costs, if the table holds a single match and a single mismatch score whose folded
        costs a wavefront can grow by (see `gap_affine.wavefront_penalties`)."""
        return wavefront_penalties(self.substitutions, self.alphabet_size(), Int(self.gaps.open), Int(self.gaps.extend))


def letter_count(alphabet: String) raises AlignmentError -> Int:
    """How many letters `alphabet` holds, refusing none, more than a staged table covers, and `-`, which a
    gapped row reads as a gap."""
    var letters = alphabet.byte_length()
    if letters == 0 or letters > MAX_ALPHABET_SIZE:
        raise AlignmentError(ErrorKind.ALPHABET_TOO_LARGE, String(letters, " letters, from 1 to ", MAX_ALPHABET_SIZE))
    if GAP_BYTE in alphabet.as_bytes():
        raise AlignmentError(ErrorKind.INVALID_SCORING, "'-' names a gap, not a letter")
    return letters


def gap_scores(opening: Int, extension: Int) raises AlignmentError -> AffineGapCosts:
    """Gaps scoring `opening + k extension` over `k` letters, in the kernels' terms: `opening + extension` for
    a gap's first letter, `extension` for each one after."""
    if opening > 0:
        raise AlignmentError(ErrorKind.INVALID_SCORING, String("a gap opening of ", opening, " rewards a gap"))
    # Scores live in 32 bits, where a gap score past them would wrap. Each part is bounded before they are
    # summed, so the sum cannot wrap an `Int`, and from `-Int32.MAX` up, so that negating one still fits.
    var lowest = -Int(Int32.MAX)
    if min(opening, extension) < lowest or extension > Int(Int32.MAX) or opening + extension < lowest:
        raise AlignmentError(
            ErrorKind.INVALID_SCORING, String("a gap opening of ", opening, " and extension of ", extension)
        )
    return AffineGapCosts.checked(Int32(opening + extension), Int32(extension))


@fieldwise_init
struct ScoreReach(ImplicitlyCopyable, TrivialRegisterPassable):
    """How far from zero a `Scoring`'s cells can get over a pair. Every kernel keeps them in 32 bits, with its
    sentinels a quarter of the way down, so they must stay above that: a gap's first letter twice over, where
    the shifted kernels start their gap layers, plus twice the dearest of a pair's score and a gap's later
    letter for each letter of both sequences, or of the shorter for a local alignment, whose cells stop at
    zero."""

    var opening: Int
    var step: Int

    @staticmethod
    def of(scoring: Scoring) -> Self:
        """The reach of `scoring`'s table and gaps."""
        var extremes = table_extremes(scoring.substitutions, scoring.alphabet_size())
        return Self(2 * -Int(scoring.gaps.open), max(max(extremes[0], -extremes[1]), -Int(scoring.gaps.extend)))

    def check(self, rows: Int, columns: Int, floored: Bool) raises AlignmentError:
        """Refuses a pair of `rows` and `columns` letters whose scores could pass 32 bits."""
        var letters = min(rows, columns) if floored else rows + columns
        if self.opening + 2 * (letters + 2) * self.step >= 1 << 28:
            raise AlignmentError(ErrorKind.INVALID_SCORING, String("scores past 32 bits over ", rows, " by ", columns))


def first_past_32_bits(scoring: Scoring, firsts: List[String], seconds: List[String], floored: Bool) -> Int:
    """The earliest pair whose scores could pass 32 bits (see `ScoreReach`), or -1 for none."""
    var reach = ScoreReach.of(scoring)
    for index in range(min(len(firsts), len(seconds))):
        try:
            reach.check(firsts[index].byte_length(), seconds[index].byte_length(), floored)
        except:
            return index
    return -1


def pair_count(firsts: List[String], seconds: List[String]) raises AlignmentError -> Int:
    """How many pairs two lists make, refusing lists of different lengths."""
    if len(firsts) == len(seconds):
        return len(firsts)
    raise AlignmentError(
        ErrorKind.LENGTH_MISMATCH, String(len(firsts), " references against ", len(seconds), " queries")
    )


# endregion Scoring

# region Pair Tape


@fieldwise_init
struct PairTape(Movable):
    """Pairs laid end to end as alphabet codes, the layout a batch kernel reads: pair `p`'s first sequence
    runs from `starts[2p]` to `starts[2p + 1]`, and its second from there to `starts[2p + 2]`."""

    var codes: List[Scalar[SymbolDType]]
    var starts: List[Scalar[OffsetDType]]

    @staticmethod
    def single(first: ImmSpan[Scalar[SymbolDType], _], second: ImmSpan[Scalar[SymbolDType], _]) -> Self:
        """One pair already in codes, a batch of one."""
        var codes = List[Scalar[SymbolDType]](capacity=len(first) + len(second))
        codes.extend(first)
        codes.extend(second)
        var starts: List[Scalar[OffsetDType]] = [0, Scalar[OffsetDType](len(first)), Scalar[OffsetDType](len(codes))]
        return Self(codes^, starts^)

    @staticmethod
    def of(
        firsts: List[String], seconds: List[String], chosen: List[Int], alphabet: String, threads: Int
    ) raises AlignmentError -> Self:
        """The `chosen` pairs in `alphabet`'s codes, in their order. Every sequence's place is laid out first,
        so that each of `threads` threads can then write its share of the pairs straight onto the tape:
        translated one sequence at a time on one thread, a batch of 500,000 short reads took three times as
        long to pack as the kernel took to score it. A letter outside the alphabet raises as `translate`
        does, for the first chosen pair that holds one."""
        var pairs = len(chosen)
        var starts = List[Scalar[OffsetDType]](capacity=2 * pairs + 1)
        var length = 0
        starts.append(0)
        for index in chosen:
            length += firsts[index].byte_length()
            starts.append(Scalar[OffsetDType](length))
            length += seconds[index].byte_length()
            starts.append(Scalar[OffsetDType](length))
        var codes = List[Scalar[SymbolDType]](capacity=length)
        codes.resize(unsafe_uninit_length=length)
        var unknown = List[Bool](length=pairs, fill=False)
        var codes_by_byte = code_table(alphabet)
        var tape = codes.unsafe_ptr()
        var places = starts.unsafe_ptr()
        var flags = unknown.unsafe_ptr()
        var shares = max(min(threads, pairs), 1)

        def write_share(
            share: Int,
        ) {
            imm firsts,
            imm seconds,
            imm chosen,
            imm codes_by_byte,
            imm pairs,
            imm shares,
            imm tape,
            imm places,
            imm flags,
        }:
            """Writes share `share` of the pairs onto the tape, flagging each that holds an unknown letter."""
            for slot in range(pairs * share // shares, pairs * (share + 1) // shares):
                var index = chosen[slot]
                var first_known = encoded_into(
                    firsts[index], codes_by_byte, tape.unsafe_offset(Int(places[unsafe_offset=2 * slot]))
                )
                var second_known = encoded_into(
                    seconds[index], codes_by_byte, tape.unsafe_offset(Int(places[unsafe_offset=2 * slot + 1]))
                )
                if not (first_known and second_known):
                    flags[unsafe_offset=slot] = True

        spread(write_share, shares, shares)
        for slot in range(pairs):
            if unknown[slot]:
                raise_unknown(firsts[chosen[slot]], seconds[chosen[slot]], alphabet)
        return Self(codes^, starts^)


@always_inline
def encoded_into(text: String, codes_by_byte: Array[UInt8, 256], target: MutPointer[Scalar[SymbolDType], _]) -> Bool:
    """Writes `text` in codes from `target` on, and whether the alphabet held every one of its letters."""
    var bytes = text.unsafe_ptr()
    var known = True
    for position in range(text.byte_length()):
        var code = codes_by_byte[Int(bytes[unsafe_offset=position])]
        known = known and code != UNKNOWN_SYMBOL
        target[unsafe_offset=position] = Scalar[SymbolDType](code)
    return known


# endregion Pair Tape

# region Host Pairs


def host_score[
    mode: AlignmentMode
](first: ImmSpan[Scalar[SymbolDType], _], second: ImmSpan[Scalar[SymbolDType], _], scoring: Scoring) -> Int32:
    """The best score on the host. A global one under a table of one match and one mismatch score tries the
    wavefront, whose work grows with the score rather than the matrix, and whose rings the default memory
    holds; the wavefront hands back a pair a sweep would finish sooner. Everything else is swept (see
    `swept_score`)."""
    var first_codes = List[UInt8](first)
    var second_codes = List[UInt8](second)
    comptime if mode == AlignmentMode.GLOBAL:
        var penalties = scoring.penalties()
        if penalties and rings_fit(penalties.value(), len(first_codes), len(second_codes), DEFAULT_MAX_MEMORY):
            var found = wavefront_score(first_codes, second_codes, penalties.value())
            if found:
                return Int32(found.value())
    return Int32(swept_score[mode](Span(first_codes), Span(second_codes), scoring))


def host_text_score[mode: AlignmentMode](first: String, second: String, scoring: Scoring) raises -> Int32:
    """`host_score` of two sequences of letters, refusing a letter outside the alphabet."""
    var first_codes = translate(first, scoring.alphabet)
    var second_codes = translate(second, scoring.alphabet)
    return host_score[mode](first_codes, second_codes, scoring)


def swept_score[mode: AlignmentMode](first: Span[UInt8, _], second: Span[UInt8, _], scoring: Scoring) -> Int:
    """The best score by a sweep over the whole matrix under any table, in lanes of 16 bits while the scores
    fit (see `tabulated_end`): a global alignment as free ends with no end free, a local one anywhere. Against
    an empty sequence a global alignment is one gap and a local one nothing."""
    var empty = len(first) == 0 or len(second) == 0
    comptime if mode == AlignmentMode.LOCAL:
        if empty:
            return 0
        return tabulated_end[ANYWHERE](first, second, scoring, EndsFree(), True)[0]
    if empty:
        var letters = len(first) + len(second)
        return 0 if letters == 0 else Int(scoring.gaps.open) + (letters - 1) * Int(scoring.gaps.extend)
    return tabulated_end[FROM_EDGE](first, second, scoring, EndsFree(), True)[0]


def host_alignment[
    mode: AlignmentMode
](
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    scoring: Scoring,
    stored_cells: Int,
) raises -> GappedAlignment:
    """An optimal alignment on the host, in memory that grows with the alignment rather than the matrix: a
    local one by its span (see `local_by_span`), a global one as `global_on_host` routes it."""
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
    """An optimal global alignment on the host.

    Under a table of one match and one mismatch score whose costs a wavefront can grow by, the wavefront
    from both ends traces it through its own fronts, split where they would outgrow `stored_cells` (see
    `gap_affine.wavefront_align`). Under any other table it is traced in a band of diagonals that holds
    every optimal path, swept by anti-diagonal with each cell's decision kept (see `vector_align`), the
    band proved by the score it yields (see `certified_alignment`); a score `known` beforehand, as a local
    alignment's span knows it, bounds the band outright (see `optimal_band`). A band whose decisions would
    pass `stored_cells`' bytes takes the linear-space path.
    """
    if len(first) == 0 or len(second) == 0:
        return serial_align[AlignmentMode.GLOBAL](
            first, second, scoring.substitutions, scoring.alphabet_size(), scoring.gaps, scoring.alphabet
        )
    var first_codes = List[UInt8](first)
    var second_codes = List[UInt8](second)
    var penalties = scoring.penalties()
    if penalties and rings_fit(penalties.value(), len(first_codes), len(second_codes), cells_bytes(stored_cells)):
        var traced = wavefront_align(
            first_codes, second_codes, penalties.value(), scoring.alphabet, fronts_within(stored_cells)
        )
        return GappedAlignment(Int32(traced[0]), traced[1], traced[2])
    var lookup = SubstitutionLookup(scoring.substitutions, scoring.alphabet_size())
    # Not `known.or_else(...)`, whose argument would sweep for the score even when it is known.
    var best: Int
    if known:
        best = known.value()
    else:
        var certified = certified_alignment(first_codes, second_codes, lookup, scoring, stored_cells)
        if certified:
            return certified.take()
        best = swept_score[AlignmentMode.GLOBAL](Span(first_codes), Span(second_codes), scoring)
    var band = optimal_band(len(first), len(second), lookup.best, scoring.gaps, best)
    var width = min(band[1], len(second)) - max(band[0], -len(first)) + 1
    # Counted with its indexes and its three diagonals of scores, which a narrow band's decisions are fewer than.
    if vector_bytes(len(first), len(second), (len(first) + 1) * width) > cells_bytes(stored_cells):
        return linear_global(first, second, scoring)
    return vector_align(
        first_codes,
        second_codes,
        lookup,
        scoring.gaps,
        scoring.substitutions,
        scoring.alphabet_size(),
        scoring.alphabet,
        band[0],
        band[1],
    )


def certified_alignment(
    first: List[UInt8], second: List[UInt8], lookup: SubstitutionLookup, scoring: Scoring, stored_cells: Int
) -> Optional[GappedAlignment]:
    """An optimal global alignment traced in bands of diagonals that grow until one proves itself, or None
    when the next would pass `stored_cells`' bytes.

    A band holding the start's and the end's diagonals and the run between them holds a path, so its
    traced score is a score some alignment earns, at most the best; every alignment scoring that much or
    more stays inside the band `optimal_band` gives for it. When that band lies within the one traced,
    so does every optimal path: the traced score is the best, and the walk, whose every decision on an
    optimal path reads only cells on optimal paths, takes the path the whole matrix's walk takes. Else
    the band grows toward the one the score proves, at most fourfold a step, so a first score far below
    the best, from a band that missed the path, costs a few narrower sweeps rather than the whole matrix.
    The first band reaches a step's lanes either side, so a pair close to its diagonal takes one sweep
    of about its own cells, where a sweep for the score alone covered the matrix.
    """
    var rows = len(first)
    var columns = len(second)
    var difference = columns - rows
    var low = max(min(0, difference) - 2 * VECTOR_WIDTH, -rows)
    var high = min(max(0, difference) + 2 * VECTOR_WIDTH, columns)
    while True:
        var width = high - low + 1
        if vector_bytes(rows, columns, (rows + 1) * width) > cells_bytes(stored_cells):
            return None
        var traced = vector_align(
            first,
            second,
            lookup,
            scoring.gaps,
            scoring.substitutions,
            scoring.alphabet_size(),
            scoring.alphabet,
            low,
            high,
        )
        var needed = optimal_band(rows, columns, lookup.best, scoring.gaps, Int(traced.score))
        var needed_low = max(needed[0], -rows)
        var needed_high = min(needed[1], columns)
        if needed_low >= low and needed_high <= high:
            return traced^
        if needed_low < low:
            low = max(needed_low, low - 3 * width // 2)
        if needed_high > high:
            high = min(needed_high, high + 3 * width // 2)


def linear_global(
    first: ImmSpan[Scalar[SymbolDType], _], second: ImmSpan[Scalar[SymbolDType], _], scoring: Scoring
) raises -> GappedAlignment:
    """An optimal global alignment in memory linear in the pair (see `gotoh.linear_path`)."""
    var rows = len(first)
    var path = RowPath(rows)
    linear_path(
        first,
        second,
        Rectangle(0, rows, 0, len(second)),
        scoring.substitutions,
        scoring.alphabet_size(),
        scoring.gaps,
        DEFAULT_LEAF_CELLS,
        path,
    )
    var score = path.score(first, second, scoring.substitutions, scoring.alphabet_size(), scoring.gaps)
    var gapped = path.gapped(first, second, scoring.alphabet, AlignmentMode.GLOBAL, 0, rows)
    return GappedAlignment(score, gapped[0], gapped[1])


def local_by_span(
    first: ImmSpan[Scalar[SymbolDType], _], second: ImmSpan[Scalar[SymbolDType], _], scoring: Scoring, stored_cells: Int
) raises -> GappedAlignment:
    """An optimal local alignment by its span, as SSW finds one. The forward sweep's best cell ends it, in
    lanes of 16 bits while the scores fit (see `tabulated_end`); a sweep back from that cell to where the
    score is earned starts it (see `reach_back`); and the letters between, aligned globally, score the same.
    Only the forward sweep covers the matrix."""
    var first_codes = List[UInt8](first)
    var second_codes = List[UInt8](second)
    var end = tabulated_end[ANYWHERE](Span(first_codes), Span(second_codes), scoring, EndsFree(), True)
    var best = end[0]
    if best <= 0:
        return GappedAlignment(0, String(), String())
    var lookup = SubstitutionLookup(scoring.substitutions, scoring.alphabet_size())
    var narrow = table_bits[FROM_EDGE](scoring, end[1], end[2]) == 16
    var start = reach_back(first_codes, second_codes, end[1], end[2], lookup, scoring.gaps, Int32(best), narrow)
    return global_on_host(first[start[0] : end[1]], second[start[1] : end[2]], scoring, stored_cells, best)


# endregion Host Pairs

# region Device Pairs


@always_inline
def batched(rows: Int, band: Int) -> Bool:
    """Whether a batch kernel takes a pair with `rows` letters in its first sequence: its carry, `band` rows
    long, holds them. A taller pair tiles over global memory."""
    return serving_space(rows, band) == Space.BANDED


def device_pair_score[
    mode: AlignmentMode
](
    scope: DeviceScope,
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    scoring: Scoring,
) raises -> Int32:
    """The best score of one pair on the device, as a batch of one if a block can carry it, else tiled."""
    if batched(len(first), band_length(scope.specs)):
        var tape = PairTape.single(first, second)
        return device_scores[mode](
            scope, tape.codes, tape.starts, scoring.substitutions, scoring.alphabet_size(), scoring.gaps
        )[0]
    return device_score[mode](scope, first, second, scoring.substitutions, scoring.alphabet_size(), scoring.gaps)


def device_pair_alignment[
    mode: AlignmentMode
](
    scope: DeviceScope,
    first: ImmSpan[Scalar[SymbolDType], _],
    second: ImmSpan[Scalar[SymbolDType], _],
    scoring: Scoring,
    stored_cells: Int,
    placement: Placement,
) raises -> GappedAlignment:
    """An optimal alignment of one pair on the device. A batch of one records every cell's decision, which
    needs a block to carry the pair and its matrix to fit both `stored_cells` and the device's own bound;
    those are separate limits, since a tall pair overruns the carry with a matrix that would fit. Any other
    pair takes the linear-space path over tiles (see `device_align.device_align`)."""
    var recorded = len(first) * len(second) <= min(stored_cells, DEVICE_STORED_CELLS)
    if recorded and batched(len(first), band_length(scope.specs)):
        var tape = PairTape.single(first, second)
        var aligned = device_alignments[mode](
            scope, tape.codes, tape.starts, scoring.substitutions, scoring.alphabet, scoring.gaps
        )
        return aligned.pop()
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


# endregion Device Pairs

# region Every Mode


def gap_costs(scoring: Scoring) -> Costs:
    """A `Scoring`'s gaps priced as `Costs` price them, in the sweep's terms: a gap scoring
    `opening + k extension` costs minus that. Its mismatch cost is a placeholder, since the table scores pairs."""
    return Costs(1, Int(scoring.gaps.extend) - Int(scoring.gaps.open), -Int(scoring.gaps.extend), -1, 0)


def table_bits[kind: Int](scoring: Scoring, rows: Int, columns: Int) -> Int:
    """The narrowest lanes that hold every score of a sweep under `scoring` (see `anti_diagonals.lane_bits`),
    bounded by the table's best pair, its dearest step, a mismatch or a gap's first letter, and how far a
    cell falls below zero (see `scored.costs_deficit`)."""
    var extremes = table_extremes(scoring.substitutions, scoring.alphabet_size())
    var substitution = -min(extremes[1], 0)
    var costs = gap_costs(scoring)
    var dearest = max(substitution, costs.opening + costs.extension)
    var deficit = costs_deficit(costs, substitution, rows, columns)
    return lane_bits[kind == ANYWHERE](max(extremes[0], 0), dearest, deficit, rows, columns)


def tabulated_end[
    kind: Int
](
    first: Span[UInt8, _], second: Span[UInt8, _], scoring: Scoring, ends: EndsFree, highest: Bool, zdrop: Int = -1
) -> Tuple[Int, Int, Int, Bool]:
    """The best score of an alignment `kind` allows under `scoring`, the letters of each sequence it ends
    after, and whether a Z-drop stopped the sweep: the sweep of `scored.sweep`, reading each pair's score
    from the table, over `first` and `second` in the alphabet's codes, in lanes along the shorter sequence as
    narrow as the scores allow."""
    var unused = List[Int]()
    return sweep[kind](
        first,
        second,
        gap_costs(scoring),
        0,
        ends,
        highest,
        table_bits[kind](scoring, len(first), len(second)),
        scoring.substitutions,
        scoring.alphabet_size(),
        zdrop,
        unused,
    )


def extension_span(
    first: List[Scalar[SymbolDType]], second: List[Scalar[SymbolDType]], scoring: Scoring, mode: Mode
) -> Tuple[Int, Int, Int, Int, Int, Bool]:
    """An extension's best stop and the span it covers, as `mode_span` gives them, and whether the Z-drop
    stopped the sweep: one sweep out from its anchor, its stop the latest of equally good ones."""
    var columns = len(first)
    var rows = len(second)
    if mode.anchor == Anchor.END:
        var reversed_first = reversed_list(Span(first))
        var reversed_second = reversed_list(Span(second))
        var stop = tabulated_end[FROM_ORIGIN](reversed_first, reversed_second, scoring, EndsFree(), True, mode.zdrop)
        return (stop[0], columns - stop[1], rows - stop[2], columns, rows, stop[3])
    var stop = tabulated_end[FROM_ORIGIN](first, second, scoring, EndsFree(), True, mode.zdrop)
    return (stop[0], 0, 0, stop[1], stop[2], stop[3])


def mode_span(
    first: List[Scalar[SymbolDType]],
    second: List[Scalar[SymbolDType]],
    scoring: Scoring,
    mode: Mode,
    started: Bool = True,
) -> Tuple[Int, Int, Int, Int, Int]:
    """The best score with free ends or as an extension, and the span it covers, placed as `Costs` place it
    under `Ties.LEFT` (see `gap_affine.free_ends_alignment`): the score, then the first and the second
    sequence's letters where it starts, then where it ends.

    Free ends sweep in from the edges for the end on the highest diagonal and, with `started`, back from
    there for the start on the highest too; without it the start stays at the origin. An extension sweeps
    out from its anchor (see `extension_span`)."""
    if mode.kind == Mode.EXTENSION:
        var stop = extension_span(first, second, scoring, mode)
        # The end bonus favours the best extension that reaches the query's far end, unless the Z-drop gave up.
        if mode.end_bonus > 0 and not stop[5]:
            var reaching = mode_span(first, second, scoring, mode.reaching_end(), started)
            if reaching[0] + mode.end_bonus > stop[0]:
                return reaching
        return (stop[0], stop[1], stop[2], stop[3], stop[4])
    var ends = EndsFree.of(mode, len(first), len(second))
    var forward = tabulated_end[FROM_EDGE](first, second, scoring, ends, True)
    if not started:
        return (forward[0], 0, 0, forward[1], forward[2])

    def back(head: List[UInt8], lead: List[UInt8], starts: EndsFree) {imm scoring} -> Tuple[Int, Int, Int, Bool]:
        """The sweep back over the reversed sequences, for the start on the lowest diagonal."""
        return tabulated_end[FROM_EDGE](Span(head), Span(lead), scoring, starts, False)

    return started_span(Span(first), Span(second), ends, forward, back)


def as_alignment(gapped: GappedAlignment, first: String, second: String, whole: Bool, eqx: Bool) -> Alignment:
    """Gotoh's gapped rows as an `Alignment`. With `whole` they span both sequences; otherwise they are a local
    alignment's, placed at the last place its letters occur in each sequence: wherever they occur, they align
    the same letters for the same score."""
    var score = Int(gapped.score)
    var cigar = gapped.cigar(eqx)
    if whole:
        return Alignment(-score, score, cigar^, 0, first.byte_length(), 0, second.byte_length())
    var first_letters = gapped.first_gapped.replace("-", "")
    var second_letters = gapped.second_gapped.replace("-", "")
    var first_start = first.rfind(first_letters) if first_letters.byte_length() > 0 else 0
    var second_start = second.rfind(second_letters) if second_letters.byte_length() > 0 else 0
    return Alignment(
        -score,
        score,
        cigar^,
        first_start,
        first_start + first_letters.byte_length(),
        second_start,
        second_start + second_letters.byte_length(),
    )


def check_scoring_mode(first: String, second: String, scoring: Scoring, mode: Mode) raises AlignmentError:
    """Refuses a mode that would add a match score of its own to the table's, and a pair whose scores could
    pass 32 bits (see `ScoreReach`)."""
    if mode.match_score > 0:
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a match score in the mode as well as in the Scoring's table")
    ScoreReach.of(scoring).check(first.byte_length(), second.byte_length(), mode.kind == Mode.SMITH_WATERMAN)


def refuse_span_modes_on_gpu(placement: Optional[Placement]) raises AlignmentError:
    """Refuses a device placement for free ends or an extension, which align only on the host."""
    if placement and placement.value().device == Device.GPU:
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a Scoring aligns only globally or locally on the GPU")


def scoring_alignment(
    first: String,
    second: String,
    scoring: Scoring,
    mode: Mode,
    placement: Optional[Placement],
    stored_cells: Int,
    eqx: Bool,
) raises -> Alignment:
    """An optimal alignment under `scoring` as `mode` asks, as an `Alignment`, whose cost is minus its score.
    Global and local alignments run on either device (see `pair_alignment`), a local one placed where its
    letters occur; free ends and extensions run on the host, their span found by sweep (see `mode_span`)
    and aligned globally."""
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
    """`scoring_alignment` through the searches of `space`, which a batch's thread reuses from pair to pair."""
    check_scoring_mode(first, second, scoring, mode)
    var on_host = not placement or placement.value().device != Device.GPU
    if on_host and mode.is_global() and first.byte_length() > 0 and second.byte_length() > 0:
        var penalties = scoring.penalties()
        if penalties and rings_fit(
            penalties.value(), first.byte_length(), second.byte_length(), cells_bytes(stored_cells)
        ):
            return wavefront_alignment(first, second, scoring, penalties.value(), stored_cells, eqx, space)
    if mode.kind == Mode.SMITH_WATERMAN:
        var gapped = pair_alignment[AlignmentMode.LOCAL](first, second, scoring, placement, stored_cells)
        return as_alignment(gapped, first, second, False, eqx)
    if mode.is_global():
        var gapped = pair_alignment[AlignmentMode.GLOBAL](first, second, scoring, placement, stored_cells)
        return as_alignment(gapped, first, second, True, eqx)
    refuse_span_modes_on_gpu(placement)
    var first_codes = translate(first, scoring.alphabet)
    var second_codes = translate(second, scoring.alphabet)
    var span = mode_span(first_codes, second_codes, scoring, mode)
    var score = span[0]
    var first_start = span[1]
    var second_start = span[2]
    var first_end = span[3]
    var second_end = span[4]
    var inner = host_alignment[AlignmentMode.GLOBAL](
        Span(first_codes)[first_start:first_end], Span(second_codes)[second_start:second_end], scoring, stored_cells
    )
    return Alignment(-score, score, inner.cigar(eqx), first_start, first_end, second_start, second_end)


def wavefront_alignment(
    first: String,
    second: String,
    scoring: Scoring,
    penalties: Penalties,
    stored_cells: Int,
    eqx: Bool,
    mut space: SearchSpace,
) raises -> Alignment:
    """A global alignment under a table of one match and one mismatch score, its CIGAR spelled straight from
    the wavefront's moves with no gapped rows between: the alignment `host_alignment` gives (see
    `gap_affine.wavefront_align`)."""
    var first_codes = translate(first, scoring.alphabet)
    var second_codes = translate(second, scoring.alphabet)
    var letters = len(first_codes) + len(second_codes)
    var moves = List[UInt8](capacity=letters)
    var cost = solve[1](
        space.forward,
        space.backward,
        Span(first_codes),
        Span(second_codes),
        penalties,
        FREE_START,
        FREE_START,
        fronts_within(stored_cells),
        moves,
    )
    var score = penalties.score(cost, letters)
    var cigar = cigar_of(first, second, moves^, cost, penalties, eqx)
    return Alignment(-score, score, cigar^, 0, len(first_codes), 0, len(second_codes))


def scoring_score(
    first: String, second: String, scoring: Scoring, mode: Mode, placement: Optional[Placement]
) raises -> Int:
    """The best score under `scoring` as `mode` asks, with no alignment traced: global and local on either
    device (see `pair_score`), free ends and extensions by their sweep on the host."""
    check_scoring_mode(first, second, scoring, mode)
    if mode.kind == Mode.SMITH_WATERMAN:
        return Int(pair_score[AlignmentMode.LOCAL](first, second, scoring, placement))
    if mode.is_global():
        return Int(pair_score[AlignmentMode.GLOBAL](first, second, scoring, placement))
    refuse_span_modes_on_gpu(placement)
    # The score alone needs no sweep back for the start.
    var first_codes = translate(first, scoring.alphabet)
    var second_codes = translate(second, scoring.alphabet)
    return mode_span(first_codes, second_codes, scoring, mode, started=False)[0]


# endregion Every Mode

# region Single Pairs


def pair_score[
    mode: AlignmentMode
](first: String, second: String, scoring: Scoring, placement: Optional[Placement] = None) raises -> Int32:
    """The best global or local score of one pair, in memory linear in the pair, on the device `placement`
    names (see `host_score`, `device_pair_score`)."""
    var resolved = placement.or_else(Placement.default())
    if resolved.device != Device.GPU:
        return host_text_score[mode](first, second, scoring)
    var first_codes = translate(first, scoring.alphabet)
    var second_codes = translate(second, scoring.alphabet)
    return device_pair_score[mode](DeviceScope(resolved.gpu_id), first_codes, second_codes, scoring)


def pair_alignment[
    mode: AlignmentMode
](
    first: String,
    second: String,
    scoring: Scoring,
    placement: Optional[Placement] = None,
    stored_cells: Int = cells_within(DEFAULT_MAX_MEMORY),
) raises -> GappedAlignment:
    """An optimal global or local alignment of one pair as two gapped rows, on the device `placement` names
    (see `host_alignment`, `device_pair_alignment`). A global alignment spans both sequences; a local one
    holds only the best-scoring stretch of each."""
    var resolved = placement.or_else(Placement.default())
    var first_codes = translate(first, scoring.alphabet)
    var second_codes = translate(second, scoring.alphabet)
    if resolved.device != Device.GPU:
        return host_alignment[mode](first_codes, second_codes, scoring, stored_cells)
    return device_pair_alignment[mode](
        DeviceScope(resolved.gpu_id), first_codes, second_codes, scoring, stored_cells, resolved
    )


# endregion Single Pairs

# region Host Lanes


def unsettle(refused: List[Bool], mut settled: List[Bool]):
    """Marks each refused pair unsettled again, so that its own call raises for its unknown letter."""
    for index in range(len(refused)):
        if refused[index]:
            settled[index] = False


def unknown_letters(firsts: List[String], seconds: List[String], alphabet: String, workers: Int) -> List[Bool]:
    """Which pairs hold a letter outside `alphabet`, checked over `workers` threads. The lanes skip those, for
    their own calls to raise as they do."""
    var pairs = len(firsts)
    var codes_by_byte = code_table(alphabet)
    var refused = List[Bool](length=pairs, fill=False)
    var flags = refused.unsafe_ptr()

    def check_share(share: Int) {imm firsts, imm seconds, imm codes_by_byte, imm pairs, imm workers, imm flags}:
        """Flags the pairs of share `share` that hold a letter outside the alphabet."""
        for index in range(pairs * share // workers, pairs * (share + 1) // workers):
            var known = True
            for letter in firsts[index].as_bytes():
                known = known and codes_by_byte[Int(letter)] != UNKNOWN_SYMBOL
            for letter in seconds[index].as_bytes():
                known = known and codes_by_byte[Int(letter)] != UNKNOWN_SYMBOL
            if not known:
                flags[unsafe_offset=index] = True

    spread(check_share, workers, workers)
    return refused^


def uniform_scores(scoring: Scoring) -> Optional[Tuple[Int, Int]]:
    """The match and the mismatch score, if the table holds one of each (see `uniform_pair`)."""
    return uniform_pair(scoring.substitutions, scoring.alphabet_size())


def laned_scores(
    firsts: List[String],
    seconds: List[String],
    alphabet: String,
    costs: LaneCosts,
    penalties: Penalties,
    threads: Int,
    scores_out: MutPointer[Int32, _],
) -> List[Bool]:
    """Writes the global score of every pair the lanes settle to `scores_out`, and returns which pairs those
    were. A pair holding a letter outside `alphabet` is left for its own call to raise."""
    var pairs = len(firsts)
    var workers = max(threads, 1)
    # Marked settled beforehand, a refused pair is one the lanes skip.
    var refused = unknown_letters(firsts, seconds, alphabet, workers)
    var settled = refused.copy()
    var found = List[Optional[Int]](length=pairs, fill=None)
    _ = lane_distances(
        pairs,
        StringTexts.of(firsts),
        StringTexts.of(seconds),
        costs,
        Band(),
        Int.MAX,
        workers,
        found.unsafe_ptr(),
        settled.unsafe_ptr(),
    )
    unsettle(refused, settled)
    for index in range(pairs):
        if settled[index]:
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
    """Writes the global alignment of every pair the lanes settle to `alignments_out`, and returns which pairs
    those were, under a table of one match and one mismatch score. Each is the alignment `scoring_alignment`
    gives, the wavefront's path by `Ties.LEFT`, traced from the decisions a group's band kept within `budget`
    bytes (see `lanes.lane_alignments`). A pair is left for its own call if it holds a letter outside the
    alphabet, if a side is empty, or if the searches might split it for memory, since each piece of a split
    would follow the tie rule on its own."""
    var pairs = len(firsts)
    var workers = max(threads, 1)
    var penalties = scoring.penalties()
    if not penalties:
        return List[Bool](length=pairs, fill=False)
    var found = penalties.value()
    var refused = unknown_letters(firsts, seconds, scoring.alphabet, workers)
    var settled = refused.copy()
    var costs = List[Optional[Int]](length=pairs, fill=None)
    var paths = List[List[UInt8]](capacity=pairs)
    for _ in range(pairs):
        paths.append(List[UInt8]())
    _ = lane_alignments(
        pairs,
        StringTexts.of(firsts),
        StringTexts.of(seconds),
        LaneCosts.of_penalties(found),
        Band(),
        Int.MAX,
        True,
        workers,
        budget,
        costs.unsafe_ptr(),
        paths.unsafe_ptr(),
        settled.unsafe_ptr(),
    )
    unsettle(refused, settled)
    var settled_ptr = settled.unsafe_ptr()
    var cost_ptr = costs.unsafe_ptr()
    var path_ptr = paths.unsafe_ptr()

    def spell_share(
        share: Int,
    ) {
        imm firsts,
        imm seconds,
        imm found,
        imm eqx,
        imm pairs,
        imm workers,
        imm budget,
        imm settled_ptr,
        imm cost_ptr,
        imm path_ptr,
        imm alignments_out,
    }:
        """Spells the CIGAR of each pair of share `share` the lanes traced, and unsettles any the searches
        could have split."""
        for index in range(pairs * share // workers, pairs * (share + 1) // workers):
            if not settled_ptr[unsafe_offset=index]:
                continue
            var cost = cost_ptr[unsafe_offset=index].value()
            var columns = firsts[index].byte_length()
            var rows = seconds[index].byte_length()
            # The searches split a pair whose fronts outgrow the budget's entries, which they cannot below this.
            if 2 * (cost + 1) * (columns + rows + 1) > budget // KEPT_BYTES:
                settled_ptr[unsafe_offset=index] = False
                continue
            var moves = List[UInt8]()
            swap(moves, path_ptr[unsafe_offset=index])
            var score = found.score(cost, columns + rows)
            var cigar = cigar_of(firsts[index], seconds[index], moves^, cost, found, eqx)
            alignments_out[unsafe_offset=index] = Alignment(-score, score, cigar^, 0, columns, 0, rows)

    spread(spell_share, workers, workers)
    return settled^


def laned_local_scores(
    firsts: List[String], seconds: List[String], scoring: Scoring, threads: Int, scores_out: MutPointer[Int32, _]
) -> List[Bool]:
    """Writes the local score of every pair the lanes settle to `scores_out`, and returns which pairs those
    were, under a table of one match and one mismatch score (see `lanes.lane_local_scores`). The lanes pad
    past a sequence's end with two bytes the alphabet lacks, so a pair holding either is left for its own
    call."""
    var pairs = len(firsts)
    var workers = max(threads, 1)
    var uniform = uniform_scores(scoring)
    if not uniform:
        return List[Bool](length=pairs, fill=False)
    var codes_by_byte = code_table(scoring.alphabet)
    var pads = List[UInt8]()
    for byte in range(256):
        if codes_by_byte[byte] == UNKNOWN_SYMBOL and len(pads) < 2:
            pads.append(UInt8(byte))
    if len(pads) < 2:
        return List[Bool](length=pairs, fill=False)
    var refused = unknown_letters(firsts, seconds, scoring.alphabet, workers)
    var settled = refused.copy()
    var match_score, mismatch_score = uniform.value()
    _ = lane_local_scores(
        pairs,
        StringTexts.of(firsts),
        StringTexts.of(seconds),
        LocalCosts.symmetric(match_score, mismatch_score, Int(scoring.gaps.open), Int(scoring.gaps.extend)),
        (pads[0], pads[1]),
        workers,
        scores_out,
        settled.unsafe_ptr(),
    )
    unsettle(refused, settled)
    return settled^


@fieldwise_init
struct CodedBatch(Movable):
    """A batch in the alphabet's codes for the lanes (see `lanes.CodeTexts`): every first sequence, then every
    second one, back to back, and which pairs hold a letter outside the alphabet, whose codes are left zero."""

    var codes: List[UInt8]
    var starts: List[Int]
    """Where sequence `i` of the firsts begins, and sequence `i` of the seconds at `pairs + i`, and the end."""
    var refused: List[Bool]

    @staticmethod
    def of(firsts: List[String], seconds: List[String], alphabet: String, workers: Int) -> Self:
        """Every pair in `alphabet`'s codes, translated over `workers` threads."""
        var pairs = len(firsts)
        var starts = List[Int](capacity=2 * pairs + 1)
        var total = 0
        for text in firsts:
            starts.append(total)
            total += text.byte_length()
        for text in seconds:
            starts.append(total)
            total += text.byte_length()
        starts.append(total)
        var codes = List[UInt8](capacity=max(total, 1))
        codes.resize(unsafe_uninit_length=total)
        var refused = List[Bool](length=pairs, fill=False)
        var codes_by_byte = code_table(alphabet)
        var code_ptr = codes.unsafe_ptr()
        var start_ptr = starts.unsafe_ptr()
        var refused_ptr = refused.unsafe_ptr()

        def translate_share(
            share: Int,
        ) {
            imm firsts,
            imm seconds,
            imm codes_by_byte,
            imm pairs,
            imm workers,
            imm code_ptr,
            imm start_ptr,
            imm refused_ptr,
        }:
            """Translates the pairs of share `share`, flagging each that holds a letter outside the alphabet."""

            @always_inline
            def write(text: String, start: Int) {imm codes_by_byte, imm code_ptr} -> Bool:
                """Writes `text` in codes from `start` on, an unknown letter as zero, and whether all were known."""
                var known = True
                var target = code_ptr.unsafe_offset(start)
                var bytes = text.as_bytes()
                for position in range(len(bytes)):
                    var code = codes_by_byte[Int(bytes[position])]
                    known = known and code != UNKNOWN_SYMBOL
                    target[unsafe_offset=position] = 0 if code == UNKNOWN_SYMBOL else code
                return known

            for index in range(pairs * share // workers, pairs * (share + 1) // workers):
                var first_known = write(firsts[index], start_ptr[unsafe_offset=index])
                var second_known = write(seconds[index], start_ptr[unsafe_offset=pairs + index])
                refused_ptr[unsafe_offset=index] = not (first_known and second_known)

        spread(translate_share, workers, workers)
        return Self(codes^, starts^, refused^)

    def firsts(self) -> CodeTexts:
        """The first sequences, as the lanes read them."""
        return CodeTexts(
            self.codes.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
            self.starts.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
        )

    def seconds(self) -> CodeTexts:
        """The second sequences, as the lanes read them."""
        return CodeTexts(
            self.codes.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin](),
            self.starts.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]().unsafe_offset(len(self.refused)),
        )


def tabled_scores[
    mode: AlignmentMode
](
    firsts: List[String], seconds: List[String], scoring: Scoring, threads: Int, scores_out: MutPointer[Int32, _]
) -> List[Bool]:
    """Writes the global or local score of every pair the lanes settle to `scores_out`, and returns which pairs
    those were, under a table of more than one mismatch score with at most `TABLE_ENTRIES` entries, the
    sequences in the alphabet's codes.

    A global score is found as a cost, folding in the table's best score as `wavefront_penalties` folds a
    uniform table's match: a pair costs twice that best less its own score, a gap's letter that best less
    twice the extension, and a gap's opening twice the extension less the opening, all divided by their
    greatest common factor."""
    var pairs = len(firsts)
    var workers = max(threads, 1)
    var letters = scoring.alphabet_size()
    var none = List[Bool](length=pairs, fill=False)
    if letters < 2 or letters * letters > TABLE_ENTRIES:
        return none^
    var best, least = table_extremes(scoring.substitutions, letters)
    var open = Int(scoring.gaps.open)
    var extend = Int(scoring.gaps.extend)
    comptime if mode == AlignmentMode.GLOBAL:
        var opening = 2 * (extend - open)
        var extension = best - 2 * extend
        if opening < 0 or extension <= 0:
            return none^
        var scale = gcd(opening, extension)
        for cell in range(letters * letters):
            scale = gcd(scale, 2 * (best - Int(scoring.substitutions[cell])))
        var table = SIMD[DType.uint8, TABLE_ENTRIES](0)
        for cell in range(letters * letters):
            var cost = 2 * (best - Int(scoring.substitutions[cell])) // scale
            if cost > 255:
                return none^
            table[cell] = UInt8(cost)
        opening //= scale
        extension //= scale
        var batch = CodedBatch.of(firsts, seconds, scoring.alphabet, workers)
        var settled = batch.refused.copy()
        var penalties = Penalties(0, opening, extension, scale, best, 0, 0)
        var costs = LaneCosts.one_piece(0, opening, extension, opening, extension).tabled(letters, table)
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
        unsettle(batch.refused, settled)
        for index in range(pairs):
            if settled[index]:
                var letters_of_pair = firsts[index].byte_length() + seconds[index].byte_length()
                scores_out[unsafe_offset=index] = Int32(penalties.score(found[index].value(), letters_of_pair))
        return settled^
    else:
        var batch = CodedBatch.of(firsts, seconds, scoring.alphabet, workers)
        var settled = batch.refused.copy()
        # The two codes past the alphabet pad the lanes.
        _ = lane_local_scores(
            pairs,
            batch.firsts(),
            batch.seconds(),
            LocalCosts.symmetric(best, least, open, extend),
            (UInt8(letters), UInt8(letters + 1)),
            workers,
            scores_out,
            settled.unsafe_ptr(),
            letters,
            shuffled_table(scoring.substitutions, letters),
        )
        unsettle(batch.refused, settled)
        return settled^


def lane_scores[
    mode: AlignmentMode
](
    firsts: List[String], seconds: List[String], scoring: Scoring, threads: Int, scores_out: MutPointer[Int32, _]
) -> List[Bool]:
    """Writes the score of every pair some lane kernel settles to `scores_out`, and returns which pairs those
    were. A table of one match and one mismatch score takes the wavefront's lanes for a global score, whose
    cost folds the match in (see `wavefront_penalties`), or the local lanes; a small table of several
    mismatch scores takes the lanes over its codes (see `tabled_scores`)."""
    if not uniform_scores(scoring):
        return tabled_scores[mode](firsts, seconds, scoring, threads, scores_out)
    comptime if mode == AlignmentMode.LOCAL:
        return laned_local_scores(firsts, seconds, scoring, threads, scores_out)
    var penalties = scoring.penalties()
    if not penalties:
        return List[Bool](length=len(firsts), fill=False)
    var found = penalties.value()
    return laned_scores(firsts, seconds, scoring.alphabet, LaneCosts.of_penalties(found), found, threads, scores_out)


# endregion Host Lanes

# region Batches

comptime CHUNKS_PER_THREAD = 8
"""How many chunks per thread a host batch is cut into, so that a thread that draws long pairs does not hold
up the rest."""


def chunk_count(pairs: Int, threads: Int) -> Int:
    """How many contiguous chunks a host batch of `pairs` is cut into for `threads` threads."""
    return min(pairs, max(threads, 1) * CHUNKS_PER_THREAD)


def batch_scores[
    mode: AlignmentMode
](firsts: List[String], seconds: List[String], scoring: Scoring, placement: Optional[Placement] = None) raises -> List[
    Int32
]:
    """The best global or local score of every pair, on the device `placement` names (see `host_batch_scores`,
    `device_batch_scores`)."""
    var resolved = placement.or_else(Placement.default())
    var pairs = pair_count(firsts, seconds)
    if pairs == 0:
        return List[Int32]()
    if resolved.device == Device.GPU:
        return device_batch_scores[mode](firsts, seconds, scoring, resolved)
    return host_batch_scores[mode](firsts, seconds, scoring, resolved.threads)


def host_batch_scores[
    mode: AlignmentMode
](firsts: List[String], seconds: List[String], scoring: Scoring, threads: Int) raises -> List[Int32]:
    """Every pair's score on the host: as many as the lanes take many at once (see `lane_scores`), and each of
    the rest whole on one thread, in contiguous chunks over `threads` threads. A pair that raises is scored
    again afterwards, in order, so the batch raises the error a serial loop would have raised first."""
    var pairs = len(firsts)
    var results = List[Int32](length=pairs, fill=0)
    var out = results.unsafe_ptr()
    var settled = lane_scores[mode](firsts, seconds, scoring, threads, out)
    var failed = List[Bool](length=pairs, fill=False)
    var settled_ptr = settled.unsafe_ptr()
    var failed_ptr = failed.unsafe_ptr()
    var chunks = chunk_count(pairs, threads)

    def score_chunk(chunk: Int) {imm}:
        """Scores the unsettled pairs of chunk `chunk`, flagging each that raises."""
        for index in range(pairs * chunk // chunks, pairs * (chunk + 1) // chunks):
            if settled_ptr[unsafe_offset=index]:
                continue
            try:
                out[unsafe_offset=index] = host_text_score[mode](firsts[index], seconds[index], scoring)
            except:
                failed_ptr[unsafe_offset=index] = True

    spread(score_chunk, chunks, max(threads, 1))
    for index in range(pairs):
        if failed[index]:
            results[index] = host_text_score[mode](firsts[index], seconds[index], scoring)
    return results^


def device_batch_scores[
    mode: AlignmentMode
](firsts: List[String], seconds: List[String], scoring: Scoring, placement: Placement) raises -> List[Int32]:
    """Every pair's score on the device. A global score first tries the band its cost proves, a thread a pair
    (see `band_groups`); the pairs that band does not prove sweep their whole matrix, several a warp wherever
    every second sequence fits a shape (see `score_groups`), and otherwise in a batch launch, a warp a pair,
    where a block can carry the pair, and tiled alone where it cannot."""
    var pairs = len(firsts)
    var scope = DeviceScope(placement.gpu_id)
    var results = List[Int32](length=pairs, fill=0)
    var remaining = List[Int](capacity=pairs)
    for index in range(pairs):
        remaining.append(index)
    comptime if mode == AlignmentMode.GLOBAL:
        var banded = banded_scores(
            scope,
            firsts,
            seconds,
            scoring.alphabet,
            scoring.substitutions,
            scoring.gaps,
            placement.threads,
            placement.gpu_id,
        )
        if banded:
            ref proven = banded.value()
            results = proven[0].copy()
            remaining.clear()
            for index in range(pairs):
                if not proven[1][index]:
                    remaining.append(index)
            if len(remaining) == 0:
                return results^
    var grouped = grouped_scores[mode](
        scope,
        firsts,
        seconds,
        remaining,
        scoring.alphabet,
        scoring.substitutions,
        scoring.gaps,
        placement.threads,
        placement.gpu_id,
    )
    if grouped:
        var scored = grouped.take()
        for slot in range(len(remaining)):
            results[remaining[slot]] = scored[slot]
        return results^
    var band = band_length(scope.specs)
    var launched = List[Int]()
    for index in remaining:
        if batched(firsts[index].byte_length(), band):
            launched.append(index)
        else:
            var first_codes = translate(firsts[index], scoring.alphabet)
            var second_codes = translate(seconds[index], scoring.alphabet)
            results[index] = device_score[mode](
                scope, first_codes, second_codes, scoring.substitutions, scoring.alphabet_size(), scoring.gaps
            )
    if len(launched) > 0:
        var tape = PairTape.of(firsts, seconds, launched, scoring.alphabet, placement.threads)
        var scored = device_scores[mode](
            scope, tape.codes, tape.starts, scoring.substitutions, scoring.alphabet_size(), scoring.gaps
        )
        for slot in range(len(launched)):
            results[launched[slot]] = scored[slot]
    return results^


def launch_cuts(firsts: List[String], seconds: List[String], launched: List[Int], largest: Int) -> List[Int]:
    """Where to cut `launched` into launches, each within the device's `largest` allocation: its pairs' recorded
    decisions, and their gapped rows at its widest pair's length each. The cuts run from 0 to `len(launched)`,
    and every launch holds at least one pair. A single launch asked 16 GB at once of a batch of 700,000 short
    reads, and as much of each pair as one long pair among many short ones needed."""
    var cuts: List[Int] = [0]
    var recorded = 0
    var widest = 1
    for slot in range(len(launched)):
        var rows = firsts[launched[slot]].byte_length()
        var columns = seconds[launched[slot]].byte_length()
        var more = launch_bytes(rows, columns)
        var wider = max(widest, rows + columns)
        var held = slot - cuts[len(cuts) - 1]
        if held > 0 and (recorded + more > largest or (held + 1) * wider > largest):
            cuts.append(slot)
            recorded = 0
            wider = max(rows + columns, 1)
        recorded += more
        widest = wider
    cuts.append(len(launched))
    return cuts^


def device_batch_alignments[
    mode: AlignmentMode
](
    firsts: List[String], seconds: List[String], scoring: Scoring, placement: Placement, stored_cells: Int
) raises -> List[GappedAlignment]:
    """An optimal global or local alignment of every pair on the device. Each pair a block can carry and whose
    decisions fit `stored_cells` joins a batch launch, cut into as many as the device's buffers need (see
    `launch_cuts`); every other pair is aligned alone (see `device_pair_alignment`)."""
    var pairs = pair_count(firsts, seconds)
    var results = List[GappedAlignment](capacity=pairs)
    for _ in range(pairs):
        results.append(GappedAlignment(0, String(), String()))
    if pairs == 0:
        return results^
    # A batch of one is a single pair, whatever way it came, and is bounded as one.
    var limit = stored_cells if pairs > 1 else min(stored_cells, DEVICE_STORED_CELLS)
    var scope = DeviceScope(placement.gpu_id)
    var band = band_length(scope.specs)
    var launched = List[Int]()
    for index in range(pairs):
        var rows = firsts[index].byte_length()
        if rows * seconds[index].byte_length() <= limit and batched(rows, band):
            launched.append(index)
            continue
        var first_codes = translate(firsts[index], scoring.alphabet)
        var second_codes = translate(seconds[index], scoring.alphabet)
        results[index] = device_pair_alignment[mode](scope, first_codes, second_codes, scoring, stored_cells, placement)
    var cuts = launch_cuts(firsts, seconds, launched, scope.specs.largest_allocation)
    for launch in range(len(cuts) - 1):
        var chosen = List[Int](launched[cuts[launch] : cuts[launch + 1]])
        if len(chosen) == 0:
            continue
        var tape = PairTape.of(firsts, seconds, chosen, scoring.alphabet, placement.threads)
        var aligned = device_alignments[mode](
            scope, tape.codes, tape.starts, scoring.substitutions, scoring.alphabet, scoring.gaps
        )
        # Moved out rather than copied, last first, as `pop` gives them.
        for slot in reversed(range(len(chosen))):
            results[chosen[slot]] = aligned.pop()
    return results^


# endregion Batches
