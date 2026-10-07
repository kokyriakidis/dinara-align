"""
The entry points a caller uses, each deciding which search serves a pair so the caller never has to.

Under `Costs`, which an alignment minimizes, `distance` finds the least cost and `align` an optimal
alignment as a CIGAR, for any `Mode`; `distances` and `alignments` take a batch over every thread.
Unit costs, global or with the query found inside or at the start of the reference, take the
bit-parallel band doubling of A*PA2 (see `edit_distance` and `edit_search`); every other case, and a
pair holding more symbols than it takes, the gap-affine wavefront from both ends (see `gap_affine`),
which gives the same alignment at the same costs. The modes with a match score maximize a score: an
extension, and a global alignment with a reward, by the wavefront; a local alignment and other free
ends with a reward by sweep (see `scored`).

Under a `Scoring`, an alphabet's substitution table and gap scores, which an alignment maximizes,
`score` and `align` take a global or a local alignment, on the host or the device (see `scoring`).
"""

from std.atomic import Atomic

from max.algorithm import parallelize

from .alignment import AlignmentMode, GappedAlignment
from .common import Placement, hardware_threads
from .edit_distance import edit_cigar, edit_distance
from .edit_search import edit_search
from .errors import AlignmentError, ErrorKind
from .cigar import cigar_matches, cigar_runs, joined_cigar
from .scored import global_rewarded, local_alignment, rewarded_alignment
from .gap_affine import (
    AffineCigar,
    EndsFree,
    cigar_within,
    extension_of,
    extension_penalties,
    outside,
    penalties_of,
    wavefront_distance,
)
from .modes import Alignment, Anchor, Band, Costs, Mode, Ties
from .scoring import (
    STORED_MATRIX_BUDGET,
    Scoring,
    align_with,
    alignments_with,
    paired_length,
    score_with,
    scores_with,
)


# region Costs


def distance(
    reference: String, query: String, costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, band: Band = Band()
) raises AlignmentError -> Int:
    """The least cost of aligning `query` to `reference` as `mode` asks, every move inside `band` (see
    `Band`), with no alignment; raises when no alignment fits the band.

    Unit costs, globally, take A*PA2's band doubling: guess a bound, sweep only the cells a path within
    it could cross, raise the guess until the answer fits under it. Other costs take the wavefront from
    both ends, keeping only a few costs' fronts, so a few rows of memory however long the pair.
    """
    var found = cost_within(reference, query, costs, mode, band, Int.MAX)
    if not found:
        raise outside(band)
    return found.value()


def distance(
    reference: String,
    query: String,
    costs: Costs = Costs.edit(),
    mode: Mode = Mode.GLOBAL,
    *,
    max_cost: Int,
    band: Band = Band(),
) raises AlignmentError -> Optional[Int]:
    """`distance`, or None when it would pass `max_cost` or no alignment fits `band`: the searches stop
    once each has grown to about half of `max_cost` without the two meeting within it, so a pair far
    over costs a fraction of its full search."""
    return cost_within(reference, query, costs, mode, band, max_cost)


def align(
    reference: String,
    query: String,
    costs: Costs = Costs.edit(),
    mode: Mode = Mode.GLOBAL,
    *,
    band: Band = Band(),
    ties: Ties = Ties.LEFT,
    extended: Bool = True,
) raises AlignmentError -> Alignment:
    """An optimal alignment of `query` to `reference` as `mode` asks, every move inside `band`, as a
    CIGAR with `=` and `X`, or with `extended` false `M` for both (see `Alignment`); raises when no
    alignment fits the band.

    Of several equally good alignments the CIGAR is always the one `ties` names (see `Ties`): by
    default every edit as far left as it goes, indels placed as minimap2 places them, or with
    `Ties.RIGHT` as far right, WFA2-lib's CIGAR byte for byte, whichever search found the cost.

    Every byte is a symbol matching only itself, so DNA in either case, or any other text, needs no
    alphabet. The wavefront's work grows with the square of the cost rather than with the matrix, and
    its memory stays bounded.
    """
    var found = aligned_within(reference, query, costs, mode, band, Int.MAX, ties, extended)
    if not found:
        raise outside(band)
    return found.take()


def align(
    reference: String,
    query: String,
    costs: Costs = Costs.edit(),
    mode: Mode = Mode.GLOBAL,
    *,
    max_cost: Int,
    band: Band = Band(),
    ties: Ties = Ties.LEFT,
    extended: Bool = True,
) raises AlignmentError -> Optional[Alignment]:
    """`align`, or None when the cost would pass `max_cost` or no alignment fits `band`, found as
    `distance` finds that, with no fronts traced. A mode with a match score, which maximizes a score,
    takes no cap."""
    return aligned_within(reference, query, costs, mode, band, max_cost, ties, extended)


def free_ends(mode: Mode, columns: Int, rows: Int) -> EndsFree:
    """A mode's free letters, none past its sequence's length."""
    return EndsFree(
        min(mode.reference_start, columns),
        min(mode.reference_end, columns),
        min(mode.query_start, rows),
        min(mode.query_end, rows),
    )


def cost_within(
    reference: String, query: String, costs: Costs, mode: Mode, band: Band, max_cost: Int
) raises AlignmentError -> Optional[Int]:
    """The least cost, or None when it would pass `max_cost` (`Int.MAX` for no cap) or none fits `band`."""
    if mode.is_scored():
        raise AlignmentError(
            ErrorKind.INVALID_ARGUMENT, "a mode with a match score maximizes a score, which `align` finds"
        )
    var columns = reference.byte_length()
    var rows = query.byte_length()
    var ends = free_ends(mode, columns, rows)
    if swept_by_bits(costs, ends, band, max_cost, columns, rows):
        try:
            var scale = costs.unit_scale()
            if ends.first_begin == 0 and ends.first_end == 0:
                return edit_distance(reference, query) * scale
            return edit_search(query, reference, ends.first_begin == 0).distance * scale
        except error:
            # More symbols than the sweep takes: the wavefront takes any.
            if error.kind != ErrorKind.UNKNOWN_SYMBOL:
                raise error
    var penalties = penalties_of(costs)
    if max_cost < 0:
        return None
    var ceiling = max_cost // penalties.scale
    var cost: Int
    if costs.pieces() == 2:
        cost = wavefront_distance[2](reference.as_bytes(), query.as_bytes(), penalties, ceiling, ends, band)
    else:
        cost = wavefront_distance[1](reference.as_bytes(), query.as_bytes(), penalties, ceiling, ends, band)
    if cost < 0:
        return None
    return cost * penalties.scale


def aligned_within(
    reference: String,
    query: String,
    costs: Costs,
    mode: Mode,
    band: Band,
    max_cost: Int,
    ties: Ties,
    extended: Bool,
) raises AlignmentError -> Optional[Alignment]:
    """An optimal alignment, or None when its cost would pass `max_cost` (`Int.MAX` for no cap) or none
    fits `band`: the least costly one, or for a mode that maximizes a score the best-scoring one."""
    if not mode.is_scored():
        return least_costly(reference, query, costs, mode, band, max_cost, ties, extended)
    if max_cost != Int.MAX:
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a mode with a match score takes no cost cap")
    _ = penalties_of(costs)
    return best_scoring(reference, query, costs, mode, band, ties, extended)


def swept_by_bits(costs: Costs, ends: EndsFree, band: Band, max_cost: Int, columns: Int, rows: Int) -> Bool:
    """Whether the bit-parallel sweep serves a pair: unit costs, no cap and no band, both sequences
    holding a letter, and free ends of none, the query inside the reference, or at its start."""
    if costs.unit_scale() == 0 or max_cost != Int.MAX or not band.covers(columns, rows) or columns == 0 or rows == 0:
        return False
    if ends.second_begin != 0 or ends.second_end != 0:
        return False
    if ends.first_begin == 0:
        return ends.first_end == 0 or ends.first_end == columns
    return ends.first_begin == columns and ends.first_end == columns


def least_costly(
    reference: String,
    query: String,
    costs: Costs,
    mode: Mode,
    band: Band,
    max_cost: Int,
    ties: Ties,
    extended: Bool,
) raises AlignmentError -> Optional[Alignment]:
    """The least costly alignment with `mode`'s free ends: by the bit-parallel sweep where it serves the
    pair (see `swept_by_bits`), else by the wavefront, the same alignment at the same costs."""
    var columns = reference.byte_length()
    var rows = query.byte_length()
    var ends = free_ends(mode, columns, rows)
    if swept_by_bits(costs, ends, band, max_cost, columns, rows):
        try:
            var scale = costs.unit_scale()
            if ends.first_begin == 0 and ends.first_end == 0:
                var whole = edit_cigar(reference, query, extended, ties)
                var cost = whole.distance * scale
                return Alignment(cost, -cost, whole.cigar, 0, columns, 0, rows)
            var hit = edit_search(query, reference, ends.first_begin == 0)
            var part = String(StringSlice(unsafe_from_utf8=reference.as_bytes()[hit.start : hit.end]))
            var found = edit_cigar(part, query, extended, ties)
            var cost = found.distance * scale
            return Alignment(cost, -cost, found.cigar, hit.start, hit.end, 0, rows)
        except error:
            if error.kind != ErrorKind.UNKNOWN_SYMBOL:
                raise error
    var penalties = penalties_of(costs)
    if max_cost < 0:
        return None
    var found: Optional[AffineCigar]
    if costs.pieces() == 2:
        found = cigar_within[2](reference, query, penalties, extended, max_cost // penalties.scale, ends, band, ties)
    else:
        found = cigar_within[1](reference, query, penalties, extended, max_cost // penalties.scale, ends, band, ties)
    if not found:
        return None
    return trimmed(found.take(), ends, columns, rows)


def best_scoring(
    reference: String, query: String, costs: Costs, mode: Mode, band: Band, ties: Ties, extended: Bool
) raises AlignmentError -> Alignment:
    """The best-scoring alignment for a mode with a match score: an extension by the wavefront from its
    anchor; a global alignment, its letters fixed so the reward folds into the costs, by the wavefront
    too; a local alignment or other free ends by sweep (see `scored`)."""
    var columns = reference.byte_length()
    var rows = query.byte_length()
    if mode.kind == Mode.EXTENSION:
        return extended_alignment(reference, query, costs, mode, band, ties, extended)
    if mode.kind == Mode.ENDS and mode.is_global():
        var whole = global_rewarded(reference, query, costs, mode.match_score, band, ties, extended)
        var score = mode.match_score * cigar_matches(reference, query, whole[1]) - whole[0]
        return Alignment(whole[0], score, whole[1], 0, columns, 0, rows)
    if mode.match_score <= 0:
        raise AlignmentError(
            ErrorKind.INVALID_SCORING, "a local alignment under Costs needs a match that earns: Mode.local"
        )
    if not band.covers(columns, rows):
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a sweep takes no band: a local alignment or free ends")
    if mode.kind == Mode.SMITH_WATERMAN:
        return local_alignment(reference, query, costs, mode.match_score, ties, extended)
    return rewarded_alignment(reference, query, costs, mode.match_score, free_ends(mode, columns, rows), ties, extended)


def extended_alignment(
    reference: String, query: String, costs: Costs, mode: Mode, band: Band, ties: Ties, extended: Bool
) raises AlignmentError -> Alignment:
    """The best extension from `mode`'s anchor (see `Mode.extension`)."""
    var two = costs.pieces() == 2
    var penalties = extension_penalties(
        mode.match_score,
        costs.mismatch,
        costs.opening,
        costs.extension,
        costs.opening2 if two else 0,
        costs.extension2 if two else 0,
    )
    var found = extension_of[2](
        reference, query, penalties, extended, mode.anchor, band, ties
    ) if two else extension_of[1](reference, query, penalties, extended, mode.anchor, band, ties)
    var columns = reference.byte_length()
    var rows = query.byte_length()
    var cost = mode.match_score * found.matches - found.score
    if mode.anchor == Anchor.END:
        return Alignment(
            cost,
            found.score,
            found.cigar,
            columns - found.first_length,
            columns,
            rows - found.second_length,
            rows,
        )
    return Alignment(cost, found.score, found.cigar, 0, found.first_length, 0, found.second_length)


def trimmed(var found: AffineCigar, ends: EndsFree, columns: Int, rows: Int) -> Alignment:
    """An alignment whose free letters, `D` and `I` runs at either end of the wavefront's CIGAR, move
    into the spans: as many as the mode frees, the rest of a run left in the CIGAR, as it is paid."""
    var runs = cigar_runs(found.cigar)
    var letters = runs[0].copy()
    var lengths = runs[1].copy()
    var reference_start = 0
    var reference_end = columns
    var query_start = 0
    var query_end = rows
    comptime DELETED = UInt8(ord("D"))
    comptime INSERTED = UInt8(ord("I"))
    var count = len(letters)
    if count > 0:
        var free = 0
        if letters[0] == DELETED:
            free = min(lengths[0], ends.first_begin)
            reference_start = free
        elif letters[0] == INSERTED:
            free = min(lengths[0], ends.second_begin)
            query_start = free
        lengths[0] -= free
        var last = count - 1
        free = 0
        if letters[last] == DELETED:
            free = min(lengths[last], ends.first_end)
            reference_end = columns - free
        elif letters[last] == INSERTED:
            free = min(lengths[last], ends.second_end)
            query_end = rows - free
        lengths[last] -= free
    return Alignment(
        found.cost, -found.cost, joined_cigar(letters, lengths), reference_start, reference_end, query_start, query_end
    )


# endregion Costs

# region Scoring


def scoring_mode(mode: Mode) raises AlignmentError -> AlignmentMode:
    """The recurrence a `Scoring` runs for `mode`: global or local alone."""
    if mode.kind == Mode.SMITH_WATERMAN:
        return AlignmentMode.LOCAL
    if mode.match_score > 0:
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a Scoring's table holds what a match earns")
    if mode.is_global():
        return AlignmentMode.GLOBAL
    raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a Scoring aligns globally or locally alone")


def score(
    reference: String, query: String, scoring: Scoring, mode: Mode = Mode.GLOBAL, placement: Optional[Placement] = None
) raises -> Int32:
    """The optimal score under `scoring`, `Mode.GLOBAL` or `Mode.LOCAL`, alone, in two rows of memory on
    either device (see `scoring.score_with`)."""
    if scoring_mode(mode) == AlignmentMode.LOCAL:
        return score_with[AlignmentMode.LOCAL](reference, query, scoring, placement)
    return score_with[AlignmentMode.GLOBAL](reference, query, scoring, placement)


def align(
    reference: String,
    query: String,
    scoring: Scoring,
    mode: Mode = Mode.GLOBAL,
    placement: Optional[Placement] = None,
    stored_budget: Int = STORED_MATRIX_BUDGET,
) raises -> GappedAlignment:
    """The optimal score under `scoring` and the two gapped rows that earn it, the reference's first,
    also read as a CIGAR (see `GappedAlignment.cigar`): both sequences whole for `Mode.GLOBAL`,
    Needleman-Wunsch, or the best-scoring window of each for `Mode.LOCAL`, Smith-Waterman."""
    if scoring_mode(mode) == AlignmentMode.LOCAL:
        return align_with[AlignmentMode.LOCAL](reference, query, scoring, placement, stored_budget)
    return align_with[AlignmentMode.GLOBAL](reference, query, scoring, placement, stored_budget)


def scores(
    references: List[String],
    queries: List[String],
    scoring: Scoring,
    mode: Mode = Mode.GLOBAL,
    placement: Optional[Placement] = None,
) raises -> List[Int32]:
    """`score` for every pair; on the device, every pair one block can carry goes out in one launch."""
    if scoring_mode(mode) == AlignmentMode.LOCAL:
        return scores_with[AlignmentMode.LOCAL](references, queries, scoring, placement)
    return scores_with[AlignmentMode.GLOBAL](references, queries, scoring, placement)


def alignments(
    references: List[String],
    queries: List[String],
    scoring: Scoring,
    mode: Mode = Mode.GLOBAL,
    placement: Optional[Placement] = None,
    stored_budget: Int = STORED_MATRIX_BUDGET,
) raises -> List[GappedAlignment]:
    """`align` for every pair; on the device, every pair both bounds admit goes out in one launch."""
    if scoring_mode(mode) == AlignmentMode.LOCAL:
        return alignments_with[AlignmentMode.LOCAL](references, queries, scoring, placement, stored_budget)
    return alignments_with[AlignmentMode.GLOBAL](references, queries, scoring, placement, stored_budget)


# endregion Scoring

# region Batches


def longest_first(references: List[String], queries: List[String]) -> List[Int]:
    """The pairs' indices, the longest pair first: taken in that order by whichever thread is free, the
    long pairs start first and the short ones fill in around them, so none is left alone at the end
    holding up the rest."""
    comptime INDEX_BITS = 25
    var pairs = len(references)
    var order = List[Int](capacity=pairs)
    if pairs >= 1 << INDEX_BITS:
        for index in range(pairs):
            order.append(index)
        return order^
    var keys = List[Int](capacity=pairs)
    for index in range(pairs):
        var length = references[index].byte_length() + queries[index].byte_length()
        keys.append((length << INDEX_BITS) | index)
    sort(keys)
    for slot in range(pairs - 1, -1, -1):
        order.append(keys[slot] & ((1 << INDEX_BITS) - 1))
    return order^


def distances(
    references: List[String],
    queries: List[String],
    costs: Costs = Costs.edit(),
    mode: Mode = Mode.GLOBAL,
    *,
    band: Band = Band(),
    threads: Optional[Int] = None,
) raises AlignmentError -> List[Int]:
    """Every pair's `distance`, the pairs spread over `threads` threads, every thread this process may
    use by default.

    The pairs are independent, so each runs on one thread start to finish, each thread taking the
    next pair of the batch, longest first, as soon as it is free (see `longest_first`). A pair that
    fails raises, after the rest, the same error a serial loop would have raised first.
    """
    var pairs = paired_length(references, queries)
    var results = List[Int](length=pairs, fill=0)
    if pairs == 0:
        return results^
    var workers = max(threads.or_else(hardware_threads()), 1)
    var out = results.unsafe_ptr()
    var failed = List[Bool](length=pairs, fill=False)
    var flags = failed.unsafe_ptr()
    var order = longest_first(references, queries)
    var taken = Atomic[Int64](0)

    def distance_worker(
        worker: Int,
    ) {mut taken, imm order, imm references, imm queries, imm out, imm flags, imm pairs, imm costs, imm mode, imm band}:
        while True:
            var dealt = Int(taken.fetch_add(1))
            if dealt >= pairs:
                return
            var index = order[dealt]
            try:
                out[unsafe_offset=index] = distance(references[index], queries[index], costs, mode, band=band)
            except:
                flags[unsafe_offset=index] = True

    parallelize(distance_worker, min(workers, pairs), min(workers, pairs))
    for index in range(pairs):
        if failed[index]:
            results[index] = distance(references[index], queries[index], costs, mode, band=band)
    return results^


def alignments(
    references: List[String],
    queries: List[String],
    costs: Costs = Costs.edit(),
    mode: Mode = Mode.GLOBAL,
    *,
    band: Band = Band(),
    ties: Ties = Ties.LEFT,
    extended: Bool = True,
    threads: Optional[Int] = None,
) raises AlignmentError -> List[Alignment]:
    """Every pair's `align`, the pairs spread over threads as `distances` spreads them."""
    var pairs = paired_length(references, queries)
    var results = List[Alignment](capacity=pairs)
    for _ in range(pairs):
        results.append(Alignment(0, 0, String(), 0, 0, 0, 0))
    if pairs == 0:
        return results^
    var workers = max(threads.or_else(hardware_threads()), 1)
    var out = results.unsafe_ptr()
    var failed = List[Bool](length=pairs, fill=False)
    var flags = failed.unsafe_ptr()
    var order = longest_first(references, queries)
    var taken = Atomic[Int64](0)

    def alignment_worker(
        worker: Int,
    ) {
        mut taken,
        imm order,
        imm references,
        imm queries,
        imm out,
        imm flags,
        imm pairs,
        imm costs,
        imm mode,
        imm band,
        imm ties,
        imm extended,
    }:
        while True:
            var dealt = Int(taken.fetch_add(1))
            if dealt >= pairs:
                return
            var index = order[dealt]
            try:
                out[unsafe_offset=index] = align(
                    references[index], queries[index], costs, mode, band=band, ties=ties, extended=extended
                )
            except:
                flags[unsafe_offset=index] = True

    parallelize(alignment_worker, min(workers, pairs), min(workers, pairs))
    for index in range(pairs):
        if failed[index]:
            results[index] = align(
                references[index], queries[index], costs, mode, band=band, ties=ties, extended=extended
            )
    return results^


# endregion Batches
