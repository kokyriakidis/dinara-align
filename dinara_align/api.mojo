# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The entry points a caller uses, each deciding which search serves a pair so the caller never has to.

Under `Costs`, which an alignment minimizes, `distance` finds the least cost and `align` an optimal
alignment as a CIGAR, for any `Mode`; `distances` and `alignments` take a batch, many pairs at once.
Every call runs on the caller's own thread unless asked for more (see `distances`), and keeps no state
between calls, so an application may call it from as many threads as it likes; an `Aligner` keeps one
thread's memory from call to call.
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
from std.bit import count_leading_zeros


from .alignment import AlignmentMode, GappedAlignment
from .common import Device, DeviceScope, Placement, next_share, spread
from .device_edit import MAX_PATTERN_WORDS, device_edit_distances
from .edit_distance import edit_cigar, edit_distance
from .edit_search import edit_search
from .errors import AlignmentError, ErrorKind
from .cigar import cigar_matches, cigar_runs, joined_cigar
from .search import Hit, local_scores_by_block
from .scored import (
    FROM_EDGE,
    LocalScores,
    end_of,
    global_rewarded,
    local_alignment,
    local_scores as local_scores_of,
    rewarded_alignment,
    swept,
)
from .lanes import LaneCosts, StringTexts, lane_alignments, lane_distances, lane_free_alignments
from .gap_affine import (
    AffineCigar,
    DEFAULT_MAX_MEMORY,
    SearchSpace,
    EndsFree,
    KEPT_BYTES,
    Penalties,
    Spanned,
    cigar_of,
    cigar_within,
    extend,
    extension_of,
    free_ends_alignment,
    rewarded_penalties,
    outside,
    penalties_of,
    wavefront_distance,
)
from .modes import Alignment, Anchor, Band, Costs, Mode, Ties
from .scoring import (
    cells_within,
    Scoring,
    alignments_with,
    as_alignment,
    paired_length,
    scores_with,
    laned_alignments,
    scoring_alignment,
    scoring_score,
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
    eqx: Bool = True,
    max_memory: Int = DEFAULT_MAX_MEMORY,
) raises AlignmentError -> Alignment:
    """An optimal alignment of `query` to `reference` as `mode` asks, every move inside `band`, as a
    CIGAR with `=` and `X`, or with `eqx` false `M` for both (see `Alignment`); raises when no
    alignment fits the band.

    Of several equally good alignments the CIGAR is always the one `ties` names (see `Ties`): by
    default every edit as far left as it goes, indels placed as minimap2 places them, or with
    `Ties.RIGHT` as far right, WFA2-lib's CIGAR byte for byte, whichever search found the cost.

    Every byte is a symbol matching only itself, so DNA in either case, or any other text, needs no
    alphabet. The wavefront's work grows with the square of the cost rather than with the matrix, and
    its memory stays bounded: the fronts it keeps for the traceback, the bulk of it, stay within
    `max_memory` bytes, past which the pair is split where an optimal path crosses and each piece
    aligned alone (see `gap_affine.solve`), the cost still the least, the tie rule followed within each
    piece. The sweeps and the bit-parallel search keep a few rows, or a band's edges, whatever the cap.
    """
    var found = aligned_within(reference, query, costs, mode, band, Int.MAX, ties, eqx, max_memory // KEPT_BYTES)
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
    eqx: Bool = True,
    max_memory: Int = DEFAULT_MAX_MEMORY,
) raises AlignmentError -> Optional[Alignment]:
    """`align`, or None when the cost would pass `max_cost` or no alignment fits `band`, found as
    `distance` finds that, with no fronts traced. A mode with a match score, which maximizes a score,
    takes no cap."""
    return aligned_within(reference, query, costs, mode, band, max_cost, ties, eqx, max_memory // KEPT_BYTES)


struct Aligner(Movable):
    """One thread's aligner: `distance` and `align` as the functions of those names give them, the memory
    their searches take kept from call to call, so a loop of calls takes none once it is warm. An
    application calling from many threads keeps one a thread; one never crosses threads, and the library
    keeps no state of its own, so any number of them work at once."""

    var space: SearchSpace

    def __init__(out self):
        """An aligner holding no memory yet: its first call takes what it needs."""
        self.space = SearchSpace()

    def distance(
        mut self,
        reference: String,
        query: String,
        costs: Costs = Costs.edit(),
        mode: Mode = Mode.GLOBAL,
        *,
        band: Band = Band(),
    ) raises AlignmentError -> Int:
        """`distance`, through this aligner's memory."""
        var found = cost_within(reference, query, costs, mode, band, Int.MAX, self.space)
        if not found:
            raise outside(band)
        return found.value()

    def distance(
        mut self,
        reference: String,
        query: String,
        costs: Costs = Costs.edit(),
        mode: Mode = Mode.GLOBAL,
        *,
        max_cost: Int,
        band: Band = Band(),
    ) raises AlignmentError -> Optional[Int]:
        """`distance` under a cap, through this aligner's memory."""
        return cost_within(reference, query, costs, mode, band, max_cost, self.space)

    def align(
        mut self,
        reference: String,
        query: String,
        costs: Costs = Costs.edit(),
        mode: Mode = Mode.GLOBAL,
        *,
        band: Band = Band(),
        ties: Ties = Ties.LEFT,
        eqx: Bool = True,
        max_memory: Int = DEFAULT_MAX_MEMORY,
    ) raises AlignmentError -> Alignment:
        """`align`, through this aligner's memory: the same alignment, the CIGAR its `ties` picks."""
        var found = aligned_within(
            reference, query, costs, mode, band, Int.MAX, ties, eqx, max_memory // KEPT_BYTES, self.space
        )
        if not found:
            raise outside(band)
        return found.take()

    def align(
        mut self,
        reference: String,
        query: String,
        costs: Costs = Costs.edit(),
        mode: Mode = Mode.GLOBAL,
        *,
        max_cost: Int,
        band: Band = Band(),
        ties: Ties = Ties.LEFT,
        eqx: Bool = True,
        max_memory: Int = DEFAULT_MAX_MEMORY,
    ) raises AlignmentError -> Optional[Alignment]:
        """`align` under a cap, through this aligner's memory."""
        return aligned_within(
            reference, query, costs, mode, band, max_cost, ties, eqx, max_memory // KEPT_BYTES, self.space
        )


def score(
    reference: String, query: String, costs: Costs, mode: Mode = Mode.GLOBAL, *, band: Band = Band()
) raises AlignmentError -> Int:
    """The best score `align` would return, with no alignment traced: for a mode with a match score its
    matches' reward less its costs, else minus the least cost, `distance`'s. A local alignment or free
    ends with a reward take the sweep alone, an extension its search alone, and a global alignment
    with a reward the wavefront's cost with the reward folded in, so each skips the traceback."""
    if not mode.is_scored():
        return -distance(reference, query, costs, mode, band=band)
    var columns = reference.byte_length()
    var rows = query.byte_length()
    var two = costs.pieces() == 2
    _ = penalties_of(costs)
    var penalties = rewarded_penalties(mode.match_score, costs)
    if mode.kind == Mode.EXTENSION:
        if not band.holds(0):
            raise outside(band)
        var drop_extension = cheapest_extension(costs)
        var at_end = mode.anchor == Anchor.END
        if mode.end_bonus > 0 and not band.covers(columns, rows):
            raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "an end bonus takes no band")
        var found = extend[2](
            reference.as_bytes(),
            query.as_bytes(),
            penalties,
            band,
            at_end,
            -1,
            mode.zdrop,
            drop_extension,
            mode.end_bonus,
        ) if two else extend[1](
            reference.as_bytes(),
            query.as_bytes(),
            penalties,
            band,
            at_end,
            -1,
            mode.zdrop,
            drop_extension,
            mode.end_bonus,
        )
        var best = penalties.score(found[0], found[1] + found[2])
        # The end bonus prefers the best extension reaching the query's far end, unless the Z-drop gave up.
        if mode.end_bonus > 0 and not found[3]:
            var reaching: Int
            if found[5] >= 0:
                reaching = penalties.score(found[4], found[5] + rows)
            elif columns > 0 and rows > 0 and penalties.reward > 0:
                # No front point reached the query's end before the search could stop: none comes close.
                return best
            else:
                # An empty side or no reward, where the search does not run: the free ends' own.
                reaching = score(reference, query, costs, mode.reaching_end(), band=band)
            if reaching + mode.end_bonus > best:
                return reaching
        return best
    if mode.kind == Mode.ENDS and mode.is_global():
        var cost = wavefront_distance[2](
            reference.as_bytes(), query.as_bytes(), penalties, Int.MAX, EndsFree(), band
        ) if two else wavefront_distance[1](
            reference.as_bytes(), query.as_bytes(), penalties, Int.MAX, EndsFree(), band
        )
        if cost < 0:
            raise outside(band)
        return penalties.score(cost, columns + rows)
    if mode.match_score <= 0:
        raise AlignmentError(
            ErrorKind.INVALID_SCORING, "a local alignment under Costs needs a match that earns: Mode.local"
        )
    if not band.covers(columns, rows):
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a sweep takes no band: a local alignment or free ends")
    if mode.kind == Mode.SMITH_WATERMAN:
        return end_of(reference.as_bytes(), query.as_bytes(), costs, mode.match_score)[0]
    return swept[FROM_EDGE](
        reference.as_bytes(), query.as_bytes(), costs, mode.match_score, free_ends(mode, columns, rows)
    )[0]


def local_scores(
    reference: String, query: String, costs: Costs, mode: Mode, *, window: Optional[Int] = None
) raises AlignmentError -> LocalScores:
    """A local alignment's best score and where it ends, with the best score of an alignment ending
    more than `window` reference letters away, as SSW's `score2` and `ref_end2` report it for a mapping
    quality (see `LocalScores`): one sweep, with no alignment traced. The window is half the query, and
    at least 15, by default, as SSW suggests."""
    if mode.kind != Mode.SMITH_WATERMAN or mode.match_score <= 0:
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "local scores for a local alignment: Mode.local")
    _ = penalties_of(costs)
    var span = window.or_else(max(query.byte_length() // 2, 15))
    return local_scores_of(reference.as_bytes(), query.as_bytes(), costs, mode.match_score, span)


def cheapest_extension(costs: Costs) -> Int:
    """The cheapest a gap grows a letter, either way, at either piece: a Z-drop's slack a diagonal, as
    KSW2 charges a long gap."""
    var cheapest = min(costs.extension, costs.deletion_extension)
    if costs.pieces() == 2:
        cheapest = min(cheapest, min(costs.extension2, costs.deletion_extension2))
    return cheapest


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
    var space = SearchSpace()
    return cost_within(reference, query, costs, mode, band, max_cost, space)


def cost_within(
    reference: String, query: String, costs: Costs, mode: Mode, band: Band, max_cost: Int, mut space: SearchSpace
) raises AlignmentError -> Optional[Int]:
    """`cost_within` through `space`'s searches, which a batch's worker keeps from pair to pair (see
    `gap_affine.SearchSpace`)."""
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
                return edit_distance(reference, query, space.edit) * scale
            return edit_search(query, reference, ends.first_begin == 0).distance * scale
        except error:
            # More symbols than the sweep takes: the wavefront takes any.
            if error.kind != ErrorKind.UNKNOWN_SYMBOL:
                raise error
    var penalties = space.penalties_for(costs)
    if max_cost < 0:
        return None
    # No cap is the usual case, and a 64-bit division per pair counts when short reads take a microsecond.
    var ceiling = Int.MAX if max_cost == Int.MAX else max_cost // penalties.scale
    var cost: Int
    if costs.pieces() == 2:
        cost = wavefront_distance[2](reference.as_bytes(), query.as_bytes(), penalties, ceiling, space, ends, band)
    else:
        cost = wavefront_distance[1](reference.as_bytes(), query.as_bytes(), penalties, ceiling, space, ends, band)
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
    eqx: Bool,
    limit: Int,
) raises AlignmentError -> Optional[Alignment]:
    """An optimal alignment, or None when its cost would pass `max_cost` (`Int.MAX` for no cap) or none
    fits `band`: the least costly one, or for a mode that maximizes a score the best-scoring one, its
    kept fronts within `limit` diagonals."""
    var space = SearchSpace()
    return aligned_within(reference, query, costs, mode, band, max_cost, ties, eqx, limit, space)


def aligned_within(
    reference: String,
    query: String,
    costs: Costs,
    mode: Mode,
    band: Band,
    max_cost: Int,
    ties: Ties,
    eqx: Bool,
    limit: Int,
    mut space: SearchSpace,
) raises AlignmentError -> Optional[Alignment]:
    """`aligned_within` through `space`'s searches, which a batch's worker keeps from pair to pair."""
    if not mode.is_scored():
        return least_costly(reference, query, costs, mode, band, max_cost, ties, eqx, limit, space)
    if max_cost != Int.MAX:
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a mode with a match score takes no cost cap")
    _ = penalties_of(costs)
    return best_scoring(reference, query, costs, mode, band, ties, eqx, limit)


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


def bits_serve(costs: Costs, mode: Mode, band: Band, max_cost: Int) -> Bool:
    """Whether the bit-parallel sweep serves a batch's pairs (see `swept_by_bits`), as it does every pair
    of a batch of unit costs with no cap and no band, globally or with the query found inside or at the
    start of the reference. Such a batch's alignments take the sweep pair by pair, which beat the lanes on
    short reads and long: 1.26 against 1.44 us a 150 bp read, 21.9 against 23.5 a 1 kbp read at 10%, on
    the Skylake-X. Its distances take the lanes, which beat it on short reads, 0.23 against 0.40 us."""
    if costs.unit_scale() == 0 or max_cost != Int.MAX or not band.covers(1 << 40, 1 << 40):
        return False
    if mode.is_scored() or mode.query_start != 0 or mode.query_end != 0:
        return False
    var whole = 1 << 40
    if mode.reference_start == 0:
        return mode.reference_end == 0 or mode.reference_end >= whole
    return mode.reference_start >= whole and mode.reference_end >= whole


def least_costly(
    reference: String,
    query: String,
    costs: Costs,
    mode: Mode,
    band: Band,
    max_cost: Int,
    ties: Ties,
    eqx: Bool,
    limit: Int,
    mut space: SearchSpace,
) raises AlignmentError -> Optional[Alignment]:
    """The least costly alignment with `mode`'s free ends: by the bit-parallel sweep where it serves the
    pair (see `swept_by_bits`), else by the wavefront, the same alignment at the same costs, through
    `space`'s searches."""
    var columns = reference.byte_length()
    var rows = query.byte_length()
    var ends = free_ends(mode, columns, rows)
    if swept_by_bits(costs, ends, band, max_cost, columns, rows):
        try:
            var scale = costs.unit_scale()
            if ends.first_begin == 0 and ends.first_end == 0:
                var whole = edit_cigar(reference, query, eqx, ties)
                var cost = whole.distance * scale
                return Alignment(cost, -cost, whole.cigar, 0, columns, 0, rows)
            var hit = edit_search(query, reference, ends.first_begin == 0, ties)
            var part = String(StringSlice(unsafe_from_utf8=reference.as_bytes()[hit.start : hit.end]))
            var found = edit_cigar(part, query, eqx, ties)
            var cost = found.distance * scale
            return Alignment(cost, -cost, found.cigar, hit.start, hit.end, 0, rows)
        except error:
            if error.kind != ErrorKind.UNKNOWN_SYMBOL:
                raise error
    var penalties = space.penalties_for(costs)
    if max_cost < 0:
        return None
    var ceiling = Int.MAX if max_cost == Int.MAX else max_cost // penalties.scale
    if mode.is_global():
        var found: Optional[AffineCigar]
        if costs.pieces() == 2:
            found = cigar_within[2](reference, query, penalties, eqx, ceiling, space, band, ties, limit)
        else:
            found = cigar_within[1](reference, query, penalties, eqx, ceiling, space, band, ties, limit)
        if not found:
            return None
        var cost = found.value().cost
        return Alignment(cost, -cost, found.take().cigar, 0, columns, 0, rows)
    var spanned: Optional[Spanned]
    if costs.pieces() == 2:
        spanned = free_ends_alignment[2](
            space.forward2, space.backward2, reference, query, penalties, eqx, ceiling, ends, band, ties, limit
        )
    else:
        spanned = free_ends_alignment[1](
            space.forward, space.backward, reference, query, penalties, eqx, ceiling, ends, band, ties, limit
        )
    if not spanned:
        return None
    ref found = spanned.value()
    return Alignment(
        found.cost,
        -found.cost,
        found.cigar,
        found.first_start,
        found.first_end,
        found.second_start,
        found.second_end,
    )


def best_scoring(
    reference: String, query: String, costs: Costs, mode: Mode, band: Band, ties: Ties, eqx: Bool, limit: Int
) raises AlignmentError -> Alignment:
    """The best-scoring alignment for a mode with a match score: an extension by the wavefront from its
    anchor; a global alignment, its letters fixed so the reward folds into the costs, by the wavefront
    too; a local alignment or other free ends by sweep (see `scored`)."""
    var columns = reference.byte_length()
    var rows = query.byte_length()
    if mode.kind == Mode.EXTENSION:
        return extended_alignment(reference, query, costs, mode, band, ties, eqx, limit)
    if mode.kind == Mode.ENDS and mode.is_global():
        var whole = global_rewarded(reference, query, costs, mode.match_score, band, ties, eqx, limit)
        var score = mode.match_score * cigar_matches(reference, query, whole[1]) - whole[0]
        return Alignment(whole[0], score, whole[1], 0, columns, 0, rows)
    if mode.match_score <= 0:
        raise AlignmentError(
            ErrorKind.INVALID_SCORING, "a local alignment under Costs needs a match that earns: Mode.local"
        )
    if not band.covers(columns, rows):
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a sweep takes no band: a local alignment or free ends")
    if mode.kind == Mode.SMITH_WATERMAN:
        return local_alignment(reference, query, costs, mode.match_score, ties, eqx, limit)
    return rewarded_alignment(
        reference, query, costs, mode.match_score, free_ends(mode, columns, rows), ties, eqx, limit
    )


def extended_alignment(
    reference: String, query: String, costs: Costs, mode: Mode, band: Band, ties: Ties, eqx: Bool, limit: Int
) raises AlignmentError -> Alignment:
    """The best extension from `mode`'s anchor (see `Mode.extension`)."""
    if mode.end_bonus > 0 and not band.covers(reference.byte_length(), query.byte_length()):
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "an end bonus takes no band")
    var two = costs.pieces() == 2
    var penalties = rewarded_penalties(mode.match_score, costs)
    # The Z-drop's slack a diagonal is the cheapest extension, as KSW2 charges a long gap.
    var drop_extension = cheapest_extension(costs)
    var found = extension_of[2](
        reference, query, penalties, eqx, mode.anchor, band, ties, -1, limit, mode.zdrop, drop_extension, mode.end_bonus
    ) if two else extension_of[1](
        reference, query, penalties, eqx, mode.anchor, band, ties, -1, limit, mode.zdrop, drop_extension, mode.end_bonus
    )
    # The search weighs the end bonus itself (see `extension_of`), but for an empty side or no reward,
    # where it does not run: there the free ends' own alignment through the query is weighed against it.
    var searched = reference.byte_length() > 0 and query.byte_length() > 0 and penalties.reward > 0
    if mode.end_bonus > 0 and not found.dropped and not searched:
        var reaching = aligned_within(reference, query, costs, mode.reaching_end(), band, Int.MAX, ties, eqx, limit)
        if reaching and reaching.value().score + mode.end_bonus > found.score:
            return reaching.take()
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


# endregion Costs

# region Scoring


def score(
    reference: String,
    query: String,
    scoring: Scoring,
    mode: Mode = Mode.GLOBAL,
    *,
    placement: Optional[Placement] = None,
) raises -> Int:
    """The optimal score under `scoring`, with no alignment traced: `Mode.GLOBAL` and `Mode.LOCAL` in two
    rows of memory on either device (see `scoring.score_with`), free ends and extensions by sweep on the
    host. The table holds what a match earns, so a mode's own match score must be zero:
    `Mode.extension(0)` for an extension."""
    return scoring_score(reference, query, scoring, mode, placement)


def align(
    reference: String,
    query: String,
    scoring: Scoring,
    mode: Mode = Mode.GLOBAL,
    *,
    placement: Optional[Placement] = None,
    max_memory: Int = DEFAULT_MAX_MEMORY,
    eqx: Bool = True,
) raises -> Alignment:
    """An optimal alignment under `scoring`, as `Costs` give one (see `Alignment`), its `cost` minus its
    score: both sequences whole for `Mode.GLOBAL`, Needleman-Wunsch, the best-scoring window of each
    for `Mode.LOCAL`, Smith-Waterman, on either device; free ends and extensions, with Z-drop as KSW2
    gauges it, on the host, their span by sweep and the letters between aligned globally (see
    `scoring.scoring_alignment`). Its rows come back with `Alignment.gapped`. Of equally good
    alignments, Gotoh's walk picks the CIGAR (see `alignment.reconstruct`), not `Ties`. A traceback whose
    stored matrix would pass `max_memory` bytes recurses in linear space instead (see
    `scoring.cells_within`)."""
    return scoring_alignment(reference, query, scoring, mode, placement, cells_within(max_memory), eqx)


def scores(
    references: List[String],
    queries: List[String],
    scoring: Scoring,
    mode: Mode = Mode.GLOBAL,
    *,
    placement: Optional[Placement] = None,
) raises -> List[Int]:
    """`score` for every pair; on the device, every pair one block can carry goes out in one launch."""
    var found: List[Int32]
    var kind = mode.kind
    if kind == Mode.SMITH_WATERMAN:
        found = scores_with[AlignmentMode.LOCAL](references, queries, scoring, placement)
    elif mode.is_global() and mode.match_score == 0:
        found = scores_with[AlignmentMode.GLOBAL](references, queries, scoring, placement)
    else:
        var pairs = paired_length(references, queries)
        var results = List[Int](capacity=pairs)
        for index in range(pairs):
            results.append(score(references[index], queries[index], scoring, mode, placement=placement))
        return results^
    var results = List[Int](capacity=len(found))
    for value in found:
        results.append(Int(value))
    return results^


def alignments(
    references: List[String],
    queries: List[String],
    scoring: Scoring,
    mode: Mode = Mode.GLOBAL,
    *,
    placement: Optional[Placement] = None,
    max_memory: Int = DEFAULT_MAX_MEMORY,
    eqx: Bool = True,
) raises -> List[Alignment]:
    """`align` for every pair; on the device, every pair both bounds admit goes out in one launch."""
    var pairs = paired_length(references, queries)
    var stored_cells = cells_within(max_memory)
    var resolved = placement.or_else(Placement.default())
    if resolved.device != Device.GPU:
        # On the host each pair as `align` takes it, the pairs over the threads asked for.
        var results = List[Alignment](capacity=pairs)
        for _ in range(pairs):
            results.append(Alignment(0, 0, String(), 0, 0, 0, 0))
        var out = results.unsafe_ptr()
        var failed = List[Bool](length=pairs, fill=False)
        var flags = failed.unsafe_ptr()
        var single = Placement.on_cpu(1)
        var workers = max(resolved.threads, 1)
        var chunks = max(min(pairs, workers * 8), 1)
        # Global alignments under a table of one match and one mismatch score: many pairs at once in the
        # lanes (see `scoring.laned_alignments`); the rest one at a time.
        var settled = List[Bool](length=pairs, fill=False)
        if mode.is_global() and mode.match_score == 0:
            settled = laned_alignments(references, queries, scoring, eqx, workers, max_memory, out)
        var settled_ptr = settled.unsafe_ptr()

        def align_chunk(
            chunk: Int,
        ) {
            imm references,
            imm queries,
            imm scoring,
            imm mode,
            imm single,
            imm stored_cells,
            imm eqx,
            imm out,
            imm flags,
            imm pairs,
            imm chunks,
            imm settled_ptr,
        }:
            """Aligns chunk `chunk` of the pairs, flagging any pair that raised for the serial retry below."""
            var space = SearchSpace()
            for index in range(pairs * chunk // chunks, pairs * (chunk + 1) // chunks):
                if settled_ptr[unsafe_offset=index]:
                    continue
                try:
                    out[unsafe_offset=index] = scoring_alignment(
                        references[index], queries[index], scoring, mode, single, stored_cells, eqx, space
                    )
                except:
                    flags[unsafe_offset=index] = True

        spread(align_chunk, chunks, workers)
        # A pair that failed raises here, the same error a serial loop would have raised first.
        for index in range(pairs):
            if failed[index]:
                results[index] = scoring_alignment(
                    references[index], queries[index], scoring, mode, single, stored_cells, eqx
                )
        return results^
    var gapped: List[GappedAlignment]
    if mode.kind == Mode.SMITH_WATERMAN and mode.match_score == 0:
        gapped = alignments_with[AlignmentMode.LOCAL](references, queries, scoring, placement, stored_cells)
    elif mode.is_global() and mode.match_score == 0:
        gapped = alignments_with[AlignmentMode.GLOBAL](references, queries, scoring, placement, stored_cells)
    else:
        var results = List[Alignment](capacity=pairs)
        for index in range(pairs):
            results.append(
                align(
                    references[index],
                    queries[index],
                    scoring,
                    mode,
                    placement=placement,
                    max_memory=max_memory,
                    eqx=eqx,
                )
            )
        return results^
    # Each pair's rows read as its CIGAR, on the host threads the placement asks for.
    var results = List[Alignment](capacity=pairs)
    for _ in range(pairs):
        results.append(Alignment(0, 0, String(), 0, 0, 0, 0))
    var out = results.unsafe_ptr()
    var whole = mode.is_global()
    var workers = max(resolved.threads, 1)
    var chunks = max(min(pairs, workers * 8), 1)

    def convert(
        chunk: Int,
    ) {imm gapped, imm references, imm queries, imm out, imm whole, imm eqx, imm pairs, imm chunks}:
        """Reads chunk `chunk` of the gapped pairs as `Alignment`s, each with its CIGAR and spans."""
        for index in range(pairs * chunk // chunks, pairs * (chunk + 1) // chunks):
            out[unsafe_offset=index] = as_alignment(gapped[index], references[index], queries[index], whole, eqx)

    spread(convert, chunks, workers)
    return results^


# endregion Scoring

# region Batches


def longest_first(references: List[String], queries: List[String], workers: Int) -> List[Int]:
    """The pairs' indices, the longest pair first: taken in that order by whichever thread is free, the
    long pairs start first and the short ones fill in around them, so none is left alone at the end
    holding up the rest.

    Balancing a batch needs no finer order than lengths within an eighth of each other, so the pairs
    are dealt into such classes, no sort, within a class in the batch's own order. Each of `workers`
    counts its stretch of the batch's classes and then places its stretch's pairs: on one thread the
    order of short reads took as long as a fifteenth of the work ten threads did on them."""
    comptime CLASSES = 8 * 62
    var pairs = len(references)
    # A stretch holds at least as many pairs as there are classes, so summing its counts costs no more
    # than counting them did.
    var stretches = max(min(workers, pairs // CLASSES), 1)
    var classes = List[Int](capacity=pairs)
    classes.resize(unsafe_uninit_length=pairs)
    var order = List[Int](capacity=pairs)
    order.resize(unsafe_uninit_length=pairs)
    var starts = List[Int](length=stretches * CLASSES, fill=0)
    var class_ptr = classes.unsafe_ptr()
    var order_ptr = order.unsafe_ptr()
    var start_ptr = starts.unsafe_ptr()

    def count(stretch: Int) {imm references, imm queries, imm pairs, imm stretches, imm class_ptr, imm start_ptr}:
        """Each pair of stretch `stretch` its class, and the stretch's count of each."""
        var counts = start_ptr.unsafe_offset(stretch * CLASSES)
        for index in range(pairs * stretch // stretches, pairs * (stretch + 1) // stretches):
            var length = references[index].byte_length() + queries[index].byte_length()
            # Lengths below 16 are classes of their own; past that, a power of two and its top three
            # bits under the leading one, longest the lowest class.
            var shift = max(60 - Int(count_leading_zeros(UInt64(length))), 0)
            var found = CLASSES - 1 - (8 * shift + (length >> shift))
            class_ptr[unsafe_offset=index] = found
            counts[unsafe_offset=found] += 1

    def place(stretch: Int) {imm pairs, imm stretches, imm class_ptr, imm order_ptr, imm start_ptr}:
        """Each pair of stretch `stretch` at the next place its class has for the stretch."""
        var next = start_ptr.unsafe_offset(stretch * CLASSES)
        for index in range(pairs * stretch // stretches, pairs * (stretch + 1) // stretches):
            var found = class_ptr[unsafe_offset=index]
            order_ptr[unsafe_offset=next[unsafe_offset=found]] = index
            next[unsafe_offset=found] += 1

    spread(count, stretches, stretches)
    # Each stretch's first place in each class: the classes in turn, each stretch's pairs in turn.
    var placed = 0
    for found in range(CLASSES):
        for stretch in range(stretches):
            var counted = starts[stretch * CLASSES + found]
            starts[stretch * CLASSES + found] = placed
            placed += counted
    spread(place, stretches, stretches)
    return order^


def distances(
    references: List[String],
    queries: List[String],
    costs: Costs = Costs.edit(),
    mode: Mode = Mode.GLOBAL,
    *,
    band: Band = Band(),
    threads: Optional[Int] = None,
    placement: Optional[Placement] = None,
) raises -> List[Int]:
    """Every pair's `distance`, on the caller's own thread by default, or spread over `threads` threads
    when asked: an application that calls the library from threads of its own spreads its work itself,
    and it alone knows how many its machine can spare.

    The pairs are independent, so each runs on one thread start to finish, each thread taking the
    next pair of the batch, longest first, as soon as it is free (see `longest_first`). A pair that
    fails raises, after the rest, the same error a serial loop would have raised first.

    On the GPU, `placement`, unit costs align globally, a thread a pair by Myers' bit-vectors (see
    `device_edit`), every pair whose shorter sequence fits a thread's 4,096 letters; the rest, and
    every other cost or mode, on the host.
    """
    if placement and placement.value().device == Device.GPU:
        if costs.unit_scale() == 0 or not mode.is_global() or not band.covers(1 << 40, 1 << 40):
            raise AlignmentError(
                ErrorKind.INVALID_ARGUMENT, "on the GPU, distances at unit costs, globally, with no band"
            )
        return gpu_distances(references, queries, costs.unit_scale(), placement.value())
    var found = capped_distances(references, queries, costs, mode, band, Int.MAX, threads)
    var results = List[Int](capacity=len(found))
    for index in range(len(found)):
        results.append(found[index].value())
    return results^


def distances(
    references: List[String],
    queries: List[String],
    costs: Costs = Costs.edit(),
    mode: Mode = Mode.GLOBAL,
    *,
    max_cost: Int,
    band: Band = Band(),
    threads: Optional[Int] = None,
) raises AlignmentError -> List[Optional[Int]]:
    """Every pair's `distance` under `max_cost`, None for a pair past it or with no alignment inside
    `band`, the pairs spread over threads as the uncapped `distances` spreads them: a batch of
    candidates filtered by cost, the far ones costing a fraction of their full search."""
    return capped_distances(references, queries, costs, mode, band, max_cost, threads)


def gpu_distances(
    references: List[String], queries: List[String], scale: Int, placement: Placement
) raises -> List[Int]:
    """Every pair's unit-cost distance times `scale`, on the device where its shorter sequence fits a
    thread (see `device_edit`), on the host otherwise."""
    var pairs = paired_length(references, queries)
    var results = List[Int](length=pairs, fill=0)
    var patterns = List[String]()
    var texts = List[String]()
    var on_device = List[Int]()
    var on_host = List[Int]()
    comptime LIMIT = 64 * MAX_PATTERN_WORDS
    for index in range(pairs):
        var reference = references[index].byte_length()
        var query = queries[index].byte_length()
        if min(reference, query) == 0:
            results[index] = max(reference, query) * scale
        elif min(reference, query) <= LIMIT:
            var shorter_query = query <= reference
            patterns.append(queries[index] if shorter_query else references[index])
            texts.append(references[index] if shorter_query else queries[index])
            on_device.append(index)
        else:
            on_host.append(index)
    if len(on_device) > 0:
        var found = device_edit_distances(DeviceScope(placement.gpu_id), patterns, texts, placement.threads)
        for slot in range(len(on_device)):
            results[on_device[slot]] = found[slot] * scale
    for index in on_host:
        results[index] = distance(references[index], queries[index], Costs.edit()) * scale
    return results^


def capped_distances(
    references: List[String],
    queries: List[String],
    costs: Costs,
    mode: Mode,
    band: Band,
    max_cost: Int,
    threads: Optional[Int],
) raises AlignmentError -> List[Optional[Int]]:
    """Every pair's `cost_within`, on every thread asked for, longest first. With no cap, a pair no
    alignment inside `band` fits raises as a failed pair does, in the batch's order."""
    var pairs = paired_length(references, queries)
    var results = List[Optional[Int]](length=pairs, fill=None)
    if pairs == 0:
        return results^
    var workers = max(threads.or_else(1), 1)
    var out = results.unsafe_ptr()
    var failed = List[Bool](length=pairs, fill=False)
    var flags = failed.unsafe_ptr()
    # Costs with no reward, globally or with free ends: as many pairs at once as a register holds lanes
    # (see `lanes`); the pairs too long for its 16 bits, a pair with an empty side and free ends, and
    # every pair of a mode with a reward, one at a time.
    var settled = List[Bool](length=pairs, fill=False)
    var settled_ptr = settled.unsafe_ptr()
    var lane_costs = LaneCosts.of(costs, mode, free_ends=True)
    var reference_texts = StringTexts(references.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]())
    var query_texts = StringTexts(queries.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]())
    if lane_costs:
        # Costs the searches would refuse raise here, as a pair's search would raise them.
        _ = penalties_of(costs)
    if (
        lane_costs
        and lane_distances(
            pairs,
            reference_texts,
            query_texts,
            lane_costs.value(),
            band,
            max_cost,
            workers,
            out,
            settled_ptr,
            mode,
        )
        == pairs
    ):
        for index in range(pairs):
            if not results[index] and max_cost == Int.MAX:
                raise outside(band)
        return results^
    var order = longest_first(references, queries, workers)
    var taken = Atomic[Int64](0)

    def distance_worker(
        worker: Int,
    ) {
        mut taken,
        imm order,
        imm references,
        imm queries,
        imm out,
        imm flags,
        imm settled_ptr,
        imm pairs,
        imm costs,
        imm mode,
        imm band,
        imm max_cost,
        imm workers,
    }:
        """Takes the next pairs in `order` until none is left, storing each one's capped cost or flagging that
        it raised, its searches kept from pair to pair; a pair the lanes settled it leaves."""
        var space = SearchSpace()
        var last = 0
        while True:
            var share = next_share(taken, pairs, workers, last)
            if share[0] >= pairs:
                return
            for dealt in range(share[0], share[1]):
                var index = order[dealt]
                if settled_ptr[unsafe_offset=index]:
                    continue
                try:
                    out[unsafe_offset=index] = cost_within(
                        references[index], queries[index], costs, mode, band, max_cost, space
                    )
                except:
                    flags[unsafe_offset=index] = True

    spread(distance_worker, min(workers, pairs), min(workers, pairs))
    for index in range(pairs):
        if failed[index]:
            results[index] = cost_within(references[index], queries[index], costs, mode, band, max_cost)
        if not results[index] and max_cost == Int.MAX:
            raise outside(band)
    return results^


def alignments(
    references: List[String],
    queries: List[String],
    costs: Costs = Costs.edit(),
    mode: Mode = Mode.GLOBAL,
    *,
    band: Band = Band(),
    ties: Ties = Ties.LEFT,
    eqx: Bool = True,
    threads: Optional[Int] = None,
    max_memory: Int = DEFAULT_MAX_MEMORY,
) raises AlignmentError -> List[Alignment]:
    """Every pair's `align`, the pairs spread over threads as `distances` spreads them, each thread's
    kept fronts within `max_memory` bytes."""
    var found = capped_alignments(
        references, queries, costs, mode, band, Int.MAX, ties, eqx, threads, max_memory // KEPT_BYTES
    )
    var results = List[Alignment](capacity=len(found))
    for index in range(len(found)):
        results.append(found[index].take())
    return results^


def alignments(
    references: List[String],
    queries: List[String],
    costs: Costs = Costs.edit(),
    mode: Mode = Mode.GLOBAL,
    *,
    max_cost: Int,
    band: Band = Band(),
    ties: Ties = Ties.LEFT,
    eqx: Bool = True,
    threads: Optional[Int] = None,
    max_memory: Int = DEFAULT_MAX_MEMORY,
) raises AlignmentError -> List[Optional[Alignment]]:
    """Every pair's `align` under `max_cost`, None for a pair past it or with no alignment inside
    `band`, the pairs spread over threads as `distances` spreads them."""
    return capped_alignments(
        references, queries, costs, mode, band, max_cost, ties, eqx, threads, max_memory // KEPT_BYTES
    )


def capped_alignments(
    references: List[String],
    queries: List[String],
    costs: Costs,
    mode: Mode,
    band: Band,
    max_cost: Int,
    ties: Ties,
    eqx: Bool,
    threads: Optional[Int],
    limit: Int,
) raises AlignmentError -> List[Optional[Alignment]]:
    """Every pair's `aligned_within`, on every thread asked for, longest first. With no cap, a pair no
    alignment inside `band` fits raises as a failed pair does, in the batch's order."""
    var pairs = paired_length(references, queries)
    var results = List[Optional[Alignment]](capacity=pairs)
    for _ in range(pairs):
        results.append(None)
    if pairs == 0:
        return results^
    var workers = max(threads.or_else(1), 1)
    var out = results.unsafe_ptr()
    var failed = List[Bool](length=pairs, fill=False)
    var flags = failed.unsafe_ptr()
    # Costs alone: as many pairs at once as a register holds lanes, each traced from the flags its band
    # kept (see `lanes`), with free ends the span the tie rule picks first; every other pair, and every pair
    # of a mode that scores, one at a time.
    var settled = List[Bool](length=pairs, fill=False)
    var settled_ptr = settled.unsafe_ptr()
    var laned = List[Optional[Int]](length=pairs, fill=None)
    var laned_ptr = laned.unsafe_ptr()
    var paths = List[List[UInt8]](capacity=pairs)
    for _ in range(pairs):
        paths.append(List[UInt8]())
    var path_ptr = paths.unsafe_ptr()
    # Each laned pair's span: the reference's first letter and the one past its last, then the query's.
    var spans = List[Int](length=4 * pairs, fill=0)
    var span_ptr = spans.unsafe_ptr()
    # Unit costs the bit-parallel sweep serves take it pair by pair, faster than the lanes on short reads and
    # long alike: 1.26 against 1.44 us a 150 bp read, 21.9 against 23.5 a 1 kbp read at 10%, on the Skylake-X.
    # Free ends take the lanes only with no band, which would bound the span's search too.
    var unbanded = band.covers(1 << 40, 1 << 40)
    var lane_costs = LaneCosts.of(costs, mode, unbanded) if not bits_serve(costs, mode, band, max_cost) else None
    var penalties: Optional[Penalties] = None
    if lane_costs and not mode.is_global():
        penalties = penalties_of(costs)
        lane_free_alignments(
            references,
            queries,
            lane_costs.value(),
            mode,
            max_cost,
            ties == Ties.RIGHT,
            workers,
            limit * KEPT_BYTES,
            laned_ptr,
            path_ptr,
            span_ptr,
            settled_ptr,
        )
    elif lane_costs:
        # Costs the searches would refuse raise here, as a pair's search would raise them.
        penalties = penalties_of(costs)
        for index in range(pairs):
            span_ptr[unsafe_offset=4 * index + 1] = references[index].byte_length()
            span_ptr[unsafe_offset=4 * index + 3] = queries[index].byte_length()
        _ = lane_alignments(
            pairs,
            StringTexts(references.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]()),
            StringTexts(queries.unsafe_ptr().unsafe_origin_cast[ImmUntrackedOrigin]()),
            lane_costs.value(),
            band,
            max_cost,
            ties == Ties.LEFT,
            workers,
            limit * KEPT_BYTES,
            laned_ptr,
            path_ptr,
            settled_ptr,
        )
    var order = longest_first(references, queries, workers)
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
        imm max_cost,
        imm ties,
        imm eqx,
        imm limit,
        imm workers,
        imm settled_ptr,
        imm laned_ptr,
        imm path_ptr,
        imm span_ptr,
        imm penalties,
    }:
        """Takes the next pairs in `order` until none is left, storing each one's capped alignment or flagging
        that it raised; a pair the lanes settled has its CIGAR spelled from their path."""
        var last = 0
        var space = SearchSpace()
        while True:
            var share = next_share(taken, pairs, workers, last)
            if share[0] >= pairs:
                return
            for dealt in range(share[0], share[1]):
                var index = order[dealt]
                if settled_ptr[unsafe_offset=index] and not laned_ptr[unsafe_offset=index]:
                    continue
                if settled_ptr[unsafe_offset=index]:
                    var cost = laned_ptr[unsafe_offset=index].value()
                    var first_start = span_ptr[unsafe_offset=4 * index]
                    var first_end = span_ptr[unsafe_offset=4 * index + 1]
                    var second_start = span_ptr[unsafe_offset=4 * index + 2]
                    var second_end = span_ptr[unsafe_offset=4 * index + 3]
                    ref scaled = penalties.value()
                    # The lanes follow the tie rule across the whole span; the searches would too unless they
                    # split it, which they never do while both sides' fronts, each at most a diagonal of the
                    # matrix at every cost, stay within `limit`. A pair they might split takes them.
                    var letters = first_end - first_start + second_end - second_start
                    if 2 * (cost // scaled.scale + 1) * (letters + 1) <= limit:
                        var moves = List[UInt8]()
                        swap(moves, path_ptr[unsafe_offset=index])
                        var whole = (
                            first_start == 0
                            and second_start == 0
                            and first_end == references[index].byte_length()
                            and second_end == queries[index].byte_length()
                        )
                        var cigar: String
                        if whole:
                            cigar = cigar_of(
                                references[index], queries[index], moves^, cost // scaled.scale, scaled, eqx
                            )
                        else:
                            cigar = cigar_of(
                                String(
                                    StringSlice(unsafe_from_utf8=references[index].as_bytes()[first_start:first_end])
                                ),
                                String(
                                    StringSlice(unsafe_from_utf8=queries[index].as_bytes()[second_start:second_end])
                                ),
                                moves^,
                                cost // scaled.scale,
                                scaled,
                                eqx,
                            )
                        out[unsafe_offset=index] = Alignment(
                            cost, -cost, cigar^, first_start, first_end, second_start, second_end
                        )
                        continue
                try:
                    out[unsafe_offset=index] = aligned_within(
                        references[index], queries[index], costs, mode, band, max_cost, ties, eqx, limit, space
                    )
                except:
                    flags[unsafe_offset=index] = True

    spread(alignment_worker, min(workers, pairs), min(workers, pairs))
    for index in range(pairs):
        if failed[index]:
            results[index] = aligned_within(
                references[index], queries[index], costs, mode, band, max_cost, ties, eqx, limit
            )
        if not results[index] and max_cost == Int.MAX:
            raise outside(band)
    return results^


def search(
    references: List[String],
    query: String,
    costs: Costs = Costs.edit(),
    mode: Mode = Mode.GLOBAL,
    *,
    best: Optional[Int] = None,
    max_cost: Optional[Int] = None,
    aligned: Bool = False,
    ties: Ties = Ties.LEFT,
    threads: Optional[Int] = None,
) raises AlignmentError -> List[Hit]:
    """The query against every reference, a database search: each reference's `Hit`, its score, the
    best first, ties by the references' order; with `best` that many alone, with `max_cost` (a mode
    with no reward) those within it alone, and with `aligned` each kept hit's alignment too.

    A local alignment scores a block of references at once, one to a SIMD lane, as SWIPE does (see
    `local_scores_by_block`); a mode with no reward takes `distances`, under the cap when there is one;
    any other mode each pair's `score`. Every kept hit is then aligned on its own, when asked for, by
    `align`."""
    var count = len(references)
    var workers = max(threads.or_else(1), 1)
    var scores = List[Int](length=count, fill=Int.MIN)
    if mode.kind == Mode.SMITH_WATERMAN and mode.match_score > 0:
        _ = penalties_of(costs)
        scores = local_scores_by_block(query, references, costs, mode.match_score, workers)
    elif not mode.is_scored():
        var queries = List[String](length=count, fill=query)
        var found = capped_distances(
            references, queries, costs, mode, Band(), max_cost.or_else(Int.MAX), Optional[Int](workers)
        )
        for index in range(count):
            if found[index]:
                scores[index] = -found[index].value()
    else:
        if max_cost:
            raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a mode with a match score takes no cost cap")
        var out = scores.unsafe_ptr()
        var failed = List[Bool](length=count, fill=False)
        var flags = failed.unsafe_ptr()
        var taken = Atomic[Int64](0)

        def work(slot: Int) {mut taken, imm}:
            """Scores the next unclaimed references against the query until none is left, flagging any that
            raised."""
            var last = 0
            while True:
                var share = next_share(taken, count, workers, last)
                if share[0] >= count:
                    return
                for index in range(share[0], share[1]):
                    try:
                        out[unsafe_offset=index] = score(references[index], query, costs, mode)
                    except:
                        flags[unsafe_offset=index] = True

        spread(work, min(workers, max(count, 1)), min(workers, max(count, 1)))
        for index in range(count):
            if failed[index]:
                scores[index] = score(references[index], query, costs, mode)
    # The best first, ties by the references' order; a reference past the cap is no hit.
    var order = List[Int](capacity=count)
    for index in range(count):
        if scores[index] != Int.MIN:
            order.append(index)

    def ahead(left: Int, right: Int) {imm scores} -> Bool:
        """Whether reference `left` ranks before `right`: the higher score, then the earlier reference."""
        if scores[left] != scores[right]:
            return scores[left] > scores[right]
        return left < right

    sort(order, ahead)
    var kept = min(len(order), best.or_else(len(order)))
    var hits = List[Hit](capacity=kept)
    for rank in range(kept):
        var index = order[rank]
        var alignment = Optional[Alignment]()
        if aligned:
            alignment = align(references[index], query, costs, mode, ties=ties)
        hits.append(Hit(index, scores[index], alignment^))
    return hits^


# endregion Batches
