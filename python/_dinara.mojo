"""
The Python extension under the `dinara_align` package (see `dinara_align/__init__.py`): each function
takes the sequences as `str` and the costs, the mode and the options as tuples of integers, in the C
API's fields (see `c/dinara.h`), and returns plain Python values. The package builds those tuples
and wraps what comes back; this module is not meant to be called directly.

    pixi run build-python   # build/python: the package, this module and the runtime libraries it loads
"""

from std.os import abort
from std.python import Python, PythonObject
from std.python.bindings import PythonModuleBuilder

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
    Scoring,
    Ties,
    align,
    alignments,
    distance,
    distances,
    local_scores,
    score,
)

comptime C_ENDS_FREE = 0
comptime C_EXTENSION = 1
comptime C_LOCAL = 2
comptime C_OVERLAP = 3


@export
def PyInit__dinara() abi("C") -> PythonObject:
    try:
        var module = PythonModuleBuilder("_dinara")
        module.def_function[py_distance]("distance")
        module.def_function[py_align]("align")
        module.def_function[py_score]("score")
        module.def_function[py_local_scores]("local_scores")
        module.def_function[py_distances]("distances")
        module.def_function[py_alignments]("alignments")
        module.def_function[py_scoring_align]("scoring_align")
        module.def_function[py_scoring_score]("scoring_score")
        return module.finalize()
    except error:
        abort(String("dinara_align: the extension failed to load: ", error))


def ints(fields: PythonObject) raises -> List[Int]:
    """A tuple of Python integers."""
    var out = List[Int]()
    for field in fields:
        out.append(Int(py=field))
    return out^


def costs_of(fields: PythonObject) raises -> Costs:
    """`(mismatch, opening, extension, opening2, extension2, deletion_opening, deletion_extension,
    deletion_opening2, deletion_extension2)`, as `dinara_costs`."""
    var at = ints(fields)
    var costs: Costs
    if at[3] < 0:
        costs = Costs.affine(at[0], at[1], at[2])
    else:
        costs = Costs.two_piece(at[0], at[1], at[2], at[3], at[4])
    if at[6] > 0:
        costs = costs.with_deletions(at[5], at[6], at[7], at[8])
    return costs


def mode_of(fields: PythonObject) raises -> Mode:
    """`(kind, reference_start, reference_end, query_start, query_end, match_score, anchor, zdrop)`,
    as `dinara_mode`; for a `Scoring`, kind `C_LOCAL` with a match score of zero is `Mode.LOCAL`."""
    var at = ints(fields)
    if at[0] == C_EXTENSION:
        var anchor = Anchor.END if at[6] != 0 else Anchor.START
        if at[7] >= 0:
            return Mode.extension(at[5], anchor, zdrop=at[7])
        return Mode.extension(at[5], anchor)
    if at[0] == C_LOCAL:
        return Mode.LOCAL if at[5] == 0 else Mode.local(at[5])
    if at[0] == C_OVERLAP:
        return Mode.overlap(at[5])
    if at[0] != C_ENDS_FREE:
        raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "an unknown mode")
    return Mode.ends_free(
        reference_start=at[1], reference_end=at[2], query_start=at[3], query_end=at[4], match_score=at[5]
    )


@fieldwise_init
struct Options(ImplicitlyCopyable):
    """`(band_low, band_high, max_cost, extended, right_ties, max_memory)`: a negative cap for none, a
    memory of zero or less for the default."""

    var band: Band
    var max_cost: Int
    var extended: Bool
    var ties: Ties
    var max_memory: Int


def options_of(fields: PythonObject) raises -> Options:
    var at = ints(fields)
    comptime EDGE = 1 << 60
    var band = Band(max(at[0], -EDGE), min(at[1], EDGE))
    var memory = at[5] if at[5] > 0 else DEFAULT_MAX_MEMORY
    return Options(band, at[2], at[3] != 0, Ties.RIGHT if at[4] != 0 else Ties.LEFT, memory)


def alignment_tuple(found: Alignment) raises -> PythonObject:
    """`(cost, score, cigar, reference_start, reference_end, query_start, query_end)`."""
    return Python.tuple(
        PythonObject(found.cost),
        PythonObject(found.score),
        PythonObject(found.cigar),
        PythonObject(found.reference_start),
        PythonObject(found.reference_end),
        PythonObject(found.query_start),
        PythonObject(found.query_end),
    )


def py_distance(
    reference: PythonObject, query: PythonObject, costs: PythonObject, mode: PythonObject, options: PythonObject
) raises -> PythonObject:
    """The least cost, or None past a cap of zero or more."""
    var asked = options_of(options)
    var first = String(py=reference)
    var second = String(py=query)
    if asked.max_cost < 0:
        return PythonObject(distance(first, second, costs_of(costs), mode_of(mode), band=asked.band))
    var found = distance(first, second, costs_of(costs), mode_of(mode), max_cost=asked.max_cost, band=asked.band)
    if not found:
        return Python.none()
    return PythonObject(found.value())


def py_align(
    reference: PythonObject, query: PythonObject, costs: PythonObject, mode: PythonObject, options: PythonObject
) raises -> PythonObject:
    """An optimal alignment as a tuple (see `alignment_tuple`), or None past a cap of zero or more."""
    var asked = options_of(options)
    var first = String(py=reference)
    var second = String(py=query)
    if asked.max_cost < 0:
        return alignment_tuple(
            align(
                first,
                second,
                costs_of(costs),
                mode_of(mode),
                band=asked.band,
                ties=asked.ties,
                extended=asked.extended,
                max_memory=asked.max_memory,
            )
        )
    var found = align(
        first,
        second,
        costs_of(costs),
        mode_of(mode),
        max_cost=asked.max_cost,
        band=asked.band,
        ties=asked.ties,
        extended=asked.extended,
        max_memory=asked.max_memory,
    )
    if not found:
        return Python.none()
    return alignment_tuple(found.value())


def py_score(
    reference: PythonObject, query: PythonObject, costs: PythonObject, mode: PythonObject, options: PythonObject
) raises -> PythonObject:
    """The score `align` would return, with no alignment traced."""
    var asked = options_of(options)
    return PythonObject(score(String(py=reference), String(py=query), costs_of(costs), mode_of(mode), band=asked.band))


def py_local_scores(
    reference: PythonObject, query: PythonObject, costs: PythonObject, mode: PythonObject, window: PythonObject
) raises -> PythonObject:
    """`(score, reference_end, query_end, second_score, second_reference_end)`; a negative window for
    SSW's default."""
    var span = Int(py=window)
    var found: LocalScores
    if span >= 0:
        found = local_scores(String(py=reference), String(py=query), costs_of(costs), mode_of(mode), window=span)
    else:
        found = local_scores(String(py=reference), String(py=query), costs_of(costs), mode_of(mode))
    return Python.tuple(
        PythonObject(found.score),
        PythonObject(found.reference_end),
        PythonObject(found.query_end),
        PythonObject(found.second_score),
        PythonObject(found.second_reference_end),
    )


def strings(items: PythonObject) raises -> List[String]:
    var out = List[String]()
    for item in items:
        out.append(String(py=item))
    return out^


def threads_of(count: PythonObject) raises -> Optional[Int]:
    var threads = Int(py=count)
    return Optional[Int](threads) if threads > 0 else None


def py_distances(
    references: PythonObject,
    queries: PythonObject,
    costs: PythonObject,
    mode: PythonObject,
    options: PythonObject,
    threads: PythonObject,
) raises -> PythonObject:
    """Every pair's least cost, None past a cap, over `threads` threads, every one for zero."""
    var asked = options_of(options)
    var out = Python.list()
    if asked.max_cost < 0:
        for value in distances(
            strings(references),
            strings(queries),
            costs_of(costs),
            mode_of(mode),
            band=asked.band,
            threads=threads_of(threads),
        ):
            out.append(PythonObject(value))
        return out
    for value in distances(
        strings(references),
        strings(queries),
        costs_of(costs),
        mode_of(mode),
        max_cost=asked.max_cost,
        band=asked.band,
        threads=threads_of(threads),
    ):
        out.append(PythonObject(value.value()) if value else Python.none())
    return out


def py_alignments(
    references: PythonObject,
    queries: PythonObject,
    costs: PythonObject,
    mode: PythonObject,
    options: PythonObject,
    threads: PythonObject,
) raises -> PythonObject:
    """Every pair's alignment as a tuple, None past a cap, over `threads` threads."""
    var asked = options_of(options)
    var out = Python.list()
    if asked.max_cost < 0:
        for found in alignments(
            strings(references),
            strings(queries),
            costs_of(costs),
            mode_of(mode),
            band=asked.band,
            ties=asked.ties,
            extended=asked.extended,
            threads=threads_of(threads),
            max_memory=asked.max_memory,
        ):
            out.append(alignment_tuple(found))
        return out
    for found in alignments(
        strings(references),
        strings(queries),
        costs_of(costs),
        mode_of(mode),
        max_cost=asked.max_cost,
        band=asked.band,
        ties=asked.ties,
        extended=asked.extended,
        threads=threads_of(threads),
        max_memory=asked.max_memory,
    ):
        out.append(alignment_tuple(found.value()) if found else Python.none())
    return out


def scoring_of(table: PythonObject) raises -> Scoring:
    """`(alphabet, cells, opening, extension)`: a `Scoring` over the alphabet's letters, its table row
    by row, a gap of `k` letters scoring `opening + k extension`."""
    var cells = List[Int8]()
    for cell in table[1]:
        cells.append(Int8(Int(py=cell)))
    return Scoring.tabulated(String(py=table[0]), cells^, Int(py=table[2]), Int(py=table[3]))


def py_scoring_align(
    reference: PythonObject, query: PythonObject, table: PythonObject, mode: PythonObject, extended: PythonObject
) raises -> PythonObject:
    """An optimal alignment under a `Scoring` (see `scoring_of`), as a tuple."""
    return alignment_tuple(
        align(
            String(py=reference),
            String(py=query),
            scoring_of(table),
            mode_of(mode),
            extended=Bool(py=extended),
        )
    )


def py_scoring_score(
    reference: PythonObject, query: PythonObject, table: PythonObject, mode: PythonObject
) raises -> PythonObject:
    """The best score under a `Scoring`, with no alignment traced."""
    return PythonObject(score(String(py=reference), String(py=query), scoring_of(table), mode_of(mode)))
