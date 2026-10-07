"""
Differential fuzzing: random pairs, costs, modes, bands, caps and tie rules, every answer held to the
full matrix (see `oracle`).

    pixi run fuzz [iterations] [seed]

Each trial checks that `align`'s alignment is one the mode allows, stays on the band, prices to its own
cost and score, and earns the matrix's optimum; that `distance` agrees, under a cap too, and raises
where no alignment fits the band; that the two tie rules mirror each other over reversed sequences;
and, every so often, that a batch is its pairs one at a time. A failing trial prints what reproduces
it and the run exits nonzero.
"""

from std.random import random_float64, random_ui64, seed
from std.sys import argv

from dinara_align import (
    DEFAULT_MAX_MEMORY,
    Alignment,
    AlignmentError,
    Anchor,
    Band,
    Costs,
    Mode,
    Ties,
    align,
    alignments,
    distance,
    score,
)

from oracle import ENDS, EXTENSION, LOCAL, Model, check_alignment, optimum, reversed_bytes, rule_span


def draw(low: Int, high: Int) -> Int:
    return Int(random_ui64(UInt64(low), UInt64(high)))


def chance(probability: Float64) -> Bool:
    return random_float64() < probability


def letters(length: Int, alphabet: String) -> String:
    var symbols = alphabet.as_bytes()
    var out = List[UInt8](capacity=length)
    for _ in range(length):
        out.append(symbols[draw(0, len(symbols) - 1)])
    return String(unsafe_from_utf8=out^)


def mutated(text: String, rate: Float64, longest_gap: Int, alphabet: String) -> String:
    var bytes = text.as_bytes()
    var symbols = alphabet.as_bytes()
    var out = List[UInt8]()
    var index = 0
    while index < len(bytes):
        var roll = random_float64()
        if roll < rate / 3:
            out.append(symbols[draw(0, len(symbols) - 1)])
            index += 1
        elif roll < 2 * rate / 3:
            index += draw(1, longest_gap)
        elif roll < rate:
            for _ in range(draw(1, longest_gap)):
                out.append(symbols[draw(0, len(symbols) - 1)])
        else:
            out.append(bytes[index])
            index += 1
    return String(unsafe_from_utf8=out^)


def sequences() -> Tuple[String, String]:
    """A pair: unrelated, one a mutation of the other, or a shared core inside unrelated flanks."""
    var alphabet = "ACGT" if chance(0.8) else ("ACGTN" if chance(0.5) else "ACGTRYKMSWN")
    var longest = [0, 3, 12, 60, 200, 700][draw(0, 5)]
    var reference = letters(draw(0, longest), alphabet)
    var shape = draw(0, 2)
    if shape == 0:
        return (reference, letters(draw(0, longest), alphabet))
    var rate = [0.0, 0.02, 0.1, 0.3, 0.6][draw(0, 4)]
    var query = mutated(reference, rate, draw(1, 12), alphabet)
    if shape == 2:
        query = letters(draw(0, 20), alphabet) + query + letters(draw(0, 20), alphabet)
        reference = letters(draw(0, 20), alphabet) + reference + letters(draw(0, 20), alphabet)
    if chance(0.5):
        return (query, reference)
    return (reference, query)


def random_costs() raises AlignmentError -> Costs:
    var kind = draw(0, 3)
    if kind == 0:
        return Costs.edit()
    var mismatch = draw(1, 9)
    if kind == 1:
        return Costs.linear(mismatch, draw(1, 5))
    var opening = draw(0, 12)
    var extension = draw(1, 5)
    if kind == 2:
        return Costs.affine(mismatch, opening, extension)
    return Costs.two_piece(mismatch, opening, extension, opening + draw(1, 30), draw(1, extension))


def allowance(length: Int) -> Int:
    var choice = draw(0, 3)
    if choice == 0:
        return 0
    if choice == 1:
        return 1 << 60
    return draw(0, length + 1)


def random_mode(columns: Int, rows: Int) raises AlignmentError -> Mode:
    var kind = draw(0, 10)
    var reward = draw(1, 4)
    if kind == 0:
        return Mode.GLOBAL
    if kind == 1:
        return Mode.INFIX
    if kind == 2:
        return Mode.PREFIX if chance(0.5) else Mode.SUFFIX
    if kind == 3:
        return Mode.REFERENCE_IN_QUERY
    if kind == 4 or kind == 5:
        var free = Mode.ends_free(
            reference_start=allowance(columns),
            reference_end=allowance(columns),
            query_start=allowance(rows),
            query_end=allowance(rows),
        )
        return free.with_match_score(reward) if kind == 5 else free
    if kind == 6:
        return Mode.INFIX.with_match_score(reward) if chance(0.5) else Mode.GLOBAL.with_match_score(reward)
    if kind == 7:
        var anchor = Anchor.END if chance(0.5) else Anchor.START
        if chance(0.3):
            return Mode.extension(draw(0, 4), anchor, zdrop=draw(0, 60))
        return Mode.extension(draw(0, 4), anchor)
    if kind == 8:
        return Mode.local(reward)
    return Mode.overlap(reward)


def model_of(costs: Costs, mode: Mode, band: Band) -> Model:
    var pieces = List[Tuple[Int, Int]]()
    pieces.append((costs.opening, costs.extension))
    if costs.pieces() == 2:
        pieces.append((costs.opening2, costs.extension2))
    var kind = ENDS
    if mode.kind == Mode.EXTENSION:
        kind = EXTENSION
    elif mode.kind == Mode.SMITH_WATERMAN:
        kind = LOCAL
    return Model(
        Model.uniform(mode.match_score, costs.mismatch),
        pieces.copy(),
        pieces^,
        kind,
        mode.reference_start,
        mode.reference_end,
        mode.query_start,
        mode.query_end,
        mode.anchor == Anchor.END,
        band.low,
        band.high,
    )


def sweeps(mode: Mode) -> Bool:
    """Whether the mode takes the sweep, which takes no band: a local alignment, or free ends with a
    reward."""
    if mode.kind == Mode.SMITH_WATERMAN:
        return True
    return mode.kind == Mode.ENDS and mode.match_score > 0 and not mode.is_global()


@fieldwise_init
struct Case(Copyable, Movable, Writable):
    var reference: String
    var query: String
    var costs: Costs
    var mode: Mode
    var band: Band
    var ties: Ties
    var extended: Bool
    var memory: Int
    """The kept fronts' budget in bytes: the default, or a few hundred, which splits every pair."""


def check(trial: Case) raises:
    """Every property of one trial; raises on the first that fails."""
    var model = model_of(trial.costs, trial.mode, trial.band)
    var best = optimum(model, trial.reference, trial.query)
    var found: Alignment
    try:
        found = align(
            trial.reference,
            trial.query,
            trial.costs,
            trial.mode,
            band=trial.band,
            ties=trial.ties,
            extended=trial.extended,
        )
    except error:
        if not best:
            return
        raise Error(String("align raised ", error, " where the matrix found ", best.value()))
    if not best:
        raise Error("align found an alignment where the matrix found none inside the band")
    var earned = check_alignment(
        model,
        trial.reference,
        trial.query,
        found.cigar,
        found.reference_start,
        found.reference_end,
        found.query_start,
        found.query_end,
    )
    if earned.cost != found.cost:
        raise Error(String("cost ", found.cost, " but the CIGAR costs ", earned.cost))
    var expected_score = earned.score if trial.mode.is_scored() else -earned.cost
    if found.score != expected_score:
        raise Error(String("score ", found.score, " but the CIGAR earns ", expected_score))
    var scored = score(trial.reference, trial.query, trial.costs, trial.mode, band=trial.band)
    if scored != found.score:
        raise Error(String("score ", scored, ", align's ", found.score))
    if trial.mode.zdrop >= 0:
        # A Z-drop may give up short of the best: an extension that earns no more than it.
        if earned.score > best.value():
            raise Error(String("a Z-dropped extension earns ", earned.score, " past the best, ", best.value()))
        return
    if earned.score != best.value():
        raise Error(String("the alignment earns ", earned.score, ", the matrix's best is ", best.value()))
    # Free ends: the span the tie rule names.
    if model.kind == ENDS and not trial.mode.is_global():
        var span = rule_span(model, trial.reference, trial.query, trial.ties == Ties.LEFT)
        if (
            span[0] != found.reference_start
            or span[1] != found.reference_end
            or span[2] != found.query_start
            or span[3] != found.query_end
        ):
            raise Error(
                String(
                    "span ",
                    found.reference_start,
                    "..",
                    found.reference_end,
                    " x ",
                    found.query_start,
                    "..",
                    found.query_end,
                    " where the rule names ",
                    span[0],
                    "..",
                    span[1],
                    " x ",
                    span[2],
                    "..",
                    span[3],
                )
            )
    if not trial.mode.is_scored():
        var cost = distance(trial.reference, trial.query, trial.costs, trial.mode, band=trial.band)
        if cost != found.cost:
            raise Error(String("distance ", cost, ", align ", found.cost))
        var cap = draw(0, found.cost + 3)
        var capped = distance(trial.reference, trial.query, trial.costs, trial.mode, max_cost=cap, band=trial.band)
        if Bool(capped) != (found.cost <= cap) or (capped and capped.value() != found.cost):
            raise Error(String("distance under cap ", cap, " gave ", capped.or_else(-1), " for cost ", found.cost))
        var aligned = align(
            trial.reference,
            trial.query,
            trial.costs,
            trial.mode,
            max_cost=cap,
            band=trial.band,
            ties=trial.ties,
            extended=trial.extended,
            max_memory=trial.memory,
        )
        if Bool(aligned) != (found.cost <= cap):
            raise Error(String("align under cap ", cap, " for cost ", found.cost))
        if aligned and aligned.value().cost != found.cost:
            raise Error(String("align under a cap cost ", aligned.value().cost, ", without ", found.cost))
        # A split follows the tie rule within its pieces alone, so only an unsplit pair spells the same.
        if (
            aligned
            and trial.memory == DEFAULT_MAX_MEMORY
            and (
                aligned.value().cigar != found.cigar
                or aligned.value().reference_start != found.reference_start
                or aligned.value().query_start != found.query_start
            )
        ):
            raise Error(
                String(
                    "align under a cap spelled ",
                    aligned.value().cigar,
                    " from ",
                    aligned.value().reference_start,
                    ",",
                    aligned.value().query_start,
                    "; without, ",
                    found.cigar,
                    " from ",
                    found.reference_start,
                    ",",
                    found.query_start,
                )
            )
    # Ties: the left rule is the right one over both sequences reversed, read backwards, span and all;
    # a split follows it within its pieces alone.
    if trial.memory != DEFAULT_MAX_MEMORY:
        return
    if trial.mode.kind == Mode.SMITH_WATERMAN:
        var n = trial.reference.byte_length()
        var m = trial.query.byte_length()
        var other = Ties.RIGHT if trial.ties == Ties.LEFT else Ties.LEFT
        var mirrored = align(
            reversed_bytes(trial.reference),
            reversed_bytes(trial.query),
            trial.costs,
            trial.mode,
            ties=other,
            extended=trial.extended,
        )
        var empty = found.score == 0 and mirrored.score == 0
        if not empty and (
            reversed_cigar(mirrored.cigar) != found.cigar
            or n - mirrored.reference_end != found.reference_start
            or m - mirrored.query_end != found.query_start
        ):
            raise Error(String("mirrored local spelled ", mirrored.cigar, " against ", found.cigar))
    if trial.mode.is_global() and trial.band.covers(trial.reference.byte_length(), trial.query.byte_length()):
        var other = Ties.RIGHT if trial.ties == Ties.LEFT else Ties.LEFT
        var mirrored = align(
            reversed_bytes(trial.reference),
            reversed_bytes(trial.query),
            trial.costs,
            trial.mode,
            ties=other,
            extended=trial.extended,
        )
        if reversed_cigar(mirrored.cigar) != found.cigar:
            raise Error(String("mirrored ties spelled ", mirrored.cigar, " against ", found.cigar))


def reversed_cigar(cigar: String) -> String:
    var runs = List[String]()
    var start = 0
    var bytes = cigar.as_bytes()
    for index in range(len(bytes)):
        if bytes[index] < UInt8(ord("0")) or bytes[index] > UInt8(ord("9")):
            runs.append(String(StringSlice(unsafe_from_utf8=bytes[start : index + 1])))
            start = index + 1
    var out = String()
    for index in range(len(runs) - 1, -1, -1):
        out += runs[index]
    return out


def random_case() raises -> Case:
    var pair = sequences()
    var columns = pair[0].byte_length()
    var rows = pair[1].byte_length()
    var costs = random_costs()
    var mode = random_mode(columns, rows)
    var band = Band()
    if not sweeps(mode) and chance(0.3):
        var low = -draw(0, rows + 2)
        var high = draw(0, columns + 2)
        if chance(0.2):
            low = draw(-3, 3)
            high = low + draw(0, 4)
        band = Band(low, high)
    var memory = DEFAULT_MAX_MEMORY if chance(0.7) else draw(0, 2000)
    return Case(pair[0], pair[1], costs, mode, band, Ties.RIGHT if chance(0.5) else Ties.LEFT, chance(0.8), memory)


def main() raises:
    var arguments = argv()
    var iterations = Int(String(arguments[1])) if len(arguments) > 1 else 1000
    var start = Int(String(arguments[2])) if len(arguments) > 2 else 1
    var failures = 0
    var batch = List[Case]()
    for iteration in range(iterations):
        seed(start * 1_000_003 + iteration)
        var trial = random_case()
        try:
            check(trial)
        except error:
            failures += 1
            print("FAIL seed", start, "iteration", iteration, ":", error)
            print("  reference", trial.reference)
            print("  query    ", trial.query)
            print("  costs", trial.costs, "mode", trial.mode, "band", trial.band, "ties", trial.ties)
            if failures >= 10:
                break
        if not sweeps(trial.mode) and trial.band.covers(trial.reference.byte_length(), trial.query.byte_length()):
            if len(batch) == 0 or (batch[0].costs == trial.costs and batch[0].mode == trial.mode):
                batch.append(trial.copy())
        if len(batch) >= 16:
            var references = List[String]()
            var queries = List[String]()
            for item in batch:
                references.append(item.reference)
                queries.append(item.query)
            var together = alignments(references, queries, batch[0].costs, batch[0].mode, threads=3)
            for index in range(len(batch)):
                var alone = align(references[index], queries[index], batch[0].costs, batch[0].mode)
                if together[index].cigar != alone.cigar or together[index].cost != alone.cost:
                    failures += 1
                    print("FAIL batch at seed", start, "iteration", iteration, ": pair", index)
            batch.clear()
    if failures > 0:
        print(failures, "of", iterations, "cases failed")
        raise Error("fuzzing failed")
    print("fuzz:", iterations, "cases from seed", start, "agree with the full matrix")
