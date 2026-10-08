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
    Scoring,
    Ties,
    align,
    alignments,
    distance,
    score,
)

from oracle import ENDS, EXTENSION, LOCAL, Model, check_alignment, optimum, reversed_bytes, rule_span


def draw(low: Int, high: Int) -> Int:
    """A random integer in `low ..= high`."""
    return Int(random_ui64(UInt64(low), UInt64(high)))


def chance(probability: Float64) -> Bool:
    """True with probability `probability`."""
    return random_float64() < probability


def letters(length: Int, alphabet: String) -> String:
    """`length` random letters drawn from `alphabet`."""
    var symbols = alphabet.as_bytes()
    var out = List[UInt8](capacity=length)
    for _ in range(length):
        out.append(symbols[draw(0, len(symbols) - 1)])
    return String(unsafe_from_utf8=out^)


def mutated(text: String, rate: Float64, longest_gap: Int, alphabet: String) -> String:
    """`text` with about a `rate` of its letters substituted, or starting a deletion or followed by an
    insertion of up to `longest_gap` letters, a third each."""
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
    """Random `Costs` of any kind, a quarter of them with deletions priced apart."""
    var costs = symmetric_costs()
    if chance(0.25):
        # Deletions of their own, one piece or two, as bwa's -O del,ins.
        var opening = draw(0, 12)
        var extension = draw(1, 5)
        if chance(0.4):
            return costs.with_deletions(opening, extension, opening + draw(1, 30), draw(1, extension))
        return costs.with_deletions(opening, extension)
    return costs


def symmetric_costs() raises AlignmentError -> Costs:
    """Random unit, linear, affine or two-piece costs, insertions and deletions alike."""
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
    """A free-end allowance for a sequence of `length` letters: none, unbounded, or up to one past it."""
    var choice = draw(0, 3)
    if choice == 0:
        return 0
    if choice == 1:
        return 1 << 60
    return draw(0, length + 1)


def random_mode(columns: Int, rows: Int) raises AlignmentError -> Mode:
    """A random mode for a reference of `columns` letters and a query of `rows`, any kind, rewarded or not."""
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
    """The oracle's `Model` of `costs` under `mode` within `band`."""
    var pieces = List[Tuple[Int, Int]]()
    pieces.append((costs.opening, costs.extension))
    if costs.pieces() == 2:
        pieces.append((costs.opening2, costs.extension2))
    var deletions = List[Tuple[Int, Int]]()
    deletions.append((costs.deletion_opening, costs.deletion_extension))
    if costs.pieces() == 2:
        deletions.append((costs.deletion_opening2, costs.deletion_extension2))
    var kind = ENDS
    if mode.kind == Mode.EXTENSION:
        kind = EXTENSION
    elif mode.kind == Mode.SMITH_WATERMAN:
        kind = LOCAL
    return Model(
        Model.uniform(mode.match_score, costs.mismatch),
        deletions^,
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
    """One trial: a pair and everything it is aligned under."""

    var reference: String
    var query: String
    var costs: Costs
    var mode: Mode
    var band: Band
    var ties: Ties
    var eqx: Bool
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
            eqx=trial.eqx,
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
            eqx=trial.eqx,
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
            eqx=trial.eqx,
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
            eqx=trial.eqx,
        )
        if reversed_cigar(mirrored.cigar) != found.cigar:
            raise Error(String("mirrored ties spelled ", mirrored.cigar, " against ", found.cigar))


def reversed_cigar(cigar: String) -> String:
    """`cigar` with its runs in reverse order, each run kept whole: the CIGAR of both sequences reversed."""
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
    """A random trial; a band only for modes that take one, and a small memory budget three times in ten."""
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


def random_scoring(alphabet: String) raises -> Scoring:
    """A table over `alphabet`, mostly rewarding equal letters and charging unequal ones, at times a
    uniform one, and affine gap scores, a free extension among them."""
    var size = alphabet.byte_length()
    var table = List[Int8](length=size * size, fill=0)
    var uniform = chance(0.3)
    var reward = draw(0, 5)
    var penalty = draw(1, 6)
    for row in range(size):
        for column in range(size):
            if uniform:
                table[row * size + column] = Int8(reward if row == column else -penalty)
            elif row == column:
                table[row * size + column] = Int8(draw(0, 6))
            else:
                table[row * size + column] = Int8(draw(0, 9) - 7)
    var extension = -draw(0, 3) if chance(0.8) else 0
    var opening = -draw(0, 10)
    if opening + extension > extension or (extension == 0 and opening == 0):
        opening = -1
    return Scoring.tabulated(alphabet, table^, opening, extension)


def scoring_mode(columns: Int, rows: Int) raises AlignmentError -> Mode:
    """A random mode a `Scoring` takes for a reference of `columns` letters and a query of `rows`: no match
    score, the table's own rewards counting instead."""
    var kind = draw(0, 6)
    if kind == 0:
        return Mode.GLOBAL
    if kind == 1:
        return Mode.LOCAL
    if kind == 2:
        return [Mode.INFIX, Mode.PREFIX, Mode.SUFFIX, Mode.REFERENCE_IN_QUERY][draw(0, 3)]
    if kind == 3 or kind == 4:
        return Mode.ends_free(
            reference_start=allowance(columns),
            reference_end=allowance(columns),
            query_start=allowance(rows),
            query_end=allowance(rows),
        )
    var anchor = Anchor.END if chance(0.5) else Anchor.START
    if chance(0.3):
        return Mode.extension(0, anchor, zdrop=draw(0, 40))
    return Mode.extension(0, anchor)


def check_scoring(reference: String, query: String, scoring: Scoring, mode: Mode, eqx: Bool) raises:
    """A `Scoring`'s alignment in any mode: one the mode allows, earning its own score, the optimum,
    and `score` agreeing."""
    var size = scoring.alphabet_size()
    var letters = scoring.alphabet.as_bytes()
    var table = List[Int](length=256 * 256, fill=-1000)
    for row in range(size):
        for column in range(size):
            table[Int(letters[row]) * 256 + Int(letters[column])] = Int(scoring.substitutions[row * size + column])
    # The scores charge `open` for a gap's first letter and `extend` for each after it; as the model's
    # costs, a gap of `k` letters is `(extend - open) + k * -extend`.
    var opening = Int(scoring.gaps.extend - scoring.gaps.open)
    var extension = Int(-scoring.gaps.extend)
    var pieces: List[Tuple[Int, Int]] = [(opening, extension)]
    var kind = ENDS
    if mode.kind == Mode.EXTENSION:
        kind = EXTENSION
    elif mode.kind == Mode.SMITH_WATERMAN:
        kind = LOCAL
    var model = Model(
        table^,
        pieces.copy(),
        pieces^,
        kind,
        mode.reference_start,
        mode.reference_end,
        mode.query_start,
        mode.query_end,
        mode.anchor == Anchor.END,
        -(1 << 60),
        1 << 60,
    )
    var best = optimum(model, reference, query).value()
    var found = align(reference, query, scoring, mode, eqx=eqx)
    var earned = check_alignment(
        model,
        reference,
        query,
        found.cigar,
        found.reference_start,
        found.reference_end,
        found.query_start,
        found.query_end,
    )
    if earned.score != found.score:
        raise Error(String("a Scoring's score ", found.score, " but its CIGAR earns ", earned.score))
    if mode.zdrop >= 0:
        if earned.score > best:
            raise Error(String("a Z-dropped extension earns ", earned.score, " past the best, ", best))
    elif earned.score != best:
        raise Error(String("a Scoring's alignment earns ", earned.score, ", the matrix's best is ", best))
    var scored = score(reference, query, scoring, mode)
    if scored != found.score:
        raise Error(String("a Scoring's score ", scored, ", its alignment's ", found.score))
    if model.kind == ENDS and not mode.is_global():
        var span = rule_span(model, reference, query, True)
        if span[0] != found.reference_start or span[1] != found.reference_end or span[2] != found.query_start:
            raise Error("a Scoring's span is not the rule's")


def main() raises:
    """Runs `iterations` trials from `seed` (the arguments, 1000 and 1 by default), each with a `Scoring`'s
    turn, and every 16 trials alike a batch; raises when any failed."""
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
        # A Scoring's turn, over its own alphabet.
        var alphabet = "ACGT" if chance(0.7) else "ACGTN"
        var scoring = random_scoring(alphabet)
        var pair = (letters(draw(0, 60), alphabet), String())
        pair[1] = mutated(pair[0], [0.0, 0.1, 0.4][draw(0, 2)], 4, alphabet) if chance(0.7) else letters(
            draw(0, 60), alphabet
        )
        var mode = scoring_mode(pair[0].byte_length(), pair[1].byte_length())
        try:
            check_scoring(pair[0], pair[1], scoring, mode, chance(0.8))
        except error:
            failures += 1
            print("FAIL seed", start, "iteration", iteration, "(Scoring):", error)
            print("  reference", pair[0])
            print("  query    ", pair[1])
            var cells = String()
            for value in scoring.substitutions:
                cells += String(Int(value), ",")
            print("  table", scoring.alphabet, cells, scoring.gaps.open, scoring.gaps.extend, "mode", mode)
        # Trials sharing the first's costs and mode, with no band in the way, gather into a batch.
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
