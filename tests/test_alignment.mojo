"""
Tests for the alignment package, run with `pixi run test`.

Self-consistency alone would let two kernels agree while both were wrong, so most properties are
checked against something that shares no code with the dynamic programming:

- an exhaustive enumerator, which lists every alignment of a short pair and keeps the best;
- a rescorer, which prices a returned pair of gapped strings under the affine rule directly;
- DNA answers small enough to derive by hand, pinned here as literals.

Every property runs on the host. The device tests run only where an accelerator answers a real
alignment, and compare the device to the host pair by pair.
"""

from std.random import random_float64, random_ui64, seed
from std.testing import TestSuite, assert_equal, assert_false, assert_raises, assert_true

from dinara_align import (
    AlignmentMode,
    AlignmentResult,
    DNA_ALPHABET,
    affine_cigar,
    affine_cigars,
    affine_distance,
    Placement,
    Scoring,
    align,
    alignments,
    colorize,
    edit_alignment,
    edit_alignments,
    edit_cigar,
    edit_distance,
    edit_distances,
    edit_search,
    edit_search_alignment,
    levenshtein_alignment,
    needleman_wunsch_gotoh_alignment,
    needleman_wunsch_gotoh_score,
    score,
    scores,
    smith_waterman_gotoh_alignment,
    smith_waterman_gotoh_score,
)
from dinara_align.seeds import SEED_COLUMNS
from dinara_align.gap_affine import wavefront_align, wavefront_penalties
from dinara_align.vector_score import vector_score

comptime GLOBAL = AlignmentMode.GLOBAL
comptime LOCAL = AlignmentMode.LOCAL
comptime REPETITIONS = 20
"""Random draws per randomized test."""

# region Oracles


def random_sequence(shortest: Int, longest: Int, alphabet: String) -> String:
    """One sequence of a random length in `[shortest, longest]` over `alphabet`."""
    var letters = alphabet.as_bytes()
    var length = Int(random_ui64(UInt64(shortest), UInt64(longest)))
    var drawn = List[UInt8](capacity=length)
    for _ in range(length):
        drawn.append(letters[Int(random_ui64(0, UInt64(len(letters) - 1)))])
    return String(unsafe_from_utf8=drawn)


def expensive_gap() raises -> Scoring:
    """Uniform scores with an opening dearer than a mismatch, the regime most tests hold fixed."""
    return Scoring.uniform(5, -4, -20, -1)


def scoring_regimes() raises -> List[Scoring]:
    """One representative per regime: an expensive gap, a cheap one, unit costs, a free extension, and the DNA default.
    """
    var regimes = List[Scoring]()
    regimes.append(Scoring.uniform(5, -4, -20, -1))
    regimes.append(Scoring.uniform(2, -1, -2, -1))
    regimes.append(Scoring.uniform(0, -1, -1, -1))
    regimes.append(Scoring.uniform(1, -1, -5, 0))
    regimes.append(Scoring.dna())
    return regimes^


def substitution(scoring: Scoring, left: UInt8, right: UInt8) -> Int:
    """The table entry for two letters, looked up by position in the alphabet."""
    var letters = scoring.alphabet.as_bytes()
    var row = 0
    var column = 0
    for index in range(len(letters)):
        if letters[index] == left:
            row = index
        if letters[index] == right:
            column = index
    return Int(scoring.substitutions[row * scoring.alphabet_size() + column])


def rescore(first_gapped: String, second_gapped: String, scoring: Scoring) -> Int:
    """Prices a gapped pair under the affine rule, sharing no code with the dynamic programming."""
    comptime GAP = UInt8(ord("-"))
    var top = first_gapped.as_bytes()
    var bottom = second_gapped.as_bytes()
    var total = 0
    var in_first = False
    var in_second = False
    for column in range(len(top)):
        if top[column] == GAP and bottom[column] == GAP:
            continue
        if top[column] == GAP:
            total += Int(scoring.gaps.extend) if in_first else Int(scoring.gaps.open)
            in_first = True
            in_second = False
        elif bottom[column] == GAP:
            total += Int(scoring.gaps.extend) if in_second else Int(scoring.gaps.open)
            in_first = False
            in_second = True
        else:
            total += substitution(scoring, top[column], bottom[column])
            in_first = False
            in_second = False
    return total


def best_extension(
    first: List[UInt8], second: List[UInt8], i: Int, j: Int, top: List[UInt8], bottom: List[UInt8], scoring: Scoring
) -> Int:
    """The best rescored alignment completing `top` and `bottom` from position `(i, j)`."""
    comptime GAP = UInt8(ord("-"))
    if i == len(first) and j == len(second):
        return rescore(String(unsafe_from_utf8=top), String(unsafe_from_utf8=bottom), scoring)
    var best = Int.MIN
    if i < len(first) and j < len(second):
        var t = top.copy()
        var b = bottom.copy()
        t.append(first[i])
        b.append(second[j])
        best = max(best, best_extension(first, second, i + 1, j + 1, t, b, scoring))
    if i < len(first):
        var t = top.copy()
        var b = bottom.copy()
        t.append(first[i])
        b.append(GAP)
        best = max(best, best_extension(first, second, i + 1, j, t, b, scoring))
    if j < len(second):
        var t = top.copy()
        var b = bottom.copy()
        t.append(GAP)
        b.append(second[j])
        best = max(best, best_extension(first, second, i, j + 1, t, b, scoring))
    return best


def best_enumerated(first: List[UInt8], second: List[UInt8], scoring: Scoring) -> Int:
    """The best global score over every alignment, by the three-way edit recursion with no table.

    Each alignment is built out in full and rescored, so a recurrence that is self-consistently
    wrong on every backend at once still disagrees with this.
    """
    return best_extension(first, second, 0, 0, List[UInt8](), List[UInt8](), scoring)


def slice_bytes(text: List[UInt8], start: Int, stop: Int) -> List[UInt8]:
    """The bytes in `[start, stop)`."""
    var piece = List[UInt8](capacity=stop - start)
    for index in range(start, stop):
        piece.append(text[index])
    return piece^


def brute_optimum[mode: AlignmentMode](first: String, second: String, scoring: Scoring) -> Int:
    """Global: the best enumerated alignment. Local: the best over every pair of substrings, or zero."""
    var top = List[UInt8](first.as_bytes())
    var bottom = List[UInt8](second.as_bytes())
    comptime if mode == AlignmentMode.GLOBAL:
        return best_enumerated(top, bottom, scoring)
    var best = 0
    for start in range(len(top)):
        for stop in range(start + 1, len(top) + 1):
            for low in range(len(bottom)):
                for high in range(low + 1, len(bottom) + 1):
                    best = max(
                        best, best_enumerated(slice_bytes(top, start, stop), slice_bytes(bottom, low, high), scoring)
                    )
    return best


def assert_well_formed[
    mode: AlignmentMode
](first: String, second: String, produced: AlignmentResult, scoring: Scoring) raises:
    """Both rows have one length, each rebuilds its input or a piece of it, and they earn their score."""
    assert_equal(produced.first_gapped.byte_length(), produced.second_gapped.byte_length())
    var core_first = produced.first_gapped.replace("-", "")
    var core_second = produced.second_gapped.replace("-", "")
    comptime if mode == AlignmentMode.GLOBAL:
        assert_equal(core_first, first)
        assert_equal(core_second, second)
    else:
        assert_true(core_first in first, "a local row is not a substring of its input")
        assert_true(core_second in second, "a local row is not a substring of its input")
    assert_equal(rescore(produced.first_gapped, produced.second_gapped, scoring), Int(produced.score))


def gpu_available() raises -> Bool:
    """Whether an accelerator serves a real alignment here, which a successful import does not prove."""
    var scoring = Scoring.dna()
    try:
        _ = score[GLOBAL]("AC", "CA", scoring, Placement.on_gpu(0, 1))
        return True
    except:
        return False


# endregion Oracles

# region Known Answers


def test_dna_default_is_minimap2() raises:
    """Match 2, mismatch -4, and minimap2's `-O4 -E2` gap, which charges the first gapped base six."""
    var dna = Scoring.dna()
    assert_equal(dna.alphabet, "ACGT")
    for left in [UInt8(ord("A")), UInt8(ord("C")), UInt8(ord("G")), UInt8(ord("T"))]:
        for right in [UInt8(ord("A")), UInt8(ord("C")), UInt8(ord("G")), UInt8(ord("T"))]:
            assert_equal(substitution(dna, left, right), 2 if left == right else -4)
    assert_equal(Int(dna.gaps.open), -6)
    assert_equal(Int(dna.gaps.extend), -2)


def test_hand_computed_global() raises:
    """Pairs whose optimum can be worked out on paper, strings included."""
    var dna = Scoring.dna()
    var same = needleman_wunsch_gotoh_alignment("ACGTACGT", "ACGTACGT", dna)
    assert_equal(same.first_gapped, "ACGTACGT")
    assert_equal(same.score, 16)

    # One substitution costs -4, which beats opening two gaps at -6 each.
    var substituted = needleman_wunsch_gotoh_alignment("ACGTACGT", "ACGTTCGT", dna)
    assert_equal(substituted.first_gapped, "ACGTACGT")
    assert_equal(substituted.second_gapped, "ACGTTCGT")
    assert_equal(substituted.score, 7 * 2 - 4)

    # A three-base deletion is one run: -6 for its first base and -2 for each of the other two.
    var deleted = needleman_wunsch_gotoh_alignment("ACGTTGCAGGGCATGACGT", "ACGTTGCACATGACGT", dna)
    assert_equal(deleted.first_gapped, "ACGTTGCAGGGCATGACGT")
    assert_equal(deleted.second_gapped, "ACGTTGCA---CATGACGT")
    assert_equal(deleted.score, 16 * 2 - 6 - 2 * 2)
    assert_equal(needleman_wunsch_gotoh_score("ACGTTGCAGGGCATGACGT", "ACGTTGCACATGACGT", dna), 22)

    # Against nothing, the whole sequence is one gap run.
    assert_equal(needleman_wunsch_gotoh_score("AAAA", "", dna), -6 - 3 * 2)
    assert_equal(needleman_wunsch_gotoh_score("", "", dna), 0)


def test_hand_computed_local() raises:
    """The shared core of two otherwise unrelated sequences, trimmed at both ends."""
    var aligned = smith_waterman_gotoh_alignment("GGGGACGTACGTGGGG", "CCCCACGTACGTCCCC", Scoring.dna())
    assert_equal(aligned.first_gapped, "ACGTACGT")
    assert_equal(aligned.second_gapped, "ACGTACGT")
    assert_equal(aligned.score, 16)
    assert_equal(smith_waterman_gotoh_score("GGGGACGTACGTGGGG", "CCCCACGTACGTCCCC", Scoring.dna()), 16)


def test_levenshtein_known_distances() raises:
    """Edit distances that can be checked by eye, with rows that rebuild both inputs."""
    var cases: List[Tuple[String, String, Int]] = [
        ("ACGT", "ACGT", 0),
        ("ACGT", "AGT", 1),
        ("ACGT", "", 4),
        ("", "ACG", 3),
        ("AAAA", "TTTT", 4),
        # A rotation by one: every position mismatches, so a deletion and an insertion are cheaper.
        ("ACGTACGT", "TACGTACG", 2),
    ]
    for example in cases:
        var aligned = levenshtein_alignment(example[0], example[1])
        assert_equal(Int(aligned.score), example[2])
        assert_equal(aligned.first_gapped.replace("-", ""), example[0])
        assert_equal(aligned.second_gapped.replace("-", ""), example[1])


# endregion Known Answers

# region Exhaustive Oracle


def check_against_enumeration[mode: AlignmentMode]() raises:
    """Every pair over a two-letter alphabet with combined length at most five."""
    var scoring = expensive_gap()
    var words = List[String]()
    words.append("")
    var frontier: List[String] = [""]
    for _ in range(3):
        var grown = List[String]()
        for word in frontier:
            grown.append(word + "A")
            grown.append(word + "C")
        words.extend(grown.copy())
        frontier = grown^
    for first in words:
        for second in words:
            if first.byte_length() + second.byte_length() > 5:
                continue
            var expected = brute_optimum[mode](first, second, scoring)
            assert_equal(Int(score[mode](first, second, scoring, Placement.on_cpu(1))), expected)
            var produced = align[mode](first, second, scoring, Placement.on_cpu(1))
            assert_equal(Int(produced.score), expected)


def test_global_matches_enumeration() raises:
    check_against_enumeration[GLOBAL]()


def test_local_matches_enumeration() raises:
    check_against_enumeration[LOCAL]()


# endregion Exhaustive Oracle

# region Properties


def check_well_formed[mode: AlignmentMode]() raises:
    """Every path is well formed, realizes its own score, and agrees with the score-only kernel."""
    seed(1)
    for scoring in scoring_regimes():
        for _ in range(REPETITIONS):
            var first = random_sequence(5, 25, DNA_ALPHABET)
            var second = random_sequence(5, 25, DNA_ALPHABET)
            var produced = align[mode](first, second, scoring)
            assert_well_formed[mode](first, second, produced, scoring)
            assert_equal(score[mode](first, second, scoring), produced.score)


def test_global_output_is_well_formed() raises:
    check_well_formed[GLOBAL]()


def test_local_output_is_well_formed() raises:
    check_well_formed[LOCAL]()


def check_linear_matches_stored[mode: AlignmentMode]() raises:
    """Both traceback strategies reach the same score, each with a path that earns it.

    Ties may land differently under a divide-and-conquer join, so the strings need not match.
    """
    seed(2)
    for scoring in scoring_regimes():
        for _ in range(3):
            var first = random_sequence(140, 260, DNA_ALPHABET)
            var second = random_sequence(140, 260, DNA_ALPHABET)
            var stored = align[mode](first, second, scoring, Placement.default(), 10**12)
            var linear = align[mode](first, second, scoring, Placement.default(), 0)
            assert_equal(stored.score, linear.score)
            assert_equal(stored.score, score[mode](first, second, scoring))
            assert_well_formed[mode](first, second, stored, scoring)
            assert_well_formed[mode](first, second, linear, scoring)


def test_global_linear_matches_stored() raises:
    check_linear_matches_stored[GLOBAL]()


def test_local_linear_matches_stored() raises:
    check_linear_matches_stored[LOCAL]()


def test_linear_space_carries_a_long_pair() raises:
    """Three thousand by three thousand with no stored matrix at all."""
    seed(3)
    var scoring = expensive_gap()
    var first = random_sequence(3000, 3000, DNA_ALPHABET)
    var second = random_sequence(3000, 3000, DNA_ALPHABET)
    assert_well_formed[GLOBAL](first, second, align[GLOBAL](first, second, scoring, Placement.default(), 0), scoring)
    assert_well_formed[LOCAL](first, second, align[LOCAL](first, second, scoring, Placement.default(), 0), scoring)


def test_symmetry() raises:
    """Swapping the arguments does not change the score."""
    seed(4)
    var scoring = expensive_gap()
    for _ in range(REPETITIONS):
        var first = random_sequence(5, 25, DNA_ALPHABET)
        var second = random_sequence(5, 25, DNA_ALPHABET)
        assert_equal(score[GLOBAL](first, second, scoring), score[GLOBAL](second, first, scoring))
        assert_equal(score[LOCAL](first, second, scoring), score[LOCAL](second, first, scoring))


def test_levenshtein_is_the_unit_cost_limit() raises:
    """At unit costs the global recurrence is the negated edit distance."""
    seed(5)
    var unit = Scoring.edit_distance()
    for _ in range(REPETITIONS):
        var first = random_sequence(3, 15, DNA_ALPHABET)
        var second = random_sequence(3, 15, DNA_ALPHABET)
        assert_equal(-score[GLOBAL](first, second, unit), levenshtein_alignment(first, second).score)


def test_local_never_scores_below_global() raises:
    """A global path is also a local candidate, and the empty window is always available."""
    seed(6)
    var scoring = expensive_gap()
    for _ in range(REPETITIONS):
        var first = random_sequence(5, 25, DNA_ALPHABET)
        var second = random_sequence(5, 25, DNA_ALPHABET)
        var local = score[LOCAL](first, second, scoring)
        assert_true(local >= 0)
        assert_true(local >= score[GLOBAL](first, second, scoring))


def test_optimum_falls_as_gaps_get_harsher() raises:
    """A harsher opening lowers every feasible alignment, so the maximum cannot rise."""
    seed(7)
    for _ in range(REPETITIONS):
        var first = random_sequence(5, 25, DNA_ALPHABET)
        var second = random_sequence(5, 25, DNA_ALPHABET)
        var previous = Int32.MAX
        for opening in [-2, -5, -10, -20, -40]:
            var current = score[GLOBAL](first, second, Scoring.uniform(5, -4, opening, -1))
            assert_true(current <= previous, "a harsher gap raised the optimum")
            previous = current


def test_free_extension_ignores_gap_width() raises:
    """A gap that costs nothing to extend costs the same however wide it is."""
    seed(8)
    var free_extension = Scoring.uniform(5, -10, -1, 0, "ACGTN")
    for _ in range(REPETITIONS):
        var first = random_sequence(5, 15, "ACGT")
        var second = random_sequence(5, 15, "ACGT")
        var bytes = second.as_bytes()
        var cut = len(bytes) // 2
        var head = String(unsafe_from_utf8=slice_bytes(List[UInt8](bytes), 0, cut))
        var tail = String(unsafe_from_utf8=slice_bytes(List[UInt8](bytes), cut, len(bytes)))
        var reference = score[GLOBAL](first, head + "N" + tail, free_extension)
        for width in range(2, 6):
            assert_equal(score[GLOBAL](first, head + "N" * width + tail, free_extension), reference)


def test_table_matches_the_uniform_costs_it_spells() raises:
    """A tabulated diagonal scores exactly what the uniform record scores."""
    seed(9)
    var table = List[Int8](length=16, fill=-1)
    for index in range(4):
        table[index * 4 + index] = 2
    var tabulated = Scoring.tabulated("ACGT", table^, -5, -1)
    var uniform = Scoring.uniform(2, -1, -5, -1, "ACGT")
    for _ in range(REPETITIONS):
        var first = random_sequence(10, 40, "ACGT")
        var second = random_sequence(10, 40, "ACGT")
        assert_equal(score[GLOBAL](first, second, tabulated), score[GLOBAL](first, second, uniform))
        assert_equal(score[LOCAL](first, second, tabulated), score[LOCAL](first, second, uniform))


def test_batch_matches_single_pairs() raises:
    """The batched entry points reproduce the single-pair ones, and an empty batch answers empty."""
    seed(10)
    var scoring = expensive_gap()
    var firsts = List[String]()
    var seconds = List[String]()
    for _ in range(24):
        firsts.append(random_sequence(5, 40, DNA_ALPHABET))
        seconds.append(random_sequence(5, 40, DNA_ALPHABET))
    var batch_scores = scores[GLOBAL](firsts, seconds, scoring)
    var batch_alignments = alignments[LOCAL](firsts, seconds, scoring)
    for index in range(len(firsts)):
        assert_equal(batch_scores[index], score[GLOBAL](firsts[index], seconds[index], scoring))
        var single = align[LOCAL](firsts[index], seconds[index], scoring)
        assert_equal(batch_alignments[index].first_gapped, single.first_gapped)
        assert_equal(batch_alignments[index].second_gapped, single.second_gapped)
        assert_equal(batch_alignments[index].score, single.score)
    assert_equal(len(scores[GLOBAL](List[String](), List[String](), scoring)), 0)


def test_bit_parallel_edit_distance_matches_the_full_matrix() raises:
    """The bit-parallel sweep returns the distance the cell-by-cell traceback reports.

    Lengths straddle a word (63, 64, 65), a four-word block (255, 256, 257) and the eight columns the
    staggered block needs, so every edge path runs: the triangles, the word left over after the last
    block, and a last word only partly filled.
    """
    seed(12)
    var lengths: List[Int] = [0, 1, 2, 7, 8, 9, 63, 64, 65, 127, 128, 129, 255, 256, 257, 300, 513]
    for first_length in lengths:
        for second_length in lengths:
            var first = random_sequence(first_length, first_length, DNA_ALPHABET)
            var second = random_sequence(second_length, second_length, DNA_ALPHABET)
            assert_equal(edit_distance(first, second), Int(levenshtein_alignment(first, second).score))
    # Similar pairs too, where long runs of matches carry through whole words.
    for _ in range(REPETITIONS):
        var first = random_sequence(500, 900, DNA_ALPHABET)
        var second = first
        var bytes = List[UInt8](second.as_bytes())
        for _ in range(10):
            bytes[Int(random_ui64(0, UInt64(len(bytes) - 1)))] = UInt8(ord("A"))
        second = String(unsafe_from_utf8=bytes)
        assert_equal(edit_distance(first, second), Int(levenshtein_alignment(first, second).score))
    with assert_raises(contains="symbols"):
        _ = edit_distance("ACGT", "ACGTNRYKM")


def test_bit_parallel_symbols_past_acgt() raises:
    """Bytes past `ACGT` are symbols of their own, each matching only itself, up to four of them.

    `N` alone and four extra symbols, short pairs and long ones whose band sweeps with seeds,
    close and divergent, and pairs where only one side holds them, either side; every distance is the global
    affine alignment's at unit costs over the same alphabet, which shares no code with the sweep,
    and every alignment rebuilds both inputs and rescores to the distance.
    """
    seed(19)
    var alphabets: List[String] = ["ACGTN", "ACGTNRYK"]
    for alphabet in alphabets:
        var unit = Scoring.edit_distance(alphabet)
        for length in [1, 2, 64, 65, 700, 3000, SEED_COLUMNS + 616]:
            for rate in [0.0, 0.05, 0.2]:
                var first = random_sequence(length, length, alphabet)
                var second = mutate(first, rate)
                var plain = random_sequence(length, length, DNA_ALPHABET)
                var pairs: List[Tuple[String, String]] = [
                    (first, second),
                    (second, first),
                    (first, plain),
                    (plain, first),
                ]
                # One long pair an alphabet: past `SEED_COLUMNS`, where bases alone would take seeds.
                if length > SEED_COLUMNS:
                    if rate != 0.05:
                        continue
                    pairs = [(first, second)]
                for pair in pairs:
                    var expected = -Int(score[GLOBAL](pair[0], pair[1], unit))
                    assert_equal(edit_distance(pair[0], pair[1]), expected)
                    var aligned = edit_alignment(pair[0], pair[1])
                    assert_equal(Int(aligned.score), expected)
                    assert_equal(aligned.first_gapped.replace("-", ""), pair[0])
                    assert_equal(aligned.second_gapped.replace("-", ""), pair[1])
                    assert_equal(rescore(aligned.first_gapped, aligned.second_gapped, unit), -expected)


def test_bit_parallel_whole_matrix_is_exact() raises:
    """Unrelated pairs, which no band narrows, are swept whole and keep the exact distance.

    Shapes are deliberately ragged, rows that do not fill a block of words, and each distance is the
    global affine alignment's at unit costs, which shares no code with the sweep.
    """
    seed(13)
    var unit = Scoring.edit_distance()
    var shapes: List[Tuple[Int, Int]] = [(9000, 8100), (8100, 9001), (12345, 6789), (20000, 3300), (4100, 30000)]
    for shape in shapes:
        var first = random_sequence(shape[0], shape[0], DNA_ALPHABET)
        var second = random_sequence(shape[1], shape[1], DNA_ALPHABET)
        assert_equal(edit_distance(first, second), -Int(score[GLOBAL](first, second, unit)))
    var first = random_sequence(9000, 9000, DNA_ALPHABET)
    var bytes = List[UInt8](first.as_bytes())
    for _ in range(900):
        bytes[Int(random_ui64(0, UInt64(len(bytes) - 1)))] = UInt8(ord("G"))
    var second = String(unsafe_from_utf8=bytes)
    assert_equal(edit_distance(first, second), Int(levenshtein_alignment(first, second).score))


def mutate(text: String, rate: Float64) -> String:
    """Substitutions, deletions and insertions in equal thirds, at `rate` edits per base."""
    var letters = String(DNA_ALPHABET).as_bytes()
    var out = List[UInt8]()
    for byte in text.as_bytes():
        var roll = random_float64()
        if roll < rate / 3:
            out.append(letters[Int(random_ui64(0, 3))])
        elif roll < 2 * rate / 3:
            continue
        elif roll < rate:
            out.append(byte)
            out.append(letters[Int(random_ui64(0, 3))])
        else:
            out.append(byte)
    return String(unsafe_from_utf8=out)


def test_bit_parallel_band_doubling_is_exact() raises:
    """Band doubling with pruning returns the full matrix's distance on pairs it actually narrows.

    Similar pairs at several divergences take the banded rounds; a block cut out of the middle
    gives a large length difference either way; unrelated pairs exhaust the band and fall back to
    the whole matrix. Every one is checked against the cell-by-cell distance.
    """
    seed(14)
    for length in [2000, 5000]:
        for rate in [0.0, 0.001, 0.01, 0.05, 0.2]:
            var first = random_sequence(length, length, DNA_ALPHABET)
            var second = mutate(first, rate)
            var bytes = List[UInt8](first.as_bytes())
            var cut = String(unsafe_from_utf8=slice_bytes(bytes, 0, length // 3)) + String(
                unsafe_from_utf8=slice_bytes(bytes, length // 2, length)
            )
            var pairs: List[Tuple[String, String]] = [(first, second), (second, first), (first, cut), (cut, first)]
            for pair in pairs:
                assert_equal(edit_distance(pair[0], pair[1]), Int(levenshtein_alignment(pair[0], pair[1]).score))
    var unrelated_first = random_sequence(3000, 3000, DNA_ALPHABET)
    var unrelated_second = random_sequence(2800, 2800, DNA_ALPHABET)
    assert_equal(
        edit_distance(unrelated_first, unrelated_second),
        Int(levenshtein_alignment(unrelated_first, unrelated_second).score),
    )


def test_edit_alignment_is_an_optimal_alignment() raises:
    """The bit-parallel traceback returns rows that rebuild both inputs and cost exactly the distance.

    Lengths from empty to a few kbp, divergences from identical to half the bases edited, and a large length difference either way; every distance is
    the cell-by-cell one, and every pair of rows is rescored independently at unit cost.
    """
    seed(15)
    for length in [0, 1, 64, 65, 700, 3000, 6000]:
        for rate in [0.0, 0.01, 0.05, 0.2, 0.5]:
            var first = random_sequence(length, length, DNA_ALPHABET)
            var second = mutate(first, rate)
            var pairs: List[Tuple[String, String]] = [(first, second), (second, first)]
            if length >= 700:
                var bytes = List[UInt8](first.as_bytes())
                var cut = String(unsafe_from_utf8=slice_bytes(bytes, 0, length // 3)) + String(
                    unsafe_from_utf8=slice_bytes(bytes, length // 2, length)
                )
                pairs.append((first, cut))
                pairs.append((cut, first))
            for pair in pairs:
                var aligned = edit_alignment(pair[0], pair[1])
                assert_equal(Int(aligned.score), Int(levenshtein_alignment(pair[0], pair[1]).score))
                assert_equal(aligned.first_gapped.replace("-", ""), pair[0])
                assert_equal(aligned.second_gapped.replace("-", ""), pair[1])
                assert_equal(
                    rescore(aligned.first_gapped, aligned.second_gapped, Scoring.edit_distance()), -Int(aligned.score)
                )


def semi_global_distance(pattern: String, text: String, prefix: Bool) -> Int:
    """The textbook dynamic program for the pattern against any part of the text, or with `prefix`
    against a prefix of it: the top row free, or the global border, and the least of the last row."""
    var p = pattern.as_bytes()
    var t = text.as_bytes()
    var previous = List[Int](length=len(t) + 1, fill=0)
    for j in range(len(t) + 1):
        previous[j] = j if prefix else 0
    for i in range(1, len(p) + 1):
        var current = List[Int](length=len(t) + 1, fill=i)
        for j in range(1, len(t) + 1):
            var diagonal = previous[j - 1] + (0 if p[i - 1] == t[j - 1] else 1)
            current[j] = min(diagonal, min(previous[j] + 1, current[j - 1] + 1))
        previous = current^
    var best = previous[0]
    for value in previous:
        best = min(best, value)
    return best


def test_edit_search_matches_the_dynamic_program() raises:
    """A pattern found inside a text, or at its start, at the distance the textbook dynamic program
    gives, and aligned to the part of the text reported, which the alignment rebuilds.

    Patterns across word boundaries, from empty to longer than the text, planted mutated copies of a
    piece of the text and unrelated ones, with an `N` among the symbols too.
    """
    seed(23)
    var unit = Scoring.edit_distance("ACGTN")
    for pattern_length in [0, 1, 63, 64, 65, 129, 300]:
        for text_length in [0, 1, 200, 700]:
            for rate in [0.0, 0.05, 0.2]:
                for prefix in [False, True]:
                    var text = random_sequence(text_length, text_length, "ACGTN")
                    var pattern = random_sequence(pattern_length, pattern_length, DNA_ALPHABET)
                    if 0 < pattern_length <= text_length:
                        var start = Int(random_ui64(0, UInt64(text_length - pattern_length)))
                        var bytes = List[UInt8](text.as_bytes())
                        pattern = mutate(
                            String(unsafe_from_utf8=slice_bytes(bytes, start, start + pattern_length)), rate
                        )
                    var expected = semi_global_distance(pattern, text, prefix)
                    var hit = edit_search(pattern, text, prefix)
                    assert_equal(hit.distance, expected)
                    if prefix:
                        assert_equal(hit.start, 0)
                    var found = edit_search_alignment(pattern, text, prefix)
                    var aligned = found[1].copy()
                    var part = String(StringSlice(unsafe_from_utf8=text.as_bytes()[hit.start : hit.end]))
                    assert_equal(Int(aligned.score), expected)
                    assert_equal(aligned.first_gapped.replace("-", ""), part)
                    assert_equal(aligned.second_gapped.replace("-", ""), pattern)
                    assert_equal(rescore(aligned.first_gapped, aligned.second_gapped, unit), -expected)


def test_edit_batches_match_single_pairs() raises:
    """A batch spread over threads answers every pair as the single call does, in order.

    Empty pairs, short and long ones, close and divergent, and one long enough for seeds, on one
    thread and on several; every alignment rebuilds both inputs. A bad base or sides of different
    lengths are refused, as a serial loop would refuse them.
    """
    seed(17)
    var firsts = List[String]()
    var seconds = List[String]()
    firsts.append(String())
    seconds.append(String())
    firsts.append(String("ACGT"))
    seconds.append(String())
    for length in [1, 100, 700, 3000, 6000, 20000]:
        for rate in [0.0, 0.05, 0.15, 0.3]:
            var first = random_sequence(length, length, DNA_ALPHABET)
            firsts.append(first)
            seconds.append(mutate(first, rate))
    for threads in [1, 8]:
        var distances = edit_distances(firsts, seconds, threads)
        var aligned = edit_alignments(firsts, seconds, threads)
        assert_equal(len(distances), len(firsts))
        assert_equal(len(aligned), len(firsts))
        for index in range(len(firsts)):
            var expected = edit_distance(firsts[index], seconds[index])
            assert_equal(distances[index], expected)
            assert_equal(Int(aligned[index].score), expected)
            assert_equal(aligned[index].first_gapped.replace("-", ""), firsts[index])
            assert_equal(aligned[index].second_gapped.replace("-", ""), seconds[index])
    var bad_firsts: List[String] = ["ACGT", "ACGT", "ACGTNRYKM"]
    var bad_seconds: List[String] = ["ACGA", "ACG", "ACGT"]
    with assert_raises():
        _ = edit_distances(bad_firsts, bad_seconds, 8)
    with assert_raises():
        _ = edit_alignments(bad_firsts, bad_seconds, 8)
    var short: List[String] = ["ACGT"]
    with assert_raises():
        _ = edit_distances(bad_firsts, short)


def test_seeded_bands_are_exact() raises:
    """Pairs long enough for the seed heuristic, exact seeds and inexact, keep the exact distance.

    From `SEED_COLUMNS`, 16 kbp but for AVX-512's 86, a band prunes with seeds; once few exact seeds
    chain, past about one edit in fifteen bases, they are rebuilt to match within one edit.
    Divergences either side of that, with errors also gathered into a burst at one end as real reads
    carry them; every distance is the global wavefront's at unit costs, which shares no code with the
    band, and every alignment is rescored independently.
    """
    seed(16)
    var unit = Scoring.edit_distance()
    var shortest = max(20000, SEED_COLUMNS + 2000)
    for rate in [0.05, 0.1, 0.15, 0.2, 0.3]:
        var first = random_sequence(shortest, shortest + 4000, DNA_ALPHABET)
        var second = mutate(first, rate)
        var bytes = List[UInt8](second.as_bytes())
        var tail = len(bytes) - len(bytes) // 20
        var burst = String(unsafe_from_utf8=slice_bytes(bytes, 0, tail)) + mutate(
            String(unsafe_from_utf8=slice_bytes(bytes, tail, len(bytes))), 0.5
        )
        var pairs: List[Tuple[String, String]] = [(first, second), (second, first), (first, burst)]
        for pair in pairs:
            var expected = -Int(score[GLOBAL](pair[0], pair[1], unit))
            assert_equal(edit_distance(pair[0], pair[1]), expected)
            var aligned = edit_alignment(pair[0], pair[1])
            assert_equal(Int(aligned.score), expected)
            assert_equal(aligned.first_gapped.replace("-", ""), pair[0])
            assert_equal(aligned.second_gapped.replace("-", ""), pair[1])
            assert_equal(rescore(aligned.first_gapped, aligned.second_gapped, unit), -expected)


def rows_from_cigar(first: String, second: String, cigar: String) raises -> Tuple[String, String, Int]:
    """The gapped rows a CIGAR spells over its two sequences, and its count of edits; raises on an
    entry that claims a match between differing bases, or a mismatch between equal ones."""
    var top = List[UInt8]()
    var bottom = List[UInt8]()
    var a = first.as_bytes()
    var b = second.as_bytes()
    var column = 0
    var row = 0
    var edits = 0
    var length = 0
    for byte in cigar.as_bytes():
        if byte >= UInt8(ord("0")) and byte <= UInt8(ord("9")):
            length = length * 10 + Int(byte - UInt8(ord("0")))
            continue
        assert_true(length > 0)
        for _ in range(length):
            if byte == UInt8(ord("D")):
                top.append(a[column])
                bottom.append(UInt8(ord("-")))
                column += 1
                edits += 1
            elif byte == UInt8(ord("I")):
                top.append(UInt8(ord("-")))
                bottom.append(b[row])
                row += 1
                edits += 1
            else:
                if byte == UInt8(ord("=")):
                    assert_equal(a[column], b[row])
                elif byte == UInt8(ord("X")):
                    assert_true(a[column] != b[row])
                    edits += 1
                else:
                    assert_equal(byte, UInt8(ord("M")))
                top.append(a[column])
                bottom.append(b[row])
                column += 1
                row += 1
        length = 0
    assert_equal(length, 0)
    assert_equal(column, len(a))
    assert_equal(row, len(b))
    return (String(unsafe_from_utf8=top), String(unsafe_from_utf8=bottom), edits)


def test_edit_cigar_spells_the_alignment() raises:
    """`edit_cigar` spells the very alignment `edit_alignment` writes out, run by run.

    Pairs settled by one diagonal front, by two, and by a band with seeds, with symbols past `ACGT`,
    and with empty sides; the CIGAR rebuilds both gapped rows, every `=` joins equal bases and every
    `X` differing ones, its substitutions and gaps number the distance, and with `M` for both the
    same runs merge.
    """
    seed(29)
    var pairs: List[Tuple[String, String]] = [("", ""), ("", "ACG"), ("ACG", ""), ("ACGT", "ACGT")]
    var alphabets: List[String] = [DNA_ALPHABET, "ACGTN"]
    for alphabet in alphabets:
        for length in [1, 40, 700, 3000, 20000]:
            for rate in [0.0, 0.02, 0.15, 0.4]:
                var first = random_sequence(length, length, alphabet)
                pairs.append((first, mutate(first, rate)))
    for pair in pairs:
        var aligned = edit_alignment(pair[0], pair[1])
        var spelled = edit_cigar(pair[0], pair[1])
        assert_equal(spelled.distance, Int(aligned.score))
        var rows = rows_from_cigar(pair[0], pair[1], spelled.cigar)
        assert_equal(rows[0], aligned.first_gapped)
        assert_equal(rows[1], aligned.second_gapped)
        assert_equal(rows[2], spelled.distance)
        var plain = edit_cigar(pair[0], pair[1], extended=False)
        assert_equal(plain.distance, spelled.distance)
        var plain_rows = rows_from_cigar(pair[0], pair[1], plain.cigar)
        assert_equal(plain_rows[0], aligned.first_gapped)
        assert_equal(plain_rows[1], aligned.second_gapped)
        assert_false("=" in plain.cigar or "X" in plain.cigar)


def sprinkle(text: String, rate: Float64, symbol: String) -> String:
    """`text` with each byte replaced by `symbol` at `rate`."""
    var out = List[UInt8]()
    var replacement = symbol.as_bytes()[0]
    for byte in text.as_bytes():
        out.append(replacement if random_float64() < rate else byte)
    return String(unsafe_from_utf8=out)


def test_seeded_bands_fold_symbols_past_acgt() raises:
    """Long pairs with `N` among the bases take seeds, and keep the exact distance.

    A seed holding an `N` goes uncounted, and the second sequence's `N` read as a base, which only adds
    matches: scattered `N`, a run of them in both sequences, and `N` on either side only, against bases,
    so tiles sweep on two planes, on the rows' mask, or on the third plane in full.
    Every distance is the global wavefront's at unit costs over `ACGTN`, which shares no code with the
    band, and every alignment is rescored independently.
    """
    seed(23)
    var unit = Scoring.edit_distance("ACGTN")
    var shortest = max(20000, SEED_COLUMNS + 2000)
    for rate in [0.05, 0.15]:
        var bases = random_sequence(shortest, shortest + 4000, DNA_ALPHABET)
        var bytes = List[UInt8](bases.as_bytes())
        var run = len(bytes) * 2 // 5
        var first = sprinkle(
            String(unsafe_from_utf8=slice_bytes(bytes, 0, run))
            + String("N") * 300
            + String(unsafe_from_utf8=slice_bytes(bytes, run + 300, len(bytes))),
            0.002,
            "N",
        )
        var second = mutate(first, rate)
        var pairs: List[Tuple[String, String]] = [
            (first, second),
            (second, first),
            (first, sprinkle(second, 0.002, "N")),
            (first, mutate(bases, rate)),
            (mutate(bases, rate), first),
        ]
        for pair in pairs:
            var expected = -Int(score[GLOBAL](pair[0], pair[1], unit))
            assert_equal(edit_distance(pair[0], pair[1]), expected)
            var aligned = edit_alignment(pair[0], pair[1])
            assert_equal(Int(aligned.score), expected)
            assert_equal(aligned.first_gapped.replace("-", ""), pair[0])
            assert_equal(aligned.second_gapped.replace("-", ""), pair[1])
            assert_equal(rescore(aligned.first_gapped, aligned.second_gapped, unit), -expected)


# endregion Properties

# region Refusals


def test_refuses_what_it_cannot_do() raises:
    """Every contradiction is refused rather than quietly answered."""
    var dna = Scoring.dna()
    with assert_raises(contains="outside the alphabet"):
        _ = score[GLOBAL]("ACGT", "ACGN", dna)
    with assert_raises(contains="outside the alphabet"):
        _ = score[GLOBAL]("ACGT", "acgt", dna)
    with assert_raises(contains="cannot be served"):
        _ = Scoring.uniform(5, -4, -1, -20)
    with assert_raises(contains="cannot be served"):
        _ = Scoring.uniform(5, -4, -1, 1)
    with assert_raises(contains="cannot be served"):
        _ = Scoring.uniform(500, -4)
    with assert_raises(contains="cannot be served"):
        _ = Scoring.tabulated("AC", List[Int8](length=9, fill=0))
    with assert_raises(contains="do not"):
        _ = scores[GLOBAL](["AC", "CA"], ["AC"], dna)
    with assert_raises(contains="ASCII"):
        _ = levenshtein_alignment("naïve", "naive")


def test_colouring_keeps_the_rows_it_paints() raises:
    """Painting wraps every column in escapes and changes nothing else; ragged rows are refused."""
    var painted = colorize("ACGTA", "AC-TT")
    for row in [painted[0], painted[1]]:
        var plain = row
        for escape in ["\x1b[32m", "\x1b[31m", "\x1b[37m", "\x1b[0m"]:
            plain = plain.replace(escape, "")
        assert_true(plain == "ACGTA" or plain == "AC-TT")
    assert_true(painted[0].startswith("\x1b[32mA"), "a match is not painted green")
    assert_true("\x1b[37m-" in painted[1], "a gap is not painted white")
    assert_true("\x1b[31mT" in painted[1], "a mismatch is not painted red")
    with assert_raises(contains="do not"):
        _ = colorize("ACGTA", "AC")


def mutated(text: String, rate: Float64, longest_gap: Int) -> String:
    """`text` with substitutions, insertions and deletions at `rate` in all, each gap up to `longest_gap` long."""
    var letters = text.as_bytes()
    var bases = DNA_ALPHABET.as_bytes()
    var out = List[UInt8]()
    var index = 0
    while index < len(letters):
        var draw = random_float64()
        if draw < rate / 3:
            out.append(bases[Int(random_ui64(0, 3))])
            index += 1
        elif draw < 2 * rate / 3:
            index += Int(random_ui64(1, UInt64(longest_gap)))
        elif draw < rate:
            for _ in range(Int(random_ui64(1, UInt64(longest_gap)))):
                out.append(bases[Int(random_ui64(0, 3))])
        else:
            out.append(letters[index])
            index += 1
    return String(unsafe_from_utf8=out)


def dna_codes(text: String) -> List[UInt8]:
    """Each base's position in `DNA_ALPHABET`."""
    var codes = List[UInt8]()
    var bases = DNA_ALPHABET.as_bytes()
    for byte in text.as_bytes():
        for code in range(4):
            if bases[code] == byte:
                codes.append(UInt8(code))
    return codes^


def test_wavefront_matches_the_full_sweep() raises:
    """The two-ended wavefront's global score and alignment are the full sweep's optimum.

    Pairs of every shape the meeting must handle: close and far, long gaps either way, one side much
    longer than the other, and a single letter; under WFA's costs and a dear opening. The scores are
    the full Gotoh sweep's, which shares no code with the wavefront.
    """
    seed(23)
    var regimes = List[Scoring]()
    regimes.append(Scoring.uniform(0, -4, -8, -2))
    regimes.append(Scoring.uniform(5, -4, -20, -1))
    regimes.append(Scoring.uniform(1, -3, -6, -1))
    var host = Placement.on_cpu(1)
    for scoring in regimes:
        var gaps = scoring.gaps
        for trial in range(120):
            var first = random_sequence(1, 400, DNA_ALPHABET)
            var rate = [0.0, 0.02, 0.1, 0.3][trial % 4]
            var second = mutated(first, rate, [1, 3, 40][trial % 3])
            if trial % 11 == 0:
                second = random_sequence(1, 3, DNA_ALPHABET)
            if second.byte_length() == 0:
                second = "A"
            var expected = vector_score[GLOBAL](
                dna_codes(first), dna_codes(second), Int(scoring.substitutions[0]), Int(scoring.substitutions[1]), gaps
            )
            var produced = align[GLOBAL](first, second, scoring, host)
            assert_equal(produced.score, expected)
            assert_well_formed[GLOBAL](first, second, produced, scoring)
            assert_equal(score[GLOBAL](first, second, scoring, host), expected)


def test_wavefront_splits_a_pair_too_large_to_keep() raises:
    """A pair whose fronts pass the limit is split where an optimal path crosses, recursively, and the
    pieces' alignment is still optimal: limits of no entries at all, a few and some, so splits fall
    between moves and inside gaps of either sequence, with pieces that must end or begin in them."""
    seed(29)
    var regimes = List[Scoring]()
    regimes.append(Scoring.uniform(0, -4, -8, -2))
    regimes.append(Scoring.uniform(5, -4, -20, -1))
    for scoring in regimes:
        var penalties = wavefront_penalties(
            scoring.substitutions, scoring.alphabet_size(), Int(scoring.gaps.open), Int(scoring.gaps.extend)
        ).value()
        for trial in range(60):
            var first = random_sequence(1, 300, DNA_ALPHABET)
            var rate = [0.02, 0.1, 0.3][trial % 3]
            var second = mutated(first, rate, [1, 5, 60][trial % 3])
            if second.byte_length() == 0:
                second = "C"
            var expected = vector_score[GLOBAL](
                dna_codes(first),
                dna_codes(second),
                Int(scoring.substitutions[0]),
                Int(scoring.substitutions[1]),
                scoring.gaps,
            )
            for limit in [0, 64, 4096]:
                var traced = wavefront_align(dna_codes(first), dna_codes(second), penalties, DNA_ALPHABET, limit)
                var produced = AlignmentResult(Int32(traced[0]), traced[1], traced[2])
                assert_equal(produced.score, expected)
                assert_well_formed[GLOBAL](first, second, produced, scoring)


def test_affine_cigar_spells_an_optimal_alignment() raises:
    """`affine_cigar`'s cost is the full sweep's optimum at WFA's costs, and its CIGAR spells an
    alignment of both sequences that costs exactly that, gap runs and all."""
    seed(31)
    for costs in [(4, 6, 2), (1, 0, 1), (3, 10, 1)]:
        var x = costs[0]
        var o = costs[1]
        var e = costs[2]
        var scoring = Scoring.uniform(0, -x, -(o + e), -e)
        for trial in range(40):
            var first = random_sequence(0, 300, DNA_ALPHABET)
            var second = mutated(first, [0.0, 0.05, 0.2][trial % 3], [1, 4, 30][trial % 3])
            if trial % 13 == 0:
                second = String()
            var found = affine_cigar(first, second, x, o, e)
            var expected = -Int(vector_score[GLOBAL](dna_codes(first), dna_codes(second), 0, -x, scoring.gaps))
            assert_equal(found.cost, expected)
            var rows = rows_from_cigar(first, second, found.cigar)
            assert_equal(rescore(rows[0], rows[1], scoring), -expected)
    # Three substitutions, 12, undercut the two single gaps the edit distance takes, 16.
    var known = affine_cigar("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 4, 6, 2)
    assert_equal(known.cost, 12)
    assert_equal(known.cigar, "4=3X6=")
    assert_equal(affine_cigar("acgu", "acgu", 4, 6, 2).cigar, "4=")
    with assert_raises(contains="must cost"):
        _ = affine_cigar("A", "C", 0, 6, 2)

    # A batch over threads gives every pair what the single call gives it.
    var firsts = List[String]()
    var seconds = List[String]()
    for trial in range(30):
        var first = random_sequence(0, 2000, DNA_ALPHABET)
        firsts.append(first)
        seconds.append(mutated(first, [0.01, 0.1, 0.3][trial % 3], [1, 8, 50][trial % 3]))
    var batch = affine_cigars(firsts, seconds, 4, 6, 2, threads=4)
    for index in range(len(firsts)):
        var single = affine_cigar(firsts[index], seconds[index], 4, 6, 2)
        assert_equal(batch[index].cost, single.cost)
        assert_equal(batch[index].cigar, single.cigar)


def test_affine_distance_and_its_cap() raises:
    """`affine_distance` is `affine_cigar`'s cost, and a cap of `max_cost` returns both up to the
    optimum and neither a unit below it: pairs close and far, empty sides, a cap below zero."""
    seed(37)
    for costs in [(4, 6, 2), (1, 0, 1), (3, 10, 1)]:
        var x = costs[0]
        var o = costs[1]
        var e = costs[2]
        for trial in range(36):
            var first = random_sequence(0, 400, DNA_ALPHABET)
            var second = mutated(first, [0.0, 0.03, 0.15, 0.4][trial % 4], [1, 5, 40][trial % 3])
            if trial % 9 == 0:
                second = String()
            var cost = affine_cigar(first, second, x, o, e).cost
            assert_equal(affine_distance(first, second, x, o, e), cost)
            assert_equal(affine_distance(first, second, x, o, e, max_cost=cost).value(), cost)
            assert_equal(affine_distance(first, second, x, o, e, max_cost=cost + 1000).value(), cost)
            var within = affine_cigar(first, second, x, o, e, max_cost=cost)
            assert_true(Bool(within), "a cap at the optimum refused it")
            assert_equal(within.value().cost, cost)
            assert_equal(within.value().cigar, affine_cigar(first, second, x, o, e).cigar)
            if cost > 0:
                assert_false(Bool(affine_distance(first, second, x, o, e, max_cost=cost - 1)), "under the cap")
                assert_false(Bool(affine_cigar(first, second, x, o, e, max_cost=cost - 1)), "under the cap")
    assert_false(Bool(affine_distance("", "", 4, 6, 2, max_cost=-1)), "a cap below zero")
    assert_equal(affine_distance("", "", 4, 6, 2, max_cost=0).value(), 0)


# endregion Refusals

# region Device


def test_device_matches_host() raises:
    """Every device path returns the host's score, with an alignment that earns it.

    The host traces a global alignment under a uniform table through the wavefront's fronts, which
    may resolve a tie between optimal paths differently from the device's Gotoh walk, so the strings
    need not match. Covers the banded batch, the tiled sweep for a pair too tall for one block, the
    linear-space device traceback, and empty sides. Skipped where no accelerator answers.
    """
    if not gpu_available():
        print("    skipped: no accelerator serves a real alignment here")
        return
    seed(11)
    var device = Placement.on_gpu(0, 4)
    var host = Placement.on_cpu(4)
    for scoring in scoring_regimes():
        var firsts = List[String]()
        var seconds = List[String]()
        for _ in range(REPETITIONS):
            firsts.append(random_sequence(5, 60, DNA_ALPHABET))
            seconds.append(random_sequence(5, 60, DNA_ALPHABET))
        comptime for mode in [GLOBAL, LOCAL]:
            var on_device = alignments[mode](firsts, seconds, scoring, device)
            var device_scores = scores[mode](firsts, seconds, scoring, device)
            for index in range(len(firsts)):
                var expected = align[mode](firsts[index], seconds[index], scoring, host)
                assert_equal(on_device[index].score, expected.score)
                assert_well_formed[mode](firsts[index], seconds[index], on_device[index], scoring)
                assert_well_formed[mode](firsts[index], seconds[index], expected, scoring)
                assert_equal(device_scores[index], expected.score)

    var scoring = expensive_gap()
    var tall = "A" * 70_000
    var short = "ACGT" * 4
    comptime for mode in [GLOBAL, LOCAL]:
        assert_equal(score[mode](tall, short, scoring, device), score[mode](tall, short, scoring, host))
        assert_equal(score[mode](tall, "", scoring, device), score[mode](tall, "", scoring, host))
        assert_equal(score[mode]("", short, scoring, device), score[mode]("", short, scoring, host))
        var long_first = random_sequence(2000, 2000, DNA_ALPHABET)
        var long_second = random_sequence(2000, 2000, DNA_ALPHABET)
        var linear = align[mode](long_first, long_second, scoring, device, 0)
        assert_well_formed[mode](long_first, long_second, linear, scoring)
        assert_equal(linear.score, score[mode](long_first, long_second, scoring, host))


# endregion Device


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
