# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
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
    Aligner,
    Alignment,
    AlignmentError,
    Anchor,
    Band,
    Costs,
    DEFAULT_MAX_MEMORY,
    Mode,
    Placement,
    Scoring,
    Ties,
    align,
    alignments,
    distance,
    distances,
    local_scores,
    score,
    scores,
    search,
)
from dinara_align.alignment import (
    AffineGapCosts,
    AlignmentMode,
    GapRun,
    GappedAlignment,
    SweepHalf,
    colorize,
    serial_align,
    vector_sweep_bands,
)
from dinara_align.scoring import DNA_ALPHABET, tabulated_end
from dinara_align.vector_score import reach_back
from dinara_align.cigar import cigar_runs
from dinara_align.cigar import reversed_text as reversed_bytes
from dinara_align.edit_distance import edit_distance as bit_parallel_distance
from dinara_align.scored import ANYWHERE, best_end, end_of
from dinara_align.seeds import SEED_COLUMNS
from dinara_align.bit_parallel import Profile
from dinara_align.diagonal import DiagonalFronts, diagonal_transition, trace_diagonals
from dinara_align.traceback import EditPath, cigar_string
from dinara_align.gap_affine import (
    ALIGNED,
    EndsFree,
    extension_of,
    extension_penalties,
    free_ends_alignment,
    traced_extension,
    FIRST_GAP,
    FREE_START,
    SECOND_GAP,
    Penalties,
    Wavefront,
    affine2p_penalties,
    affine_penalties,
    solve,
    trace,
    wavefront_align,
    wavefront_penalties,
)
from dinara_align.substitutions import SubstitutionLookup

comptime GLOBAL = Mode.GLOBAL
comptime LOCAL = Mode.LOCAL
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
    return Scoring.uniform(5, -4, -19, -1)


def scoring_regimes() raises -> List[Scoring]:
    """One representative per regime: an expensive gap, a cheap one, unit costs, a free extension, and the DNA default.
    """
    var regimes = List[Scoring]()
    regimes.append(Scoring.uniform(5, -4, -19, -1))
    regimes.append(Scoring.uniform(2, -1, -1, -1))
    regimes.append(Scoring.uniform(0, -1, 0, -1))
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


def brute_optimum(mode: Mode, first: String, second: String, scoring: Scoring) -> Int:
    """Global: the best enumerated alignment. Local: the best over every pair of substrings, or zero."""
    var top = List[UInt8](first.as_bytes())
    var bottom = List[UInt8](second.as_bytes())
    if mode == GLOBAL:
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


def gapped_rows(found: Alignment, first: String, second: String) -> GappedAlignment:
    """An alignment's two gapped rows and its score, as Gotoh's rows were."""
    var rows = found.gapped(first, second)
    return GappedAlignment(Int32(found.score), rows[0], rows[1])


def assert_well_formed(mode: Mode, first: String, second: String, found: Alignment, scoring: Scoring) raises:
    """Both rows have one length, each rebuilds its input or a piece of it, and they earn their score;
    the spans and the CIGAR read the rows back, and the cost is minus the score."""
    assert_well_formed(mode, first, second, gapped_rows(found, first, second), scoring)
    assert_equal(found.cost, -found.score)
    if mode == GLOBAL:
        assert_equal(found.reference_start, 0)
        assert_equal(found.reference_end, first.byte_length())


def assert_well_formed(mode: Mode, first: String, second: String, produced: GappedAlignment, scoring: Scoring) raises:
    """Both rows have one length, each rebuilds its input or a piece of it, and they earn their score."""
    assert_equal(produced.first_gapped.byte_length(), produced.second_gapped.byte_length())
    var core_first = produced.first_gapped.replace("-", "")
    var core_second = produced.second_gapped.replace("-", "")
    if mode == GLOBAL:
        assert_equal(core_first, first)
        assert_equal(core_second, second)
    else:
        assert_true(core_first in first, "a local row is not a substring of its input")
        assert_true(core_second in second, "a local row is not a substring of its input")
    assert_equal(rescore(produced.first_gapped, produced.second_gapped, scoring), Int(produced.score))
    # The rows read back as a CIGAR spell the same rows, over the letters they align.
    var spelled = rows_from_cigar(core_first, core_second, produced.cigar())
    assert_equal(spelled[0], produced.first_gapped)
    assert_equal(spelled[1], produced.second_gapped)


def gpu_available() raises -> Bool:
    """Whether an accelerator serves a real alignment here, which a successful import does not prove."""
    var scoring = Scoring.dna()
    try:
        _ = score("AC", "CA", scoring, GLOBAL, placement=Placement.on_gpu(0, 1))
        return True
    except:
        return False


def unit_distance(first: String, second: String) raises -> Int:
    """The edit distance by Gotoh's recurrence at unit costs over the pair's own letters, which shares no
    code with the bit-parallel sweep or the wavefront."""
    var seen = List[Bool](length=256, fill=False)
    var letters = List[UInt8]()
    for text in [first, second]:
        for byte in text.as_bytes():
            if not seen[Int(byte)]:
                seen[Int(byte)] = True
                letters.append(byte)
    if len(letters) == 0:
        letters.append(UInt8(ord("A")))
    return -Int(score(first, second, Scoring.edit_distance(String(unsafe_from_utf8=letters^)), GLOBAL))


def edit_rows(first: String, second: String, ties: Ties = Ties.LEFT) raises -> GappedAlignment:
    """`align` at unit costs as two gapped rows, its score the distance."""
    var aligned = align(first, second, ties=ties)
    var rows = aligned.gapped(first, second)
    return GappedAlignment(Int32(aligned.cost), rows[0], rows[1])


# endregion Oracles

# region Known Answers


def test_dna_default_is_minimap2() raises:
    """Match 2, mismatch -4, and minimap2's `-O4 -E2` gap, `-(4 + 2k)`: its first gapped base six."""
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
    var same = align("ACGTACGT", "ACGTACGT", dna)
    assert_equal(same.gapped("ACGTACGT", "ACGTACGT")[0], "ACGTACGT")
    assert_equal(same.score, 16)

    # One substitution costs -4, which beats opening two gaps at -6 each.
    var substituted = align("ACGTACGT", "ACGTTCGT", dna)
    assert_equal(substituted.gapped("ACGTACGT", "ACGTTCGT")[0], "ACGTACGT")
    assert_equal(substituted.gapped("ACGTACGT", "ACGTTCGT")[1], "ACGTTCGT")
    assert_equal(substituted.score, 7 * 2 - 4)

    # A three-base deletion is one gap, scoring -(4 + 2 * 3).
    var deleted = align("ACGTTGCAGGGCATGACGT", "ACGTTGCACATGACGT", dna)
    assert_equal(deleted.gapped("ACGTTGCAGGGCATGACGT", "ACGTTGCACATGACGT")[0], "ACGTTGCAGGGCATGACGT")
    assert_equal(deleted.gapped("ACGTTGCAGGGCATGACGT", "ACGTTGCACATGACGT")[1], "ACGTTGCA---CATGACGT")
    assert_equal(deleted.cigar, "8=3D8=")
    assert_equal(substituted.cigar, "4=1X3=")
    assert_equal(align("ACGTACGT", "ACGTTCGT", dna, eqx=False).cigar, "8M")
    assert_equal(deleted.score, 16 * 2 - 6 - 2 * 2)
    assert_equal(score("ACGTTGCAGGGCATGACGT", "ACGTTGCACATGACGT", dna), 22)

    # Against nothing, the whole sequence is one gap run.
    assert_equal(score("AAAA", "", dna), -6 - 3 * 2)
    assert_equal(score("", "", dna), 0)


def test_hand_computed_local() raises:
    """The shared core of two otherwise unrelated sequences, trimmed at both ends."""
    var aligned = align("GGGGACGTACGTGGGG", "CCCCACGTACGTCCCC", Scoring.dna(), LOCAL)
    assert_equal(aligned.cigar, "8=")
    assert_equal(aligned.reference_start, 4)
    assert_equal(aligned.query_end, 12)
    assert_equal(aligned.score, 16)
    assert_equal(score("GGGGACGTACGTGGGG", "CCCCACGTACGTCCCC", Scoring.dna(), LOCAL), 16)


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
        var aligned = align(example[0], example[1])
        assert_equal(aligned.cost, example[2])
        var rows = aligned.gapped(example[0], example[1])
        assert_equal(rows[0].replace("-", ""), example[0])
        assert_equal(rows[1].replace("-", ""), example[1])


# endregion Known Answers

# region Exhaustive Oracle


def check_against_enumeration(mode: Mode) raises:
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
            var expected = brute_optimum(mode, first, second, scoring)
            assert_equal(Int(score(first, second, scoring, mode, placement=Placement.on_cpu(1))), expected)
            var produced = align(first, second, scoring, mode, placement=Placement.on_cpu(1))
            assert_equal(Int(produced.score), expected)


def test_global_matches_enumeration() raises:
    """A global alignment's score is every tiny pair's optimum by enumeration."""
    check_against_enumeration(GLOBAL)


def test_local_matches_enumeration() raises:
    """A local alignment's score is every tiny pair's optimum by enumeration."""
    check_against_enumeration(LOCAL)


# endregion Exhaustive Oracle

# region Properties


def check_well_formed(mode: Mode) raises:
    """Every path is well formed, realizes its own score, and agrees with the score-only kernel."""
    seed(1)
    for scoring in scoring_regimes():
        for _ in range(REPETITIONS):
            var first = random_sequence(5, 25, DNA_ALPHABET)
            var second = random_sequence(5, 25, DNA_ALPHABET)
            var produced = align(first, second, scoring, mode)
            assert_well_formed(mode, first, second, produced, scoring)
            assert_equal(score(first, second, scoring, mode), produced.score)


def test_global_output_is_well_formed() raises:
    check_well_formed(GLOBAL)


def test_local_output_is_well_formed() raises:
    check_well_formed(LOCAL)


def check_linear_matches_stored(mode: Mode) raises:
    """Both traceback strategies reach the same score, each with a path that earns it.

    Ties may land differently under a divide-and-conquer join, so the strings need not match.
    """
    seed(2)
    for scoring in scoring_regimes():
        for _ in range(3):
            var first = random_sequence(140, 260, DNA_ALPHABET)
            var second = random_sequence(140, 260, DNA_ALPHABET)
            var stored = align(first, second, scoring, mode, max_memory=10**12)
            var linear = align(first, second, scoring, mode, max_memory=0)
            assert_equal(stored.score, linear.score)
            assert_equal(stored.score, score(first, second, scoring, mode))
            assert_well_formed(mode, first, second, stored, scoring)
            assert_well_formed(mode, first, second, linear, scoring)


def test_global_linear_matches_stored() raises:
    """A global alignment traced in linear space scores what one traced from a stored matrix does."""
    check_linear_matches_stored(GLOBAL)


def test_local_linear_matches_stored() raises:
    """A local alignment traced in linear space scores what one traced from a stored matrix does."""
    check_linear_matches_stored(LOCAL)


def test_linear_space_carries_a_long_pair() raises:
    """Three thousand by three thousand with no stored matrix at all."""
    seed(3)
    var scoring = expensive_gap()
    var first = random_sequence(3000, 3000, DNA_ALPHABET)
    var second = random_sequence(3000, 3000, DNA_ALPHABET)
    assert_well_formed(GLOBAL, first, second, align(first, second, scoring, GLOBAL, max_memory=0), scoring)
    assert_well_formed(LOCAL, first, second, align(first, second, scoring, LOCAL, max_memory=0), scoring)


def test_symmetry() raises:
    """Swapping the arguments does not change the score."""
    seed(4)
    var scoring = expensive_gap()
    for _ in range(REPETITIONS):
        var first = random_sequence(5, 25, DNA_ALPHABET)
        var second = random_sequence(5, 25, DNA_ALPHABET)
        assert_equal(score(first, second, scoring, GLOBAL), score(second, first, scoring, GLOBAL))
        assert_equal(score(first, second, scoring, LOCAL), score(second, first, scoring, LOCAL))


def test_levenshtein_is_the_unit_cost_limit() raises:
    """At unit costs the global recurrence is the negated edit distance."""
    seed(5)
    var unit = Scoring.edit_distance()
    for _ in range(REPETITIONS):
        var first = random_sequence(3, 15, DNA_ALPHABET)
        var second = random_sequence(3, 15, DNA_ALPHABET)
        assert_equal(-Int(score(first, second, unit, GLOBAL)), distance(first, second))


def test_local_never_scores_below_global() raises:
    """A global path is also a local candidate, and the empty window is always available."""
    seed(6)
    var scoring = expensive_gap()
    for _ in range(REPETITIONS):
        var first = random_sequence(5, 25, DNA_ALPHABET)
        var second = random_sequence(5, 25, DNA_ALPHABET)
        var local = score(first, second, scoring, LOCAL)
        assert_true(local >= 0)
        assert_true(local >= score(first, second, scoring, GLOBAL))


def test_optimum_falls_as_gaps_get_harsher() raises:
    """A harsher opening lowers every feasible alignment, so the maximum cannot rise."""
    seed(7)
    for _ in range(REPETITIONS):
        var first = random_sequence(5, 25, DNA_ALPHABET)
        var second = random_sequence(5, 25, DNA_ALPHABET)
        var previous = Int.MAX
        for opening in [-2, -5, -10, -20, -40]:
            var current = score(first, second, Scoring.uniform(5, -4, opening, -1), GLOBAL)
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
        var reference = score(first, head + "N" + tail, free_extension, GLOBAL)
        for width in range(2, 6):
            assert_equal(score(first, head + "N" * width + tail, free_extension, GLOBAL), reference)


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
        assert_equal(score(first, second, tabulated, GLOBAL), score(first, second, uniform, GLOBAL))
        assert_equal(score(first, second, tabulated, LOCAL), score(first, second, uniform, LOCAL))


def test_batch_matches_single_pairs() raises:
    """The batched entry points reproduce the single-pair ones, and an empty batch answers empty."""
    seed(10)
    var scoring = expensive_gap()
    var firsts = List[String]()
    var seconds = List[String]()
    for _ in range(24):
        firsts.append(random_sequence(5, 40, DNA_ALPHABET))
        seconds.append(random_sequence(5, 40, DNA_ALPHABET))
    var batch_scores = scores(firsts, seconds, scoring, GLOBAL)
    var batch_alignments = alignments(firsts, seconds, scoring, LOCAL)
    for index in range(len(firsts)):
        assert_equal(batch_scores[index], score(firsts[index], seconds[index], scoring, GLOBAL))
        var single = align(firsts[index], seconds[index], scoring, LOCAL)
        assert_equal(batch_alignments[index].cigar, single.cigar)
        assert_equal(batch_alignments[index].reference_start, single.reference_start)
        assert_equal(batch_alignments[index].score, single.score)
    assert_equal(len(scores(List[String](), List[String](), scoring, GLOBAL)), 0)


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
            assert_equal(distance(first, second), unit_distance(first, second))
    # Similar pairs too, where long runs of matches carry through whole words.
    for _ in range(REPETITIONS):
        var first = random_sequence(500, 900, DNA_ALPHABET)
        var second = first
        var bytes = List[UInt8](second.as_bytes())
        for _ in range(10):
            bytes[Int(random_ui64(0, UInt64(len(bytes) - 1)))] = UInt8(ord("A"))
        second = String(unsafe_from_utf8=bytes)
        assert_equal(distance(first, second), unit_distance(first, second))
    with assert_raises(contains="symbols"):
        _ = bit_parallel_distance("ACGT", "ACGTNRYKM")
    # Past what the sweep takes, the wavefront serves the pair at the same costs.
    assert_equal(distance("ACGT", "ACGTNRYKM"), 5)


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
                    var expected = -Int(score(pair[0], pair[1], unit, GLOBAL))
                    assert_equal(distance(pair[0], pair[1]), expected)
                    var aligned = edit_rows(pair[0], pair[1])
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
        assert_equal(distance(first, second), -Int(score(first, second, unit, GLOBAL)))
    var first = random_sequence(9000, 9000, DNA_ALPHABET)
    var bytes = List[UInt8](first.as_bytes())
    for _ in range(900):
        bytes[Int(random_ui64(0, UInt64(len(bytes) - 1)))] = UInt8(ord("G"))
    var second = String(unsafe_from_utf8=bytes)
    assert_equal(distance(first, second), unit_distance(first, second))


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
                assert_equal(distance(pair[0], pair[1]), unit_distance(pair[0], pair[1]))
    var unrelated_first = random_sequence(3000, 3000, DNA_ALPHABET)
    var unrelated_second = random_sequence(2800, 2800, DNA_ALPHABET)
    assert_equal(
        distance(unrelated_first, unrelated_second),
        unit_distance(unrelated_first, unrelated_second),
    )


def test_unit_cost_alignment_is_optimal() raises:
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
                var aligned = edit_rows(pair[0], pair[1])
                assert_equal(Int(aligned.score), unit_distance(pair[0], pair[1]))
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


def test_infix_and_prefix_match_the_dynamic_program() raises:
    """A pattern found inside a text, or at its start, at the distance the textbook dynamic program
    gives, and aligned to the part of the text reported, which the alignment rebuilds.

    Patterns across word boundaries, from empty to longer than the text, planted mutated copies of a
    piece of the text and unrelated ones, with an `N` among the symbols too; a long text under a long
    pattern moves the prefix search's band far down.
    """
    seed(23)
    var unit = Scoring.edit_distance("ACGTN")
    for pattern_length in [0, 1, 63, 64, 65, 129, 300, 1000]:
        for text_length in [0, 1, 200, 700, 2500]:
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
                    var mode = Mode.PREFIX if prefix else Mode.INFIX
                    assert_equal(distance(text, pattern, Costs.edit(), mode), expected)
                    var found = align(text, pattern, Costs.edit(), mode)
                    assert_equal(found.cost, expected)
                    if prefix:
                        assert_equal(found.reference_start, 0)
                    assert_equal(found.query_start, 0)
                    assert_equal(found.query_end, pattern.byte_length())
                    var rows = found.gapped(text, pattern)
                    var part = String(
                        StringSlice(unsafe_from_utf8=text.as_bytes()[found.reference_start : found.reference_end])
                    )
                    assert_equal(rows[0].replace("-", ""), part)
                    assert_equal(rows[1].replace("-", ""), pattern)
                    assert_equal(rescore(rows[0], rows[1], unit), -expected)


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
        var found = distances(firsts, seconds, threads=threads)
        var aligned = alignments(firsts, seconds, threads=threads)
        assert_equal(len(found), len(firsts))
        assert_equal(len(aligned), len(firsts))
        for index in range(len(firsts)):
            var expected = distance(firsts[index], seconds[index])
            assert_equal(found[index], expected)
            assert_equal(aligned[index].cost, expected)
            assert_equal(aligned[index].cigar, align(firsts[index], seconds[index]).cigar)
            var rows = aligned[index].gapped(firsts[index], seconds[index])
            assert_equal(rows[0].replace("-", ""), firsts[index])
            assert_equal(rows[1].replace("-", ""), seconds[index])
    # A band no alignment of one pair fits fails the batch, as a serial loop would; so do sides of
    # different lengths.
    var bad_firsts: List[String] = ["ACGT", "ACGT", "ACGTACGTACGT"]
    var bad_seconds: List[String] = ["ACGA", "ACG", "ACGT"]
    with assert_raises(contains="band"):
        _ = distances(bad_firsts, bad_seconds, band=Band.around(2), threads=8)
    with assert_raises(contains="band"):
        _ = alignments(bad_firsts, bad_seconds, band=Band.around(2), threads=8)
    var short: List[String] = ["ACGT"]
    with assert_raises():
        _ = distances(bad_firsts, short)


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
            var expected = -Int(score(pair[0], pair[1], unit, GLOBAL))
            assert_equal(distance(pair[0], pair[1]), expected)
            var aligned = edit_rows(pair[0], pair[1])
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


def test_unit_cost_cigar_spells_the_alignment() raises:
    """At unit costs the CIGAR spells an optimal alignment, run by run.

    Pairs settled by one diagonal front, by two, and by a band with seeds, with symbols past `ACGT`,
    and with empty sides; the CIGAR rebuilds both sequences, every `=` joins equal bases and every `X`
    differing ones, its substitutions and gaps number the distance, and with `M` for both the same
    runs merge.
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
        var spelled = align(pair[0], pair[1])
        assert_equal(spelled.cost, distance(pair[0], pair[1]))
        var rows = rows_from_cigar(pair[0], pair[1], spelled.cigar)
        assert_equal(rows[2], spelled.cost)
        var gapped = spelled.gapped(pair[0], pair[1])
        assert_equal(rows[0], gapped[0])
        assert_equal(rows[1], gapped[1])
        var plain = align(pair[0], pair[1], eqx=False)
        assert_equal(plain.cost, spelled.cost)
        var plain_rows = rows_from_cigar(pair[0], pair[1], plain.cigar)
        assert_equal(plain_rows[0], rows[0])
        assert_equal(plain_rows[1], rows[1])
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
            var expected = -Int(score(pair[0], pair[1], unit, GLOBAL))
            assert_equal(distance(pair[0], pair[1]), expected)
            var aligned = edit_rows(pair[0], pair[1])
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
        _ = score("ACGT", "ACGN", dna, GLOBAL)
    with assert_raises(contains="outside the alphabet"):
        _ = score("ACGT", "acgt", dna, GLOBAL)
    with assert_raises(contains="cannot be served"):
        _ = Scoring.uniform(5, -4, 1, -2)
    with assert_raises(contains="cannot be served"):
        _ = Scoring.uniform(5, -4, -1, 1)
    with assert_raises(contains="cannot be served"):
        _ = Scoring.uniform(500, -4)
    with assert_raises(contains="cannot be served"):
        _ = Scoring.tabulated("AC", List[Int8](length=9, fill=0))
    with assert_raises(contains="do not"):
        _ = scores(["AC", "CA"], ["AC"], dna, GLOBAL)
    var rewarded = Mode.INFIX.with_match_score(2)
    with assert_raises(contains="holds what a match earns"):
        _ = score("ACGT", "ACG", dna, rewarded)
    with assert_raises(contains="cannot be served"):
        _ = Costs.affine(0, 6, 2)
    with assert_raises(contains="cannot be served"):
        _ = Costs.two_piece(4, 6, 2, 24, 0)
    with assert_raises(contains="cannot be served"):
        _ = Mode.extension(-1)
    with assert_raises(contains="rejected"):
        _ = Mode.ends_free(reference_start=-1)
    with assert_raises(contains="score"):
        _ = distance("ACGT", "ACG", Costs.edit(), Mode.local(1))
    with assert_raises(contains="earns"):
        _ = align("ACGT", "ACG", Costs.edit(), Mode.LOCAL)
    with assert_raises(contains="earns"):
        _ = Mode.local(0)
    with assert_raises(contains="earns"):
        _ = Mode.overlap(0)
    with assert_raises(contains="score"):
        _ = distance("ACGT", "ACG", Costs.edit(), Mode.overlap(1))
    with assert_raises(contains="score"):
        _ = distance("ACGT", "ACG", Costs.edit(), Mode.INFIX.with_match_score(1))
    with assert_raises(contains="free ends alone"):
        _ = Mode.local(1).with_match_score(2)
    var overlap = Mode.overlap(1)
    with assert_raises(contains="table"):
        _ = score("ACGT", "ACG", dna, overlap)
    with assert_raises(contains="score"):
        _ = distance("ACGT", "ACG", Costs.edit(), Mode.extension(1))
    with assert_raises(contains="cap"):
        _ = align("ACGT", "ACG", Costs.edit(), Mode.extension(1), max_cost=3)


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
    regimes.append(Scoring.uniform(5, -4, -19, -1))
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
            var expected = serial_align[AlignmentMode.GLOBAL](
                dna_codes(first), dna_codes(second), scoring.substitutions, 4, gaps, DNA_ALPHABET
            ).score
            var produced = align(first, second, scoring, GLOBAL, placement=host)
            assert_equal(produced.score, Int(expected))
            assert_well_formed(GLOBAL, first, second, produced, scoring)
            assert_equal(score(first, second, scoring, GLOBAL, placement=host), Int(expected))


def test_wavefront_splits_a_pair_too_large_to_keep() raises:
    """A pair whose fronts pass the limit is split where an optimal path crosses, recursively, and the
    pieces' alignment is still optimal: limits of no entries at all, a few and some, so splits fall
    between moves and inside gaps of either sequence, with pieces that must end or begin in them."""
    seed(29)
    var regimes = List[Scoring]()
    regimes.append(Scoring.uniform(0, -4, -8, -2))
    regimes.append(Scoring.uniform(5, -4, -19, -1))
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
            var expected = serial_align[AlignmentMode.GLOBAL](
                dna_codes(first), dna_codes(second), scoring.substitutions, 4, scoring.gaps, DNA_ALPHABET
            ).score
            for limit in [0, 64, 4096]:
                var traced = wavefront_align(dna_codes(first), dna_codes(second), penalties, DNA_ALPHABET, limit)
                var produced = GappedAlignment(Int32(traced[0]), traced[1], traced[2])
                assert_equal(produced.score, expected)
                assert_well_formed(GLOBAL, first, second, produced, scoring)


def test_affine_cigar_spells_an_optimal_alignment() raises:
    """`align`'s cost under `Costs.affine` is the full sweep's optimum at WFA's costs, and its CIGAR spells an
    alignment of both sequences that costs exactly that, gap runs and all."""
    seed(31)
    for costs in [(4, 6, 2), (1, 0, 1), (3, 10, 1)]:
        var x = costs[0]
        var o = costs[1]
        var e = costs[2]
        var scoring = Scoring.uniform(0, -x, -o, -e)
        for trial in range(40):
            var first = random_sequence(0, 300, DNA_ALPHABET)
            var second = mutated(first, [0.0, 0.05, 0.2][trial % 3], [1, 4, 30][trial % 3])
            if trial % 13 == 0:
                second = String()
            var found = align(first, second, Costs.affine(x, o, e))
            var expected = -Int(
                serial_align[AlignmentMode.GLOBAL](
                    dna_codes(first), dna_codes(second), scoring.substitutions, 4, scoring.gaps, DNA_ALPHABET
                ).score
            )
            assert_equal(found.cost, expected)
            var rows = rows_from_cigar(first, second, found.cigar)
            assert_equal(rescore(rows[0], rows[1], scoring), -expected)
    # Three substitutions, 12, undercut the two single gaps the edit distance takes, 16.
    var known = align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", Costs.affine(4, 6, 2))
    assert_equal(known.cost, 12)
    assert_equal(known.cigar, "4=3X6=")
    assert_equal(align("acgu", "acgu", Costs.affine(4, 6, 2)).cigar, "4=")
    with assert_raises(contains="must cost"):
        _ = align("A", "C", Costs.affine(0, 6, 2))

    # A batch over threads gives every pair what the single call gives it.
    var firsts = List[String]()
    var seconds = List[String]()
    for trial in range(30):
        var first = random_sequence(0, 2000, DNA_ALPHABET)
        firsts.append(first)
        seconds.append(mutated(first, [0.01, 0.1, 0.3][trial % 3], [1, 8, 50][trial % 3]))
    var batch = alignments(firsts, seconds, Costs.affine(4, 6, 2), threads=4)
    for index in range(len(firsts)):
        var single = align(firsts[index], seconds[index], Costs.affine(4, 6, 2))
        assert_equal(batch[index].cost, single.cost)
        assert_equal(batch[index].cigar, single.cigar)


def test_affine_distance_and_its_cap() raises:
    """`distance` under `Costs.affine` is `align`'s cost, and a cap of `max_cost` returns both up to the
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
            var cost = align(first, second, Costs.affine(x, o, e)).cost
            assert_equal(distance(first, second, Costs.affine(x, o, e)), cost)
            assert_equal(distance(first, second, Costs.affine(x, o, e), max_cost=cost).value(), cost)
            assert_equal(distance(first, second, Costs.affine(x, o, e), max_cost=cost + 1000).value(), cost)
            var within = align(first, second, Costs.affine(x, o, e), max_cost=cost)
            assert_true(Bool(within), "a cap at the optimum refused it")
            assert_equal(within.value().cost, cost)
            assert_equal(within.value().cigar, align(first, second, Costs.affine(x, o, e)).cigar)
            if cost > 0:
                assert_false(Bool(distance(first, second, Costs.affine(x, o, e), max_cost=cost - 1)), "under the cap")
                assert_false(Bool(align(first, second, Costs.affine(x, o, e), max_cost=cost - 1)), "under the cap")
    assert_false(Bool(distance("", "", Costs.affine(4, 6, 2), max_cost=-1)), "a cap below zero")
    assert_equal(distance("", "", Costs.affine(4, 6, 2), max_cost=0).value(), 0)


def ends_mode(ends: EndsFree) raises AlignmentError -> Mode:
    """The mode freeing what `ends` frees."""
    return Mode.ends_free(
        reference_start=ends.first_begin,
        reference_end=ends.first_end,
        query_start=ends.second_begin,
        query_end=ends.second_end,
    )


def whole(found: Alignment, first: String, second: String) -> Alignment:
    """`found` with its free letters back in its CIGAR as `D` and `I` runs, spanning both sequences, so
    it prices and walks as a global alignment's would."""
    var leading = String()
    if found.reference_start > 0:
        leading = String(found.reference_start, "D")
    elif found.query_start > 0:
        leading = String(found.query_start, "I")
    var trailing = String()
    if found.reference_end < first.byte_length():
        trailing = String(first.byte_length() - found.reference_end, "D")
    elif found.query_end < second.byte_length():
        trailing = String(second.byte_length() - found.query_end, "I")
    return Alignment(
        found.cost,
        found.score,
        merged_cigar(leading + found.cigar + trailing),
        0,
        first.byte_length(),
        0,
        second.byte_length(),
    )


def merged_cigar(cigar: String) -> String:
    """A CIGAR with neighbouring runs of one letter joined."""
    var letters = List[UInt8]()
    var lengths = List[Int]()
    var length = 0
    for byte in cigar.as_bytes():
        if byte >= UInt8(ord("0")) and byte <= UInt8(ord("9")):
            length = length * 10 + Int(byte - UInt8(ord("0")))
            continue
        if len(letters) > 0 and letters[len(letters) - 1] == byte:
            lengths[len(lengths) - 1] += length
        else:
            letters.append(byte)
            lengths.append(length)
        length = 0
    var out = String()
    for index in range(len(letters)):
        out += String(lengths[index], chr(Int(letters[index])))
    return out


def matches_in(cigar: String) -> Int:
    """The `=` letters of a CIGAR."""
    var total = 0
    var length = 0
    for byte in cigar.as_bytes():
        if byte >= UInt8(ord("0")) and byte <= UInt8(ord("9")):
            length = length * 10 + Int(byte - UInt8(ord("0")))
            continue
        if byte == UInt8(ord("=")):
            total += length
        length = 0
    return total


def ends_free_optimum(
    first: String,
    second: String,
    x: Int,
    o: Int,
    e: Int,
    ends: EndsFree,
    o2: Int = -1,
    e2: Int = 0,
    band: Band = Band(),
) -> Int:
    """The least gap-affine cost of aligning two sequences with the letters `ends` allows at each end
    left unaligned for nothing, by Gotoh's recurrence over the whole matrix: starts free along the
    first row and column up to the leading allowances, the best end along the last row and column
    within the trailing ones. With `o2` not negative, a gap costs the less of its two pieces, each
    with a layer either way of its own. Cells off `band`'s diagonals, `i - j`, are never reached; `1 <<
    40` when no end is. Shares no code with the wavefront."""
    comptime HIGH = 1 << 40
    var a = first.as_bytes()
    var b = second.as_bytes()
    var n = len(a)
    var m = len(b)
    var width = m + 1
    var best = List[Int](length=(n + 1) * width, fill=HIGH)
    var across = List[Int](length=(n + 1) * width, fill=HIGH)
    var down = List[Int](length=(n + 1) * width, fill=HIGH)
    var across2 = List[Int](length=(n + 1) * width, fill=HIGH)
    var down2 = List[Int](length=(n + 1) * width, fill=HIGH)
    for i in range(n + 1):
        for j in range(m + 1):
            var at = i * width + j
            if not band.holds(i - j):
                continue
            if (j == 0 and i <= ends.first_begin) or (i == 0 and j <= ends.second_begin):
                best[at] = 0
                continue
            if i > 0:
                across[at] = min(best[at - width] + o + e, across[at - width] + e)
            if j > 0:
                down[at] = min(best[at - 1] + o + e, down[at - 1] + e)
            var value = min(across[at], down[at])
            if o2 >= 0:
                if i > 0:
                    across2[at] = min(best[at - width] + o2 + e2, across2[at - width] + e2)
                if j > 0:
                    down2[at] = min(best[at - 1] + o2 + e2, down2[at - 1] + e2)
                value = min(value, min(across2[at], down2[at]))
            if i > 0 and j > 0:
                value = min(value, best[at - width - 1] + (0 if a[i - 1] == b[j - 1] else x))
            best[at] = value
    var answer = HIGH
    for j in range(m + 1):
        if m - j <= ends.second_end:
            answer = min(answer, best[n * width + j])
    for i in range(n + 1):
        if n - i <= ends.first_end:
            answer = min(answer, best[i * width + m])
    return answer


def ends_free_price(cigar: String, x: Int, o: Int, e: Int, ends: EndsFree, o2: Int = -1, e2: Int = 0) raises -> Int:
    """What a CIGAR's alignment costs with its leading and trailing gap runs free up to `ends`, each
    run at the cheaper piece with `o2` not negative."""
    var kinds = List[UInt8]()
    var lengths = List[Int]()
    var length = 0
    for byte in cigar.as_bytes():
        if byte >= UInt8(ord("0")) and byte <= UInt8(ord("9")):
            length = length * 10 + Int(byte - UInt8(ord("0")))
            continue
        kinds.append(byte)
        lengths.append(length)
        length = 0
    var total = 0
    for index in range(len(kinds)):
        var kind = kinds[index]
        var run = lengths[index]
        if kind == UInt8(ord("X")):
            total += x * run
        elif kind == UInt8(ord("D")) or kind == UInt8(ord("I")):
            var deletion = kind == UInt8(ord("D"))
            var free = 0
            if index == 0:
                free += ends.first_begin if deletion else ends.second_begin
            if index == len(kinds) - 1:
                free += ends.first_end if deletion else ends.second_end
            var paid = max(run - free, 0)
            if paid > 0:
                total += min(o + e * paid, o2 + e2 * paid) if o2 >= 0 else o + e * paid
    return total


def test_affine_ends_free_matches_the_full_matrix() raises:
    """With ends free, the cost is the whole matrix's under the same allowances, the CIGAR spells an
    alignment of both sequences that costs it, its free runs as `D` and `I`, the cap and the split of
    a pair too large to keep both honour the allowances."""
    seed(41)
    for costs in [(4, 6, 2), (1, 0, 1), (3, 10, 1)]:
        var x = costs[0]
        var o = costs[1]
        var e = costs[2]
        var penalties = affine_penalties(x, o, e)
        for trial in range(60):
            var core = random_sequence(1, 200, DNA_ALPHABET)
            var first = random_sequence(0, 40, DNA_ALPHABET) + core + random_sequence(0, 40, DNA_ALPHABET)
            var second = mutated(core, [0.0, 0.03, 0.15][trial % 3], [1, 5, 20][trial % 3])
            if trial % 5 == 0:
                first, second = second, first
            var sizes = [0, 3, 40, 1000]
            var ends = EndsFree(
                sizes[trial % 4], sizes[(trial // 4) % 4], sizes[(trial // 2) % 4], sizes[(trial // 3) % 4]
            )
            var expected = ends_free_optimum(first, second, x, o, e, ends)
            assert_equal(distance(first, second, Costs.affine(x, o, e), ends_mode(ends)), expected)
            var found = whole(align(first, second, Costs.affine(x, o, e), ends_mode(ends)), first, second)
            assert_equal(found.cost, expected)
            var rows = rows_from_cigar(first, second, found.cigar)
            assert_equal(rows[0].replace("-", ""), first)
            assert_equal(rows[1].replace("-", ""), second)
            assert_equal(ends_free_price(found.cigar, x, o, e, ends), expected)
            var capped = align(first, second, Costs.affine(x, o, e), ends_mode(ends), max_cost=expected)
            assert_equal(capped.value().cost, expected)
            if expected > 0:
                assert_false(
                    Bool(distance(first, second, Costs.affine(x, o, e), ends_mode(ends), max_cost=expected - 1))
                )
            for limit in [0, 64]:
                var split = (
                    free_ends_alignment[1](first, second, penalties, True, Int.MAX, ends, Band(), Ties.LEFT, limit)
                    .value()
                    .copy()
                )
                assert_equal(split.cost, expected)
                var spanned = Alignment(
                    split.cost,
                    -split.cost,
                    split.cigar,
                    split.first_start,
                    split.first_end,
                    split.second_start,
                    split.second_end,
                )
                assert_equal(ends_free_price(whole(spanned, first, second).cigar, x, o, e, ends), expected)
    var placed = whole(
        align("TTTTACGTACGTTTTT", "ACGTACGT", Costs.affine(4, 6, 2), ends_mode(EndsFree(16, 16, 0, 0))),
        "TTTTACGTACGTTTTT",
        "ACGTACGT",
    )
    assert_equal(placed.cost, 0)
    assert_equal(placed.cigar, "4D8=4D")


def stays_inside(cigar: String, ends: EndsFree, band: Band) raises -> Bool:
    """Whether every cell of a CIGAR's path from where it starts to where it stops lies on `band`'s
    diagonals, its letters left unaligned for nothing at either end, as `ends_free_price` counts them,
    before its start and after its stop."""
    var kinds = List[UInt8]()
    var lengths = List[Int]()
    var length = 0
    for byte in cigar.as_bytes():
        if byte >= UInt8(ord("0")) and byte <= UInt8(ord("9")):
            length = length * 10 + Int(byte - UInt8(ord("0")))
            continue
        kinds.append(byte)
        lengths.append(length)
        length = 0
    var cells = List[Int]()
    var diagonal = 0
    var skipped = 0
    var dropped = 0
    for index in range(len(kinds)):
        var deletion = kinds[index] == UInt8(ord("D"))
        var gap = deletion or kinds[index] == UInt8(ord("I"))
        if gap and index == 0:
            skipped = min(lengths[index], ends.first_begin if deletion else ends.second_begin)
        if gap and index == len(kinds) - 1 and len(kinds) > 1:
            dropped = min(lengths[index], ends.first_end if deletion else ends.second_end)
        for _ in range(lengths[index]):
            if kinds[index] == UInt8(ord("D")):
                diagonal += 1
            elif kinds[index] == UInt8(ord("I")):
                diagonal -= 1
            cells.append(diagonal)
    # The start cell, after the letters skipped, and every cell up to the stop.
    var path = [0] if skipped == 0 else List[Int]()
    for index in range(max(skipped - 1, 0), len(cells) - dropped):
        path.append(cells[index])
    for cell in path:
        if not band.holds(cell):
            return False
    return True


def extension_optimum(
    first: String, second: String, a: Int, x: Int, o: Int, e: Int, o2: Int, e2: Int, band: Band
) -> Int:
    """The best score of an alignment from both sequences' first letters to any cell, by Gotoh's
    recurrence over scores: a match earning `a`, a substitution costing `x` and a gap the less of its
    pieces, the second's only with `o2` not negative, every cell inside `band`. Aligning nothing scores
    zero. Shares no code with the wavefront."""
    comptime LOW = -(1 << 40)
    var p = first.as_bytes()
    var q = second.as_bytes()
    var n = len(p)
    var m = len(q)
    var width = m + 1
    var best = List[Int](length=(n + 1) * width, fill=LOW)
    var layers = List[List[Int]]()
    for _ in range(4):
        layers.append(List[Int](length=(n + 1) * width, fill=LOW))
    var answer = 0
    for i in range(n + 1):
        for j in range(m + 1):
            var at = i * width + j
            if not band.holds(i - j):
                continue
            if i == 0 and j == 0:
                best[at] = 0
                continue
            var value = LOW
            for piece in range(2 if o2 >= 0 else 1):
                var opening = o if piece == 0 else o2
                var extension = e if piece == 0 else e2
                if i > 0:
                    layers[2 * piece][at] = max(
                        best[at - width] - opening - extension, layers[2 * piece][at - width] - extension
                    )
                if j > 0:
                    layers[2 * piece + 1][at] = max(
                        best[at - 1] - opening - extension, layers[2 * piece + 1][at - 1] - extension
                    )
                value = max(value, max(layers[2 * piece][at], layers[2 * piece + 1][at]))
            if i > 0 and j > 0:
                value = max(value, best[at - width - 1] + (a if p[i - 1] == q[j - 1] else -x))
            best[at] = value
            answer = max(answer, value)
    return answer


def local_optimum(first: String, second: String, a: Int, x: Int, o: Int, e: Int, o2: Int, e2: Int) -> Int:
    """The best score of an alignment of any part of each sequence, Smith-Waterman's recurrence over
    the whole matrix, every cell floored at zero; shares no code with the library."""
    comptime LOW = -(1 << 40)
    var p = first.as_bytes()
    var q = second.as_bytes()
    var n = len(p)
    var m = len(q)
    var width = m + 1
    var best = List[Int](length=(n + 1) * width, fill=0)
    var layers = List[List[Int]]()
    for _ in range(4):
        layers.append(List[Int](length=(n + 1) * width, fill=LOW))
    var answer = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            var at = i * width + j
            var value = 0
            for piece in range(2 if o2 >= 0 else 1):
                var opening = o if piece == 0 else o2
                var extension = e if piece == 0 else e2
                layers[2 * piece][at] = max(
                    best[at - width] - opening - extension, layers[2 * piece][at - width] - extension
                )
                layers[2 * piece + 1][at] = max(
                    best[at - 1] - opening - extension, layers[2 * piece + 1][at - 1] - extension
                )
                value = max(value, max(layers[2 * piece][at], layers[2 * piece + 1][at]))
            value = max(value, best[at - width - 1] + (a if p[i - 1] == q[j - 1] else -x))
            best[at] = value
            answer = max(answer, value)
    return answer


def rewarded_optimum(
    first: String, second: String, a: Int, x: Int, o: Int, e: Int, o2: Int, e2: Int, ends: EndsFree
) -> Int:
    """The best score, a match earning `a`, of an alignment with the letters `ends` allows free at either
    end: it starts on the first row or column, free within those letters and past them paying a gap,
    and ends on the last row or column within the letters free there. Gotoh's recurrence over the whole
    matrix; shares no code with the library."""
    comptime LOW = -(1 << 40)
    var p = first.as_bytes()
    var q = second.as_bytes()
    var n = len(p)
    var m = len(q)
    var width = m + 1

    def edge(letters: Int, free: Int) {imm o, imm e, imm o2, imm e2} -> Int:
        """The score of reaching the first row or column `letters` in: nothing within `free`, a gap past it."""
        if letters <= free:
            return 0
        var paid = o + e * (letters - free)
        if o2 >= 0:
            paid = min(paid, o2 + e2 * (letters - free))
        return -paid

    var best = List[Int](length=(n + 1) * width, fill=0)
    var layers = List[List[Int]]()
    for _ in range(4):
        layers.append(List[Int](length=(n + 1) * width, fill=LOW))
    for i in range(n + 1):
        best[i * width] = edge(i, ends.first_begin)
    for j in range(m + 1):
        best[j] = edge(j, ends.second_begin)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            var at = i * width + j
            var value = LOW
            for piece in range(2 if o2 >= 0 else 1):
                var opening = o if piece == 0 else o2
                var extension = e if piece == 0 else e2
                layers[2 * piece][at] = max(
                    best[at - width] - opening - extension, layers[2 * piece][at - width] - extension
                )
                layers[2 * piece + 1][at] = max(
                    best[at - 1] - opening - extension, layers[2 * piece + 1][at - 1] - extension
                )
                value = max(value, max(layers[2 * piece][at], layers[2 * piece + 1][at]))
            value = max(value, best[at - width - 1] + (a if p[i - 1] == q[j - 1] else -x))
            best[at] = value
    var answer = LOW
    for j in range(m + 1):
        if m - j <= ends.second_end:
            answer = max(answer, best[n * width + j])
    for i in range(n + 1):
        if n - i <= ends.first_end:
            answer = max(answer, best[i * width + m])
    return answer


def extension_price(cigar: String, a: Int, x: Int, o: Int, e: Int, o2: Int, e2: Int) -> Int:
    """What a CIGAR's alignment scores, a match earning `a`."""
    var total = 0
    var length = 0
    for byte in cigar.as_bytes():
        if byte >= UInt8(ord("0")) and byte <= UInt8(ord("9")):
            length = length * 10 + Int(byte - UInt8(ord("0")))
            continue
        if byte == UInt8(ord("=")):
            total += a * length
        elif byte == UInt8(ord("X")):
            total -= x * length
        else:
            total -= min(o + e * length, o2 + e2 * length) if o2 >= 0 else o + e * length
        length = 0
    return total


def reversed_text(text: String) -> String:
    """`text`'s bytes in reverse order."""
    return reversed_bytes(text.as_bytes())


def reversed_cigar(cigar: String) -> String:
    """A CIGAR's entries in the other order, the alignment read from its far end."""
    var entries = List[String]()
    var start = 0
    var bytes = cigar.as_bytes()
    for index in range(len(bytes)):
        if bytes[index] < UInt8(ord("0")) or bytes[index] > UInt8(ord("9")):
            entries.append(String(cigar[byte = start : index + 1]))
            start = index + 1
    var out = String()
    for index in range(len(entries) - 1, -1, -1):
        out += entries[index]
    return out


def cigar_of_moves(first: String, second: String, moves: List[UInt8]) -> String:
    """The CIGAR of moves `solve` appended right to left: 0 two letters aligned, 1 the first's alone, 2 the second's."""
    var a = first.as_bytes()
    var b = second.as_bytes()
    var out = String()
    var last = UInt8(0)
    var run = 0
    var i = 0
    var j = 0
    for index in range(len(moves) - 1, -1, -1):
        var move = moves[index]
        var op: UInt8
        if move == 0:
            op = UInt8(ord("=")) if a[i] == b[j] else UInt8(ord("X"))
            i += 1
            j += 1
        elif move == 1:
            op = UInt8(ord("D"))
            i += 1
        else:
            op = UInt8(ord("I"))
            j += 1
        if op != last and run > 0:
            out += String(run, chr(Int(last)))
            run = 0
        last = op
        run += 1
    if run > 0:
        out += String(run, chr(Int(last)))
    return out


def test_two_piece_matches_the_full_matrix() raises:
    """Two-piece gap costs, a gap the cheaper of two affine costs: the cost is the whole matrix's with
    a layer either way per piece, the CIGAR spells an alignment that costs it, and the cap, the ends
    free and the split of a pair too large to keep all agree, the split's crossings inside a gap of
    either piece included."""
    seed(43)
    for costs in [(4, 6, 2, 24, 1), (3, 10, 1, 2, 3), (2, 0, 3, 7, 1), (5, 3, 2, 3, 2)]:
        var x = costs[0]
        var o = costs[1]
        var e = costs[2]
        var o2 = costs[3]
        var e2 = costs[4]
        var penalties = affine2p_penalties(x, o, e, o2, e2)
        for trial in range(60):
            var core = random_sequence(1, 200, DNA_ALPHABET)
            var first = core
            var ends = EndsFree()
            if trial % 2 == 1:
                first = random_sequence(0, 40, DNA_ALPHABET) + core + random_sequence(0, 40, DNA_ALPHABET)
                var sizes = [0, 3, 40, 1000]
                ends = EndsFree(
                    sizes[trial % 4], sizes[(trial // 4) % 4], sizes[(trial // 2) % 4], sizes[(trial // 3) % 4]
                )
            var second = mutated(core, [0.0, 0.03, 0.15][trial % 3], [1, 5, 40][(trial // 3) % 3])
            if trial % 5 == 0:
                first, second = second, first
            var expected = ends_free_optimum(first, second, x, o, e, ends, o2, e2)
            assert_equal(distance(first, second, Costs.two_piece(x, o, e, o2, e2), ends_mode(ends)), expected)
            var found = whole(align(first, second, Costs.two_piece(x, o, e, o2, e2), ends_mode(ends)), first, second)
            assert_equal(found.cost, expected)
            var rows = rows_from_cigar(first, second, found.cigar)
            assert_equal(rows[0].replace("-", ""), first)
            assert_equal(rows[1].replace("-", ""), second)
            assert_equal(ends_free_price(found.cigar, x, o, e, ends, o2, e2), expected)
            var capped = align(first, second, Costs.two_piece(x, o, e, o2, e2), ends_mode(ends), max_cost=expected)
            assert_equal(capped.value().cost, expected)
            if expected > 0:
                assert_false(
                    Bool(
                        distance(
                            first, second, Costs.two_piece(x, o, e, o2, e2), ends_mode(ends), max_cost=expected - 1
                        )
                    )
                )
            for limit in [0, 64]:
                var split = (
                    free_ends_alignment[2](first, second, penalties, True, Int.MAX, ends, Band(), Ties.LEFT, limit)
                    .value()
                    .copy()
                )
                assert_equal(split.cost, expected)
                var spanned = Alignment(
                    split.cost,
                    -split.cost,
                    split.cigar,
                    split.first_start,
                    split.first_end,
                    split.second_start,
                    split.second_end,
                )
                assert_equal(ends_free_price(whole(spanned, first, second).cigar, x, o, e, ends, o2, e2), expected)
    # A long gap at the second piece, 24 + 30, and a short one at the first, 6 + 2: 62, where the first
    # piece alone charges 74 and the second 79.
    var gapped = "GATTACAGCTTGCA" + "C" * 30 + "TGGACCATGAGTCA" + "TTGACCAGTCGATC"
    var plain = "GATTACAGCTTGCA" + "TGGACCATGAGTCA" + "G" + "TTGACCAGTCGATC"
    var both = align(gapped, plain, Costs.two_piece(4, 6, 2, 24, 1))
    assert_equal(both.cost, 62)
    assert_equal(both.cigar, "14=30D14=1I14=")
    assert_equal(distance(gapped, plain, Costs.affine(4, 6, 2)), 74)
    assert_equal(distance(gapped, plain, Costs.affine(4, 24, 1)), 79)
    assert_equal(distance("", "ACGT", Costs.two_piece(4, 6, 2, 24, 1)), 14)
    assert_equal(distance("", "A" * 30, Costs.two_piece(4, 6, 2, 24, 1)), 54)
    with assert_raises():
        _ = distance("A", "C", Costs.two_piece(4, 6, 2, 24, 0))


def test_affine_band_matches_the_full_matrix() raises:
    """With a band, the cost is the whole matrix's over its diagonals alone, for one gap piece or two,
    global or with ends free: the CIGAR costs it and stays inside, the cap agrees, a band no alignment
    fits raises or comes back empty, and the split of a pair too large to keep, whose pieces see the
    band from their own origins, finds the same cost."""
    seed(47)
    for costs in [(4, 6, 2, -1, 0), (1, 0, 1, -1, 0), (4, 6, 2, 24, 1), (3, 2, 3, 10, 1)]:
        var x = costs[0]
        var o = costs[1]
        var e = costs[2]
        var o2 = costs[3]
        var e2 = costs[4]
        var two = o2 >= 0
        var penalties = affine2p_penalties(x, o, e, o2, e2) if two else affine_penalties(x, o, e)
        for trial in range(80):
            var core = random_sequence(1, 200, DNA_ALPHABET)
            var first = core
            var ends = EndsFree()
            if trial % 3 == 2:
                first = random_sequence(0, 30, DNA_ALPHABET) + core + random_sequence(0, 30, DNA_ALPHABET)
                var sizes = [0, 3, 30, 1000]
                ends = EndsFree(
                    sizes[trial % 4], sizes[(trial // 4) % 4], sizes[(trial // 2) % 4], sizes[(trial // 3) % 4]
                )
            var second = mutated(core, [0.03, 0.15, 0.3][trial % 3], [1, 5, 30][(trial // 3) % 3])
            if trial % 5 == 0:
                first, second = second, first
            var width = Int(random_ui64(0, 40))
            var shift = Int(random_ui64(0, 20)) - 10
            var band = Band(-width + shift, width + shift)
            var expected = ends_free_optimum(first, second, x, o, e, ends, o2 if two else -1, e2, band)
            if expected >= 1 << 40:
                with assert_raises(contains="band"):
                    if two:
                        _ = distance(first, second, Costs.two_piece(x, o, e, o2, e2), ends_mode(ends), band=band)
                    else:
                        _ = distance(first, second, Costs.affine(x, o, e), ends_mode(ends), band=band)
                with assert_raises(contains="band"):
                    if two:
                        _ = whole(
                            align(first, second, Costs.two_piece(x, o, e, o2, e2), ends_mode(ends), band=band),
                            first,
                            second,
                        )
                    else:
                        _ = whole(
                            align(first, second, Costs.affine(x, o, e), ends_mode(ends), band=band), first, second
                        )
                assert_false(
                    Bool(align(first, second, Costs.affine(x, o, e), ends_mode(ends), max_cost=1 << 30, band=band))
                )
                continue
            var least: Int
            var cigar: String
            var cost: Int
            if two:
                least = distance(first, second, Costs.two_piece(x, o, e, o2, e2), ends_mode(ends), band=band)
                var found = whole(
                    align(first, second, Costs.two_piece(x, o, e, o2, e2), ends_mode(ends), band=band), first, second
                )
                cost = found.cost
                cigar = found.cigar
            else:
                least = distance(first, second, Costs.affine(x, o, e), ends_mode(ends), band=band)
                var found = whole(
                    align(first, second, Costs.affine(x, o, e), ends_mode(ends), band=band), first, second
                )
                cost = found.cost
                cigar = found.cigar
            assert_equal(least, expected)
            assert_equal(cost, expected)
            var rows = rows_from_cigar(first, second, cigar)
            assert_equal(rows[0].replace("-", ""), first)
            assert_equal(rows[1].replace("-", ""), second)
            assert_equal(ends_free_price(cigar, x, o, e, ends, o2 if two else -1, e2), expected)
            assert_true(stays_inside(cigar, ends, band), String("the CIGAR leaves the band: ", cigar))
            if not two:
                assert_equal(
                    align(first, second, Costs.affine(x, o, e), ends_mode(ends), max_cost=expected, band=band)
                    .value()
                    .cost,
                    expected,
                )
                if expected > 0:
                    assert_false(
                        Bool(
                            distance(
                                first, second, Costs.affine(x, o, e), ends_mode(ends), max_cost=expected - 1, band=band
                            )
                        )
                    )
            for limit in [0, 64]:
                var split = free_ends_alignment[2](
                    first, second, penalties, True, Int.MAX, ends, band, Ties.LEFT, limit
                ) if two else free_ends_alignment[1](
                    first, second, penalties, True, Int.MAX, ends, band, Ties.LEFT, limit
                )
                assert_equal(split.value().cost, expected)
                ref piece = split.value()
                var spelled = whole(
                    Alignment(
                        piece.cost,
                        -piece.cost,
                        piece.cigar,
                        piece.first_start,
                        piece.first_end,
                        piece.second_start,
                        piece.second_end,
                    ),
                    first,
                    second,
                ).cigar
                assert_equal(ends_free_price(spelled, x, o, e, ends, o2 if two else -1, e2), expected)
                assert_true(stays_inside(spelled, ends, band), String("the split's CIGAR leaves the band: ", spelled))
    # A band of one diagonal admits only the diagonal itself: four substitutions, where a gap either
    # way would be cheaper.
    assert_equal(align("ACGTACGT", "ACGAACGT", Costs.affine(4, 6, 2), band=Band(0, 0)).cigar, "3=1X4=")
    assert_equal(align("AAAACCCC", "CCCCAAAA", Costs.affine(4, 6, 2), band=Band.around(0)).cost, 32)
    with assert_raises(contains="band"):
        _ = distance("ACGT", "AC", Costs.affine(4, 6, 2), band=Band.around(1))


def single_search_cigar[
    pieces: Int
](first: String, second: String, penalties: Penalties, ends: EndsFree, band: Band) -> Tuple[Int, String]:
    """WFA2-lib's alignment as it finds it: one search from the origin, unpruned, every cost kept,
    until an end the trailing allowances admit is reached, the first on the highest diagonal, traced
    back from there. The reference `Ties.RIGHT` must match however the library found the cost."""
    var columns = first.byte_length()
    var rows = second.byte_length()
    var search = Wavefront[pieces](
        first.as_bytes(),
        second.as_bytes(),
        penalties,
        FREE_START,
        True,
        False,
        ends.first_begin,
        ends.second_begin,
        band,
    )
    while True:
        var cost = search.cost
        var history = Pointer(to=search.history)
        for diagonal in range(history[].highs[cost], history[].lows[cost] - 1, -1):
            var column = history[].column(cost, diagonal)
            if column < 0:
                continue
            var row = column - diagonal
            var trailing = -1
            var along_first = True
            if row >= rows and columns - column <= ends.first_end:
                trailing = columns - column
            elif column >= columns and rows - row <= ends.second_end:
                trailing = rows - row
                along_first = False
            if trailing >= 0:
                var moves = List[UInt8]()
                for _ in range(trailing):
                    moves.append(UInt8(FIRST_GAP) if along_first else UInt8(SECOND_GAP))
                trace(history[], penalties, ALIGNED, cost, diagonal, column, moves)
                return (cost * penalties.scale, cigar_of_moves(first, second, moves))
        search.advance[True]()


def cost_matrix(
    first: String, second: String, x: Int, o: Int, e: Int, o2: Int, e2: Int, starts: EndsFree, band: Band
) -> List[Int]:
    """Every cell's least gap-affine cost from a start `starts` frees, by Gotoh's recurrence over the
    whole matrix, `(len(second) + 1)` cells a row of the first's letters; `1 << 40` off `band`."""
    comptime HIGH = 1 << 40
    var a = first.as_bytes()
    var b = second.as_bytes()
    var n = len(a)
    var m = len(b)
    var width = m + 1
    var best = List[Int](length=(n + 1) * width, fill=HIGH)
    var layers = List[List[Int]]()
    for _ in range(4):
        layers.append(List[Int](length=(n + 1) * width, fill=HIGH))
    for i in range(n + 1):
        for j in range(m + 1):
            var at = i * width + j
            if not band.holds(i - j):
                continue
            if (j == 0 and i <= starts.first_begin) or (i == 0 and j <= starts.second_begin):
                best[at] = 0
                continue
            var value = HIGH
            for piece in range(2 if o2 >= 0 else 1):
                var opening = o if piece == 0 else o2
                var extension = e if piece == 0 else e2
                if i > 0:
                    layers[2 * piece][at] = min(
                        best[at - width] + opening + extension, layers[2 * piece][at - width] + extension
                    )
                if j > 0:
                    layers[2 * piece + 1][at] = min(
                        best[at - 1] + opening + extension, layers[2 * piece + 1][at - 1] + extension
                    )
                value = min(value, min(layers[2 * piece][at], layers[2 * piece + 1][at]))
            if i > 0 and j > 0:
                value = min(value, best[at - width - 1] + (0 if a[i - 1] == b[j - 1] else x))
            best[at] = value
    return best^


def rule_span(
    first: String, second: String, x: Int, o: Int, e: Int, o2: Int, e2: Int, ends: EndsFree, band: Band
) -> Tuple[Int, Int, Int, Int]:
    """The span `Ties.LEFT` names for free ends, from the whole matrix: of the ends an optimal alignment
    reaches, the one on the highest diagonal, and of the starts an optimal alignment ending there
    leaves from, the one on the highest diagonal too. Start column and row, end column and row."""
    var n = first.byte_length()
    var m = second.byte_length()
    var width = m + 1
    var forward = cost_matrix(first, second, x, o, e, o2, e2, ends, band)
    var least = 1 << 40
    var end = (0, 0)
    for i in range(n + 1):
        for j in range(m + 1):
            if not ((i == n and m - j <= ends.second_end) or (j == m and n - i <= ends.first_end)):
                continue
            var value = forward[i * width + j]
            if value < least or (value == least and i - j > end[0] - end[1]):
                least = value
                end = (i, j)
    # Back from that end alone over both reversed, to the starts.
    var head = reversed_text(String(StringSlice(unsafe_from_utf8=first.as_bytes()[: end[0]])))
    var lead = reversed_text(String(StringSlice(unsafe_from_utf8=second.as_bytes()[: end[1]])))
    var backward = cost_matrix(head, lead, x, o, e, o2, e2, EndsFree(), band.mirrored(end[0] - end[1]))
    var start = (-1, -1)
    var back_width = end[1] + 1
    for i in range(end[0] + 1):
        for j in range(end[1] + 1):
            # The forward cell `(end[0] - i, end[1] - j)`, a start if on an edge within the allowance.
            var column = end[0] - i
            var row = end[1] - j
            if not ((row == 0 and column <= ends.first_begin) or (column == 0 and row <= ends.second_begin)):
                continue
            if backward[i * back_width + j] != least:
                continue
            if start[0] < 0 or column - row > start[0] - start[1]:
                start = (column, row)
    return (start[0], start[1], end[0], end[1])


def test_ties_follow_a_fixed_rule() raises:
    """Of several equally good alignments the CIGAR is always the one `Ties` names, however the two
    searches found the cost: `Ties.RIGHT` is one search traced back from the far end, as WFA2-lib's,
    and `Ties.LEFT` the same rule run from the start, the right rule's CIGAR of both sequences
    reversed, read backwards. With free ends the span comes first, for the left rule the end and then
    the start on the highest diagonal an optimum allows (see `rule_span`). For one gap piece or two, global, with ends
    free, inside a band."""
    seed(59)
    for costs in [(4, 6, 2, -1, 0), (1, 0, 1, -1, 0), (3, 1, 4, -1, 0), (4, 6, 2, 24, 1), (2, 2, 3, 9, 1)]:
        var x = costs[0]
        var o = costs[1]
        var e = costs[2]
        var o2 = costs[3]
        var e2 = costs[4]
        var two = o2 >= 0
        var penalties = affine2p_penalties(x, o, e, o2, e2) if two else affine_penalties(x, o, e)
        for trial in range(60):
            # Repeats make ties: gaps that could sit anywhere in a run, substitutions that trade for gaps.
            var unit = random_sequence(1, 4, DNA_ALPHABET)
            var core = random_sequence(1, 60, DNA_ALPHABET)
            for _ in range(Int(random_ui64(0, 6))):
                core += unit + random_sequence(0, 20, DNA_ALPHABET)
            var first = core
            var second = mutated(core, [0.03, 0.1, 0.25][trial % 3], [1, 4, 12][(trial // 3) % 3])
            if second.byte_length() == 0:
                second = "A"
            var ends = EndsFree()
            if trial % 3 == 1:
                var sizes = [0, 2, 9, 1000]
                ends = EndsFree(
                    sizes[trial % 4], sizes[(trial // 4) % 4], sizes[(trial // 2) % 4], sizes[(trial // 3) % 4]
                )
            var band = Band() if trial % 4 != 3 else Band.around(Int(random_ui64(0, 40)))
            var target = first.byte_length() - second.byte_length()
            if not band.holds(0) or not band.holds(target):
                band = Band()
            var reference = single_search_cigar[2](
                first, second, penalties, ends, band
            ) if two else single_search_cigar[1](first, second, penalties, ends, band)
            var right: String
            var left: String
            var mirrored: String
            var flipped = EndsFree(ends.first_end, ends.first_begin, ends.second_end, ends.second_begin)
            var back = band.mirrored(target)
            if two:
                right = whole(
                    align(first, second, Costs.two_piece(x, o, e, o2, e2), ends_mode(ends), band=band, ties=Ties.RIGHT),
                    first,
                    second,
                ).cigar
                left = whole(
                    align(first, second, Costs.two_piece(x, o, e, o2, e2), ends_mode(ends), band=band), first, second
                ).cigar
                mirrored = whole(
                    align(
                        reversed_text(first),
                        reversed_text(second),
                        Costs.two_piece(x, o, e, o2, e2),
                        ends_mode(flipped),
                        band=back,
                        ties=Ties.RIGHT,
                    ),
                    reversed_text(first),
                    reversed_text(second),
                ).cigar
            else:
                right = whole(
                    align(first, second, Costs.affine(x, o, e), ends_mode(ends), band=band, ties=Ties.RIGHT),
                    first,
                    second,
                ).cigar
                left = whole(
                    align(first, second, Costs.affine(x, o, e), ends_mode(ends), band=band), first, second
                ).cigar
                mirrored = whole(
                    align(
                        reversed_text(first),
                        reversed_text(second),
                        Costs.affine(x, o, e),
                        ends_mode(flipped),
                        band=back,
                        ties=Ties.RIGHT,
                    ),
                    reversed_text(first),
                    reversed_text(second),
                ).cigar
            # With free ends the span comes first, by the rule, and the letters between are a global
            # alignment by the same rule: the left one WFA2-lib's trace over both reversed, read back.
            var span = rule_span(first, second, x, o, e, o2, e2, ends, band)
            var part = reversed_text(String(StringSlice(unsafe_from_utf8=first.as_bytes()[span[0] : span[2]])))
            var piece = reversed_text(String(StringSlice(unsafe_from_utf8=second.as_bytes()[span[1] : span[3]])))
            var inner_band = band.shifted(span[0] - span[1]).mirrored(part.byte_length() - piece.byte_length())
            var inner = single_search_cigar[2](
                part, piece, penalties, EndsFree(), inner_band
            ) if two else single_search_cigar[1](part, piece, penalties, EndsFree(), inner_band)
            var expected = whole(
                Alignment(inner[0], -inner[0], reversed_cigar(inner[1]), span[0], span[2], span[1], span[3]),
                first,
                second,
            ).cigar
            assert_equal(left, expected, String("the left rule's span, trial ", trial))
            assert_equal(left, reversed_cigar(mirrored), String("the left rule, trial ", trial))
            if ends.first_begin + ends.first_end + ends.second_begin + ends.second_end > 0:
                continue
            assert_equal(right, reference[1], String("WFA2-lib's rule, trial ", trial))
            # The split of a pair too large to keep follows the rule within its pieces; one that fits
            # whole, whatever the limit, follows it throughout.
            for limit in [1 << 20, 1 << 30]:
                var moves = List[UInt8]()
                if two:
                    _ = solve[2](
                        first.as_bytes(),
                        second.as_bytes(),
                        penalties,
                        FREE_START,
                        FREE_START,
                        limit,
                        moves,
                        True,
                        Int.MAX,
                        band,
                        Ties.RIGHT,
                    )
                else:
                    _ = solve[1](
                        first.as_bytes(),
                        second.as_bytes(),
                        penalties,
                        FREE_START,
                        FREE_START,
                        limit,
                        moves,
                        True,
                        Int.MAX,
                        band,
                        Ties.RIGHT,
                    )
                assert_equal(cigar_of_moves(first, second, moves), reference[1])
    # A gap in a run of repeats sits at its left end by default, at its right end for WFA2-lib's rule.
    assert_equal(align("ACGTTTTACG", "ACGTTTACG", Costs.affine(4, 6, 2)).cigar, "3=1D6=")
    assert_equal(align("ACGTTTTACG", "ACGTTTACG", Costs.affine(4, 6, 2), ties=Ties.RIGHT).cigar, "6=1D3=")


def test_edit_ties_follow_a_fixed_rule() raises:
    """At unit costs the alignment is the one `Ties` names whichever search finds the distance: near
    pairs settled from one end, moderate ones from both, divergent ones by band doubling, and under a
    cap the wavefront. `Ties.RIGHT` is one front from the start traced back from the corner, WFA2-lib's,
    and `Ties.LEFT` the right rule's CIGAR of both sequences reversed, read backwards."""
    seed(61)
    for trial in range(90):
        var unit = random_sequence(1, 4, DNA_ALPHABET)
        var core = random_sequence(20, [300, 3000, 12000][trial % 3], DNA_ALPHABET)
        for _ in range(Int(random_ui64(0, 8))):
            core += unit * 2 + random_sequence(0, 30, DNA_ALPHABET)
        var first = core
        var second = mutated(core, [0.01, 0.06, 0.15, 0.3][(trial // 3) % 4], [1, 3, 12][(trial // 12) % 3])
        if second.byte_length() == 0:
            second = "A"
        var profile = Profile(first, second)
        var fronts = DiagonalFronts()
        var reference = diagonal_transition(profile, 1 << 30, fronts)
        var right = align(first, second, ties=Ties.RIGHT)
        if reference.distance >= 0:
            var moves = List[UInt8]()
            trace_diagonals(profile, fronts, reference.distance, moves)
            var path = EditPath(moves^, List[UInt8](), first.byte_length(), second.byte_length(), reference.distance)
            assert_equal(right.cost, reference.distance)
            assert_equal(right.cigar, cigar_string(first, second, path, True), String("the right rule, trial ", trial))
        var left = align(first, second)
        var mirrored = align(reversed_text(first), reversed_text(second), ties=Ties.RIGHT)
        assert_equal(left.cost, right.cost)
        assert_equal(left.cigar, reversed_cigar(mirrored.cigar), String("the left rule, trial ", trial))
        var capped = align(first, second, max_cost=right.cost, ties=Ties.RIGHT)
        assert_equal(capped.value().cigar, right.cigar, String("the wavefront's right rule, trial ", trial))
        assert_equal(align(first, second, max_cost=left.cost).value().cigar, left.cigar)
    # A gap in a run of repeats: at its left end by default, at its right end under WFA2-lib's rule.
    assert_equal(align("ACGTTTTACG", "ACGTTTACG").cigar, "3=1D6=")
    assert_equal(align("ACGTTTTACG", "ACGTTTACG", ties=Ties.RIGHT).cigar, "6=1D3=")


def test_affine_extension_matches_the_full_matrix() raises:
    """An extension fixed at the start or at the end scores the best of every alignment from that end
    to any cell, for one gap piece or two, with a band or without: its CIGAR spells an alignment of
    the letters it covers from that end, scores what it claims, and stays inside the band."""
    seed(53)
    for scores in [(1, 4, 6, 2, -1, 0), (2, 4, 4, 2, -1, 0), (2, 4, 4, 2, 24, 1), (3, 1, 0, 2, 5, 1)]:
        var a = scores[0]
        var x = scores[1]
        var o = scores[2]
        var e = scores[3]
        var o2 = scores[4]
        var e2 = scores[5]
        var two = o2 >= 0
        for trial in range(80):
            var core = random_sequence(0, 200, DNA_ALPHABET)
            var first = core + random_sequence(0, 60, DNA_ALPHABET)
            var second = mutated(core, [0.0, 0.05, 0.2][trial % 3], [1, 5, 30][(trial // 3) % 3]) + random_sequence(
                0, 60, DNA_ALPHABET
            )
            if trial % 7 == 0:
                first, second = second, first
            var band = Band()
            if trial % 2 == 1:
                band = Band.around(Int(random_ui64(0, 30)))
            for anchor in [Anchor.START, Anchor.END]:
                var at_end = anchor == Anchor.END
                var left = reversed_text(first) if at_end else first
                var right = reversed_text(second) if at_end else second
                var expected = extension_optimum(left, right, a, x, o, e, o2, e2, band)
                var found = align(
                    left if at_end else first,
                    right if at_end else second,
                    Costs.two_piece(x, o, e, o2, e2),
                    Mode.extension(a),
                    band=band,
                ) if two else align(
                    left if at_end else first,
                    right if at_end else second,
                    Costs.affine(x, o, e),
                    Mode.extension(a),
                    band=band,
                )
                # Asked from the end of the reversed texts, which the oracle reads from their start.
                if at_end:
                    found = align(
                        first, second, Costs.two_piece(x, o, e, o2, e2), Mode.extension(a, Anchor.END), band=band
                    ) if two else align(first, second, Costs.affine(x, o, e), Mode.extension(a, Anchor.END), band=band)
                assert_equal(found.score, expected)
                var covered_first = String(first[byte = found.reference_start : found.reference_end])
                var covered_second = String(second[byte = found.query_start : found.query_end])
                if at_end:
                    assert_equal(found.reference_end, first.byte_length())
                    assert_equal(found.query_end, second.byte_length())
                else:
                    assert_equal(found.reference_start, 0)
                    assert_equal(found.query_start, 0)
                assert_equal(found.cost, a * matches_in(found.cigar) - found.score)
                if covered_first.byte_length() + covered_second.byte_length() > 0:
                    var rows = rows_from_cigar(covered_first, covered_second, found.cigar)
                    assert_equal(rows[0].replace("-", ""), covered_first)
                    assert_equal(rows[1].replace("-", ""), covered_second)
                else:
                    assert_equal(found.cigar, "")
                assert_equal(extension_price(found.cigar, a, x, o, e, o2, e2), expected)
                var walked = reversed_cigar(found.cigar) if at_end else found.cigar
                assert_true(
                    stays_inside(walked, EndsFree(), band), String("the extension leaves the band: ", found.cigar)
                )
    # A read that matches its reference for twelve letters, then not at all: the extension stops there.
    var read = align(
        "ACGTTGCAAGGC" + "TTTTTTTTTT", "ACGTTGCAAGGC" + "GAGAGAGAGA", Costs.affine(4, 6, 2), Mode.extension(1)
    )
    assert_equal(read.score, 12)
    assert_equal(read.cigar, "12=")
    var back = align(
        "TTTTTTTTTT" + "ACGTTGCAAGGC",
        "GAGAGAGAGA" + "ACGTTGCAAGGC",
        Costs.affine(4, 6, 2),
        Mode.extension(1, Anchor.END),
    )
    assert_equal(back.score, 12)
    assert_equal(back.reference_end - back.reference_start, 12)
    assert_equal(back.reference_start, 10)
    assert_equal(align("ACGT", "TGCA", Costs.affine(4, 6, 2), Mode.extension(0)).cigar, "")


def test_every_mode_matches_the_full_matrix() raises:
    """Every cost model in every mode is the whole matrix's optimum: unit costs, linear, affine and
    two-piece gaps, globally, a query inside, at the start or at the end of a reference, with random
    free ends, and as an extension. The CIGAR spans what it claims and prices to the cost with its free
    letters put back, and the capped search, the wavefront's whatever the costs, agrees with the
    uncapped one, the bit-parallel sweep's at unit costs."""
    seed(71)
    var models: List[Costs] = [Costs.edit(), Costs.linear(2, 3), Costs.affine(4, 6, 2), Costs.two_piece(4, 6, 2, 24, 1)]
    for costs in models:
        var x = costs.mismatch
        var o = costs.opening
        var e = costs.extension
        var o2 = costs.opening2
        var e2 = costs.extension2
        for trial in range(40):
            var core = random_sequence(1, 150, DNA_ALPHABET)
            var reference = random_sequence(0, 40, DNA_ALPHABET) + core + random_sequence(0, 40, DNA_ALPHABET)
            var query = mutated(core, [0.0, 0.05, 0.2][trial % 3], [1, 4, 20][(trial // 3) % 3])
            var sizes = [0, 2, 30, 1000]
            var modes: List[Mode] = [
                Mode.GLOBAL,
                Mode.INFIX,
                Mode.PREFIX,
                Mode.SUFFIX,
                Mode.ends_free(
                    reference_start=sizes[trial % 4],
                    reference_end=sizes[(trial // 4) % 4],
                    query_start=sizes[(trial // 2) % 4],
                    query_end=sizes[(trial // 3) % 4],
                ),
            ]
            for mode in modes:
                var ends = EndsFree(
                    min(mode.reference_start, reference.byte_length()),
                    min(mode.reference_end, reference.byte_length()),
                    min(mode.query_start, query.byte_length()),
                    min(mode.query_end, query.byte_length()),
                )
                var expected = ends_free_optimum(reference, query, x, o, e, ends, o2, e2)
                assert_equal(distance(reference, query, costs, mode), expected)
                var found = align(reference, query, costs, mode)
                assert_equal(found.cost, expected)
                assert_equal(found.score, -expected)
                var spelled = rows_from_cigar(
                    String(reference[byte = found.reference_start : found.reference_end]),
                    String(query[byte = found.query_start : found.query_end]),
                    found.cigar,
                )
                assert_equal(spelled[0].replace("-", "").byte_length(), found.reference_end - found.reference_start)
                var full = whole(found, reference, query)
                assert_equal(ends_free_price(full.cigar, x, o, e, ends, o2, e2), expected)
                assert_equal(distance(reference, query, costs, mode, max_cost=expected).value(), expected)
                var capped = align(reference, query, costs, mode, max_cost=expected).value().copy()
                assert_equal(capped.cost, expected)
                if expected > 0:
                    assert_false(Bool(align(reference, query, costs, mode, max_cost=expected - 1)))
            var best = extension_optimum(reference, query, 1, x, o, e, o2, e2, Band())
            var eqx = align(reference, query, costs, Mode.extension(1))
            assert_equal(eqx.score, best)
            assert_equal(extension_price(eqx.cigar, 1, x, o, e, o2, e2), best)
            assert_equal(eqx.cost, matches_in(eqx.cigar) - best)
            for ties in [Ties.LEFT, Ties.RIGHT]:
                var local = align(reference, query, costs, Mode.local(2), ties=ties)
                var top = local_optimum(reference, query, 2, x, o, e, o2, e2)
                assert_equal(local.score, top)
                assert_equal(extension_price(local.cigar, 2, x, o, e, o2, e2), top)
                assert_equal(local.cost, 2 * matches_in(local.cigar) - top)
                if top > 0:
                    var part = rows_from_cigar(
                        String(reference[byte = local.reference_start : local.reference_end]),
                        String(query[byte = local.query_start : local.query_end]),
                        local.cigar,
                    )
                    assert_true(part[0].byte_length() > 0)
                else:
                    assert_equal(local.cigar, "")
                    assert_equal(local.reference_end, 0)
            # Every mode of free ends again with a match earning 2, and the overlap, all four free.
            var scored: List[Mode] = [Mode.overlap(2)]
            for mode in modes:
                scored.append(mode.with_match_score(2))
            for mode in scored:
                var ends = EndsFree(
                    min(mode.reference_start, reference.byte_length()),
                    min(mode.reference_end, reference.byte_length()),
                    min(mode.query_start, query.byte_length()),
                    min(mode.query_end, query.byte_length()),
                )
                var most = rewarded_optimum(reference, query, 2, x, o, e, o2, e2, ends)
                for ties in [Ties.LEFT, Ties.RIGHT]:
                    var found = align(reference, query, costs, mode, ties=ties)
                    assert_equal(found.score, most)
                    assert_equal(extension_price(found.cigar, 2, x, o, e, o2, e2), most)
                    assert_equal(found.cost, 2 * matches_in(found.cigar) - most)
                    # It starts within the letters free at the start, and ends having consumed one
                    # sequence, the other's rest within the letters free at the end.
                    assert_true(found.reference_start <= ends.first_begin and found.query_start <= ends.second_begin)
                    assert_true(found.reference_start == 0 or found.query_start == 0)
                    var reference_rest = reference.byte_length() - found.reference_end
                    var query_rest = query.byte_length() - found.query_end
                    assert_true(
                        (reference_rest == 0 and query_rest <= ends.second_end)
                        or (query_rest == 0 and reference_rest <= ends.first_end),
                        String(
                            "spans ",
                            found.reference_start,
                            "..",
                            found.reference_end,
                            " of ",
                            reference.byte_length(),
                            ", ",
                            found.query_start,
                            "..",
                            found.query_end,
                            " of ",
                            query.byte_length(),
                            ", ends ",
                            ends,
                            ", cigar ",
                            found.cigar,
                        ),
                    )
                    if found.cigar.byte_length() > 0:
                        _ = rows_from_cigar(
                            String(reference[byte = found.reference_start : found.reference_end]),
                            String(query[byte = found.query_start : found.query_end]),
                            found.cigar,
                        )
            var inside = align(query, reference, costs, Mode.REFERENCE_IN_QUERY)
            var flipped = align(reference, query, costs, Mode.INFIX)
            assert_equal(inside.cost, flipped.cost)
    # Unit costs inside a reference: the sweep and the wavefront, under a cap, find the same span.
    var placed = align("TTTTACGTACGTTTTT", "ACGTCGT", Costs.edit(), Mode.INFIX)
    assert_equal(placed.cost, 1)
    assert_equal(placed.cigar, "4=1D3=")
    assert_equal(placed.reference_start, 4)
    assert_equal(placed.reference_end, 12)


def test_local_sweeps_agree_in_either_width() raises:
    """The local end found in 16-bit lanes is the one found in 32-bit lanes, and with the lanes along
    the reference the one found with them along the query, for one gap piece or two, pairs close and
    far, short and long."""
    seed(73)
    for costs in [Costs.affine(4, 6, 2), Costs.two_piece(4, 6, 2, 24, 1), Costs.edit()]:
        for trial in range(30):
            var core = random_sequence(1, 600, DNA_ALPHABET)
            var reference = random_sequence(0, 300, DNA_ALPHABET) + core + random_sequence(0, 300, DNA_ALPHABET)
            var query = mutated(core, [0.0, 0.05, 0.2, 0.5][trial % 4], [1, 4, 20][(trial // 4) % 3])
            var a = reference.as_bytes()
            var b = query.as_bytes()
            # Lanes along either sequence, in either width, find the same end.
            var found = List[Tuple[Int, Int, Int, Bool]]()
            if costs.pieces() == 1:
                found.append(best_end[1, DType.int16, 32, False](a, b, costs, 2))
                found.append(best_end[1, DType.int16, 32, True](a, b, costs, 2))
                found.append(best_end[1, DType.int32, 16, False](a, b, costs, 2))
                found.append(best_end[1, DType.int32, 16, True](a, b, costs, 2))
            else:
                found.append(best_end[2, DType.int16, 32, False](a, b, costs, 2))
                found.append(best_end[2, DType.int16, 32, True](a, b, costs, 2))
                found.append(best_end[2, DType.int32, 16, False](a, b, costs, 2))
                found.append(best_end[2, DType.int32, 16, True](a, b, costs, 2))
            for index in range(1, 4):
                assert_equal(found[index][0], found[0][0])
                assert_equal(found[index][1], found[0][1])
                assert_equal(found[index][2], found[0][2])


def test_gotoh_sweeps_agree_in_either_width() raises:
    """A half's last row swept in 16-bit lanes is the one swept in 32-bit lanes, from either end and with
    a deletion run entering open or not, and a local alignment's start found back from its end is the
    same in either, under uniform tables and one telling transitions from transversions."""
    seed(29)
    var regimes = scoring_regimes()
    var transitions: List[Int8] = [2, -3, -1, -3, -3, 2, -3, -1, -1, -3, 2, -3, -3, -1, -3, 2]
    regimes.append(Scoring("ACGT", transitions^, AffineGapCosts(-5, -2)))
    for scoring in regimes:
        var lookup = SubstitutionLookup(scoring.substitutions, scoring.alphabet_size())
        for trial in range(12):
            var first = List[UInt8]()
            var second = List[UInt8]()
            for _ in range(Int(random_ui64(1, 200))):
                first.append(UInt8(random_ui64(0, 3)))
            for _ in range(Int(random_ui64(1, 200))):
                second.append(UInt8(random_ui64(0, 3)))
            var run = GapRun.EXTENDS if trial % 2 == 1 else GapRun.OPENS
            var columns = len(second) + 1
            var narrow_scores = List[Int32](length=columns, fill=0)
            var narrow_deletes = List[Int32](length=columns, fill=0)
            var wide_scores = List[Int32](length=columns, fill=0)
            var wide_deletes = List[Int32](length=columns, fill=0)
            for reverse in [False, True]:
                if reverse:
                    vector_sweep_bands[SweepHalf.REVERSE, DType.int16, 32](
                        first,
                        second,
                        0,
                        len(first),
                        0,
                        len(second),
                        run,
                        lookup,
                        scoring.gaps,
                        narrow_scores,
                        narrow_deletes,
                    )
                    vector_sweep_bands[SweepHalf.REVERSE, DType.int32, 16](
                        first,
                        second,
                        0,
                        len(first),
                        0,
                        len(second),
                        run,
                        lookup,
                        scoring.gaps,
                        wide_scores,
                        wide_deletes,
                    )
                else:
                    vector_sweep_bands[SweepHalf.FORWARD, DType.int16, 32](
                        first,
                        second,
                        0,
                        len(first),
                        0,
                        len(second),
                        run,
                        lookup,
                        scoring.gaps,
                        narrow_scores,
                        narrow_deletes,
                    )
                    vector_sweep_bands[SweepHalf.FORWARD, DType.int32, 16](
                        first,
                        second,
                        0,
                        len(first),
                        0,
                        len(second),
                        run,
                        lookup,
                        scoring.gaps,
                        wide_scores,
                        wide_deletes,
                    )
                for column in range(columns):
                    assert_equal(narrow_scores[column], wide_scores[column])
                    assert_equal(narrow_deletes[column], wide_deletes[column])
            var end = tabulated_end[ANYWHERE](Span(first), Span(second), scoring, EndsFree(), True)
            if end[0] > 0:
                var narrow = reach_back(first, second, end[1], end[2], lookup, scoring.gaps, Int32(end[0]), True)
                var wide = reach_back(first, second, end[1], end[2], lookup, scoring.gaps, Int32(end[0]), False)
                assert_equal(narrow[0], wide[0])
                assert_equal(narrow[1], wide[1])


def test_extreme_arguments_are_held_or_refused() raises:
    """An integer type's limits as a band's bounds, a Z-drop or an end bonus mean no bound at all, never an
    overflow; an alphabet naming `-`, the gapped rows' gap, and a gap score past 32 bits are refused."""
    var reference = "TTACGTACGTTTGCAGG"
    var query = "ACGTCGTTTTGCA"
    var costs = Costs.affine(4, 6, 2)
    var free = align(reference, query, costs)
    for band in [Band.around(Int.MAX), Band(Int.MIN, 1 << 40), Band(Int.MIN, Int.MAX)]:
        var banded = align(reference, query, costs, band=band)
        assert_equal(banded.cost, free.cost)
        assert_equal(banded.cigar, free.cigar)
    # A band holding no diagonal holds no alignment, at either limit.
    for band in [Band(Int.MIN, Int.MIN), Band(Int.MAX, Int.MAX), Band(0, Int.MIN)]:
        with assert_raises():
            _ = align("", "ACG", Costs.edit(), band=band)
        with assert_raises():
            _ = distance("", "ACG", Costs.edit(), band=band)
    var start = "ACGTTGCAAGGCGAGATTGCAAGGCATTACG"
    var other = "ACGTTGCAAGGCCTCTTTGCAAGGCATTACG"
    for drop in [1 << 62, Int.MAX]:
        assert_equal(
            align(start, other, costs, Mode.extension(3, zdrop=drop)).score,
            align(start, other, costs, Mode.extension(3)).score,
        )
        assert_equal(
            score(start, other, Scoring.dna(), Mode.extension(0, zdrop=drop)),
            score(start, other, Scoring.dna(), Mode.extension(0)),
        )
    var read = "ACGTTGCAAGGCGAGAACGT"
    var window = "ACGTTGCAAGGCTTTTACGT"
    var reaching = align(window, read, costs, Mode.extension(1, end_bonus=50))
    for bonus in [1 << 61, 1 << 62, Int.MAX]:
        var found = align(window, read, costs, Mode.extension(1, end_bonus=bonus))
        assert_equal(found.score, reaching.score)
        assert_equal(found.cigar, reaching.cigar)
        assert_equal(
            score(window, read, Scoring.dna(), Mode.extension(0, end_bonus=bonus)),
            score(window, read, Scoring.dna(), Mode.extension(0, end_bonus=50)),
        )
    with assert_raises():
        _ = Scoring.uniform(2, -3, -5, -1, "ACGT-")
    with assert_raises():
        _ = Scoring.uniform(2, -4, -(1 << 40), -1)


def test_lanes_leave_what_they_cannot_hold() raises:
    """A pair longer than the lanes' 16-bit coordinates, and a banded pair whose cost reaches 16 bits' far
    value, are left to their own searches, which the batch then agrees with."""
    var query = "ACGTTGCAACGTGGCATTACGATCGATCGGATCCATGCAAGT"
    var long_reference = query + String("A") * 65536
    var references: List[String] = [long_reference, query]
    var queries: List[String] = [query, long_reference]
    var found = distances(references, queries, Costs.edit())
    for index in range(2):
        assert_equal(found[index], distance(references[index], queries[index], Costs.edit()))
    var ungapped = Band.around(0)
    var costs = Costs.affine(300, 0, 1)
    var all_a: List[String] = [String("A") * 55]
    var all_c: List[String] = [String("C") * 55]
    var single = distance(all_a[0], all_c[0], costs, band=ungapped)
    assert_equal(distances(all_a, all_c, costs, band=ungapped)[0], single)
    assert_false(Bool(distances(all_a, all_c, costs, band=ungapped, max_cost=single - 100)[0]))
    assert_equal(
        distances(all_a, all_c, costs, Mode.PREFIX, band=ungapped)[0],
        distance(all_a[0], all_c[0], costs, Mode.PREFIX, band=ungapped),
    )


def test_end_bonus_reaches_the_end_when_it_pays() raises:
    """An extension's end bonus, KSW2's: the read is aligned to its end once the bonus passes what
    stopping short gains, and not a point before, the comparison strict, as KSW2's is; the alignment
    reaching the end scores its own score, the free ends' best. From either anchor, under `Costs` and a
    `Scoring`, and refused under a band narrower than the pair."""
    var reference = "ACGTTGCAAGGCTTACGATCAGGCGGGGGGGG"
    var query = "ACGTTGCAAGGCTTACGATCAGGCTTTTATTT"
    var costs = Costs.affine(4, 6, 2)
    for from_end in [False, True]:
        var anchor = Anchor.END if from_end else Anchor.START
        var first = reversed_text(reference) if from_end else reference
        var second = reversed_text(query) if from_end else query
        var plain = Mode.extension(1, anchor)
        var stop = align(first, second, costs, plain)
        var reaching = score(first, second, costs, plain.reaching_end())
        assert_true(reaching < stop.score)
        var short = align(first, second, costs, Mode.extension(1, anchor, end_bonus=stop.score - reaching))
        assert_equal(short.score, stop.score)
        assert_equal(short.cigar, stop.cigar)
        var whole = align(first, second, costs, Mode.extension(1, anchor, end_bonus=stop.score - reaching + 1))
        assert_equal(whole.score, reaching)
        assert_equal(whole.query_end - whole.query_start, second.byte_length())
        assert_equal(
            score(first, second, costs, Mode.extension(1, anchor, end_bonus=stop.score - reaching + 1)), reaching
        )
        with assert_raises(contains="end bonus"):
            _ = align(first, second, costs, Mode.extension(1, anchor, end_bonus=5), band=Band.around(2))
        var scoring = Scoring.dna()
        var table_stop = align(first, second, scoring, Mode.extension(0, anchor))
        var table_reaching = score(first, second, scoring, Mode.extension(0, anchor).reaching_end())
        assert_true(table_reaching < table_stop.score)
        var bonus = table_stop.score - table_reaching
        assert_equal(align(first, second, scoring, Mode.extension(0, anchor, end_bonus=bonus)).score, table_stop.score)
        var table_whole = align(first, second, scoring, Mode.extension(0, anchor, end_bonus=bonus + 1))
        assert_equal(table_whole.score, table_reaching)
        assert_equal(table_whole.query_end - table_whole.query_start, second.byte_length())


def test_traced_extension_is_the_searched_one() raises:
    """The extension back from an end, traced as it searched, is the one the extension's own search and
    split finds over both sequences reversed, read backwards: the same stop, cost and CIGAR."""
    seed(79)
    for costs in [Costs.affine(4, 6, 2), Costs.two_piece(4, 6, 2, 24, 1), Costs.edit()]:
        var two = costs.pieces() == 2
        var penalties = extension_penalties(
            2,
            costs.mismatch,
            costs.opening,
            costs.extension,
            costs.opening2 if two else 0,
            costs.extension2 if two else 0,
        )
        for trial in range(40):
            var core = random_sequence(1, 300, DNA_ALPHABET)
            var head = random_sequence(0, 100, DNA_ALPHABET) + core
            var lead = random_sequence(0, 100, DNA_ALPHABET) + mutated(core, [0.0, 0.05, 0.2, 0.4][trial % 4], 8)
            var found = end_of(head.as_bytes(), lead.as_bytes(), costs, 2)
            if found[0] == 0:
                continue
            var a = String(StringSlice(unsafe_from_utf8=head.as_bytes()[: found[1]]))
            var b = String(StringSlice(unsafe_from_utf8=lead.as_bytes()[: found[2]]))
            var traced = traced_extension[2](a, b, penalties, True, found[0]) if two else traced_extension[1](
                a, b, penalties, True, found[0]
            )
            var searched = extension_of[2](
                reversed_text(a), reversed_text(b), penalties, True, Anchor.START, Band(), Ties.RIGHT, found[0]
            ) if two else extension_of[1](
                reversed_text(a), reversed_text(b), penalties, True, Anchor.START, Band(), Ties.RIGHT, found[0]
            )
            assert_equal(traced.value().score, found[0])
            assert_equal(searched.score, found[0])
            assert_equal(traced.value().first_length, searched.first_length)
            assert_equal(traced.value().second_length, searched.second_length)
            assert_equal(traced.value().matches, searched.matches)
            assert_equal(traced.value().cigar, searched.cigar)


# endregion Refusals

# region Output


def test_sam_fields_describe_the_alignment() raises:
    """The clipped CIGAR, the edit count, the identity and the `MD` string of known alignments, and on
    random ones the edit count is the unit cost of the alignment's own edits and the `MD` string
    rebuilds the reference's aligned part from the query and the CIGAR."""
    var costs = Costs.affine(4, 6, 2)
    var core = align("GGGGACGTACGTGGGG", "CCCCACGTTCGTCC", costs, Mode.local(2))
    assert_equal(core.cigar, "4=1X3=")
    assert_equal(core.clipped_cigar(14), "4S4=1X3=2S")
    assert_equal(core.clipped_cigar(14, hard=True), "4H4=1X3=2H")
    assert_equal(core.edit_distance("GGGGACGTACGTGGGG", "CCCCACGTTCGTCC"), 1)
    assert_equal(core.mismatch_string("GGGGACGTACGTGGGG", "CCCCACGTTCGTCC"), "4A3")
    assert_equal(core.identity("GGGGACGTACGTGGGG", "CCCCACGTTCGTCC"), 7.0 / 8.0)
    var gapped = align("ACGTTTGCAAC", "ACGTGCATAC", Costs.edit(), eqx=False)
    assert_equal(gapped.cigar, "3M2D4M1I2M")
    assert_equal(gapped.clipped_cigar(10), "3M2D4M1I2M")
    assert_equal(gapped.mismatch_string("ACGTTTGCAAC", "ACGTGCATAC"), "3^TT6")
    var counted = gapped.counts("ACGTTTGCAAC", "ACGTGCATAC")
    assert_equal(counted.matches, 9)
    assert_equal(counted.mismatches, 0)
    assert_equal(counted.deleted, 2)
    assert_equal(counted.inserted, 1)
    assert_equal(gapped.edit_distance("ACGTTTGCAAC", "ACGTGCATAC"), 3)
    var empty = Alignment(0, 0, String(), 3, 3, 2, 2)
    assert_equal(empty.identity("ACGT", "AC"), 0.0)
    assert_equal(empty.mismatch_string("ACGT", "AC"), "0")
    assert_equal(empty.clipped_cigar(2), "2S")
    seed(83)
    for trial in range(200):
        var reference = random_sequence(0, 80, DNA_ALPHABET)
        var query = mutated(reference, 0.2, 4) if trial % 2 == 0 else random_sequence(0, 80, DNA_ALPHABET)
        var found = align(reference, query, costs, Mode.INFIX if trial % 3 == 0 else Mode.GLOBAL, eqx=trial % 5 != 0)
        var edits = found.edit_distance(reference, query)
        var part = String(
            StringSlice(unsafe_from_utf8=reference.as_bytes()[found.reference_start : found.reference_end])
        )
        var piece = String(StringSlice(unsafe_from_utf8=query.as_bytes()[found.query_start : found.query_end]))
        assert_true(edits >= distance(part, piece))
        assert_equal(rebuilt_reference(piece, found.cigar, found.mismatch_string(reference, query)), part)


def gotoh_cost(first: String, second: String, mismatch: Int, opening: Int, extension: Int) -> Int:
    """The least global cost by Gotoh's full matrix, a gap of `k` letters `opening + k extension`: the
    test's own, sharing nothing with the library."""
    comptime FAR = 1 << 50
    var a = first.as_bytes()
    var b = second.as_bytes()
    var n = len(a)
    var m = len(b)
    var h = List[List[Int]]()
    var d = List[List[Int]]()
    var g = List[List[Int]]()
    for _ in range(n + 1):
        h.append(List[Int](length=m + 1, fill=FAR))
        d.append(List[Int](length=m + 1, fill=FAR))
        g.append(List[Int](length=m + 1, fill=FAR))
    h[0][0] = 0
    for i in range(1, n + 1):
        d[i][0] = opening + i * extension
        h[i][0] = d[i][0]
    for j in range(1, m + 1):
        g[0][j] = opening + j * extension
        h[0][j] = g[0][j]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d[i][j] = min(h[i - 1][j] + opening + extension, d[i - 1][j] + extension)
            g[i][j] = min(h[i][j - 1] + opening + extension, g[i][j - 1] + extension)
            var diagonal = h[i - 1][j - 1] + (0 if a[i - 1] == b[j - 1] else mismatch)
            h[i][j] = min(diagonal, min(d[i][j], g[i][j]))
    return h[n][m]


def test_dear_gaps_cost_what_the_matrix_says() raises:
    """Gap costs far dearer than a mismatch leave most costs reachable by no path, which a search must
    not pay for: the distance and the score of pairs up to 40 letters, unrelated and close, empty
    sides among them, are the full matrix's, at an opening of 100,000 and at ordinary costs alike."""
    seed(71)
    for costs in [(4, 100_000, 10_000), (3, 50_000, 1), (1, 2, 1)]:
        var mismatch = costs[0]
        var opening = costs[1]
        var extension = costs[2]
        var scoring = Scoring.uniform(0, -mismatch, -opening, -extension)
        for trial in range(25):
            var first = random_sequence(0, 40, DNA_ALPHABET)
            var second = random_sequence(0, 40, DNA_ALPHABET) if trial % 2 == 0 else mutated(first, 0.15, 4)
            var expected = gotoh_cost(first, second, mismatch, opening, extension)
            assert_equal(distance(first, second, Costs.affine(mismatch, opening, extension)), expected)
            assert_equal(score(first, second, scoring, GLOBAL), -expected)


def test_lane_batches_match_single_pairs() raises:
    """A batch whose pairs go many at once into the lanes of a register (see `lanes`) gives each the
    cost a call of its own gives, under a band, under a cap, both and neither: unit, affine, linear and
    two-piece costs and deletions priced apart, pairs of every length up to a few hundred, close, divergent and
    unrelated, empty sides, and one too long for 16 bits among them, which goes the other way."""
    seed(29)
    var references = List[String]()
    var queries = List[String]()
    for trial in range(150):
        var reference = random_sequence(0, 300, DNA_ALPHABET)
        references.append(reference)
        if trial % 5 == 0:
            queries.append(random_sequence(0, 300, DNA_ALPHABET))
        else:
            queries.append(mutated(reference, [0.0, 0.01, 0.05, 0.2][trial % 4], 12))
    references.append(random_sequence(9000, 9000, DNA_ALPHABET))
    queries.append(mutated(references[len(references) - 1], 0.02, 6))
    var all_costs: List[Costs] = [
        Costs.edit(),
        Costs.affine(4, 6, 2),
        Costs.affine(1, 2, 1).with_deletions(3, 2),
        Costs.linear(2, 3),
        Costs.two_piece(4, 6, 2, 24, 1),
    ]
    var bands: List[Band] = [Band(), Band.around(6), Band(-3, 40)]
    for costs in all_costs:
        var uncapped = distances(references, queries, costs, threads=3)
        for index in range(len(references)):
            assert_equal(uncapped[index], distance(references[index], queries[index], costs))
        # A cap past every cost here leaves a pair no alignment inside its band fits as None, where no
        # cap would raise.
        for band in bands:
            for cap in [1 << 40, 8, 60]:
                var found = distances(references, queries, costs, max_cost=cap, band=band, threads=3)
                for index in range(len(references)):
                    var expected = distance(references[index], queries[index], costs, max_cost=cap, band=band)
                    assert_equal(found[index].or_else(-1), expected.or_else(-1))


def test_lane_free_ends_match_single_pairs() raises:
    """A batch of distances with free ends, which goes many pairs at once into the lanes with each lane's
    own free letters (see `lanes.LaneEnds`), gives each pair the cost a call of its own gives: reads in
    windows, prefixes and suffixes, a read holding the reference, a few letters free at each end and an
    overlap's, under unit, affine and two-piece costs, a band and a cap, empty sides among the pairs."""
    seed(47)
    var references = List[String]()
    var queries = List[String]()
    for trial in range(120):
        var reference = random_sequence(0, 260, DNA_ALPHABET)
        references.append(reference)
        if trial % 6 == 0:
            queries.append(random_sequence(0, 200, DNA_ALPHABET))
        else:
            var bytes = reference.as_bytes()
            var start = Int(random_ui64(0, UInt64(len(bytes) // 4)))
            var stop = len(bytes) - Int(random_ui64(0, UInt64(len(bytes) // 4)))
            var piece = String(StringSlice(unsafe_from_utf8=bytes[start : max(start, stop)]))
            queries.append(mutated(piece, [0.0, 0.02, 0.1][trial % 3], 6))
    var modes: List[Mode] = [
        Mode.INFIX,
        Mode.PREFIX,
        Mode.SUFFIX,
        Mode.REFERENCE_IN_QUERY,
        Mode.ends_free(reference_start=8, reference_end=8, query_start=3, query_end=3),
        Mode.ends_free(reference_end=40, query_start=40),
    ]
    var all_costs: List[Costs] = [Costs.edit(), Costs.affine(4, 6, 2), Costs.two_piece(4, 6, 2, 24, 1)]
    for costs in all_costs:
        for mode in modes:
            var found = distances(references, queries, costs, mode, threads=3)
            for index in range(len(references)):
                assert_equal(found[index], distance(references[index], queries[index], costs, mode))
            for band in [Band(), Band.around(20)]:
                var capped = distances(references, queries, costs, mode, max_cost=30, band=band, threads=3)
                for index in range(len(references)):
                    var expected = distance(references[index], queries[index], costs, mode, max_cost=30, band=band)
                    assert_equal(capped[index].or_else(-1), expected.or_else(-1))


def test_lane_free_alignments_match_single_pairs() raises:
    """A batch of alignments with free ends, which goes many pairs at once into the lanes, the end, then
    the start, then the span's own alignment (see `lanes.lane_free_alignments`), gives each pair the
    alignment a call of its own gives, its span and its CIGAR the ones either tie rule picks: reads in
    windows, prefixes and suffixes, a read holding the reference, a few letters free at each end and an
    overlap's, under unit, affine, linear and two-piece costs, with a cap and too little memory for a
    group's flags, empty sides and unrelated pairs among them."""
    seed(53)
    var references = List[String]()
    var queries = List[String]()
    for trial in range(120):
        var reference = random_sequence(0, 260, DNA_ALPHABET)
        references.append(reference)
        if trial % 6 == 0:
            queries.append(random_sequence(0, 200, DNA_ALPHABET))
        else:
            var bytes = reference.as_bytes()
            var start = Int(random_ui64(0, UInt64(len(bytes) // 4)))
            var stop = len(bytes) - Int(random_ui64(0, UInt64(len(bytes) // 4)))
            var piece = String(StringSlice(unsafe_from_utf8=bytes[start : max(start, stop)]))
            queries.append(mutated(piece, [0.0, 0.02, 0.1][trial % 3], 6))
    var modes: List[Mode] = [
        Mode.INFIX,
        Mode.PREFIX,
        Mode.SUFFIX,
        Mode.REFERENCE_IN_QUERY,
        Mode.ends_free(reference_start=8, reference_end=8, query_start=3, query_end=3),
        Mode.ends_free(reference_end=40, query_start=40),
    ]
    var all_costs: List[Costs] = [
        Costs.edit(),
        Costs.affine(4, 6, 2),
        Costs.linear(2, 3),
        Costs.two_piece(4, 6, 2, 24, 1),
    ]
    for costs in all_costs:
        for mode in modes:
            for ties in [Ties.LEFT, Ties.RIGHT]:
                for memory in [DEFAULT_MAX_MEMORY, 20000]:
                    var whole = alignments(references, queries, costs, mode, ties=ties, threads=3, max_memory=memory)
                    for index in range(len(references)):
                        var single = align(references[index], queries[index], costs, mode, ties=ties, max_memory=memory)
                        assert_equal(whole[index].cost, single.cost)
                        assert_equal(whole[index].cigar, single.cigar)
                        assert_equal(whole[index].reference_start, single.reference_start)
                        assert_equal(whole[index].reference_end, single.reference_end)
                        assert_equal(whole[index].query_start, single.query_start)
                        assert_equal(whole[index].query_end, single.query_end)
                var capped = alignments(references, queries, costs, mode, max_cost=30, ties=ties, threads=3)
                for index in range(len(references)):
                    var expected = align(references[index], queries[index], costs, mode, max_cost=30, ties=ties)
                    assert_equal(Bool(capped[index]), Bool(expected))
                    if expected:
                        assert_equal(capped[index].value().cost, expected.value().cost)
                        assert_equal(capped[index].value().cigar, expected.value().cigar)
                        assert_equal(capped[index].value().reference_start, expected.value().reference_start)
                        assert_equal(capped[index].value().query_start, expected.value().query_start)


def test_scoring_batches_match_single_pairs_in_every_mode() raises:
    """A batch of scores under a `Scoring`, with free ends or as an extension, spread over the threads asked
    for, gives each pair the score a call of its own gives; and a mode a single pair refuses, a local
    alignment with its own match score, the batch refuses too."""
    seed(59)
    var references = List[String]()
    var queries = List[String]()
    for trial in range(60):
        var reference = random_sequence(0, 220, DNA_ALPHABET)
        references.append(reference)
        queries.append(
            mutated(reference, [0.0, 0.05, 0.2][trial % 3], 6) if trial % 5
            != 0 else random_sequence(0, 150, DNA_ALPHABET)
        )
    var scoring = Scoring.dna()
    var modes: List[Mode] = [
        Mode.INFIX,
        Mode.PREFIX,
        Mode.ends_free(reference_end=40, query_start=40),
        Mode.extension(0),
        Mode.extension(0, Anchor.END, zdrop=20),
    ]
    for mode in modes:
        var found = scores(references, queries, scoring, mode, placement=Placement.on_cpu(3))
        for index in range(len(references)):
            assert_equal(found[index], score(references[index], queries[index], scoring, mode))
    var local = Mode.local(2)
    with assert_raises():
        _ = score(references[1], queries[1], scoring, local)
    with assert_raises():
        _ = scores(references, queries, scoring, local)


def test_unit_suffixes_take_the_bits_as_the_wavefront_aligns_them() raises:
    """A suffix at unit costs, which the bit-parallel search serves as a prefix of both sequences reversed,
    gives the cost, the span and the CIGAR the wavefront gives, which a cap sends it to, under either tie
    rule, one pair at a time and in a batch."""
    seed(67)
    var references = List[String]()
    var queries = List[String]()
    for trial in range(80):
        var reference = random_sequence(1, 400, DNA_ALPHABET)
        references.append(reference)
        var bytes = reference.as_bytes()
        var start = Int(random_ui64(0, UInt64(len(bytes))))
        var tail = String(StringSlice(unsafe_from_utf8=bytes[start:]))
        var query = mutated(tail, [0.0, 0.03, 0.15][trial % 3], 4) if trial % 7 != 0 else random_sequence(
            1, 80, DNA_ALPHABET
        )
        if query.byte_length() == 0:
            query = "ACGT"
        queries.append(query)
    var unit = Costs.edit()
    for ties in [Ties.LEFT, Ties.RIGHT]:
        var batch = alignments(references, queries, unit, Mode.SUFFIX, ties=ties, threads=3)
        var costs = distances(references, queries, unit, Mode.SUFFIX, threads=3)
        for index in range(len(references)):
            var swept = align(references[index], queries[index], unit, Mode.SUFFIX, ties=ties)
            var capped = align(references[index], queries[index], unit, Mode.SUFFIX, max_cost=1 << 40, ties=ties)
            var searched = capped.take()
            assert_equal(swept.cost, searched.cost)
            assert_equal(swept.cigar, searched.cigar)
            assert_equal(swept.reference_start, searched.reference_start)
            assert_equal(swept.reference_end, references[index].byte_length())
            assert_equal(distance(references[index], queries[index], unit, Mode.SUFFIX), searched.cost)
            assert_equal(batch[index].cigar, searched.cigar)
            assert_equal(batch[index].reference_start, searched.reference_start)
            assert_equal(costs[index], searched.cost)


def test_long_prefixes_find_their_distance() raises:
    """A read of a few kbp at up to a fifth divergence placed at a reference's start costs at unit costs
    what the wavefront finds, its CIGAR spending it: the bit-parallel prefix search gives a try up at a
    checkpoint once its climb projects past the bound (see `edit_search.banded_last_row`), and a later
    try must still find the least cost."""
    seed(61)
    for trial in range(24):
        var length = [1000, 2000, 3000][trial % 3]
        var source = random_sequence(length + 500, length + 500, DNA_ALPHABET)
        var read = mutated(
            String(StringSlice(unsafe_from_utf8=source.as_bytes()[0:length])), [0.05, 0.1, 0.2][trial % 3], 6
        )
        var swept = distance(source, read, Costs.edit(), Mode.PREFIX)
        # A cap the batch never reaches sends the pair to the wavefront instead (see `api.swept_by_bits`).
        var searched = distance(source, read, Costs.edit(), Mode.PREFIX, max_cost=1 << 40)
        assert_equal(swept, searched.value())
        var aligned = align(source, read, Costs.edit(), Mode.PREFIX)
        assert_equal(aligned.cost, swept)
        assert_equal(aligned.reference_start, 0)


def test_aligner_matches_the_functions() raises:
    """One `Aligner` reused across every call, its memory kept from pair to pair, gives each pair what
    the functions give: distances and alignments under unit, affine and two-piece costs, globally and
    with free ends, under either tie rule, a cap and a band, short pairs after long and empty sides."""
    seed(59)
    var aligner = Aligner()
    var all_costs: List[Costs] = [Costs.edit(), Costs.affine(4, 6, 2), Costs.two_piece(4, 6, 2, 24, 1)]
    var modes: List[Mode] = [Mode.GLOBAL, Mode.INFIX, Mode.PREFIX]
    for trial in range(120):
        var reference = random_sequence(0, [40, 600, 12, 250][trial % 4], DNA_ALPHABET)
        var query = mutated(reference, [0.0, 0.02, 0.1][trial % 3], 8) if trial % 7 else random_sequence(
            0, 300, DNA_ALPHABET
        )
        var costs = all_costs[trial % 3]
        var mode = modes[(trial // 3) % 3]
        var ties = Ties.LEFT if trial % 2 else Ties.RIGHT
        assert_equal(aligner.distance(reference, query, costs, mode), distance(reference, query, costs, mode))
        var mine = aligner.align(reference, query, costs, mode, ties=ties)
        var theirs = align(reference, query, costs, mode, ties=ties)
        assert_equal(mine.cost, theirs.cost)
        assert_equal(mine.cigar, theirs.cigar)
        assert_equal(mine.reference_start, theirs.reference_start)
        var capped = aligner.align(reference, query, costs, mode, max_cost=40, band=Band.around(30), ties=ties)
        var expected = align(reference, query, costs, mode, max_cost=40, band=Band.around(30), ties=ties)
        assert_equal(Bool(capped), Bool(expected))
        if expected:
            assert_equal(capped.value().cigar, expected.value().cigar)
        assert_equal(
            aligner.distance(reference, query, costs, mode, max_cost=40).or_else(-1),
            distance(reference, query, costs, mode, max_cost=40).or_else(-1),
        )


def test_lane_tables_match_single_pairs() raises:
    """A batch of global or local scores under a table of more than one mismatch score, small enough for
    a register, which goes many pairs at once into the lanes over the alphabet's codes (see
    `scoring.tabled_scores`), gives each pair the score a call of its own gives: a transition and a
    transversion table, close and unrelated pairs, empty sides, and a letter outside the alphabet raising
    as it does alone."""
    seed(53)
    var cells = List[Int8]()
    for row in range(4):
        for column in range(4):
            cells.append(Int8(2 if row == column else (-2 if (row + column) % 2 == 0 else -4)))
    var tables: List[Scoring] = [
        Scoring.tabulated("ACGT", cells.copy(), -4, -2),
        Scoring.tabulated("ACGT", cells.copy(), -10, -1),
    ]
    var firsts = List[String]()
    var seconds = List[String]()
    for trial in range(120):
        var first = random_sequence(0, 300, DNA_ALPHABET)
        firsts.append(first)
        seconds.append(random_sequence(0, 300, DNA_ALPHABET) if trial % 6 == 0 else mutated(first, 0.05, 8))
    for scoring in tables:
        for mode in [GLOBAL, LOCAL]:
            var found = scores(firsts, seconds, scoring, mode, placement=Placement.on_cpu(3))
            for index in range(len(firsts)):
                assert_equal(found[index], score(firsts[index], seconds[index], scoring, mode))
    var odd_firsts: List[String] = ["ACGT", "ACNT"]
    var odd_seconds: List[String] = ["ACGT", "ACGT"]
    var raised = False
    try:
        _ = scores(odd_firsts, odd_seconds, tables[0], GLOBAL, placement=Placement.on_cpu(2))
    except:
        raised = True
    assert_true(raised)


def test_lane_alignments_match_single_pairs() raises:
    """A batch of global alignments, which goes many pairs at once into the lanes and traces each from the
    flags its band kept (see `lanes.traced`), gives each pair the alignment a call of its own gives, its
    CIGAR the one either tie rule picks: under a band, under a cap, both and neither, and with too little
    memory for a group's flags, which leaves its pairs to the searches; unit, affine, linear and
    two-piece costs and deletions priced apart, pairs of every length up to a few hundred, close,
    divergent and unrelated, empty sides, and one long pair among them."""
    seed(43)
    var references = List[String]()
    var queries = List[String]()
    for trial in range(150):
        var reference = random_sequence(0, 300, DNA_ALPHABET)
        references.append(reference)
        if trial % 5 == 0:
            queries.append(random_sequence(0, 300, DNA_ALPHABET))
        else:
            queries.append(mutated(reference, [0.0, 0.01, 0.05, 0.2][trial % 4], 12))
    references.append(random_sequence(3000, 3000, DNA_ALPHABET))
    queries.append(mutated(references[len(references) - 1], 0.02, 6))
    var all_costs: List[Costs] = [
        Costs.edit(),
        Costs.affine(4, 6, 2),
        Costs.affine(1, 2, 1).with_deletions(3, 2),
        Costs.linear(2, 3),
        Costs.two_piece(4, 6, 2, 24, 1),
    ]
    var bands: List[Band] = [Band(), Band.around(6), Band(-3, 40)]
    for costs in all_costs:
        for ties in [Ties.LEFT, Ties.RIGHT]:
            for memory in [DEFAULT_MAX_MEMORY, 200000, 20000]:
                var whole = alignments(references, queries, costs, ties=ties, threads=3, max_memory=memory)
                for index in range(len(references)):
                    var single = align(references[index], queries[index], costs, ties=ties, max_memory=memory)
                    assert_equal(whole[index].cost, single.cost)
                    assert_equal(whole[index].cigar, single.cigar)
            for band in bands:
                for cap in [1 << 40, 8, 60]:
                    var found = alignments(references, queries, costs, max_cost=cap, band=band, ties=ties, threads=3)
                    for index in range(len(references)):
                        var expected = align(
                            references[index], queries[index], costs, max_cost=cap, band=band, ties=ties
                        )
                        assert_equal(Bool(found[index]), Bool(expected))
                        if expected:
                            assert_equal(found[index].value().cost, expected.value().cost)
                            assert_equal(found[index].value().cigar, expected.value().cigar)


def test_lane_scores_match_single_pairs() raises:
    """A batch of global scores or alignments under a table of one match and one mismatch score, which
    goes many pairs at once into the lanes (see `lanes`), gives each pair the score and the CIGAR a call
    of its own gives: minimap2's
    scores, a dear gap and a cheap one, unit costs, pairs of every length up to a few hundred, close and
    unrelated, empty sides among them; and a letter outside the alphabet raises as it does alone."""
    seed(37)
    var firsts = List[String]()
    var seconds = List[String]()
    for trial in range(120):
        var first = random_sequence(0, 300, DNA_ALPHABET)
        firsts.append(first)
        seconds.append(random_sequence(0, 300, DNA_ALPHABET) if trial % 6 == 0 else mutated(first, 0.05, 8))
    var all_scoring: List[Scoring] = [
        Scoring.dna(),
        Scoring.uniform(1, -3, -12, -1),
        Scoring.uniform(2, -2, -2, -1),
        Scoring.edit_distance(),
    ]
    for scoring in all_scoring:
        var found = scores(firsts, seconds, scoring, GLOBAL, placement=Placement.on_cpu(3))
        for index in range(len(firsts)):
            assert_equal(found[index], score(firsts[index], seconds[index], scoring, GLOBAL))
        # Local scores too, the whole matrix swept many pairs at once.
        var local = scores(firsts, seconds, scoring, LOCAL, placement=Placement.on_cpu(3))
        for index in range(len(firsts)):
            assert_equal(local[index], score(firsts[index], seconds[index], scoring, LOCAL))
        # Their alignments too, each traced from its lanes' flags as the wavefront's tie rule picks it.
        var aligned = alignments(firsts, seconds, scoring, GLOBAL, placement=Placement.on_cpu(3))
        for index in range(len(firsts)):
            var single = align(firsts[index], seconds[index], scoring, GLOBAL, placement=Placement.on_cpu(1))
            assert_equal(aligned[index].score, single.score)
            assert_equal(aligned[index].cigar, single.cigar)
    var odd_firsts: List[String] = ["ACGT", "ACXT"]
    var odd_seconds: List[String] = ["ACGT", "ACGT"]
    var dna = Scoring.dna()
    var raised = False
    try:
        _ = scores(odd_firsts, odd_seconds, dna, GLOBAL, placement=Placement.on_cpu(2))
    except:
        raised = True
    assert_true(raised)


def test_batches_reuse_their_searches_cleanly() raises:
    """A batch's worker keeps its searches from pair to pair (see `gap_affine.SearchSpace`), so a pair
    must never see the last one's state: long pairs after short and short after long, identical and
    unrelated, an empty side, under one and two gap pieces and deletions priced apart, globally and with
    free ends, each pair's batch cost the one a call of its own gives, and each alignment the one a call
    of its own gives under either tie rule, with its kept fronts whole or split for want of memory."""
    seed(91)
    var references = List[String]()
    var queries = List[String]()
    for trial in range(90):
        var length = [0, 3, 40, 600, 12, 250][trial % 6]
        var reference = random_sequence(length, length, DNA_ALPHABET)
        references.append(reference)
        if trial % 7 == 0:
            queries.append(random_sequence(0, 300, DNA_ALPHABET))
        else:
            queries.append(mutated(reference, [0.0, 0.01, 0.1, 0.4][trial % 4], 9))
    var all_costs: List[Costs] = [
        Costs.affine(4, 6, 2),
        Costs.two_piece(4, 6, 2, 24, 1),
        Costs.affine(1, 2, 1).with_deletions(3, 2),
    ]
    for costs in all_costs:
        for mode in [Mode.GLOBAL, Mode.INFIX]:
            for threads in [1, 4]:
                var found = distances(references, queries, costs, mode, threads=threads)
                for index in range(len(references)):
                    assert_equal(found[index], distance(references[index], queries[index], costs, mode))
                for ties in [Ties.LEFT, Ties.RIGHT]:
                    for memory in [DEFAULT_MAX_MEMORY, 3000]:
                        var aligned = alignments(
                            references, queries, costs, mode, ties=ties, threads=threads, max_memory=memory
                        )
                        for index in range(len(references)):
                            var single = align(
                                references[index], queries[index], costs, mode, ties=ties, max_memory=memory
                            )
                            assert_equal(aligned[index].cost, single.cost)
                            assert_equal(aligned[index].cigar, single.cigar)
                            assert_equal(aligned[index].reference_start, single.reference_start)


def test_capped_batches_match_single_pairs() raises:
    """A batch under a cost cap is each pair's capped `distance` and `align`, None past the cap, on any
    number of threads; with no cap a pair outside the band raises."""
    seed(89)
    var references = List[String]()
    var queries = List[String]()
    for trial in range(60):
        var reference = random_sequence(0, 200, DNA_ALPHABET)
        references.append(reference)
        queries.append(mutated(reference, [0.0, 0.05, 0.3][trial % 3], 6))
    var costs = Costs.affine(4, 6, 2)
    for threads in [1, 3]:
        var found = distances(references, queries, costs, max_cost=40, threads=threads)
        var aligned = alignments(references, queries, costs, Mode.INFIX, max_cost=40, threads=threads)
        for index in range(len(references)):
            var expected = distance(references[index], queries[index], costs, max_cost=40)
            assert_equal(Bool(found[index]), Bool(expected))
            if expected:
                assert_equal(found[index].value(), expected.value())
            var single = align(references[index], queries[index], costs, Mode.INFIX, max_cost=40)
            assert_equal(Bool(aligned[index]), Bool(single))
            if single:
                assert_equal(aligned[index].value().cigar, single.value().cigar)
                assert_equal(aligned[index].value().reference_start, single.value().reference_start)
    var outside: List[String] = ["ACGT", "ACGTACGT"]
    var short: List[String] = ["ACGT", "AC"]
    with assert_raises():
        _ = distances(outside, short, costs, band=Band.around(1))
    var capped = distances(outside, short, costs, max_cost=100, band=Band.around(1))
    assert_equal(capped[0].value(), 0)
    assert_false(Bool(capped[1]))


def test_memory_budget_keeps_the_cost() raises:
    """Under a budget too small to keep any front, every mode splits its pairs again and again and
    still finds an alignment of the least cost, or the best score; the default budget keeps them whole."""
    seed(97)
    var costs = Costs.affine(4, 6, 2)
    var modes: List[Mode] = [Mode.GLOBAL, Mode.INFIX, Mode.SUFFIX, Mode.extension(2), Mode.local(2), Mode.overlap(2)]
    for trial in range(12):
        var reference = random_sequence(300, 1500, DNA_ALPHABET)
        var query = mutated(reference, [0.02, 0.1, 0.3][trial % 3], 8)
        for mode in modes:
            var whole_pair = align(reference, query, costs, mode)
            for budget in [0, 5000]:
                var split = align(reference, query, costs, mode, max_memory=budget)
                assert_equal(
                    split.cost if not mode.is_scored() else split.score,
                    whole_pair.cost if not mode.is_scored() else whole_pair.score,
                )
                var part = String(
                    StringSlice(unsafe_from_utf8=reference.as_bytes()[split.reference_start : split.reference_end])
                )
                var piece = String(StringSlice(unsafe_from_utf8=query.as_bytes()[split.query_start : split.query_end]))
                _ = rows_from_cigar(part, piece, split.cigar)


def test_zdrop_gives_up_on_noise() raises:
    """An extension with a Z-drop stops where a read turns to noise, while the exact one crosses the
    noise to the matching stretch past it; a Z-drop larger than any fall is the exact extension. Every
    stop is an alignment that earns its score."""
    seed(101)
    var costs = Costs.affine(4, 6, 2)
    var core = random_sequence(200, 200, DNA_ALPHABET)
    var tail = random_sequence(300, 300, DNA_ALPHABET)
    var reference = core + random_sequence(60, 60, DNA_ALPHABET) + tail
    var query = core + random_sequence(60, 60, DNA_ALPHABET) + tail
    var exact = align(reference, query, costs, Mode.extension(1))
    assert_true(exact.reference_end > 500)
    var dropped = align(reference, query, costs, Mode.extension(1, zdrop=50))
    assert_equal(dropped.reference_end, 200)
    assert_equal(dropped.score, 200)
    assert_equal(dropped.cigar, "200=")
    var back = align(reversed_text(reference), reversed_text(query), costs, Mode.extension(1, Anchor.END, zdrop=50))
    assert_equal(back.reference_start, reference.byte_length() - 200)
    for trial in range(40):
        var a = random_sequence(0, 400, DNA_ALPHABET)
        var b = mutated(a, [0.02, 0.1, 0.3, 0.6][trial % 4], 6) + random_sequence(0, 50, DNA_ALPHABET)
        var whole_way = align(a, b, costs, Mode.extension(2))
        var generous = align(a, b, costs, Mode.extension(2, zdrop=1 << 30))
        assert_equal(generous.score, whole_way.score)
        assert_equal(generous.cigar, whole_way.cigar)
        for drop in [0, 10, 100]:
            var found = align(a, b, costs, Mode.extension(2, zdrop=drop))
            assert_true(found.score <= whole_way.score)
            assert_equal(extension_price(found.cigar, 2, 4, 6, 2, -1, 0), found.score)
    with assert_raises():
        _ = Mode.extension(1, zdrop=-5)


def column_maxima(first: String, second: String, a: Int, x: Int, o: Int, e: Int, o2: Int, e2: Int) -> List[Int]:
    """Each reference column's best Smith-Waterman score, of every cell with that many reference
    letters, over the whole matrix; shares no code with the library."""
    comptime LOW = -(1 << 40)
    var p = first.as_bytes()
    var q = second.as_bytes()
    var n = len(p)
    var m = len(q)
    var width = m + 1
    var best = List[Int](length=(n + 1) * width, fill=0)
    var layers = List[List[Int]]()
    for _ in range(4):
        layers.append(List[Int](length=(n + 1) * width, fill=LOW))
    var out = List[Int](length=n + 1, fill=0)
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            var at = i * width + j
            var value = 0
            for piece in range(2 if o2 >= 0 else 1):
                var opening = o if piece == 0 else o2
                var extension = e if piece == 0 else e2
                layers[2 * piece][at] = max(
                    best[at - width] - opening - extension, layers[2 * piece][at - width] - extension
                )
                layers[2 * piece + 1][at] = max(
                    best[at - 1] - opening - extension, layers[2 * piece + 1][at - 1] - extension
                )
                value = max(value, max(layers[2 * piece][at], layers[2 * piece + 1][at]))
            value = max(value, best[at - width - 1] + (a if p[i - 1] == q[j - 1] else -x))
            best[at] = value
            out[i] = max(out[i], value)
    return out^


def test_scores_without_alignments() raises:
    """`score` is `align`'s score for every mode, and `local_scores` the local optimum with SSW's second
    best: the best column more than the window from the best end, the first on a tie."""
    seed(103)
    var costs = Costs.affine(4, 6, 2)
    var modes: List[Mode] = [
        Mode.GLOBAL,
        Mode.INFIX,
        Mode.extension(1),
        Mode.extension(2, Anchor.END),
        Mode.local(2),
        Mode.overlap(2),
        Mode.INFIX.with_match_score(2),
        Mode.GLOBAL.with_match_score(1),
    ]
    for trial in range(40):
        var core = random_sequence(1, 200, DNA_ALPHABET)
        var reference = random_sequence(0, 100, DNA_ALPHABET) + core + random_sequence(0, 100, DNA_ALPHABET)
        var query = mutated(core, [0.0, 0.05, 0.2][trial % 3], 6)
        if trial % 4 == 0:
            # A second copy elsewhere, for a second best worth its name.
            reference += random_sequence(30, 60, DNA_ALPHABET) + mutated(core, 0.1, 4)
        for mode in modes:
            assert_equal(score(reference, query, costs, mode), align(reference, query, costs, mode).score)
        var found = local_scores(reference, query, costs, Mode.local(2))
        var maxima = column_maxima(reference, query, 2, 4, 6, 2, -1, 0)
        assert_equal(found.score, local_optimum(reference, query, 2, 4, 6, 2, -1, 0))
        var window = max(query.byte_length() // 2, 15)
        var second = 0
        var second_end = 0
        for letters in range(len(maxima)):
            if abs(letters - found.reference_end) > window and maxima[letters] > second:
                second = maxima[letters]
                second_end = letters
        assert_equal(found.second_score, second)
        assert_equal(found.second_reference_end, second_end)
        var aligned = align(reference, query, costs, Mode.local(2), ties=Ties.LEFT)
        assert_equal(found.reference_end, aligned.reference_end)
        assert_equal(found.query_end, aligned.query_end)
        # The shorter sequence as the reference: the lanes run along it, the columns kept by row.
        var swapped = local_scores(query, reference, costs, Mode.local(2), window=10)
        var swapped_maxima = column_maxima(query, reference, 2, 4, 6, 2, -1, 0)
        var swapped_second = 0
        for letters in range(len(swapped_maxima)):
            if abs(letters - swapped.reference_end) > 10:
                swapped_second = max(swapped_second, swapped_maxima[letters])
        assert_equal(swapped.second_score, swapped_second)
    with assert_raises():
        _ = local_scores("ACGT", "ACGT", costs, Mode.GLOBAL)


def test_deletions_cost_their_own() raises:
    """Deletions and insertions priced apart, as bwa's `-O del,ins`: the dearer side is avoided where
    the other serves, every mode prices its CIGAR by side, and costs equal both ways are the costs
    they always were."""
    var dear_deletions = Costs.affine(4, 2, 1).with_deletions(20, 5)
    assert_equal(align("ACGTACGT", "ACGACGT", dear_deletions).cost, 25)
    assert_equal(align("ACGACGT", "ACGTACGT", dear_deletions).cigar, "3=1I4=")
    assert_equal(distance("ACGACGT", "ACGTACGT", dear_deletions), 3)
    assert_true(dear_deletions.symmetric() == False)
    assert_true(Costs.affine(4, 6, 2).with_deletions(6, 2).symmetric())
    assert_equal(dear_deletions.unit_scale(), 0)
    # Two pieces on one side alone: the other counts its one piece twice.
    var long_deletions = Costs.affine(4, 6, 2).with_deletions(6, 2, 24, 1)
    assert_equal(long_deletions.pieces(), 2)
    assert_equal(long_deletions.gap(30, True), 54)
    assert_equal(long_deletions.gap(30, False), 66)
    var gapped = "GATTACAGCTTGCA" + "C" * 30 + "TGGACCATGAGTCA"
    var plain = "GATTACAGCTTGCA" + "TGGACCATGAGTCA"
    assert_equal(distance(gapped, plain, long_deletions), 54)
    assert_equal(distance(plain, gapped, long_deletions), 66)
    with assert_raises():
        _ = Costs.affine(4, 6, 2).with_deletions(6, 0)


def test_search_ranks_every_reference() raises:
    """A search's every hit scores what `score` gives its pair, the best first, ties by the
    references' order; `best` keeps that many, a cap drops what passes it, and `aligned` aligns the
    kept hits. Local scores come from groups of references one to a lane, of every length, a group
    shorter than a lane's width among them, in 16 bits and, where a reward passes them, in 32."""
    seed(107)
    var query = random_sequence(80, 160, DNA_ALPHABET)
    var references = List[String]()
    for index in range(77):
        var kind = index % 4
        if kind == 0:
            references.append(random_sequence(0, 300, DNA_ALPHABET))
        elif kind == 1:
            references.append(
                random_sequence(0, 100, DNA_ALPHABET) + mutated(query, 0.1, 4) + random_sequence(0, 100, DNA_ALPHABET)
            )
        elif kind == 2:
            references.append(mutated(query, 0.3, 6))
        else:
            references.append(random_sequence(0, 20, DNA_ALPHABET))
    for costs in [Costs.affine(4, 6, 2), Costs.two_piece(4, 6, 2, 24, 1), Costs.affine(2, 3, 1).with_deletions(9, 2)]:
        for mode in [Mode.local(2), Mode.INFIX, Mode.GLOBAL, Mode.overlap(1)]:
            var hits = search(references, query, costs, mode)
            assert_equal(len(hits), len(references))
            for rank in range(len(hits)):
                assert_equal(hits[rank].score, score(references[hits[rank].index], query, costs, mode))
                if rank > 0:
                    var before = hits[rank - 1].score
                    assert_true(
                        before > hits[rank].score
                        or (before == hits[rank].score and hits[rank - 1].index < hits[rank].index)
                    )
            var top = search(references, query, costs, mode, best=5, aligned=True, threads=3)
            assert_equal(len(top), 5)
            for rank in range(5):
                assert_equal(top[rank].index, hits[rank].index)
                assert_equal(top[rank].alignment.value().score, top[rank].score)
    var capped = search(references, query, Costs.affine(4, 6, 2), Mode.INFIX, max_cost=60)
    for hit in capped:
        assert_true(-hit.score <= 60)
    var whole = search(references, query, Costs.affine(4, 6, 2), Mode.INFIX)
    var within = 0
    for hit in whole:
        if -hit.score <= 60:
            within += 1
    assert_equal(len(capped), within)
    assert_equal(len(search(List[String](), query, Costs.affine(4, 6, 2), Mode.local(2))), 0)
    # A reward 16 bits cannot hold over the query, which the lanes score in 32.
    var wide = Mode.local(300)
    for hit in search(references, query, Costs.affine(400, 600, 200), wide, threads=2):
        assert_equal(hit.score, score(references[hit.index], query, Costs.affine(400, 600, 200), wide))


def rebuilt_reference(query: String, cigar: String, md: String) raises -> String:
    """The reference's aligned part from the query's, the CIGAR and the `MD` string, as SAM readers
    rebuild it: the query's letters through matches and substitutions, those `MD` names replaced, and
    `MD`'s deleted letters spliced in."""
    # The `MD` string as a queue of events: matched letters, a substituted letter, or deleted letters.
    var counts = List[Int]()
    var letters = List[String]()
    var deleted = List[Bool]()
    var bytes = md.as_bytes()
    var index = 0
    while index < len(bytes):
        var number = 0
        while index < len(bytes) and bytes[index] >= UInt8(ord("0")) and bytes[index] <= UInt8(ord("9")):
            number = number * 10 + Int(bytes[index] - UInt8(ord("0")))
            index += 1
        counts.append(number)
        if index == len(bytes):
            letters.append(String())
            deleted.append(False)
            break
        var gap = bytes[index] == UInt8(ord("^"))
        if gap:
            index += 1
        var start = index
        while index < len(bytes) and (bytes[index] < UInt8(ord("0")) or bytes[index] > UInt8(ord("9"))):
            index += 1
            if not gap:
                break
        letters.append(String(StringSlice(unsafe_from_utf8=bytes[start:index])))
        deleted.append(gap)
    var out = String()
    var row = 0
    var event = 0
    var left = counts[0]
    var operations = cigar_runs(cigar)
    for slot in range(len(operations[0])):
        var operation = operations[0][slot]
        var length = operations[1][slot]
        if operation == UInt8(ord("I")):
            row += length
        elif operation == UInt8(ord("D")):
            assert_true(deleted[event])
            assert_equal(letters[event].byte_length(), length)
            out += letters[event]
            event += 1
            left = counts[event]
        else:
            for _ in range(length):
                if left > 0:
                    out += String(StringSlice(unsafe_from_utf8=query.as_bytes()[row : row + 1]))
                    left -= 1
                else:
                    assert_false(deleted[event])
                    out += letters[event]
                    event += 1
                    left = counts[event]
                row += 1
    return out


# endregion Output

# region Device


def test_device_edit_distances() raises:
    """On the GPU, a batch's unit-cost distances are the host's, patterns of one word and of many,
    symbols past ACGT, empty sides, and a pair too long for a thread, which the host takes; other
    costs are refused there. Skipped where no accelerator answers."""
    if not gpu_available():
        print("    skipped: no accelerator serves a real alignment here")
        return
    seed(109)
    var references = List[String]()
    var queries = List[String]()
    for trial in range(400):
        var reference = random_sequence(0, [60, 300, 1500][trial % 3], "ACGTN" if trial % 5 == 0 else DNA_ALPHABET)
        references.append(reference)
        queries.append(mutated(reference, [0.0, 0.1, 0.4][trial % 3], 6))
    references.append(random_sequence(5000, 5000, DNA_ALPHABET))
    queries.append(mutated(references[len(references) - 1], 0.05, 4))
    var device = Placement.on_gpu(0, 4)
    var found = distances(references, queries, placement=device)
    var expected = distances(references, queries)
    for index in range(len(found)):
        assert_equal(found[index], expected[index])
    assert_equal(distances(references, queries, Costs.linear(2, 2), placement=device)[3], 2 * expected[3])
    var affine = Costs.affine(4, 6, 2)
    with assert_raises(contains="unit costs"):
        _ = distances(references, queries, affine, placement=device)


def test_device_grouped_scores_match_host() raises:
    """A batch scored several pairs a warp gives every pair the host's score, whichever shape the batch's
    widest second sequence picks (see `score_groups.SHAPES`), the widest past every shape scored a warp
    a pair: pairs of every length up to the width, empty sides among them, the rows longer or shorter
    than the columns, close pairs and unrelated ones. Skipped where no accelerator answers."""
    if not gpu_available():
        print("    skipped: no accelerator serves a real alignment here")
        return
    seed(23)
    var device = Placement.on_gpu(0, 4)
    var host = Placement.on_cpu(4)
    for scoring in scoring_regimes():
        for widest in [1, 31, 48, 64, 95, 150, 200, 256, 380, 512, 700]:
            var firsts = List[String]()
            var seconds = List[String]()
            for trial in range(70):
                var second = random_sequence(0 if trial % 9 == 0 else 1, widest, DNA_ALPHABET)
                if trial == 0:
                    second = random_sequence(widest, widest, DNA_ALPHABET)
                seconds.append(second)
                firsts.append(
                    random_sequence(0, widest + 40, DNA_ALPHABET) if trial % 4 == 0 else mutated(second, 0.08, 5)
                )
            comptime for mode in [GLOBAL, LOCAL]:
                var found = scores(firsts, seconds, scoring, mode, placement=device)
                for index in range(len(firsts)):
                    assert_equal(found[index], score(firsts[index], seconds[index], scoring, mode, placement=host))


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
            var on_device = alignments(firsts, seconds, scoring, mode, placement=device)
            var device_scores = scores(firsts, seconds, scoring, mode, placement=device)
            for index in range(len(firsts)):
                var expected = align(firsts[index], seconds[index], scoring, mode, placement=host)
                assert_equal(on_device[index].score, expected.score)
                assert_well_formed(mode, firsts[index], seconds[index], on_device[index], scoring)
                assert_well_formed(mode, firsts[index], seconds[index], expected, scoring)
                assert_equal(device_scores[index], expected.score)

    var scoring = expensive_gap()
    var tall = "A" * 70_000
    var short = "ACGT" * 4
    comptime for mode in [GLOBAL, LOCAL]:
        assert_equal(
            score(tall, short, scoring, mode, placement=device), score(tall, short, scoring, mode, placement=host)
        )
        assert_equal(score(tall, "", scoring, mode, placement=device), score(tall, "", scoring, mode, placement=host))
        assert_equal(score("", short, scoring, mode, placement=device), score("", short, scoring, mode, placement=host))
        var long_first = random_sequence(2000, 2000, DNA_ALPHABET)
        var long_second = random_sequence(2000, 2000, DNA_ALPHABET)
        var linear = align(long_first, long_second, scoring, mode, placement=device, max_memory=0)
        assert_well_formed(mode, long_first, long_second, linear, scoring)
        assert_equal(linear.score, score(long_first, long_second, scoring, mode, placement=host))


# endregion Device


def main() raises:
    """Runs every `test_` function here."""
    TestSuite.discover_tests[__functions_in_module()]().run()
