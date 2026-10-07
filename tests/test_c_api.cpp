// The C API through its C++ wrapper (see c/dinara.h), and once through C's own structs: known distances
// and CIGARs, empty sequences, symbols past ACGT, every mode, long pairs whose CIGAR spells the distance,
// four threads at once, and batches.
//
//     pixi run test-c
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <thread>
#include <vector>

#include "dinara.h"

static int failures = 0;

#define CHECK(condition)                                                       \
    do {                                                                       \
        if (!(condition)) {                                                    \
            std::fprintf(stderr, "%s:%d: failed: %s\n", __FILE__, __LINE__, #condition); \
            ++failures;                                                        \
        }                                                                      \
    } while (0)

// The substitutions and gaps a CIGAR counts, after checking it consumes both sequences exactly.
static int64_t edits(const std::string &cigar, size_t first_length, size_t second_length) {
    int64_t total = 0, length = 0;
    size_t columns = 0, rows = 0;
    for (char symbol : cigar) {
        if (symbol >= '0' && symbol <= '9') {
            length = length * 10 + (symbol - '0');
            continue;
        }
        if (symbol == 'D') columns += length, total += length;
        else if (symbol == 'I') rows += length, total += length;
        else columns += length, rows += length, total += symbol == 'X' ? length : 0;
        length = 0;
    }
    CHECK(columns == first_length && rows == second_length);
    return total;
}

static std::string mutated(const std::string &text, double rate, std::mt19937 &random) {
    static const char bases[] = "ACGT";
    std::uniform_real_distribution<double> roll(0.0, 1.0);
    std::uniform_int_distribution<int> base(0, 3);
    std::string out;
    for (char symbol : text) {
        double draw = roll(random);
        if (draw < rate / 3) out += bases[base(random)];
        else if (draw < 2 * rate / 3) continue;
        else if (draw < rate) out += symbol, out += bases[base(random)];
        else out += symbol;
    }
    return out;
}

int main() {
    using dinara::Costs;
    using dinara::Mode;
    CHECK(dinara::distance("ACGTACGTTTGCA", "ACGTCGTTTTGCA") == 2);
    dinara::Alignment aligned = dinara::align("ACGTACGTTTGCA", "ACGTCGTTTTGCA");
    CHECK(aligned.cost == 2 && aligned.score == -2 && aligned.cigar == "4=1D2=1I6=");
    CHECK(aligned.reference_start == 0 && aligned.reference_end == 13 && aligned.query_end == 13);
    CHECK(dinara::align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", Costs::edit(), Mode::global(), {}, dinara::Ties::left, false)
              .cigar == "4M1D2M1I6M");

    CHECK(dinara::align("", "ACG").cigar == "3I");
    CHECK(dinara::align("ACG", "").cigar == "3D");
    dinara::Alignment empty = dinara::align("", "");
    CHECK(empty.cost == 0 && empty.cigar.empty());

    CHECK(dinara::align("ACGTNNAC", "ACGTNNAC").cigar == "8=");
    CHECK(dinara::distance("ACGTN", "ACGTA") == 1);
    // More symbols than the bit-parallel sweep takes: the wavefront takes any byte.
    CHECK(dinara::distance("NRYKMSW", "NRYKMSA") == 1);
    bool refused = false;
    try {
        dinara::distance("\xfe", "A");
    } catch (const dinara::UnsupportedSymbols &) {
        refused = true;
    }
    CHECK(refused);

    Costs affine = Costs::affine(4, 6, 2);
    dinara::Alignment cost12 = dinara::align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", affine);
    CHECK(cost12.cost == 12 && cost12.cigar == "4=3X6=");
    CHECK(dinara::align("ACGT", "", affine).cigar == "4D");
    CHECK(dinara::align("ACGT", "", affine).cost == 14);
    bool invalid = false;
    try {
        dinara::align("A", "C", Costs::affine(0, 6, 2));
    } catch (const std::invalid_argument &) {
        invalid = true;
    }
    CHECK(invalid);
    CHECK(dinara::distance("ACGTACGTTTGCA", "ACGTCGTTTTGCA", affine) == 12);
    CHECK(dinara::distance_within("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 12, affine) == 12);
    CHECK(!dinara::distance_within("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 11, affine));
    CHECK(!dinara::align_within("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 11, affine));
    CHECK(dinara::align_within("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 12, affine)->cigar == "4=3X6=");
    // A read placed inside a reference: its span, and the CIGAR over it alone.
    dinara::Alignment placed = dinara::align("TTTTACGTACGTTTTT", "ACGTACGT", affine, Mode::infix());
    CHECK(placed.cost == 0 && placed.cigar == "8=" && placed.reference_start == 4 && placed.reference_end == 12);
    dinara::Alignment found = dinara::align("TTTTACGTACGTTTTT", "ACGTCGT", Costs::edit(), Mode::infix());
    CHECK(found.cost == 1 && found.reference_start == 4);
    CHECK(dinara::distance("TTTTACGTACGTTTTT", "ACGTACGT", affine) > 0);
    CHECK(dinara::distance("TTTTACGTACGTTTTT", "ACGTACGT", affine, Mode::infix()) == 0);
    CHECK(dinara::distance("ACGTACGTTTTT", "ACGTACGT", affine, Mode::prefix()) == 0);
    CHECK(dinara::distance("TTTTACGTACGT", "ACGTACGT", affine, Mode::suffix()) == 0);
    CHECK(dinara::distance("TTTTACGTACGT", "ACGTACGT", Costs::edit(), Mode::suffix()) == 0);
    dinara::Alignment overlap = dinara::align("TTTTTACGTACGT", "ACGTACGTGGGGG", affine, Mode::ends_free(5, 0, 0, 5));
    CHECK(overlap.cost == 0 && overlap.cigar == "8=" && overlap.reference_start == 5 && overlap.query_end == 8);
    // Two-piece gap costs: a long gap at the second piece, 24 + 30, and a short one at the first, 6 + 2.
    std::string gapped = "GATTACAGCTTGCA" + std::string(30, 'C') + "TGGACCATGAGTCATTGACCAGTCGATC";
    std::string plain = "GATTACAGCTTGCATGGACCATGAGTCAGTTGACCAGTCGATC";
    Costs two_piece = Costs::two_piece(4, 6, 2, 24, 1);
    dinara::Alignment two = dinara::align(gapped, plain, two_piece);
    CHECK(two.cost == 62 && two.cigar == "14=30D14=1I14=");
    CHECK(dinara::distance(gapped, plain, affine) == 74);
    CHECK(dinara::distance(gapped, plain, two_piece) == 62);
    CHECK(dinara::distance_within(gapped, plain, 62, two_piece) == 62);
    CHECK(!dinara::distance_within(gapped, plain, 61, two_piece));
    CHECK(!dinara::align_within(gapped, plain, 61, two_piece));
    CHECK(dinara::align("TTTTACGTACGTTTTT", "ACGTACGT", two_piece, Mode::infix()).cost == 0);
    // Deletions dearer than insertions: a letter missing from the query costs a deletion's 20 + 5, one
    // extra in it an insertion's 2 + 1.
    Costs dear_deletions = Costs::affine(4, 2, 1).with_deletions(20, 5);
    CHECK(dinara::align("ACGTACGT", "ACGACGT", dear_deletions).cost == 25);
    CHECK(dinara::align("ACGACGT", "ACGTACGT", dear_deletions).cost == 3);
    CHECK(dinara::align("ACGACGT", "ACGTACGT", dear_deletions).cigar == "3=1I4=");
    // A band of the one diagonal: substitutions only, where a gap either way would be cheaper.
    CHECK(dinara::align("ACGTACGT", "ACGAACGT", affine, Mode::global(), dinara::Band{0, 0}).cigar == "3=1X4=");
    CHECK(dinara::distance("AAAACCCC", "CCCCAAAA", affine, Mode::global(), dinara::Band::around(0)) == 32);
    CHECK(dinara::distance("AAAACCCC", "CCCCAAAA", affine) < 32);
    CHECK(dinara::distance("AAAACCCC", "CCCCAAAA", Costs::edit(), Mode::global(), dinara::Band::around(0)) == 8);
    CHECK(!dinara::distance_within("ACGT", "AC", 100, affine, Mode::global(), dinara::Band::around(1)));
    bool outside = false;
    try {
        dinara::distance("ACGT", "AC", affine, Mode::global(), dinara::Band::around(1));
    } catch (const dinara::OutsideBand &) {
        outside = true;
    }
    CHECK(outside);
    // A gap in a run of repeats: at its left end by default, as minimap2 places it, at its right end
    // under WFA2-lib's rule.
    CHECK(dinara::align("ACGTTTTACG", "ACGTTTACG", affine).cigar == "3=1D6=");
    CHECK(dinara::align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", Costs::edit(), Mode::global(), {}, dinara::Ties::right)
              .cigar == "4=1D5=1I3=");
    CHECK(dinara::align("ACGTTTTACG", "ACGTTTACG", affine, Mode::global(), {}, dinara::Ties::right).cigar ==
          "6=1D3=");
    // An extension stops where the read stops matching, from either end.
    dinara::Alignment right = dinara::align("ACGTTGCAAGGCTTTTTTTTTT", "ACGTTGCAAGGCGAGAGAGAGA", affine, Mode::extension(1));
    CHECK(right.score == 12 && right.cost == 0 && right.cigar == "12=" && right.reference_end == 12 &&
          right.query_end == 12);
    dinara::Alignment left = dinara::align("TTTTTTTTTTACGTTGCAAGGC", "GAGAGAGAGAACGTTGCAAGGC", affine,
                                           Mode::extension(1, dinara::Anchor::end));
    CHECK(left.score == 12 && left.cigar == "12=" && left.reference_start == 10 && left.query_start == 10);
    CHECK(dinara::align("ACGTTGCAAGGCTTTT", "ACGTTGCAAGGCGAGA", two_piece, Mode::extension(1)).score == 12);
    // A Z-drop gives up at the noise: the matching stretch past it is never reached.
    std::string seed_core(40, 'A'), noise_a = "CGTCGTCGTGCTGCTGACGT", noise_b = "TGCATGCATTGACGTACGTG",
                                    rest = std::string(30, 'C') + std::string(30, 'G') + std::string(30, 'T');
    CHECK(dinara::align(seed_core + noise_a + rest, seed_core + noise_b + rest, affine, Mode::extension(1, dinara::Anchor::start, 10))
              .reference_end == 40);
    CHECK(dinara::align(seed_core + noise_a + rest, seed_core + noise_b + rest, affine, Mode::extension(1)).reference_end >
          40);
    // A local alignment: the shared core, whichever ends surround it.
    dinara::Alignment local = dinara::align("GGGGACGTACGTGGGG", "CCCCACGTACGTCCCC", affine, Mode::local(2));
    CHECK(local.score == 16 && local.cigar == "8=" && local.reference_start == 4 && local.query_start == 4);
    CHECK(dinara::align("ACGT", "ACGT", Costs::edit(), Mode::local(1)).score == 4);
    // A score alone, with no alignment, and SSW's second best: the later copy of the core is the best end.
    CHECK(dinara::score("GGGGACGTACGTGGGG", "CCCCACGTACGTCCCC", affine, Mode::local(2)) == 16);
    CHECK(dinara::score("ACGTACGTTTGCA", "ACGTCGTTTTGCA", affine) == -12);
    dinara_local_scores twice = dinara::local_scores("ACGTACGTAC" + std::string(20, 'T') + "ACGTACGTAC", "ACGTACGTAC",
                                                     affine, Mode::local(2), 5);
    CHECK(twice.score == 20 && twice.reference_end == 40 && twice.second_score == 20 && twice.second_reference_end == 10);
    // An overlap of two reads, the first's suffix on the second's prefix.
    dinara::Alignment over = dinara::align("TTTTTACGTACGT", "ACGTACGTGGGGG", affine, Mode::overlap(2));
    CHECK(over.score == 16 && over.cigar == "8=" && over.reference_start == 5 && over.query_end == 8);
    // A read placed in a window as a mapper scores it, a match earning 2: 7 matches, one gap of 8.
    dinara::Alignment placed_scored =
        dinara::align("TTTTACGTACGTTTTT", "ACGTCGT", affine, Mode::infix().with_match_score(2));
    CHECK(placed_scored.score == 6 && placed_scored.cigar == "4=1D3=" && placed_scored.reference_start == 4);
    // The whole reference inside the query: one gap where it lies at 4, or two mismatches at 8, both 8.
    dinara::Alignment inside = dinara::align("ACGTCGT", "TTTTACGTACGTTTTT", affine, Mode::reference_in_query());
    CHECK(inside.cost == 8 && inside.reference_start == 0 && inside.reference_end == 7);
    bool no_cost = false;
    try {
        dinara::distance("ACGT", "ACGT", affine, Mode::extension(1));
    } catch (const std::invalid_argument &) {
        no_cost = true;
    }
    CHECK(no_cost);

    // The C structs themselves, null for every default.
    dinara_alignment raw{};
    CHECK(dinara_align("ACGTACGTTTGCA", 13, "ACGTCGTTTTGCA", 13, nullptr, nullptr, nullptr, &raw) == 0);
    CHECK(raw.cost == 2 && std::string(raw.cigar) == "4=1D2=1I6=" && raw.cigar_length == 10);
    dinara_free(raw.cigar);
    dinara_costs c_affine{4, 6, 2, -1, 0, 0, 0, -1, 0};
    dinara_mode c_infix{DINARA_ENDS_FREE, DINARA_ALL, DINARA_ALL, 0, 0, 0, 0, 0};
    CHECK(dinara_distance("TTTTACGTACGTTTTT", 16, "ACGTACGT", 8, &c_affine, &c_infix, nullptr) == 0);
    dinara_options capped{INT64_MIN, INT64_MAX, 11, 1, 0, 0};
    CHECK(dinara_distance("ACGTACGTTTGCA", 13, "ACGTCGTTTTGCA", 13, &c_affine, nullptr, &capped) == DINARA_ABOVE_MAX);
    CHECK(dinara_align("ACGTACGTTTGCA", 13, "ACGTCGTTTTGCA", 13, &c_affine, nullptr, &capped, &raw) == DINARA_ABOVE_MAX);
    dinara_costs c_free{0, 6, 2, -1, 0, 0, 0, -1, 0};
    CHECK(dinara_distance("A", 1, "C", 1, &c_free, nullptr, nullptr) == DINARA_INVALID_COSTS);

    std::mt19937 random(7);
    std::uniform_int_distribution<int> base(0, 3);
    std::vector<std::pair<std::string, std::string>> pairs;
    for (double rate : {0.01, 0.05, 0.15, 0.3}) {
        std::string first;
        for (int index = 0; index < 20000; ++index) first += "ACGT"[base(random)];
        pairs.emplace_back(first, mutated(first, rate, random));
    }
    std::vector<int64_t> distances;
    for (auto &pair : pairs) {
        dinara::Alignment long_aligned = dinara::align(pair.first, pair.second);
        CHECK(long_aligned.cost == dinara::distance(pair.first, pair.second));
        CHECK(edits(long_aligned.cigar, pair.first.size(), pair.second.size()) == long_aligned.cost);
        distances.push_back(long_aligned.cost);
        // Unit costs under a cap or a band take the wavefront, and the same CIGAR.
        CHECK(dinara::align_within(pair.first, pair.second, long_aligned.cost, Costs::edit())->cigar ==
              long_aligned.cigar);
        CHECK(dinara::distance(pair.first, pair.second, Costs::linear(2, 2)) == 2 * long_aligned.cost);
        // A memory budget of a few kilobytes splits the pair again and again, at the same cost.
        dinara::Alignment small = dinara::align(pair.first, pair.second, Costs::affine(4, 6, 2), Mode::global(), {},
                                                dinara::Ties::left, true, 4096);
        CHECK(small.cost == dinara::distance(pair.first, pair.second, Costs::affine(4, 6, 2)));
        CHECK(!dinara::distance_within(pair.first, pair.second, long_aligned.cost - 1));
    }

    // No state between calls: four threads aligning the same pairs agree with the single thread.
    std::vector<std::vector<int64_t>> seen(4);
    std::vector<std::thread> threads;
    for (int worker = 0; worker < 4; ++worker)
        threads.emplace_back([&, worker] {
            for (auto &pair : pairs) seen[worker].push_back(dinara::align(pair.first, pair.second).cost);
        });
    for (auto &thread : threads) thread.join();
    for (auto &worker : seen) CHECK(worker == distances);

    // A batch over every thread, and over two, agrees with the pairs one at a time.
    std::vector<std::string_view> firsts, seconds;
    for (auto &pair : pairs) firsts.push_back(pair.first), seconds.push_back(pair.second);
    CHECK(dinara::distances(firsts, seconds) == distances);
    CHECK(dinara::distances(firsts, seconds, Costs::edit(), Mode::global(), {}, 2) == distances);
    std::vector<dinara::Alignment> batch = dinara::alignments(firsts, seconds);
    for (size_t index = 0; index < batch.size(); ++index) {
        CHECK(batch[index].cost == distances[index]);
        CHECK(batch[index].cigar == dinara::align(pairs[index].first, pairs[index].second).cigar);
    }
    // Under a cap, the pairs past it come back empty, the rest as they were.
    int64_t cap = distances[1];
    std::vector<std::optional<int64_t>> within = dinara::distances_within(firsts, seconds, cap);
    std::vector<std::optional<dinara::Alignment>> capped_aligned = dinara::alignments_within(firsts, seconds, cap);
    for (size_t index = 0; index < pairs.size(); ++index) {
        CHECK(within[index].has_value() == (distances[index] <= cap));
        CHECK(capped_aligned[index].has_value() == (distances[index] <= cap));
        if (within[index]) CHECK(*within[index] == distances[index]);
        if (capped_aligned[index]) CHECK(capped_aligned[index]->cigar == batch[index].cigar);
    }
    // A pair's own failure stays its own: a sentinel byte fails that pair alone.
    std::vector<std::string_view> odd_firsts{"ACGT", "\xfe", "ACGA"}, odd_seconds{"ACGT", "A", "ACGT"};
    std::vector<int64_t> codes(3);
    dinara::detail::Batch odd(odd_firsts, odd_seconds);
    CHECK(dinara_distances(3, odd.references.data(), odd.reference_lengths.data(), odd.queries.data(),
                           odd.query_lengths.data(), nullptr, nullptr, nullptr, 0, codes.data()) == 0);
    CHECK(codes[0] == 0 && codes[1] == DINARA_UNSUPPORTED_SYMBOLS && codes[2] == 1);
    bool batch_refused = false;
    try {
        dinara::distances(odd_firsts, odd_seconds);
    } catch (const dinara::UnsupportedSymbols &) {
        batch_refused = true;
    }
    CHECK(batch_refused);
    CHECK(dinara_distances(3, odd.references.data(), odd.reference_lengths.data(), odd.queries.data(),
                           odd.query_lengths.data(), &c_free, nullptr, nullptr, 0, codes.data()) ==
          DINARA_INVALID_COSTS);

    if (failures) {
        std::fprintf(stderr, "%d checks failed\n", failures);
        return 1;
    }
    std::printf("C API: every check passed\n");
    return 0;
}
