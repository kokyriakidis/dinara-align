// The C API through its C++ wrapper (see c/dinara.h): known distances and CIGARs, empty sequences, symbols
// past ACGT, too many of them, long pairs whose CIGAR spells the distance, the affine cost, and four
// threads at once.
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
    CHECK(dinara::edit_distance("ACGTACGTTTGCA", "ACGTCGTTTTGCA") == 2);
    dinara::Alignment aligned = dinara::edit_cigar("ACGTACGTTTGCA", "ACGTCGTTTTGCA");
    CHECK(aligned.distance == 2 && aligned.cigar == "4=1D2=1I6=");
    CHECK(dinara::edit_cigar("ACGTACGTTTGCA", "ACGTCGTTTTGCA", false).cigar == "4M1D2M1I6M");

    CHECK(dinara::edit_cigar("", "ACG").cigar == "3I");
    CHECK(dinara::edit_cigar("ACG", "").cigar == "3D");
    dinara::Alignment empty = dinara::edit_cigar("", "");
    CHECK(empty.distance == 0 && empty.cigar.empty());

    CHECK(dinara::edit_cigar("ACGTNNAC", "ACGTNNAC").cigar == "8=");
    CHECK(dinara::edit_distance("ACGTN", "ACGTA") == 1);
    bool refused = false;
    try {
        dinara::edit_distance("NRYKM", "A");
    } catch (const dinara::UnsupportedSymbols &) {
        refused = true;
    }
    CHECK(refused);

    dinara::AffineAlignment affine = dinara::affine_cigar("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 4, 6, 2);
    CHECK(affine.cost == 12 && affine.cigar == "4=3X6=");
    CHECK(dinara::affine_cigar("ACGT", "", 4, 6, 2).cigar == "4D");
    CHECK(dinara::affine_cigar("ACGT", "", 4, 6, 2).cost == 14);
    bool invalid = false;
    try {
        dinara::affine_cigar("A", "C", 0, 6, 2);
    } catch (const std::invalid_argument &) {
        invalid = true;
    }
    CHECK(invalid);
    CHECK(dinara::affine_distance("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 4, 6, 2) == 12);
    CHECK(dinara::affine_distance_within("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 4, 6, 2, 12) == 12);
    CHECK(!dinara::affine_distance_within("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 4, 6, 2, 11));
    CHECK(!dinara::affine_cigar_within("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 4, 6, 2, 11));
    CHECK(dinara::affine_cigar_within("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 4, 6, 2, 12)->cigar == "4=3X6=");
    // A read placed inside a reference, the reference's ends free.
    dinara::EndsFree inside{16, 16, 0, 0};
    dinara::AffineAlignment placed = dinara::affine_cigar("TTTTACGTACGTTTTT", "ACGTACGT", 4, 6, 2, true, inside);
    CHECK(placed.cost == 0 && placed.cigar == "4D8=4D");
    CHECK(dinara::affine_distance("TTTTACGTACGTTTTT", "ACGTACGT", 4, 6, 2) > 0);
    CHECK(dinara::affine_distance("TTTTACGTACGTTTTT", "ACGTACGT", 4, 6, 2, inside) == 0);
    // Two-piece gap costs: a long gap at the second piece, 24 + 30, and a short one at the first, 6 + 2.
    std::string gapped = "GATTACAGCTTGCA" + std::string(30, 'C') + "TGGACCATGAGTCATTGACCAGTCGATC";
    std::string plain = "GATTACAGCTTGCATGGACCATGAGTCAGTTGACCAGTCGATC";
    dinara::SecondPiece cheap_long{24, 1};
    dinara::AffineAlignment two = dinara::affine2p_cigar(gapped, plain, 4, 6, 2, cheap_long);
    CHECK(two.cost == 62 && two.cigar == "14=30D14=1I14=");
    CHECK(dinara::affine_distance(gapped, plain, 4, 6, 2) == 74);
    CHECK(dinara::affine2p_distance(gapped, plain, 4, 6, 2, cheap_long) == 62);
    CHECK(dinara::affine2p_distance_within(gapped, plain, 4, 6, 2, cheap_long, 62) == 62);
    CHECK(!dinara::affine2p_distance_within(gapped, plain, 4, 6, 2, cheap_long, 61));
    CHECK(!dinara::affine2p_cigar_within(gapped, plain, 4, 6, 2, cheap_long, 61));
    CHECK(dinara::affine2p_cigar("TTTTACGTACGTTTTT", "ACGTACGT", 4, 6, 2, cheap_long, true, inside).cost == 0);
    // A band of the one diagonal: substitutions only, where a gap either way would be cheaper.
    CHECK(dinara::affine_cigar("ACGTACGT", "ACGAACGT", 4, 6, 2, true, {}, dinara::Band{0, 0}).cigar == "3=1X4=");
    CHECK(dinara::affine_distance("AAAACCCC", "CCCCAAAA", 4, 6, 2, {}, dinara::Band::around(0)) == 32);
    CHECK(dinara::affine_distance("AAAACCCC", "CCCCAAAA", 4, 6, 2) < 32);
    CHECK(!dinara::affine_distance_within("ACGT", "AC", 4, 6, 2, 100, {}, dinara::Band::around(1)));
    bool outside = false;
    try {
        dinara::affine_distance("ACGT", "AC", 4, 6, 2, {}, dinara::Band::around(1));
    } catch (const dinara::OutsideBand &) {
        outside = true;
    }
    CHECK(outside);
    // A gap in a run of repeats: at its left end by default, as minimap2 places it, at its right end
    // under WFA2-lib's rule.
    CHECK(dinara::affine_cigar("ACGTTTTACG", "ACGTTTACG", 4, 6, 2).cigar == "3=1D6=");
    CHECK(dinara::edit_cigar("ACGTACGTTTGCA", "ACGTCGTTTTGCA", true, dinara::Ties::right).cigar == "4=1D5=1I3=");
    CHECK(dinara::affine_cigar("ACGTTTTACG", "ACGTTTACG", 4, 6, 2, true, {}, {}, dinara::Ties::right).cigar ==
          "6=1D3=");
    // An extension stops where the read stops matching, from either end.
    dinara::Extension right = dinara::affine_extension("ACGTTGCAAGGCTTTTTTTTTT", "ACGTTGCAAGGCGAGAGAGAGA", 1, 4, 6, 2);
    CHECK(right.score == 12 && right.cigar == "12=" && right.first_length == 12 && right.second_length == 12);
    dinara::Extension left =
        dinara::affine_extension("TTTTTTTTTTACGTTGCAAGGC", "GAGAGAGAGAACGTTGCAAGGC", 1, 4, 6, 2, dinara::Anchor::end);
    CHECK(left.score == 12 && left.cigar == "12=");
    CHECK(dinara::affine2p_extension("ACGTTGCAAGGCTTTT", "ACGTTGCAAGGCGAGA", 1, 4, 6, 2, cheap_long).score == 12);

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
        dinara::Alignment long_aligned = dinara::edit_cigar(pair.first, pair.second);
        CHECK(long_aligned.distance == dinara::edit_distance(pair.first, pair.second));
        CHECK(edits(long_aligned.cigar, pair.first.size(), pair.second.size()) == long_aligned.distance);
        distances.push_back(long_aligned.distance);
        // At unit costs the affine cost is the edit distance plus an opening of zero.
        CHECK(dinara::affine_cigar(pair.first, pair.second, 1, 0, 1).cost == long_aligned.distance);
        CHECK(dinara::affine_distance(pair.first, pair.second, 1, 0, 1) == long_aligned.distance);
        CHECK(!dinara::affine_distance_within(pair.first, pair.second, 1, 0, 1, long_aligned.distance - 1));
    }

    // No state between calls: four threads aligning the same pairs agree with the single thread.
    std::vector<std::vector<int64_t>> seen(4);
    std::vector<std::thread> threads;
    for (int worker = 0; worker < 4; ++worker)
        threads.emplace_back([&, worker] {
            for (auto &pair : pairs) seen[worker].push_back(dinara::edit_cigar(pair.first, pair.second).distance);
        });
    for (auto &thread : threads) thread.join();
    for (auto &worker : seen) CHECK(worker == distances);

    if (failures) {
        std::fprintf(stderr, "%d checks failed\n", failures);
        return 1;
    }
    std::printf("C API: every check passed\n");
    return 0;
}
