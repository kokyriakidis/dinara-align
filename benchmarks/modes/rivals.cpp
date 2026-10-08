/* This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
 * MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/. */

// The rivals of `mode_bench.py`, each in the mode dinara-align runs on the same workload (see there), each
// with its CIGAR:
//
//     rivals <workload file>
//
// Edlib (HW, SHW) and WFA2-lib (ends-free) for free ends, KSW2 for the extension and two-piece gaps, WFA2-lib
// for two-piece gaps, parasail and SSW for a DNA substitution table. Each pair of the file,
// `name<TAB>reference<TAB>query`, is aligned by each tool that offers the workload's mode; a tool's time is
// the faster of two passes over the file, its mean per pair, and its answer the sum and position-weighted sum
// of each pair's cost, or score where the mode rewards matches. Prints
// `tool<TAB>workload<TAB>seconds<TAB>answer` rows.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <functional>
#include <string>
#include <vector>

#include "edlib.h"
extern "C" {
#include "parasail.h"
#include "ssw.h"
#include "wavefront/wavefront_align.h"
#ifdef WITH_KSW2
#include "ksw2.h"
#endif
}

namespace {

double now() {
    timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec / 1e9;
}

struct Pair {
    std::string reference, query;
};

// The workload's name and its pairs, `name<TAB>reference<TAB>query` a line.
std::vector<Pair> load(const char *path, std::string &name) {
    std::vector<Pair> pairs;
    FILE *in = fopen(path, "r");
    char *line = nullptr;
    size_t capacity = 0;
    ssize_t length;
    while ((length = getline(&line, &capacity, in)) > 0) {
        if (line[length - 1] == '\n') line[--length] = 0;
        char *first = strchr(line, '\t');
        if (!first) continue;
        char *second = strchr(first + 1, '\t');
        *first = *second = 0;
        name = line;
        pairs.push_back({first + 1, second + 1});
    }
    free(line);
    fclose(in);
    return pairs;
}

// Times `aligner` over the pairs, the faster of two passes, and prints its row.
void time_tool(const char *tool, const std::string &workload, const std::vector<Pair> &pairs,
               const std::function<long(const Pair &)> &aligner) {
    double best = 1e300;
    long total = 0, weighted = 0;
    for (int pass = 0; pass < 2; ++pass) {
        total = weighted = 0;
        double started = now();
        for (size_t index = 0; index < pairs.size(); ++index) {
            long found = aligner(pairs[index]);
            total += found;
            weighted += static_cast<long>(index + 1) * found;
        }
        best = std::min(best, now() - started);
    }
    printf("%s\t%s\t%.9g\t%ld:%ld\n", tool, workload.c_str(), best / pairs.size(), total, weighted);
}

// Edlib's edit distance of the query against the reference, in `mode`, its path traced.
long edlib(const Pair &pair, EdlibAlignMode mode) {
    EdlibAlignResult found = edlibAlign(pair.query.data(), pair.query.size(), pair.reference.data(),
                                        pair.reference.size(), edlibNewAlignConfig(-1, mode, EDLIB_TASK_PATH, nullptr, 0));
    long distance = found.editDistance;
    edlibFreeAlignResult(found);
    return distance;
}

// WFA2-lib's cost of the query, the pattern, against the reference, the text, its CIGAR computed: the text's
// ends free when `infix` or its end alone when `prefix`, else end to end.
long wfa(wavefront_aligner_t *aligner, const Pair &pair, bool infix, bool prefix) {
    int text = static_cast<int>(pair.reference.size());
    if (infix || prefix) {
        wavefront_aligner_set_alignment_free_ends(aligner, 0, 0, infix ? text : 0, text);
    } else {
        wavefront_aligner_set_alignment_end_to_end(aligner);
    }
    wavefront_align(aligner, pair.query.data(), static_cast<int>(pair.query.size()), pair.reference.data(), text);
    // WFA2-lib scores a cost as its negative.
    return labs(aligner->cigar->score);
}

// An exact WFA2-lib aligner keeping every front, for its CIGAR, under `metric` at WFA's (4, 6, 2), or two-piece
// (4, 6, 2, 24, 1).
wavefront_aligner_t *wfa_aligner(distance_metric_t metric) {
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = metric;
    attributes.alignment_scope = compute_alignment;
    // Exact: WFA2-lib's defaults turn on its WF-adaptive heuristic, which drops lagging diagonals.
    attributes.heuristic.strategy = wf_heuristic_none;
    attributes.affine_penalties = {0, 4, 6, 2};
    attributes.affine2p_penalties = {0, 4, 6, 2, 24, 1};
    return wavefront_aligner_new(&attributes);
}

// Bases as the codes SSW and KSW2 read: A, C, G, T, then anything else.
std::vector<int8_t> codes(const std::string &text) {
    std::vector<int8_t> out(text.size());
    for (size_t i = 0; i < text.size(); ++i) {
        char c = text[i];
        out[i] = c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : 4;
    }
    return out;
}

#ifdef WITH_KSW2
// A five-letter KSW2 table: `match` on the diagonal, `mismatch` elsewhere, zero against the fifth letter.
std::vector<int8_t> ksw2_table(int8_t match, int8_t mismatch) {
    std::vector<int8_t> table(25, 0);
    for (int row = 0; row < 4; ++row)
        for (int column = 0; column < 4; ++column) table[row * 5 + column] = row == column ? match : mismatch;
    return table;
}

// KSW2's best extension score from both sequences' starts, a match 2, a mismatch -4, a gap of `k` letters
// `6 + 2k`, no band and no Z-drop: the exact best stop, as dinara-align's `Mode.extension(2)` finds it.
long ksw2_extension(const Pair &pair, const std::vector<int8_t> &table) {
    std::vector<int8_t> query = codes(pair.query), target = codes(pair.reference);
    ksw_extz_t found;
    memset(&found, 0, sizeof(found));
    ksw_extz2_sse(nullptr, static_cast<int>(query.size()), reinterpret_cast<uint8_t *>(query.data()),
                  static_cast<int>(target.size()), reinterpret_cast<uint8_t *>(target.data()), 5, table.data(), 6, 2,
                  -1, -1, 0, KSW_EZ_EXTZ_ONLY, &found);
    free(found.cigar);
    return found.max;
}

// KSW2's global cost at two-piece gaps, a gap of `k` letters the less of `6 + 2k` and `24 + k`, from its
// score with no reward for a match.
long ksw2_two_piece(const Pair &pair, const std::vector<int8_t> &table) {
    std::vector<int8_t> query = codes(pair.query), target = codes(pair.reference);
    ksw_extz_t found;
    memset(&found, 0, sizeof(found));
    ksw_extd2_sse(nullptr, static_cast<int>(query.size()), reinterpret_cast<uint8_t *>(query.data()),
                  static_cast<int>(target.size()), reinterpret_cast<uint8_t *>(target.data()), 5, table.data(), 6, 2,
                  24, 1, -1, -1, 0, 0, &found);
    free(found.cigar);
    return -found.score;
}
#endif

// parasail's score of the query against the reference, globally or locally, its CIGAR read back.
long parasail(const Pair &pair, const parasail_matrix_t *matrix, int open, int extend, bool local) {
    const char *query = pair.query.data(), *reference = pair.reference.data();
    int qlen = static_cast<int>(pair.query.size()), rlen = static_cast<int>(pair.reference.size());
    parasail_result_t *result = local ? parasail_sw_trace_striped_sat(query, qlen, reference, rlen, open, extend, matrix)
                                      : parasail_nw_trace_striped_sat(query, qlen, reference, rlen, open, extend, matrix);
    parasail_cigar_t *cigar = parasail_result_get_cigar(result, query, qlen, reference, rlen, matrix);
    long score = parasail_result_get_score(result);
    parasail_cigar_free(cigar);
    parasail_result_free(result);
    return score;
}

// SSW's best local score, with its CIGAR, over `size` codes that `code_of` gives each letter.
long ssw(const Pair &pair, const std::vector<int8_t> &table, int size, const std::function<int8_t(char)> &code_of,
         int open, int extend) {
    std::vector<int8_t> query(pair.query.size()), reference(pair.reference.size());
    for (size_t i = 0; i < query.size(); ++i) query[i] = code_of(pair.query[i]);
    for (size_t i = 0; i < reference.size(); ++i) reference[i] = code_of(pair.reference[i]);
    int qlen = static_cast<int>(query.size());
    s_profile *profile = ssw_init(query.data(), qlen, table.data(), size, 2);
    s_align *found = ssw_align(profile, reference.data(), static_cast<int>(reference.size()), open, extend, 1, 0, 0,
                               qlen / 2 < 15 ? 15 : qlen / 2);
    long score = found->score1;
    align_destroy(found);
    init_destroy(profile);
    return score;
}

}  // namespace

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "usage: rivals <workload file>\n");
        return 1;
    }
    std::string workload;
    std::vector<Pair> pairs = load(argv[1], workload);
    bool infix = workload.rfind("infix", 0) == 0, prefix = workload == "prefix-edit";
    if (workload.rfind("infix-edit", 0) == 0 || prefix) {
        time_tool("Edlib", workload, pairs, [&](const Pair &pair) { return edlib(pair, prefix ? EDLIB_MODE_SHW : EDLIB_MODE_HW); });
        wavefront_aligner_t *aligner = wfa_aligner(edit);
        time_tool("WFA2-lib", workload, pairs, [&](const Pair &pair) { return wfa(aligner, pair, infix, prefix); });
        wavefront_aligner_delete(aligner);
    } else if (workload == "infix-affine" || workload == "two-piece") {
        wavefront_aligner_t *aligner = wfa_aligner(workload == "two-piece" ? gap_affine_2p : gap_affine);
        time_tool("WFA2-lib", workload, pairs, [&](const Pair &pair) { return wfa(aligner, pair, infix, false); });
        wavefront_aligner_delete(aligner);
#ifdef WITH_KSW2
        if (workload == "two-piece") {
            std::vector<int8_t> table = ksw2_table(0, -4);
            time_tool("KSW2", workload, pairs, [&](const Pair &pair) { return ksw2_two_piece(pair, table); });
        }
#endif
    } else if (workload == "extension") {
#ifdef WITH_KSW2
        std::vector<int8_t> table = ksw2_table(2, -4);
        time_tool("KSW2", workload, pairs, [&](const Pair &pair) { return ksw2_extension(pair, table); });
#endif
    } else if (workload.rfind("table", 0) == 0) {
        // A match 2, a transition -2, a transversion -4, a gap of `k` letters `4 + 2k`: an opening of 6 and an
        // extension of 2 as parasail and SSW count them.
        parasail_matrix_t *matrix = parasail_matrix_copy(parasail_matrix_create("ACGT", 2, -4));
        std::vector<int8_t> table(25, 0);
        for (int row = 0; row < 4; ++row) {
            for (int column = 0; column < 4; ++column) {
                int value = row == column ? 2 : ((row + column) % 2 == 0 ? -2 : -4);
                parasail_matrix_set_value(matrix, row, column, value);
                table[row * 5 + column] = static_cast<int8_t>(value);
            }
        }
        bool local = workload == "table-local";
        time_tool("parasail", workload, pairs, [&](const Pair &pair) { return parasail(pair, matrix, 6, 2, local); });
        if (local) {
            auto code_of = [](char c) -> int8_t { return c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : 4; };
            time_tool("SSW", workload, pairs, [&](const Pair &pair) { return ssw(pair, table, 5, code_of, 6, 2); });
        }
        parasail_matrix_free(matrix);
    }
    return 0;
}
