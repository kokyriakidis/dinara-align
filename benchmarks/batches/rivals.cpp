/* This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
 * MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/. */

// The CPU rivals of `batch_bench.py`, each scoring a whole batch of global alignments on every thread through
// OpenMP, one aligner a thread, scores alone, as Accelign's short-read case study runs them:
//
//     batch-rivals <pairs file> <workload>
//
// `illumina-affine`: WFA2-lib exact, KSW2's `extz2` with no band or Z-drop, and parasail's striped `nw`, at a
// match 0, a mismatch 1 and a gap of `k` letters `2 + k`. `illumina-edit`: WFA2-lib exact and Edlib's NW,
// the edit distance alone. Each tool's time is the faster of two passes, and its answer the sum and
// position-weighted sum of the costs.
// Prints `tool<TAB>workload<TAB>seconds<TAB>answer` rows.
#include <omp.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "edlib.h"
extern "C" {
#include "parasail.h"
#include "wavefront/wavefront_align.h"
#ifdef WITH_KSW2
#include "ksw2.h"
#endif
}

namespace {

struct Pair {
    std::string reference, query;
};

// The pairs, `name<TAB>reference<TAB>query` a line.
std::vector<Pair> load(const char *path) {
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
        pairs.push_back({first + 1, second + 1});
    }
    free(line);
    fclose(in);
    return pairs;
}

// Times one tool over the batch on every thread, the faster of two passes, and prints its row. `score`
// takes a pair and the thread's own state, which `make` builds once a thread and `drop` frees.
template <class State, class Make, class Score, class Drop>
void time_batch(const char *tool, const char *workload, const std::vector<Pair> &pairs, Make make, Score score,
                Drop drop) {
    double best = 1e300;
    long total = 0, weighted = 0;
    for (int pass = 0; pass < 2; ++pass) {
        long sum = 0, positioned = 0;
        double started = omp_get_wtime();
#pragma omp parallel reduction(+ : sum, positioned)
        {
            State state = make();
#pragma omp for schedule(dynamic, 256)
            for (long index = 0; index < static_cast<long>(pairs.size()); ++index) {
                long cost = score(pairs[index], state);
                sum += cost;
                positioned += (index + 1) * cost;
            }
            drop(state);
        }
        best = std::min(best, omp_get_wtime() - started);
        total = sum;
        weighted = positioned;
    }
    printf("%s\t%s\t%.9g\t%ld:%ld\n", tool, workload, best, total, weighted);
}

// Bases as the codes KSW2 reads.
void codes(const std::string &text, std::vector<uint8_t> &out) {
    out.resize(text.size());
    for (size_t i = 0; i < text.size(); ++i) {
        char c = text[i];
        out[i] = c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : 4;
    }
}

}  // namespace

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "usage: batch-rivals <pairs file> <workload>\n");
        return 1;
    }
    std::vector<Pair> pairs = load(argv[1]);
    std::string workload = argv[2];
    if (workload == "illumina-affine") {
        // WFA2-lib exact, its score alone: its defaults turn on the WF-adaptive heuristic, so it is turned off.
        time_batch<wavefront_aligner_t *>(
            "WFA2-lib", argv[2], pairs,
            [] {
                wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
                attributes.distance_metric = gap_affine;
                attributes.alignment_scope = compute_score;
                attributes.heuristic.strategy = wf_heuristic_none;
                attributes.affine_penalties = {0, 1, 2, 1};
                return wavefront_aligner_new(&attributes);
            },
            [](const Pair &pair, wavefront_aligner_t *aligner) {
                wavefront_align(aligner, pair.query.data(), static_cast<int>(pair.query.size()), pair.reference.data(),
                                static_cast<int>(pair.reference.size()));
                // WFA2-lib scores a cost as its negative.
                return labs(aligner->cigar->score);
            },
            [](wavefront_aligner_t *aligner) { wavefront_aligner_delete(aligner); });
#ifdef WITH_KSW2
        // KSW2's global score alone, with no band and no Z-drop: a match 0, a mismatch -1, a gap of `l`
        // letters `-(2 + l)`.
        struct Buffers {
            std::vector<uint8_t> query, target;
        };
        static int8_t table[25];
        for (int row = 0; row < 5; ++row)
            for (int column = 0; column < 5; ++column) table[row * 5 + column] = row == column && row < 4 ? 0 : -1;
        time_batch<Buffers *>(
            "KSW2", argv[2], pairs, [] { return new Buffers(); },
            [](const Pair &pair, Buffers *buffers) {
                codes(pair.query, buffers->query);
                codes(pair.reference, buffers->target);
                ksw_extz_t found;
                memset(&found, 0, sizeof(found));
                ksw_extz2_sse(nullptr, static_cast<int>(buffers->query.size()), buffers->query.data(),
                              static_cast<int>(buffers->target.size()), buffers->target.data(), 5, table, 2, 1, -1, -1,
                              0, KSW_EZ_SCORE_ONLY, &found);
                return static_cast<long>(-found.score);
            },
            [](Buffers *buffers) { delete buffers; });
#endif
        // parasail's striped NW, its score alone: an opening of 3 and an extension of 1, as it counts a gap.
        parasail_matrix_t *matrix = parasail_matrix_create("ACGT", 0, -1);
        time_batch<int>(
            "parasail", argv[2], pairs, [] { return 0; },
            [matrix](const Pair &pair, int) {
                parasail_result_t *result =
                    parasail_nw_striped_sat(pair.query.data(), static_cast<int>(pair.query.size()), pair.reference.data(),
                                            static_cast<int>(pair.reference.size()), 3, 1, matrix);
                long cost = -parasail_result_get_score(result);
                parasail_result_free(result);
                return cost;
            },
            [](int) {});
        parasail_matrix_free(matrix);
    } else if (workload == "illumina-edit") {
        // WFA2-lib exact at unit costs, its score alone, its WF-adaptive heuristic off as above.
        time_batch<wavefront_aligner_t *>(
            "WFA2-lib", argv[2], pairs,
            [] {
                wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
                attributes.distance_metric = edit;
                attributes.alignment_scope = compute_score;
                attributes.heuristic.strategy = wf_heuristic_none;
                return wavefront_aligner_new(&attributes);
            },
            [](const Pair &pair, wavefront_aligner_t *aligner) {
                wavefront_align(aligner, pair.query.data(), static_cast<int>(pair.query.size()), pair.reference.data(),
                                static_cast<int>(pair.reference.size()));
                return labs(aligner->cigar->score);
            },
            [](wavefront_aligner_t *aligner) { wavefront_aligner_delete(aligner); });
        // Edlib's global edit distance alone, with its own band, as Accelign's study ran it.
        time_batch<int>(
            "Edlib", argv[2], pairs, [] { return 0; },
            [](const Pair &pair, int) {
                EdlibAlignResult found =
                    edlibAlign(pair.query.data(), static_cast<int>(pair.query.size()), pair.reference.data(),
                               static_cast<int>(pair.reference.size()),
                               edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, nullptr, 0));
                long distance = found.editDistance;
                edlibFreeAlignResult(found);
                return distance;
            },
            [](int) {});
    }
    return 0;
}
