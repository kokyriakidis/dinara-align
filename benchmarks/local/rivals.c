// The C rivals of `local_bench.py`: SSW, parasail and abPOA, local alignment and, where they have them,
// overlap alignment and a query placed whole in a reference scored with a reward, each with its CIGAR, at dinara-align's Mode.local(2) and Mode.overlap(2) under
// Costs.affine(4, 6, 2): a match 2, a mismatch -4, a gap of k letters 6 + 2k, which SSW and parasail
// take as an opening of 8 and an extension of 2, and abPOA as 6 and 2.
//
//     rivals <workload file>
//
// Each pair of the file, `name<TAB>reference<TAB>query`, is aligned by each tool; a tool's time is the
// faster of two passes over the file, its mean per pair, and its answer the sum and position-weighted
// sum of its scores. Prints `tool<TAB>workload<TAB>task<TAB>seconds<TAB>answer` rows.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "abpoa.h"
#include "parasail.h"
#include "ssw.h"

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + t.tv_nsec / 1e9;
}

static int8_t code(char c) { return c == 'A' ? 0 : c == 'C' ? 1 : c == 'G' ? 2 : c == 'T' ? 3 : 4; }

typedef struct {
    int count;
    char **names, **references, **queries;
} pairs_t;

static pairs_t load(const char *path) {
    FILE *in = fopen(path, "r");
    pairs_t pairs = {0, NULL, NULL, NULL};
    size_t room = 0;
    char *line = NULL;
    size_t capacity = 0;
    ssize_t length;
    while ((length = getline(&line, &capacity, in)) > 0) {
        if (line[length - 1] == '\n') line[--length] = 0;
        char *first = strchr(line, '\t');
        if (!first) continue;
        char *second = strchr(first + 1, '\t');
        *first = 0;
        *second = 0;
        if ((size_t)pairs.count == room) {
            room = room ? 2 * room : 64;
            pairs.names = realloc(pairs.names, room * sizeof(char *));
            pairs.references = realloc(pairs.references, room * sizeof(char *));
            pairs.queries = realloc(pairs.queries, room * sizeof(char *));
        }
        pairs.names[pairs.count] = strdup(line);
        pairs.references[pairs.count] = strdup(first + 1);
        pairs.queries[pairs.count] = strdup(second + 1);
        ++pairs.count;
    }
    free(line);
    fclose(in);
    return pairs;
}

typedef long (*aligner_t)(const char *reference, const char *query);

static int8_t ssw_matrix[25];
static parasail_matrix_t *parasail_dna;
static abpoa_t *ab;
static abpoa_para_t *abpt;

static long ssw_local(const char *reference, const char *query) {
    int rlen = strlen(reference), qlen = strlen(query);
    int8_t *r = malloc(rlen), *q = malloc(qlen);
    for (int i = 0; i < rlen; ++i) r[i] = code(reference[i]);
    for (int i = 0; i < qlen; ++i) q[i] = code(query[i]);
    s_profile *profile = ssw_init(q, qlen, ssw_matrix, 5, 2);
    s_align *found = ssw_align(profile, r, rlen, 8, 2, 1, 0, 0, qlen / 2 < 15 ? 15 : qlen / 2);
    long score = found->score1;
    align_destroy(found);
    init_destroy(profile);
    free(r);
    free(q);
    return score;
}

static long parasail_local(const char *reference, const char *query) {
    int rlen = strlen(reference), qlen = strlen(query);
    parasail_result_t *result = parasail_sw_trace_striped_sat(query, qlen, reference, rlen, 8, 2, parasail_dna);
    parasail_cigar_t *cigar = parasail_result_get_cigar(result, query, qlen, reference, rlen, parasail_dna);
    long score = parasail_result_get_score(result);
    parasail_cigar_free(cigar);
    parasail_result_free(result);
    return score;
}

static long parasail_infix(const char *reference, const char *query) {
    int rlen = strlen(reference), qlen = strlen(query);
    parasail_result_t *result = parasail_sg_dx_trace_striped_sat(query, qlen, reference, rlen, 8, 2, parasail_dna);
    parasail_cigar_t *cigar = parasail_result_get_cigar(result, query, qlen, reference, rlen, parasail_dna);
    long score = parasail_result_get_score(result);
    parasail_cigar_free(cigar);
    parasail_result_free(result);
    return score;
}

static long parasail_overlap(const char *reference, const char *query) {
    int rlen = strlen(reference), qlen = strlen(query);
    parasail_result_t *result = parasail_sg_trace_striped_sat(query, qlen, reference, rlen, 8, 2, parasail_dna);
    parasail_cigar_t *cigar = parasail_result_get_cigar(result, query, qlen, reference, rlen, parasail_dna);
    long score = parasail_result_get_score(result);
    parasail_cigar_free(cigar);
    parasail_result_free(result);
    return score;
}

// abPOA aligns to a graph: the reference goes in as the graph's one sequence, then the query is
// aligned to it, both timed, as any pairwise use of abPOA pays for both.
static long abpoa_local(const char *reference, const char *query) {
    int rlen = strlen(reference), qlen = strlen(query);
    uint8_t *r = malloc(rlen), *q = malloc(qlen);
    for (int i = 0; i < rlen; ++i) r[i] = code(reference[i]);
    for (int i = 0; i < qlen; ++i) q[i] = code(query[i]);
    abpoa_reset(ab, abpt, rlen);
    abpoa_res_t res;
    res.graph_cigar = 0;
    res.n_cigar = 0;
    abpoa_add_graph_alignment(ab, abpt, r, NULL, rlen, NULL, res, 0, 2, 1);
    res.graph_cigar = 0;
    res.n_cigar = 0;
    abpoa_align_sequence_to_graph(ab, abpt, q, qlen, &res);
    long score = res.best_score;
    if (res.n_cigar) free(res.graph_cigar);
    free(r);
    free(q);
    return score;
}

static void run(const char *tool, const char *task, aligner_t aligner, pairs_t *pairs) {
    double best = 1e300;
    long total = 0, weighted = 0;
    for (int pass = 0; pass < 2; ++pass) {
        total = weighted = 0;
        double started = now();
        for (int index = 0; index < pairs->count; ++index) {
            long score = aligner(pairs->references[index], pairs->queries[index]);
            total += score;
            weighted += (index + 1) * score;
        }
        double spent = now() - started;
        if (spent < best) best = spent;
    }
    printf("%s\t%s\t%s\t%.9f\t%ld:%ld\n", tool, pairs->names[0], task, best / pairs->count, total, weighted);
    fflush(stdout);
}

int main(int argc, char **argv) {
    pairs_t pairs = load(argv[1]);
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j) ssw_matrix[i * 5 + j] = (i < 4 && i == j) ? 2 : -4;
    parasail_dna = parasail_matrix_create("ACGT", 2, -4);
    ab = abpoa_init();
    abpt = abpoa_init_para();
    abpt->align_mode = ABPOA_LOCAL_MODE;
    abpt->match = 2;
    abpt->mismatch = 4;
    abpt->gap_open1 = 6;
    abpt->gap_ext1 = 2;
    abpt->gap_open2 = 0;
    abpt->gap_ext2 = 0;
    abpt->ret_cigar = 1;
    abpt->disable_seeding = 1;
    abpt->progressive_poa = 0;
    abpt->out_msa = 0;
    abpt->out_cons = 0;
    abpoa_post_set_para(abpt);
    if (strstr(pairs.names[0], "overlap")) {
        run("parasail", "overlap", parasail_overlap, &pairs);
    } else if (strstr(pairs.names[0], "infix")) {
        run("parasail", "infix", parasail_infix, &pairs);
    } else {
        run("SSW", "local", ssw_local, &pairs);
        run("parasail", "local", parasail_local, &pairs);
        run("abPOA", "local", abpoa_local, &pairs);
    }
    parasail_matrix_free(parasail_dna);
    abpoa_free(ab);
    abpoa_free_para(abpt);
    return 0;
}
