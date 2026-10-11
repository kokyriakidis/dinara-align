// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
// MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// BSAlign (Shao and Ruan, Bioinformatics 2024) as the benchmarks' runners run every aligner:
//
//     runner seq bsalign-edit <budget seconds> <file.seq>...
//     runner seq bsalign:<x,o,e> <budget seconds> <file.seq>...
//
// Each pair of a `.seq` file, a `>` line and a `<` line, is aligned globally once with its traceback until
// the budget is spent; one row a pair, flushed as it finishes: the tool, the file, the seconds, the cost
// and the pair's growth of peak resident memory. `bsalign-edit` is its striped bit-vector edit distance,
// `bsalign:x,o,e` its striped 8-bit difference recurrence at WFA's mismatch, opening and extension costs
// (a gap of length l costs o + e l, as BSAlign counts it too). Both run with the band off, a band as wide
// as the query, so both are exact. Only the alignment call is timed.

// `dna.h` holds the two-bit code table, `base_bit_table`, that `bsalign.h` leaves out.
#include <dna.h>
#include <bsalign.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>

int verbose = 0;

static double now(void) {
    struct timespec moment;
    clock_gettime(CLOCK_MONOTONIC, &moment);
    return (double)moment.tv_sec + (double)moment.tv_nsec * 1e-9;
}

// The process's peak resident memory so far, in bytes: `getrusage` counts kilobytes on Linux.
static long peak_resident(void) {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    return (long)usage.ru_maxrss * 1024;
}

// The next line of `file` without its newline, in a buffer `read_line` grows; NULL at the end.
static char *read_line(FILE *file, char **buffer, size_t *capacity) {
    ssize_t length = getline(buffer, capacity, file);
    if (length < 0) return NULL;
    while (length > 0 && ((*buffer)[length - 1] == '\n' || (*buffer)[length - 1] == '\r')) (*buffer)[--length] = 0;
    return *buffer;
}

// `letters` as BSAlign's two-bit codes, in `codes`, grown to fit.
static u4i encode(const char *letters, u1v *codes) {
    u4i length = (u4i)strlen(letters);
    clear_and_encap_u1v(codes, length);
    for (u4i index = 0; index < length; index++) codes->buffer[index] = base_bit_table[(int)letters[index]];
    codes->size = length;
    return length;
}

int main(int count, char **arguments) {
    if (count < 5 || strcmp(arguments[1], "seq") != 0) {
        fprintf(stderr, "usage: runner seq bsalign-edit|bsalign:<x,o,e> <budget seconds> <file.seq>...\n");
        return 2;
    }
    const char *tool = arguments[2];
    int edit = strcmp(tool, "bsalign-edit") == 0, mismatch = 0, open = 0, extend = 0;
    if (!edit && sscanf(tool, "bsalign:%d,%d,%d", &mismatch, &open, &extend) != 3) {
        fprintf(stderr, "unknown tool %s\n", tool);
        return 2;
    }
    double budget = atof(arguments[3]);
    b1i matrix[16];
    banded_striped_epi8_seqalign_set_score_matrix(matrix, 0, -mismatch);
    b1v *memory = adv_init_b1v(1024, 0, WORDSIZE, 0);
    u4v *cigars = init_u4v(64);
    u1v *pattern_codes = init_u1v(1024), *text_codes = init_u1v(1024);
    // A short spin first, so the scheduler has moved this process onto a fast core.
    for (double started = now(); now() - started < 0.05;) {
    }
    double spent = 0;
    char *pattern = NULL, *text = NULL;
    size_t pattern_capacity = 0, text_capacity = 0;
    for (int index = 4; index < count && spent < budget; index++) {
        FILE *file = fopen(arguments[index], "r");
        if (!file) {
            perror(arguments[index]);
            return 1;
        }
        while (spent < budget) {
            if (!read_line(file, &pattern, &pattern_capacity)) break;
            u4i pattern_length = encode(pattern + 1, pattern_codes);
            if (!read_line(file, &text, &text_capacity)) break;
            u4i text_length = encode(text + 1, text_codes);
            long before = peak_resident();
            double started = now();
            seqalign_result_t result;
            if (edit) {
                result = striped_seqedit_pairwise(pattern_codes->buffer, pattern_length, text_codes->buffer, text_length,
                                                  SEQALIGN_MODE_GLOBAL, 0, memory, cigars, 0);
            } else {
                result = banded_striped_epi8_seqalign_pairwise(
                    pattern_codes->buffer, pattern_length, text_codes->buffer, text_length, memory, cigars,
                    SEQALIGN_MODE_GLOBAL, roundup_times(pattern_length, 16), matrix, -open, -extend, 0, 0, 0);
            }
            double seconds = now() - started;
            long growth = peak_resident() - before;
            spent += seconds;
            // The cost from the alignment's own counts, so a score's sign convention cannot mislead.
            long cost = edit ? (long)result.mis + result.ins + result.del : -(long)result.score;
            printf("%s\t%s\t%.9f\t%ld\t%ld\n", tool, arguments[index], seconds, cost, growth);
            fflush(stdout);
        }
        fclose(file);
    }
    free(pattern);
    free(text);
    return 0;
}
