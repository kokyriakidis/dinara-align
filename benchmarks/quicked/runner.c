// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
// MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// QuickEd (Doblas et al., Bioinformatics 2025) as the benchmarks' runners run every aligner:
//
//     runner seq quicked <budget seconds> <file.seq>...
//
// Each pair of a `.seq` file, a `>` line and a `<` line, is aligned once with its traceback at edit
// distance, by QuickEd's default bound-and-align, until the budget is spent; one row a pair, flushed as
// it finishes: the tool, the file, the seconds, the cost and the pair's growth of peak resident memory.
// As QuickEd's own `align_benchmark` does, each pair gets an aligner of its own and only the alignment
// call is timed.

#include <quicked.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>
#include <time.h>

static double now(void) {
    struct timespec moment;
    clock_gettime(CLOCK_MONOTONIC, &moment);
    return (double)moment.tv_sec + (double)moment.tv_nsec * 1e-9;
}

// The process's peak resident memory so far, in bytes: `getrusage` counts kilobytes on Linux and bytes
// on macOS.
static long peak_resident(void) {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
#ifdef __APPLE__
    return (long)usage.ru_maxrss;
#else
    return (long)usage.ru_maxrss * 1024;
#endif
}

// The next line of `file` without its newline, in a buffer `read_line` grows; NULL at the end.
static char *read_line(FILE *file, char **buffer, size_t *capacity) {
    ssize_t length = getline(buffer, capacity, file);
    if (length < 0) return NULL;
    while (length > 0 && ((*buffer)[length - 1] == '\n' || (*buffer)[length - 1] == '\r')) (*buffer)[--length] = 0;
    return *buffer;
}

int main(int count, char **arguments) {
    if (count < 5 || strcmp(arguments[1], "seq") != 0 || strcmp(arguments[2], "quicked") != 0) {
        fprintf(stderr, "usage: runner seq quicked <budget seconds> <file.seq>...\n");
        return 2;
    }
    double budget = atof(arguments[3]);
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
            // `read_line` reuses its buffer, so the pattern moves out before the text is read.
            char *first = strdup(pattern + 1);
            if (!read_line(file, &text, &text_capacity)) {
                free(first);
                break;
            }
            quicked_params_t params = quicked_default_params();
            quicked_aligner_t aligner;
            quicked_status_t status = quicked_new(&aligner, &params);
            if (quicked_check_error(status)) {
                fprintf(stderr, "%s\n", quicked_status_msg(status));
                return 1;
            }
            long before = peak_resident();
            double started = now();
            status = quicked_align(&aligner, first, (int)strlen(first), text + 1, (int)strlen(text + 1));
            double seconds = now() - started;
            long growth = peak_resident() - before;
            if (quicked_check_error(status)) {
                fprintf(stderr, "%s\n", quicked_status_msg(status));
                return 1;
            }
            spent += seconds;
            printf("quicked\t%s\t%.9f\t%d\t%ld\n", arguments[index], seconds, aligner.score, growth);
            fflush(stdout);
            quicked_free(&aligner);
            free(first);
        }
        fclose(file);
    }
    free(pattern);
    free(text);
    return 0;
}
