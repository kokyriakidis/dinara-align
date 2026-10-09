/* This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
 * MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/. */

/*
 * dinara-align from C and C++: the least cost of aligning a query to a reference, and an optimal
 * alignment as a CIGAR, under unit costs (the edit distance), gap-affine or two-piece gap-affine costs,
 * globally, with free ends (a read inside a reference, at its start or end, an overlap), as an
 * extension from one end, locally, within a band or under a cost cap. Every result is exact. Build the library
 * with `pixi run build-c [target-cpu]`, which leaves it in build/c beside the Mojo runtime libraries it
 * loads and this header; link with `-Lbuild/c -ldinara` and put build/c on the program's library search
 * path (an rpath, say). build/c also holds the libstdc++ the runtime loads, which a C++ program then
 * shares: it is GCC 15's, as new as any compiler the program is likely built with. Built for an
 * explicit `target-cpu`, the library runs only on CPUs with that instruction set; the default is the
 * platform's oldest (x86-64, apple-m1).
 *
 * Sequences are bytes, each matching only itself, save 0xFE and 0xFF, which UTF-8 never holds. The
 * functions keep no state between calls, so they may run on many threads at once. Failures come back
 * as a negative result, one of the `DINARA_` codes below.
 *
 * The CIGAR reads the first sequence as the reference: `=` a match, `X` a substitution (or `M` for
 * either), `D` a reference letter alone, `I` a query letter alone. It spans the aligned parts alone,
 * `reference[reference_start:reference_end]` against `query[query_start:query_end]`; letters a mode
 * leaves unaligned for nothing lie outside them.
 */
#ifndef DINARA_H
#define DINARA_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* A 0xFE or 0xFF byte in either sequence. */
#define DINARA_UNSUPPORTED_SYMBOLS (-1)
/* The CIGAR's memory could not be allocated, or costs this dear would keep more fronts than the memory
 * allowed (`max_memory`, or 80 MB for a distance or a score). */
#define DINARA_OUT_OF_MEMORY (-2)
/* Costs no alignment can be searched by: a mismatch or extension of zero, or a negative cost. */
#define DINARA_INVALID_COSTS (-3)
/* Every alignment costs more than the `max_cost` asked for, or under a cap none fits the band. */
#define DINARA_ABOVE_MAX (-4)
/* No alignment stays inside the band of diagonals asked for. */
#define DINARA_OUTSIDE_BAND (-5)
/* A mode that cannot serve what was asked: the least cost of an extension, a local alignment or an
 * overlap, which maximize a score, a cap on any of them, or a band on the last two. */
#define DINARA_INVALID_MODE (-6)
/* A sequence's length below zero. */
#define DINARA_INVALID_LENGTH (-7)

/*
 * What each edit costs: a substitution `mismatch`, a gap of `k` letters `opening + k * extension`, or
 * with `opening2` zero or more the less of that and `opening2 + k * extension2`. With
 * `deletion_extension` above zero, a deletion, a run of reference letters alone, costs by the
 * `deletion_` fields instead, as bwa's `-O del,ins`, and the others are an insertion's; zero leaves
 * deletions costing what insertions do. A null pointer is unit costs, the edit distance,
 * {1, 0, 1, -1, 0}.
 */
typedef struct {
    int64_t mismatch;            /* A substitution's cost. */
    int64_t opening;             /* A gap's cost before its letters. */
    int64_t extension;           /* Each gapped letter's cost. */
    int64_t opening2;            /* The second piece's opening; negative for one piece. */
    int64_t extension2;          /* Each gapped letter's cost on the second piece. */
    int64_t deletion_opening;    /* A deletion's `opening`, with `deletion_extension` above zero. */
    int64_t deletion_extension;  /* A deletion's `extension`; zero for deletions costing as insertions. */
    int64_t deletion_opening2;   /* A deletion's `opening2`; negative for one piece. */
    int64_t deletion_extension2; /* A deletion's `extension2`. */
} dinara_costs;

/* A `dinara_mode`'s kind: free ends, a global alignment when none is free. */
#define DINARA_ENDS_FREE 0
/* A `dinara_mode`'s kind: an extension from both sequences' starts, or ends. */
#define DINARA_EXTENSION 1
/* A `dinara_mode`'s kind: a local alignment, Smith-Waterman. */
#define DINARA_LOCAL 2
/* As many free letters as any sequence has. */
#define DINARA_ALL INT64_MAX

/*
 * Which alignments count. DINARA_ENDS_FREE: up to so many letters at each end of each sequence left
 * unaligned for nothing; all zero is a global alignment, the reference's two at DINARA_ALL a query
 * placed anywhere in it (infix), its end alone a prefix. With `match_score` zero the alignment has the
 * least cost; above zero, the best score, a match earning it, as parasail's and hyalite's semi-global
 * modes count it. DINARA_EXTENSION: fixed at both sequences' starts, or with `anchor` nonzero their
 * ends, and free to stop anywhere, a match earning `match_score`: a seed's extension, as KSW2's, with
 * its Z-drop when `zdrop` is above zero, and its end bonus when `end_bonus` is. DINARA_LOCAL: any part of each, a match earning
 * `match_score`, Smith-Waterman, as abPOA's local mode. An overlap, parasail's `sg`, is free ends with
 * every count DINARA_ALL and a match earning. A null pointer is a global alignment.
 */
typedef struct {
    int64_t kind;            /* DINARA_ENDS_FREE, DINARA_EXTENSION or DINARA_LOCAL. */
    int64_t reference_start; /* Free ends: the reference's letters free before the alignment. */
    int64_t reference_end;   /* Free ends: the reference's letters free after it. */
    int64_t query_start;     /* Free ends: the query's letters free before it. */
    int64_t query_end;       /* Free ends: the query's letters free after it. */
    int64_t match_score;     /* What a match earns; zero for free ends' least cost. */
    int64_t anchor;          /* An extension: zero fixes it at the starts, nonzero at the ends. */
    int64_t zdrop;           /* An extension: its Z-drop when above zero, else none. */
    int64_t end_bonus;       /* An extension: what reaching the query's far end earns, KSW2's end bonus. */
} dinara_mode;

/*
 * Every move stays on the diagonals `band_low ..= band_high`, a diagonal being the reference's letters
 * aligned or skipped less the query's, from the alignment's origin: the exact optimum over the
 * alignments inside; INT64_MIN and INT64_MAX for no band, KSW2's band of width `w` `-w ..= w`. A
 * `max_cost` of zero or more caps the cost, which the search proves past after about half of it.
 * `eqx` nonzero writes `=` and `X`, else `M`; of equally good alignments the CIGAR places indels
 * left, as minimap2 does, or with `right_ties` nonzero right, WFA2-lib's CIGAR byte for byte. A
 * `max_memory` above zero bounds the bytes of fronts an alignment keeps for its traceback, about 80 MB
 * by default; past it the pair is split, its cost still the least and the tie rule followed within
 * each piece. A null pointer is no band, no cap, `=` and `X`, indels left, the default memory.
 */
typedef struct {
    int64_t band_low;   /* The lowest diagonal a move may reach; INT64_MIN for no bound. */
    int64_t band_high;  /* The highest diagonal a move may reach; INT64_MAX for no bound. */
    int64_t max_cost;   /* The cost cap when zero or more; negative for none. */
    int64_t eqx;        /* Nonzero writes `=` and `X`, zero `M`. */
    int64_t right_ties; /* Nonzero places tied indels right, zero left. */
    int64_t max_memory; /* The bytes of fronts kept for the traceback; zero or less for the default. */
} dinara_options;

/*
 * An optimal alignment: its cost, its score (an extension's, a local alignment's or an overlap's
 * matches' reward less the cost, else minus the cost), the spans it aligns, and its CIGAR,
 * NUL-terminated and `cigar_length` bytes long, which the caller frees with `dinara_free`.
 */
typedef struct {
    int64_t cost;            /* The alignment's cost. */
    int64_t score;           /* Its matches' reward less its cost, or minus its cost with no reward. */
    int64_t reference_start; /* The reference's first aligned letter. */
    int64_t reference_end;   /* One past the reference's last aligned letter. */
    int64_t query_start;     /* The query's first aligned letter. */
    int64_t query_end;       /* One past the query's last aligned letter. */
    char *cigar;             /* The CIGAR, NUL-terminated, for `dinara_free`. */
    int64_t cigar_length;    /* The CIGAR's bytes, its NUL not counted. */
} dinara_alignment;

/* The least cost of aligning the query to the reference, or a DINARA_ code. */
int64_t dinara_distance(const char *reference, int64_t reference_length, const char *query, int64_t query_length,
                        const dinara_costs *costs, const dinara_mode *mode, const dinara_options *options);

/* An optimal alignment into `*alignment`; zero, or a DINARA_ code and then no CIGAR to free. */
int64_t dinara_align(const char *reference, int64_t reference_length, const char *query, int64_t query_length,
                     const dinara_costs *costs, const dinara_mode *mode, const dinara_options *options,
                     dinara_alignment *alignment);

/*
 * One thread's aligner: `dinara_distance` and `dinara_align` through memory kept from call to call, so a
 * loop of calls on one thread takes none once it is warm. Every function here runs on the caller's own
 * thread and keeps no state of its own, so an application calls them from as many threads as it likes,
 * an aligner a thread; an aligner is never used by two threads at once.
 */
typedef struct dinara_aligner dinara_aligner;

/* A new aligner, or null when memory runs out; free it with `dinara_aligner_free`. */
dinara_aligner *dinara_aligner_new(void);

/* Frees an aligner and the memory it kept; null is nothing to free. */
void dinara_aligner_free(dinara_aligner *aligner);

/* `dinara_distance` through `aligner`. */
int64_t dinara_aligner_distance(dinara_aligner *aligner, const char *reference, int64_t reference_length,
                                const char *query, int64_t query_length, const dinara_costs *costs,
                                const dinara_mode *mode, const dinara_options *options);

/* `dinara_align` through `aligner`. */
int64_t dinara_aligner_align(dinara_aligner *aligner, const char *reference, int64_t reference_length,
                             const char *query, int64_t query_length, const dinara_costs *costs,
                             const dinara_mode *mode, const dinara_options *options, dinara_alignment *alignment);

/* The best score `dinara_align` would return, with no alignment traced, into `*score`: for a mode with
 * a match score its matches' reward less its costs, else minus the least cost. Zero, or a DINARA_ code.
 * The options' cap, `eqx`, `right_ties` and memory change nothing. */
int64_t dinara_score(const char *reference, int64_t reference_length, const char *query, int64_t query_length,
                     const dinara_costs *costs, const dinara_mode *mode, const dinara_options *options, int64_t *score);

/* A local alignment's best score and the reference's and query's letters up to its end, and the best
 * score of an alignment ending more than a window of reference letters away, as SSW's `score2` and
 * `ref_end2`; zero at zero when there is none. */
typedef struct {
    int64_t score;                /* The best local alignment's score. */
    int64_t reference_end;        /* The reference's letters up to that alignment's end. */
    int64_t query_end;            /* The query's letters up to that alignment's end. */
    int64_t second_score;         /* The best score ending more than a window away, as SSW's `score2`. */
    int64_t second_reference_end; /* The reference's letters up to its end, as SSW's `ref_end2`. */
} dinara_local_scores_result;

/* `dinara_local_scores_result` for a DINARA_LOCAL mode, the window `window` letters, or for a negative
 * one half the query and at least 15, as SSW suggests; one sweep, no alignment traced. Zero, or a code. */
int64_t dinara_local_scores(const char *reference, int64_t reference_length, const char *query,
                            int64_t query_length, const dinara_costs *costs, const dinara_mode *mode,
                            int64_t window, dinara_local_scores_result *scores);

/* The query against `count` references, a database search: the hits, best first, ties by order, their
 * places into `indices` and their scores (minus their costs with no reward) into `scores`, each with
 * room for `count`; returns how many, or a DINARA_ code. `best` above zero keeps that many, and the
 * options' `max_cost` drops what passes it; a band is refused, DINARA_INVALID_MODE. A local search
 * scores many references at once, one to a SIMD lane. */
int64_t dinara_search(int64_t count, const char *const *references, const int64_t *reference_lengths,
                      const char *query, int64_t query_length, const dinara_costs *costs, const dinara_mode *mode,
                      const dinara_options *options, int64_t best, int64_t threads, int64_t *indices,
                      int64_t *scores);

/* Frees a CIGAR that `dinara_align` or `dinara_alignments` returned; null is nothing to free. */
void dinara_free(char *cigar);

/*
 * A batch: pair `i` is `references[i]` of `reference_lengths[i]` bytes against `queries[i]` of
 * `query_lengths[i]`, all under the same costs, mode and options, spread over `threads` threads, or on
 * the caller's own thread alone for zero: the library starts threads only when asked, its caller being
 * the one that knows how to spread its work. Each thread takes the next pair as soon as it is free. Both
 * return zero, or a DINARA_ code that fails the whole batch (costs or a mode no pair can take), and
 * then write nothing.
 */

/* Every pair's `dinara_distance` into `results[i]`: its least cost, or its own DINARA_ code. */
int64_t dinara_distances(int64_t pairs, const char *const *references, const int64_t *reference_lengths,
                         const char *const *queries, const int64_t *query_lengths, const dinara_costs *costs,
                         const dinara_mode *mode, const dinara_options *options, int64_t threads, int64_t *results);

/* Every pair's `dinara_align` into `alignments[i]`, its result into `statuses[i]`: zero and a CIGAR
 * the caller frees, or the pair's DINARA_ code and no CIGAR. */
int64_t dinara_alignments(int64_t pairs, const char *const *references, const int64_t *reference_lengths,
                          const char *const *queries, const int64_t *query_lengths, const dinara_costs *costs,
                          const dinara_mode *mode, const dinara_options *options, int64_t threads,
                          dinara_alignment *alignments, int64_t *statuses);

#ifdef __cplusplus
}

#include <new>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace dinara {

/* Raised on a 0xFE or 0xFF byte. */
struct UnsupportedSymbols : std::invalid_argument {
    UnsupportedSymbols() : std::invalid_argument("dinara: a 0xFE or 0xFF byte") {}
};

/* Raised when no alignment stays inside the band asked for. */
struct OutsideBand : std::invalid_argument {
    OutsideBand() : std::invalid_argument("dinara: no alignment stays inside the band") {}
};

/* What each edit costs (see `dinara_costs`). */
struct Costs {
    /* The fields of `dinara_costs`, unit costs by default. */
    int64_t mismatch = 1, opening = 0, extension = 1, opening2 = -1, extension2 = 0;
    int64_t deletion_opening = 0, deletion_extension = 0, deletion_opening2 = -1, deletion_extension2 = 0;
    /* Unit costs: the edit distance. */
    static Costs edit() { return {}; }
    /* A substitution `mismatch` and every gapped letter `gap`. */
    static Costs linear(int64_t mismatch, int64_t gap) { return {mismatch, 0, gap, -1, 0}; }
    /* A substitution `mismatch`, a gap of `k` letters `opening + k * extension`. */
    static Costs affine(int64_t mismatch, int64_t opening, int64_t extension) {
        return {mismatch, opening, extension, -1, 0};
    }
    /* A gap of `k` letters the less of `opening + k * extension` and `opening2 + k * extension2`. */
    static Costs two_piece(int64_t mismatch, int64_t opening, int64_t extension, int64_t opening2,
                           int64_t extension2) {
        return {mismatch, opening, extension, opening2, extension2};
    }
    /* These costs with deletions of their own, the others an insertion's: bwa's `-O6,5 -E1,2` is
     * `affine(4, 5, 2).with_deletions(6, 1)`. */
    Costs with_deletions(int64_t del_opening, int64_t del_extension, int64_t del_opening2 = -1,
                         int64_t del_extension2 = 0) const {
        Costs out = *this;
        out.deletion_opening = del_opening;
        out.deletion_extension = del_extension;
        out.deletion_opening2 = del_opening2;
        out.deletion_extension2 = del_extension2;
        return out;
    }
};

/* Which end of both sequences an extension is fixed at. */
enum class Anchor { start, end };

/* Which alignments count (see `dinara_mode`). */
struct Mode {
    /* The `dinara_mode` the C functions take, a global alignment by default. */
    dinara_mode fields{DINARA_ENDS_FREE, 0, 0, 0, 0, 0, 0, 0};
    /* Both sequences end to end. */
    static Mode global() { return {}; }
    /* The whole query against wherever in the reference it fits best. */
    static Mode infix() { return ends_free(DINARA_ALL, DINARA_ALL, 0, 0); }
    /* The whole query against the reference's best prefix. */
    static Mode prefix() { return ends_free(0, DINARA_ALL, 0, 0); }
    /* The whole query against the reference's best suffix. */
    static Mode suffix() { return ends_free(DINARA_ALL, 0, 0, 0); }
    /* Up to so many letters at each end of each sequence left unaligned for nothing; the query's both ends
     * free place the whole reference inside it. */
    static Mode ends_free(int64_t reference_start, int64_t reference_end, int64_t query_start, int64_t query_end) {
        return {{DINARA_ENDS_FREE, reference_start, reference_end, query_start, query_end, 0, 0, 0}};
    }
    /* These free ends with a match earning `match_score`: the best score rather than the least cost. */
    Mode with_match_score(int64_t match_score) const {
        Mode scored = *this;
        scored.fields.match_score = match_score;
        return scored;
    }
    /* The best-scoring alignment from one end, a match earning `match_score`; with `zdrop` above zero,
     * given up once it falls that far below its best, as KSW2's Z-drop gives up; with `end_bonus` above
     * zero, aligned to the query's far end when that scores within the bonus of the best stop. */
    static Mode extension(int64_t match_score, Anchor anchor = Anchor::start, int64_t zdrop = 0,
                          int64_t end_bonus = 0) {
        return {{DINARA_EXTENSION, 0, 0, 0, 0, match_score, anchor == Anchor::end ? 1 : 0, zdrop, end_bonus}};
    }
    /* The best-scoring alignment of any part of each, a match earning `match_score`: Smith-Waterman; under a
     * table's scores, its own rewards, `local()`. */
    static Mode local(int64_t match_score = 0) { return {{DINARA_LOCAL, 0, 0, 0, 0, match_score, 0, 0}}; }
    /* The best-scoring alignment with every end gap free, a match earning `match_score`: an overlap. */
    static Mode overlap(int64_t match_score) {
        // With costs alone every end free lets the empty alignment win.
        if (match_score <= 0) throw std::invalid_argument("dinara: an overlap needs a match that earns");
        return ends_free(DINARA_ALL, DINARA_ALL, DINARA_ALL, DINARA_ALL).with_match_score(match_score);
    }
};

/* The diagonals every move stays on, `low ..= high`; the default is every diagonal. */
struct Band {
    int64_t low = INT64_MIN;  /* The lowest diagonal. */
    int64_t high = INT64_MAX; /* The highest diagonal. */
    /* KSW2's band of width `w`: at most `w` diagonals from the origin's, either way. */
    static Band around(int64_t width) { return Band{-width, width}; }
};

/* Which of several equally good alignments a CIGAR spells: indels placed left, as minimap2 places them,
 * or right, WFA2-lib's CIGAR byte for byte. */
enum class Ties { left, right };

/* An optimal alignment (see `dinara_alignment`). */
struct Alignment {
    int64_t cost;       /* The alignment's cost. */
    int64_t score;      /* Its matches' reward less its cost, or minus its cost with no reward. */
    std::string cigar;  /* The CIGAR, its memory the string's own. */
    /* The aligned spans, `reference[reference_start:reference_end]` and `query[query_start:query_end]`. */
    int64_t reference_start, reference_end, query_start, query_end;
};

/* The wrapper's internals, called by the functions below it. */
namespace detail {
/* Raises the exception a DINARA_ code stands for; a result of zero or more, or DINARA_ABOVE_MAX, passes. */
inline void check(int64_t result) {
    if (result == DINARA_INVALID_COSTS) throw std::invalid_argument("dinara: a mismatch and an extension must cost");
    if (result == DINARA_INVALID_MODE) throw std::invalid_argument("dinara: a mode these costs cannot serve");
    if (result == DINARA_OUTSIDE_BAND) throw OutsideBand();
    if (result == DINARA_OUT_OF_MEMORY) throw std::bad_alloc();
    if (result < 0 && result != DINARA_ABOVE_MAX) throw UnsupportedSymbols();
}

/* The `dinara_costs` the C functions take for `costs`. */
inline dinara_costs c_costs(const Costs &costs) {
    return {costs.mismatch,          costs.opening,          costs.extension,
            costs.opening2,          costs.extension2,       costs.deletion_opening,
            costs.deletion_extension, costs.deletion_opening2, costs.deletion_extension2};
}

/* The least cost, or DINARA_ABOVE_MAX past a `max_cost` of zero or more; raises on any other code. */
inline int64_t distance(std::string_view reference, std::string_view query, const Costs &costs, const Mode &mode,
                        Band band, int64_t max_cost, dinara_aligner *aligner = nullptr) {
    dinara_costs c = c_costs(costs);
    dinara_options options{band.low, band.high, max_cost, 1, 0, 0};
    int64_t cost =
        aligner ? dinara_aligner_distance(aligner, reference.data(), static_cast<int64_t>(reference.size()),
                                          query.data(), static_cast<int64_t>(query.size()), &c, &mode.fields, &options)
                : dinara_distance(reference.data(), static_cast<int64_t>(reference.size()), query.data(),
                                  static_cast<int64_t>(query.size()), &c, &mode.fields, &options);
    check(cost);
    return cost;
}

/* An optimal alignment with its CIGAR copied and the C one freed, or nothing past a `max_cost` of zero
 * or more; raises on any other code. */
inline std::optional<Alignment> align(std::string_view reference, std::string_view query, const Costs &costs,
                                      const Mode &mode, Band band, int64_t max_cost, Ties ties, bool eqx,
                                      int64_t max_memory, dinara_aligner *aligner = nullptr) {
    dinara_costs c = c_costs(costs);
    dinara_options options{band.low, band.high, max_cost, eqx ? 1 : 0, ties == Ties::right ? 1 : 0, max_memory};
    dinara_alignment found{};
    int64_t status =
        aligner ? dinara_aligner_align(aligner, reference.data(), static_cast<int64_t>(reference.size()),
                                       query.data(), static_cast<int64_t>(query.size()), &c, &mode.fields, &options,
                                       &found)
                : dinara_align(reference.data(), static_cast<int64_t>(reference.size()), query.data(),
                               static_cast<int64_t>(query.size()), &c, &mode.fields, &options, &found);
    check(status);
    if (status == DINARA_ABOVE_MAX) return std::nullopt;
    Alignment result{found.cost,          found.score,         std::string(found.cigar, static_cast<size_t>(found.cigar_length)),
                     found.reference_start, found.reference_end, found.query_start, found.query_end};
    dinara_free(found.cigar);
    return result;
}

/* A batch's sequences as C takes them: each sequence's first byte and its length. */
struct Batch {
    std::vector<const char *> references, queries;          /* Each pair's first bytes. */
    std::vector<int64_t> reference_lengths, query_lengths; /* Each pair's lengths. */

    /* Pair `i` is `firsts[i]` against `seconds[i]`; the two must be as long. */
    Batch(const std::vector<std::string_view> &firsts, const std::vector<std::string_view> &seconds) {
        if (firsts.size() != seconds.size()) throw std::invalid_argument("dinara: a batch's sides differ in length");
        for (size_t index = 0; index < firsts.size(); ++index) {
            references.push_back(firsts[index].data());
            reference_lengths.push_back(static_cast<int64_t>(firsts[index].size()));
            queries.push_back(seconds[index].data());
            query_lengths.push_back(static_cast<int64_t>(seconds[index].size()));
        }
    }
    /* The number of pairs. */
    int64_t size() const { return static_cast<int64_t>(references.size()); }
};

/* Every pair's least cost, or its own DINARA_ code, which the caller checks; raises on a code failing the
 * whole batch. */
inline std::vector<int64_t> distances(const Batch &batch, const Costs &costs, const Mode &mode, Band band,
                                      int64_t max_cost, int threads) {
    dinara_costs c = c_costs(costs);
    dinara_options options{band.low, band.high, max_cost, 1, 0, 0};
    std::vector<int64_t> results(batch.references.size());
    check(dinara_distances(batch.size(), batch.references.data(), batch.reference_lengths.data(),
                           batch.queries.data(), batch.query_lengths.data(), &c, &mode.fields, &options, threads,
                           results.data()));
    return results;
}

/* Every pair's optimal alignment, or nothing for a pair past `max_cost`; with every CIGAR copied and
 * freed, raises as the first pair failing otherwise would. */
inline std::vector<std::optional<Alignment>> alignments(const Batch &batch, const Costs &costs, const Mode &mode,
                                                        Band band, int64_t max_cost, Ties ties, bool eqx,
                                                        int threads, int64_t max_memory) {
    dinara_costs c = c_costs(costs);
    dinara_options options{band.low, band.high, max_cost, eqx ? 1 : 0, ties == Ties::right ? 1 : 0, max_memory};
    std::vector<dinara_alignment> found(batch.references.size());
    std::vector<int64_t> statuses(batch.references.size());
    check(dinara_alignments(batch.size(), batch.references.data(), batch.reference_lengths.data(),
                            batch.queries.data(), batch.query_lengths.data(), &c, &mode.fields, &options, threads,
                            found.data(), statuses.data()));
    std::vector<std::optional<Alignment>> results;
    int64_t failed = 0;
    for (size_t index = 0; index < found.size(); ++index) {
        if (statuses[index] != 0) {
            if (statuses[index] != DINARA_ABOVE_MAX && !failed) failed = statuses[index];
            results.emplace_back();
            continue;
        }
        results.push_back(Alignment{found[index].cost, found[index].score,
                                    std::string(found[index].cigar, static_cast<size_t>(found[index].cigar_length)),
                                    found[index].reference_start, found[index].reference_end,
                                    found[index].query_start, found[index].query_end});
        dinara_free(found[index].cigar);
    }
    // Every CIGAR freed first, the first failure is raised as a single pair's would be.
    if (failed) check(failed);
    return results;
}
}  // namespace detail

/* The best score `align` would return, with no alignment traced. */
inline int64_t score(std::string_view reference, std::string_view query, const Costs &costs = Costs::edit(),
                     const Mode &mode = Mode::global(), Band band = {}) {
    dinara_costs c = detail::c_costs(costs);
    dinara_options options{band.low, band.high, -1, 1, 0, 0};
    int64_t found = 0;
    detail::check(dinara_score(reference.data(), static_cast<int64_t>(reference.size()), query.data(),
                               static_cast<int64_t>(query.size()), &c, &mode.fields, &options, &found));
    return found;
}

/* A local alignment's best score and end, and SSW's second best, more than `window` reference letters
 * away; a negative window is half the query and at least 15. */
inline dinara_local_scores_result local_scores(std::string_view reference, std::string_view query, const Costs &costs,
                                               const Mode &mode, int64_t window = -1) {
    dinara_costs c = detail::c_costs(costs);
    dinara_local_scores_result found{};
    detail::check(dinara_local_scores(reference.data(), static_cast<int64_t>(reference.size()), query.data(),
                                      static_cast<int64_t>(query.size()), &c, &mode.fields, window, &found));
    return found;
}

/* One hit of a `search`: the reference's place in the list, and its score. */
struct Hit {
    int64_t index; /* The reference's place in the list searched. */
    int64_t score; /* Its score, minus its cost with no reward. */
};

/* The query against every reference, the best first; `best` above zero keeps that many. */
inline std::vector<Hit> search(const std::vector<std::string_view> &references, std::string_view query,
                               const Costs &costs = Costs::edit(), const Mode &mode = Mode::global(),
                               int64_t best = 0, int threads = 0) {
    std::vector<const char *> texts;
    std::vector<int64_t> lengths;
    for (auto reference : references) texts.push_back(reference.data()), lengths.push_back(static_cast<int64_t>(reference.size()));
    std::vector<int64_t> indices(references.size()), scores(references.size());
    dinara_costs c = detail::c_costs(costs);
    int64_t found = dinara_search(static_cast<int64_t>(references.size()), texts.data(), lengths.data(), query.data(),
                                  static_cast<int64_t>(query.size()), &c, &mode.fields, nullptr, best, threads,
                                  indices.data(), scores.data());
    detail::check(found);
    std::vector<Hit> hits;
    for (int64_t rank = 0; rank < found; ++rank) hits.push_back({indices[rank], scores[rank]});
    return hits;
}

/* The least cost of aligning the query to the reference. */
inline int64_t distance(std::string_view reference, std::string_view query, const Costs &costs = Costs::edit(),
                        const Mode &mode = Mode::global(), Band band = {}) {
    return detail::distance(reference, query, costs, mode, band, -1);
}

/* The least cost, or nothing when it passes `max_cost`. */
inline std::optional<int64_t> distance_within(std::string_view reference, std::string_view query, int64_t max_cost,
                                              const Costs &costs = Costs::edit(), const Mode &mode = Mode::global(),
                                              Band band = {}) {
    if (max_cost < 0) return std::nullopt;
    int64_t cost = detail::distance(reference, query, costs, mode, band, max_cost);
    if (cost == DINARA_ABOVE_MAX) return std::nullopt;
    return cost;
}

/* An optimal alignment, its kept fronts within `max_memory` bytes when above zero. */
inline Alignment align(std::string_view reference, std::string_view query, const Costs &costs = Costs::edit(),
                       const Mode &mode = Mode::global(), Band band = {}, Ties ties = Ties::left,
                       bool eqx = true, int64_t max_memory = 0) {
    return *detail::align(reference, query, costs, mode, band, -1, ties, eqx, max_memory);
}

/* An optimal alignment, or nothing when its cost passes `max_cost`. */
inline std::optional<Alignment> align_within(std::string_view reference, std::string_view query, int64_t max_cost,
                                             const Costs &costs = Costs::edit(), const Mode &mode = Mode::global(),
                                             Band band = {}, Ties ties = Ties::left, bool eqx = true,
                                             int64_t max_memory = 0) {
    if (max_cost < 0) return std::nullopt;
    return detail::align(reference, query, costs, mode, band, max_cost, ties, eqx, max_memory);
}

/* Every pair's least cost, `references[i]` against `queries[i]`, over `threads` threads (zero: the
 * caller's own alone); raises as the first failing pair's `distance` would. */
inline std::vector<int64_t> distances(const std::vector<std::string_view> &references,
                                      const std::vector<std::string_view> &queries,
                                      const Costs &costs = Costs::edit(), const Mode &mode = Mode::global(),
                                      Band band = {}, int threads = 0) {
    std::vector<int64_t> results = detail::distances(detail::Batch(references, queries), costs, mode, band, -1, threads);
    for (int64_t result : results) detail::check(result);
    return results;
}

/* Every pair's least cost, or nothing for a pair past `max_cost`. */
inline std::vector<std::optional<int64_t>> distances_within(const std::vector<std::string_view> &references,
                                                            const std::vector<std::string_view> &queries,
                                                            int64_t max_cost, const Costs &costs = Costs::edit(),
                                                            const Mode &mode = Mode::global(), Band band = {},
                                                            int threads = 0) {
    std::vector<std::optional<int64_t>> results;
    if (max_cost < 0) return std::vector<std::optional<int64_t>>(references.size());
    for (int64_t result :
         detail::distances(detail::Batch(references, queries), costs, mode, band, max_cost, threads)) {
        detail::check(result);
        results.push_back(result == DINARA_ABOVE_MAX ? std::nullopt : std::optional<int64_t>(result));
    }
    return results;
}

/* Every pair's optimal alignment over `threads` threads; raises as the first failing pair's `align` would. */
inline std::vector<Alignment> alignments(const std::vector<std::string_view> &references,
                                         const std::vector<std::string_view> &queries,
                                         const Costs &costs = Costs::edit(), const Mode &mode = Mode::global(),
                                         Band band = {}, Ties ties = Ties::left, bool eqx = true,
                                         int threads = 0, int64_t max_memory = 0) {
    std::vector<Alignment> results;
    for (auto &found : detail::alignments(detail::Batch(references, queries), costs, mode, band, -1, ties, eqx,
                                          threads, max_memory))
        results.push_back(std::move(*found));
    return results;
}

/* Every pair's optimal alignment, or nothing for a pair whose cost passes `max_cost`. */
inline std::vector<std::optional<Alignment>> alignments_within(const std::vector<std::string_view> &references,
                                                               const std::vector<std::string_view> &queries,
                                                               int64_t max_cost, const Costs &costs = Costs::edit(),
                                                               const Mode &mode = Mode::global(), Band band = {},
                                                               Ties ties = Ties::left, bool eqx = true,
                                                               int threads = 0, int64_t max_memory = 0) {
    if (max_cost < 0) return std::vector<std::optional<Alignment>>(references.size());
    return detail::alignments(detail::Batch(references, queries), costs, mode, band, max_cost, ties, eqx,
                              threads, max_memory);
}

/* One thread's aligner: `distance` and `align` through memory kept from call to call, freed with it.
 * Movable, not copyable; one a thread, as the functions themselves keep no state between calls. */
class Aligner {
  public:
    Aligner() : handle_(dinara_aligner_new()) {
        if (!handle_) throw std::bad_alloc();
    }
    ~Aligner() { dinara_aligner_free(handle_); }
    Aligner(const Aligner &) = delete;
    Aligner &operator=(const Aligner &) = delete;
    Aligner(Aligner &&other) noexcept : handle_(other.handle_) { other.handle_ = nullptr; }
    Aligner &operator=(Aligner &&other) noexcept {
        std::swap(handle_, other.handle_);
        return *this;
    }

    /* `dinara::distance` through this aligner. */
    int64_t distance(std::string_view reference, std::string_view query, const Costs &costs = Costs::edit(),
                     const Mode &mode = Mode::global(), Band band = {}) {
        return detail::distance(reference, query, costs, mode, band, -1, handle_);
    }

    /* `dinara::distance_within` through this aligner. */
    std::optional<int64_t> distance_within(std::string_view reference, std::string_view query, int64_t max_cost,
                                           const Costs &costs = Costs::edit(), const Mode &mode = Mode::global(),
                                           Band band = {}) {
        if (max_cost < 0) return std::nullopt;
        int64_t cost = detail::distance(reference, query, costs, mode, band, max_cost, handle_);
        if (cost == DINARA_ABOVE_MAX) return std::nullopt;
        return cost;
    }

    /* `dinara::align` through this aligner. */
    Alignment align(std::string_view reference, std::string_view query, const Costs &costs = Costs::edit(),
                    const Mode &mode = Mode::global(), Band band = {}, Ties ties = Ties::left, bool eqx = true,
                    int64_t max_memory = 0) {
        return *detail::align(reference, query, costs, mode, band, -1, ties, eqx, max_memory, handle_);
    }

    /* `dinara::align_within` through this aligner. */
    std::optional<Alignment> align_within(std::string_view reference, std::string_view query, int64_t max_cost,
                                          const Costs &costs = Costs::edit(), const Mode &mode = Mode::global(),
                                          Band band = {}, Ties ties = Ties::left, bool eqx = true,
                                          int64_t max_memory = 0) {
        if (max_cost < 0) return std::nullopt;
        return detail::align(reference, query, costs, mode, band, max_cost, ties, eqx, max_memory, handle_);
    }

  private:
    dinara_aligner *handle_;
};

}  // namespace dinara
#endif

#endif /* DINARA_H */
