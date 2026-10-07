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
/* The CIGAR's memory could not be allocated. */
#define DINARA_OUT_OF_MEMORY (-2)
/* Costs no alignment can be searched by: a mismatch or extension of zero, or a negative cost. */
#define DINARA_INVALID_COSTS (-3)
/* Every alignment costs more than the `max_cost` asked for, or under a cap none fits the band. */
#define DINARA_ABOVE_MAX (-4)
/* No alignment stays inside the band of diagonals asked for. */
#define DINARA_OUTSIDE_BAND (-5)
/* A mode that cannot serve what was asked: the least cost of an extension or a local alignment, which
 * maximize a score, a cap on either, or a band on a local alignment. */
#define DINARA_INVALID_MODE (-6)

/*
 * What each edit costs: a substitution `mismatch`, a gap of `k` letters `opening + k * extension`, or
 * with `opening2` zero or more the less of that and `opening2 + k * extension2`. A null pointer is unit
 * costs, the edit distance, {1, 0, 1, -1, 0}.
 */
typedef struct {
    int64_t mismatch;
    int64_t opening;
    int64_t extension;
    int64_t opening2;
    int64_t extension2;
} dinara_costs;

#define DINARA_ENDS_FREE 0
#define DINARA_EXTENSION 1
#define DINARA_LOCAL 2
/* As many free letters as any sequence has. */
#define DINARA_ALL INT64_MAX

/*
 * Which alignments count. DINARA_ENDS_FREE: up to so many letters at each end of each sequence left
 * unaligned for nothing; all zero is a global alignment, the reference's two at DINARA_ALL a query
 * placed anywhere in it (infix), its end alone a prefix. DINARA_EXTENSION: fixed at both sequences'
 * starts, or with `anchor` nonzero their ends, and free to stop anywhere, a match earning
 * `match_score`: a seed's extension, as KSW2's without Z-drop. DINARA_LOCAL: any part of each, a match
 * earning `match_score`, Smith-Waterman, as abPOA's local mode. A null pointer is a global alignment.
 */
typedef struct {
    int64_t kind;
    int64_t reference_start;
    int64_t reference_end;
    int64_t query_start;
    int64_t query_end;
    int64_t match_score;
    int64_t anchor;
} dinara_mode;

/*
 * Every move stays on the diagonals `band_low ..= band_high`, a diagonal being the reference's letters
 * aligned or skipped less the query's, from the alignment's origin: the exact optimum over the
 * alignments inside; INT64_MIN and INT64_MAX for no band, KSW2's band of width `w` `-w ..= w`. A
 * `max_cost` of zero or more caps the cost, which the search proves past after about half of it.
 * `extended` nonzero writes `=` and `X`, else `M`; of equally good alignments the CIGAR places indels
 * left, as minimap2 does, or with `right_ties` nonzero right, WFA2-lib's CIGAR byte for byte. A null
 * pointer is no band, no cap, `=` and `X`, indels left.
 */
typedef struct {
    int64_t band_low;
    int64_t band_high;
    int64_t max_cost;
    int64_t extended;
    int64_t right_ties;
} dinara_options;

/*
 * An optimal alignment: its cost, its score (an extension's or a local alignment's matches' reward
 * less the cost, else minus the cost), the spans it aligns, and its CIGAR, NUL-terminated and `cigar_length` bytes long, which the
 * caller frees with `dinara_free`.
 */
typedef struct {
    int64_t cost;
    int64_t score;
    int64_t reference_start;
    int64_t reference_end;
    int64_t query_start;
    int64_t query_end;
    char *cigar;
    int64_t cigar_length;
} dinara_alignment;

/* The least cost of aligning the query to the reference, or a DINARA_ code. */
int64_t dinara_distance(const char *reference, int64_t reference_length, const char *query, int64_t query_length,
                        const dinara_costs *costs, const dinara_mode *mode, const dinara_options *options);

/* An optimal alignment into `*alignment`; zero, or a DINARA_ code and then no CIGAR to free. */
int64_t dinara_align(const char *reference, int64_t reference_length, const char *query, int64_t query_length,
                     const dinara_costs *costs, const dinara_mode *mode, const dinara_options *options,
                     dinara_alignment *alignment);

/* Frees a CIGAR that `dinara_align` returned; null is nothing to free. */
void dinara_free(char *cigar);

#ifdef __cplusplus
}

#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

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
    int64_t mismatch = 1, opening = 0, extension = 1, opening2 = -1, extension2 = 0;
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
};

/* Which end of both sequences an extension is fixed at. */
enum class Anchor { start, end };

/* Which alignments count (see `dinara_mode`). */
struct Mode {
    dinara_mode fields{DINARA_ENDS_FREE, 0, 0, 0, 0, 0, 0};
    /* Both sequences end to end. */
    static Mode global() { return {}; }
    /* The whole query against wherever in the reference it fits best. */
    static Mode infix() { return ends_free(DINARA_ALL, DINARA_ALL, 0, 0); }
    /* The whole query against the reference's best prefix. */
    static Mode prefix() { return ends_free(0, DINARA_ALL, 0, 0); }
    /* The whole query against the reference's best suffix. */
    static Mode suffix() { return ends_free(DINARA_ALL, 0, 0, 0); }
    /* Up to so many letters at each end of each sequence left unaligned for nothing. */
    static Mode ends_free(int64_t reference_start, int64_t reference_end, int64_t query_start, int64_t query_end) {
        return {{DINARA_ENDS_FREE, reference_start, reference_end, query_start, query_end, 0, 0}};
    }
    /* The best-scoring alignment from one end, a match earning `match_score`. */
    static Mode extension(int64_t match_score, Anchor anchor = Anchor::start) {
        return {{DINARA_EXTENSION, 0, 0, 0, 0, match_score, anchor == Anchor::end ? 1 : 0}};
    }
    /* The best-scoring alignment of any part of each, a match earning `match_score`: Smith-Waterman. */
    static Mode local(int64_t match_score) { return {{DINARA_LOCAL, 0, 0, 0, 0, match_score, 0}}; }
};

/* The diagonals every move stays on, `low ..= high`; the default is every diagonal. */
struct Band {
    int64_t low = INT64_MIN;
    int64_t high = INT64_MAX;
    /* KSW2's band of width `w`: at most `w` diagonals from the origin's, either way. */
    static Band around(int64_t width) { return Band{-width, width}; }
};

/* Which of several equally good alignments a CIGAR spells: indels placed left, as minimap2 places them,
 * or right, WFA2-lib's CIGAR byte for byte. */
enum class Ties { left, right };

/* An optimal alignment (see `dinara_alignment`). */
struct Alignment {
    int64_t cost;
    int64_t score;
    std::string cigar;
    int64_t reference_start, reference_end, query_start, query_end;
};

namespace detail {
inline void check(int64_t result) {
    if (result == DINARA_INVALID_COSTS) throw std::invalid_argument("dinara: a mismatch and an extension must cost");
    if (result == DINARA_INVALID_MODE) throw std::invalid_argument("dinara: a mode these costs cannot serve");
    if (result == DINARA_OUTSIDE_BAND) throw OutsideBand();
    if (result == DINARA_OUT_OF_MEMORY) throw std::bad_alloc();
    if (result < 0 && result != DINARA_ABOVE_MAX) throw UnsupportedSymbols();
}

inline dinara_costs c_costs(const Costs &costs) {
    return {costs.mismatch, costs.opening, costs.extension, costs.opening2, costs.extension2};
}

inline int64_t distance(std::string_view reference, std::string_view query, const Costs &costs, const Mode &mode,
                        Band band, int64_t max_cost) {
    dinara_costs c = c_costs(costs);
    dinara_options options{band.low, band.high, max_cost, 1, 0};
    int64_t cost = dinara_distance(reference.data(), static_cast<int64_t>(reference.size()), query.data(),
                                   static_cast<int64_t>(query.size()), &c, &mode.fields, &options);
    check(cost);
    return cost;
}

inline std::optional<Alignment> align(std::string_view reference, std::string_view query, const Costs &costs,
                                      const Mode &mode, Band band, int64_t max_cost, Ties ties, bool extended) {
    dinara_costs c = c_costs(costs);
    dinara_options options{band.low, band.high, max_cost, extended ? 1 : 0, ties == Ties::right ? 1 : 0};
    dinara_alignment found{};
    int64_t status = dinara_align(reference.data(), static_cast<int64_t>(reference.size()), query.data(),
                                  static_cast<int64_t>(query.size()), &c, &mode.fields, &options, &found);
    check(status);
    if (status == DINARA_ABOVE_MAX) return std::nullopt;
    Alignment result{found.cost,          found.score,         std::string(found.cigar, static_cast<size_t>(found.cigar_length)),
                     found.reference_start, found.reference_end, found.query_start, found.query_end};
    dinara_free(found.cigar);
    return result;
}
}  // namespace detail

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

/* An optimal alignment. */
inline Alignment align(std::string_view reference, std::string_view query, const Costs &costs = Costs::edit(),
                       const Mode &mode = Mode::global(), Band band = {}, Ties ties = Ties::left,
                       bool extended = true) {
    return *detail::align(reference, query, costs, mode, band, -1, ties, extended);
}

/* An optimal alignment, or nothing when its cost passes `max_cost`. */
inline std::optional<Alignment> align_within(std::string_view reference, std::string_view query, int64_t max_cost,
                                             const Costs &costs = Costs::edit(), const Mode &mode = Mode::global(),
                                             Band band = {}, Ties ties = Ties::left, bool extended = true) {
    if (max_cost < 0) return std::nullopt;
    return detail::align(reference, query, costs, mode, band, max_cost, ties, extended);
}

}  // namespace dinara
#endif

#endif /* DINARA_H */
