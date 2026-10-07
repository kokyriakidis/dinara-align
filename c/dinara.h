/*
 * dinara-align from C and C++: the unit-cost edit distance between two sequences, and the least
 * gap-affine cost as WFA counts it, global, ends-free or banded, each with an optimal alignment as a
 * CIGAR, and the best extension from one end. Every result is exact. Build the library with `pixi run build-c [target-cpu]`, which leaves it in
 * build/c beside the Mojo runtime libraries it loads and this header; link with `-Lbuild/c -ldinara`
 * and put build/c on the program's library search path (an rpath, say). build/c also holds the
 * libstdc++ the runtime loads, which a C++ program then shares: it is GCC 15's, as new as any compiler
 * the program is likely built with. Built for an explicit `target-cpu`, the library runs only on CPUs
 * with that instruction set; the default is the platform's oldest (x86-64, apple-m1).
 *
 * Sequences are bytes. For the edit distance: `A`, `C`, `G` and `T`, and up to four other symbols
 * between the two, each matching only itself. For the affine cost, every byte matches only itself,
 * save 0xFE and 0xFF, which UTF-8 never holds. The functions keep no state between calls, so they may run on many threads at
 * once. Failures come back as a negative result, one of the `DINARA_` codes below.
 *
 * The CIGAR takes the first sequence as the reference: `=` a match, `X` a substitution (or `M` for
 * either), `D` a base of the first sequence alone, `I` a base of the second alone.
 */
#ifndef DINARA_H
#define DINARA_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* More than four symbols past ACGT between the two sequences, or a 0xFE or 0xFF byte for the affine cost. */
#define DINARA_UNSUPPORTED_SYMBOLS (-1)
/* Affine costs no alignment can be searched by: a mismatch or extension of zero, or a negative cost. */
#define DINARA_INVALID_COSTS (-3)
/* Every alignment costs more than the `max_cost` asked for. */
#define DINARA_ABOVE_MAX (-4)
/* No alignment stays inside the band of diagonals asked for. */
#define DINARA_OUTSIDE_BAND (-5)

/* The global edit distance, or a DINARA_ code. */
int64_t dinara_edit_distance(const char *first, int64_t first_length, const char *second, int64_t second_length);

/*
 * The global edit distance, or a DINARA_ code, and an optimal alignment's CIGAR: `=` and `X` when
 * `extended` is nonzero, else `M`. On success `*cigar` holds the CIGAR, NUL-terminated and
 * `*cigar_length` bytes long, which the caller frees with `dinara_free`.
 */
int64_t dinara_edit_cigar(const char *first, int64_t first_length, const char *second, int64_t second_length,
                          int extended, char **cigar, int64_t *cigar_length);

/*
 * The least global cost under gap-affine costs as WFA counts them, a substitution `mismatch` and a
 * gap of `k` letters `opening + k * extension`, with no alignment, or a DINARA_ code: DINARA_ABOVE_MAX
 * when the cost passes a `max_cost` of zero or more, which the search proves after about half of it.
 * A negative `max_cost` caps nothing. The four `free` counts are the letters at each end of each
 * sequence that may go unaligned for nothing, as WFA2-lib's ends-free alignment counts them: all zero
 * for a global alignment, the first sequence's two at its length to place the second anywhere in it.
 * Every move stays on the diagonals `band_low ..= band_high`, a diagonal being the letters of the first
 * sequence aligned or skipped less those of the second: the exact optimum over the alignments inside,
 * DINARA_OUTSIDE_BAND when none is (DINARA_ABOVE_MAX under a cap); INT64_MIN and INT64_MAX for no band.
 * KSW2's band of width `w` is `-w ..= w`. By a wavefront from both ends keeping a few fronts.
 */
int64_t dinara_affine_distance(const char *first, int64_t first_length, const char *second, int64_t second_length,
                               int64_t mismatch, int64_t opening, int64_t extension, int64_t max_cost,
                               int64_t first_begin_free, int64_t first_end_free, int64_t second_begin_free,
                               int64_t second_end_free, int64_t band_low, int64_t band_high);

/*
 * `dinara_affine_distance`'s cost, or a DINARA_ code, and an optimal alignment's CIGAR, returned as
 * `dinara_edit_cigar` returns it; no CIGAR for DINARA_ABOVE_MAX; free letters as `D` and `I` runs.
 * Its time grows with the square of the cost, and its memory stays bounded.
 */
int64_t dinara_affine_cigar(const char *first, int64_t first_length, const char *second, int64_t second_length,
                            int64_t mismatch, int64_t opening, int64_t extension, int64_t max_cost,
                            int64_t first_begin_free, int64_t first_end_free, int64_t second_begin_free,
                            int64_t second_end_free, int64_t band_low, int64_t band_high, int extended, char **cigar,
                            int64_t *cigar_length);

/*
 * `dinara_affine_distance` under two-piece gap-affine costs, WFA's gap-affine-2p: a gap of `k`
 * letters costs the less of `opening1 + k * extension1` and `opening2 + k * extension2`.
 */
int64_t dinara_affine2p_distance(const char *first, int64_t first_length, const char *second, int64_t second_length,
                                 int64_t mismatch, int64_t opening1, int64_t extension1, int64_t opening2,
                                 int64_t extension2, int64_t max_cost, int64_t first_begin_free,
                                 int64_t first_end_free, int64_t second_begin_free, int64_t second_end_free,
                                 int64_t band_low, int64_t band_high);

/* `dinara_affine_cigar` under two-piece gap-affine costs, as `dinara_affine2p_distance` counts them. */
int64_t dinara_affine2p_cigar(const char *first, int64_t first_length, const char *second, int64_t second_length,
                              int64_t mismatch, int64_t opening1, int64_t extension1, int64_t opening2,
                              int64_t extension2, int64_t max_cost, int64_t first_begin_free, int64_t first_end_free,
                              int64_t second_begin_free, int64_t second_end_free, int64_t band_low,
                              int64_t band_high, int extended, char **cigar, int64_t *cigar_length);

/*
 * The best score of an alignment fixed at both sequences' starts, or with `at_end` nonzero their ends,
 * and free to stop anywhere: a seed's extension, as KSW2's extension without Z-drop, exact. A match
 * earns `match_score` (zero or more), the costs as for `dinara_affine_cigar`; aligning nothing scores
 * zero. The letters of each sequence covered from that end go to `*first_covered` and
 * `*second_covered`, and the CIGAR over them is returned as `dinara_edit_cigar` returns it. The band
 * counts diagonals from the anchor. A negative result is a DINARA_ code.
 */
int64_t dinara_affine_extension(const char *first, int64_t first_length, const char *second, int64_t second_length,
                                int64_t match_score, int64_t mismatch, int64_t opening, int64_t extension, int at_end,
                                int64_t band_low, int64_t band_high, int extended, char **cigar,
                                int64_t *cigar_length, int64_t *first_covered, int64_t *second_covered);

/* `dinara_affine_extension` under two-piece gap-affine costs. */
int64_t dinara_affine2p_extension(const char *first, int64_t first_length, const char *second,
                                  int64_t second_length, int64_t match_score, int64_t mismatch, int64_t opening1,
                                  int64_t extension1, int64_t opening2, int64_t extension2, int at_end,
                                  int64_t band_low, int64_t band_high, int extended, char **cigar,
                                  int64_t *cigar_length, int64_t *first_covered, int64_t *second_covered);

/* Frees a CIGAR that any function here returned. */
void dinara_free(char *cigar);

#ifdef __cplusplus
}

#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace dinara {

/* The global edit distance and an optimal alignment as a CIGAR. */
struct Alignment {
    int64_t distance;
    std::string cigar;
};

/* Raised on sequences holding more than four symbols past ACGT between them. */
struct UnsupportedSymbols : std::invalid_argument {
    UnsupportedSymbols() : std::invalid_argument("dinara: more than four symbols past ACGT") {}
};

inline int64_t edit_distance(std::string_view first, std::string_view second) {
    int64_t distance = dinara_edit_distance(first.data(), static_cast<int64_t>(first.size()), second.data(),
                                            static_cast<int64_t>(second.size()));
    if (distance < 0) throw UnsupportedSymbols();
    return distance;
}

/* The least gap-affine cost and an optimal alignment as a CIGAR. */
struct AffineAlignment {
    int64_t cost;
    std::string cigar;
};

/* `extended` writes `=` and `X`; false writes `M` for both. */
inline Alignment edit_cigar(std::string_view first, std::string_view second, bool extended = true) {
    char *text = nullptr;
    int64_t length = 0;
    int64_t distance = dinara_edit_cigar(first.data(), static_cast<int64_t>(first.size()), second.data(),
                                         static_cast<int64_t>(second.size()), extended ? 1 : 0, &text, &length);
    if (distance < 0) throw UnsupportedSymbols();
    Alignment result{distance, std::string(text, static_cast<size_t>(length))};
    dinara_free(text);
    return result;
}

/* Letters at each end of each sequence an alignment may leave unaligned for nothing; all zero is global. */
struct EndsFree {
    int64_t first_begin = 0;
    int64_t first_end = 0;
    int64_t second_begin = 0;
    int64_t second_end = 0;
};

/* A second gap piece: a gap of `k` letters then costs the less of the two pieces' `opening + k * extension`. */
struct SecondPiece {
    int64_t opening;
    int64_t extension;
};

/* The diagonals every move stays on, `low ..= high`; the default is every diagonal. */
struct Band {
    int64_t low = INT64_MIN;
    int64_t high = INT64_MAX;
    /* KSW2's band of width `w`: at most `w` diagonals from the origin's, either way. */
    static Band around(int64_t width) { return Band{-width, width}; }
};

/* Raised when no alignment stays inside the band asked for. */
struct OutsideBand : std::invalid_argument {
    OutsideBand() : std::invalid_argument("dinara: no alignment stays inside the band") {}
};

/* Which end of both sequences an extension is fixed at. */
enum class Anchor { start, end };

/* The best extension's score, the letters of each sequence it covers from its anchor, and its CIGAR. */
struct Extension {
    int64_t score;
    int64_t first_length;
    int64_t second_length;
    std::string cigar;
};

namespace detail {
inline void check_affine(int64_t result) {
    if (result == DINARA_INVALID_COSTS) throw std::invalid_argument("dinara: a mismatch and an extension must cost");
    if (result == DINARA_OUTSIDE_BAND) throw OutsideBand();
    if (result < 0 && result != DINARA_ABOVE_MAX) throw UnsupportedSymbols();
}

inline int64_t affine_distance(std::string_view first, std::string_view second, int64_t mismatch, int64_t opening,
                               int64_t extension, std::optional<SecondPiece> piece, int64_t max_cost, EndsFree ends,
                               Band band) {
    int64_t cost =
        piece ? dinara_affine2p_distance(first.data(), static_cast<int64_t>(first.size()), second.data(),
                                         static_cast<int64_t>(second.size()), mismatch, opening, extension,
                                         piece->opening, piece->extension, max_cost, ends.first_begin, ends.first_end,
                                         ends.second_begin, ends.second_end, band.low, band.high)
              : dinara_affine_distance(first.data(), static_cast<int64_t>(first.size()), second.data(),
                                       static_cast<int64_t>(second.size()), mismatch, opening, extension, max_cost,
                                       ends.first_begin, ends.first_end, ends.second_begin, ends.second_end, band.low,
                                       band.high);
    check_affine(cost);
    return cost;
}

inline std::optional<AffineAlignment> affine_cigar(std::string_view first, std::string_view second, int64_t mismatch,
                                                   int64_t opening, int64_t extension, std::optional<SecondPiece> piece,
                                                   int64_t max_cost, bool extended, EndsFree ends, Band band) {
    char *text = nullptr;
    int64_t length = 0;
    int64_t cost =
        piece ? dinara_affine2p_cigar(first.data(), static_cast<int64_t>(first.size()), second.data(),
                                      static_cast<int64_t>(second.size()), mismatch, opening, extension,
                                      piece->opening, piece->extension, max_cost, ends.first_begin, ends.first_end,
                                      ends.second_begin, ends.second_end, band.low, band.high, extended ? 1 : 0, &text,
                                      &length)
              : dinara_affine_cigar(first.data(), static_cast<int64_t>(first.size()), second.data(),
                                    static_cast<int64_t>(second.size()), mismatch, opening, extension, max_cost,
                                    ends.first_begin, ends.first_end, ends.second_begin, ends.second_end, band.low,
                                    band.high, extended ? 1 : 0, &text, &length);
    check_affine(cost);
    if (cost == DINARA_ABOVE_MAX) return std::nullopt;
    AffineAlignment result{cost, std::string(text, static_cast<size_t>(length))};
    dinara_free(text);
    return result;
}

inline Extension extension(std::string_view first, std::string_view second, int64_t match_score, int64_t mismatch,
                           int64_t opening, int64_t extension, std::optional<SecondPiece> piece, Anchor anchor,
                           Band band, bool extended) {
    char *text = nullptr;
    int64_t length = 0, first_covered = 0, second_covered = 0;
    int at_end = anchor == Anchor::end ? 1 : 0;
    int64_t score =
        piece ? dinara_affine2p_extension(first.data(), static_cast<int64_t>(first.size()), second.data(),
                                          static_cast<int64_t>(second.size()), match_score, mismatch, opening,
                                          extension, piece->opening, piece->extension, at_end, band.low, band.high,
                                          extended ? 1 : 0, &text, &length, &first_covered, &second_covered)
              : dinara_affine_extension(first.data(), static_cast<int64_t>(first.size()), second.data(),
                                        static_cast<int64_t>(second.size()), match_score, mismatch, opening, extension,
                                        at_end, band.low, band.high, extended ? 1 : 0, &text, &length, &first_covered,
                                        &second_covered);
    check_affine(score);
    Extension result{score, first_covered, second_covered, std::string(text, static_cast<size_t>(length))};
    dinara_free(text);
    return result;
}
}  // namespace detail

/* Gap-affine costs as WFA counts them: a substitution `mismatch`, a gap of `k` letters `opening + k * extension`. */
inline int64_t affine_distance(std::string_view first, std::string_view second, int64_t mismatch, int64_t opening,
                               int64_t extension, EndsFree ends = {}, Band band = {}) {
    return detail::affine_distance(first, second, mismatch, opening, extension, std::nullopt, -1, ends, band);
}

/* The cost, or nothing when it passes `max_cost`, which the search proves after about half of it. */
inline std::optional<int64_t> affine_distance_within(std::string_view first, std::string_view second,
                                                     int64_t mismatch, int64_t opening, int64_t extension,
                                                     int64_t max_cost, EndsFree ends = {}, Band band = {}) {
    if (max_cost < 0) return std::nullopt;
    int64_t cost =
        detail::affine_distance(first, second, mismatch, opening, extension, std::nullopt, max_cost, ends, band);
    if (cost == DINARA_ABOVE_MAX) return std::nullopt;
    return cost;
}

/* The cost and an optimal alignment's CIGAR. */
inline AffineAlignment affine_cigar(std::string_view first, std::string_view second, int64_t mismatch, int64_t opening,
                                    int64_t extension, bool extended = true, EndsFree ends = {}, Band band = {}) {
    return *detail::affine_cigar(first, second, mismatch, opening, extension, std::nullopt, -1, extended, ends, band);
}

/* The cost and an optimal alignment's CIGAR, or nothing when the cost passes `max_cost`. */
inline std::optional<AffineAlignment> affine_cigar_within(std::string_view first, std::string_view second,
                                                          int64_t mismatch, int64_t opening, int64_t extension,
                                                          int64_t max_cost, bool extended = true,
                                                          EndsFree ends = {}, Band band = {}) {
    if (max_cost < 0) return std::nullopt;
    return detail::affine_cigar(first, second, mismatch, opening, extension, std::nullopt, max_cost, extended, ends,
                                band);
}

/* Two-piece gap-affine costs, WFA's gap-affine-2p: each of the above with a second gap piece. */
inline int64_t affine2p_distance(std::string_view first, std::string_view second, int64_t mismatch, int64_t opening,
                                 int64_t extension, SecondPiece piece, EndsFree ends = {}, Band band = {}) {
    return detail::affine_distance(first, second, mismatch, opening, extension, piece, -1, ends, band);
}

inline std::optional<int64_t> affine2p_distance_within(std::string_view first, std::string_view second,
                                                       int64_t mismatch, int64_t opening, int64_t extension,
                                                       SecondPiece piece, int64_t max_cost, EndsFree ends = {},
                                                       Band band = {}) {
    if (max_cost < 0) return std::nullopt;
    int64_t cost = detail::affine_distance(first, second, mismatch, opening, extension, piece, max_cost, ends, band);
    if (cost == DINARA_ABOVE_MAX) return std::nullopt;
    return cost;
}

inline AffineAlignment affine2p_cigar(std::string_view first, std::string_view second, int64_t mismatch,
                                      int64_t opening, int64_t extension, SecondPiece piece, bool extended = true,
                                      EndsFree ends = {}, Band band = {}) {
    return *detail::affine_cigar(first, second, mismatch, opening, extension, piece, -1, extended, ends, band);
}

inline std::optional<AffineAlignment> affine2p_cigar_within(std::string_view first, std::string_view second,
                                                            int64_t mismatch, int64_t opening, int64_t extension,
                                                            SecondPiece piece, int64_t max_cost, bool extended = true,
                                                            EndsFree ends = {}, Band band = {}) {
    if (max_cost < 0) return std::nullopt;
    return detail::affine_cigar(first, second, mismatch, opening, extension, piece, max_cost, extended, ends, band);
}

/*
 * The best-scoring alignment fixed at one end of both sequences and free to stop anywhere, a match
 * earning `match_score`: a seed's extension, as KSW2's without Z-drop, exact.
 */
inline Extension affine_extension(std::string_view first, std::string_view second, int64_t match_score,
                                  int64_t mismatch, int64_t opening, int64_t extension, Anchor anchor = Anchor::start,
                                  Band band = {}, bool extended = true) {
    return detail::extension(first, second, match_score, mismatch, opening, extension, std::nullopt, anchor, band,
                             extended);
}

inline Extension affine2p_extension(std::string_view first, std::string_view second, int64_t match_score,
                                    int64_t mismatch, int64_t opening, int64_t extension, SecondPiece piece,
                                    Anchor anchor = Anchor::start, Band band = {}, bool extended = true) {
    return detail::extension(first, second, match_score, mismatch, opening, extension, piece, anchor, band, extended);
}

}  // namespace dinara
#endif

#endif /* DINARA_H */
