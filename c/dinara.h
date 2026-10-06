/*
 * dinara-align from C and C++: the unit-cost edit distance between two sequences, and the least
 * gap-affine cost as WFA counts it, each with an optimal alignment as a CIGAR. Build the library with `pixi run build-c [target-cpu]`, which leaves it in
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
 * gap of `k` letters `opening + k * extension`, or a DINARA_ code, and an optimal alignment's CIGAR,
 * returned as `dinara_edit_cigar` returns it. By a wavefront from both ends: its time grows with the
 * square of the cost, and its memory stays bounded.
 */
int64_t dinara_affine_cigar(const char *first, int64_t first_length, const char *second, int64_t second_length,
                            int64_t mismatch, int64_t opening, int64_t extension, int extended, char **cigar,
                            int64_t *cigar_length);

/* Frees a CIGAR that `dinara_edit_cigar` or `dinara_affine_cigar` returned. */
void dinara_free(char *cigar);

#ifdef __cplusplus
}

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

/* Gap-affine costs as WFA counts them: a substitution `mismatch`, a gap of `k` letters `opening + k * extension`. */
inline AffineAlignment affine_cigar(std::string_view first, std::string_view second, int64_t mismatch, int64_t opening,
                                    int64_t extension, bool extended = true) {
    char *text = nullptr;
    int64_t length = 0;
    int64_t cost = dinara_affine_cigar(first.data(), static_cast<int64_t>(first.size()), second.data(),
                                       static_cast<int64_t>(second.size()), mismatch, opening, extension,
                                       extended ? 1 : 0, &text, &length);
    if (cost == DINARA_INVALID_COSTS) throw std::invalid_argument("dinara: a mismatch and an extension must cost");
    if (cost < 0) throw UnsupportedSymbols();
    AffineAlignment result{cost, std::string(text, static_cast<size_t>(length))};
    dinara_free(text);
    return result;
}

}  // namespace dinara
#endif

#endif /* DINARA_H */
