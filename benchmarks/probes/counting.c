// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
// MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
//
// Counts a process's allocations and bulk copies, preloaded under it on Linux (`LD_PRELOAD`), for
// `benchmarks/profile.py probes`: every call of malloc, calloc, realloc, the aligned allocators and free,
// and of the Mojo runtime's own allocator, `KGEN_CompilerRT_AlignedAlloc` and its free, which Mojo's lists
// and strings take and which reaches no malloc, with the bytes asked for, and every call of memcpy, memmove and memset that reaches the C library, with
// the bytes moved. Copies the compiler inlined, short ones of a known size, never reach it and are not
// counted. The runner built with probes finds `dinara_counts` by `dlsym` and reads it before and after
// each pair; a process without this library finds nothing and reports no counts.
//
//     cc -O2 -shared -fPIC benchmarks/probes/counting.c -o counting.so -ldl
#define _GNU_SOURCE
#include <dlfcn.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

enum {
    ALLOCATIONS, ALLOCATED, FREES, REALLOCATIONS, COPIES, COPIED, MOVES, MOVED, SETS, SET, COUNTS
};

// Indexed by the enum above; relaxed atomics, since a batch may allocate on several threads.
uint64_t dinara_counts[COUNTS];

#define COUNT(slot, amount) __atomic_fetch_add(&dinara_counts[slot], (uint64_t)(amount), __ATOMIC_RELAXED)

static void *(*real_malloc)(size_t);
static void *(*real_calloc)(size_t, size_t);
static void *(*real_realloc)(void *, size_t);
static void (*real_free)(void *);
static int (*real_posix_memalign)(void **, size_t, size_t);
static void *(*real_aligned_alloc)(size_t, size_t);
static void *(*real_memalign)(size_t, size_t);
static void *(*real_memcpy)(void *, const void *, size_t);
static void *(*real_memmove)(void *, const void *, size_t);
static void *(*real_memset)(void *, int, size_t);
static void *(*real_kgen_alloc)(size_t, size_t);
static void (*real_kgen_free)(void *);

// dlsym may itself call calloc before `real_calloc` is known: those few calls take this arena.
static char early[8192];
static size_t early_used;
static int resolving;

static void resolve(void) {
    if (real_malloc || resolving)
        return;
    resolving = 1;
    real_malloc = dlsym(RTLD_NEXT, "malloc");
    real_calloc = dlsym(RTLD_NEXT, "calloc");
    real_realloc = dlsym(RTLD_NEXT, "realloc");
    real_free = dlsym(RTLD_NEXT, "free");
    real_posix_memalign = dlsym(RTLD_NEXT, "posix_memalign");
    real_aligned_alloc = dlsym(RTLD_NEXT, "aligned_alloc");
    real_memalign = dlsym(RTLD_NEXT, "memalign");
    real_memcpy = dlsym(RTLD_NEXT, "memcpy");
    real_memmove = dlsym(RTLD_NEXT, "memmove");
    real_memset = dlsym(RTLD_NEXT, "memset");
    real_kgen_alloc = dlsym(RTLD_NEXT, "KGEN_CompilerRT_AlignedAlloc");
    real_kgen_free = dlsym(RTLD_NEXT, "KGEN_CompilerRT_AlignedFree");
    resolving = 0;
}

static int is_early(void *pointer) {
    return (char *)pointer >= early && (char *)pointer < early + sizeof early;
}

void *malloc(size_t size) {
    resolve();
    COUNT(ALLOCATIONS, 1);
    COUNT(ALLOCATED, size);
    return real_malloc(size);
}

void *calloc(size_t count, size_t size) {
    if (!real_calloc) {
        if (resolving) {
            size_t bytes = (count * size + 15) & ~(size_t)15;
            if (early_used + bytes > sizeof early)
                return NULL;
            void *pointer = early + early_used;
            early_used += bytes;
            return pointer;  // static memory, already zero
        }
        resolve();
    }
    COUNT(ALLOCATIONS, 1);
    COUNT(ALLOCATED, count * size);
    return real_calloc(count, size);
}

void *realloc(void *pointer, size_t size) {
    resolve();
    COUNT(REALLOCATIONS, 1);
    COUNT(ALLOCATED, size);
    return real_realloc(pointer, size);
}

void free(void *pointer) {
    if (!pointer || is_early(pointer))
        return;
    resolve();
    COUNT(FREES, 1);
    real_free(pointer);
}

int posix_memalign(void **out, size_t alignment, size_t size) {
    resolve();
    COUNT(ALLOCATIONS, 1);
    COUNT(ALLOCATED, size);
    return real_posix_memalign(out, alignment, size);
}

void *aligned_alloc(size_t alignment, size_t size) {
    resolve();
    COUNT(ALLOCATIONS, 1);
    COUNT(ALLOCATED, size);
    return real_aligned_alloc(alignment, size);
}

void *memalign(size_t alignment, size_t size) {
    resolve();
    COUNT(ALLOCATIONS, 1);
    COUNT(ALLOCATED, size);
    return real_memalign(alignment, size);
}

// The Mojo runtime's allocator, as the standard library calls it: an alignment, then the bytes.
void *KGEN_CompilerRT_AlignedAlloc(size_t alignment, size_t size) {
    resolve();
    COUNT(ALLOCATIONS, 1);
    COUNT(ALLOCATED, size);
    return real_kgen_alloc(alignment, size);
}

void KGEN_CompilerRT_AlignedFree(void *pointer) {
    resolve();
    if (pointer)
        COUNT(FREES, 1);
    real_kgen_free(pointer);
}

void *memcpy(void *to, const void *from, size_t size) {
    if (!real_memcpy) {
        char *a = to;
        const char *b = from;
        for (size_t i = 0; i < size; i++)
            a[i] = b[i];
        if (!resolving)
            resolve();
        return to;
    }
    COUNT(COPIES, 1);
    COUNT(COPIED, size);
    return real_memcpy(to, from, size);
}

void *memmove(void *to, const void *from, size_t size) {
    if (!real_memmove) {
        resolve();
        if (!real_memmove) {
            char *a = to;
            const char *b = from;
            if (a < b)
                for (size_t i = 0; i < size; i++)
                    a[i] = b[i];
            else
                for (size_t i = size; i > 0; i--)
                    a[i - 1] = b[i - 1];
            return to;
        }
    }
    COUNT(MOVES, 1);
    COUNT(MOVED, size);
    return real_memmove(to, from, size);
}

void *memset(void *to, int value, size_t size) {
    if (!real_memset) {
        resolve();
        if (!real_memset) {
            char *a = to;
            for (size_t i = 0; i < size; i++)
                a[i] = (char)value;
            return to;
        }
    }
    COUNT(SETS, 1);
    COUNT(SET, size);
    return real_memset(to, value, size);
}
