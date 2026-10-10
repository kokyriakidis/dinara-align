# dinara-align

**Exact pairwise sequence alignment, written in Mojo, for the CPU and the GPU.**

dinara-align finds the optimal alignment of two sequences: the least edit distance, the least
gap-affine cost or the best score. It aligns globally, with free ends, as a seed extension or locally,
and returns the alignment as a CIGAR. Every answer is exact: no heuristic decides the result, and
equally good alignments are broken by one fixed rule. The library is available from Mojo, C, C++
and Python, and as a command-line aligner.

[![CI](https://github.com/kokyriakidis/dinara-align/actions/workflows/ci.yml/badge.svg)](https://github.com/kokyriakidis/dinara-align/actions/workflows/ci.yml)
[![License: MPL 2.0](https://img.shields.io/badge/license-MPL%202.0-blue.svg)](LICENSE)
![Platforms](https://img.shields.io/badge/platforms-Linux%20x86--64%20%7C%20Linux%20AArch64%20%7C%20macOS%20Apple%20silicon-lightgrey)

## Highlights

- **Fast.** In our benchmarks it is the fastest exact aligner on every workload measured, against
  A\*PA2, A\*PA, Edlib, WFA2-lib (WFA and BiWFA), KSW2, parasail, SSW, abPOA and hyalite
  ([results](benchmarks/README.md)).
- **Exact.** Bands, cost caps and memory limits restrict *which* alignments count, never whether the
  best of them is found. The test suite and a fuzzer hold the results to a full dynamic-programming
  matrix, and to WFA2-lib's own regression set, scores and CIGARs alike.
- **Deterministic.** Of equally good alignments, one rule picks the CIGAR: indels placed left, as
  minimap2 places them, or right, WFA2-lib's backtrace byte for byte. The CIGAR does not depend on the
  band, the cap or which algorithm found the cost.
- **Every common cost model and mode.** Unit costs, linear, gap-affine and two-piece gap-affine costs,
  deletions priced apart from insertions, and substitution tables; global, infix, prefix, suffix,
  overlap, any free ends, extension with Z-drop and end bonus, and local alignment.
- **Batches and GPUs.** Batches of short pairs are aligned many at a time in SIMD lanes. Global and
  local scores and alignments under a substitution table, and global unit-cost distances, also run on
  NVIDIA GPUs and Apple silicon.
- **Bounded memory.** Every alignment keeps its traceback within a memory limit (`max_memory`, 80 MB
  by default), splitting a large problem where an optimal path crosses, as BiWFA does.

## How it works

| problem | method |
| :-- | :-- |
| unit costs (edit distance) | A\*PA2's bit-parallel band doubling, pruned by its gap-chaining seed heuristic on long pairs, after a diagonal transition that settles near-identical pairs |
| gap-affine and two-piece costs | a gap-affine wavefront from both ends at once, after WFA and BiWFA |
| local alignment, overlaps, rewarded free ends | an anti-diagonal sweep in 16-bit SIMD lanes while the scores fit, then an extension back from the best end |
| batches of short pairs | inter-sequence SIMD: one pair a lane, over a band of diagonals each pair's own cost proves |
| GPU scores | a band of sixteen diagonals a thread, the rest scored whole, several pairs a warp |

## Performance

Mean time per alignment on one core of an Intel Core i9-7900X at a fixed 3.3 GHz, every tool built
for the machine's instruction set and every tool's answers checked against every other's. A selection
of the rows in [benchmarks/README.md](benchmarks/README.md), which has the full tables, the datasets
and how to reproduce them.

| workload | dinara-align | next fastest exact aligner |
| :-- | --: | :-- |
| ONT reads, mean 0.8 kbp (ont-1k), edit distance | **22 µs** | WFA2-lib, 31 µs |
| ONT reads, mean 3.6 kbp (ont-10k), edit distance | **165 µs** | A\*PA2-simple, 248 µs |
| ONT reads, mean 9.5 kbp (ont-50k), edit distance | **599 µs** | A\*PA2-simple, 940 µs |
| SARS-CoV-2 genomes, 30 kbp, edit distance | **282 µs** | A\*PA2-simple, 728 µs |
| 100 kbp pairs at 1% divergence, edit distance | **1.3 ms** | WFA2-lib, 3.38 ms |
| ONT reads, mean 3.6 kbp (ont-10k), gap-affine (4, 6, 2) | **1.69 ms** | WFA2-lib, 4.84 ms |
| local alignment, 1 kbp read in a 10 kbp window | **1.47 ms** | SSW, 3.19 ms |
| overlap of two 2 kbp reads | **2.22 ms** | parasail, 15.3 ms |
| extension with an end bonus | **157 µs** | KSW2, 1.23 ms |
| 500,000 short reads, gap-affine scores, one core | **230 ms** | WFA2-lib, 1.06 s |

## Installation

### Prebuilt releases

Each [release](https://github.com/kokyriakidis/dinara-align/releases) provides, for Linux x86-64,
Linux AArch64 and macOS on Apple silicon, three archives that run without Mojo or pixi:

| archive | contents |
| :-- | :-- |
| `dinara-align-<version>-<platform>-cli.tar.gz` | the command-line aligner, `cli/dinara-align` |
| `dinara-align-<version>-<platform>-python.tar.gz` | the Python package, `python/dinara_align` |
| `dinara-align-<version>-<platform>-c.tar.gz` | the C library and its header, `c/libdinara` and `c/dinara.h` |

Each archive carries the Mojo runtime libraries it needs beside it. On Linux they need glibc 2.35 or
later (Ubuntu 22.04, Debian 12, RHEL 10). The builds target each platform's baseline CPU.

```bash
tar -xzf dinara-align-v0.1.0-linux-64-cli.tar.gz && ./cli/dinara-align --help
tar -xzf dinara-align-v0.1.0-linux-64-python.tar.gz && PYTHONPATH=$PWD/python python3 -c "import dinara_align"
```

### From source

Building needs [pixi](https://pixi.sh), which installs Mojo and everything else into the project.

```bash
git clone https://github.com/kokyriakidis/dinara-align
cd dinara-align
pixi run test            # build and run the test suite
pixi run build-cli       # build/cli/dinara-align
pixi run build-python    # build/python/dinara_align
pixi run build-wheel     # build/dist/*.whl, for pip install
pixi run build-c         # build/c: libdinara and dinara.h
```

Each build script takes an optional target CPU (for example `pixi run build-cli x86-64-v3`) for a
build tuned to newer machines, which then runs only on CPUs with that instruction set.

## Quick start

### Command line

```bash
$ dinara-align -r ACGTACGTTTGCA -q ACGTCGTTTTGCA
query	reference	strand	cost	score	reference_start	reference_end	query_start	query_end	cigar
query	reference	+	2	-2	0	13	0	13	4=1D2=1I6=

# Reads placed in a reference at gap-affine costs, as SAM
$ dinara-align --costs affine:4,6,2 --mode infix --format sam reference.fa reads.fq
```

Inputs are FASTA or FASTQ files, a tab-separated file of pairs (`--pairs`), or two sequences given
whole (`-r`, `-q`). Output is a table, SAM or PAF; `dinara-align --help` lists every option.

### Python

```python
import dinara_align as da

da.distance("ACGTACGTTTGCA", "ACGTCGTTTTGCA")                          # 2

costs = da.Costs.affine(4, 6, 2)                                        # mismatch 4, gap 6 + 2k
found = da.align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", costs)
found.cost, found.cigar                                                  # (12, '4=3X6=')

placed = da.align("TTTTACGTACGTTTTT", "ACGTACGT", costs, da.Mode.INFIX)
placed.cigar, placed.reference_start, placed.reference_end               # ('8=', 4, 12)

da.align("GGGGACGTACGTGGGG", "CCCCACGTACGTCCCC", costs, da.Mode.local(2)).score   # 16

da.distances(["ACGTACGT", "TTGCA"], ["ACGACGT", "TTGGCA"], threads=4)  # [1, 1]
```

### C and C++

```c
#include "dinara.h"

int64_t edits = dinara_distance("ACGTACGTTTGCA", 13, "ACGTCGTTTTGCA", 13, NULL, NULL, NULL);  /* 2 */

dinara_costs affine = {4, 6, 2, -1, 0, 0, 0, -1, 0};
dinara_alignment found;
if (dinara_align("ACGTACGTTTGCA", 13, "ACGTCGTTTTGCA", 13, &affine, NULL, NULL, &found) == 0) {
    /* found.cost == 12, found.cigar == "4=3X6=" */
    dinara_free(found.cigar);
}
```

```cpp
#include "dinara.h"

auto found = dinara::align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", dinara::Costs::affine(4, 6, 2));
// found.cost == 12, found.cigar == "4=3X6="
```

Link with `-L build/c -ldinara` and put `build/c` on the program's library search path. Every
function returns a negative `DINARA_` code on failure; the header documents each one.

### Mojo

```mojo
from dinara_align import Costs, Mode, align, distance

def main() raises:
    print(distance("ACGTACGTTTGCA", "ACGTCGTTTTGCA"))  # 2
    var found = align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", Costs.affine(4, 6, 2))
    print(found.cost, found.cigar)  # 12 4=3X6=
    var core = align("GGGGACGTACGTGGGG", "CCCCACGTACGTCCCC", Costs.affine(4, 6, 2), Mode.local(2))
    print(core.score)  # 16
```

Build a program against the source tree with `mojo build -I path/to/dinara-align program.mojo`, or
inside this repository with `pixi run mojo build -I . program.mojo`.

## Alignment modes

A mode says which ends of the two sequences an alignment must reach. There are three kinds, free
ends, extension and local, and the named free-end modes are presets of `ends_free`.

| mode | reference | query | also known as |
| :-- | :-- | :-- | :-- |
| `GLOBAL` | whole | whole | end to end, Needleman-Wunsch, Edlib's NW |
| `INFIX` | any part | whole | semi-global, glocal, Edlib's HW |
| `PREFIX`, `SUFFIX` | a prefix, a suffix | whole | Edlib's SHW |
| `overlap(m)` | a prefix or suffix | a suffix or prefix | parasail's `sg` |
| `ends_free(...)` | as many letters free at each end as asked | likewise | WFA2-lib's ends-free |
| `extension(m, anchor)` | from one end, as far as it pays | from the same end | KSW2's extension, with Z-drop and end bonus |
| `local(m)` | any part | any part | Smith-Waterman |

Free ends minimize the cost alone, as Edlib and WFA2-lib count it; `mode.with_match_score(m)` makes a
match earn `m` and maximizes the score instead, as parasail and mappers score a read in a window.

## Cost models

| costs | a gap of `k` letters costs |
| :-- | :-- |
| `Costs.edit()` | `k`, every edit 1: the edit distance (the default) |
| `Costs.linear(x, g)` | `g k`, a mismatch `x` |
| `Costs.affine(x, o, e)` | `o + e k`, as WFA2-lib and minimap2 count it |
| `Costs.two_piece(x, o, e, o2, e2)` | `min(o + e k, o2 + e2 k)`, minimap2's `-O o,o2 -E e,e2` |
| `costs.with_deletions(o, e, ...)` | deletions priced apart from insertions, bwa's `-O del,ins` |
| `Scoring.tabulated(...)`, `Scoring.dna()` | any alphabet's substitution table with affine gap scores |

Within a band of diagonals (`band=Band.around(w)`, KSW2's `w`) and under a cost cap (`max_cost=`)
the result is still the exact optimum of the alignments allowed.

## Documentation

- [docs/API.md](docs/API.md): the reference for everything the Mojo package exports.
- [c/dinara.h](c/dinara.h): the C and C++ API, every function and error code documented.
- The Python package's docstrings, and `dinara-align --help` for the command line.
- [benchmarks/README.md](benchmarks/README.md): the benchmarks, their setup and full results.
- [CHANGELOG.md](CHANGELOG.md): what each release changes.

## Testing

```bash
pixi run test          # the Mojo test suite
pixi run fuzz          # random pairs, costs, modes, bands and caps, held to a full matrix
pixi run test-wfa      # WFA2-lib's regression set, scores and CIGARs
pixi run test-c        # the C and C++ API
pixi run test-python   # the Python package
pixi run test-cli      # the command line
```

## Acknowledgements

The unit-cost edit distance is ported from `pa-bitpacking` in
[A\*PA](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) by Ragnar Groot Koerkamp and
Pesho Ivanov, and follows A\*PA2's band doubling and seed heuristic. The gap-affine search follows WFA
and BiWFA by Santiago Marco-Sola and colleagues. If you use dinara-align, please also cite the work
it builds on:

- R. Groot Koerkamp. *A\*PA2: up to 19× faster exact global alignment.* WABI 2024.
- R. Groot Koerkamp and P. Ivanov. *Exact global alignment using A\* with chaining seed heuristic and
  match pruning.* Bioinformatics, 2024.
- S. Marco-Sola, J. C. Moure, M. Moreto and A. Espinosa. *Fast gap-affine pairwise alignment using the
  wavefront algorithm.* Bioinformatics, 2021.
- S. Marco-Sola, J. M. Eizenga, A. Guarracino, B. Paten, E. Garrison and M. Moreto. *Optimal gap-affine
  alignment in O(s) space.* Bioinformatics, 2023.

## License

dinara-align is licensed under the [Mozilla Public License 2.0](LICENSE). [NOTICE](NOTICE) credits
the work it builds on.
