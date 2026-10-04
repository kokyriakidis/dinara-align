# Benchmarks

Times dinara-align against other exact global DNA aligners on the same inputs, and refuses to report a time for any workload where the tools disagree on the answer.

```bash
pixi run bench                    # 1 kbp and 10 kbp workloads, about a minute
pixi run bench --full             # adds the 100 kbp pairs, about four and a half minutes
pixi run bench --install-rust     # also installs the nightly Rust A*PA needs, under .cache/
```

The table prints to the terminal and lands in `.cache/results/results.md`, with every raw row in `results.tsv`.

## What Runs

Every rival is cloned at a pinned commit into `.cache/` and built there; nothing of theirs is vendored into this repository.

| Tool                                                                            | Commit    | License    | What it is run on                                                          |
| :------------------------------------------------------------------------------ | :-------- | :--------- | :------------------------------------------------------------------------- |
| dinara-align                                                                    | this tree | Apache-2.0 | the read batches and the affine pairs, on the CPU and on the GPU where one answers |
| dinara-align (bit-parallel, …)                                                  | this tree | MPL-2.0    | the edit-distance workloads: `edit_distance` and `edit_alignment`, ported from A\*PA2-simple, on one thread and on all of them |
| [hyalite](https://github.com/Psy-Fer/hyalite)                                   | `0189bcb` | MIT        | the read batches and the affine pairs, global mode (`Mode::Nw`), one CPU thread with its NEON or AVX2 kernels |
| [A\*PA, A\*PA2](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) | `bf2e14e` | MPL-2.0    | the edit-distance workloads, because it computes edit distance and nothing else |

A\*PA runs four of its aligners: A\*PA2-full (its default), A\*PA2-simple, the original A\*PA, and `astarpa2_nw`, a bit-parallel full-matrix Needleman-Wunsch that is the closest like-for-like to a full sweep.
The three A\*PA2 aligners run twice, with traceback and without, so each has a score row and an alignment row; the original A\*PA has no cost-only entry point and reports its alignment alone.
`a*pa2-nw` against `dinara-align (bit-parallel, 1 thread)` is the same algorithm in Rust and in Mojo; the all-threads column is the same kernel tiled across cores, which A\*PA does not do.
A\*PA2-simple's score-only mode is slower than its alignment mode on long pairs, so read its best row, not its score row, as its speed.

The edit-distance pairs go to the edit-distance aligners alone.
dinara-align's affine-gap kernels and hyalite would answer them with a full quadratic sweep, over a minute per 100 kbp pair, which compares nothing the affine pairs do not and made a full run take forty minutes.

Many rows run in microseconds, where a single cold call measures the allocator and the scheduler as much as the aligner.
Every runner therefore times each measurement the same way, warm and in-process: a 200 ms spin first so the process sits on a fast core, then one call to size the batches, then twenty batches of about 10 ms each, reporting the fastest batch's average; a call of 0.1 s or more is timed once, where noise is small beside it.
So one run of each runner suffices, and `run.py` makes one unless `--repeat` asks for more.

## Workloads

All global alignment over `ACGT`, generated from a fixed seed by `run.py`.
The second sequence of every pair is a mutation of the first, by substitutions, deletions and insertions in equal thirds.

| Workload                    | Pairs                                                   | Scoring                                             |
| :-------------------------- | :------------------------------------------------------ | :-------------------------------------------------- |
| `reads-150bp`               | 10,000 pairs of 150 bp at 2% divergence, as one batch   | `Scoring.dna()`: minimap2's match 2, mismatch −4, `-O4 -E2` |
| `reads-1kbp`                | 1,000 pairs of 1 kbp at 10% divergence, as one batch    | the same                                            |
| `affine-{1k,10k,100k}`      | one pair each at 5% divergence; 100k only with `--full` | the same                                            |
| `edit-{1k,10k,100k}-{1,5,15}%` | one pair each at that divergence                     | edit distance: unit costs                           |

The affine pairs use a single divergence because a full sweep takes the same time however similar the sequences are.
Edit distance sweeps three, because A\*PA's time depends on it.
dinara-align writes its default scoring to the data directory and hyalite reads it from there, so both score with the same numbers; both charge a gap's first base the opening penalty and each further base the extension.

## Results

Measured with `pixi run bench --full`, which took four and a half minutes, on an Apple M2 (4 performance and 4 efficiency cores, 24 GB) under macOS 27.0.1, with Mojo 1.1.0 and Xcode 27.0.
Hyalite built with Rust 1.92, A\*PA with the nightly its repository pins.
A dash marks a task the tool does not offer, or a workload it is not run on.

| workload | task | dinara-align (cpu) | dinara-align (gpu) | dinara-align (bit-parallel, 1 thread) | dinara-align (bit-parallel, 8 threads) | hyalite | a*pa2-full | a*pa2-simple | a*pa2-nw | a*pa | agree |
| :-- | :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| reads-150bp | score | 642 ms | 15.6 ms | — | — | 279 ms | — | — | — | — | ✓ |
| reads-150bp | alignment | 1.05 s | 139 ms | — | — | 930 ms | — | — | — | — | ✓ |
| reads-1kbp | score | 2.98 s | 39.2 ms | — | — | 1.23 s | — | — | — | — | ✓ |
| reads-1kbp | alignment | 4.55 s | 211 ms | — | — | 5.69 s | — | — | — | — | ✓ |
| affine-1k | score | 2.93 ms | 5.15 ms | — | — | 1.24 ms | — | — | — | — | ✓ |
| affine-1k | alignment | 4.38 ms | 9.16 ms | — | — | 4.3 ms | — | — | — | — | ✓ |
| affine-10k | score | 305 ms | 18.7 ms | — | — | 209 ms | — | — | — | — | ✓ |
| affine-10k | alignment | 729 ms | 33.7 ms | — | — | 788 ms | — | — | — | — | ✓ |
| affine-100k | score | 31.6 s | 362 ms | — | — | 24.9 s | — | — | — | — | ✓ |
| affine-100k | alignment | 76.6 s | 723 ms | — | — | 81 s | — | — | — | — | ✓ |
| edit-1k-1% | score | — | — | 2 µs | 2 µs | — | 37 µs | 20 µs | 22 µs | — | ✓ |
| edit-1k-1% | alignment | — | — | 3 µs | 3 µs | — | 39 µs | 18 µs | 57 µs | 149 µs | ✓ |
| edit-1k-5% | score | — | — | 10 µs | 11 µs | — | 37 µs | 20 µs | 22 µs | — | ✓ |
| edit-1k-5% | alignment | — | — | 11 µs | 11 µs | — | 40 µs | 18 µs | 58 µs | 144 µs | ✓ |
| edit-1k-15% | score | — | — | 10 µs | 10 µs | — | 31 µs | 20 µs | 22 µs | — | ✓ |
| edit-1k-15% | alignment | — | — | 21 µs | 21 µs | — | 48 µs | 31 µs | 60 µs | 508 µs | ✓ |
| edit-10k-1% | score | — | — | 76 µs | 73 µs | — | 489 µs | 931 µs | 1.69 ms | — | ✓ |
| edit-10k-1% | alignment | — | — | 54 µs | 55 µs | — | 504 µs | 153 µs | 2.97 ms | 1.61 ms | ✓ |
| edit-10k-5% | score | — | — | 98 µs | 92 µs | — | 502 µs | 1.63 ms | 1.67 ms | — | ✓ |
| edit-10k-5% | alignment | — | — | 156 µs | 145 µs | — | 537 µs | 267 µs | 2.97 ms | 1.49 ms | ✓ |
| edit-10k-15% | score | — | — | 204 µs | 164 µs | — | 867 µs | 2.79 ms | 1.72 ms | — | ✓ |
| edit-10k-15% | alignment | — | — | 336 µs | 244 µs | — | 1.09 ms | 770 µs | 3.25 ms | 39.1 ms | ✓ |
| edit-100k-1% | score | — | — | 1.67 ms | 1.2 ms | — | 4.28 ms | 189 ms | 162 ms | — | ✓ |
| edit-100k-1% | alignment | — | — | 2.02 ms | 1.5 ms | — | 4.43 ms | 3.32 ms | 821 ms | 19.6 ms | ✓ |
| edit-100k-5% | score | — | — | 4.53 ms | 3.55 ms | — | 4.6 ms | 274 ms | 162 ms | — | ✓ |
| edit-100k-5% | alignment | — | — | 5.3 ms | 3.95 ms | — | 5 ms | 16.1 ms | 752 ms | 17.8 ms | ✓ |
| edit-100k-15% | score | — | — | 11.4 ms | 8.63 ms | — | 21.7 ms | 242 ms | 163 ms | — | ✓ |
| edit-100k-15% | alignment | — | — | 13 ms | 9.63 ms | — | 23.8 ms | 25.3 ms | 485 ms | 7.66 s | ✓ |

## Reading the Numbers

- **Every answer agreed**, so each row compares tools answering the same question. The answers are the optimal scores and distances: a sum and a position-weighted sum per workload, so a reordered batch cannot pass.
- **On edit distance, the bit-parallel path beats A\*PA2-simple on every workload, on one thread**, distance and alignment alike: 3 against 18 µs aligning 1 kbp at 1%, 21 against 31 µs at 1 kbp and 15%, 54 against 153 µs at 10 kbp and 1%, 13 against 25 ms at 100 kbp and 15%. Its score rows beat every A\*PA2 aligner's score row by more, since A\*PA2's cost-only mode is its slower one.
- **A\*PA2-full is still faster aligning 100 kbp at 5%**: 5.0 ms, against 5.3 ms on one thread and 3.95 ms on all of them. It is not a faster kernel: per cell it is about half as fast. Its gap-chaining seed heuristic, built from k-mer matches between the sequences in about 1.5 ms, prunes the band to about a third of the cells the gap heuristic leaves, which is the heuristic dinara-align ports from A\*PA2-simple.
- **Near-identical pairs finish before any band.** dinara-align first runs diagonal transition, as WFA does, while it stays cheaper than a band would be; at 1% divergence it finishes there, alignment included, which is where A\*PA2's own write-up places its weak spot.
- **The other rivals run on one CPU thread.** dinara-align's affine CPU batch also runs one pair at a time today, and its linear-space traceback does not yet speed up with more threads, so those CPU columns compare single-threaded code throughout. The bit-parallel all-threads column sweeps a long pair from both ends at once, one direction per thread.
- **dinara-align's GPU times include the overheads a caller pays**: opening the device context on every call, copying sequences over, and copying results back. That is why a single short pair is slower on the GPU than on the CPU.
- **hyalite's traceback budget is 1 GiB.** It keeps the whole matrix when it fits and switches to a checkpointed sweep above that, as the 100 kbp pairs do. dinara-align switches to its linear-space recursion above six million cells on the host and one million on the device.
- **A\*PA answers only edit distance**, and its time depends on how similar the sequences are, as does the bit-parallel path's; the affine columns do not.
- **One machine, one fixed pair per workload.** On fresh pairs from two other seeds, every pair from 1 to 6 kbp at 2 to 30% divergence beat both A\*PA2 aligners on distance and alignment, some by only a few percent; on 100 kbp pairs at 2 to 10%, A\*PA2-full stayed ahead on one thread. Rerun on your own hardware before quoting any of this.
