# Benchmarks

Times dinara-align against other exact global DNA aligners on the same inputs, and refuses to report a time for any workload where the tools disagree on the answer.

```bash
pixi run bench                    # 1 kbp and 10 kbp workloads, best of three, about two minutes
pixi run bench --full             # adds the 100 kbp pairs, about forty minutes
pixi run bench --install-rust     # also installs the nightly Rust A*PA needs, under .cache/
```

The table prints to the terminal and lands in `.cache/results/results.md`, with every raw row in `results.tsv`.

## What Runs

Every rival is cloned at a pinned commit into `.cache/` and built there; nothing of theirs is vendored into this repository.

| Tool                                                                            | Commit    | License    | What it is run on                                                          |
| :------------------------------------------------------------------------------ | :-------- | :--------- | :------------------------------------------------------------------------- |
| dinara-align                                                                    | this tree | Apache-2.0 | every workload, on the CPU and on the GPU where one answers                |
| dinara-align (bit-parallel, …)                                                  | this tree | MPL-2.0    | the edit-distance workloads: `edit_distance` and `edit_alignment`, ported from A\*PA2-simple, on one thread and on all of them |
| [hyalite](https://github.com/Psy-Fer/hyalite)                                   | `0189bcb` | MIT        | every workload, global mode (`Mode::Nw`), one CPU thread with its NEON or AVX2 kernels |
| [A\*PA, A\*PA2](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) | `bf2e14e` | MPL-2.0    | the edit-distance workloads, because it computes edit distance and nothing else |

A\*PA runs four of its aligners: A\*PA2-full (its default), A\*PA2-simple, the original A\*PA, and `astarpa2_nw`, a bit-parallel full-matrix Needleman-Wunsch that is the closest like-for-like to a full sweep.
The three A\*PA2 aligners run twice, with traceback and without, so each has a score row and an alignment row; the original A\*PA has no cost-only entry point and reports its alignment alone.
`a*pa2-nw` against `dinara-align (bit-parallel, 1 thread)` is the same algorithm in Rust and in Mojo; the all-threads column is the same kernel tiled across cores, which A\*PA does not do.
A\*PA2-simple's score-only mode is slower than its alignment mode on long pairs, so read its best row, not its score row, as its speed.

The edit-distance rows of A\*PA and of dinara-align's bit-parallel columns run in microseconds, where a single cold call measures the allocator and the scheduler as much as the aligner.
Both runners therefore time those rows the same way, warm and in-process: a 200 ms spin first so the process sits on a fast core, then one call to size the batches, then twenty batches of about 10 ms each, reporting the fastest batch's average; a call of 0.1 s or more is timed once.
Every other row is one call per run, the fastest of three runs.

## Workloads

All global alignment over `ACGT`, generated from a fixed seed by `run.py`.
The second sequence of every pair is a mutation of the first, by substitutions, deletions and insertions in equal thirds.

| Workload                    | Pairs                                                   | Scoring                                             |
| :-------------------------- | :------------------------------------------------------ | :-------------------------------------------------- |
| `reads-150bp`               | 10,000 pairs of 150 bp at 2% divergence, as one batch   | `Scoring.dna()`: minimap2's match 2, mismatch −4, `-O4 -E2` |
| `reads-1kbp`                | 1,000 pairs of 1 kbp at 10% divergence, as one batch    | the same                                            |
| `affine-{1k,10k,100k}`      | one pair each at 5% divergence; 100k only with `--full` | the same                                            |
| `edit-{1k,10k,100k}-{1,5,15}%` | one pair each at that divergence                     | edit distance: match 0, mismatch and gaps −1        |

The affine pairs use a single divergence because a full sweep takes the same time however similar the sequences are.
Edit distance sweeps three, because A\*PA's time depends on it.
dinara-align writes its default scoring to the data directory and hyalite reads it from there, so both score with the same numbers; both charge a gap's first base the opening penalty and each further base the extension.

## Results

Measured with `pixi run bench --full` at commit `641819a`, on an Apple M2 (4 performance and 4 efficiency cores, 24 GB) under macOS 27.0.1, with Mojo 1.1.0 and Xcode 27.0.
Hyalite built with Rust 1.92, A\*PA with the nightly its repository pins.
A dash marks a task the tool does not offer.

| workload | task | dinara-align (cpu) | dinara-align (gpu) | dinara-align (bit-parallel, 1 thread) | dinara-align (bit-parallel, 8 threads) | hyalite | a*pa2-full | a*pa2-simple | a*pa2-nw | a*pa | agree |
| :-- | :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| reads-150bp | score | 644 ms | 31.8 ms | — | — | 277 ms | — | — | — | — | ✓ |
| reads-150bp | alignment | 1.04 s | 114 ms | — | — | 909 ms | — | — | — | — | ✓ |
| reads-1kbp | score | 2.99 s | 39.9 ms | — | — | 1.24 s | — | — | — | — | ✓ |
| reads-1kbp | alignment | 4.51 s | 184 ms | — | — | 4.81 s | — | — | — | — | ✓ |
| affine-1k | score | 2.94 ms | 5.47 ms | — | — | 1.24 ms | — | — | — | — | ✓ |
| affine-1k | alignment | 4.39 ms | 8.89 ms | — | — | 3.89 ms | — | — | — | — | ✓ |
| affine-10k | score | 305 ms | 18.7 ms | — | — | 205 ms | — | — | — | — | ✓ |
| affine-10k | alignment | 725 ms | 35.2 ms | — | — | 767 ms | — | — | — | — | ✓ |
| affine-100k | score | 30.8 s | 357 ms | — | — | 19.7 s | — | — | — | — | ✓ |
| affine-100k | alignment | 73.7 s | 717 ms | — | — | 77.7 s | — | — | — | — | ✓ |
| edit-1k-1% | score | 2.95 ms | 5.64 ms | 8 µs | 8 µs | 1.1 ms | 37 µs | 20 µs | 22 µs | — | ✓ |
| edit-1k-1% | alignment | 4.55 ms | 9.49 ms | 15 µs | 15 µs | 5.05 ms | 38 µs | 18 µs | 57 µs | 149 µs | ✓ |
| edit-1k-5% | score | 2.96 ms | 5.37 ms | 9 µs | 9 µs | 1.11 ms | 37 µs | 20 µs | 22 µs | — | ✓ |
| edit-1k-5% | alignment | 4.5 ms | 9.31 ms | 17 µs | 17 µs | 4.98 ms | 40 µs | 18 µs | 58 µs | 139 µs | ✓ |
| edit-1k-15% | score | 3 ms | 5.38 ms | 16 µs | 16 µs | 1.13 ms | 31 µs | 20 µs | 22 µs | — | ✓ |
| edit-1k-15% | alignment | 4.59 ms | 11.6 ms | 29 µs | 29 µs | 5.69 ms | 48 µs | 31 µs | 59 µs | 505 µs | ✓ |
| edit-10k-1% | score | 307 ms | 18.4 ms | 78 µs | 77 µs | 112 ms | 485 µs | 928 µs | 1.69 ms | — | ✓ |
| edit-10k-1% | alignment | 728 ms | 34.6 ms | 140 µs | 129 µs | 879 ms | 496 µs | 152 µs | 2.87 ms | 1.6 ms | ✓ |
| edit-10k-5% | score | 308 ms | 17.7 ms | 158 µs | 152 µs | 112 ms | 497 µs | 1.62 ms | 1.67 ms | — | ✓ |
| edit-10k-5% | alignment | 728 ms | 33.9 ms | 242 µs | 205 µs | 869 ms | 535 µs | 264 µs | 2.9 ms | 1.49 ms | ✓ |
| edit-10k-15% | score | 308 ms | 18.2 ms | 232 µs | 206 µs | 112 ms | 863 µs | 2.78 ms | 1.69 ms | — | ✓ |
| edit-10k-15% | alignment | 732 ms | 34 ms | 377 µs | 304 µs | 826 ms | 1.08 ms | 759 µs | 3.13 ms | 39 ms | ✓ |
| edit-100k-1% | score | 30.8 s | 356 ms | 2.38 ms | 1.53 ms | 17.7 s | 4.26 ms | 182 ms | 161 ms | — | ✓ |
| edit-100k-1% | alignment | 73.3 s | 714 ms | 3.01 ms | 1.97 ms | 77.8 s | 4.39 ms | 3.32 ms | 404 ms | 19.5 ms | ✓ |
| edit-100k-5% | score | 30.7 s | 360 ms | 6.36 ms | 4.48 ms | 17.6 s | 4.58 ms | 272 ms | 161 ms | — | ✓ |
| edit-100k-5% | alignment | 73.5 s | 723 ms | 7.28 ms | 5.08 ms | 77.9 s | 4.96 ms | 16.1 ms | 335 ms | 17.7 ms | ✓ |
| edit-100k-15% | score | 30.7 s | 357 ms | 16.3 ms | 10.3 ms | 17.6 s | 21.7 ms | 241 ms | 161 ms | — | ✓ |
| edit-100k-15% | alignment | 73.6 s | 718 ms | 18.1 ms | 11.4 ms | 77.6 s | 23.6 ms | 25 ms | 333 ms | 6.44 s | ✓ |

## Reading the Numbers

- **Every answer agreed**, so each row compares tools answering the same question. The answers are the optimal scores and distances: a sum and a position-weighted sum per workload, so a reordered batch cannot pass.
- **On edit distance, the bit-parallel path beats A\*PA2-simple on every workload, on one thread**, distance and alignment alike: 15 against 18 µs aligning 1 kbp at 1%, 140 against 152 µs at 10 kbp and 1%, 18.1 against 25 ms at 100 kbp and 15%. Its score rows beat every A\*PA2 aligner's score row by more, since A\*PA2's cost-only mode is its slower one.
- **A\*PA2-full is still faster at 100 kbp and 5%**: 4.96 ms aligning, against 7.28 ms on one thread and 5.08 ms on all of them. Its gap-chaining seed heuristic prunes far more of a long, moderately divergent pair than the gap heuristic dinara-align ports from A\*PA2-simple.
- **The other rivals run on one CPU thread.** dinara-align's affine CPU batch also runs one pair at a time today, and its linear-space traceback does not yet speed up with more threads, so those CPU columns compare single-threaded code throughout. The bit-parallel all-threads column sweeps a long pair from both ends at once, one direction per thread.
- **dinara-align's GPU times include the overheads a caller pays**: opening the device context on every call, copying sequences over, and copying results back. That is why a single short pair is slower on the GPU than on the CPU.
- **hyalite's traceback budget is 1 GiB.** It keeps the whole matrix when it fits and switches to a checkpointed sweep above that, as the 100 kbp pairs do. dinara-align switches to its linear-space recursion above six million cells on the host and one million on the device.
- **A\*PA answers only edit distance**, and its time depends on how similar the sequences are, as does the bit-parallel path's; the affine columns do not. On near-identical DNA both skip almost the whole matrix.
- **One machine.** Rerun on your own hardware before quoting any of this.
