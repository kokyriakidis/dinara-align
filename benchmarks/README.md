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
| dinara-align (bit-parallel, …)                                                  | this tree | MPL-2.0    | the edit-distance workloads, score only: `edit_distance`, ported from A\*PA, on one thread and on all of them |
| [hyalite](https://github.com/Psy-Fer/hyalite)                                   | `0189bcb` | MIT        | every workload, global mode (`Mode::Nw`), one CPU thread with its NEON or AVX2 kernels |
| [A\*PA, A\*PA2](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) | `bf2e14e` | MPL-2.0    | the edit-distance workloads, because it computes edit distance and nothing else |

A\*PA runs four of its aligners: A\*PA2-full (its default), A\*PA2-simple, the original A\*PA, and `astarpa2_nw`, a bit-parallel full-matrix Needleman-Wunsch that is the closest like-for-like to a full sweep.
The three A\*PA2 aligners run twice, with traceback and without, so each has a score row and an alignment row; the original A\*PA has no cost-only entry point and reports its alignment alone.
`a*pa2-nw` against `dinara-align (bit-parallel, 1 thread)` is the same algorithm in Rust and in Mojo; the all-threads column is the same kernel tiled across cores, which A\*PA does not do.
A\*PA2-simple's score-only mode is slower than its alignment mode on long pairs, so read its best row, not its score row, as its speed.

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

Measured with `pixi run bench --full`, best of three, on an Apple M2 (4 performance and 4 efficiency cores, 24 GB) under macOS 27.0.1, with Mojo 1.1.0 and Xcode 27.0.
Hyalite built with Rust 1.92, A\*PA with the nightly its repository pins.
A dash marks a task the tool does not offer.

| workload | task | dinara-align (cpu) | dinara-align (gpu) | hyalite | a*pa2-full | a*pa2-simple | a*pa | a*pa2-nw | agree |
| :-- | :-- | --: | --: | --: | --: | --: | --: | --: | :-: |
| reads-150bp | score | 641 ms | 26.7 ms | 278 ms | — | — | — | — | ✓ |
| reads-150bp | alignment | 1.05 s | 99.1 ms | 905 ms | — | — | — | — | ✓ |
| reads-1kbp | score | 2.98 s | 39.8 ms | 1.23 s | — | — | — | — | ✓ |
| reads-1kbp | alignment | 4.5 s | 134 ms | 5.58 s | — | — | — | — | ✓ |
| affine-1k | score | 2.93 ms | 5.43 ms | 1.24 ms | — | — | — | — | ✓ |
| affine-1k | alignment | 4.39 ms | 8.99 ms | 5.28 ms | — | — | — | — | ✓ |
| affine-10k | score | 304 ms | 18.1 ms | 205 ms | — | — | — | — | ✓ |
| affine-10k | alignment | 730 ms | 33.9 ms | 942 ms | — | — | — | — | ✓ |
| affine-100k | score | 30.7 s | 356 ms | 19.7 s | — | — | — | — | ✓ |
| affine-100k | alignment | 73.3 s | 717 ms | 77.6 s | — | — | — | — | ✓ |
| edit-1k-1% | score | 2.95 ms | 5.39 ms | 1.1 ms | — | — | — | — | ✓ |
| edit-1k-1% | alignment | 4.49 ms | 9.12 ms | 5.01 ms | 90 µs | 29 µs | 246 µs | 93 µs | ✓ |
| edit-1k-5% | score | 2.97 ms | 5.18 ms | 1.12 ms | — | — | — | — | ✓ |
| edit-1k-5% | alignment | 4.48 ms | 8.92 ms | 4.99 ms | 67 µs | 31 µs | 211 µs | 73 µs | ✓ |
| edit-1k-15% | score | 3 ms | 5.15 ms | 1.13 ms | — | — | — | — | ✓ |
| edit-1k-15% | alignment | 4.54 ms | 10.5 ms | 4.92 ms | 82 µs | 47 µs | 606 µs | 77 µs | ✓ |
| edit-10k-1% | score | 305 ms | 18.4 ms | 112 ms | — | — | — | — | ✓ |
| edit-10k-1% | alignment | 729 ms | 34.5 ms | 826 ms | 593 µs | 166 µs | 1.81 ms | 5.03 ms | ✓ |
| edit-10k-5% | score | 306 ms | 17.8 ms | 112 ms | — | — | — | — | ✓ |
| edit-10k-5% | alignment | 725 ms | 34.5 ms | 819 ms | 627 µs | 294 µs | 1.68 ms | 4.02 ms | ✓ |
| edit-10k-15% | score | 305 ms | 18.4 ms | 112 ms | — | — | — | — | ✓ |
| edit-10k-15% | alignment | 732 ms | 34.2 ms | 951 ms | 1.18 ms | 808 µs | 40.2 ms | 3.23 ms | ✓ |
| edit-100k-1% | score | 31.1 s | 361 ms | 17.7 s | — | — | — | — | ✓ |
| edit-100k-1% | alignment | 73.3 s | 717 ms | 77.9 s | 4.7 ms | 3.36 ms | 22.1 ms | 404 ms | ✓ |
| edit-100k-5% | score | 31.1 s | 361 ms | 17.7 s | — | — | — | — | ✓ |
| edit-100k-5% | alignment | 73.3 s | 716 ms | 77.8 s | 6.23 ms | 16.3 ms | 20.6 ms | 339 ms | ✓ |
| edit-100k-15% | score | 30.9 s | 361 ms | 17.6 s | — | — | — | — | ✓ |
| edit-100k-15% | alignment | 73.3 s | 717 ms | 77.8 s | 24.4 ms | 25.2 ms | 6.41 s | 336 ms | ✓ |

## Reading the Numbers

- **Every answer agreed**, so each row compares tools answering the same question. The answers are the optimal scores and distances: a sum and a position-weighted sum per workload, so a reordered batch cannot pass.
- **The rivals run on one CPU thread.** dinara-align's CPU batch also runs one pair at a time today, and its linear-space traceback does not yet speed up with more threads, so the CPU columns compare single-threaded code throughout.
- **dinara-align's GPU times include the overheads a caller pays**: opening the device context on every call, copying sequences over, and copying results back. That is why a single short pair is slower on the GPU than on the CPU.
- **hyalite's traceback budget is 1 GiB.** It keeps the whole matrix when it fits and switches to a checkpointed sweep above that, as the 100 kbp pairs do. dinara-align switches to its linear-space recursion above six million cells on the host and one million on the device.
- **A\*PA answers only edit distance**, and its time depends on how similar the sequences are, where every other column here does not. On near-identical DNA it skips almost the whole matrix.
- **One machine, best of three.** Rerun on your own hardware before quoting any of this.
