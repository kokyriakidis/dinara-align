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
| [Edlib](https://github.com/Martinsos/edlib), [WFA2-lib](https://github.com/smarco/WFA2-lib) through [pa-bench](https://github.com/pairwise-alignment/pa-bench) | `af7a50d` | MIT | the edit-distance workloads: the other exact aligners A\*PA2's evaluation compares against, with its parameters |

A\*PA runs four of its aligners: A\*PA2-full (its default), A\*PA2-simple, the original A\*PA, and `astarpa2_nw`, a bit-parallel full-matrix Needleman-Wunsch that is the closest like-for-like to a full sweep.
The three A\*PA2 aligners run twice, with traceback and without, so each has a score row and an alignment row; the original A\*PA has no cost-only entry point and reports its alignment alone.
`a*pa2-nw` against `dinara-align (bit-parallel, 1 thread)` is the same algorithm in Rust and in Mojo; the all-threads column is the same kernel tiled across cores, which A\*PA does not do.
A\*PA2-simple's score-only mode is slower than its alignment mode on long pairs, so read its best row, not its score row, as its speed.

pa-bench's `pa-wrapper` builds the rest of A\*PA2's comparison, each as its evaluation configured it: Edlib (band doubling over Myers' bit-vectors), BiWFA (WFA2-lib's bidirectional wavefront, its lowest-memory mode), and WFA (WFA2-lib keeping every wavefront, its fastest alignment).
Every column computes the optimal alignment, which is the contest here: WFA-adaptive, which drops diagonals lagging behind and may miss the optimum, and Block Aligner, which bounds its band, are left out, as A\*PA2's evaluation counts them approximate.
KSW2's kernels are SSE only and do not build on ARM, and TripleAccel's quadratic time would add minutes per 100 kbp pair.

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

The affine pairs use a single divergence, 5%, where a full sweep takes the same time however similar the sequences are and dinara-align's global score instead runs a wavefront whose time grows with the score (see below).
Edit distance sweeps three, because A\*PA's time depends on it.
dinara-align writes its default scoring to the data directory and hyalite reads it from there, so both score with the same numbers; both charge a gap's first base the opening penalty and each further base the extension.

## Results

Measured with `pixi run bench --full`, which took about five minutes, on an Apple M2 (4 performance and 4 efficiency cores, 24 GB) under macOS 27.0.1, with Mojo 1.1.0 and Xcode 27.0.
Hyalite built with Rust 1.92, A\*PA with the nightly its repository pins.
A dash marks a task the tool does not offer, or a workload it is not run on.

| workload | task | dinara-align (cpu) | dinara-align (gpu) | dinara-align (bit-parallel, 1 thread) | dinara-align (bit-parallel, 8 threads) | hyalite | a*pa2-full | a*pa2-simple | a*pa2-nw | a*pa | edlib | biwfa | wfa | agree |
| :-- | :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| reads-150bp | score | 13 ms | 15.6 ms | — | — | 280 ms | — | — | — | — | — | — | — | ✓ |
| reads-150bp | alignment | 283 ms | 137 ms | — | — | 919 ms | — | — | — | — | — | — | — | ✓ |
| reads-1kbp | score | 59.4 ms | 39.3 ms | — | — | 1.25 s | — | — | — | — | — | — | — | ✓ |
| reads-1kbp | alignment | 1.24 s | 266 ms | — | — | 5.05 s | — | — | — | — | — | — | — | ✓ |
| affine-1k | score | 109 µs | 5.49 ms | — | — | 1.24 ms | — | — | — | — | — | — | — | ✓ |
| affine-1k | alignment | 4.46 ms | 9.41 ms | — | — | 4.26 ms | — | — | — | — | — | — | — | ✓ |
| affine-10k | score | 5.04 ms | 18.4 ms | — | — | 215 ms | — | — | — | — | — | — | — | ✓ |
| affine-10k | alignment | 737 ms | 33.7 ms | — | — | 870 ms | — | — | — | — | — | — | — | ✓ |
| affine-100k | score | 438 ms | 366 ms | — | — | 20.3 s | — | — | — | — | — | — | — | ✓ |
| affine-100k | alignment | 77.3 s | 739 ms | — | — | 79.3 s | — | — | — | — | — | — | — | ✓ |
| edit-1k-1% | score | — | — | 1 µs | 1 µs | — | 37 µs | 20 µs | 22 µs | — | 13 µs | 1 µs | 1 µs | ✓ |
| edit-1k-1% | alignment | — | — | 2 µs | 2 µs | — | 39 µs | 18 µs | 57 µs | 152 µs | 48 µs | 4 µs | 3 µs | ✓ |
| edit-1k-5% | score | — | — | 3 µs | 3 µs | — | 37 µs | 20 µs | 22 µs | — | 12 µs | 3 µs | 3 µs | ✓ |
| edit-1k-5% | alignment | — | — | 5 µs | 5 µs | — | 40 µs | 18 µs | 58 µs | 140 µs | 50 µs | 8 µs | 5 µs | ✓ |
| edit-1k-15% | score | — | — | 11 µs | 11 µs | — | 31 µs | 20 µs | 22 µs | — | 32 µs | 15 µs | 15 µs | ✓ |
| edit-1k-15% | alignment | — | — | 19 µs | 19 µs | — | 48 µs | 31 µs | 60 µs | 506 µs | 74 µs | 30 µs | 24 µs | ✓ |
| edit-10k-1% | score | — | — | 11 µs | 11 µs | — | 492 µs | 932 µs | 1.69 ms | — | 240 µs | 14 µs | 14 µs | ✓ |
| edit-10k-1% | alignment | — | — | 23 µs | 23 µs | — | 501 µs | 152 µs | 2.89 ms | 1.6 ms | 1.48 ms | 45 µs | 35 µs | ✓ |
| edit-10k-5% | score | — | — | 106 µs | 92 µs | — | 500 µs | 1.63 ms | 1.68 ms | — | 372 µs | 170 µs | 171 µs | ✓ |
| edit-10k-5% | alignment | — | — | 136 µs | 122 µs | — | 541 µs | 263 µs | 2.98 ms | 1.5 ms | 1.85 ms | 372 µs | 288 µs | ✓ |
| edit-10k-15% | score | — | — | 209 µs | 163 µs | — | 870 µs | 2.79 ms | 1.69 ms | — | 1.18 ms | 1.16 ms | 1.16 ms | ✓ |
| edit-10k-15% | alignment | — | — | 299 µs | 227 µs | — | 1.1 ms | 767 µs | 3.19 ms | 39.5 ms | 3.34 ms | 2.42 ms | 1.95 ms | ✓ |
| edit-100k-1% | score | — | — | 605 µs | 428 µs | — | 4.28 ms | 185 ms | 163 ms | — | 5.64 ms | 692 µs | 692 µs | ✓ |
| edit-100k-1% | alignment | — | — | 907 µs | 672 µs | — | 4.42 ms | 3.33 ms | 900 ms | 19.9 ms | 32.3 ms | 1.53 ms | 1.49 ms | ✓ |
| edit-100k-5% | score | — | — | 2.38 ms | 2.36 ms | — | 4.59 ms | 275 ms | 164 ms | — | 45.4 ms | 16.4 ms | 16 ms | ✓ |
| edit-100k-5% | alignment | — | — | 2.79 ms | 2.79 ms | — | 5 ms | 16.2 ms | 854 ms | 18.2 ms | 106 ms | 31.9 ms | 25.8 ms | ✓ |
| edit-100k-15% | score | — | — | 11.1 ms | 9.26 ms | — | 21.9 ms | 245 ms | 184 ms | — | 71.9 ms | 126 ms | 126 ms | ✓ |
| edit-100k-15% | alignment | — | — | 12.4 ms | 10.1 ms | — | 23.9 ms | 25.4 ms | 726 ms | 7.27 s | 208 ms | 265 ms | 218 ms | ✓ |

## Reading the Numbers

- **Every answer agreed**, so each row compares tools answering the same question. The answers are the optimal scores and distances: a sum and a position-weighted sum per workload, so a reordered batch cannot pass.
- **On edit distance, dinara-align on one thread beats or matches every exact aligner on every workload**, distance and alignment alike: 23 against WFA's 35 µs and A\*PA2-simple's 152 µs aligning 10 kbp at 1%, 136 against A\*PA2-simple's 263 µs at 10 kbp and 5%, 2.8 against A\*PA2-full's 5.0 ms at 100 kbp and 5%, 12.4 against A\*PA2-full's 23.9 ms at 100 kbp and 15%, and 605 against WFA's 692 µs scoring 100 kbp at 1%. The ties are the 1 kbp pairs at 1 and 5%, a few microseconds each, where BiWFA and WFA run the same diagonal transition.
- **Near-identical pairs run diagonal transition**, as WFA does, before any band: one front keeping its history for a close alignment, and two from both ends otherwise, as BiWFA scores, keeping both histories for an alignment and tracing it back through each from where they met. It stays on while it costs less than the band would, which at 1 to 3% divergence covers most pairs up to tens of kbp.
- **A band re-aims its first bound as it sweeps.** The bound starts from a projection of the first few edits, which strays by up to half either way. At an eighth, a quarter and half of the columns, the band projects the distance from its own climb, hundreds of edits in, lowering the bound when it was set too high and giving the round up early when it was set too low. Only a distance within the final bound is accepted, so the answer stays exact.
- **Long, moderately divergent pairs prune with A\*PA2-full's seed heuristic.** From about 1,500 projected edits up to one in seven bases, the first sequence is cut into 12-base seeds, their exact matches in the second are found by a rolling hash, and matches that cannot shorten any path over the next 14 seeds are dropped, as A\*PA's local pruning drops them. The band then keeps a row only while its score plus the seeds still ahead, less the longest chain of matches it can still reach, fits the bound. At 2 to 3% divergence that bound at the origin lands within 1% of the distance, so the band starts just past it and finishes in one round: the 100 kbp pair at 5% drops from 5.3 to 3.0 ms. A\*PA2-full also prunes matches between rounds, which dinara-align leaves out, as its rounds now rarely number more than two.
- **Global affine scores run a wavefront.** Under a table of one match and one mismatch score, as `Scoring.dna()` is, the match reward folds away for a global alignment (Eizenga and Lindquist), and WFA's three-layer wavefront then searches by cost: 438 ms against hyalite's 20.3 s at 100 kbp, 109 µs against 1.24 ms at 1 kbp. A pair whose projected wavefront would cost more than a full sweep is handed to the sweep, which runs sixteen cells at a time by anti-diagonal and also serves local scores: 39 ms against hyalite's 181 to 198 ms on a 10 kbp pair at 30% divergence or unrelated.
- **Affine alignments still sweep the whole matrix**, cell by cell on the host, with the tie-breaking the GPU reproduces exactly. hyalite's striped kernel edges it at 1 kbp, 4.26 against 4.46 ms; dinara-align leads at 10 and 100 kbp, and its GPU leads everywhere it is not dominated by setup.
- **The other rivals run on one CPU thread.** dinara-align's CPU batches spread their pairs over every thread, one pair a thread, which the reads rows show; the single affine pairs run on one. The bit-parallel all-threads column sweeps a long unseeded pair from both ends at once, one direction per thread, and runs the two-ended diagonal transition's fronts on two threads once they are long: 428 against one thread's 605 µs scoring 100 kbp at 1%.
- **dinara-align's GPU times include the overheads a caller pays**: opening the device context on every call, copying sequences over, and copying results back. That is why a single short pair is slower on the GPU than on the CPU.
- **hyalite's traceback budget is 1 GiB.** It keeps the whole matrix when it fits and switches to a checkpointed sweep above that, as the 100 kbp pairs do. dinara-align switches to its linear-space recursion above six million cells on the host and one million on the device.
- **A\*PA answers only edit distance**, and its time depends on how similar the sequences are, as does the bit-parallel path's and now dinara-align's global affine score; the affine alignments and hyalite's columns do not.
- **One machine, one fixed pair per workload.** On fresh pairs from two other seeds, from 1 to 100 kbp at 2 to 30% divergence, dinara-align on one thread beat or matched every exact aligner here, A\*PA2, A\*PA, Edlib, BiWFA and WFA, on every pair, distance and alignment alike. Rerun on your own hardware before quoting any of this.
