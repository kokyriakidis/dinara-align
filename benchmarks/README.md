# Benchmarks

Times dinara-align against other exact global DNA aligners on the same inputs, and refuses to report a time for any workload where the tools disagree on the answer.

```bash
pixi run bench                    # 1 kbp and 10 kbp workloads, about a minute
pixi run bench --full             # adds the 100 kbp pairs, about four and a half minutes
pixi run bench --install-rust     # also installs the nightly Rust A*PA needs, under .cache/
```

The table prints to the terminal and lands in `.cache/results/results.md`, with every raw row in `results.tsv`.

```bash
pixi run bench-astarpa2           # A*PA2's own evaluation datasets, about fifteen seconds warm
pixi run bench-astarpa2 --fresh   # the rivals too, rather than their kept results, about three minutes
```

That table lands in `.cache/results/astarpa2.md`; see [A\*PA2's Evaluation](#apa2s-evaluation).

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

Measured with `pixi run bench --full`, which took about three minutes, on an Apple M2 (4 performance and 4 efficiency cores, 24 GB) under macOS 27.0.1, with Mojo 1.1.0 and Xcode 27.0.
Hyalite built with Rust 1.92, A\*PA with the nightly its repository pins.
A dash marks a task the tool does not offer, or a workload it is not run on.

| workload | task | dinara-align (cpu) | dinara-align (gpu) | dinara-align (bit-parallel, 1 thread) | dinara-align (bit-parallel, 8 threads) | hyalite | a*pa2-full | a*pa2-simple | a*pa2-nw | a*pa | edlib | biwfa | wfa | agree |
| :-- | :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| reads-150bp | score | 15.1 ms | 15.6 ms | — | — | 300 ms | — | — | — | — | — | — | — | ✓ |
| reads-150bp | alignment | 37.6 ms | 147 ms | — | — | 950 ms | — | — | — | — | — | — | — | ✓ |
| reads-1kbp | score | 67.6 ms | 39 ms | — | — | 1.26 s | — | — | — | — | — | — | — | ✓ |
| reads-1kbp | alignment | 103 ms | 271 ms | — | — | 5.42 s | — | — | — | — | — | — | — | ✓ |
| affine-1k | score | 126 µs | 5.24 ms | — | — | 1.24 ms | — | — | — | — | — | — | — | ✓ |
| affine-1k | alignment | 217 µs | 8.7 ms | — | — | 4.02 ms | — | — | — | — | — | — | — | ✓ |
| affine-10k | score | 5.16 ms | 18.6 ms | — | — | 206 ms | — | — | — | — | — | — | — | ✓ |
| affine-10k | alignment | 70 ms | 34 ms | — | — | 917 ms | — | — | — | — | — | — | — | ✓ |
| affine-100k | score | 430 ms | 360 ms | — | — | 20.3 s | — | — | — | — | — | — | — | ✓ |
| affine-100k | alignment | 8.96 s | 719 ms | — | — | 84.4 s | — | — | — | — | — | — | — | ✓ |
| edit-1k-1% | score | — | — | 1 µs | 1 µs | — | 37 µs | 20 µs | 22 µs | — | 13 µs | 1 µs | 1 µs | ✓ |
| edit-1k-1% | alignment | — | — | 2 µs | 2 µs | — | 39 µs | 18 µs | 57 µs | 150 µs | 47 µs | 4 µs | 3 µs | ✓ |
| edit-1k-5% | score | — | — | 3 µs | 3 µs | — | 37 µs | 20 µs | 22 µs | — | 12 µs | 3 µs | 3 µs | ✓ |
| edit-1k-5% | alignment | — | — | 5 µs | 5 µs | — | 40 µs | 18 µs | 58 µs | 143 µs | 49 µs | 8 µs | 6 µs | ✓ |
| edit-1k-15% | score | — | — | 12 µs | 13 µs | — | 32 µs | 20 µs | 22 µs | — | 32 µs | 15 µs | 15 µs | ✓ |
| edit-1k-15% | alignment | — | — | 20 µs | 20 µs | — | 48 µs | 31 µs | 60 µs | 507 µs | 74 µs | 30 µs | 24 µs | ✓ |
| edit-10k-1% | score | — | — | 11 µs | 11 µs | — | 491 µs | 932 µs | 1.69 ms | — | 243 µs | 14 µs | 15 µs | ✓ |
| edit-10k-1% | alignment | — | — | 24 µs | 23 µs | — | 501 µs | 153 µs | 2.93 ms | 1.6 ms | 1.53 ms | 47 µs | 35 µs | ✓ |
| edit-10k-5% | score | — | — | 115 µs | 94 µs | — | 497 µs | 1.63 ms | 1.68 ms | — | 387 µs | 175 µs | 172 µs | ✓ |
| edit-10k-5% | alignment | — | — | 150 µs | 139 µs | — | 540 µs | 266 µs | 3.05 ms | 1.5 ms | 1.94 ms | 372 µs | 290 µs | ✓ |
| edit-10k-15% | score | — | — | 222 µs | 165 µs | — | 880 µs | 2.8 ms | 1.69 ms | — | 1.19 ms | 1.18 ms | 1.17 ms | ✓ |
| edit-10k-15% | alignment | — | — | 321 µs | 236 µs | — | 1.11 ms | 773 µs | 3.28 ms | 39.1 ms | 3.38 ms | 2.43 ms | 1.97 ms | ✓ |
| edit-100k-1% | score | — | — | 634 µs | 432 µs | — | 4.27 ms | 183 ms | 232 ms | — | 5.66 ms | 689 µs | 687 µs | ✓ |
| edit-100k-1% | alignment | — | — | 913 µs | 675 µs | — | 4.38 ms | 3.33 ms | 851 ms | 19.6 ms | 32.3 ms | 1.53 ms | 1.51 ms | ✓ |
| edit-100k-5% | score | — | — | 2.49 ms | 2.55 ms | — | 4.61 ms | 280 ms | 162 ms | — | 46.2 ms | 16.9 ms | 16.8 ms | ✓ |
| edit-100k-5% | alignment | — | — | 2.97 ms | 2.94 ms | — | 4.99 ms | 16.2 ms | 757 ms | 17.7 ms | 105 ms | 32.5 ms | 25.9 ms | ✓ |
| edit-100k-15% | score | — | — | 6.99 ms | 7.71 ms | — | 21.8 ms | 247 ms | 163 ms | — | 72.2 ms | 127 ms | 139 ms | ✓ |
| edit-100k-15% | alignment | — | — | 8.16 ms | 8.43 ms | — | 23.8 ms | 25.2 ms | 469 ms | 6.92 s | 235 ms | 266 ms | 274 ms | ✓ |

## Reading the Numbers

- **Every answer agreed**, so each row compares tools answering the same question. The answers are the optimal scores and distances: a sum and a position-weighted sum per workload, so a reordered batch cannot pass.
- **On edit distance, dinara-align on one thread beats or matches every exact aligner on every workload**, distance and alignment alike: 24 against WFA's 35 µs and A\*PA2-simple's 153 µs aligning 10 kbp at 1%, 150 against A\*PA2-simple's 266 µs at 10 kbp and 5%, 3.0 against A\*PA2-full's 5.0 ms at 100 kbp and 5%, 8.2 against A\*PA2-full's 23.8 ms at 100 kbp and 15%, and 634 against WFA's 687 µs scoring 100 kbp at 1%. The ties are the 1 kbp pairs at 1 and 5%, a few microseconds each, where BiWFA and WFA run the same diagonal transition.
- **Near-identical pairs run diagonal transition**, as WFA does, before any band: one front keeping its history for a close alignment, and two from both ends otherwise, as BiWFA scores, keeping both histories for an alignment and tracing it back through each from where they met. It stays on while it costs less than the band would, which at 1 to 3% divergence covers most pairs up to tens of kbp.
- **A band re-aims its first bound as it sweeps.** The bound starts from a projection of the first few edits, which strays by up to half either way. At an eighth, a quarter and half of the columns, the band projects the distance from its own climb, hundreds of edits in, lowering the bound when it was set too high and giving the round up early when it was set too low. Only a distance within the final bound is accepted, so the answer stays exact.
- **Long, moderately divergent pairs prune with A\*PA2-full's seed heuristic.** From about 1,500 projected edits up to one in seven bases, the first sequence is cut into 12-base seeds, their exact matches in the second are found by a rolling hash, and matches that cannot shorten any path over the next 14 seeds are dropped, as A\*PA's local pruning drops them. The band then keeps a row only while its score plus the seeds still ahead, less the longest chain of matches it can still reach, fits the bound. At 2 to 3% divergence that bound at the origin lands within 1% of the distance, so the band starts just past it and finishes in one round: the 100 kbp pair at 5% drops from 5.3 to 3.0 ms. A\*PA2-full also prunes matches between rounds, which dinara-align leaves out, as its rounds now rarely number more than two.
- **Long, divergent pairs match seeds within one edit**, as A\*PA's `r = 2` does. Past about one edit in fifteen bases most exact 12-base seeds are broken and the bound at the origin falls far short: 8,277 for a distance of 12,196 on a 100 kbp pair at 15%. Seeds of 16 bases that may match with one edit, charging two edits where they match nowhere, put it at 11,304. Each match is found by the half of it that matches exactly, a check on its quarters turns most chance lookups away, and an exact match's one-edit neighbours are left out, which keeps the bound a lower bound. They are used from 64 kbp, where the band they save outgrows their setup, when the projection says one edit in ten bases or under 40% of the exact seeds chain: the 100 kbp alignment at 15% drops from 12.4 to 8.2 ms.
- **Global affine scores run a wavefront.** Under a table of one match and one mismatch score, as `Scoring.dna()` is, the match reward folds away for a global alignment (Eizenga and Lindquist), and WFA's three-layer wavefront then searches by cost: 430 ms against hyalite's 20.3 s at 100 kbp, 126 µs against 1.24 ms at 1 kbp. A pair whose projected wavefront would cost more than a full sweep is handed to the sweep, which runs sixteen cells at a time by anti-diagonal and also serves local scores: 39 ms against hyalite's 181 to 198 ms on a 10 kbp pair at 30% divergence or unrelated.
- **Affine alignments sweep the whole matrix sixteen cells at a time**, by anti-diagonal under the same uniform table, with the same values in every cell as the cell-by-cell recurrence and so the same gapped strings the GPU returns. A pair up to six million cells stores its three layers and walks them back; a larger one splits on rows in linear space, Myers and Miller's way, its sweeps vectorized alike: 217 µs against hyalite's 4.02 ms at 1 kbp, 8.96 against 84.4 s at 100 kbp. The GPU still leads on the 10 and 100 kbp alignments.
- **The other rivals run on one CPU thread.** dinara-align's CPU batches spread their pairs over every thread, one pair a thread, which the reads rows show; the single affine pairs run on one. The bit-parallel all-threads column sweeps a long unseeded pair from both ends at once, one direction per thread, and runs the two-ended diagonal transition's fronts on two threads once they are long: 432 against one thread's 634 µs scoring 100 kbp at 1%.
- **dinara-align's GPU times include the overheads a caller pays**: opening the device context on every call, copying sequences over, and copying results back. That is why a single short pair is slower on the GPU than on the CPU.
- **hyalite's traceback budget is 1 GiB.** It keeps the whole matrix when it fits and switches to a checkpointed sweep above that, as the 100 kbp pairs do. dinara-align switches to its linear-space recursion above six million cells on the host and one million on the device.
- **A\*PA answers only edit distance**, and its time depends on how similar the sequences are, as does the bit-parallel path's and now dinara-align's global affine score; the affine alignments and hyalite's columns do not.
- **One machine, one fixed pair per workload.** On fresh pairs from two other seeds, from 1 to 100 kbp at 2 to 30% divergence, dinara-align on one thread beat or matched every exact aligner here, A\*PA2, A\*PA, Edlib, BiWFA and WFA, on every pair, distance and alignment alike. Rerun on your own hardware before quoting any of this.

## A\*PA2's Evaluation

`pa_bench.py` runs the exact aligners of A\*PA2's evaluation, with its parameters, on its own datasets: Oxford Nanopore reads and SARS-CoV-2 genomes from pa-bench's release, and the uniform-error pairs pa-generate regenerates exactly from the evaluation's seed and lock.
Each dataset is a fixed shuffled sample of about 2 Mbp and at least four pairs, the same for every tool; each tool aligns its pairs with traceback, once each, as pa-bench times them, within a budget of five seconds, and every cost is checked against every other tool's and against the costs A\*PA2's published results recorded.
The rivals are pinned, so their results are kept and reused until `--fresh`; a second table gives each tool's peak resident memory.
dinara-align also runs on all eight threads two ways: one pair at a time across them, and as a batch, `edit_alignments` over the whole sample, one pair a thread, whose column is the batch's wall-clock time over its pairs, a throughput where the others are each one pair's latency.

Measured as above, dinara-align's columns the fastest of three warm runs; a count marks a tool its budget stopped partway, and more than the budget one that finished no pair:

| dataset | pairs | mean length | dinara-align (bit-parallel, 1 thread) | dinara-align (bit-parallel, 8 threads) | dinara-align (bit-parallel, batch, 8 threads) | a*pa2-full | a*pa2-simple | a*pa | edlib | biwfa | agree |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| ont-1k | 1221 of 12477 | 0.818 kbp | 27 µs | 28 µs | 7 µs | 64 µs | 40 µs | 498 µs | 71 µs | 27 µs | ✓ |
| ont-10k | 277 of 5000 | 3.6 kbp | 162 µs | 150 µs | 47 µs | 293 µs | 205 µs | 8.54 ms | 705 µs | 367 µs | ✓ |
| ont-50k | 104 of 10000 | 9.52 kbp | 735 µs | 716 µs | 264 µs | 999 µs | 958 µs | 93 ms (54/104) | 4.62 ms | 3.49 ms | ✓ |
| ont-500k | 4 of 50 | 638 kbp | 144 ms | 145 ms | 101 ms | 280 ms | 805 ms | > 5 s | > 5 s | > 5 s | ✓ |
| ont-500k-genvar | 4 of 48 | 659 kbp | 207 ms | 207 ms | 88.6 ms | 269 ms | 627 ms | > 5 s | 4.63 s (1/4) | 4.18 s (1/4) | ✓ |
| sars-cov-2 | 33 of 10000 | 29.6 kbp | 242 µs | 255 µs | 114 µs | 1.39 ms | 700 µs | 3.42 ms | 6.95 ms | 460 µs | ✓ |
| Uniform-t10000000-n3000-e0.05 | 333 of 3333 | 3 kbp | 32 µs | 40 µs | 10 µs | 149 µs | 65 µs | 177 µs | 344 µs | 51 µs | ✓ |
| Uniform-t10000000-n10000-e0.05 | 99 of 1000 | 10 kbp | 198 µs | 199 µs | 60 µs | 581 µs | 298 µs | 569 µs | 1.92 ms | 290 µs | ✓ |
| Uniform-t10000000-n30000-e0.05 | 33 of 333 | 30 kbp | 909 µs | 909 µs | 274 µs | 1.53 ms | 1.86 ms | 1.61 ms | 11.7 ms | 2.49 ms | ✓ |
| Uniform-t10000000-n100000-e0.05 | 10 of 100 | 100 kbp | 2.94 ms | 2.92 ms | 1.05 ms | 5 ms | 18 ms | 6.85 ms | 107 ms | 32.6 ms | ✓ |
| Uniform-t10000000-n300000-e0.05 | 4 of 33 | 300 kbp | 12.5 ms | 12.6 ms | 3.76 ms | 16.3 ms | 90 ms | 19.1 ms | 627 ms | 245 ms | ✓ |
| Uniform-t10000000-n1000000-e0.05 | 4 of 10 | 1e+03 kbp | 41.2 ms | 40.8 ms | 12.2 ms | 75 ms | 1.03 s | 69.9 ms | > 5 s | 2.43 s (2/4) | ✓ |
| Uniform-t10000000-n3000-e0.15 | 333 of 3333 | 3 kbp | 81 µs | 89 µs | 22 µs | 166 µs | 135 µs | 1.08 ms | 423 µs | 190 µs | ✓ |
| Uniform-t10000000-n10000-e0.15 | 99 of 1000 | 10 kbp | 351 µs | 304 µs | 97 µs | 800 µs | 757 µs | 4.06 ms | 3.35 ms | 1.9 ms | ✓ |
| Uniform-t10000000-n30000-e0.15 | 33 of 333 | 30 kbp | 2.1 ms | 1.86 ms | 585 µs | 3.75 ms | 2.62 ms | 14.2 ms | 19.3 ms | 15.8 ms | ✓ |
| Uniform-t10000000-n100000-e0.15 | 9 of 100 | 100 kbp | 8.93 ms | 8.94 ms | 3.11 ms | 14.5 ms | 27.2 ms | 53.1 ms | 199 ms | 196 ms | ✓ |
| Uniform-t10000000-n300000-e0.15 | 4 of 33 | 300 kbp | 46.4 ms | 46.7 ms | 13.1 ms | 134 ms | 345 ms | 224 ms | 2.11 s (2/4) | 1.72 s (3/4) | ✓ |
| Uniform-t10000000-n1000000-e0.15 | 4 of 10 | 1e+03 kbp | 388 ms | 387 ms | 107 ms | 1.75 s (3/4) | 1.85 s (3/4) | 951 ms | > 5 s | > 5 s | ✓ |

Four long reads are few: on all 50 ont-500k reads dinara-align on one thread took 7.5 s and A\*PA2-full 8.4 s, and on all 48 ont-500k-genvar reads 11.8 s against A\*PA2-full's 11.6 s, A\*PA2-full the faster on 31 of them, most of them between 4 and 8% divergent.
