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
| ont-1k | 1221 of 12477 | 0.818 kbp | 27 µs | 27 µs | 7 µs | 64 µs | 40 µs | 498 µs | 71 µs | 27 µs | ✓ |
| ont-10k | 277 of 5000 | 3.6 kbp | 153 µs | 146 µs | 42 µs | 293 µs | 205 µs | 8.54 ms | 705 µs | 367 µs | ✓ |
| ont-50k | 104 of 10000 | 9.52 kbp | 685 µs | 677 µs | 184 µs | 999 µs | 958 µs | 93 ms (54/104) | 4.62 ms | 3.49 ms | ✓ |
| ont-500k | 4 of 50 | 638 kbp | 140 ms | 104 ms | 100 ms | 280 ms | 805 ms | > 5 s | > 5 s | > 5 s | ✓ |
| ont-500k-genvar | 4 of 48 | 659 kbp | 204 ms | 143 ms | 87.2 ms | 269 ms | 627 ms | > 5 s | 4.63 s (1/4) | 4.18 s (1/4) | ✓ |
| sars-cov-2 | 33 of 10000 | 29.6 kbp | 236 µs | 246 µs | 106 µs | 1.39 ms | 700 µs | 3.42 ms | 6.95 ms | 460 µs | ✓ |
| Uniform-t10000000-n3000-e0.05 | 333 of 3333 | 3 kbp | 33 µs | 40 µs | 10 µs | 149 µs | 65 µs | 177 µs | 344 µs | 51 µs | ✓ |
| Uniform-t10000000-n10000-e0.05 | 99 of 1000 | 10 kbp | 203 µs | 202 µs | 59 µs | 581 µs | 298 µs | 569 µs | 1.92 ms | 290 µs | ✓ |
| Uniform-t10000000-n30000-e0.05 | 33 of 333 | 30 kbp | 908 µs | 922 µs | 277 µs | 1.53 ms | 1.86 ms | 1.61 ms | 11.7 ms | 2.49 ms | ✓ |
| Uniform-t10000000-n100000-e0.05 | 10 of 100 | 100 kbp | 3.01 ms | 2.87 ms | 1.01 ms | 5 ms | 18 ms | 6.85 ms | 107 ms | 32.6 ms | ✓ |
| Uniform-t10000000-n300000-e0.05 | 4 of 33 | 300 kbp | 12.6 ms | 11.1 ms | 3.7 ms | 16.3 ms | 90 ms | 19.1 ms | 627 ms | 245 ms | ✓ |
| Uniform-t10000000-n1000000-e0.05 | 4 of 10 | 1e+03 kbp | 41.6 ms | 36.6 ms | 12.2 ms | 75 ms | 1.03 s | 69.9 ms | > 5 s | 2.43 s (2/4) | ✓ |
| Uniform-t10000000-n3000-e0.15 | 333 of 3333 | 3 kbp | 84 µs | 96 µs | 22 µs | 166 µs | 135 µs | 1.08 ms | 423 µs | 190 µs | ✓ |
| Uniform-t10000000-n10000-e0.15 | 99 of 1000 | 10 kbp | 352 µs | 305 µs | 98 µs | 800 µs | 757 µs | 4.06 ms | 3.35 ms | 1.9 ms | ✓ |
| Uniform-t10000000-n30000-e0.15 | 33 of 333 | 30 kbp | 2.1 ms | 1.86 ms | 612 µs | 3.75 ms | 2.62 ms | 14.2 ms | 19.3 ms | 15.8 ms | ✓ |
| Uniform-t10000000-n100000-e0.15 | 9 of 100 | 100 kbp | 9.06 ms | 8.63 ms | 3.17 ms | 14.5 ms | 27.2 ms | 53.1 ms | 199 ms | 196 ms | ✓ |
| Uniform-t10000000-n300000-e0.15 | 4 of 33 | 300 kbp | 46.2 ms | 31.6 ms | 13.2 ms | 134 ms | 345 ms | 224 ms | 2.11 s (2/4) | 1.72 s (3/4) | ✓ |
| Uniform-t10000000-n1000000-e0.15 | 4 of 10 | 1e+03 kbp | 385 ms | 190 ms | 104 ms | 1.75 s (3/4) | 1.85 s (3/4) | 951 ms | > 5 s | > 5 s | ✓ |

Four long reads are few: on all 50 ont-500k reads dinara-align took 7.5 s on one thread and 6.5 s on eight, and A\*PA2-full 8.4 s; on all 48 ont-500k-genvar reads 11.8 s on one thread, level with A\*PA2-full's 11.6 s (A\*PA2-full the faster on 31 of them, most between 4 and 8% divergent), and 8.4 s on eight, where pairs from 200 kbp whose seeds match within one edit or chain poorly sweep their bands in stripes across the threads.

## A\*PA2's Results, Redone

A\*PA2's own results section ([curiouscoding.nl/posts/astarpa2](https://curiouscoding.nl/posts/astarpa2/#results)), redone with dinara-align beside the exact aligners it compares, on a machine set up as A\*PA2's was: its real datasets, its uniform pairs swept over divergence and over length, and in place of its ablation, dinara-align's history on the long reads.

```bash
pixi run results-astarpa2 collect   # forty minutes the first time; after that only dinara-align runs, a few minutes
pixi run results-astarpa2 tables    # the tables below, into .cache/results/astarpa2-results.md
```

### Setup

As A\*PA2's evaluation runs them: one single-threaded job at a time, every pair aligned once with its traceback, the time the average wall clock per alignment, reading the data left out.

The machine is an Intel Core i9-7900X under Ubuntu 26.04 (Linux 7.0), with Mojo 1.1.0, set up as A\*PA2's i7-10750H was: every core fixed at 3.3 GHz, turbo boost and hyper-threading off, and each collection pinned to one core with `taskset`.
Two things differ: the jobs ran at normal priority, not A\*PA2's niceness −20, and a job of the machine's owner ran on another core during part of the rivals' runs, which may have slowed those a little.

Each dataset is a fixed shuffled sample, the same pairs for every aligner: about 2 Mbp of sequence and at least four pairs, and 20 Mbp, fifteen reads, of the two 500 kbp read sets.
Each aligner may spend 20 seconds on a sample; one that has not finished stops there, its numbers covering the pairs it finished, counted beside them, and one that finished none shows a dash.
Every cost is checked against every other aligner's on the pairs both aligned, and against the costs A\*PA2's published results recorded; a disagreement fails the run, and none occurred.
The rivals, and dinara-align's older commits, never change, so each is run on a sample once and kept; later runs time today's dinara-align alone.
The same tables measured on the Apple M2 above, without any of this pinning, close the section.

### Real datasets

Mean time per alignment, median in brackets:

| dataset | pairs | mean length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: |
| ont-1k | 1221 | 0.818 kbp | 38 µs (31 µs) | 90 µs (101 µs) | 57 µs (64 µs) | 603 µs (544 µs) | 136 µs (133 µs) | 44 µs (46 µs) |
| ont-10k | 277 | 3.6 kbp | 236 µs (152 µs) | 464 µs (309 µs) | 340 µs (233 µs) | 15.4 ms (3.22 ms) | 1.07 ms (582 µs) | 746 µs (282 µs) |
| ont-50k | 104 | 9.52 kbp | 1.04 ms (303 µs) | 1.66 ms (600 µs) | 1.54 ms (452 µs) | 262 ms (11.3 ms), 78/104 | 5.82 ms (1.71 ms) | 7.3 ms (862 µs) |
| ont-500k | 15 | 638 kbp | 171 ms (120 ms) | 303 ms (150 ms) | 1.15 s (795 ms) | — | 4.72 s (2.88 s), 4/15 | 19.8 s (19.8 s), 1/15 |
| ont-500k-genvar | 15 | 651 kbp | 230 ms (188 ms) | 367 ms (264 ms) | 1.33 s (1.09 s) | — | 4.94 s (4.91 s), 4/15 | 6.35 s (7.55 s), 3/15 |
| sars-cov-2 | 33 | 29.6 kbp | 399 µs (229 µs) | 2.43 ms (2.51 ms) | 1.14 ms (1.1 ms) | 6.54 ms (2.03 ms) | 8.17 ms (7.78 ms) | 888 µs (438 µs) |

- **dinara-align has the lowest mean and the lowest median on every dataset.** On the short reads, the case A\*PA2's post leaves to BiWFA, 38 against BiWFA's 44 µs, median 31 against 46; on ont-10k and ont-50k it takes about two thirds of A\*PA2-simple's time, and on the SARS-CoV-2 genomes 399 against BiWFA's 888 µs.
- **On the 500 kbp reads it takes 56 to 63% of A\*PA2-full's time:** 171 against 303 ms on ont-500k and 230 against 367 ms on ont-500k-genvar, medians 120 against 150 and 188 against 264.
- **Edlib and BiWFA finish only some of the long reads within 20 seconds, and A\*PA none;** their means cover the reads they finished.

### Divergence, 100 kbp pairs

Mean time per alignment, A\*PA with whichever of `r = 1` and `r = 2` was the faster at each divergence:

| divergence | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: |
| 0% | 147 µs | 8.02 ms | 2.57 ms | 6.72 ms | 42.5 ms | 621 µs |
| 1% | 1.65 ms | 9.4 ms | 5.98 ms | 8.2 ms | 57.6 ms | 3.28 ms |
| 2% | 5.48 ms | 9.66 ms | 8.97 ms | 8.75 ms | 70.7 ms | 10.2 ms |
| 3% | 5.12 ms | 9.73 ms | 16.4 ms | 8.44 ms | 94.4 ms | 21.3 ms |
| 4% | 5.2 ms | 9.93 ms | 14.2 ms | 10.4 ms | 97.1 ms | 36 ms |
| 5% | 5.7 ms | 10.4 ms | 31.7 ms | 13 ms | 141 ms | 54.4 ms |
| 6% | 11.3 ms | 10.8 ms | 28.7 ms | 18.9 ms | 143 ms | 74.3 ms |
| 7% | 11 ms | 10.7 ms | 26.7 ms | 39.5 ms | 146 ms | 97.9 ms |
| 8% | 11.3 ms | 15.4 ms | 24.6 ms | 46.8 ms | 150 ms | 124 ms |
| 9% | 11.3 ms | 17 ms | 22.9 ms | 47 ms | 153 ms | 154 ms |
| 10% | 11.8 ms | 14.7 ms | 59.2 ms | 53 ms | 234 ms | 186 ms |
| 11% | 11.8 ms | 21.5 ms | 56.6 ms | 70.5 ms | 237 ms | 221 ms |
| 12% | 12.2 ms | 18.3 ms | 54.6 ms | 60.7 ms | 239 ms | 255 ms |
| 13% | 11.8 ms | 30.7 ms | 52.3 ms | 66.1 ms | 242 ms | 293 ms |
| 14% | 13.1 ms | 28.1 ms | 51.1 ms | 76.9 ms | 245 ms | 332 ms |
| 15% | 13.2 ms | 25.8 ms | 49.6 ms | 111 ms | 248 ms | 374 ms |

- **Near-identical pairs finish in the diagonal transition:** 146 µs against BiWFA's 621 at 0%, 1.68 against 3.28 ms at 1%.
- **From 2% to 5% the exact seeds' bound lands close to the distance,** and dinara-align runs in 5.2 to 6.0 ms, A\*PA2-full in 9.7 to 10.4.
- **At 6% and 7% A\*PA2-full is the faster,** 10.8 and 10.7 ms against 11.6 and 11.3: at 6% only 37% of the exact seeds chain, under the 40% at which dinara-align builds inexact seeds in their place, and the second seed setup is not yet repaid.
- **From 8% on dinara-align stays between 11.6 and 13.2 ms,** as inexact seeds keep the bound close, while A\*PA2-full climbs to 26 to 31 ms by 13 to 15%; at 15% A\*PA2-simple takes 50 ms, A\*PA 111, Edlib 248 and BiWFA 374.

### Length

Mean time per alignment at 5% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: |
| 3 kbp | 49 µs | 214 µs | 91 µs | 261 µs | 829 µs | 82 µs |
| 10 kbp | 300 µs | 919 µs | 465 µs | 919 µs | 3.83 ms | 554 µs |
| 30 kbp | 1.66 ms | 2.86 ms | 3.14 ms | 3.04 ms | 14.2 ms | 4.9 ms |
| 100 kbp | 5.48 ms | 9.61 ms | 31.1 ms | 11.9 ms | 139 ms | 53.4 ms |
| 300 kbp | 21.5 ms | 32.5 ms | 133 ms | 38.9 ms | 712 ms | 463 ms |
| 1 Mbp | 73.4 ms | 148 ms | 1.84 s | 160 ms | 8.17 s (2/4) | 5.38 s (3/4) |

At 15% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: |
| 3 kbp | 119 µs | 243 µs | 212 µs | 1.52 ms | 884 µs | 380 µs |
| 10 kbp | 425 µs | 1.19 ms | 1.31 ms | 6.61 ms | 5.66 ms | 3.92 ms |
| 30 kbp | 2.42 ms | 6.33 ms | 4.59 ms | 24.5 ms | 23.2 ms | 34.1 ms |
| 100 kbp | 13.3 ms | 27.2 ms | 49.6 ms | 111 ms | 248 ms | 368 ms |
| 300 kbp | 55.2 ms | 237 ms | 619 ms | 514 ms | 2.38 s | 3.37 s |
| 1 Mbp | 361 ms | 2.76 s | 2.89 s | 2.16 s | 17.7 s (1/4) | — |

- **dinara-align is the fastest at every length and both divergences.** At 5% it aligns a 1 Mbp pair in 75.6 ms against A\*PA2-full's 148 and A\*PA's 160; at 15% in 362 ms against A\*PA's 2.16 s and A\*PA2-full's 2.76 s.
- **At 15% its band still grows faster than the length,** but less than A\*PA2-full's: from 300 kbp to 1 Mbp, 3.3 times the length, dinara-align takes 6.5 times as long and A\*PA2-full 11.6 times, where A\*PA, whose pruning makes its heuristic nearly exact, takes 4.2 times.

### dinara-align's history on the long reads

A\*PA2 measures each of its methods by adding them one at a time. dinara-align's methods came in commits that also changed other things, so here instead are the same fifteen reads of each 500 kbp set aligned by the package as of the commit that added each method, a 60-second budget each. The commits before the x86 build fix take that fix's one line, the inline assembly's register constraint, and nothing else:

| commit | ont-500k | ont-500k-genvar |
| :-- | --: | --: |
| `641819a` port | 1.32 s (481 ms) | 1.19 s (854 ms) |
| `5a6ee58` + diagonal transition | 1.24 s (872 ms) | 8.54 s (10.1 s), 7/15 |
| `3851282` + seed heuristic | 1.04 s (713 ms) | 8.54 s (10.1 s), 7/15 |
| `28c1241` + local pruning | 762 ms (169 ms) | 8.1 s (10.2 s), 7/15 |
| `c53edc6` + two-ended search | 472 ms (161 ms) | 4.33 s (5.41 s), 13/15 |
| `b60ece0` + real-read fixes | 456 ms (160 ms) | 752 ms (499 ms) |
| `501a6bb` + retries aimed | 162 ms (94.4 ms) | 199 ms (157 ms) |
| `ef0e73a` + inexact seeds | 181 ms (131 ms) | 244 ms (204 ms) |
| today | 171 ms (120 ms) | 230 ms (188 ms) |

- **On ont-500k the time falls from 1.32 s at the port of A\*PA2-simple's band doubling to 162 ms once the retries are aimed.**
- **On ont-500k-genvar the commits between the port and `b60ece0` were several times slower,** 4.3 to 8.5 s a read and up to eight reads past the budget: real reads gather their errors at the ends, and the projections those commits aimed by ran many times over the distance, so they swept bands far wider than needed. `b60ece0` trusted a projection only when two of them agreed and otherwise grew the band from what was known, bringing genvar back to 752 ms, and aiming the retries to 199 ms.
- **Inexact seeds make these reads slower here,** 162 to 181 ms and 199 to 244 ms, where on the M2 they made them faster, 171 to 131 and 262 to 197: likely this machine sweeps the band faster beside the seeds' setup, so the point past which inexact seeds pay lies further out than the one tuned on the M2. They still win on the uniform pairs past 8%.

### Memory

Peak resident memory, measured with `pixi run bench-astarpa2` on the M2 (on Linux a runner inherits the harness's own peak when it is forked, which the measurement does not yet undo):

- dinara-align starts at 12 MB, the Mojo runtime, where the Rust aligners start at 2 MB, and stays under 20 MB on pairs up to 100 kbp.
- On the 500 kbp reads it peaks at 70 to 74 MB against A\*PA2-full's 116 to 146 MB, and on a 1 Mbp pair at 15% at 96 MB against 192 MB; A\*PA reaches over 800 MB.

### Beyond one thread

All of the above is one thread, as A\*PA2's evaluation measures. dinara-align also runs one pair across threads and batches of pairs across threads, which none of these aligners does; see [A\*PA2's Evaluation](#apa2s-evaluation) above for those columns, measured on the M2: on eight threads a 1 Mbp pair at 15% aligns in 190 ms rather than 385, and a batch of ont-1k reads at 7 µs a read.

### On the Apple M2

The same collection on the M2 above, unpinned, its frequency free, 2 to 6% between repeated runs. The order is mostly the same; there A\*PA2-full has the lower median on ont-500k, 75.9 against 90.9 ms, dinara-align and BiWFA tie on ont-1k, 27 against 28 µs, and inexact seeds speed the long reads up rather than down.

<details>
<summary>The M2's tables</summary>

Real datasets, mean time per alignment, median in brackets:

| dataset | pairs | mean length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: |
| ont-1k | 1221 | 0.818 kbp | 27 µs (22 µs) | 64 µs (70 µs) | 40 µs (46 µs) | 389 µs (357 µs) | 73 µs (78 µs) | 28 µs (29 µs) |
| ont-10k | 277 | 3.6 kbp | 152 µs (97 µs) | 293 µs (204 µs) | 209 µs (141 µs) | 8.8 ms (2.02 ms) | 708 µs (372 µs) | 369 µs (151 µs) |
| ont-50k | 104 | 9.52 kbp | 681 µs (208 µs) | 1.01 ms (397 µs) | 890 µs (278 µs) | 121 ms (5.53 ms) | 4.8 ms (1.24 ms) | 3.5 ms (427 µs) |
| ont-500k | 15 | 638 kbp | 129 ms (90.9 ms) | 160 ms (75.9 ms) | 635 ms (443 ms) | — | 4.16 s (2.53 s), 4/15 | 2.03 s (891 ms), 10/15 |
| ont-500k-genvar | 15 | 651 kbp | 192 ms (134 ms) | 199 ms (133 ms) | 745 ms (609 ms) | — | 4.38 s (4.33 s), 4/15 | 2.59 s (2.75 s), 8/15 |
| sars-cov-2 | 33 | 29.6 kbp | 241 µs (148 µs) | 1.38 ms (1.38 ms) | 705 µs (661 µs) | 4.12 ms (1.1 ms) | 6.97 ms (6.78 ms) | 464 µs (216 µs) |

Divergence, 100 kbp pairs:

| divergence | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: |
| 0% | 200 µs | 3.83 ms | 1.5 ms | 4.16 ms | 22.5 ms | 344 µs |
| 1% | 872 µs | 4.57 ms | 3.55 ms | 3.94 ms | 32.3 ms | 1.46 ms |
| 2% | 3.15 ms | 4.6 ms | 5.13 ms | 4.14 ms | 43.6 ms | 4.93 ms |
| 3% | 2.74 ms | 4.71 ms | 9.23 ms | 4.44 ms | 64.8 ms | 9.94 ms |
| 4% | 2.8 ms | 4.87 ms | 9 ms | 5.03 ms | 66.8 ms | 16.2 ms |
| 5% | 3.07 ms | 5 ms | 17.6 ms | 6.46 ms | 103 ms | 24.6 ms |
| 6% | 6.29 ms | 5.25 ms | 15.9 ms | 9.85 ms | 105 ms | 35.2 ms |
| 7% | 6.52 ms | 5 ms | 14.9 ms | 20.2 ms | 108 ms | 49 ms |
| 8% | 6.52 ms | 7.58 ms | 14 ms | 20.4 ms | 112 ms | 59.7 ms |
| 9% | 6.76 ms | 9.71 ms | 12.7 ms | 21.2 ms | 115 ms | 74.7 ms |
| 10% | 6.89 ms | 8.89 ms | 32.8 ms | 23.1 ms | 189 ms | 90.5 ms |
| 11% | 6.97 ms | 12.1 ms | 31.4 ms | 27.2 ms | 189 ms | 108 ms |
| 12% | 7.39 ms | 10 ms | 30.4 ms | 28.5 ms | 192 ms | 128 ms |
| 13% | 7.12 ms | 17 ms | 29.5 ms | 33 ms | 192 ms | 143 ms |
| 14% | 8.32 ms | 16.1 ms | 28.6 ms | 40.4 ms | 199 ms | 175 ms |
| 15% | 9.01 ms | 14.4 ms | 27.6 ms | 53.4 ms | 199 ms | 196 ms |

Length at 5% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: |
| 3 kbp | 38 µs | 150 µs | 66 µs | 178 µs | 346 µs | 51 µs |
| 10 kbp | 205 µs | 582 µs | 306 µs | 546 µs | 1.96 ms | 295 µs |
| 30 kbp | 913 µs | 1.53 ms | 1.85 ms | 1.61 ms | 11.7 ms | 2.3 ms |
| 100 kbp | 2.94 ms | 4.96 ms | 17.6 ms | 6.51 ms | 103 ms | 24.6 ms |
| 300 kbp | 12.5 ms | 17.3 ms | 74.5 ms | 18 ms | 619 ms | 235 ms |
| 1 Mbp | 41.1 ms | 71.4 ms | 1.04 s | 73 ms | 7.26 s (2/4) | 2.45 s |

Length at 15% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: |
| 3 kbp | 84 µs | 170 µs | 135 µs | 999 µs | 411 µs | 192 µs |
| 10 kbp | 351 µs | 804 µs | 760 µs | 3.91 ms | 3.38 ms | 1.87 ms |
| 30 kbp | 2.11 ms | 3.89 ms | 2.63 ms | 14.3 ms | 19.3 ms | 16.4 ms |
| 100 kbp | 8.91 ms | 14.4 ms | 27.3 ms | 52.7 ms | 200 ms | 196 ms |
| 300 kbp | 46.3 ms | 134 ms | 348 ms | 227 ms | 2.1 s | 1.65 s |
| 1 Mbp | 395 ms | 1.56 s | 1.6 s | 949 ms | 15.7 s (1/4) | 17.3 s (1/4) |

dinara-align's history:

| commit | ont-500k | ont-500k-genvar |
| :-- | --: | --: |
| `641819a` port | 800 ms (380 ms) | 791 ms (546 ms) |
| `5a6ee58` + diagonal transition | 821 ms (676 ms) | 5.42 s (6.98 s), 10/15 |
| `3851282` + seed heuristic | 640 ms (439 ms) | 5.52 s (6.88 s), 10/15 |
| `28c1241` + local pruning | 453 ms (105 ms) | 5.44 s (7.02 s), 10/15 |
| `c53edc6` + two-ended search | 288 ms (103 ms) | 3.53 s (3.88 s) |
| `b60ece0` + real-read fixes | 312 ms (120 ms) | 629 ms (411 ms) |
| `501a6bb` + retries aimed | 171 ms (63.1 ms) | 262 ms (212 ms) |
| `ef0e73a` + inexact seeds | 131 ms (87.1 ms) | 197 ms (154 ms) |
| today | 129 ms (90.9 ms) | 192 ms (134 ms) |

</details>
