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

| workload | task | dinara-align (cpu) | dinara-align (gpu) | dinara-align (bit-parallel, 1 thread) | hyalite | a*pa2-full | a*pa2-simple | a*pa2-nw | a*pa | edlib | biwfa | wfa | agree |
| :-- | :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| reads-150bp | score | 15.1 ms | 15.6 ms | — | 300 ms | — | — | — | — | — | — | — | ✓ |
| reads-150bp | alignment | 37.6 ms | 147 ms | — | 950 ms | — | — | — | — | — | — | — | ✓ |
| reads-1kbp | score | 67.6 ms | 39 ms | — | 1.26 s | — | — | — | — | — | — | — | ✓ |
| reads-1kbp | alignment | 103 ms | 271 ms | — | 5.42 s | — | — | — | — | — | — | — | ✓ |
| affine-1k | score | 126 µs | 5.24 ms | — | 1.24 ms | — | — | — | — | — | — | — | ✓ |
| affine-1k | alignment | 217 µs | 8.7 ms | — | 4.02 ms | — | — | — | — | — | — | — | ✓ |
| affine-10k | score | 5.16 ms | 18.6 ms | — | 206 ms | — | — | — | — | — | — | — | ✓ |
| affine-10k | alignment | 70 ms | 34 ms | — | 917 ms | — | — | — | — | — | — | — | ✓ |
| affine-100k | score | 430 ms | 360 ms | — | 20.3 s | — | — | — | — | — | — | — | ✓ |
| affine-100k | alignment | 8.96 s | 719 ms | — | 84.4 s | — | — | — | — | — | — | — | ✓ |
| edit-1k-1% | score | — | — | 1 µs | — | 37 µs | 20 µs | 22 µs | — | 13 µs | 1 µs | 1 µs | ✓ |
| edit-1k-1% | alignment | — | — | 2 µs | — | 39 µs | 18 µs | 57 µs | 150 µs | 47 µs | 4 µs | 3 µs | ✓ |
| edit-1k-5% | score | — | — | 3 µs | — | 37 µs | 20 µs | 22 µs | — | 12 µs | 3 µs | 3 µs | ✓ |
| edit-1k-5% | alignment | — | — | 5 µs | — | 40 µs | 18 µs | 58 µs | 143 µs | 49 µs | 8 µs | 6 µs | ✓ |
| edit-1k-15% | score | — | — | 12 µs | — | 32 µs | 20 µs | 22 µs | — | 32 µs | 15 µs | 15 µs | ✓ |
| edit-1k-15% | alignment | — | — | 20 µs | — | 48 µs | 31 µs | 60 µs | 507 µs | 74 µs | 30 µs | 24 µs | ✓ |
| edit-10k-1% | score | — | — | 11 µs | — | 491 µs | 932 µs | 1.69 ms | — | 243 µs | 14 µs | 15 µs | ✓ |
| edit-10k-1% | alignment | — | — | 24 µs | — | 501 µs | 153 µs | 2.93 ms | 1.6 ms | 1.53 ms | 47 µs | 35 µs | ✓ |
| edit-10k-5% | score | — | — | 115 µs | — | 497 µs | 1.63 ms | 1.68 ms | — | 387 µs | 175 µs | 172 µs | ✓ |
| edit-10k-5% | alignment | — | — | 150 µs | — | 540 µs | 266 µs | 3.05 ms | 1.5 ms | 1.94 ms | 372 µs | 290 µs | ✓ |
| edit-10k-15% | score | — | — | 222 µs | — | 880 µs | 2.8 ms | 1.69 ms | — | 1.19 ms | 1.18 ms | 1.17 ms | ✓ |
| edit-10k-15% | alignment | — | — | 321 µs | — | 1.11 ms | 773 µs | 3.28 ms | 39.1 ms | 3.38 ms | 2.43 ms | 1.97 ms | ✓ |
| edit-100k-1% | score | — | — | 634 µs | — | 4.27 ms | 183 ms | 232 ms | — | 5.66 ms | 689 µs | 687 µs | ✓ |
| edit-100k-1% | alignment | — | — | 913 µs | — | 4.38 ms | 3.33 ms | 851 ms | 19.6 ms | 32.3 ms | 1.53 ms | 1.51 ms | ✓ |
| edit-100k-5% | score | — | — | 2.49 ms | — | 4.61 ms | 280 ms | 162 ms | — | 46.2 ms | 16.9 ms | 16.8 ms | ✓ |
| edit-100k-5% | alignment | — | — | 2.97 ms | — | 4.99 ms | 16.2 ms | 757 ms | 17.7 ms | 105 ms | 32.5 ms | 25.9 ms | ✓ |
| edit-100k-15% | score | — | — | 6.99 ms | — | 21.8 ms | 247 ms | 163 ms | — | 72.2 ms | 127 ms | 139 ms | ✓ |
| edit-100k-15% | alignment | — | — | 8.16 ms | — | 23.8 ms | 25.2 ms | 469 ms | 6.92 s | 235 ms | 266 ms | 274 ms | ✓ |

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
dinara-align also runs as a batch on all eight threads, `edit_alignments` over the whole sample, one pair a thread, whose column is the batch's wall-clock time over its pairs, a throughput where the others are each one pair's latency.

Measured as above, dinara-align's columns the fastest of three warm runs; a count marks a tool its budget stopped partway, and more than the budget one that finished no pair:

| dataset | pairs | mean length | dinara-align (bit-parallel, 1 thread) | dinara-align (bit-parallel, batch, 8 threads) | a*pa2-full | a*pa2-simple | a*pa | edlib | biwfa | agree |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| ont-1k | 1221 of 12477 | 0.818 kbp | 26 µs | 9 µs | 64 µs | 40 µs | 389 µs | 73 µs | 28 µs | ✓ |
| ont-10k | 277 of 5000 | 3.6 kbp | 153 µs | 41 µs | 293 µs | 209 µs | 8.8 ms | 708 µs | 369 µs | ✓ |
| ont-50k | 104 of 10000 | 9.52 kbp | 664 µs | 180 µs | 1.01 ms | 890 µs | 118 ms (43/104) | 4.8 ms | 3.5 ms | ✓ |
| ont-500k | 4 of 50 | 638 kbp | 134 ms | 93.9 ms | 280 ms | 805 ms | > 5 s | > 5 s | > 5 s | ✓ |
| ont-500k-genvar | 4 of 48 | 659 kbp | 183 ms | 84.4 ms | 269 ms | 627 ms | > 5 s | 4.63 s (1/4) | 4.18 s (1/4) | ✓ |
| sars-cov-2 | 33 of 10000 | 29.6 kbp | 235 µs | 107 µs | 1.38 ms | 705 µs | 4.12 ms | 6.97 ms | 464 µs | ✓ |
| Uniform-t10000000-n3000-e0.05 | 333 of 3333 | 3 kbp | 39 µs | 10 µs | 150 µs | 66 µs | 178 µs | 346 µs | 51 µs | ✓ |
| Uniform-t10000000-n10000-e0.05 | 99 of 1000 | 10 kbp | 203 µs | 61 µs | 582 µs | 306 µs | 546 µs | 1.96 ms | 295 µs | ✓ |
| Uniform-t10000000-n30000-e0.05 | 33 of 333 | 30 kbp | 880 µs | 268 µs | 1.53 ms | 1.85 ms | 1.61 ms | 11.7 ms | 2.3 ms | ✓ |
| Uniform-t10000000-n100000-e0.05 | 10 of 100 | 100 kbp | 2.87 ms | 979 µs | 4.96 ms | 17.6 ms | 6.51 ms | 103 ms | 24.6 ms | ✓ |
| Uniform-t10000000-n300000-e0.05 | 4 of 33 | 300 kbp | 12.3 ms | 3.63 ms | 17.3 ms | 74.5 ms | 18 ms | 619 ms | 235 ms | ✓ |
| Uniform-t10000000-n1000000-e0.05 | 4 of 10 | 1e+03 kbp | 41.7 ms | 11.8 ms | 71.4 ms | 1.04 s | 73 ms | > 5 s | 2.46 s (2/4) | ✓ |
| Uniform-t10000000-n3000-e0.15 | 333 of 3333 | 3 kbp | 83 µs | 22 µs | 170 µs | 135 µs | 999 µs | 411 µs | 192 µs | ✓ |
| Uniform-t10000000-n10000-e0.15 | 99 of 1000 | 10 kbp | 352 µs | 98 µs | 804 µs | 760 µs | 3.91 ms | 3.38 ms | 1.87 ms | ✓ |
| Uniform-t10000000-n30000-e0.15 | 33 of 333 | 30 kbp | 2.08 ms | 609 µs | 3.89 ms | 2.63 ms | 14.3 ms | 19.3 ms | 16.4 ms | ✓ |
| Uniform-t10000000-n100000-e0.15 | 9 of 100 | 100 kbp | 8.53 ms | 3.01 ms | 14.4 ms | 27.3 ms | 52.7 ms | 200 ms | 196 ms | ✓ |
| Uniform-t10000000-n300000-e0.15 | 4 of 33 | 300 kbp | 45.1 ms | 12.6 ms | 134 ms | 348 ms | 227 ms | 2.1 s (2/4) | 1.64 s (3/4) | ✓ |
| Uniform-t10000000-n1000000-e0.15 | 4 of 10 | 1e+03 kbp | 384 ms | 105 ms | 1.56 s (3/4) | 1.61 s (3/4) | 949 ms | > 5 s | > 5 s | ✓ |

Four long reads are few: on all 50 ont-500k reads dinara-align took 6.5 s, and A\*PA2-full 8.4 s (the faster on 16 of them); on all 48 ont-500k-genvar reads 10.4 s against A\*PA2-full's 11.5 s (the faster on 21 of them).

## A\*PA2's Results, Redone

A\*PA2's own results section ([curiouscoding.nl/posts/astarpa2](https://curiouscoding.nl/posts/astarpa2/#results)), redone with dinara-align beside the exact aligners it compares, on a machine set up as A\*PA2's was: its real datasets, its uniform pairs swept over divergence and over length, and in place of its ablation, dinara-align's history on the long reads.

```bash
pixi run results-astarpa2 collect   # forty minutes the first time; after that only dinara-align runs, a few minutes
pixi run results-astarpa2 tables    # the tables below, into .cache/results/astarpa2-results.md
```

### Setup

As A\*PA2's evaluation runs them: one single-threaded job at a time, every pair aligned once with its traceback, the time the average wall clock per alignment, reading the data left out.

The machine is an Intel Core i9-7900X under Ubuntu 26.04 (Linux 7.0), with Mojo 1.1.0, set up as A\*PA2's i7-10750H was: every core fixed at 3.3 GHz, turbo boost and hyper-threading off, and each collection pinned to one core with `taskset`.
One thing differs: the jobs ran at normal priority, not A\*PA2's niceness −20, on an otherwise idle machine.

Each dataset is a fixed shuffled sample, the same pairs for every aligner: about 2 Mbp of sequence and at least four pairs, and 20 Mbp, fifteen reads, of the two 500 kbp read sets.
Each aligner may spend 20 seconds on a sample; one that has not finished stops there, its numbers covering the pairs it finished, counted beside them, and one that finished none shows a dash.
Every cost is checked against every other aligner's on the pairs both aligned, and against the costs A\*PA2's published results recorded; a disagreement fails the run, and none occurred.
The rivals, and dinara-align's older commits, never change, so each is run on a sample once and kept; later runs time today's dinara-align alone.
The same tables measured on the Apple M2 above, without any of this pinning, close the section.

### Real datasets

Mean time per alignment, median in brackets:

| dataset | pairs | mean length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: |
| ont-1k | 1221 | 0.818 kbp | 38 µs (32 µs) | 91 µs (101 µs) | 57 µs (65 µs) | 604 µs (546 µs) | 136 µs (135 µs) | 45 µs (47 µs) |
| ont-10k | 277 | 3.6 kbp | 233 µs (149 µs) | 466 µs (309 µs) | 342 µs (231 µs) | 15.2 ms (3.23 ms) | 1.07 ms (582 µs) | 755 µs (283 µs) |
| ont-50k | 104 | 9.52 kbp | 987 µs (303 µs) | 1.67 ms (603 µs) | 1.56 ms (457 µs) | 187 ms (10.4 ms), 81/104 | 5.85 ms (1.73 ms) | 7.43 ms (872 µs) |
| ont-500k | 15 | 638 kbp | 141 ms (105 ms) | 283 ms (135 ms) | 1.14 s (792 ms) | — | 4.71 s (2.88 s), 4/15 | 19.7 s (19.7 s), 1/15 |
| ont-500k-genvar | 15 | 651 kbp | 187 ms (137 ms) | 356 ms (248 ms) | 1.33 s (1.08 s) | — | 4.93 s (4.9 s), 4/15 | 6.34 s (7.53 s), 3/15 |
| sars-cov-2 | 33 | 29.6 kbp | 395 µs (236 µs) | 2.43 ms (2.52 ms) | 1.15 ms (1.1 ms) | 6.41 ms (2.03 ms) | 8.25 ms (7.87 ms) | 897 µs (446 µs) |

- **dinara-align has the lowest mean and the lowest median on every dataset.** On the short reads, the case A\*PA2's post leaves to BiWFA, 38 against BiWFA's 45 µs, median 32 against 47; on ont-10k and ont-50k it takes about two thirds of A\*PA2-simple's time, and on the SARS-CoV-2 genomes 395 against BiWFA's 897 µs.
- **On the 500 kbp reads it takes about half of A\*PA2-full's time:** 141 against 283 ms on ont-500k and 187 against 356 ms on ont-500k-genvar, medians 105 against 135 and 137 against 248.
- **Edlib and BiWFA finish only some of the long reads within 20 seconds, and A\*PA none;** their means cover the reads they finished.

### Divergence, 100 kbp pairs

Mean time per alignment, A\*PA with whichever of `r = 1` and `r = 2` was the faster at each divergence:

| divergence | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: |
| 0% | 149 µs | 7.45 ms | 2.52 ms | 6.38 ms | 41.9 ms | 605 µs |
| 1% | 1.65 ms | 9 ms | 5.89 ms | 7.52 ms | 56.9 ms | 3.22 ms |
| 2% | 5.56 ms | 8.97 ms | 8.82 ms | 8.19 ms | 70.1 ms | 9.97 ms |
| 3% | 4.92 ms | 9.29 ms | 16.2 ms | 7.84 ms | 93.7 ms | 21.1 ms |
| 4% | 4.99 ms | 9.31 ms | 14 ms | 9.54 ms | 96.1 ms | 35.8 ms |
| 5% | 5.39 ms | 9.63 ms | 31.4 ms | 11.9 ms | 140 ms | 54.2 ms |
| 6% | 5.94 ms | 10.4 ms | 28.4 ms | 17.1 ms | 142 ms | 74 ms |
| 7% | 6.52 ms | 10.3 ms | 26.5 ms | 36.1 ms | 145 ms | 97.5 ms |
| 8% | 6.5 ms | 14.7 ms | 24.4 ms | 39.7 ms | 149 ms | 124 ms |
| 9% | 10.5 ms | 16.4 ms | 22.7 ms | 40.5 ms | 153 ms | 153 ms |
| 10% | 10.7 ms | 14.2 ms | 58.6 ms | 43.9 ms | 233 ms | 186 ms |
| 11% | 10.6 ms | 20.6 ms | 56.2 ms | 48.4 ms | 236 ms | 220 ms |
| 12% | 10.7 ms | 17.6 ms | 54.3 ms | 52.2 ms | 238 ms | 255 ms |
| 13% | 10.4 ms | 30.5 ms | 52.2 ms | 58.7 ms | 242 ms | 293 ms |
| 14% | 11.5 ms | 27.9 ms | 50.9 ms | 69.3 ms | 244 ms | 326 ms |
| 15% | 11.7 ms | 25.4 ms | 49.4 ms | 96.9 ms | 248 ms | 374 ms |

- **dinara-align is the fastest at every divergence.**
- **Near-identical pairs finish in the diagonal transition:** 149 µs against BiWFA's 605 at 0%, 1.65 against 3.22 ms at 1%.
- **From 2% to 8% the exact seeds' bound lands close to the distance:** dinara-align runs in 4.9 to 6.5 ms, the next aligner, A\*PA at 2 and 3% and A\*PA2-full beyond, in 7.8 to 14.7.
- **From 9% on, where fewer than a fifth of the exact seeds chain,** dinara-align rebuilds them to match within one edit and stays between 10.4 and 11.7 ms, its narrowest lead 10.7 against A\*PA2-full's 14.2 at 10%, while A\*PA2-full climbs to 25 to 31 ms by 13 to 15%; at 15% A\*PA2-simple takes 49 ms, A\*PA 97, Edlib 248 and BiWFA 374. On the M2, whose band is slower beside its seeds' setup, the rebuild starts at 6% (see the history below).

### Length

Mean time per alignment at 5% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: |
| 3 kbp | 49 µs | 216 µs | 91 µs | 261 µs | 829 µs | 83 µs |
| 10 kbp | 300 µs | 926 µs | 468 µs | 923 µs | 3.85 ms | 557 µs |
| 30 kbp | 1.55 ms | 2.86 ms | 3.15 ms | 3.05 ms | 14.3 ms | 4.98 ms |
| 100 kbp | 5.1 ms | 9.58 ms | 31.4 ms | 12.4 ms | 140 ms | 56.1 ms |
| 300 kbp | 20.6 ms | 32.5 ms | 135 ms | 39.3 ms | 714 ms | 472 ms |
| 1 Mbp | 70.2 ms | 152 ms | 1.83 s | 168 ms | 8.16 s (2/4) | 5.28 s (3/4) |

At 15% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: |
| 3 kbp | 119 µs | 243 µs | 212 µs | 1.52 ms | 883 µs | 377 µs |
| 10 kbp | 421 µs | 1.2 ms | 1.31 ms | 6.29 ms | 5.69 ms | 3.99 ms |
| 30 kbp | 2.33 ms | 6.32 ms | 4.62 ms | 22.7 ms | 23.4 ms | 34.6 ms |
| 100 kbp | 11.8 ms | 26.3 ms | 49.2 ms | 97.5 ms | 248 ms | 375 ms |
| 300 kbp | 50.7 ms | 237 ms | 624 ms | 453 ms | 2.38 s | 3.37 s |
| 1 Mbp | 344 ms | 2.76 s | 2.88 s | 2.01 s | 17.7 s (1/4) | — |

- **dinara-align is the fastest at every length and both divergences.** At 5% it aligns a 1 Mbp pair in 70.2 ms against A\*PA2-full's 152 and A\*PA's 168; at 15% in 344 ms against A\*PA's 2.01 s and A\*PA2-full's 2.76 s.
- **At 15% its band still grows faster than the length,** but less than A\*PA2-full's: from 300 kbp to 1 Mbp, 3.3 times the length, dinara-align takes 6.8 times as long and A\*PA2-full 11.6 times, where A\*PA, whose pruning makes its heuristic nearly exact, takes 4.4 times.

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
| today | 141 ms (105 ms) | 187 ms (137 ms) |

- **On ont-500k the time falls from 1.32 s at the port of A\*PA2-simple's band doubling to 162 ms once the retries are aimed.**
- **On ont-500k-genvar the commits between the port and `b60ece0` were several times slower,** 4.3 to 8.5 s a read and up to eight reads past the budget: real reads gather their errors at the ends, and the projections those commits aimed by ran many times over the distance, so they swept bands far wider than needed. `b60ece0` trusted a projection only when two of them agreed and otherwise grew the band from what was known, bringing genvar back to 752 ms, and aiming the retries to 199 ms.
- **Inexact seeds made these reads slower here at first,** 162 to 181 ms and 199 to 244 ms, where on the M2 they made them faster, 171 to 131 and 262 to 197. This machine's seeds' setup ran twice as slow as the M2's while its band, on AVX-512, runs as fast, so inexact seeds paid only on more divergent pairs; the cutoff now follows the vector width, 20% of the exact seeds chained here against the M2's 40%.
- **A leaner seeds' setup since** (`fc0e5e4`: a branch-free one-edit test, only the seeds whose matches can reach the end, a leaner local pruning) brings both below their best before inexact seeds: 141 ms on ont-500k and 187 ms on ont-500k-genvar.

### Memory

Peak resident memory on the same machine, from `pixi run bench-astarpa2`, on its smaller samples of four pairs a dataset (on Linux each runner runs under GNU time, so that it does not start from the harness's own peak):

- dinara-align starts at 11 MB, the Mojo runtime, where the Rust aligners start at 3 MB, and stays under 20 MB on pairs up to 100 kbp.
- On the 500 kbp reads it peaks at 68 to 71 MB; A\*PA2-full at 89 MB on ont-500k, more, and at 56 MB on ont-500k-genvar, less.
- On 1 Mbp pairs at 15% it peaks at 100 MB, against A\*PA2-full's 145 MB and A\*PA2-simple's 98; A\*PA passes 183 MB before its budget stops it.

### Beyond one thread

All of the above is one thread, as A\*PA2's evaluation measures, and one pair always runs on one thread. dinara-align also aligns batches of pairs across threads, one pair a thread, which the evaluation runs none of the others doing. On this machine's ten cores, from `pixi run bench-astarpa2`, a batch takes 6 µs a read on ont-1k, 34 µs on ont-10k and 74 ms on four ont-500k-genvar reads, against 45 µs, 342 µs and 376 ms for the fastest single-threaded aligner on each.

### On the Apple M2

The same collection on the M2 above, unpinned, its frequency free, 2 to 6% between repeated runs. The order is mostly the same, with one exception: on 100 kbp pairs at 6 and 7% divergence A\*PA2-full is the faster, 5.25 and 5.0 against 6.3 ms, where the M2's cutoff, 40% of the exact seeds chained, already rebuilds the seeds to match within one edit. Elsewhere dinara-align leads, narrowly on ont-1k, 25 against BiWFA's 28 µs, and inexact seeds speed the long reads up rather than down.

<details>
<summary>The M2's tables</summary>

Real datasets, mean time per alignment, median in brackets:

| dataset | pairs | mean length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: |
| ont-1k | 1221 | 0.818 kbp | 25 µs (22 µs) | 64 µs (70 µs) | 40 µs (46 µs) | 389 µs (357 µs) | 73 µs (78 µs) | 28 µs (29 µs) |
| ont-10k | 277 | 3.6 kbp | 151 µs (96 µs) | 293 µs (204 µs) | 209 µs (141 µs) | 8.8 ms (2.02 ms) | 708 µs (372 µs) | 369 µs (151 µs) |
| ont-50k | 104 | 9.52 kbp | 663 µs (203 µs) | 1.01 ms (397 µs) | 890 µs (278 µs) | 121 ms (5.53 ms) | 4.8 ms (1.24 ms) | 3.5 ms (427 µs) |
| ont-500k | 15 | 638 kbp | 114 ms (63.7 ms) | 160 ms (75.9 ms) | 635 ms (443 ms) | — | 4.16 s (2.53 s), 4/15 | 2.03 s (891 ms), 10/15 |
| ont-500k-genvar | 15 | 651 kbp | 169 ms (113 ms) | 199 ms (133 ms) | 745 ms (609 ms) | — | 4.38 s (4.33 s), 4/15 | 2.59 s (2.75 s), 8/15 |
| sars-cov-2 | 33 | 29.6 kbp | 238 µs (147 µs) | 1.38 ms (1.38 ms) | 705 µs (661 µs) | 4.12 ms (1.1 ms) | 6.97 ms (6.78 ms) | 464 µs (216 µs) |

Divergence, 100 kbp pairs:

| divergence | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: |
| 0% | 90 µs | 3.83 ms | 1.5 ms | 4.16 ms | 22.5 ms | 344 µs |
| 1% | 867 µs | 4.57 ms | 3.55 ms | 3.94 ms | 32.3 ms | 1.46 ms |
| 2% | 3.05 ms | 4.6 ms | 5.13 ms | 4.14 ms | 43.6 ms | 4.93 ms |
| 3% | 2.66 ms | 4.71 ms | 9.23 ms | 4.44 ms | 64.8 ms | 9.94 ms |
| 4% | 2.73 ms | 4.87 ms | 9 ms | 5.03 ms | 66.8 ms | 16.2 ms |
| 5% | 3.07 ms | 5 ms | 17.6 ms | 6.46 ms | 103 ms | 24.6 ms |
| 6% | 6.3 ms | 5.25 ms | 15.9 ms | 9.85 ms | 105 ms | 35.2 ms |
| 7% | 6.25 ms | 5 ms | 14.9 ms | 20.2 ms | 108 ms | 49 ms |
| 8% | 6.36 ms | 7.58 ms | 14 ms | 20.4 ms | 112 ms | 59.7 ms |
| 9% | 6.4 ms | 9.71 ms | 12.7 ms | 21.2 ms | 115 ms | 74.7 ms |
| 10% | 6.47 ms | 8.89 ms | 32.8 ms | 23.1 ms | 189 ms | 90.5 ms |
| 11% | 6.55 ms | 12.1 ms | 31.4 ms | 27.2 ms | 189 ms | 108 ms |
| 12% | 6.83 ms | 10 ms | 30.4 ms | 28.5 ms | 192 ms | 128 ms |
| 13% | 6.77 ms | 17 ms | 29.5 ms | 33 ms | 192 ms | 143 ms |
| 14% | 7.89 ms | 16.1 ms | 28.6 ms | 40.4 ms | 199 ms | 175 ms |
| 15% | 8.6 ms | 14.4 ms | 27.6 ms | 53.4 ms | 199 ms | 196 ms |

Length at 5% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: |
| 3 kbp | 32 µs | 150 µs | 66 µs | 178 µs | 346 µs | 51 µs |
| 10 kbp | 206 µs | 582 µs | 306 µs | 546 µs | 1.96 ms | 295 µs |
| 30 kbp | 883 µs | 1.53 ms | 1.85 ms | 1.61 ms | 11.7 ms | 2.3 ms |
| 100 kbp | 2.86 ms | 4.96 ms | 17.6 ms | 6.51 ms | 103 ms | 24.6 ms |
| 300 kbp | 12.2 ms | 17.3 ms | 74.5 ms | 18 ms | 619 ms | 235 ms |
| 1 Mbp | 40.5 ms | 71.4 ms | 1.04 s | 73 ms | 7.26 s (2/4) | 2.45 s |

Length at 15% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA |
| :-- | --: | --: | --: | --: | --: | --: |
| 3 kbp | 82 µs | 170 µs | 135 µs | 999 µs | 411 µs | 192 µs |
| 10 kbp | 362 µs | 804 µs | 760 µs | 3.91 ms | 3.38 ms | 1.87 ms |
| 30 kbp | 2.1 ms | 3.89 ms | 2.63 ms | 14.3 ms | 19.3 ms | 16.4 ms |
| 100 kbp | 8.56 ms | 14.4 ms | 27.3 ms | 52.7 ms | 200 ms | 196 ms |
| 300 kbp | 45.2 ms | 134 ms | 348 ms | 227 ms | 2.1 s | 1.65 s |
| 1 Mbp | 383 ms | 1.56 s | 1.6 s | 949 ms | 15.7 s (1/4) | 17.3 s (1/4) |

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
| today | 114 ms (63.7 ms) | 169 ms (113 ms) |

</details>
