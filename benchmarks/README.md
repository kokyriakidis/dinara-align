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
| dinara-align (bit-parallel, …)                                                  | this tree | MPL-2.0    | the edit-distance workloads: `distance` and `align` at unit costs, ported from A\*PA2-simple, on one thread and on all of them |
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

Measured with `pixi run bench --full`, which took about three minutes, on an Apple M2 (4 performance and 4 efficiency cores, 24 GB) under macOS 27.0.1, with Mojo 1.1.0 and Xcode 27.0, beside other work at a load of about three.
Hyalite built with Rust 1.92, A\*PA with the nightly its repository pins.
A dash marks a task the tool does not offer, or a workload it is not run on.

| workload | task | dinara-align (cpu) | dinara-align (gpu) | dinara-align (bit-parallel, 1 thread) | hyalite | a*pa2-full | a*pa2-simple | a*pa2-nw | a*pa | edlib | biwfa | wfa | agree |
| :-- | :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| reads-150bp | score | 10.8 ms | 15.7 ms | — | 280 ms | — | — | — | — | — | — | — | ✓ |
| reads-150bp | alignment | 22.8 ms | 117 ms | — | 924 ms | — | — | — | — | — | — | — | ✓ |
| reads-1kbp | score | 31.2 ms | 39.6 ms | — | 1.23 s | — | — | — | — | — | — | — | ✓ |
| reads-1kbp | alignment | 41.8 ms | 220 ms | — | 4.8 s | — | — | — | — | — | — | — | ✓ |
| affine-1k | score | 61 µs | 5.24 ms | — | 1.24 ms | — | — | — | — | — | — | — | ✓ |
| affine-1k | alignment | 88 µs | 9.27 ms | — | 3.85 ms | — | — | — | — | — | — | — | ✓ |
| affine-10k | score | 1.81 ms | 18.2 ms | — | 206 ms | — | — | — | — | — | — | — | ✓ |
| affine-10k | alignment | 2.24 ms | 33.2 ms | — | 792 ms | — | — | — | — | — | — | — | ✓ |
| affine-100k | score | 139 ms | 368 ms | — | 19.8 s | — | — | — | — | — | — | — | ✓ |
| affine-100k | alignment | 267 ms | 747 ms | — | 81.2 s | — | — | — | — | — | — | — | ✓ |
| edit-1k-1% | score | — | — | 1 µs | — | 37 µs | 20 µs | 22 µs | — | 13 µs | 1 µs | 1 µs | ✓ |
| edit-1k-1% | alignment | — | — | 2 µs | — | 39 µs | 18 µs | 57 µs | 153 µs | 47 µs | 4 µs | 3 µs | ✓ |
| edit-1k-5% | score | — | — | 3 µs | — | 37 µs | 20 µs | 22 µs | — | 12 µs | 3 µs | 3 µs | ✓ |
| edit-1k-5% | alignment | — | — | 5 µs | — | 40 µs | 18 µs | 58 µs | 144 µs | 49 µs | 8 µs | 6 µs | ✓ |
| edit-1k-15% | score | — | — | 14 µs | — | 31 µs | 20 µs | 22 µs | — | 33 µs | 15 µs | 15 µs | ✓ |
| edit-1k-15% | alignment | — | — | 22 µs | — | 48 µs | 31 µs | 60 µs | 507 µs | 74 µs | 30 µs | 24 µs | ✓ |
| edit-10k-1% | score | — | — | 10 µs | — | 488 µs | 930 µs | 1.69 ms | — | 240 µs | 15 µs | 15 µs | ✓ |
| edit-10k-1% | alignment | — | — | 22 µs | — | 502 µs | 153 µs | 2.89 ms | 1.61 ms | 1.47 ms | 45 µs | 35 µs | ✓ |
| edit-10k-5% | score | — | — | 103 µs | — | 507 µs | 1.63 ms | 1.67 ms | — | 374 µs | 169 µs | 169 µs | ✓ |
| edit-10k-5% | alignment | — | — | 137 µs | — | 544 µs | 264 µs | 2.94 ms | 1.49 ms | 1.85 ms | 374 µs | 287 µs | ✓ |
| edit-10k-15% | score | — | — | 208 µs | — | 876 µs | 2.78 ms | 1.69 ms | — | 1.18 ms | 1.15 ms | 1.16 ms | ✓ |
| edit-10k-15% | alignment | — | — | 299 µs | — | 1.1 ms | 763 µs | 3.13 ms | 39.2 ms | 3.35 ms | 2.44 ms | 1.95 ms | ✓ |
| edit-100k-1% | score | — | — | 598 µs | — | 4.26 ms | 187 ms | 162 ms | — | 5.61 ms | 686 µs | 685 µs | ✓ |
| edit-100k-1% | alignment | — | — | 909 µs | — | 4.39 ms | 3.33 ms | 912 ms | 19.6 ms | 32.1 ms | 1.52 ms | 1.49 ms | ✓ |
| edit-100k-5% | score | — | — | 1.94 ms | — | 4.58 ms | 277 ms | 163 ms | — | 45.3 ms | 16.4 ms | 16.2 ms | ✓ |
| edit-100k-5% | alignment | — | — | 2.4 ms | — | 4.97 ms | 16.2 ms | 826 ms | 17.9 ms | 104 ms | 31 ms | 25.7 ms | ✓ |
| edit-100k-15% | score | — | — | 6.59 ms | — | 21.8 ms | 242 ms | 161 ms | — | 71.2 ms | 120 ms | 120 ms | ✓ |
| edit-100k-15% | alignment | — | — | 7.77 ms | — | 23.7 ms | 25.1 ms | 732 ms | 7.09 s | 205 ms | 252 ms | 221 ms | ✓ |

## Reading the Numbers

- **Every answer agreed**, so each row compares tools answering the same question. The answers are the optimal scores and distances: a sum and a position-weighted sum per workload, so a reordered batch cannot pass.
- **On edit distance, dinara-align on one thread beats or matches every exact aligner on every workload**, distance and alignment alike: 22 against WFA's 35 µs and A\*PA2-simple's 153 µs aligning 10 kbp at 1%, 137 against A\*PA2-simple's 264 µs at 10 kbp and 5%, 2.4 against A\*PA2-full's 5.0 ms at 100 kbp and 5%, 7.8 against A\*PA2-full's 23.7 ms at 100 kbp and 15%, and 598 against WFA's 685 µs scoring 100 kbp at 1%. The ties are the 1 kbp pairs at 1 and 5%, a few microseconds each, where BiWFA and WFA run the same diagonal transition.
- **Near-identical pairs run diagonal transition**, as WFA does, before any band: one front keeping its history for a close alignment, and two from both ends otherwise, as BiWFA scores, keeping both histories for an alignment and tracing it back through each from where they met. It stays on while it costs less than the band would, which at 1 to 3% divergence covers most pairs up to tens of kbp.
- **A band re-aims its first bound as it sweeps.** The bound starts from a projection of the first few edits, which strays by up to half either way. At an eighth, a quarter and half of the columns, the band projects the distance from its own climb, hundreds of edits in, lowering the bound when it was set too high and giving the round up early when it was set too low. Only a distance within the final bound is accepted, so the answer stays exact.
- **Long, moderately divergent pairs prune with A\*PA2-full's seed heuristic.** From about 1,500 projected edits up to one in seven bases, the first sequence is cut into 12-base seeds, their exact matches in the second are found by a rolling hash, and matches that cannot shorten any path over the next 14 seeds are dropped, as A\*PA's local pruning drops them. The band then keeps a row only while its score plus the seeds still ahead, less the longest chain of matches it can still reach, fits the bound. At 2 to 3% divergence that bound at the origin lands within 1% of the distance, so the band starts just past it and finishes in one round: the 100 kbp pair at 5% drops from 5.3 to 3.0 ms. A\*PA2-full also prunes matches between rounds, which dinara-align leaves out, as its rounds now rarely number more than two. On AVX-512, whose band runs fast beside the seeds' setup, seeds are used only from 86 kbp, where they start to pay there.
- **Long, divergent pairs match seeds within one edit**, as A\*PA's `r = 2` does. Past about one edit in fifteen bases most exact 12-base seeds are broken and the bound at the origin falls far short: 8,277 for a distance of 12,196 on a 100 kbp pair at 15%. Seeds of 16 bases that may match with one edit, charging two edits where they match nowhere, put it at 11,304. Each match is found by the half of it that matches exactly, a check on its quarters turns most chance lookups away, and an exact match's one-edit neighbours are left out, which keeps the bound a lower bound. They are used from 64 kbp, where the band they save outgrows their setup, when the projection says one edit in ten bases or under 40% of the exact seeds chain: the 100 kbp alignment at 15% drops from 12.4 to 8.2 ms.
- **Global affine scores run a wavefront from both ends.** Under a table of one match and one mismatch score, as `Scoring.dna()` is, the match reward folds away for a global alignment (Eizenga and Lindquist), and WFA's three-layer wavefront then searches by cost, from the origin and from the corner at once until they meet: 139 ms against hyalite's 19.8 s at 100 kbp, 61 µs against 1.24 ms at 1 kbp. A pair whose projected wavefront would cost more than a full sweep is handed to the sweep, which runs sixteen cells at a time by anti-diagonal and also serves local scores: 39 ms against hyalite's 181 to 198 ms on a 10 kbp pair at 30% divergence or unrelated.
- **Global affine alignments are traced back through the same wavefront's fronts**, five bytes a diagonal kept from each end and the path walked from where they met, a pair whose fronts would pass 80 MB split where an optimal path crosses: 88 µs against hyalite's 3.85 ms at 1 kbp, 2.24 ms against 792 ms at 10 kbp and 267 ms against 81.2 s at 100 kbp, ahead of dinara-align's own GPU sweep, 33.2 and 747 ms. A tie between optimal paths may resolve differently from the GPU's, the score the same. Local alignments, and tables of more than one mismatch score, sweep the whole matrix sixteen cells at a time by anti-diagonal, storing the three layers up to six million cells and splitting on rows in linear space beyond, Myers and Miller's way.
- **The other rivals run on one CPU thread.** dinara-align's CPU batches spread their pairs over every thread, one pair a thread, which the reads rows show; the single affine pairs, and every edit-distance pair, run on one.
- **dinara-align's GPU times include the overheads a caller pays**: opening the device context on every call, copying sequences over, and copying results back. That is why a single short pair is slower on the GPU than on the CPU.
- **hyalite's traceback budget is 1 GiB.** It keeps the whole matrix when it fits and switches to a checkpointed sweep above that, as the 100 kbp pairs do. dinara-align switches to its linear-space recursion above six million cells on the host and one million on the device.
- **A\*PA answers only edit distance**, and its time depends on how similar the sequences are, as do the bit-parallel path's and dinara-align's global affine score and alignment; the local alignments and hyalite's columns do not.
- **One machine, one fixed pair per workload.** On fresh pairs from two other seeds, from 1 to 100 kbp at 2 to 30% divergence, dinara-align on one thread beat or matched every exact aligner here, A\*PA2, A\*PA, Edlib, BiWFA and WFA, on every pair, distance and alignment alike. Rerun on your own hardware before quoting any of this.

## A\*PA2's Evaluation

`pa_bench.py` runs the exact aligners of A\*PA2's evaluation, with its parameters, and WFA beside BiWFA, on its own datasets: Oxford Nanopore reads and SARS-CoV-2 genomes from pa-bench's release, and the uniform-error pairs pa-generate regenerates exactly from the evaluation's seed and lock.
Each dataset is a fixed shuffled sample of about 2 Mbp and at least four pairs, the same for every tool; each tool aligns its pairs with traceback, once each, as pa-bench times them, within a budget of five seconds, and every cost is checked against every other tool's and against the costs A\*PA2's published results recorded.
The rivals are pinned, so their results are kept and reused until `--fresh`; a second table gives each tool's peak resident memory.
dinara-align also runs as a batch on all eight threads, `alignments` over the whole sample, one pair a thread, whose column is the batch's wall-clock time over its pairs, a throughput where the others are each one pair's latency.

Measured as above, dinara-align's columns the fastest of three warm runs; a count marks a tool its budget, or WFA's 8 GiB cap, stopped partway, and more than the budget one that finished no pair:

| dataset | pairs | mean length | dinara-align (bit-parallel, 1 thread) | dinara-align (bit-parallel, batch, 8 threads) | a*pa2-full | a*pa2-simple | a*pa | edlib | biwfa | wfa | agree |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| ont-1k | 1221 of 12477 | 0.818 kbp | 23 µs | 8 µs | 64 µs | 40 µs | 389 µs | 70 µs | 27 µs | 21 µs | ✓ |
| ont-10k | 277 of 5000 | 3.6 kbp | 149 µs | 43 µs | 293 µs | 209 µs | 8.8 ms | 715 µs | 369 µs | 316 µs | ✓ |
| ont-50k | 104 of 10000 | 9.52 kbp | 685 µs | 191 µs | 1.01 ms | 890 µs | 118 ms (43/104) | 4.64 ms | 3.51 ms | 3.06 ms | ✓ |
| ont-500k | 4 of 50 | 638 kbp | 142 ms | 98.4 ms | 280 ms | 805 ms | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| ont-500k-genvar | 4 of 48 | 659 kbp | 184 ms | 83.4 ms | 269 ms | 627 ms | > 5 s | 4.35 s (1/4) | 4.23 s (1/4) | > 5 s | ✓ |
| sars-cov-2 | 33 of 10000 | 29.6 kbp | 234 µs | 109 µs | 1.38 ms | 705 µs | 4.12 ms | 7.34 ms | 481 µs | 478 µs | ✓ |
| Uniform-t10000000-n3000-e0.05 | 333 of 3333 | 3 kbp | 37 µs | 10 µs | 150 µs | 66 µs | 178 µs | 349 µs | 49 µs | 41 µs | ✓ |
| Uniform-t10000000-n10000-e0.05 | 99 of 1000 | 10 kbp | 196 µs | 58 µs | 582 µs | 306 µs | 546 µs | 1.95 ms | 291 µs | 250 µs | ✓ |
| Uniform-t10000000-n30000-e0.05 | 33 of 333 | 30 kbp | 873 µs | 266 µs | 1.53 ms | 1.85 ms | 1.61 ms | 12.3 ms | 2.31 ms | 2.04 ms | ✓ |
| Uniform-t10000000-n100000-e0.05 | 10 of 100 | 100 kbp | 2.86 ms | 957 µs | 4.96 ms | 17.6 ms | 6.51 ms | 104 ms | 24.7 ms | 21.7 ms | ✓ |
| Uniform-t10000000-n300000-e0.05 | 4 of 33 | 300 kbp | 12.1 ms | 3.57 ms | 17.3 ms | 74.5 ms | 18 ms | 617 ms | 234 ms | 203 ms | ✓ |
| Uniform-t10000000-n1000000-e0.05 | 4 of 10 | 1e+03 kbp | 40.3 ms | 11.7 ms | 71.4 ms | 1.04 s | 73 ms | > 5 s | 2.44 s (2/4) | 2.37 s (2/4) | ✓ |
| Uniform-t10000000-n3000-e0.15 | 333 of 3333 | 3 kbp | 79 µs | 24 µs | 170 µs | 135 µs | 999 µs | 411 µs | 189 µs | 165 µs | ✓ |
| Uniform-t10000000-n10000-e0.15 | 99 of 1000 | 10 kbp | 348 µs | 103 µs | 804 µs | 760 µs | 3.91 ms | 3.36 ms | 1.87 ms | 1.57 ms | ✓ |
| Uniform-t10000000-n30000-e0.15 | 33 of 333 | 30 kbp | 2.06 ms | 591 µs | 3.89 ms | 2.63 ms | 14.3 ms | 19.4 ms | 16 ms | 13.1 ms | ✓ |
| Uniform-t10000000-n100000-e0.15 | 9 of 100 | 100 kbp | 8.51 ms | 3.04 ms | 14.4 ms | 27.3 ms | 52.7 ms | 200 ms | 197 ms | 145 ms | ✓ |
| Uniform-t10000000-n300000-e0.15 | 4 of 33 | 300 kbp | 46 ms | 12.6 ms | 134 ms | 348 ms | 227 ms | 2.11 s (2/4) | 1.67 s (3/4) | 1.34 s | ✓ |
| Uniform-t10000000-n1000000-e0.15 | 4 of 10 | 1e+03 kbp | 381 ms | 104 ms | 1.56 s (3/4) | 1.61 s (3/4) | 949 ms | > 5 s | > 5 s | > 5 s | ✓ |

Four long reads are few: on all 50 ont-500k reads dinara-align took 6.5 s, and A\*PA2-full 8.4 s (the faster on 16 of them); on all 48 ont-500k-genvar reads 10.4 s against A\*PA2-full's 11.5 s (the faster on 21 of them).

## Affine Costs

`pa_bench.py --affine x,o,e` runs the same samples at gap-affine costs as WFA counts them, a mismatch `x` and a gap of `k` letters `o + k e`, against the exact aligners that take them: WFA2-lib keeping every front (WFA), its lowest-memory mode (BiWFA), and KSW2's banded SSE kernel with band doubling, on x86-64 alone.
dinara-align answers with `align` under `Costs.affine`, a CIGAR as the others hand back, from a wavefront grown from both ends at once, which keeps five bytes a diagonal for its traceback and splits a pair whose fronts would pass 80 MB where an optimal path crosses, as BiWFA does.

```bash
pixi run bench-astarpa2 --affine 4,6,2
```

At WFA's costs (4, 6, 2) on the i9-7900X above, pinned to one core, each tool five seconds a sample; a count marks a tool its budget stopped partway, and more than the budget one that finished no pair:

| dataset | pairs | mean length | dinara-align | WFA | BiWFA | KSW2 | agree |
| :-- | --: | --: | --: | --: | --: | --: | :-: |
| ont-1k | 1221 of 12477 | 0.818 kbp | 111 µs | 187 µs | 305 µs | 1.11 ms | ✓ |
| ont-10k | 277 of 5000 | 3.6 kbp | 1.8 ms | 4.81 ms | 5.85 ms | 30.8 ms (163/277) | ✓ |
| ont-50k | 104 of 10000 | 9.52 kbp | 24.8 ms | 69.1 ms (82/104) | 57.6 ms (90/104) | 363 ms (16/104) | ✓ |
| ont-500k | 4 of 50 | 638 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| ont-500k-genvar | 4 of 48 | 659 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| sars-cov-2 | 33 of 10000 | 29.6 kbp | 785 µs | 2.06 ms | 1.77 ms | 88.8 ms | ✓ |
| Uniform-t10000000-n3000-e0.05 | 333 of 3333 | 3 kbp | 298 µs | 580 µs | 916 µs | 6.11 ms | ✓ |
| Uniform-t10000000-n10000-e0.05 | 99 of 1000 | 10 kbp | 2.51 ms | 6.35 ms | 8.02 ms | 97.2 ms (52/99) | ✓ |
| Uniform-t10000000-n30000-e0.05 | 33 of 333 | 30 kbp | 20.9 ms | 56 ms | 62.6 ms | 1.35 s (4/33) | ✓ |
| Uniform-t10000000-n100000-e0.05 | 10 of 100 | 100 kbp | 346 ms | 760 ms (7/10) | 642 ms (8/10) | > 5 s | ✓ |
| Uniform-t10000000-n300000-e0.05 | 4 of 33 | 300 kbp | 3.5 s (1/4) | > 5 s | > 5 s | > 5 s | ✓ |
| Uniform-t10000000-n1000000-e0.05 | 4 of 10 | 1e+03 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| Uniform-t10000000-n3000-e0.15 | 333 of 3333 | 3 kbp | 1.42 ms | 3.75 ms | 5.07 ms | 18.7 ms (267/333) | ✓ |
| Uniform-t10000000-n10000-e0.15 | 99 of 1000 | 10 kbp | 13.5 ms | 38.7 ms | 44.8 ms | 203 ms (25/99) | ✓ |
| Uniform-t10000000-n30000-e0.15 | 33 of 333 | 30 kbp | 174 ms (29/33) | 386 ms (13/33) | 372 ms (14/33) | 2.48 s (2/33) | ✓ |
| Uniform-t10000000-n100000-e0.15 | 9 of 100 | 100 kbp | 2.35 s (2/9) | > 5 s | 4.17 s (1/9) | > 5 s | ✓ |
| Uniform-t10000000-n300000-e0.15 | 4 of 33 | 300 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| Uniform-t10000000-n1000000-e0.15 | 4 of 10 | 1e+03 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |

- **dinara-align is the fastest on every dataset any tool finished,** 1.7 to 2.9 times faster than WFA: 111 against 187 µs on ont-1k, 1.8 against 4.81 ms on ont-10k, 785 µs against 2.06 ms on the SARS-CoV-2 genomes, and 13.5 against 38.7 ms on 10 kbp at 15%. BiWFA beats WFA only on the genomes and where WFA's budget runs out.
- **It alone finishes every ont-50k read,** in 24.8 ms on average, where WFA finishes 82 of the 104 and BiWFA 90 within the budget, and it finishes more of the 100 and 300 kbp pairs than anyone.
- **KSW2's band doubles from eight diagonals against a loose bound,** so it times out on most reads past 1 kbp.
- **No aligner finishes a 500 kbp read or a 1 Mbp pair in five seconds:** the wavefront's work grows with the square of the cost, and gap-affine costs have no seed heuristic here to prune it (see TODO.md).

Growth of peak memory aligning each pair, median / largest, as A\*PA2's Table 10 measures it:

| dataset | dinara-align | WFA | BiWFA | KSW2 |
| :-- | --: | --: | --: | --: |
| ont-1k | 0.0 / 0.5 MB | 0.0 / 3.0 MB | 0.0 / 1.5 MB | 0.0 / 0.9 MB |
| ont-10k | 0.0 / 21 MB | 0.0 / 97 MB | 0.0 / 2.1 MB | 0.0 / 58 MB |
| ont-50k | 0.0 / 84 MB | 0.0 / 1420 MB | 0.0 / 4.4 MB | 0.0 / 317 MB |
| sars-cov-2 | 0.0 / 11 MB | 0.0 / 51 MB | 0.0 / 2.5 MB | 0.0 / 211 MB |
| Uniform-t10000000-n3000-e0.05 | 0.0 / 1.2 MB | 0.0 / 5.2 MB | 0.0 / 1.8 MB | 0.0 / 3.1 MB |
| Uniform-t10000000-n10000-e0.05 | 0.0 / 6.4 MB | 0.0 / 29 MB | 0.0 / 1.8 MB | 0.0 / 40 MB |
| Uniform-t10000000-n30000-e0.05 | 0.0 / 45 MB | 0.0 / 212 MB | 0.0 / 3.2 MB | 31 / 455 MB |
| Uniform-t10000000-n100000-e0.05 | 0.0 / 85 MB | 12 / 2413 MB | 0.0 / 9.1 MB | — |
| Uniform-t10000000-n300000-e0.05 | 99 / 99 MB | — | — | — |
| Uniform-t10000000-n3000-e0.15 | 0.0 / 4.0 MB | 0.0 / 18 MB | 0.0 / 1.1 MB | 0.0 / 12 MB |
| Uniform-t10000000-n10000-e0.15 | 0.0 / 34 MB | 0.0 / 160 MB | 0.0 / 3.1 MB | 0.0 / 79 MB |
| Uniform-t10000000-n30000-e0.15 | 0.0 / 104 MB | 0.5 / 1416 MB | 0.0 / 5.5 MB | 795 / 795 MB |
| Uniform-t10000000-n100000-e0.15 | 84 / 84 MB | — | 16 / 16 MB | — |

dinara-align keeps under a tenth of WFA's memory on the largest pairs (85 against 2413 MB at 100 kbp), as its traceback keeps a byte of flags and the column of one front a diagonal, from fronts each half as long; BiWFA, keeping only its last few fronts, stays smallest.

## A\*PA2's Results, Redone

A\*PA2's own results section ([curiouscoding.nl/posts/astarpa2](https://curiouscoding.nl/posts/astarpa2/#results)), redone with dinara-align beside the exact aligners it compares, and WFA, on a machine set up as A\*PA2's was: its real datasets, its uniform pairs swept over divergence and over length, and in place of its ablation, dinara-align's history on the long reads.

```bash
pixi run results-astarpa2 collect   # forty minutes the first time; after that only dinara-align runs, a few minutes
pixi run results-astarpa2 tables    # the tables below, into .cache/results/astarpa2-results.md
```

### Setup

As A\*PA2's evaluation runs them: one single-threaded job at a time, every pair aligned once with its traceback, which every aligner hands back as a CIGAR, the time the average wall clock per alignment, reading the data left out.

The machine is an Intel Core i9-7900X under Ubuntu 26.04 (Linux 7.0), with Mojo 1.1.0, set up as A\*PA2's i7-10750H was: every core fixed at 3.3 GHz, turbo boost and hyper-threading off, and each collection pinned to one core with `taskset`.
One thing differs: the jobs ran at normal priority, not A\*PA2's niceness −20, on an otherwise idle machine.

Each dataset is a fixed shuffled sample, the same pairs for every aligner: about 2 Mbp of sequence and at least four pairs, and 20 Mbp, fifteen reads, of the two 500 kbp read sets.
Each aligner may spend 20 seconds on a sample; one that has not finished stops there, its numbers covering the pairs it finished, counted beside them, and one that finished none shows a dash.
WFA, which A\*PA2's evaluation leaves out, is WFA2-lib keeping every front for its traceback, BiWFA's faster but hungrier twin: its memory grows with the square of the distance, so it also stops at 8 GiB, as its budget would stop it.
On the real datasets, as in A\*PA2's results, two approximate aligners run beside the exact ones, with the evaluation's parameters: WFA-adaptive, WFA2-lib's lowest-memory mode dropping lagging diagonals at its defaults (10, 50, 10), and Block Aligner, with blocks from 0.1 to 1% of the input and, as it takes only affine costs, a gap's opening one more. As pa-bench does, each is held against an exact alignment at its own costs, the exact distance for WFA-adaptive and for Block Aligner BiWFA at the same affine costs, and the share it aligned optimally goes beside its time; a share over fewer pairs than it finished counts only those its exact reference finished within the budget.
Every cost is checked against every other aligner's on the pairs both aligned, and against the costs A\*PA2's published results recorded; a disagreement fails the run, and none occurred.
The rivals, and dinara-align's older commits, never change, so each is run on a sample once and kept; later runs time today's dinara-align alone.
The same tables measured on the Apple M2 above, without any of this pinning, close the section.

### Real datasets

Mean time per alignment, median in brackets:

| dataset | pairs | mean length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA | WFA-adaptive* | Block Aligner* |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| ont-1k | 1221 | 0.818 kbp | 22 µs (22 µs) | 91 µs (101 µs) | 57 µs (65 µs) | 602 µs (547 µs) | 137 µs (135 µs) | 44 µs (46 µs) | 31 µs (32 µs) | 41 µs (42 µs), 93% optimal | 43 µs (49 µs), 85% optimal |
| ont-10k | 277 | 3.6 kbp | 168 µs (115 µs) | 463 µs (306 µs) | 340 µs (232 µs) | 15 ms (3.28 ms) | 1.08 ms (582 µs) | 752 µs (282 µs) | 602 µs (224 µs) | 386 µs (173 µs), 61% optimal | 208 µs (151 µs), 61% optimal |
| ont-50k | 104 | 9.52 kbp | 600 µs (244 µs) | 1.65 ms (599 µs) | 1.55 ms (450 µs) | 197 ms (10.7 ms), 81/104 | 5.86 ms (1.72 ms) | 7.42 ms (873 µs) | 6.43 ms (753 µs) | 2.33 ms (459 µs), 52% optimal | 676 µs (327 µs), 61% optimal |
| ont-500k | 15 | 638 kbp | 99.5 ms (81.7 ms) | 286 ms (135 ms) | 1.14 s (791 ms) | — | 4.71 s (2.88 s), 4/15 | 19.7 s (19.7 s), 1/15 | — | 973 ms (618 ms), 53% optimal | 736 ms (724 ms) |
| ont-500k-genvar | 15 | 651 kbp | 135 ms (105 ms) | 364 ms (256 ms) | 1.33 s (1.08 s) | — | 4.94 s (4.9 s), 4/15 | 6.33 s (7.53 s), 3/15 | — | 533 ms (273 ms), 7% optimal | 875 ms (679 ms), 0% optimal of 1 |
| sars-cov-2 | 33 | 29.6 kbp | 294 µs (146 µs) | 2.39 ms (2.47 ms) | 1.14 ms (1.09 ms) | 6.31 ms (1.99 ms) | 8.24 ms (7.85 ms) | 906 µs (450 µs) | 1.01 ms (368 µs) | 614 µs (355 µs), 97% optimal | 2.63 ms (2.39 ms), 30% optimal |

- **Among the exact aligners, dinara-align has the lowest mean and median on every dataset.** On the short reads, the case A\*PA2's post leaves to BiWFA, 22 against WFA's 31 µs and BiWFA's 44, median 22 against 32, since a short pair no longer hands a band, which covers its whole matrix there, what the diagonal transition does faster; on ont-10k and ont-50k it takes half and two fifths of A\*PA2-simple's time, and on the SARS-CoV-2 genomes 294 µs against WFA's 1.01 ms.
- **On the 500 kbp reads it takes less than half of A\*PA2-full's time:** 99.5 against 286 ms on ont-500k and 135 against 364 ms on ont-500k-genvar, medians 81.7 against 135 and 105 against 256.
- **Edlib finishes only some of the long reads within 20 seconds, BiWFA one ont-500k read and three genvar reads, and A\*PA and WFA none;** their means cover the reads they finished.
- **The approximate aligners trade optimality for speed, and are still slower than dinara-align on every dataset.** The closest, Block Aligner, takes 208 and 676 µs on ont-10k and ont-50k against dinara-align's 168 and 600, aligning 61% of those reads optimally; WFA-adaptive aligns 7% of the genvar reads optimally. A\*PA2 reports the same shares: Block Aligner 85% on ont-1k, WFA-adaptive 93%.

### Divergence, 100 kbp pairs

Mean time per alignment, A\*PA with whichever of `r = 1` and `r = 2` was the faster at each divergence:

| divergence | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| 0% | 84 µs | 7.74 ms | 2.55 ms | 6.42 ms | 42.2 ms | 633 µs | 330 µs |
| 1% | 1.38 ms | 9.06 ms | 5.91 ms | 7.61 ms | 57.2 ms | 3.25 ms | 3.46 ms |
| 2% | 4.49 ms | 9.23 ms | 8.83 ms | 8.23 ms | 70.5 ms | 10.1 ms | 9.77 ms |
| 3% | 4.55 ms | 9.25 ms | 16.3 ms | 8.02 ms | 93.8 ms | 21.1 ms | 19.3 ms |
| 4% | 4.64 ms | 9.52 ms | 14 ms | 9.55 ms | 96.2 ms | 35.8 ms | 31.9 ms |
| 5% | 4.8 ms | 9.66 ms | 31.5 ms | 11.9 ms | 140 ms | 54.2 ms | 47.1 ms |
| 6% | 6.34 ms | 10.3 ms | 28.5 ms | 17.9 ms | 142 ms | 73.9 ms | 64.7 ms |
| 7% | 6.4 ms | 10.4 ms | 26.6 ms | 35.9 ms | 145 ms | 97.2 ms | 84.5 ms |
| 8% | 6.58 ms | 14.8 ms | 24.4 ms | 41.1 ms | 149 ms | 124 ms | 108 ms |
| 9% | 9.58 ms | 16.3 ms | 22.7 ms | 42.6 ms | 153 ms | 153 ms | 130 ms |
| 10% | 10.3 ms | 14.2 ms | 58.6 ms | 45.6 ms | 232 ms | 182 ms | 154 ms |
| 11% | 10.8 ms | 20.9 ms | 56.4 ms | 51.9 ms | 236 ms | 220 ms | 185 ms |
| 12% | 10.6 ms | 17.9 ms | 54.5 ms | 56.2 ms | 238 ms | 254 ms | 214 ms |
| 13% | 10.5 ms | 30.6 ms | 52.2 ms | 65.2 ms | 242 ms | 292 ms | 245 ms |
| 14% | 11.9 ms | 28.2 ms | 51.2 ms | 77.3 ms | 245 ms | 331 ms | 274 ms |
| 15% | 11.4 ms | 25.7 ms | 49.5 ms | 109 ms | 248 ms | 373 ms | 307 ms |

- **dinara-align is the fastest at every divergence.**
- **Near-identical pairs finish in the diagonal transition:** 84 µs against WFA's 330 at 0%, 1.38 against BiWFA's 3.25 ms at 1%.
- **From 2% to 8% the exact seeds' bound lands close to the distance:** dinara-align runs in 4.5 to 6.6 ms, the next aligner, A\*PA at 2 and 3% and A\*PA2-full beyond, in 8.0 to 14.8.
- **From 9% on, where fewer than a fifth of the exact seeds chain,** dinara-align rebuilds them to match within one edit and stays between 9.6 and 11.9 ms, its narrowest lead 10.3 against A\*PA2-full's 14.2 at 10%, while A\*PA2-full climbs to 26 to 31 ms by 13 to 15%; at 15% A\*PA2-simple takes 49 ms, A\*PA 109, Edlib 248, WFA 307 and BiWFA 373. On the M2, whose band is slower beside its seeds' setup, real reads rebuild earlier, but pairs like these, their errors spread, rebuild at 9% there too.

### Length

Mean time per alignment at 5% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| 3 kbp | 43 µs | 212 µs | 90 µs | 255 µs | 831 µs | 82 µs | 58 µs |
| 10 kbp | 295 µs | 913 µs | 463 µs | 905 µs | 3.86 ms | 558 µs | 441 µs |
| 30 kbp | 953 µs | 2.82 ms | 3.14 ms | 3.01 ms | 14.2 ms | 4.99 ms | 4.47 ms |
| 100 kbp | 4.82 ms | 9.45 ms | 31.2 ms | 11.6 ms | 139 ms | 54.4 ms | 46.7 ms |
| 300 kbp | 19.6 ms | 30.9 ms | 134 ms | 35.6 ms | 714 ms | 471 ms | 457 ms |
| 1 Mbp | 61.4 ms | 144 ms | 1.85 s | 154 ms | 8.17 s (2/4) | 5.38 s (3/4) | 5.59 s (3/4) |

At 15% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| 3 kbp | 118 µs | 241 µs | 211 µs | 1.52 ms | 886 µs | 375 µs | 293 µs |
| 10 kbp | 415 µs | 1.18 ms | 1.3 ms | 6.35 ms | 5.7 ms | 3.99 ms | 3.06 ms |
| 30 kbp | 1.73 ms | 6.28 ms | 4.59 ms | 22.7 ms | 23.3 ms | 34.7 ms | 27.1 ms |
| 100 kbp | 11.1 ms | 26.1 ms | 49 ms | 94.4 ms | 249 ms | 375 ms | 318 ms |
| 300 kbp | 44.7 ms | 236 ms | 624 ms | 447 ms | 2.38 s | 3.3 s | 3.05 s |
| 1 Mbp | 261 ms | 2.76 s | 2.89 s | 2.38 s | 17.7 s (1/4) | — | — |

- **dinara-align is the fastest at every length and both divergences.** At 5% it aligns a 1 Mbp pair in 61.4 ms against A\*PA2-full's 144 and A\*PA's 154; at 15% in 261 ms against A\*PA's 2.38 s and A\*PA2-full's 2.76 s.
- **At 15% its band still grows faster than the length,** but less than A\*PA2-full's: from 300 kbp to 1 Mbp, 3.3 times the length, dinara-align takes 5.8 times as long and A\*PA2-full 11.7 times, where A\*PA, whose pruning makes its heuristic nearly exact, takes 5.3 times.

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
| today | 99.5 ms (81.7 ms) | 135 ms (105 ms) |

- **On ont-500k the time falls from 1.32 s at the port of A\*PA2-simple's band doubling to 162 ms once the retries are aimed.**
- **On ont-500k-genvar the commits between the port and `b60ece0` were several times slower,** 4.3 to 8.5 s a read and up to eight reads past the budget: real reads gather their errors at the ends, and the projections those commits aimed by ran many times over the distance, so they swept bands far wider than needed. `b60ece0` trusted a projection only when two of them agreed and otherwise grew the band from what was known, bringing genvar back to 752 ms, and aiming the retries to 199 ms.
- **Inexact seeds made these reads slower here at first,** 162 to 181 ms and 199 to 244 ms, where on the M2 they made them faster, 171 to 131 and 262 to 197. This machine's seeds' setup ran twice as slow as the M2's while its band, on AVX-512, runs as fast, so inexact seeds paid only on more divergent pairs; the cutoff now follows the vector width, 20% of the exact seeds chained here against the M2's 40%.
- **A leaner seeds' setup since** (`fc0e5e4`: a branch-free one-edit test, only the seeds whose matches can reach the end, a leaner local pruning) brings both below their best before inexact seeds, and with retries aimed closer, AVX-512 sweeping two groups at a time (`697fef1`) and exact seeds' matching filtered (`bcdc7dc`) since, they take 99.5 ms on ont-500k and 135 ms on ont-500k-genvar, a CIGAR included (`3152d27`); the diagonal transition's slides by AVX-512 gathers (`bad5ba5`) took the short reads and low-divergence pairs down 7 to 19% since. Every CIGAR now follows a fixed rule for ties (`ea094ba`), the band's retraced a tile at a time by a forward search (`8f9c92f`).

### Memory

From `pixi run bench-astarpa2` on the same machine, over its samples of each dataset (the 1221 reads of ont-1k, four pairs of each long set). Each cell is the most memory aligning one pair added to the most the process had held before it, the measure of A\*PA2's Table 10, with the process's whole peak resident memory in brackets; the lowest of the aligners that finished every pair is bold, and † marks one that finished only some within its budget, a dash none:

| dataset | dinara-align | A\*PA2-full | A\*PA2-simple | A\*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| ***Real datasets*** |  |  |  |  |  |  |  |
| ont-1k | 0.4 MB (13 MB) | 0.3 MB (5 MB) | 0.3 MB (5 MB) | 0.2 MB (6 MB) | 0.5 MB (7 MB) | **0.1 MB** (7 MB) | 0.5 MB (8 MB) |
| ont-10k | 1.0 MB (15 MB) | 0.7 MB (6 MB) | **0.4 MB** (6 MB) | 20 MB (39 MB) | 0.5 MB (8 MB) | 0.9 MB (8 MB) | 6.1 MB (38 MB) |
| ont-50k | 1.9 MB (18 MB) | 1.2 MB (8 MB) | **0.6 MB** (7 MB) | 83 MB (144 MB) † | 0.8 MB (8 MB) | 1.6 MB (9 MB) | 71 MB (187 MB) |
| ont-500k | **48 MB** (64 MB) | 81 MB (89 MB) | 82 MB (90 MB) † | — | — | — | — |
| ont-500k-genvar | 41 MB (68 MB) | 49 MB (57 MB) | **32 MB** (55 MB) | — | 4.6 MB (?) † | — | — |
| sars-cov-2 | 1.2 MB (17 MB) | 1.9 MB (7 MB) | **0.7 MB** (6 MB) | 13 MB (19 MB) | 1.2 MB (8 MB) | 1.5 MB (9 MB) | 15 MB (26 MB) |
| ***Uniform pairs, 5% divergence*** |  |  |  |  |  |  |  |
| 3 kbp, 5% | 0.4 MB (13 MB) | 0.4 MB (5 MB) | 0.3 MB (5 MB) | **0.1 MB** (5 MB) | 0.8 MB (7 MB) | 0.4 MB (7 MB) | 0.6 MB (7 MB) |
| 10 kbp, 5% | 1.5 MB (14 MB) | 0.6 MB (6 MB) | **0.4 MB** (5 MB) | 0.5 MB (5 MB) | 0.7 MB (7 MB) | 1.3 MB (8 MB) | 2.0 MB (9 MB) |
| 30 kbp, 5% | 0.9 MB (14 MB) | 1.8 MB (7 MB) | 0.8 MB (6 MB) | 1.2 MB (6 MB) | **0.6 MB** (7 MB) | 1.4 MB (8 MB) | 7.8 MB (42 MB) |
| 100 kbp, 5% | 3.3 MB (16 MB) | 5.5 MB (11 MB) | 2.5 MB (8 MB) | 4.3 MB (9 MB) | **1.4 MB** (8 MB) | 3.0 MB (10 MB) | 74 MB (100 MB) |
| 300 kbp, 5% | 11 MB (26 MB) | 16 MB (21 MB) | 8.1 MB (13 MB) | 13 MB (20 MB) | **2.4 MB** (9 MB) | 7.8 MB (17 MB) | 648 MB (673 MB) |
| 1 Mbp, 5% | **30 MB** (56 MB) | 53 MB (65 MB) | 58 MB (69 MB) † | 42 MB (63 MB) | — | — | — |
| ***Uniform pairs, 15% divergence*** |  |  |  |  |  |  |  |
| 3 kbp, 15% | 0.4 MB (13 MB) | 0.4 MB (5 MB) | **0.3 MB** (5 MB) | 0.6 MB (6 MB) | 0.8 MB (7 MB) | 0.9 MB (7 MB) | 1.5 MB (8 MB) |
| 10 kbp, 15% | 0.5 MB (13 MB) | 0.8 MB (6 MB) | **0.4 MB** (6 MB) | 2.2 MB (7 MB) | 0.7 MB (7 MB) | 0.9 MB (7 MB) | 6.4 MB (38 MB) |
| 30 kbp, 15% | 1.1 MB (14 MB) | 1.8 MB (7 MB) | **0.9 MB** (6 MB) | 4.6 MB (10 MB) | **0.9 MB** (7 MB) | 2.1 MB (9 MB) | 48 MB (78 MB) |
| 100 kbp, 15% | 3.6 MB (17 MB) | 5.8 MB (10 MB) | 3.0 MB (8 MB) | 18 MB (24 MB) | **1.4 MB** (8 MB) | 3.8 MB (11 MB) | 511 MB (551 MB) |
| 300 kbp, 15% | **13 MB** (29 MB) | 21 MB (27 MB) | 19 MB (25 MB) | 38 MB (50 MB) | 2.2 MB (?) † | 8.9 MB (?) † | 4.5 GB (?) † |
| 1 Mbp, 15% | **52 MB** (81 MB) | 133 MB (145 MB) † | 86 MB (97 MB) † | 149 MB (?) † | — | — | — |
| *runtime alone (one 8 bp pair)* | *10 MB* | *3 MB* | *3 MB* | *3 MB* | *5 MB* | *5 MB* | *5 MB* |

- **The memory an alignment adds is what compares aligners**: each runner reads `getrusage` before and after each pair, outside its timer, so the runtime and the sample read in count for nothing. The whole peak counts them, 13 to 15 MB for dinara-align's runner on short pairs, most of it the Mojo runtime's 10 MB (a Mojo program printing one line peaks at 9.3 MB), against 5 to 8 MB for the Rust ones, so there the brackets show runtimes, not aligners. Each runner reads its sample in one call of the file's size; read through a growing buffer instead, the memory it left behind went to the alignments unseen, and dinara-align's growth on ont-500k measured 32 MB rather than 48.
- **On pairs up to 30 kbp every aligner adds about 2 MB at most**, but WFA, keeping every front, and A\*PA; at 100 kbp Edlib adds the least, 1.4 MB, dinara-align 3.3 to 3.6 and A\*PA2 2.5 to 5.8.
- **On the long reads dinara-align adds the least but on ont-500k-genvar**: 48 MB on ont-500k against A\*PA2's 81 and 82, and on 1 Mbp pairs 30 MB at 5% against 42 to 58 and 52 MB at 15% against 86 to 149, A\*PA2's on the two pairs they finished; on ont-500k-genvar 41 MB, between A\*PA2-simple's 32 and A\*PA2-full's 49. A\*PA2 reports the same shape on whole datasets: on reads over 500 kbp A\*PA2-full adds 30 MB in median and 82 at most.
- **WFA's grows with the square of the distance**, 511 MB at 100 kbp and 4.5 GB at 300 kbp at 15%, where its memory cap stops it.

### Beyond one thread

All of the above is one thread, as A\*PA2's evaluation measures, and one pair always runs on one thread. dinara-align also aligns batches of pairs across threads, one pair a thread, which the evaluation runs none of the others doing. On this machine's ten cores, from `pixi run bench-astarpa2`, a batch takes 7 µs a read on ont-1k, 37 µs on ont-10k and 76 ms on four ont-500k-genvar reads, against 31 µs, 342 µs and 376 ms for the fastest single-threaded aligner on each.

### On the Apple M2

The same collection on the M2 above, unpinned, its frequency free, 2 to 6% between repeated runs, dinara-align's alignment still written as two gapped rows, within a few percent of its CIGAR either way. The order is the same: dinara-align is the fastest exact aligner on every row, on ont-1k 20 against WFA's 21 µs, and on 100 kbp pairs at 6 and 7% 3.57 and 3.87 against A\*PA2-full's 5.25 and 5.0 ms, since pairs whose errors spread along them keep their exact seeds down to a fifth chained. On every real dataset it is faster than both approximate aligners, and inexact seeds speed the long reads up rather than down.

<details>
<summary>The M2's tables</summary>

Real datasets, mean time per alignment, median in brackets:

| dataset | pairs | mean length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA | WFA-adaptive* | Block Aligner* |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| ont-1k | 1221 | 0.818 kbp | 20 µs (19 µs) | 64 µs (70 µs) | 40 µs (46 µs) | 389 µs (357 µs) | 70 µs (77 µs) | 27 µs (28 µs) | 21 µs (21 µs) | 25 µs (25 µs), 93% optimal | 54 µs (61 µs), 85% optimal |
| ont-10k | 277 | 3.6 kbp | 143 µs (92 µs) | 293 µs (204 µs) | 209 µs (141 µs) | 8.8 ms (2.02 ms) | 705 µs (368 µs) | 367 µs (150 µs) | 324 µs (133 µs) | 186 µs (92 µs), 61% optimal | 235 µs (176 µs), 61% optimal |
| ont-50k | 104 | 9.52 kbp | 600 µs (188 µs) | 1.01 ms (397 µs) | 890 µs (278 µs) | 121 ms (5.53 ms) | 4.61 ms (1.13 ms) | 3.52 ms (435 µs) | 3.12 ms (368 µs) | 1.07 ms (225 µs), 52% optimal | 1.03 ms (463 µs), 61% optimal |
| ont-500k | 15 | 638 kbp | 110 ms (59.7 ms) | 160 ms (75.9 ms) | 635 ms (443 ms) | — | 4.18 s (2.54 s), 4/15 | 2.22 s (1.07 s), 9/15 | — | 524 ms (343 ms), 53% optimal | 875 ms (879 ms) |
| ont-500k-genvar | 15 | 651 kbp | 142 ms (95.4 ms) | 199 ms (133 ms) | 745 ms (609 ms) | — | 4.44 s (4.45 s), 4/15 | 2.95 s (2.97 s), 7/15 | 4.6 s (4.6 s), 4/15 | 263 ms (124 ms), 7% optimal | 1.06 s (856 ms), 33% optimal of 3 |
| sars-cov-2 | 33 | 29.6 kbp | 219 µs (143 µs) | 1.38 ms (1.38 ms) | 705 µs (661 µs) | 4.12 ms (1.1 ms) | 7.17 ms (6.89 ms) | 471 µs (211 µs) | 457 µs (181 µs) | 310 µs (210 µs), 97% optimal | 2.83 ms (2.68 ms), 30% optimal |

Divergence, 100 kbp pairs:

| divergence | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| 0% | 86 µs | 3.83 ms | 1.5 ms | 4.16 ms | 22 ms | 344 µs | 261 µs |
| 1% | 906 µs | 4.57 ms | 3.55 ms | 3.94 ms | 32.5 ms | 1.47 ms | 1.64 ms |
| 2% | 3.16 ms | 4.6 ms | 5.13 ms | 4.14 ms | 43.7 ms | 5.11 ms | 4.73 ms |
| 3% | 2.2 ms | 4.71 ms | 9.23 ms | 4.44 ms | 63.3 ms | 10.2 ms | 9.19 ms |
| 4% | 2.26 ms | 4.87 ms | 9 ms | 5.03 ms | 66.2 ms | 16.6 ms | 15 ms |
| 5% | 2.58 ms | 5 ms | 17.6 ms | 6.46 ms | 105 ms | 26.8 ms | 22.1 ms |
| 6% | 3.57 ms | 5.25 ms | 15.9 ms | 9.85 ms | 106 ms | 36.1 ms | 30.9 ms |
| 7% | 3.87 ms | 5 ms | 14.9 ms | 20.2 ms | 111 ms | 47.7 ms | 40.7 ms |
| 8% | 4.13 ms | 7.58 ms | 14 ms | 20.4 ms | 112 ms | 63.7 ms | 50.3 ms |
| 9% | 5.74 ms | 9.71 ms | 12.7 ms | 21.2 ms | 116 ms | 77.2 ms | 62.5 ms |
| 10% | 5.83 ms | 8.89 ms | 32.8 ms | 23.1 ms | 190 ms | 93.1 ms | 74 ms |
| 11% | 5.89 ms | 12.1 ms | 31.4 ms | 27.2 ms | 189 ms | 111 ms | 87.3 ms |
| 12% | 6.05 ms | 10 ms | 30.4 ms | 28.5 ms | 192 ms | 133 ms | 104 ms |
| 13% | 6.31 ms | 17 ms | 29.5 ms | 33 ms | 194 ms | 148 ms | 118 ms |
| 14% | 7.35 ms | 16.1 ms | 28.6 ms | 40.4 ms | 197 ms | 182 ms | 132 ms |
| 15% | 7.55 ms | 14.4 ms | 27.6 ms | 53.4 ms | 200 ms | 206 ms | 164 ms |

Length at 5% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| 3 kbp | 36 µs | 150 µs | 66 µs | 178 µs | 346 µs | 51 µs | 40 µs |
| 10 kbp | 195 µs | 582 µs | 306 µs | 546 µs | 1.94 ms | 289 µs | 255 µs |
| 30 kbp | 712 µs | 1.53 ms | 1.85 ms | 1.61 ms | 11.8 ms | 2.31 ms | 2.1 ms |
| 100 kbp | 2.39 ms | 4.96 ms | 17.6 ms | 6.51 ms | 104 ms | 25.3 ms | 23.2 ms |
| 300 kbp | 9.8 ms | 17.3 ms | 74.5 ms | 18 ms | 626 ms | 250 ms | 200 ms |
| 1 Mbp | 33.1 ms | 71.4 ms | 1.04 s | 73 ms | 7.33 s (2/4) | 2.65 s | 2.58 s |

Length at 15% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| 3 kbp | 80 µs | 170 µs | 135 µs | 999 µs | 407 µs | 189 µs | 167 µs |
| 10 kbp | 344 µs | 804 µs | 760 µs | 3.91 ms | 3.35 ms | 1.87 ms | 1.55 ms |
| 30 kbp | 1.8 ms | 3.89 ms | 2.63 ms | 14.3 ms | 19.2 ms | 15.8 ms | 13.1 ms |
| 100 kbp | 7.52 ms | 14.4 ms | 27.3 ms | 52.7 ms | 199 ms | 201 ms | 143 ms |
| 300 kbp | 37.7 ms | 134 ms | 348 ms | 227 ms | 2.1 s | 1.72 s | 1.29 s |
| 1 Mbp | 306 ms | 1.56 s | 1.6 s | 949 ms | 15.7 s (1/4) | 17.8 s (1/4) | — |

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
| today | 110 ms (59.7 ms) | 142 ms (95.4 ms) |

</details>
