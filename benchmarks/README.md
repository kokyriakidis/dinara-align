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
| dinara-align                                                                    | this tree | MPL-2.0    | the read batches and the affine pairs, on the CPU and on the GPU where one answers |
| dinara-align (bit-parallel, …)                                                  | this tree | MPL-2.0    | the edit-distance workloads: `distance` and `align` at unit costs, ported from A\*PA2-simple, on one thread and on all of them |
| [hyalite](https://github.com/Psy-Fer/hyalite)                                   | `0189bcb` | MIT        | the read batches and the affine pairs, global mode (`Mode::Nw`), one CPU thread with its NEON or AVX2 kernels |
| [A\*PA, A\*PA2](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) | `bf2e14e` | MPL-2.0    | the edit-distance workloads, because it computes edit distance and nothing else |
| [Edlib](https://github.com/Martinsos/edlib), [WFA2-lib](https://github.com/smarco/WFA2-lib) through [pa-bench](https://github.com/pairwise-alignment/pa-bench) | `af7a50d` | MIT | the edit-distance workloads: the other exact aligners A\*PA2's evaluation compares against, with its parameters |
| [Edlib](https://github.com/Martinsos/edlib) v1.2.7, [KSW2](https://github.com/lh3/ksw2), [WFA2-lib](https://github.com/smarco/WFA2-lib) v2.3.6, [parasail](https://github.com/jeffdaily/parasail), [SSW](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) | `ec2310e`, `289609b`, `bcf473a`, `fb985ee`, `a66636b` | MIT | the other modes: free ends, seed extension, two-piece gaps and substitution tables (see Other Modes) |

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

Measured with `pixi run bench --full`, which took about five minutes, on the Intel Core i9-7900X of A\*PA2's results below (ten cores at a fixed 3.3 GHz, turbo boost and hyper-threading off, 91 GB) with an NVIDIA GeForce RTX 2070, unpinned on the otherwise idle machine so dinara-align's batches use every core, every tool built for the machine's own instruction set, AVX-512 included (`--cpu native`, the default), dinara-align at `6dcb6fa`.
Hyalite built with Rust 1.92, A\*PA with the nightly its repository pins, both with `-C target-cpu=native`.
A dash marks a task the tool does not offer, or a workload it is not run on.

| workload | task | dinara-align (cpu) | dinara-align (gpu) | dinara-align (bit-parallel, 1 thread) | hyalite | a*pa2-full | a*pa2-simple | a*pa2-nw | a*pa | edlib | biwfa | wfa | agree |
| :-- | :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| reads-150bp | score | 12.3 ms | 7.27 ms | — | 282 ms | — | — | — | — | — | — | — | ✓ |
| reads-150bp | alignment | 15.8 ms | 18.7 ms | — | 2.09 s | — | — | — | — | — | — | — | ✓ |
| reads-1kbp | score | 28.5 ms | 8.59 ms | — | 1.11 s | — | — | — | — | — | — | — | ✓ |
| reads-1kbp | alignment | 42.2 ms | 16.6 ms | — | 14.8 s | — | — | — | — | — | — | — | ✓ |
| affine-1k | score | 126 µs | 899 µs | — | 1.12 ms | — | — | — | — | — | — | — | ✓ |
| affine-1k | alignment | 200 µs | 1.59 ms | — | 14.5 ms | — | — | — | — | — | — | — | ✓ |
| affine-10k | score | 2.92 ms | 6.46 ms | — | 142 ms | — | — | — | — | — | — | — | ✓ |
| affine-10k | alignment | 4.01 ms | 9.42 ms | — | 1.83 s | — | — | — | — | — | — | — | ✓ |
| affine-100k | score | 244 ms | 82.2 ms | — | 15.5 s | — | — | — | — | — | — | — | ✓ |
| affine-100k | alignment | 495 ms | 177 ms | — | 186 s | — | — | — | — | — | — | — | ✓ |
| edit-1k-1% | score | — | — | 2 µs | — | 43 µs | 17 µs | 19 µs | — | 27 µs | 3 µs | 3 µs | ✓ |
| edit-1k-1% | alignment | — | — | 2 µs | — | 45 µs | 18 µs | 205 µs | 307 µs | 107 µs | 6 µs | 4 µs | ✓ |
| edit-1k-5% | score | — | — | 5 µs | — | 43 µs | 17 µs | 19 µs | — | 26 µs | 8 µs | 8 µs | ✓ |
| edit-1k-5% | alignment | — | — | 6 µs | — | 49 µs | 20 µs | 198 µs | 289 µs | 110 µs | 17 µs | 10 µs | ✓ |
| edit-1k-15% | score | — | — | 12 µs | — | 33 µs | 17 µs | 19 µs | — | 65 µs | 41 µs | 41 µs | ✓ |
| edit-1k-15% | alignment | — | — | 31 µs | — | 64 µs | 45 µs | 205 µs | 875 µs | 159 µs | 73 µs | 48 µs | ✓ |
| edit-10k-1% | score | — | — | 17 µs | — | 802 µs | 782 µs | 1.42 ms | — | 388 µs | 39 µs | 39 µs | ✓ |
| edit-10k-1% | alignment | — | — | 22 µs | — | 827 µs | 146 µs | 17.5 ms | 2.81 ms | 1.97 ms | 82 µs | 50 µs | ✓ |
| edit-10k-5% | score | — | — | 138 µs | — | 584 µs | 1.37 ms | 1.39 ms | — | 599 µs | 460 µs | 461 µs | ✓ |
| edit-10k-5% | alignment | — | — | 276 µs | — | 659 µs | 282 µs | 17.4 ms | 2.56 ms | 2.56 ms | 919 µs | 568 µs | ✓ |
| edit-10k-15% | score | — | — | 169 µs | — | 859 µs | 2.36 ms | 1.42 ms | — | 1.67 ms | 3.11 ms | 3.11 ms | ✓ |
| edit-10k-15% | alignment | — | — | 351 µs | — | 1.26 ms | 886 µs | 17.8 ms | 80.1 ms | 4.54 ms | 6.31 ms | 3.7 ms | ✓ |
| edit-100k-1% | score | — | — | 764 µs | — | 5.33 ms | 151 ms | 134 ms | — | 8.35 ms | 1.72 ms | 1.68 ms | ✓ |
| edit-100k-1% | alignment | — | — | 1.07 ms | — | 5.53 ms | 3.1 ms | 1.81 s | 38 ms | 39.3 ms | 3.44 ms | 2.39 ms | ✓ |
| edit-100k-5% | score | — | — | 3.8 ms | — | 5.66 ms | 226 ms | 134 ms | — | 57.7 ms | 37.4 ms | 37.4 ms | ✓ |
| edit-100k-5% | alignment | — | — | 4.54 ms | — | 6.35 ms | 14.1 ms | 1.81 s | 33.7 ms | 129 ms | 77.1 ms | 51.7 ms | ✓ |
| edit-100k-15% | score | — | — | 8.26 ms | — | 19.8 ms | 200 ms | 134 ms | — | 88.2 ms | 297 ms | 298 ms | ✓ |
| edit-100k-15% | alignment | — | — | 14 ms | — | 23.2 ms | 22.9 ms | 1.81 s | 13.9 s | 251 ms | 608 ms | 684 ms | ✓ |

## Reading the Numbers

- **Every answer agreed**, so each row compares tools answering the same question. The answers are the optimal scores and distances: a sum and a position-weighted sum per workload, so a reordered batch cannot pass.
- **On edit distance, dinara-align on one thread beats or matches every exact aligner on every workload**, distance and alignment alike: 22 against WFA's 50 µs and A\*PA2-simple's 146 µs aligning 10 kbp at 1%, 276 against A\*PA2-simple's 282 µs at 10 kbp and 5%, a tie, 4.54 against A\*PA2-full's 6.35 ms at 100 kbp and 5%, 14 against A\*PA2-simple's 22.9 ms at 100 kbp and 15%, 764 µs against WFA's 1.68 ms scoring 100 kbp at 1%, and 12 against A\*PA2-simple's 17 µs scoring 1 kbp at 15%, where a short pair's diagonal transition gives way to the whole matrix once it would cost more. Its alignment follows a fixed rule for ties (see `Ties`), traced a tile at a time by a forward search, which costs about twice what its alignment of divergent long pairs took before that rule. The ties are the 1 kbp pairs at 1 and 5%, a few microseconds each, where BiWFA and WFA run the same diagonal transition.
- **Near-identical pairs run diagonal transition**, as WFA does, before any band: one front keeping its history for a close alignment, and two from both ends otherwise, as BiWFA scores, keeping both histories for an alignment and tracing it back through each from where they met. It stays on while it costs less than the band would, which at 1 to 3% divergence covers most pairs up to tens of kbp.
- **A band re-aims its first bound as it sweeps.** The bound starts from a projection of the first few edits, which strays by up to half either way. At an eighth, a quarter and half of the columns, the band projects the distance from its own climb, hundreds of edits in, lowering the bound when it was set too high and giving the round up early when it was set too low. Only a distance within the final bound is accepted, so the answer stays exact.
- **Long, moderately divergent pairs prune with A\*PA2-full's seed heuristic.** From about 1,500 projected edits up to one in seven bases, the first sequence is cut into 12-base seeds, their exact matches in the second are found by a rolling hash, and matches that cannot shorten any path over the next 14 seeds are dropped, as A\*PA's local pruning drops them. The band then keeps a row only while its score plus the seeds still ahead, less the longest chain of matches it can still reach, fits the bound. At 2 to 3% divergence that bound at the origin lands within 1% of the distance, so the band starts just past it and finishes in one round. A\*PA2-full also prunes matches between rounds, which dinara-align leaves out, as its rounds now rarely number more than two. On AVX-512, whose band runs fast beside the seeds' setup, seeds are used only from 86 kbp, where they start to pay there.
- **Long, divergent pairs match seeds within one edit**, as A\*PA's `r = 2` does. Past about one edit in fifteen bases most exact 12-base seeds are broken and the bound at the origin falls far short: 8,277 for a distance of 12,196 on a 100 kbp pair at 15%. Seeds of 16 bases that may match with one edit, charging two edits where they match nowhere, put it at 11,304. Each match is found by the half of it that matches exactly, a check on its quarters turns most chance lookups away, and an exact match's one-edit neighbours are left out, which keeps the bound a lower bound. They are used from 64 kbp, where the band they save outgrows their setup, when the projection says one edit in ten bases or under 40% of the exact seeds chain.
- **Global affine scores run a wavefront from both ends.** Under a table of one match and one mismatch score, as `Scoring.dna()` is, the match reward folds away for a global alignment (Eizenga and Lindquist), and WFA's three-layer wavefront then searches by cost, from the origin and from the corner at once until they meet: 244 ms against hyalite's 15.5 s at 100 kbp, 126 µs against 1.12 ms at 1 kbp. A pair whose projected wavefront would cost more than a full sweep is handed to the sweep, which runs sixteen cells at a time by anti-diagonal and also serves local scores.
- **Global affine alignments are traced back through the same wavefront's fronts**, five bytes a diagonal kept from each end and the path walked from where they met, a pair whose fronts would pass 80 MB split where an optimal path crosses: 200 µs against hyalite's 14.5 ms at 1 kbp, 4.01 ms against 1.83 s at 10 kbp and 495 ms against 186 s at 100 kbp; dinara-align's own GPU sweep is slower on the 1 and 10 kbp pairs, 1.59 and 9.42 ms, and faster on the 100 kbp one, 177 ms. A tie between optimal paths may resolve differently from the GPU's, the score the same. Under any other table the matrix is swept sixteen cells at a time by anti-diagonal, each pair's score read by comparison, by byte shuffle from a register for a table of up to sixteen entries, or gathered from memory for a larger one: a global alignment stores only the band of diagonals its score bounds, a local one is found by its span and aligned globally (see Other Modes), and past `max_memory`, 80 MB by default, a global one splits on rows in linear space, Myers and Miller's way.
- **The other rivals run on one CPU thread.** dinara-align's CPU batches spread their pairs over every thread, one pair a thread, which the reads rows show; the single affine pairs, and every edit-distance pair, run on one.
- **dinara-align's GPU times include the overheads a caller pays**: opening the device context on every call, copying sequences over, and copying results back. That is why a single short pair is slower on the GPU than on the CPU, while the reads batches, every pair in one launch, run fastest there: 7.27 against 12.3 ms scoring the 150 bp reads.
- **hyalite's traceback budget is 1 GiB.** It keeps the whole matrix when it fits and switches to a checkpointed sweep above that, as the 100 kbp pairs do. dinara-align switches to its linear-space recursion above six million cells on the host and one million on the device.
- **A\*PA answers only edit distance**, and its time depends on how similar the sequences are, as do the bit-parallel path's and dinara-align's global affine score and alignment; the local alignments and hyalite's columns do not.
- **One machine, one fixed pair per workload, one run.** The rows of microseconds move by a few percent between runs; read a difference of that size as a tie. Rerun on your own hardware before quoting any of this.

## A\*PA2's Evaluation

`pa_bench.py` runs the exact aligners of A\*PA2's evaluation, with its parameters, and WFA beside BiWFA, on its own datasets: Oxford Nanopore reads and SARS-CoV-2 genomes from pa-bench's release, and the uniform-error pairs pa-generate regenerates exactly from the evaluation's seed and lock.
Each dataset is a fixed shuffled sample of about 2 Mbp and at least four pairs, the same for every tool; each tool aligns its pairs with traceback, once each, as pa-bench times them, within a budget of five seconds, and every cost is checked against every other tool's and against the costs A\*PA2's published results recorded.
The rivals are pinned, so their results are kept and reused until `--fresh`; a second table gives each tool's peak resident memory.
dinara-align also runs as a batch on all ten threads, `alignments` over the whole sample, one pair a thread, whose column is the batch's wall-clock time over its pairs, a throughput where the others are each one pair's latency.

Measured on the i9-7900X of A\*PA2's results below, unpinned on the otherwise idle machine, every tool built for its own instruction set, dinara-align at `2673fe3`; a count marks a tool its budget, or WFA's 8 GiB cap, stopped partway, and more than the budget one that finished no pair:

| dataset | pairs | mean length | dinara-align (bit-parallel, 1 thread) | dinara-align (bit-parallel, batch, 10 threads) | a*pa2-full | a*pa2-simple | a*pa | edlib | biwfa | wfa | agree |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| ont-1k | 1221 of 12477 | 0.818 kbp | 22 µs | 3 µs | 81 µs | 48 µs | 616 µs | 142 µs | 47 µs | 31 µs | ✓ |
| ont-10k | 277 of 5000 | 3.6 kbp | 164 µs | 24 µs | 393 µs | 248 µs | 15.5 ms | 1.12 ms | 885 µs | 598 µs | ✓ |
| ont-50k | 104 of 10000 | 9.52 kbp | 590 µs | 87 µs | 1.28 ms | 940 µs | 221 ms (25/104) | 6.15 ms | 9.1 ms | 6.36 ms | ✓ |
| ont-500k | 4 of 50 | 638 kbp | 116 ms | 64.9 ms | 283 ms | 686 ms | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| ont-500k-genvar | 4 of 48 | 659 kbp | 143 ms | 48.4 ms | 225 ms | 511 ms | > 5 s | 5.2 s (1/4) | > 5 s | > 5 s | ✓ |
| sars-cov-2 | 33 of 10000 | 29.6 kbp | 279 µs | 96 µs | 2.01 ms | 728 µs | 6.56 ms | 8.63 ms | 998 µs | 988 µs | ✓ |
| Uniform-t10000000-n3000-e0.05 | 333 of 3333 | 3 kbp | 41 µs | 7 µs | 175 µs | 61 µs | 267 µs | 857 µs | 88 µs | 58 µs | ✓ |
| Uniform-t10000000-n10000-e0.05 | 99 of 1000 | 10 kbp | 285 µs | 44 µs | 787 µs | 292 µs | 937 µs | 3.97 ms | 615 µs | 439 µs | ✓ |
| Uniform-t10000000-n30000-e0.05 | 33 of 333 | 30 kbp | 947 µs | 151 µs | 2.45 ms | 1.72 ms | 3.11 ms | 14.9 ms | 5.82 ms | 4.31 ms | ✓ |
| Uniform-t10000000-n100000-e0.05 | 10 of 100 | 100 kbp | 4.82 ms | 1.12 ms | 8.31 ms | 15.5 ms | 12.3 ms | 147 ms | 62.9 ms | 46.5 ms | ✓ |
| Uniform-t10000000-n300000-e0.05 | 4 of 33 | 300 kbp | 19.6 ms | 6.72 ms | 27.5 ms | 65.5 ms | 38.5 ms | 758 ms | 536 ms | 454 ms | ✓ |
| Uniform-t10000000-n1000000-e0.05 | 4 of 10 | 1e+03 kbp | 61.9 ms | 19.3 ms | 125 ms | 862 ms | 158 ms | > 5 s | > 5 s | > 5 s | ✓ |
| Uniform-t10000000-n3000-e0.15 | 333 of 3333 | 3 kbp | 114 µs | 16 µs | 202 µs | 153 µs | 1.56 ms | 910 µs | 417 µs | 311 µs | ✓ |
| Uniform-t10000000-n10000-e0.15 | 99 of 1000 | 10 kbp | 407 µs | 55 µs | 916 µs | 786 µs | 6.55 ms | 5.9 ms | 4.86 ms | 2.98 ms | ✓ |
| Uniform-t10000000-n30000-e0.15 | 33 of 333 | 30 kbp | 1.69 ms | 222 µs | 4.07 ms | 2.62 ms | 23.5 ms | 24.7 ms | 42.2 ms | 26.8 ms | ✓ |
| Uniform-t10000000-n100000-e0.15 | 9 of 100 | 100 kbp | 11.1 ms | 1.58 ms | 16.4 ms | 24.5 ms | 99.8 ms | 263 ms | 454 ms | 311 ms | ✓ |
| Uniform-t10000000-n300000-e0.15 | 4 of 33 | 300 kbp | 44.8 ms | 14.2 ms | 125 ms | 293 ms | 474 ms | 2.54 s (2/4) | 4.06 s (1/4) | 4.65 s (1/4) | ✓ |
| Uniform-t10000000-n1000000-e0.15 | 4 of 10 | 1e+03 kbp | 265 ms | 71.1 ms | 1.33 s | 1.35 s | 2.09 s (2/4) | > 5 s | > 5 s | > 5 s | ✓ |


## Affine Costs

`pa_bench.py --affine x,o,e` runs the same samples at gap-affine costs as WFA counts them, a mismatch `x` and a gap of `k` letters `o + k e`, against the exact aligners that take them: WFA2-lib keeping every front (WFA), its lowest-memory mode (BiWFA), and KSW2's banded SSE kernel with band doubling, on x86-64 alone.
dinara-align answers with `align` under `Costs.affine`, a CIGAR as the others hand back, from a wavefront grown from both ends at once, which keeps five bytes a diagonal for its traceback and splits a pair whose fronts would pass 80 MB where an optimal path crosses, as BiWFA does.

```bash
pixi run bench-astarpa2 --affine 4,6,2
```

At WFA's costs (4, 6, 2) on the i9-7900X above, pinned to one core, every tool built for its own instruction set, dinara-align at `2673fe3`, each tool five seconds a sample; a count marks a tool its budget stopped partway, and more than the budget one that finished no pair:

| dataset | pairs | mean length | dinara-align | WFA | BiWFA | KSW2 | agree |
| :-- | --: | --: | --: | --: | --: | --: | :-: |
| ont-1k | 1221 of 12477 | 0.818 kbp | 139 µs | 190 µs | 329 µs | 1.13 ms | ✓ |
| ont-10k | 277 of 5000 | 3.6 kbp | 1.84 ms | 4.84 ms | 6.96 ms | 31.1 ms (161/277) | ✓ |
| ont-50k | 104 of 10000 | 9.52 kbp | 22.5 ms | 69.2 ms (82/104) | 69.4 ms (82/104) | 366 ms (16/104) | ✓ |
| ont-500k | 4 of 50 | 638 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| ont-500k-genvar | 4 of 48 | 659 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| sars-cov-2 | 33 of 10000 | 29.6 kbp | 877 µs | 2.08 ms | 1.98 ms | 89.6 ms | ✓ |
| Uniform-t10000000-n3000-e0.05 | 333 of 3333 | 3 kbp | 347 µs | 589 µs | 1.01 ms | 6.27 ms | ✓ |
| Uniform-t10000000-n10000-e0.05 | 99 of 1000 | 10 kbp | 2.37 ms | 6.42 ms | 9.38 ms | 97.9 ms (52/99) | ✓ |
| Uniform-t10000000-n30000-e0.05 | 33 of 333 | 30 kbp | 18.8 ms | 56.8 ms | 75.4 ms | 1.36 s (4/33) | ✓ |
| Uniform-t10000000-n100000-e0.05 | 10 of 100 | 100 kbp | 310 ms | 765 ms (7/10) | 780 ms (7/10) | > 5 s | ✓ |
| Uniform-t10000000-n300000-e0.05 | 4 of 33 | 300 kbp | 3.15 s (1/4) | > 5 s | > 5 s | > 5 s | ✓ |
| Uniform-t10000000-n1000000-e0.05 | 4 of 10 | 1e+03 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| Uniform-t10000000-n3000-e0.15 | 333 of 3333 | 3 kbp | 1.42 ms | 3.64 ms | 5.9 ms | 19 ms (263/333) | ✓ |
| Uniform-t10000000-n10000-e0.15 | 99 of 1000 | 10 kbp | 12.8 ms | 39.2 ms | 54.3 ms (93/99) | 204 ms (25/99) | ✓ |
| Uniform-t10000000-n30000-e0.15 | 33 of 333 | 30 kbp | 156 ms (32/33) | 391 ms (13/33) | 461 ms (11/33) | 2.48 s (2/33) | ✓ |
| Uniform-t10000000-n100000-e0.15 | 9 of 100 | 100 kbp | 2.09 s (2/9) | > 5 s | 5.16 s (1/9) | > 5 s | ✓ |
| Uniform-t10000000-n300000-e0.15 | 4 of 33 | 300 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| Uniform-t10000000-n1000000-e0.15 | 4 of 10 | 1e+03 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |

- **dinara-align is the fastest on every dataset any tool finished,** 1.4 to 3.1 times faster than WFA: 139 against 190 µs on ont-1k, 1.84 against 4.84 ms on ont-10k, 877 µs against 2.08 ms on the SARS-CoV-2 genomes, and 12.8 against 39.2 ms on 10 kbp at 15%. Its CIGAR follows a fixed rule for ties (see `Ties`), its search grown on to the optimum and traced back from the far end, which costs the short pairs most: 139 µs on ont-1k where it took 111 before that rule. BiWFA beats WFA only on the genomes and where WFA's budget runs out.
- **It alone finishes every ont-50k read,** in 22.5 ms on average, where WFA and BiWFA each finish 82 of the 104 within the budget, and it finishes more of the 100 and 300 kbp pairs than anyone.
- **KSW2's band doubles from eight diagonals against a loose bound,** so it times out on most reads past 1 kbp.
- **No aligner finishes a 500 kbp read or a 1 Mbp pair in five seconds:** the wavefront's work grows with the square of the cost, and gap-affine costs have no seed heuristic here to prune it (see TODO.md).

Growth of peak memory aligning each pair, median / largest, as A\*PA2's Table 10 measures it:

| dataset | dinara-align | WFA | BiWFA | KSW2 |
| :-- | --: | --: | --: | --: |
| ont-1k | 0.0 / 1.0 MB | 0.0 / 3.0 MB | 0.0 / 1.5 MB | 0.0 / 0.9 MB |
| ont-10k | 0.0 / 22 MB | 0.0 / 97 MB | 0.0 / 2.2 MB | 0.0 / 58 MB |
| ont-50k | 0.0 / 85 MB | 0.0 / 1420 MB | 0.0 / 4.4 MB | 0.0 / 317 MB |
| sars-cov-2 | 0.0 / 11 MB | 0.0 / 51 MB | 0.0 / 2.4 MB | 0.0 / 211 MB |
| Uniform-t10000000-n3000-e0.05 | 0.0 / 1.9 MB | 0.0 / 5.4 MB | 0.0 / 1.8 MB | 0.0 / 3.2 MB |
| Uniform-t10000000-n10000-e0.05 | 0.0 / 7.1 MB | 0.0 / 29 MB | 0.0 / 1.8 MB | 0.0 / 40 MB |
| Uniform-t10000000-n30000-e0.05 | 0.0 / 46 MB | 0.0 / 212 MB | 0.0 / 2.9 MB | 31 / 455 MB |
| Uniform-t10000000-n100000-e0.05 | 0.0 / 86 MB | 12 / 2413 MB | 0.0 / 9.1 MB | — |
| Uniform-t10000000-n300000-e0.05 | 97 / 97 MB | — | — | — |
| Uniform-t10000000-n3000-e0.15 | 0.0 / 4.5 MB | 0.0 / 18 MB | 0.0 / 1.4 MB | 0.0 / 12 MB |
| Uniform-t10000000-n10000-e0.15 | 0.0 / 35 MB | 0.0 / 160 MB | 0.0 / 3.2 MB | 0.0 / 79 MB |
| Uniform-t10000000-n30000-e0.15 | 0.0 / 101 MB | 0.5 / 1416 MB | 0.0 / 5.5 MB | 796 / 796 MB |
| Uniform-t10000000-n100000-e0.15 | 85 / 85 MB | — | 16 / 16 MB | — |

dinara-align keeps under a tenth of WFA's memory on the largest pairs (86 against 2413 MB at 100 kbp), as its traceback keeps a byte of flags and the column of one front a diagonal, from fronts each half as long; BiWFA, keeping only its last few fronts, stays smallest.

## Local and Overlap Alignment

Local alignment, Smith-Waterman, overlap alignment, every end gap free, and a read placed whole in a window with a reward, each with its CIGAR, against the aligners that offer them:

```bash
pixi run bench-local   # builds SSW, parasail and abPOA the first time; then a few minutes
```

Every tool scores a match 2, a mismatch -4 and a gap of `k` letters `6 + 2k` (dinara-align's `Mode.local(2)` and `Mode.overlap(2)` under `Costs.affine(4, 6, 2)`), and every tool's scores must agree on every pair, or the run fails; none disagreed.
A tool's time is the faster of two passes over a workload, its mean per pair, on one thread of the Skylake-X, pinned (see A\*PA2's results below for the machine), every tool built for its own instruction set, AVX-512 included, at `b132767`.
abPOA aligns to a graph, so its time includes adding the reference to one, as any pairwise use of it pays; SSW and abPOA have no overlap or infix mode.
The infix row scores with a reward, `Mode.INFIX.with_match_score(2)`, as parasail's `sg_dx` and hyalite's HW do; without one, dinara-align's `Mode.INFIX` minimizes the costs, Edlib's way, on the bit-parallel sweep at unit costs or the wavefront at any others.

| workload | dinara-align | SSW | parasail | abPOA | hyalite |
| :-- | --: | --: | --: | --: | --: |
| local, 600 short noisy pairs (20 to 700 bp) | **36 µs** | 52 µs | 67 µs | 175 µs | 845 µs |
| local, 1 kbp read at 10% in a 10 kbp window | **1.62 ms** | 3.21 ms | 4.77 ms | 16.8 ms | 162 ms |
| local, 10 kbp at 5% against 12 kbp | **46.5 ms** | 85.1 ms | 288 ms | 205 ms | 2.29 s |
| overlap, 2 kbp reads overlapping by 0.5 to 1.5 kbp | **2.29 ms** | — | 15.3 ms | — | 42.7 ms |
| infix with a reward, 1 kbp read at 10% placed whole in a 3 kbp window | **1.46 ms** | — | 14.9 ms | — | 45.7 ms |

- **dinara-align is the fastest on every workload**, 1.4 to 2 times SSW on local alignment, and 7 and 10 times parasail on overlaps and scored infixes.
- **Its sweep runs by anti-diagonal in 16-bit lanes along the shorter sequence,** while the scores fit, and finds only the best end; an extension back from that end, stopping once it earns the sweep's score, gives the alignment, traced as it searched. SSW sweeps striped (Farrar), whose lazy pass runs most of every column when a long alignment scores high, as on the 10 kbp pairs.

## Other Modes

Free ends, seed extension, two-piece gaps and substitution tables of more than one mismatch score, each against the aligners that offer the same mode, on DNA:

```bash
pixi run bench-modes   # builds Edlib, WFA2-lib's library and KSW2 the first time; then a few minutes
```

Every tool aligns every pair with its CIGAR on one thread of the Skylake-X, pinned, every tool built for its own instruction set, AVX-512 included, dinara-align at `9767ca4`; a tool's time is the faster of two passes over a workload, its mean per pair, and every tool's costs or scores must agree on every pair, or the run fails; none disagreed.
WFA2-lib runs exact, keeping every front, its WF-adaptive heuristic off.
KSW2 extends with no Z-drop: it gauges one by anti-diagonal and dinara-align as WFA2-lib does, a cost at a time, so the two would stop at different places and their times compare different work.
Deletions priced apart from insertions are left out: no rival offers them.

| workload | dinara-align | Edlib | WFA2-lib | KSW2 | parasail | SSW |
| :-- | --: | --: | --: | --: | --: | --: |
| infix at unit costs, 1 kbp read at 10% placed whole in a 3 kbp window | **143 µs** | 288 µs | 1.06 ms | — | — | — |
| infix at unit costs, 10 kbp read at 5% in a 30 kbp window | **3.31 ms** | 13.8 ms | 55.6 ms | — | — | — |
| prefix at unit costs, 1 kbp read at 10% against 2 kbp of reference | **45 µs** | 129 µs | 47 µs | — | — | — |
| infix at WFA's (4, 6, 2), the 1 kbp reads above | **3.10 ms** | — | 4.11 ms | — | — | — |
| extension, a match 2 at (4, 6, 2), 300 to 1,400 bases at 5% then noise | **779 µs** | — | — | 2.09 ms | — | — |
| two-piece gaps (4, 6, 2, 24, 1), 5 kbp at 5% with three indels of 100 to 400 bp | **9.47 ms** | — | 15.8 ms | 46.8 ms | — | — |
| a table, match 2, transition -2, transversion -4, global, 1 kbp at 10% | **557 µs** | — | — | — | 2.96 ms | — |
| the same table, local, 1 kbp read at 10% in a 10 kbp window | **2.62 ms** | — | — | — | 5.16 ms | 3.32 ms |

- **dinara-align is the fastest on every workload**: 2 and 4 times Edlib on infixes, 2.7 times KSW2 on extension, 1.7 times WFA2-lib on two-piece gaps and 5.3 times parasail on a table globally. The prefix search ties WFA2-lib's diagonal transition, 45 against 47 µs.
- **The prefix search sweeps only the band its bound allows**, which Edlib's prefix mode does not: the band's top moves down with the diagonal, no column past the read's length plus the best end so far is swept, and a try whose band has already emptied stops. The part of the reference it finds is then aligned globally by the edit-distance traceback.
- **A table of more than one mismatch score sweeps sixteen cells at a time**, by anti-diagonal, each pair's score read by byte shuffle from a register when the table has at most sixteen entries, as DNA's four letters do, and gathered from memory otherwise. A local alignment is found as SSW finds one: a 16-bit sweep for its end, one anchored there for the first cell that earns its score, its start, and the span between aligned globally. A global one stores only the band of diagonals its score bounds. Both used to sweep a cell at a time: 13.3 ms globally and 133 ms locally.

## A\*PA2's Results, Redone

A\*PA2's own results section ([curiouscoding.nl/posts/astarpa2](https://curiouscoding.nl/posts/astarpa2/#results)), redone with dinara-align beside the exact aligners it compares, and WFA, on a machine set up as A\*PA2's was: its real datasets, its uniform pairs swept over divergence and over length, and in place of its ablation, dinara-align's history on the long reads.

```bash
pixi run results-astarpa2 collect   # forty minutes the first time; after that only dinara-align runs, a few minutes
pixi run results-astarpa2 tables    # the tables below, into .cache/results/astarpa2-results.md
```

### Setup

As A\*PA2's evaluation runs them: one single-threaded job at a time, every pair aligned once with its traceback, which every aligner hands back as a CIGAR, the time the average wall clock per alignment, reading the data left out.

The machine is an Intel Core i9-7900X under Ubuntu 26.04 (Linux 7.0), with Mojo 1.1.0, set up as A\*PA2's i7-10750H was: every core fixed at 3.3 GHz, turbo boost and hyper-threading off, and each collection pinned to one core with `taskset`.
One thing differs: the jobs ran at normal priority, not A\*PA2's niceness −20, on an otherwise idle machine. Every tool is built for this machine's own instruction set, AVX-512 included (`--cpu native`, the default), as A\*PA2's own repository builds it; dinara-align at `b132767`.

Each dataset is a fixed shuffled sample, the same pairs for every aligner: about 2 Mbp of sequence and at least four pairs, and 20 Mbp, fifteen reads, of the two 500 kbp read sets.
Each aligner may spend 20 seconds on a sample; one that has not finished stops there, its numbers covering the pairs it finished, counted beside them, and one that finished none shows a dash.
WFA, which A\*PA2's evaluation leaves out, is WFA2-lib keeping every front for its traceback, BiWFA's faster but hungrier twin: its memory grows with the square of the distance, so it also stops at 8 GiB, as its budget would stop it.
On the real datasets, as in A\*PA2's results, two approximate aligners run beside the exact ones, with the evaluation's parameters: WFA-adaptive, WFA2-lib's lowest-memory mode dropping lagging diagonals at its defaults (10, 50, 10), and Block Aligner, with blocks from 0.1 to 1% of the input and, as it takes only affine costs, a gap's opening one more. As pa-bench does, each is held against an exact alignment at its own costs, the exact distance for WFA-adaptive and for Block Aligner BiWFA at the same affine costs, and the share it aligned optimally goes beside its time; a share over fewer pairs than it finished counts only those its exact reference finished within the budget.
Every cost is checked against every other aligner's on the pairs both aligned, and against the costs A\*PA2's published results recorded; a disagreement fails the run, and none occurred.
The rivals, and dinara-align's older commits, never change, so each is run on a sample once and kept; later runs time today's dinara-align alone.

### Real datasets

Mean time per alignment, median in brackets:

| dataset | pairs | mean length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA | WFA-adaptive* | Block Aligner* |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| ont-1k | 1221 | 0.818 kbp | 22 µs (22 µs) | 81 µs (90 µs) | 48 µs (54 µs) | 616 µs (561 µs) | 142 µs (141 µs) | 47 µs (49 µs) | 31 µs (32 µs) | 42 µs (43 µs), 93% optimal | 41 µs (47 µs), 85% optimal |
| ont-10k | 277 | 3.6 kbp | 165 µs (112 µs) | 393 µs (266 µs) | 248 µs (179 µs) | 15.5 ms (3.29 ms) | 1.12 ms (605 µs) | 885 µs (315 µs) | 598 µs (224 µs) | 418 µs (183 µs), 61% optimal | 202 µs (145 µs), 61% optimal |
| ont-50k | 104 | 9.52 kbp | 599 µs (245 µs) | 1.28 ms (515 µs) | 940 µs (306 µs) | 194 ms (10.6 ms), 81/104 | 6.15 ms (1.8 ms) | 9.1 ms (1.02 ms) | 6.36 ms (741 µs) | 2.61 ms (491 µs), 52% optimal | 654 µs (326 µs), 61% optimal |
| ont-500k | 15 | 638 kbp | 100 ms (82.4 ms) | 177 ms (116 ms) | 539 ms (372 ms) | — | 5.02 s (3.06 s), 4/15 | — | — | 1.11 s (678 ms), 53% optimal | 735 ms (721 ms) |
| ont-500k-genvar | 15 | 651 kbp | 136 ms (105 ms) | 211 ms (181 ms) | 625 ms (509 ms) | — | 5.26 s (5.21 s), 3/15 | 9.62 s (9.62 s), 2/15 | — | 610 ms (303 ms), 7% optimal | 847 ms (657 ms) |
| sars-cov-2 | 33 | 29.6 kbp | 282 µs (136 µs) | 2.01 ms (2.08 ms) | 728 µs (662 µs) | 6.56 ms (2.08 ms) | 8.63 ms (8.23 ms) | 998 µs (455 µs) | 988 µs (359 µs) | 639 µs (353 µs), 97% optimal | 2.48 ms (2.26 ms), 30% optimal |

- **Among the exact aligners, dinara-align has the lowest mean and median on every dataset.** On the short reads, the case A\*PA2's post leaves to BiWFA, 22 against WFA's 31 µs and BiWFA's 47, median 22 against 32, since a short pair no longer hands a band, which covers its whole matrix there, what the diagonal transition does faster; on ont-10k and ont-50k it takes about two thirds of A\*PA2-simple's time, 165 against 248 µs and 599 against 940 µs, and on the SARS-CoV-2 genomes 282 µs against A\*PA2-simple's 728 and WFA's 988.
- **On the 500 kbp reads it takes about three fifths of A\*PA2-full's time:** 100 against 177 ms on ont-500k and 136 against 211 ms on ont-500k-genvar, medians 82.4 against 116 and 105 against 181.
- **Edlib finishes only some of the long reads within 20 seconds, BiWFA two genvar reads and no ont-500k read, and A\*PA and WFA none;** their means cover the reads they finished.
- **The approximate aligners trade optimality for speed, and are still slower than dinara-align on every dataset.** The closest, Block Aligner, takes 202 and 654 µs on ont-10k and ont-50k against dinara-align's 165 and 599, aligning 61% of those reads optimally; WFA-adaptive aligns 7% of the genvar reads optimally. A\*PA2 reports the same shares: Block Aligner 85% on ont-1k, WFA-adaptive 93%.

### Divergence, 100 kbp pairs

Mean time per alignment, A\*PA with whichever of `r = 1` and `r = 2` was the faster at each divergence:

| divergence | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| 0% | 81 µs | 6.65 ms | 1.55 ms | 6.56 ms | 43.9 ms | 614 µs | 305 µs |
| 1% | 1.3 ms | 7.46 ms | 3.4 ms | 7.83 ms | 58.9 ms | 3.43 ms | 3.38 ms |
| 2% | 4.23 ms | 7.57 ms | 4.84 ms | 8.43 ms | 72.9 ms | 11.2 ms | 9.67 ms |
| 3% | 4.56 ms | 7.77 ms | 8.37 ms | 8.2 ms | 97.8 ms | 24 ms | 19.3 ms |
| 4% | 4.65 ms | 8.07 ms | 8.48 ms | 10.2 ms | 101 ms | 40.9 ms | 31.6 ms |
| 5% | 4.81 ms | 8.2 ms | 15.5 ms | 12.1 ms | 147 ms | 63.2 ms | 47.3 ms |
| 6% | 6.27 ms | 8.85 ms | 14.2 ms | 17.6 ms | 149 ms | 86.4 ms | 65.1 ms |
| 7% | 6.37 ms | 8.89 ms | 13.3 ms | 35.8 ms | 153 ms | 114 ms | 85.5 ms |
| 8% | 6.58 ms | 11.8 ms | 12.4 ms | 41.2 ms | 157 ms | 147 ms | 108 ms |
| 9% | 9.63 ms | 12 ms | 11.6 ms | 41.3 ms | 160 ms | 181 ms | 130 ms |
| 10% | 10.3 ms | 10.5 ms | 28.3 ms | 45.7 ms | 247 ms | 221 ms | 157 ms |
| 11% | 10.8 ms | 13.7 ms | 27.3 ms | 48.8 ms | 249 ms | 264 ms | 184 ms |
| 12% | 10.6 ms | 11.9 ms | 26.4 ms | 55.2 ms | 251 ms | 307 ms | 214 ms |
| 13% | 10.5 ms | 18 ms | 25.7 ms | 60.5 ms | 255 ms | 352 ms | 241 ms |
| 14% | 11.8 ms | 16.7 ms | 25.2 ms | 72 ms | 259 ms | 404 ms | 276 ms |
| 15% | 11.4 ms | 15.6 ms | 24.6 ms | 101 ms | 262 ms | 454 ms | 310 ms |

- **dinara-align is the fastest at every divergence.**
- **Near-identical pairs finish in the diagonal transition:** 81 µs against WFA's 305 at 0%, 1.3 against WFA's 3.38 ms at 1%.
- **From 2% to 8% the exact seeds' bound lands close to the distance:** dinara-align runs in 4.2 to 6.6 ms, the next aligner, A\*PA2-simple at 2% and A\*PA2-full beyond, in 4.8 to 11.8.
- **From 9% on, where fewer than a fifth of the exact seeds chain,** dinara-align rebuilds them to match within one edit and stays between 9.6 and 11.8 ms, its narrowest lead 10.3 against A\*PA2-full's 10.5 at 10%, a tie, while A\*PA2-full climbs to 16 to 18 ms by 13 to 15%; at 15% A\*PA2-simple takes 24.6 ms, A\*PA 101, Edlib 262, WFA 310 and BiWFA 454.

### Length

Mean time per alignment at 5% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| 3 kbp | 42 µs | 175 µs | 61 µs | 267 µs | 857 µs | 88 µs | 58 µs |
| 10 kbp | 287 µs | 787 µs | 292 µs | 937 µs | 3.97 ms | 615 µs | 439 µs |
| 30 kbp | 962 µs | 2.45 ms | 1.72 ms | 3.11 ms | 14.9 ms | 5.82 ms | 4.31 ms |
| 100 kbp | 4.79 ms | 8.31 ms | 15.5 ms | 12.3 ms | 147 ms | 62.9 ms | 46.5 ms |
| 300 kbp | 19.5 ms | 27.5 ms | 65.5 ms | 38.5 ms | 758 ms | 536 ms | 454 ms |
| 1 Mbp | 61.8 ms | 125 ms | 862 ms | 158 ms | 8.71 s (2/4) | 6.02 s (3/4) | 5.59 s (3/4) |

At 15% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| 3 kbp | 115 µs | 202 µs | 153 µs | 1.56 ms | 910 µs | 417 µs | 311 µs |
| 10 kbp | 410 µs | 916 µs | 786 µs | 6.55 ms | 5.9 ms | 4.86 ms | 2.98 ms |
| 30 kbp | 1.71 ms | 4.07 ms | 2.62 ms | 23.5 ms | 24.7 ms | 42.2 ms | 26.8 ms |
| 100 kbp | 11.1 ms | 16.4 ms | 24.5 ms | 99.8 ms | 263 ms | 454 ms | 311 ms |
| 300 kbp | 44.6 ms | 125 ms | 293 ms | 474 ms | 2.53 s | 4.04 s | 3.05 s |
| 1 Mbp | 265 ms | 1.33 s | 1.35 s | 2.13 s | 18.9 s (1/4) | — | — |

- **dinara-align is the fastest at every length and both divergences.** At 5% it aligns a 1 Mbp pair in 61.8 ms against A\*PA2-full's 125 and A\*PA's 158; at 15% in 265 ms against A\*PA2-full's 1.33 s and A\*PA's 2.13 s. Its narrowest leads are at 10 kbp, 287 against A\*PA2-simple's 292 µs at 5%, and at 100 kbp and 15%, 11.1 against A\*PA2-full's 16.4 ms.
- **At 15% its band still grows faster than the length,** but less than A\*PA2-full's: from 300 kbp to 1 Mbp, 3.3 times the length, dinara-align takes 5.9 times as long and A\*PA2-full 10.6 times, where A\*PA, whose pruning makes its heuristic nearly exact, takes 4.5 times.

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
| today | 100 ms (82.4 ms) | 136 ms (105 ms) |

- **On ont-500k the time falls from 1.32 s at the port of A\*PA2-simple's band doubling to 162 ms once the retries are aimed.**
- **On ont-500k-genvar the commits between the port and `b60ece0` were several times slower,** 4.3 to 8.5 s a read and up to eight reads past the budget: real reads gather their errors at the ends, and the projections those commits aimed by ran many times over the distance, so they swept bands far wider than needed. `b60ece0` trusted a projection only when two of them agreed and otherwise grew the band from what was known, bringing genvar back to 752 ms, and aiming the retries to 199 ms.
- **Inexact seeds made these reads slower here at first,** 162 to 181 ms and 199 to 244 ms: this machine's band, on AVX-512, runs fast beside its seeds' setup, so inexact seeds paid only on more divergent pairs; the cutoff now follows the vector width, 20% of the exact seeds chained on AVX-512 against 40% on narrower vectors.
- **A leaner seeds' setup since** (`fc0e5e4`: a branch-free one-edit test, only the seeds whose matches can reach the end, a leaner local pruning) brings both below their best before inexact seeds, and with retries aimed closer, AVX-512 sweeping two groups at a time (`697fef1`) and exact seeds' matching filtered (`bcdc7dc`) since, they take 100 ms on ont-500k and 136 ms on ont-500k-genvar, a CIGAR included (`3152d27`); the diagonal transition's slides by AVX-512 gathers (`bad5ba5`) took the short reads and low-divergence pairs down 7 to 19% since. Every CIGAR now follows a fixed rule for ties (`ea094ba`), the band's retraced a tile at a time by a forward search (`8f9c92f`).

### Memory

From `pixi run bench-astarpa2` on the same machine, over its samples of each dataset (the 1221 reads of ont-1k, four pairs of each long set). Each cell is the most memory aligning one pair added to the most the process had held before it, the measure of A\*PA2's Table 10, with the process's whole peak resident memory in brackets; the lowest of the aligners that finished every pair is bold, and † marks one that finished only some within its budget, a dash none:

| dataset | dinara-align | A\*PA2-full | A\*PA2-simple | A\*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| ***Real datasets*** |  |  |  |  |  |  |  |
| ont-1k | 0.4 MB (15 MB) | 0.2 MB (5 MB) | 0.4 MB (5 MB) | 0.2 MB (6 MB) | 0.5 MB (8 MB) | **0.1 MB** (7 MB) | 0.4 MB (8 MB) |
| ont-10k | 1.1 MB (17 MB) | 0.6 MB (6 MB) | **0.5 MB** (6 MB) | 20 MB (40 MB) | **0.5 MB** (8 MB) | 1.0 MB (8 MB) | 6.1 MB (39 MB) |
| ont-50k | 2.0 MB (20 MB) | 1.5 MB (8 MB) | **0.7 MB** (7 MB) | 83 MB (—) † | 1.0 MB (8 MB) | 1.4 MB (9 MB) | 71 MB (187 MB) |
| ont-500k | **50 MB** (67 MB) | 81 MB (89 MB) | 82 MB (90 MB) | — | — | — | — |
| ont-500k-genvar | 35 MB (59 MB) | 48 MB (56 MB) | **32 MB** (55 MB) | — | 4.7 MB (14 MB) † | — | — |
| sars-cov-2 | 2.2 MB (18 MB) | 1.9 MB (7 MB) | **0.6 MB** (6 MB) | 13 MB (19 MB) | 1.4 MB (8 MB) | 1.6 MB (10 MB) | 15 MB (26 MB) |
| ***Uniform pairs, 5% divergence*** |  |  |  |  |  |  |  |
| 3 kbp, 5% | 0.7 MB (15 MB) | 0.4 MB (5 MB) | **0.2 MB** (5 MB) | 0.3 MB (5 MB) | 0.9 MB (7 MB) | 0.4 MB (7 MB) | 0.8 MB (7 MB) |
| 10 kbp, 5% | 1.7 MB (15 MB) | 0.6 MB (6 MB) | **0.4 MB** (5 MB) | 0.5 MB (6 MB) | 0.9 MB (7 MB) | 1.1 MB (8 MB) | 2.0 MB (9 MB) |
| 30 kbp, 5% | 1.5 MB (15 MB) | 1.9 MB (7 MB) | **0.8 MB** (6 MB) | 1.4 MB (7 MB) | **0.8 MB** (7 MB) | 1.5 MB (8 MB) | 7.8 MB (42 MB) |
| 100 kbp, 5% | 3.6 MB (18 MB) | 5.6 MB (11 MB) | 2.5 MB (8 MB) | 4.3 MB (10 MB) | **1.5 MB** (8 MB) | 3.1 MB (10 MB) | 74 MB (100 MB) |
| 300 kbp, 5% | 11 MB (28 MB) | 16 MB (22 MB) | 8.0 MB (13 MB) | 13 MB (20 MB) | **2.3 MB** (9 MB) | 7.8 MB (17 MB) | 648 MB (673 MB) |
| 1 Mbp, 5% | **30 MB** (58 MB) | 53 MB (65 MB) | 58 MB (69 MB) | 42 MB (63 MB) | — | — | — |
| ***Uniform pairs, 15% divergence*** |  |  |  |  |  |  |  |
| 3 kbp, 15% | 0.6 MB (14 MB) | 0.4 MB (5 MB) | **0.2 MB** (5 MB) | 0.6 MB (6 MB) | 0.9 MB (7 MB) | 1.0 MB (8 MB) | 1.5 MB (8 MB) |
| 10 kbp, 15% | 0.7 MB (15 MB) | 0.8 MB (6 MB) | **0.4 MB** (6 MB) | 2.2 MB (8 MB) | 0.8 MB (7 MB) | 1.0 MB (8 MB) | 6.4 MB (39 MB) |
| 30 kbp, 15% | 1.3 MB (15 MB) | 1.8 MB (7 MB) | **0.9 MB** (6 MB) | 4.6 MB (10 MB) | **0.9 MB** (7 MB) | 1.9 MB (9 MB) | 48 MB (79 MB) |
| 100 kbp, 15% | 3.9 MB (18 MB) | 5.7 MB (11 MB) | 3.1 MB (8 MB) | 18 MB (24 MB) | **1.6 MB** (8 MB) | 3.6 MB (11 MB) | 511 MB (551 MB) |
| 300 kbp, 15% | **12 MB** (29 MB) | 21 MB (27 MB) | 19 MB (25 MB) | 38 MB (50 MB) | 2.4 MB (9 MB) † | 8.9 MB (≥ 20 MB) † | 4632 MB (≥ 4699 MB) † |
| 1 Mbp, 15% | **52 MB** (85 MB) | 134 MB (145 MB) | 86 MB (98 MB) | 149 MB (≥ 183 MB) † | — | — | — |
| *runtime alone (one 8 bp pair)* | *12 MB* | *3 MB* | *3 MB* | *3 MB* | *5 MB* | *5 MB* | *5 MB* |

- **The memory an alignment adds is what compares aligners**: each runner reads `getrusage` before and after each pair, outside its timer, so the runtime and the sample read in count for nothing. The whole peak counts them, 14 to 15 MB for dinara-align's runner on short pairs, most of it the Mojo runtime's 12 MB, against 5 to 8 MB for the Rust ones, so there the brackets show runtimes, not aligners. Each runner reads its sample in one call of the file's size; read through a growing buffer instead, the memory it left behind went to the alignments unseen, and dinara-align's growth on ont-500k measured 32 MB rather than about 50.
- **On pairs up to 30 kbp every aligner adds about 2 MB at most**, but WFA, keeping every front, and A\*PA; at 100 kbp Edlib adds the least, 1.5 to 1.6 MB, dinara-align 3.6 to 3.9 and A\*PA2 2.5 to 5.7.
- **On the long reads dinara-align adds the least but on ont-500k-genvar**: 50 MB on ont-500k against A\*PA2's 81 and 82, and on 1 Mbp pairs 30 MB at 5% against 42 to 58 and 52 MB at 15% against 86 to 149; on ont-500k-genvar 35 MB, between A\*PA2-simple's 32 and A\*PA2-full's 48. A\*PA2 reports the same shape on whole datasets: on reads over 500 kbp A\*PA2-full adds 30 MB in median and 82 at most.
- **WFA's grows with the square of the distance**, 511 MB at 100 kbp and 4.6 GB at 300 kbp at 15%, where its memory cap stops it.

### Beyond one thread

All of the above is one thread, as A\*PA2's evaluation measures, and one pair always runs on one thread. dinara-align also aligns batches of pairs across threads, one pair a thread, which the evaluation runs none of the others doing. On this machine's ten cores, from `pixi run bench-astarpa2`, a batch takes 3 µs a read on ont-1k, 24 µs on ont-10k and 48.4 ms on four ont-500k-genvar reads, against 31 µs, 248 µs and 225 ms for the fastest other aligner on one thread.
