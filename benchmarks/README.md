# Benchmarks

Times dinara-align against other exact global DNA aligners on the same inputs, and refuses to report a time for any workload where the tools disagree on the answer.

```bash
pixi run bench                    # 1 kbp and 10 kbp workloads, about a minute the first time, then seconds
pixi run bench --full             # adds the 100 kbp pairs, about four and a half minutes the first time
pixi run bench --install-rust     # also installs the nightly Rust A*PA needs, under .cache/
pixi run bench --remeasure        # times the pinned rivals again rather than replaying their kept rows
```

The table prints to the terminal and lands in `.cache/results/results.md`, with every raw row in `results.tsv`.

```bash
pixi run bench-astarpa2           # A*PA2's own evaluation datasets, about fifteen seconds warm
pixi run bench-astarpa2 --fresh   # the rivals too, rather than their kept results, about three minutes
```

That table lands in `.cache/results/astarpa2.md`; see [A\*PA2's Evaluation](#apa2s-evaluation).

## What Runs

Every rival is cloned at a pinned commit into `.cache/` and built there; nothing of theirs is vendored into this repository.
A pinned rival's times do not change, so each harness times every rival once on a machine and keeps its rows, and a later run replays them and times dinara-align alone, in seconds; a new pin, runner, workload or CPU, or `--remeasure`, times the rival again. Every answer is still checked against every other.

| Tool                                                                            | Commit    | License    | What it is run on                                                          |
| :------------------------------------------------------------------------------ | :-------- | :--------- | :------------------------------------------------------------------------- |
| dinara-align                                                                    | this tree | MPL-2.0    | the read batches and the affine pairs, on the CPU and on the GPU where one answers |
| dinara-align (bit-parallel, …)                                                  | this tree | MPL-2.0    | the edit-distance workloads: `distance` and `align` at unit costs, ported from A\*PA2-simple, on one thread, pair by pair and as a batch |
| [hyalite](https://github.com/Psy-Fer/hyalite)                                   | `0189bcb` | MIT        | the read batches and the affine pairs, global mode (`Mode::Nw`), one CPU thread with its NEON or AVX2 kernels |
| [A\*PA, A\*PA2](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) | `bf2e14e` | MPL-2.0    | the edit-distance workloads, because it computes edit distance and nothing else |
| [Edlib](https://github.com/Martinsos/edlib), [WFA2-lib](https://github.com/smarco/WFA2-lib) through [pa-bench](https://github.com/pairwise-alignment/pa-bench) | `af7a50d` | MIT | the edit-distance workloads: the other exact aligners A\*PA2's evaluation compares against, with its parameters |
| [Edlib](https://github.com/Martinsos/edlib) v1.2.7, [KSW2](https://github.com/lh3/ksw2), [WFA2-lib](https://github.com/smarco/WFA2-lib) v2.3.6, [parasail](https://github.com/jeffdaily/parasail), [SSW](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) | `ec2310e`, `289609b`, `bcf473a`, `fb985ee`, `a66636b` | MIT | the other modes: free ends, seed extension, two-piece gaps and substitution tables (see Other Modes) |

A\*PA runs four of its aligners: A\*PA2-full (its default), A\*PA2-simple, the original A\*PA, and `astarpa2_nw`, a bit-parallel full-matrix Needleman-Wunsch that is the closest like-for-like to a full sweep.
The three A\*PA2 aligners run twice, with traceback and without, so each has a score row and an alignment row; the original A\*PA has no cost-only entry point and reports its alignment alone.
`a*pa2-nw` against `dinara-align (bit-parallel, 1 thread)` is the same algorithm in Rust and in Mojo.
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

Measured with `pixi run bench --full`, which took about five minutes, on the Intel Core i9-7900X of A\*PA2's results below (ten cores at a fixed 3.3 GHz, turbo boost and hyper-threading off, 91 GB) with an NVIDIA GeForce RTX 2070, on the otherwise idle machine, every column on one core, dinara-align's CPU column as every rival's, since a library runs on its caller's thread (see Short-Read Batches for every core), every tool built for the machine's own instruction set, AVX-512 included (`--cpu native`, the default), dinara-align and its GPU column at `19e5f18`, its edit-distance column at `ebe290f`, every tool the fastest of three runs, the rivals measured again beside it (`--repeat 3`).
Hyalite built with Rust 1.92, A\*PA with the nightly its repository pins, both with `-C target-cpu=native`.
A dash marks a task the tool does not offer, or a workload it is not run on.

| workload | task | dinara-align (cpu) | dinara-align (gpu) | dinara-align (bit-parallel, 1 thread) | hyalite | a*pa2-full | a*pa2-simple | a*pa2-nw | a*pa | edlib | biwfa | wfa | agree |
| :-- | :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| reads-150bp | score | 5.26 ms | 1.31 ms | — | 284 ms | — | — | — | — | — | — | — | ✓ |
| reads-150bp | alignment | 16.4 ms | 14 ms | — | 1.61 s | — | — | — | — | — | — | — | ✓ |
| reads-1kbp | score | 33.5 ms | 6.69 ms | — | 1.07 s | — | — | — | — | — | — | — | ✓ |
| reads-1kbp | alignment | 250 ms | 14.2 ms | — | 12.8 s | — | — | — | — | — | — | — | ✓ |
| affine-1k | score | 81 µs | 888 µs | — | 1.08 ms | — | — | — | — | — | — | — | ✓ |
| affine-1k | alignment | 141 µs | 1.69 ms | — | 12.3 ms | — | — | — | — | — | — | — | ✓ |
| affine-10k | score | 2.54 ms | 5.55 ms | — | 159 ms | — | — | — | — | — | — | — | ✓ |
| affine-10k | alignment | 3.56 ms | 7.51 ms | — | 1.55 s | — | — | — | — | — | — | — | ✓ |
| affine-100k | score | 239 ms | 72.3 ms | — | 16.3 s | — | — | — | — | — | — | — | ✓ |
| affine-100k | alignment | 459 ms | 152 ms | — | 163 s | — | — | — | — | — | — | — | ✓ |
| edit-1k-1% | score | — | — | 2 µs | — | 43 µs | 17 µs | 19 µs | — | 25 µs | 3 µs | 3 µs | ✓ |
| edit-1k-1% | alignment | — | — | 2 µs | — | 45 µs | 17 µs | 205 µs | 306 µs | 91 µs | 6 µs | 4 µs | ✓ |
| edit-1k-5% | score | — | — | 4 µs | — | 43 µs | 17 µs | 19 µs | — | 24 µs | 7 µs | 7 µs | ✓ |
| edit-1k-5% | alignment | — | — | 7 µs | — | 48 µs | 19 µs | 202 µs | 288 µs | 94 µs | 16 µs | 10 µs | ✓ |
| edit-1k-15% | score | — | — | 12 µs | — | 33 µs | 17 µs | 19 µs | — | 63 µs | 35 µs | 35 µs | ✓ |
| edit-1k-15% | alignment | — | — | 31 µs | — | 64 µs | 45 µs | 202 µs | 876 µs | 143 µs | 66 µs | 47 µs | ✓ |
| edit-10k-1% | score | — | — | 14 µs | — | 789 µs | 777 µs | 1.41 ms | — | 373 µs | 35 µs | 35 µs | ✓ |
| edit-10k-1% | alignment | — | — | 22 µs | — | 824 µs | 146 µs | 17.6 ms | 2.81 ms | 1.8 ms | 80 µs | 52 µs | ✓ |
| edit-10k-5% | score | — | — | 136 µs | — | 585 µs | 1.37 ms | 1.39 ms | — | 586 µs | 373 µs | 374 µs | ✓ |
| edit-10k-5% | alignment | — | — | 276 µs | — | 663 µs | 281 µs | 17.3 ms | 2.55 ms | 2.41 ms | 796 µs | 571 µs | ✓ |
| edit-10k-15% | score | — | — | 167 µs | — | 854 µs | 2.36 ms | 1.41 ms | — | 1.65 ms | 2.54 ms | 2.54 ms | ✓ |
| edit-10k-15% | alignment | — | — | 354 µs | — | 1.25 ms | 886 µs | 17.8 ms | 79.4 ms | 4.38 ms | 5.25 ms | 3.67 ms | ✓ |
| edit-100k-1% | score | — | — | 734 µs | — | 5.34 ms | 151 ms | 133 ms | — | 8.23 ms | 1.56 ms | 1.56 ms | ✓ |
| edit-100k-1% | alignment | — | — | 1.07 ms | — | 5.54 ms | 3.09 ms | 1.82 s | 37.5 ms | 37.6 ms | 3.29 ms | 2.41 ms | ✓ |
| edit-100k-5% | score | — | — | 3.86 ms | — | 5.64 ms | 226 ms | 134 ms | — | 57.3 ms | 33.1 ms | 33.1 ms | ✓ |
| edit-100k-5% | alignment | — | — | 4.62 ms | — | 6.35 ms | 14.1 ms | 1.81 s | 33.5 ms | 127 ms | 67.2 ms | 51.4 ms | ✓ |
| edit-100k-15% | score | — | — | 8.25 ms | — | 19.7 ms | 200 ms | 134 ms | — | 87.8 ms | 245 ms | 245 ms | ✓ |
| edit-100k-15% | alignment | — | — | 14 ms | — | 23.1 ms | 22.9 ms | 1.81 s | 13.7 s | 249 ms | 493 ms | 681 ms | ✓ |

## Reading the Numbers

- **Every answer agreed**, so each row compares tools answering the same question. The answers are the optimal scores and distances: a sum and a position-weighted sum per workload, so a reordered batch cannot pass.
- **On edit distance, dinara-align on one thread beats or matches every exact aligner on every workload**, distance and alignment alike: 22 against WFA's 52 µs and A\*PA2-simple's 146 µs aligning 10 kbp at 1%, 276 against A\*PA2-simple's 281 µs at 10 kbp and 5%, a tie, 4.62 against A\*PA2-full's 6.35 ms at 100 kbp and 5%, 14 against A\*PA2-simple's 22.9 ms at 100 kbp and 15%, 734 µs against WFA's 1.56 ms scoring 100 kbp at 1%, and 12 against A\*PA2-simple's 17 µs scoring 1 kbp at 15%, where a short pair's diagonal transition gives way to the whole matrix once it would cost more. Its alignment follows a fixed rule for ties (see `Ties`), traced a tile at a time by a forward search, which costs about twice what its alignment of divergent long pairs took before that rule. The ties are the 1 kbp pairs at 1 and 5%, a few microseconds each, where BiWFA and WFA run the same diagonal transition.
- **Near-identical pairs run diagonal transition**, as WFA does, before any band: one front keeping its history for a close alignment, and two from both ends otherwise, as BiWFA scores, keeping both histories for an alignment and tracing it back through each from where they met. It stays on while it costs less than the band would, which at 1 to 3% divergence covers most pairs up to tens of kbp.
- **A band re-aims its first bound as it sweeps.** The bound starts from a projection of the first few edits, which strays by up to half either way. At an eighth, a quarter and half of the columns, the band projects the distance from its own climb, hundreds of edits in, lowering the bound when it was set too high and giving the round up early when it was set too low. Only a distance within the final bound is accepted, so the answer stays exact.
- **Long, moderately divergent pairs prune with A\*PA2-full's seed heuristic.** From about 1,500 projected edits up to one in seven bases, the first sequence is cut into 12-base seeds, their exact matches in the second are found by a rolling hash, and matches that cannot shorten any path over the next 14 seeds are dropped, as A\*PA's local pruning drops them. The band then keeps a row only while its score plus the seeds still ahead, less the longest chain of matches it can still reach, fits the bound. At 2 to 3% divergence that bound at the origin lands within 1% of the distance, so the band starts just past it and finishes in one round. A\*PA2-full also prunes matches between rounds, which dinara-align leaves out, as its rounds now rarely number more than two. On AVX-512, whose band runs fast beside the seeds' setup, seeds are used only from 86 kbp, where they start to pay there.
- **Long, divergent pairs match seeds within one edit**, as A\*PA's `r = 2` does. Past about one edit in fifteen bases most exact 12-base seeds are broken and the bound at the origin falls far short: 8,277 for a distance of 12,196 on a 100 kbp pair at 15%. Seeds of 16 bases that may match with one edit, charging two edits where they match nowhere, put it at 11,304. Each match is found by the half of it that matches exactly, a check on its quarters turns most chance lookups away, and an exact match's one-edit neighbours are left out, which keeps the bound a lower bound. They are used from 64 kbp, where the band they save outgrows their setup, when the projection says one edit in ten bases or under 40% of the exact seeds chain.
- **Global affine scores run a wavefront from both ends.** Under a table of one match and one mismatch score, as `Scoring.dna()` is, the match reward folds away for a global alignment (Eizenga and Lindquist), and WFA's three-layer wavefront then searches by cost, from the origin and from the corner at once until they meet: 239 ms against hyalite's 16.3 s at 100 kbp, 80 µs against 1.08 ms at 1 kbp. A pair whose projected wavefront would cost more than a full sweep is handed to the sweep, which runs sixteen cells at a time by anti-diagonal and also serves local scores.
- **Global affine alignments are traced back through the same wavefront's fronts**, five bytes a diagonal kept from each end and the path walked from where they met, a pair whose fronts would pass 80 MB split where an optimal path crosses: 141 µs against hyalite's 12.3 ms at 1 kbp, 3.56 ms against 1.55 s at 10 kbp and 460 ms against 163 s at 100 kbp; dinara-align's own GPU sweep is slower on the 1 and 10 kbp pairs, 1.69 and 7.5 ms, and faster on the 100 kbp one, 153 ms. A tie between optimal paths may resolve differently from the GPU's, the score the same. Under any other table the matrix is swept sixteen cells at a time by anti-diagonal, each pair's score read by comparison, by byte shuffle from a register for a table of up to sixteen entries, or gathered from memory for a larger one: a global alignment stores only the band of diagonals its score bounds, a local one is found by its span and aligned globally (see Other Modes), and past `max_memory`, 80 MB by default, a global one splits on rows in linear space, Myers and Miller's way.
- **Every CPU column is one core.** dinara-align, a library, runs on its caller's thread and starts no threads of its own: an application spreads its calls over its threads, as the command line does (see Short-Read Batches, on every core). The reads rows are each one call over the whole batch, its global scores and alignments many pairs at once, one a lane of a SIMD register over a band each pair's cost proves, each alignment traced from the flags its band kept: 5.25 ms scoring the 150 bp reads and 16.5 ms aligning them, against hyalite's 284 ms and 1.61 s.
- **dinara-align's GPU times include the overheads a caller pays**: opening the device context on every call, copying sequences over, and copying results back. That is why a single short pair is slower on the GPU than on the CPU, while the reads batches, every pair in one launch, several to a warp (see Short-Read Batches), run fastest there, the device against one core: 1.24 against 5.25 ms scoring the 150 bp reads, 6.03 against 33.5 ms the 1 kbp ones, and 14 against 249 ms aligning those. The 150 bp reads' scores go first over a band each pair's cost proves (see Short-Read Batches), 1.69 ms before it; the 1 kbp ones, too long for the whole matrices a warp holds and too divergent for the band to prove many, pay for its pass, 5.31 ms before it.
- **hyalite's traceback budget is 1 GiB.** It keeps the whole matrix when it fits and switches to a checkpointed sweep above that, as the 100 kbp pairs do. dinara-align switches to its linear-space recursion above six million cells on the host and one million on the device.
- **A\*PA answers only edit distance**, and its time depends on how similar the sequences are, as do the bit-parallel path's and dinara-align's global affine score and alignment; the local alignments and hyalite's columns do not.
- **One machine, one fixed pair per workload, one run.** The rows of microseconds move by a few percent between runs; read a difference of that size as a tie. Rerun on your own hardware before quoting any of this.

## A\*PA2's Evaluation

`pa_bench.py` runs the exact aligners of A\*PA2's evaluation, with its parameters, and WFA beside BiWFA, on its own datasets: Oxford Nanopore reads and SARS-CoV-2 genomes from pa-bench's release, and the uniform-error pairs pa-generate regenerates exactly from the evaluation's seed and lock.
Each dataset is a fixed shuffled sample of about 2 Mbp and at least four pairs, the same for every tool; each tool aligns its pairs with traceback, once each, as pa-bench times them, within a budget of five seconds, and every cost is checked against every other tool's and against the costs A\*PA2's published results recorded.
The rivals are pinned, so their results are kept and reused until `--fresh`; a second table gives each tool's peak resident memory.
dinara-align also runs as a batch, `alignments` over the whole sample in one call on the same one thread, whose column is the call's time over its pairs, a throughput where the others are each one pair's latency: at unit costs a batch takes the bit-parallel sweep pair by pair, as single calls do.

Measured on the i9-7900X of A\*PA2's results below, unpinned on the otherwise idle machine, every tool built for its own instruction set, dinara-align at `ebe290f`; a count marks a tool its budget, or WFA's 8 GiB cap, stopped partway, and more than the budget one that finished no pair:

| dataset | pairs | mean length | dinara-align (bit-parallel, 1 thread) | dinara-align (bit-parallel, batch, 1 thread) | a*pa2-full | a*pa2-simple | a*pa | edlib | biwfa | wfa | agree |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| ont-1k | 1221 of 12477 | 0.818 kbp | 23 µs | 22 µs | 81 µs | 48 µs | 612 µs | 126 µs | 44 µs | 31 µs | ✓ |
| ont-10k | 277 of 5000 | 3.6 kbp | 164 µs | 158 µs | 391 µs | 247 µs | 15.3 ms | 1.1 ms | 753 µs | 595 µs | ✓ |
| ont-50k | 104 of 10000 | 9.52 kbp | 586 µs | 564 µs | 1.27 ms | 944 µs | 218 ms (25/104) | 6.1 ms | 7.43 ms | 6.37 ms | ✓ |
| ont-500k | 4 of 50 | 638 kbp | 117 ms | 117 ms | 283 ms | 687 ms | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| ont-500k-genvar | 4 of 48 | 659 kbp | 156 ms | 156 ms | 225 ms | 511 ms | > 5 s | 5.18 s (1/4) | > 5 s | > 5 s | ✓ |
| sars-cov-2 | 33 of 10000 | 29.6 kbp | 283 µs | 255 µs | 2 ms | 728 µs | 6.49 ms | 8.15 ms | 895 µs | 983 µs | ✓ |
| Uniform-t10000000-n3000-e0.05 | 333 of 3333 | 3 kbp | 42 µs | 43 µs | 174 µs | 60 µs | 266 µs | 813 µs | 82 µs | 57 µs | ✓ |
| Uniform-t10000000-n10000-e0.05 | 99 of 1000 | 10 kbp | 285 µs | 284 µs | 788 µs | 292 µs | 937 µs | 3.8 ms | 557 µs | 440 µs | ✓ |
| Uniform-t10000000-n30000-e0.05 | 33 of 333 | 30 kbp | 927 µs | 913 µs | 2.45 ms | 1.71 ms | 3.1 ms | 14.4 ms | 4.98 ms | 4.27 ms | ✓ |
| Uniform-t10000000-n100000-e0.05 | 10 of 100 | 100 kbp | 4.87 ms | 4.76 ms | 8.22 ms | 15.5 ms | 12.1 ms | 145 ms | 54.5 ms | 45.8 ms | ✓ |
| Uniform-t10000000-n300000-e0.05 | 4 of 33 | 300 kbp | 20 ms | 19.3 ms | 27.6 ms | 65 ms | 36.9 ms | 749 ms | 471 ms | 450 ms | ✓ |
| Uniform-t10000000-n1000000-e0.05 | 4 of 10 | 1e+03 kbp | 63 ms | 63 ms | 125 ms | 861 ms | 152 ms | > 5 s | 5.39 s (1/4) | > 5 s | ✓ |
| Uniform-t10000000-n3000-e0.15 | 333 of 3333 | 3 kbp | 114 µs | 115 µs | 202 µs | 152 µs | 1.55 ms | 873 µs | 376 µs | 289 µs | ✓ |
| Uniform-t10000000-n10000-e0.15 | 99 of 1000 | 10 kbp | 405 µs | 408 µs | 914 µs | 785 µs | 6.47 ms | 5.73 ms | 3.98 ms | 2.94 ms | ✓ |
| Uniform-t10000000-n30000-e0.15 | 33 of 333 | 30 kbp | 1.7 ms | 1.65 ms | 4.07 ms | 2.62 ms | 23.4 ms | 24 ms | 34.1 ms | 26.7 ms | ✓ |
| Uniform-t10000000-n100000-e0.15 | 9 of 100 | 100 kbp | 11.5 ms | 11.4 ms | 16.3 ms | 24.4 ms | 96.4 ms | 260 ms | 375 ms | 307 ms | ✓ |
| Uniform-t10000000-n300000-e0.15 | 4 of 33 | 300 kbp | 46 ms | 45.6 ms | 125 ms | 293 ms | 464 ms | 2.53 s (2/4) | 3.38 s (1/4) | 4.65 s (1/4) | ✓ |
| Uniform-t10000000-n1000000-e0.15 | 4 of 10 | 1e+03 kbp | 266 ms | 265 ms | 1.33 s | 1.34 s | 2.03 s (2/4) | > 5 s | > 5 s | > 5 s | ✓ |


## Affine Costs

`pa_bench.py --affine x,o,e` runs the same samples at gap-affine costs as WFA counts them, a mismatch `x` and a gap of `k` letters `o + k e`, against the exact aligners that take them: WFA2-lib keeping every front (WFA), its lowest-memory mode (BiWFA), and KSW2's banded SSE kernel with band doubling, on x86-64 alone.
dinara-align answers with `align` under `Costs.affine`, a CIGAR as the others hand back, from a wavefront grown from both ends at once, which keeps five bytes a diagonal for its traceback and splits a pair whose fronts would pass 80 MB where an optimal path crosses, as BiWFA does.

```bash
pixi run bench-astarpa2 --affine 4,6,2
```

At WFA's costs (4, 6, 2) on the i9-7900X above, pinned to one core, every tool built for its own instruction set, dinara-align at `19e5f18`, each tool five seconds a sample; a count marks a tool its budget stopped partway, and more than the budget one that finished no pair:

| dataset | pairs | mean length | dinara-align | WFA | BiWFA | KSW2 | agree |
| :-- | --: | --: | --: | --: | --: | --: | :-: |
| ont-1k | 1221 of 12477 | 0.818 kbp | 112 µs | 186 µs | 304 µs | 1.1 ms | ✓ |
| ont-10k | 277 of 5000 | 3.6 kbp | 1.68 ms | 4.8 ms | 5.84 ms | 30.8 ms (163/277) | ✓ |
| ont-50k | 104 of 10000 | 9.52 kbp | 21.2 ms | 68.4 ms (82/104) | 57.3 ms (90/104) | 367 ms (16/104) | ✓ |
| ont-500k | 4 of 50 | 638 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| ont-500k-genvar | 4 of 48 | 659 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| sars-cov-2 | 33 of 10000 | 29.6 kbp | 755 µs | 2.06 ms | 1.79 ms | 88.4 ms | ✓ |
| Uniform-t10000000-n3000-e0.05 | 333 of 3333 | 3 kbp | 283 µs | 579 µs | 913 µs | 6.12 ms | ✓ |
| Uniform-t10000000-n10000-e0.05 | 99 of 1000 | 10 kbp | 2.17 ms | 6.36 ms | 8.02 ms | 97.1 ms (52/99) | ✓ |
| Uniform-t10000000-n30000-e0.05 | 33 of 333 | 30 kbp | 18.1 ms | 56 ms | 62.7 ms | 1.35 s (4/33) | ✓ |
| Uniform-t10000000-n100000-e0.05 | 10 of 100 | 100 kbp | 304 ms | 754 ms (7/10) | 643 ms (8/10) | > 5 s | ✓ |
| Uniform-t10000000-n300000-e0.05 | 4 of 33 | 300 kbp | 3.07 s (1/4) | > 5 s | > 5 s | > 5 s | ✓ |
| Uniform-t10000000-n1000000-e0.05 | 4 of 10 | 1e+03 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| Uniform-t10000000-n3000-e0.15 | 333 of 3333 | 3 kbp | 1.26 ms | 3.6 ms | 5.08 ms | 18.7 ms (268/333) | ✓ |
| Uniform-t10000000-n10000-e0.15 | 99 of 1000 | 10 kbp | 12.2 ms | 38.5 ms | 44.8 ms | 203 ms (25/99) | ✓ |
| Uniform-t10000000-n30000-e0.15 | 33 of 333 | 30 kbp | 151 ms | 385 ms (13/33) | 372 ms (14/33) | 2.48 s (2/33) | ✓ |
| Uniform-t10000000-n100000-e0.15 | 9 of 100 | 100 kbp | 2.02 s (2/9) | > 5 s | 4.19 s (1/9) | > 5 s | ✓ |
| Uniform-t10000000-n300000-e0.15 | 4 of 33 | 300 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |
| Uniform-t10000000-n1000000-e0.15 | 4 of 10 | 1e+03 kbp | > 5 s | > 5 s | > 5 s | > 5 s | ✓ |

- **dinara-align is the fastest on every dataset any tool finished,** 1.7 to 3.2 times faster than WFA: 112 against 186 µs on ont-1k, 1.68 against 4.8 ms on ont-10k, 755 µs against 2.06 ms on the SARS-CoV-2 genomes, and 12.2 against 38.5 ms on 10 kbp at 15%. Its CIGAR follows a fixed rule for ties (see `Ties`), its search grown on to the optimum and traced back from the far end. BiWFA beats WFA only on the genomes and where WFA's budget runs out.
- **It alone finishes every ont-50k read,** in 21.2 ms on average, where WFA finishes 82 of the 104 within the budget and BiWFA 90, and it finishes more of the 100 and 300 kbp pairs than anyone.
- **KSW2's band doubles from eight diagonals against a loose bound,** so it times out on most reads past 1 kbp.
- **No aligner finishes a 500 kbp read or a 1 Mbp pair in five seconds:** the wavefront's work grows with the square of the cost, and gap-affine costs have no seed heuristic here to prune it (see TODO.md).

Growth of peak memory aligning each pair, median / largest, as A\*PA2's Table 10 measures it:

| dataset | dinara-align | WFA | BiWFA | KSW2 |
| :-- | --: | --: | --: | --: |
| ont-1k | 0.0 / 1.1 MB | 0.0 / 2.9 MB | 0.0 / 1.5 MB | 0.0 / 0.8 MB |
| ont-10k | 0.0 / 22 MB | 0.0 / 97 MB | 0.0 / 2.4 MB | 0.0 / 58 MB |
| ont-50k | 0.0 / 85 MB | 0.0 / 1420 MB | 0.0 / 4.6 MB | 0.0 / 317 MB |
| sars-cov-2 | 0.0 / 11 MB | 0.0 / 51 MB | 0.0 / 2.5 MB | 0.0 / 211 MB |
| Uniform-t10000000-n3000-e0.05 | 0.0 / 1.7 MB | 0.0 / 5.4 MB | 0.0 / 1.6 MB | 0.0 / 3.1 MB |
| Uniform-t10000000-n10000-e0.05 | 0.0 / 7.2 MB | 0.0 / 29 MB | 0.0 / 1.8 MB | 0.0 / 40 MB |
| Uniform-t10000000-n30000-e0.05 | 0.0 / 46 MB | 0.0 / 212 MB | 0.0 / 3.0 MB | 31 / 455 MB |
| Uniform-t10000000-n100000-e0.05 | 0.0 / 90 MB | 12 / 2413 MB | 0.0 / 9.0 MB | — |
| Uniform-t10000000-n300000-e0.05 | 98 / 98 MB | — | — | — |
| Uniform-t10000000-n3000-e0.15 | 0.0 / 4.8 MB | 0.0 / 18 MB | 0.0 / 1.3 MB | 0.0 / 12 MB |
| Uniform-t10000000-n10000-e0.15 | 0.0 / 35 MB | 0.0 / 160 MB | 0.0 / 3.2 MB | 0.0 / 79 MB |
| Uniform-t10000000-n30000-e0.15 | 0.0 / 89 MB | 0.5 / 1417 MB | 0.0 / 5.6 MB | 795 / 795 MB |
| Uniform-t10000000-n100000-e0.15 | 97 / 97 MB | — | 16 / 16 MB | — |

dinara-align keeps under a tenth of WFA's memory on the largest pairs (90 against 2413 MB at 100 kbp), as its traceback keeps a byte of flags and the column of one front a diagonal, from fronts each half as long; BiWFA, keeping only its last few fronts, stays smallest.

## Local and Overlap Alignment

Local alignment, Smith-Waterman, overlap alignment, every end gap free, and a read placed whole in a window with a reward, each with its CIGAR, against the aligners that offer them:

```bash
pixi run bench-local   # builds and times SSW, parasail and abPOA the first time, a few minutes; then seconds
```

Every tool scores a match 2, a mismatch -4 and a gap of `k` letters `6 + 2k` (dinara-align's `Mode.local(2)` and `Mode.overlap(2)` under `Costs.affine(4, 6, 2)`), and every tool's scores must agree on every pair, or the run fails; none disagreed.
A tool's time is the faster of two passes over a workload, its mean per pair, on one thread of the Skylake-X, pinned (see A\*PA2's results below for the machine), every tool built for its own instruction set, AVX-512 included, dinara-align at `19e5f18`.
abPOA aligns to a graph, so its time includes adding the reference to one, as any pairwise use of it pays; SSW and abPOA have no overlap or infix mode.
The infix row scores with a reward, `Mode.INFIX.with_match_score(2)`, as parasail's `sg_dx` and hyalite's HW do; without one, dinara-align's `Mode.INFIX` minimizes the costs, Edlib's way, on the bit-parallel sweep at unit costs or the wavefront at any others.

| workload | dinara-align | SSW | parasail | abPOA | hyalite |
| :-- | --: | --: | --: | --: | --: |
| local, 600 short noisy pairs (20 to 700 bp) | **32 µs** | 52 µs | 66 µs | 173 µs | 843 µs |
| local, 1 kbp read at 10% in a 10 kbp window | **1.40 ms** | 3.19 ms | 4.70 ms | 16.7 ms | 160 ms |
| local, 10 kbp at 5% against 12 kbp | **50.5 ms** | 85.0 ms | 285 ms | 207 ms | 2.29 s |
| overlap, 2 kbp reads overlapping by 0.5 to 1.5 kbp | **2.15 ms** | — | 15.3 ms | — | 42.5 ms |
| infix with a reward, 1 kbp read at 10% placed whole in a 3 kbp window | **1.40 ms** | — | 14.8 ms | — | 45.4 ms |

- **dinara-align is the fastest on every workload**, 1.6 to 2.3 times SSW on local alignment, and 7 and 11 times parasail on overlaps and scored infixes.
- **Its sweep runs by anti-diagonal in 16-bit lanes along the shorter sequence,** while the scores fit, and finds only the best end; an extension back from that end, stopping once it earns the sweep's score, gives the alignment, traced as it searched. SSW sweeps striped (Farrar), whose lazy pass runs most of every column when a long alignment scores high, as on the 10 kbp pairs.

## Other Modes

Free ends, seed extension, two-piece gaps and substitution tables of more than one mismatch score, each against the aligners that offer the same mode, on DNA:

```bash
pixi run bench-modes   # builds and times Edlib, WFA2-lib, KSW2, parasail and SSW the first time; then seconds
```

Every tool aligns every pair with its CIGAR on one thread of the Skylake-X, pinned, every tool built for its own instruction set, AVX-512 included, dinara-align at `19e5f18`; a tool's time is the faster of two passes over a workload, its mean per pair, and every tool's costs or scores must agree on every pair, or the run fails; none disagreed.
WFA2-lib runs exact, keeping every front, its WF-adaptive heuristic off.
KSW2 extends with no Z-drop: it gauges one by anti-diagonal and dinara-align as WFA2-lib does, a cost at a time, so the two would stop at different places and their times compare different work.
Deletions priced apart from insertions are left out: no rival offers them.

| workload | dinara-align | Edlib | WFA2-lib | KSW2 | parasail | SSW |
| :-- | --: | --: | --: | --: | --: | --: |
| infix at unit costs, 1 kbp read at 10% placed whole in a 3 kbp window | **135 µs** | 270 µs | 986 µs | — | — | — |
| infix at unit costs, 10 kbp read at 5% in a 30 kbp window | **3.14 ms** | 13.4 ms | 53.7 ms | — | — | — |
| prefix at unit costs, 1 kbp read at 10% against 2 kbp of reference | **43 µs** | 114 µs | 47 µs | — | — | — |
| infix at WFA's (4, 6, 2), the 1 kbp reads above | **3.00 ms** | — | 4.10 ms | — | — | — |
| extension, a match 2 at (4, 6, 2), 300 to 1,400 bases at 5% then noise | **726 µs** | — | — | 2.07 ms | — | — |
| the same with an end bonus of 50, 0 to 80 bases of noise at the end | **158 µs** | — | — | 1.23 ms | — | — |
| two-piece gaps (4, 6, 2, 24, 1), 5 kbp at 5% with three indels of 100 to 400 bp | **8.44 ms** | — | 15.9 ms | 46.6 ms | — | — |
| a table, match 2, transition -2, transversion -4, global, 1 kbp at 10% | **542 µs** | — | — | — | 2.94 ms | — |
| the same table, local, 1 kbp read at 10% in a 10 kbp window | **2.45 ms** | — | — | — | 5.16 ms | 3.27 ms |

- **dinara-align is the fastest on every workload**: 2.0 and 4.3 times Edlib on infixes, 2.9 times KSW2 on extension and 7.8 times with an end bonus, 1.9 times WFA2-lib on two-piece gaps and 5.4 times parasail on a table globally; on an infix at WFA2-lib's own costs 3.00 against 4.10 ms. The prefix search runs 43 against WFA2-lib's 47 µs: it gives a doomed try up once its climb projects past the bound, as the global band does, where its first try used to sweep two thirds of the way.
- **An end bonus aligns a read to its end when that scores within the bonus of the best stop**, as KSW2's `end_bonus` and BWA-MEM's clipping penalty decide, and both tools choose the same alignment on every read. The extension's own search weighs it, keeping the best of its points on the read's last row as it grows and going on only while one could still win, so a read ending in a short stretch of noise costs less than a plain extension of a long one.
- **Some rows move by up to a third between commits that do not touch them.** The infix at WFA's costs takes 3.00 ms at `16732cf` and 3.92 at `2f8b935`, which changed only the prefix search: both run the same 15.29 billion instructions and 2.17 billion branches, but the slower takes 8.9 billion of its micro-ops from the legacy decoder against the faster's 4.8, its loop placed where the micro-op cache misses it, in 6.5 billion cycles against 5.0. The Skylake-X penalises a jump that crosses a 32-byte boundary (Intel's JCC erratum fix). `mojo build` passes no option to LLVM's assembler, and its documentation names none; assembling Mojo's own `--emit asm` output with LLVM's `-mbranches-within-32B-boundaries` and linking it as `mojo build` does runs this row in 2.95 ms at both commits, which shows the cause, but the builds here do not, as that would depend on Mojo's internal output and on the machine's own LLVM. Later x86 cores and Apple silicon are not affected. Read a difference of that size between commits as noise; the global table's 527 µs at `16732cf` and 630 at `2f8b935` are another, and so is the 11% the 150 bp reads' batch score lost at `348e3a4`, a commit that changed only the GPU's path.
- **The prefix search sweeps only the band its bound allows**, which Edlib's prefix mode does not: the band's top moves down with the diagonal, no column past the read's length plus the best end so far is swept, and a try whose band has already emptied stops. The part of the reference it finds is then aligned globally by the edit-distance traceback.
- **A table of more than one mismatch score sweeps sixteen cells at a time**, by anti-diagonal, each pair's score read by byte shuffle from a register when the table has at most sixteen entries, as DNA's four letters do, and gathered from memory otherwise. A local alignment is found as SSW finds one: a 16-bit sweep for its end, one anchored there for the first cell that earns its score, its start, and the span between aligned globally. A global one stores only the band of diagonals its score bounds. Both used to sweep a cell at a time: 13.3 ms globally and 133 ms locally.

## Short-Read Batches

Accelign's short-read case study (Kallenborn et al., BMC Bioinformatics 2026) scores millions of Illumina reads against the reference sections BWA placed them in, on every thread: each pair's global score alone, a match 0, a mismatch −1 and a gap of `k` letters `−(2 + k)`, dinara-align's `Costs.affine(1, 2, 1)`.
Its reads are not ours to ship, so the pairs here are drawn to the same shape from a fixed seed: 500,000 reads of 148 bp at 1% error against sections of their source, most the read's own 148 bp, a fifth with up to 22 bases fewer or 30 more at either end.

```bash
pixi run bench-batch   # builds the rivals the first time; then a minute or two
```

Every tool scores the whole batch twice: on one core, as a library runs on its caller's thread, and on all ten cores of the Skylake-X, each tool's calls spread by its own driver, the rivals' through OpenMP, one aligner a thread, and dinara-align's by the benchmark program, the batch cut before the clock starts into forty pieces, each one `distances` call, each thread taking the next piece as it finishes; the library itself starts no threads. A time is the faster of two passes; every tool's costs must agree, summed and position-weighted, or the run fails. WFA2-lib runs exact, its score alone, its WF-adaptive heuristic off. dinara-align on the GPU is `scores` on the RTX 2070, the faster of five calls once three have brought the device's clocks up, each the whole call: packing the batch, copying it over, scoring it and copying back. dinara-align and its GPU column at `19e5f18`, every rival as measured beside it at `16732cf` on one core and on ten.

| workload | threads | dinara-align | dinara-align GPU | WFA2-lib | KSW2 | parasail | Edlib |
| :-- | :-- | --: | --: | --: | --: | --: | --: |
| affine, `distances(..., Costs.affine(1, 2, 1))`, on the GPU `scores` at the same scores; KSW2's `extz2` and parasail's striped `nw`, scores alone | one core | **233 ms** | — | 1.06 s | 10.9 s | 21.7 s | — |
| unit costs, `distances(...)`; WFA2-lib exact, score alone; Edlib's NW distance | one core | **212 ms** | — | 450 ms | — | — | 2 s |
| affine, `distances(..., Costs.affine(1, 2, 1))`, on the GPU `scores` at the same scores; KSW2's `extz2` and parasail's striped `nw`, scores alone | 10 cores | **26 ms** | 32 ms | 115 ms | 1.18 s | 2.34 s | — |
| unit costs, `distances(...)`; WFA2-lib exact, score alone; Edlib's NW distance | 10 cores | **26 ms** | — | 48 ms | — | — | 223 ms |

- **On the CPU a batch of global distances runs many pairs at once**, one a 16-bit lane of a register, 32 under AVX-512, the inter-sequence vectorization SeqAn's and parasail's batch modes sweep whole matrices with, here over a band of diagonals each pair's own cost proves (`lanes`): a path visiting a diagonal past the end diagonal and the main one runs a gap out to it and one back, so it costs at least their two openings and extensions, and a cost no dearer than every path off the band is the least. The first band is the end and main diagonals and one more either side; the pairs it does not prove, a seventh here, are scored again over every diagonal a cheaper path could visit. Every pair goes first into a byte whose additions saturate, 64 to an AVX-512 register, as SSW and parasail do, a cost under 255 exact, and a pair whose byte saturates into 16 bits after. Their letters go side by side eight positions at a time, a lane's eight one load, the block turned in registers. One core scores a pair in 0.47 µs at the affine costs and 0.42 at unit costs, WFA2-lib in 2.12 and 0.90; the second pass's pairs go by the width of their bands, so a group's lanes need bands alike, and a group stops once every lane's row has passed what a byte or the cap holds, no row costing less than the one before.
- **On the GPU the batch takes 32 ms**, 44 before the band below, most of it the host packing the batch, ten threads at it, which no faster device would shorten. Each pair goes first over a band of sixteen diagonals, a thread a pair, its scores, deletions and insertions a register a diagonal, swept anti-diagonal by anti-diagonal, letters read two bits a base off the tape in the words they crossed the bus in; the band's score stands where the pair's own cost proves it, as on the CPU, and the thread gives up as soon as the cheapest cell of its anti-diagonal costs more than every path off the band, a path's cost only growing. Five pairs in six are proved; the rest go on their chunk's list and are scored whole on the device from the same tape, all while the host packs the next chunk: a warp scores several at once, a group of lanes to each, every lane holding a run of its pair's columns in registers and handing its last on by shuffle, the shape Accelign's paper describes, chosen for the batch from the compiled ones by the fewest cells it sweeps, not tuned for a device. The cells of a global alignment are held shifted by their anti-diagonal times the gap extension, which takes Gotoh's recurrence from five additions a cell to three, what Hopper's fused add-and-max instructions do in hardware, done by algebra on any device. Tried on the RTX 2070 and dropped as slower: wider bands, a group of lanes to a pair handing its edges on by shuffle, 33,000 noisy reads over 128 diagonals taking 10 ms; scores as floats; each column's substitutions held in a register; and row letters handed on by shuffle. Reads the band seldom proves cost a little more than whole matrices alone, packed for a band that then proves little: 150 bp reads at 15% errors 34 ms against 31.5, 1 kbp reads at 5% 253 against 236; at 2% the 250 bp ones take 29 ms against 64. Accelign's own case study ran on a Grace Hopper GH200 with 16-bit arithmetic and DPX instructions, which the RTX 2070 has neither of.
- **Alignments with free ends go into the lanes too**, the span first, as the rule for ties picks it one pair at a time: the end, of the cells an optimal alignment may stop on the one furthest along the reference, from a sweep proving its cost with every optimal path inside its band; then the start, by a sweep back from that end over the reversed prefixes; then the span's own global alignment, traced as above. On one core of the Skylake-X, 150 bp reads at 3% in 200 bp windows align under affine costs (4, 6, 2) in 8.0 µs a pair as an infix where one at a time took 26.6, in 12.1 against 29.3 as a prefix, and in 20.0 against 50.3 with 60 letters free at the window's end and the read's start; at unit costs that overlap takes 8.3 µs against 32.4, while an infix or a prefix at unit costs keeps the bit-parallel sweep, faster still.
- **Pairs the lanes do not take go one at a time**: a pair whose dearest path would not fit 16 bits, with free ends one with an empty side or whose alignment ends on the first row or column, and with free ends under a band. For those the two-ended wavefront and the bit-parallel sweep remain, each worker keeping its searches' memory from pair to pair, the forward search alone taking a cheap pair as far as its ring holds, and every worker's memory kept two cache lines clear of another's, since x86's L2 fetches lines in pairs.

## A\*PA2's Results, Redone

A\*PA2's own results section ([curiouscoding.nl/posts/astarpa2](https://curiouscoding.nl/posts/astarpa2/#results)), redone with dinara-align beside the exact aligners it compares, and WFA, on a machine set up as A\*PA2's was: its real datasets, its uniform pairs swept over divergence and over length, and in place of its ablation, dinara-align's history on the long reads.

```bash
pixi run results-astarpa2 collect   # forty minutes the first time; after that only dinara-align runs, a few minutes
pixi run results-astarpa2 tables    # the tables below, into .cache/results/astarpa2-results.md
```

### Setup

As A\*PA2's evaluation runs them: one single-threaded job at a time, every pair aligned once with its traceback, which every aligner hands back as a CIGAR, the time the average wall clock per alignment, reading the data left out.

The machine is an Intel Core i9-7900X under Ubuntu 26.04 (Linux 7.0), with Mojo 1.1.0, set up as A\*PA2's i7-10750H was: every core fixed at 3.3 GHz, turbo boost and hyper-threading off, and each collection pinned to one core with `taskset`.
One thing differs: the jobs ran at normal priority, not A\*PA2's niceness −20, on an otherwise idle machine. Every tool is built for this machine's own instruction set, AVX-512 included (`--cpu native`, the default), as A\*PA2's own repository builds it; dinara-align at `ebe290f`.

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
| ont-1k | 1221 | 0.818 kbp | 23 µs (23 µs) | 81 µs (89 µs) | 48 µs (54 µs) | 612 µs (556 µs) | 126 µs (120 µs) | 44 µs (46 µs) | 31 µs (32 µs) | 41 µs (42 µs), 93% optimal | 41 µs (47 µs), 85% optimal |
| ont-10k | 277 | 3.6 kbp | 163 µs (111 µs) | 391 µs (264 µs) | 247 µs (177 µs) | 15.3 ms (3.29 ms) | 1.1 ms (579 µs) | 753 µs (284 µs) | 595 µs (223 µs) | 385 µs (172 µs), 61% optimal | 198 µs (141 µs), 61% optimal |
| ont-50k | 104 | 9.52 kbp | 587 µs (242 µs) | 1.27 ms (513 µs) | 944 µs (307 µs) | 190 ms (10.5 ms), 81/104 | 6.1 ms (2.02 ms) | 7.43 ms (875 µs) | 6.37 ms (722 µs) | 2.32 ms (459 µs), 52% optimal | 657 µs (319 µs), 61% optimal |
| ont-500k | 15 | 638 kbp | 104 ms (83.8 ms) | 176 ms (113 ms) | 537 ms (371 ms) | — | 4.99 s (3.04 s), 4/15 | 19.8 s (19.8 s), 1/15 | — | 989 ms (629 ms), 53% optimal | 728 ms (714 ms) |
| ont-500k-genvar | 15 | 651 kbp | 144 ms (114 ms) | 210 ms (177 ms) | 624 ms (506 ms) | — | 5.23 s (5.2 s), 4/15 | 6.36 s (7.57 s), 3/15 | — | 533 ms (273 ms), 7% optimal | 836 ms (647 ms), 0% optimal of 1 |
| sars-cov-2 | 33 | 29.6 kbp | 279 µs (188 µs) | 2 ms (2.07 ms) | 728 µs (662 µs) | 6.49 ms (2.07 ms) | 8.15 ms (7.72 ms) | 895 µs (446 µs) | 983 µs (359 µs) | 611 µs (353 µs), 97% optimal | 2.51 ms (2.3 ms), 30% optimal |

- **Among the exact aligners, dinara-align has the lowest mean and median on every dataset.** On the short reads, the case A\*PA2's post leaves to BiWFA, 23 against WFA's 31 µs and BiWFA's 44, median 23 against 32, since a short pair no longer hands a band, which covers its whole matrix there, what the diagonal transition does faster; on ont-10k and ont-50k it takes about two thirds of A\*PA2-simple's time, 163 against 247 µs and 587 against 944 µs, and on the SARS-CoV-2 genomes 279 µs against A\*PA2-simple's 728 and WFA's 983.
- **On the 500 kbp reads it takes three fifths to seven tenths of A\*PA2-full's time:** 104 against 176 ms on ont-500k and 144 against 210 ms on ont-500k-genvar, medians 83.8 against 113 and 114 against 177.
- **Edlib and BiWFA finish only some of the long reads within 20 seconds, BiWFA three genvar reads and one ont-500k read, and A\*PA and WFA none;** their means cover the reads they finished.
- **The approximate aligners trade optimality for speed, and are still slower than dinara-align on every dataset.** The closest, Block Aligner, takes 198 and 657 µs on ont-10k and ont-50k against dinara-align's 163 and 587, aligning 61% of those reads optimally; WFA-adaptive aligns 7% of the genvar reads optimally. A\*PA2 reports the same shares: Block Aligner 85% on ont-1k, WFA-adaptive 93%.

### Divergence, 100 kbp pairs

Mean time per alignment, A\*PA with whichever of `r = 1` and `r = 2` was the faster at each divergence:

| divergence | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| 0% | 75 µs | 6.58 ms | 1.54 ms | 6.58 ms | 41.7 ms | 630 µs | 330 µs |
| 1% | 1.33 ms | 7.46 ms | 3.38 ms | 7.75 ms | 57 ms | 3.3 ms | 3.38 ms |
| 2% | 4.41 ms | 7.52 ms | 4.84 ms | 8.37 ms | 71.3 ms | 10.2 ms | 9.72 ms |
| 3% | 4.68 ms | 7.71 ms | 8.37 ms | 8.08 ms | 96 ms | 21.1 ms | 19.1 ms |
| 4% | 4.75 ms | 7.92 ms | 7.39 ms | 9.74 ms | 98.9 ms | 35.8 ms | 31.4 ms |
| 5% | 4.89 ms | 8.11 ms | 15.5 ms | 12.1 ms | 145 ms | 54.1 ms | 46.9 ms |
| 6% | 6.3 ms | 8.76 ms | 14.1 ms | 17.5 ms | 147 ms | 74 ms | 64.1 ms |
| 7% | 6.53 ms | 8.8 ms | 13.2 ms | 34.8 ms | 150 ms | 97.4 ms | 84 ms |
| 8% | 6.79 ms | 11.7 ms | 12.3 ms | 40.4 ms | 153 ms | 124 ms | 107 ms |
| 9% | 10.4 ms | 11.9 ms | 11.5 ms | 40.5 ms | 159 ms | 153 ms | 129 ms |
| 10% | 11 ms | 10.5 ms | 28.4 ms | 44.4 ms | 244 ms | 186 ms | 156 ms |
| 11% | 11.5 ms | 13.6 ms | 27.2 ms | 47.7 ms | 247 ms | 220 ms | 184 ms |
| 12% | 11.3 ms | 11.8 ms | 26.4 ms | 53.4 ms | 249 ms | 254 ms | 212 ms |
| 13% | 11.2 ms | 17.9 ms | 25.5 ms | 60.2 ms | 253 ms | 292 ms | 239 ms |
| 14% | 12.3 ms | 16.7 ms | 25 ms | 71 ms | 256 ms | 332 ms | 273 ms |
| 15% | 11.8 ms | 15.6 ms | 24.4 ms | 99.1 ms | 259 ms | 374 ms | 304 ms |

- **dinara-align is the fastest at every divergence but 10%,** where A\*PA2-full is: 10.5 against 11 ms.
- **Near-identical pairs finish in the diagonal transition:** 75 µs against WFA's 330 at 0%, 1.33 against BiWFA's 3.3 ms at 1%.
- **From 2% to 8% the exact seeds' bound lands close to the distance:** dinara-align runs in 4.4 to 6.8 ms, the next aligner, A\*PA2-simple at 2% and A\*PA2-full beyond, in 4.8 to 11.7.
- **From 9% on, where fewer than a fifth of the exact seeds chain,** dinara-align rebuilds them to match within one edit and stays between 10.4 and 12.3 ms, while A\*PA2-full runs from 10.5 at 10% to 16 to 18 ms by 13 to 15%; at 15% A\*PA2-simple takes 24.4 ms, A\*PA 99, Edlib 259, WFA 304 and BiWFA 374. These seeds take an exact match's neighbours a base shorter and longer, without which the bound could overestimate (`11b24d2`), from the exact match once it is kept, with no scan or search of their own (`ebe290f`).

### Length

Mean time per alignment at 5% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| 3 kbp | 42 µs | 174 µs | 60 µs | 266 µs | 813 µs | 82 µs | 57 µs |
| 10 kbp | 284 µs | 788 µs | 292 µs | 937 µs | 3.8 ms | 557 µs | 440 µs |
| 30 kbp | 941 µs | 2.45 ms | 1.71 ms | 3.1 ms | 14.4 ms | 4.98 ms | 4.27 ms |
| 100 kbp | 4.85 ms | 8.22 ms | 15.5 ms | 12.1 ms | 145 ms | 54.5 ms | 45.8 ms |
| 300 kbp | 19.4 ms | 27.6 ms | 65 ms | 36.9 ms | 749 ms | 471 ms | 450 ms |
| 1 Mbp | 63.3 ms | 125 ms | 861 ms | 152 ms | 8.49 s (2/4) | 5.39 s (3/4) | 5.21 s |

At 15% divergence:

| length | dinara-align | A*PA2-full | A*PA2-simple | A*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| 3 kbp | 114 µs | 202 µs | 152 µs | 1.55 ms | 873 µs | 376 µs | 289 µs |
| 10 kbp | 405 µs | 914 µs | 785 µs | 6.47 ms | 5.73 ms | 3.98 ms | 2.94 ms |
| 30 kbp | 1.68 ms | 4.07 ms | 2.62 ms | 23.4 ms | 24 ms | 34.1 ms | 26.7 ms |
| 100 kbp | 11.5 ms | 16.3 ms | 24.4 ms | 96.4 ms | 260 ms | 375 ms | 307 ms |
| 300 kbp | 46 ms | 125 ms | 293 ms | 464 ms | 2.52 s | 3.37 s | 3.04 s |
| 1 Mbp | 265 ms | 1.33 s | 1.34 s | 2.07 s | 18.8 s (1/4) | — | — |

- **dinara-align is the fastest at every length and both divergences.** At 5% it aligns a 1 Mbp pair in 63.3 ms against A\*PA2-full's 125 and A\*PA's 152; at 15% in 265 ms against A\*PA2-full's 1.33 s and A\*PA's 2.07 s. Its narrowest leads are at 10 kbp, 284 against A\*PA2-simple's 292 µs at 5%, and at 100 kbp and 15%, 11.5 against A\*PA2-full's 16.3 ms.
- **At 15% its band still grows faster than the length,** but less than A\*PA2-full's: from 300 kbp to 1 Mbp, 3.3 times the length, dinara-align takes 5.8 times as long and A\*PA2-full 10.6 times, where A\*PA, whose pruning makes its heuristic nearly exact, takes 4.5 times.

### dinara-align's history on the long reads

A\*PA2 measures each of its methods by adding them one at a time. dinara-align's methods came in commits that also changed other things, so here instead are the same fifteen reads of each 500 kbp set aligned by the package as of the commit that added each method, a 60-second budget each. The commits before the x86 build fix take that fix's one line, the inline assembly's register constraint, and nothing else:

| commit | ont-500k | ont-500k-genvar |
| :-- | --: | --: |
| `495017d` port | 1.32 s (481 ms) | 1.19 s (855 ms) |
| `1907e33` + diagonal transition | 1.19 s (845 ms) | 8.26 s (9.83 s), 7/15 |
| `3f548c4` + seed heuristic | 1.04 s (714 ms) | 8.57 s (10.2 s), 7/15 |
| `8800593` + local pruning | 764 ms (168 ms) | 8.12 s (10.2 s), 7/15 |
| `16af9c8` + two-ended search | 472 ms (161 ms) | 4.33 s (5.4 s), 13/15 |
| `95d14ed` + real-read fixes | 456 ms (160 ms) | 754 ms (494 ms) |
| `2b038c8` + retries aimed | 163 ms (94.5 ms) | 200 ms (156 ms) |
| `117018a` + inexact seeds | 179 ms (127 ms) | 240 ms (200 ms) |
| today | 104 ms (83.8 ms) | 144 ms (114 ms) |

- **On ont-500k the time falls from 1.32 s at the port of A\*PA2-simple's band doubling to 163 ms once the retries are aimed.**
- **On ont-500k-genvar the commits between the port and `95d14ed` were several times slower,** 4.3 to 8.5 s a read and up to eight reads past the budget: real reads gather their errors at the ends, and the projections those commits aimed by ran many times over the distance, so they swept bands far wider than needed. `95d14ed` trusted a projection only when two of them agreed and otherwise grew the band from what was known, bringing genvar back to 753 ms, and aiming the retries to 200 ms.
- **Inexact seeds made these reads slower here at first,** 163 to 179 ms and 200 to 240 ms: this machine's band, on AVX-512, runs fast beside its seeds' setup, so inexact seeds paid only on more divergent pairs; the cutoff now follows the vector width, 20% of the exact seeds chained on AVX-512 against 40% on narrower vectors.
- **A leaner seeds' setup since** (`7435edd`: a branch-free one-edit test, only the seeds whose matches can reach the end, a leaner local pruning) brings both below their best before inexact seeds, and with retries aimed closer, AVX-512 sweeping two groups at a time (`475684c`) and exact seeds' matching filtered (`ff26cea`) since, they took 100 ms on ont-500k and 136 ms on ont-500k-genvar, a CIGAR included (`9bf5295`), and 112 and 165 ms once the inexact seeds kept an exact match's neighbours, which the bound needs to stay a lower bound (`11b24d2`), and take 104 and 144 ms since an exact match brings its neighbours itself (`ebe290f`); the diagonal transition's slides by AVX-512 gathers (`45697e0`) took the short reads and low-divergence pairs down 7 to 19% since. Every CIGAR now follows a fixed rule for ties (`f95d819`), the band's retraced a tile at a time by a forward search (`5a0e514`).

### Memory

From `pixi run bench-astarpa2` on the same machine, over its samples of each dataset (the 1221 reads of ont-1k, four pairs of each long set). Each cell is the most memory aligning one pair added to the most the process had held before it, the measure of A\*PA2's Table 10, with the process's whole peak resident memory in brackets; the lowest of the aligners that finished every pair is bold, and † marks one that finished only some within its budget, a dash none:

| dataset | dinara-align | A\*PA2-full | A\*PA2-simple | A\*PA | Edlib | BiWFA | WFA |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| ***Real datasets*** |  |  |  |  |  |  |  |
| ont-1k | 0.8 MB (16 MB) | 0.2 MB (5 MB) | 0.4 MB (5 MB) | 0.2 MB (6 MB) | 0.5 MB (8 MB) | **0.1 MB** (7 MB) | 0.4 MB (8 MB) |
| ont-10k | 1.0 MB (17 MB) | 0.6 MB (6 MB) | **0.5 MB** (6 MB) | 20 MB (40 MB) | **0.5 MB** (8 MB) | 1.0 MB (8 MB) | 6.1 MB (39 MB) |
| ont-50k | 1.9 MB (21 MB) | 1.5 MB (8 MB) | **0.7 MB** (7 MB) | 83 MB (—) † | 1.0 MB (8 MB) | 1.4 MB (9 MB) | 71 MB (187 MB) |
| ont-500k | **51 MB** (69 MB) | 81 MB (89 MB) | 82 MB (90 MB) | — | — | — | — |
| ont-500k-genvar | 37 MB (65 MB) | 48 MB (56 MB) | **32 MB** (55 MB) | — | 4.7 MB (14 MB) † | — | — |
| sars-cov-2 | 2.2 MB (19 MB) | 1.9 MB (7 MB) | **0.6 MB** (6 MB) | 13 MB (19 MB) | 1.4 MB (8 MB) | 1.6 MB (10 MB) | 15 MB (26 MB) |
| ***Uniform pairs, 5% divergence*** |  |  |  |  |  |  |  |
| 3 kbp, 5% | 1.4 MB (16 MB) | 0.4 MB (5 MB) | **0.2 MB** (5 MB) | 0.3 MB (5 MB) | 0.9 MB (7 MB) | 0.4 MB (7 MB) | 0.8 MB (7 MB) |
| 10 kbp, 5% | 1.9 MB (16 MB) | 0.6 MB (6 MB) | **0.4 MB** (5 MB) | 0.5 MB (6 MB) | 0.9 MB (7 MB) | 1.1 MB (8 MB) | 2.0 MB (9 MB) |
| 30 kbp, 5% | 1.6 MB (16 MB) | 1.9 MB (7 MB) | **0.8 MB** (6 MB) | 1.4 MB (7 MB) | **0.8 MB** (7 MB) | 1.5 MB (8 MB) | 7.8 MB (42 MB) |
| 100 kbp, 5% | 4.2 MB (19 MB) | 5.6 MB (11 MB) | 2.5 MB (8 MB) | 4.3 MB (10 MB) | **1.5 MB** (8 MB) | 3.1 MB (10 MB) | 74 MB (100 MB) |
| 300 kbp, 5% | 13 MB (32 MB) | 16 MB (22 MB) | 8.0 MB (13 MB) | 13 MB (20 MB) | **2.3 MB** (9 MB) | 7.8 MB (17 MB) | 648 MB (673 MB) |
| 1 Mbp, 5% | **31 MB** (65 MB) | 53 MB (65 MB) | 58 MB (69 MB) | 42 MB (63 MB) | — | — | — |
| ***Uniform pairs, 15% divergence*** |  |  |  |  |  |  |  |
| 3 kbp, 15% | 1.1 MB (15 MB) | 0.4 MB (5 MB) | **0.2 MB** (5 MB) | 0.6 MB (6 MB) | 0.9 MB (7 MB) | 1.0 MB (8 MB) | 1.5 MB (8 MB) |
| 10 kbp, 15% | 1.1 MB (15 MB) | 0.8 MB (6 MB) | **0.4 MB** (6 MB) | 2.2 MB (8 MB) | 0.8 MB (7 MB) | 1.0 MB (8 MB) | 6.4 MB (39 MB) |
| 30 kbp, 15% | 1.6 MB (16 MB) | 1.8 MB (7 MB) | **0.9 MB** (6 MB) | 4.6 MB (10 MB) | **0.9 MB** (7 MB) | 1.9 MB (9 MB) | 48 MB (79 MB) |
| 100 kbp, 15% | 5.1 MB (19 MB) | 5.7 MB (11 MB) | 3.1 MB (8 MB) | 18 MB (24 MB) | **1.6 MB** (8 MB) | 3.6 MB (11 MB) | 511 MB (551 MB) |
| 300 kbp, 15% | **13 MB** (32 MB) | 21 MB (27 MB) | 19 MB (25 MB) | 38 MB (50 MB) | 2.4 MB (9 MB) † | 8.9 MB (≥ 20 MB) † | 4632 MB (≥ 4699 MB) † |
| 1 Mbp, 15% | **55 MB** (86 MB) | 134 MB (145 MB) | 86 MB (98 MB) | 149 MB (≥ 183 MB) † | — | — | — |
| *runtime alone (one 8 bp pair)* | *13 MB* | *3 MB* | *3 MB* | *3 MB* | *5 MB* | *5 MB* | *5 MB* |

- **The memory an alignment adds is what compares aligners**: each runner reads `getrusage` before and after each pair, outside its timer, so the runtime and the sample read in count for nothing. The whole peak counts them, 15 to 16 MB for dinara-align's runner on short pairs, most of it the Mojo runtime's 13 MB, against 5 to 8 MB for the Rust ones, so there the brackets show runtimes, not aligners. Each runner reads its sample in one call of the file's size; read through a growing buffer instead, the memory it left behind went to the alignments unseen, and dinara-align's growth on ont-500k measured 32 MB rather than about 50.
- **On pairs up to 30 kbp every aligner adds about 2 MB at most**, but WFA, keeping every front, and A\*PA; at 100 kbp Edlib adds the least, 1.5 to 1.6 MB, dinara-align 4.2 to 5.1 and A\*PA2 2.5 to 5.7.
- **On the long reads dinara-align adds the least but on ont-500k-genvar**: 51 MB on ont-500k against A\*PA2's 81 and 82, and on 1 Mbp pairs 31 MB at 5% against 42 to 58 and 55 MB at 15% against 86 to 149; on ont-500k-genvar 37 MB, between A\*PA2-simple's 32 and A\*PA2-full's 48. A\*PA2 reports the same shape on whole datasets: on reads over 500 kbp A\*PA2-full adds 30 MB in median and 82 at most.
- **WFA's grows with the square of the distance**, 511 MB at 100 kbp and 4.6 GB at 300 kbp at 15%, where its memory cap stops it.

### Beyond one thread

All of the above is one thread, as A\*PA2's evaluation measures, and one pair always runs on one thread. dinara-align, a library, runs every call on its caller's thread and keeps no state between calls, so an application spreads its pairs over its own threads, as `dinara-align`'s command line does over every thread it may use; Short-Read Batches times that spread on all ten cores, 26 ms for 500,000 reads at the affine costs where WFA2-lib takes 115.
