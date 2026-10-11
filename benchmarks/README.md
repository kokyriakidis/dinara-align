# Benchmarks

Times dinara-align against the exact aligners of the two latest evaluations, A\*PA2's (Groot Koerkamp,
Bioinformatics 2025) and QuickEd's (Doblas et al., Bioinformatics 2025), on their datasets, and fails
any run in which two tools disagree on an answer.

```bash
pixi run bench                       # unit costs, every dataset, 5 s a tool a dataset
pixi run bench --affine 4,6,2        # gap-affine costs as WFA counts them: mismatch, opening, extension
pixi run bench --only ont l10000.    # the datasets whose names hold one of these
pixi run bench --budget 30           # the budget the tables below ran with
pixi run bench --fresh               # the rivals too, rather than their kept results
pixi run bench-modes                 # local, overlap, infix, prefix, extension, two-piece, tables
pixi run bench-batch                 # 500,000 short reads, on one core and on every core
pixi run bench-ablation              # dinara-align with each technique switched off in turn
```

Each writes its tables to the terminal and to `.cache/results/`: `bench.md` and `bench-affine-4-6-2.md`
with every table below and more (the totals over each whole dataset among them), `bench-scaling.json`
for the paper's charts, `modes.md` and `batch-results.md`.

## What Runs

Every rival is cloned at a pinned commit into `.cache/` and built there for the host's instruction set;
nothing of theirs is vendored. A pinned rival's results do not change, so each runs on a sample once and
its rows are kept, keyed by its pin, its runner's source and the CPU; a rerun times dinara-align alone,
in seconds to minutes, and still checks every answer against every other.

| tool | commit | run on |
| :-- | :-- | :-- |
| dinara-align | this tree | everything, pair by pair, and as a batch, every pair of a sample in one call on the same thread |
| [A\*PA, A\*PA2-simple, A\*PA2-full](https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner) | `bf2e14e` | unit costs, with the settings A\*PA2's evaluation gave each kind of data |
| [Edlib, BiWFA, WFA](https://github.com/pairwise-alignment/pa-bench) through pa-bench | `af7a50d` | unit costs, as A\*PA2's evaluation configured them |
| [QuickEd](https://github.com/maxdoblas/QuickEd) | `2b1d44b` | unit costs, its bound-and-align (x86-64) |
| [BSAlign](https://github.com/ruanjue/bsalign) | `1707aab` | unit and affine costs, its striped kernels with the band off, so exact (x86-64) |
| [WFA2-lib](https://github.com/smarco/WFA2-lib) v2.3.6 | `bcf473a` | affine costs (WFA, BiWFA), the other modes, the short reads |
| [KSW2](https://github.com/lh3/ksw2) | `289609b` | unit and affine costs with band doubling, extension, two-piece gaps, the short reads (x86-64) |
| [Edlib](https://github.com/Martinsos/edlib) v1.2.7 | `ec2310e` | infix and prefix at unit costs, the short reads |
| [parasail](https://github.com/jeffdaily/parasail) | `fb985ee` | local, overlap and scored infix alignment, substitution tables, the short reads |
| [SSW](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library) | `a66636b` | local alignment |

**The main benchmark** (`bench.py`) samples each dataset: its pairs in a fixed random order, as many as
fit 2 Mbp of sequence and at least four, the same sample for every tool. Each tool aligns the sample's
pairs with their traceback, once each, on one thread, until its budget runs out; a pair still running a
second past it is stopped and the pairs before it count. One run gives both evaluations' views: A\*PA2's
mean time per alignment and memory growth per pair, and QuickEd's total over the whole dataset, the mean
times its pairs, and peak resident memory. Every cost is checked against every other tool's on the pairs
both aligned, and on A\*PA2's datasets against the costs its published results recorded. WFA, KSW2 and
BSAlign, whose memory grows with the matrix or the square of the distance, run under an 8 GB cap; *t.o.*
marks a tool its budget stopped before a pair, *o.o.m.* one the cap stopped, as QuickEd's tables do.

**The other suites** time many short calls, where one cold call measures the allocator as much as the
aligner: each tool's time is the faster of two passes over a workload, its mean per pair, in process.

All tables were measured on an Intel Core i9-7900X (Skylake-X, AVX-512), ten cores at a fixed 3.3 GHz,
turbo boost and hyper-threading off, 91 GB, Ubuntu 26.04, on the otherwise idle machine, every tool on one
core pinned with `taskset` unless a row says otherwise, dinara-align at `dc980d7`, the affine and short-read tables at `289bcf0`, which nothing since touches; Rust nightly
2026-08-14, GCC 15.2, Mojo 1.1.0.

## Datasets

A\*PA2's real datasets, as pa-bench releases them: Oxford Nanopore reads of about 1, 10, 50 and 500 kbp,
the 500 kbp ones also against a reference with genetic variation (QuickEd's MiniION-a and MiniION-b), and
SARS-CoV-2 genomes. QuickEd's simulated ones, from Zenodo: random sequences of 10 kbp, 100 kbp, 500 kbp and
1 Mbp with 1, 5, 10 and 20% uniformly spread errors, 100 Mbp a file, downloaded the first time (3.2 GB);
they cover the lengths and divergences of A\*PA2's synthetic sets. Over each sample, as QuickEd's Table S1:

| dataset | pairs | length, min / mean / max (bp) | divergence, min / mean / max (%) |
| :-- | --: | --: | --: |
| ont-1k | 1221 of 12477 | 42 / 830 / 1,089 | 0 / 11 / 21 |
| ont-10k | 277 of 5000 | 342 / 3,643 / 9,791 | 4 / 13 / 20 |
| ont-50k | 104 of 10000 | 360 / 9,633 / 42,626 | 4 / 12 / 19 |
| ont-500k | 4 of 50 | 553,834 / 642,465 / 683,520 | 4 / 7 / 13 |
| ont-500k-genvar | 4 of 48 | 594,436 / 662,980 / 861,994 | 5 / 7 / 10 |
| sars-cov-2 | 33 of 10000 | 29,244 / 29,746 / 29,864 | 0 / 1 / 7 |
| l10000.n10000.e1 | 100 of 10000 | 10,000 / 10,003 / 10,023 | 1 / 1 / 1 |
| l10000.n10000.e5 | 99 of 10000 | 10,000 / 10,008 / 10,045 | 5 / 5 / 5 |
| l10000.n10000.e10 | 100 of 10000 | 10,000 / 10,008 / 10,056 | 9 / 9 / 9 |
| l10000.n10000.e20 | 99 of 10000 | 10,000 / 10,014 / 10,095 | 17 / 17 / 17 |
| l100000.n1000.e1 | 9 of 1000 | 100,000 / 100,013 / 100,035 | 1 / 1 / 1 |
| l100000.n1000.e5 | 9 of 1000 | 100,000 / 100,028 / 100,123 | 5 / 5 / 5 |
| l100000.n1000.e10 | 10 of 1000 | 100,000 / 100,039 / 100,120 | 9 / 9 / 9 |
| l100000.n1000.e20 | 10 of 1000 | 100,000 / 100,040 / 100,155 | 17 / 17 / 17 |
| l500000.n200.e1 | 4 of 200 | 500,000 / 500,038 / 500,103 | 1 / 1 / 1 |
| l500000.n200.e5 | 4 of 200 | 500,016 / 500,088 / 500,178 | 5 / 5 / 5 |
| l500000.n200.e10 | 4 of 200 | 500,000 / 500,052 / 500,105 | 9 / 9 / 9 |
| l500000.n200.e20 | 4 of 200 | 500,000 / 500,055 / 500,220 | 17 / 17 / 17 |
| l1000000.n100.e1 | 4 of 100 | 1,000,000 / 1,000,034 / 1,000,107 | 1 / 1 / 1 |
| l1000000.n100.e5 | 4 of 100 | 1,000,185 / 1,000,301 / 1,000,551 | 5 / 5 / 5 |
| l1000000.n100.e10 | 4 of 100 | 1,000,000 / 1,000,054 / 1,000,120 | 9 / 9 / 9 |
| l1000000.n100.e20 | 4 of 100 | 1,000,000 / 1,000,166 / 1,000,355 | 17 / 17 / 17 |

## Unit Costs

Edit distance with its alignment, a 30 s budget a tool a dataset. Mean time per alignment:

| dataset | pairs | dinara-align (bit-parallel, 1 thread) | dinara-align (bit-parallel, batch, 1 thread) | a*pa2-full | a*pa2-simple | a*pa | edlib | biwfa | wfa | ksw2 | quicked | bsalign | agree |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | :-: |
| ont-1k | 1221 of 12477 | 23 µs | 22 µs | 81 µs | 48 µs | 612 µs | 124 µs | 44 µs | 31 µs | 539 µs | 75 µs | 119 µs | ✓ |
| ont-10k | 277 of 5000 | 157 µs | 151 µs | 391 µs | 247 µs | 15.3 ms | 1.09 ms | 741 µs | 599 µs | 13.1 ms | 688 µs | 1.3 ms | ✓ |
| ont-50k | 104 of 10000 | 584 µs | 565 µs | 1.27 ms | 944 µs | 263 ms | 6.14 ms | 7.43 ms | 6.38 ms | 167 ms | 2.31 ms | 12.3 ms | ✓ |
| ont-500k | 4 of 50 | 119 ms | 117 ms | 283 ms | 687 ms | t.o. (30 s) | 5 s | 6.33 s | o.o.m. | o.o.m. | 847 ms | t.o. (30 s) | ✓ |
| ont-500k-genvar | 4 of 48 | 156 ms | 155 ms | 225 ms | 511 ms | t.o. (30 s) | 5.23 s | 6.2 s | o.o.m. | o.o.m. | 1.39 s | t.o. (30 s) | ✓ |
| sars-cov-2 | 33 of 10000 | 279 µs | 260 µs | 2 ms | 728 µs | 6.49 ms | 8.18 ms | 909 µs | 987 µs | 68.5 ms | 3.9 ms | 47.1 ms | ✓ |
| l10000.n10000.e1 | 100 of 10000 | 33 µs | 34 µs | 729 µs | 160 µs | 564 µs | 3.22 ms | 92 µs | 66 µs | 6.31 ms | 611 µs | 5.55 ms | ✓ |
| l10000.n10000.e5 | 99 of 10000 | 269 µs | 267 µs | 797 µs | 288 µs | 1.08 ms | 3.79 ms | 642 µs | 508 µs | 21.4 ms | 1.25 ms | 5.59 ms | ✓ |
| l10000.n10000.e10 | 100 of 10000 | 355 µs | 356 µs | 726 µs | 450 µs | 3.26 ms | 4.56 ms | 2.26 ms | 1.81 ms | 40.1 ms | 1.99 ms | 5.61 ms | ✓ |
| l10000.n10000.e20 | 99 of 10000 | 523 µs | 515 µs | 1.44 ms | 1.02 ms | 205 ms | 5.93 ms | 7.33 ms | 5.32 ms | 96.7 ms | 1.55 ms | 5.64 ms | ✓ |
| l100000.n1000.e1 | 9 of 1000 | 1.49 ms | 1.4 ms | 7.53 ms | 3.26 ms | 7.76 ms | 57.6 ms | 3.71 ms | 4.05 ms | 609 ms | 12.5 ms | 641 ms | ✓ |
| l100000.n1000.e5 | 9 of 1000 | 4.99 ms | 4.93 ms | 8.45 ms | 14.9 ms | 14.4 ms | 146 ms | 63.3 ms | 55 ms | 4.72 s (6/9) | 24 ms | 643 ms | ✓ |
| l100000.n1000.e10 | 10 of 1000 | 7.22 ms | 7.17 ms | 13.9 ms | 27.6 ms | 48.3 ms | 247 ms | 218 ms | 182 ms | 9.19 s (3/10) | 39 ms | 631 ms | ✓ |
| l100000.n1000.e20 | 10 of 1000 | 15.8 ms | 15.8 ms | 36.7 ms | 56.8 ms | t.o. (30 s) | 446 ms | 695 ms | 551 ms | 17.3 s (1/10) | 62.8 ms | 630 ms | ✓ |
| l500000.n200.e1 | 4 of 200 | 32.8 ms | 32.4 ms | 42.7 ms | 72.3 ms | 41.2 ms | 670 ms | 80.8 ms | 74.8 ms | 24.1 s (1/4) | 124 ms | t.o. (30 s) | ✓ |
| l500000.n200.e5 | 4 of 200 | 31.3 ms | 30.6 ms | 45 ms | 214 ms | 83.1 ms | 2.28 s | 1.56 s | 1.48 s | o.o.m. | 401 ms | t.o. (30 s) | ✓ |
| l500000.n200.e10 | 4 of 200 | 67.5 ms | 66.9 ms | 119 ms | 425 ms | 406 ms | 4.38 s | 5.45 s | 5.05 s | o.o.m. | 711 ms | t.o. (30 s) | ✓ |
| l500000.n200.e20 | 4 of 200 | 253 ms | 251 ms | 545 ms | 901 ms | t.o. (30 s) | 8.52 s (3/4) | 18 s (1/4) | o.o.m. | o.o.m. | 1.24 s | t.o. (30 s) | ✓ |
| l1000000.n100.e1 | 4 of 100 | 59.2 ms | 56.9 ms | 89.1 ms | 253 ms | 97.3 ms | 2.38 s | 366 ms | 391 ms | o.o.m. | 398 ms | t.o. (30 s) | ✓ |
| l1000000.n100.e5 | 4 of 100 | 63.9 ms | 64.4 ms | 118 ms | 818 ms | 183 ms | 8.67 s (3/4) | 6.35 s | o.o.m. | o.o.m. | 1.44 s | t.o. (30 s) | ✓ |
| l1000000.n100.e10 | 4 of 100 | 183 ms | 183 ms | 369 ms | 1.65 s | 978 ms | 17.3 s (1/4) | 22.2 s (1/4) | o.o.m. | o.o.m. | 2.65 s | t.o. (30 s) | ✓ |
| l1000000.n100.e20 | 4 of 100 | 897 ms | 897 ms | 2.02 s | 3.52 s | t.o. (30 s) | t.o. (30 s) | t.o. (30 s) | o.o.m. | o.o.m. | 4.8 s | t.o. (30 s) | ✓ |

Each rival's time over dinara-align's on the pairs both aligned, across the datasets, as QuickEd's text
summarizes its own:

| rival | least | geometric mean | most | datasets |
| :-- | --: | --: | --: | --: |
| a*pa2-full | 1.3× | 2.5× | 22.2× | 22 |
| a*pa2-simple | 1.1× | 3.3× | 12.8× | 22 |
| a*pa | 1.3× | 10.8× | 450.6× | 17 |
| edlib | 5.4× | 29.0× | 132.1× | 21 |
| biwfa | 1.9× | 13.8× | 110.7× | 21 |
| wfa | 1.4× | 7.2× | 74.8× | 16 |
| ksw2 | 23.5× | 246.8× | 1179.1× | 13 |
| quicked | 3.0× | 6.8× | 22.6× | 22 |
| bsalign | 5.2× | 40.2× | 429.4× | 12 |

Growth of peak memory aligning a pair, median and largest, as A\*PA2's Table 10 measures it:

| dataset | dinara-align (bit-parallel, 1 thread) | dinara-align (bit-parallel, batch, 1 thread) | a*pa2-full | a*pa2-simple | a*pa | edlib | biwfa | wfa | ksw2 | quicked | bsalign |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| ont-1k | 0.0 / 0.5 MB | — | 0.0 / 0.3 MB | 0.0 / 0.3 MB | 0.0 / 0.2 MB | 0.0 / 0.5 MB | 0.0 / 0.3 MB | 0.0 / 0.5 MB | 0.0 / 0.4 MB | 0.0 / 0.1 MB | 0.0 / 0.2 MB |
| ont-10k | 0.0 / 0.9 MB | — | 0.0 / 0.7 MB | 0.0 / 0.4 MB | 0.0 / 20 MB | 0.0 / 0.5 MB | 0.0 / 1.0 MB | 0.0 / 6.2 MB | 0.0 / 29 MB | 0.0 / 2.5 MB | 0.0 / 13 MB |
| ont-50k | 0.0 / 2.0 MB | — | 0.0 / 1.3 MB | 0.0 / 0.6 MB | 0.0 / 410 MB | 0.0 / 0.8 MB | 0.0 / 1.6 MB | 0.0 / 71 MB | 0.0 / 167 MB | 0.0 / 16 MB | 0.0 / 163 MB |
| ont-500k | 0.0 / 55 MB | — | 0.0 / 82 MB | 0.0 / 82 MB | — | 0.3 / 5.3 MB | 4.5 / 15 MB | — | — | 0.0 / 22 MB | — |
| ont-500k-genvar | 8.7 / 37 MB | — | 0.0 / 48 MB | 15 / 32 MB | — | 2.2 / 4.8 MB | 6.0 / 12 MB | — | — | 3.0 / 20 MB | — |
| sars-cov-2 | 0.0 / 2.2 MB | — | 0.0 / 1.8 MB | 0.0 / 0.8 MB | 0.0 / 13 MB | 0.0 / 1.2 MB | 0.0 / 1.6 MB | 0.0 / 15 MB | 0.0 / 97 MB | 0.0 / 7.9 MB | 0.0 / 211 MB |
| l10000.n10000.e1 | 0.0 / 0.6 MB | — | 0.0 / 0.9 MB | 0.0 / 0.5 MB | 0.0 / 0.6 MB | 0.0 / 0.6 MB | 0.0 / 0.5 MB | 0.0 / 0.7 MB | 0.0 / 3.2 MB | 0.0 / 0.5 MB | 0.0 / 24 MB |
| l10000.n10000.e5 | 0.0 / 1.4 MB | — | 0.0 / 0.7 MB | 0.0 / 0.4 MB | 0.0 / 0.6 MB | 0.0 / 0.7 MB | 0.0 / 1.4 MB | 0.0 / 2.1 MB | 0.0 / 11 MB | 0.0 / 1.5 MB | 0.0 / 24 MB |
| l10000.n10000.e10 | 0.0 / 1.3 MB | — | 0.0 / 0.8 MB | 0.0 / 0.5 MB | 0.0 / 2.4 MB | 0.0 / 0.6 MB | 0.0 / 1.4 MB | 0.0 / 4.3 MB | 0.0 / 20 MB | 0.0 / 2.6 MB | 0.0 / 24 MB |
| l10000.n10000.e20 | 0.0 / 0.9 MB | — | 0.0 / 0.8 MB | 0.0 / 0.4 MB | 0.0 / 28 MB | 0.0 / 0.6 MB | 0.0 / 1.4 MB | 0.0 / 11 MB | 0.0 / 40 MB | 0.0 / 4.5 MB | 0.0 / 24 MB |
| l100000.n1000.e1 | 0.0 / 5.3 MB | — | 0.0 / 5.9 MB | 0.0 / 1.9 MB | 0.0 / 4.1 MB | 0.0 / 1.4 MB | 0.0 / 4.0 MB | 3.4 / 5.9 MB | 0.0 / 203 MB | 0.0 / 12 MB | 0.0 / 2388 MB |
| l100000.n1000.e5 | 0.0 / 3.9 MB | — | 0.0 / 5.5 MB | 0.0 / 2.6 MB | 0.0 / 4.5 MB | 0.0 / 1.4 MB | 0.0 / 2.9 MB | 3.8 / 88 MB | 0.1 / 1555 MB | 0.0 / 11 MB | 0.0 / 2392 MB |
| l100000.n1000.e10 | 0.1 / 4.5 MB | — | 0.0 / 5.5 MB | 0.0 / 3.3 MB | 0.0 / 18 MB | 0.0 / 1.4 MB | 0.0 / 3.4 MB | 5.1 / 301 MB | 22 / 2986 MB | 0.0 / 17 MB | 0.0 / 2386 MB |
| l100000.n1000.e20 | 0.1 / 5.9 MB | — | 0.0 / 7.3 MB | 0.0 / 4.7 MB | — | 0.0 / 1.4 MB | 0.0 / 3.6 MB | 2.1 / 955 MB | 5476 / 5476 MB | 0.0 / 8.8 MB | 0.0 / 2392 MB |
| l500000.n200.e1 | 6.8 / 24 MB | — | 1.0 / 28 MB | 0.1 / 11 MB | 0.0 / 22 MB | 0.7 / 3.4 MB | 3.1 / 13 MB | 8.9 / 100 MB | 7835 / 7835 MB | 0.0 / 18 MB | — |
| l500000.n200.e5 | 4.5 / 16 MB | — | 0.7 / 26 MB | 0.2 / 19 MB | 2.8 / 23 MB | 0.0 / 4.0 MB | 2.8 / 12 MB | 22 / 2111 MB | — | 0.0 / 19 MB | — |
| l500000.n200.e10 | 1.1 / 21 MB | — | 0.2 / 30 MB | 0.1 / 30 MB | 5.4 / 76 MB | 0.1 / 4.0 MB | 2.8 / 14 MB | 123 / 7416 MB | — | 0.4 / 13 MB | — |
| l500000.n200.e20 | 18 / 53 MB | — | 0.1 / 59 MB | 0.3 / 53 MB | — | 0.1 / 4.1 MB | 15 / 15 MB | — | — | 0.6 / 16 MB | — |
| l1000000.n100.e1 | 6.9 / 40 MB | — | 2.2 / 57 MB | 0.2 / 28 MB | 0.0 / 43 MB | 0.0 / 7.6 MB | 6.8 / 17 MB | 8.8 / 382 MB | — | 0.0 / 25 MB | — |
| l1000000.n100.e5 | 4.4 / 32 MB | — | 2.2 / 53 MB | 0.2 / 57 MB | 5.5 / 46 MB | 0.1 / 7.9 MB | 6.7 / 16 MB | — | — | 0.9 / 25 MB | — |
| l1000000.n100.e10 | 4.9 / 42 MB | — | 0.6 / 67 MB | 0.6 / 100 MB | 13 / 151 MB | 7.9 / 7.9 MB | 19 / 19 MB | — | — | 0.9 / 24 MB | — |
| l1000000.n100.e20 | 20 / 125 MB | — | 0.4 / 184 MB | 2.9 / 192 MB | — | — | — | — | — | 0.1 / 26 MB | — |

- **dinara-align is the fastest on every one of the 22 datasets**, 1.1 to 2.6 times the next fastest exact
  aligner: 23 against WFA's 31 µs on ont-1k, 584 against A\*PA2-simple's 944 µs on ont-50k, 119 against
  A\*PA2-full's 283 ms on the 500 kbp reads, 279 against 728 µs on the SARS-CoV-2 genomes, and 897 ms
  against 2.02 s on the 1 Mbp pairs at 20%. The closest is QuickEd's 10 kbp pairs at 5%, 269 against
  A\*PA2-simple's 288 µs. Against QuickEd itself, the latest of them, 3.0 to 22.6 times, 6.8 in the
  geometric mean.
- **Near-identical pairs run diagonal transition**, as WFA does, before any band: one front keeping its
  history for a close pair, two from both ends otherwise, each kept for the traceback. It stays on while
  it costs less than the band would, which at 1 to 3% covers most pairs up to tens of kbp.
- **A band re-aims its first bound as it sweeps.** The bound starts from a projection of the first edits;
  at an eighth, a quarter and half of the columns the band projects the distance from its own climb,
  lowering a bound set too high and giving up early on one set too low. Only a distance within the final
  bound is accepted, so the answer stays exact.
- **Long pairs prune with A\*PA2's gap-chaining seed heuristic**, 12-base seeds matched exactly, or past
  about one edit in fifteen bases 16-base seeds matched within one edit, as A\*PA's `r = 2`, which keep
  the bound at the origin close to the distance on divergent pairs where most exact seeds break.
- **Memory stays low and bounded.** Up to 50 kbp an alignment adds about 2 MB at most. The band keeps its
  tiles' edges for the traceback within `max_memory`, 80 MB by default, and past it keeps some and
  sweeps the stretches between them again, the same CIGAR: on the 1 Mbp pairs at 20% it adds at most
  125 MB, where it added 167 before the cap, against A\*PA2's 184 and 192; the rest of that is the
  profile and the seeds, which grow with the length. QuickEd, splitting Hirschberg's way past 16 MB,
  stays near 25 MB.

## Affine Costs

Gap-affine alignment at WFA's costs (4, 6, 2), a mismatch 4 and a gap of `k` letters `6 + 2k`, against the
exact aligners that take them, a 30 s budget a tool a dataset:

| dataset | pairs | dinara-align (1 thread) | dinara-align (batch, 1 thread) | wfa | biwfa | ksw2 | bsalign | agree |
| :-- | --: | --: | --: | --: | --: | --: | --: | :-: |
| ont-1k | 1221 of 12477 | 112 µs | 105 µs | 186 µs | 304 µs | 1.1 ms | 336 µs | ✓ |
| ont-10k | 277 of 5000 | 1.71 ms | 1.61 ms | 4.8 ms | 5.84 ms | 29.8 ms | 7.84 ms | ✓ |
| ont-50k | 104 of 10000 | 22 ms | 21.7 ms | 63.5 ms | 54.9 ms | 343 ms (81/104) | 87 ms | ✓ |
| ont-500k | 4 of 50 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| ont-500k-genvar | 4 of 48 | 21.7 s (1/4) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| sars-cov-2 | 33 of 10000 | 737 µs | 679 µs | 2.06 ms | 1.79 ms | 88.4 ms | 341 ms | ✓ |
| l10000.n10000.e1 | 100 of 10000 | 208 µs | 219 µs | 461 µs | 650 µs | 21.3 ms | 37 ms | ✓ |
| l10000.n10000.e5 | 99 of 10000 | 2.32 ms | 2.33 ms | 7.05 ms | 8.72 ms | 97.5 ms | 37 ms | ✓ |
| l10000.n10000.e10 | 100 of 10000 | 7.18 ms | 7.07 ms | 22.7 ms | 27.3 ms | 204 ms | 37.3 ms | ✓ |
| l10000.n10000.e20 | 99 of 10000 | 20.6 ms | 20.2 ms | 65.2 ms | 75.2 ms | 374 ms (81/99) | 37.4 ms | ✓ |
| l100000.n1000.e1 | 9 of 1000 | 14 ms | 12.1 ms | 39.7 ms | 39 ms | 2.44 s | o.o.m. | ✓ |
| l100000.n1000.e5 | 9 of 1000 | 337 ms | 332 ms | 811 ms | 715 ms | 9.25 s (3/9) | o.o.m. | ✓ |
| l100000.n1000.e10 | 10 of 1000 | 1.19 s | 1.16 s | o.o.m. | 2.43 s | 17.5 s (1/10) | o.o.m. | ✓ |
| l100000.n1000.e20 | 10 of 1000 | 3.62 s (8/10) | t.o. (30 s) | o.o.m. | 7.34 s (4/10) | o.o.m. | o.o.m. | ✓ |
| l500000.n200.e1 | 4 of 200 | 490 ms | 477 ms | 1.29 s | 942 ms | o.o.m. | o.o.m. | ✓ |
| l500000.n200.e5 | 4 of 200 | 9.92 s (3/4) | t.o. (30 s) | o.o.m. | 19.4 s (1/4) | o.o.m. | o.o.m. | ✓ |
| l500000.n200.e10 | 4 of 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| l500000.n200.e20 | 4 of 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| l1000000.n100.e1 | 4 of 100 | 2.32 s | 2.31 s | o.o.m. | 3.92 s | o.o.m. | o.o.m. | ✓ |
| l1000000.n100.e5 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| l1000000.n100.e10 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| l1000000.n100.e20 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |

| rival | least | geometric mean | most | datasets |
| :-- | --: | --: | --: | --: |
| wfa | 1.7× | 2.7× | 3.2× | 11 |
| biwfa | 1.7× | 2.6× | 3.8× | 15 |
| ksw2 | 9.8× | 33.2× | 174.4× | 11 |
| bsalign | 1.8× | 12.7× | 462.3× | 8 |

- **dinara-align is the fastest on every dataset any tool finished**, 1.7 to 3.8 times WFA and BiWFA, 2.7
  and 2.6 in the geometric mean: 112 against WFA's 186 µs on ont-1k, 1.71 against 4.8 ms on ont-10k, 22.0
  against BiWFA's 54.9 ms on ont-50k, 737 µs against 1.79 ms on the SARS-CoV-2 genomes, and 20.6 against
  BSAlign's 37.4 ms on 10 kbp pairs at 20%, where BSAlign's 8-bit difference sweep, whose time does not
  depend on the divergence, comes closest.
- **It finishes more of the long pairs than anyone**: 8 of the 10 100 kbp pairs at 20% where BiWFA
  finishes 4, three 500 kbp pairs at 5% to BiWFA's one, and one ont-500k-genvar read. No tool finishes a
  1 Mbp pair past 1% in 30 s: the wavefront's work grows with the square of the cost, and gap-affine costs
  have no seed heuristic here to prune it.
- **A pair far apart is handed to a full sweep** in difference recurrences (Suzuki and Kasahara), 8-bit
  lanes at any length, once the wavefront projects it would finish later, the same alignment byte for
  byte; the wavefront gives a pair up only once its work has reached a 24th of the sweep's, as ONT reads,
  their errors gathered at the ends, otherwise looked far dearer from their first edits than they were.
- **Memory**: the wavefront keeps five bytes a diagonal of its fronts for the traceback, within 80 MB, and
  a pair whose fronts would pass that is split where an optimal path crosses, as BiWFA does. The rest of
  its growth, 112 MB on the 100 kbp pairs at 20% and 144 MB on a 500 kbp read, is the searches' rings,
  which grow with the cost rather than its square: a 100 kbp pair at 20% grows 115, 73, 42 and 33 MB under
  `max_memory` of 80, 40, 10 and 2 MB, its fronts filling the cap and its rings about 31 MB beside them.
  BiWFA, keeping only its last few fronts, stays smallest.

## Other Modes

Local, overlap and infix alignment, free ends, seed extension, two-piece gaps and substitution tables,
each against the aligners that offer the same mode (`bench-modes`). Every tool aligns every pair with its
CIGAR, and every tool's costs or scores must agree on every pair, or the run fails.

| workload | pairs | dinara-align | rivals |
| :-- | :-- | :-- | :-- |
| infix-edit | 100 reads of 1 kbp at 10%, each placed whole in a 3 kbp window | `Costs.edit()`, `Mode.INFIX` | Edlib's HW, WFA2-lib's ends-free edit |
| infix-edit-long | 20 reads of 10 kbp at 5% in 30 kbp windows | the same | the same |
| prefix-edit | 100 reads of 1 kbp at 10% against 2 kbp of reference | `Costs.edit()`, `Mode.PREFIX` | Edlib's SHW, WFA2-lib's ends-free edit |
| infix-affine | the infix-edit pairs | `Costs.affine(4, 6, 2)`, `Mode.INFIX` | WFA2-lib's ends-free gap-affine |
| extension | 200 reads following 2 kbp of reference for 300 to 1,400 bases at 5%, then noise | `Costs.affine(4, 6, 2)`, `Mode.extension(2)` | KSW2's `extz2`, no Z-drop |
| extension-bonus | the same with 0 to 80 bases of noise at the end, an end bonus of 50 | `end_bonus=50` | KSW2's `extz2` with an end bonus of 50 |
| two-piece | 50 pairs of 5 kbp at 5% with three indels of 100 to 400 bp | `Costs.two_piece(4, 6, 2, 24, 1)` | KSW2's `extd2`, WFA2-lib's two-piece |
| table-global | 100 pairs of 1 kbp at 10% | a `Scoring` with match 2, transition −2, transversion −4 | parasail |
| table-local | 20 reads of 1 kbp at 10% in 10 kbp windows | the same, `Mode.local()` | parasail, SSW |
| local-short | 600 short noisy pairs with planted cores of 20 to 300 bp | `Costs.affine(4, 6, 2)`, `Mode.local(2)` | SSW, parasail's `sw` |
| local-window | 20 reads of 1 kbp at 10% in 10 kbp windows | the same | the same |
| local-long | 5 pairs, 10 kbp at 5% against 12 kbp | the same | the same |
| overlap-reads | 100 pairs of 2 kbp reads overlapping by 0.5 to 1.5 kbp, both at 5% | `Costs.affine(4, 6, 2)`, `Mode.overlap(2)` | parasail's `sg` |
| infix-reads | 100 reads of 1 kbp at 10% placed whole in 3 kbp windows | `Mode.INFIX.with_match_score(2)` | parasail's `sg_dx` |

| workload | dinara-align | Edlib | WFA2-lib | KSW2 | parasail | SSW |
| :-- | --: | --: | --: | --: | --: | --: |
| infix-edit | **132 µs** | 273 µs | 1.01 ms | — | — | — |
| infix-edit-long | **2.01 ms** | 13.6 ms | 54 ms | — | — | — |
| prefix-edit | **44 µs** | 115 µs | 47 µs | — | — | — |
| infix-affine | **3.93 ms** | — | 3.94 ms | — | — | — |
| extension | **734 µs** | — | — | 2.06 ms | — | — |
| extension-bonus | **166 µs** | — | — | 1.23 ms | — | — |
| two-piece | **8.3 ms** | — | 15.9 ms | 45.3 ms | — | — |
| table-global | **198 µs** | — | — | — | 2.92 ms | — |
| table-local | **1.91 ms** | — | — | — | 5.63 ms | 3.27 ms |
| local-short | **32 µs** | — | — | — | 66 µs | 51 µs |
| local-window | **1.36 ms** | — | — | — | 4.73 ms | 3.15 ms |
| local-long | **17.5 ms** | — | — | — | 292 ms | 85.1 ms |
| overlap-reads | **724 µs** | — | — | — | 15.4 ms | — |
| infix-reads | **693 µs** | — | — | — | 15.2 ms | — |

- **dinara-align is the fastest on every workload**: twice Edlib on short infixes and 6.8 times
  on long ones, 2.8 times KSW2 on extension and 7.4 times with an end bonus, 1.9 times WFA2-lib on two-piece
  gaps, 14.7 times parasail on a table globally, and 4.9 times SSW on the long local pairs.
- **An end bonus aligns a read to its end when that scores within the bonus of the best stop**, as KSW2's
  `end_bonus` and BWA-MEM's clipping penalty decide, and both choose the same alignment on every read.
- **Local alignment sweeps by anti-diagonal in 16-bit lanes** while the scores fit, a strip of rows at a
  time so the rows stay in the first-level cache, and finds only the best end; an extension back from that
  end, stopping once it earns the sweep's score, gives the alignment. SSW sweeps striped (Farrar), whose
  lazy pass runs most of every column when a long alignment scores high, as on the 10 kbp pairs.
- **A table of more than one mismatch score sweeps sixteen cells at a time**, each pair's score read by
  byte shuffle from a register for up to sixteen entries, gathered from memory otherwise; a global
  alignment is traced in bands that prove themselves, a local one found by its end and start and the span
  between aligned globally.
- **Some rows move by up to a third between builds that do not touch them**, where a loop lands against
  the Skylake-X's micro-op cache (Intel's JCC erratum fix). The infix at WFA2-lib's costs, 3.93 ms here,
  takes 3.01 or 3.92 ms from build to build, the same 6.13 billion instructions either way and 2.0 or 2.6
  billion cycles; read a difference of that size between builds as noise.

## Short-Read Batches

Accelign's short-read case study (Kallenborn et al., BMC Bioinformatics 2026) scores millions of Illumina
reads against the reference sections BWA placed them in: each pair's global score alone, a mismatch 1 and a
gap of `k` letters `2 + k`, `Costs.affine(1, 2, 1)`, and here unit costs too. Its reads are not ours to
ship, so the pairs are drawn to the same shape from a fixed seed: 500,000 reads of 148 bp at 1% error
against sections of their source, most the read's own 148 bp, a fifth with up to 22 bases fewer or 30 more
at either end (`bench-batch`).

Every tool scores the whole batch on one core, as a library runs on its caller's thread, and on all ten,
each tool's calls spread by its own driver: the rivals' through OpenMP, one aligner a thread, and
dinara-align's by the benchmark program, the batch cut into forty pieces, each one `distances` call, each
thread taking the next piece as it finishes; the library itself starts no threads. WFA2-lib runs exact,
its score alone. dinara-align on the GPU is `scores` on an NVIDIA GeForce RTX 2070, the whole call:
packing the batch, copying it over, scoring it and copying back.

| workload | threads | dinara-align | dinara-align GPU | WFA2-lib | KSW2 | parasail | Edlib |
| :-- | :-- | --: | --: | --: | --: | --: | --: |
| illumina-affine | one core | **230 ms** | — | 1062 ms | 10.9 s | 21.7 s | — |
| illumina-edit | one core | **212 ms** | — | 450 ms | — | — | 2002 ms |
| illumina-affine | 10 cores | **27 ms** | 36 ms | 116 ms | 1211 ms | 2321 ms | — |
| illumina-edit | 10 cores | **29 ms** | — | 64 ms | — | — | 216 ms |

- **On the CPU a batch of global distances runs many pairs at once**, one a 16-bit lane of a register, 32
  under AVX-512, over a band of diagonals each pair's own cost proves: a path visiting a diagonal past the
  end diagonal and the main one runs a gap out to it and one back, so it costs at least their two openings
  and extensions, and a cost no dearer than every path off the band is the least. The pairs the first band
  does not prove are scored again over every diagonal a cheaper path could visit. Every pair goes first
  into a byte whose additions saturate, 64 to an AVX-512 register, a cost under 255 exact, and a pair whose
  byte saturates into 16 bits after. One core scores a pair in 0.46 µs at the affine costs and 0.42 at
  unit costs, WFA2-lib in 2.12 and 0.90; switching the certified band off triples the time (see Ablation).
- **On the GPU the batch takes 36 ms**, most of it the host packing the batch: each pair goes first over a
  band of sixteen diagonals, a thread a pair, and the band's score stands where the pair's own cost proves
  it, as on the CPU; the rest are scored whole on the device while the host packs the next chunk.

## Ablation

What each technique contributes: dinara-align built with one switched off at a time against the build
with none, everything else alike (`bench-ablation`, switches in `dinara_align/ablation.mojo`). Every
configuration aligns samples of the main benchmark's datasets with traceback, the 500 kbp reads ten times
the default sample, and the certified-band switch scores the short-read batch. The configurations
alternate within each of three rounds, so a drift in the machine falls on all of them, and the table gives
each one's median time, in ms, against the baseline's; every configuration must return the baseline's
costs. The Skylake-X penalises jumps that cross a 32-byte boundary, so where a rebuild places a loop moves
a time by several percent; the alignment builds keep their jumps clear of those boundaries
(`-mbranches-within-32B-boundaries`).

| switch | what it switches off |
| :-- | :-- |
| regroup | Myers' step as A\*PA2 writes it, rather than regrouped for a shorter chain (`bit_parallel.advance`) |
| paired | one group of eight words a sweep on AVX-512 rather than two interleaved |
| gather | diagonal transition's slides a scalar compare a lane rather than gathered (`slides`) |
| agreement | every projected bound trusted, whether or not the two projections agree (`diagonal.trusted_projection`) |
| neighbours | an exact inexact-seed match's neighbours found by the scan and pruned each alone, rather than built from the match |
| certified | a batch swept over each group's whole matrices rather than a band its pairs' costs certify (`lanes`) |

| dataset | baseline (ms) | regroup | paired | gather | agreement | neighbours | certified |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| ont-500k | 96.33 | +1.5% | +5.0% | +0.5% | +0.0% | +9.6% |  |
| ont-500k-genvar | 133.3 | +1.8% | +5.0% | +0.4% | -0.2% | +16.7% |  |
| ont-1k | 0.02082 | +0.4% | +0.4% | +15.6% | +0.3% | -0.1% |  |
| ont-10k | 0.1511 | -1.2% | +0.6% | +3.7% | +17.3% | -0.0% |  |
| sars-cov-2 | 0.1947 | +0.4% | +1.5% | +14.3% | +6.0% | -0.0% |  |
| l100000.n1000.e1 | 1.104 | +0.5% | +0.1% | +26.7% | +0.4% | +0.3% |  |
| l100000.n1000.e5 | 4.451 | -0.1% | -0.1% | +0.4% | +0.3% | +0.2% |  |
| l100000.n1000.e10 | 9.728 | +0.0% | +0.2% | +1.3% | +0.3% | +11.5% |  |
| l100000.n1000.e20 | 14.57 | +2.7% | +7.8% | +4.2% | +3.7% | +2.9% |  |
| illumina-affine | 228.2 |  |  |  |  |  | +192.0% |
| illumina-edit | 212.1 |  |  |  |  |  | +214.4% |

Every technique pays where it was built to: the seed neighbours 10 to 17% on the ultra-long reads and 12%
on 100 kbp pairs at 10%, gathered diagonal transition 14 to 27% where diagonal transition does most of
the work, the agreement rule 17% on the 10 kbp reads, interleaved groups 5 to 8% on the long and divergent
pairs, the regrouped step 1.5 to 2.7%, and the certified band three times on the short reads; switching one
off never made a dataset more than 1.2% faster. Measured at `ffbf46e`: the alignments here are at unit
costs, which neither the divergent-pair handover nor the later changes to it touch, and their bands' edges
fit under the memory cap.

## Reading the Numbers

- **One machine.** Every number here comes from one Skylake-X. Rows of microseconds move by a few percent
  between runs, and on this core by up to a third between builds that do not touch the code involved,
  where a loop lands against the micro-op cache (Intel's JCC erratum fix); read a difference of that size
  as a tie, and rerun on your own hardware before quoting any of this.
- **One thread, the caller's.** dinara-align is a library: it runs on its caller's thread and starts none
  of its own, so every column but the short reads' ten-core rows is one core.
- **Samples, not whole datasets.** The totals over each dataset in `bench.md` multiply a sample's mean by
  the dataset's pairs, marked ≈ where the sample is not the whole dataset; QuickEd's tables align every
  pair.
- **Memory growth leaves the input and the runtime out**, as A\*PA2's Table 10 does; peak resident
  memory, as QuickEd reports it, includes them, and dinara-align's runtime holds about 13 MB where the
  others' hold 2 to 8.
