# Benchmarks

Times dinara-align against the exact aligners of the two latest evaluations, A\*PA2's (Groot Koerkamp,
Bioinformatics 2025) and QuickEd's (Doblas et al., Bioinformatics 2025), on their datasets, and fails
any run in which two tools disagree on an answer.

```bash
pixi run bench                       # unit costs, every dataset, 5 s a tool a dataset
pixi run bench --affine 4,6,2        # gap-affine costs as WFA counts them: mismatch, opening, extension
pixi run bench --only ont l10000.    # the datasets whose names hold one of these
pixi run bench --budget 30           # the budget the tables below ran with
pixi run bench --repeats 3           # every tool three times a dataset, the rounds alternating
pixi run bench --fresh               # the rivals too, rather than their kept results
pixi run bench-modes                 # local, overlap, infix, prefix, extension, two-piece, tables
pixi run bench-batch                 # 500,000 short reads, on one core and on every core
pixi run bench-ablation              # dinara-align with each technique switched off in turn
pixi run profile all                 # where the time goes and where a rival wins (see Profiling)
```

Each writes its tables to the terminal and to `.cache/results/`: `bench.md` and `bench-affine-4-6-2.md`
with every table below and more (the totals over each whole dataset among them), `bench-exact.json` and
`bench-affine-4-6-2-exact.json` with every tool's unrounded figures, from which the paper's tables are
written, `bench-scaling.json` for its charts, `modes.md` and `batch-results.md`.

The tables below are the run the paper reports, kept in `results/2026-10-11/`: every tool timed again in
one session (`--budget 30 --fresh --repeats 3`), a pair's time its median over three alternating rounds
and a tool credited with the pairs it finished in every round; the modes and the batches are each the
median of three runs.

## What Runs

Every rival is cloned at a pinned commit into `.cache/` and built there for the host's instruction set;
nothing of theirs is vendored. A pinned rival's results do not change, so each runs on a sample once and
its rows are kept, keyed by its pin, its runner's source and the CPU; a rerun times dinara-align alone,
in seconds to minutes, and still checks every answer against every other.

| tool | commit | run on |
| :-- | :-- | :-- |
| dinara-align | this tree | everything, pair by pair through one `Aligner` a file, as every rival's runner keeps one aligner a file, and as a batch, every pair of a sample in one call on the same thread |
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
marks a tool its budget stopped before a pair, *o.o.m.* one the cap stopped, as QuickEd's tables do, and
*failed* one that crashed: BSAlign's edit distance on every pair of 500 kbp or more, with or without the cap.

**The other suites** time many short calls, where one cold call measures the allocator as much as the
aligner: each tool's time is the faster of two passes over a workload, its mean per pair, in process.

All tables were measured on an Intel Core i9-7900X (Skylake-X, AVX-512), ten cores at a fixed 3.3 GHz,
turbo boost and hyper-threading off, 91 GB, Ubuntu 26.04, on the otherwise idle machine, every tool on one
core pinned with `taskset` unless a row says otherwise, dinara-align at `87d365f`
for every table but the ablation (`dc980d7`); Rust nightly
2026-08-14, GCC 15.2, Mojo 1.1.0.

## Datasets

A\*PA2's real datasets, as pa-bench releases them: Oxford Nanopore reads of about 1, 10, 50 and 500 kbp,
the 500 kbp ones of CHM13 and, in ont-500k-genvar, of another individual, which adds genetic variation
(QuickEd's MiniION-a and MiniION-b), and
SARS-CoV-2 genomes. Real Illumina, PacBio HiFi and PromethION reads, the technologies QuickEd adds
(`giab.py`): Genome in a Bottle's HG002 reads, its 2x250 bp Illumina reads of chr20:40.00-40.02 Mbp, its
Sequel II HiFi reads of chr20:40.0-40.2 Mbp and UCSC's ultra-long PromethION reads of chr12:80.0-80.5 Mbp
of 10 kbp or more (ont-promethion, as QuickEd's PromethION set), streamed from GIAB's GRCh38 alignments, each read mapped by minimap2 to both haplotypes of HG002's own assembly, the T2T HG002 Q100
v1.1, and paired with the stretch it covers, so a pair differs by sequencing errors alone; reads clipped
past a tenth are left out. The assembly's chromosome 20 is cut from `DINARA_HG002` if that names a local
copy of `hg002v1.1.fasta`, else streamed from its S3 release; samtools and minimap2 come through
`pixi exec`. QuickEd's simulated ones, from Zenodo: random sequences of 10 kbp, 100 kbp, 500 kbp and
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
| illumina-250 | 4032 of 5992 | 77 / 248 / 251 | 0 / 0 / 10 |
| pacbio-hifi | 69 of 890 | 10,956 / 14,372 / 22,816 | 0 / 0 / 1 |
| ont-promethion | 22 of 584 | 10,362 / 44,631 / 155,664 | 4 / 10 / 20 |
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
| ont-1k | 1221 of 12477 | 20.8 µs | 21.4 µs | 81 µs | 48.1 µs | 613 µs | 124 µs | 44.2 µs | 31.3 µs | 538 µs | 73.2 µs | 119 µs | ✓ |
| ont-10k | 277 of 5000 | 154 µs | 151 µs | 393 µs | 247 µs | 15.4 ms | 1.11 ms | 752 µs | 599 µs | 13 ms | 679 µs | 1.3 ms | ✓ |
| ont-50k | 104 of 10000 | 582 µs | 562 µs | 1.28 ms | 941 µs | 263 ms | 6.11 ms | 7.43 ms | 6.35 ms | 169 ms | 2.31 ms | 12.3 ms | ✓ |
| ont-500k | 4 of 50 | 118 ms | 118 ms | 284 ms | 688 ms | t.o. (30 s) | 5 s | 6.37 s | o.o.m. | o.o.m. | 848 ms | failed | ✓ |
| ont-500k-genvar | 4 of 48 | 156 ms | 156 ms | 225 ms | 512 ms | t.o. (30 s) | 5.23 s | 6.18 s | o.o.m. | o.o.m. | 1.39 s | failed | ✓ |
| sars-cov-2 | 33 of 10000 | 277 µs | 255 µs | 2.01 ms | 727 µs | 6.55 ms | 8.18 ms | 906 µs | 979 µs | 68.6 ms | 3.86 ms | 47 ms | ✓ |
| illumina-250 | 4032 of 5992 | 0.593 µs | 0.799 µs | 14.5 µs | 4.11 µs | 13.7 µs | 23.4 µs | 2.06 µs | 1.12 µs | 12.4 µs | 16.2 µs | 20.5 µs | ✓ |
| pacbio-hifi | 69 of 890 | 12.9 µs | 15.7 µs | 980 µs | 259 µs | 1.06 ms | 3.41 ms | 86.5 µs | 66 µs | 3.09 ms | 951 µs | 12.2 ms | ✓ |
| ont-promethion | 22 of 584 | 3.99 ms | 3.93 ms | 7.02 ms | 6.76 ms | 929 ms | 73.9 ms | 71.1 ms | 79.7 ms | 1.81 s (15/22) | 17.1 ms | 181 ms (failed after 7/22) | ✓ |
| l10000.n10000.e1 | 100 of 10000 | 25.8 µs | 28.1 µs | 738 µs | 160 µs | 566 µs | 3.25 ms | 92.3 µs | 62.3 µs | 6.31 ms | 614 µs | 5.55 ms | ✓ |
| l10000.n10000.e5 | 99 of 10000 | 268 µs | 266 µs | 795 µs | 288 µs | 1.08 ms | 3.81 ms | 637 µs | 513 µs | 21.4 ms | 1.26 ms | 5.58 ms | ✓ |
| l10000.n10000.e10 | 100 of 10000 | 348 µs | 358 µs | 724 µs | 444 µs | 3.23 ms | 4.57 ms | 2.26 ms | 1.84 ms | 40.1 ms | 2.01 ms | 5.61 ms | ✓ |
| l10000.n10000.e20 | 99 of 10000 | 511 µs | 516 µs | 1.44 ms | 1.01 ms | 204 ms | 5.96 ms | 7.31 ms | 5.34 ms | 98 ms | 1.54 ms | 5.64 ms | ✓ |
| l100000.n1000.e1 | 9 of 1000 | 1.47 ms | 1.42 ms | 7.53 ms | 3.26 ms | 7.76 ms | 58.1 ms | 3.72 ms | 4.08 ms | 617 ms | 12.5 ms | 642 ms | ✓ |
| l100000.n1000.e5 | 9 of 1000 | 5.1 ms | 5.09 ms | 8.35 ms | 14.9 ms | 14.1 ms | 147 ms | 63.3 ms | 54.9 ms | 4.77 s (6/9) | 23.9 ms | 642 ms | ✓ |
| l100000.n1000.e10 | 10 of 1000 | 7.24 ms | 7.16 ms | 13.9 ms | 27.7 ms | 48.4 ms | 247 ms | 218 ms | 182 ms | 9.24 s (3/10) | 39 ms | 629 ms | ✓ |
| l100000.n1000.e20 | 10 of 1000 | 15.8 ms | 15.8 ms | 36.8 ms | 56.8 ms | t.o. (30 s) | 445 ms | 709 ms | 553 ms | 17.5 s (1/10) | 62.6 ms | 630 ms | ✓ |
| l500000.n200.e1 | 4 of 200 | 33 ms | 35.4 ms | 43.1 ms | 72.5 ms | 41.2 ms | 669 ms | 82.1 ms | 75.8 ms | 24.4 s (1/4) | 124 ms | failed | ✓ |
| l500000.n200.e5 | 4 of 200 | 31.2 ms | 31.1 ms | 44.7 ms | 213 ms | 81.7 ms | 2.29 s | 1.55 s | 1.48 s | o.o.m. | 394 ms | failed | ✓ |
| l500000.n200.e10 | 4 of 200 | 67.8 ms | 68.6 ms | 119 ms | 423 ms | 407 ms | 4.38 s | 5.45 s | 5.02 s | o.o.m. | 699 ms | failed | ✓ |
| l500000.n200.e20 | 4 of 200 | 253 ms | 251 ms | 544 ms | 898 ms | t.o. (30 s) | 8.53 s (3/4) | 18 s (1/4) | o.o.m. | o.o.m. | 1.24 s | failed | ✓ |
| l1000000.n100.e1 | 4 of 100 | 61.1 ms | 57.1 ms | 89.2 ms | 253 ms | 96.6 ms | 2.37 s | 365 ms | 392 ms | o.o.m. | 393 ms | failed | ✓ |
| l1000000.n100.e5 | 4 of 100 | 65 ms | 64.3 ms | 116 ms | 819 ms | 183 ms | 8.83 s (3/4) | 6.35 s | o.o.m. | o.o.m. | 1.44 s | failed | ✓ |
| l1000000.n100.e10 | 4 of 100 | 185 ms | 185 ms | 368 ms | 1.65 s | 974 ms | 17.3 s (1/4) | 22.3 s (1/4) | o.o.m. | o.o.m. | 2.65 s | failed | ✓ |
| l1000000.n100.e20 | 4 of 100 | 890 ms | 898 ms | 2.02 s | 3.51 s | t.o. (30 s) | t.o. (30 s) | t.o. (30 s) | o.o.m. | o.o.m. | 4.8 s | failed | ✓ |

Each rival's time over dinara-align's on the pairs both aligned, across the datasets, as QuickEd's text
summarizes its own:

| rival | least | geometric mean | most | datasets |
| :-- | --: | --: | --: | --: |
| a*pa2-full | 1.3× | 3.1× | 76.1× | 25 |
| a*pa2-simple | 1.1× | 3.6× | 20.1× | 25 |
| a*pa | 1.2× | 14.6× | 452.1× | 20 |
| edlib | 6.0× | 32.1× | 264.6× | 24 |
| biwfa | 2.1× | 13.0× | 109.0× | 24 |
| wfa | 1.5× | 7.1× | 74.0× | 19 |
| ksw2 | 20.8× | 229.0× | 1193.1× | 16 |
| quicked | 3.0× | 7.8× | 73.8× | 25 |
| bsalign | 5.7× | 52.7× | 943.9× | 15 |

Growth of peak memory aligning a pair, median and largest, as A\*PA2's Table 10 measures it:

| dataset | dinara-align (bit-parallel, 1 thread) | dinara-align (bit-parallel, batch, 1 thread) | a*pa2-full | a*pa2-simple | a*pa | edlib | biwfa | wfa | ksw2 | quicked | bsalign |
| :-- | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: | --: |
| ont-1k | 0.0 / 0.1 MB | — | 0.0 / 0.4 MB | 0.0 / 0.3 MB | 0.0 / 0.4 MB | 0.0 / 0.5 MB | 0.0 / 0.3 MB | 0.0 / 0.5 MB | 0.0 / 0.4 MB | 0.0 / 0.1 MB | 0.0 / 0.2 MB |
| ont-10k | 0.0 / 0.6 MB | — | 0.0 / 0.7 MB | 0.0 / 0.4 MB | 0.0 / 20 MB | 0.0 / 0.5 MB | 0.0 / 1.0 MB | 0.0 / 6.3 MB | 0.0 / 29 MB | 0.0 / 2.5 MB | 0.0 / 13 MB |
| ont-50k | 0.0 / 1.8 MB | — | 0.0 / 1.4 MB | 0.0 / 0.8 MB | 0.0 / 410 MB | 0.0 / 0.8 MB | 0.0 / 1.7 MB | 0.0 / 71 MB | 0.0 / 167 MB | 0.0 / 16 MB | 0.0 / 163 MB |
| ont-500k | 0.0 / 48 MB | — | 0.0 / 81 MB | 0.0 / 82 MB | — | 0.2 / 5.3 MB | 4.5 / 15 MB | — | — | 0.0 / 22 MB | — |
| ont-500k-genvar | 6.6 / 39 MB | — | 0.0 / 49 MB | 15 / 32 MB | — | 2.4 / 4.6 MB | 6.0 / 12 MB | — | — | 4.0 / 19 MB | — |
| sars-cov-2 | 0.0 / 2.4 MB | — | 0.0 / 2.1 MB | 0.0 / 0.8 MB | 0.0 / 13 MB | 0.0 / 1.2 MB | 0.0 / 1.7 MB | 0.0 / 15 MB | 0.0 / 97 MB | 0.0 / 7.9 MB | 0.0 / 211 MB |
| illumina-250 | 0.0 / 0.1 MB | — | 0.0 / 0.2 MB | 0.0 / 0.3 MB | 0.0 / 0.1 MB | 0.0 / 0.1 MB | 0.0 / 0.3 MB | 0.0 / 0.1 MB | 0.0 / 0.1 MB | 0.0 / 0.0 MB | 0.0 / 0.1 MB |
| pacbio-hifi | 0.0 / 0.2 MB | — | 0.0 / 0.9 MB | 0.0 / 0.5 MB | 0.0 / 0.6 MB | 0.0 / 0.8 MB | 0.0 / 0.5 MB | 0.0 / 0.6 MB | 0.0 / 6.4 MB | 0.0 / 1.0 MB | 0.0 / 32 MB |
| ont-promethion | 0.0 / 3.9 MB | — | 0.0 / 4.6 MB | 0.0 / 2.0 MB | 0.0 / 338 MB | 0.0 / 1.0 MB | 0.0 / 2.6 MB | 0.1 / 639 MB | 0.0 / 3587 MB | 0.0 / 8.5 MB | 0.0 / 1288 MB |
| l10000.n10000.e1 | 0.0 / 0.1 MB | — | 0.0 / 0.8 MB | 0.0 / 0.4 MB | 0.0 / 0.6 MB | 0.0 / 0.6 MB | 0.0 / 0.6 MB | 0.0 / 0.7 MB | 0.0 / 3.1 MB | 0.0 / 0.5 MB | 0.0 / 24 MB |
| l10000.n10000.e5 | 0.0 / 0.8 MB | — | 0.0 / 0.7 MB | 0.0 / 0.5 MB | 0.0 / 0.6 MB | 0.0 / 0.7 MB | 0.0 / 1.4 MB | 0.0 / 2.1 MB | 0.0 / 11 MB | 0.0 / 1.5 MB | 0.0 / 24 MB |
| l10000.n10000.e10 | 0.0 / 0.2 MB | — | 0.0 / 1.0 MB | 0.0 / 0.5 MB | 0.0 / 2.4 MB | 0.0 / 0.6 MB | 0.0 / 1.4 MB | 0.0 / 4.3 MB | 0.0 / 20 MB | 0.0 / 2.6 MB | 0.0 / 24 MB |
| l10000.n10000.e20 | 0.0 / 0.5 MB | — | 0.0 / 0.8 MB | 0.0 / 0.4 MB | 0.0 / 28 MB | 0.0 / 0.7 MB | 0.0 / 1.4 MB | 0.0 / 11 MB | 0.0 / 40 MB | 0.0 / 4.5 MB | 0.0 / 24 MB |
| l100000.n1000.e1 | 0.0 / 4.9 MB | — | 0.0 / 6.0 MB | 0.0 / 2.0 MB | 0.0 / 4.1 MB | 0.0 / 1.4 MB | 0.0 / 4.0 MB | 3.4 / 5.9 MB | 0.0 / 203 MB | 0.0 / 12 MB | 0.0 / 2388 MB |
| l100000.n1000.e5 | 0.0 / 3.9 MB | — | 0.0 / 5.5 MB | 0.0 / 2.6 MB | 0.0 / 4.5 MB | 0.0 / 1.4 MB | 0.0 / 3.0 MB | 3.8 / 88 MB | 0.1 / 1555 MB | 0.0 / 11 MB | 0.0 / 2392 MB |
| l100000.n1000.e10 | 0.0 / 4.0 MB | — | 0.0 / 5.6 MB | 0.0 / 3.4 MB | 0.0 / 18 MB | 0.0 / 1.4 MB | 0.0 / 3.4 MB | 5.1 / 301 MB | 22 / 2986 MB | 0.0 / 17 MB | 0.0 / 2386 MB |
| l100000.n1000.e20 | 0.0 / 5.3 MB | — | 0.0 / 7.2 MB | 0.0 / 4.6 MB | — | 0.0 / 1.4 MB | 0.0 / 3.7 MB | 2.1 / 955 MB | 5476 / 5476 MB | 0.0 / 8.8 MB | 0.0 / 2392 MB |
| l500000.n200.e1 | 9.1 / 23 MB | — | 1.0 / 29 MB | 0.1 / 11 MB | 0.0 / 22 MB | 0.7 / 3.4 MB | 3.1 / 13 MB | 8.9 / 100 MB | 7835 / 7835 MB | 0.0 / 18 MB | — |
| l500000.n200.e5 | 4.4 / 16 MB | — | 0.7 / 26 MB | 0.2 / 19 MB | 2.9 / 23 MB | 0.0 / 4.0 MB | 2.8 / 12 MB | 22 / 2111 MB | — | 0.0 / 19 MB | — |
| l500000.n200.e10 | 3.2 / 18 MB | — | 0.2 / 30 MB | 0.1 / 30 MB | 5.4 / 76 MB | 0.1 / 4.0 MB | 2.8 / 14 MB | 123 / 7416 MB | — | 0.4 / 13 MB | — |
| l500000.n200.e20 | 14 / 53 MB | — | 0.1 / 59 MB | 0.3 / 53 MB | — | 0.0 / 4.0 MB | 15 / 15 MB | — | — | 0.6 / 15 MB | — |
| l1000000.n100.e1 | 11 / 37 MB | — | 2.0 / 57 MB | 0.2 / 28 MB | 0.0 / 43 MB | 0.1 / 7.6 MB | 6.8 / 17 MB | 8.8 / 382 MB | — | 0.0 / 25 MB | — |
| l1000000.n100.e5 | 11 / 30 MB | — | 2.2 / 53 MB | 0.2 / 57 MB | 5.6 / 46 MB | 0.0 / 7.9 MB | 6.7 / 16 MB | — | — | 0.9 / 25 MB | — |
| l1000000.n100.e10 | 4.3 / 39 MB | — | 0.6 / 67 MB | 0.6 / 99 MB | 13 / 151 MB | 7.9 / 7.9 MB | 19 / 19 MB | — | — | 0.9 / 24 MB | — |
| l1000000.n100.e20 | 14 / 123 MB | — | 0.4 / 184 MB | 2.9 / 192 MB | — | — | — | — | — | 0.1 / 26 MB | — |

- **dinara-align is the fastest on every one of the 25 datasets**, 1.1 to 5.1 times the next fastest exact
  aligner: 20.8 against WFA's 31.3 µs on ont-1k, 582 against A\*PA2-simple's 941 µs on ont-50k, 3.99
  against 6.76 ms on the PromethION reads, 118 against A\*PA2-full's 284 ms on the 500 kbp reads, 277
  against 727 µs on the SARS-CoV-2 genomes, 0.593 against WFA's 1.12 µs on the Illumina reads and 12.9
  against 66.0 µs on the HiFi reads, and 890 ms against 2.02 s on the 1 Mbp pairs at 20%. The closest is
  QuickEd's 10 kbp pairs at 5%, 268 against A\*PA2-simple's 288 µs. Against QuickEd itself, the latest of
  them, 3.0 to 74 times, 7.8 in the geometric mean.
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
  123 MB, where it added 167 before the cap, against A\*PA2's 184 and 192; the rest of that is the
  profile and the seeds, which grow with the length. QuickEd, splitting Hirschberg's way past 16 MB,
  stays near 25 MB.

## Affine Costs

Gap-affine alignment at WFA's costs (4, 6, 2), a mismatch 4 and a gap of `k` letters `6 + 2k`, against the
exact aligners that take them, a 30 s budget a tool a dataset:

| dataset | pairs | dinara-align (1 thread) | dinara-align (batch, 1 thread) | wfa | biwfa | ksw2 | bsalign | agree |
| :-- | --: | --: | --: | --: | --: | --: | --: | :-: |
| ont-1k | 1221 of 12477 | 105 µs | 106 µs | 189 µs | 305 µs | 1.1 ms | 336 µs | ✓ |
| ont-10k | 277 of 5000 | 1.59 ms | 1.62 ms | 4.86 ms | 5.85 ms | 29.9 ms | 7.83 ms | ✓ |
| ont-50k | 104 of 10000 | 21.7 ms | 21.8 ms | 63.4 ms | 54.9 ms | 341 ms (81/104) | 86.8 ms | ✓ |
| ont-500k | 4 of 50 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| ont-500k-genvar | 4 of 48 | 27.3 s (1/4) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| sars-cov-2 | 33 of 10000 | 698 µs | 672 µs | 2.06 ms | 1.78 ms | 88.6 ms | 341 ms | ✓ |
| illumina-250 | 4032 of 5992 | 1.51 µs | 1.81 µs | 1.54 µs | 2.82 µs | 12.8 µs | 53.8 µs | ✓ |
| pacbio-hifi | 69 of 890 | 92.9 µs | 812 µs | 224 µs | 232 µs | 9.74 ms | 86.1 ms | ✓ |
| ont-promethion | 22 of 584 | 286 ms | 265 ms | 286 ms (o.o.m. after 7/22) | 576 ms | 2.35 s (7/22) | o.o.m. | ✓ |
| l10000.n10000.e1 | 100 of 10000 | 195 µs | 218 µs | 409 µs | 649 µs | 21.2 ms | 37 ms | ✓ |
| l10000.n10000.e5 | 99 of 10000 | 2.37 ms | 2.37 ms | 7.13 ms | 8.73 ms | 97.6 ms | 37.2 ms | ✓ |
| l10000.n10000.e10 | 100 of 10000 | 7.08 ms | 7.07 ms | 22.8 ms | 27.2 ms | 204 ms | 37.2 ms | ✓ |
| l10000.n10000.e20 | 99 of 10000 | 19.9 ms | 20.3 ms | 65.5 ms | 75.2 ms | 375 ms (80/99) | 37.5 ms | ✓ |
| l100000.n1000.e1 | 9 of 1000 | 11.9 ms | 12.1 ms | 39.9 ms | 39 ms | 2.44 s | o.o.m. | ✓ |
| l100000.n1000.e5 | 9 of 1000 | 357 ms | 340 ms | 810 ms | 715 ms | 9.22 s (3/9) | o.o.m. | ✓ |
| l100000.n1000.e10 | 10 of 1000 | 1.37 s | 1.26 s | o.o.m. | 2.43 s | 17.4 s (1/10) | o.o.m. | ✓ |
| l100000.n1000.e20 | 10 of 1000 | 3.94 s (7/10) | t.o. (30 s) | o.o.m. | 7.35 s (4/10) | o.o.m. | o.o.m. | ✓ |
| l500000.n200.e1 | 4 of 200 | 499 ms | 479 ms | 1.29 s | 942 ms | o.o.m. | o.o.m. | ✓ |
| l500000.n200.e5 | 4 of 200 | 11.4 s (2/4) | t.o. (30 s) | o.o.m. | 19.4 s (1/4) | o.o.m. | o.o.m. | ✓ |
| l500000.n200.e10 | 4 of 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| l500000.n200.e20 | 4 of 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| l1000000.n100.e1 | 4 of 100 | 2.44 s | 2.37 s | o.o.m. | 3.92 s | o.o.m. | o.o.m. | ✓ |
| l1000000.n100.e5 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| l1000000.n100.e10 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| l1000000.n100.e20 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |

| rival | least | geometric mean | most | datasets |
| :-- | --: | --: | --: | --: |
| wfa | 1.0× | 2.6× | 3.4× | 14 |
| biwfa | 1.6× | 2.5× | 3.8× | 18 |
| ksw2 | 8.5× | 33.1× | 205.2× | 14 |
| bsalign | 1.9× | 22.3× | 926.5× | 10 |

- **dinara-align is the fastest on every dataset any tool finished, the Illumina reads a tie**: 105 against
  WFA's 189 µs on ont-1k, 1.62 against 4.86 ms on ont-10k, 21.7 against BiWFA's 54.9 ms on ont-50k, 286
  against BiWFA's 576 ms on the PromethION reads, where WFA ran out of memory after 7 of 22, 698 µs against
  1.78 ms on the SARS-CoV-2 genomes, 92.9 against 225 µs on the HiFi reads, and 19.9 against BSAlign's
  37.5 ms on 10 kbp pairs at 20%, where BSAlign's 8-bit difference sweep, whose time does not depend on the
  divergence, comes close. On the Illumina reads it takes 1.51 against WFA's 1.54 µs, within the spread of
  its rounds: a pair that close is aligned by one search, as WFA aligns it, grown alone while it steps no
  more diagonals than the two searches would, its fronts measured as it grows.
- **The batch call is slow on HiFi reads at gap-affine costs**, 812 µs a pair against 93 one at a time,
  for a reason not yet found (see `TODO.md`).
- **It finishes more of the long pairs than anyone**: 7 of the 10 100 kbp pairs at 20% in every round where
  BiWFA finishes 4, two 500 kbp pairs at 5% to BiWFA's one, and one ont-500k-genvar read. No tool finishes a
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
| infix-edit | **130 µs** | 268 µs | 1.00 ms | — | — | — |
| infix-edit-long | **2.01 ms** | 13.5 ms | 54.0 ms | — | — | — |
| prefix-edit | **42.4 µs** | 114 µs | 47.0 µs | — | — | — |
| infix-affine | **3.92 ms** | — | 3.96 ms | — | — | — |
| extension | **726 µs** | — | — | 2.06 ms | — | — |
| extension-bonus | **166 µs** | — | — | 1.24 ms | — | — |
| two-piece | **8.25 ms** | — | 15.8 ms | 46.9 ms | — | — |
| table-global | **196 µs** | — | — | — | 2.94 ms | — |
| table-local | **1.91 ms** | — | — | — | 5.16 ms | 3.28 ms |
| local-short | **32.0 µs** | — | — | — | 67.2 µs | 51.2 µs |
| local-window | **1.36 ms** | — | — | — | 4.75 ms | 3.15 ms |
| local-long | **17.4 ms** | — | — | — | 293 ms | 84.7 ms |
| overlap-reads | **703 µs** | — | — | — | 15.3 ms | — |
| infix-reads | **687 µs** | — | — | — | 15.1 ms | — |

- **dinara-align is the fastest on every workload, or tied**: twice Edlib on short infixes and 6.7 times
  on long ones, 2.8 times KSW2 on extension and 7.5 times with an end bonus, 1.9 times WFA2-lib on two-piece
  gaps, 15 times parasail on a table globally, and 4.9 times SSW on the long local pairs; the infix at
  WFA2-lib's costs, 3.92 against 3.96 ms, is a tie.
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
  the Skylake-X's micro-op cache (Intel's JCC erratum fix). The infix at WFA2-lib's costs, 3.92 ms here,
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
| illumina-affine | one core | **231 ms** | — | 1.06 s | 10.9 s | 21.7 s | — |
| illumina-edit | one core | **212 ms** | — | 450 ms | — | — | 2.00 s |
| illumina-affine | 10 cores | **29.0 ms** | 35.0 ms | 113 ms | 1.20 s | 2.35 s | — |
| illumina-edit | 10 cores | **28.0 ms** | — | 49.0 ms | — | — | 209 ms |

- **On the CPU a batch of global distances runs many pairs at once**, one a 16-bit lane of a register, 32
  under AVX-512, over a band of diagonals each pair's own cost proves: a path visiting a diagonal past the
  end diagonal and the main one runs a gap out to it and one back, so it costs at least their two openings
  and extensions, and a cost no dearer than every path off the band is the least. The pairs the first band
  does not prove are scored again over every diagonal a cheaper path could visit. Every pair goes first
  into a byte whose additions saturate, 64 to an AVX-512 register, a cost under 255 exact, and a pair whose
  byte saturates into 16 bits after. One core scores a pair in 0.46 µs at the affine costs and 0.42 at
  unit costs, WFA2-lib in 2.13 and 0.90; switching the certified band off triples the time (see Ablation).
- **On the GPU the batch takes 35 ms** (the median of three runs, which varied by half), most of it the host packing the batch: each pair goes first over a
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
| ont-500k | 98.8 | +1.1% | +4.5% | +0.5% | -0.2% | +8.9% |  |
| ont-500k-genvar | 136.0 | +1.9% | +5.0% | +0.4% | 0.0% | +16.4% |  |
| ont-1k | 0.021 | +0.2% | +0.3% | +14.9% | 0.0% | +0.7% |  |
| ont-10k | 0.145 | -0.6% | +1.2% | +6.4% | +22.2% | +0.1% |  |
| sars-cov-2 | 0.197 | -0.9% | +0.2% | +13.4% | +4.9% | -0.5% |  |
| l100000.n1000.e1 | 1.10 | -1.4% | +0.4% | +26.5% | -0.1% | -1.0% |  |
| l100000.n1000.e5 | 4.61 | -0.6% | -0.2% | +0.1% | 0.0% | +0.2% |  |
| l100000.n1000.e10 | 6.72 | +0.3% | +3.4% | +1.9% | +0.5% | -0.2% |  |
| l100000.n1000.e20 | 15.0 | +4.1% | +11.5% | +3.9% | -2.9% | 0.0% |  |
| illumina-affine | 232 |  |  |  |  |  | +188% |
| illumina-edit | 213 |  |  |  |  |  | +215% |

Every technique pays where it was built to: the seed neighbours 9 and 16% on the ultra-long reads, gathered
diagonal transition 6 to 27% where diagonal transition does most of the work, the agreement rule 22% on the
ont-10k reads, interleaved groups 4.5 to 11.5% on the long and divergent pairs, the regrouped step up to
4.1%, and the certified band about three times on the short reads; switching one off made a dataset at most
2.9% faster. Measured at `dc980d7`, as the paper's Table S2 gives it: the alignments here are at unit costs,
which neither the divergent-pair handover nor the later changes to it touch, and their bands' edges fit under
the memory cap.

## Profiling

`pixi run profile` pinpoints where dinara-align's time goes and where a rival beats it, on the samples
`bench` times (`profile.py`):

```bash
pixi run profile duel [--affine 4,6,2]   # every pair against every rival's time for it
pixi run profile probes --only ont-10k   # per pair: engine, phases, work, counters, allocations, copies
pixi run profile perf --only ont-10k     # Linux: counters and hot functions, files and lines; who copies
pixi run profile compare main~1          # this checkout against a commit, pair by pair, with intervals
pixi run profile all                     # duel, probes and perf in one report
```

- **duel** sets dinara-align's median time per pair, over `--rounds` rounds, beside each rival's kept time
  for the same pair. A dataset's mean can win while pairs lose: it lists the margin over the fastest
  rival, the pairs some rival aligned faster, those faster than dinara-align's fastest round (beyond its
  noise), the share of time those losses cost, and the worst of them with length, cost and divergence.
- **probes** builds `probed.mojo` with `-D DINARA_PROBES`, which compiles in the library's probes
  (`dinara_align/probes.mojo`; a normal build folds every one away): the engine that settled each pair,
  each phase's time (encoding, diagonal transition from one end and two, the fronts grown, the seed
  heuristic, the band's rounds and traceback, the CIGAR, a full sweep, the wavefront, the difference
  sweep) and its work (diagonals stepped, bounds tried against the distance, tiles and cells swept, the
  share swept by rounds that failed, seeds built, matched and kept, tiles swept again or recomputed for
  the traceback, wavefront steps, meeting checks, splits). Around each pair it reads the processor's
  counters through `perf_event_open` (instructions, cycles, branches and misses, L1d, last-level and TLB
  misses) and, with `probes/counting.c` preloaded on Linux, every allocation, `memcpy`, `memmove` and
  `memset` that reaches the C library, with the bytes, and the page faults. The report gives each
  dataset's time by engine and phase, the work against what the pair needed, the processor's and the
  allocator's view per pair after the first, the pairs whose probes match a pattern that costs time
  (rounds wasted, a band wider than needed, seeds that took longer than the band, tiles recomputed,
  allocations in a warm aligner, copies beyond the sequences' size), and the slowest pairs in full.
- **perf** samples the runner built with line tables: cycles, IPC, branch and cache misses per sample,
  the share of samples by function, by source file (one engine each) and by source line, inlined code
  charged to its own line, and the samples in copy and allocation routines, with their callers where the
  processor keeps a last-branch record.
- **compare** builds the same runner against two commits, or a commit and the checkout as it stands,
  times both in rounds that alternate, holds every pair's costs equal, and gives each dataset's ratio of
  mean times with a 95% bootstrap interval over pairs and rounds: faster or slower only when the interval
  excludes 1. On x86 both builds keep jumps clear of 32-byte boundaries, as the ablation's do.

Every run pins itself to one core (`--core`, default 2) and stops if that core or the machine is busy,
unless given `--force`. Reports go to `.cache/results/profile-*.md`, each pair's figures to a JSON beside it.

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
