Every tool built for `native` (Rust `-C target-cpu=native`, C and C++ `-march=native`, Mojo its default, the host), run on Intel(R) Core(TM) i9-7900X CPU @ 3.30GHz.

Rivals at pa-bench `af7a50d`, astar-pairwise-aligner `bf2e14e`, WFA2-lib `bcf473a`, SSW `a66636b`, parasail `fb985ee`, edlib `ec2310e`, ksw2 `289609b`, QuickEd `2b1d44b`, bsalign `1707aab`. Compilers: cc (Ubuntu 15.2.0-16ubuntu1) 15.2.0; rustc 1.99.0-nightly (d453bdd8f 2026-08-14); Mojo 1.1.0 (8189361e). One thread each, pinned on Linux by the caller.

Datasets, over each sample (QuickEd's Table S1):

| dataset | pairs | length, min / mean / max (bp) | divergence, min / mean / max (%) |
| :-- | --: | --: | --: |
| ont-1k | 1221 of 12477 | 42 / 830 / 1,089 | — |
| ont-10k | 277 of 5000 | 342 / 3,643 / 9,791 | — |
| ont-50k | 104 of 10000 | 360 / 9,633 / 42,626 | — |
| ont-500k | 4 of 50 | 553,834 / 642,465 / 683,520 | — |
| ont-500k-genvar | 4 of 48 | 594,436 / 662,980 / 861,994 | — |
| sars-cov-2 | 33 of 10000 | 29,244 / 29,746 / 29,864 | — |
| illumina-250 | 4032 of 5992 | 77 / 248 / 251 | — |
| pacbio-hifi | 69 of 890 | 10,956 / 14,372 / 22,816 | — |
| ont-promethion | 22 of 584 | 10,362 / 44,631 / 155,664 | — |
| l10000.n10000.e1 | 100 of 10000 | 10,000 / 10,003 / 10,023 | — |
| l10000.n10000.e5 | 99 of 10000 | 10,000 / 10,008 / 10,045 | — |
| l10000.n10000.e10 | 100 of 10000 | 10,000 / 10,008 / 10,056 | — |
| l10000.n10000.e20 | 99 of 10000 | 10,000 / 10,014 / 10,095 | — |
| l100000.n1000.e1 | 9 of 1000 | 100,000 / 100,013 / 100,035 | — |
| l100000.n1000.e5 | 9 of 1000 | 100,000 / 100,028 / 100,123 | — |
| l100000.n1000.e10 | 10 of 1000 | 100,000 / 100,039 / 100,120 | — |
| l100000.n1000.e20 | 10 of 1000 | 100,000 / 100,040 / 100,155 | — |
| l500000.n200.e1 | 4 of 200 | 500,000 / 500,038 / 500,103 | — |
| l500000.n200.e5 | 4 of 200 | 500,016 / 500,088 / 500,178 | — |
| l500000.n200.e10 | 4 of 200 | 500,000 / 500,052 / 500,105 | — |
| l500000.n200.e20 | 4 of 200 | 500,000 / 500,055 / 500,220 | — |
| l1000000.n100.e1 | 4 of 100 | 1,000,000 / 1,000,034 / 1,000,107 | — |
| l1000000.n100.e5 | 4 of 100 | 1,000,185 / 1,000,301 / 1,000,551 | — |
| l1000000.n100.e10 | 4 of 100 | 1,000,000 / 1,000,054 / 1,000,120 | — |
| l1000000.n100.e20 | 4 of 100 | 1,000,000 / 1,000,166 / 1,000,355 | — |

Mean time per alignment over each dataset's sample, at most 30 s a tool a dataset, the pairs
a tool aligned before it was stopped counted; t.o. a tool its budget stopped before a pair, o.o.m. one the
8 GB memory cap stopped (A*PA2's view):

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

Median time per alignment over the same pairs (A*PA2's supplement):

| dataset | pairs | dinara-align (1 thread) | dinara-align (batch, 1 thread) | wfa | biwfa | ksw2 | bsalign |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| ont-1k | 1221 of 12477 | 112 µs | 106 µs | 201 µs | 279 µs | 1.04 ms | 393 µs |
| ont-10k | 277 of 5000 | 607 µs | 1.62 ms | 1.89 ms | 2.29 ms | 10.1 ms | 3.09 ms |
| ont-50k | 104 of 10000 | 1.77 ms | 21.8 ms | 5.2 ms | 6.77 ms | 44.8 ms (81/104) | 14.1 ms |
| ont-500k | 4 of 50 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| ont-500k-genvar | 4 of 48 | 27.3 s (1/4) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| sars-cov-2 | 33 of 10000 | 310 µs | 672 µs | 681 µs | 942 µs | 64 ms | 315 ms |
| illumina-250 | 4032 of 5992 | 1.12 µs | 1.81 µs | 1.21 µs | 1.99 µs | 12.1 µs | 54.2 µs |
| pacbio-hifi | 69 of 890 | 36.5 µs | 812 µs | 68.4 µs | 102 µs | 4.54 ms | 68.4 ms |
| ont-promethion | 22 of 584 | 63.7 ms | 265 ms | 102 ms (o.o.m. after 7/22) | 156 ms | 1.82 s (7/22) | o.o.m. |
| l10000.n10000.e1 | 100 of 10000 | 189 µs | 218 µs | 382 µs | 641 µs | 21.1 ms | 36.1 ms |
| l10000.n10000.e5 | 99 of 10000 | 2.33 ms | 2.37 ms | 6.92 ms | 8.72 ms | 97.1 ms | 36.1 ms |
| l10000.n10000.e10 | 100 of 10000 | 6.98 ms | 7.07 ms | 22.4 ms | 27.3 ms | 204 ms | 36.2 ms |
| l10000.n10000.e20 | 99 of 10000 | 19.5 ms | 20.3 ms | 64 ms | 75.2 ms | 374 ms (80/99) | 36.4 ms |
| l100000.n1000.e1 | 9 of 1000 | 10.2 ms | 12.1 ms | 33.6 ms | 38.6 ms | 2.43 s | o.o.m. |
| l100000.n1000.e5 | 9 of 1000 | 350 ms | 340 ms | 666 ms | 710 ms | 9.2 s (3/9) | o.o.m. |
| l100000.n1000.e10 | 10 of 1000 | 1.36 s | 1.26 s | o.o.m. | 2.44 s | 17.4 s (1/10) | o.o.m. |
| l100000.n1000.e20 | 10 of 1000 | 3.93 s (7/10) | t.o. (30 s) | o.o.m. | 7.34 s (4/10) | o.o.m. | o.o.m. |
| l500000.n200.e1 | 4 of 200 | 490 ms | 479 ms | 956 ms | 939 ms | o.o.m. | o.o.m. |
| l500000.n200.e5 | 4 of 200 | 11.4 s (2/4) | t.o. (30 s) | o.o.m. | 19.4 s (1/4) | o.o.m. | o.o.m. |
| l500000.n200.e10 | 4 of 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l500000.n200.e20 | 4 of 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l1000000.n100.e1 | 4 of 100 | 2.44 s | 2.37 s | o.o.m. | 3.91 s | o.o.m. | o.o.m. |
| l1000000.n100.e5 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l1000000.n100.e10 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l1000000.n100.e20 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |

Total over each whole dataset, its mean time per pair times its pairs, ≈ where the sample is not the
whole (QuickEd's Table 2 and S3):

| dataset | pairs | dinara-align (1 thread) | dinara-align (batch, 1 thread) | wfa | biwfa | ksw2 | bsalign |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| ont-1k | 12477 | ≈ 1.31 s | ≈ 1.32 s | ≈ 2.35 s | ≈ 3.81 s | ≈ 13.8 s | ≈ 4.19 s |
| ont-10k | 5000 | ≈ 7.94 s | ≈ 8.08 s | ≈ 24.3 s | ≈ 29.2 s | ≈ 150 s | ≈ 39.1 s |
| ont-50k | 10000 | ≈ 217 s | ≈ 218 s | ≈ 634 s | ≈ 549 s | ≈ 3.41e+03 s | ≈ 868 s |
| ont-500k | 50 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| ont-500k-genvar | 48 | ≈ 1.31e+03 s | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| sars-cov-2 | 10000 | ≈ 6.98 s | ≈ 6.72 s | ≈ 20.6 s | ≈ 17.8 s | ≈ 886 s | ≈ 3.41e+03 s |
| illumina-250 | 5992 | ≈ 9.03 ms | ≈ 10.8 ms | ≈ 9.25 ms | ≈ 16.9 ms | ≈ 76.7 ms | ≈ 322 ms |
| pacbio-hifi | 890 | ≈ 82.7 ms | ≈ 722 ms | ≈ 200 ms | ≈ 206 ms | ≈ 8.67 s | ≈ 76.6 s |
| ont-promethion | 584 | ≈ 167 s | ≈ 155 s | o.o.m. | ≈ 336 s | ≈ 1.38e+03 s | o.o.m. |
| l10000.n10000.e1 | 10000 | ≈ 1.95 s | ≈ 2.18 s | ≈ 4.09 s | ≈ 6.49 s | ≈ 212 s | ≈ 370 s |
| l10000.n10000.e5 | 10000 | ≈ 23.7 s | ≈ 23.7 s | ≈ 71.3 s | ≈ 87.3 s | ≈ 976 s | ≈ 372 s |
| l10000.n10000.e10 | 10000 | ≈ 70.8 s | ≈ 70.7 s | ≈ 228 s | ≈ 272 s | ≈ 2.04e+03 s | ≈ 372 s |
| l10000.n10000.e20 | 10000 | ≈ 199 s | ≈ 203 s | ≈ 655 s | ≈ 752 s | ≈ 3.75e+03 s | ≈ 375 s |
| l100000.n1000.e1 | 1000 | ≈ 11.9 s | ≈ 12.1 s | ≈ 39.9 s | ≈ 39 s | ≈ 2.44e+03 s | o.o.m. |
| l100000.n1000.e5 | 1000 | ≈ 357 s | ≈ 340 s | ≈ 810 s | ≈ 715 s | ≈ 9.22e+03 s | o.o.m. |
| l100000.n1000.e10 | 1000 | ≈ 1.37e+03 s | ≈ 1.26e+03 s | o.o.m. | ≈ 2.43e+03 s | ≈ 1.74e+04 s | o.o.m. |
| l100000.n1000.e20 | 1000 | ≈ 3.94e+03 s | t.o. (30 s) | o.o.m. | ≈ 7.35e+03 s | o.o.m. | o.o.m. |
| l500000.n200.e1 | 200 | ≈ 99.7 s | ≈ 95.9 s | ≈ 257 s | ≈ 188 s | o.o.m. | o.o.m. |
| l500000.n200.e5 | 200 | ≈ 2.28e+03 s | t.o. (30 s) | o.o.m. | ≈ 3.88e+03 s | o.o.m. | o.o.m. |
| l500000.n200.e10 | 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l500000.n200.e20 | 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l1000000.n100.e1 | 100 | ≈ 244 s | ≈ 237 s | o.o.m. | ≈ 392 s | o.o.m. | o.o.m. |
| l1000000.n100.e5 | 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l1000000.n100.e10 | 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l1000000.n100.e20 | 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |

Each rival's time over dinara-align (1 thread)'s, on the pairs both aligned:

| dataset | wfa | biwfa | ksw2 | bsalign |
| :-- | --: | --: | --: | --: |
| ont-1k | 1.8× | 2.9× | 10.5× | 3.2× |
| ont-10k | 3.1× | 3.7× | 18.8× | 4.9× |
| ont-50k | 2.9× | 2.5× | 19.4× | 4.0× |
| ont-500k | — | — | — | — |
| ont-500k-genvar | — | — | — | — |
| sars-cov-2 | 3.0× | 2.5× | 127× | 489× |
| illumina-250 | 1.0× | 1.9× | 8.5× | 35.7× |
| pacbio-hifi | 2.4× | 2.5× | 105× | 926× |
| ont-promethion | 3.4× | 2.0× | 28.2× | — |
| l10000.n10000.e1 | 2.1× | 3.3× | 109× | 189× |
| l10000.n10000.e5 | 3.0× | 3.7× | 41.2× | 15.7× |
| l10000.n10000.e10 | 3.2× | 3.8× | 28.8× | 5.3× |
| l10000.n10000.e20 | 3.3× | 3.8× | 18.8× | 1.9× |
| l100000.n1000.e1 | 3.4× | 3.3× | 205× | — |
| l100000.n1000.e5 | 2.3× | 2.0× | 25.2× | — |
| l100000.n1000.e10 | — | 1.8× | 12.5× | — |
| l100000.n1000.e20 | — | 1.9× | — | — |
| l500000.n200.e1 | 2.6× | 1.9× | — | — |
| l500000.n200.e5 | — | 1.7× | — | — |
| l500000.n200.e10 | — | — | — | — |
| l500000.n200.e20 | — | — | — | — |
| l1000000.n100.e1 | — | 1.6× | — | — |
| l1000000.n100.e5 | — | — | — | — |
| l1000000.n100.e10 | — | — | — | — |
| l1000000.n100.e20 | — | — | — | — |

The same across datasets, as QuickEd's text summarizes its own:

| rival | least | geometric mean | most | datasets |
| :-- | --: | --: | --: | --: |
| wfa | 1.0× | 2.6× | 3.4× | 14 |
| biwfa | 1.6× | 2.5× | 3.8× | 18 |
| ksw2 | 8.5× | 33.1× | 205.2× | 14 |
| bsalign | 1.9× | 22.3× | 926.5× | 10 |

Growth of peak memory aligning a pair, median / largest, input and runtime left out (A*PA2's Table 10):

| dataset | dinara-align (1 thread) | dinara-align (batch, 1 thread) | wfa | biwfa | ksw2 | bsalign |
| :-- | --: | --: | --: | --: | --: | --: |
| ont-1k | 0.0 / 0.5 MB | — | 0.0 / 3.0 MB | 0.0 / 1.6 MB | 0.0 / 0.9 MB | 0.0 / 1.9 MB |
| ont-10k | 0.0 / 22 MB | — | 0.0 / 97 MB | 0.0 / 2.4 MB | 0.0 / 58 MB | 0.0 / 104 MB |
| ont-50k | 0.0 / 85 MB | — | 0.0 / 1420 MB | 0.0 / 4.5 MB | 0.0 / 317 MB | 0.0 / 1303 MB |
| ont-500k | — | — | — | — | — | — |
| ont-500k-genvar | 144 / 144 MB | — | — | — | — | — |
| sars-cov-2 | 0.0 / 11 MB | — | 0.0 / 51 MB | 0.0 / 2.5 MB | 0.0 / 211 MB | 0.0 / 1682 MB |
| illumina-250 | 0.0 / 0.3 MB | — | 0.0 / 0.1 MB | 0.0 / 0.1 MB | 0.0 / 0.0 MB | 0.0 / 0.2 MB |
| pacbio-hifi | 0.0 / 1.8 MB | — | 0.0 / 7.1 MB | 0.0 / 1.5 MB | 0.0 / 26 MB | 0.0 / 258 MB |
| ont-promethion | 0.0 / 89 MB | — | 1.8 / 1768 MB | 0.0 / 8.9 MB | 0.0 / 2155 MB | — |
| l10000.n10000.e1 | 0.0 / 0.7 MB | — | 0.0 / 4.2 MB | 0.0 / 1.5 MB | 0.0 / 11 MB | 0.0 / 192 MB |
| l10000.n10000.e5 | 0.0 / 7.3 MB | — | 0.0 / 30 MB | 0.0 / 2.0 MB | 0.0 / 40 MB | 0.0 / 192 MB |
| l10000.n10000.e10 | 0.0 / 21 MB | — | 0.0 / 94 MB | 0.0 / 2.4 MB | 0.0 / 79 MB | 0.0 / 193 MB |
| l10000.n10000.e20 | 0.0 / 60 MB | — | 0.0 / 282 MB | 0.0 / 2.8 MB | 0.0 / 140 MB | 0.0 / 191 MB |
| l100000.n1000.e1 | 0.0 / 29 MB | — | 2.5 / 129 MB | 0.0 / 4.4 MB | 0.0 / 789 MB | — |
| l100000.n1000.e5 | 0.0 / 89 MB | — | 4.9 / 2671 MB | 0.0 / 9.3 MB | 29 / 2990 MB | — |
| l100000.n1000.e10 | 0.0 / 96 MB | — | — | 0.0 / 15 MB | 5468 / 5468 MB | — |
| l100000.n1000.e20 | 0.0 / 112 MB | — | — | 0.5 / 22 MB | — | — |
| l500000.n200.e1 | 1.0 / 91 MB | — | 88 / 3062 MB | 3.9 / 17 MB | — | — |
| l500000.n200.e5 | 114 / 114 MB | — | — | 33 / 33 MB | — | — |
| l500000.n200.e10 | — | — | — | — | — | — |
| l500000.n200.e20 | — | — | — | — | — | — |
| l1000000.n100.e1 | 7.5 / 101 MB | — | — | 6.7 / 26 MB | — | — |
| l1000000.n100.e5 | — | — | — | — | — | — |
| l1000000.n100.e10 | — | — | — | — | — | — |
| l1000000.n100.e20 | — | — | — | — | — | — |

Peak resident memory of each runner, its sample included, a run stopped mid-pair at least this
(QuickEd's Table 2 and S3):

| dataset | dinara-align (1 thread) | dinara-align (batch, 1 thread) | wfa | biwfa | ksw2 | bsalign |
| :-- | --: | --: | --: | --: | --: | --: |
| ont-1k | 15 MB | 17 MB | 11 MB | 8 MB | 8 MB | 4 MB |
| ont-10k | 49 MB | 50 MB | 163 MB | 10 MB | 114 MB | 184 MB |
| ont-50k | 114 MB | 117 MB | 3427 MB | 18 MB | — | 3365 MB |
| ont-500k | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| ont-500k-genvar | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| sars-cov-2 | 26 MB | 27 MB | 65 MB | 11 MB | 265 MB | 1701 MB |
| illumina-250 | 13 MB | 17 MB | 7 MB | 7 MB | 7 MB | 2 MB |
| pacbio-hifi | 16 MB | 28 MB | 19 MB | 9 MB | 57 MB | 998 MB |
| ont-promethion | 156 MB | 156 MB | o.o.m. | 32 MB | — | o.o.m. |
| l10000.n10000.e1 | 14 MB | 17 MB | 11 MB | 8 MB | 17 MB | 195 MB |
| l10000.n10000.e5 | 21 MB | 22 MB | 58 MB | 9 MB | 67 MB | 195 MB |
| l10000.n10000.e10 | 36 MB | 37 MB | 122 MB | 9 MB | 106 MB | 195 MB |
| l10000.n10000.e20 | 83 MB | 84 MB | 315 MB | 9 MB | ≥ 168 MB | 196 MB |
| l100000.n1000.e1 | 42 MB | 42 MB | 151 MB | 11 MB | 826 MB | o.o.m. |
| l100000.n1000.e5 | 101 MB | 101 MB | 2756 MB | 17 MB | — | o.o.m. |
| l100000.n1000.e10 | 116 MB | 116 MB | o.o.m. | 22 MB | — | o.o.m. |
| l100000.n1000.e20 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| l500000.n200.e1 | 107 MB | 107 MB | 3180 MB | 32 MB | o.o.m. | o.o.m. |
| l500000.n200.e5 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| l500000.n200.e10 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| l500000.n200.e20 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| l1000000.n100.e1 | 128 MB | 128 MB | o.o.m. | 45 MB | o.o.m. | o.o.m. |
| l1000000.n100.e5 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| l1000000.n100.e10 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| l1000000.n100.e20 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
