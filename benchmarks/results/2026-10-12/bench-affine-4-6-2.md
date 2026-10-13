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
| ont-1k | 1221 of 12477 | 102 µs | 103 µs | 189 µs | 306 µs | 1.11 ms | 342 µs | ✓ |
| ont-10k | 277 of 5000 | 1.57 ms | 1.6 ms | 4.86 ms | 5.86 ms | 30 ms | 7.85 ms | ✓ |
| ont-50k | 104 of 10000 | 21.3 ms | 21.4 ms | 63.4 ms | 54.9 ms | 340 ms (81/104) | 86.8 ms | ✓ |
| ont-500k | 4 of 50 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| ont-500k-genvar | 4 of 48 | 21.5 s (1/4) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| sars-cov-2 | 33 of 10000 | 674 µs | 650 µs | 2.06 ms | 1.78 ms | 88.6 ms | 341 ms | ✓ |
| illumina-250 | 4032 of 5992 | 1.16 µs | 1.88 µs | 1.57 µs | 2.82 µs | 12.8 µs | 54.3 µs | ✓ |
| pacbio-hifi | 69 of 890 | 78.7 µs | 810 µs | 227 µs | 233 µs | 9.73 ms | 85.9 ms | ✓ |
| ont-promethion | 22 of 584 | 261 ms | 263 ms | 285 ms (o.o.m. after 7/22) | 576 ms | 2.35 s (7/22) | o.o.m. | ✓ |
| l10000.n10000.e1 | 100 of 10000 | 186 µs | 209 µs | 408 µs | 648 µs | 21.3 ms | 37 ms | ✓ |
| l10000.n10000.e5 | 99 of 10000 | 2.29 ms | 2.32 ms | 7.1 ms | 8.72 ms | 97.8 ms | 37 ms | ✓ |
| l10000.n10000.e10 | 100 of 10000 | 6.96 ms | 6.97 ms | 22.8 ms | 27.3 ms | 204 ms | 37.6 ms | ✓ |
| l10000.n10000.e20 | 99 of 10000 | 19.7 ms | 19.8 ms | 65.5 ms | 75.2 ms | 375 ms (81/99) | 37.6 ms | ✓ |
| l100000.n1000.e1 | 9 of 1000 | 11.7 ms | 11.8 ms | 40 ms | 39.1 ms | 2.43 s | o.o.m. | ✓ |
| l100000.n1000.e5 | 9 of 1000 | 327 ms | 329 ms | 810 ms | 715 ms | 9.22 s (3/9) | o.o.m. | ✓ |
| l100000.n1000.e10 | 10 of 1000 | 1.17 s | 1.17 s | o.o.m. | 2.43 s | 17.4 s (1/10) | o.o.m. | ✓ |
| l100000.n1000.e20 | 10 of 1000 | 3.56 s (8/10) | t.o. (30 s) | o.o.m. | 7.35 s (4/10) | o.o.m. | o.o.m. | ✓ |
| l500000.n200.e1 | 4 of 200 | 463 ms | 465 ms | 1.29 s | 943 ms | o.o.m. | o.o.m. | ✓ |
| l500000.n200.e5 | 4 of 200 | 9.71 s (3/4) | t.o. (30 s) | o.o.m. | 19.4 s (1/4) | o.o.m. | o.o.m. | ✓ |
| l500000.n200.e10 | 4 of 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| l500000.n200.e20 | 4 of 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| l1000000.n100.e1 | 4 of 100 | 2.23 s | 2.22 s | o.o.m. | 3.91 s | o.o.m. | o.o.m. | ✓ |
| l1000000.n100.e5 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| l1000000.n100.e10 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |
| l1000000.n100.e20 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. | ✓ |

Median time per alignment over the same pairs (A*PA2's supplement):

| dataset | pairs | dinara-align (1 thread) | dinara-align (batch, 1 thread) | wfa | biwfa | ksw2 | bsalign |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| ont-1k | 1221 of 12477 | 109 µs | 103 µs | 201 µs | 279 µs | 1.04 ms | 400 µs |
| ont-10k | 277 of 5000 | 599 µs | 1.6 ms | 1.89 ms | 2.29 ms | 10.1 ms | 3.09 ms |
| ont-50k | 104 of 10000 | 1.72 ms | 21.4 ms | 5.45 ms | 6.92 ms | 44.8 ms (81/104) | 13.9 ms |
| ont-500k | 4 of 50 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| ont-500k-genvar | 4 of 48 | 21.5 s (1/4) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| sars-cov-2 | 33 of 10000 | 292 µs | 650 µs | 682 µs | 947 µs | 64 ms | 315 ms |
| illumina-250 | 4032 of 5992 | 0.817 µs | 1.88 µs | 1.21 µs | 1.97 µs | 12.1 µs | 54.2 µs |
| pacbio-hifi | 69 of 890 | 25.8 µs | 810 µs | 70.8 µs | 103 µs | 4.54 ms | 68.3 ms |
| ont-promethion | 22 of 584 | 62.5 ms | 263 ms | 101 ms (o.o.m. after 7/22) | 155 ms | 1.82 s (7/22) | o.o.m. |
| l10000.n10000.e1 | 100 of 10000 | 180 µs | 209 µs | 382 µs | 642 µs | 21.1 ms | 35.9 ms |
| l10000.n10000.e5 | 99 of 10000 | 2.26 ms | 2.32 ms | 6.87 ms | 8.72 ms | 97.3 ms | 35.9 ms |
| l10000.n10000.e10 | 100 of 10000 | 6.87 ms | 6.97 ms | 22.4 ms | 27.4 ms | 204 ms | 36.3 ms |
| l10000.n10000.e20 | 99 of 10000 | 19.2 ms | 19.8 ms | 64 ms | 75.1 ms | 374 ms (81/99) | 36.4 ms |
| l100000.n1000.e1 | 9 of 1000 | 10.1 ms | 11.8 ms | 33.7 ms | 38.7 ms | 2.43 s | o.o.m. |
| l100000.n1000.e5 | 9 of 1000 | 321 ms | 329 ms | 667 ms | 710 ms | 9.2 s (3/9) | o.o.m. |
| l100000.n1000.e10 | 10 of 1000 | 1.17 s | 1.17 s | o.o.m. | 2.44 s | 17.4 s (1/10) | o.o.m. |
| l100000.n1000.e20 | 10 of 1000 | 3.55 s (8/10) | t.o. (30 s) | o.o.m. | 7.34 s (4/10) | o.o.m. | o.o.m. |
| l500000.n200.e1 | 4 of 200 | 455 ms | 465 ms | 957 ms | 940 ms | o.o.m. | o.o.m. |
| l500000.n200.e5 | 4 of 200 | 9.7 s (3/4) | t.o. (30 s) | o.o.m. | 19.4 s (1/4) | o.o.m. | o.o.m. |
| l500000.n200.e10 | 4 of 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l500000.n200.e20 | 4 of 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l1000000.n100.e1 | 4 of 100 | 2.23 s | 2.22 s | o.o.m. | 3.91 s | o.o.m. | o.o.m. |
| l1000000.n100.e5 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l1000000.n100.e10 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l1000000.n100.e20 | 4 of 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |

Total over each whole dataset, its mean time per pair times its pairs, ≈ where the sample is not the
whole (QuickEd's Table 2 and S3):

| dataset | pairs | dinara-align (1 thread) | dinara-align (batch, 1 thread) | wfa | biwfa | ksw2 | bsalign |
| :-- | --: | --: | --: | --: | --: | --: | --: |
| ont-1k | 12477 | ≈ 1.27 s | ≈ 1.29 s | ≈ 2.36 s | ≈ 3.81 s | ≈ 13.8 s | ≈ 4.27 s |
| ont-10k | 5000 | ≈ 7.87 s | ≈ 7.99 s | ≈ 24.3 s | ≈ 29.3 s | ≈ 150 s | ≈ 39.3 s |
| ont-50k | 10000 | ≈ 213 s | ≈ 214 s | ≈ 634 s | ≈ 549 s | ≈ 3.4e+03 s | ≈ 868 s |
| ont-500k | 50 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| ont-500k-genvar | 48 | ≈ 1.03e+03 s | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| sars-cov-2 | 10000 | ≈ 6.74 s | ≈ 6.5 s | ≈ 20.6 s | ≈ 17.8 s | ≈ 886 s | ≈ 3.41e+03 s |
| illumina-250 | 5992 | ≈ 6.97 ms | ≈ 11.3 ms | ≈ 9.43 ms | ≈ 16.9 ms | ≈ 76.8 ms | ≈ 325 ms |
| pacbio-hifi | 890 | ≈ 70 ms | ≈ 721 ms | ≈ 202 ms | ≈ 208 ms | ≈ 8.66 s | ≈ 76.4 s |
| ont-promethion | 584 | ≈ 153 s | ≈ 154 s | o.o.m. | ≈ 336 s | ≈ 1.37e+03 s | o.o.m. |
| l10000.n10000.e1 | 10000 | ≈ 1.86 s | ≈ 2.09 s | ≈ 4.08 s | ≈ 6.48 s | ≈ 213 s | ≈ 370 s |
| l10000.n10000.e5 | 10000 | ≈ 22.9 s | ≈ 23.2 s | ≈ 71 s | ≈ 87.2 s | ≈ 978 s | ≈ 370 s |
| l10000.n10000.e10 | 10000 | ≈ 69.6 s | ≈ 69.7 s | ≈ 228 s | ≈ 273 s | ≈ 2.04e+03 s | ≈ 376 s |
| l10000.n10000.e20 | 10000 | ≈ 197 s | ≈ 198 s | ≈ 655 s | ≈ 752 s | ≈ 3.75e+03 s | ≈ 376 s |
| l100000.n1000.e1 | 1000 | ≈ 11.7 s | ≈ 11.8 s | ≈ 40 s | ≈ 39.1 s | ≈ 2.43e+03 s | o.o.m. |
| l100000.n1000.e5 | 1000 | ≈ 327 s | ≈ 329 s | ≈ 810 s | ≈ 715 s | ≈ 9.22e+03 s | o.o.m. |
| l100000.n1000.e10 | 1000 | ≈ 1.17e+03 s | ≈ 1.17e+03 s | o.o.m. | ≈ 2.43e+03 s | ≈ 1.74e+04 s | o.o.m. |
| l100000.n1000.e20 | 1000 | ≈ 3.56e+03 s | t.o. (30 s) | o.o.m. | ≈ 7.35e+03 s | o.o.m. | o.o.m. |
| l500000.n200.e1 | 200 | ≈ 92.7 s | ≈ 93.1 s | ≈ 257 s | ≈ 189 s | o.o.m. | o.o.m. |
| l500000.n200.e5 | 200 | ≈ 1.94e+03 s | t.o. (30 s) | o.o.m. | ≈ 3.88e+03 s | o.o.m. | o.o.m. |
| l500000.n200.e10 | 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l500000.n200.e20 | 200 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l1000000.n100.e1 | 100 | ≈ 223 s | ≈ 222 s | o.o.m. | ≈ 391 s | o.o.m. | o.o.m. |
| l1000000.n100.e5 | 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l1000000.n100.e10 | 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |
| l1000000.n100.e20 | 100 | t.o. (30 s) | t.o. (30 s) | o.o.m. | t.o. (30 s) | o.o.m. | o.o.m. |

Each rival's time over dinara-align (1 thread)'s, on the pairs both aligned:

| dataset | wfa | biwfa | ksw2 | bsalign |
| :-- | --: | --: | --: | --: |
| ont-1k | 1.9× | 3.0× | 10.9× | 3.4× |
| ont-10k | 3.1× | 3.7× | 19.0× | 5.0× |
| ont-50k | 3.0× | 2.6× | 19.5× | 4.1× |
| ont-500k | — | — | — | — |
| ont-500k-genvar | — | — | — | — |
| sars-cov-2 | 3.1× | 2.6× | 131× | 506× |
| illumina-250 | 1.4× | 2.4× | 11.0× | 46.6× |
| pacbio-hifi | 2.9× | 3.0× | 124× | 1091× |
| ont-promethion | 3.5× | 2.2× | 29.2× | — |
| l10000.n10000.e1 | 2.2× | 3.5× | 114× | 199× |
| l10000.n10000.e5 | 3.1× | 3.8× | 42.7× | 16.2× |
| l10000.n10000.e10 | 3.3× | 3.9× | 29.3× | 5.4× |
| l10000.n10000.e20 | 3.3× | 3.8× | 19.0× | 1.9× |
| l100000.n1000.e1 | 3.4× | 3.3× | 208× | — |
| l100000.n1000.e5 | 2.5× | 2.2× | 27.3× | — |
| l100000.n1000.e10 | — | 2.1× | 14.4× | — |
| l100000.n1000.e20 | — | 2.1× | — | — |
| l500000.n200.e1 | 2.8× | 2.0× | — | — |
| l500000.n200.e5 | — | 2.0× | — | — |
| l500000.n200.e10 | — | — | — | — |
| l500000.n200.e20 | — | — | — | — |
| l1000000.n100.e1 | — | 1.8× | — | — |
| l1000000.n100.e5 | — | — | — | — |
| l1000000.n100.e10 | — | — | — | — |
| l1000000.n100.e20 | — | — | — | — |

The same across datasets, as QuickEd's text summarizes its own:

| rival | least | geometric mean | most | datasets |
| :-- | --: | --: | --: | --: |
| wfa | 1.4× | 2.7× | 3.5× | 14 |
| biwfa | 1.8× | 2.7× | 3.9× | 18 |
| ksw2 | 10.9× | 35.3× | 208.1× | 14 |
| bsalign | 1.9× | 23.8× | 1091.2× | 10 |

Growth of peak memory aligning a pair, median / largest, input and runtime left out (A*PA2's Table 10):

| dataset | dinara-align (1 thread) | dinara-align (batch, 1 thread) | wfa | biwfa | ksw2 | bsalign |
| :-- | --: | --: | --: | --: | --: | --: |
| ont-1k | 0.0 / 0.4 MB | — | 0.0 / 3.0 MB | 0.0 / 1.7 MB | 0.0 / 0.9 MB | 0.0 / 1.9 MB |
| ont-10k | 0.0 / 22 MB | — | 0.0 / 97 MB | 0.0 / 2.2 MB | 0.0 / 58 MB | 0.0 / 104 MB |
| ont-50k | 0.0 / 85 MB | — | 0.0 / 1420 MB | 0.0 / 4.5 MB | 0.0 / 317 MB | 0.0 / 1303 MB |
| ont-500k | — | — | — | — | — | — |
| ont-500k-genvar | 143 / 143 MB | — | — | — | — | — |
| sars-cov-2 | 0.0 / 11 MB | — | 0.0 / 51 MB | 0.0 / 2.5 MB | 0.0 / 211 MB | 0.0 / 1682 MB |
| illumina-250 | 0.0 / 0.1 MB | — | 0.0 / 0.1 MB | 0.0 / 0.1 MB | 0.0 / 0.0 MB | 0.0 / 0.2 MB |
| pacbio-hifi | 0.0 / 1.6 MB | — | 0.0 / 7.1 MB | 0.0 / 1.5 MB | 0.0 / 26 MB | 0.0 / 258 MB |
| ont-promethion | 0.0 / 89 MB | — | 1.8 / 1768 MB | 0.0 / 8.9 MB | 0.0 / 2155 MB | — |
| l10000.n10000.e1 | 0.0 / 0.5 MB | — | 0.0 / 4.0 MB | 0.0 / 1.5 MB | 0.0 / 11 MB | 0.0 / 192 MB |
| l10000.n10000.e5 | 0.0 / 6.8 MB | — | 0.0 / 30 MB | 0.0 / 2.0 MB | 0.0 / 40 MB | 0.0 / 192 MB |
| l10000.n10000.e10 | 0.0 / 21 MB | — | 0.0 / 94 MB | 0.0 / 2.2 MB | 0.0 / 79 MB | 0.0 / 193 MB |
| l10000.n10000.e20 | 0.0 / 60 MB | — | 0.0 / 282 MB | 0.0 / 2.8 MB | 0.0 / 140 MB | 0.0 / 191 MB |
| l100000.n1000.e1 | 0.0 / 29 MB | — | 2.5 / 129 MB | 0.0 / 4.2 MB | 0.0 / 789 MB | — |
| l100000.n1000.e5 | 0.0 / 88 MB | — | 4.9 / 2671 MB | 0.0 / 9.1 MB | 29 / 2990 MB | — |
| l100000.n1000.e10 | 0.0 / 96 MB | — | — | 0.0 / 15 MB | 5468 / 5468 MB | — |
| l100000.n1000.e20 | 0.0 / 111 MB | — | — | 0.5 / 22 MB | — | — |
| l500000.n200.e1 | 0.5 / 91 MB | — | 88 / 3062 MB | 3.9 / 17 MB | — | — |
| l500000.n200.e5 | 0.0 / 113 MB | — | — | 33 / 33 MB | — | — |
| l500000.n200.e10 | — | — | — | — | — | — |
| l500000.n200.e20 | — | — | — | — | — | — |
| l1000000.n100.e1 | 7.6 / 100 MB | — | — | 6.7 / 26 MB | — | — |
| l1000000.n100.e5 | — | — | — | — | — | — |
| l1000000.n100.e10 | — | — | — | — | — | — |
| l1000000.n100.e20 | — | — | — | — | — | — |

Peak resident memory of each runner, its sample included, a run stopped mid-pair at least this
(QuickEd's Table 2 and S3):

| dataset | dinara-align (1 thread) | dinara-align (batch, 1 thread) | wfa | biwfa | ksw2 | bsalign |
| :-- | --: | --: | --: | --: | --: | --: |
| ont-1k | 15 MB | 18 MB | 11 MB | 8 MB | 8 MB | 4 MB |
| ont-10k | 49 MB | 51 MB | 163 MB | 10 MB | 114 MB | 184 MB |
| ont-50k | 114 MB | 117 MB | 3427 MB | 18 MB | — | 3365 MB |
| ont-500k | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| ont-500k-genvar | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| sars-cov-2 | 26 MB | 27 MB | 65 MB | 11 MB | 265 MB | 1701 MB |
| illumina-250 | 13 MB | 17 MB | 7 MB | 7 MB | 7 MB | 2 MB |
| pacbio-hifi | 16 MB | 28 MB | 19 MB | 9 MB | 57 MB | 998 MB |
| ont-promethion | 156 MB | 157 MB | o.o.m. | 32 MB | — | o.o.m. |
| l10000.n10000.e1 | 14 MB | 18 MB | 11 MB | 8 MB | 17 MB | 195 MB |
| l10000.n10000.e5 | 20 MB | 22 MB | 58 MB | 8 MB | 66 MB | 195 MB |
| l10000.n10000.e10 | 36 MB | 37 MB | 122 MB | 9 MB | 106 MB | 195 MB |
| l10000.n10000.e20 | 83 MB | 85 MB | 315 MB | 9 MB | ≥ 168 MB | 196 MB |
| l100000.n1000.e1 | 42 MB | 42 MB | 152 MB | 11 MB | 826 MB | o.o.m. |
| l100000.n1000.e5 | 101 MB | 101 MB | 2756 MB | 17 MB | — | o.o.m. |
| l100000.n1000.e10 | 116 MB | 116 MB | o.o.m. | 22 MB | — | o.o.m. |
| l100000.n1000.e20 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| l500000.n200.e1 | 106 MB | 106 MB | 3180 MB | 32 MB | o.o.m. | o.o.m. |
| l500000.n200.e5 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| l500000.n200.e10 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| l500000.n200.e20 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| l1000000.n100.e1 | 128 MB | 128 MB | o.o.m. | 45 MB | o.o.m. | o.o.m. |
| l1000000.n100.e5 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| l1000000.n100.e10 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
| l1000000.n100.e20 | — | — | o.o.m. | — | o.o.m. | o.o.m. |
