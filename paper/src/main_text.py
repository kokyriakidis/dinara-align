"""The main paper, written once; render.py turns it into paper.html and dinara-align.tex.

Inline markup: **bold**, *italic*, `code`, [@Key; @Key] parenthetical citation, [@!Key] narrative citation,
\\( \\) inline maths, \\[ \\] display maths, [#tab:x] / [#fig:x] cross-references.
Blocks: ("sec", title), ("sub", title), ("p", text), ("list", [items]) numbered, ("table", key), ("figure", key).
"""

TITLE = 'dinara-align: hardware-aware algorithms for exact pairwise alignment'
# Authors in order, each with an email; both are corresponding authors.
AUTHORS = [("Konstantinos Kyriakidis", "kkyriaki@ucsc.edu"), ("Benedict Paten", "bpaten@ucsc.edu")]
AFFILIATION = "UC Santa Cruz Genomics Institute, Santa Cruz, CA, USA"

ABSTRACT = [('Motivation',
  'Exact pairwise alignment of long sequencing reads remains computationally demanding. Band doubling over '
  'bit-parallel blocks, the wavefront algorithm and A* with seed heuristics have greatly reduced the number of cells '
  'an exact aligner must compute. A large part of the remaining running time depends on how the computation maps '
  'onto the processor: loop-carried dependency chains, memory-bound table lookups, and bounds and thresholds that '
  'are poorly matched to the input or the hardware.'),
 ('Results',
  'We present dinara-align, an exact pairwise aligner built on algorithmic changes that address these bottlenecks: a '
  "regrouping of Myers' bit-parallel recurrence that shortens its loop-carried dependency chain, band bounds "
  'accepted only when independent estimates agree, diagonal transition vectorised across diagonals, a direct '
  'construction of the neighbours of exact seed matches that provably preserves local pruning, and per-pair '
  'optimality certificates for banded inter-sequence batches. On an Intel Skylake-X processor, dinara-align had the '
  'lowest mean running time among the exact aligners tested in 33 of the 34 configurations of the A*PA2 benchmark, '
  'aligning ultra-long nanopore reads in 104 ms on average, compared with 176 ms for A*PA2-full. It was 1.7 to 3.2 '
  'times faster than WFA at gap-affine costs and computed global scores for 500,000 short-read pairs on one core 4.5 '
  'times faster than WFA2-lib.'),
 ('Availability and implementation',
  'dinara-align is free software under the Mozilla Public License 2.0, available at '
  'https://github.com/kokyriakidis/dinara-align with C, C++, Python and command-line interfaces.'),
 ('Contact', 'kkyriaki@ucsc.edu, bpaten@ucsc.edu'),
 ('Supplementary information', 'Supplementary data are available at *Bioinformatics* online.')]

BODY = [
('sec', 'Introduction'),
('p',
 'Pairwise alignment determines the minimum-cost series of substitutions, insertions and deletions that transforms '
 'one sequence into another. The classical dynamic programming algorithm computes an \\(n \\times m\\) matrix '
 '[@Gotoh1982], which is practical for short reads but not for the hundreds of kilobases spanned by current nanopore '
 'reads. Heuristic methods restrict the computed region and lose the guarantee of optimality; exact methods instead '
 'avoid only those cells that provably cannot lie on an optimal path.'),
('p',
 'Three ideas underlie current exact aligners. *Band doubling* computes only the diagonals that a path within a '
 'given cost bound can reach, and doubles the bound until the optimal cost lies within it [@Ukkonen1985]. '
 '*Bit-parallelism* encodes a column of 64 cells in one machine word [@Myers1999], and Edlib combines it with band '
 'doubling [@Sosic2017]. *Diagonal transition* extends the furthest-reaching paths one cost at a time, following '
 'runs of matches at no cost [@Ukkonen1985; @Myers1986]; WFA generalised it to gap-affine costs [@MarcoSola2021] and '
 'BiWFA reduced its memory to linear in the distance [@MarcoSola2023]. A*PA introduced A* search with a gap-chaining '
 'seed heuristic [@GrootKoerkamp2024a], and A*PA2 combined this heuristic with band doubling over SIMD blocks of '
 'bit-parallel words, aligning ultra-long reads up to 19 times faster than Edlib and BiWFA [@GrootKoerkamp2024b]. A '
 'second line of work accelerates the full dynamic program on vector hardware, using anti-diagonal [@Wozniak1997], '
 'striped [@Farrar2007; @Zhao2013] and prefix-scan [@Daily2016] layouts, difference recurrences [@Suzuki2018; '
 '@Li2018], adaptive blocks [@Liu2023] and inter-sequence parallelism [@Rognes2011; @Rahn2018]; these methods '
 'compute the full matrix or a fixed band, so their running time grows with its size rather than with the distance.'),
('p',
 'In the exact methods, a large part of the remaining running time is determined by how the computation maps onto '
 'the processor. The bit-parallel sweep carries a dependency from each column to the next, so its speed is set by '
 'the latency of that chain rather than by its number of operations. Band doubling spends whole rounds when its '
 'bound is far from the distance. The seed heuristic is built by scalar lookups in tables larger than the caches. '
 "Batches of short pairs are vectorised across pairs, but each pair's full matrix is computed. Here we present "
 'dinara-align, an exact aligner whose algorithms address these bottlenecks. Our contributions are:'),
('list',
 ["**A shorter critical path for Myers' recurrence.** Two Boolean identities remove two operations from the "
  'loop-carried dependency chain of the bit-parallel step; on AVX-512, two staggered groups of words are '
  'interleaved, extending the instruction-level parallelism of A*PA2.',
  '**Band bounds from agreeing estimates.** A band-doubling bound is extrapolated from partial progress only when '
  'two independent estimates agree.',
  '**Gathered diagonal transition.** Diagonal transition and the gap-affine wavefront extend eight diagonals per '
  'step using gather instructions.',
  '**Direct neighbour construction for the seed heuristic.** The inexact neighbours of an exact seed match are '
  'inserted without search, and we prove that this preserves the outcome of local pruning (Lemma 1).',
  '**Certified bands for inter-sequence batches.** Each short pair is aligned in one vector lane within a narrow '
  'band whose optimality is certified per lane by a lower bound on the cost of any path that leaves it (Proposition '
  '1), adapting the speculate-and-verify approach of the SeedEx accelerator [@Fujiki2020] to software batches.']),
('sec', 'Methods'),
('sub', 'Problem definition and overview'),
('p',
 'Let \\(A=a_1\\dots a_n\\) be the reference and \\(B=b_1\\dots b_m\\) the query, and let cell \\((i,j)\\) lie on '
 'diagonal \\(j-i\\). Under gap-affine costs, a mismatch costs \\(x\\) and a gap of \\(k\\) letters costs '
 '\\(g(k)=o+ke\\); unit costs correspond to \\(x=e=1\\) and \\(o=0\\), and two-piece costs take the minimum of two '
 'such functions [@Li2018]. A substitution matrix instead assigns each pair of letters a score to be maximised. We '
 'consider global, ends-free, extension and local alignment. An alignment is *exact* if it is optimal under the '
 'requested costs, mode, band and cost limit; among optimal alignments, a fixed rule places indels leftmost or '
 'rightmost, so that the result does not depend on which algorithm computed the cost.'),
('p',
 'Each problem is assigned to an algorithm ([#tab:engines]). For unit costs we follow A*PA2 [@GrootKoerkamp2024b]: '
 'band doubling over tiles of 256 columns of bit-parallel words, pruned by the gap-chaining seed heuristic '
 '[@GrootKoerkamp2024a], with a diagonal-transition traceback inside each tile; near-identical pairs are aligned by '
 'diagonal transition alone. Gap-affine costs use the bidirectional wavefront of WFA and BiWFA [@MarcoSola2021; '
 '@MarcoSola2023]. A short diagonal-transition search estimates the distance, and the estimate, the length and the '
 'cost model select the algorithm. The thresholds were set from measurements on each instruction set and are fixed '
 'at compile time (Supplementary Section S1).'),
('table', 'engines'),
('sub', "A shorter critical path for Myers' recurrence"),
('p',
 "As in Edlib and A*PA2 [@Sosic2017; @GrootKoerkamp2024b], the query is packed 64 bases per word, and Myers' "
 'recurrence [@Myers1999] advances a word by one column. With match mask \\(E\\), vertical difference masks \\(P_v, '
 "M_v\\) and the incoming negative horizontal difference folded into \\(E' = E \\lor M_h^{\\mathrm{in}}\\), A*PA2 "
 'computes'),
('p',
 "\\[ c = (E' \\land P_v) + P_v, \\qquad X_h = (c \\oplus P_v) \\lor E', \\qquad P_h = M_v \\lor \\lnot(X_h \\lor "
 'P_v), \\qquad M_h = P_v \\land X_h. \\]'),
('p',
 'The horizontal masks \\(P_h\\) and \\(M_h\\) determine the vertical masks of the next column, so the path from '
 '\\(c\\) to them lies on the loop-carried chain. Since \\((c \\oplus P_v) \\lor P_v = c \\lor P_v\\) and \\(P_v \\land (c \\oplus P_v) = P_v \\land '
 '\\lnot c\\), the same masks are obtained as'),
('p',
 "\\[ P_h = M_v \\lor \\lnot\\bigl(c \\lor (E' \\lor P_v)\\bigr), \\qquad M_h = (P_v \\land \\lnot c) \\lor (P_v "
 "\\land E'), \\]"),
('p',
 "where \\(E' \\lor P_v\\) and \\(P_v \\land E'\\) do not depend on \\(c\\) and are computed off the critical path. "
 'This removes two operations from the loop-carried chain and, for a single word, reduces the cost from 14 to 8 '
 'cycles per column on the Skylake-X. Myers [@Myers1999], Edlib and A*PA2 use the original form, and we are not '
 'aware of a previous description of the regrouped one.'),
('p',
 'Words occupy the lanes of a vector, staggered by one column, so that the horizontal difference leaving one word '
 'enters the next through a lane rotation [@Wozniak1997; @GrootKoerkamp2024b]. Every column then waits for the '
 'rotation of the previous one, and A*PA2 processes two four-lane vectors together to overlap these waits. On '
 'AVX-512, where a vector holds eight words, dinara-align interleaves two eight-word groups, the lower lagging ten '
 'columns behind the upper, so that the instructions of each group fill the latency of the other; with narrower '
 'vectors a single group is used. Bands of two or three words, in which every operation lies on the critical path, are '
 'computed in general-purpose registers, whose operations have lower latency. On AVX-512 the sweep requires about '
 '2.4 cycles per word and column.'),
('sub', 'Band bounds from agreeing estimates'),
('p',
 'Band doubling returns the exact distance for any initial bound, but its running time depends on how close that '
 'bound is to the distance. A partial computation that has reached cost \\(s\\) at anti-diagonal \\(r\\) suggests '
 'the extrapolated distance \\(\\hat d = s(n+m)/r\\). This extrapolation is reliable when edits are spread evenly '
 'along the pair and can be too high by orders of magnitude when they are concentrated, as near the ends of nanopore '
 'reads. dinara-align therefore computes two extrapolations, \\(\\hat d_1\\) from the first eight edits of a '
 'diagonal-transition search and \\(\\hat d_2\\) from the same search continued until its budget, and sets the first '
 'bound from them only if \\(\\max(\\hat d_1,\\hat d_2) \\le \\tfrac43 \\min(\\hat d_1,\\hat d_2)\\). Otherwise, as '
 'in A*PA2, the bound starts 256 above the known lower bound and doubles. Within a round, the distance is '
 'extrapolated again after one eighth, one quarter and one half of the columns; the bound is lowered if the '
 'extrapolation lies below it, and the round is abandoned if the extrapolation exceeds it. Because these rules '
 'change only the sequence of bounds, the result remains exact.'),
('sub', 'Gathered diagonal transition'),
('p',
 'Diagonal transition extends each diagonal of a wavefront along its run of matches. dinara-align compares eight '
 'letters at a time by the exclusive-or of two 64-bit words and a count of trailing zeros. On AVX-512, eight '
 'diagonals are extended together: their letters are loaded with gather instructions, and the wavefront is padded '
 'with unreachable diagonals so that every step processes full vectors under branch-free masks. The gather masks are '
 'derived from the lanes of each group, which avoids a dependency between consecutive gathers that would otherwise '
 'serialise them. The same vectorised step serves near-identical pairs, computed from one or both ends '
 '[@MarcoSola2023], the traceback inside each tile, restricted to the diagonals that can still reach the traced '
 'cell, and the three-layer gap-affine wavefront.'),
('sub', 'Seed heuristic: direct neighbour construction'),
('p',
 'The reference is divided into disjoint seeds of 12 bases, matched exactly, or of 16 bases, matched with at most '
 'one edit [@GrootKoerkamp2024a]. A 16-base seed that matches with at most one edit matches one of its two 8-base '
 'halves exactly, so each half indexes a table. Since these tables exceed the caches, the lookups for 16 query rows '
 'are issued together, so that their memory accesses overlap, and only seeds whose matches could still chain to the '
 'end of the pair are considered.'),
('p',
 'An exact match \\(\\mu\\) of a seed at query rows \\([i, i+16)\\) always has four neighbouring matches with one '
 'edit: two share its start, at rows \\([i, i+15)\\) and \\([i, i+17)\\), and two share its end, at \\([i+1, '
 'i+16)\\) and \\([i-1, i+16)\\). Because chaining is all-or-nothing, a neighbour one diagonal away may chain where '
 '\\(\\mu\\) cannot, and the neighbours are required for the heuristic to remain a lower bound. Rather than locating '
 'them by search and testing each for local pruning, dinara-align inserts them whenever \\(\\mu\\) is retained, '
 'which is justified by the following lemma. Local pruning retains a match \\(\\nu\\) of a seed only if \\(c(\\nu) + '
 'h(e(\\nu)) < \\tau\\), where \\(c(\\nu)\\) is its cost, \\(e(\\nu)\\) its end cell, \\(h(u)\\) the least cost of a '
 'path from \\(u\\) across the following seeds, and \\(\\tau\\) a threshold that depends only on the seeds.'),
('p', '**Lemma 1.** If an exact match \\(\\mu\\) fails local pruning, so does each of its four neighbours.'),
('p',
 '*Proof.* A neighbour \\(\\nu\\) that shares the end of \\(\\mu\\) has \\(e(\\nu)=e(\\mu)\\) and '
 '\\(c(\\nu)=1>c(\\mu)=0\\). A neighbour that shares the start ends on an adjacent diagonal, and one insertion or '
 'deletion from \\(e(\\mu)\\) reaches a cell on that diagonal at or beyond \\(e(\\nu)\\). The least remaining cost '
 'does not increase along a diagonal [@Ukkonen1985], so \\(h(e(\\mu)) \\le 1 + h(e(\\nu))\\) and '
 '\\(c(\\mu)+h(e(\\mu)) \\le c(\\nu)+h(e(\\nu))\\). In both cases the pruning test value of \\(\\nu\\) is at least '
 'that of \\(\\mu\\). \\(\\square\\)'),
('p',
 'Neighbours sharing the start therefore raise the layer of that start by the chain length available from their '
 'ends, and neighbours sharing the end are placed one layer below it, without a search or a pruning test of their '
 'own.'),
('sub', 'Certified bands for inter-sequence batches'),
('p',
 'Batches of short pairs are aligned with one pair per vector lane [@Rognes2011; @Rahn2018], first in saturating '
 '8-bit lanes and, for pairs that saturate, in 16-bit lanes [@Zhao2013; @Daily2016]. Rather than computing each '
 'matrix in full, a lane computes the band of diagonals \\([d_\\ell, d_h]\\) containing the main diagonal 0 and the '
 'end diagonal \\(\\delta = m - n\\), together with one further diagonal on either side. Its optimality is '
 'certified with the following bound, a cost form of the out-of-band bounds of [@!Gibrat2018] and '
 '[@!Fujiki2020], stated for gap-affine and two-piece costs and with a strict form for tie-breaking.'),
('p',
 '**Proposition 1.** Let \\(C\\) be the least cost of an alignment within \\([d_\\ell, d_h]\\), and let \\(L(d) = '
 'g(|d|) + g(|d-\\delta|)\\). If \\(C \\le \\min\\bigl(L(d_h+1), L(d_\\ell-1)\\bigr)\\), then \\(C\\) is the optimal '
 'cost; if the inequality is strict, every optimal alignment lies within the band.'),
('p',
 '*Proof.* An alignment that leaves the band visits diagonal \\(d_h+1\\) or \\(d_\\ell-1\\). Consider \\(d=d_h+1 > '
 '\\max(0,\\delta)\\); the case below the band is symmetric. Starting on diagonal 0, the alignment contains '
 'insertions totalling at least \\(d\\) letters to reach \\(d\\), and deletions totalling at least \\(d-\\delta\\) '
 'letters to end on \\(\\delta\\). These form at least one gap of each kind, and since \\(g\\) is subadditive (for '
 'two-piece costs, with the usual convention that the longer piece has the lower extension cost), their total cost '
 'is at least \\(g(d)+g(d-\\delta)=L(d)\\). Every alignment outside the band therefore costs at least '
 '\\(\\min\\bigl(L(d_h+1), L(d_\\ell-1)\\bigr)\\). \\(\\square\\)'),
('p',
 'Bounds of this kind underlie band doubling [@Ukkonen1985] and have been used to verify narrow-band results in '
 'software [@Gibrat2018] and in hardware [@Fujiki2020]; here the certificate is evaluated independently in each '
 'vector lane of an inter-sequence batch. For alignments, the strict form guarantees that the alignment selected by '
 'the tie rule lies within the band. Pairs that fail the certificate, '
 'about one in seven in our short-read benchmark, are recomputed over every diagonal that a cheaper path could '
 'visit, grouped by the width of that band so that the lanes of a group finish together.'),
('sub', 'Gap-affine and score-based alignment'),
('p',
 'A scoring scheme with match score \\(a\\), mismatch score \\(b\\) and gap scores \\(-(o+ke)\\) can be converted '
 'into costs for global alignment [@Eizenga2022]. Every alignment covers \\(n+m\\) letters, so its score \\(S\\) and '
 'a gap-affine cost \\(C\\) satisfy'),
('p',
 '\\[ S = \\tfrac12\\bigl(a(n+m) - C\\bigr), \\qquad C = 2(a-b)\\,X + \\sum_{\\text{gaps}} \\bigl(2o + '
 'k(2e+a)\\bigr), \\]'),
('p',
 'where \\(X\\) is the number of mismatches, so the wavefront computes score-maximising alignments. The cost also '
 'bounds the band in which the alignment must be stored. Each gapped letter contributes at least \\(2e+a\\) to '
 '\\(C\\), and an alignment that deviates by \\(t\\) diagonals from both the main and the end diagonal contains at '
 'least \\(2t\\) gapped letters; hence only diagonals within \\(C/(2(2e+a))\\) of them can contain an optimal '
 'alignment. For a 1 kbp pair at 5% divergence, this stores 137 diagonals instead of 1,000. Under a general '
 'substitution matrix, anti-diagonals are computed in 16-bit vector lanes while the scores provably fit '
 '[@Wozniak1997]. Local alignment first locates the end of the optimal alignment with a forward sweep, as in SSW '
 '[@Zhao2013], and then extends a wavefront backwards from that end, stopping at the first cell that attains the '
 'maximal score; this avoids the correction pass that striped layouts may need when long alignments score highly '
 '[@Farrar2007].'),
('sub', 'Validation'),
('p',
 'A differential fuzzer compares results for random pairs, cost models, modes, bands, cost limits and tie-breaking '
 'rules with a full-matrix reference implementation that shares no code with the library, checking validity, '
 'optimality, the symmetry of the two tie-breaking rules and agreement between batched and single calls. The '
 'WFA2-lib regression set is reproduced exactly, including all CIGAR strings under the rightmost tie-breaking rule. '
 'The seed heuristic is compared with the true remaining cost at every cell of random pairs; this check found and '
 'confirmed the correction of an implementation error that caused occasional overestimates (Supplementary Section '
 'S1).'),
('sec', 'Results'),
('sub', 'Experimental setup'),
('p',
 'All measurements were taken on an Intel Core i9-7900X (Skylake-X, AVX-512) with ten cores fixed at 3.3 GHz, turbo '
 'boost and hyper-threading disabled, and 91 GB of memory, running Ubuntu 26.04. dinara-align is implemented in Mojo '
 '[@Modular] and was compiled with Mojo 1.1.0 at commit `ebe290f`; all tools were compiled for the native '
 'instruction set. The exact aligners compared are A*PA2-full and A*PA2-simple [@GrootKoerkamp2024b], A*PA '
 '[@GrootKoerkamp2024a], Edlib [@Sosic2017], WFA and BiWFA [@MarcoSola2021; @MarcoSola2023], KSW2 [@Li2018; '
 '@Suzuki2018], parasail [@Daily2016], SSW [@Zhao2013], abPOA [@Gao2021] and hyalite [@hyalite]. Every cost was '
 'compared with those of every other tool and with the costs recorded in the published A*PA2 results; no '
 'discrepancies occurred.'),
('p',
 'For unit costs we followed the A*PA2 protocol [@GrootKoerkamp2024b]: one single-threaded job at a time, pinned to '
 'one core, each pair aligned once with traceback, and a time limit of 20 s per sample, after which only the '
 'completed pairs are reported. The datasets are those of A*PA2: nanopore reads of about 1, 10 and 50 kbp, '
 'ultra-long reads over 500 kbp with and without genetic variation, SARS-CoV-2 genomes, and uniform-error pairs '
 'regenerated from the original random seed. We report the mean wall-clock time per alignment, excluding input, and '
 'the increase in peak resident memory during an alignment, measured with `getrusage`. As in A*PA2, each sample was '
 'run once; the one configuration in which dinara-align was not the fastest was repeated five times. Other workloads '
 'were timed in-process after a warm-up, as the fastest of twenty batches of about 10 ms, or once for calls of 0.1 s '
 'or longer; short-read batches as the faster of two passes.'),
('sub', 'Edit distance'),
('p',
 'On the real datasets, dinara-align had the lowest mean and median time among the exact aligners ([#tab:real]; '
 'Supplementary Table S4). It aligned the 1 kbp reads in 23 µs, compared with 31 µs for WFA, required about two '
 'thirds of the time of A*PA2-simple on the 10 and 50 kbp reads, and aligned the ultra-long reads in 104 and 144 ms, '
 'compared with 176 and 210 ms for A*PA2-full. The approximate aligners WFA-adaptive [@MarcoSola2021] and Block '
 'Aligner [@Liu2023] were slower than dinara-align on every dataset.'),
('table', 'real'),
('p',
 'On uniform 100 kbp pairs ([#fig:sweeps]a), dinara-align required 4.4 to 6.8 ms between 2 and 8% divergence and '
 '10.4 to 12.3 ms from 9 to 15%, where A*PA2-full required up to 18 ms. A*PA2-full was faster only at 10% divergence '
 '(10.5 versus 11.0 ms; in five repeated runs, 10.43–10.49 versus 11.02–11.06 ms). With increasing length '
 '([#fig:sweeps]b), dinara-align aligned 1 Mbp pairs at 15% divergence in 265 ms, compared with 1.33 s for '
 'A*PA2-full and 2.07 s for A*PA.'),
('figure', 'sweeps'),
('p',
 'On the ultra-long reads, an alignment increased peak memory by 51 MB without and 37 MB with genetic variation, '
 'compared with 81 and 48 MB for A*PA2-full and 82 and 32 MB for A*PA2-simple; on 1 Mbp pairs at 15% divergence, by '
 '55 MB, compared with 86 to 149 MB for A*PA and A*PA2. On shorter inputs dinara-align generally used more memory '
 'than the most economical tool, for example 13 MB against 2.3 MB for Edlib on 300 kbp pairs at 5% divergence, and '
 'its runtime occupies about 13 MB, compared with 3 to 5 MB for the other tools.'),
('sub', 'Gap-affine costs, other modes and batches'),
('p',
 'At the WFA costs (4, 6, 2), dinara-align was 1.7 to 3.2 times faster than WFA on every dataset that any tool '
 'completed, and was the only aligner to complete all 104 reads of the 50 kbp dataset within the time limit '
 '([#tab:affine]). In each of the 14 workloads covering local, overlap, infix, prefix, extension, two-piece and '
 'substitution-matrix alignment, it was faster than the fastest alternative, by factors from 1.1 to 10.6 '
 '(Supplementary Table S2).'),
('table', 'affine'),
('p',
 'Following the short-read case study of Accelign [@Kallenborn2026], we computed global scores for 500,000 pairs of '
 '148 bp reads with 1% error against reference segments of similar length ([#tab:batch]). Under gap-affine costs '
 '(mismatch 1, gap \\(2+k\\)), dinara-align required 0.47 µs per pair on one core, compared with 2.12 µs for '
 'WFA2-lib, and 26 ms for the whole batch on ten cores, compared with 115 ms.'),
('table', 'batch'),
('sub', 'Effect of each technique'),
('p',
 'To measure what each technique contributes, dinara-align was built with one technique switched off at a time, '
 'everything else unchanged, and every build aligned the A*PA2 samples with traceback and scored the short-read '
 'batch. The builds alternated within each of three rounds, pinned to one core, and we report the median against '
 'the build with every technique on; the baseline varied by at most 1.4% between rounds, and every build returned '
 'the same costs on every dataset ([#tab:ablation]). Constructing the neighbours of exact seed matches directly '
 'contributed most on long, divergent pairs: without it, the ultra-long reads took 11% and 19% longer and the '
 '100 kbp pairs at 10 and 15% divergence 13% longer. Interleaving two groups of bit-parallel words and the regrouped '
 'recurrence each contributed 4 to 6% on the ultra-long reads and on the most divergent uniform pairs. Gathered '
 'diagonal transition contributed most where diagonal transition does most of the work, on the 1 kbp reads, the '
 'SARS-CoV-2 genomes and near-identical 100 kbp pairs (16 to 25%). The agreement rule for extrapolated bounds '
 'mattered on the 10 kbp reads (17%) and the SARS-CoV-2 genomes (5%), and not on the ultra-long reads, where the '
 'seed heuristic sets the first bound. Without the certified band, the batch took 2.9 and 3.1 times as long at '
 'affine and unit costs. Effects below about 6% are of the order of the variation between builds that we attribute '
 'to code alignment (Section 4): a second, independent set of builds reproduced every effect above 15% but gave 1% '
 'for the regrouped recurrence and 3% for the paired groups on the ultra-long reads, and 6 to 8% for the seed '
 'neighbours on the 100 kbp pairs.'),
('table', 'ablation'),
('sec', 'Discussion'),
('p',
 'The algorithms in dinara-align build on those of earlier exact aligners; the gains come from changes to how these '
 'algorithms use the processor. Two observations stand out. First, on long pairs the techniques that reduce work '
 'outside the inner loop contributed more than those that speed it up: constructing seed neighbours directly saved '
 'up to 19%, against 4 to 6% each for the two changes to the bit-parallel sweep ([#tab:ablation]). Second, as the '
 'inner loops become faster, the scalar, memory-bound construction of the seed heuristic dominates: on 100 kbp pairs at 10% divergence it accounts for 58% of the '
 'running time and is the main cost in the one configuration in which A*PA2-full remained faster. Neither the '
 'regrouped recurrence nor the per-lane band certificate depends on the implementation, and both apply directly to '
 'other bit-parallel and inter-sequence aligners.'),
('p',
 'The study has several limitations. First, the evaluation covers a single x86 processor with AVX-512. The AVX2 '
 'configuration was measured only for its inexact-seed threshold, built for Haswell on the same machine (20 and 30% '
 'equal within 0.5%, 40% up to 6% slower), and the NEON configuration, which the library also supports, is not '
 'evaluated here. Second, some thresholds were chosen using workloads that also appear in the '
 'evaluation; in particular, the AVX-512 threshold for inexact seeds was lowered from 40 to 20% after measurements '
 'on the uniform 100 kbp pairs of [#fig:sweeps]a, so results on those pairs may overstate performance on unseen '
 'data. Third, as in the A*PA2 protocol, most configurations were run once, and on the Skylake-X some running times '
 'varied by up to one third between builds that did not change the code involved, which we attribute to code '
 'alignment effects in the instruction decoder. Fourth, the ablation switches off one technique at a time against '
 'the full configuration, so it does not separate interactions between techniques. Finally, release builds target the '
 'baseline processor of each platform, so the AVX-512 code paths require a build for the host processor.'),
('p',
 'Future work includes speeding up seed construction, reusing computation across band-doubling rounds as in the '
 'incremental doubling of A*PA2, measurements on further processors, and GPU implementations of the same algorithms.'),
]

BACK = [
    ("Funding", "This work received no specific funding."),
    ("Conflict of interest", "None declared."),
    ("Data availability", "The software, benchmark harness and the commits of every compared tool are available at "
     "https://github.com/kokyriakidis/dinara-align. The benchmarks are reproduced by `pixi run bench-astarpa2`, "
     "`results-astarpa2`, `bench`, `bench-local`, `bench-modes`, `bench-batch` and `bench-ablation`."),
]

TABLES = {
"engines": {
  "caption": "Algorithms in dinara-align and the changes introduced in this work.",
  "header": ["Problem", "Algorithm", "Changes (section)"],
  "align": "lll",
  "wrap": [1, 2],
  "rows": [
    ["Unit costs, long pairs", "Band doubling over bit-parallel tiles with the gap-chaining seed heuristic (A*PA2)", "Regrouped recurrence, staggered groups (2.2); agreeing bounds (2.3); neighbour construction (2.5)"],
    ["Unit costs, similar pairs", "Diagonal transition (WFA, BiWFA)", "Gathered extension (2.4)"],
    ["Gap-affine, two-piece", "Bidirectional wavefront (WFA, BiWFA)", "Gathered extension (2.4); score-derived band (2.7)"],
    ["Substitution matrices, local", "Anti-diagonal vector sweep; end-first local alignment (SSW)", "Backward wavefront from the best end (2.7)"],
    ["Batches of short pairs", "Inter-sequence vectorisation (SWIPE, SeqAn, parasail)", "Certified per-lane bands (2.6)"],
  ]},
"ablation": {
  "caption": "Effect of each technique: median time with all techniques on (ms per alignment; per batch of 500,000 pairs for the short reads), and the change when one technique is switched off, over three alternating rounds on the Skylake-X. Every configuration returned the same costs.",
  "header": ["Dataset", "All on", "Regrouped recurrence", "Paired groups", "Gathers", "Agreement rule", "Seed neighbours", "Certified band"],
  "align": "lrrrrrrr",
  "rows": [
    ["ont-500k", "97.8", "+4.2%", "+6.1%", "+2.5%", "−0.4%", "**+11.4%**", ""],
    ["ont-500k-genvar", "135.7", "+4.3%", "+6.5%", "+2.3%", "−0.3%", "**+19.0%**", ""],
    ["ont-1k", "0.021", "+1.0%", "+0.2%", "**+17.2%**", "−0.1%", "+0.4%", ""],
    ["ont-10k", "0.152", "0.0%", "+1.4%", "+7.0%", "**+17.3%**", "+1.0%", ""],
    ["SARS-CoV-2", "0.198", "−0.7%", "−0.5%", "**+15.9%**", "+4.5%", "−0.9%", ""],
    ["100 kbp, 1%", "0.978", "−0.1%", "+0.1%", "**+24.5%**", "−0.9%", "+0.7%", ""],
    ["100 kbp, 5%", "4.52", "+0.3%", "+0.4%", "+2.2%", "0.0%", "−0.1%", ""],
    ["100 kbp, 10%", "10.3", "+1.9%", "+2.6%", "+3.9%", "+0.1%", "**+13.0%**", ""],
    ["100 kbp, 15%", "10.85", "+3.6%", "+4.5%", "+6.8%", "0.0%", "**+13.3%**", ""],
    ["Short reads, affine", "232", "", "", "", "", "", "**+187%**"],
    ["Short reads, unit", "212", "", "", "", "", "", "**+214%**"],
  ]},
"real": {
  "caption": "Unit-cost alignment of the A*PA2 real datasets: mean time per alignment, traceback included. Superscripts give the number of pairs completed within 20 s where not all were; a dash marks a tool that completed none.",
  "header": ["Dataset", "dinara-align", "A*PA2-full", "A*PA2-simple", "A*PA", "Edlib", "BiWFA", "WFA"],
  "align": "lrrrrrrr",
  "rows": [
    ["ont-1k (0.8 kbp)", "**23 µs**", "81 µs", "48 µs", "612 µs", "126 µs", "44 µs", "31 µs"],
    ["ont-10k (3.6 kbp)", "**163 µs**", "391 µs", "247 µs", "15.3 ms", "1.1 ms", "753 µs", "595 µs"],
    ["ont-50k (9.5 kbp)", "**587 µs**", "1.27 ms", "944 µs", "190 ms^81/104", "6.1 ms", "7.43 ms", "6.37 ms"],
    ["ont-500k (638 kbp)", "**104 ms**", "176 ms", "537 ms", "–", "4.99 s^4/15", "19.8 s^1/15", "–"],
    ["genvar (651 kbp)", "**144 ms**", "210 ms", "624 ms", "–", "5.23 s^4/15", "6.36 s^3/15", "–"],
    ["SARS-CoV-2 (30 kbp)", "**279 µs**", "2.00 ms", "728 µs", "6.49 ms", "8.15 ms", "895 µs", "983 µs"],
  ]},
"affine": {
  "caption": "Gap-affine alignment at costs (4, 6, 2): mean time per alignment with a 5 s budget per sample. Superscripts give the number of pairs completed where not all were.",
  "header": ["Dataset", "dinara-align", "WFA", "BiWFA", "KSW2"],
  "align": "lrrrr",
  "rows": [
    ["ont-1k", "**112 µs**", "186 µs", "304 µs", "1.1 ms"],
    ["ont-10k", "**1.68 ms**", "4.8 ms", "5.84 ms", "30.8 ms^163/277"],
    ["ont-50k", "**21.2 ms**", "68.4 ms^82/104", "57.3 ms^90/104", "367 ms^16/104"],
    ["SARS-CoV-2", "**755 µs**", "2.06 ms", "1.79 ms", "88.4 ms"],
    ["10 kbp, 5%", "**2.17 ms**", "6.36 ms", "8.02 ms", "97.1 ms^52/99"],
    ["10 kbp, 15%", "**12.2 ms**", "38.5 ms", "44.8 ms", "203 ms^25/99"],
  ]},
"batch": {
  "caption": "Global scores of 500,000 short-read pairs; each time is the whole call.",
  "header": ["Costs", "Threads", "dinara-align", "WFA2-lib", "KSW2", "parasail", "Edlib"],
  "align": "llrrrrr",
  "rows": [
    ["Affine (1, 2, 1)", "1", "**233 ms**", "1.06 s", "10.9 s", "21.7 s", "–"],
    ["Unit", "1", "**212 ms**", "450 ms", "–", "–", "2.0 s"],
    ["Affine (1, 2, 1)", "10", "**26 ms**", "115 ms", "1.18 s", "2.34 s", "–"],
    ["Unit", "10", "**26 ms**", "48 ms", "–", "–", "223 ms"],
  ]},
}

FIGURES = {
"sweeps": "Unit-cost alignment of uniform pairs, mean time per alignment on a logarithmic scale. (a) 100 kbp pairs by divergence; dinara-align had the lowest mean time at every divergence except 10%. (b) Pairs at 15% divergence by length; lines end where a tool completed no pair within 20 s, and hollow points mark tools that completed only some pairs.",
}

REFS = {
 "Gotoh1982": ("Gotoh", "1982"), "Ukkonen1985": ("Ukkonen", "1985"), "Myers1999": ("Myers", "1999"),
 "Sosic2017": ("Šošić and Šikić", "2017"), "Myers1986": ("Myers", "1986"), "MarcoSola2021": ("Marco-Sola et al.", "2021"),
 "MarcoSola2023": ("Marco-Sola et al.", "2023"), "GrootKoerkamp2024a": ("Groot Koerkamp and Ivanov", "2024"),
 "GrootKoerkamp2024b": ("Groot Koerkamp", "2024"), "Wozniak1997": ("Wozniak", "1997"), "Farrar2007": ("Farrar", "2007"),
 "Zhao2013": ("Zhao et al.", "2013"), "Suzuki2018": ("Suzuki and Kasahara", "2018"), "Li2018": ("Li", "2018"),
 "Daily2016": ("Daily", "2016"), "Liu2023": ("Liu and Steinegger", "2023"), "Rognes2011": ("Rognes", "2011"),
 "Rahn2018": ("Rahn et al.", "2018"), "Kallenborn2026": ("Kallenborn et al.", "2026"), "Modular": ("Modular Inc.", "2026"),
 "Lattner2021": ("Lattner et al.", "2021"), "Eizenga2022": ("Eizenga and Paten", "2022"), "hyalite": ("Ferguson", "2026"),
 "Gao2021": ("Gao et al.", "2021"), "Gibrat2018": ("Gibrat", "2018"), "Fujiki2020": ("Fujiki et al.", "2020"),
}
# The bibliography keys these map to in references.bib.
BIBKEYS = {
 "Gotoh1982": "Gotoh1982", "Ukkonen1985": "Ukkonen1985", "Myers1999": "Myers1999bitvector", "Sosic2017": "Sosic2017edlib",
 "Myers1986": "Myers1986ond", "MarcoSola2021": "MarcoSola2021wfa", "MarcoSola2023": "MarcoSola2023biwfa",
 "GrootKoerkamp2024a": "GrootKoerkamp2024astarpa", "GrootKoerkamp2024b": "GrootKoerkamp2024astarpa2", "Wozniak1997": "Wozniak1997",
 "Farrar2007": "Farrar2007striped", "Zhao2013": "Zhao2013ssw", "Suzuki2018": "Suzuki2018diff", "Li2018": "Li2018minimap2",
 "Daily2016": "Daily2016parasail", "Liu2023": "Liu2023blockaligner", "Rognes2011": "Rognes2011swipe", "Rahn2018": "Rahn2018seqan",
 "Kallenborn2026": "Kallenborn2026accelign", "Modular": "Modular2025mojo", "Lattner2021": "Lattner2021mlir",
 "Eizenga2022": "Eizenga2022wfalm", "hyalite": "hyalite", "Gao2021": "Gao2021abpoa", "Gibrat2018": "Gibrat2018band", "Fujiki2020": "Fujiki2020seedex",
}
