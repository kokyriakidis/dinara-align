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
  'direct construction of the neighbours of exact seed matches that provably preserves local pruning, diagonal '
  'transition vectorised across diagonals, band bounds accepted only when independent estimates agree, and per-pair '
  'optimality certificates for banded inter-sequence batches, together with a bit-parallel sweep arranged for '
  'instruction latency. On an Intel Skylake-X processor, dinara-align had the '
  'lowest mean running time among the exact aligners tested on all 24 datasets of the A*PA2 and QuickEd benchmarks, '
  'real Illumina and PacBio HiFi reads among them, aligning ultra-long nanopore reads in 118 ms on average, '
  'compared with 283 ms for A*PA2-full, and HiFi reads 5.2 times faster than WFA, and was 8.0 times faster than '
  'QuickEd in the geometric mean. It was 1.1 to 3.8 times faster than WFA and BiWFA at gap-affine costs '
  'and computed global scores for 500,000 short-read pairs on one core 4.6 times faster than WFA2-lib.'),
 ('Availability and implementation',
  'dinara-align is free software under the Mozilla Public License 2.0, available at '
  'https://github.com/kokyriakidis/dinara-align with C, C++, Python and command-line interfaces.'),
 ('Contact', 'kkyriaki@ucsc.edu, bpaten@ucsc.edu'),
 ('Supplementary information', 'Supplementary data are available at *Bioinformatics* online.')]

BODY = [
('sec', 'Introduction'),
('p',
 'Pairwise alignment is a building block of read mapping, genome assembly and variant calling. It determines the'
 ' minimum-cost series of substitutions, insertions and deletions that transforms one sequence into another, and'
 ' the classical dynamic programming algorithm computes it over an \\(n \\times m\\) matrix [@Gotoh1982]. That is '
 'practical for short reads but not for current long reads: nanopore reads span tens to hundreds of kilobases, '
 'often with 5 to 15% errors, and one pair of 500 kbp reads fills a matrix of \\(2.5 \\times 10^{11}\\) cells.'),
('p',
 'Heuristic aligners compute only part of the matrix, a band, a window or the cells before a drop in score, and '
 'lose the guarantee of optimality. Exact aligners instead skip only those cells that provably cannot lie on an '
 'optimal path, so that their running time depends on the distance between the sequences rather than on the '
 'product of their lengths.'),
('p',
 'Three ideas underlie current exact aligners. *Band doubling* computes only the diagonals that a path within a '
 'given cost bound can reach, and doubles the bound until the optimal cost lies within it [@Ukkonen1985]. '
 '*Bit-parallelism* encodes a column of 64 cells in one machine word [@Myers1999], and Edlib combines it with '
 'band doubling [@Sosic2017]. *Diagonal transition* extends the furthest-reaching paths one cost at a time, '
 'following runs of matches at no cost [@Ukkonen1985; @Myers1986]; WFA generalised it to gap-affine costs '
 '[@MarcoSola2021] and BiWFA reduced its memory to linear in the distance [@MarcoSola2023]. A*PA introduced A* '
 'search with a gap-chaining seed heuristic [@GrootKoerkamp2024a], and A*PA2 combined this heuristic with band '
 'doubling over SIMD blocks of bit-parallel words, aligning ultra-long reads up to 19 times faster than Edlib '
 'and BiWFA [@GrootKoerkamp2024b]. QuickEd first bounds the distance with windowed and banded heuristics and '
 'then aligns within the bound [@Doblas2025]. A second line of work accelerates the full dynamic program on '
 'vector hardware, using anti-diagonal [@Wozniak1997], striped [@Farrar2007; @Zhao2013] and prefix-scan '
 '[@Daily2016] layouts, difference recurrences [@Suzuki2018; @Li2018], adaptive blocks [@Liu2023] and '
 'inter-sequence parallelism [@Rognes2011; @Rahn2018]; these methods compute the full matrix or a fixed band, so'
 ' their running time grows with its size rather than with the distance.'),
('p',
 'In the exact methods, a large part of the remaining running time is determined by how the computation maps '
 'onto the processor. The seed heuristic is built by scalar lookups in tables larger than the caches. Diagonal '
 'transition extends one diagonal at a time. Band doubling spends whole rounds when its bound is far from the '
 "distance. Batches of short pairs are vectorised across pairs, but each pair's full matrix is computed. And the"
 ' bit-parallel sweep carries a dependency from each column to the next, so its speed is set by the latency of '
 'that chain rather than by its number of operations. The choice among these algorithms rests on thresholds '
 'measured on earlier versions of the code, and the boundaries that a band keeps for its traceback take memory '
 'that grows with the distance. Here we present dinara-align, an exact aligner whose algorithms address these '
 'bottlenecks. Our contributions are:'),
('list',
 ['**Direct neighbour construction for the seed heuristic.** The inexact neighbours of an exact seed match are '
  'inserted without search, and we prove that this preserves the outcome of local pruning (Lemma 1).',
 '**Gathered diagonal transition.** Diagonal transition and the gap-affine wavefront extend eight diagonals per '
  'step using gather instructions.',
 '**Band bounds from agreeing estimates.** A band-doubling bound is extrapolated from partial progress only when'
  ' two independent estimates agree.',
 '**Certified bands for inter-sequence batches.** Each short pair is aligned in one vector lane within a narrow '
  'band whose optimality is certified per lane by a lower bound on the cost of any path that leaves it '
  '(Proposition 1), adapting the speculate-and-verify approach of the SeedEx accelerator [@Fujiki2020] to '
  'software batches.',
 '**A latency-aware bit-parallel sweep.** On AVX-512, two staggered groups of words are interleaved, extending '
  'the instruction-level parallelism of A*PA2, and two Boolean identities shorten the loop-carried chain of '
  "Myers' step; the two are worth 4.5 to 5% and 1 to 2% on the ultra-long reads, and 11.5% and 4% on 100 kbp "
  'pairs at 20% errors.',
 '**A traceback within a memory limit that keeps the tie rule.** When the tile boundaries a band keeps for its '
  'traceback would exceed the memory limit, only some are kept, and the stretches between them are computed again'
  ' from the recorded state of the band, so the alignment is the same byte for byte.',
 '**Thresholds measured against the code.** Every threshold that selects an algorithm or bounds a search was '
  'measured again on the final code, which exposed thresholds that had cost up to a third of the running time.']),
('p',
 'Section 2 reviews the algorithms dinara-align builds on, Section 3 describes our changes, Section 4 compares '
 'it with the exact aligners of the A*PA2 and QuickEd benchmarks, and Section 5 discusses the results and their '
 'limitations.'),
('sec', 'Background'),
('sub', 'Alignment and costs'),
('p',
 'Let \\(A=a_1\\dots a_n\\) be the reference and \\(B=b_1\\dots b_m\\) the query, and let cell \\((i,j)\\) lie on '
 'diagonal \\(j-i\\). Under gap-affine costs, a mismatch costs \\(x\\) and a gap of \\(k\\) letters costs '
 '\\(g(k)=o+ke\\); unit costs correspond to \\(x=e=1\\) and \\(o=0\\), and two-piece costs take the minimum of two '
 'such functions [@Li2018]. A substitution matrix instead assigns each pair of letters a score to be maximised. '
 'We consider global, ends-free, extension and local alignment. An alignment is *exact* if it is optimal under '
 'the requested costs, mode, band and cost limit; among optimal alignments, a fixed rule places indels leftmost '
 'or rightmost, so that, under costs, the result does not depend on which algorithm computed it (Section 3.1). '
 'At unit costs, the distance \\(D(i,j)\\) between the prefixes \\(a_1\\dots a_i\\) and \\(b_1\\dots b_j\\) satisfies'),
('p',
 '\\[ D(i,j) = \\min\\bigl(D(i-1,j)+1,\\; D(i,j-1)+1, \\qquad D(i-1,j-1)+[a_i \\ne b_j]\\bigr), \\]'),
('p',
 'with \\(D(i,0)=i\\) and \\(D(0,j)=j\\), and the distance is \\(D(n,m)\\); gap-affine costs need three such layers '
 '[@Gotoh1982].'),
('sub', 'Band doubling and bit-parallel computation'),
('p',
 'An alignment of cost at most \\(t\\) visits only diagonals \\(d\\) with \\(|d| + |d - \\delta| \\le t\\), where '
 '\\(\\delta = m - n\\) is the diagonal on which it ends. Band doubling computes only those cells, for \\(t = t_0, '
 '2t_0, \\dots\\), until \\(D(n,m) \\le t\\), in \\(O(ns)\\) time for distance \\(s\\) [@Ukkonen1985]. Myers [@Myers1999]'
 ' encodes the vertical differences of 64 consecutive cells of a column, each \\(-1\\), 0 or \\(+1\\), in two '
 'machine words \\(P_v\\) and \\(M_v\\), and advances them by one column with about twenty bitwise operations. Edlib'
 ' combines this encoding with band doubling [@Sosic2017]. A*PA2 computes the band in tiles of 256 columns, '
 'several words at once in vector registers, and drops the rows of a tile whose score already rules them out '
 '[@GrootKoerkamp2024b].'),
('sub', 'Diagonal transition'),
('p',
 'Diagonal transition stores, for each cost \\(s\\) and diagonal \\(d\\), the furthest cell on \\(d\\) that an '
 'alignment of cost \\(s\\) reaches, and extends it along matches at no cost; the distance is the first \\(s\\) '
 'whose wavefront reaches \\((n, m)\\), after \\(O(n + s^2)\\) steps [@Ukkonen1985; @Myers1986]. WFA generalises it '
 'to gap-affine costs with three wavefronts per cost [@MarcoSola2021], and BiWFA runs it from both ends until '
 'the wavefronts meet, which keeps its memory linear in \\(s\\) [@MarcoSola2023].'),
('sub', 'The seed heuristic'),
('p',
 'A*PA divides the reference into disjoint seeds of \\(k\\) bases and finds their matches in the query, exactly or'
 ' with one edit [@GrootKoerkamp2024a]. Its gap-chaining seed heuristic bounds the cost from a cell to the end '
 'of the alignment by the seeds ahead that no chain of matches from that cell covers, counting the gaps between '
 'chained matches, and never overestimates. Local pruning discards a match \\(\\nu\\) unless \\(c(\\nu) + h(e(\\nu)) <'
 ' \\tau\\), where \\(c(\\nu)\\) is its cost, \\(e(\\nu)\\) its end cell, \\(h(u)\\) the least cost of a path from \\(u\\) '
 'across the following seeds, and \\(\\tau\\) a threshold that depends only on the seeds. A*PA2 uses the heuristic '
 'to drop the rows of each band whose score plus the heuristic exceeds the bound [@GrootKoerkamp2024b].'),
('sec', 'Methods'),
('sub', 'Overview'),
('p',
 'Each problem is assigned to an algorithm ([#tab:engines]; [#fig:overview]). For unit costs we follow A*PA2 '
 '[@GrootKoerkamp2024b]: band doubling over tiles of 256 columns of bit-parallel words, pruned by the '
 'gap-chaining seed heuristic [@GrootKoerkamp2024a], with a diagonal-transition traceback inside each tile; '
 'near-identical pairs are aligned by diagonal transition alone. Gap-affine costs use the bidirectional '
 'wavefront of WFA and BiWFA [@MarcoSola2021; @MarcoSola2023]. A short diagonal-transition search estimates the '
 'distance, and the estimate, the length and the cost model select the algorithm. The thresholds are fixed at '
 'compile time for each instruction set (Section 3.9).'),
('p',
 'Under costs, every engine selects the same alignment among equally good ones. All follow the backtrace order '
 'of WFA2-lib: from the end of the alignment, an edit is taken into a cell whenever one is optimal, a '
 'substitution before a deletion before an insertion, and a gap is extended before it is opened. The tile '
 'traceback applies this order cell by cell (Section 3.7), the batches and the difference sweep store per-cell '
 'flags under it, and the wavefronts take the furthest-reaching source, which makes the same choice. Leftmost '
 'placement of indels is obtained by running this rightmost rule on both sequences reversed and reversing the '
 'result. Two exceptions remain: a gap-affine pair whose wavefronts would exceed the memory limit is split where'
 ' an optimal path crosses and follows the rule within each piece, and alignments under a substitution table '
 "follow the traceback order of Gotoh's recurrence."),
('figure', 'overview'),
('table', 'engines'),
('sub', "A shorter critical path for Myers' recurrence"),
('p',
 "As in Edlib and A*PA2 [@Sosic2017; @GrootKoerkamp2024b], the query is packed 64 bases per word, and Myers' "
 'recurrence [@Myers1999] advances a word by one column. With match mask \\(E\\), vertical difference masks \\(P_v,'
 " M_v\\) and the incoming negative horizontal difference folded into \\(E' = E \\lor M_h^{\\mathrm{in}}\\), A*PA2 "
 'computes'),
('p',
 "\\[ c = (E' \\land P_v) + P_v, \\qquad X_h = (c \\oplus P_v) \\lor E', \\qquad P_h = M_v \\lor \\lnot(X_h \\lor P_v), "
 '\\qquad M_h = P_v \\land X_h. \\]'),
('p',
 'The horizontal masks \\(P_h\\) and \\(M_h\\) determine the vertical masks of the next column, so the path from '
 '\\(c\\) to them lies on the loop-carried chain. Since \\((c \\oplus P_v) \\lor P_v = c \\lor P_v\\) and \\(P_v \\land '
 '(c \\oplus P_v) = P_v \\land \\lnot c\\), the same masks are obtained as'),
('p',
 "\\[ P_h = M_v \\lor \\lnot\\bigl(c \\lor (E' \\lor P_v)\\bigr), \\qquad M_h = (P_v \\land \\lnot c) \\lor (P_v \\land E'),"
 ' \\]'),
('p',
 "where \\(E' \\lor P_v\\) and \\(P_v \\land E'\\) do not depend on \\(c\\) and are computed off the critical path. This"
 ' removes two operations from the loop-carried chain and, for a single word, reduces the cost from 14 to 8 '
 'cycles per column on the Skylake-X. Myers [@Myers1999], Edlib and A*PA2 use the original form, and we are not '
 'aware of a previous description of the regrouped one.'),
('p',
 'Words occupy the lanes of a vector, staggered by one column, so that the horizontal difference leaving one '
 'word enters the next through a lane rotation [@Wozniak1997; @GrootKoerkamp2024b]. Every column then waits for '
 'the rotation of the previous one, and A*PA2 processes two four-lane vectors together to overlap these waits. '
 'On AVX-512, where a vector holds eight words, dinara-align interleaves two eight-word groups, the lower '
 'lagging ten columns behind the upper, so that the instructions of each group fill the latency of the other; '
 'with narrower vectors a single group is used. Bands of two or three words, in which every operation lies on '
 'the critical path, are computed in general-purpose registers, whose operations have lower latency. On AVX-512 '
 'the sweep requires about 2.4 cycles per word and column.'),
('sub', 'Band bounds from agreeing estimates'),
('p',
 'Band doubling returns the exact distance for any initial bound, but its running time depends on how close that'
 ' bound is to the distance. A partial computation that has reached cost \\(s\\) at anti-diagonal \\(r\\) suggests '
 'the extrapolated distance \\(\\hat d = s(n+m)/r\\). This extrapolation is reliable when edits are spread evenly '
 'along the pair and can be too high by orders of magnitude when they are concentrated, as near the ends of '
 'nanopore reads. dinara-align therefore computes two extrapolations, \\(\\hat d_1\\) from the first eight edits of'
 ' a diagonal-transition search and \\(\\hat d_2\\) from the same search continued until its budget, and sets the '
 'first bound from them only if \\(\\max(\\hat d_1,\\hat d_2) \\le \\tfrac43 \\min(\\hat d_1,\\hat d_2)\\). Otherwise, as '
 'in A*PA2, the bound starts 256 above the known lower bound, and each retry at most doubles it. A failed round '
 'sets the next bound from where it stopped, its progress extrapolated to the whole width, rather than doubling.'
 ' Within a round, the distance is extrapolated again after one eighth, one quarter and one half of the columns;'
 ' the bound is lowered if the extrapolation lies below it, and the round is abandoned if the extrapolation '
 'exceeds it. Because these rules change only the sequence of bounds, the result remains exact.'),
('p',
 'Within a round, the band reads the scores down the right edge of each tile and drops the rows at its top and '
 'bottom whose score plus the heuristic exceeds the bound, as A*PA2 does. A score changes by at most one per '
 'row, and the heuristic by at most its climb \\(c\\), one for exact seeds and two for inexact ones, so a row that'
 ' exceeds the bound by \\(x\\) rules out the next \\(\\lceil x/(1+c) \\rceil\\) rows without reading them. The bottom'
 ' of the next tile is bounded by the deepest row kept, since going \\(k\\) rows past it costs at least \\(k\\) '
 'more. With seeds, when at least one seed in 32 chains from the origin, the first bound starts 64 above the '
 'heuristic at the origin, which on long reads lies within a few percent of the distance.'),
('sub', 'Diagonal transition for near-identical pairs'),
('p',
 'Diagonal transition extends each diagonal of a wavefront along its run of matches. dinara-align compares eight'
 ' letters at a time by the exclusive-or of two 64-bit words and a count of trailing zeros. On AVX-512, eight '
 'diagonals are extended together: their letters are loaded with gather instructions, and the wavefront is '
 'padded with unreachable diagonals so that every step processes full vectors under branch-free masks. The '
 'gather masks are derived from the lanes of each group, which avoids a dependency between consecutive gathers '
 'that would otherwise serialise them. The same vectorised step serves near-identical pairs, computed from one '
 'or both ends [@MarcoSola2023], the traceback inside each tile, restricted to the diagonals that can still '
 'reach the traced cell, and the three-layer gap-affine wavefront.'),
('p',
 'A pair whose diagonal transition finishes within its budget is aligned from its wavefronts alone. From one end'
 ' the wavefronts are kept whole. From both ends, as in BiWFA, the forward wavefronts are grown on to the proven'
 ' distance, each cost kept only on the diagonals on which the backward wavefronts show that an optimal path '
 'passes, so that the extra costs grow few diagonals. The backtrace then takes, at each cost, the '
 'furthest-reaching source, which is the order of Section 3.1, and the CIGAR is written run by run from the '
 'wavefronts.'),
('sub', 'Seed heuristic: direct neighbour construction'),
('p',
 'The reference is divided into disjoint seeds of 12 bases, matched exactly, or of 16 bases, matched with at '
 'most one edit [@GrootKoerkamp2024a]. A 16-base seed that matches with at most one edit matches one of its two '
 '8-base halves exactly, so each half indexes a table. Since these tables exceed the caches, the lookups for 16 '
 'query rows are issued together, so that their memory accesses overlap, and only seeds whose matches could '
 'still chain to the end of the pair are considered.'),
('p',
 'An exact match \\(\\mu\\) of a seed at query rows \\([i, i+16)\\) always has four neighbouring matches with one '
 'edit: two share its start, at rows \\([i, i+15)\\) and \\([i, i+17)\\), and two share its end, at \\([i+1, i+16)\\) '
 'and \\([i-1, i+16)\\). Because chaining is all-or-nothing, a neighbour one diagonal away may chain where \\(\\mu\\)'
 ' cannot, and the neighbours are required for the heuristic to remain a lower bound. Rather than locating them '
 'by search and testing each for local pruning, dinara-align inserts them whenever \\(\\mu\\) is retained, which is'
 ' justified by the following lemma. Local pruning retains a match \\(\\nu\\) only if \\(c(\\nu) + h(e(\\nu)) < \\tau\\)'
 ' (Section 2.4).'),
('p',
 '**Lemma 1.** If an exact match \\(\\mu\\) fails local pruning, so does each of its four neighbours.'),
('p',
 '*Proof.* A neighbour \\(\\nu\\) that shares the end of \\(\\mu\\) has \\(e(\\nu)=e(\\mu)\\) and \\(c(\\nu)=1>c(\\mu)=0\\). A'
 ' neighbour that shares the start ends on an adjacent diagonal, and one insertion or deletion from \\(e(\\mu)\\) '
 'reaches a cell on that diagonal at or beyond \\(e(\\nu)\\). The least remaining cost does not increase along a '
 'diagonal [@Ukkonen1985], so \\(h(e(\\mu)) \\le 1 + h(e(\\nu))\\) and \\(c(\\mu)+h(e(\\mu)) \\le c(\\nu)+h(e(\\nu))\\). In '
 'both cases the pruning test value of \\(\\nu\\) is at least that of \\(\\mu\\). \\(\\square\\)'),
('p',
 'Neighbours sharing the start therefore raise the layer of that start by the chain length available from their '
 'ends, and neighbours sharing the end are placed one layer below it, without a search or a pruning test of '
 'their own.'),
('sub', 'Certified bands for inter-sequence batches'),
('p',
 'Batches of short pairs are aligned with one pair per vector lane [@Rognes2011; @Rahn2018], first in saturating'
 ' 8-bit lanes and, for the scores of pairs that saturate, in 16-bit lanes [@Zhao2013; @Daily2016]; alignments '
 'stay in 8-bit lanes. Rather than computing each matrix in full, a lane computes the band of diagonals '
 '\\([d_\\ell, d_h]\\) containing the main diagonal 0 and the end diagonal \\(\\delta = m - n\\), together with one '
 'further diagonal on either side. Its optimality is certified with the following bound, a cost form of the '
 'out-of-band bounds of [@!Gibrat2018] and [@!Fujiki2020], stated for gap-affine and two-piece costs and with a '
 'strict form for tie-breaking.'),
('p',
 '**Proposition 1.** Let \\(C\\) be the least cost of an alignment within \\([d_\\ell, d_h]\\), and let \\(L(d) = '
 'g(|d|) + g(|d-\\delta|)\\). If \\(C \\le \\min\\bigl(L(d_h+1), L(d_\\ell-1)\\bigr)\\), then \\(C\\) is the optimal cost; '
 'if the inequality is strict, every optimal alignment lies within the band.'),
('p',
 '*Proof.* An alignment that leaves the band visits diagonal \\(d_h+1\\) or \\(d_\\ell-1\\). Consider \\(d=d_h+1 > '
 '\\max(0,\\delta)\\); the case below the band is symmetric. Starting on diagonal 0, the alignment contains '
 'insertions totalling at least \\(d\\) letters to reach \\(d\\), and deletions totalling at least \\(d-\\delta\\) '
 'letters to end on \\(\\delta\\). These form at least one gap of each kind, and since \\(g\\) is subadditive (for '
 'two-piece costs, with the usual convention that the longer piece has the lower extension cost), their total '
 'cost is at least \\(g(d)+g(d-\\delta)=L(d)\\). Every alignment outside the band therefore costs at least '
 '\\(\\min\\bigl(L(d_h+1), L(d_\\ell-1)\\bigr)\\). \\(\\square\\)'),
('p',
 'Bounds of this kind underlie band doubling [@Ukkonen1985] and have been used to verify narrow-band results in '
 'software [@Gibrat2018] and in hardware [@Fujiki2020]; here the certificate is evaluated independently in each '
 'vector lane of an inter-sequence batch. For alignments, the strict form guarantees that the alignment selected'
 ' by the tie rule lies within the band. Pairs that fail the certificate, about one in seven in our short-read '
 'benchmark, are recomputed over every diagonal that a cheaper path could visit, grouped by the width of that '
 'band so that the lanes of a group finish together.'),
('sub', 'Traceback within a memory limit'),
('p',
 'The traceback of a band runs from its last tile back, and traces each tile from the boundary recorded on its '
 'left: the vertical differences down that boundary, within the band. On a 1 Mbp pair at 20% errors these '
 'boundaries take up to 175 MB, so the memory limit of an alignment, 80 MB by default, applies to them. While '
 'they fit, every boundary is kept; beyond the limit only every \\(k\\)-th is, \\(k\\) doubling as the round '
 'proceeds, with room left for \\(k\\) of the widest boundaries. Each tile also records the state that the band '
 'carries into it: its first and last words, the score at its top, the deepest row kept and its score, and the '
 'bound. The traceback computes a stretch of tiles again from the last kept boundary before a missing one. Every'
 ' step of a band depends only on this state and on the scores it computes, so each boundary is recovered as the'
 ' round left it, and the alignment is the same byte for byte. Splitting the pair instead, as QuickEd does '
 '[@Doblas2025] and as our gap-affine wavefront does past its own limit, aligns the parts independently and '
 'resolves ties within each. On the 1 Mbp pairs at 20% errors, the growth of peak memory fell from 175 to 133 MB'
 ' on the worst pair, and the four pairs took 9% longer; pairs whose boundaries fit are unchanged.'),
('sub', 'Gap-affine and score-based alignment'),
('p',
 'A scoring scheme with match score \\(a\\), mismatch score \\(b\\) and gap scores \\(-(o+ke)\\) can be converted into'
 ' costs for global alignment [@Eizenga2022]. Every alignment covers \\(n+m\\) letters, so its score \\(S\\) and a '
 'gap-affine cost \\(C\\) satisfy'),
('p',
 '\\[ S = \\tfrac12\\bigl(a(n+m) - C\\bigr), \\qquad C = 2(a-b)\\,X + \\sum_{\\text{gaps}} \\bigl(2o + k(2e+a)\\bigr), \\]'),
('p',
 'where \\(X\\) is the number of mismatches, so the wavefront computes score-maximising alignments. The cost also '
 'bounds the band in which the alignment must be stored. Each gapped letter contributes at least \\(2e+a\\) to '
 '\\(C\\), and an alignment that deviates by \\(t\\) diagonals from both the main and the end diagonal contains at '
 'least \\(2t\\) gapped letters; hence only diagonals within \\(C/(2(2e+a))\\) of them can contain an optimal '
 'alignment. For a 1 kbp pair at 5% divergence, this stores 137 diagonals instead of 1,000. Local alignment '
 'first locates the end of the optimal alignment with a forward sweep, as in SSW [@Zhao2013], and then extends a'
 ' wavefront backwards from that end, stopping at the first cell that attains the maximal score; this avoids the'
 ' correction pass that striped layouts may need when long alignments score highly [@Farrar2007].'),
('p',
 'The bidirectional wavefront keeps, for every cost and diagonal, the furthest column and one byte naming the '
 'source of each layer, five bytes per diagonal, and traces the alignment through these fronts from where the '
 'two searches met; a pair whose fronts would exceed the memory limit is split where an optimal path crosses, as'
 ' in BiWFA. A pair that ends within the costs one search can hold is finished by that search alone, before '
 'the other starts: for a distance the forward search, which saved a third of the time of short reads at 1% '
 'errors; for an alignment the search the traceback reads, traced straight back from the end. Two searches '
 'must each pass half the cost and a margin, the widest opening and the dearest move, before they may stop, '
 'so a pair of cost d costs them about 2W((d + margin)/2), W(c) being the diagonals a search steps to cost c, '
 'and one search W(d); that search therefore grows alone while the W it has measured stays within twice its '
 'value at (c + margin)/2, which hands over at about 2.4 margins where fronts widen with the cost and never '
 'where the sequences bound them, and the handover loses no meeting. On 250 bp Illumina reads at gap-affine '
 'costs, starting with one search took an alignment from 2.41 to 1.52 µs, and growing it by measured work '
 'saved a further 5%; and a cost '
 'without traceback stops as soon as its best lies within one of the costs both searches have checked against '
 'the start of the other. A global pair whose projected work exceeds that of a full sweep in difference '
 'recurrences [@Suzuki2018], in 8-bit lanes, is handed to that sweep once the wavefront has spent one '
 'twenty-fourth of it; the sweep stores the predecessor flags of Section 3.1, so the alignment is the same, and '
 'unrelated pairs of 5 kbp aligned in 4.8 ms instead of 32.9 ms.'),
('p',
 'Under a general substitution matrix, anti-diagonals are computed in 16-bit lanes while the scores provably fit'
 ' [@Wozniak1997], a strip of rows of at most 16 KB at a time, so that the rows a strip reads stay in the '
 'first-level cache; this took local alignment of 10 kbp against 12 kbp from 50.0 to 17.4 ms. A global alignment'
 ' is traced in a band of diagonals that grows, at most fourfold per step, until the score it yields proves that'
 ' every optimal path lies inside it, so a pair close to its diagonal takes one sweep of about its own cells '
 'rather than a sweep of the matrix for its score: 1 kbp pairs at 10% divergence aligned in about 200 µs instead'
 ' of 544 µs. A band whose decisions would exceed the memory limit keeps every \\(k\\)-th anti-diagonal and '
 'computes each stretch again for the traceback.'),
('sub', 'Thresholds measured against the code'),
('p',
 'The choice of algorithm and the bounds of each search rest on 41 thresholds, from the step budgets of the '
 'diagonal transition to the seed gates and the limits of the traceback, each set from measurements when it was '
 'introduced, on the code of the time. We measured all of them again on the final code: dinara-align was built '
 'with each threshold halved and with each doubled, and the 83 builds were timed on 17 unit-cost and 8 '
 'gap-affine samples against the unchanged build, which was timed again after every six. About thirty thresholds'
 ' moved no sample by more than 3% in either direction, and most of the others traded one dataset against '
 'another. One was stale: the length from which seeds match with one edit, which had cost 100 kbp pairs at 10% '
 'errors a third of their running time. Two further stale choices were found by examining the slowest results '
 'directly: a narrower tile for narrow bands, kept from the first version, which had cost 10 kbp pairs at 5% '
 'errors a fifth of their time, and a projection that handed nanopore reads to a full sweep at gap-affine costs,'
 ' which had slowed them by half (Supplementary Section S1).'),
('sub', 'Validation'),
('p',
 'A differential fuzzer compares every mode, cost model, band, cost limit and tie rule on random pairs with a '
 'full-matrix implementation that shares no code with the library, and the WFA2-lib regression set is reproduced'
 ' exactly, every CIGAR included (Supplementary Section S1).'),
('sec', 'Results'),
('sub', 'Experimental setup'),
('p',
 'All measurements were taken on an Intel Core i9-7900X (Skylake-X, AVX-512), its ten cores fixed at 3.3 GHz '
 'with turbo boost and hyper-threading disabled, with dinara-align compiled by Mojo 1.1.0 [@Modular] at commit '
 '`036e12b` (gap-affine costs at `a6d9166`, the other alignment modes and batches at `289bcf0`) and every tool '
 'compiled for the native instruction set. We combined the protocols of A*PA2 [@GrootKoerkamp2024b] and QuickEd '
 '[@Doblas2025]: the real datasets of A*PA2, nanopore reads of about 1, 10, 50 and 500 kbp, the last also with '
 'genetic variation, and SARS-CoV-2 genomes; real Illumina 2×250 bp and PacBio HiFi reads, the two technologies '
 'QuickEd adds, here reads of HG002 from Genome in a Bottle [@Zook2016] covering one region of chromosome 20, '
 "each mapped by minimap2 [@Li2018] to both haplotypes of HG002's own assembly [@Hansen2025] and paired with the "
 'stretch it covers, so that a pair differs by sequencing errors alone; and the simulated datasets of QuickEd, '
 'random sequences of 10 kbp to 1 Mbp with 1, 5, 10 and 20% errors, each sampled to about 2 Mbp and at least '
 'four pairs (Supplementary Table S5). Each tool aligned each pair once, with traceback, on one pinned core, '
 "within 30 s per sample, through one aligner object per sample, as every tool's benchmark runner keeps one. The "
 'exact aligners compared are A*PA2-full and A*PA2-simple, A*PA [@GrootKoerkamp2024a], Edlib [@Sosic2017], WFA '
 'and BiWFA [@MarcoSola2021; @MarcoSola2023], QuickEd, BSAlign [@Shao2024], KSW2 [@Li2018], parasail '
 '[@Daily2016] and SSW [@Zhao2013], each on the workloads it supports; every cost agreed with those of every '
 'other tool and with the published A*PA2 results.'),
('sub', 'Speed and memory'),
('p',
 'dinara-align had the lowest mean time per alignment of the exact aligners on every workload ([#tab:unit]; '
 '[#tab:affine]). On the nanopore reads it was 1.4 to 2.4 times faster than the fastest alternative, on the '
 'Illumina and HiFi reads 2.0 and 5.2 times faster than WFA, the fastest of the others there, and on the '
 'simulated pairs 1.1 to 2.6 times ([#fig:sweeps]): 100 kbp pairs at 1, 5, 10 and 20% errors took 1.5, 5.1, 7.3 '
 'and 15.7 ms, against 7.5, 8.4, 13.9 and 36.7 ms for A*PA2-full and 12.5, 24.0, 39.0 and 62.8 ms for QuickEd, '
 'and its narrowest margin was on 10 kbp pairs at 5%, 0.267 against 0.288 ms for A*PA2-simple. Over all 24 '
 'datasets QuickEd took 3.0 to 75 times as long, 8.0 times in the geometric mean. At gap-affine costs '
 'dinara-align was 1.1 to 3.8 times faster than WFA and BiWFA wherever a tool completed a dataset, its '
 'narrowest margin on the Illumina reads, 1.45 against 1.55 µs for WFA, and it '
'completed more of the long pairs than any other tool; in the other modes it was 1.0 to 22 times faster than '
 'the fastest tool offering the same mode, the one tie lying within the variation between builds. Global scores '
 'of 500,000 short-read pairs took 230 ms on one core and 27 ms on ten at (1, 2, 1), 4.6 and 4.3 times faster '
 'than WFA2-lib.'),
('table', 'unit'),
('table', 'affine'),
('figure', 'sweeps'),
('p',
 'An alignment of an ultra-long read increased peak memory by at most 48 MB, against 82 MB for A*PA2-full, and '
 'of a 1 Mbp pair at 20% errors by at most 123 MB, against 184 MB, its traceback held within the default limit '
 "of 80 MB (Section 3.7; Supplementary Table S6); QuickEd stayed near 25 MB. On short reads dinara-align's peak "
 "is higher than the others' ([#tab:unit]), its runtime alone taking about 13 MB."),
('sub', 'Effect of each technique'),
('p',
 'To measure what each technique contributes, dinara-align was built with one switched off at a time and every '
 'build aligned samples of the A*PA2 real datasets and of the 100 kbp simulated pairs, and scored the short-read'
 ' batch, the builds alternating within three rounds on one pinned core, assembled with branches kept clear of '
 '32-byte boundaries, which the Skylake-X penalises ([#tab:ablation]); every build returned the same costs. '
 'Constructing seed neighbours directly contributed most on the ultra-long reads, 9% and 16%; 100 kbp pairs, '
 'which take exact seeds on AVX-512, do not use it. Interleaving two groups of bit-parallel words contributed '
 '4.5 to 5% on the ultra-long reads and 11.5% on 100 kbp pairs at 20% errors, whose bands are widest; the '
 "regrouped recurrence contributed 1 to 4%, as bands narrow enough to be bound by a single word's step are a "
 'small part of the work. Gathered diagonal transition contributed 13 to 27% where diagonal transition does most'
 ' of the work, and the agreement rule 22% on the 10 kbp reads. Without the certified band, the short-read batch'
 ' took 2.9 and 3.2 times as long.'),
('table', 'ablation'),
('sec', 'Discussion'),
('p',
 'The algorithms in dinara-align build on those of earlier exact aligners; the gains come from changes to how '
 'these algorithms use the processor. Two observations stand out. First, on long pairs the techniques that '
 'reduce work outside the inner loop contributed more than those that speed it up: constructing seed neighbours '
 'directly saved up to 16%, against 4.5 to 11.5% for interleaving bit-parallel groups and 1 to 4% for the '
 'regrouped recurrence ([#tab:ablation]). Second, thresholds that encode a cost model go stale as the code '
 'around them changes, and measuring them all again belongs in the evaluation of any aligner that chooses its '
 'algorithm by fixed thresholds (Section 3.9). Neither the regrouped recurrence nor the per-lane band '
 'certificate depends on the implementation, and both apply directly to other bit-parallel and inter-sequence '
 'aligners.'),
('p',
 'The study has several limitations:'),
('list',
 ['The evaluation covers a single x86 processor with AVX-512; the AVX2 configuration was measured only for its '
  'inexact-seed threshold, and the NEON configuration is not evaluated here.',
 'Some thresholds were chosen using workloads that also appear in the evaluation; in particular, the AVX-512 '
  'thresholds for inexact seeds, the share of chained seeds (lowered from 40 to 20%) and the length from which '
  'they are used (raised from 64 to 128 kbp), were set after measurements on uniform 100 kbp pairs like those of '
  '[#fig:sweeps]a, so results on those pairs may overstate performance on unseen data.',
 'As in the A*PA2 protocol, each pair was aligned once, and on the Skylake-X some running times varied by up to '
  'one third between builds that did not change the code involved, which we attribute to code alignment effects '
  'in the instruction decoder; the ablation controls for this, the comparison with other tools does not.',
 'The Illumina and HiFi datasets come from one region of chromosome 20, 20 and 200 kbp long, and keep only reads '
  'minimap2 maps to HG002 with at most a tenth clipped, so they leave out reads from repeats and structural '
  'variants.',
 'The ablation switches off one technique at a time against the full configuration, so it does not separate '
  'interactions between techniques.',
 'Release builds target the baseline processor of each platform, so the AVX-512 code paths require a build for '
  'the host processor.']),
('p',
 'Future work includes:'),
('list',
 ['faster construction of the seed heuristic, whose scalar lookups dominate on long divergent pairs;',
 'reusing computation across band-doubling rounds, as in the incremental doubling of A*PA2;',
 'a seed heuristic for gap-affine costs, at which no tool completed a 1 Mbp pair beyond 1% errors;',
 'measurements on further processors, and an evaluation of the GPU kernels the library provides for batch '
  'scores and alignments (Supplementary Section S1).']),
]

BACK = [
    ("Funding", "This work received no specific funding."),
    ("Conflict of interest", "None declared."),
    ("Data availability", "The software, benchmark harness and the commits of every compared tool are available at "
     "https://github.com/kokyriakidis/dinara-align. The benchmarks are reproduced by `pixi run bench` (unit and "
     "`--affine` costs), `bench-modes`, `bench-batch` and `bench-ablation`."),
]

TABLES = {
"engines": {
  "caption": "Algorithms in dinara-align and the changes introduced in this work.",
  "header": ["Problem", "Algorithm", "Changes (section)"],
  "align": "lll",
  "wrap": [1, 2],
  "rows": [
    ["Unit costs, long pairs", "Band doubling over bit-parallel tiles with the gap-chaining seed heuristic (A*PA2)", "Regrouped recurrence, staggered groups (3.2); agreeing bounds (3.3); neighbour construction (3.5); traceback within a memory limit (3.7)"],
    ["Unit costs, similar pairs", "Diagonal transition (WFA, BiWFA)", "Gathered extension (3.4)"],
    ["Gap-affine, two-piece", "Bidirectional wavefront (WFA, BiWFA); full sweep in difference recurrences for divergent pairs (Suzuki and Kasahara)", "Gathered extension (3.4); score-derived band (3.8); handover gated on work spent (S1)"],
    ["Substitution matrices, local", "Anti-diagonal vector sweep; end-first local alignment (SSW)", "Backward wavefront from the best end (3.8)"],
    ["Batches of short pairs", "Inter-sequence vectorisation (SWIPE, SeqAn, parasail)", "Certified per-lane bands (3.6)"],
  ]},
"ablation": {
  "caption": "Effect of each technique: median time with all techniques on (ms per alignment; per batch of 500,000 pairs for the short reads), and the change when one technique is switched off, over three alternating rounds on the Skylake-X at commit `dc980d7`, with branches kept clear of 32-byte boundaries. Every configuration returned the same costs.",
  "header": ["Dataset", "All on", "Regrouped recurrence", "Paired groups", "Gathers", "Agreement rule", "Seed neighbours", "Certified band"],
  "align": "lrrrrrrr",
  "rows": [
    ["ont-500k", "98.8", "+1.1%", "**+4.5%**", "+0.5%", "−0.2%", "**+8.9%**", ""],
    ["ont-500k-genvar", "136.0", "+1.9%", "**+5.0%**", "+0.4%", "0.0%", "**+16.4%**", ""],
    ["ont-1k", "0.021", "+0.2%", "+0.3%", "**+14.9%**", "0.0%", "+0.7%", ""],
    ["ont-10k", "0.145", "−0.6%", "+1.2%", "**+6.4%**", "**+22.2%**", "+0.1%", ""],
    ["SARS-CoV-2", "0.197", "−0.9%", "+0.2%", "**+13.4%**", "**+4.9%**", "−0.5%", ""],
    ["100 kbp, 1%", "1.10", "−1.4%", "+0.4%", "**+26.5%**", "−0.1%", "−1.0%", ""],
    ["100 kbp, 5%", "4.61", "−0.6%", "−0.2%", "+0.1%", "0.0%", "+0.2%", ""],
    ["100 kbp, 10%", "6.72", "+0.3%", "+3.4%", "+1.9%", "+0.5%", "−0.2%", ""],
    ["100 kbp, 20%", "15.0", "**+4.1%**", "**+11.5%**", "+3.9%", "−2.9%", "0.0%", ""],
    ["Short reads, affine", "232", "", "", "", "", "", "**+188%**"],
    ["Short reads, unit", "213", "", "", "", "", "", "**+215%**"],
  ]},
"unit": {
  "fit": True,
  "caption": "Unit-cost alignment of the real datasets by every exact aligner: mean time per alignment in ms, traceback included, and below it the peak resident memory of the run in MB, as QuickEd reports it. t.o. marks a tool that completed no pair within 30 s and o.o.m. one the 8 GB limit stopped. The simulated datasets are in Fig. 2 and Supplementary Table S4.",
  "header": ["Dataset", "dinara-align", "A*PA2-full", "A*PA2-simple", "A*PA", "QuickEd", "Edlib", "BiWFA", "WFA", "BSAlign", "KSW2"],
  "align": "lrrrrrrrrrr",
  "rows": [
    ["*Mean time per alignment*", "", "", "", "", "", "", "", "", "", ""],
    ["ont-1k (830 bp)", "**0.021**", "0.081", "0.048", "0.612", "0.075", "0.124", "0.044", "0.031", "0.119", "0.539"],
    ["ont-10k (3.6 kbp)", "**0.156**", "0.391", "0.247", "15.3", "0.688", "1.09", "0.741", "0.599", "1.30", "13.1"],
    ["ont-50k (9.6 kbp)", "**0.580**", "1.27", "0.944", "263", "2.31", "6.14", "7.43", "6.38", "12.3", "167"],
    ["ont-500k (642 kbp)", "**118**", "283", "687", "t.o.", "847", "5,000", "6,330", "o.o.m.", "t.o.", "o.o.m."],
    ["ont-500k-genvar (663 kbp)", "**156**", "225", "511", "t.o.", "1,390", "5,230", "6,200", "o.o.m.", "t.o.", "o.o.m."],
    ["SARS-CoV-2 (30 kbp)", "**0.286**", "2.00", "0.728", "6.49", "3.90", "8.18", "0.909", "0.987", "47.1", "68.5"],
    ["Illumina (248 bp)", "**0.000576**", "0.014", "0.00410", "0.014", "0.016", "0.023", "0.00204", "0.00112", "0.021", "0.012"],
    ["PacBio HiFi (14 kbp)", "**0.013**", "0.977", "0.259", "1.05", "0.950", "3.35", "0.086", "0.066", "12.2", "3.09"],
    ["*Peak resident memory (MB)*", "", "", "", "", "", "", "", "", "", ""],
    ["ont-1k", "14", "5", "5", "6", "2", "7", "7", "8", "2", "7"],
    ["ont-10k", "16", "6", "6", "39", "6", "8", "8", "39", "25", "72"],
    ["ont-50k", "19", "8", "7", "858", "18", "8", "10", "187", "423", "650"],
    ["ont-500k", "64", "89", "90", "t.o.", "26", "15", "29", "o.o.m.", "t.o.", "o.o.m."],
    ["ont-500k-genvar", "64", "57", "55", "t.o.", "27", "16", "27", "o.o.m.", "t.o.", "o.o.m."],
    ["SARS-CoV-2", "18", "7", "6", "19", "15", "8", "9", "26", "215", "151"],
    ["Illumina", "13", "5", "5", "5", "2", "7", "7", "7", "2", "7"],
    ["PacBio HiFi", "13", "7", "6", "6", "4", "8", "8", "9", "127", "21"],
  ]},
"affine": {
  "caption": "Gap-affine alignment at costs (4, 6, 2) by every exact aligner that takes them: mean time per alignment in ms, a 30 s budget per sample. Superscripts give the pairs completed where not all were; t.o. marks a tool that completed none and o.o.m. one the 8 GB limit stopped. Every dataset is in Supplementary Table S3.",
  "header": ["Dataset", "dinara-align", "WFA", "BiWFA", "BSAlign", "KSW2"],
  "align": "lrrrrr",
  "rows": [
    ["ont-1k", "**0.105**", "0.188", "0.308", "0.336", "1.10"],
    ["ont-10k", "**1.60**", "4.89", "5.86", "7.84", "30.0"],
    ["ont-50k", "**21.5**", "63.4", "54.9", "87.0", "341^81/104"],
    ["ont-500k", "t.o.", "o.o.m.", "t.o.", "o.o.m.", "o.o.m."],
    ["ont-500k-genvar", "**22,000^1/4**", "o.o.m.", "t.o.", "o.o.m.", "o.o.m."],
    ["SARS-CoV-2", "**0.710**", "2.05", "1.78", "341", "87.9"],
    ["Illumina", "**0.00145**", "0.00155", "0.00281", "0.054", "0.013"],
    ["PacBio HiFi", "**0.095**", "0.227", "0.231", "86.0", "9.72"],
    ["10 kbp, 1%", "**0.199**", "0.416", "0.646", "37.0", "21.3"],
    ["10 kbp, 5%", "**2.32**", "7.14", "8.72", "37.0", "96.1"],
    ["10 kbp, 10%", "**7.06**", "23.1", "27.1", "37.3", "202"],
    ["10 kbp, 20%", "**20.1**", "66.0", "75.1", "37.4", "375^81/99"],
    ["100 kbp, 1%", "**11.9**", "39.9", "39.0", "o.o.m.", "2,420"],
    ["100 kbp, 5%", "**331**", "808", "714", "o.o.m.", "9,170^3/9"],
  ]},
}

FIGURES = {
"overview": "How dinara-align assigns a pair to an algorithm. At unit costs, diagonal transition from one end and then from both settles near-identical pairs while it costs less than a band would; other pairs go to band doubling over bit-parallel tiles, pruned by the seed heuristic on long pairs, and are traced back a tile at a time within the memory limit. At gap-affine costs, the bidirectional wavefront hands a pair to a full sweep in difference recurrences once its projected work exceeds the sweep's, the alignment being the same either way. Batches of short pairs are aligned one pair per vector lane within a band whose optimality each lane certifies.",
"sweeps": "Unit-cost alignment of QuickEd's simulated pairs, mean time per alignment on a logarithmic scale. (a) 100 kbp pairs by error rate. (b) Pairs at 20% errors by length (5% in the interactive version); lines end where a tool completed no pair within 30 s, and hollow points mark tools that completed only some pairs. dinara-align had the lowest mean time at every length and error rate.",
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
 "Gao2021": ("Gao et al.", "2021"), "Doblas2025": ("Doblas et al.", "2025"), "Shao2024": ("Shao and Ruan", "2024"), "Gibrat2018": ("Gibrat", "2018"), "Fujiki2020": ("Fujiki et al.", "2020"),
 "Zook2016": ("Zook et al.", "2016"), "Hansen2025": ("Hansen et al.", "2025"),
}
# The bibliography keys these map to in references.bib.
BIBKEYS = {
 "Gotoh1982": "Gotoh1982", "Ukkonen1985": "Ukkonen1985", "Myers1999": "Myers1999bitvector", "Sosic2017": "Sosic2017edlib",
 "Myers1986": "Myers1986ond", "MarcoSola2021": "MarcoSola2021wfa", "MarcoSola2023": "MarcoSola2023biwfa",
 "GrootKoerkamp2024a": "GrootKoerkamp2024astarpa", "GrootKoerkamp2024b": "GrootKoerkamp2024astarpa2", "Wozniak1997": "Wozniak1997",
 "Farrar2007": "Farrar2007striped", "Zhao2013": "Zhao2013ssw", "Suzuki2018": "Suzuki2018diff", "Li2018": "Li2018minimap2",
 "Daily2016": "Daily2016parasail", "Liu2023": "Liu2023blockaligner", "Rognes2011": "Rognes2011swipe", "Rahn2018": "Rahn2018seqan",
 "Kallenborn2026": "Kallenborn2026accelign", "Modular": "Modular2025mojo", "Lattner2021": "Lattner2021mlir",
 "Eizenga2022": "Eizenga2022wfalm", "hyalite": "hyalite", "Gao2021": "Gao2021abpoa", "Doblas2025": "Doblas2025quicked", "Shao2024": "Shao2024bsalign", "Gibrat2018": "Gibrat2018band", "Fujiki2020": "Fujiki2020seedex",
 "Zook2016": "Zook2016giab", "Hansen2025": "Hansen2025hg002",
}
