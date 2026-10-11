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
  'Pairwise alignment underlies read mapping, genome assembly, variant calling and pangenome construction. '
  'Heuristic aligners can miss or misplace alignments with large indels or clustered errors, whereas exact '
  'aligners always return an optimal alignment, a reproducible result and a reference for evaluating heuristics. Yet no single exact aligner is fastest across workloads: diagonal '
  'transition wins on short and near-identical pairs, seeded band '
  'doubling on long nanopore reads, and SIMD methods on batches of short pairs.'),
 ('Results',
  'We present dinara-align, an exact aligner that chooses its algorithm per pair, making every choice by a rule '
  'that either provably preserves optimality or bounds the work it can waste. A gap-affine alignment starts with '
  'one search and adds a second only when its measured work exceeds what two would need, which we prove loses no '
  'optimal alignment; searches abandon a projection only after spending a fixed share of the alternative\'s cost; '
  'short pairs are batched in narrow bands certified per vector lane; and inexact seed neighbours are inserted '
  'without search, provably losing none that pruning keeps. Under costs, every engine returns the same alignment '
  'among equally optimal ones. On an Intel Skylake-X processor, dinara-align was the fastest of the exact aligners compared on '
  'every dataset, from 250 bp reads to 1 Mbp pairs and from 0 to 20% divergence: up to 2.4 '
  'times faster than A*PA2 on ultra-long nanopore reads, 2.0 and 5.2 times faster than WFA on Illumina and HiFi '
  'reads, and 1.1 to 3.8 times faster than WFA and BiWFA at gap-affine costs.'),
 ('Availability and implementation',
  'dinara-align is free software under the Mozilla Public License 2.0, available at '
  'https://github.com/kokyriakidis/dinara-align with C, C++, Python and command-line interfaces.'),
 ('Contact', 'kkyriaki@ucsc.edu, bpaten@ucsc.edu'),
 ('Supplementary information', 'Supplementary data are available at *Bioinformatics* online.')]

BODY = [
('sec', 'Introduction'),
('p',
 'Pairwise alignment finds the least costly series of substitutions, insertions and deletions that transforms '
 'one sequence into another. It is a core component of genome analysis, applied to read mapping, variant '
 'detection, de novo genome assembly and multiple sequence alignment [@Prjibelski2019]: read mappers align each '
 'read to the reference window its seeds select and extend alignments from seed matches, assemblers align reads '
 'to one another to find their overlaps [@Li2018], and assemblies, haplotypes and viral genomes are compared '
 'base by base.'),
('p',
 'Read lengths have grown from hundreds of bases to hundreds of kilobases, and nanopore reads still carry 5 to '
 '15% errors. The classical dynamic program computes an optimal alignment over an \\(n \\times m\\) matrix '
 '[@Gotoh1982], practical for short reads but not for long ones: one pair of 500 kbp reads fills \\(2.5 \\times '
 '10^{11}\\) cells, and the total cost grows with both the number of reads and their length.'),
('p',
 'Heuristic aligners therefore compute only part of the matrix, a fixed band, a series of windows, or the cells '
 'before the score drops by more than a threshold (X-drop, and Z-drop in minimap2 [@Li2018]), and lose the '
 'guarantee of optimality: their result depends on these cut-offs, and alignments with large indels or clustered '
 'errors can be missed or misplaced. Exact aligners instead skip only those cells that provably cannot lie on an '
 'optimal path, so that their running time depends on the distance between the sequences rather than on the '
 'product of their lengths.'),
('p',
 'Three ideas underlie current exact aligners. *Band doubling* computes only the diagonals a path within a cost '
 'bound can reach, doubling the bound until the optimum lies within it [@Ukkonen1985]; *bit-parallelism* packs '
 '64 cells of a column into a machine word [@Myers1999], and Edlib combines the two [@Sosic2017]. *Diagonal '
 'transition* extends the furthest-reaching paths one cost at a time, following matches at no cost '
 '[@Ukkonen1985; @Myers1986]; WFA generalised it to gap-affine costs [@MarcoSola2021] and BiWFA made its memory '
 'linear in the distance [@MarcoSola2023]. A*PA added A* search with a gap-chaining seed heuristic '
 '[@GrootKoerkamp2024a], and A*PA2 combined this heuristic with band doubling over SIMD blocks, aligning '
 'ultra-long reads up to 19 times faster than Edlib and BiWFA [@GrootKoerkamp2024b]; QuickEd bounds the distance '
 'heuristically and then aligns within the bound [@Doblas2025]. A second line of work vectorises the full '
 'dynamic program, with anti-diagonal [@Wozniak1997], striped [@Farrar2007; @Zhao2013] and prefix-scan '
 '[@Daily2016] layouts, difference recurrences [@Suzuki2018; @Li2018], adaptive blocks [@Liu2023] and '
 'inter-sequence parallelism [@Rognes2011; @Rahn2018], at a cost that grows with the matrix or the band rather '
 'than with the distance.'),
('p',
 'Each of these methods is fastest on part of the range of workloads: diagonal transition, whose work grows with '
 'the square of the distance, on short and near-identical pairs; seeded band doubling, whose work grows with the '
 'length, on long divergent reads; and inter-sequence vectorisation on batches of short pairs. An aligner that '
 'is fastest across workloads must therefore choose among them for each pair, and an exact aligner must make '
 'every such choice, which algorithm, which bound and when to abandon a search, without giving up optimality. '
 'dinara-align makes each choice by a rule that either provably preserves optimality or bounds the work the '
 'choice can waste, and implements each algorithm for the latency and memory behaviour of current processors. '
 'Our contributions are:'),
('list',
 ['**One wavefront search or two, chosen by measured work** (Section 3.9). A gap-affine alignment starts with one '
  'search and adds a second only when its measured work exceeds what two would need, without losing any optimal '
  'alignment (Proposition 2).',
 "**Projections used with guarantees** (Sections 3.3 and 3.10). A projected distance sets a band's first bound "
  'only when two independent estimates agree, and a search abandons a projection only after spending a fixed '
  "share of the alternative's cost.",
 '**Certified bands for inter-sequence batches** (Section 3.6). Each short pair is aligned in a vector lane '
  'within a narrow band that the lane certifies optimal (Proposition 1), adapting the speculate-and-verify '
  'approach of SeedEx [@Fujiki2020] to software.',
 '**Direct neighbour construction for the seed heuristic** (Section 3.5). Inexact neighbours of exact seed '
  'matches are inserted without search, and every neighbour that local pruning would keep is provably among them '
  '(Lemma 1), so the heuristic remains a lower bound.',
 '**One alignment whichever engine computes it** (Sections 3.1 and 3.7). Under costs, every engine follows the '
  'same rule among equally optimal alignments, and a unit-cost band past the memory limit recomputes stretches of '
  'its traceback from recorded state rather than splitting the pair, so its alignment is the same under any '
  'memory limit.',
 '**Engines arranged for the processor** (Sections 3.2, 3.4 and 3.11): gathered diagonal transition, an '
  'interleaved bit-parallel sweep with a shorter dependency chain, and fewer constant costs per pair on short '
  'reads.',
 '**One library for every mode and cost model** ([#tab:engines]): global, ends-free, extension and local '
  'alignment under unit, gap-affine, two-piece and substitution-matrix costs, for single pairs and batches, '
  'through C, C++, Python and command-line interfaces.']),
('p',
 'We evaluate dinara-align on the A*PA2 and QuickEd benchmarks, to which we add real Illumina and PacBio HiFi '
 'reads paired with the genome they came from. Sections 2 to 5 give the background, the methods, the results and '
 'the discussion.'),
('sec', 'Background'),
('sub', 'Alignment and costs'),
('p',
 'Let \\(A=a_1\\dots a_n\\) be the reference and \\(B=b_1\\dots b_m\\) the query, and let cell \\((i,j)\\) lie on '
 'diagonal \\(j-i\\). Under gap-affine costs a mismatch costs \\(x\\) and a gap of \\(k\\) letters \\(g(k)=o+ke\\); unit '
 'costs have \\(x=e=1\\) and \\(o=0\\), two-piece costs take the minimum of two such functions [@Li2018], and a '
 'substitution matrix instead scores each pair of letters. Global alignment compares two whole sequences, such '
 'as two assemblies or a read and the reference stretch it maps to; ends-free alignment places a read inside a '
 'longer window or finds the overlap of two reads; extension grows an alignment from a seed match, as read '
 'mappers do [@Li2018]; and local alignment finds the best-matching parts of two sequences. An alignment is '
 '*exact* if it is optimal under the requested costs, mode, band and cost limit. At unit costs the distance '
 '\\(D(i,j)\\) between the prefixes \\(a_1\\dots a_i\\) and \\(b_1\\dots b_j\\) satisfies'),
('p',
 '\\[ D(i,j) = \\min\\bigl(D(i-1,j)+1,\\; D(i,j-1)+1, \\qquad D(i-1,j-1)+[a_i \\ne b_j]\\bigr), \\]'),
('p',
 'with \\(D(i,0)=i\\) and \\(D(0,j)=j\\), and the distance is \\(D(n,m)\\); gap-affine costs need three such layers '
 '[@Gotoh1982].'),
('sub', 'Band doubling and bit-parallel computation'),
('p',
 'An alignment of cost at most \\(t\\) visits only diagonals \\(d\\) with \\(|d| + |d - \\delta| \\le t\\), where '
 '\\(\\delta = m - n\\); band doubling computes only those cells, doubling \\(t\\) until \\(D(n,m) \\le t\\) '
 "[@Ukkonen1985]. Myers' recurrence advances 64 cells of a column per machine word [@Myers1999], and A*PA2 "
 'sweeps the band in tiles of 256 columns of vector registers, dropping rows that a score already rules out '
 '[@GrootKoerkamp2024b].'),
('sub', 'Diagonal transition'),
('p',
 'Diagonal transition keeps, for each cost \\(s\\) and diagonal, the furthest cell an alignment of cost \\(s\\) '
 'reaches, extended along matches at no cost, in \\(O(n + s^2)\\) steps [@Ukkonen1985; @Myers1986]; WFA '
 'generalises it to gap-affine costs with three wavefronts per cost [@MarcoSola2021], and BiWFA runs it from '
 'both ends [@MarcoSola2023].'),
('sub', 'The seed heuristic'),
('p',
 'A*PA divides the reference into disjoint seeds of \\(k\\) bases and finds their matches in the query, exactly or '
 'with one edit [@GrootKoerkamp2024a]. Its gap-chaining seed heuristic bounds the cost from a cell to the end of '
 'the alignment by the seeds ahead that no chain of matches from that cell covers, counting the gaps between '
 'chained matches, and never overestimates. Local pruning discards a match \\(\\nu\\) unless \\(c(\\nu) + h(e(\\nu)) < '
 '\\tau\\), where \\(c(\\nu)\\) is its cost, \\(e(\\nu)\\) its end cell, \\(h(u)\\) the least cost of a path from \\(u\\) '
 'across the following seeds, and \\(\\tau\\) a threshold that depends only on the seeds. A*PA2 uses the heuristic '
 'to drop the rows of each band whose score plus the heuristic exceeds the bound [@GrootKoerkamp2024b].'),
('sec', 'Methods'),
('sub', 'Overview'),
('p',
 'Each problem is assigned to an algorithm ([#tab:engines]; [#fig:overview]). Unit costs follow A*PA2 '
 '[@GrootKoerkamp2024b]: band doubling over 256-column tiles of bit-parallel words, pruned by the gap-chaining '
 'seed heuristic [@GrootKoerkamp2024a], with near-identical pairs aligned by diagonal transition alone. '
 'Gap-affine costs use the wavefront of WFA and BiWFA [@MarcoSola2021; @MarcoSola2023]. A short '
 'diagonal-transition search estimates the distance, which with the length and the cost model selects the '
 'algorithm.'),
('p',
 'Under costs every engine selects the same alignment among equally good ones, following the backtrace order of '
 'WFA2-lib: from the end, an edit is taken into a cell whenever one is optimal, a substitution before a deletion '
 'before an insertion, and a gap is extended before it is opened. The tile traceback applies this order cell by '
 'cell (Section 3.7), the batches and the difference sweep store per-cell flags under it, and the wavefronts '
 'take the furthest-reaching source, which makes the same choice; leftmost indels come from the same rule '
 'applied to both sequences reversed. Two exceptions remain: a gap-affine pair whose wavefronts would exceed the '
 'memory limit is split and follows the rule within each piece, and alignments under a substitution matrix '
 "follow the traceback order of Gotoh's recurrence."),
('figure', 'overview'),
('table', 'engines'),
('sub', "A shorter critical path for Myers' recurrence"),
('p',
 "As in Edlib and A*PA2 [@Sosic2017; @GrootKoerkamp2024b], the query is packed 64 bases per word, and Myers' "
 'recurrence [@Myers1999] advances a word by one column. With match mask \\(E\\), vertical difference masks \\(P_v, '
 "M_v\\) and the incoming negative horizontal difference folded into \\(E' = E \\lor M_h^{\\mathrm{in}}\\), A*PA2 "
 'computes'),
('p',
 "\\[ c = (E' \\land P_v) + P_v, \\qquad X_h = (c \\oplus P_v) \\lor E', \\qquad P_h = M_v \\lor \\lnot(X_h \\lor P_v), "
 '\\qquad M_h = P_v \\land X_h. \\]'),
('p',
 'The horizontal masks \\(P_h\\) and \\(M_h\\) determine the vertical masks of the next column, so the path from '
 '\\(c\\) to them lies on the loop-carried chain. Since \\((c \\oplus P_v) \\lor P_v = c \\lor P_v\\) and \\(P_v \\land '
 '(c \\oplus P_v) = P_v \\land \\lnot c\\), the same masks are obtained as'),
('p',
 "\\[ P_h = M_v \\lor \\lnot\\bigl(c \\lor (E' \\lor P_v)\\bigr), \\qquad M_h = (P_v \\land \\lnot c) \\lor (P_v \\land E'), "
 '\\]'),
('p',
 "where \\(E' \\lor P_v\\) and \\(P_v \\land E'\\) do not depend on \\(c\\) and are computed off the critical path. This "
 'removes two operations from the loop-carried chain and, for a single word, reduces the cost from 14 to 8 '
 'cycles per column on the Skylake-X. Myers [@Myers1999], Edlib and A*PA2 use the original form, and we are not '
 'aware of a previous description of the regrouped one.'),
('p',
 'Words occupy the lanes of a vector, staggered by one column, so that the horizontal difference leaving one '
 'word enters the next through a lane rotation [@Wozniak1997; @GrootKoerkamp2024b], and each column waits for '
 'the rotation of the previous one. A*PA2 overlaps these waits with two four-lane vectors; on AVX-512, where a '
 'vector holds eight words, dinara-align interleaves two eight-word groups, the lower ten columns behind the '
 'upper, so that each group fills the latency of the other, and bands of two or three words run in '
 'general-purpose registers, whose operations have lower latency. The sweep takes about 2.4 cycles per word and '
 'column.'),
('sub', 'Band bounds from agreeing estimates'),
('p',
 'Band doubling is exact for any initial bound, but its running time depends on how close that bound is to the '
 'distance. A computation that has reached cost \\(s\\) at anti-diagonal \\(r\\) extrapolates the distance \\(\\hat d '
 '= s(n+m)/r\\), which is reliable when edits are spread evenly and too high by orders of magnitude when they '
 'gather, as at the ends of nanopore reads. dinara-align therefore sets the first bound from two extrapolations, '
 '\\(\\hat d_1\\) from the first eight edits of a diagonal-transition search and \\(\\hat d_2\\) from the same search '
 'at its budget, only if \\(\\max(\\hat d_1,\\hat d_2) \\le \\tfrac43 \\min(\\hat d_1,\\hat d_2)\\); otherwise, as in '
 'A*PA2, the bound starts 256 above the known lower bound and each retry at most doubles it. A failed round aims '
 'the next bound from where it stopped, its progress extrapolated to the whole width, and within a round the '
 'distance is extrapolated again after one eighth, one quarter and one half of the columns, lowering the bound '
 'or abandoning the round. These rules change only the sequence of bounds, so the result stays exact.'),
('p',
 'Within a round, the band reads the scores down the right edge of each tile and drops the rows at its top and '
 'bottom whose score plus the heuristic exceeds the bound, as A*PA2 does. A score changes by at most one per '
 'row, and the heuristic by at most its climb \\(c\\), one for exact seeds and two for inexact ones, so a row that '
 'exceeds the bound by \\(x\\) rules out the next \\(\\lceil x/(1+c) \\rceil\\) rows without reading them. The bottom '
 'of the next tile is bounded by the deepest row kept, since going \\(k\\) rows past it costs at least \\(k\\) more. '
 'With seeds, when at least one seed in 32 chains from the origin, the first bound starts 64 above the heuristic '
 'at the origin, which on long reads lies within a few percent of the distance.'),
('sub', 'Diagonal transition for near-identical pairs'),
('p',
 'Diagonal transition extends each diagonal of a wavefront along its run of matches. dinara-align compares eight '
 'letters at a time by the exclusive-or of two 64-bit words and a count of trailing zeros. On AVX-512, eight '
 'diagonals are extended together: their letters are loaded with gather instructions, and the wavefront is '
 'padded with unreachable diagonals so that every step processes full vectors under branch-free masks. The '
 'gather masks are derived from the lanes of each group, which avoids a dependency between consecutive gathers '
 'that would otherwise serialise them. The same vectorised step serves near-identical pairs, computed from one '
 'or both ends [@MarcoSola2023], the traceback inside each tile, restricted to the diagonals that can still '
 'reach the traced cell, and the three-layer gap-affine wavefront.'),
('p',
 'A pair whose diagonal transition finishes within its budget is aligned from its wavefronts alone. From one end '
 'the wavefronts are kept whole. From both ends, as in BiWFA, the forward wavefronts are grown on to the proven '
 'distance, each cost kept only on the diagonals on which the backward wavefronts show that an optimal path '
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
 'An exact match \\(\\mu\\) of a seed at query rows \\([i, i+16)\\) has four neighbouring matches with one edit, '
 'fewer at the ends of the query: two share its start, at rows \\([i, i+15)\\) and \\([i, i+17)\\), and two share '
 'its end, at \\([i+1, i+16)\\) and \\([i-1, i+16)\\). Because chaining is all-or-nothing, a neighbour one diagonal '
 'away may chain where \\(\\mu\\) cannot, and the neighbours are required for the heuristic to remain a lower '
 'bound. Rather than locating them by search and testing each for local pruning, dinara-align inserts them '
 'whenever \\(\\mu\\) is retained, which is justified by the following lemma. Local pruning retains a match \\(\\nu\\) '
 'only if \\(c(\\nu) + h(e(\\nu)) < \\tau\\) (Section 2.4).'),
('p',
 '**Lemma 1.** If an exact match \\(\\mu\\) fails local pruning, so does each of its four neighbours.'),
('p',
 '*Proof.* Each neighbour matches the same seed as \\(\\mu\\), so it is tested against the same \\(\\tau\\). A '
 'neighbour \\(\\nu\\) that shares the end of \\(\\mu\\) has \\(e(\\nu)=e(\\mu)\\) and \\(c(\\nu)=1>c(\\mu)=0\\). A neighbour '
 'that shares the start ends on an adjacent diagonal, and one insertion or deletion from \\(e(\\mu)\\) reaches a '
 'cell on that diagonal at or beyond \\(e(\\nu)\\). The least cost from a cell to beyond the following seeds does '
 "not increase further along its diagonal, by Ukkonen's lemma applied to the reversed sequences [@Ukkonen1985], "
 'so \\(h(e(\\mu)) \\le 1 + h(e(\\nu))\\) and \\(c(\\mu)+h(e(\\mu)) \\le c(\\nu)+h(e(\\nu))\\). In both cases the pruning '
 'test value of \\(\\nu\\) is at least that of \\(\\mu\\). \\(\\square\\)'),
('p',
 'Every neighbour that local pruning would keep is therefore inserted; one it would have discarded can only '
 'lower the heuristic, which remains a lower bound, so the alignment stays exact. Neighbours sharing the start '
 'raise the layer of that start by the chain length available from their ends, and neighbours sharing the end '
 'are placed one layer below it, without a search or a pruning test of their own.'),
('sub', 'Certified bands for inter-sequence batches'),
('p',
 'Batches of short pairs are aligned with one pair per vector lane [@Rognes2011; @Rahn2018], first in saturating '
 '8-bit lanes and, for the scores of pairs that saturate, in 16-bit lanes [@Zhao2013; @Daily2016]; alignments '
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
 'letters to end on \\(\\delta\\). These form at least one gap of each kind, and since \\(g\\) is subadditive, as is '
 'the minimum of two such functions for two-piece costs, their total cost is at least \\(g(d)+g(d-\\delta)=L(d)\\). '
 'Every alignment outside the band therefore costs at least \\(\\min\\bigl(L(d_h+1), L(d_\\ell-1)\\bigr)\\). '
 '\\(\\square\\)'),
('p',
 'Such bounds underlie band doubling [@Ukkonen1985] and verify narrow-band results in software [@Gibrat2018] and '
 'hardware [@Fujiki2020]; here each vector lane of an inter-sequence batch evaluates its own certificate, and '
 'the strict form guarantees that the alignment the tie rule selects lies within the band. Pairs that fail, '
 'about one in seven in our short-read benchmark, are recomputed over every diagonal a cheaper path could visit, '
 'grouped by the width of that band so that the lanes of a group finish together.'),
('sub', 'Traceback within a memory limit'),
('p',
 'The traceback of a band traces each tile from the boundary recorded on its left, the vertical differences down '
 'that boundary. On a 1 Mbp pair at 20% errors these boundaries take up to 175 MB, so they are held to the '
 "alignment's memory limit, 80 MB by default: beyond it only every \\(k\\)-th is kept, \\(k\\) doubling as the round "
 'proceeds. Each tile also records the state the band carries into it, its first and last words, the score at '
 'its top, the deepest row kept and its score, and the bound, and the traceback computes a stretch of tiles '
 'again from the last kept boundary. Every step of a band depends only on this state and on the scores it '
 'computes, so each boundary is recovered as the round left it and the alignment is the same byte for byte; '
 'splitting the pair instead, as QuickEd does [@Doblas2025], resolves ties within each part. On the 1 Mbp pairs '
 'at 20% errors, the growth of peak memory fell from 175 to 133 MB on the worst pair, for 9% more time; pairs '
 'whose boundaries fit are unchanged.'),
('sub', 'Gap-affine and score-based alignment'),
('p',
 'A scoring scheme with match score \\(a\\), mismatch score \\(b\\) and gap scores \\(-(o+ke)\\) can be converted into '
 'costs for global alignment [@Eizenga2022]. Every alignment covers \\(n+m\\) letters, so its score \\(S\\) and a '
 'gap-affine cost \\(C\\) satisfy'),
('p',
 '\\[ S = \\tfrac12\\bigl(a(n+m) - C\\bigr), \\qquad C = 2(a-b)\\,X + \\sum_{\\text{gaps}} \\bigl(2o + k(2e+a)\\bigr), \\]'),
('p',
 'where \\(X\\) is the number of mismatches, so the wavefront computes score-maximising alignments. The cost also '
 'bounds the band in which the alignment must be stored: each gapped letter contributes at least \\(2e+a\\) to '
 '\\(C\\), and an alignment that deviates by \\(t\\) diagonals from both the main and the end diagonal contains at '
 'least \\(2t\\) gapped letters, so only diagonals within \\(C/(2(2e+a))\\) of them are stored, 137 instead of 1,000 '
 'for a 1 kbp pair at 5% divergence. Local alignment locates the end of the optimal alignment with a forward '
 'sweep, as in SSW [@Zhao2013], and extends a wavefront back from it to the first cell attaining the maximal '
 'score, which avoids the correction pass that striped layouts may need [@Farrar2007].'),
('p',
 'The wavefront keeps five bytes per diagonal for the traceback, the furthest column and one byte naming each '
 "layer's source, and a pair whose fronts would exceed the memory limit is split where an optimal path crosses, "
 'as in BiWFA. Whether one search or two align a pair is decided while it is aligned (Section 3.9); a cost '
 'without traceback stops once its best lies within one of the costs both searches have checked against the '
 'start of the other; and a global pair projected to cost more than a full sweep in 8-bit difference recurrences '
 '[@Suzuki2018] is handed to that sweep (Section 3.10), which keeps the predecessor flags of Section 3.1, so '
 'that unrelated 5 kbp pairs aligned in 4.8 ms instead of 32.9 ms.'),
('p',
 'Under a substitution matrix, anti-diagonals are computed in 16-bit lanes while the scores provably fit '
 '[@Wozniak1997], in strips of rows of at most 16 KB that stay in the first-level cache, which took local '
 'alignment of 10 kbp against 12 kbp from 50.0 to 17.4 ms. A global alignment is traced in a band of diagonals '
 'that grows at most fourfold per step until its score proves that every optimal path lies inside it, so a pair '
 'close to its diagonal costs about its own cells: 1 kbp pairs at 10% divergence aligned in about 200 µs instead '
 'of 544 µs.'),
('sub', 'One search or two at gap-affine costs'),
('p',
 'WFA grows the wavefront from one end of the pair and keeps every front for the traceback [@MarcoSola2021]; '
 'BiWFA grows it from both ends and keeps only the latest fronts [@MarcoSola2023]. For an alignment, two '
 'searches must also prove where an optimal path crosses between them, so each must pass half the cost by a '
 'margin \\(m = o + w\\), where \\(o\\) is the widest gap opening and \\(w\\) the dearest single move: every optimal '
 'path has a cell whose costs from the two ends differ by at most \\(w\\), and a meeting inside a gap pays its '
 'opening once for both halves. With \\(W(c)\\) the number of diagonals a search steps to reach cost \\(c\\), a pair '
 'of cost \\(d\\) costs two searches about \\(2W((d+m)/2)\\) plus the setup of the second, and one search \\(W(d)\\). '
 'Where the fronts widen with the cost, \\(W(c)\\) grows with \\(c^2\\) and one search is cheaper while \\(d < '
 '(1+\\sqrt2)\\,m\\); where the ends of the sequences bound the fronts, it is always cheaper. On 250 bp Illumina '
 'reads most pairs cost less than \\(m\\), 14 at costs (4, 6, 2), and two searches took 2.8 µs per read against '
 '1.5 µs for WFA.'),
('p',
 'dinara-align therefore begins with the one search whose fronts the traceback reads, the backward search for '
 'the default tie rule, keeps its fronts, and traces the alignment back as soon as it reaches the far corner. '
 'Every cost \\(c\\) passed without reaching it shows that \\(d > c\\), and the search measures \\(W\\) as it grows: '
 'it continues alone while \\(c < m\\), costs the second search would pass in any case, and then while \\(W(c) \\le '
 '2W(\\lfloor (c+m)/2 \\rfloor)\\), that is, while it has stepped no more diagonals than two searches would need '
 'for the cheapest cost still possible. It also stops at half the memory limit or when a band leaves it no '
 'reachable cell. No threshold is tuned: the search hands over at about \\(2.4m\\) where the fronts widen with the '
 'cost, and never where the sequences bound them, aligning as WFA does.'),
('p',
 '**Proposition 2.** If the first search stops at cost \\(B\\) without reaching the far corner, and the second '
 'search, growing from cost 0, checks each new front against the last \\(w+1\\) fronts of the first, costs \\(B-w\\) '
 'to \\(B\\), while its own cost is at most \\(B\\), then the two searches miss no optimal alignment.'),
('p',
 '*Proof.* Let \\(P\\) be an optimal path, of cost \\(d\\). Having checked at every cost up to \\(B\\) that it had not '
 'reached the far corner, where the second search starts, the first search shows that \\(d > B\\). Along \\(P\\), '
 "the cost \\(g\\) from the first search's start never decreases and the cost \\(f\\) from the second's never "
 'increases, so \\(f - g\\) never increases. While its cost is at most \\(B\\), the second search checks each front '
 'against costs \\(B-w\\) to \\(B\\); afterwards the two alternate as in BiWFA, where a cell with \\(|f-g| \\le w\\) is '
 'checked once both searches have passed it. Since \\(g\\) rises in steps of at most \\(w\\) and \\(d > B\\), some '
 'cell of \\(P\\) has \\(g \\in [B-w, B]\\). If \\(f - g < -w\\) there, then \\(f < B - w\\), and the cell is checked in '
 'the first phase. Otherwise, as \\(f - g\\) never increases, \\(P\\) has a cell with \\(|f-g| \\le w\\) and \\(g \\ge '
 'B-w\\): either \\(g \\le B\\) and \\(f \\le B + w\\), checked against the costs the first search still holds, or \\(g '
 '> B\\), a cell of the alternating phase. In the first case \\(f < (d+o-w)/2\\), as \\(f + g \\le d + o\\), and \\(g\\) '
 'is already held; in the second both costs are at most \\((d+o+w)/2\\). The stopping rule ends the searches only '
 'once both have passed half the best cost found plus the margin, so the cell is checked first. \\(\\square\\)'),
('p',
 'An alignment traced from the unpruned fronts of the single search is the one the two searches would have '
 'traced after meeting: the traceback steps only from a cell to a predecessor on an optimal path, and the fronts '
 'that the two-search path prunes away hold no such cell, so the tie rule is unchanged. A distance, which needs '
 'no traceback, grows the forward search alone through the costs its ring can hold. Starting with one search '
 'took a 250 bp Illumina read at (4, 6, 2) from 2.41 to 1.52 µs, and the measured rule saved a further 5%, a '
 'third on the pairs costing 9 to 32; nanopore and HiFi reads were unchanged.'),
('sub', 'Giving up on a projection, with bounded waste'),
('p',
 'Two hand-overs rest on projections: diagonal transition gives way to band doubling, and the gap-affine '
 'wavefront to the difference sweep, once continuing is projected to cost more than the alternative. '
 'Extrapolated from the edits seen so far, a projection is reliable when edits are spread evenly and can exceed '
 'the distance a hundredfold when they gather at the ends, as on nanopore reads with adapters or poly-A tails, '
 'handing a pair that the search would finish in microseconds to a band over the whole matrix. dinara-align '
 'gives up on a projection only once the work already spent has reached one twenty-fourth of the cheapest '
 'alternative that could still be needed: for diagonal transition, the band whose bound is the cost already '
 'reached, a lower bound on any band because the distance is at least that cost; for the wavefront, the full '
 'sweep. Waiting thus costs a search that gives up in the end at most about 4% of what the alternative costs in '
 'any case, while a search that would have finished is no longer abandoned on a projection from a handful of '
 'edits. For an alignment, diagonal transition also widens a projection from \\(k\\) edits by \\(1.5/\\sqrt{k}\\), '
 'about its relative error. On nanopore reads at gap-affine costs this took the 10 kbp reads from 2.33 to 1.55 '
 'ms and the 50 kbp reads from 30.5 to 21.0 ms (Supplementary Table S1).'),
('sub', 'Constant costs on short reads'),
('p',
 'A 250 bp read aligns in well under a microsecond, so costs that do not depend on the pair dominate. An aligner '
 'object keeps every buffer its searches use, the encoded sequences, the rings and the fronts traced back, from '
 "pair to pair, one aligner per thread, as each competitor's benchmark runner keeps its own. Both sequences are "
 'encoded, and reversed for the leftmost tie rule, in one vectorised pass of 32 bytes with an overlapping tail, '
 'and results pass between internal calls as plain values rather than optional ones, which the compiler moves a '
 'byte at a time. Together these took a unit-cost alignment of a 250 bp Illumina read through one aligner from '
 '0.97 to about 0.40 µs (Supplementary Table S1).'),
('sub', 'Thresholds measured against the code'),
('p',
 'The choice of algorithm and the bounds of each search rest on 41 thresholds, each set from measurements on the '
 'code of its time. We built dinara-align with each halved and with each doubled and timed the 83 builds against '
 'the unchanged one (Supplementary Section S1): about thirty moved no sample by more than 3%, and most others '
 'traded one dataset against another, but three choices had gone stale, among them the length from which seeds '
 'match with one edit, which had cost 100 kbp pairs at 10% errors a third of their running time.'),
('sub', 'Validation'),
('p',
 'A differential fuzzer compares every mode, cost model, band, cost limit and tie rule on random pairs with a '
 'full-matrix implementation that shares no code with the library, and the WFA2-lib regression set is reproduced '
 'exactly, every CIGAR included (Supplementary Section S1).'),
('sec', 'Results'),
('sub', 'Experimental setup'),
('p',
 'All measurements were taken on an Intel Core i9-7900X (Skylake-X, AVX-512), its ten cores at a fixed 3.3 GHz '
 'with turbo boost and hyper-threading disabled, with dinara-align compiled by Mojo 1.1.0 [@Modular] at commit '
 '`036e12b` (gap-affine costs at `a6d9166`, the other modes and batches at `289bcf0`) and every tool compiled '
 'for the native instruction set. We combined the protocols of A*PA2 [@GrootKoerkamp2024b] and QuickEd '
 "[@Doblas2025]. The real datasets are A*PA2's nanopore reads of about 1, 10, 50 and 500 kbp, the 500 kbp ones "
 'of CHM13 against its own reference and also of a second person against it, which adds genetic variation and '
 'gaps of up to 42 kbp; its SARS-CoV-2 genomes; and Illumina 2×250 bp and PacBio HiFi reads, the two '
 'technologies QuickEd adds, here reads of HG002 from Genome in a Bottle [@Zook2016] in one region of chromosome '
 "20, each mapped by minimap2 [@Li2018] to both haplotypes of HG002's own assembly [@Hansen2025] and paired with "
 "the stretch it covers, so that a pair differs by sequencing errors alone. The simulated datasets are QuickEd's "
 'random sequences of 10 kbp to 1 Mbp whose pairs differ by single-base substitutions, insertions and deletions '
 'spread uniformly at 1, 5, 10 and 20%, modelling no sequencing technology. Each dataset was sampled to about 2 '
 'Mbp and at least four pairs (Supplementary Table S6), and each tool aligned each pair once, with traceback, on '
 "one pinned core, within 30 s per sample, through one aligner object per sample, as every tool's benchmark "
 'runner keeps one. The exact aligners compared are A*PA2-full and A*PA2-simple, A*PA [@GrootKoerkamp2024a], '
 'Edlib [@Sosic2017], WFA and BiWFA [@MarcoSola2021; @MarcoSola2023], QuickEd, BSAlign [@Shao2024], KSW2 '
 '[@Li2018], parasail [@Daily2016] and SSW [@Zhao2013], each on the workloads it supports; every cost agreed '
 'with those of every other tool and with the published A*PA2 results.'),
('sub', 'Speed and memory'),
('p',
 'dinara-align had the lowest mean time per alignment of the exact aligners on every dataset. On the real data '
 '([#tab:unit]) it was 1.4 to 2.4 times faster than the fastest alternative on the nanopore reads and 2.5 times '
 'on the SARS-CoV-2 genomes, and 2.0 and 5.2 times faster than WFA, the fastest of the others there, on the '
 'Illumina and HiFi reads. On the simulated data ([#fig:sweeps]; Supplementary Table S5) it was 1.1 to 2.6 times '
 'faster: 100 kbp pairs at 1, 5, 10 and 20% errors took 1.5, 5.1, 7.3 and 15.7 ms, against 7.5, 8.4, 13.9 and '
 '36.7 ms for A*PA2-full and 12.5, 24.0, 39.0 and 62.8 ms for QuickEd, and its narrowest margin was on 10 kbp '
 'pairs at 5%, 0.267 against 0.288 ms for A*PA2-simple. Over all 24 datasets QuickEd took 3.0 to 75 times as '
 'long, 8.0 times in the geometric mean. At gap-affine costs dinara-align was 1.1 to 3.8 times faster than WFA '
 'and BiWFA wherever a tool completed a dataset ([#tab:affine]; [#tab:affine_sim]), its narrowest margin on the '
 'Illumina reads, 1.45 against 1.55 µs for WFA, and it completed more of the long simulated pairs than any other '
 'tool; in the other modes it was 1.0 to 22 times faster than the fastest tool offering the same mode, the one '
 'tie lying within the variation between builds. Global scores of 500,000 short-read pairs took 230 ms on one '
 'core and 27 ms on ten at (1, 2, 1), 4.6 and 4.3 times faster than WFA2-lib.'),
('table', 'unit'),
('table', 'affine'),
('table', 'affine_sim'),
('figure', 'sweeps'),
('p',
 'An alignment of an ultra-long read increased peak memory by at most 48 MB, against 82 MB for A*PA2-full, and '
 'of a 1 Mbp pair at 20% errors by at most 123 MB, against 184 MB, its traceback held within the default limit '
 "of 80 MB (Section 3.7; Supplementary Table S7); QuickEd stayed near 25 MB. On short reads dinara-align's peak "
 "is higher than the others' ([#tab:unit]), its runtime alone taking about 13 MB."),
('sec', 'Discussion'),
('p',
 'No single algorithm is fastest across the workloads evaluated here. dinara-align is fastest on all of them '
 'because it chooses among algorithms for each pair, and because every choice preserves optimality or bounds its '
 'waste, choosing never costs the guarantee or changes the alignment. Three observations stand out. First, on '
 'long pairs the techniques that reduce work outside the inner loop contributed more than those that speed it '
 'up: constructing seed neighbours directly saved up to 16%, interleaving bit-parallel groups 4.5 to 11.5% and '
 'the regrouped recurrence 1 to 4% (Supplementary Table S2). Second, on short reads the costs that do not depend '
 'on the pair dominate: before the changes of Section 3.11, a 250 bp alignment at gap-affine costs took a median '
 '4.8 µs through a fresh call and 2.2 µs through a reused aligner, so a benchmark should call every tool as its '
 'own runner does. Third, thresholds that encode a cost model go stale as the code changes, and measuring them '
 'again belongs in the evaluation of any aligner that selects its algorithm by thresholds (Section 3.12). The '
 'measured-work rule, the regrouped recurrence and the per-lane certificate do not depend on this implementation '
 'and apply directly to other wavefront, bit-parallel and inter-sequence aligners.'),
('p',
 'The study has several limitations:'),
('list',
 ['The evaluation covers one x86 processor with AVX-512; the AVX2 and NEON configurations are not evaluated here.',
 'Some thresholds, in particular the AVX-512 gates for inexact seeds, were set after measurements on uniform 100 '
  'kbp pairs like those of [#fig:sweeps]a, so results on those pairs may overstate performance on unseen data.',
 'As in the A*PA2 protocol, each pair was aligned once; on the Skylake-X some running times varied by up to a '
  'third between builds that did not change the code involved, which the ablation controls for (Supplementary '
  'Table S2) and the comparison with other tools does not.',
 'The Illumina and HiFi reads come from one region of chromosome 20 and exclude reads that minimap2 leaves more '
  'than a tenth clipped, such as those spanning structural variants.',
 'The ablation switches off one technique at a time, so it does not separate their interactions.',
 'Release builds target the baseline processor of each platform, so the AVX-512 code paths require a build for '
  'the host processor.']),
('p',
 'Future work includes:'),
('list',
 ['a cheap upper bound on the distance to cap band doubling: a first bound equal to the distance would save up to '
  '27% on long divergent pairs, but a bound 15% too high loses up to 23% on ultra-long reads, so only a bound '
  'within a few percent of the distance would pay;',
 'reusing computation across band-doubling rounds, as in the incremental doubling of A*PA2, once its exactness '
  'for these tiles is proven;',
 'faster construction of the seed heuristic, whose scalar lookups dominate on long divergent pairs, and a seed '
  'heuristic for gap-affine costs, at which no tool completed a 1 Mbp pair beyond 1% errors;',
 'measurements on further processors, and an evaluation of the GPU kernels the library provides for batch scores '
  'and alignments (Supplementary Section S1).']),
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
    ["Unit costs, long pairs", "Band doubling over bit-parallel tiles with the gap-chaining seed heuristic (A*PA2)", "Regrouped recurrence, staggered groups (3.2); agreeing bounds (3.3); neighbour construction (3.5); traceback within a memory limit (3.7); constant costs per pair (3.11)"],
    ["Unit costs, similar pairs", "Diagonal transition (WFA, BiWFA)", "Gathered extension (3.4); giving up on a projection with bounded waste (3.10)"],
    ["Gap-affine, two-piece", "Bidirectional wavefront (WFA, BiWFA); full sweep in difference recurrences for divergent pairs (Suzuki and Kasahara)", "Gathered extension (3.4); score-derived band (3.8); one search or two (3.9); handover with bounded waste (3.10)"],
    ["Substitution matrices, local", "Anti-diagonal vector sweep; end-first local alignment (SSW)", "Backward wavefront from the best end (3.8)"],
    ["Batches of short pairs", "Inter-sequence vectorisation (SWIPE, SeqAn, parasail)", "Certified per-lane bands (3.6)"],
  ]},
"unit": {
  "fit": True,
  "caption": "Unit-cost alignment of the real datasets by every exact aligner: mean time per alignment in ms, traceback included, and below it the peak resident memory of the run in MB, as QuickEd reports it. t.o. marks a tool that completed no pair within 30 s and o.o.m. one the 8 GB limit stopped. The simulated datasets are in Fig. 2 and Supplementary Table S5.",
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
"affine_sim": {
  "caption": "Gap-affine alignment of QuickEd's simulated datasets, random sequences with uniformly spread substitutions, insertions and deletions, at costs (4, 6, 2) by every exact aligner that takes them: mean time per alignment in ms, a 30 s budget per sample. Brackets give the pairs completed where not all were; t.o. marks a tool that completed none and o.o.m. one the 8 GB limit stopped. Datasets no tool completed are left out; every dataset is in Supplementary Table S4.",
  "header": ["Dataset", "dinara-align", "WFA", "BiWFA", "BSAlign", "KSW2"],
  "align": "lrrrrr",
  "rows": [
    ["10 kbp, 1%", "**0.199**", "0.416", "0.646", "37.0", "21.3"],
    ["10 kbp, 5%", "**2.32**", "7.14", "8.72", "37.0", "96.1"],
    ["10 kbp, 10%", "**7.06**", "23.1", "27.1", "37.3", "202"],
    ["10 kbp, 20%", "**20.1**", "66.0", "75.1", "37.4", "375 [81/99]"],
    ["100 kbp, 1%", "**11.9**", "39.9", "39.0", "o.o.m.", "2,420"],
    ["100 kbp, 5%", "**331**", "808", "714", "o.o.m.", "9,170 [3/9]"],
    ["100 kbp, 10%", "**1,160**", "o.o.m.", "2,430", "o.o.m.", "17,300 [1/10]"],
    ["100 kbp, 20%", "**3,610 [8/10]**", "o.o.m.", "7,360 [4/10]", "o.o.m.", "o.o.m."],
    ["500 kbp, 1%", "**477**", "1,290", "942", "o.o.m.", "o.o.m."],
    ["500 kbp, 5%", "**9,910 [3/4]**", "o.o.m.", "19,400 [1/4]", "o.o.m.", "o.o.m."],
    ["1 Mbp, 1%", "**2,300**", "o.o.m.", "3,940", "o.o.m.", "o.o.m."],
  ]},
"affine": {
  "caption": "Gap-affine alignment of the real datasets at costs (4, 6, 2) by every exact aligner that takes them: mean time per alignment in ms, a 30 s budget per sample. Brackets give the pairs completed where not all were; t.o. marks a tool that completed none and o.o.m. one the 8 GB limit stopped. Every dataset is in Supplementary Table S4.",
  "header": ["Dataset", "dinara-align", "WFA", "BiWFA", "BSAlign", "KSW2"],
  "align": "lrrrrr",
  "rows": [
    ["ont-1k", "**0.105**", "0.188", "0.308", "0.336", "1.10"],
    ["ont-10k", "**1.60**", "4.89", "5.86", "7.84", "30.0"],
    ["ont-50k", "**21.5**", "63.4", "54.9", "87.0", "341 [81/104]"],
    ["ont-500k", "t.o.", "o.o.m.", "t.o.", "o.o.m.", "o.o.m."],
    ["ont-500k-genvar", "**22,000 [1/4]**", "o.o.m.", "t.o.", "o.o.m.", "o.o.m."],
    ["SARS-CoV-2", "**0.710**", "2.05", "1.78", "341", "87.9"],
    ["Illumina", "**0.00145**", "0.00155", "0.00281", "0.054", "0.013"],
    ["PacBio HiFi", "**0.095**", "0.227", "0.231", "86.0", "9.72"],
  ]},
}

FIGURES = {
"overview": "How dinara-align assigns a pair to an algorithm, and the rule behind each hand-over. At unit costs, diagonal transition from one end and then from both settles near-identical pairs, giving way to band doubling over bit-parallel tiles only once its projected work exceeds the band's and it has spent one twenty-fourth of the cheapest band (Section 3.10); the tiles are traced back within the memory limit. At gap-affine costs one search grows alone, its fronts kept, and is traced back at once if it reaches the end; a second search joins once one would step more diagonals than two would need (Section 3.9), and a pair projected to cost more than a full sweep in difference recurrences is handed to that sweep after one twenty-fourth of its work, the alignment the same either way. Batches of short pairs are aligned one pair per vector lane within a band whose optimality each lane certifies (Section 3.6).",

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
 "Zook2016": ("Zook et al.", "2016"), "Hansen2025": ("Hansen et al.", "2025"), "Prjibelski2019": ("Prjibelski et al.", "2019"),
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
 "Zook2016": "Zook2016giab", "Hansen2025": "Hansen2025hg002", "Prjibelski2019": "Prjibelski2019sequence",
}
