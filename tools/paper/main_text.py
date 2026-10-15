"""The main paper, written once; render.py turns it into paper.html and dinara-align.tex.

Inline markup: **bold**, *italic*, `code`, [@Key; @Key] parenthetical citation, [@!Key] narrative citation,
\\( \\) inline maths, \\[ \\] display maths, [#tab:x] / [#fig:x] cross-references.
Blocks: ("sec", title), ("sub", title), ("p", text), ("list", [items]) numbered, ("table", key), ("figure", key).
"""

TITLE = 'dinara-align: fast, hardware-aware exact pairwise alignment that adapts its algorithm to each pair'
# Authors in order, each with an email; both are corresponding authors.
AUTHORS = [("Konstantinos Kyriakidis", "kkyriaki@ucsc.edu"), ("Benedict Paten", "bpaten@ucsc.edu")]
AFFILIATION = "UC Santa Cruz Genomics Institute, Santa Cruz, CA, USA"
KEYWORDS = ["pairwise alignment", "edit distance", "gap-affine alignment", "band doubling", "diagonal transition", "SIMD"]

ABSTRACT = [('Motivation',
  "Exact pairwise alignment returns an optimal alignment under the chosen costs, so its results never depend on a "
  "heuristic's cut-offs, and it is the reference against which heuristic aligners are judged. Yet no exact "
  'algorithm is fastest on every workload: diagonal transition wins on near-identical pairs, bit-parallel band '
  'doubling on long divergent reads, and inter-sequence vectorisation on batches of short pairs.'),
 ('Results',
  'We present dinara-align, an exact aligner that chooses the algorithm for each pair as it runs, without giving '
  'up optimality. Three rules make these choices safe. A switch suggested by a runtime estimate is taken only '
  'once the current search has spent a small fixed share of the least work the alternative could need, so '
  'estimates inflated by errors clustered at read ends cannot abandon searches that would have finished. A '
  'gap-affine alignment starts a second wavefront search only when measured work shows that two searches would '
  'be cheaper than one, a rule we prove preserves optimality. Short pairs in batches are aligned in narrow bands '
  'whose optimality each vector lane certifies. On one Intel Skylake-X core and 25 datasets, from 250 bp reads to '
  '1 Mbp pairs, dinara-align was the fastest exact aligner at unit costs, up to 2.4 times faster than '
  'A*PA2 on ultra-long nanopore reads and 2.5 and 5.7 times faster than WFA on Illumina and HiFi '
  'reads, and at gap-affine costs the fastest on every dataset any tool completed, 1.3 times faster than '
  'WFA on Illumina reads.'),
 ('Availability and implementation',
  'dinara-align is free software under the Mozilla Public License 2.0, available at '
  'https://github.com/kokyriakidis/dinara-align with C, C++, Python and command-line interfaces.'),
]

BODY = [
('sec', 'Introduction'),
('p',
 'Pairwise alignment finds the least costly series of substitutions, insertions and deletions that transforms '
 'one sequence into another. It underlies read mapping, variant detection, genome assembly and multiple sequence '
 'alignment [@Prjibelski2019]: read mappers align each read to the reference window its seeds select, assemblers '
 'compare reads with one another to find their overlaps [@Li2018], and assemblies, haplotypes and viral genomes '
 'are compared base by base.'),
('p',
 'Read lengths have grown from hundreds of bases to hundreds of kilobases, and the gap-affine dynamic program '
 '[@Gotoh1982] fills an \\(n \\times m\\) matrix, \\(2.5 \\times 10^{11}\\) cells for one pair of 500 kbp reads. '
 'Heuristic aligners therefore compute only part of it: a fixed or adaptive band [@Liu2023], a series of '
 "windows, or the cells before the score drops by more than a threshold, as in X-drop or minimap2's Z-drop "
 '[@Li2018]. Their result then depends on these cut-offs, and alignments with large indels or clustered errors '
 'can be missed or misplaced. Exact aligners return an optimal alignment under the chosen costs whatever the '
 'data, which also makes them the ground truth against which heuristics are evaluated, and the fastest of them '
 'skip only cells that provably cannot lie on an optimal path, so that their running time grows with the '
 'distance as well as the length rather than with the product of the lengths.'),
('p',
 'Three ideas underlie current exact aligners. *Band doubling* computes only the diagonals a path within a cost '
 'bound can reach, doubling the bound until the optimum lies within it [@Ukkonen1985]; *bit-parallelism* encodes '
 'the score differences of 64 cells of a column in two machine words [@Myers1999], and Edlib combines the two '
 '[@Sosic2017]. *Diagonal transition* extends the furthest-reaching paths one cost at a time [@Ukkonen1985; '
 '@Myers1986]; WFA extended it to gap-affine costs [@MarcoSola2021] and BiWFA made its memory linear in the '
 'distance [@MarcoSola2023]. A*PA guided an A* search by a gap-chaining seed heuristic [@GrootKoerkamp2024a]; '
 'A*PA2 combined that heuristic with band doubling over SIMD blocks, aligning ultra-long reads up to 19 times '
 'faster than Edlib and BiWFA [@GrootKoerkamp2024b]; and QuickEd aligns within a bound it finds heuristically '
 '[@Doblas2025]. Full-matrix vectorisation, in anti-diagonal [@Wozniak1997], striped [@Farrar2007; @Zhao2013], '
 'prefix-scan [@Daily2016], difference-recurrence [@Suzuki2018; @Li2018] and inter-sequence [@Rognes2011; '
 '@Rahn2018] layouts, takes time in proportion to the matrix or the band instead.'),
('p',
 'Each of these methods is fastest on part of the range of workloads: diagonal transition, whose expected work '
 'grows with the square of the distance, on short and near-identical pairs; bit-parallel band doubling, whose '
 'work grows with the length times the distance, on long divergent reads, with the seed heuristic on ultra-long '
 'ones; and inter-sequence vectorisation on batches of short pairs. An aligner that is fastest across workloads '
 'must therefore choose among them for each pair, and an exact aligner must make every such choice (which '
 'algorithm, which bound, when to abandon a search) without giving up optimality. In dinara-align every engine '
 'is exact, so no choice changes the optimal cost, and a choice that rests on a runtime estimate can delay the '
 "work by at most a fixed small share of the alternative's. Our contributions are:"),
('list',
 ['**One wavefront search or two, chosen by measured work** (Section 3.9): a gap-affine alignment adds a second '
  'search only when its measured work exceeds what two would need, and still returns an optimal alignment '
  '(Proposition 2).',
 "**Runtime estimates that cannot derail a search** (Sections 3.3 and 3.10): a projected distance sets a band's "
  'first bound only when two projections agree, and a search hands over to an algorithm projected to be cheaper '
  "only after spending one twenty-fourth of that algorithm's least possible work.",
 '**Certified bands for inter-sequence batches** (Section 3.6): each pair is aligned in a vector lane within a '
  'narrow band that the lane certifies optimal under gap-affine and two-piece costs (Proposition 1), bringing the '
  'speculate-and-verify approach of narrow-band aligners [@Gibrat2018; @Fujiki2020] to every lane of a batch.',
 '**Direct neighbour construction for the seed heuristic** (Section 3.5): inexact neighbours of exact seed '
  'matches are inserted without search, provably including every one that local pruning keeps (Lemma 1); we also '
  "prove that A*PA2's pruning keeps the heuristic a lower bound (Lemma S1).",
 '**One alignment whichever engine computes it** (Sections 3.1 and 3.7): at unit, gap-affine and two-piece costs '
  'every engine breaks ties by the same rule, and a unit-cost band recomputes stretches of its traceback rather '
  'than splitting the pair, so its alignment does not depend on the memory limit; only a gap-affine pair split to '
  'fit that limit follows the rule within each piece.',
 '**Hardware-aware engines** (Sections 3.2, 3.4 and 3.11): diagonal transition with gathered AVX-512 loads, '
  'interleaved bit-parallel groups with a shorter dependency chain, and lower per-pair overhead on short reads.',
 '**One library** ([#tab:engines]; Supplementary Table S3) for global, ends-free, extension and local alignment '
  'under unit, gap-affine, two-piece and substitution-matrix costs, for single pairs and batches, through C, C++, '
  'Python and command-line interfaces.']),
('p',
 'We evaluate dinara-align on the 25 datasets of the A*PA2 and QuickEd benchmarks, nanopore, Illumina, PacBio '
 "HiFi and PromethION reads, SARS-CoV-2 genomes and simulated pairs of up to 1 Mbp, with QuickEd's Illumina, "
 "HiFi and PromethION sets built from HG002 reads, each paired with the stretch of HG002's own assembly it "
 'covers. It was the fastest exact aligner on all of them at unit costs, and at gap-affine costs on every '
 'dataset any tool completed.'),
('sec', 'Background'),
('sub', 'Alignment and costs'),
('p',
 'Let \\(A=a_1\\dots a_n\\) be the reference and \\(B=b_1\\dots b_m\\) the query, and let cell \\((i,j)\\) lie on '
 'diagonal \\(j-i\\). Under gap-affine costs \\((x, o, e)\\), for example (4, 6, 2), a match costs 0, a mismatch '
 '\\(x\\) and a gap of \\(k\\) letters \\(g(k)=o+ke\\); unit costs have \\(x=e=1\\) and \\(o=0\\), two-piece costs take '
 'the minimum of two such functions [@Li2018], and a substitution matrix instead scores each pair of letters. '
 'Global alignment compares two whole sequences; ends-free alignment places a read inside a longer window or '
 'finds the overlap of two reads; extension grows an alignment from a seed match [@Li2018]; and local alignment '
 'finds the best-matching parts of two sequences. Extension, local and, where asked, ends-free alignment instead '
 'maximise a score that rewards matches (Section 3.8). An alignment is *optimal* if no alignment has a lower '
 'cost, or a higher score, under the requested costs and mode; a caller may further restrict alignments to a '
 'band of diagonals or a maximum cost, and an exact aligner then returns an optimal alignment among those or '
 'reports that none exists. At unit costs the distance \\(D(i,j)\\) between the prefixes \\(a_1\\dots a_i\\) and '
 '\\(b_1\\dots b_j\\) satisfies'),
('p',
 '\\[ D(i,j) = \\min\\bigl(D(i-1,j)+1,\\; D(i,j-1)+1, \\qquad D(i-1,j-1)+[a_i \\ne b_j]\\bigr), \\]'),
('p',
 'where \\([a_i \\ne b_j]\\) is 1 if the letters differ and 0 otherwise, with \\(D(i,0)=i\\) and \\(D(0,j)=j\\), and '
 'the distance is \\(D(n,m)\\); gap-affine costs need three such layers [@Gotoh1982].'),
('sub', 'Band doubling and bit-parallel computation'),
('p',
 'At unit costs, an alignment of cost at most \\(t\\) visits only diagonals \\(d\\) with \\(|d| + |d - \\delta| \\le '
 't\\), where \\(\\delta = m - n\\); band doubling computes only those cells, doubling \\(t\\) until \\(D(n,m) \\le t\\) '
 "[@Ukkonen1985]. Myers' recurrence advances the score differences of 64 cells of a column, held in two machine "
 'words, by one column at a time [@Myers1999], and A*PA2 computes the band in tiles of 256 columns with SIMD '
 'vectors [@GrootKoerkamp2024b].'),
('sub', 'Diagonal transition'),
('p',
 'Diagonal transition keeps, for each cost \\(s\\) and diagonal, the furthest cell an alignment of cost \\(s\\) '
 'reaches, extended along matches at no cost, in \\(O(ns)\\) time in the worst case and \\(O(n + s^2)\\) expected '
 'time on random sequences [@Ukkonen1985; @Myers1986]; WFA generalises it to gap-affine costs, where the '
 'wavefront, or front, of each cost has three layers [@MarcoSola2021], and BiWFA runs it from both ends '
 '[@MarcoSola2023].'),
('sub', 'The seed heuristic'),
('p',
 'A*PA divides the reference into disjoint seeds of fixed length and finds their matches in the query, exactly '
 'or with one edit [@GrootKoerkamp2024a]. The potential of a seed is the cost \\(r\\) of crossing it unmatched, 1 '
 'for exact seeds and 2 for seeds matched with one edit, and a match crosses a seed for less than \\(r\\). The '
 'gap-chaining seed heuristic charges a chain of matches from a cell the cost of each match and, for each gap '
 'between consecutive matches and from the last to the end, the larger of the potential of the seeds that start '
 'in the gap and the number of diagonals it shifts. Its value, the least charge of any chain, never exceeds the '
 'remaining cost, and A*PA keeps the match starts in layers by how much potential the best chain from each saves '
 "[@GrootKoerkamp2024a]. A*PA2 prunes matches before the alignment, in place of A*PA's pruning during the "
 'search: a match \\(\\nu\\) with cost \\(c(\\nu)\\) and end cell \\(\\mathrm{end}(\\nu)\\) is kept only if a path from '
 '\\(\\mathrm{end}(\\nu)\\), starting at cost \\(c(\\nu)\\), crosses the next \\(p = 14\\) seeds, counting that of '
 '\\(\\nu\\), while its cost stays below the potential of the seeds it has crossed [@GrootKoerkamp2024b]. A '
 'diagonal transition from \\(\\mathrm{end}(\\nu)\\) tests this, dropping each furthest-reaching point whose cost '
 'reaches that potential; a further point on a diagonal has crossed at least as much potential, so no such path '
 'is missed. A*PA2 uses the heuristic to drop the rows of each band whose score plus the heuristic exceeds the '
 'bound [@GrootKoerkamp2024b].'),
('sec', 'Methods'),
('sub', 'Overview'),
('p',
 'Each pair is assigned to an algorithm ([#tab:engines]; [#fig:overview]). Unit costs follow A*PA2 '
 '[@GrootKoerkamp2024b]: band doubling over 256-column tiles of bit-parallel words, pruned on long pairs by the '
 'gap-chaining seed heuristic [@GrootKoerkamp2024a], with near-identical pairs aligned by diagonal transition '
 'alone. Gap-affine costs use the wavefront of WFA and BiWFA [@MarcoSola2021; @MarcoSola2023]. At unit costs, '
 "diagonal transition runs first and hands over once its projected remaining work exceeds the band's and it has "
 "spent a twenty-fourth of the band's least work (Section 3.10); its projection, the length and the cost model "
 'then choose the next algorithm and its parameters.'),
('p',
 'At unit, gap-affine and two-piece costs, every engine selects the same alignment among equally good ones, by '
 'the backtrace order of WFA2-lib: from the end, an edit is preferred to a match whenever one is optimal, a '
 'substitution before a deletion (a letter of \\(A\\) alone) before an insertion, the second gap piece before the '
 "first, and an extension before an opening. This gives WFA2-lib's CIGARs exactly; the default applies it to both "
 'sequences reversed, which places every edit as far left as possible, as variant callers normalise indels. Every '
 'engine takes the furthest-reaching source or stores per-cell flags under the rule. A gap-affine pair split for '
 'the memory limit follows the rule within each piece, and alignments under a substitution matrix follow the '
 "traceback order of Gotoh's recurrence."),
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
 "where \\(E' \\lor P_v\\) and \\(P_v \\land E'\\) do not depend on \\(c\\) and are computed off the critical path. The "
 'vertical update is regrouped in the same way, its negated part computed from \\(E\\), \\(M_v\\) and the incoming '
 'horizontal difference before \\(P_h\\) is known. Together the rewrites shorten the loop-carried chain of a '
 'single word from 14 to 8 cycles per column. [@!Myers1999], Edlib and A*PA2 use the original form, and we are '
 'not aware of a previous description of the regrouped one.'),
('p',
 'Words occupy the lanes of a vector, staggered by one column, so that the horizontal difference leaving one word '
 'enters the next by a lane rotation [@Wozniak1997; @GrootKoerkamp2024b], on which each column waits. On AVX-512 '
 'dinara-align interleaves two eight-word groups, ten columns apart, so that each fills the latency of the other, '
 'and bands of two or three words run in general-purpose registers.'),
('sub', 'Band bounds from agreeing estimates'),
('p',
 'Band doubling is exact for any initial bound, but its running time depends on how close that bound is to the '
 'distance. A computation that has reached cost \\(s\\) at anti-diagonal \\(\\rho\\) estimates the distance as \\(\\hat '
 'D = s(n+m)/\\rho\\), which is reliable when edits are spread evenly and too high by orders of magnitude when they '
 'gather, as at the ends of nanopore reads. dinara-align takes two estimates, \\(\\hat D_1\\) from a one-ended '
 'diagonal transition stopped after eight edits and \\(\\hat D_2\\) from the two-ended search where it handed over. '
 'When at least one seed in 32 chains from the origin, the first bound starts 64 above the seed heuristic there. '
 'Otherwise it is aimed just above \\(\\hat D_1\\) when the two estimates agree within a third, and failing that, as '
 'in A*PA2, it starts at most 256 above the known lower bound and at most doubles until a failed round has '
 'crossed an eighth of the columns. A failed round aims the next bound from where it stopped, and the first round '
 'estimates again after an eighth, a quarter and a half of the columns. A distance is accepted only within the '
 'bound its round ended on, so these rules change only the sequence of bounds and the result stays exact.'),
('p',
 'Within a round, the band drops the rows at the top and bottom of each tile whose score plus the heuristic '
 'exceeds the bound, as A*PA2 does. A score changes by at most one per row and the heuristic by at most its climb '
 '\\(\\gamma\\), so a row that exceeds the bound by \\(y\\) also rules out the next \\(\\lceil y/(1+\\gamma) \\rceil - 1\\) '
 'rows, which are skipped unread.'),
('sub', 'Diagonal transition for near-identical pairs'),
('p',
 'Diagonal transition extends each diagonal along its run of matches, eight letters at a time by the exclusive-or '
 'of two 64-bit words, and on AVX-512 eight diagonals together, their letters gathered under branch-free masks '
 '(Supplementary Section S1.1). The same step serves near-identical pairs, from one or both ends '
 '[@MarcoSola2023], the traceback inside each tile, and the gap-affine front. A pair finished from both ends has '
 'its forward wavefronts grown on to the proven distance, kept only on the diagonals an optimal path passes, and '
 'traced from them (Section 3.1).'),
('sub', 'Seed heuristic: direct neighbour construction'),
('p',
 'The reference is divided into disjoint seeds of 12 bases, matched exactly, or of 16 bases, matched with at most '
 'one edit, each of whose 8-base halves indexes a table [@GrootKoerkamp2024a]; the lookups for 16 query rows are '
 'issued together, so that their cache misses overlap.'),
('p',
 'An exact match \\(\\mu\\) of a 16-base seed at query rows \\([j, j+16)\\) has four neighbouring matches with one '
 'edit: two share its start, at rows \\([j, j+15)\\) and \\([j, j+17)\\), and two share its end, at \\([j+1, j+16)\\) '
 'and \\([j-1, j+16)\\). Chaining is all-or-nothing, so a neighbour may chain where \\(\\mu\\) cannot, and the '
 'neighbours are needed for the heuristic to remain a lower bound. Rather than finding each by search and testing '
 'it for local pruning, dinara-align inserts them whenever \\(\\mu\\) is kept, as the following lemma justifies; the '
 'search is that of local pruning (Section 2.4).'),
('p',
 '**Lemma 1.** If the search from an exact match \\(\\mu\\) finds no path across the next \\(p\\) seeds, it finds '
 'none from any of its four neighbours.'),
('p',
 '*Proof idea.* The search from \\(\\mu\\) dominates that from each neighbour: at every cost and on every diagonal '
 'it reaches at least as far, and a further point has crossed at least as much potential, so it is dropped only '
 'if the nearer one is too (proof in Supplementary Section S6).'),
('p',
 'Every neighbour that the search would keep is therefore inserted. As in A*PA2, a match whose search runs into a '
 'match already kept on its diagonal is kept too. A match is dropped only if its search would find no path, which '
 'keeps the heuristic a lower bound (Lemma S1), and each further match kept can only lower it, so the alignment '
 'stays exact.'),
('sub', 'Certified bands for inter-sequence batches'),
('p',
 'Batches of short pairs are aligned with one pair per vector lane [@Rognes2011; @Rahn2018] in saturating 8-bit '
 'lanes [@Zhao2013; @Daily2016], a saturated score computed again in 16-bit lanes. Rather than the full matrix, '
 'the lanes of a group compute the band of diagonals \\([d_\\ell, d_h]\\) containing the main diagonal 0 and each '
 "pair's end diagonal \\(\\delta\\), with one further diagonal on either side; pairs are grouped by end diagonal, so "
 'the band stays narrow. Its optimality is certified by a cost form of the out-of-band bounds of [@!Gibrat2018] '
 'and [@!Fujiki2020], with a strict form for ties.'),
('p',
 '**Proposition 1.** Let \\(C\\) be the least cost of an alignment within \\([d_\\ell, d_h]\\), and let \\(L(d) = '
 'g(|d|) + g(|d-\\delta|)\\), or \\(\\infty\\) for a diagonal outside the matrix. If \\(C \\le \\min\\bigl(L(d_h+1), '
 'L(d_\\ell-1)\\bigr)\\), then \\(C\\) is the optimal cost; if the inequality is strict, every optimal alignment lies '
 'within the band.'),
('p',
 '*Proof idea.* An alignment that leaves the band must open a gap reaching the diagonal beyond it and another '
 'returning to the end diagonal, which together cost at least \\(L(d)\\) for that diagonal \\(d\\) (proof in '
 'Supplementary Section S6).'),
('p',
 'Each lane evaluates its own certificate. When every optimal alignment lies in the band, its in-band scores are '
 'exact and the tie rule chooses as in the full matrix; pairs that fail are recomputed over every diagonal that a '
 'path at most as cheap could visit, grouped by the width of that band.'),
('sub', 'Traceback within a memory limit'),
('p',
 'The traceback of a band traces each tile from the boundary recorded on its left. Boundaries are kept within the '
 "alignment's memory limit, 80 MiB by default: beyond it only one in two is kept, then one in four, and so on, "
 'and each tile records the state the band carries into it, from which the traceback recomputes a stretch of '
 "tiles from the last kept boundary (Supplementary Section S1.1). A band's steps depend only on that state, so "
 'each boundary is recovered as the round left it and the alignment is the same byte for byte; splitting the pair '
 'instead, as QuickEd does [@Doblas2025], resolves ties within each part.'),
('sub', 'Gap-affine and score-based alignment'),
('p',
 'A scoring scheme with match score \\(\\alpha\\), mismatch score \\(\\beta\\) and gap scores \\(-(o+ke)\\) can be '
 'converted into costs for global alignment [@Eizenga2022]. Every alignment covers \\(n+m\\) letters, so its score '
 '\\(S\\) and a gap-affine cost \\(C\\) satisfy'),
('p',
 '\\[ S = \\tfrac12\\bigl(\\alpha(n+m) - C\\bigr), \\qquad C = 2(\\alpha-\\beta)\\,X + \\sum_{\\text{gaps}} \\bigl(2o + '
 'k(2e+\\alpha)\\bigr), \\]'),
('p',
 'where \\(X\\) is the number of mismatches, so the wavefront computes score-maximising alignments. Local '
 'alignment locates the end of the optimal alignment with a forward sweep, as in SSW [@Zhao2013], and extends a '
 'wavefront back from it to the first cell attaining the maximal score, which avoids the correction pass that '
 'striped layouts may need [@Farrar2007].'),
('p',
 'The wavefront keeps five bytes per diagonal for the traceback, and a pair whose fronts would exceed the memory '
 'limit is split where an optimal path crosses, as in BiWFA (Section 3.1). A global pair projected to take more '
 'work than a full sweep in 8-bit difference recurrences [@Suzuki2018] is handed to that sweep (Section 3.10), '
 'which keeps the predecessor flags of Section 3.1.'),
('p',
 'Under a substitution matrix, anti-diagonals [@Wozniak1997] are computed in 16-bit lanes while the scores '
 'provably fit, in strips of rows that stay in the first-level cache, and a global alignment is traced in a band '
 'of diagonals that grows until its own score proves that every optimal path lies inside it (Supplementary '
 'Section S1.4).'),
('sub', 'One search or two at gap-affine costs'),
('p',
 'WFA grows the wavefront from one end of the pair and keeps every front for the traceback [@MarcoSola2021]; '
 'BiWFA grows it from both ends and keeps only the latest fronts [@MarcoSola2023]. Two searches must also prove '
 'where an optimal path crosses between them, so they stop only once \\(2\\min(s_1, s_2) \\ge C^* - 1 + \\theta\\), '
 'where \\(s_1, s_2\\) are their costs, \\(C^*\\) the cheapest meeting found and \\(\\theta = o + w\\) a margin: \\(o\\) '
 'is the largest gap-open cost and \\(w = \\max(x, o+e)\\) the largest cost of a single edit, over both pieces for '
 'two-piece costs. Every optimal path has a cell whose costs from the two ends differ by at most \\(w\\), and a '
 'meeting inside a gap pays its opening once for both halves; BiWFA stops by the same margin. A score-only '
 'computation, which needs no meeting cell, stops sooner: once its best is at most one more than the highest cost '
 "to which either search has checked the other's start. With \\(W(s)\\) the number of diagonals a search steps to "
 'reach cost \\(s\\), two searches need about \\(2W((C_{\\mathrm{opt}}+\\theta)/2)\\) steps for a pair of cost '
 '\\(C_{\\mathrm{opt}}\\), plus the setup of the second, and one search \\(W(C_{\\mathrm{opt}})\\). Where the fronts '
 'widen with the cost, \\(W(s)\\) grows with \\(s^2\\) and one search is cheaper while \\(C_{\\mathrm{opt}} < '
 '(1+\\sqrt2)\\,\\theta\\); where the ends of the sequences bound the fronts, it is always cheaper. Most 250 bp '
 'Illumina reads cost less than \\(\\theta\\), 14 at (4, 6, 2).'),
('p',
 'dinara-align therefore starts with a single search, the one whose fronts the traceback reads (the backward '
 'search under the default tie rule), keeps its fronts and traces the alignment as soon as the search reaches '
 'the far corner. Every cost \\(s\\) passed without reaching it shows that \\(C_{\\mathrm{opt}} > s\\), and the '
 'search measures \\(W\\) as it grows: it continues alone while \\(s < \\theta\\), costs the second search would pass '
 'in any case, and then while \\(W(s) \\le 2W(\\lfloor (s+\\theta)/2 \\rfloor)\\), that is, while it has stepped no '
 'more diagonals than two searches would need for the cheapest cost still possible. It also hands over at half '
 'the memory limit or when a band leaves it no reachable cell. The rule has no tuned constant: it hands over at '
 'about \\(2.4\\theta\\) where the fronts widen with the cost, and never where the sequences bound the fronts from '
 'early costs on, where it aligns exactly as WFA does.'),
('p',
 '**Proposition 2.** Suppose the first search stops at cost \\(\\sigma\\) without reaching the far corner, and the '
 'second search, growing from cost 0, checks each new front against the last \\(w+1\\) fronts of the first, costs '
 '\\(\\sigma-w\\) to \\(\\sigma\\), while its own cost is at most \\(\\sigma\\), after which the two alternate as in '
 "BiWFA, each new front checked against the other side's last \\(w+1\\). When they stop, \\(C^*\\) is the optimal "
 'cost and its meeting lies on an optimal path.'),
('p',
 '*Proof idea.* On every optimal path the costs from the two starts change by at most \\(w\\) a move, so some cell '
 'has them within \\(w\\) of each other, and the two phases check every such cell before the stopping rule can end '
 'the searches (proof in Supplementary Section S6).'),
('p',
 'An alignment traced from the unpruned fronts of the single search is the one the two searches would have '
 'traced: the traceback steps only from a cell to a predecessor on an optimal path, and the fronts that two '
 'searches prune hold no such cell, so the tie rule is unchanged (Supplementary Table S1 gives the effect on '
 'short reads).'),
('sub', 'Hand-overs that wait on a runtime estimate'),
('p',
 'Two hand-overs rest on projections of the work still needed: diagonal transition hands over to band doubling '
 "once its remaining work is projected to exceed a band's, and the gap-affine wavefront to the difference sweep "
 "once its projected total exceeds the sweep's. A projection can exceed the distance a hundredfold when edits "
 'gather at the ends, as on nanopore reads with adapters or poly-A tails, and hand a pair that the search would '
 'finish in microseconds to a band over the whole matrix. dinara-align hands over on a projection only once the '
 'work already spent has reached one twenty-fourth of the cheapest alternative that could still be needed: for '
 'diagonal transition, the band whose bound equals the cost already reached, which no band can undercut; for the '
 'wavefront, the full sweep. Once the projection favours the alternative, this delays the hand-over by at most '
 "about 4% of the work the alternative needs in any case. An alignment's one-ended diagonal transition also waits "
 "until the projected work exceeds the band's by a factor \\(1 + 1.5/\\sqrt{q}\\) after \\(q\\) edits, 1.5 times the "
 'relative error of a projection from \\(q\\) edits. A pair of at most 2,048 columns hands over on no projection: '
 'its band covers the whole matrix at a cost known in advance, so its searches give up only once their steps '
 'would have paid for that band, and it costs at most twice the cheaper of the two (Supplementary Table S1).'),
('sub', 'Per-pair overhead on short reads'),
('p',
 'A 250 bp read aligns in well under a microsecond, so work that does not depend on the pair dominates. An '
 'aligner object keeps every buffer its searches use from pair to pair, one aligner per thread, as each '
 "competitor's benchmark runner keeps its own; both sequences are encoded, and reversed for the leftmost tie "
 'rule, in one vectorised pass; and an aligner warms its buffers on two built-in pairs when it is made, about 80 '
 'µs once, while a single call takes one that does not (Supplementary Table S1).'),
('sub', 'Thresholds measured against the code'),
('p',
 'The choice of algorithm and the bounds of each search rest on 41 thresholds. We built dinara-align with each '
 'halved and with each doubled and timed the 82 builds against the unchanged one: about thirty moved no sample by '
 'more than 3%, most others traded one dataset against another, and one had become outdated, the length from '
 'which seeds match with one edit (Supplementary Section S1.1).'),
('sub', 'Validation'),
('p',
 'A differential fuzzer checks every mode, cost model, band, cost limit and memory limit on random pairs against '
 'a full-matrix implementation that shares no code with the library, for validity, optimality and the tie rule. A '
 "unit-cost alignment is the same byte for byte under any memory limit, and all 305 pairs of WFA2-lib's "
 'regression set give its scores in every mode it tests and its CIGARs under its tie rule (Supplementary Section '
 'S1.7); at the benchmarked commit 20,000 random cases and the regression set passed.'),
('sec', 'Results'),
('sub', 'Experimental setup'),
('p',
 'All measurements were taken on an Intel Core i9-7900X (Skylake-X, AVX-512), its ten cores at a fixed 3.3 GHz '
 'with turbo boost and hyper-threading disabled, with dinara-align compiled by Mojo 1.1.0 [@Modular] at commit '
 "`87d365f` and every tool compiled for the native instruction set except BSAlign's gap-affine kernel, which its "
 'own build keeps to SSE4.2. We combined the protocols of A*PA2 [@GrootKoerkamp2024b] and QuickEd [@Doblas2025].'),
('p',
 "The real datasets come from both. From A*PA2's release [@PaBench] we use the nanopore sets ont-1k, ont-10k, "
 'ont-50k and ont-500k, named by length class (mean 0.83, 3.6, 9.6 and 642 kbp), CHM13 reads paired with the '
 'CHM13 reference; ont-500k-genvar, ultra-long reads of another individual against the same reference, which '
 'adds genetic variation and gaps of up to 42 kbp; and 30 kbp SARS-CoV-2 genomes. We call the two 500k sets '
 "ultra-long. For QuickEd's Illumina 2×250 bp, PacBio HiFi and nanopore PromethION technologies we use HG002 "
 "reads from Genome in a Bottle [@Zook2016]: Illumina and HiFi reads of one region of chromosome 20, and UCSC's "
 "ultra-long PromethION reads [@Shafin2020] of one region of chromosome 12, of 10 kbp or more as in QuickEd's "
 'set (ont-promethion; 45 kbp on average in the sample, up to 156 kbp). Each read was picked with mapping '
 "quality at least 20 on GRCh38, mapped by minimap2 [@Li2018] to both haplotypes of HG002's own assembly "
 '[@Hansen2025] and paired with the stretch it covers, its soft-clipped ends trimmed, so that pairs differ '
 "mainly by sequencing errors. The simulated datasets are QuickEd's random sequences [@QuickEdData] of 10 kbp to "
 '1 Mbp whose pairs differ by single-base substitutions, insertions and deletions spread uniformly at nominal '
 'error rates of 1, 5, 10 and 20%, modelling no sequencing technology; the 20% sets measure 17% divergence '
 '(Supplementary Table S6).'),
('p',
 'Each dataset was sampled to about 2 Mbp of sequence, both sides counted, and at least four pairs. Every tool '
 'aligned each pair with traceback, on one pinned core, through one aligner object per sample as its own '
 'benchmark runner does, with a 30 s budget per sample, in three rounds; a '
 "pair's time is its median over the rounds, and a tool is credited only with the pairs it completed in every "
 'round. The other tools were timed in rounds that alternated between them, and dinara-align at its final version '
 'a day later, on the same machine, builds and samples, every answer checked against theirs again. The rounds '
 "agreed closely: the median spread between a tool's fastest and slowest round was 0.5% at unit and 0.6% at "
 "gap-affine costs, at most 24%, on dinara-align's batch alignments of the HiFi reads. Timings quoted "
 'in Section 3 were recorded when each change was made, against the version before it (Supplementary Table S1). '
 'The exact aligners compared are A*PA2-full (with the seed heuristic and pruning) and A*PA2-simple (band '
 'doubling alone), A*PA [@GrootKoerkamp2024a], Edlib [@Sosic2017], WFA and BiWFA, which denote WFA2-lib in its '
 'high-memory and bidirectional modes [@MarcoSola2021; @MarcoSola2023], QuickEd, BSAlign [@Shao2024], KSW2 '
 '[@Li2018], parasail [@Daily2016] and SSW [@Zhao2013], each on the workloads it supports. On the pairs any two '
 'tools completed, their costs agreed, and at unit costs they agreed with the published A*PA2 results on its '
 'datasets. The other alignment modes were timed in-process on simulated workloads, and the batches on 500,000 '
 'simulated pairs of 148 bp at 1% error, each the median of three runs (Supplementary Tables S3 and S8).'),
('sub', 'Speed and memory'),
('p',
 'At unit costs dinara-align had the lowest mean time per alignment of the exact aligners on all 25 datasets. On '
 'the real data ([#tab:unit]) it was 1.5 to 2.4 times faster than the fastest alternative on the nanopore reads '
 '(1.6 times on PromethION) and 3.1 times on the SARS-CoV-2 genomes, and 2.5 and 5.7 times faster than WFA, the '
 'fastest of the others there, on the Illumina and HiFi reads. On the simulated data ([#fig:sweeps]; '
 'Supplementary Table S5) it was 1.4 to 2.9 times faster than the fastest alternative: 100 kbp pairs at 1, 5, 10 '
 'and 20% errors took 1.1, 4.9, 7.0 and 15.4 ms, against 7.5, 8.4, 13.9 and 36.8 ms for A*PA2-full and 12.5, '
 '23.9, 39.0 and 62.6 ms for QuickEd. Its narrowest margin was on 10 kbp pairs at 10%, 0.315 against 0.444 ms for '
 'A*PA2-simple. Over all 25 datasets QuickEd took 3.5 to 82 times as long, 8.7 times in the geometric mean.'),
('p',
 'At gap-affine costs dinara-align was 1.3 to 3.3 times faster than the faster of WFA and BiWFA on the 16 '
 'datasets all three completed ([#tab:affine]; Supplementary Table S4), the lowest on the Illumina reads, 1.16 '
 'against 1.54 µs for WFA. It aligned every PromethION pair, in 261 ms '
 'against 576 ms for BiWFA, where WFA ran out of memory after 7 of 22, and completed more of the simulated pairs '
 'of 100 kbp or more than any other tool, 47 against 41 for BiWFA; but it completed no ont-500k read and one of '
 'the four ont-500k-genvar reads. In the other modes it was as fast as, or up to 22 times faster than, the '
 'fastest tool offering the same mode, infix alignment at (4, 6, 2) being a tie (Supplementary Table S3). Global '
 'scores of 500,000 short-read pairs took 229 ms on one core and 30 ms on ten at (1, 2, 1), 4.6 and 3.7 times '
 'faster than WFA2-lib (Supplementary Table S8).'),
('table', 'unit'),
('table', 'affine'),
('table', 'affine_sim'),
('figure', 'sweeps'),
('p',
 'An alignment of an ultra-long read increased peak memory by at most 54 MiB, against 81 MiB for A*PA2-full, and '
 'of a 1 Mbp pair at 20% errors by at most 131 MiB, against 184 MiB (Supplementary Table S7). This growth '
 'includes the seed heuristic and the band; the stored tile boundaries stayed within the default limit of 80 MiB '
 "(Section 3.7). QuickEd stayed near 25 MiB. On short reads dinara-align's peak was higher than the others' "
 '([#tab:unit]): the Mojo runtime and the sample, read whole, take about 13 MiB.'),
('sec', 'Discussion'),
('p',
 'dinara-align was the fastest exact aligner compared on every dataset at unit costs, and at gap-affine costs on '
 'every dataset that any tool completed. This is consistent with choosing '
 'an algorithm for each pair: no single engine was fastest across these workloads, and since every engine is '
 'exact, choosing never forfeits the guarantee or changes the alignment. Three observations stand out. First, '
 'each technique paid where its engine does the work (Supplementary Table S2): on the ultra-long reads, '
 'switching off direct neighbour construction, which reduces work outside the inner loop, slowed alignment by 9 '
 'and 16%, more than either change to the loop itself (at most 5%), while the interleaved groups paid most on '
 'the widest bands and gathered extension where diagonal transition does most of the work. Second, on short '
 'reads the work that does not depend on the pair dominates: a 250 bp alignment at gap-affine costs once took '
 '5.2 µs through a fresh call and 2.4 µs through a reused aligner (Supplementary Table S1), so a benchmark '
 'should call every tool as its own runner does. Third, thresholds that encode a cost model become outdated as '
 'the code changes, and measuring them again belongs in the evaluation of any aligner that selects its algorithm '
 'by thresholds. Gap-affine alignment of ultra-long reads remains open: no tool completed an ont-500k read '
 'within 30 s. Where memory matters more than time, Edlib, QuickEd and BiWFA (15 to 29 MiB on ultra-long reads) '
 "use less than dinara-align's 71 to 81 MiB. The measured-work rule, the regrouped recurrence and the per-lane "
 'certificate do not depend on this implementation and apply to other wavefront, bit-parallel and inter-sequence '
 'aligners. And because its alignments are exact and the same whichever engine computes them, dinara-align can '
 'serve as the reference against which heuristic aligners are measured, at the scale of current long reads.'),
('sub', 'Limitations'),
('list',
 ['The evaluation covers one x86 processor with AVX-512, and single-thread latency except for the batches; the '
  'AVX2 and NEON configurations are not evaluated here.',
 'The evaluated samples also guided development: several changes were motivated and checked on the Illumina, '
  'HiFi and ont-500k-genvar samples, and some thresholds, in particular the AVX-512 gates for inexact seeds, were '
  'set on uniform 100 kbp pairs like those of [#fig:sweeps]a, so results on unseen data may be less favourable.',
 'The ultra-long and the largest simulated sets were sampled to four pairs. The three rounds were timed from one '
  'build, while on the Skylake-X running times varied by up to a third between builds that left the code involved '
  'unchanged (Supplementary Table S2), so margins of less than about 10% may not hold across builds. dinara-align '
  "was timed a day after the other tools, on the same machine; replayed then, the other tools' times stayed within "
  '3% of their own rounds.',
 'The Illumina, HiFi and PromethION reads come from one region of one chromosome each and exclude reads that '
  'minimap2 leaves more than a tenth clipped, such as those spanning structural variants.',
 'The ablation switches off one technique at a time, so it does not separate their interactions.']),
('sub', 'Future work'),
('list',
 ['A cheap upper bound on the distance. With the true distance supplied by an oracle, starting band doubling at '
  'it saves up to 27% on 100 kbp pairs at 20% errors, but starting 15% above it loses up to 23% on ultra-long '
  'reads, whereas capping the bound at an upper bound within 5% of the distance saves up to 27% on 1 Mbp pairs at '
  '20% (Supplementary Section S5).',
 'Reusing computation across band-doubling rounds, as in the incremental doubling of A*PA2, once its exactness '
  'for these tiles is proven.',
 'Faster construction of the seed heuristic, whose scalar lookups dominate on long divergent pairs, and a seed '
  'heuristic for gap-affine costs, at which no tool completed a 1 Mbp pair beyond 1% errors.',
 'Measurements on further processors, and an evaluation of the GPU kernels the library provides for batches '
  '(Supplementary Section S1.5).']),
]

BACK = [
    ("Author contributions", "Konstantinos Kyriakidis (Conceptualization [lead], Methodology [lead], Software [lead], "
     "Validation [lead], Formal analysis [lead], Investigation [lead], Data curation [lead], Visualization [lead], "
     "Writing—original draft [lead]), Benedict Paten (Supervision [lead], Writing—review & editing [equal])."),
    ("Supplementary data", "Supplementary data are available at *Bioinformatics* online."),
    ("Conflict of interest", "None declared."),
    ("Funding", "This work received no specific funding."),
    ("Data availability", "The software, benchmark harness and the commits of every compared tool are available at "
     "https://github.com/kokyriakidis/dinara-align. The benchmarks are reproduced by `pixi run bench` (unit and "
     "`--affine` costs), `bench-modes`, `bench-batch` and `bench-ablation`, and the raw results behind every table "
     "are in `benchmarks/results/2026-10-11`. Release builds target the baseline "
     "processor of each platform; the AVX-512 code paths require a build for the host processor."),
]

TABLES = {
"engines": {
  "caption": "Algorithms in dinara-align and the changes introduced in this work.",
  "header": ["Problem", "Algorithm", "Changes (section)"],
  "align": "lll",
  "wrap": [1, 2],
  "rows": [
    ["Unit costs, divergent pairs", "Band doubling over bit-parallel tiles with the gap-chaining seed heuristic (A*PA2)", "Regrouped recurrence, interleaved groups (3.2); agreeing estimates (3.3); direct neighbours (3.5); traceback within a memory limit (3.7)"],
    ["Unit costs, near-identical pairs", "Diagonal transition (WFA, BiWFA)", "Gathered extension (3.4); estimate-delayed hand-over (3.10); per-pair overhead (3.11)"],
    ["Gap-affine, two-piece", "Wavefront, one search or two (WFA, BiWFA); full sweep in difference recurrences for divergent pairs (Suzuki and Kasahara, 2018)", "Gathered extension (3.4); one search or two (3.9); estimate-delayed hand-over (3.10)"],
    ["Substitution matrices", "Anti-diagonal vector sweep (Wozniak, 1997)", "Score-certified band, cache-sized strips (3.8)"],
    ["Local alignment", "End found by a forward sweep (Zhao et al., 2013)", "Backward wavefront from the best end (3.8)"],
    ["Batches of short pairs", "Inter-sequence vectorisation (Rognes, 2011; Rahn et al., 2018; Daily, 2016)", "Certified per-lane bands (3.6)"],
  ]},
"unit": {
  "fit": True,
  "caption": "Unit-cost alignment of the real datasets by every exact aligner: mean time per alignment in ms, traceback included, and below it the peak resident memory of the run in MiB, as QuickEd reports it. t.o. marks a tool that completed no pair within 30 s, o.o.m. one the 8 GiB limit stopped, and crash one that crashed: BSAlign's edit distance crashes on pairs of about 150 kbp or more, with or without the limit, on the eighth PromethION pair among them. Times are medians of three rounds: the other tools' alternated between them, and dinara-align's were run a day later at its final version, on the same machine, builds and samples. The simulated datasets are in Fig. 2 and Supplementary Table S5. Bold marks the lowest mean of the tools that completed the most pairs.",
  "header": ["Dataset", "dinara-align", "A*PA2-full", "A*PA2-simple", "A*PA", "QuickEd", "Edlib", "BiWFA", "WFA", "BSAlign", "KSW2"],
  "align": "lrrrrrrrrrr",
  "rows": [
    ["*Mean time per alignment*", "", "", "", "", "", "", "", "", "", ""],
    ["ont-1k (830 bp)", "**0.0170**", "0.0810", "0.0481", "0.613", "0.0732", "0.124", "0.0442", "0.0313", "0.119", "0.538"],
    ["ont-10k (3.6 kbp)", "**0.143**", "0.393", "0.247", "15.4", "0.679", "1.11", "0.752", "0.599", "1.30", "13.0"],
    ["ont-50k (9.6 kbp)", "**0.535**", "1.28", "0.941", "263", "2.31", "6.11", "7.43", "6.35", "12.3", "169"],
    ["ont-promethion (45 kbp)", "**4.11**", "7.02", "6.76", "929", "17.1", "73.9", "71.1", "79.7", "181 [7/22]", "1,810 [15/22]"],
    ["ont-500k (642 kbp)", "**116**", "284", "688", "t.o.", "848", "5,000", "6,370", "o.o.m.", "crash", "o.o.m."],
    ["ont-500k-genvar (663 kbp)", "**154**", "225", "512", "t.o.", "1,390", "5,230", "6,180", "o.o.m.", "crash", "o.o.m."],
    ["SARS-CoV-2 (30 kbp)", "**0.237**", "2.01", "0.727", "6.55", "3.86", "8.18", "0.906", "0.979", "47.0", "68.6"],
    ["Illumina (248 bp)", "**0.000441**", "0.0145", "0.00411", "0.0137", "0.0162", "0.0234", "0.00206", "0.00112", "0.0205", "0.0124"],
    ["PacBio HiFi (14 kbp)", "**0.0116**", "0.980", "0.259", "1.06", "0.951", "3.41", "0.0865", "0.0660", "12.2", "3.09"],
    ["*Peak resident memory (MiB)*", "", "", "", "", "", "", "", "", "", ""],
    ["ont-1k", "14", "5", "5", "6", "2", "7", "7", "8", "2", "7"],
    ["ont-10k", "16", "6", "6", "39", "6", "8", "8", "39", "25", "72"],
    ["ont-50k", "19", "8", "7", "858", "18", "8", "9", "187", "423", "650"],
    ["ont-promethion", "32", "14", "9", "455", "19", "8", "14", "1069", "crash", "—"],
    ["ont-500k", "71", "89", "90", "t.o.", "26", "15", "29", "o.o.m.", "crash", "o.o.m."],
    ["ont-500k-genvar", "81", "56", "55", "t.o.", "27", "17", "27", "o.o.m.", "crash", "o.o.m."],
    ["SARS-CoV-2", "18", "7", "6", "19", "15", "8", "9", "26", "215", "151"],
    ["Illumina", "12", "5", "5", "5", "2", "7", "7", "7", "2", "7"],
    ["PacBio HiFi", "13", "7", "6", "6", "4", "8", "8", "9", "127", "21"],
  ]},
"affine_sim": {
  "caption": "Gap-affine alignment of QuickEd's simulated datasets, random sequences with uniformly spread substitutions, insertions and deletions, at costs (4, 6, 2) by every exact aligner that takes them: mean time per alignment in ms, a 30 s budget per sample. Brackets give the pairs completed where not all were; t.o. marks a tool that completed none and o.o.m. one the 8 GiB limit stopped. Datasets no tool completed are left out; every dataset is in Supplementary Table S4. Bold marks the lowest mean of the tools that completed the most pairs.",
  "header": ["Dataset", "dinara-align", "WFA", "BiWFA", "BSAlign", "KSW2"],
  "align": "lrrrrr",
  "rows": [
    ["10 kbp, 1%", "**0.186**", "0.409", "0.649", "37.0", "21.2"],
    ["10 kbp, 5%", "**2.29**", "7.13", "8.73", "37.2", "97.6"],
    ["10 kbp, 10%", "**6.96**", "22.8", "27.2", "37.2", "204"],
    ["10 kbp, 20%", "**19.7**", "65.5", "75.2", "37.5", "375 [80/99]"],
    ["100 kbp, 1%", "**11.7**", "39.9", "39.0", "o.o.m.", "2,440"],
    ["100 kbp, 5%", "**327**", "810", "715", "o.o.m.", "9,220 [3/9]"],
    ["100 kbp, 10%", "**1,170**", "o.o.m.", "2,430", "o.o.m.", "17,400 [1/10]"],
    ["100 kbp, 20%", "**3,560 [8/10]**", "o.o.m.", "7,350 [4/10]", "o.o.m.", "o.o.m."],
    ["500 kbp, 1%", "**463**", "1,290", "942", "o.o.m.", "o.o.m."],
    ["500 kbp, 5%", "**9,710 [3/4]**", "o.o.m.", "19,400 [1/4]", "o.o.m.", "o.o.m."],
    ["1 Mbp, 1%", "**2,230**", "o.o.m.", "3,920", "o.o.m.", "o.o.m."],
  ]},
"affine": {
  "caption": "Gap-affine alignment of the real datasets at costs (4, 6, 2) by every exact aligner that takes them: mean time per alignment in ms, a 30 s budget per sample, the median of three rounds, the other tools' alternating between them and dinara-align's a day later at its final version, on the same machine, builds and samples. Brackets give the pairs completed where not all were; t.o. marks a tool that completed none and o.o.m. one the 8 GiB limit stopped. Every dataset is in Supplementary Table S4. Bold marks the lowest mean of the tools that completed the most pairs.",
  "header": ["Dataset", "dinara-align", "WFA", "BiWFA", "BSAlign", "KSW2"],
  "align": "lrrrrr",
  "rows": [
    ["ont-1k", "**0.102**", "0.189", "0.305", "0.336", "1.10"],
    ["ont-10k", "**1.57**", "4.86", "5.85", "7.83", "29.9"],
    ["ont-50k", "**21.3**", "63.4", "54.9", "86.8", "341 [81/104]"],
    ["ont-promethion", "**261**", "286 [7/22]", "576", "o.o.m.", "2,350 [7/22]"],
    ["ont-500k", "t.o.", "o.o.m.", "t.o.", "o.o.m.", "o.o.m."],
    ["ont-500k-genvar", "21,500 [1/4]", "o.o.m.", "t.o.", "o.o.m.", "o.o.m."],
    ["SARS-CoV-2", "**0.674**", "2.06", "1.78", "341", "88.6"],
    ["Illumina", "**0.00116**", "0.00154", "0.00282", "0.0538", "0.0128"],
    ["PacBio HiFi", "**0.0787**", "0.224", "0.232", "86.1", "9.74"],
  ]},
}

FIGURES = {
"overview": "How dinara-align assigns a pair to an algorithm, and the rule behind each hand-over. At unit costs, diagonal transition from one end and then from both aligns near-identical pairs, handing over to band doubling over bit-parallel tiles only once its projected work exceeds the band's and it has spent one twenty-fourth of the cheapest band (Section 3.10); the tiles are traced back within the memory limit. At gap-affine costs one search grows alone, its fronts kept, and is traced back at once if it reaches the end; a second search joins once one would step more diagonals than two would need (Section 3.9), and a global pair projected to take more work than a full sweep in difference recurrences is handed to that sweep after one twenty-fourth of its work, the alignment the same either way. Batches of short pairs are aligned one pair per vector lane within a band whose optimality each lane certifies (Section 3.6). The figure shows alignments; a unit-cost distance skips the one-ended search.",

"sweeps": "Unit-cost alignment of QuickEd's simulated pairs, mean time per alignment on a logarithmic scale. (a) 100 kbp pairs by error rate. (b) Pairs at 20% errors by length; lines end where a tool completed no pair within 30 s or the 8 GiB limit stopped it, and hollow points mark tools that completed only some pairs. BSAlign and KSW2 are left out for legibility; their times are in the repository's benchmark tables. dinara-align had the lowest mean time at every length and error rate.",
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
 "Zook2016": ("Zook et al.", "2016"), "Shafin2020": ("Shafin et al.", "2020"), "QuickEdData": ("Doblas et al.", "2024"), "PaBench": ("Groot Koerkamp", "2024"), "Hansen2025": ("Hansen et al.", "2026"), "Prjibelski2019": ("Prjibelski et al.", "2019"),
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
 "Zook2016": "Zook2016giab", "Shafin2020": "Shafin2020shasta", "QuickEdData": "Doblas2024datasets", "PaBench": "GrootKoerkamp2024pabench", "Hansen2025": "Hansen2025hg002", "Prjibelski2019": "Prjibelski2019sequence",
}
