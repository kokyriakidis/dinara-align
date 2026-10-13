"""The supplementary material, in the same markup as main_text.py."""
from change_record import SUPPLEMENT as _OLD

TITLE = "Supplementary material for “dinara-align: fast, hardware-aware exact pairwise alignment that adapts its algorithm to each pair”"

# Table S1's rows are kept with the change record they came from; only labels that read informally are reworded.
_changes = next(b for b in _OLD if b[0] == "table")
_LABELS = {"Prefix search gives up a doomed try": "Prefix search abandons a try that cannot succeed"}
for _row in _changes[3]:
    _row[1] = _LABELS.get(_row[1], _row[1])
# The changes since the change record, each measured on the Skylake-X when it was made.
_changes[3].extend([
    ["a329af5", "Table sweeps in 16-bit lanes where the straight path bounds their scores", "table-global 544 → 352 µs", "Skylake-X"],
    ["0395d7e", "Banded traceback keeps each cell's decision in a byte", "table-global 352 → 253 µs", "Skylake-X"],
    ["c5aee28", "Table global alignment traced in bands that prove themselves", "table-global 243 → 157 µs", "M2"],
    ["691f43a", "Divergent gap-affine pairs handed to a full sweep in difference recurrences", "unrelated 1, 2, 5 kbp at (4, 6, 2) 1.50, 5.3, 32.9 → 0.28, 0.87, 4.8 ms", "Skylake-X"],
    ["64be0ec", "Anti-diagonal sweeps a 16 KiB strip of rows at a time", "local, 10 kbp against 12 kbp 50.0 → 17.4 ms", "Skylake-X"],
    ["ffbf46e", "Banded traceback, local end search and linear-space halves in place", "divergent table pairs 18 to 37% faster", "Skylake-X"],
    ["94bf719", "Unit-cost traceback within the memory limit, stretches swept again", "1 Mbp at 20%: peak growth 175 → 133 MiB a pair; 3.35 → 3.64 s for four (with 5742ef7, 92deb39)", "Skylake-X"],
    ["289bcf0", "Gap-affine hand-over only after one twenty-fourth of the sweep's work", "ont-10k at (4, 6, 2) 2.33 → 1.56 ms; ont-50k 30.5 → 21.0 ms", "Skylake-X"],
    ["ea3f09e", "256-column tiles for every band", "10 kbp at 5% 318 → 254 µs; ont-10k 161 → 147 µs", "Skylake-X"],
    ["dc980d7", "Inexact seeds from 128 kbp on AVX-512", "100 kbp at 10% 10.3 → 6.8 ms", "Skylake-X"],
    ["7aab823", "A space holds its profile directly, not in an Optional, and an aligner's traceback reuses it", "250 bp Illumina, through one aligner 0.97 → 0.80 µs; distance 1.21 → 0.68 µs", "Skylake-X"],
    ["90e0040", "An aligner keeps the fronts it traces back", "250 bp Illumina 0.82 → 0.65 µs", "Skylake-X"],
    ["32e78aa", "A gap-affine alignment begins with the one search it traces", "250 bp Illumina at (4, 6, 2) through one aligner 2.41 → 1.52 µs (a fresh call 5.16 µs before)", "Skylake-X"],
    ["7d38a99", "A reversed pair coded back to front in the coding pass", "250 bp Illumina 0.649 → 0.583 µs", "Skylake-X"],
    ["4558186", "Alignments pass between internal calls without an Optional", "250 bp Illumina 0.585 → 0.420 µs; at gap-affine costs 0.995 → 0.840 µs", "Skylake-X"],
    ["036e12b", "A profile coded 32 bytes at a time", "250 bp Illumina 0.421 → 0.400 µs; PacBio HiFi 5.7 → 4.95 µs", "Skylake-X"],
    ["a6d9166", "The first gap-affine search grows alone while it steps no more than two searches would", "250 bp Illumina at (4, 6, 2) 1.26 → 1.20 µs; pairs costing 9 to 32 a third less", "Skylake-X"],
    ["02a8455", "Alignments reuse the two-ended search's memory and write kept fronts in place; CIGARs written into the string; the diagonal transition abandons a distance it cannot reach once it has spent its share (with c8ae809 to 4f9aa01)", "250 bp Illumina 0.582 → 0.444 µs; ont-1k 20.4 → 16.8 µs; 100 kbp at 1% 1.43 → 1.13 ms; 500 kbp at 1% 32.3 → 24.1 ms", "Skylake-X"],
    ["04233cf", "An alignment's hand-over charges the band its recording and traceback", "10 kbp at 5% 250 → 193 µs", "Skylake-X"],
    ["383bbc5", "A distance's hand-over charges the band as measured, 0.6 times the earlier model", "distances: 1 Mbp at 1% 47.0 → 43.0 ms; 10 kbp at 10% 195 → 186 µs", "Skylake-X"],
    ["04906ec", "Gap-affine traceback writes runs of moves sixteen at a time into the alignment's own buffer (with df44eaf)", "250 bp Illumina at (4, 6, 2) 1.36 → 1.14 µs; PacBio HiFi 86.6 → 77.4 µs", "Skylake-X"],
    ["d9dc492", "An untrusted band follows a failed round's estimate once the round crossed one eighth of the columns", "10 kbp at 10% 341 → 310 µs, at 20% 497 → 437 µs; ont-50k 557 → 528 µs", "Skylake-X"],
    ["5a500f8", "A CIGAR's diagonal runs split by a 64-base mismatch mask; moves already in order written where they lie", "ont-10k 141 → 139 µs; 10 kbp at 20% 435 → 428 µs", "Skylake-X"],
])

BODY = [
('p',
 'This supplement gives the implementation details, per-change measurements and complete results behind the main '
 'text. Section S1 describes implementation details and how the library was validated (Section S1.7); Section S2 '
 'the measured effect of each change and of switching each technique off; Section S3 the complete results at '
 'unit and gap-affine costs, the other alignment modes, the samples, memory and batches; Sections S4 and S5 '
 'robustness measures and approaches evaluated and not adopted; and Section S6 the proofs of Lemma 1 and '
 'Propositions 1 and 2 of the main text, and of Lemma S1. Unless stated otherwise, the measurements in Sections '
 'S1, S4 and S5 and in Table S1 were taken when each change was introduced, on the machine indicated, against '
 'the preceding version; Table S2 was measured at commit `dc980d7`. The results of Section 4 of the main text '
 'and of Section S3 were obtained on the Skylake-X at commit `87d365f`, every tool timed in three alternating '
 'rounds. The machines are an Intel Core i9-7900X (Skylake-X, AVX-512) and an Apple M2 (NEON), and in Section '
 'S1.5 an NVIDIA GeForce RTX 2070. Memory is given in MiB; references are listed in the main text.'),
('sec', 'Implementation details'),
('sub', 'Unit costs'),
('p',
 '**Engine selection.** The diagonal transition receives a step budget equal to the cost of the band it would '
 'replace, in steps: for a band over \\(n\\) columns at an estimated distance \\(\\hat D\\), \\(n(a + \\hat D/b)\\) '
 'steps, the band growing taller with the distance. Both terms were measured on the Skylake-X against bands run '
 'to completion: for a distance \\(a = 3.6\\) and \\(b = 138\\), and for an alignment, whose band also records '
 'its tiles and is traced back, \\(a = 8\\) and \\(b = 52\\). Charging alignments the cost of a distance had '
 'abandoned pairs that the search finishes in half the time of the band, and 10 kbp pairs at 5% divergence took '
 '250 instead of 193 µs. From 128 kbp on '
 'AVX-512 and 64 kbp with narrower vectors, inexact seeds replace exact seeds when fewer than 20, 30 or 40% of '
 'the exact seeds chain from the origin on AVX-512, AVX2 and NEON, respectively: the wider the vector unit, the '
 'cheaper the bit-parallel band relative to the scalar construction of inexact seeds. All thresholds are fixed '
 'at compile time for each instruction set, so that none is dispatched at run time.'),
('p',
 '**Register and memory layout of the band.** An empty inline-assembly barrier prevents the compiler from '
 'reassociating the regrouped recurrence back to its original form. Vertical differences remain in registers '
 "across a word's columns, and the horizontal differences of successive columns are stored at distinct "
 'addresses, so that no load waits on a preceding store. Bands of five to seven words combine a four-lane vector '
 'with scalar words.'),
('p',
 '**Symbols other than ACGT.** `N` and other non-ACGT symbols are encoded in a third bit plane. Each tile uses '
 'only the planes its symbols require: two when neither its columns nor any row contains such a symbol, the '
 'third plane of the rows as a mask when only rows contain one, and all three otherwise. Bases are coded 32 '
 'bytes at a time, and symbols and column planes 16 bytes at a time; only chunks that contain such a symbol are '
 'processed byte by byte. This made reads containing `N` 7 to 9% faster on both machines. A seed that contains '
 '`N` is not counted, so that pairs containing `N` can still use seeds; this made the ultra-long reads 2.5 to 3 '
 'times faster on the M2 and 1.4 to 1.7 times faster on the Skylake-X.'),
('p',
 '**Tile width.** Tiles are 256 columns wide, as in A*PA2, whatever the height of the band. Earlier versions '
 'used 64-column tiles in bands bounded below 1,024, a rule kept from the first version but never measured: '
 'narrower tiles reduce the slope a band adds across each tile, but each tile is pruned, recorded and traced '
 'back separately. With one width, 10 kbp pairs at 5% divergence aligned in 254 µs instead of 318, the ont-10k '
 'reads in 147 instead of 161 µs, and every other dataset but one in the same time or up to 6.5% faster; 10 kbp '
 'pairs at 1%, settled mostly by the diagonal transition, were within 3%.'),
('p',
 "**Traceback within a memory limit.** The band keeps every tile boundary while they fit within the alignment's "
 'memory limit (main text, Section 3.7), and otherwise every \\(k\\)-th, \\(k\\) doubling as the round proceeds, '
 'leaving room for one stretch of the widest boundaries. The state each tile records is its words, the score at '
 'its top, the deepest kept row and its score, the bound and the checkpoints passed. The traceback, from the '
 'last tile back, sweeps a stretch again from the last kept boundary before a tile it lacks; kept boundaries '
 'form a stack in tile order and are released as they are traced. The boundary buffers grow by doubling up to '
 'half the limit and then to the whole of it at once, since doubling further would hold the old buffer and its '
 'copy together, up to twice the limit. On the four 1 Mbp pairs at 20% errors, measured one per process, an '
 'alignment grew the peak memory by 133, 98, 124 and 113 MiB instead of 175, 99, 143 and 146 MiB, in 3.64 s '
 'instead of 3.35 s for the four, and pairs whose boundaries fit were unchanged. With a limit of zero, '
 'alignments took about twice as long.'),
('p',
 '**Retries.** Section 3.3 of the main text sets the first bound. A retry raises the bound by at least a quarter '
 'and, unless the two estimates agreed or the failed round crossed one eighth of the columns, at most doubles it, '
 'and a failed round grows its bound from the one it ended on, which a checkpoint may have lowered, rather than '
 'the one it started on. Before these rules, a 30 kbp SARS-CoV-2 pair at distance 49 had been extrapolated to '
 '14,448 and computed over the full matrix. A round that dies early on a read whose edits gather at its start '
 'extrapolates from those edits alone, so its estimate is used only once the round has crossed an eighth of the '
 'columns, the share a seeded round must also cross; doubling regardless, 10 kbp pairs at 10% divergence whose '
 'estimates disagreed had spent 37% of their band computation in failed rounds, and the 10 kbp pairs at 10 and '
 '20% took 341 and 497 µs instead of 310 and 437.'),
('p',
 'When they were introduced, these rules reduced the mean time on the ont-500k-genvar reads from 4.33 s, with 13 '
 'of 15 reads completed within the time limit, to 754 ms. At commit `dc980d7`, where Table S2 was measured, the '
 'seed heuristic sets the first bound on long reads, and switching the agreement test off changes the ultra-long '
 'reads by less than 1% but slows the ont-10k reads by 22% (Section S2.1).'),
('p',
 '**The bottom of the next tile.** The rows of a tile below the deepest row kept in the previous tile are '
 'bounded before the tile is swept: a cell \\(j\\) rows below the diagonal through the deepest kept cell costs at '
 "least \\(j\\) more than that cell's score, so the band stops where that plus the heuristic exceeds the bound, as "
 'A*PA2 bounds a block.'),
('p',
 '**Retries with seeds.** With seeds, a retry aims one quarter of the extrapolated climb above the heuristic at '
 'the origin beyond the estimate; with one half, a round had failed after 90% of the columns of a 527 kbp read. '
 'A seeded round that fails within the first eighth of its columns does not extrapolate; instead, its margin '
 'over the heuristic bound at the origin increases fourfold, as in A*PA2. Previously, an 884 kbp read had failed '
 'at column 1,280 and been extrapolated to 125,000 for a distance of 75,000. Only the first round revises its '
 'bound at the checkpoints; when later rounds also lowered their bounds, they failed late, four times on one 897 '
 'kbp read.'),
('p',
 "**Traceback.** Each tile is traced by a wavefront search from the traced cell back to the tile's left edge; a "
 'tile with more edits than that search is allowed is swept again over a window of rows above the traced cell, '
 'starting at four words and doubling, as in A*PA2. The top boundary of the window is initialised to +1, so that '
 'its scores correspond to real paths. A 600 kbp nanopore pair at 13% divergence was aligned in 617 ms instead '
 'of 1.38 s, of which 840 ms had been spent recomputing tiles of 900 words. Recomputation buffers are reused '
 'without initialisation, and tile boundaries are allocated once and written through pointers; loading the '
 'boundaries had accounted for a quarter of the traceback time. The CIGAR string is produced by reversing moves '
 '16 at a time and splitting each diagonal run into matches and substitutions from a mask of the mismatches of '
 "64 bases at a time, whose trailing zeros give each piece's length, which costs 1 to 3% of the alignment time."),
('p',
 '**Diagonal transition.** Wavefronts are padded with unreached diagonals and advanced eight diagonals at a time '
 'with branch-free validity masks, in a separate function; a step costs about 1.3 ns, compared with 5 ns before. '
 'The two-ended search tests the overlap of the two wavefronts with vector comparisons and requires about half '
 'as many steps. Reducing the step budget (Engine selection, above) from 12 to 5 steps per column reduced the '
 'time for 10 kbp pairs at 5% divergence from 241 to 144 µs. The choice between one and two wavefronts was '
 'fitted on 1,300 pairs of 1 to 30 kbp. For pairs of up to 2,048 columns, the band had been set wide across the '
 'whole matrix, and 42% of such pairs were assigned to a band although the diagonal transition would have been '
 'faster; such a band is now charged an additional 30,000 steps, which made the ont-1k reads 25% faster on the '
 'Skylake-X. A short noisy pair abandons the diagonal search in favour of the full matrix (a 1 kbp pair at 15% '
 'divergence, 19.2 to 11.8 µs), and the first bound of a short pair is 1.7 times its extrapolation, at most 128 '
 'above it.'),
('p',
 '**Gather masks.** A gather consumes its mask register, and compilers generate the required all-ones mask with '
 '`kxnor`, which depends on the previous value of that register; consecutive gathers were therefore serialised '
 "at about 21 cycles per group. Deriving the mask from the group's own lanes removes this dependency. Gathered "
 'extension made the ont-1k reads, the SARS-CoV-2 genomes and 100 kbp pairs at 1 to 2% divergence 7 to 19% '
 'faster (Table S1, commit `45697e0`).'),
('p',
 '**Measured thresholds.** A full bit-parallel sweep of a 64-row word costs 35%, 45% and 77% of a '
 'diagonal-transition step on AVX-512, AVX2 and NEON, respectively, measured on 1 kbp pairs at 15% divergence '
 '(the AVX2 value on the Skylake-X with code compiled for Haswell). On AVX-512, seeds reduce the running time '
 'only from about 86 kbp: below that length, the Skylake-X aligned real reads of 16 to 64 kbp 14% faster without '
 'seeds, and uniform 30 kbp pairs at 5% divergence 30% faster, whereas at 100 kbp seeds reduced the time by 7 to '
 '24%. With AVX2 and NEON, seeds are used from 16 kbp. For pairs whose two extrapolations agree, as with uniform '
 'errors, the inexact-seed threshold is 20% on every platform, and an extrapolated error rate of at least one '
 'edit in ten bases selects inexact seeds directly from 128 kbp on AVX-512 and from 64 kbp with narrower '
 'vectors: once every band used 256-column tiles, 100 kbp pairs at 10% aligned in 6.8 ms on exact seeds against '
 '10.3 ms on inexact ones on the Skylake-X, while 500 kbp pairs at 10 and 20% ran 7 to 10% slower without '
 'inexact seeds.'),
('p',
 '**Every threshold measured again.** Each of the 41 thresholds that select an algorithm or bound a search, from '
 "the step budgets of the diagonal transition to the seed gates and the traceback's limits, was halved and "
 'doubled in turn, and each of the 82 builds was timed on 17 unit-cost samples, with and without traceback, or '
 'on 8 gap-affine ones, with the unchanged build timed after every six; its spread was 1.2% in the median. About '
 'thirty thresholds moved no sample by more than 3% either way, and most of the others traded one dataset '
 'against another, such as the first bound of an untrusted band, or the agreement and give-up rules for '
 'extrapolations, and were left as they were. One had become outdated: the length from which inexact seeds are '
 'used (Measured thresholds, above). Two more outdated choices were found by examining the slowest results: '
 '64-column tiles for narrow bands (Tile width, above) and the projection that handed nanopore reads to the full '
 'sweep (Divergent pairs at gap-affine costs, Section S1.3).'),
('p',
 '**Interleaving.** llvm-mca estimates the interleaved AVX-512 loop at 24 cycles per column, against a '
 'throughput bound of 11, so the loop remains limited by the lane rotation. A third group still fits in the 32 '
 'AVX-512 registers, and llvm-mca predicted 15% fewer cycles per group, but the measured time changed by less '
 'than 1.5%.'),
('sub', 'Seed heuristic'),
('p',
 '**Batched lookups.** Query rows are processed in batches of 16: the bucket bounds of all rows are loaded '
 'first, so that the loads overlap, and their entries are prefetched; each entry stores a seed together with its '
 'code, so that a single load provides both. A lookup tests each candidate window with a branch-free expression. '
 'Together with considering only seeds that can still chain to the end, these steps reduced the inexact setup on '
 'the ont-500k-genvar reads from 5.02 to 2.90 s in total. On 100 kbp pairs at 10% divergence, inexact-seed '
 'construction accounted for 58% of the running time on the Skylake-X; there exact seeds now win (Measured '
 'thresholds, above).'),
('p',
 '**Admissibility check and correction.** Comparing the heuristic with the true remaining cost at every cell of '
 'random pairs showed that an earlier version of our implementation overestimated the remaining cost by one on '
 'optimal paths in 1 to 2% of pairs. The cause was a filter that discarded matches whose gap to the end of the '
 'pair exceeded the potential of the intervening seeds: a match ending one diagonal beyond this limit was '
 'discarded, although a path through it costs only one more than the potential, and the neighbour that would '
 'otherwise bound such paths could be removed by local pruning. Such a match now chains to the end at the cost '
 'of its excess: an exact match of an inexact seed, which saves 2, one diagonal beyond the limit saves 1 and '
 'starts in layer 1, and no further match can follow it. After the correction, no cell exceeded the true cost in '
 '2,200 random pairs, compared with 34,450 cells in 400 pairs before. The reported distances agreed with those '
 'of the other aligners in all benchmarks, both before and after the correction. The overestimate arose from our '
 'filter; we do not claim that it applies to A*PA or A*PA2.'),
('p',
 '**Exact-seed filter.** A 32-bit filter per seed, placed before a half-full open-addressing hash table, '
 'resolves most lookups with a single predictable test. On the M2, matching a sample of 320 reads took 32 ms '
 'instead of 87 ms, and 30 kbp pairs at 5% divergence became 15% faster.'),
('p',
 '**Local pruning.** The diagonal-transition search used for local pruning updates its wavefronts in place, with '
 'an unreachable diagonal on either side replacing three bounds checks per diagonal.'),
('sub', 'Gap-affine costs'),
('p',
 "**Divergent pairs at gap-affine costs.** The wavefront's work grows with the cost times the length, so on "
 "unrelated pairs it is slower than an 8-bit difference sweep such as BSAlign's. A global pair with the whole "
 'matrix open is therefore handed to a full sweep in difference recurrences (Suzuki and Kasahara, 2018), in '
 "8-bit lanes at any length, once its projected work exceeds the sweep's, at seven sweep cells per wavefront "
 "step, the rate measured at 10 kbp. The sweep keeps each cell's predecessor flag under the tie rule of the "
 'wavefront, so the alignment is identical byte for byte, and for the traceback it keeps every \\(k\\)-th '
 'anti-diagonal, \\(k\\) the square root of the anti-diagonals times the bytes a cell takes, and sweeps each '
 'stretch again. Each anti-diagonal is swept in place from the last row up. On unrelated pairs of 1, 2 and 5 kbp '
 'at (4, 6, 2) alignment took 0.28, 0.87 and 4.8 ms instead of 1.50, 5.3 and 32.9 ms. The projection '
 'extrapolates the cost from how far the two wavefronts have advanced, which on nanopore reads, whose errors '
 'gather at the ends, ran up to 2.8 times too high in cost and 8 times in work, and handed reads to sweeps that '
 'took up to eight times longer than their wavefront: the ont-10k reads took 2.33 ms instead of 1.56. As for the '
 "diagonal transition (main text, Section 3.10), the projection is now trusted only once the wavefront's work "
 "has reached one twenty-fourth of the sweep's."),
('p',
 '**Compact traceback.** For the traceback of the wavefront, dinara-align stores, for every cost and diagonal, '
 'the furthest-reaching column in four bytes and one byte of flags recording the predecessor of every layer: '
 'five bytes per diagonal, compared with twelve in the high-memory mode of WFA. When the stored wavefronts would '
 'exceed the memory limit (80 MiB by default), the pair is split at a cell through which an optimal path passes. '
 'Under a general substitution matrix, substitution scores are obtained by byte shuffles (`pshufb`, `tbl`) for '
 'matrices of up to 16 entries and by gathers otherwise.'),
('p',
 '**Gathers.** In the gap-affine wavefront, AVX-512 gathers reduced the time for the ont-1k reads from 149 to '
 "129 µs. Deriving each gather's mask from the group's own lanes, with one reduction after the loop rather than "
 'one per group, reduced it further to 111 µs, and that for the ont-10k reads from 2.34 to 1.84 ms. Extension '
 'along matches was moved into the wavefront step, with columns stored through scalar loads and stores rather '
 'than reinserted into a vector, and retained wavefronts are held in fixed blocks that are never moved; memory '
 'moves had accounted for 7.5% of profiling samples (Table S1).'),
('p',
 '**Memory layout.** The wavefront updates the bounds and furthest reaches of its slots at every cost. When '
 "these were separate arrays of a few words, the allocator placed them next to another thread's data, so that "
 'the cache line moved between cores at every write; ten threads on the Skylake-X then required twice the cycles '
 'per pair of one thread, with 40 cache lines fetched from another core per pair. Each such array now occupies '
 'cache lines that no other allocation shares. A worker retains its profile, reversed sequences and search '
 'buffers between pairs, so that a warm worker performs no allocation; per-pair allocations had accounted for '
 'one third of the time of a short read, and caused contention in the allocator with many threads. Padding each '
 'row of a wavefront buffer by one cache line, to avoid 4 KiB aliasing between loads and earlier stores, removed '
 '1.5 billion blocked loads over ten 100 kbp pairs. Separating the working memory of threads by two cache lines '
 'rather than one reduced the cache lines fetched from other cores from 3.4 million to 0.3 million over 2 '
 'million reads on ten threads.'),
('sub', 'Other modes'),
('p',
 '**Extension.** Extension alignment uses a single wavefront search from the anchor, with the match score folded '
 'into the costs, and records the best end point on the furthest anti-diagonal of each wavefront. The search '
 'stops once a proven bound shows that no higher cost can yield a better score, and the prefix up to the best '
 'end point is then aligned globally. With an end bonus, the search also records the best point on the last row '
 'of the query and continues only while such an alignment can still come within the bonus; this reduced the time '
 'for an extension with an end bonus from 1.02 ms to 186 µs.'),
('p',
 '**Ends-free alignment.** Under gap costs with free ends, the forward search starts at cost zero on every '
 'diagonal on which a free start lies. The span is determined first: one search from the free side finds the end '
 'on the highest diagonal, and a narrow search backwards finds the start; the span is then aligned globally. On '
 'x86, this aligns a read at the start of a reference, an overlap or a suffix 1.4 to 3.5 times faster, and a '
 'read anywhere within the reference 0 to 14% slower, than the previous approach.'),
('p',
 '**Infix and prefix alignment at unit costs.** Infix alignment uses the bit-parallel sweep with the cut-off of '
 'Edlib, with a first bound of 64; a failed try extrapolates the distance from the deepest row its band kept '
 'within the bound and the retry aims past it, which on the M2 took a 10 kbp infix from 3.39 to 1.99 ms; words '
 'entering the band are initialised to +1, and the last word advances one column at a time. Prefix alignment '
 'sweeps only the band permitted by its bound: the top of the band moves down with the diagonal, no column '
 'beyond the read length plus the best end so far is computed, an attempt whose band becomes empty stops, and at '
 'one eighth, one quarter and one half of the columns an attempt that can no longer succeed is abandoned. This '
 'reduced the time for a prefix distance from 23.7 to 20.5 µs on the Skylake-X.'),
('p',
 '**A band derived from the score.** A table with one match and one mismatch score is converted to costs and '
 'aligned by the wavefront (main text, Section 3.8). Under any other substitution matrix, a global alignment is '
 'traced in a band of diagonals that grows at most fourfold per step until its own score proves that it holds '
 "every optimal path, each cell's decision kept in one byte. With \\(\\alpha_{\\max}\\) the table's largest score, "
 'each gapped letter adds at least \\(2e+\\alpha_{\\max}\\) to \\(C = \\alpha_{\\max}(n+m) - 2S\\), and a path that '
 'strays \\(t\\) diagonals beyond those between the main and the end diagonal \\(\\delta\\) has at least \\(|\\delta| + '
 '2t\\) gapped letters, so the score of any alignment confines every optimal path to \\(t \\le '
 '(C/(2e+\\alpha_{\\max}) - |\\delta|)/2\\). A pair close to its diagonal thus costs about its own cells: with '
 'sweeps in 16-bit lanes and one-byte traceback decisions, 1 kbp pairs at 10% divergence under a '
 'transition-transversion table aligned in about 200 µs instead of 544 µs on the Skylake-X (Table S1). '
 'Anti-diagonals are swept in strips of rows of at most 16 KiB, which stay in the first-level cache; local '
 'alignment of 10 kbp against 12 kbp took 17.4 ms instead of 50.0 ms. The backward sweep from the end of a local '
 'alignment and each half of a linear-space alignment also use 16-bit lanes while the scores fit, which made 1 '
 'kbp local alignments 1.2 times and linear-space global alignments 1.4 times faster on the M2.'),
('sub', 'Batches'),
('p',
 '**Letter layout.** The pairs of a group are interleaved so that one 64-bit load reads eight letters of a lane, '
 'and eight lanes are transposed in registers by two shuffles.'),
('p',
 '**Alignment in 8-bit lanes.** Batch alignments remain in 8-bit lanes, 64 pairs per AVX-512 register, and a '
 'pair whose cost exceeds the 8-bit range is passed to a single call. With 16-bit lanes, half as many lanes '
 'would sweep the wide band that such a pair requires, which on the Skylake-X took as long as the searches on 1 '
 'kbp reads at 10% divergence, and twice as long on NEON.'),
('p',
 "**Cache-efficient traceback.** Storing the traceback flags of a group's whole band would exceed the caches on "
 'long pairs and touch new memory pages for every group; with ten threads, this doubled the time for 1 kbp pairs '
 'relative to the searches alone. The sweep therefore stores one row of every layer every \\(s\\) rows and '
 'recomputes a segment to obtain its flags during the traceback. With \\(R\\) rows, \\(s\\) is the square root of '
 '\\(R\\) times the bytes per cell of a stored row divided by those of a flag, which balances the stored rows '
 'against the flags of one segment and minimises memory use.'),
('p',
 '**Lane groups.** A group of 64 bytes occupies one register on AVX-512 and several on narrower vector units; '
 'with the 128-bit registers of the M2, 1 kbp pairs at 10% divergence required 73 µs per pair under gap-affine '
 'costs with four registers, compared with 132 µs with one. Local scores use 32 lanes of 16 bits or 16 lanes of '
 '32 bits, which on 128-bit registers spreads a group over four independent registers; on the M2, this scored '
 '20,000 references against a read in 83 ms, compared with 160 ms with one register. Pairs whose first band does '
 'not satisfy the bound are grouped by power-of-two classes of the width of their second band, and a group stops '
 'once every lane has exceeded its cost limit.'),
('p',
 '**Work shares.** Batch calls asked for more than one thread spread their pairs over them, longest first: the '
 'pairs are dealt into classes of lengths within an eighth of each other, with no sort, since sorting the short '
 'reads on one thread had taken a fifteenth of the work ten threads did on them. Each worker then takes a share '
 "of half of what is left divided among the workers, as OpenMP's guided schedule does, but never more than twice "
 'its last share, starting from one; a first share of a twentieth of the batch had left one worker holding the '
 "longest pairs while the others finished, adding half the batch's time on the Skylake-X. Building the packed "
 'batch is split across the threads in the same way; on one thread it had taken three times as long as the GPU '
 'took to score it.'),
('p',
 '**Free ends in lanes.** Ends-free alignments take the lanes too, their span found first as the tie rule picks '
 'it for a single pair: first the end, the cell furthest along the reference on which an optimal alignment may '
 'stop, found by a sweep whose band provably holds every optimal path; then the start, by a sweep back from that '
 'end over the reversed prefixes; then the global alignment of the span. On one core of the Skylake-X, 150 bp '
 'reads at 3% errors in 200 bp windows aligned at (4, 6, 2) in 8.0 µs per pair as an infix, against 26.6 µs one '
 'pair at a time, in 12.1 against 29.3 µs as a prefix, and in 20.0 against 50.3 µs with 60 letters free at the '
 'end of the window and the start of the read; at unit costs that overlap took 8.3 against 32.4 µs.'),
('p',
 '**GPU scores.** Batch scores and alignments under a substitution table, and unit-cost distances, also run on '
 'NVIDIA GPUs and Apple silicon. Each pair is first scored over a band of sixteen diagonals by one thread, its '
 "scores and gap layers one register per diagonal, swept anti-diagonal by anti-diagonal; the band's score stands "
 "where the pair's own score certifies it, by the bound of Proposition 1 of the main text with the table's least "
 'gap penalties, and the thread stops once the cheapest cell of its anti-diagonal costs more than every path off '
 'the band. The remaining pairs, about one in six of the short-read batch, are scored over the whole matrix '
 "while the host packs the next chunk, several to a warp, each lane holding a run of its pair's columns in "
 'registers and passing its last column on by shuffle. The cells of a global alignment are held shifted by their '
 "anti-diagonal times the gap extension, which reduces Gotoh's recurrence from five additions per cell to three. "
 'On an NVIDIA GeForce RTX 2070, the 500,000 short-read pairs were scored in 35 ms at (1, 2, 1), the median of '
 'three runs that varied by half, the whole call included and most of it the host packing the batch, against 29 '
 'ms on ten cores of the Skylake-X.'),
('sub', 'Code generation'),
('p',
 '**Avoiding division by a variable.** A 64-bit division by a variable takes tens of cycles on x86. The '
 'potential of the heuristic, required at every row that the band prunes, involves a division by the seed '
 'length, so each seed length has its own code path with a constant divisor, which compiles to a multiplication; '
 'batches likewise avoid a division per pair.'),
('p',
 '**Kernels compiled separately.** When inlined next to each other, or selected per band round, the bit-parallel '
 'kernels for plain bases and for additional symbols slowed the plain kernel by 1 to 4%; when inlined into the '
 'seed construction, the local-pruning search slowed the construction loop by a few per cent. Separate seeded '
 'and unseeded versions of the band keep the code of the heuristic out of unseeded pairs, which it had slowed by '
 '3 to 7%.'),
('p',
 '**Operand order.** Under a substitution matrix, the recurrence first takes the maximum of the two gap layers, '
 'so that the substitution score, the operand that arrives last, waits on only one maximum; a score obtained by '
 'comparison is available early and is used first.'),
('p',
 '**Branches.** On a 1 Mbp pair at 15% divergence, branch mispredictions account for about 7% of the cycles, '
 'mostly at short loop exits, compared with 2 to 3% in the bit-parallel sweep. A fully branch-free seed scan '
 'would spend about 25 million cycles per pair to save about 34 million, a gain of 1 to 2%, and was not adopted.'),
('sub', 'Validation'),
('p',
 '**Differential fuzzing.** `pixi run fuzz` draws random pairs with costs, modes, bands, cost limits and memory '
 'limits, three trials in ten with a small memory limit and every sixteenth as part of a batch, and checks each '
 'answer against a full-matrix implementation (`tests/oracle.mojo`) that shares no code with the library: the '
 'alignment must be valid and its cost or score optimal. An alignment under a cost limit must cost the same as '
 'without one, and the tie rule is checked at the default memory limit, since a pair split to fit a smaller one '
 'follows the rule within each piece: there the alignment under a cost limit must be the same as without one, '
 "and the default rule must give WFA2-lib's alignment of both sequences reversed, read backwards. A unit test "
 'checks that a unit-cost alignment is the same byte for byte under any memory limit. At commit `98fa91b`, whose '
 'library is the one Section 4 of the main text measures, 20,000 random cases (seed 20261011) passed, and so did '
 "WFA2-lib's regression set."),
('p',
 "**WFA2-lib's regression set.** `pixi run test-wfa` aligns the 305 pairs that WFA2-lib ships with its tests in "
 "every mode it tests, and compares each score with WFA2-lib's and each CIGAR, under WFA2-lib's tie rule, with "
 'the one WFA2-lib reports.'),
('sec', 'Measured effect of each change'),
('table', 'changes'),
('sub', 'Effect of each technique'),
('p',
 'To measure what each technique contributes to the final code, dinara-align was built with one technique '
 'switched off at a time, and every build aligned samples of the A*PA2 real datasets and of the 100 kbp '
 'simulated pairs and scored the short-read batch, the builds alternating within three rounds on one pinned core '
 'and assembled with branches kept clear of 32-byte boundaries, which the Skylake-X penalises (Table S2); every '
 'build returned the same costs. Without direct neighbour construction (main text, Section 3.5), the ultra-long '
 'reads took 9% and 16% longer; 100 kbp pairs, which take exact seeds on AVX-512, do not use it. Without the '
 'interleaved groups (Section 3.2), the ultra-long reads took 4.5 to 5% longer, 100 kbp pairs at 10% and 20% '
 'errors, whose bands are widest, 3.4% and 11.5%, and the other samples within 1.2%. The regrouped recurrence '
 "changed the samples by −1.4 to +4.1%, as bands narrow enough to be bound by a single word's step are a small "
 'part of the work. Without gathered extension (Section 3.4), the samples on which diagonal transition does most '
 'of the work took 6.4 to 26.5% longer, and without the agreeing estimates (Section 3.3), the ont-10k reads took '
 '22% and the SARS-CoV-2 genomes 4.9% longer. Without the certified bands (Section 3.6), the short-read batch '
 'took 2.9 and 3.2 times as long.'),
('table', 'ablation'),
('sec', 'Additional results'),
('p',
 'Table S3 gives the results for the other alignment modes, each compared on one pinned thread of the Skylake-X '
 'with the aligners that support the same mode. Table S4 gives the complete gap-affine comparison at costs (4, '
 '6, 2). Table S5 gives the mean and median time per alignment at unit costs on every dataset, Table S6 the '
 'samples, Table S7 the memory and Table S8 the short-read batch. The approximate aligners of the A*PA2 '
 'benchmark, WFA-adaptive and Block Aligner, are not evaluated, as they are not exact. At unit costs, KSW2 and '
 "BSAlign, which ran out of memory or crashed on every dataset of 500 kbp and more, are in the benchmark's own "
 'tables (benchmarks/README.md in the repository); Table S4 includes them.'),
('table', 'modes'),
('table', 'affine_full'),
('table', 'unit_full'),
('table', 'datasets'),
('table', 'memory'),
('table', 'batch'),
('sec', 'Robustness'),
('p',
 'Three changes address inputs on which the running time had been disproportionate:'),
('bullets', [
 '**Low-complexity sequence.** Wavefront layers that spread into a staircase pattern are handled explicitly, and '
 'the pruning search is bounded; a poly-A sequence of 32,000 bases aligned against one of 16,000 bases took 0.15 '
 's instead of 17 s.',
 '**Repetitive seeds.** Seeds are abandoned when they exceed 64 matches per seed on average; setting up a 200 '
 'kbp tandem repeat took 0.045 s and 34 MiB instead of 73 s and 760 MiB.',
 '**Unreachable wavefronts.** A wavefront that reaches no cell is not tested further; one such case took 13 ms '
 'instead of 1.5 s.',
]),
('sec', 'Approaches evaluated and not adopted'),
('bullets', [
 '**Multi-threading within a pair.** Striped bands over eight threads aligned 1 Mbp pairs at 15% divergence in '
 '575 ms, compared with 1,122 ms on one thread. The implementation was removed because a single alignment runs '
 'on the calling thread, and distributing pairs over threads, which batch calls do when asked, yields a larger '
 'gain.',
 '**One sixteen-lane group instead of two eight-lane groups.** Equal on long reads, and up to 2% slower on short '
 'divergent pairs.',
 '**A Bloom filter for inexact keys.** It halved the number of candidates, but reduced the alignment time by '
 'only 1.5 to 2%.',
 '**Shorter look-ahead, or no local pruning.** Without pruning, the ultra-long reads with and without genetic '
 'variation were 21% and 72% slower, respectively.',
 '**A* search within the diagonal transition, and an upper bound from beam search.** Each required several tuned '
 'parameters for small gains, and A* within the diagonal transition was 3 to 5 times slower on real reads.',
 '**Traceback during the sweep.** On the ultra-long reads and on 1 Mbp pairs at up to 15% errors, the stored '
 'tile boundaries occupied 0.6 to 13 MiB, compared with 10 to 17 MiB for the sequences and bit planes, so little '
 'memory could be saved; the traceback within a memory limit (Section S1.1) now bounds the worst case.',
 "**A tighter first bound.** With each pair's distance supplied by an oracle, starting band doubling at it saved "
 'up to 27% (100 kbp pairs at 20% errors), but starting 15% above it lost 18 to 23% on the ultra-long reads; '
 'capping the bound at an upper bound 5% above the distance saved up to 27% (1 Mbp pairs at 20%), and one 15% '
 'above about half that or nothing. No upper bound cheap and tight enough was found: beam search cost more than '
 'it saved, and a batch edit-distance certificate proved 578 of 70,519 pairs.',
]),
('sec', 'Proofs'),
('p',
 'This section gives the proofs of the statements in Sections 3.5, 3.6 and 3.9 of the main text, whose notation '
 'it uses, and a lemma on the pruning that Section 3.5 relies on.'),
('p',
 '**Lemma 1.** If the search from an exact match \\(\\mu\\) finds no path across the next \\(p\\) seeds, it finds '
 'none from any of its four neighbours.'),
('p',
 '*Proof.* Each neighbour \\(\\nu\\) matches the same seed as \\(\\mu\\), so its search crosses the same seeds against '
 'the same potentials. If \\(\\nu\\) shares the end of \\(\\mu\\), its search starts from \\(\\mathrm{end}(\\mu)\\) at '
 'cost 1 instead of 0. If \\(\\nu\\) shares the start, \\(\\mathrm{end}(\\nu)\\) lies on an adjacent diagonal, and '
 'after one edit the search from \\(\\mu\\) holds a point on that diagonal at or beyond \\(\\mathrm{end}(\\nu)\\), at '
 'cost \\(1 = c(\\nu)\\), which is kept because one edit is below the potential 2 of the seed itself. In both '
 'cases, at each cost and on each diagonal, the search from \\(\\mu\\) reaches at least as far as the search from '
 '\\(\\nu\\), by induction on the cost: an edit takes the furthest of three points, a run of matching bases from a '
 'point goes on from any later point on it, and a further point has crossed at least as much potential, so it is '
 'dropped only if the nearer one is too. Hence if the search from \\(\\nu\\) crosses the last of the seeds, so does '
 'that from \\(\\mu\\). \\(\\square\\)'),
('p',
 '**Proposition 1.** Let \\(C\\) be the least cost of an alignment within \\([d_\\ell, d_h]\\), and let \\(L(d) = '
 'g(|d|) + g(|d-\\delta|)\\), or \\(\\infty\\) for a diagonal outside the matrix. If \\(C \\le \\min\\bigl(L(d_h+1), '
 'L(d_\\ell-1)\\bigr)\\), then \\(C\\) is the optimal cost; if the inequality is strict, every optimal alignment lies '
 'within the band.'),
('p',
 '*Proof.* An alignment that leaves the band visits diagonal \\(d_h+1\\) or \\(d_\\ell-1\\). Consider \\(d=d_h+1 > '
 '\\max(0,\\delta)\\); the case below the band is symmetric. Starting on diagonal 0, the alignment contains '
 'insertions totalling at least \\(d\\) letters to reach \\(d\\), and deletions totalling at least \\(d-\\delta\\) '
 'letters to end on \\(\\delta\\). Since \\(g\\) is nondecreasing and subadditive for \\(o, e \\ge 0\\), as is the '
 'minimum of two such functions for two-piece costs, gaps totalling at least \\(k\\) letters cost at least '
 '\\(g(k)\\), so their total cost is at least \\(g(d)+g(d-\\delta)=L(d)\\). Every alignment outside the band '
 'therefore costs at least \\(\\min\\bigl(L(d_h+1), L(d_\\ell-1)\\bigr)\\), so \\(C\\) at most this minimum is optimal, '
 'and \\(C\\) below it leaves every optimal alignment inside the band. \\(\\square\\)'),
('p',
 '**Proposition 2.** Suppose the first search stops at cost \\(\\sigma\\) without reaching the far corner, and the '
 'second search, growing from cost 0, checks each new front against the last \\(w+1\\) fronts of the first, costs '
 '\\(\\sigma-w\\) to \\(\\sigma\\), while its own cost is at most \\(\\sigma\\), after which the two alternate as in '
 "BiWFA, each new front checked against the other side's last \\(w+1\\). When they stop, once \\(2\\min(s_1, s_2) "
 '\\ge C^* - 1 + \\theta\\), \\(C^*\\) is the optimal cost and its meeting lies on an optimal path.'),
('p',
 '*Proof.* Let \\(P\\) be an optimal path, of cost \\(C_{\\mathrm{opt}}\\), and suppose the searches stop with \\(C^* '
 '> C_{\\mathrm{opt}}\\). As costs are integers, \\(C^* - 1 \\ge C_{\\mathrm{opt}}\\), so both searches have passed '
 "\\((C_{\\mathrm{opt}}+\\theta)/2\\) when they stop. Write \\(u_1\\) and \\(u_2\\) for a cell's costs from the starts "
 'of the first and the second search. Along \\(P\\), \\(u_1\\) never decreases and \\(u_2\\) never increases, \\(u_1 + '
 'u_2 = C_{\\mathrm{opt}}\\) between moves and \\(u_1 + u_2 \\le C_{\\mathrm{opt}} + o\\) inside a gap, and one move '
 'changes each by at most \\(w\\). Having checked at each cost up to \\(\\sigma\\) that it had not reached the far '
 "corner, which is meeting the second search's front of cost 0, the first search shows that \\(C_{\\mathrm{opt}} > "
 '\\sigma\\), so \\(P\\) has a cell with \\(\\sigma - w < u_1 \\le \\sigma\\). (i) If \\(u_2 < u_1 - w\\) there, then \\(u_2 '
 '< \\sigma - w\\) and \\(2u_2 < C_{\\mathrm{opt}} + o - w < C_{\\mathrm{opt}} + \\theta\\), so the second search '
 'reaches the cell in its first phase and checks it against cost \\(u_1\\), which the first search holds. (ii) '
 'Otherwise, since \\(u_2 - u_1\\) falls to \\(-C_{\\mathrm{opt}}\\) at the far corner in steps of at most \\(2w\\), '
 '\\(P\\) has a cell at or after this one with \\(|u_1 - u_2| \\le w\\) and \\(u_1 > \\sigma - w\\). Both of its costs '
 'are at most \\((C_{\\mathrm{opt}} + o + w)/2 = (C_{\\mathrm{opt}} + \\theta)/2\\), so both searches reach it; it is '
 'checked in the first phase if both costs are at most \\(\\sigma\\), and otherwise in the alternation, where each '
 "new front meets the other side's last \\(w + 1\\). Either way a meeting of cost \\(C_{\\mathrm{opt}}\\) is found, "
 'contradicting \\(C^* > C_{\\mathrm{opt}}\\). \\(\\square\\)'),
('p',
 "**Lemma S1.** Let \\(M\\) hold every match of a seed with cost below its potential \\(r\\), and let \\(M'\\) be "
 '\\(M\\) without the matches whose local-pruning search (Section 2.4 of the main text) finds no path across the '
 "next \\(p\\) seeds. The gap-chaining seed heuristic over \\(M'\\) does not exceed the cost of any path from its "
 'cell to the end of the pair.'),
('p',
 '*Proof.* The heuristic charges a chain of matches the cost of each match and, for each gap between consecutive '
 'matches and from the last to the end, the larger of the potential of the seeds that start in the gap and the '
 'number of diagonals the gap shifts (Groot Koerkamp and Ivanov, 2024). Let \\(\\pi\\) be a path from cell \\(u\\) to '
 'the end, and sweep it from left to right from a column \\(x\\), at first that of \\(u\\). Let \\(s\\) be the next '
 'seed that starts at or after \\(x\\); no seed starts between \\(x\\) and \\(s\\). If the piece of \\(\\pi\\) across '
 '\\(s\\) costs at least \\(r\\), it pays the potential of \\(s\\), and the sweep continues from the end of \\(s\\). '
 "Otherwise the crossing is a match \\(\\mu \\in M\\). If \\(\\mu \\in M'\\), it joins the chain and the sweep continues "
 'from its end. If not, \\(\\pi\\) continued from \\(\\mathrm{end}(\\mu)\\) at cost \\(c(\\mu)\\) is a path that the '
 'search did not find, so its cost reaches the potential of the seeds it has crossed at some column \\(y\\) no '
 'later than the end of the \\(p\\)-th seed: the piece of \\(\\pi\\) from the start of \\(s\\) to \\(y\\) pays the '
 'potential of the seeds that start in it, and the sweep continues from \\(y\\). The matches joined appear along '
 "\\(\\pi\\) in order, so they form a chain in \\(M'\\). On each joined match \\(\\pi\\) pays its cost. Each gap is cut "
 'into pieces, each of which pays at least the potential of the seeds that start in it and at least the number '
 'of diagonals it shifts. Summed over the pieces, \\(\\pi\\) pays in the gap at least the potential of the seeds '
 'that start in it and at least the diagonals it shifts, so at least the larger of the two. Hence \\(\\pi\\) costs '
 'at least what the heuristic charges this chain, and so at least the heuristic. \\(\\square\\)'),
('p',
 'Besides local pruning, the implementation drops only matches whose gap to the end of the pair costs at least '
 'what they save, which no chain can use to lower the heuristic; a seed left uncounted, as in a long tandem '
 "repeat, has potential 0, which the proof covers. It keeps every other match of \\(M'\\) and some that local "
 'pruning would drop (Section 3.5 of the main text); a heuristic over more matches is no larger, so it too never '
 'exceeds the remaining cost. The cell-by-cell check of Section S1.2 runs with local pruning enabled.'),
]

TABLES = {
"ablation": {
  "label": "S2",
  "fit": True,
  "caption": "Effect of each technique: median time with all techniques on (ms per alignment; per batch of 500,000 pairs for the short reads), and the change when one technique is switched off, over three alternating rounds on the Skylake-X at commit `dc980d7`, with branches kept clear of 32-byte boundaries. Every configuration returned the same costs.",
  "header": ["Dataset", "All on", "Regrouped recurrence", "Interleaved groups", "Gathered extension", "Agreeing estimates", "Direct neighbours", "Certified bands"],
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
    ["Short reads, gap-affine", "232", "", "", "", "", "", "**+188%**"],
    ["Short reads, unit", "213", "", "", "", "", "", "**+215%**"],
  ]},
"modes": {
  "label": "S3",
  "caption": "Other modes: mean time per pair in ms on one pinned thread, the median of three runs, each the faster of two passes. Every tool's costs or scores agree on every pair.",
  "header": ["Workload", "dinara-align", "Next fastest"],
  "align": "lrl",
  "wrap": [0],
  "rows": [
    ["Infix at unit costs, 1 kbp read at 10% in a 3 kbp window", "**0.130**", "Edlib 0.268"],
    ["Infix at unit costs, 10 kbp read at 5% in a 30 kbp window", "**2.01**", "Edlib 13.5"],
    ["Prefix at unit costs, 1 kbp read against 2 kbp", "**0.0424**", "WFA2-lib 0.0470"],
    ["Infix at (4, 6, 2), 1 kbp read in a 3 kbp window", "**3.92**", "WFA2-lib 3.96"],
    ["Extension, match 2 at (4, 6, 2)", "**0.726**", "KSW2 2.06"],
    ["Extension with an end bonus of 50", "**0.166**", "KSW2 1.24"],
    ["Two-piece gaps (4, 6, 2, 24, 1), 5 kbp at 5%", "**8.25**", "WFA2-lib 15.8"],
    ["Table, global, 1 kbp at 10%", "**0.196**", "parasail 2.94"],
    ["Table, local, 1 kbp read in a 10 kbp window", "**1.91**", "SSW 3.28"],
    ["Local, 600 short noisy pairs (20–700 bp)", "**0.0320**", "SSW 0.0512"],
    ["Local, 1 kbp read at 10% in a 10 kbp window", "**1.36**", "SSW 3.15"],
    ["Local, 10 kbp at 5% against 12 kbp", "**17.4**", "SSW 84.7"],
    ["Overlap of 2 kbp reads by 0.5–1.5 kbp", "**0.703**", "parasail 15.3"],
    ["Infix with a reward, 1 kbp read in a 3 kbp window", "**0.687**", "parasail 15.1"],
  ]},
"affine_full": {
  "label": "S4",
  "fit": True,
  "caption": "Gap-affine alignment at costs (4, 6, 2), all datasets: mean time per alignment in ms, with a 30 s budget per sample, the median of three alternating rounds. Brackets give the pairs completed where not all were; t.o. marks a tool that completed none within 30 s and o.o.m. one the 8 GiB limit stopped. Bold marks the lowest mean of the tools that completed the most pairs.",
  "header": ["Dataset", "dinara-align", "WFA", "BiWFA", "BSAlign", "KSW2"],
  "align": "lrrrrr",
  "rows": [
    ["*Real data*", "", "", "", "", ""],
    ["ont-1k", "**0.105**", "0.189", "0.305", "0.336", "1.10"],
    ["ont-10k", "**1.59**", "4.86", "5.85", "7.83", "29.9"],
    ["ont-50k", "**21.7**", "63.4", "54.9", "86.8", "341 [81/104]"],
    ["ont-promethion", "**286**", "286 [7/22]", "576", "o.o.m.", "2,350 [7/22]"],
    ["ont-500k", "t.o.", "o.o.m.", "t.o.", "o.o.m.", "o.o.m."],
    ["ont-500k-genvar", "27,300 [1/4]", "o.o.m.", "t.o.", "o.o.m.", "o.o.m."],
    ["SARS-CoV-2", "**0.698**", "2.06", "1.78", "341", "88.6"],
    ["Illumina", "**0.00151**", "0.00154", "0.00282", "0.0538", "0.0128"],
    ["PacBio HiFi", "**0.0929**", "0.224", "0.232", "86.1", "9.74"],
    ["*Simulated: random sequences, uniform errors (QuickEd)*", "", "", "", "", ""],
    ["10 kbp, 1%", "**0.195**", "0.409", "0.649", "37.0", "21.2"],
    ["10 kbp, 5%", "**2.37**", "7.13", "8.73", "37.2", "97.6"],
    ["10 kbp, 10%", "**7.08**", "22.8", "27.2", "37.2", "204"],
    ["10 kbp, 20%", "**19.9**", "65.5", "75.2", "37.5", "375 [80/99]"],
    ["100 kbp, 1%", "**11.9**", "39.9", "39.0", "o.o.m.", "2,440"],
    ["100 kbp, 5%", "**357**", "810", "715", "o.o.m.", "9,220 [3/9]"],
    ["100 kbp, 10%", "**1,370**", "o.o.m.", "2,430", "o.o.m.", "17,400 [1/10]"],
    ["100 kbp, 20%", "**3,940 [7/10]**", "o.o.m.", "7,350 [4/10]", "o.o.m.", "o.o.m."],
    ["500 kbp, 1%", "**499**", "1,290", "942", "o.o.m.", "o.o.m."],
    ["500 kbp, 5%", "**11,400 [2/4]**", "o.o.m.", "19,400 [1/4]", "o.o.m.", "o.o.m."],
    ["500 kbp, 10%", "t.o.", "o.o.m.", "t.o.", "o.o.m.", "o.o.m."],
    ["500 kbp, 20%", "t.o.", "o.o.m.", "t.o.", "o.o.m.", "o.o.m."],
    ["1 Mbp, 1%", "**2,440**", "o.o.m.", "3,920", "o.o.m.", "o.o.m."],
    ["1 Mbp, 5%", "t.o.", "o.o.m.", "t.o.", "o.o.m.", "o.o.m."],
    ["1 Mbp, 10%", "t.o.", "o.o.m.", "t.o.", "o.o.m.", "o.o.m."],
    ["1 Mbp, 20%", "t.o.", "o.o.m.", "t.o.", "o.o.m.", "o.o.m."],
  ]},
"unit_full": {
  "label": "S5",
  "fit": True,
  "caption": "Unit costs, all datasets: mean time per alignment in ms, with the median in parentheses, traceback included, a 30 s budget per sample, the median of three alternating rounds. Brackets give the pairs completed where not all were, the median over those; t.o. marks a tool that completed none within 30 s and o.o.m. one the 8 GiB limit stopped. Bold marks the lowest mean of the tools that completed the most pairs.",
  "header": ["Dataset", "dinara-align", "A*PA2-full", "A*PA2-simple", "A*PA", "Edlib", "BiWFA", "WFA", "QuickEd"],
  "align": "lrrrrrrrr",
  "rows": [
    ["*Real data*", "", "", "", "", "", "", "", ""],
    ["ont-1k", "**0.0208 (0.0212)**", "0.0810 (0.0897)", "0.0481 (0.0542)", "0.613 (0.557)", "0.124 (0.117)", "0.0442 (0.0461)", "0.0313 (0.0322)", "0.0732 (0.0723)"],
    ["ont-10k", "**0.154 (0.105)**", "0.393 (0.267)", "0.247 (0.179)", "15.4 (3.30)", "1.11 (0.576)", "0.752 (0.283)", "0.599 (0.224)", "0.679 (0.449)"],
    ["ont-50k", "**0.582 (0.221)**", "1.28 (0.515)", "0.941 (0.305)", "263 (8.99)", "6.11 (2.02)", "7.43 (0.872)", "6.35 (0.739)", "2.31 (1.05)"],
    ["ont-promethion", "**3.99 (1.99)**", "7.02 (4.12)", "6.76 (2.94)", "929 (31.8)", "73.9 (26.7)", "71.1 (19.0)", "79.7 (15.5)", "17.1 (9.16)"],
    ["ont-500k", "**118 (73.9)**", "284 (111)", "688 (399)", "t.o.", "5,000 (3,040)", "6,370 (2,230)", "o.o.m.", "848 (586)"],
    ["ont-500k-genvar", "**156 (156)**", "225 (199)", "512 (565)", "t.o.", "5,230 (5,200)", "6,180 (6,600)", "o.o.m.", "1,390 (1,320)"],
    ["SARS-CoV-2", "**0.277 (0.130)**", "2.01 (2.08)", "0.727 (0.664)", "6.55 (2.08)", "8.18 (7.77)", "0.906 (0.453)", "0.979 (0.367)", "3.86 (4.31)"],
    ["Illumina", "**0.000593 (0.000496)**", "0.0145 (0.0143)", "0.00411 (0.00406)", "0.0137 (0.0132)", "0.0234 (0.0233)", "0.00206 (0.00182)", "0.00112 (0.000987)", "0.0162 (0.0154)"],
    ["PacBio HiFi", "**0.0129 (0.00533)**", "0.980 (0.929)", "0.259 (0.239)", "1.06 (0.966)", "3.41 (3.19)", "0.0865 (0.0686)", "0.0660 (0.0474)", "0.951 (0.883)"],
    ["*Simulated: random sequences, uniform errors (QuickEd)*", "", "", "", "", "", "", "", ""],
    ["10 kbp, 1%", "**0.0258 (0.0228)**", "0.738 (0.733)", "0.160 (0.156)", "0.566 (0.562)", "3.25 (3.25)", "0.0923 (0.0882)", "0.0623 (0.0576)", "0.614 (0.614)"],
    ["10 kbp, 5%", "**0.268 (0.242)**", "0.795 (0.790)", "0.288 (0.284)", "1.08 (1.07)", "3.81 (3.82)", "0.637 (0.629)", "0.513 (0.499)", "1.26 (1.25)"],
    ["10 kbp, 10%", "**0.348 (0.346)**", "0.724 (0.712)", "0.444 (0.438)", "3.23 (3.21)", "4.57 (4.56)", "2.26 (2.26)", "1.84 (1.67)", "2.01 (2.03)"],
    ["10 kbp, 20%", "**0.511 (0.467)**", "1.44 (1.44)", "1.01 (1.01)", "204 (202)", "5.96 (5.96)", "7.31 (7.31)", "5.34 (5.11)", "1.54 (1.57)"],
    ["100 kbp, 1%", "**1.47 (1.13)**", "7.53 (7.24)", "3.26 (3.15)", "7.76 (7.47)", "58.1 (58.0)", "3.72 (3.45)", "4.08 (4.32)", "12.5 (12.9)"],
    ["100 kbp, 5%", "**5.10 (4.82)**", "8.35 (8.06)", "14.9 (14.7)", "14.1 (14.4)", "147 (147)", "63.3 (63.2)", "54.9 (50.9)", "23.9 (23.7)"],
    ["100 kbp, 10%", "**7.24 (7.04)**", "13.9 (13.2)", "27.7 (27.7)", "48.4 (46.6)", "247 (247)", "218 (218)", "182 (169)", "39.0 (38.9)"],
    ["100 kbp, 20%", "**15.8 (15.6)**", "36.8 (36.2)", "56.8 (56.7)", "t.o.", "445 (446)", "709 (709)", "553 (510)", "62.6 (62.6)"],
    ["500 kbp, 1%", "**33.0 (31.6)**", "43.1 (41.7)", "72.5 (72.4)", "41.2 (39.6)", "669 (670)", "82.1 (81.4)", "75.8 (66.1)", "124 (124)"],
    ["500 kbp, 5%", "**31.2 (29.9)**", "44.7 (43.9)", "213 (212)", "81.7 (80.7)", "2,290 (2,290)", "1,550 (1,550)", "1,480 (1,250)", "394 (394)"],
    ["500 kbp, 10%", "**67.8 (63.8)**", "119 (118)", "423 (423)", "407 (406)", "4,380 (4,380)", "5,450 (5,450)", "5,020 (4,240)", "699 (699)"],
    ["500 kbp, 20%", "**253 (252)**", "544 (543)", "898 (893)", "t.o.", "8,530 [3/4] (8,530)", "18,000 [1/4] (18,000)", "o.o.m.", "1,240 (1,240)"],
    ["1 Mbp, 1%", "**61.1 (60.3)**", "89.2 (87.8)", "253 (253)", "96.6 (93.1)", "2,370 (2,380)", "365 (364)", "392 (352)", "393 (393)"],
    ["1 Mbp, 5%", "**65.0 (62.9)**", "116 (114)", "819 (816)", "183 (185)", "8,830 [3/4] (8,830)", "6,350 (6,350)", "o.o.m.", "1,440 (1,440)"],
    ["1 Mbp, 10%", "**185 (184)**", "368 (365)", "1,650 (1,640)", "974 (971)", "17,300 [1/4] (17,300)", "22,300 [1/4] (22,300)", "o.o.m.", "2,650 (2,650)"],
    ["1 Mbp, 20%", "**890 (687)**", "2,020 (2,010)", "3,510 (3,490)", "t.o.", "t.o.", "t.o.", "o.o.m.", "4,800 (4,800)"],
  ]},
"datasets": {
  "label": "S6",
  "caption": "The samples: pairs in each sample of the dataset, their lengths and their divergence, the edit distance over the longer sequence, as QuickEd's Table S1 describes its datasets.",
  "header": ["Dataset", "Pairs", "Length, min / mean / max (bp)", "Divergence, min / mean / max (%)"],
  "align": "lrrr",
  "rows": [
    ["*Real data*", "", "", ""],
    ["ont-1k", "1221 of 12477", "42 / 830 / 1,089", "0 / 11 / 21"],
    ["ont-10k", "277 of 5000", "342 / 3,643 / 9,791", "4 / 13 / 20"],
    ["ont-50k", "104 of 10000", "360 / 9,633 / 42,626", "4 / 12 / 19"],
    ["ont-promethion", "22 of 584", "10,362 / 44,631 / 155,664", "4 / 10 / 20"],
    ["ont-500k", "4 of 50", "553,834 / 642,465 / 683,520", "4 / 7 / 13"],
    ["ont-500k-genvar", "4 of 48", "594,436 / 662,980 / 861,994", "5 / 7 / 10"],
    ["SARS-CoV-2", "33 of 10000", "29,244 / 29,746 / 29,864", "0 / 1 / 7"],
    ["Illumina", "4032 of 5992", "77 / 248 / 251", "0 / 0 / 10"],
    ["PacBio HiFi", "69 of 890", "10,956 / 14,372 / 22,816", "0 / 0 / 1"],
    ["*Simulated: random sequences, uniform errors (QuickEd)*", "", "", ""],
    ["10 kbp, 1%", "100 of 10000", "10,000 / 10,003 / 10,023", "1 / 1 / 1"],
    ["10 kbp, 5%", "99 of 10000", "10,000 / 10,008 / 10,045", "5 / 5 / 5"],
    ["10 kbp, 10%", "100 of 10000", "10,000 / 10,008 / 10,056", "9 / 9 / 9"],
    ["10 kbp, 20%", "99 of 10000", "10,000 / 10,014 / 10,095", "17 / 17 / 17"],
    ["100 kbp, 1%", "9 of 1000", "100,000 / 100,013 / 100,035", "1 / 1 / 1"],
    ["100 kbp, 5%", "9 of 1000", "100,000 / 100,028 / 100,123", "5 / 5 / 5"],
    ["100 kbp, 10%", "10 of 1000", "100,000 / 100,039 / 100,120", "9 / 9 / 9"],
    ["100 kbp, 20%", "10 of 1000", "100,000 / 100,040 / 100,155", "17 / 17 / 17"],
    ["500 kbp, 1%", "4 of 200", "500,000 / 500,038 / 500,103", "1 / 1 / 1"],
    ["500 kbp, 5%", "4 of 200", "500,016 / 500,088 / 500,178", "5 / 5 / 5"],
    ["500 kbp, 10%", "4 of 200", "500,000 / 500,052 / 500,105", "9 / 9 / 9"],
    ["500 kbp, 20%", "4 of 200", "500,000 / 500,055 / 500,220", "17 / 17 / 17"],
    ["1 Mbp, 1%", "4 of 100", "1,000,000 / 1,000,034 / 1,000,107", "1 / 1 / 1"],
    ["1 Mbp, 5%", "4 of 100", "1,000,185 / 1,000,301 / 1,000,551", "5 / 5 / 5"],
    ["1 Mbp, 10%", "4 of 100", "1,000,000 / 1,000,054 / 1,000,120", "9 / 9 / 9"],
    ["1 Mbp, 20%", "4 of 100", "1,000,000 / 1,000,166 / 1,000,355", "17 / 17 / 17"],
  ]},
"memory": {
  "label": "S7",
  "caption": "Memory at unit costs: the largest growth of peak resident memory aligning one pair of the sample, in MiB, as in A*PA2, and the peak resident memory of the whole run in parentheses, as in QuickEd, runtime and input included. Brackets give the pairs completed where not all were, with no whole-run peak; t.o. marks a tool that completed no pair within 30 s.",
  "header": ["Dataset", "dinara-align", "A*PA2-full", "A*PA2-simple", "Edlib", "BiWFA", "QuickEd"],
  "align": "lrrrrrr",
  "rows": [
    ["*Real data*", "", "", "", "", "", ""],
    ["ont-500k", "48 (64)", "81 (89)", "82 (90)", "5.3 (15)", "15 (29)", "22 (26)"],
    ["ont-500k-genvar", "39 (64)", "49 (56)", "32 (55)", "4.6 (17)", "12 (27)", "19 (27)"],
    ["*Simulated: random sequences, uniform errors (QuickEd)*", "", "", "", "", "", ""],
    ["100 kbp, 20%", "5.3 (19)", "7.2 (12)", "4.6 (9)", "1.4 (8)", "3.7 (11)", "8.8 (12)"],
    ["500 kbp, 20%", "53 (82)", "59 (66)", "53 (61)", "4.0 [3/4]", "15 [1/4]", "15 (19)"],
    ["1 Mbp, 20%", "123 (178)", "184 (196)", "192 (205)", "t.o.", "t.o.", "26 (31)"],
  ]},
"batch": {
  "label": "S8",
  "caption": "Global scores of 500,000 short-read pairs: each time is the whole call, in ms, the median of three runs.",
  "header": ["Costs", "Threads", "dinara-align", "WFA2-lib", "KSW2", "parasail", "Edlib"],
  "align": "llrrrrr",
  "rows": [
    ["Affine (1, 2, 1)", "1", "**231**", "1,060", "10,900", "21,700", "–"],
    ["Unit", "1", "**212**", "450", "–", "–", "2,000"],
    ["Affine (1, 2, 1)", "10", "**29.0**", "113", "1,200", "2,340", "–"],
    ["Unit", "10", "**28.0**", "49.0", "–", "–", "209"],
  ]},
"changes": {
  "label": "S1",
  "caption": "Changes behind the results and their measured effects, as recorded when each was made. Each row times what its commit timed (an in-process loop through one aligner, the benchmark runner, a distance or an alignment), so rows for one dataset need not agree with each other or with Section 4 of the main text. Machines: the i9-7900X (Skylake-X, AVX-512) and the Apple M2 (NEON); “not recorded” where the commit does not say.",
  "header": _changes[2],
  "align": "llll",
  "wrap": [1, 2],
  "rows": _changes[3],
},
}

REFS = {}
BIBKEYS = {}
