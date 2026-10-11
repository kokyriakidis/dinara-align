"""The supplementary material, in the same markup as main_text.py."""
from change_record import SUPPLEMENT as _OLD

TITLE = "Supplementary material for “dinara-align: hardware-aware algorithms for exact pairwise alignment”"

# Table S1's rows are kept with the change record they came from; only labels that read informally are reworded.
_changes = next(b for b in _OLD if b[0] == "table")
_LABELS = {"Prefix search gives up a doomed try": "Prefix search abandons a try that cannot succeed"}
for _row in _changes[3]:
    _row[1] = _LABELS.get(_row[1], _row[1])
# The changes since the change record, each measured on the Skylake-X when it was made.
_changes[3].extend([
    ["c5aee28", "Table global alignment traced in bands that prove themselves", "table-global 544 → about 200 µs", "Skylake-X"],
    ["691f43a", "Divergent gap-affine pairs handed to a full sweep in difference recurrences", "unrelated 1, 2, 5 kbp at (4, 6, 2) 1.50, 5.3, 32.9 → 0.28, 0.87, 4.8 ms", "Skylake-X"],
    ["64be0ec", "Anti-diagonal sweeps a 16 KB strip of rows at a time", "local, 10 kbp against 12 kbp 50.0 → 17.4 ms", "Skylake-X"],
    ["ffbf46e", "Banded traceback, local end search and linear-space halves in place", "divergent table pairs 18 to 37% faster", "Skylake-X"],
    ["94bf719", "Unit-cost traceback within the memory limit, stretches swept again", "1 Mbp at 20%: peak growth 175 → 133 MB a pair; 3.35 → 3.64 s for four", "Skylake-X"],
    ["289bcf0", "Gap-affine handover only after a 24th of the sweep's work", "ont-10k at (4, 6, 2) 2.33 → 1.55 ms; ont-50k 30.5 → 21.0 ms", "Skylake-X"],
    ["ea3f09e", "256-column tiles for every band", "10 kbp at 5% 318 → 254 µs; ont-10k 161 → 147 µs", "Skylake-X"],
    ["dc980d7", "Inexact seeds from 128 kbp on AVX-512", "100 kbp at 10% 10.3 → 6.8 ms", "Skylake-X"],
    ["7aab823", "A space holds its profile directly, not in an Optional, and an aligner's traceback reuses it", "250 bp Illumina, through one aligner 0.97 → 0.80 µs; distance 1.21 → 0.68 µs", "Skylake-X"],
    ["90e0040", "An aligner keeps the fronts it traces back", "250 bp Illumina 0.82 → 0.65 µs", "Skylake-X"],
    ["32e78aa", "A gap-affine alignment begins with the one search it traces", "250 bp Illumina at (4, 6, 2) 2.41 → 1.52 µs", "Skylake-X"],
    ["7d38a99", "A reversed pair coded back to front in the coding pass", "250 bp Illumina 0.649 → 0.583 µs", "Skylake-X"],
    ["4558186", "Alignments pass between internal calls without an Optional", "250 bp Illumina 0.585 → 0.420 µs; at (4, 6, 2) 0.995 → 0.840 µs", "Skylake-X"],
    ["036e12b", "A profile coded 32 bytes at a time", "250 bp Illumina 0.421 → 0.400 µs; PacBio HiFi 5.7 → 4.95 µs", "Skylake-X"],
    ["a6d9166", "The first gap-affine search grows alone while it steps no more than two searches would", "250 bp Illumina at (4, 6, 2) 1.26 → 1.20 µs; pairs costing 9 to 32 a third less", "Skylake-X"],
])

BODY = [
("p", "This supplement provides implementation details, per-change measurements and additional results that the main "
 "text summarises. Section S1 describes techniques not covered in the main text, Section S2 lists the measured "
 "effect of each change, Section S3 gives the complete results for gap-affine costs and the other alignment modes, "
 "and Sections S4 and S5 describe robustness measures and approaches that were evaluated and not adopted. Unless "
 "stated otherwise, the measurements in Sections S1, S2, S4 and S5 were taken when each change was introduced, on "
 "the machine indicated, and compare the change with the preceding version; the results in the main text and in "
 "Section S3 were obtained on the Skylake-X at commit `036e12b`, gap-affine costs at `a6d9166`, and the other "
 "alignment modes and batches at `289bcf0`. Two machines were used: an Intel Core i9-7900X (Skylake-X, AVX-512) and an "
 "Apple M2 (NEON)."),

("sec", "Implementation details"),
("sub", "Unit costs"),
("p", "**Engine selection.** The diagonal transition receives a step budget proportional to the cost of the band it "
 "would replace; on the Skylake-X a step costs about 1.3 ns and a band column about 6 ns. Inexact seeds replace exact "
 "seeds when fewer than 20, 30 or 40% of exact seeds chain on AVX-512, AVX2 and NEON, respectively: the wider the "
 "vector unit, the cheaper the bit-parallel band relative to the scalar construction of inexact seeds. All "
 "thresholds are fixed at compile time for each instruction set, so that no dispatch occurs at run time."),
("p", "**Divergent pairs at gap-affine costs.** The wavefront's work grows with the cost times the length, so on "
 "unrelated 2 kbp pairs it took 4.6 ms where the 8-bit difference sweep of BSAlign took 1.5 ms. A global pair with "
 "the whole matrix open is therefore handed to a full sweep in difference recurrences (Suzuki and Kasahara, 2018), in 8-bit lanes "
 "at any length, once its projected work exceeds the sweep's, at seven sweep cells per wavefront step, the rate "
 "measured at 10 kbp. The sweep keeps each cell's predecessor flag under the tie rule of the wavefront, so the "
 "alignment is identical byte for byte, and for the traceback it keeps every \\(k\\)-th anti-diagonal, \\(k\\) the "
 "square root of the anti-diagonals times the bytes a cell takes, and sweeps each stretch again. Each anti-diagonal "
 "is swept in place from the last row up. On unrelated pairs of 1, 2 and 5 kbp at (4, 6, 2) alignment took 0.28, "
 "0.87 and 4.8 ms instead of 1.50, 5.3 and 32.9 ms. The projection extrapolates the cost from how far the two "
 "wavefronts have advanced, which on nanopore reads, whose errors gather at the ends, ran up to 2.8 times too high "
 "in cost and 8 times in work, and handed reads to sweeps that took up to eight times longer than their wavefront: "
 "the 10 kbp reads took 2.33 ms instead of 1.55. As for the diagonal transition, the projection is now trusted only "
 "once the wavefront's work has reached one twenty-fourth of the sweep's."),
("p", "**Register and memory layout of the sweep.** An empty inline-assembly barrier prevents the compiler from "
 "reassociating the regrouped recurrence back to its original form. Vertical differences remain in registers across "
 "a word's columns, and the horizontal differences of successive columns are stored at distinct addresses, so that "
 "no load waits on a preceding store. Bands of five to seven words combine a four-lane vector with scalar words."),
("p", "**Symbols other than ACGT.** `N` and other non-ACGT symbols are encoded in a third bit plane. Each tile uses "
 "only the planes its symbols require: two when neither its columns nor any row contains such a symbol, the third "
 "plane of the rows as a mask when only rows contain one, and all three otherwise. Sequences are encoded 16 bytes or "
 "eight rows at a time, and only chunks that contain such a symbol are processed byte by byte. This made reads "
 "containing `N` 7 to 9% faster on both machines. A seed that contains `N` is not counted, so that pairs containing "
 "`N` can still use seeds; this made the ultra-long reads 2.5 to 3 times faster on the M2 and 1.4 to 1.7 times faster "
 "on the Skylake-X."),
("p", "**Tile width.** Tiles are 256 columns wide, as in A*PA2, whatever the height of the band. Earlier versions "
 "used tiles of 64 columns in bands bounded below 1,024, a rule kept from the first version on an argument never "
 "measured: narrower tiles cut the slope a band adds across each tile, but each tile is pruned, recorded and traced "
 "back on its own. With one width, 10 kbp pairs at 5% divergence aligned in 254 µs instead of 318, the 10 kbp "
 "nanopore reads in 147 instead of 161 µs, and every other dataset but one in the same time or up to 6.5% faster; "
 "10 kbp pairs at 1%, settled mostly by the diagonal transition, were within 3%."),
("p", "**Traceback within a memory limit.** The traceback reads the left boundary of every tile, which the band "
 "records: on 1 Mbp pairs at 20% errors, up to 175 MB. The band keeps every boundary while they fit within the "
 "memory limit of the alignment, 80 MB by default, and otherwise only every \\(k\\)-th, \\(k\\) doubling as the "
 "round proceeds, leaving room for one stretch of the widest boundaries. Each tile also records the state the band "
 "carries into it: its words, the score at its top, the deepest kept row and its score, and the bound. The "
 "traceback, from the last tile back, sweeps a stretch again from the last kept boundary before a tile it lacks; "
 "every step of the sweep depends only on that state and the scores it computes, so each boundary is recovered "
 "exactly and the alignment is the same byte for byte, unlike the Hirschberg split of QuickEd (Doblas et al., 2025), which breaks ties "
 "within each half. Kept boundaries form a stack in tile order and are released as they are traced. The boundary "
 "buffers grow by doubling up to half the limit and then to the whole of it at once, since doubling further would "
 "hold the old buffer and its copy together, up to twice the limit. On the four 1 Mbp pairs at 20% errors, measured "
 "one per process, an alignment grew the peak memory by 133, 98, 124 and 113 MB instead of 175, 99, 143 and 146 MB, "
 "in 3.64 s instead of 3.35 s for the four, and pairs whose boundaries fit were unchanged. With no memory at all, "
 "alignments took about twice the time."),
("p", "**Accepting an extrapolation.** The diagonal-transition estimate is abandoned once it has spent one "
 "twenty-fourth of the cost of the cheapest band. An extrapolation sets the band only when the single-wavefront "
 "extrapolation from the first edits and the extrapolation of the longer search agree to within one third; "
 "otherwise the first bound is set 256 above the known lower bound, and each retry at most doubles it. A later "
 "checkpoint never lowers the bound below the previous bound plus one quarter, and a failed round continues from the "
 "bound at which it ended rather than the one at which it started. Before these rules, a 30 kbp SARS-CoV-2 pair at "
 "distance 49 had been extrapolated to 14,448 and computed over the full matrix."),
("p", "When they were introduced, these rules reduced the mean time on the ultra-long reads with genetic variation "
 "from 4.33 s, with 13 of 15 reads completed within the time limit, to 754 ms. At the evaluated commit the seed "
 "heuristic sets the first bound on long reads, and switching the agreement test off changes the ultra-long reads "
 "by less than 1% but slows the 10 kbp reads by 22% (main text, Effect of each technique)."),
("p", "**Retries with seeds.** When seeds are used, a retry sets the bound one quarter of the extrapolated increase "
 "beyond the estimate; with one half, a round had failed after 90% of the columns of a 527 kbp read. A seeded round "
 "that fails within the first eighth of its columns does not extrapolate; instead, its margin over the heuristic "
 "bound at the origin increases fourfold, as in A*PA2. Previously, an 884 kbp read had failed at column 1,280 and "
 "been extrapolated to 125,000 for a distance of 75,000. Only the first round revises its bound at the "
 "checkpoints; when later rounds also lowered their bounds, they failed late, four times on one 897 kbp read."),
("p", "**Traceback.** A tile on which the forward search exceeds its limit is recomputed over a window of rows above "
 "the traced cell, starting at four words and doubling. The top boundary of the window is initialised to +1, so that "
 "its scores correspond to real paths. A 600 kbp nanopore pair at 13% divergence was aligned in 617 ms instead of "
 "1.38 s, of which 840 ms had been spent recomputing tiles of 900 words. Recomputation buffers are reused without "
 "initialisation, and tile boundaries are allocated once and written through pointers; loading the boundaries had "
 "accounted for a quarter of the traceback time. The CIGAR string is produced by reversing moves 16 at a time and "
 "splitting diagonal runs eight bases at a time, which costs 1 to 3% of the alignment time."),
("p", "**Diagonal transition.** Wavefronts are padded with unreached diagonals and advanced eight diagonals at a time "
 "with branch-free validity masks, in a separate function; a step costs about 1.3 ns, compared with 5 ns before. The "
 "two-ended search tests the overlap of the two wavefronts with vector comparisons and requires about half as many "
 "steps. Since a step costs about 1.3 ns and a band column about 6 ns, the diagonal transition is given a step budget "
 "proportional to the band it would replace; reducing this budget from 12 to 5 steps per column reduced the time "
 "for 10 kbp pairs at 5% divergence from 241 to 144 µs. The choice between one and two wavefronts was fitted on "
 "1,300 pairs of 1 to 30 kbp. For pairs of up to 2,048 columns, the band had been set wide across the whole matrix, "
 "and 42% of such pairs were assigned to a band although the diagonal transition would have been faster; such a "
 "band is now charged an additional 30,000 steps, which made the 1 kbp reads 25% faster on the Skylake-X. A short "
 "noisy pair abandons the diagonal search in favour of the full matrix (a 1 kbp pair at 15% divergence, 19.2 to "
 "11.8 µs), and the first bound of a short pair is 1.7 times its extrapolation, at most 128 above it."),
("p", "**Gather masks.** A gather consumes its mask register, and compilers generate the required all-ones mask with "
 "`kxnor`, which depends on the previous value of that register; consecutive gathers were therefore serialised at "
 "about 21 cycles per group. Deriving the mask from the group's own lanes removes this dependency. Gathered "
 "extension made the 1 kbp reads, the SARS-CoV-2 genomes and 100 kbp pairs at 1 to 2% divergence 7 to 19% faster "
 "(Table S1, commit 45697e0)."),
("p", "**Measured thresholds.** A full bit-parallel sweep of a 64-row word costs 35%, 45% and 77% of a "
 "diagonal-transition step on AVX-512, AVX2 and NEON, respectively, measured on 1 kbp pairs at 15% divergence (the "
 "AVX2 value on the Skylake-X with code compiled for Haswell). On AVX-512, seeds reduce the running time only from "
 "about 86 kbp: below that length, the Skylake-X aligned real reads of 16 to 64 kbp 14% faster without seeds, and "
 "uniform 30 kbp pairs at 5% divergence 30% faster, whereas at 100 kbp seeds reduced the time by 7 to 24%. With AVX2 "
 "and NEON, seeds are used from 16 kbp. For pairs whose two extrapolations agree, as with uniform errors, the "
 "inexact-seed threshold is 20% on every platform, and an extrapolated error rate of at least one edit in ten bases "
 "selects inexact seeds directly from 128 kbp on AVX-512 and from 64 kbp with narrower vectors: once every band "
 "used 256-column tiles, 100 kbp pairs at 10% aligned in 6.8 ms on exact seeds against 10.3 ms on inexact ones on "
 "the Skylake-X, while 500 kbp pairs at 10 and 20% ran 7 to 10% slower without inexact seeds."),
("p", "**Every threshold measured again.** Each of the 41 thresholds that select an algorithm or bound a search, "
 "from the step budgets of the diagonal transition to the seed gates and the traceback's limits, was halved and "
 "doubled in turn, and each of the 83 builds was timed on 17 unit-cost samples, with and without traceback, or on 8 "
 "gap-affine ones, with the unchanged build timed after every six; its spread was 1.2% in the median. About thirty "
 "thresholds moved no sample by more than 3% either way, and most of the others traded one dataset against "
 "another, such as the first bound of an untrusted band, or the agreement and give-up rules for extrapolations, "
 "and were left as they were. One was stale: the length from which inexact seeds are used, above."),
("p", "**Interleaving.** llvm-mca estimates the interleaved AVX-512 loop at 24 cycles per column, against a "
 "throughput bound of 11, so the loop remains limited by the lane rotation. A third group still fits in the 32 "
 "AVX-512 registers, and llvm-mca predicted 15% fewer cycles per group, but the measured time changed by less than "
 "1.5%."),

("sub", "Seed heuristic"),
("p", "**Batched lookups.** Query rows are processed in batches of 16: the bucket bounds of all rows are loaded "
 "first, so that the loads overlap, and their entries are prefetched; each entry stores a seed together with its "
 "code, so that a single load provides both. A lookup tests each candidate window with a branch-free expression. "
 "Together with considering only seeds that can still chain to the end, these steps reduced the inexact setup on "
 "the ultra-long reads with genetic variation from 5.02 to 2.90 s in total. On 100 kbp pairs at 10% divergence, "
 "inexact-seed construction accounted for 58% of the running time on the Skylake-X; there exact seeds now win "
 "(Measured thresholds, above)."),
("p", "**Admissibility check and correction.** Comparing the heuristic with the true remaining cost at every cell of "
 "random pairs showed that an earlier version of our implementation overestimated the remaining cost by one on "
 "optimal paths in 1 to 2% of pairs. The cause was a filter that discarded matches whose gap to the end of the pair "
 "exceeded the potential of the intervening seeds: a match ending one diagonal beyond this limit was discarded, "
 "although a path through it costs only one more than the potential, and the neighbour that would otherwise bound "
 "such paths could be removed by local pruning. Such a match now chains to the end at the cost of its excess, "
 "entering layer \\(2-1=1\\); no further match can follow it. After the correction, no cell exceeded the true "
 "cost in 2,200 random pairs, compared with 34,450 cells in 400 pairs before. The reported distances agreed with "
 "those of the other aligners in all benchmarks, both before and after the correction. The overestimate arose from "
 "our filter; we do not claim that it applies to A*PA or A*PA2."),
("p", "**Exact-seed filter.** A 32-bit filter per seed, placed before a half-full open-addressing hash table, "
 "resolves most lookups with a single predictable test. On the M2, matching the 320 mid-length reads took 32 ms "
 "instead of 87 ms, and 30 kbp pairs at 5% divergence became 15% faster."),
("p", "**Local pruning.** The diagonal-transition search used for local pruning updates its wavefronts in place, with "
 "an unreachable diagonal on either side replacing three bounds checks per diagonal."),

("sub", "Gap-affine costs"),
("p", "**Compact traceback.** For the traceback of the bidirectional wavefront, dinara-align stores, for every cost "
 "and diagonal, the furthest-reaching column and one byte recording the predecessor of each layer: five bytes per "
 "diagonal, compared with twelve in the high-memory mode of WFA. When the stored wavefronts would exceed the memory "
 "limit (80 MB by default), the pair is split at a cell through which an optimal path passes. Under a general "
 "substitution matrix, substitution scores are obtained by byte shuffles (`pshufb`, `tbl`) for matrices of up to "
 "16 entries and by gathers otherwise."),
("p", "**Gathers.** In the gap-affine wavefront, AVX-512 gathers reduced the time for 1 kbp reads from 149 to 129 µs. "
 "Deriving each gather's mask from the group's own lanes, with one reduction after the loop rather than one per "
 "group, reduced it further to 111 µs, and the time for 10 kbp reads from 2.34 to 1.84 ms. Extension along matches "
 "was moved into the wavefront step, with columns stored through scalar loads and stores rather than reinserted "
 "into a vector, and retained wavefronts are held in fixed blocks that are never moved; memory moves had accounted "
 "for 7.5% of profiling samples (Table S1)."),
("p", "**Memory layout.** The wavefront updates the bounds and furthest reaches of its slots at every cost. When "
 "these were separate arrays of a few words, the allocator placed them next to another thread's data, so that the "
 "cache line moved between cores at every write; ten threads on the Skylake-X then required twice the cycles per "
 "pair of one thread, with 40 cache lines fetched from another core per pair. Each such array now occupies cache "
 "lines that no other allocation shares. A worker retains its profile, reversed sequences and search buffers between "
 "pairs, so that a warm worker performs no allocation; per-pair allocations had accounted for one third of the time "
 "of a short read, and caused contention in the allocator with many threads. Padding each row of a wavefront buffer "
 "by one cache line, to avoid 4 KB aliasing between loads and earlier stores, removed 1.5 billion blocked loads over "
 "ten 100 kbp pairs. Separating the working memory of threads by two cache lines rather than one reduced the cache "
 "lines fetched from other cores from 3.4 million to 0.3 million over 2 million reads on ten threads."),

("sub", "Other modes"),
("p", "**Extension.** Extension alignment uses a single wavefront search from the anchor, with the match score folded "
 "into the costs, and records the best end point on the furthest anti-diagonal of each wavefront. The search stops "
 "once a proven bound shows that no higher cost can yield a better score, and the prefix up to the best end point is "
 "then aligned globally. With an end bonus, the search also records the best point on the last row of the query and "
 "continues only while such an alignment can still come within the bonus; this reduced the time for an extension "
 "with an end bonus from 1.02 ms to 186 µs."),
("p", "**Ends-free alignment.** Under gap costs with free ends, the forward search starts at cost zero on every "
 "diagonal on which a free start lies. The span is determined first: one search from the free side finds the end "
 "on the highest diagonal, and a narrow search backwards finds the start; the span is then aligned globally. On "
 "x86, this aligns a read at the start of a reference, an overlap or a suffix 1.4 to 3.5 times faster, and a read "
 "anywhere within the reference 0 to 14% slower, than the previous approach. Keeping the search state outside an "
 "optional wrapper while the search runs, so that its fields remain in registers, reduced the time for an infix "
 "alignment at (4, 6, 2) from 3.92 to 3.00 ms."),
("p", "**Infix and prefix alignment at unit costs.** Infix alignment uses the bit-parallel sweep with the cut-off of "
 "Edlib, with a first bound of 64; a failed try extrapolates the distance from the deepest row its band kept "
 "within the bound and the retry aims past it, which on the M2 took a 10 kbp infix from 3.39 to 1.99 ms; words "
 "entering the band are initialised to +1, and the last word advances "
 "one column at a time. Prefix alignment sweeps only the band permitted by its bound: the top of the band moves down "
 "with the diagonal, no column beyond the read length plus the best end so far is computed, an attempt whose band "
 "becomes empty stops, and at one eighth, one quarter and one half of the columns an attempt that can no longer "
 "succeed is abandoned. This reduced the time for a prefix distance from 23.7 to 20.5 µs on the Skylake-X."),
("p", "**A band derived from the score.** A global alignment under a single match score and a single mismatch score "
 "stores only the diagonals permitted by its score. In the converted costs, every gapped letter contributes at "
 "least \\(2e+a\\) to the total \\(C=a(n+m)-2S\\), and a path that deviates by \\(t\\) diagonals from both the main "
 "and the end diagonal contains at least \\(2t\\) gapped letters, so only diagonals within \\(C/(2(2e+a))\\) of "
 "these can contain an optimal path. For a 1 kbp pair at 5% divergence, 137 diagonals are stored instead of 1,000. "
 "Under any other substitution matrix, the alignment is traced in a band that grows until the score it yields "
 "proves that it holds every optimal path (main text, Section 3.8), each cell's decision kept in one byte. The backward sweep "
 "from the end of a local alignment and each half of a linear-space alignment also use 16-bit lanes while the scores "
 "fit, which made 1 kbp local alignments 1.2 times and linear-space global alignments 1.4 times faster on the M2."),

("sub", "Batches"),
("p", "**Letter layout.** The pairs of a group are interleaved so that one 64-bit load reads eight letters of a "
 "lane, and eight lanes are transposed in registers by two shuffles."),
("p", "**Alignment in 8-bit lanes.** Batch alignments remain in 8-bit lanes, 64 pairs per AVX-512 register, and a "
 "pair whose cost exceeds the 8-bit range is passed to a single call. With 16-bit lanes, half as many lanes would "
 "sweep the wide band that such a pair requires, which on the Skylake-X took as long as the searches on 1 kbp reads "
 "at 10% divergence, and twice as long on NEON."),
("p", "**Cache-efficient traceback.** Storing the traceback flags of a group's whole band would exceed the caches on "
 "long pairs and touch new memory pages for every group; with ten threads, this doubled the time for 1 kbp pairs "
 "relative to the searches alone. The sweep therefore stores one row of every layer every \\(s\\) rows and "
 "recomputes a segment to obtain its flags during the traceback. With \\(R\\) rows, \\(s\\) is the square root of "
 "\\(R\\) times the bytes per cell of a stored row divided by those of a flag, which balances the stored rows "
 "against the flags of one segment and minimises memory use."),
("p", "**Lane groups.** A group of 64 bytes occupies one register on AVX-512 and several on narrower vector units; "
 "with the 128-bit registers of the M2, 1 kbp pairs at 10% divergence required 73 µs per pair under affine costs "
 "with four registers, compared with 132 µs with one. Local scores use 32 lanes of 16 bits or 16 lanes of 32 bits, "
 "four independent registers' worth of work; on the M2, this scored 20,000 references against a read in 83 ms, "
 "compared with 160 ms with one register. Pairs whose first band does not satisfy the bound are grouped by the exact "
 "width of their second band, and a group stops once every lane has exceeded its cost limit."),

("p", "**Work shares.** Batch calls spread their pairs over threads, longest first: the pairs are dealt into "
 "classes of lengths within an eighth of each other, with no sort, since sorting the short reads on one thread had "
 "taken a fifteenth of the work ten threads did on them. Each worker then takes a share of half of what is left "
 "divided among the workers, as OpenMP's guided schedule does, but never more than twice its last share, starting "
 "from one; a first share of a twentieth of the batch had left one worker holding the longest pairs while the "
 "others finished, adding half the batch's time on the Skylake-X. Building the packed batch is split across the "
 "threads in the same way; on one thread it had taken three times as long as the GPU took to score it."),
("p", "**Free ends in lanes.** Ends-free alignments take the lanes too, their span found first as the tie rule "
 "picks it for a single pair: the end, of the cells an optimal alignment may stop on the one furthest along the "
 "reference, from a sweep that proves its cost with every optimal path inside its band; then the start, by a sweep "
 "back from that end over the reversed prefixes; then the global alignment of the span. On one core of the Skylake-X, "
 "150 bp reads at 3% errors in 200 bp windows aligned at (4, 6, 2) in 8.0 µs per pair as an infix, against 26.6 µs "
 "one pair at a time, in 12.1 against 29.3 µs as a prefix, and in 20.0 against 50.3 µs with 60 letters free at the "
 "end of the window and the start of the read; at unit costs that overlap took 8.3 against 32.4 µs."),
("p", "**GPU scores.** Batch scores and alignments under a substitution table, and unit-cost distances, also run on "
 "NVIDIA GPUs and Apple silicon. Each pair is first scored over a band of sixteen diagonals by one thread, its "
 "scores and gap layers one register per diagonal, swept anti-diagonal by anti-diagonal; the band's score stands "
 "where the pair's own cost certifies it (Proposition 1 of the main text), and the thread stops once the cheapest "
 "cell of its anti-diagonal costs more than every path off the band. The remaining pairs, about one in six of the "
 "short-read batch, are scored over the whole matrix while the host packs the next chunk, several to a warp, each "
 "lane holding a run of its pair's columns in registers and passing its last column on by shuffle. The cells of a "
 "global alignment are held shifted by their anti-diagonal times the gap extension, which reduces Gotoh's "
 "recurrence from five additions per cell to three. On an NVIDIA GeForce RTX 2070, the 500,000 short-read pairs "
 "were scored in 36 ms at (1, 2, 1), the whole call included and most of it the host packing the batch, against "
 "27 ms on ten cores of the Skylake-X."),

("sub", "Code generation"),
("p", "**Avoiding division by a variable.** A 64-bit division by a variable takes tens of cycles on x86. The "
 "potential of the heuristic, required at every row that the band prunes, involves a division by the seed length, "
 "so each seed length has its own code path with a constant divisor, which compiles to a multiplication; batches "
 "likewise avoid a division per pair."),
("p", "**Kernels compiled separately.** When inlined next to each other, or selected per band round, the "
 "bit-parallel kernels for plain bases and for additional symbols slowed the plain kernel by 1 to 4%; when inlined "
 "into the seed construction, the local-pruning search slowed the construction loop by a few percent. Separate "
 "seeded and unseeded versions of the band keep the code of the heuristic out of unseeded pairs, which it had slowed "
 "by 3 to 7%."),
("p", "**Operand order.** Under a substitution matrix, the recurrence first takes the maximum of the two gap layers, "
 "so that the substitution score, the operand that arrives last, waits on only one maximum; a score obtained by "
 "comparison is available early and is used first. This was measured on the M2."),
("p", "**Branches.** On a 1 Mbp pair at 15% divergence, branch mispredictions account for about 7% of the cycles, "
 "mostly at short loop exits, compared with 2 to 3% in the bit-parallel sweep. A fully branch-free seed scan would "
 "spend about 25 million cycles per pair to save about 34 million, a gain of 1 to 2%, and was not adopted."),

("sec", "Measured effect of each change"),
("table", "changes"),

("sec", "Additional results"),
("p", "Table S2 gives the results for the other alignment modes, each compared on one pinned thread of the Skylake-X "
 "with the aligners that support the same mode. Table S3 gives the complete gap-affine comparison at costs "
 "(4, 6, 2). Table S4 gives the mean and median time per alignment at unit costs on every dataset, Table S5 the "
 "samples, Table S6 the memory and Table S7 the short-read batch. The approximate aligners of the A*PA2 benchmark, WFA-adaptive and Block Aligner, are no longer "
 "run. KSW2 and BSAlign, which ran out of time or memory on most datasets of 500 kbp and more, are in the "
 "benchmark's own tables (benchmarks/README.md in the repository)."),
("table", "modes"),
("table", "affine_full"),
("table", "unit_full"),
("table", "datasets"),
("table", "memory"),
("table", "batch"),

("sec", "Robustness"),
("p", "Three changes address inputs on which the running time had been disproportionate:"),
("bullets", [
 "**Low-complexity sequence.** Wavefront layers that spread into a staircase pattern are handled explicitly, and the "
 "pruning search is bounded; a poly-A sequence of 32,000 bases aligned against one of 16,000 bases took 0.15 s "
 "instead of 17 s.",
 "**Repetitive seeds.** Seeds are abandoned when they exceed 64 matches per seed on average; setting up a 200 kbp "
 "tandem repeat took 0.045 s and 34 MB instead of 73 s and 760 MB.",
 "**Unreachable wavefronts.** A wavefront that reaches no cell is not tested further; one such case took 13 ms "
 "instead of 1.5 s.",
]),

("sec", "Approaches evaluated and not adopted"),
("bullets", [
 "**Multi-threading within a pair.** Striped bands over eight threads aligned 1 Mbp pairs at 15% divergence in "
 "575 ms, compared with 1,122 ms on one thread. The implementation was removed because the library runs on the "
 "calling thread and distributing pairs over threads yields a larger gain.",
 "**A third bit-parallel group on AVX-512.** The measured time changed by less than 1.5%.",
 "**One sixteen-lane group instead of two eight-lane groups.** Equal on long reads, and up to 2% slower on short "
 "divergent pairs.",
 "**A fully branch-free seed scan.** It would save 1 to 2% at a cost of about 25 million cycles per pair.",
 "**A Bloom filter for inexact keys.** It halved the number of candidates, but reduced the alignment time by only "
 "1.5 to 2%.",
 "**Shorter look-ahead, or no local pruning.** Without pruning, the ultra-long reads with and without genetic "
 "variation were 21% and 72% slower, respectively.",
 "**A* search within the diagonal transition, and an upper bound from beam search.** Each required several tuned "
 "parameters for small gains.",
 "**Traceback during the sweep.** The stored tile boundaries occupy 0.6 to 13 MB, compared with 10 to 17 MB for the "
 "sequences and bit planes, so little memory could be saved.",
]),
]

TABLES = {
"modes": {
  "label": "S2",
  "caption": "Other modes: mean time per pair in ms on one pinned thread, the faster of two passes. Every tool's costs or scores agree on every pair.",
  "header": ["Workload", "dinara-align", "Next fastest"],
  "align": "lrl",
  "wrap": [0],
  "rows": [
    ["Infix at unit costs, 1 kbp read at 10% in a 3 kbp window", "**0.132**", "Edlib 0.273"],
    ["Infix at unit costs, 10 kbp read at 5% in a 30 kbp window", "**2.01**", "Edlib 13.6"],
    ["Prefix at unit costs, 1 kbp read against 2 kbp", "**0.044**", "WFA2-lib 0.047"],
    ["Infix at (4, 6, 2), 1 kbp read in a 3 kbp window", "**3.93**", "WFA2-lib 3.94"],
    ["Extension, match 2 at (4, 6, 2)", "**0.734**", "KSW2 2.06"],
    ["Extension with an end bonus of 50", "**0.166**", "KSW2 1.23"],
    ["Two-piece gaps (4, 6, 2, 24, 1), 5 kbp at 5%", "**8.30**", "WFA2-lib 15.9"],
    ["Table, global, 1 kbp at 10%", "**0.198**", "parasail 2.92"],
    ["Table, local, 1 kbp read in a 10 kbp window", "**1.91**", "SSW 3.27"],
    ["Local, 600 short noisy pairs (20–700 bp)", "**0.032**", "SSW 0.051"],
    ["Local, 1 kbp read at 10% in a 10 kbp window", "**1.36**", "SSW 3.15"],
    ["Local, 10 kbp at 5% against 12 kbp", "**17.5**", "SSW 85.1"],
    ["Overlap of 2 kbp reads by 0.5–1.5 kbp", "**0.724**", "parasail 15.4"],
    ["Infix with a reward, 1 kbp read in a 3 kbp window", "**0.693**", "parasail 15.2"],
  ]},
"affine_full": {
  "label": "S3",
  "fit": True,
  "caption": "Gap-affine alignment at costs (4, 6, 2), all datasets: mean time per alignment in ms, with a 30 s budget per sample. Superscripts give the pairs completed where not all were; a dash marks a tool that completed none within 30 s, o.o.m. one the 8 GB limit stopped.",
  "header": ["Dataset", "dinara-align", "WFA", "BiWFA", "KSW2", "BSAlign"],
  "align": "lrrrrr",
  "rows": [
    ["ont-1k", "**0.105**", "0.188", "0.308", "1.10", "0.336"],
    ["ont-10k", "**1.60**", "4.89", "5.86", "30.0", "7.84"],
    ["ont-50k", "**21.5**", "63.4", "54.9", "341^81/104", "87.0"],
    ["ont-500k", "–", "o.o.m.", "–", "o.o.m.", "o.o.m."],
    ["ont-500k-genvar", "**22,000^1/4**", "o.o.m.", "–", "o.o.m.", "o.o.m."],
    ["SARS-CoV-2", "**0.710**", "2.05", "1.78", "87.9", "341"],
    ["Illumina", "**0.00145**", "0.00155", "0.00281", "0.013", "0.054"],
    ["PacBio HiFi", "**0.095**", "0.227", "0.231", "9.72", "86.0"],
    ["10 kbp, 1%", "**0.199**", "0.416", "0.646", "21.3", "37.0"],
    ["10 kbp, 5%", "**2.32**", "7.14", "8.72", "96.1", "37.0"],
    ["10 kbp, 10%", "**7.06**", "23.1", "27.1", "202", "37.3"],
    ["10 kbp, 20%", "**20.1**", "66.0", "75.1", "375^81/99", "37.4"],
    ["100 kbp, 1%", "**11.9**", "39.9", "39.0", "2,420", "o.o.m."],
    ["100 kbp, 5%", "**331**", "808", "714", "9,170^3/9", "o.o.m."],
    ["100 kbp, 10%", "**1,160**", "o.o.m.", "2,430", "17,300^1/10", "o.o.m."],
    ["100 kbp, 20%", "**3,610^8/10**", "o.o.m.", "7,360^4/10", "o.o.m.", "o.o.m."],
    ["500 kbp, 1%", "**477**", "1,290", "942", "o.o.m.", "o.o.m."],
    ["500 kbp, 5%", "**9,910^3/4**", "o.o.m.", "19,400^1/4", "o.o.m.", "o.o.m."],
    ["500 kbp, 10%", "–", "o.o.m.", "–", "o.o.m.", "o.o.m."],
    ["500 kbp, 20%", "–", "o.o.m.", "–", "o.o.m.", "o.o.m."],
    ["1 Mbp, 1%", "**2,300**", "o.o.m.", "3,940", "o.o.m.", "o.o.m."],
    ["1 Mbp, 5%", "–", "o.o.m.", "–", "o.o.m.", "o.o.m."],
    ["1 Mbp, 10%", "–", "o.o.m.", "–", "o.o.m.", "o.o.m."],
    ["1 Mbp, 20%", "–", "o.o.m.", "–", "o.o.m.", "o.o.m."],
  ]},
"unit_full": {
  "label": "S4",
  "fit": True,
  "caption": "Unit costs, all datasets: mean time per alignment in ms, with the median in parentheses, traceback included, a 30 s budget per sample. Superscripts give the pairs completed where not all were, the median over those; a dash marks a tool that completed none, o.o.m. one the 8 GB limit stopped.",
  "header": ["Dataset", "dinara-align", "A*PA2-full", "A*PA2-simple", "A*PA", "Edlib", "BiWFA", "WFA", "QuickEd"],
  "align": "lrrrrrrrr",
  "rows": [
    ["ont-1k", "**0.021 (0.021)**", "0.081 (0.089)", "0.048 (0.054)", "0.612 (0.556)", "0.124 (0.118)", "0.044 (0.046)", "0.031 (0.032)", "0.075 (0.073)"],
    ["ont-10k", "**0.156 (0.109)**", "0.391 (0.264)", "0.247 (0.177)", "15.3 (3.29)", "1.09 (0.569)", "0.741 (0.278)", "0.599 (0.225)", "0.688 (0.451)"],
    ["ont-50k", "**0.580 (0.229)**", "1.27 (0.513)", "0.944 (0.307)", "263 (9.26)", "6.14 (2.05)", "7.43 (0.872)", "6.38 (0.726)", "2.31 (1.05)"],
    ["ont-500k", "**118 (73.6)**", "283 (110)", "687 (398)", "–", "5,000 (3,050)", "6,330 (2,190)", "o.o.m.", "847 (585)"],
    ["ont-500k-genvar", "**156 (156)**", "225 (196)", "511 (564)", "–", "5,230 (5,200)", "6,200 (6,630)", "o.o.m.", "1,390 (1,320)"],
    ["SARS-CoV-2", "**0.286 (0.146)**", "2.00 (2.07)", "0.728 (0.662)", "6.49 (2.07)", "8.18 (7.74)", "0.909 (0.456)", "0.987 (0.362)", "3.90 (4.27)"],
    ["Illumina", "**0.000576 (0.000482)**", "0.014 (0.014)", "0.00410 (0.00405)", "0.014 (0.013)", "0.023 (0.023)", "0.00204 (0.00179)", "0.00112 (0.000999)", "0.016 (0.015)"],
    ["PacBio HiFi", "**0.013 (0.00527)**", "0.977 (0.921)", "0.259 (0.237)", "1.05 (0.962)", "3.35 (3.14)", "0.086 (0.069)", "0.066 (0.047)", "0.950 (0.881)"],
    ["10 kbp, 1%", "**0.026 (0.023)**", "0.729 (0.724)", "0.160 (0.156)", "0.564 (0.560)", "3.22 (3.21)", "0.092 (0.088)", "0.066 (0.058)", "0.611 (0.609)"],
    ["10 kbp, 5%", "**0.267 (0.242)**", "0.797 (0.794)", "0.288 (0.285)", "1.08 (1.07)", "3.79 (3.80)", "0.642 (0.635)", "0.508 (0.495)", "1.25 (1.24)"],
    ["10 kbp, 10%", "**0.355 (0.361)**", "0.726 (0.714)", "0.450 (0.439)", "3.26 (3.23)", "4.56 (4.55)", "2.26 (2.26)", "1.81 (1.66)", "1.99 (2.01)"],
    ["10 kbp, 20%", "**0.512 (0.468)**", "1.44 (1.44)", "1.02 (1.01)", "205 (203)", "5.93 (5.93)", "7.33 (7.32)", "5.32 (5.11)", "1.55 (1.56)"],
    ["100 kbp, 1%", "**1.50 (1.11)**", "7.53 (7.27)", "3.26 (3.14)", "7.76 (7.53)", "57.6 (57.4)", "3.71 (3.47)", "4.05 (4.28)", "12.5 (12.8)"],
    ["100 kbp, 5%", "**5.10 (4.77)**", "8.45 (8.11)", "14.9 (14.7)", "14.4 (14.8)", "146 (146)", "63.3 (63.2)", "55.0 (51.2)", "24.0 (23.8)"],
    ["100 kbp, 10%", "**7.28 (7.02)**", "13.9 (13.3)", "27.6 (27.5)", "48.3 (45.8)", "247 (247)", "218 (219)", "182 (170)", "39.0 (38.9)"],
    ["100 kbp, 20%", "**15.7 (15.6)**", "36.7 (36.0)", "56.8 (56.6)", "–", "446 (446)", "695 (695)", "551 (510)", "62.8 (62.7)"],
    ["500 kbp, 1%", "**33.5 (32.3)**", "42.7 (41.3)", "72.3 (72.2)", "41.2 (39.2)", "670 (670)", "80.8 (79.4)", "74.8 (63.8)", "124 (124)"],
    ["500 kbp, 5%", "**31.4 (30.2)**", "45.0 (44.4)", "214 (213)", "83.1 (82.5)", "2,280 (2,280)", "1,560 (1,560)", "1,480 (1,250)", "401 (400)"],
    ["500 kbp, 10%", "**68.4 (64.3)**", "119 (117)", "425 (424)", "406 (406)", "4,380 (4,380)", "5,450 (5,450)", "5,050 (4,250)", "711 (711)"],
    ["500 kbp, 20%", "**255 (257)**", "545 (545)", "901 (896)", "–", "8,520^3/4 (8,520)", "18,000^1/4 (18,000)", "o.o.m.", "1,240 (1,240)"],
    ["1 Mbp, 1%", "**59.7 (56.8)**", "89.1 (87.0)", "253 (253)", "97.3 (93.3)", "2,380 (2,380)", "366 (364)", "391 (353)", "398 (398)"],
    ["1 Mbp, 5%", "**65.0 (62.8)**", "118 (117)", "818 (815)", "183 (182)", "8,670^3/4 (8,670)", "6,350 (6,350)", "o.o.m.", "1,440 (1,440)"],
    ["1 Mbp, 10%", "**183 (183)**", "369 (366)", "1,650 (1,640)", "978 (980)", "17,300^1/4 (17,300)", "22,200^1/4 (22,200)", "o.o.m.", "2,650 (2,650)"],
    ["1 Mbp, 20%", "**892 (688)**", "2,020 (2,010)", "3,520 (3,500)", "–", "–", "–", "o.o.m.", "4,800 (4,800)"],
  ]},
"datasets": {
  "label": "S5",
  "caption": "The samples: pairs in each sample of the dataset, their lengths and their divergence, the edit distance over the longer sequence, as QuickEd's Table S1 describes its datasets.",
  "header": ["Dataset", "Pairs", "Length, min / mean / max (bp)", "Divergence, min / mean / max (%)"],
  "align": "lrrr",
  "rows": [
    ["ont-1k", "1221 of 12477", "42 / 830 / 1,089", "0 / 11 / 21"],
    ["ont-10k", "277 of 5000", "342 / 3,643 / 9,791", "4 / 13 / 20"],
    ["ont-50k", "104 of 10000", "360 / 9,633 / 42,626", "4 / 12 / 19"],
    ["ont-500k", "4 of 50", "553,834 / 642,465 / 683,520", "4 / 7 / 13"],
    ["ont-500k-genvar", "4 of 48", "594,436 / 662,980 / 861,994", "5 / 7 / 10"],
    ["SARS-CoV-2", "33 of 10000", "29,244 / 29,746 / 29,864", "0 / 1 / 7"],
    ["Illumina", "4032 of 5992", "77 / 248 / 251", "0 / 0 / 10"],
    ["PacBio HiFi", "69 of 890", "10,956 / 14,372 / 22,816", "0 / 0 / 1"],
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
  "label": "S6",
  "caption": "Memory at unit costs: the largest growth of peak resident memory aligning one pair of the sample, in MB, as in A*PA2, and the peak resident memory of the whole run in parentheses, as in QuickEd, runtime and input included. A dash marks a tool that completed no pair within 30 s.",
  "header": ["Dataset", "dinara-align", "A*PA2-full", "A*PA2-simple", "Edlib", "BiWFA", "QuickEd"],
  "align": "lrrrrrr",
  "rows": [
    ["ont-500k", "48 (64)", "82 (89)", "82 (90)", "5.3 (15)", "15 (29)", "22 (26)"],
    ["ont-500k-genvar", "39 (64)", "48 (57)", "32 (55)", "4.8 (16)", "12 (27)", "20 (27)"],
    ["100 kbp, 20%", "5.2 (18)", "7.3 (12)", "4.7 (9)", "1.4 (8)", "3.6 (11)", "8.8 (12)"],
    ["500 kbp, 20%", "53 (82)", "59 (66)", "53 (61)", "4.1^3/4", "15^1/4", "16 (20)"],
    ["1 Mbp, 20%", "123 (178)", "184 (196)", "192 (205)", "–", "–", "26 (31)"],
  ]},
"batch": {
  "label": "S7",
  "caption": "Global scores of 500,000 short-read pairs: each time is the whole call, in ms.",
  "header": ["Costs", "Threads", "dinara-align", "WFA2-lib", "KSW2", "parasail", "Edlib"],
  "align": "llrrrrr",
  "rows": [
    ["Affine (1, 2, 1)", "1", "**230**", "1,060", "10,900", "21,700", "–"],
    ["Unit", "1", "**212**", "450", "–", "–", "2,000"],
    ["Affine (1, 2, 1)", "10", "**27.0**", "116", "1,210", "2,320", "–"],
    ["Unit", "10", "**29.0**", "64.0", "–", "–", "216"],
  ]},
"changes": {
  "label": "S1",
  "caption": "Changes behind the results and their measured effects, as recorded when each was made. Machines: the i9-7900X (Skylake-X, AVX-512) and the Apple M2 (NEON); “not recorded” where the commit does not say.",
  "header": _changes[2],
  "align": "llll",
  "wrap": [1, 2],
  "rows": _changes[3],
},
}

REFS = {}
BIBKEYS = {}
