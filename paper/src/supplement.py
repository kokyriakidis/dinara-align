"""The supplementary material, in the same markup as main_text.py."""
from change_record import SUPPLEMENT as _OLD

TITLE = "Supplementary material for “dinara-align: exact pairwise alignment specialised to the target processor”"

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
])

BODY = [
("p", "This supplement provides implementation details, per-change measurements and additional results that the main "
 "text summarises. Section S1 describes techniques not covered in the main text, Section S2 lists the measured "
 "effect of each change, Section S3 gives the complete results for gap-affine costs and the other alignment modes, "
 "and Sections S4 and S5 describe robustness measures and approaches that were evaluated and not adopted. Unless "
 "stated otherwise, the measurements in Sections S1, S2, S4 and S5 were taken when each change was introduced, on "
 "the machine indicated, and compare the change with the preceding version; the results in the main text and in "
 "Section S3 were obtained on the Skylake-X at commit `dc980d7`, and for gap-affine costs and batches at `289bcf0`, "
 "whose code those workloads share. Two machines were used: an Intel Core i9-7900X (Skylake-X, AVX-512) and an "
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
 "Edlib, with a bound doubling from 64; words entering the band are initialised to +1, and the last word advances "
 "one column at a time. Prefix alignment sweeps only the band permitted by its bound: the top of the band moves down "
 "with the diagonal, no column beyond the read length plus the best end so far is computed, an attempt whose band "
 "becomes empty stops, and at one eighth, one quarter and one half of the columns an attempt that can no longer "
 "succeed is abandoned. This reduced the time for a prefix distance from 23.7 to 20.5 µs on the Skylake-X."),
("p", "**A band derived from the score.** A global alignment under a single match score and a single mismatch score "
 "stores only the diagonals permitted by its score. In the converted costs, every gapped letter contributes at "
 "least \\(2e+a\\) to the total \\(C=a(n+m)-2S\\), and a path that deviates by \\(t\\) diagonals from both the main "
 "and the end diagonal contains at least \\(2t\\) gapped letters, so only diagonals within \\(C/(2(2e+a))\\) of "
 "these can contain an optimal path. For a 1 kbp pair at 5% divergence, 137 diagonals are stored instead of 1,000. "
 "Under any other substitution matrix, the score of the sweep bounds the band in the same way. The backward sweep "
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
 "(4, 6, 2). Table S4 gives the mean and median time per alignment at unit costs on every dataset, and Table S5 "
 "the samples. The approximate aligners of the A*PA2 benchmark, WFA-adaptive and Block Aligner, are no longer "
 "run. KSW2 and BSAlign, which ran out of time or memory on most datasets of 500 kbp and more, are in the "
 "benchmark's own tables (benchmarks/README.md in the repository)."),
("table", "modes"),
("table", "affine_full"),
("table", "unit_full"),
("table", "datasets"),

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
  "caption": "Other modes: mean time per pair on one pinned thread, the faster of two passes. Every tool's costs or scores agree on every pair.",
  "header": ["Workload", "dinara-align", "Next fastest"],
  "align": "lrl",
  "wrap": [0],
  "rows": [
    ["Infix at unit costs, 1 kbp read at 10% in a 3 kbp window", "**132 µs**", "Edlib 273 µs"],
    ["Infix at unit costs, 10 kbp read at 5% in a 30 kbp window", "**2.01 ms**", "Edlib 13.6 ms"],
    ["Prefix at unit costs, 1 kbp read against 2 kbp", "**44 µs**", "WFA2-lib 47 µs"],
    ["Infix at (4, 6, 2), 1 kbp read in a 3 kbp window", "**3.93 ms**", "WFA2-lib 3.94 ms"],
    ["Extension, match 2 at (4, 6, 2)", "**734 µs**", "KSW2 2.06 ms"],
    ["Extension with an end bonus of 50", "**166 µs**", "KSW2 1.23 ms"],
    ["Two-piece gaps (4, 6, 2, 24, 1), 5 kbp at 5%", "**8.3 ms**", "WFA2-lib 15.9 ms"],
    ["Table, global, 1 kbp at 10%", "**198 µs**", "parasail 2.92 ms"],
    ["Table, local, 1 kbp read in a 10 kbp window", "**1.91 ms**", "SSW 3.27 ms"],
    ["Local, 600 short noisy pairs (20–700 bp)", "**32 µs**", "SSW 51 µs"],
    ["Local, 1 kbp read at 10% in a 10 kbp window", "**1.36 ms**", "SSW 3.15 ms"],
    ["Local, 10 kbp at 5% against 12 kbp", "**17.5 ms**", "SSW 85.1 ms"],
    ["Overlap of 2 kbp reads by 0.5–1.5 kbp", "**724 µs**", "parasail 15.4 ms"],
    ["Infix with a reward, 1 kbp read in a 3 kbp window", "**693 µs**", "parasail 15.2 ms"],
  ]},
"affine_full": {
  "label": "S3",
  "fit": True,
  "caption": "Gap-affine alignment at costs (4, 6, 2), all datasets: mean time per alignment with a 30 s budget per sample. Superscripts give the pairs completed where not all were; a dash marks a tool that completed none within 30 s, o.o.m. one the 8 GB limit stopped.",
  "header": ["Dataset", "dinara-align", "WFA", "BiWFA", "KSW2", "BSAlign"],
  "align": "lrrrrr",
  "rows": [
    ["ont-1k", "**112 µs**", "186 µs", "304 µs", "1.1 ms", "336 µs"],
    ["ont-10k", "**1.71 ms**", "4.8 ms", "5.84 ms", "29.8 ms", "7.84 ms"],
    ["ont-50k", "**22 ms**", "63.5 ms", "54.9 ms", "343 ms^81/104", "87 ms"],
    ["ont-500k", "–", "o.o.m.", "–", "o.o.m.", "o.o.m."],
    ["ont-500k-genvar", "**21.7 s^1/4**", "o.o.m.", "–", "o.o.m.", "o.o.m."],
    ["SARS-CoV-2", "**737 µs**", "2.06 ms", "1.79 ms", "88.4 ms", "341 ms"],
    ["10 kbp, 1%", "**208 µs**", "461 µs", "650 µs", "21.3 ms", "37 ms"],
    ["10 kbp, 5%", "**2.32 ms**", "7.05 ms", "8.72 ms", "97.5 ms", "37 ms"],
    ["10 kbp, 10%", "**7.18 ms**", "22.7 ms", "27.3 ms", "204 ms", "37.3 ms"],
    ["10 kbp, 20%", "**20.6 ms**", "65.2 ms", "75.2 ms", "374 ms^81/99", "37.4 ms"],
    ["100 kbp, 1%", "**14 ms**", "39.7 ms", "39 ms", "2.44 s", "o.o.m."],
    ["100 kbp, 5%", "**337 ms**", "811 ms", "715 ms", "9.25 s^3/9", "o.o.m."],
    ["100 kbp, 10%", "**1.19 s**", "o.o.m.", "2.43 s", "17.5 s^1/10", "o.o.m."],
    ["100 kbp, 20%", "**3.62 s^8/10**", "o.o.m.", "7.34 s^4/10", "o.o.m.", "o.o.m."],
    ["500 kbp, 1%", "**490 ms**", "1.29 s", "942 ms", "o.o.m.", "o.o.m."],
    ["500 kbp, 5%", "**9.92 s^3/4**", "o.o.m.", "19.4 s^1/4", "o.o.m.", "o.o.m."],
    ["500 kbp, 10%", "–", "o.o.m.", "–", "o.o.m.", "o.o.m."],
    ["500 kbp, 20%", "–", "o.o.m.", "–", "o.o.m.", "o.o.m."],
    ["1 Mbp, 1%", "**2.32 s**", "o.o.m.", "3.92 s", "o.o.m.", "o.o.m."],
    ["1 Mbp, 5%", "–", "o.o.m.", "–", "o.o.m.", "o.o.m."],
    ["1 Mbp, 10%", "–", "o.o.m.", "–", "o.o.m.", "o.o.m."],
    ["1 Mbp, 20%", "–", "o.o.m.", "–", "o.o.m.", "o.o.m."],
  ]},
"unit_full": {
  "label": "S4",
  "fit": True,
  "caption": "Unit costs, all datasets: mean time per alignment, with the median in parentheses, traceback included, a 30 s budget per sample. Superscripts give the pairs completed where not all were, the median over those; a dash marks a tool that completed none, o.o.m. one the 8 GB limit stopped.",
  "header": ["Dataset", "dinara-align", "A*PA2-full", "A*PA2-simple", "A*PA", "Edlib", "BiWFA", "WFA", "QuickEd"],
  "align": "lrrrrrrrr",
  "rows": [
    ["ont-1k", "**23 µs (23 µs)**", "81 µs (89 µs)", "48 µs (54 µs)", "612 µs (556 µs)", "124 µs (118 µs)", "44 µs (46 µs)", "31 µs (32 µs)", "75 µs (73 µs)"],
    ["ont-10k", "**157 µs (107 µs)**", "391 µs (264 µs)", "247 µs (177 µs)", "15.3 ms (3.29 ms)", "1.09 ms (569 µs)", "741 µs (278 µs)", "599 µs (225 µs)", "688 µs (451 µs)"],
    ["ont-50k", "**590 µs (224 µs)**", "1.27 ms (513 µs)", "944 µs (307 µs)", "263 ms (9.26 ms)", "6.14 ms (2.05 ms)", "7.43 ms (872 µs)", "6.38 ms (726 µs)", "2.31 ms (1.05 ms)"],
    ["ont-500k", "**118 ms (71.1 ms)**", "283 ms (110 ms)", "687 ms (398 ms)", "–", "5 s (3.05 s)", "6.33 s (2.19 s)", "o.o.m.", "847 ms (585 ms)"],
    ["ont-500k-genvar", "**156 ms (154 ms)**", "225 ms (196 ms)", "511 ms (564 ms)", "–", "5.23 s (5.2 s)", "6.2 s (6.63 s)", "o.o.m.", "1.39 s (1.32 s)"],
    ["SARS-CoV-2", "**285 µs (169 µs)**", "2 ms (2.07 ms)", "728 µs (662 µs)", "6.49 ms (2.07 ms)", "8.18 ms (7.74 ms)", "909 µs (456 µs)", "987 µs (362 µs)", "3.9 ms (4.27 ms)"],
    ["10 kbp, 1%", "**33 µs (27 µs)**", "729 µs (724 µs)", "160 µs (156 µs)", "564 µs (560 µs)", "3.22 ms (3.21 ms)", "92 µs (88 µs)", "66 µs (58 µs)", "611 µs (609 µs)"],
    ["10 kbp, 5%", "**269 µs (244 µs)**", "797 µs (794 µs)", "288 µs (285 µs)", "1.08 ms (1.07 ms)", "3.79 ms (3.8 ms)", "642 µs (635 µs)", "508 µs (495 µs)", "1.25 ms (1.24 ms)"],
    ["10 kbp, 10%", "**354 µs (357 µs)**", "726 µs (714 µs)", "450 µs (439 µs)", "3.26 ms (3.23 ms)", "4.56 ms (4.55 ms)", "2.26 ms (2.26 ms)", "1.81 ms (1.66 ms)", "1.99 ms (2.01 ms)"],
    ["10 kbp, 20%", "**519 µs (481 µs)**", "1.44 ms (1.44 ms)", "1.02 ms (1.01 ms)", "205 ms (203 ms)", "5.93 ms (5.93 ms)", "7.33 ms (7.32 ms)", "5.32 ms (5.11 ms)", "1.55 ms (1.56 ms)"],
    ["100 kbp, 1%", "**1.5 ms (1.16 ms)**", "7.53 ms (7.27 ms)", "3.26 ms (3.14 ms)", "7.76 ms (7.53 ms)", "57.6 ms (57.4 ms)", "3.71 ms (3.47 ms)", "4.05 ms (4.28 ms)", "12.5 ms (12.8 ms)"],
    ["100 kbp, 5%", "**5.02 ms (4.8 ms)**", "8.45 ms (8.11 ms)", "14.9 ms (14.7 ms)", "14.4 ms (14.8 ms)", "146 ms (146 ms)", "63.3 ms (63.2 ms)", "55 ms (51.2 ms)", "24 ms (23.8 ms)"],
    ["100 kbp, 10%", "**7.21 ms (7 ms)**", "13.9 ms (13.3 ms)", "27.6 ms (27.5 ms)", "48.3 ms (45.8 ms)", "247 ms (247 ms)", "218 ms (219 ms)", "182 ms (170 ms)", "39 ms (38.9 ms)"],
    ["100 kbp, 20%", "**16 ms (15.8 ms)**", "36.7 ms (36 ms)", "56.8 ms (56.6 ms)", "–", "446 ms (446 ms)", "695 ms (695 ms)", "551 ms (510 ms)", "62.8 ms (62.7 ms)"],
    ["500 kbp, 1%", "**33.9 ms (31.8 ms)**", "42.7 ms (41.3 ms)", "72.3 ms (72.2 ms)", "41.2 ms (39.2 ms)", "670 ms (670 ms)", "80.8 ms (79.4 ms)", "74.8 ms (63.8 ms)", "124 ms (124 ms)"],
    ["500 kbp, 5%", "**31 ms (29.8 ms)**", "45 ms (44.4 ms)", "214 ms (213 ms)", "83.1 ms (82.5 ms)", "2.28 s (2.28 s)", "1.56 s (1.56 s)", "1.48 s (1.25 s)", "401 ms (400 ms)"],
    ["500 kbp, 10%", "**67.4 ms (63.4 ms)**", "119 ms (117 ms)", "425 ms (424 ms)", "406 ms (406 ms)", "4.38 s (4.38 s)", "5.45 s (5.45 s)", "5.05 s (4.25 s)", "711 ms (711 ms)"],
    ["500 kbp, 20%", "**254 ms (254 ms)**", "545 ms (545 ms)", "901 ms (896 ms)", "–", "8.52 s^3/4 (8.52 s)", "18 s^1/4 (18 s)", "o.o.m.", "1.24 s (1.24 s)"],
    ["1 Mbp, 1%", "**60.2 ms (56.5 ms)**", "89.1 ms (87 ms)", "253 ms (253 ms)", "97.3 ms (93.3 ms)", "2.38 s (2.38 s)", "366 ms (364 ms)", "391 ms (353 ms)", "398 ms (398 ms)"],
    ["1 Mbp, 5%", "**64.1 ms (61.6 ms)**", "118 ms (117 ms)", "818 ms (815 ms)", "183 ms (182 ms)", "8.67 s^3/4 (8.67 s)", "6.35 s (6.35 s)", "o.o.m.", "1.44 s (1.44 s)"],
    ["1 Mbp, 10%", "**184 ms (185 ms)**", "369 ms (366 ms)", "1.65 s (1.64 s)", "978 ms (980 ms)", "17.3 s^1/4 (17.3 s)", "22.2 s^1/4 (22.2 s)", "o.o.m.", "2.65 s (2.65 s)"],
    ["1 Mbp, 20%", "**901 ms (701 ms)**", "2.02 s (2.01 s)", "3.52 s (3.5 s)", "–", "–", "–", "o.o.m.", "4.8 s (4.8 s)"],
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
