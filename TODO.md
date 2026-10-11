# To do

`[x]` done, `[~]` partly done, `[-]` tried or weighed and left, with why; `[ ]` open.

Where the time still goes in a single-threaded alignment, profiled on the Skylake-X (3.3 GHz,
pinned) over all 48 reads of ont-500k-genvar, about 10 s in all.

- [x] **Band rounds that fail (about 1 s of the band's 6.5 s).** A retry now aims a quarter of the
  estimated climb past its estimate, not half: genvar 2% and ont-500k 7.5% faster on the Skylake-X.
  Lowering later rounds' bounds at checkpoints, as the first round does, never ended on some reads. The sweep itself runs near the
  hardware's limit for Myers' recurrence, about 2.4 cycles a word-column on AVX-512, so the band
  only gets faster by sweeping fewer cells. On genvar 16% of the band's word-columns go to rounds
  that fail and are retried (7% on ont-500k), three to six rounds a pair. What is left of them is
  an upper bound to aim at, tried and left below (see "An upper bound on the distance").
- [x] **Wider vectors on AVX-512.** A sweep there now runs two eight-lane groups at a time, one under
  the other in one loop: genvar 7% and ont-500k 5% faster on the Skylake-X, the sample sets within
  noise (measured best of five beside the owner's jobs). One sixteen-lane group did about as well on
  the long reads but cost up to 2% on short divergent pairs; on the M2 both were slower, so NEON and
  AVX2 keep one group.
- [-] **Three eight-lane groups at a time on AVX-512.** llvm-mca puts the paired sweep's steady loop at
  24 cycles a column against a throughput bound of 11, latency-bound on each column's lane rotation,
  and a third group, which still fits AVX-512's 32 registers unspilled, at 30 cycles for three: 15%
  less a group. Measured on the Skylake-X it changed nothing within 1.5% on the long reads and pairs
  (2026-10-06): the bands it needs, 24 words and more, are rare once pruning narrows them. Left.
- [-] **Mispredicted branches in the inexact seed scan.** perf on the Skylake-X (2026-10-06): branch
  misses cost about 7% of a 1 Mbp pair's cycles at 15% and half are `inexact_matches`, the rest in
  `worth_keeping`, `score` and `exact_matches`; the bit-parallel sweep has 2 to 3%. They are short
  loops' exits, not the one-edit test, which is branch-free: about one seed a bucket, so each row's
  loops run 0 to 2 times. A branch-free scan, eight entries a bucket in one vector and the test on
  all eight, would spend about 25M cycles a pair to save about 34M: 1 to 2% on long pairs, nothing
  on short reads. Left for that.
- [-] **Seed setup.** Still 34 to 45% of a long pair's time on the Skylake-X and 9 to 32% on the M2
  (2026-10-06, `9c656ec`): on genvar 157 ms of inexact matching, 94 of local pruning and 19 of layers
  in 698; on 100 kbp pairs at 15% 27, 16 and 1 in 104; on ont-500k 64, 69 and 10 in 468.
  - Inexact matching: the half tables' lookups cost 5.5 ns a row, and the rest is the candidates,
    0.74 a row on genvar, 17% of them within one edit: about 15 ns a row plus 43 a candidate there,
    against 3 and 28 on the M2. A chance candidate shares an 8-base half with the window, likely
    with 41,000 seeds a read over 65,536 halves. Keying on two exact parts, a 12-base prefix and a
    quarter at the end, would cut them some thirtyfold for about eight hashed lookups a row in place
    of two direct ones: perhaps 5 to 10% on genvar and divergent 100 kbp pairs. Tried and left
    (2026-10-09, M2): those keys as a gate before the half tables, a Bloom filter of two bits a key
    and 16 bits a key in all, then one cache line a row keyed by the row's half, found the same
    matches and halved the candidates (genvar 0.74 to 0.37 a row), but the matches themselves are
    dense there, 0.13 a row, and the gate's eight tests a row cost about 11 ns, nearly what the
    candidates did: inexact matching 12 to 15% faster, the alignments 1.5 to 2%. Hashed tables on the
    same keys would cost more a lookup than the gate's tests.
  - Tried and left: no local pruning (genvar 21% and ont-500k 72% slower: it repays itself many
    times), a lookahead of 8 seeds for inexact matches (neutral to 4% slower), and of 6 to 10 for all
    (mid-length reads up to 6% faster on the M2, ont-500k 17% and 100 kbp at 5% 30% slower), and
    scan batches of 32 rows (neutral).
  - Exact matching already checks a 32-bit-a-seed filter before its hash table (87 to 32 ms on the
    M2's mid-length reads); exact pruning, 600 thousand matches of which 174 thousand are kept, is
    now the larger half there. Going further than this means indexing the second sequence instead
    of the seeds, at more memory.
  - Profiled again on the Skylake-X at `2b52c06` (2026-10-10): 70.9 ms a genvar read, inexact
    matching 47%, local pruning 26%, the layer loop 10%, inserting and querying the layers 10%, bounds
    checks about 4% across them. A quarter of inexact matching sits in `first_reachable`'s short skip
    loops, their samples on the instruction past the exit: the exit mispredicted, 0 to 2 steps a
    lookup, as the branch profile above found. Fetching the bucket bounds' lines a batch ahead and
    the entries a batch ahead of their tests, a two-stage pipeline, found the same matches and was
    1 to 2% slower there and neutral on the M2: the entries were not late. Local pruning has no hot
    spot, its time spread over the fronts' steps and `leftmost`'s lookups. Left.
- [x] **Traceback (about 0.9 s).** Retracing the final round's tiles from their recorded left
  edges. Since this was profiled, each tile is traced by a forward search over one window of
  diagonals, eight at a time, its recompute reusing its buffers (`6108666`, `5a0e514`, `38837f2`,
  `d2492b3`, `5e76145`): on the M2 (2026-10-09) `forward_segment` is 65 of 3,900 samples aligning
  genvar and 122 of 5,600 on ont-500k, about 2%, the recompute too rare to show.
- [x] **The inexact seeds' neighbour windows (12 to 21% on the long reads).** An exact match's windows a
  base shorter and longer were left out of the inexact seeds, and without them the heuristic overestimated,
  by an edit a seed on crafted pairs; keeping them (2026-10-09) made genvar and ont-500k 12% slower on the
  M2 and 12 and 21% on the i9-7900X. They are always real matches, so the scan leaves them out again and an
  exact match, once kept, puts its four into the layers itself, with no scan and no search of their own;
  local pruning drops none that its exact match passes (2026-10-10). The M2 aligns ont-500k and genvar 4
  and 6% faster than with them scanned, 100 kbp pairs at 10 to 15% 4 to 6%. Checking the bound at every
  cell of random pairs against the true cost found it still over by one edit on optimal paths of 1 to 2%
  of pairs: an exact match whose end lay a diagonal past the pair's end was dropped, and local pruning
  dropped the neighbour that stood in for it. Such a match now chains to the end for its one edit of
  excess (`chained_layer`), and 2,200 random pairs, both seed kinds, read under the cost at every cell.
  No distance had come out wrong: an overestimate of one costs the band a round, not the answer.
  On the i9-7900X (`ebe290f`) ont-500k and genvar align in 104 and 144 ms where they took 112 and 165,
  and 100 kbp pairs at 9 to 15% 8 to 10% faster, A\*PA2-full now ahead at 10% alone (10.5 against 11 ms).
- [x] **The M2's two lost rows.** 100 kbp pairs at 6 and 7% divergence lost to A*PA2-full there, the
  M2's 40% cutoff rebuilding inexact seeds that did not pay on spread errors. A pair whose two
  projections agree, as spread errors' do, now rebuilds only below 20%: 3.75 and 3.92 ms against
  A*PA2-full's 5.25 and 5.0, at a cost of 1 to 2.5% on the M2's real long reads, still well ahead.

## Consolidation

What the library computes twice, from four audits on 2026-10-08 (DP sweeps, cost types, batch plumbing,
traceback and text helpers). Each item shares one implementation, every copy's optimizations kept, its
answers held to the copies it replaces by the tests, and a hot path benchmarked before and after.

### Phase 1: copies that no longer agree

- [x] One batch core, `api.distances_in_lanes`, `alignments_in_lanes`, `alignment_from_lanes`, `each_pair` and
  `longest_first` over `Texts`, for the API's batches and the C API's: C now deals pairs longest first and
  takes the same memory budget. Bytes 0xFE and 0xFF stay C's to refuse: a `String` never holds them.
- [x] `Scoring` alignments honour `max_memory` (`fronts_within`), not `HISTORY_LIMIT`.
- [x] The GPU distances' host fallback over the threads asked for. Left: `search` copying the query a
  reference, about a millisecond for 20,000, which sharing would need a per-pair path over any texts for.
- [x] Latent: `cigar_of`'s room with a free mismatch, `band_flags` reading a table, one `looked_up` for
  `byte_lookup` (which fell short under 16 lanes), `joined_cigar` merging neighbours.

### Phase 2: shared leaf helpers

- [x] Byte reversal: `cigar.reverse_bytes`, `reversed_into`, `reversed_list`, `append_reversed`,
  `reversed_text` in place of six copies.
- [x] CIGARs: one `CigarWriter` (in `cigar`) for `GappedAlignment.cigar`, `joined_cigar` and the tracebacks;
  one walk, `cigar_counts`, for `cigar_matches`, `Alignment.counts` and `costs_of_cigar`.
- [x] Move codes: `gap_affine`'s defined from `bit_parallel`'s, the lanes' mask and layers from
  `gap_affine`'s (`layer_bit`, `gap_layer`).
- [x] Constants in `common`: `UNREACHED`, `NEGATIVE_INFINITY` for `NOWHERE`, the sentinel bytes; `Band.covers_any`,
  `Band.clamped` and `ANY_LENGTH` for the bindings' and the API's own spellings.
- [x] `common.code_table` for five tables, `raise_unknown` for three loops, `EndsFree.of` for three clamps,
  `Alignment.mirrored`, `matches_along`, `AffineGapCosts.run`, `text_of`, `StringTexts.of` throughout.
- [x] Tests: `reversed_text` from the library. Kept apart on purpose: the fuzzer's and the tests'
  `reversed_cigar`, `rescore`, `priced`, `cigar_of_moves`, which check the library rather than share it.

### Phase 3: one cost model

- [x] `affine_penalties`, `affine2p_penalties`, `extension_penalties` through the `Costs` factories and
  `penalties_of`/`rewarded_penalties`, and since removed, the tests calling those directly; one uniform-table check (`substitutions.uniform_pair`); `Scoring.penalties`
  for five calls; `LaneCosts.of_penalties` for three; `table_extremes` for the table scans.
- [x] One 16-bit fit rule, `scored.fits_16_bits`, behind `narrow_enough` and `table_fits`; `Costs.dearest_step`
  and `Costs.cheapest_extension`, which the sweep's Z-drop now takes as the wavefront's does.
- Kept: `Scoring.gaps` an `AffineGapCosts`, a public field; `scaled_penalties`' raw positivity and
  `wavefront_penalties`' folded one, two rules on purpose: a table may score a mismatch above zero.

### Phase 4: engines

- [x] One span rule, `scored.started_span`, for `Costs`' and `Scoring`'s free ends; one sweep dispatcher,
  `scored.sweep`, for `swept`, `local_scores` and `tabulated_end`; one front scan, `gap_affine.front_best`,
  for the extension's three.
- [x] The cost lanes 64 bytes a group, several registers where the CPU's are narrower: on the M2 1 kbp
  distances 1.8 times as fast.
- [x] One Gotoh recurrence, `anti_diagonals.gotoh_lanes`, for `gotoh_cell` and the three vector sweeps.
- [x] One anti-diagonal sweep, `anti_diagonals`: its letters, its rolling diagonals in registers
  (`DiagonalCells`), its step and its 16-bit rule, for `swept_cells`, `reach_back` and `vector_sweep_bands`,
  each keeping its borders and what it watches for; `vector_align`, which keeps its band whole, shares the
  letters and the recurrence. `reach_back` and the halves now take 16-bit lanes while they fit: on the M2 a
  `Scoring`'s 1 kbp local alignment 1.2 times as fast, its linear-space global one 1.4 times. The score
  sweeps came out 2 to 5% slower on the M2 at 1 kbp, in 32-bit lanes, their inner loops instruction for
  instruction the old ones.
- [x] One row fill and one walk, now in `gotoh.mojo`, `fill_rows` and `walk`, for `serial_align`,
  `solve_rectangle`, `sweep_bands` (over two rolling rows) and `reconstruct`; read unchecked, the
  linear-space global alignment 1.2 to 1.4 times as fast.
- [x] One source rule, `diagonal.furthest_source`, for `best_source` and `grow_to`; one checkpoint rule,
  `band.checkpoints_passed` and `check_margin`, and one score down a run of words, `Frontier.climb`, for
  `Band` and `edit_search`.
- Kept apart: `step_front` from `grow_to`, which grows only the diagonals the backward fronts keep and
  slides only those, and from `forward_segment`, whose substitutions stop at the tile's end; `edit_search`
  from `Band`, a free start, the last row read every column and a bound by words rather than rows.

## Paper

- [x] **A controlled ablation.** Build switches turn each technique off alone (`dinara_align/ablation.mojo`,
  `-D ABLATE_...=1`): the regrouped recurrence, the second AVX-512 group, gathered slides, the agreement
  rule, neighbours put in from their exact match, and the batch's certified band. Every switched build
  passes the test suite and returns the same costs on every dataset. `pixi run bench-ablation` builds all
  of them and times them alternately, pinned. Plain builds moved the smaller effects by several percent
  from one set of builds to the next (4 and 6% against 1 and 3% for the recurrence and the pairing), the
  Skylake-X's penalty on jumps crossing 32-byte boundaries; the runner now assembles Mojo's own assembly
  with `clang -mbranches-within-32B-boundaries` and links it as `mojo build` does (`mitigated_build`), and
  a baseline built so from another source file ran within 0.7% of it on every sample. So built
  (2026-10-10, three rounds, the baseline within 1.2% between rounds): switching off the neighbours cost
  the long reads 9 and 16% and 100 kbp pairs at 10 to 15% 9 to 10%; the gathers 13 to 24% on the 1 kbp
  reads, SARS-CoV-2 and pairs at 1%; the agreement rule 18% on ont-10k and 6% on SARS-CoV-2, nothing on
  the long reads, where the seeds set the first bound; the paired groups 4 to 5% on the long reads;
  the regrouped recurrence 1 to 2% (the runner's own builds, one round, agreed within a point), its 14 to 8 cycles a word mattering only to the few narrow bands; the
  certified band 2.9 and 3.1 times on the short-read batch (its harness, carrying GPU kernels, built as
  usual). In the paper's Results ("Effect of each technique").

## From the paper's literature check

Ideas from the work the paper cites, each checked (2026-10-10, M2 unless noted).

- [-] **A batch band certified by its edges' costs** (SeedEx's check of a narrow band's boundary,
  Fujiki et al., MICRO 2020). Of 500,000 short-read pairs 70,519 fail the band's certificate at affine
  costs (1, 2, 1) and 49,001 at unit costs, and the second pass takes 58 of 192 ms and 43 of 174 ms;
  25,004 and 3,335 of the failing pairs had their optimum inside the first band, so a tighter bound
  could have spared them. Kept the least cost of a step off each edge of the band during the sweep
  and bounded every path off it by that and a gap back to the end diagonal: no pair more was proven,
  and the sweep took 4 to 6% longer. In costs the bound cannot beat `off_band`: the path that leaves at
  the origin, along row zero or down column zero, is one of the edge steps and costs exactly
  `off_band`'s opening and extensions. SeedEx's check gains in scores, where a path that leaves early
  forgoes its matches' rewards; in costs a closer bound needs a lower bound on what remains past the
  step, which costs about what the second pass does. Left.
- [-] **Difference recurrences for the sweeps under a table** (Suzuki and Kasahara 2018, as KSW2
  runs them), for 8-bit lanes regardless of length. They take global and extension alignment only: a
  local alignment's clamp at zero has no form in differences, and local alignment is where SSW comes
  closest (1.3 to 2.3 times slower). Global alignment under a table already runs 5.4 times parasail's
  speed. Not built.
- [-] **8-bit saturating lanes for the local sweep** (SSW's and parasail's first pass). Exact, since a
  local score never goes below zero, and twice the cells a register, but only while every score stays
  under 255: at a match score of 2, alignments under about 125 bases. A 150 bp read already scores about
  300, and the local workloads score in the thousands, so each would saturate and sweep again in 16
  bits. Not built.
- [-] **A striped bit-parallel layout** (BSAlign, Shao and Ruan 2024), which might take the lane
  rotation off each column's chain. Its edit-distance mode is described only in its supplement; it
  reports 2.1 times Edlib's speed over 1 to 100 kbp, where dinara-align's unit-cost engine runs 6.7
  times Edlib's speed on ont-10k and 22 times on 100 kbp pairs at 15%. A striped layout also needs a
  correction pass for the bits crossing between segments, each column, the chain it means to shorten.
  A*PA2 weighed it and left it too. Not built.
- [-] **The aligners since A*PA2** (surveyed 2026-10-10). QuickEd (Doblas et al., Bioinformatics 2025):
  an upper bound from overlapping windows of bit-parallel Myers (128 and 640 wide), then one tiled pass
  with Ukkonen's cutoff; its bound runs 1 to 1.5% over on ultra-long and PromethION reads but takes
  about a quarter of QuickEd's own time, and QuickEd itself runs about as fast as A*PA2 (1.03 times
  A*PA2-simple, 0.94 times A*PA2-full). A bound only spares the rounds that fail, at most 16% of the
  band's word-columns on genvar and 7% on ont-500k with the bound free; see "An upper bound on the
  distance" below. Sassy and Sassy2 (Beeloo and Groot Koerkamp, 2025 to 2026): bit-vectors along the text,
  four text chunks in the SIMD lanes, stopping once every cell passes k; for searching patterns up to
  about 1,000 bp in long texts, a search mode dinara-align does not have. TALCO (Walia et al., HPCA 2024)
  and miniwfa (Li) save traceback memory, not time; the GPU, FPGA and in-memory designs (WFA-GPU,
  eWFA, GeneTEK, Scrooge, GenASM, RAPIDx) carry no idea to the CPU paths beyond what is here. Nothing
  promised a gain past the builds' noise. Left.
- [x] **The infix search's retries aimed, as the global band's are.** Measured beside Sassy (2026-10-10,
  M2): a free-start search failing its first bound of 64 only doubled, a 1 kbp read at 10% taking two
  sweeps and a 10 kbp read at 5% four or five, where one at the right bound costs half (35 against 67
  µs, 1.29 against 2.88 ms). A failed try now projects the distance from the deepest row its band kept
  within the bound, edits spread along the pattern putting that row near `bound * rows / distance`, and
  the retry aims past it; the sweep for the start, which knows the distance, begins there. The 10 kbp
  infix aligns in 1.99 ms where it took 3.39, the 1 kbp one and the prefix mode 3 to 4% faster; tests,
  400,000 fuzz cases and WFA2-lib's regression set agree.
- [ ] **Short patterns in long texts, the text split across lanes** (Sassy, Beeloo and Groot Koerkamp,
  Bioinformatics 2026). Sassy itself is slower here than dinara-align's infix (M2, search only): a 1
  kbp read at 10% in 3 kbp 204 µs with k doubling and 121 with k given, against 80; a 10 kbp read at 5%
  in 30 kbp 13.4 and 7.8 ms against 1.8. Its speed is a small k: given the distance, a 30 bp pattern
  at 5% in 100 kbp takes 68 µs and a 100 bp one in 10 kbp 12, where dinara-align takes 667 and 97, and
  with k doubled from 64 Sassy takes 775 and 99. A pattern of one word runs the infix sweep's scalar
  last-word step a column at a time, about 6.7 ns a column. Sassy's idea that carries over: the same
  pattern word in every lane, each lane a stretch of the text overlapping the next by twice the
  pattern's length so every match lies whole in one, exact; with the first bound from the pattern's
  length rather than 64. Not built.
- [-] **The other aligners since 2023, read for ideas** (2026-10-10). Sassy2 (2026) puts many short
  patterns in SIMD lanes behind a suffix filter: at its reported 6 Gbp/s a pattern a thread, the inexact
  seeds' 40,000 16-mers against a 650 kbp read would take about 4 s, where the hashed halves take about 30
  ms. SeqMatcher (J. Supercomputing 2025) packs sequences to two bits with AVX-512: building the profile
  is 1.5% of a SARS-CoV-2 alignment and does not show on the 1 kbp reads (perf, Skylake-X), so packing has
  nothing to win; its banded mode is a fixed threshold, 3 to 10% less accurate. BSAlign's striped sweep
  with an active F loop claims 2 times other SIMD aligners; measured since, see "Rivals' implementations"
  below. FILTR (2026) chooses a wavefront or an
  anti-diagonal schedule by divergence, as the quarter-matrix rule does; Medlib's (2025) threshold mode
  is `max_cost`'s early stop; hashed longest-common-extension queries (Ding et al., ESA 2023) are
  probabilistic, so not taken; composition and q-gram lower bounds (certified-alignment) are far looser
  than the seeds' bound at the origin; Theseus (2026) is sequence-to-graph. Left.
- [-] **Incremental doubling** (A*PA2, Section 3.8): a round after the first skips the rows the last
  round fixed, from horizontal differences stored along a fixed row in each tile. A*PA2 reports it only
  with three other methods, together 3 times faster. Measured what it could skip (2026-10-10, M2): the
  word-columns a retried round sweeps inside the last round's kept rows at both edges of a tile are 8.8%
  of the band's on genvar (58 of 126 rounds reuse any), 11.0% on ont-500k, 0.6% on 100 kbp pairs at 15%
  and none at 10%: with the band about 65% of a long pair's time, at most 6 to 7% there before the cost
  of storing every tile's edge and a row inside it, nothing elsewhere. Left, for robustness first: the
  skipped block's bottom row must be exact at every column of the tile, not only at the edges where
  pruning reads the scores, and A*PA2 asserts it without a proof that carries over to these tiles; were
  it wrong, a round could accept a distance too high. It would also change every bit-parallel kernel's
  top word, which reads `+1` from above. Worth taking up only with that proof in hand.

## Rivals' implementations of ideas dinara-align has, measured

Each rival built at a pinned commit and run as the benchmarks run every aligner (2026-10-10, Skylake-X
unless noted; BSAlign and QuickEd are x86-only).

- [x] **BSAlign** (Shao and Ruan, Bioinformatics 2024), its striped kernels with the band off, so exact:
  now a column of `pa_bench.py`. Its edit distance fills the whole matrix, 5 times slower than
  dinara-align on ont-1k, 8 on ont-10k and 170 on SARS-CoV-2; its affine alignment (4, 6, 2), 8-bit
  difference recurrences with an active F loop, 3 to 4.6 times slower on the ONT sets and 3 times faster
  than KSW2. Its kernel fills about 0.38 ns a cell with the traceback, where the Gotoh sweep under a
  table took 0.55 ns for the score alone in 32-bit lanes, which led to four changes, each holding CIGARs
  byte for byte: the 16-bit bound charges the straight path, not every letter the dearest move; the
  unreachable cell sits 4096 above 16 bits' least, not a quarter of the way down; the banded traceback
  keeps each cell's decision, a byte, not three 32-bit scores; and the band proves itself from its own
  score rather than after a sweep of the whole matrix for it. table-global went from 544 to 253 us on
  the Skylake-X with the first and third, and to about 157 on the M2 with all four.
- [x] **QuickEd** (Doblas et al., Bioinformatics 2025), bound-and-align beside the band doubling: on its
  own simulated sets, 10,000 pairs of 10 kbp take dinara-align 0.3, 2.9 and 3.7 s at 1, 5 and 10% where
  QuickEd takes 6.1, 12.6 and 19.9; see `quicked_bench.py`.
- [x] **Parallel output-sensitive edit distance** (Ding et al., ESA 2023), its BFS-Hash and BFS-SA on one
  thread, distance only against dinara-align's with its traceback: 40 to 150 times slower on 100 kbp
  and 1 Mbp pairs at 1 and 5%, its design point, and about 100 times on divergent pairs.
- [x] **Divergent pairs at affine costs.** The wavefront's work grows with the cost times the length,
  so where BSAlign's 8-bit matrix took 1.5 ms for unrelated 2 kbp pairs the wavefront took 4.6. Now the
  wavefront hands a pair over once its projected work passes a full sweep in difference recurrences
  (Suzuki and Kasahara 2018), 8-bit lanes at any length, that keeps each cell's flag by the lanes' rule
  for the wavefront's ties, so the alignment is the same byte for byte (see `differences`). Its traceback
  keeps every few anti-diagonals and sweeps each stretch again, memory the square root of the matrix's:
  a byte a cell faulted in fresh pages that cost more than the sweep. Each diagonal is swept in place
  from the last row up, which kept 5 and 10 kbp pairs at 0.18 to 0.24 ns a cell where two diagonals took
  0.36 and 0.45. On the Skylake-X, `align` at (4, 6, 2): unrelated 1, 2 and 5 kbp pairs 1.50, 5.3 and 32.9
  ms to 0.28, 0.87 and 4.8; at 30% 0.57, 2.2 and 16.4 to 0.28, 0.88 and 4.7; `distance` 2 to 14 times
  faster; pairs under 15% apart unchanged. The exchange rate is `CELLS_PER_STEP`, measured.

## From A*PA2's discussion

Its limitations and future work (curiouscoding.nl/posts/astarpa2/#discussion), and where dinara-align
stands on each.

- [x] **Symbols beyond `ACGT`.** The edit distance takes `ACGT` and up to four other bytes, `N`
  among them, each matching only itself, through a third bit plane its own kernels read, at no cost
  to bases alone. Such a pair now takes seeds too: a seed holding an `N` goes uncounted, matching
  nowhere and charging nothing, and the second sequence reads its `N` as a base for the seeds alone,
  which only adds matches. With one `N` a read, scattered ones or a 100- or 5000-base run in both,
  the 500 kbp reads align 2.5 to 3 times faster on the M2 and 1.4 to 1.7 times on the Skylake-X,
  and 100 kbp pairs up to twice as fast on the M2, at the same costs and bases alone unchanged
  within noise. Folding the first sequence's `N` too was as fast on short runs but flooded a long
  run with matches, a run of `A` against a run of `A`, up to 800 times slower.
- [x] **Seeds on ont-50k's reads with an `N`, on x86.** They cost about 6% there, and not for the
  `N`: on AVX-512 seeds lost on every pair below about 85 kbp, real reads of 16 to 64 kbp aligning
  18% faster without them and uniform 30 kbp pairs at 5% 38% faster. Seeds now start at 86 kbp on AVX-512 and nowhere
  shorter on the projection; AVX2 broke even on those reads and the M2 gained 9%, so both keep the
  16 kbp gate. Exact seeds' matching also checks a filter of 32 bits a seed before its half-full hash
  table, so nearly every window that matches nothing costs one predictable test: matching took 87 ms
  of the M2's 320 mid-length reads and takes 32. On the Skylake-X, mid-length reads 17% faster, with
  an `N` 14%, 30 kbp pairs 18 to 40%, 100 kbp at 5% 10%, ont-500k 3%.
- [x] **The cost of an `N`.** ont-50k's reads with one `N` took 12% longer than without on the
  Skylake-X, and not mostly in the sweep: coding such a pair went byte by byte, five times as slow,
  and its third plane was built a bit at a time, four times the planes' cost. Both now go sixteen
  bytes or eight rows at a time, only chunks holding a symbol past `ACGT` byte by byte, and each tile
  matches on no more planes than its own symbols need (see `Profile.symbols`): two where neither its
  columns nor any row holds one, the rows' third plane as a mask where only rows do. Reads with an
  `N` 7 to 9% faster on both machines, 5 to 7% over clean ones on ont-50k on the Skylake-X.
- [x] **Semi-global and open-ended alignment.** Every cost model takes every mode (`Mode`): global, a
  query inside a reference (Edlib's infix), at its start or end, any free ends as WFA's ends-free mode
  allows, and an extension from either end; the CIGAR spans the aligned parts, their bounds beside it.
  At unit costs the infix and prefix modes sweep a band as Edlib does, a bound doubled and only the
  rows some score within it can still reach (Ukkonen's cutoff): 10 kbp in 100 kbp in 13 ms on the M2,
  100 kbp in 1 Mbp in 0.75 s, against 21 ms and 2 s for the whole matrix; the other modes, a band or a
  cap take the wavefront at the same costs. Still open: a first bound from an estimate instead of 64,
  and a seed heuristic for the search.
- [-] **Traceback while sweeping, to save memory** (as TALCO does). Feasible, and left (2026-10-06):
  every few tiles the traces from the band's top and bottom kept rows, deterministic and unable to
  cross without meeting, would settle the alignment behind where they meet and free the tile edges
  there. But the edges are the smaller part of an alignment's memory: 0.6 to 13 MB a pair against 10 to
  17 MB of sequences and bit planes (1.1 to 13.3 MB on ont-500k, 0.6 at 1 Mbp and 5%, 13 at 1 Mbp and
  15%), and the convergence traces would cost a share of the traceback's 5 to 9% again. The larger
  items: the first sequence's two 64-bit masks a base, 16 bytes, 17 MB at 1 Mbp, which the kernel could
  expand from the codes at some cost to the sweep. The peak is set during the band, by the planes,
  the tile edges and the seeds' layers together: the diagonal transition's fronts are gone by then,
  as Mojo destroys a value after its last use, and the moves copied reversed to write a CIGAR come
  after the edges are freed, so writing them in place left the peak where it was.
- [-] **A\* on the diagonal transition** at low divergence. dinara-align runs a plain diagonal
  transition first and falls back to the band; it is fast below 2% but uses no heuristic there.
  Tried and left (2026-10-09, M2): a front trimmed at both ends wherever its score plus the exact
  seeds' heuristic at its furthest cell passes a bound, the bound from the origin's heuristic
  doubling its margin from 16, is exact (the heuristic never overestimates, and a cell further down
  a diagonal never costs more to finish from) and keeps the tie rule, since every cell an optimal
  path passes keeps its value. On uniform 100 kbp pairs at 2 and 3% it took 0.93 and 1.08 ms
  against the two-ended search's 2.02 and 2.27, but 0.9 against 0.55 at 1%, 1.8 against 2.0 at 5%
  and 8.6 against 3.4 at 8%, where the fronts widen and rounds repeat; and on real reads, whose
  errors gather, 3 to 5 times slower: ont-10k 0.30 against 0.08 ms, ont-50k 2.2 against 0.41,
  sars-cov-2 0.33 against 0.13. Two thirds of its time at 2% is the seeds' setup, 0.7 ms, and most
  of the rest the trimming's heuristic queries, about 40 ns each. A gate taking it only for uniform
  1.5 to 4% would rest on constants fitted to one machine, for synthetic pairs alone.
- [-] **An upper bound on the distance** to keep bounds from overshooting. Tried and left
  (2026-10-06): a beam, a band of 128 to 1024 rows kept around the lowest score down each tile's edge,
  whose corner is a real alignment's cost. It found the distance itself on most pairs, and one round
  at exactly the distance would save 16 to 34% of the time on long and mid-length reads on the M2,
  the failed rounds and the last round's overshoot. But the beam costs about a band round, so it pays
  only after a failed round and only with seeds; even then, on the Skylake-X it took 5% off ont-500k,
  added 1.5% to genvar, whose insertions it loses, and nothing elsewhere, while the M2 gained 4 to 12%
  on seeded sets. Four tuned parts for that was judged not worth it. A cheaper bound, or one that
  follows large indels, might still be.
- [-] **Affine costs with the seed heuristic.** The affine aligners are exact but sweep without a
  heuristic; A*PA2's method for affine costs needs a gap-chaining seed heuristic of its own. Weighed
  and left (2026-10-09): the affine aligners are wavefronts, so a heuristic would trim their fronts
  as above, and the unit-cost trial there paid only on uniform pairs at 2 to 3% and lost on every
  real dataset. An affine heuristic is also looser: a seed broken by one edit may have cost a
  mismatch, so it charges `min(x, o + e)` while the edit may cost more, and its fronts would trim
  less.
- [x] **Low divergence (below 2%)** was A*PA2's weak spot against BiWFA; the diagonal transition
  before any band now beats BiWFA and WFA there (146 µs against 302 and 605 at 0%, 100 kbp).
- [x] **The seeds' setup** was A*PA2's other limitation. Exact matching is filtered, seeds holding an
  `N` are handled, and AVX-512 builds skip seeds below 86 kbp where they do not pay, but it is still a
  third or more of a long read's alignment on the Skylake-X (see above), where every idea listed
  has now been tried.
