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
- [x] **Traceback (about 0.9 s).** Retracing the final round's tiles from their recorded left
  edges. Since this was profiled, each tile is traced by a forward search over one window of
  diagonals, eight at a time, its recompute reusing its buffers (`6108666`, `5a0e514`, `38837f2`,
  `d2492b3`, `5e76145`): on the M2 (2026-10-09) `forward_segment` is 65 of 3,900 samples aligning
  genvar and 122 of 5,600 on ont-500k, about 2%, the recompute too rare to show.
- [ ] **The inexact seeds' neighbour windows (12 to 21% on the long reads).** An exact match's windows a base
  shorter and longer were left out of the inexact seeds, and without them the heuristic overestimated,
  by an edit a seed on crafted pairs, though no distance came out wrong; keeping them (2026-10-09) made
  genvar and ont-500k 12% slower on the M2, and 100 kbp pairs at 15% 8%, every cost the same; on
  the i9-7900X (2026-10-09) ont-500k 12% and genvar 21% (100 and 136 ms to 112 and 165), and the
  100 kbp pairs at 9 to 15% 10 to 20%, A\*PA2-full now ahead at 10 and 12%. A
  neighbour only matters where it chains on and the exact match cannot, a diagonal off, so most could
  be dropped again by a test that keeps the bound exact.
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
