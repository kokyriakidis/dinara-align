# To do

`[ ]` open; `[-]` tried or weighed and left, with why, so it is not tried again. Finished work is in the
git history and `CHANGELOG.md`.

Where the time still goes in a single-threaded alignment, profiled on the Skylake-X (3.3 GHz,
pinned) over all 48 reads of ont-500k-genvar, about 10 s in all.

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

## Kept apart on purpose

Copies the consolidation of 2026-10-08 left separate, and why.

- Kept apart: the fuzzer's and the tests' `reversed_cigar`, `rescore`, `priced` and `cigar_of_moves`,
  which check the library rather than share its code.
- Kept: `Scoring.gaps` an `AffineGapCosts`, a public field; `scaled_penalties`' raw positivity and
  `wavefront_penalties`' folded one, two rules on purpose: a table may score a mismatch above zero.
- Kept apart: `step_front` from `grow_to`, which grows only the diagonals the backward fronts keep and
  slides only those, and from `forward_segment`, whose substitutions stop at the tile's end; `edit_search`
  from `Band`, a free start, the last row read every column and a bound by words rather than rows.

## Open

- [ ] **ont-500k-genvar and the inexact seeds.** With no seeds matching within one edit below 1 Mbp, the
  genvar reads aligned 25% faster (146 to 110 ms) and ont-500k's distance 30% slower: `INEXACT_CHAINED`
  rebuilds the seeds inexact on reads where the exact ones would win. Find what tells the two apart
  before touching the gate; a cutoff moved for one dataset is a tuned knob.
- [ ] **The batch on HiFi reads at gap-affine costs.** `alignments` takes 806 µs a pair on the pacbio-hifi
  sample at (4, 6, 2), where one `Aligner` takes 95 µs a pair; at unit costs the batch keeps pace (16.5
  against 12.6 µs). Profile where the batch's time goes on 14 kbp pairs before changing anything.
- [ ] **Gap-affine pairs costing 33 to 64 on short reads.** On the Illumina sample at (4, 6, 2) these 40
  pairs take 10.6 µs to WFA's 9.6, where the first search hands over to two (see
  `gap_affine.traced_alone`); its model of the two searches leaves out their setup and the far side's
  grow-on, so it hands over a little early. Count those in from what the searches measure, not a factor.
- [ ] **The gap-affine step's instructions a cycle on long divergent pairs.** On 100 kbp pairs at 20%
  (2026-10-14) a pair steps 2.26 billion diagonals in 14.9 billion instructions, 4.7 times fewer than BiWFA's
  70.4 billion, at 1.28 instructions a cycle to BiWFA's 2.41. Top-down: core-bound 40%, memory-bound 20%; the
  step's two gathers (`slides.gathered_slides`) take 38% of its samples. Its L1 misses (precisely sampled,
  nine in ten on its source rows' loads, one in eighty past L2) are L2 hits the hardware already streams:
  fetching the rows or the letters ahead lost (below). The gathers are the cheapest slide this CPU has: with
  `-D ABLATE_GATHER`, a scalar compare a lane, 100 kbp pairs at 20% took 7.17 s to 3.56, at 5% 655 ms to 327,
  10 kbp at 10% 12.3 ms to 7.0, ont-10k 2.65 ms to 1.58. What is left is fewer gathers, not cheaper ones:
  a lane group whose lanes all mismatch at their first letter gathers nothing more than one load would show.
- [ ] **Which nanopore reads the diagonal transition should finish.** Regret (`profile.py regret`, 2026-10-12)
  on ont-10k: run to its end, the diagonal transition beats the default on 114 of 277 reads, 3.3 of 38.7 ms,
  nearly all under 8% divergence (4.8 against 6.5 ms there); forced on every read it takes 1.53 times as
  long. Two signals tried and left below (the search's own projection, the projection from recent
  progress) did not pick those reads out. Compare per-pair probes of the winners and the losers before
  writing a third rule.

## Tried in the optimization rounds of 2026-10-12 to 14

Each measured on the Skylake-X with `profile.py compare` and `regret`; none kept.

- [-] **Aiming the band's first bound at the search's own projection** (`-D DINARA_FORCE=aim-search` still
  times it). The two-ended search's projection, from every edit it reached, lies far closer to the
  distance than the one front's from 8 edits (median 0.98 against 0.85 of it on 10 kbp pairs at 10%,
  1.5 against 3.9 on ont-10k), and forced on every pair it saved 3 to 8% on most sets, but cost 4% on
  ont-10k distances: where it runs 2 to 4 times high the band starts that much too wide, while a bound
  too low costs only a cheap failed round. Safe only below the projection by a factor fitted to the
  sample, a tuned knob. Left.
- [-] **A CIGAR entry's digits packed into a word and stored at once** (5a500f8, reverted f2dbca5).
  `CigarWriter.flush` then stayed out of line and its four divisions cost a third more instructions on
  close pairs, whose long runs the digit loop predicts well: 10 kbp pairs at 1% 1.5% slower, nothing
  gained on divergent ones. The mismatch-mask split of diagonal runs from the same change was kept.
- [-] **Charging a recording search for the history memory it has not held before** (fa34391, reverted
  98e5724). Measured on one 50 kbp PromethION read, its two-ended search spent 5 of 8.7 ms faulting in
  11 MB of fresh pages where the band took 1.95 ms, and the charge fixed that read. But faults are paid
  once an aligner, not once a pair: a fresh aligner holds no history, so the first search is charged and
  gives up, the memory never grows, and every later pair is charged too; 10 kbp pairs at 5% went to the
  band and took 23% longer. Model only costs every pair pays. That read's real cause is the alignment's
  band model (`ALIGNMENT_COST`, fitted on uniform pairs) overestimating its band about twice.
- [-] **Prefetching the letters a gap-affine lane group two on will slide over** (bfb06c0 and b292774,
  reverted), from its front a mismatch back, two prefetches a group. At every width: long divergent pairs at
  (4, 6, 2) 3 to 4% faster (100 kbp at 5%, 1 Mbp at 1%, PromethION) but ont-1k 1.9% and Illumina 0.6%
  slower. Only when a cost's rings outgrow L1 (32 KB, 32 bytes a diagonal): pairs from 100 kbp still 2 to
  4% faster, but 10 kbp pairs at 5 to 20% 1.1 to 2.4% and ont-10k 1.6% slower, their fronts past L1 too. A
  width that keeps only the gains would be fitted to these sets. Left.
- [-] **Prefetching the gap-affine step's source rows** (da955aa, reverted), five rows 64 diagonals ahead
  when a cost's rings outgrow L1, the loads precise sampling put nine in ten of its L1 misses on: 10 kbp
  pairs at 5 to 20% 3.6 to 7.6% slower, 100 kbp at 1% 6%, ont-10k 3.6%, only 100 kbp at 20% 2.4% faster.
  The hardware prefetchers already stream the rows from L2; five more fetches a group only add to the
  loop. Left.
- [-] **Giving up only when the projection from the second half of the edits agrees** (1db41cc, reverted
  a495e85). Meant for nanopore reads whose errors gather at their ends, where the two-ended fronts start:
  no change on ont-10k or ont-50k, so their winners are not told apart that way, and recording the progress
  at every check made short pairs' distances 2 to 4% slower (Illumina, HiFi, 10 kbp at 1%). Left.

## Tighter bounds: the audit of 2026-10-11

Every bound checked against A*PA2's code and paper, QuickEd, Edlib, KSW2, BSAlign, WFA2-lib, SeedEx,
Gibrat 2018 and POASTA, and its ceiling measured on the Skylake-X with an oracle (the pair's own
distance) in a scratch build, before anything is built:

- At parity, nothing to gain: BiWFA's stop (`wavefront_bialign.c`: F + G >= best + o + w, ours exactly;
  checking against every kept front instead of the ring cannot tighten it, costs jumping by w a move);
  A*PA2's first bound and growth, kept rows, tile top and bottom, GCSH, pre-pruning, tile width; the
  O(NP) band shape of Edlib and QuickEd. Band-edge DP values never beat Proposition 1's L(d): the edge
  cell on row 0 is a pure gap.
- [ ] **The band's first bound.** A first bound equal to the distance would save 27% on 100 kbp at 20%,
  22% ont-50k, 19% ont-500k, 17% 10 kbp at 20%, 14% 100 kbp at 10%, 13% genvar, 5-7% ont-10k, HiFi,
  1 Mbp at 5%. Failed rounds take 18-36% of the time on those (checkpoint lowering causes little: 6% on
  ont-10k, ~0 elsewhere); on 1 Mbp at 20% no round fails and the loss is a first bound aimed high. A
  loose bound is not robust: starting at 1.15 d loses 18-23% on the 500 kbp reads. Robust use is an
  upper bound U as a cap (`best = U`, so a retry at U must succeed and no round widens): with U = 1.05 d
  it saves 27% on 1 Mbp at 20%, 12% ont-50k, 7% 100 kbp at 20%, 6% ont-10k, 5% HiFi, ~0 elsewhere; with
  1.15 d about half of that or none. QuickEd's windowed bound deviates 1-6% on ONT, 5% HiFi, up to 15%
  (its paper's Table 1; to check). Worth it only for a U within ~5% that costs far less than it saves.
  A beam as U was tried and left (see "An upper bound on the distance").
- [ ] **Gap-heuristic pruning of the near-identical diagonal transition** (POASTA, A*PA's `dt.rs`): with
  U = d, keeping only |k - delta| <= U - s saves 24% on 10 kbp at 1%, 13% ont-1k, 8% HiFi. Without U it
  needs A* order, whose fronts depend on their own neighbours towards delta and would serialize the
  gathered extension; the same holds for the one-search gap-affine phase. Find an exact U or a vector
  form first.
- [ ] **Incremental doubling** (A*PA2-full): ceiling 6-7% on long ONT reads, ~0 elsewhere; blocked on an
  exactness proof (see "From the paper's literature check").
- [-] **The batch certificate, strengthened by edit distance** (after SeedEx's edit-distance check): a
  path off the band has an insertion gap and a deletion gap, and its mismatches and gapped letters are at
  least ED, so it costs at least max(L, 2o + min(x, e) ED). Measured and dropped: of 500,000 Illumina
  pairs at (1, 2, 1), 70,519 are swept again, 59 of about 300 ms; the bound proves 578 of them. 65% of
  them hold a dearer cost in the first band than their least and need the wider band whatever proves
  it; the other 35% have their least cost already, but more than two gaps, so their cost passes ED plus
  two openings. A perfect certificate would save at most about 7% of the batch.
- [ ] Small and exact: the score-based storage band can subtract the two gaps' opening surcharge
  (`vector_score.optimal_band`); Edlib's guaranteed mid-round bound, score + max(rows, columns left).

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
  with an active F loop claims 2 times other SIMD aligners; measured since as a benchmark column, its
  affine alignment 3 to 4.6 times slower than dinara-align's on the ONT sets and its edit distance 5 to 170
  times. FILTR (2026) chooses a wavefront or an
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

## Memory locality

Sweeps whose working set outgrew the first-level cache, or whose traceback faulted in fresh pages, found by
the difference sweep's case (2026-10-10, Skylake-X).

- [-] **The wavefront's kept fronts, fresh memory each call.** `align` builds its searches anew each call,
  as the library keeps no state of its own, so its kept fronts come from the allocator each time. Measured
  on fresh 5 kbp pairs 15% apart from a fresh process: the first three calls 8.4, 6.2 and 4.3 ms, then 3.5,
  as much as one `Aligner` reused takes; the allocator keeps the blocks once it has given them out, so it is
  a process's warm-up, not a call's cost. Left.

## From A*PA2's discussion

Its limitations and future work (curiouscoding.nl/posts/astarpa2/#discussion), and where dinara-align
stands on each.

- [ ] **The infix and prefix searches' first bound and seeds.** At unit costs they sweep a band as Edlib
  does, a bound doubled from 64 with Ukkonen's cutoff; still open: a first bound from an estimate instead
  of 64, and a seed heuristic for the search.
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
