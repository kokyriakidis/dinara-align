# To do

Where the time still goes in a single-threaded alignment, profiled on the Skylake-X (3.3 GHz,
pinned) over all 48 reads of ont-500k-genvar, about 10 s in all.

- [~] **Band rounds that fail (about 1 s of the band's 6.5 s).** A retry now aims a quarter of the
  estimated climb past its estimate, not half: genvar 2% and ont-500k 7.5% faster on the Skylake-X.
  Lowering later rounds' bounds at checkpoints, as the first round does, never ended on some reads. The sweep itself runs near the
  hardware's limit for Myers' recurrence, about 2.4 cycles a word-column on AVX-512, so the band
  only gets faster by sweeping fewer cells. On genvar 16% of the band's word-columns go to rounds
  that fail and are retried (7% on ont-500k), three to six rounds a pair.
- [x] **Wider vectors on AVX-512.** A sweep there now runs two eight-lane groups at a time, one under
  the other in one loop: genvar 7% and ont-500k 5% faster on the Skylake-X, the sample sets within
  noise (measured best of five beside the owner's jobs). One sixteen-lane group did about as well on
  the long reads but cost up to 2% on short divergent pairs; on the M2 both were slower, so NEON and
  AVX2 keep one group.
- [ ] **Seed setup.** Still 34 to 45% of a long pair's time on the Skylake-X and 9 to 32% on the M2
  (2026-10-06, `197656a`): on genvar 157 ms of inexact matching, 94 of local pruning and 19 of layers
  in 698; on 100 kbp pairs at 15% 27, 16 and 1 in 104; on ont-500k 64, 69 and 10 in 468.
  - Inexact matching: the half tables' lookups cost 5.5 ns a row, and the rest is the candidates,
    0.74 a row on genvar, 17% of them within one edit: about 15 ns a row plus 43 a candidate there,
    against 3 and 28 on the M2. A chance candidate shares an 8-base half with the window, likely
    with 41,000 seeds a read over 65,536 halves. Keying on two exact parts, a 12-base prefix and a
    quarter at the end, would cut them some thirtyfold for about eight hashed lookups a row in place
    of two direct ones: perhaps 5 to 10% on genvar and divergent 100 kbp pairs, untried.
  - Tried and left: no local pruning (genvar 21% and ont-500k 72% slower: it repays itself many
    times), a lookahead of 8 seeds for inexact matches (neutral to 4% slower), and of 6 to 10 for all
    (mid-length reads up to 6% faster on the M2, ont-500k 17% and 100 kbp at 5% 30% slower), and
    scan batches of 32 rows (neutral).
  - Exact matching already checks a 32-bit-a-seed filter before its hash table (87 to 32 ms on the
    M2's mid-length reads); exact pruning, 600 thousand matches of which 174 thousand are kept, is
    now the larger half there. Going further than this means indexing the second sequence instead
    of the seeds, at more memory.
- [ ] **Traceback (about 0.9 s).** Retracing the final round's tiles from their recorded left
  edges.
- [x] **The M2's two lost rows.** 100 kbp pairs at 6 and 7% divergence lost to A*PA2-full there, the
  M2's 40% cutoff rebuilding inexact seeds that did not pay on spread errors. A pair whose two
  projections agree, as spread errors' do, now rebuilds only below 20%: 3.75 and 3.92 ms against
  A*PA2-full's 5.25 and 5.0, at a cost of 1 to 2.5% on the M2's real long reads, still well ahead.

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
- [~] **Semi-global and open-ended alignment.** `edit_search` finds a pattern inside a text, Edlib's
  infix mode, or at its start, its prefix mode, and `edit_search_alignment` aligns it there, exact.
  It sweeps a band as Edlib does, a bound doubled and only the rows some score within it can still
  reach (Ukkonen's cutoff): 10 kbp in 100 kbp in 13 ms on the M2, 100 kbp in 1 Mbp in 0.75 s, against
  21 ms and 2 s for the whole matrix. Unrelated text scores about half an edit a base when the start
  is free, so a bound of `k` still reaches some `2k` rows, which limits the cutoff. Still open: a first
  bound from an estimate instead of 64, a seed heuristic for the search, and a free end for the pattern
  too, as WFA's ends-free mode allows.
- [ ] **Traceback while sweeping, to save memory** (as TALCO does). Feasible, and left (2026-10-06):
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
- [ ] **A\* on the diagonal transition** at low divergence. dinara-align runs a plain diagonal
  transition first and falls back to the band; it is fast below 2% but uses no heuristic there.
- [ ] **An upper bound on the distance** to keep bounds from overshooting. Tried and left
  (2026-10-06): a beam, a band of 128 to 1024 rows kept around the lowest score down each tile's edge,
  whose corner is a real alignment's cost. It found the distance itself on most pairs, and one round
  at exactly the distance would save 16 to 34% of the time on long and mid-length reads on the M2,
  the failed rounds and the last round's overshoot. But the beam costs about a band round, so it pays
  only after a failed round and only with seeds; even then, on the Skylake-X it took 5% off ont-500k,
  added 1.5% to genvar, whose insertions it loses, and nothing elsewhere, while the M2 gained 4 to 12%
  on seeded sets. Four tuned parts for that was judged not worth it. A cheaper bound, or one that
  follows large indels, might still be.
- [ ] **Affine costs with the seed heuristic.** The affine aligners are exact but sweep without a
  heuristic; A*PA2's method for affine costs needs a gap-chaining seed heuristic of its own.
- [x] **Low divergence (below 2%)** was A*PA2's weak spot against BiWFA; the diagonal transition
  before any band now beats BiWFA and WFA there (146 µs against 302 and 605 at 0%, 100 kbp).
- [~] **The seeds' setup** was A*PA2's other limitation. Exact matching is filtered, seeds holding an
  `N` are handled, and AVX-512 builds skip seeds below 86 kbp where they do not pay, but it is still a
  third or more of a long read's alignment on the Skylake-X (see above).
