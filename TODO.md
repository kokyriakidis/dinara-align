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
- [ ] **Seed setup (about 2.5 s).** Inexact seeds: the scan of the second sequence (about 1.4 s
  when forced on every read), local pruning (about 0.8 s) and layer insertion (about 0.4 s).
  The scan is bound by loop overhead and memory, not one hot spot; going further means indexing
  the second sequence instead of the seeds, at more memory.
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
- [ ] **The third plane's cost on x86.** A read with one `N` sweeps about 13% slower than without on
  ont-50k on the Skylake-X. When only one sequence holds symbols past `ACGT` they match nothing, so a
  per-row or per-column mask could stand in for the third plane.
- [~] **Semi-global and open-ended alignment.** `edit_search` finds a pattern inside a text, Edlib's
  infix mode, or at its start, its prefix mode, and `edit_search_alignment` aligns it there, exact.
  It sweeps a band as Edlib does, a bound doubled and only the rows some score within it can still
  reach (Ukkonen's cutoff): 10 kbp in 100 kbp in 13 ms on the M2, 100 kbp in 1 Mbp in 0.75 s, against
  21 ms and 2 s for the whole matrix. Unrelated text scores about half an edit a base when the start
  is free, so a bound of `k` still reaches some `2k` rows, which limits the cutoff. Still open: a first
  bound from an estimate instead of 64, a seed heuristic for the search, and a free end for the pattern
  too, as WFA's ends-free mode allows.
- [ ] **Traceback while sweeping, to save memory** (as TALCO does). The band records every tile's left
  edge and traces back after; memory is level with A*PA2-full's.
- [ ] **A\* on the diagonal transition** at low divergence. dinara-align runs a plain diagonal
  transition first and falls back to the band; it is fast below 2% but uses no heuristic there.
- [ ] **An upper bound on the distance** to keep bounds from overshooting. Bounds are aimed from
  projections and checkpoints rather than doubled, but 16% of genvar's band still goes to failed
  rounds (see above); a bound from a quick approximate alignment is untried.
- [ ] **Affine costs with the seed heuristic.** The affine aligners are exact but sweep without a
  heuristic; A*PA2's method for affine costs needs a gap-chaining seed heuristic of its own.
- [x] **Low divergence (below 2%)** was A*PA2's weak spot against BiWFA; the diagonal transition
  before any band now beats BiWFA and WFA there (146 µs against 302 and 605 at 0%, 100 kbp).
- [~] **The seeds' setup** was A*PA2's other limitation; it is a quarter of the time it was, but still
  a quarter of a long read's alignment (see above).
