# To do

Where the time still goes in a single-threaded alignment, profiled on the Skylake-X (3.3 GHz,
pinned) over all 48 reads of ont-500k-genvar, about 10 s in all.

- [ ] **Band rounds that fail (about 1 s of the band's 6.5 s).** The sweep itself runs near the
  hardware's limit for Myers' recurrence, about 2.4 cycles a word-column on AVX-512, so the band
  only gets faster by sweeping fewer cells. On genvar 16% of the band's word-columns go to rounds
  that fail and are retried (7% on ont-500k), three to six rounds a pair.
- [ ] **Seed setup (about 2.5 s).** Inexact seeds: the scan of the second sequence (about 1.4 s
  when forced on every read), local pruning (about 0.8 s) and layer insertion (about 0.4 s).
  The scan is bound by loop overhead and memory, not one hot spot; going further means indexing
  the second sequence instead of the seeds, at more memory.
- [ ] **Traceback (about 0.9 s).** Retracing the final round's tiles from their recorded left
  edges.
- [ ] **The M2's two lost rows.** 100 kbp pairs at 6 and 7% divergence lose to A*PA2-full there:
  the M2's cutoff, 40% of exact seeds chained, rebuilds inexact seeds from 6%, and was tuned when
  inexact setup cost twice what it does now. A sweep of `INEXACT_CHAINED` on the M2 may win them back.

## From A*PA2's discussion

Its limitations and future work (curiouscoding.nl/posts/astarpa2/#discussion), and where dinara-align
stands on each.

- [x] **Symbols beyond `ACGT`.** The edit distance takes `ACGT` and up to four other bytes, `N`
  among them, each matching only itself, through a third bit plane its own kernels read, at no cost
  to bases alone. Still open: such a pair runs without seeds, which pack two bits a base, so a long
  divergent pair with an `N` is slower; seeds that skip any window holding a symbol past `ACGT`, and
  a potential that charges nothing for those seeds, would bring them back.
- [~] **Semi-global and open-ended alignment.** `edit_search` finds a pattern inside a text, Edlib's
  infix mode, or at its start, its prefix mode, and `edit_search_alignment` aligns it there, exact.
  It sweeps the whole matrix, `n * m / 64` word-columns: 10 kbp in 100 kbp in 21 ms on the M2, 100 kbp
  in 1 Mbp in 2 s. Still open: Edlib's band for it, a bound `k` doubled and only the words some score
  within `k` can still reach swept (Ukkonen's cutoff), `k * n / 64` and about fourteen times less on
  the 100 kbp read; and a free end for the pattern too, as WFA's ends-free mode allows.
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
