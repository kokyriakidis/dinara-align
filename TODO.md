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
