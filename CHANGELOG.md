# Changelog

Every release's changes, newest first. A tag `vX.Y.Z` publishes the section of that version as the
release's notes (see `.github/workflows/release.yml`).

## [Unreleased]

The first release: exact pairwise alignment in every mode, from Mojo, C, C++, Python and the command
line.

### Alignment

- One API for every cost model and mode: `distance`, `align` and `score` take `Costs` (unit, linear,
  gap-affine, two-piece gap-affine, deletions priced apart from insertions) and a `Mode` (global,
  infix, prefix, suffix, the reference inside the query, any free ends, extension from either end
  with an optional Z-drop, local, overlap, and free ends rewarding every match).
- Unit costs by A\*PA2's bit-parallel band doubling with its seed heuristic; every other cost model by
  a gap-affine wavefront from both ends after BiWFA; local alignment, overlaps and rewarded free
  ends by an anti-diagonal sweep in 16-bit lanes; all exact, every alignment's memory bounded
  (`max_memory`).
- A `Scoring`, any alphabet's substitution table with affine gap scores, in every mode on the CPU,
  and globally or locally on the GPU.
- One rule for ties (`Ties`): of equally good alignments, indels placed left by default, as minimap2
  places them, or right, WFA2-lib's CIGAR byte for byte; with free ends the span first, by the same
  rule, whichever search found the cost, under any band or cap.
- Bands of diagonals, cost caps (`max_cost`), and batches over every thread (`distances`,
  `alignments`), capped too.
- What a SAM record holds: soft- or hard-clipped CIGARs, `NM`, `MD`, identity
  (`Alignment.clipped_cigar`, `edit_distance`, `mismatch_string`, `identity`), and SSW's second-best
  local score for a mapping quality (`local_scores`).

### Interfaces

- A C API with a C++ wrapper (`c/dinara.h`, `pixi run build-c`), batches included.
- A Python package (`pixi run build-python`, `pixi run build-wheel`), any CPython 3.
- A command-line aligner, `dinara-align` (`pixi run build-cli`): FASTA, FASTQ or pairs in, a table,
  SAM or PAF out.

### Quality

- Every answer held to a full-matrix oracle by a differential fuzzer (`pixi run fuzz`), over pairs,
  costs, tables, modes, bands, caps, tie rules and memory budgets.
- WFA2-lib's regression set (`pixi run test-wfa`), byte for byte under `Ties.RIGHT`.
- Continuous integration on x86-64 and AArch64 Linux and Apple silicon.
- Benchmarks build every tool for the same CPU, named in each table.
