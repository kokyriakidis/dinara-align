# API reference

Everything `dinara_align` exports, from its docstrings; regenerate with `pixi run api-docs`.

Exact pairwise alignment of DNA, or any text, on the CPU and the GPU.

Every call takes a reference and a query, `Costs` that price each edit and a `Mode` that says which
ends of the two the alignment must reach, and returns the least cost, `distance`, or an optimal
alignment as a CIGAR, `align`. Unit costs, the edit distance, are the default: A*PA2's bit-parallel
band doubling with its seed heuristic finds them; any other costs, and every other mode, a gap-affine
wavefront from both ends, after WFA. Every answer is exact.

```mojo
from dinara_align import Anchor, Band, Costs, Mode, Ties, align, alignments, distance, distances

var edits = distance("ACGTACGTTTGCA", "ACGTCGTTTTGCA")  # 2
var aligned = align("ACGTACGTTTGCA", "ACGTCGTTTTGCA")  # cost 2, cigar "4=1D2=1I6="
var costs = Costs.affine(4, 6, 2)  # a mismatch 4, a gap of k letters 6 + 2k: WFA2-lib's defaults
var affine = align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", costs)  # cost 12, cigar "4=3X6="
# A read placed wherever it fits best in a reference: the reference's letters outside are free.
var placed = align("TTTTACGTACGTTTTT", "ACGTACGT", costs, Mode.INFIX)  # cost 0, "8=", reference 4..12
# Two-piece gap costs, as minimap2's -O and -E take two values each: a gap of k letters the less of
# 6 + 2k and 24 + k.
var long_gap = Costs.two_piece(4, 6, 2, 24, 1)
# Deletions priced apart from insertions, as bwa's -O del,ins: a run of k reference letters 6 + k.
var lopsided = Costs.affine(4, 5, 2).with_deletions(6, 1)
# Of equally good alignments a fixed rule picks one: indels placed left, as minimap2 places them, or
# right, WFA2-lib's CIGARs byte for byte.
var left = align("ACGTTTTACG", "ACGTTTACG", costs)  # "3=1D6="
var right = align("ACGTTTTACG", "ACGTTTACG", costs, ties=Ties.RIGHT)  # "6=1D3="
# Exact within a band of diagonals, KSW2's `w`, or under a cost cap: None when it would pass it.
var banded = align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", costs, band=Band.around(2))
var capped = distance("ACGTACGTTTGCA", "ACGTCGTTTTGCA", costs, max_cost=10)  # None: it costs 12
# A seed's extension, fixed at one end and stopping where it scores best, a match earning 1.
var onward = align("ACGTTGCAAGGCTTTT", "ACGTTGCAAGGCGAGA", costs, Mode.extension(1))  # score 12, "12="
var back = align("TTTTACGTTGCAAGGC", "GAGAACGTTGCAAGGC", costs, Mode.extension(1, Anchor.END))
# The best-scoring part of each, Smith-Waterman, a match earning 2.
var core = align("GGGGACGTACGTGGGG", "CCCCACGTACGTCCCC", costs, Mode.local(2))  # score 16, "8=", 4..12
# A read placed in a window as a mapper scores it, a match earning 2, rather than at the least cost.
var mapped = align("TTTTACGTACGTTTTT", "ACGTCGT", costs, Mode.INFIX.with_match_score(2))  # score 6
# Two reads overlapping, every end gap free: the first's suffix on the second's prefix.
var joined = align("TTTTTACGTACGT", "ACGTACGTGGGGG", costs, Mode.overlap(2))  # score 16, "8=", 5..13
# What a SAM record holds: the CIGAR with the query's unaligned letters soft-clipped, and the NM and
# MD tags.
var sam_cigar = core.clipped_cigar(16)  # "4S8=4S"
var edits = core.edit_distance("GGGGACGTACGTGGGG", "CCCCACGTACGTCCCC")  # NM: 0
var md = core.mismatch_string("GGGGACGTACGTGGGG", "CCCCACGTACGTCCCC")  # MD: "8"
# A batch, many pairs at once on the caller's thread, and one under a cap, None for a pair past it.
var references: List[String] = ["ACGTACGT", "TTGCA"]
var queries: List[String] = ["ACGACGT", "TTGGCA"]
var batch = distances(references, queries)  # [1, 1]
var near = distances(references, queries, costs, max_cost=7)  # [None, None]: a gap of one costs 8
# One thread's aligner, its memory kept from call to call: a loop of single pairs, one a thread.
var aligner = Aligner()
for index in range(len(references)):
    var each = aligner.align(references[index], queries[index], costs)  # as `align` gives it
```

| mode | the reference | the query |
| :-- | :-- | :-- |
| `Mode.GLOBAL` | whole | whole |
| `Mode.INFIX` | any part | whole |
| `Mode.PREFIX`, `Mode.SUFFIX` | a prefix, a suffix | whole |
| `Mode.ends_free(...)` | as many letters free at either end as asked | likewise |
| `Mode.extension(match_score, anchor)` | from one end, as far as pays | from the same end |
| `Mode.local(match_score)`, `Mode.local()` under a `Scoring` | any part | any part |
| `Mode.overlap(match_score)` | a prefix or suffix | a suffix or prefix, or whole |

Free ends minimize the costs alone, as Edlib and WFA2-lib count them; `mode.with_match_score(a)`
rewards every match instead, as parasail's and hyalite's semi-global modes do. The named free ends are
presets of `ends_free`: the reference whole inside the query is `Mode.ends_free(query_start=n,
query_end=n)` for any `n` past the query's length.

A `Scoring`, an alphabet's substitution table and gap scores, which an alignment maximizes, aligns
globally or locally by Gotoh's Needleman-Wunsch or Smith-Waterman, with the initialization corrections
Flouri et al. found missing from the 1982 paper, on the CPU or the GPU, its traceback in linear memory
when the matrix is large; and with free ends or as an extension on the CPU, its span found by sweep
and aligned globally:

```mojo
from dinara_align import Mode, Scoring, align, score

var scoring = Scoring.dna()  # minimap2's: match 2, mismatch -4, a gap of k letters -(4 + 2k)
var found = align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", scoring)  # an Alignment, its cost minus its score
var rows = found.gapped("ACGTACGTTTGCA", "ACGTCGTTTTGCA")  # the two gapped rows
var best = score("TTTTACGTACGTTTTT", "ACGTACGT", scoring, Mode.local())  # 16
# Any table, in every mode on the CPU: a read placed in a window, or a seed's extension.
var placed = align("TTTTACGTACGTTTTT", "ACGTACGT", scoring, Mode.INFIX)  # score 16, reference 4..12
```

The edit distance is ported from A*PA by Ragnar Groot Koerkamp and Pesho Ivanov (see NOTICE).

## Functions

### `align`

```mojo
def align(reference: String, query: String, costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, band: Band = Band(), ties: Ties = Ties.LEFT, eqx: Bool = True, max_memory: Int = Int(83886080)) -> Alignment
```

An optimal alignment of `query` to `reference` as `mode` asks, every move inside `band`, as a CIGAR with `=` and `X`, or with `eqx` false `M` for both (see `Alignment`); raises when no alignment fits the band.

Of several equally good alignments the CIGAR is always the one `ties` names (see `Ties`): by
default every edit as far left as it goes, indels placed as minimap2 places them, or with
`Ties.RIGHT` as far right, WFA2-lib's CIGAR byte for byte, whichever search found the cost.

Every byte is a symbol matching only itself, so DNA in either case, or any other text, needs no
alphabet. The wavefront's work grows with the square of the cost rather than with the matrix, and
its memory stays bounded: the fronts it keeps for the traceback, the bulk of it, stay within
`max_memory` bytes, past which the pair is split where an optimal path crosses and each piece
aligned alone (see `gap_affine.solve`), the cost still the least, the tie rule followed within each
piece. The sweeps and the bit-parallel search keep a few rows, or a band's edges, whatever the cap.

```mojo
def align(reference: String, query: String, costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, max_cost: Int, band: Band = Band(), ties: Ties = Ties.LEFT, eqx: Bool = True, max_memory: Int = Int(83886080)) -> Optional[Alignment]
```

`align`, or None when the cost would pass `max_cost` or no alignment fits `band`, found as `distance` finds that, with no fronts traced. A mode with a match score, which maximizes a score, takes no cap.

```mojo
def align(reference: String, query: String, scoring: Scoring, mode: Mode = Mode.GLOBAL, *, placement: Optional[Placement] = None, max_memory: Int = Int(83886080), eqx: Bool = True) -> Alignment
```

An optimal alignment under `scoring`, as `Costs` give one (see `Alignment`), its `cost` minus its score: both sequences whole for `Mode.GLOBAL`, Needleman-Wunsch, the best-scoring window of each for `Mode.local()`, Smith-Waterman, on either device; free ends and extensions, with Z-drop as KSW2 gauges it, on the host, their span by sweep and the letters between aligned globally (see `scoring.scoring_alignment`). Its rows come back with `Alignment.gapped`. Of equally good alignments, Gotoh's walk picks the CIGAR (see `alignment.reconstruct`), not `Ties`. A traceback whose stored matrix would pass `max_memory` bytes recurses in linear space instead (see `scoring.cells_within`).

### `alignments`

```mojo
def alignments(references: List[String], queries: List[String], scoring: Scoring, mode: Mode = Mode.GLOBAL, *, placement: Optional[Placement] = None, max_memory: Int = Int(83886080), eqx: Bool = True) -> List[Alignment]
```

`align` for every pair; on the device, every pair both bounds admit goes out in one launch.

```mojo
def alignments(references: List[String], queries: List[String], costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, band: Band = Band(), ties: Ties = Ties.LEFT, eqx: Bool = True, threads: Optional[Int] = None, max_memory: Int = Int(83886080)) -> List[Alignment]
```

Every pair's `align`, the pairs spread over threads as `distances` spreads them, each thread's kept fronts within `max_memory` bytes.

```mojo
def alignments(references: List[String], queries: List[String], costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, max_cost: Int, band: Band = Band(), ties: Ties = Ties.LEFT, eqx: Bool = True, threads: Optional[Int] = None, max_memory: Int = Int(83886080)) -> List[Optional[Alignment]]
```

Every pair's `align` under `max_cost`, None for a pair past it or with no alignment inside `band`, the pairs spread over threads as `distances` spreads them.

### `distance`

```mojo
def distance(reference: String, query: String, costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, band: Band = Band()) -> Int
```

The least cost of aligning `query` to `reference` as `mode` asks, every move inside `band` (see `Band`), with no alignment; raises when no alignment fits the band.

Unit costs, globally, take A*PA2's band doubling: guess a bound, sweep only the cells a path within
it could cross, raise the guess until the answer fits under it. Other costs take the wavefront from
both ends, keeping only a few costs' fronts, so a few rows of memory however long the pair.

```mojo
def distance(reference: String, query: String, costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, max_cost: Int, band: Band = Band()) -> Optional[Int]
```

`distance`, or None when it would pass `max_cost` or no alignment fits `band`: the searches stop once each has grown to about half of `max_cost` without the two meeting within it, so a pair far over costs a fraction of its full search.

### `distances`

```mojo
def distances(references: List[String], queries: List[String], costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, band: Band = Band(), threads: Optional[Int] = None, placement: Optional[Placement] = None) -> List[Int]
```

Every pair's `distance`, on the caller's own thread by default, or spread over `threads` threads when asked: an application that calls the library from threads of its own spreads its work itself, and it alone knows how many its machine can spare.

The pairs are independent, so each runs on one thread start to finish, each thread taking the
next pair of the batch, longest first, as soon as it is free (see `longest_first`). A pair that
fails raises, after the rest, the same error a serial loop would have raised first.

On the GPU, `placement`, unit costs align globally, a thread a pair by Myers' bit-vectors (see
`device_edit`), every pair whose shorter sequence fits a thread's 4,096 letters; the rest, and
every other cost or mode, on the host.

```mojo
def distances(references: List[String], queries: List[String], costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, max_cost: Int, band: Band = Band(), threads: Optional[Int] = None) -> List[Optional[Int]]
```

Every pair's `distance` under `max_cost`, None for a pair past it or with no alignment inside `band`, the pairs spread over threads as the uncapped `distances` spreads them: a batch of candidates filtered by cost, the far ones costing a fraction of their full search.

### `local_scores`

```mojo
def local_scores(reference: String, query: String, costs: Costs, mode: Mode, *, window: Optional[Int] = None) -> LocalScores
```

A local alignment's best score and where it ends, with the best score of an alignment ending more than `window` reference letters away, as SSW's `score2` and `ref_end2` report it for a mapping quality (see `LocalScores`): one sweep, with no alignment traced. The window is half the query, and at least 15, by default, as SSW suggests.

### `score`

```mojo
def score(reference: String, query: String, costs: Costs, mode: Mode = Mode.GLOBAL, *, band: Band = Band()) -> Int
```

The best score `align` would return, with no alignment traced: for a mode with a match score its matches' reward less its costs, else minus the least cost, `distance`'s. A local alignment or free ends with a reward take the sweep alone, an extension its search alone, and a global alignment with a reward the wavefront's cost with the reward folded in, so each skips the traceback.

```mojo
def score(reference: String, query: String, scoring: Scoring, mode: Mode = Mode.GLOBAL, *, placement: Optional[Placement] = None) -> Int
```

The optimal score under `scoring`, with no alignment traced: `Mode.GLOBAL` and `Mode.local()` in two rows of memory on either device (see `scoring.score_with`), free ends and extensions by sweep on the host. The table holds what a match earns, so a mode's own match score must be zero: `Mode.extension(0)` for an extension.

### `scores`

```mojo
def scores(references: List[String], queries: List[String], scoring: Scoring, mode: Mode = Mode.GLOBAL, *, placement: Optional[Placement] = None) -> List[Int]
```

`score` for every pair; on the device, every pair one block can carry goes out in one launch.

### `search`

```mojo
def search(references: List[String], query: String, costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, best: Optional[Int] = None, max_cost: Optional[Int] = None, aligned: Bool = False, ties: Ties = Ties.LEFT, threads: Optional[Int] = None) -> List[Hit]
```

The query against every reference, a database search: each reference's `Hit`, its score, the best first, ties by the references' order; with `best` that many alone, with `max_cost` (a mode with no reward) those within it alone, and with `aligned` each kept hit's alignment too.

A local alignment scores a group of references at once, one to a SIMD lane, as SWIPE does (see
`local_scores_by_lane`); a mode with no reward takes `distances`, under the cap when there is one;
any other mode each pair's `score`. Every kept hit is then aligned on its own, when asked for, by
`align`.

## Types

### `Aligner`

```mojo
struct Aligner
```

One thread's aligner: `distance` and `align` as the functions of those names give them, the memory their searches take kept from call to call, so a loop of calls takes none once it is warm. An application calling from many threads keeps one a thread; one never crosses threads, and the library keeps no state of its own, so any number of them work at once.

| field | type | |
| :-- | :-- | :-- |
| `space` | `SearchSpace` |  |

#### `__init__`

```mojo
def Aligner.__init__(out self)
```

An aligner holding no memory yet: its first call takes what it needs.

#### `distance`

```mojo
def distance(mut self, reference: String, query: String, costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, band: Band = Band()) -> Int
```

`distance`, through this aligner's memory.

```mojo
def distance(mut self, reference: String, query: String, costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, max_cost: Int, band: Band = Band()) -> Optional[Int]
```

`distance` under a cap, through this aligner's memory.

#### `align`

```mojo
def align(mut self, reference: String, query: String, costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, band: Band = Band(), ties: Ties = Ties.LEFT, eqx: Bool = True, max_memory: Int = Int(83886080)) -> Alignment
```

`align`, through this aligner's memory: the same alignment, the CIGAR its `ties` picks.

```mojo
def align(mut self, reference: String, query: String, costs: Costs = Costs.edit(), mode: Mode = Mode.GLOBAL, *, max_cost: Int, band: Band = Band(), ties: Ties = Ties.LEFT, eqx: Bool = True, max_memory: Int = Int(83886080)) -> Optional[Alignment]
```

`align` under a cap, through this aligner's memory.

### `Device`

```mojo
struct Device
```

The host's cores or a GPU.

| field | type | |
| :-- | :-- | :-- |
| `kind` | `UInt8` |  |

- `Device.CPU` = `Device(UInt8(0))`: 
- `Device.GPU` = `Device(UInt8(1))`: 

### `Placement`

```mojo
struct Placement
```

Where a call runs, which GPU, and how many of the host's threads it may take. Built by `on_cpu` or `on_gpu`, so a call on the host never carries a GPU's number.

| field | type | |
| :-- | :-- | :-- |
| `device` | `Device` |  |
| `gpu_id` | `Int` | The GPU's number, zero on the host. |
| `threads` | `Int` | Host threads the call may take, at least one and no more than the process can run. |

#### `__init__`

```mojo
def Placement.__init__(device: Device, gpu_id: Int, threads: Int) -> Self
```

`device`, the GPU numbered `gpu_id` on a GPU, and `threads` held to what the process runs.

#### `on_cpu`

```mojo
def Placement.on_cpu(threads: Int) -> Self
```

The host, over up to `threads` threads.

#### `on_gpu`

```mojo
def Placement.on_gpu(gpu_id: Int, threads: Int) -> Self
```

The GPU numbered `gpu_id`, the host's part of the work, its packing, over up to `threads` threads.

#### `default`

```mojo
def Placement.default() -> Self
```

The host, on the caller's own thread: an application spreads its calls over its threads itself, and asks for more here only when it wants this call spread too.

### `AlignmentError`

```mojo
struct AlignmentError
```

A failed call: which way it failed, and what it failed on, the value or the input to look at.

| field | type | |
| :-- | :-- | :-- |
| `kind` | `ErrorKind` |  |
| `detail` | `String` |  |

#### `write_to`

```mojo
def write_to(self, mut writer: T)
```

`dinara-align: <phrase> [<detail>]`, the form the C API, the Python package and the command line pass on.

### `ErrorKind`

```mojo
struct ErrorKind
```

The way a call failed, which callers branch on: the C API turns each into its own code.

| field | type | |
| :-- | :-- | :-- |
| `id` | `UInt8` |  |

- `ErrorKind.UNKNOWN_SYMBOL` = `ErrorKind(UInt8(1))`: A sequence holds a letter its `Scoring`'s alphabet lacks, or a byte no sequence may hold.
- `ErrorKind.ALPHABET_TOO_LARGE` = `ErrorKind(UInt8(2))`: An alphabet with more letters than a substitution table holds.
- `ErrorKind.SEQUENCE_TOO_LONG` = `ErrorKind(UInt8(3))`: The input needs more memory than it may take: past `max_memory`, or past what the device can allocate at once.
- `ErrorKind.SCRATCH_TOO_SMALL` = `ErrorKind(UInt8(4))`: Device memory sized for a smaller problem than the one launched on it.
- `ErrorKind.LENGTH_MISMATCH` = `ErrorKind(UInt8(5))`: Two lists, or two gapped rows, of unequal length where each item pairs with one of the other.
- `ErrorKind.INVALID_SCORING` = `ErrorKind(UInt8(6))`: Costs or scores no search can use: an edit that costs nothing, a reward that costs, a value past the range its arithmetic holds.
- `ErrorKind.INVALID_ARGUMENT` = `ErrorKind(UInt8(7))`: An argument outside what the call takes: a mode, a band, a cap or a count it refuses.
- `ErrorKind.OUTSIDE_BAND` = `ErrorKind(UInt8(8))`: Every alignment leaves the band of diagonals asked for.

#### `phrase`

```mojo
def phrase(self) -> StaticString
```

The failure in a few words, the start of every message of this kind.

#### `write_to`

```mojo
def write_to(self, mut writer: T)
```

The kind's phrase.

### `AlignedCounts`

```mojo
struct AlignedCounts
```

An alignment's columns by kind (see `cigar_counts`).

| field | type | |
| :-- | :-- | :-- |
| `matches` | `Int` |  |
| `mismatches` | `Int` |  |
| `deleted` | `Int` | Reference letters against a gap, `D`. |
| `inserted` | `Int` | Query letters against a gap, `I`. |

### `Alignment`

```mojo
struct Alignment
```

An optimal alignment: of `reference[reference_start:reference_end]` against `query[query_start:query_end]`, as a CIGAR over those letters alone, `=` a match and `X` a substitution (or `M` for either), `D` a reference letter alone and `I` a query letter alone, each run its length then its letter. A global alignment spans both; letters a mode leaves unaligned for nothing lie outside the spans, and a gap past them stays in the CIGAR, as it is paid.

`cost` is what the CIGAR's edits cost. `score` is an extension's: the matches' reward less the cost;
outside an extension, with no reward, it is minus the cost.

| field | type | |
| :-- | :-- | :-- |
| `cost` | `Int` |  |
| `score` | `Int` |  |
| `cigar` | `String` |  |
| `reference_start` | `Int` |  |
| `reference_end` | `Int` |  |
| `query_start` | `Int` |  |
| `query_end` | `Int` |  |

#### `gapped`

```mojo
def gapped(self, reference: String, query: String) -> Tuple[String, String]
```

The aligned parts of both sequences as two rows of one length, `-` against each gapped letter.

#### `clipped_cigar`

```mojo
def clipped_cigar(self, query_length: Int, *, hard: Bool = False) -> String
```

The CIGAR as a SAM record writes it: the query's letters outside the span clipped, soft (`S`), their letters kept in the record's sequence, or with `hard` hard (`H`), dropped from it. The record's position is `reference_start + 1`; the reference's letters outside the span need no operation.

#### `mirrored`

```mojo
def mirrored(self, reference_length: Int, query_length: Int) -> Self
```

This alignment of both sequences reversed, `reference_length` and `query_length` letters, turned back: its CIGAR's runs in the other order and its spans counted from the other ends.

#### `counts`

```mojo
def counts(self, reference: String, query: String) -> AlignedCounts
```

How many letters the alignment pairs equal and unequal, and leaves gapped either way, `M` runs compared letter by letter.

#### `edit_distance`

```mojo
def edit_distance(self, reference: String, query: String) -> Int
```

The alignment's edits, SAM's `NM` tag: its substitutions and gapped letters.

#### `identity`

```mojo
def identity(self, reference: String, query: String) -> Float64
```

The matches over the alignment's columns, BLAST's identity: one for an exact match, zero for an empty alignment.

#### `mismatch_string`

```mojo
def mismatch_string(self, reference: String, query: String) -> String
```

SAM's `MD` tag: the reference's letters the alignment does not match, each substitution's letter after the matches before it and each deletion's after a `^`, so the reference can be rebuilt from the query and the CIGAR; insertions leave no mark.

### `Anchor`

```mojo
struct Anchor
```

Which end of both sequences an extension is fixed at (see `Mode.extension`).

| field | type | |
| :-- | :-- | :-- |
| `identifier` | `UInt8` |  |

- `Anchor.START` = `Anchor(UInt8(0))`: Both sequences' first letters: the alignment runs right from there, a seed's right extension.
- `Anchor.END` = `Anchor(UInt8(1))`: Both sequences' last letters: the alignment runs left from there, a seed's left extension.

### `Band`

```mojo
struct Band
```

The diagonals an alignment may use, `low ..= high`: a cell's diagonal is the reference's letters aligned or skipped up to it less the query's, counted from the alignment's fixed origin, so every move from the origin on stays inside. A global alignment starts on diagonal zero and ends on the reference's length less the query's, which the band must hold.

The alignment found is the optimum over every path inside the band, exact, not a heuristic: only
the diagonals outside are never searched. KSW2's band of width `w` is `Band.around(w)`; WFA2-lib's
static band counts diagonals the other way, its `min_k ..= max_k` being `Band(-max_k, -min_k)`.

| field | type | |
| :-- | :-- | :-- |
| `low` | `Int` | The lowest diagonal the alignment may use. |
| `high` | `Int` | The highest diagonal the alignment may use. |

#### `__init__`

```mojo
def Band.__init__() -> Self
```

No band: every diagonal.

```mojo
def Band.__init__(low: Int, high: Int) -> Self
```

The diagonals `low ..= high`, each bound held within `UNBOUNDED` of the origin's, past any pair's diagonals, so an integer type's limits stand for no bound and nothing overflows.

#### `around`

```mojo
def Band.around(width: Int) -> Self
```

The diagonals at most `width` from the origin's, either way.

#### `holds`

```mojo
def holds(self, diagonal: Int) -> Bool
```

Whether `diagonal` lies inside the band.

#### `covers`

```mojo
def covers(self, columns: Int, rows: Int) -> Bool
```

Whether every diagonal of a `columns` by `rows` matrix lies inside: no band at all for it.

#### `covers_any`

```mojo
def covers_any(self) -> Bool
```

Whether every diagonal of any pair a batch could hold lies inside, sequences of up to `1 << 40` letters: no band for any pair.

#### `shifted`

```mojo
def shifted(self, origin: Int) -> Self
```

The band as seen from a cell on diagonal `origin`, the start of a piece after a split.

#### `mirrored`

```mojo
def mirrored(self, target: Int) -> Self
```

The band in the reversed sequences of a pair whose end lies on diagonal `target`, where diagonal `k` reads `target - k`.

### `Costs`

```mojo
struct Costs
```

What each edit costs, as WFA and minimap2 count it: a substitution `mismatch`, and a gap of `k` letters `opening + k extension`, or with a second piece the less of that and `opening2 + k extension2`. A deletion, a run of reference letters alone, may cost otherwise than an insertion, a run of query letters alone, as bwa's `-O del,ins` (see `with_deletions`). Every alignment minimizes the total, save an extension (see `Mode.extension`).

Built through `edit`, `linear`, `affine` or `two_piece`, which refuse costs no search can run by.

| field | type | |
| :-- | :-- | :-- |
| `mismatch` | `Int` |  |
| `opening` | `Int` | An insertion's opening, and a deletion's unless `with_deletions` set its own. |
| `extension` | `Int` |  |
| `opening2` | `Int` | The second gap piece's opening, -1 when there is none. |
| `extension2` | `Int` |  |
| `deletion_opening` | `Int` | A deletion's own opening, its extension and its second piece's after it: `opening` and the rest unless `with_deletions` set them. |
| `deletion_extension` | `Int` |  |
| `deletion_opening2` | `Int` |  |
| `deletion_extension2` | `Int` |  |

#### `__init__`

```mojo
def Costs.__init__(mismatch: Int, opening: Int, extension: Int, opening2: Int, extension2: Int) -> Self
```

Costs whose deletions cost what insertions do, trusted as given: the factories check them.

#### `edit`

```mojo
def Costs.edit() -> Self
```

Unit costs: the edit (Levenshtein) distance, a substitution, an insertion and a deletion one each. A global alignment, or a query found inside or at the start of the reference, takes the bit-parallel band doubling of A*PA2; the other modes the wavefront, at the same costs.

#### `linear`

```mojo
def Costs.linear(mismatch: Int, gap: Int) -> Self
```

A substitution `mismatch` and every gapped letter `gap`, with no opening.

#### `affine`

```mojo
def Costs.affine(mismatch: Int, opening: Int, extension: Int) -> Self
```

Gap-affine costs: a substitution `mismatch`, a gap of `k` letters `opening + k extension`. minimap2's `-B4 -O4 -E2` is `affine(4, 4, 2)`, WFA2-lib's default `affine(4, 6, 2)`.

#### `two_piece`

```mojo
def Costs.two_piece(mismatch: Int, opening: Int, extension: Int, opening2: Int, extension2: Int) -> Self
```

Two-piece gap-affine costs, minimap2's `-O4,24 -E2,1` as `two_piece(4, 4, 2, 24, 1)`: a gap of `k` letters the less of `opening + k extension` and `opening2 + k extension2`, one piece usually cheap to open and the other cheap to extend, so a long gap costs less than one piece charges.

#### `with_deletions`

```mojo
def with_deletions(self, opening: Int, extension: Int, opening2: Int = Int(-1), extension2: Int = Int(0)) -> Self
```

These costs with a deletion, a run of `k` reference letters alone, costing `opening + k extension`, or the less of that and `opening2 + k extension2`, and an insertion as before: bwa's `-O6,5 -E1,2` is `affine(4, 5, 2).with_deletions(6, 1)`, its insertions first. Either side may have a second piece the other lacks: the one without counts its one piece twice.

#### `pieces`

```mojo
def pieces(self) -> Int
```

How many gap pieces the costs have, either side.

#### `symmetric`

```mojo
def symmetric(self) -> Bool
```

Whether a deletion costs what an insertion does.

#### `unit_scale`

```mojo
def unit_scale(self) -> Int
```

The factor these costs are of unit costs, zero when they are not: a pair's edit distance times it is then their least cost, which the bit-parallel search finds.

#### `cheapest_extension`

```mojo
def cheapest_extension(self) -> Int
```

The cheapest a gap grows a letter, either way, at either piece: a Z-drop's slack a diagonal, as KSW2 charges a long gap.

#### `dearest_step`

```mojo
def dearest_step(self) -> Int
```

The dearest single move: a substitution, or a gap's first letter at either piece either way.

#### `gap`

```mojo
def gap(self, letters: Int, deleted: Bool) -> Int
```

What a gap of `letters` letters costs, a deletion or an insertion, at its cheaper piece.

### `Mode`

```mojo
struct Mode
```

Which alignments of the two sequences count: how many letters at each end of each may be left unaligned for nothing, or an extension from one end, or for a `Scoring` a local alignment.

| mode | reference | query | also called |
| :-- | :-- | :-- | :-- |
| `GLOBAL` | whole | whole | end to end, Needleman-Wunsch, Edlib's NW |
| `INFIX` | any part | whole | semi-global, glocal, Edlib's HW |
| `PREFIX` | a prefix | whole | Edlib's SHW |
| `SUFFIX` | a suffix | whole | |
| `ends_free(...)` | as asked | as asked | WFA2-lib's ends-free; the query's ends free place the reference inside it |
| `extension(...)` | from one end | from the same end | KSW2's extension, with or without Z-drop |
| `local(...)` | any part | any part | Smith-Waterman, abPOA's local mode |
| `overlap(...)` | a prefix or suffix | a suffix or prefix | semi-global, parasail's `sg`, hyalite's OV |

Free ends minimize the costs alone, as Edlib and WFA2-lib count them, unless a match earns
something (see `with_match_score`), as parasail's and hyalite's do. Three kinds lie beneath, as WFA2-lib's
ends-free and extension and abPOA's local are: the named free ends are presets of `ends_free`.

| field | type | |
| :-- | :-- | :-- |
| `kind` | `UInt8` |  |
| `reference_start` | `Int` | Letters at the reference's start that may go unaligned for nothing; past them a gap is paid. |
| `reference_end` | `Int` |  |
| `query_start` | `Int` |  |
| `query_end` | `Int` |  |
| `match_score` | `Int` | What a match earns: in an extension or a local alignment, which maximize a score, and with free ends when asked (see `with_match_score`); zero elsewhere, the costs alone minimized. |
| `anchor` | `Anchor` |  |
| `zdrop` | `Int` | An extension's Z-drop, -1 for none (see `extension`). |
| `end_bonus` | `Int` | What an extension reaching the query's far end earns over its score, -1 for none (see `extension`). |

- `Mode.ENDS` = `UInt8(0)`: The kind of a global alignment and of every one with free ends.
- `Mode.EXTENSION` = `UInt8(1)`: The kind of an extension from one end (see `extension`).
- `Mode.SMITH_WATERMAN` = `UInt8(2)`: The kind of a local alignment (see `local`).
- `Mode.GLOBAL` = `Mode(Mode.ENDS, Int(0), Int(0), Int(0), Int(0), Int(0), Anchor.START, Int(-1), Int(-1))`: Both sequences end to end.
- `Mode.INFIX` = `Mode(Mode.ENDS, Int(1152921504606846976), Int(1152921504606846976), Int(0), Int(0), Int(0), Anchor.START, Int(-1), Int(-1))`: The whole query against wherever in the reference it fits best: a read placed in a window.
- `Mode.PREFIX` = `Mode(Mode.ENDS, Int(0), Int(1152921504606846976), Int(0), Int(0), Int(0), Anchor.START, Int(-1), Int(-1))`: The whole query against the reference's best prefix.
- `Mode.SUFFIX` = `Mode(Mode.ENDS, Int(1152921504606846976), Int(0), Int(0), Int(0), Int(0), Anchor.START, Int(-1), Int(-1))`: The whole query against the reference's best suffix.

#### `ends_free`

```mojo
def Mode.ends_free(*, reference_start: Int = Int(0), reference_end: Int = Int(0), query_start: Int = Int(0), query_end: Int = Int(0)) -> Self
```

Up to so many letters at each end of each sequence left unaligned for nothing, as WFA2-lib's ends-free alignment counts them; all zero is `GLOBAL`. An overlap of two reads frees one's start and the other's end, and the query's both ends place the whole reference inside it. With costs alone, freeing both ends of both lets the empty alignment win, at no cost; a match score makes the alignment the best-scoring one instead (see `with_match_score`).

#### `with_match_score`

```mojo
def with_match_score(self, match_score: Int) -> Self
```

These free ends with every match earning `match_score`: the best-scoring alignment, the reward less the costs, as parasail's and hyalite's semi-global modes count it, rather than the least costly one. With some letters left free the two can differ, as a reward pays for aligning letters a cost alone would leave out. `Mode.INFIX.with_match_score(2)` is a read placed in a window as a mapper scores it. Its time grows with the matrix, as `local`'s does, but for a global alignment, whose letters are all aligned, so the reward folds into the costs and the wavefront finds it.

#### `extension`

```mojo
def Mode.extension(match_score: Int, anchor: Anchor = Anchor.START, *, zdrop: Optional[Int] = None, end_bonus: Optional[Int] = None) -> Self
```

The best-scoring alignment fixed at one end of both sequences, `anchor`, and free to stop anywhere: a read mapper's seed extension. A match earns `match_score` and every edit costs what `Costs` charges; aligning nothing scores zero. A reward is what makes stopping a choice: with costs alone, aligning nothing would always win.

Exact by default: the best stop of all. With `zdrop`, minimap2's `-z` and KSW2's Z-drop, the
search gives up once every alignment it is growing scores more than `zdrop`, plus a gap
extension a diagonal between them, below the best so far, as WFA2-lib's Z-drop gauges it a cost
at a time, and the best stop it found stands: a heuristic, faster on a seed whose read turns to
noise, which may miss a better stop past a divergent stretch.

With `end_bonus`, KSW2's and minimap2's `--end-bonus`, BWA-MEM's clipping penalty, an extension
that reaches the query's far end is preferred whenever its score plus the bonus passes the best
stop's: the read is aligned to its end unless stopping short gains more than the bonus. The one
reaching the end is the best of those that do, with the reference's far end free, and its score
is its own, the bonus only choosing it. As in KSW2, a search that the Z-drop gave up never
reaches the end. A bonus of zero changes nothing: no alignment reaching the end scores more
than the best stop.

#### `reaching_end`

```mojo
def reaching_end(self) -> Self
```

The free ends an extension that reaches the query's far end takes, for its end bonus: from the same anchor, the whole query against the reference, whose far end is free.

#### `local`

```mojo
def Mode.local(match_score: Int = Int(0)) -> Self
```

The best-scoring alignment of any part of the reference against any part of the query, Smith-Waterman: under `Costs` a match earns `match_score`, above zero, and every edit costs what they charge; under a `Scoring`, `Mode.local()`, its table says what each pair earns. It is every end free, with a reward: with costs alone, aligning nothing would always win. Its time grows with the matrix, as every local aligner's does (see `scored`).

#### `overlap`

```mojo
def Mode.overlap(match_score: Int) -> Self
```

The best-scoring alignment with every end gap free, semi-global: it starts on either sequence's first letter and ends on either's last, so one may overhang the other at each end, an overlap of two reads, or one may lie inside the other. A match earns `match_score`, as with all four ends free and costs alone the empty alignment would win. Its time grows with the matrix, as `local`'s does.

#### `is_global`

```mojo
def is_global(self) -> Bool
```

Whether these are free ends with none free: both sequences end to end, whatever a match earns.

#### `is_scored`

```mojo
def is_scored(self) -> Bool
```

Whether the alignment maximizes a score: an extension, a local alignment, or free ends with a match that earns.

### `Ties`

```mojo
struct Ties
```

Which of several equally good alignments a CIGAR spells. Both rules are WFA2-lib's backtrace: at each step back the edit that reached furthest, ties going to a substitution, then a letter of the reference alone, then one of the query, the second gap piece before the first and a gap's extension before its opening; they differ in the end it runs from. So the CIGAR is the same however the search found the cost, whatever the band or cap, and whatever the costs' common factor.

With free ends the span comes first. `LEFT` decides it from the end back, as it places edits: the
end on the highest diagonal an equally good alignment reaches, the most reference letters less query
letters, the furthest along the reference for a read placed in it, then the start on the highest
diagonal of those ending there; `RIGHT` is that over both sequences reversed, the start on the
lowest diagonal, then the end. The letters between are then aligned globally by the rule. A local
alignment likewise ends as late as an equally good one allows and starts as late too under `LEFT`,
and under `RIGHT` starts and ends as early (see `scored`).

| field | type | |
| :-- | :-- | :-- |
| `identifier` | `UInt8` |  |

- `Ties.LEFT` = `Ties(UInt8(0))`: Every edit as early as an equally good alignment allows, gaps shifted left through repeats: the rule run from the start over both sequences reversed, as KSW2 places gaps by default and as variant callers normalize indels.
- `Ties.RIGHT` = `Ties(UInt8(1))`: Every edit as late as it allows, gaps shifted right: WFA2-lib's own CIGARs, byte for byte.

### `LocalScores`

```mojo
struct LocalScores
```

A local alignment's best score and where it ends, and the best score of an alignment ending elsewhere, as SSW reports them for a mapping quality: `second_score` the best of any cell whose reference letters lie more than the window from `reference_end`, at `second_reference_end`, the first such column; zero, at zero, when there is none.

| field | type | |
| :-- | :-- | :-- |
| `score` | `Int` |  |
| `reference_end` | `Int` |  |
| `query_end` | `Int` |  |
| `second_score` | `Int` |  |
| `second_reference_end` | `Int` |  |

### `Hit`

```mojo
struct Hit
```

One reference's result in a search: its place in the list searched, its best score (minus its least cost for a mode with no reward), and its alignment when asked for.

| field | type | |
| :-- | :-- | :-- |
| `index` | `Int` |  |
| `score` | `Int` |  |
| `alignment` | `Optional[Alignment]` |  |

### `Scoring`

```mojo
struct Scoring
```

An alphabet, the substitution table it indexes, and the affine gap model, which travel together.

Built through `dna`, `edit_distance`, `uniform` or `tabulated` rather than field by field, so a
table whose shape disagrees with its alphabet cannot be expressed. A gap of `k` letters scores
`opening + k extension`, both scores zero or less, as `Costs` counts a gap's cost.

| field | type | |
| :-- | :-- | :-- |
| `alphabet` | `String` | The letters a sequence may hold, in the order the table is indexed by. |
| `substitutions` | `List[Int8]` | Row-major, one row per letter of `alphabet`. |
| `gaps` | `AffineGapCosts` | The gap scores as the kernels take them: a gap's first letter scores `open`, each further one `extend`. |

#### `__init__`

```mojo
def Scoring.__init__(out self, var alphabet: String, var substitutions: List[Int8], gaps: AffineGapCosts)
```

Trusts its arguments; the factories are the checked way in.

#### `dna`

```mojo
def Scoring.dna() -> Self
```

Minimap2's scoring over `ACGT`: match 2, mismatch -4, a gap of `k` letters `-(4 + 2k)`.

#### `edit_distance`

```mojo
def Scoring.edit_distance(alphabet: String = DNA_ALPHABET) -> Self
```

Unit costs, under which a global score is the negated Levenshtein distance.

#### `uniform`

```mojo
def Scoring.uniform(match_score: Int, mismatch_score: Int, opening: Int = Int(-4), extension: Int = Int(-2), alphabet: String = DNA_ALPHABET) -> Self
```

One score for equal letters and one for unequal, over any alphabet.

#### `tabulated`

```mojo
def Scoring.tabulated(alphabet: String, var substitutions: List[Int8], opening: Int = Int(-4), extension: Int = Int(-2)) -> Self
```

A caller's own table, refused unless it is square in the alphabet that indexes it.

#### `penalties`

```mojo
def penalties(self) -> Optional[Penalties]
```

The wavefront's costs for this table, if it holds one match and one mismatch score whose folded costs a wavefront can grow by (see `gap_affine.wavefront_penalties`).

#### `alphabet_size`

```mojo
def alphabet_size(self) -> Int
```

Letters the table is indexed by, which is its stride.

## Constants

- `DEFAULT_MAX_MEMORY` = `83886080`: The bytes of kept fronts an alignment may hold by default, about 80 MB (see `HISTORY_LIMIT`).
