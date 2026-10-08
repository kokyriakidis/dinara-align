"""
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
# A batch, over every thread, and one under a cap, None for a pair past it.
var references: List[String] = ["ACGTACGT", "TTGCA"]
var queries: List[String] = ["ACGACGT", "TTGGCA"]
var batch = distances(references, queries)  # [1, 1]
var near = distances(references, queries, costs, max_cost=7)  # [None, None]: a gap of one costs 8
```

| mode | the reference | the query |
| :-- | :-- | :-- |
| `Mode.GLOBAL` | whole | whole |
| `Mode.INFIX` | any part | whole |
| `Mode.PREFIX`, `Mode.SUFFIX` | a prefix, a suffix | whole |
| `Mode.ends_free(...)` | as many letters free at either end as asked | likewise |
| `Mode.extension(match_score, anchor)` | from one end, as far as pays | from the same end |
| `Mode.REFERENCE_IN_QUERY` | whole | any part |
| `Mode.local(match_score)` | any part | any part |
| `Mode.overlap(match_score)` | a prefix or suffix | a suffix or prefix, or whole |

Free ends minimize the costs alone, as Edlib and WFA2-lib count them; `mode.with_match_score(a)`
rewards every match instead, as parasail's and hyalite's semi-global modes do.

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
var best = score("TTTTACGTACGTTTTT", "ACGTACGT", scoring, Mode.LOCAL)  # 16
# Any table, in every mode on the CPU: a read placed in a window, or a seed's extension.
var placed = align("TTTTACGTACGTTTTT", "ACGTACGT", scoring, Mode.INFIX)  # score 16, reference 4..12
```

Ported from AffineGaps by Ash Vardanian, https://github.com/unum-science/AffineGaps, alignment only;
the edit distance from A*PA by Ragnar Groot Koerkamp and Pesho Ivanov (see NOTICE).
"""

from .api import align, alignments, distance, distances, local_scores, score, scores, search
from .common import Device, Placement
from .errors import AlignmentError, ErrorKind
from .gap_affine import DEFAULT_MAX_MEMORY
from .modes import AlignedCounts, Alignment, Anchor, Band, Costs, Mode, Ties
from .scored import LocalScores
from .search import Hit
from .scoring import Scoring
