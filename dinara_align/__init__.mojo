"""
Exact affine-gap DNA alignment with traceback, on the CPU and the GPU.

Gotoh's Needleman-Wunsch and Smith-Waterman, with the initialization corrections Flouri et al.
found missing from the 1982 paper, and Levenshtein as their unit-cost limit. The alignment itself
is reconstructed, not only the score, in linear memory through a Hirschberg recursion with a
Myers-Miller affine join, and on the GPU the reconstruction runs on the device too.

```mojo
from dinara_align import Scoring, needleman_wunsch_gotoh_alignment

var aligned = needleman_wunsch_gotoh_alignment("ACGTACGTTTGCA", "ACGTCGTTTTGCA", Scoring.dna())
print(aligned.first_gapped)
print(aligned.second_gapped)
print(aligned.score)
```

Unit-cost edit distance has its own bit-parallel path, A*PA2's band doubling with its seed heuristic,
for one pair on one thread, or a batch spread over every thread a pair at a time:

```mojo
from dinara_align import edit_alignment, edit_alignments, edit_cigar, edit_distance, edit_distances

var distance = edit_distance("ACGTACGTTTGCA", "ACGTCGTTTTGCA")
var aligned = edit_alignment("ACGTACGTTTGCA", "ACGTCGTTTTGCA")
var cigar = edit_cigar("ACGTACGTTTGCA", "ACGTCGTTTTGCA").cigar  # "4=1D2=1I6=", the first the reference
var firsts: List[String] = ["ACGTACGT", "TTGCA"]
var seconds: List[String] = ["ACGACGT", "TTGGCA"]
var distances = edit_distances(firsts, seconds)
```

Gap-affine costs as WFA counts them, a mismatch and a gap's opening and extension, have a wavefront
from both ends whose work grows with the square of the cost; every byte is a symbol of its own:

```mojo
from dinara_align import Anchor, Band, EndsFree, Ties, affine2p_cigar, affine_cigar, affine_cigars, affine_distance
from dinara_align import affine_extension

var aligned = affine_cigar("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 4, 6, 2)  # cost 12, CIGAR "4=3X6="
var firsts: List[String] = ["ACGTACGT", "TTGCA"]
var seconds: List[String] = ["ACGACGT", "TTGGCA"]
var batch = affine_cigars(firsts, seconds, 4, 6, 2)  # every pair, over every thread
var cost = affine_distance("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 4, 6, 2)  # 12, with no alignment
var capped = affine_cigar("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 4, 6, 2, max_cost=10)  # None: it costs 12
# A read placed anywhere in a reference: the reference's ends free, as WFA2-lib's ends-free mode.
var placed = affine_cigar("TTTTACGTACGTTTTT", "ACGTACGT", 4, 6, 2, ends_free=EndsFree(16, 16, 0, 0))  # 4D8=4D
# Two-piece gap costs, WFA's gap-affine-2p: a gap of k letters costs min(6 + 2k, 24 + k).
var long_gap = affine2p_cigar("ACGTACGTTTTTTTTTTTTTTTTTTTTTTTTTACGTACGT", "ACGTACGTACGTACGT", 4, 6, 2, 24, 1)
# Of equally good alignments a fixed rule picks one: indels placed left, as minimap2 places them, or
# right, WFA2-lib's CIGARs byte for byte.
var left = affine_cigar("ACGTTTTACG", "ACGTTTACG", 4, 6, 2)  # "3=1D6="
var right = affine_cigar("ACGTTTTACG", "ACGTTTACG", 4, 6, 2, ties=Ties.RIGHT)  # "6=1D3="
# Exact within a band of diagonals, KSW2's `w`: the best alignment whose every move stays inside.
var banded = affine_cigar("ACGTACGTTTGCA", "ACGTCGTTTTGCA", 4, 6, 2, band=Band.around(2))
# A seed's extension, fixed at one end and stopping where it scores best, a match earning 1.
var onward = affine_extension("ACGTTGCAAGGCTTTT", "ACGTTGCAAGGCGAGA", 1, 4, 6, 2)  # score 12, "12="
var back = affine_extension("TTTTACGTTGCAAGGC", "GAGAACGTTGCAAGGC", 1, 4, 6, 2, anchor=Anchor.END)
```

A pattern can also be found inside a text, Edlib's infix mode, or at its start, its prefix mode:

```mojo
from dinara_align import edit_search, edit_search_alignment

var hit = edit_search("ACGTCG", "TTTTACGTACGTTTTT")  # hit.distance, hit.start, hit.end
var found = edit_search_alignment("ACGTCG", "TTTTACGTACGTTTTT", prefix=True)
```

Ported from AffineGaps by Ash Vardanian, https://github.com/unum-science/AffineGaps, alignment only;
the edit distance from A*PA by Ragnar Groot Koerkamp and Pesho Ivanov (see NOTICE).
"""

from .alignment import AffineGapCosts, AlignmentMode, AlignmentResult, colorize
from .api import (
    DEFAULT_GAP_EXTENSION,
    DEFAULT_GAP_OPENING,
    DEFAULT_MATCH,
    DEFAULT_MISMATCH,
    DNA_ALPHABET,
    STORED_MATRIX_BUDGET,
    Scoring,
    affine_cigars,
    align,
    alignments,
    edit_alignments,
    edit_distances,
    levenshtein_alignment,
    needleman_wunsch_gotoh_alignment,
    needleman_wunsch_gotoh_score,
    score,
    scores,
    smith_waterman_gotoh_alignment,
    smith_waterman_gotoh_score,
)
from .edit_distance import EditCigar, edit_alignment, edit_cigar, edit_distance
from .edit_search import EditHit, edit_search, edit_search_alignment
from .gap_affine import (
    AffineCigar,
    AffineExtension,
    Anchor,
    Band,
    EndsFree,
    Ties,
    affine2p_cigar,
    affine2p_distance,
    affine2p_extension,
    affine_cigar,
    affine_distance,
    affine_extension,
)
from .common import Device, DeviceScope, GpuSpecs, Placement, hardware_threads
from .errors import AlignmentError, ErrorKind
