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

Ported from AffineGaps by Ash Vardanian, https://github.com/unum-science/AffineGaps, alignment only.
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
from .edit_distance import edit_alignment, edit_distance
from .common import Device, DeviceScope, GpuSpecs, Placement, hardware_threads
from .errors import AlignmentError, ErrorKind
