# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Exact pairwise alignment of DNA, or any text, from Python: dinara-align's Mojo library.

>>> import dinara_align as da
>>> da.distance("ACGTACGTTTGCA", "ACGTCGTTTTGCA")
2
>>> found = da.align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", da.Costs.affine(4, 6, 2))
>>> found.cost, found.cigar
(12, '4=3X6=')
>>> placed = da.align("TTTTACGTACGTTTTT", "ACGTACGT", da.Costs.affine(4, 6, 2), da.Mode.INFIX)
>>> placed.cigar, placed.reference_start, placed.reference_end
('8=', 4, 12)
>>> da.align("GGGGACGTACGTGGGG", "CCCCACGTACGTCCCC", da.Costs.affine(4, 6, 2), da.Mode.local(2)).score
16

Every call takes a reference and a query as `str` (or ASCII `bytes`), `Costs` that price each edit and
a `Mode` that says which ends of the two the alignment must reach, and returns the least cost
(`distance`), the best score (`score`) or an optimal alignment (`align`). Every answer is exact. A
`Scoring`, an alphabet's substitution table with affine gap scores, takes the place of `Costs` for
tables beyond one match and one mismatch score. Batches (`distances`, `alignments`) take many pairs at
once on the caller's thread, or spread over `threads` threads when asked. No call keeps state between
calls, so any number of threads may call at once.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional, Sequence, Union

from . import _dinara

__all__ = [
    "Aligner",
    "Alignment",
    "Anchor",
    "Band",
    "Costs",
    "LocalScores",
    "Mode",
    "Scoring",
    "align",
    "alignments",
    "distance",
    "distances",
    "local_scores",
    "score",
    "search",
    "Hit",
]

UNBOUNDED = 1 << 60
"""As many free letters as any sequence has, and a band edge past every diagonal."""

Text = Union[str, bytes]
"""A sequence as `str`, or as ASCII `bytes`."""

# A mode's kind, as `dinara_mode` numbers it.
_ENDS, _EXTENSION, _LOCAL, _OVERLAP = 0, 1, 2, 3


def _text(sequence: Text) -> str:
    """`sequence` as the `str` the extension takes, `bytes` decoded as ASCII: refused unless a `str` or
    `bytes`, with `TypeError`, and unless ASCII, with `ValueError`, as positions count letters as bytes."""
    if isinstance(sequence, bytes):
        return sequence.decode("ascii")
    if not isinstance(sequence, str):
        raise TypeError(f"a sequence is a str or bytes, not {type(sequence).__name__}")
    if not sequence.isascii():
        raise ValueError("a sequence holds ASCII letters alone")
    return sequence


def _ints(*values) -> tuple:
    """`values`, refused with `TypeError` unless each is an `int`."""
    for value in values:
        if not isinstance(value, int):
            raise TypeError(f"an integer, not {type(value).__name__}")
    return values


def _call(function, *args):
    """Calls the extension, raising its refusals as `ValueError`."""
    try:
        return function(*args)
    except Exception as error:
        message = str(error)
        raise ValueError(message.removeprefix("dinara-align: ")) from None


@dataclass(frozen=True)
class Costs:
    """What each edit costs, as WFA and minimap2 count it: a substitution `mismatch`, a gap of `k`
    letters `opening + k * extension`, or with a second piece the less of that and `opening2 + k *
    extension2`. Deletions, runs of reference letters alone, may cost their own (`with_deletions`).
    Build through `edit`, `linear`, `affine` or `two_piece`."""

    mismatch: int = 1
    opening: int = 0
    extension: int = 1
    opening2: int = -1
    extension2: int = 0
    deletion_opening: int = 0
    deletion_extension: int = 0
    deletion_opening2: int = -1
    deletion_extension2: int = 0

    @staticmethod
    def edit() -> "Costs":
        """Unit costs: the edit (Levenshtein) distance."""
        return Costs()

    @staticmethod
    def linear(mismatch: int, gap: int) -> "Costs":
        """A substitution `mismatch` and every gapped letter `gap`."""
        return Costs(mismatch, 0, gap)

    @staticmethod
    def affine(mismatch: int, opening: int, extension: int) -> "Costs":
        """A gap of `k` letters `opening + k * extension`: minimap2's `-B4 -O4 -E2` is `affine(4, 4, 2)`."""
        return Costs(mismatch, opening, extension)

    @staticmethod
    def two_piece(mismatch: int, opening: int, extension: int, opening2: int, extension2: int) -> "Costs":
        """A gap the less of two affine costs, minimap2's `-O4,24 -E2,1` as `two_piece(4, 4, 2, 24, 1)`."""
        # A negative `opening2` is one piece to the extension, so it is refused here, as the library does.
        if _ints(extension2)[0] <= 0 or _ints(opening2)[0] < 0:
            raise ValueError(f"second piece {opening2}, {extension2}: its extension must cost")
        return Costs(mismatch, opening, extension, opening2, extension2)

    def with_deletions(self, opening: int, extension: int, opening2: int = -1, extension2: int = 0) -> "Costs":
        """These costs with deletions of their own, the others an insertion's, as bwa's `-O del,ins`."""
        # An extension of zero is no deletions of their own to the extension, so it is refused here.
        if _ints(extension)[0] <= 0 or _ints(opening)[0] < 0 or (_ints(opening2)[0] >= 0 and _ints(extension2)[0] <= 0):
            raise ValueError(f"deletions {opening}, {extension}: an extension must cost")
        return replace(
            self,
            deletion_opening=opening,
            deletion_extension=extension,
            deletion_opening2=opening2,
            deletion_extension2=extension2,
        )

    def _fields(self) -> tuple:
        """The costs as the extension takes them, in `dinara_costs`'s order."""
        return _ints(
            self.mismatch,
            self.opening,
            self.extension,
            self.opening2,
            self.extension2,
            self.deletion_opening,
            self.deletion_extension,
            self.deletion_opening2,
            self.deletion_extension2,
        )


class Anchor:
    """Which end of both sequences an extension is fixed at."""

    START = 0
    END = 1


@dataclass(frozen=True)
class Mode:
    """Which alignments count: `ends_free(...)` and its presets `GLOBAL`, `INFIX` (the whole query inside the
    reference), `PREFIX`, `SUFFIX` and `overlap(...)`; `extension(...)`; and `local(...)` (`local()` under
    a `Scoring`). Free ends minimize costs alone unless `with_match_score` rewards every match."""

    kind: int = _ENDS
    reference_start: int = 0
    reference_end: int = 0
    query_start: int = 0
    query_end: int = 0
    match_score: int = 0
    anchor: int = Anchor.START
    zdrop: int = -1
    end_bonus: int = -1

    GLOBAL = None  # set below
    INFIX = None
    PREFIX = None
    SUFFIX = None

    @staticmethod
    def ends_free(
        reference_start: int = 0,
        reference_end: int = 0,
        query_start: int = 0,
        query_end: int = 0,
    ) -> "Mode":
        """Up to so many letters at each end of each sequence left unaligned for nothing; the query's both
        ends free place the whole reference inside it."""
        return Mode(_ENDS, reference_start, reference_end, query_start, query_end)

    def with_match_score(self, match_score: int) -> "Mode":
        """These free ends with every match earning `match_score`: the best score, not the least cost."""
        return replace(self, match_score=match_score)

    @staticmethod
    def extension(
        match_score: int, anchor: int = Anchor.START, zdrop: Optional[int] = None, end_bonus: Optional[int] = None
    ) -> "Mode":
        """The best-scoring alignment fixed at one end of both, free to stop anywhere: a seed's
        extension, with KSW2's Z-drop when `zdrop` is given, and with `end_bonus` aligned to the query's
        far end when that scores within the bonus of the best stop, as KSW2's end bonus chooses."""
        if zdrop is not None and _ints(zdrop)[0] < 0:
            raise ValueError("a Z-drop below zero")
        if end_bonus is not None and _ints(end_bonus)[0] < 0:
            raise ValueError("an end bonus below zero")
        return Mode(
            _EXTENSION,
            match_score=match_score,
            anchor=anchor,
            zdrop=-1 if zdrop is None else zdrop,
            end_bonus=-1 if end_bonus is None else end_bonus,
        )

    @staticmethod
    def local(match_score: int = 0) -> "Mode":
        """The best-scoring alignment of any part of each, Smith-Waterman, a match earning `match_score`; under a
        `Scoring`, `local()`, its table's own rewards."""
        return Mode(_LOCAL, match_score=match_score)

    @staticmethod
    def overlap(match_score: int) -> "Mode":
        """Every end gap free, a match earning `match_score`: two reads overlapping, or one inside the other."""
        return Mode(_OVERLAP, match_score=match_score)

    def _fields(self) -> tuple:
        """The mode as the extension takes it, in `dinara_mode`'s order, free letters capped at `UNBOUNDED`."""
        _ints(self.kind, self.reference_start, self.reference_end, self.query_start, self.query_end)
        _ints(self.match_score, self.anchor, self.zdrop, self.end_bonus)
        return (
            self.kind,
            min(self.reference_start, UNBOUNDED),
            min(self.reference_end, UNBOUNDED),
            min(self.query_start, UNBOUNDED),
            min(self.query_end, UNBOUNDED),
            self.match_score,
            self.anchor,
            self.zdrop,
            self.end_bonus,
        )


Mode.GLOBAL = Mode()
Mode.INFIX = Mode(_ENDS, UNBOUNDED, UNBOUNDED, 0, 0)
Mode.PREFIX = Mode(_ENDS, 0, UNBOUNDED, 0, 0)
Mode.SUFFIX = Mode(_ENDS, UNBOUNDED, 0, 0, 0)


@dataclass(frozen=True)
class Band:
    """The diagonals every move stays on, `low ..= high`, a cell's diagonal its reference letters less
    its query letters: KSW2's band of width `w` is `Band.around(w)`."""

    low: int = -UNBOUNDED
    high: int = UNBOUNDED

    @staticmethod
    def around(width: int) -> "Band":
        """KSW2's band of width `width`: at most `width` diagonals from the origin's, either way."""
        return Band(-width, width)


@dataclass(frozen=True)
class Scoring:
    """An alphabet's substitution table, row by row in the alphabet's order, and affine gap scores, a
    gap of `k` letters scoring `opening + k * extension`, both zero or less: what an alignment
    maximizes in place of `Costs`."""

    alphabet: str
    table: tuple
    opening: int = -4
    extension: int = -2

    @staticmethod
    def dna() -> "Scoring":
        """minimap2's: match 2, mismatch -4, a gap of `k` letters -(4 + 2k)."""
        return Scoring.uniform(2, -4)

    @staticmethod
    def uniform(match: int, mismatch: int, opening: int = -4, extension: int = -2, alphabet: str = "ACGT") -> "Scoring":
        """One score for equal letters and one for unequal, over `alphabet`."""
        size = len(alphabet)
        table = tuple(match if row == column else mismatch for row in range(size) for column in range(size))
        return Scoring(alphabet, table, opening, extension)

    @staticmethod
    def tabulated(alphabet: str, table: Sequence[int], opening: int = -4, extension: int = -2) -> "Scoring":
        """A table of your own, `len(alphabet)` squared cells, row by row."""
        return Scoring(alphabet, tuple(table), opening, extension)

    def _fields(self) -> tuple:
        """`(alphabet, cells, opening, extension)`, as the extension takes a `Scoring`: each cell an `int`
        of -128 to 127, which the table holds."""
        if not isinstance(self.alphabet, str) or not self.alphabet.isascii():
            raise ValueError("an alphabet is a str of ASCII letters")
        cells = list(_ints(*self.table))
        if any(cell < -128 or cell > 127 for cell in cells):
            raise ValueError("a table's cells lie within -128 and 127")
        _ints(self.opening, self.extension)
        return (self.alphabet, cells, self.opening, self.extension)


@dataclass(frozen=True)
class Alignment:
    """An optimal alignment of `reference[reference_start:reference_end]` against
    `query[query_start:query_end]`, as a CIGAR over those letters alone: `=` a match, `X` a
    substitution (or `M` for either), `D` a reference letter alone, `I` a query letter alone. `cost` is
    what its edits cost, `score` its matches' reward less the cost (minus the cost with no reward)."""

    cost: int
    score: int
    cigar: str
    reference_start: int
    reference_end: int
    query_start: int
    query_end: int

    def runs(self) -> list:
        """The CIGAR as `(length, operation)` pairs."""
        out, length = [], 0
        for letter in self.cigar:
            if letter.isdigit():
                length = length * 10 + int(letter)
            else:
                out.append((length, letter))
                length = 0
        return out

    def gapped(self, reference: Text, query: Text) -> tuple:
        """The aligned parts of both as two rows of one length, `-` against each gapped letter."""
        reference, query = _text(reference), _text(query)
        top, bottom = [], []
        column, row = self.reference_start, self.query_start
        for length, operation in self.runs():
            if operation == "D":
                top.append(reference[column : column + length])
                bottom.append("-" * length)
                column += length
            elif operation == "I":
                top.append("-" * length)
                bottom.append(query[row : row + length])
                row += length
            else:
                top.append(reference[column : column + length])
                bottom.append(query[row : row + length])
                column += length
                row += length
        return "".join(top), "".join(bottom)

    def clipped_cigar(self, query_length: int, hard: bool = False) -> str:
        """The CIGAR as a SAM record writes it, the query's letters outside the span clipped, soft (`S`)
        or with `hard` hard (`H`); the record's position is `reference_start + 1`."""
        clip = "H" if hard else "S"
        out = f"{self.query_start}{clip}" if self.query_start > 0 else ""
        out += self.cigar
        if query_length > self.query_end:
            out += f"{query_length - self.query_end}{clip}"
        return out

    def edit_distance(self, reference: Text, query: Text) -> int:
        """SAM's `NM`: the alignment's substitutions and gapped letters."""
        top, bottom = self.gapped(reference, query)
        return sum(1 for a, b in zip(top, bottom) if a != b)

    def identity(self, reference: Text, query: Text) -> float:
        """Matches over the alignment's columns, BLAST's identity."""
        top, bottom = self.gapped(reference, query)
        return sum(1 for a, b in zip(top, bottom) if a == b) / len(top) if top else 0.0

    def mismatch_string(self, reference: Text, query: Text) -> str:
        """SAM's `MD` tag: the reference's letters the alignment does not match."""
        top, bottom = self.gapped(reference, query)
        out, matched, deleting = [], 0, False
        for a, b in zip(top, bottom):
            if b == "-":
                if not deleting:
                    out.append(f"{matched}^")
                    matched = 0
                    deleting = True
                out.append(a)
                continue
            deleting = False
            if a == "-":
                continue
            if a == b:
                matched += 1
            else:
                out.append(f"{matched}{a}")
                matched = 0
        out.append(str(matched))
        return "".join(out)


@dataclass(frozen=True)
class LocalScores:
    """A local alignment's best score and where it ends, and SSW's second best: the best of an alignment
    ending more than a window of reference letters away."""

    score: int
    reference_end: int
    query_end: int
    second_score: int
    second_reference_end: int


def _options(band: Optional[Band], max_cost: Optional[int], eqx: bool, ties: str, max_memory: Optional[int]) -> tuple:
    """The options as the extension takes them, in `dinara_options`'s order and then whether there is a cap,
    0 for the default memory. Raises `ValueError` for a `ties` other than `"left"` or `"right"`."""
    if ties not in ("left", "right"):
        raise ValueError("ties: 'left' or 'right'")
    band = band or Band()
    _ints(band.low, band.high)
    if max_memory is not None:
        _ints(max_memory)
    return (
        band.low,
        band.high,
        0 if max_cost is None else _ints(max_cost)[0],
        1 if eqx else 0,
        1 if ties == "right" else 0,
        max_memory or 0,
        0 if max_cost is None else 1,
    )


def distance(
    reference: Text,
    query: Text,
    costs: Costs = Costs(),
    mode: Mode = Mode.GLOBAL,
    *,
    band: Optional[Band] = None,
    max_cost: Optional[int] = None,
) -> Optional[int]:
    """The least cost of aligning `query` to `reference` as `mode` asks, or None past `max_cost`."""
    options = _options(band, max_cost, True, "left", None)
    arguments = (_text(reference), _text(query), costs._fields(), mode._fields(), options)
    return _call(_dinara.distance, *arguments)


def align(
    reference: Text,
    query: Text,
    costs: Union[Costs, Scoring] = Costs(),
    mode: Mode = Mode.GLOBAL,
    *,
    band: Optional[Band] = None,
    max_cost: Optional[int] = None,
    ties: str = "left",
    eqx: bool = True,
    max_memory: Optional[int] = None,
) -> Optional[Alignment]:
    """An optimal alignment as `mode` asks, or None past `max_cost`. Of equally good alignments the
    one `ties` names: indels placed left, as minimap2 places them, or `"right"`, WFA2-lib's CIGAR. The
    fronts kept for the traceback stay within `max_memory` bytes, about 80 MB by default. Under a
    `Scoring`, which takes no band, cap, tie rule or memory, Gotoh's walk picks among equals."""
    if isinstance(costs, Scoring):
        found = _call(_dinara.scoring_align, _text(reference), _text(query), costs._fields(), mode._fields(), eqx)
        return Alignment(*found)
    options = _options(band, max_cost, eqx, ties, max_memory)
    arguments = (_text(reference), _text(query), costs._fields(), mode._fields(), options)
    found = _call(_dinara.align, *arguments)
    return None if found is None else Alignment(*found)


class Aligner:
    """One thread's aligner: `distance` and `align` as the functions give them, under `Costs`, the memory
    their searches take kept from call to call, so a loop of calls takes none once it is warm. Keep one
    a thread; the functions keep no state of their own, so any number of aligners work at once."""

    def __init__(self) -> None:
        self._aligner = _dinara.Aligner()

    def distance(
        self,
        reference: Text,
        query: Text,
        costs: Costs = Costs(),
        mode: Mode = Mode.GLOBAL,
        *,
        band: Optional[Band] = None,
        max_cost: Optional[int] = None,
    ) -> Optional[int]:
        """`distance`, through this aligner's memory."""
        options = _options(band, max_cost, True, "left", None)
        arguments = (_text(reference), _text(query), costs._fields(), mode._fields(), options)
        return _call(self._aligner.distance, *arguments)

    def align(
        self,
        reference: Text,
        query: Text,
        costs: Costs = Costs(),
        mode: Mode = Mode.GLOBAL,
        *,
        band: Optional[Band] = None,
        max_cost: Optional[int] = None,
        ties: str = "left",
        eqx: bool = True,
        max_memory: Optional[int] = None,
    ) -> Optional[Alignment]:
        """`align`, through this aligner's memory: the same alignment, the CIGAR its `ties` picks."""
        options = _options(band, max_cost, eqx, ties, max_memory)
        arguments = (_text(reference), _text(query), costs._fields(), mode._fields(), options)
        found = _call(self._aligner.align, *arguments)
        return None if found is None else Alignment(*found)


def score(
    reference: Text,
    query: Text,
    costs: Union[Costs, Scoring] = Costs(),
    mode: Mode = Mode.GLOBAL,
    *,
    band: Optional[Band] = None,
) -> int:
    """The score `align` would return, with no alignment traced: the matches' reward less the costs,
    or minus the least cost with no reward; under a `Scoring`, its score."""
    if isinstance(costs, Scoring):
        return _call(_dinara.scoring_score, _text(reference), _text(query), costs._fields(), mode._fields())
    options = _options(band, None, True, "left", None)
    return _call(_dinara.score, _text(reference), _text(query), costs._fields(), mode._fields(), options)


def local_scores(
    reference: Text, query: Text, costs: Costs, mode: Mode, *, window: Optional[int] = None
) -> LocalScores:
    """For `Mode.local(...)`: the best score and its end, and SSW's second best, more than `window`
    reference letters away, half the query and at least 15 by default."""
    found = _call(
        _dinara.local_scores, _text(reference), _text(query), costs._fields(), mode._fields(), -1 if window is None else window
    )
    return LocalScores(*found)


def distances(
    references: Sequence[Text],
    queries: Sequence[Text],
    costs: Costs = Costs(),
    mode: Mode = Mode.GLOBAL,
    *,
    band: Optional[Band] = None,
    max_cost: Optional[int] = None,
    threads: int = 0,
) -> list:
    """Every pair's `distance`, the pairs spread over `threads` threads, the caller's own alone for zero."""
    options = _options(band, max_cost, True, "left", None)
    arguments = (
        [_text(item) for item in references],
        [_text(item) for item in queries],
        costs._fields(),
        mode._fields(),
    )
    return _call(_dinara.distances, *arguments, options, _ints(threads)[0])


def alignments(
    references: Sequence[Text],
    queries: Sequence[Text],
    costs: Costs = Costs(),
    mode: Mode = Mode.GLOBAL,
    *,
    band: Optional[Band] = None,
    max_cost: Optional[int] = None,
    ties: str = "left",
    eqx: bool = True,
    max_memory: Optional[int] = None,
    threads: int = 0,
) -> list:
    """Every pair's `align`, the pairs spread over `threads` threads, the caller's own alone for zero."""
    options = _options(band, max_cost, eqx, ties, max_memory)
    arguments = (
        [_text(item) for item in references],
        [_text(item) for item in queries],
        costs._fields(),
        mode._fields(),
    )
    found = _call(_dinara.alignments, *arguments, options, _ints(threads)[0])
    return [None if item is None else Alignment(*item) for item in found]


@dataclass(frozen=True)
class Hit:
    """One reference's result in a `search`: its place in the list, its best score (minus its least
    cost with no reward), and its alignment when asked for."""

    index: int
    score: int
    alignment: Optional[Alignment] = None


def search(
    references: Sequence[Text],
    query: Text,
    costs: Costs = Costs(),
    mode: Mode = Mode.GLOBAL,
    *,
    best: Optional[int] = None,
    max_cost: Optional[int] = None,
    aligned: bool = False,
    ties: str = "left",
    threads: int = 0,
) -> list:
    """The query against every reference, the best first: each reference's `Hit`; with `best` that
    many alone, with `max_cost` those within it, with `aligned` their alignments. A local search scores
    a block of references at once, one to a SIMD lane."""
    options = _options(None, max_cost, True, ties, None)
    if best is not None and _ints(best)[0] < 0:
        raise ValueError(f"a best of {best} hits")
    arguments = ([_text(item) for item in references], _text(query), costs._fields(), mode._fields())
    found = _call(_dinara.search, *arguments, options, (-1 if best is None else _ints(best)[0], aligned, _ints(threads)[0]))
    return [Hit(index, score, None if alignment is None else Alignment(*alignment)) for index, score, alignment in found]
