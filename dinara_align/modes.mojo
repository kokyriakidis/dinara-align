# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
What an alignment is asked for: `Costs`, the price of each edit, and `Mode`, which ends of the two
sequences it must reach; `Band` and `Ties` narrow and choose among the alignments, and `Alignment` is
what comes back. The first sequence is always the reference and the second the query, as a CIGAR
reads them: `D` a letter of the reference alone, `I` one of the query alone.
"""

from .cigar import cigar_runs
from .errors import AlignmentError, ErrorKind


comptime UNBOUNDED = 1 << 60
"""Past any diagonal or length, and far enough from overflow to shift by any sequence's length."""


@fieldwise_init
struct Costs(Equatable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """What each edit costs, as WFA and minimap2 count it: a substitution `mismatch`, and a gap of `k`
    letters `opening + k extension`, or with a second piece the less of that and `opening2 + k
    extension2`. Every alignment minimizes the total, save an extension (see `Mode.extension`).

    Built through `edit`, `linear`, `affine` or `two_piece`, which refuse costs no search can run by.
    """

    var mismatch: Int
    var opening: Int
    var extension: Int
    var opening2: Int
    """The second gap piece's opening, -1 when there is none."""
    var extension2: Int

    @staticmethod
    def edit() -> Self:
        """Unit costs: the edit (Levenshtein) distance, a substitution, an insertion and a deletion one each.
        A global alignment, or a query found inside or at the start of the reference, takes the
        bit-parallel band doubling of A*PA2; the other modes the wavefront, at the same costs."""
        return Self(1, 0, 1, -1, 0)

    @staticmethod
    def linear(mismatch: Int, gap: Int) raises AlignmentError -> Self:
        """A substitution `mismatch` and every gapped letter `gap`, with no opening."""
        return Self.affine(mismatch, 0, gap)

    @staticmethod
    def affine(mismatch: Int, opening: Int, extension: Int) raises AlignmentError -> Self:
        """Gap-affine costs: a substitution `mismatch`, a gap of `k` letters `opening + k extension`.
        minimap2's `-B4 -O4 -E2` is `affine(4, 4, 2)`, WFA2-lib's default `affine(4, 6, 2)`."""
        if mismatch <= 0 or extension <= 0 or opening < 0:
            raise AlignmentError(
                ErrorKind.INVALID_SCORING,
                String("costs ", mismatch, ", ", opening, ", ", extension, ": a mismatch and an extension must cost"),
            )
        return Self(mismatch, opening, extension, -1, 0)

    @staticmethod
    def two_piece(
        mismatch: Int, opening: Int, extension: Int, opening2: Int, extension2: Int
    ) raises AlignmentError -> Self:
        """Two-piece gap-affine costs, minimap2's `-O4,24 -E2,1` as `two_piece(4, 4, 2, 24, 1)`: a gap of
        `k` letters the less of `opening + k extension` and `opening2 + k extension2`, one piece usually
        cheap to open and the other cheap to extend, so a long gap costs less than one piece charges."""
        var first = Self.affine(mismatch, opening, extension)
        if extension2 <= 0 or opening2 < 0:
            raise AlignmentError(
                ErrorKind.INVALID_SCORING,
                String("second piece ", opening2, ", ", extension2, ": its extension must cost"),
            )
        return Self(first.mismatch, first.opening, first.extension, opening2, extension2)

    def pieces(self) -> Int:
        """How many gap pieces the costs have."""
        return 2 if self.opening2 >= 0 else 1

    def unit_scale(self) -> Int:
        """The factor these costs are of unit costs, zero when they are not: a pair's edit distance
        times it is then their least cost, which the bit-parallel search finds."""
        if self.opening2 >= 0 or self.opening != 0 or self.mismatch != self.extension:
            return 0
        return self.mismatch


@fieldwise_init
struct Anchor(Equatable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """Which end of both sequences an extension is fixed at (see `Mode.extension`)."""

    var identifier: UInt8
    comptime START = Self(0)
    """Both sequences' first letters: the alignment runs right from there, a seed's right extension."""
    comptime END = Self(1)
    """Both sequences' last letters: the alignment runs left from there, a seed's left extension."""


@fieldwise_init
struct Mode(Equatable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """Which alignments of the two sequences count: how many letters at each end of each may be left
    unaligned for nothing, or an extension from one end, or for a `Scoring` a local alignment.

    | mode | reference | query | also called |
    | :-- | :-- | :-- | :-- |
    | `GLOBAL` | whole | whole | end to end, Needleman-Wunsch, Edlib's NW |
    | `INFIX` | any part | whole | semi-global, glocal, Edlib's HW |
    | `PREFIX` | a prefix | whole | Edlib's SHW |
    | `SUFFIX` | a suffix | whole | |
    | `ends_free(...)` | as asked | as asked | WFA2-lib's ends-free, overlaps |
    | `extension(...)` | from one end | from the same end | KSW2's extension, without Z-drop |
    | `local(...)`, `LOCAL` | any part | any part | Smith-Waterman, abPOA's local mode |
    | `overlap(...)` | a prefix or suffix | a suffix or prefix | semi-global, parasail's `sg`, hyalite's OV |

    Free ends minimize the costs alone, as Edlib and WFA2-lib count them, unless a match earns
    something (see `with_match_score`), as parasail's and hyalite's do.
    | `REFERENCE_IN_QUERY` | whole | any part | hyalite's SHW |
    """

    var kind: UInt8
    var reference_start: Int
    """Letters at the reference's start that may go unaligned for nothing; past them a gap is paid."""
    var reference_end: Int
    var query_start: Int
    var query_end: Int
    var match_score: Int
    """What a match earns: in an extension or a local alignment, which maximize a score, and with free
    ends when asked (see `with_match_score`); zero elsewhere, the costs alone minimized."""
    var anchor: Anchor

    comptime ENDS = UInt8(0)
    comptime EXTENSION = UInt8(1)
    comptime SMITH_WATERMAN = UInt8(2)

    comptime GLOBAL = Self(Self.ENDS, 0, 0, 0, 0, 0, Anchor.START)
    """Both sequences end to end."""
    comptime INFIX = Self(Self.ENDS, UNBOUNDED, UNBOUNDED, 0, 0, 0, Anchor.START)
    """The whole query against wherever in the reference it fits best: a read placed in a window."""
    comptime PREFIX = Self(Self.ENDS, 0, UNBOUNDED, 0, 0, 0, Anchor.START)
    """The whole query against the reference's best prefix."""
    comptime SUFFIX = Self(Self.ENDS, UNBOUNDED, 0, 0, 0, 0, Anchor.START)
    """The whole query against the reference's best suffix."""
    comptime REFERENCE_IN_QUERY = Self(Self.ENDS, 0, 0, UNBOUNDED, UNBOUNDED, 0, Anchor.START)
    """`INFIX` the other way round: the whole reference against wherever in the query it fits best."""
    comptime LOCAL = Self(Self.SMITH_WATERMAN, 0, 0, 0, 0, 0, Anchor.START)
    """The best-scoring part of each under a `Scoring`, Smith-Waterman, whose table says what a match
    earns; under `Costs`, `local` names the reward."""

    @staticmethod
    def ends_free(
        *,
        reference_start: Int = 0,
        reference_end: Int = 0,
        query_start: Int = 0,
        query_end: Int = 0,
        match_score: Int = 0,
    ) raises AlignmentError -> Self:
        """Up to so many letters at each end of each sequence left unaligned for nothing, as WFA2-lib's
        ends-free alignment counts them; all zero is `GLOBAL`. An overlap of two reads frees one's start
        and the other's end. With costs alone, freeing both ends of both lets the empty alignment win,
        at no cost; a `match_score` makes the alignment the best-scoring one instead (see
        `with_match_score`)."""
        if min(min(reference_start, reference_end), min(query_start, query_end)) < 0:
            raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a free end of fewer than no letters")
        if match_score < 0:
            raise AlignmentError(ErrorKind.INVALID_SCORING, "a match that costs")
        return Self(Self.ENDS, reference_start, reference_end, query_start, query_end, match_score, Anchor.START)

    def with_match_score(self, match_score: Int) raises AlignmentError -> Self:
        """These free ends with every match earning `match_score`: the best-scoring alignment, the
        reward less the costs, as parasail's and hyalite's semi-global modes count it, rather than the
        least costly one. With some letters left free the two can differ, as a reward pays for
        aligning letters a cost alone would leave out. `Mode.INFIX.with_match_score(2)` is a read placed
        in a window as a mapper scores it. Its time grows with the matrix, as `local`'s does, but for a
        global alignment, whose letters are all aligned, so the reward folds into the costs and the
        wavefront finds it."""
        if self.kind != Self.ENDS:
            raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a match score for free ends alone")
        if match_score < 0:
            raise AlignmentError(ErrorKind.INVALID_SCORING, "a match that costs")
        return Self(
            Self.ENDS,
            self.reference_start,
            self.reference_end,
            self.query_start,
            self.query_end,
            match_score,
            Anchor.START,
        )

    @staticmethod
    def extension(match_score: Int, anchor: Anchor = Anchor.START) raises AlignmentError -> Self:
        """The best-scoring alignment fixed at one end of both sequences, `anchor`, and free to stop
        anywhere: a read mapper's seed extension. A match earns `match_score` and every edit costs what
        `Costs` charges; aligning nothing scores zero. A reward is what makes stopping a choice: with
        costs alone, aligning nothing would always win."""
        if match_score < 0:
            raise AlignmentError(ErrorKind.INVALID_SCORING, "a match that costs")
        return Self(Self.EXTENSION, 0, 0, 0, 0, match_score, anchor)

    @staticmethod
    def local(match_score: Int) raises AlignmentError -> Self:
        """The best-scoring alignment of any part of the reference against any part of the query,
        Smith-Waterman: a match earns `match_score` and every edit costs what `Costs` charges. It is
        every end free, with a reward: with costs alone, aligning nothing would always win. Its time
        grows with the matrix, as every local aligner's does (see `scored`)."""
        if match_score <= 0:
            raise AlignmentError(ErrorKind.INVALID_SCORING, "a local alignment needs a match that earns")
        return Self(Self.SMITH_WATERMAN, 0, 0, 0, 0, match_score, Anchor.START)

    @staticmethod
    def overlap(match_score: Int) raises AlignmentError -> Self:
        """The best-scoring alignment with every end gap free, semi-global: it starts on either
        sequence's first letter and ends on either's last, so one may overhang the other at each end,
        an overlap of two reads, or one may lie inside the other. A match earns `match_score`, as with
        all four ends free and costs alone the empty alignment would win. Its time grows with the
        matrix, as `local`'s does."""
        if match_score <= 0:
            raise AlignmentError(ErrorKind.INVALID_SCORING, "an overlap needs a match that earns")
        return Self(Self.ENDS, UNBOUNDED, UNBOUNDED, UNBOUNDED, UNBOUNDED, match_score, Anchor.START)

    def is_global(self) -> Bool:
        return (
            self.kind == Self.ENDS
            and (self.reference_start | self.reference_end | self.query_start | self.query_end) == 0
        )

    def is_scored(self) -> Bool:
        """Whether the alignment maximizes a score: an extension, a local alignment, or free ends with a
        match that earns."""
        return self.kind != Self.ENDS or self.match_score > 0


@fieldwise_init
struct Ties(Equatable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """Which of several equally good alignments a CIGAR spells. Both rules are WFA2-lib's backtrace: at
    each step back the edit that reached furthest, ties going to a substitution, then a letter of the
    reference alone, then one of the query, the second gap piece before the first and a gap's extension
    before its opening; they differ in the end it runs from. So the CIGAR is the same however the search
    found the cost, whatever the memory limit or band, and whatever the costs' common factor."""

    var identifier: UInt8
    comptime LEFT = Self(0)
    """Every edit as early as an equally good alignment allows, gaps shifted left through repeats: the
    rule run from the start over both sequences reversed, as KSW2 places gaps by default and as variant
    callers normalize indels."""
    comptime RIGHT = Self(1)
    """Every edit as late as it allows, gaps shifted right: WFA2-lib's own CIGARs, byte for byte."""


@fieldwise_init
struct Band(ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """The diagonals an alignment may use, `low ..= high`: a cell's diagonal is the reference's letters
    aligned or skipped up to it less the query's, counted from the alignment's fixed origin, so every
    move from the origin on stays inside. A global alignment starts on diagonal zero and ends on the
    reference's length less the query's, which the band must hold.

    The alignment found is the optimum over every path inside the band, exact, not a heuristic: only
    the diagonals outside are never searched. KSW2's band of width `w` is `Band.around(w)`; WFA2-lib's
    static band counts diagonals the other way, its `min_k ..= max_k` being `Band(-max_k, -min_k)`."""

    var low: Int
    var high: Int

    def __init__(out self):
        """No band: every diagonal."""
        self.low = -UNBOUNDED
        self.high = UNBOUNDED

    @staticmethod
    def around(width: Int) -> Band:
        """The diagonals at most `width` from the origin's, either way."""
        return Band(-width, width)

    def holds(self, diagonal: Int) -> Bool:
        return self.low <= diagonal and diagonal <= self.high

    def covers(self, columns: Int, rows: Int) -> Bool:
        """Whether every diagonal of a `columns` by `rows` matrix lies inside: no band at all for it."""
        return self.low <= -rows and self.high >= columns

    def shifted(self, origin: Int) -> Band:
        """The band as seen from a cell on diagonal `origin`, the start of a piece after a split."""
        return Band(self.low - origin, self.high - origin)

    def mirrored(self, target: Int) -> Band:
        """The band in the reversed sequences of a pair whose end lies on diagonal `target`, where
        diagonal `k` reads `target - k`."""
        return Band(target - self.high, target - self.low)


@fieldwise_init
struct Alignment(Copyable, Movable, Writable):
    """An optimal alignment: of `reference[reference_start:reference_end]` against
    `query[query_start:query_end]`, as a CIGAR over those letters alone, `=` a match and `X` a
    substitution (or `M` for either), `D` a reference letter alone and `I` a query letter alone, each run
    its length then its letter. A global alignment spans both; letters a mode leaves unaligned for
    nothing lie outside the spans, and a gap past them stays in the CIGAR, as it is paid.

    `cost` is what the CIGAR's edits cost. `score` is an extension's: the matches' reward less the cost;
    outside an extension, with no reward, it is minus the cost."""

    var cost: Int
    var score: Int
    var cigar: String
    var reference_start: Int
    var reference_end: Int
    var query_start: Int
    var query_end: Int

    def gapped(self, reference: String, query: String) -> Tuple[String, String]:
        """The aligned parts of both sequences as two rows of one length, `-` against each gapped letter."""
        var top = List[UInt8]()
        var bottom = List[UInt8]()
        var first = reference.as_bytes()
        var second = query.as_bytes()
        var column = self.reference_start
        var row = self.query_start
        comptime GAP = UInt8(ord("-"))
        var runs = cigar_runs(self.cigar)
        for index in range(len(runs[0])):
            var letter = runs[0][index]
            for _ in range(runs[1][index]):
                if letter == UInt8(ord("D")):
                    top.append(first[column])
                    bottom.append(GAP)
                    column += 1
                elif letter == UInt8(ord("I")):
                    top.append(GAP)
                    bottom.append(second[row])
                    row += 1
                else:
                    top.append(first[column])
                    bottom.append(second[row])
                    column += 1
                    row += 1
        return (String(unsafe_from_utf8=top^), String(unsafe_from_utf8=bottom^))
