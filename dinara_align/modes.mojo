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
    extension2`. A deletion, a run of reference letters alone, may cost otherwise than an insertion,
    a run of query letters alone, as bwa's `-O del,ins` (see `with_deletions`). Every alignment
    minimizes the total, save an extension (see `Mode.extension`).

    Built through `edit`, `linear`, `affine` or `two_piece`, which refuse costs no search can run by.
    """

    var mismatch: Int
    var opening: Int
    """An insertion's opening, and a deletion's unless `with_deletions` set its own."""
    var extension: Int
    var opening2: Int
    """The second gap piece's opening, -1 when there is none."""
    var extension2: Int
    var deletion_opening: Int
    """A deletion's own opening, its extension and its second piece's after it: `opening` and the
    rest unless `with_deletions` set them."""
    var deletion_extension: Int
    var deletion_opening2: Int
    var deletion_extension2: Int

    def __init__(out self, mismatch: Int, opening: Int, extension: Int, opening2: Int, extension2: Int):
        """Costs whose deletions cost what insertions do, trusted as given: the factories check them."""
        self.mismatch = mismatch
        self.opening = opening
        self.extension = extension
        self.opening2 = opening2
        self.extension2 = extension2
        self.deletion_opening = opening
        self.deletion_extension = extension
        self.deletion_opening2 = opening2
        self.deletion_extension2 = extension2

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

    def with_deletions(
        self, opening: Int, extension: Int, opening2: Int = -1, extension2: Int = 0
    ) raises AlignmentError -> Self:
        """These costs with a deletion, a run of `k` reference letters alone, costing `opening + k
        extension`, or the less of that and `opening2 + k extension2`, and an insertion as before:
        bwa's `-O6,5 -E1,2` is `affine(4, 5, 2).with_deletions(6, 1)`, its insertions first. Either side
        may have a second piece the other lacks: the one without counts its one piece twice."""
        if extension <= 0 or opening < 0 or (opening2 >= 0 and extension2 <= 0):
            raise AlignmentError(
                ErrorKind.INVALID_SCORING,
                String("deletions ", opening, ", ", extension, ": an extension must cost"),
            )
        var out = self
        out.deletion_opening = opening
        out.deletion_extension = extension
        out.deletion_opening2 = opening2 if opening2 >= 0 else opening
        out.deletion_extension2 = extension2 if opening2 >= 0 else extension
        if opening2 >= 0 and out.opening2 < 0:
            # Insertions gain the second piece deletions have, their own once more.
            out.opening2 = out.opening
            out.extension2 = out.extension
        elif opening2 < 0 and out.opening2 < 0:
            out.deletion_opening2 = -1
            out.deletion_extension2 = 0
        return out

    def pieces(self) -> Int:
        """How many gap pieces the costs have, either side."""
        return 2 if self.opening2 >= 0 else 1

    def symmetric(self) -> Bool:
        """Whether a deletion costs what an insertion does."""
        return (
            self.deletion_opening == self.opening
            and self.deletion_extension == self.extension
            and self.deletion_opening2 == self.opening2
            and self.deletion_extension2 == self.extension2
        )

    def unit_scale(self) -> Int:
        """The factor these costs are of unit costs, zero when they are not: a pair's edit distance
        times it is then their least cost, which the bit-parallel search finds."""
        if self.opening2 >= 0 or self.opening != 0 or self.mismatch != self.extension or not self.symmetric():
            return 0
        return self.mismatch

    def gap(self, letters: Int, deleted: Bool) -> Int:
        """What a gap of `letters` letters costs, a deletion or an insertion, at its cheaper piece."""
        if letters == 0:
            return 0
        var opening = self.deletion_opening if deleted else self.opening
        var extension = self.deletion_extension if deleted else self.extension
        var cost = opening + extension * letters
        if self.opening2 >= 0:
            var opening2 = self.deletion_opening2 if deleted else self.opening2
            var extension2 = self.deletion_extension2 if deleted else self.extension2
            cost = min(cost, opening2 + extension2 * letters)
        return cost


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
    | `extension(...)` | from one end | from the same end | KSW2's extension, with or without Z-drop |
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
    var zdrop: Int
    """An extension's Z-drop, -1 for none (see `extension`)."""

    comptime ENDS = UInt8(0)
    comptime EXTENSION = UInt8(1)
    comptime SMITH_WATERMAN = UInt8(2)

    comptime GLOBAL = Self(Self.ENDS, 0, 0, 0, 0, 0, Anchor.START, -1)
    """Both sequences end to end."""
    comptime INFIX = Self(Self.ENDS, UNBOUNDED, UNBOUNDED, 0, 0, 0, Anchor.START, -1)
    """The whole query against wherever in the reference it fits best: a read placed in a window."""
    comptime PREFIX = Self(Self.ENDS, 0, UNBOUNDED, 0, 0, 0, Anchor.START, -1)
    """The whole query against the reference's best prefix."""
    comptime SUFFIX = Self(Self.ENDS, UNBOUNDED, 0, 0, 0, 0, Anchor.START, -1)
    """The whole query against the reference's best suffix."""
    comptime REFERENCE_IN_QUERY = Self(Self.ENDS, 0, 0, UNBOUNDED, UNBOUNDED, 0, Anchor.START, -1)
    """`INFIX` the other way round: the whole reference against wherever in the query it fits best."""
    comptime LOCAL = Self(Self.SMITH_WATERMAN, 0, 0, 0, 0, 0, Anchor.START, -1)
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
        return Self(Self.ENDS, reference_start, reference_end, query_start, query_end, match_score, Anchor.START, -1)

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
            -1,
        )

    @staticmethod
    def extension(
        match_score: Int, anchor: Anchor = Anchor.START, *, zdrop: Optional[Int] = None
    ) raises AlignmentError -> Self:
        """The best-scoring alignment fixed at one end of both sequences, `anchor`, and free to stop
        anywhere: a read mapper's seed extension. A match earns `match_score` and every edit costs what
        `Costs` charges; aligning nothing scores zero. A reward is what makes stopping a choice: with
        costs alone, aligning nothing would always win.

        Exact by default: the best stop of all. With `zdrop`, minimap2's `-z` and KSW2's Z-drop, the
        search gives up once every alignment it is growing scores more than `zdrop`, plus a gap
        extension a diagonal between them, below the best so far, as WFA2-lib's Z-drop gauges it a cost
        at a time, and the best stop it found stands: a heuristic, faster on a seed whose read turns to
        noise, which may miss a better stop past a divergent stretch."""
        if match_score < 0:
            raise AlignmentError(ErrorKind.INVALID_SCORING, "a match that costs")
        var drop = zdrop.or_else(-1)
        if zdrop and drop < 0:
            raise AlignmentError(ErrorKind.INVALID_ARGUMENT, "a Z-drop below zero")
        return Self(Self.EXTENSION, 0, 0, 0, 0, match_score, anchor, drop)

    @staticmethod
    def local(match_score: Int) raises AlignmentError -> Self:
        """The best-scoring alignment of any part of the reference against any part of the query,
        Smith-Waterman: a match earns `match_score` and every edit costs what `Costs` charges. It is
        every end free, with a reward: with costs alone, aligning nothing would always win. Its time
        grows with the matrix, as every local aligner's does (see `scored`)."""
        if match_score <= 0:
            raise AlignmentError(ErrorKind.INVALID_SCORING, "a local alignment needs a match that earns")
        return Self(Self.SMITH_WATERMAN, 0, 0, 0, 0, match_score, Anchor.START, -1)

    @staticmethod
    def overlap(match_score: Int) raises AlignmentError -> Self:
        """The best-scoring alignment with every end gap free, semi-global: it starts on either
        sequence's first letter and ends on either's last, so one may overhang the other at each end,
        an overlap of two reads, or one may lie inside the other. A match earns `match_score`, as with
        all four ends free and costs alone the empty alignment would win. Its time grows with the
        matrix, as `local`'s does."""
        if match_score <= 0:
            raise AlignmentError(ErrorKind.INVALID_SCORING, "an overlap needs a match that earns")
        return Self(Self.ENDS, UNBOUNDED, UNBOUNDED, UNBOUNDED, UNBOUNDED, match_score, Anchor.START, -1)

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
    found the cost, whatever the band or cap, and whatever the costs' common factor.

    With free ends the span comes first. `LEFT` decides it from the end back, as it places edits: the
    end on the highest diagonal an equally good alignment reaches, the most reference letters less query
    letters, the furthest along the reference for a read placed in it, then the start on the highest
    diagonal of those ending there; `RIGHT` is that over both sequences reversed, the start on the
    lowest diagonal, then the end. The letters between are then aligned globally by the rule. A local
    alignment likewise ends as late as an equally good one allows and starts as late too under `LEFT`,
    and under `RIGHT` starts and ends as early (see `scored`)."""

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

    def clipped_cigar(self, query_length: Int, *, hard: Bool = False) -> String:
        """The CIGAR as a SAM record writes it: the query's letters outside the span clipped, soft (`S`),
        their letters kept in the record's sequence, or with `hard` hard (`H`), dropped from it. The
        record's position is `reference_start + 1`; the reference's letters outside the span need no
        operation."""
        var clip = "H" if hard else "S"
        var out = String()
        if self.query_start > 0:
            out += String(self.query_start, clip)
        out += self.cigar
        if query_length > self.query_end:
            out += String(query_length - self.query_end, clip)
        return out

    def counts(self, reference: String, query: String) -> AlignedCounts:
        """How many letters the alignment pairs equal and unequal, and leaves gapped either way, `M` runs
        compared letter by letter."""
        var first = reference.as_bytes()
        var second = query.as_bytes()
        var column = self.reference_start
        var row = self.query_start
        var counted = AlignedCounts(0, 0, 0, 0)
        var runs = cigar_runs(self.cigar)
        for index in range(len(runs[0])):
            var letter = runs[0][index]
            var length = runs[1][index]
            if letter == UInt8(ord("D")):
                counted.deleted += length
                column += length
            elif letter == UInt8(ord("I")):
                counted.inserted += length
                row += length
            else:
                for _ in range(length):
                    if first[column] == second[row]:
                        counted.matches += 1
                    else:
                        counted.mismatches += 1
                    column += 1
                    row += 1
        return counted

    def edit_distance(self, reference: String, query: String) -> Int:
        """The alignment's edits, SAM's `NM` tag: its substitutions and gapped letters."""
        var counted = self.counts(reference, query)
        return counted.mismatches + counted.deleted + counted.inserted

    def identity(self, reference: String, query: String) -> Float64:
        """The matches over the alignment's columns, BLAST's identity: one for an exact match, zero for
        an empty alignment."""
        var counted = self.counts(reference, query)
        var columns = counted.matches + counted.mismatches + counted.deleted + counted.inserted
        if columns == 0:
            return 0
        return Float64(counted.matches) / Float64(columns)

    def mismatch_string(self, reference: String, query: String) -> String:
        """SAM's `MD` tag: the reference's letters the alignment does not match, each substitution's
        letter after the matches before it and each deletion's after a `^`, so the reference can be
        rebuilt from the query and the CIGAR; insertions leave no mark."""
        var first = reference.as_bytes()
        var second = query.as_bytes()
        var column = self.reference_start
        var row = self.query_start
        var out = String()
        var matched = 0
        var runs = cigar_runs(self.cigar)
        for index in range(len(runs[0])):
            var letter = runs[0][index]
            var length = runs[1][index]
            if letter == UInt8(ord("I")):
                row += length
            elif letter == UInt8(ord("D")):
                out += String(matched, "^")
                matched = 0
                out += String(StringSlice(unsafe_from_utf8=first[column : column + length]))
                column += length
            else:
                for _ in range(length):
                    if first[column] == second[row]:
                        matched += 1
                    else:
                        out += String(matched)
                        matched = 0
                        out += String(StringSlice(unsafe_from_utf8=first[column : column + 1]))
                    column += 1
                    row += 1
        out += String(matched)
        return out


@fieldwise_init
struct AlignedCounts(Equatable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """An alignment's columns by kind (see `Alignment.counts`)."""

    var matches: Int
    var mismatches: Int
    var deleted: Int
    """Reference letters against a gap, `D`."""
    var inserted: Int
    """Query letters against a gap, `I`."""
