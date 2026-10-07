"""
The full matrix as an oracle: Gotoh's recurrence over every cell, written for clarity and sharing no
code with the library, under the most general model any entry point takes, and a check that an
`Alignment` is one the model allows and earns what the oracle says is best.

A `Model` is a score to maximize: a substitution scores what a table says, a match's reward or minus a
mismatch's cost, and a gap of `k` letters costs the least over its direction's pieces of `opening + k
extension`, a deletion (a reference letter alone, `D`) and an insertion (a query letter alone, `I`)
each with pieces of their own. The alignment's ends are what the kind allows: free letters at each end
(`ENDS`), one end fixed and the other anywhere (`EXTENSION`), or anywhere with every cell floored at
zero (`LOCAL`); every cell of the path lies on the band's diagonals, `i - j` counted from the fixed end.
"""

comptime ENDS = 0
comptime EXTENSION = 1
comptime LOCAL = 2
comptime LOW = -(1 << 40)
"""A cell no path reaches."""


@fieldwise_init
struct Model(Copyable, Movable, Writable):
    var table: List[Int]
    """What aligning byte `a` against byte `b` scores, at `a * 256 + b`."""
    var deletion: List[Tuple[Int, Int]]
    """A deletion's pieces, each an opening and an extension."""
    var insertion: List[Tuple[Int, Int]]
    var kind: Int
    var reference_start: Int
    var reference_end: Int
    var query_start: Int
    var query_end: Int
    var at_end: Bool
    """An extension fixed at both sequences' ends rather than their starts."""
    var band_low: Int
    var band_high: Int

    @staticmethod
    def uniform(reward: Int, mismatch: Int) -> List[Int]:
        """A table of one reward for equal bytes and one cost for unequal ones."""
        var table = List[Int](length=256 * 256, fill=-mismatch)
        for byte in range(256):
            table[byte * 256 + byte] = reward
        return table^

    def gap(self, letters: Int, deleted: Bool) -> Int:
        """What a gap of `letters` costs, at its direction's cheapest piece."""
        ref pieces = self.deletion if deleted else self.insertion
        var cheapest = 1 << 50
        for piece in pieces:
            cheapest = min(cheapest, piece[0] + piece[1] * letters)
        return cheapest

    def holds(self, diagonal: Int) -> Bool:
        return self.band_low <= diagonal and diagonal <= self.band_high


def reversed_bytes(text: String) -> String:
    var bytes = text.as_bytes()
    var out = List[UInt8](capacity=len(bytes))
    for index in range(len(bytes) - 1, -1, -1):
        out.append(bytes[index])
    return String(unsafe_from_utf8=out^)


def scores(model: Model, reference: String, query: String) -> List[Int]:
    """Every cell's best score from a start the model allows, `len(query) + 1` cells a row of the
    reference's letters, `LOW` where no path on the band reaches. An extension from the end is not
    taken here: its caller reverses both sequences."""
    var a = reference.as_bytes()
    var b = query.as_bytes()
    var n = len(a)
    var m = len(b)
    var width = m + 1
    var cells = (n + 1) * width
    var best = List[Int](length=cells, fill=LOW)
    # Each gap piece's layer either way: a run that may go on for an extension alone.
    var deleting = List[List[Int]]()
    var inserting = List[List[Int]]()
    for _ in model.deletion:
        deleting.append(List[Int](length=cells, fill=LOW))
    for _ in model.insertion:
        inserting.append(List[Int](length=cells, fill=LOW))
    var local = model.kind == LOCAL
    for i in range(n + 1):
        for j in range(m + 1):
            var at = i * width + j
            if not local and not model.holds(i - j):
                continue
            var value = LOW
            # A start: the origin, for free ends a cell on the first row or column within the letters
            # free there, and for a local alignment any cell. Past the free letters the edge is a gap,
            # reached by the recurrence through cells that must lie on the band too.
            if local or (i == 0 and j == 0):
                value = 0
            elif model.kind == ENDS and j == 0 and i <= model.reference_start:
                value = 0
            elif model.kind == ENDS and i == 0 and j <= model.query_start:
                value = 0
            for piece in range(len(model.deletion)):
                if i > 0:
                    var opening = model.deletion[piece][0]
                    var extension = model.deletion[piece][1]
                    deleting[piece][at] = max(
                        best[at - width] - opening - extension, deleting[piece][at - width] - extension
                    )
                    value = max(value, deleting[piece][at])
            for piece in range(len(model.insertion)):
                if j > 0:
                    var opening = model.insertion[piece][0]
                    var extension = model.insertion[piece][1]
                    inserting[piece][at] = max(best[at - 1] - opening - extension, inserting[piece][at - 1] - extension)
                    value = max(value, inserting[piece][at])
            if i > 0 and j > 0 and best[at - width - 1] > LOW:
                value = max(value, best[at - width - 1] + model.table[Int(a[i - 1]) * 256 + Int(b[j - 1])])
            best[at] = max(value, LOW)
    return best^


def is_end(model: Model, i: Int, j: Int, n: Int, m: Int) -> Bool:
    """Whether an alignment may end at cell `(i, j)`: anywhere for a local alignment or an extension,
    else on the last row or column within the letters free there."""
    if model.kind != ENDS:
        return True
    return (i == n and m - j <= model.query_end) or (j == m and n - i <= model.reference_end)


def optimum(model: Model, reference: String, query: String) -> Optional[Int]:
    """The best score the model allows, or None when no alignment stays inside the band."""
    if model.kind == EXTENSION and model.at_end:
        var forward = Model(
            model.table.copy(),
            model.deletion.copy(),
            model.insertion.copy(),
            EXTENSION,
            0,
            0,
            0,
            0,
            False,
            model.band_low,
            model.band_high,
        )
        return optimum(forward, reversed_bytes(reference), reversed_bytes(query))
    var n = reference.byte_length()
    var m = query.byte_length()
    var best = scores(model, reference, query)
    var answer = LOW
    for i in range(n + 1):
        for j in range(m + 1):
            if is_end(model, i, j, n, m):
                answer = max(answer, best[i * (m + 1) + j])
    if answer <= LOW // 2:
        return None
    return answer


def rule_span(model: Model, reference: String, query: String, left: Bool) -> Tuple[Int, Int, Int, Int]:
    """The span the tie rule names for free ends: with `left`, of the ends an optimal alignment
    reaches the one on the highest diagonal, `i - j`, and of the starts an optimal alignment ending
    there leaves from, the one on the highest diagonal too; otherwise that rule over both sequences
    reversed, read back. The reference's start and end, then the query's."""
    var n = reference.byte_length()
    var m = query.byte_length()
    if not left:
        var mirrored = Model(
            model.table.copy(),
            model.deletion.copy(),
            model.insertion.copy(),
            model.kind,
            model.reference_end,
            model.reference_start,
            model.query_end,
            model.query_start,
            False,
            (n - m) - model.band_high,
            (n - m) - model.band_low,
        )
        var span = rule_span(mirrored, reversed_bytes(reference), reversed_bytes(query), True)
        return (n - span[1], n - span[0], m - span[3], m - span[2])
    var width = m + 1
    var best = scores(model, reference, query)
    var top = LOW
    var end = (0, 0)
    for i in range(n + 1):
        for j in range(m + 1):
            if not is_end(model, i, j, n, m):
                continue
            var value = best[i * width + j]
            if value > top or (value == top and i - j > end[0] - end[1]):
                top = value
                end = (i, j)
    # Back from that end alone, over both reversed, to the starts: the end fixed, as the origin.
    var back = Model(
        model.table.copy(),
        model.deletion.copy(),
        model.insertion.copy(),
        EXTENSION,
        0,
        0,
        0,
        0,
        False,
        (end[0] - end[1]) - model.band_high,
        (end[0] - end[1]) - model.band_low,
    )
    var head = reversed_bytes(String(StringSlice(unsafe_from_utf8=reference.as_bytes()[: end[0]])))
    var lead = reversed_bytes(String(StringSlice(unsafe_from_utf8=query.as_bytes()[: end[1]])))
    var behind = scores(back, head, lead)
    var back_width = end[1] + 1
    var start = (-1, -1)
    for i in range(end[0] + 1):
        for j in range(end[1] + 1):
            var column = end[0] - i
            var row = end[1] - j
            if not ((row == 0 and column <= model.reference_start) or (column == 0 and row <= model.query_start)):
                continue
            if behind[i * back_width + j] != top:
                continue
            if start[0] < 0 or column - row > start[0] - start[1]:
                start = (column, row)
    return (start[0], end[0], start[1], end[1])


@fieldwise_init
struct Priced(ImplicitlyCopyable, Writable):
    """What an alignment's CIGAR earns under a model."""

    var score: Int
    var cost: Int
    """Its edits' costs alone."""
    var matches: Int


def priced(model: Model, reference: String, query: String, cigar: String, start: Tuple[Int, Int]) raises -> Priced:
    """Walks `cigar` from `start` over both sequences, refusing a step off either or off the band, an
    `=` over unequal letters or an `X` over equal ones, and prices it."""
    var a = reference.as_bytes()
    var b = query.as_bytes()
    var column = start[0]
    var row = start[1]
    var score = 0
    var cost = 0
    var matches = 0
    var length = 0
    # From the end, the band counts diagonals back from the corner.
    var mirrored = model.kind == EXTENSION and model.at_end
    var corner = len(a) - len(b)

    def inside(column: Int, row: Int) {imm model, imm mirrored, imm corner} -> Bool:
        if model.kind == LOCAL:
            return True
        return model.holds(corner - (column - row)) if mirrored else model.holds(column - row)

    if not inside(column, row):
        raise Error("the alignment starts off the band")
    for byte in cigar.as_bytes():
        if byte >= UInt8(ord("0")) and byte <= UInt8(ord("9")):
            length = length * 10 + Int(byte - UInt8(ord("0")))
            continue
        if length == 0:
            raise Error("an empty run")
        if byte == UInt8(ord("D")) or byte == UInt8(ord("I")):
            var deleted = byte == UInt8(ord("D"))
            var paid = model.gap(length, deleted)
            score -= paid
            cost += paid
            for _ in range(length):
                if deleted:
                    column += 1
                else:
                    row += 1
                if column > len(a) or row > len(b):
                    raise Error("a gap past a sequence's end")
                if not inside(column, row):
                    raise Error("a gap off the band")
        elif byte == UInt8(ord("=")) or byte == UInt8(ord("X")) or byte == UInt8(ord("M")):
            for _ in range(length):
                if column >= len(a) or row >= len(b):
                    raise Error("a pair past a sequence's end")
                var equal = a[column] == b[row]
                if byte == UInt8(ord("=")) and not equal:
                    raise Error("an = over unequal letters")
                if byte == UInt8(ord("X")) and equal:
                    raise Error("an X over equal letters")
                var earned = model.table[Int(a[column]) * 256 + Int(b[row])]
                score += earned
                if equal:
                    matches += 1
                else:
                    cost -= earned
                column += 1
                row += 1
                if not inside(column, row):
                    raise Error("a pair off the band")
        else:
            raise Error(String("an unknown operation ", chr(Int(byte))))
        length = 0
    if length != 0:
        raise Error("a run with no operation")
    return Priced(score, cost, matches)


def check_alignment(
    model: Model,
    reference: String,
    query: String,
    cigar: String,
    reference_start: Int,
    reference_end: Int,
    query_start: Int,
    query_end: Int,
) raises -> Priced:
    """Raises unless the alignment spells a path between its spans, starts and ends where the model
    allows, and stays on the band; returns what it earns, its free letters outside the spans earning
    nothing."""
    var n = reference.byte_length()
    var m = query.byte_length()
    if reference_start < 0 or reference_start > reference_end or reference_end > n:
        raise Error("a reference span out of order")
    if query_start < 0 or query_start > query_end or query_end > m:
        raise Error("a query span out of order")
    var found = priced(model, reference, query, cigar, (reference_start, query_start))
    # The walk must end exactly at both spans' ends.
    var columns = 0
    var rows = 0
    var length = 0
    for byte in cigar.as_bytes():
        if byte >= UInt8(ord("0")) and byte <= UInt8(ord("9")):
            length = length * 10 + Int(byte - UInt8(ord("0")))
            continue
        if byte != UInt8(ord("I")):
            columns += length
        if byte != UInt8(ord("D")):
            rows += length
        length = 0
    if reference_start + columns != reference_end or query_start + rows != query_end:
        raise Error("the CIGAR does not span the spans")
    if model.kind == EXTENSION:
        if model.at_end and (reference_end != n or query_end != m):
            raise Error("an extension from the end that does not reach it")
        if not model.at_end and (reference_start != 0 or query_start != 0):
            raise Error("an extension from the start that does not begin there")
    elif model.kind == ENDS:
        var starts = (reference_start <= model.reference_start and query_start == 0) or (
            query_start <= model.query_start and reference_start == 0
        )
        var ends = (n - reference_end <= model.reference_end and query_end == m) or (
            m - query_end <= model.query_end and reference_end == n
        )
        if not starts:
            raise Error("a start past the letters free there")
        if not ends:
            raise Error("an end past the letters free there")
    return found
