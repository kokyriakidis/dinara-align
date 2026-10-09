# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
A CIGAR read back: its runs, reversed, priced and its matches counted, as the modes put alignments
together from pieces (see `api`, `scored`). Writing one from traceback moves is the traceback's own
(see `traceback.cigar_string`).
"""


def cigar_runs(cigar: String) -> Tuple[List[UInt8], List[Int]]:
    """A CIGAR's letters and run lengths."""
    var letters = List[UInt8]()
    var lengths = List[Int]()
    var length = 0
    for byte in cigar.as_bytes():
        if byte >= UInt8(ord("0")) and byte <= UInt8(ord("9")):
            length = length * 10 + Int(byte - UInt8(ord("0")))
            continue
        letters.append(byte)
        lengths.append(length)
        length = 0
    return (letters^, lengths^)


def joined_cigar(letters: List[UInt8], lengths: List[Int]) -> String:
    """Runs written back as a CIGAR, the empty ones left out and neighbours of one move made one."""
    # A run's digits and its letter, 21 bytes at most.
    var writer = CigarWriter(21 * len(letters) + 1)
    for index in range(len(letters)):
        if lengths[index] > 0:
            writer.add(letters[index], lengths[index])
    return writer^.finish()


def reversed_cigar(cigar: String) -> String:
    """A CIGAR's runs in the other order: the alignment of both sequences reversed."""
    var runs = cigar_runs(cigar)
    var letters = List[UInt8](capacity=len(runs[0]))
    var lengths = List[Int](capacity=len(runs[1]))
    for index in range(len(runs[0]) - 1, -1, -1):
        letters.append(runs[0][index])
        lengths.append(runs[1][index])
    return joined_cigar(letters, lengths)


def text_of(bytes: ImmSpan[UInt8, _]) -> String:
    """`bytes` as a sequence's text, whatever they hold: a C caller's bytes need not be UTF-8."""
    return String(StringSlice(unsafe_from_utf8=bytes))


def reversed_text(bytes: ImmSpan[UInt8, _]) -> String:
    """`bytes` back to front."""
    return String(unsafe_from_utf8=reversed_list(bytes))


def reversed_list(bytes: ImmSpan[UInt8, _]) -> List[UInt8]:
    """`bytes` back to front."""
    var out = List[UInt8](capacity=len(bytes))
    out.resize(unsafe_uninit_length=len(bytes))
    reversed_into(out.unsafe_ptr(), bytes.unsafe_ptr(), len(bytes))
    return out^


def append_reversed(mut out: List[UInt8], bytes: ImmSpan[UInt8, _]):
    """`bytes` back to front, after `out`'s own."""
    var start = len(out)
    out.resize(unsafe_uninit_length=start + len(bytes))
    reversed_into(out.unsafe_ptr().unsafe_offset(start), bytes.unsafe_ptr(), len(bytes))


@always_inline
def reversed_into(destination: MutPointer[UInt8, _], source: ImmPointer[UInt8, _], count: Int):
    """`count` bytes from `source` into `destination` back to front, sixteen at a time; the two must not
    overlap."""
    comptime CHUNK = 16
    var index = 0
    while index + CHUNK <= count:
        destination.unsafe_offset(index).unsafe_store(
            source.unsafe_offset(count - index - CHUNK).unsafe_load[width=CHUNK]().reversed()
        )
        index += CHUNK
    while index < count:
        destination[unsafe_offset=index] = source[unsafe_offset=count - 1 - index]
        index += 1


def reverse_bytes(bytes: MutPointer[UInt8, _], count: Int):
    """Turns `count` bytes back to front in place, sixteen from either end at a time."""
    comptime CHUNK = 16
    var front = 0
    var back = count
    while back - front >= 2 * CHUNK:
        var head = bytes.unsafe_offset(front).unsafe_load[width=CHUNK]()
        var tail = bytes.unsafe_offset(back - CHUNK).unsafe_load[width=CHUNK]()
        bytes.unsafe_offset(front).unsafe_store(tail.reversed())
        bytes.unsafe_offset(back - CHUNK).unsafe_store(head.reversed())
        front += CHUNK
        back -= CHUNK
    while back - front >= 2:
        back -= 1
        var kept = bytes[unsafe_offset=front]
        bytes[unsafe_offset=front] = bytes[unsafe_offset=back]
        bytes[unsafe_offset=back] = kept
        front += 1


@fieldwise_init
struct AlignedCounts(Equatable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """An alignment's columns by kind (see `cigar_counts`)."""

    var matches: Int
    var mismatches: Int
    var deleted: Int
    """Reference letters against a gap, `D`."""
    var inserted: Int
    """Query letters against a gap, `I`."""


def cigar_counts(
    reference: ImmSpan[UInt8, _],
    query: ImmSpan[UInt8, _],
    cigar: String,
    reference_start: Int = 0,
    query_start: Int = 0,
) -> AlignedCounts:
    """How many letters a CIGAR of `reference` against `query`, from `reference_start` and `query_start`, pairs
    equal and unequal and leaves gapped either way, `M` runs compared letter by letter."""
    var column = reference_start
    var row = query_start
    var counted = AlignedCounts(0, 0, 0, 0)
    var runs = cigar_runs(cigar)
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
                if reference[column] == query[row]:
                    counted.matches += 1
                else:
                    counted.mismatches += 1
                column += 1
                row += 1
    return counted


def cigar_matches(first: String, second: String, cigar: String) -> Int:
    """The equal pairs a CIGAR of `first` against `second` aligns (see `cigar_counts`)."""
    return cigar_counts(first.as_bytes(), second.as_bytes(), cigar).matches


struct CigarWriter:
    """A CIGAR string written a run at a time into bytes reserved once, a run of the same letter as the
    last one joining it; each length's digits are written by hand, as formatting one through a `String`,
    or growing the bytes a run at a time, took longer than the gapped rows' whole copy on short reads."""

    var text: List[UInt8]
    var used: Int
    var letter: UInt8
    """The letter of the run still being added to, zero before the first."""
    var length: Int
    """The length of that run so far."""

    def __init__(out self, capacity: Int):
        """Room for `capacity` bytes, which the caller bounds: nothing past it is checked."""
        self.text = List[UInt8](capacity=capacity)
        self.text.resize(unsafe_uninit_length=capacity)
        self.used = 0
        self.letter = 0
        self.length = 0

    @inline(.always)
    def add(mut self, letter: UInt8, length: Int):
        """`length` more of `letter`, joining the run being added to when it has the same letter."""
        if letter != self.letter:
            self.flush()
            self.letter = letter
        self.length += length

    def flush(mut self):
        """Writes the run being added to, if any, as its length's digits then its letter."""
        if self.length == 0:
            return
        var digits = 1
        var power = 10
        while power <= self.length:
            digits += 1
            power *= 10
        var at = self.used
        self.used += digits + 1
        var out = self.text.unsafe_ptr()
        var rest = self.length
        for place in range(digits - 1, -1, -1):
            out[unsafe_offset=at + place] = UInt8(ord("0") + rest % 10)
            rest //= 10
        out[unsafe_offset=at + digits] = self.letter
        self.length = 0

    def finish(var self) -> String:
        """The CIGAR string, its last run written."""
        self.flush()
        self.text.resize(self.used, 0)
        return String(unsafe_from_utf8=self.text)
