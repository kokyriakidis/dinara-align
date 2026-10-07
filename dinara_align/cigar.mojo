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
    """Runs written back as a CIGAR, the empty ones left out."""
    var out = String()
    for index in range(len(letters)):
        if lengths[index] > 0:
            out += String(lengths[index], chr(Int(letters[index])))
    return out


def reversed_cigar(cigar: String) -> String:
    """A CIGAR's runs in the other order: the alignment of both sequences reversed."""
    var runs = cigar_runs(cigar)
    var letters = List[UInt8](capacity=len(runs[0]))
    var lengths = List[Int](capacity=len(runs[1]))
    for index in range(len(runs[0]) - 1, -1, -1):
        letters.append(runs[0][index])
        lengths.append(runs[1][index])
    return joined_cigar(letters, lengths)


def reversed_text(text: String) -> String:
    """A sequence back to front."""
    var bytes = text.as_bytes()
    var out = List[UInt8](capacity=len(bytes))
    for index in range(len(bytes) - 1, -1, -1):
        out.append(bytes[index])
    return String(unsafe_from_utf8=out^)


def cigar_matches(first: String, second: String, cigar: String) -> Int:
    """The equal pairs a CIGAR of `first` against `second` aligns, `M` runs compared letter by letter."""
    var a = first.as_bytes()
    var b = second.as_bytes()
    var runs = cigar_runs(cigar)
    var column = 0
    var row = 0
    var total = 0
    for index in range(len(runs[0])):
        var letter = runs[0][index]
        var length = runs[1][index]
        if letter == UInt8(ord("D")):
            column += length
        elif letter == UInt8(ord("I")):
            row += length
        else:
            for _ in range(length):
                if a[column] == b[row]:
                    total += 1
                column += 1
                row += 1
    return total


def cigar_cost(
    first: String,
    second: String,
    cigar: String,
    mismatch: Int,
    opening: Int,
    extension: Int,
    opening2: Int = -1,
    extension2: Int = 0,
) -> Int:
    """What a CIGAR of `first` against `second` costs: its substitutions, `M` runs compared letter by
    letter, and each gap run at the cheaper piece, the second only with `opening2` not negative."""
    var runs = cigar_runs(cigar)
    var mismatches = 0
    var gaps = 0
    for index in range(len(runs[0])):
        var letter = runs[0][index]
        var length = runs[1][index]
        if letter == UInt8(ord("D")) or letter == UInt8(ord("I")):
            var gap = opening + extension * length
            if opening2 >= 0:
                gap = min(gap, opening2 + extension2 * length)
            gaps += gap
    var aligned = 0
    for index in range(len(runs[0])):
        var letter = runs[0][index]
        if letter != UInt8(ord("D")) and letter != UInt8(ord("I")):
            aligned += runs[1][index]
    mismatches = aligned - cigar_matches(first, second, cigar)
    return mismatches * mismatch + gaps
