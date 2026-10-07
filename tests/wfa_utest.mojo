"""
WFA2-lib's regression set through dinara-align, run with `pixi run test-wfa` (see `scripts/test_wfa.py`).

WFA2-lib ships 305 pairs (`tests/wfa.utest.seq`) and its own results for each mode it tests
(`tests/wfa.utest.check`), a score and a CIGAR a pair. Each mode dinara-align offers is held to them:
the edit and indel distances, gap-affine costs at four penalties, and three more with a match reward.
Every score must equal WFA2-lib's, and every CIGAR must consume both sequences exactly, its `M` a match
and its `X` a mismatch, and cost the expected score when priced afresh, so a tie resolved another way
passes and a wrong alignment does not. WFA2-lib's own CIGARs are priced the same way, which checks the
pricing. Its two-piece affine and approximate heuristic modes, which dinara-align does not offer, are
left out.

The indel distance runs as gap-affine costs (2, 0, 1), where a substitution costs what a deletion and
an insertion do; its CIGARs may write one as `X`, priced as that pair.
"""

from std.sys import argv

from dinara_align import AlignmentMode, Placement, Scoring, affine_cigar, affine_distance, align, edit_cigar

comptime EDIT = 0
comptime GAP_AFFINE = 1
"""A cost: no match reward, so `affine_cigar` serves it."""
comptime REWARDED = 2
"""A match earns a reward, which `align` takes as a score."""


@fieldwise_init
struct Mode(Copyable, Movable):
    var name: String
    var kind: Int
    var reward: Int
    var mismatch: Int
    var opening: Int
    var extension: Int


def lines_of(path: String) raises -> List[String]:
    var out = List[String]()
    for line in open(path, "r").read().split("\n"):
        if line.byte_length() > 0:
            out.append(String(line))
    return out^


def rows_cigar(top: String, bottom: String) -> String:
    """A CIGAR from two gapped rows in WFA2-lib's letters: `M` a match, `X` a mismatch, `D` the top's
    letter alone, `I` the bottom's."""
    var a = top.as_bytes()
    var b = bottom.as_bytes()
    var out = String()
    var last = UInt8(0)
    var run = 0
    for index in range(len(a)):
        var op = UInt8(ord("I")) if a[index] == UInt8(ord("-")) else (
            UInt8(ord("D")) if b[index]
            == UInt8(ord("-")) else (UInt8(ord("M")) if a[index] == b[index] else UInt8(ord("X")))
        )
        if op != last and run > 0:
            out += String(run, chr(Int(last)))
            run = 0
        last = op
        run += 1
    if run > 0:
        out += String(run, chr(Int(last)))
    return out


def priced(cigar: String, first: String, second: String, mode: Mode) -> Optional[Int]:
    """The score a CIGAR earns under `mode` in WFA2-lib's form (a distance, a cost negated, or a reward
    less the costs), or None when it does not spell an alignment of the two sequences."""
    var a = first.as_bytes()
    var b = second.as_bytes()
    var i = 0
    var j = 0
    var matches = 0
    var cost = 0
    var length = 0
    for byte in cigar.as_bytes():
        if byte >= UInt8(ord("0")) and byte <= UInt8(ord("9")):
            length = length * 10 + Int(byte - UInt8(ord("0")))
            continue
        if byte == UInt8(ord("M")) or byte == UInt8(ord("X")) or byte == UInt8(ord("=")):
            for _ in range(length):
                if i >= len(a) or j >= len(b):
                    return None
                var equal = a[i] == b[j]
                if equal != (byte != UInt8(ord("X"))):
                    return None
                if equal:
                    matches += 1
                else:
                    cost += mode.mismatch
                i += 1
                j += 1
        elif byte == UInt8(ord("D")) or byte == UInt8(ord("I")):
            cost += mode.opening + mode.extension * length
            if byte == UInt8(ord("D")):
                i += length
            else:
                j += length
        else:
            return None
        length = 0
    if i != len(a) or j != len(b):
        return None
    if mode.kind == EDIT or mode.name == "indel":
        return cost
    if mode.kind == GAP_AFFINE:
        return -cost
    return mode.reward * matches - cost


def main() raises:
    var root = String(argv()[1])
    var seq = lines_of(root + "/wfa.utest.seq")
    var firsts = List[String]()
    var seconds = List[String]()
    for index in range(0, len(seq), 2):
        firsts.append(String(seq[index][byte=1:]))
        seconds.append(String(seq[index + 1][byte=1:]))
    # WFA2-lib's names and penalties, `match,mismatch,opening,extension` with a reward as a negative
    # match, as its `wfa.utest.sh` runs them; a gap of `k` costs `opening + k extension`.
    var modes: List[Mode] = [
        Mode("edit", EDIT, 0, 1, 0, 1),
        Mode("indel", GAP_AFFINE, 0, 2, 0, 1),
        Mode("affine", GAP_AFFINE, 0, 4, 6, 2),
        Mode("affine.p0", GAP_AFFINE, 0, 1, 2, 1),
        Mode("affine.p1", GAP_AFFINE, 0, 3, 1, 4),
        Mode("affine.p2", GAP_AFFINE, 0, 5, 3, 2),
        Mode("affine.p3", REWARDED, 5, 1, 2, 1),
        Mode("affine.p4", REWARDED, 2, 3, 1, 4),
        Mode("affine.p5", REWARDED, 3, 5, 3, 2),
    ]
    var failures = 0
    for mode in modes:
        var expected = lines_of(root + "/wfa.utest.check/test." + mode.name + ".alg")
        var wrong = 0
        var differing = 0
        for index in range(len(firsts)):
            var first = firsts[index]
            var second = seconds[index]
            var fields = expected[index].split("\t")
            var want = Int(String(fields[0]))
            var theirs = String(fields[1])
            if not priced(theirs, first, second, mode) or priced(theirs, first, second, mode).value() != want:
                print("   ", mode.name, "pair", index, ": WFA2-lib's own CIGAR does not price to its score")
                wrong += 1
                continue
            var score: Int
            var cigar: String
            if mode.kind == EDIT:
                var aligned = edit_cigar(first, second)
                score = aligned.distance
                cigar = aligned.cigar
            elif mode.kind == GAP_AFFINE:
                var aligned = affine_cigar(first, second, mode.mismatch, mode.opening, mode.extension)
                score = aligned.cost if mode.name == "indel" else -aligned.cost
                cigar = aligned.cigar
                if affine_distance(first, second, mode.mismatch, mode.opening, mode.extension) != aligned.cost:
                    print("   ", mode.name, "pair", index, ": affine_distance and affine_cigar disagree")
                    wrong += 1
            else:
                var scoring = Scoring.uniform(
                    mode.reward, -mode.mismatch, -(mode.opening + mode.extension), -mode.extension
                )
                var aligned = align[AlignmentMode.GLOBAL](first, second, scoring, Placement.on_cpu(1))
                score = Int(aligned.score)
                cigar = rows_cigar(aligned.first_gapped, aligned.second_gapped)
            var earned = priced(cigar, first, second, mode)
            if score != want or not earned or earned.value() != want:
                print("   ", mode.name, "pair", index, ": score", score, "want", want, "CIGAR", cigar)
                wrong += 1
            if cigar.replace("=", "M") != theirs:
                differing += 1
        print(
            mode.name,
            ":",
            len(firsts) - wrong,
            "of",
            len(firsts),
            "pairs agree;",
            differing,
            "CIGARs differ from WFA2-lib's at the same score",
        )
        failures += wrong
    if failures:
        raise Error(String(failures, " disagreements with WFA2-lib"))
    print("WFA2-lib's regression set: every score agrees, every CIGAR earns it")
