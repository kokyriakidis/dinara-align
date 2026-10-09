# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
WFA2-lib's regression set through dinara-align, run with `pixi run test-wfa` (see `scripts/test_wfa.py`).

WFA2-lib ships 305 pairs (`tests/wfa.utest.seq`) and its own results for each mode it tests
(`tests/wfa.utest.check`), a score and a CIGAR a pair. Each mode dinara-align offers is held to them:
the edit and indel distances, gap-affine costs at four penalties, three more with a match reward, and
two-piece gap-affine costs. Every score must equal WFA2-lib's, and every CIGAR must consume both
sequences exactly, its `M` a match and its `X` a mismatch, and cost the expected score when priced
afresh. WFA2-lib's own CIGARs are priced the same way, which checks the pricing.

Every gap-affine mode must also give WFA2-lib's CIGAR byte for byte under `Ties.RIGHT`, its own rule
for ties: the indel distance as costs (3, 0, 1), which give it exactly and never a substitution, and the
modes with a match reward by the costs WFA2-lib folds the reward into. Those run with no memory limit,
so no pair is split (see `gap_affine.solve`); and the edit distance's `align` with `Ties.RIGHT`,
whichever of its searches finds the distance. Its approximate heuristic modes, which dinara-align does
not offer, are left out.
"""

from std.sys import argv

from dinara_align import Band, Costs, Mode, Placement, Scoring, Ties, align, distance
from dinara_align.gap_affine import (
    EndsFree,
    FREE_START,
    Penalties,
    penalties_of,
    cigar_of,
    solve,
)

comptime EDIT = 0
"""Unit costs, the edit distance, a score WFA2-lib reports as the distance itself."""
comptime GAP_AFFINE = 1
"""A cost: no match reward, so `Costs.affine` serves it."""
comptime REWARDED = 2
"""A match earns a reward, which a `Scoring` takes as a score."""
comptime TWO_PIECE = 3
"""A gap costs the less of two affine costs, as `Costs.two_piece` counts them."""


@fieldwise_init
struct WfaMode(Copyable, Movable):
    """One of WFA2-lib's tested modes: its name and penalties."""

    var name: String
    """WFA2-lib's name for the mode, which names its results file."""
    var kind: Int
    """`EDIT`, `GAP_AFFINE`, `REWARDED` or `TWO_PIECE`."""
    var reward: Int
    """What a match earns, for `REWARDED` alone."""
    var mismatch: Int
    """A substitution's cost."""
    var opening: Int
    """A gap's cost before its letters."""
    var extension: Int
    """Each gapped letter's cost."""
    var opening2: Int
    """The second piece's costs, for `TWO_PIECE` alone."""
    var extension2: Int


def lines_of(path: String) raises -> List[String]:
    """A file's lines, empty ones dropped."""
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


def wfa_cigar[pieces: Int](first: String, second: String, penalties: Penalties) -> String:
    """The alignment `Ties.RIGHT` picks, in WFA2-lib's letters, found with no memory limit."""
    var moves = List[UInt8]()
    var cost = solve[pieces](
        first.as_bytes(),
        second.as_bytes(),
        penalties,
        FREE_START,
        FREE_START,
        Int.MAX,
        moves,
        True,
        Int.MAX,
        Band(),
        Ties.RIGHT,
    )
    return cigar_of(first, second, moves^, cost, penalties, True).replace("=", "M")


def priced(cigar: String, first: String, second: String, mode: WfaMode) -> Optional[Int]:
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
            var gap = mode.opening + mode.extension * length
            if mode.kind == TWO_PIECE:
                gap = min(gap, mode.opening2 + mode.extension2 * length)
            cost += gap
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
    if mode.kind == GAP_AFFINE or mode.kind == TWO_PIECE:
        return -cost
    return mode.reward * matches - cost


def main() raises:
    """Holds every mode to WFA2-lib's results, its tests directory the first argument; raises on any
    disagreement."""
    var root = String(argv()[1])
    var seq = lines_of(root + "/wfa.utest.seq")
    # A pair is two lines, the first sequence after a `>` and the second after a `<`.
    var firsts = List[String]()
    var seconds = List[String]()
    for index in range(0, len(seq), 2):
        firsts.append(String(seq[index][byte=1:]))
        seconds.append(String(seq[index + 1][byte=1:]))
    # WFA2-lib's names and penalties, `match,mismatch,opening,extension` with a reward as a negative
    # match, as its `wfa.utest.sh` runs them; a gap of `k` costs `opening + k extension`.
    var modes: List[WfaMode] = [
        WfaMode("edit", EDIT, 0, 1, 0, 1, 0, 0),
        WfaMode("indel", GAP_AFFINE, 0, 2, 0, 1, 0, 0),
        WfaMode("affine", GAP_AFFINE, 0, 4, 6, 2, 0, 0),
        WfaMode("affine.p0", GAP_AFFINE, 0, 1, 2, 1, 0, 0),
        WfaMode("affine.p1", GAP_AFFINE, 0, 3, 1, 4, 0, 0),
        WfaMode("affine.p2", GAP_AFFINE, 0, 5, 3, 2, 0, 0),
        WfaMode("affine.p3", REWARDED, 5, 1, 2, 1, 0, 0),
        WfaMode("affine.p4", REWARDED, 2, 3, 1, 4, 0, 0),
        WfaMode("affine.p5", REWARDED, 3, 5, 3, 2, 0, 0),
        # Its align_benchmark's default two-piece penalties, `0,4,6,2,24,1`.
        WfaMode("affine2p", TWO_PIECE, 0, 4, 6, 2, 24, 1),
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
                var aligned = align(first, second)
                var exact = align(first, second, ties=Ties.RIGHT).cigar.replace("=", "M")
                if exact != theirs:
                    print("   ", mode.name, "pair", index, ": Ties.RIGHT gives", exact, "where WFA2-lib gives", theirs)
                    wrong += 1
                score = aligned.cost
                cigar = aligned.cigar
            elif mode.kind == GAP_AFFINE:
                var costs = Costs.affine(mode.mismatch, mode.opening, mode.extension)
                var aligned = align(first, second, costs)
                score = aligned.cost if mode.name == "indel" else -aligned.cost
                cigar = aligned.cigar
                if distance(first, second, costs) != aligned.cost:
                    print("   ", mode.name, "pair", index, ": distance and align disagree")
                    wrong += 1
            elif mode.kind == TWO_PIECE:
                var costs = Costs.two_piece(mode.mismatch, mode.opening, mode.extension, mode.opening2, mode.extension2)
                var aligned = align(first, second, costs)
                score = -aligned.cost
                cigar = aligned.cigar
                if distance(first, second, costs) != aligned.cost:
                    print("   ", mode.name, "pair", index, ": distance and align disagree")
                    wrong += 1
            else:
                var scoring = Scoring.uniform(mode.reward, -mode.mismatch, -mode.opening, -mode.extension)
                var aligned = align(first, second, scoring, Mode.GLOBAL, placement=Placement.on_cpu(1))
                score = Int(aligned.score)
                var rows = aligned.gapped(first, second)
                cigar = rows_cigar(rows[0], rows[1])
            var earned = priced(cigar, first, second, mode)
            if score != want or not earned or earned.value() != want:
                print("   ", mode.name, "pair", index, ": score", score, "want", want, "CIGAR", cigar)
                wrong += 1
            if cigar.replace("=", "M") != theirs:
                differing += 1
            if mode.kind != EDIT:
                var exact: String
                if mode.kind == TWO_PIECE:
                    exact = wfa_cigar[2](
                        first,
                        second,
                        penalties_of(
                            Costs.two_piece(mode.mismatch, mode.opening, mode.extension, mode.opening2, mode.extension2)
                        ),
                    )
                elif mode.kind == REWARDED:
                    # WFA2-lib's own fold of the reward into costs (see `gap_affine`).
                    var a = mode.reward
                    exact = wfa_cigar[1](
                        first,
                        second,
                        penalties_of(Costs.affine(2 * (mode.mismatch + a), 2 * mode.opening, 2 * mode.extension + a)),
                    )
                elif mode.name == "indel":
                    exact = wfa_cigar[1](first, second, penalties_of(Costs.affine(3, 0, 1)))
                else:
                    exact = wfa_cigar[1](
                        first, second, penalties_of(Costs.affine(mode.mismatch, mode.opening, mode.extension))
                    )
                if exact != theirs:
                    print("   ", mode.name, "pair", index, ": Ties.RIGHT gives", exact, "where WFA2-lib gives", theirs)
                    wrong += 1
        print(
            mode.name,
            ":",
            len(firsts) - wrong,
            "of",
            len(firsts),
            "pairs agree;",
            differing,
            "default CIGARs differ from WFA2-lib's at the same score",
            "(gaps placed left), none under Ties.RIGHT",
        )
        failures += wrong
    if failures:
        raise Error(String(failures, " disagreements with WFA2-lib"))
    print("WFA2-lib's regression set: every score agrees, every CIGAR earns it")
