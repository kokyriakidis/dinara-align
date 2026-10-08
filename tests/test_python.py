# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""The Python package (python/dinara_align) against known answers, its own docstring's examples, and a
small full-matrix oracle of its own, sharing no code with the library.

    pixi run test-python
"""

import doctest
import random
import sys

import dinara_align as da


def edit_distance(a: str, b: str) -> int:
    """The edit distance of `a` and `b`, by the full matrix one row at a time."""
    row = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        previous, row[0] = row[0], i
        for j, y in enumerate(b, 1):
            previous, row[j] = row[j], min(row[j] + 1, row[j - 1] + 1, previous + (x != y))
    return row[-1]


def rows_cost(top: str, bottom: str, costs: da.Costs) -> int:
    """What two gapped rows cost: substitutions, and each gap run at its side's cheaper piece."""
    total, index = 0, 0
    while index < len(top):
        if top[index] == "-" or bottom[index] == "-":
            deleted = bottom[index] == "-"
            length = 0
            while index < len(top) and (bottom[index] == "-" if deleted else top[index] == "-"):
                length += 1
                index += 1
            opening = costs.deletion_opening if deleted and costs.deletion_extension else costs.opening
            extension = costs.deletion_extension if deleted and costs.deletion_extension else costs.extension
            cost = opening + extension * length
            if costs.opening2 >= 0:
                cost = min(cost, costs.opening2 + costs.extension2 * length)
            total += cost
            continue
        total += costs.mismatch if top[index] != bottom[index] else 0
        index += 1
    return total


def mutated(text: str, rate: float, rng: random.Random) -> str:
    """`text` with about a `rate` of its letters substituted, deleted or followed by an inserted letter,
    a third each."""
    out = []
    for letter in text:
        roll = rng.random()
        if roll < rate / 3:
            out.append(rng.choice("ACGT"))
        elif roll < 2 * rate / 3:
            continue
        elif roll < rate:
            out.append(letter + rng.choice("ACGT"))
        else:
            out.append(letter)
    return "".join(out)


def main() -> None:
    """Runs every check, the package's own examples first; a failed one raises."""
    failures, tried = doctest.testmod(da)
    assert failures == 0, f"{failures} of the package's examples failed"

    assert da.distance("ACGTACGTTTGCA", "ACGTCGTTTTGCA") == 2
    found = da.align("ACGTACGTTTGCA", "ACGTCGTTTTGCA")
    assert (found.cost, found.cigar) == (2, "4=1D2=1I6=")
    assert da.align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", eqx=False).cigar == "4M1D2M1I6M"
    assert da.align("ACGTTTTACG", "ACGTTTACG", da.Costs.affine(4, 6, 2)).cigar == "3=1D6="
    assert da.align("ACGTTTTACG", "ACGTTTACG", da.Costs.affine(4, 6, 2), ties="right").cigar == "6=1D3="
    assert da.align(b"ACGT", b"ACGT").cigar == "4="
    affine = da.Costs.affine(4, 6, 2)
    assert da.distance("ACGTACGTTTGCA", "ACGTCGTTTTGCA", affine, max_cost=11) is None
    assert da.align("ACGTACGTTTGCA", "ACGTCGTTTTGCA", affine, max_cost=12).cigar == "4=3X6="
    assert da.distance("AAAACCCC", "CCCCAAAA", affine, band=da.Band.around(0)) == 32
    assert da.score("ACGTACGTTTGCA", "ACGTCGTTTTGCA", affine) == -12
    two = da.Costs.two_piece(4, 6, 2, 24, 1)
    gapped = "GATTACAGCTTGCA" + "C" * 30 + "TGGACCATGAGTCA"
    assert da.distance(gapped, "GATTACAGCTTGCATGGACCATGAGTCA", two) == 54
    assert da.distance("ACGACGT", "ACGTACGT", affine.with_deletions(20, 5)) == 8
    assert da.distance("ACGTACGT", "ACGACGT", affine.with_deletions(20, 5)) == 25
    # Every mode, and the SAM fields.
    core = da.align("GGGGACGTACGTGGGG", "CCCCACGTTCGTCC", affine, da.Mode.local(2))
    assert core.cigar == "4=1X3=" and core.clipped_cigar(14) == "4S4=1X3=2S"
    assert core.mismatch_string("GGGGACGTACGTGGGG", "CCCCACGTTCGTCC") == "4A3"
    assert core.edit_distance("GGGGACGTACGTGGGG", "CCCCACGTTCGTCC") == 1
    assert da.align("TTTTTACGTACGT", "ACGTACGTGGGGG", affine, da.Mode.overlap(2)).reference_start == 5
    onward = da.align("ACGTTGCAAGGCTTTT", "ACGTTGCAAGGCGAGA", affine, da.Mode.extension(1))
    assert (onward.score, onward.cigar) == (12, "12=")
    assert da.align("ACGTTGCAAGGCTTTT", "ACGTTGCAAGGCGAGA", affine, da.Mode.extension(1, zdrop=5)).score == 12
    scores = da.local_scores("ACGTACGTAC" + "T" * 20 + "ACGTACGTAC", "ACGTACGTAC", affine, da.Mode.local(2), window=5)
    assert (scores.score, scores.second_score, scores.second_reference_end) == (20, 20, 10)
    dna = da.Scoring.dna()
    assert da.score("TTTTACGTACGTTTTT", "ACGTACGT", dna, da.Mode.LOCAL) == 16
    placed = da.align("TTTTACGTACGTTTTT", "ACGTACGT", dna, da.Mode.INFIX)
    assert (placed.score, placed.reference_start, placed.cigar) == (16, 4, "8=")
    blosum_like = da.Scoring.tabulated("ACGT", [5, -1, -2, -1, -1, 5, -3, -2, -2, -3, 5, -1, -1, -2, -1, 5], -6, -1)
    assert da.align("ACGTAC", "ACGTAC", blosum_like).score == 30
    try:
        da.align("A", "C", da.Costs.affine(0, 6, 2))
        raise AssertionError("free mismatches were not refused")
    except ValueError:
        pass
    try:
        da.distance("ACGT", "ACGT", affine, da.Mode.local(2))
        raise AssertionError("a local distance was not refused")
    except ValueError:
        pass

    # Random pairs: the edit distance against a plain dynamic program, alignments priced afresh, and
    # batches equal to their pairs one by one, on two threads and on all.
    rng = random.Random(7)
    references, queries = [], []
    for trial in range(200):
        reference = "".join(rng.choice("ACGT") for _ in range(rng.randrange(0, 120)))
        query = mutated(reference, [0.0, 0.05, 0.2, 0.5][trial % 4], rng)
        references.append(reference)
        queries.append(query)
        assert da.distance(reference, query) == edit_distance(reference, query)
        for costs in (da.Costs.edit(), affine, two, affine.with_deletions(9, 1)):
            found = da.align(reference, query, costs)
            top, bottom = found.gapped(reference, query)
            assert top.replace("-", "") == reference and bottom.replace("-", "") == query
            assert rows_cost(top, bottom, costs) == found.cost == da.distance(reference, query, costs)
    for threads in (2, 0):
        assert da.distances(references, queries, affine, threads=threads) == [
            da.distance(r, q, affine) for r, q in zip(references, queries)
        ]
        assert [a.cigar for a in da.alignments(references, queries, affine, da.Mode.INFIX, threads=threads)] == [
            da.align(r, q, affine, da.Mode.INFIX).cigar for r, q in zip(references, queries)
        ]
    capped = da.distances(references, queries, affine, max_cost=30)
    assert all((c is None) == (da.distance(r, q, affine) > 30) for c, r, q in zip(capped, references, queries))
    hits = da.search(["TTTT", "ACGTACGTACGT", "ACGTTCGTACGT"], "ACGTACGTACGT", affine, da.Mode.local(2), best=2, aligned=True)
    assert [hit.index for hit in hits] == [1, 2] and hits[0].score == 24 and hits[0].alignment.cigar == "12="
    print(f"Python package: every check passed, {tried} examples of its own among them")


if __name__ == "__main__":
    sys.exit(main())
