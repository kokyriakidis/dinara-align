#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Regenerates the ported kernels from an AffineGaps checkout, so a later upstream fix can be pulled in.

    git clone https://github.com/unum-science/AffineGaps /tmp/AffineGaps
    pixi run port /tmp/AffineGaps
    git diff  # review, then `pixi run test`

Copies `alignment.mojo`, `common.mojo` and `errors.mojo` into the `dinara_align` package, drops
everything that serves folding, cofolding, proteins, the Python bindings or the command line, and
applies the Mojo 1.1 migrations. Every other module, `api.mojo`, `modes.mojo`, `scoring.mojo` and
`__init__.mojo` among them, is this repository's own and is left alone.

Every textual edit is asserted, so an upstream rewrite fails here loudly rather than shipping half a
port. The script is the record of how this package differs from upstream.
"""

import os
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
PACKAGE = ROOT / "dinara_align"


HEADER = (
    "# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the\n"
    "# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.\n"
    "#\n"
    "# Derived from AffineGaps (https://github.com/unum-science/AffineGaps), Copyright Ash Vardanian, under the\n"
    "# Apache License, Version 2.0: see LICENSES/Apache-2.0.txt and NOTICE.\n"
)
"""The license header every generated module starts with: MPL-2.0, over code derived from Apache-2.0."""

def cut(text: str, start: str, end: str) -> str:
    """Removes everything from `start` up to, but not including, `end`."""
    begin = text.index(start)
    return text[:begin] + text[text.index(end, begin) :]


def replace(text: str, old: str, new: str) -> str:
    """A replacement that must happen, so an upstream rewrite fails here rather than passing silently."""
    assert old in text, f"upstream no longer contains: {old[:80]!r}"
    return text.replace(old, new, 1)


def migrate(text: str) -> str:
    """The changes every ported module needs: package-relative imports, the 1.1 renames, the error type."""
    text = re.sub(r"^from (errors|common|alignment) import", r"from .\1 import", text, flags=re.M)
    text = re.sub(r"\bInlineArray\b", "Array", text)
    return text.replace("AffineGapsError", "AlignmentError")


GAPPED_ROWS_END = '''    """The second sequence, gapped to the same columns."""\n'''
GAPPED_ROWS_CIGAR = '''
    def cigar(self, eqx: Bool = True) -> String:
        """The rows as a CIGAR, the first sequence the reference, as `Alignment`'s reads: `=` a match and
        `X` a substitution, or with `eqx` false `M` for both, `D` a letter of the first alone and
        `I` one of the second, each run its length then its letter."""
        comptime GAP = UInt8(ord("-"))
        var top = self.first_gapped.as_bytes()
        var bottom = self.second_gapped.as_bytes()
        # Bytes written straight, each run's digits then its letter: no string per run.
        var out = List[UInt8](capacity=64)

        def emit(mut out: List[UInt8], run: Int, letter: UInt8):
            """Appends one run to `out`: the decimal digits of `run`, most significant first, then `letter`."""
            var digits = Array[UInt8, 20](fill=0)
            var count = 0
            var value = run
            while value > 0:
                digits[count] = UInt8(ord("0")) + UInt8(value % 10)
                value //= 10
                count += 1
            for index in range(count - 1, -1, -1):
                out.append(digits[index])
            out.append(letter)

        var last = UInt8(0)
        var run = 0
        for column in range(len(top)):
            var letter: UInt8
            if top[column] == GAP:
                letter = UInt8(ord("I"))
            elif bottom[column] == GAP:
                letter = UInt8(ord("D"))
            elif not eqx:
                letter = UInt8(ord("M"))
            else:
                letter = UInt8(ord("=")) if top[column] == bottom[column] else UInt8(ord("X"))
            if letter != last and run > 0:
                emit(out, run, last)
                run = 0
            last = letter
            run += 1
        if run > 0:
            emit(out, run, last)
        return String(unsafe_from_utf8=out^)
'''
"""The gapped rows' `cigar`, which upstream has not."""


def port_alignment(text: str) -> str:
    """The kernels, with the 1.1 import moves and the references to a Python oracle that is not here."""
    text = replace(
        text,
        "transcribed from `affinegaps.py`, which stays\nthe parity oracle.",
        "transcribed from the NumPy reference in\nAffineGaps, which stays the parity oracle upstream.",
    )
    text = replace(
        text,
        "the recurrence that `affinegaps.py` holds as the oracle;",
        "the recurrence that AffineGaps' NumPy reference holds as the oracle;",
    )
    text = replace(text, "functions of `affinegaps.py`.", "functions of AffineGaps' NumPy reference.")
    text = text.replace("The Python walks a single layer", "AffineGaps' NumPy reference walks a single layer")
    text = text.replace("walking exactly as\n    the Python does.", "walking exactly as\n    that reference does.")
    text = text.replace("the Python's scan", "the reference's scan").replace("the Python scan's", "the reference scan's")
    text = text.replace("# Python scan's strict", "# reference scan's strict")
    text = replace(text, "dim for a gap, mirroring `colorize_alignment`.", "dim for a gap, as ANSI escapes.")
    text = replace(
        text,
        "comptime ALL_MODES = (AlignmentMode.GLOBAL, AlignmentMode.LOCAL)\n",
        "comptime ALL_MODES: Array[AlignmentMode, 2] = [AlignmentMode.GLOBAL, AlignmentMode.LOCAL]\n"
        '"""An array rather than a tuple: since Mojo 1.1 an imported tuple indexed in a `comptime for` yields the tuple."""\n',
    )
    text = replace(text, "from std.gpu import block_dim, block_idx, grid_dim, lane_id, thread_idx\n", "")
    text = replace(text, "from std.gpu.primitives.warp import WARP_SIZE, shuffle_down, shuffle_up, shuffle_xor\n", "")
    text = replace(
        text,
        "from max.gpu import barrier\n",
        "from max.gpu import WARP_SIZE, barrier, block_dim, block_idx, grid_dim, lane_id, thread_idx\n"
        "from max.gpu.primitives.warp import shuffle_down, shuffle_up, shuffle_xor\n",
    )
    text = replace(text, "\n            @parameter\n            def solve_leaf(slot: Int):", "\n            def solve_leaf(slot: Int) {imm}:")
    text = replace(
        text,
        "parallelize[solve_leaf](len(leaves), placement.threads)",
        "parallelize(solve_leaf, len(leaves), placement.threads)",
    )
    # DNA only: the protein defaults and the BLOSUM62 table go, and `scoring.mojo` supplies DNA presets.
    text = replace(text, "    DEFAULT_PROTEINS_ALPHABET,\n", "")
    text = replace(text, "    translate,\n", "")
    # `Alignment` is this package's CIGAR result (see `modes.mojo`); the gapped rows are named as such.
    assert "AlignmentResult" in text, "upstream no longer names AlignmentResult"
    text = re.sub(r"\bAlignmentResult\b", "GappedAlignment", text)
    # The gapped rows read as a CIGAR too, as every `Alignment` of this package is.
    text = replace(text, GAPPED_ROWS_END, GAPPED_ROWS_END + GAPPED_ROWS_CIGAR)
    text = cut(text, "comptime DEFAULT_PROTEINS_SCALE", "comptime CORNER_BYTES")
    return cut(text, "# BLOSUM62 scaled by five", "@fieldwise_init\nstruct AffineGapCosts")


def port_common(text: str) -> str:
    """Shared primitives, without the RNA alphabet and the dot-bracket bytes."""
    text = replace(text, text[: text.index('"""', 3) + 3], '''"""
Primitives the alignment kernels share with the routing layer.

Symbol codes, the device staging helpers, and the scoring records every entry point reads. Anything
that presumes a rotating band or an affine gap belongs to `alignment.mojo`.
"""''')
    text = cut(text, 'comptime DEFAULT_RNA_ALPHABET = "ACGU"', "comptime THREADS_PER_BLOCK")
    text = cut(text, "comptime PositionDType", "comptime MAX_ALPHABET_SIZE")
    text = cut(text, "comptime DEFAULT_PROTEINS_ALPHABET", "comptime THREADS_PER_BLOCK")
    text = replace(
        text,
        "Thirty-two covers the twenty-three protein letters with room to\nspare, and costs one kilobyte per block.",
        "Thirty-two holds all fifteen IUPAC nucleotide codes in both\ncases, and costs one kilobyte per block.",
    )
    text = replace(text, "from std.gpu.primitives.warp import WARP_SIZE\n", "")
    text = replace(text, "from max.gpu.host import DeviceAttribute", "from max.gpu import WARP_SIZE\nfrom max.gpu.host import DeviceAttribute")
    text = replace(text, "    @parameter\n    if CompilationTarget.is_linux():", "    comptime if CompilationTarget.is_linux():")
    text = replace(text, "            mask[index] = 0", "            mask[unsafe_offset=index] = 0")
    text = replace(text, "mask[index].reduce_bit_count()", "mask[unsafe_offset=index].reduce_bit_count()")
    text = replace(
        text,
        "from std.sys.info import CompilationTarget, num_logical_cores, size_of",
        "from std.sys.info import CompilationTarget, has_apple_gpu_accelerator, num_logical_cores, size_of",
    )
    return replace(
        text,
        '    Each is a live driver call, so they are asked together and once.\n    """\n',
        '''    Each is a live driver call, so they are asked together and once.
    """
    comptime if has_apple_gpu_accelerator():
        # Metal answers neither the per-multiprocessor budget, the opt-in ceiling nor the resident-block
        # count. A threadgroup gets one fixed allotment, so that is the whole budget with nothing reserved
        # out of it, and the thread ceiling over the one-warp blocks the sweeps launch stands in for the
        # resident count, which only steers how finely a level splits.
        return GpuSpecs(
            Int(context.get_attribute(DeviceAttribute.MAX_SHARED_MEMORY_PER_BLOCK)),
            0,
            Int(context.max_single_alloc_size()),
            Int(context.get_attribute(DeviceAttribute.MULTIPROCESSOR_COUNT)),
            Int(context.get_attribute(DeviceAttribute.MAX_THREADS_PER_BLOCK)) // WARP_SIZE,
        )
''',
    )


def port_errors(text: str) -> str:
    """The error type, without the folding table's kind, the unit-cost aligner's ASCII refusal (the
    wavefront takes any byte) and the command line's exit status, and with a band no alignment fits."""
    text = cut(text, "    comptime INCONSISTENT_TABLE", "    def write_to")
    text = cut(text, "        elif self == Self.INCONSISTENT_TABLE:", "        else:")
    text = cut(text, "    comptime NOT_ASCII", "    def write_to")
    text = cut(text, "        elif self == Self.NOT_ASCII:", "        else:")
    text = replace(
        text,
        "    def write_to",
        '    comptime OUTSIDE_BAND = Self(-9)\n    """No alignment stays inside the band of diagonals asked for."""\n\n    def write_to',
    )
    text = replace(
        text,
        "        else:\n            writer.write(\"an unnamed failure\")",
        '        elif self == Self.OUTSIDE_BAND:\n            writer.write("no alignment stays inside the band")\n'
        '        else:\n            writer.write("an unnamed failure")',
    )
    text = cut(text, "    def exit_status(self) -> Int:", "\n\n@fieldwise_init")
    return replace(text, 'writer.write("AffineGaps: ", self.kind', 'writer.write("dinara-align: ", self.kind')


def unused_literals() -> dict[str, list[int]]:
    """Asks the compiler which string literals stand alone as statements, by file and line.

    Telling those apart from field and function docstrings is what the compiler already does, so
    its warnings are the list rather than a second parser here. The test suite is compiled because
    it instantiates every entry point; the accelerator is pinned only so no vendor toolchain is needed.
    """
    outcome = subprocess.run(
        ["mojo", "build", "-I", ".", "tests/test_alignment.mojo", "-o", os.devnull, "--target-accelerator", "sm_90a"],
        cwd=ROOT,
        capture_output=True,
        text=True,
    )
    if re.search(r"\berror:", outcome.stderr):
        # Without a clean compile the warnings are missing too, and the hoist would silently do nothing.
        sys.exit(f"the test suite does not compile, so the literals cannot be located:\n{outcome.stderr}")
    found: dict[str, list[int]] = {}
    for match in re.finditer(r"([\w./]+\.mojo):(\d+):\d+: warning: 'StringLiteral", outcome.stderr):
        found.setdefault(match.group(1), []).append(int(match.group(2)))
    return found


def hoist(path: pathlib.Path, starts: list[int]) -> None:
    """Turns each listed literal into a comment above the statement it documents.

    Upstream writes these after the statement, the way an attribute docstring is written, and Mojo
    1.1 warns that the value is unused. As a comment it has to come first, so the statement is found
    by walking back to the literal's own indentation, past continuation lines and closing brackets.
    """
    lines = path.read_text().split("\n")
    for start in sorted(set(starts), reverse=True):
        index = start - 1
        stripped = lines[index].strip()
        indent = lines[index][: len(lines[index]) - len(lines[index].lstrip())]
        if stripped.count('"""') >= 2:
            end, body = index, [stripped[3 : stripped.rindex('"""')]]
        else:
            end = next(j for j in range(index + 1, len(lines)) if '"""' in lines[j])
            body = [stripped[3:]] + [lines[j].strip() for j in range(index + 1, end)]
            body.append(lines[end].strip()[: lines[end].strip().rindex('"""')])
        comment = [indent + ("# " + piece.strip() if piece.strip() else "#") for piece in body]
        while comment and comment[0].strip() == "#":
            comment.pop(0)
        while comment and comment[-1].strip() == "#":
            comment.pop()
        del lines[index : end + 1]
        statement = index - 1
        while statement > 0:
            candidate = lines[statement]
            if (
                candidate.strip()
                and not candidate.strip().startswith(")")
                and len(candidate) - len(candidate.lstrip()) == len(indent)
            ):
                break
            statement -= 1
        lines[statement:statement] = comment
    path.write_text("\n".join(lines))


def main(upstream: pathlib.Path) -> None:
    """Ports the three modules from the affine-gaps checkout at `upstream` into the package, then formats them."""
    for name, port in (("alignment", port_alignment), ("common", port_common), ("errors", port_errors)):
        source = (upstream / f"{name}.mojo").read_text()
        (PACKAGE / f"{name}.mojo").write_text(HEADER + migrate(port(source)))
    for name, starts in unused_literals().items():
        hoist(PACKAGE / pathlib.Path(name).name, starts)
    # Last, so a regenerated file differs from the committed one only where upstream changed.
    ported = [str(PACKAGE / f"{name}.mojo") for name in ("alignment", "common", "errors")]
    subprocess.run(["mojo", "format", "-q", "-l", "120", *ported], cwd=ROOT, check=True)


if __name__ == "__main__":
    main(pathlib.Path(sys.argv[1]))
