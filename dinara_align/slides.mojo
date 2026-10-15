# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
#
# The slide moved here from diagonal.mojo, ported from `pa-bitpacking` in A*PA
# (https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner, commit bf2e14e), by Ragnar Groot
# Koerkamp and Pesho Ivanov; the gathered slides are this repository's.
"""
How far matches carry a cell along its diagonal: eight letters at a time, by one scalar compare, or on
AVX-512 for a lane group of diagonals at once by two gathers (see `gathered_slides`). Both sequences
end in sentinels that differ from every letter and from each other, so a slide stops at the edge.
The diagonal transition and the affine wavefront both slide their fronts this way.
"""

from std.bit import count_trailing_zeros
from std.sys import llvm_intrinsic
from std.sys.info import CompilationTarget

from .ablation import ABLATE_GATHER

comptime LANES = 8
"""Diagonals a gathered slide takes at once: one AVX-512 vector of eight 64-bit words."""


@inline(.always)
def slide(first: ImmPointer[UInt8, _], second: ImmPointer[UInt8, _], start: Int, diagonal: Int) -> Int:
    """How far matches carry column `start` of `diagonal`, eight letters at a time, to the sentinels."""
    var column = start
    var lag = second.unsafe_offset(-diagonal)
    var mismatches = (
        first.unsafe_offset(column).unsafe_bitcast[UInt64]().unsafe_load()
        ^ lag.unsafe_offset(column).unsafe_bitcast[UInt64]().unsafe_load()
    )
    while mismatches == 0:
        column += 8
        mismatches = (
            first.unsafe_offset(column).unsafe_bitcast[UInt64]().unsafe_load()
            ^ lag.unsafe_offset(column).unsafe_bitcast[UInt64]().unsafe_load()
        )
    return column + (Int(count_trailing_zeros(mismatches)) >> 3)


comptime GATHERED_SLIDES = CompilationTarget.has_avx512f() and not ABLATE_GATHER
"""Whether a lane group's slides start with one gather of each sequence's next eight letters per lane,
AVX-512's `vpgatherqq` at byte offsets, rather than a scalar compare per lane."""


@inline(.always)
def gathered_words(
    base: ImmPointer[UInt8, _],
    offsets: SIMD[DType.int32, LANES],
    lanes: SIMD[DType.bool, LANES],
    elsewhere: SIMD[DType.uint64, LANES],
) -> SIMD[DType.uint64, LANES]:
    """The eight bytes from each of `offsets` past `base` in the `lanes` chosen, `elsewhere` in the rest,
    one AVX-512 gather by 32-bit offsets, the lane group's columns as they are, with no widening on
    the shuffle port. A lane left out reads nothing, so its offset may point anywhere."""
    return llvm_intrinsic["llvm.x86.avx512.mask.gather.dpq.512", SIMD[DType.uint64, LANES], has_side_effect=False](
        elsewhere, base, offsets, lanes, Int32(1)
    )


@inline(.always)
def gathered_slides(
    first: ImmPointer[UInt8, _],
    second: ImmPointer[UInt8, _],
    entries: SIMD[DType.int32, LANES],
    diagonals: SIMD[DType.int32, LANES],
) -> SIMD[DType.int32, LANES]:
    """`slide` of every reached lane of a group at once: both sequences' next eight letters gathered per
    lane and compared. A lane whose eight all match, rare off the path, finishes by `slide`.

    The gathers take the reached lanes as their mask, worked out from this group alone: an all-lanes
    mask, which LLVM builds with a `kxnor` that waits on the mask register's last writer, chained
    every group's gathers to the previous group's, about 21 cycles of latency a group on Skylake.
    An unreached lane gathers nothing, its two words zero and all ones, which never read as a match.
    """
    comptime Wide = SIMD[DType.uint64, LANES]
    var reached = entries.ge(0)
    var mismatches = gathered_words(first, entries, reached, Wide(0)) ^ gathered_words(
        second, entries - diagonals, reached, ~Wide(0)
    )
    var slid = entries + (count_trailing_zeros(mismatches) >> 3).cast[DType.int32]()
    var whole = mismatches.eq(0)
    if whole.reduce_or():
        comptime for lane in range(LANES):
            if whole[lane]:
                slid[lane] = Int32(slide(first, second, Int(slid[lane]), Int(diagonals[lane])))
    return reached.select(slid, entries)


@inline(.always)
def shifted_right(words: SIMD[DType.uint64, LANES], bits: SIMD[DType.uint64, LANES]) -> SIMD[DType.uint64, LANES]:
    """Each lane shifted right by its own count, zero from 64 on, as AVX-512's `vpsrlvq` defines it."""
    return llvm_intrinsic["llvm.x86.avx512.psrlv.q.512", SIMD[DType.uint64, LANES], has_side_effect=False](words, bits)


@inline(.always)
def shifted_left(words: SIMD[DType.uint64, LANES], bits: SIMD[DType.uint64, LANES]) -> SIMD[DType.uint64, LANES]:
    """Each lane shifted left by its own count, zero from 64 on, as AVX-512's `vpsllvq` defines it."""
    return llvm_intrinsic["llvm.x86.avx512.psllv.q.512", SIMD[DType.uint64, LANES], has_side_effect=False](words, bits)


@inline(.always)
def windowed_words(
    base: ImmPointer[UInt8, _],
    offsets: SIMD[DType.int32, LANES],
    lanes: SIMD[DType.bool, LANES],
    elsewhere: SIMD[DType.uint64, LANES],
) -> SIMD[DType.uint64, LANES]:
    """`gathered_words`, from two loads rather than a gather when every chosen lane's offset lies within four
    of lane 0's: the sixteen bytes from four before it hold each lane's eight, which a shift a lane takes out.
    The bytes must stay readable four before lane 0's offset and twelve past it."""
    comptime Wide = SIMD[DType.uint64, LANES]
    var anchor = Int(offsets[0])
    var start = anchor - 4
    var within = lanes.select((offsets - Int32(start)).cast[DType.uint32](), SIMD[DType.uint32, LANES](0))
    if lanes[0] and within.le(UInt32(8)).reduce_and():
        var low = Wide(base.unsafe_offset(start).unsafe_bitcast[UInt64]().unsafe_load())
        var high = Wide(base.unsafe_offset(start + 8).unsafe_bitcast[UInt64]().unsafe_load())
        var bits = within.cast[DType.uint64]() * 8
        return lanes.select(shifted_right(low, bits) | shifted_left(high, Wide(64) - bits), elsewhere)
    return gathered_words(base, offsets, lanes, elsewhere)


@inline(.always)
def windowed_slides(
    first: ImmPointer[UInt8, _],
    second: ImmPointer[UInt8, _],
    entries: SIMD[DType.int32, LANES],
    diagonals: SIMD[DType.int32, LANES],
) -> SIMD[DType.int32, LANES]:
    """`gathered_slides` with `windowed_words` for each sequence: a lane group on neighbouring diagonals
    mostly reads within a few letters of each, where two loads and a shift a lane cost less than a gather.
    Both sequences must be readable four letters before every entry and twelve past it."""
    comptime Wide = SIMD[DType.uint64, LANES]
    var reached = entries.ge(0)
    var mismatches = windowed_words(first, entries, reached, Wide(0)) ^ windowed_words(
        second, entries - diagonals, reached, ~Wide(0)
    )
    var slid = entries + (count_trailing_zeros(mismatches) >> 3).cast[DType.int32]()
    var whole = mismatches.eq(0)
    if whole.reduce_or():
        comptime for lane in range(LANES):
            if whole[lane]:
                slid[lane] = Int32(slide(first, second, Int(slid[lane]), Int(diagonals[lane])))
    return reached.select(slid, entries)
