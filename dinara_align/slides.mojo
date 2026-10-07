"""
How far matches carry a cell along its diagonal: eight letters at a time, by one scalar compare, or on
AVX-512 for a lane group of diagonals at once by two gathers (see `gathered_slides`). Both sequences
end in sentinels that differ from every letter and from each other, so a slide stops at the edge.
The diagonal transition and the affine wavefront both slide their fronts this way.
"""

from std.bit import count_trailing_zeros
from std.sys import llvm_intrinsic
from std.sys.info import CompilationTarget

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


comptime GATHERED_SLIDES = CompilationTarget.has_avx512f()
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
