# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Build-time switches that each turn one technique off, for measuring what it contributes: a build
with `-D ABLATE_REGROUP` and the like against one without, everything else alike. A normal build
defines none of them, and each is then a constant false the compiler folds away. Every switch leaves
the answers exact; only the time changes.
"""

from std.sys import is_defined


comptime ABLATE_REGROUP = is_defined["ABLATE_REGROUP"]()
"""Myers' step as A*PA2 writes it, rather than regrouped for a shorter chain (see `bit_parallel.advance`)."""


comptime ABLATE_PAIRED = is_defined["ABLATE_PAIRED"]()
"""One group of eight words a sweep on AVX-512 rather than two interleaved (see `bit_parallel`'s
`PAIRED_SWEEP`); the thresholds tuned on the paired sweep stay as they are."""


comptime ABLATE_GATHER = is_defined["ABLATE_GATHER"]()
"""Diagonal transition's slides a scalar compare a lane rather than gathered (see `slides.GATHERED_SLIDES`)."""


comptime ABLATE_AGREEMENT = is_defined["ABLATE_AGREEMENT"]()
"""Every projection trusted, whether or not the two agree (see `diagonal.trusted_projection`)."""


comptime ABLATE_NEIGHBOURS = is_defined["ABLATE_NEIGHBOURS"]()
"""An exact inexact-seed match's four neighbours found by the scan and pruned each on its own, as before
`ebe290f`, rather than put in the layers from the exact match (see `seeds.SeedHeuristic.add_neighbours`)."""


comptime ABLATE_CERTIFIED = is_defined["ABLATE_CERTIFIED"]()
"""A batch's distances swept over each group's whole matrices rather than a band its pairs' costs certify
(see `lanes.lane_stage`)."""
