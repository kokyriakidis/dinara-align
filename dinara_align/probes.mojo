# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Probes for profiling, compiled in only by a build with `-D DINARA_PROBES` (see `benchmarks/profile.py`):
per pair, the engine that settled it, the time each phase took and the work each did. A normal build
defines no `DINARA_PROBES`, every probe is then a constant branch the compiler folds away, and a
`Tally` holds an empty list it never touches.

A `Tally` lives in each `Aligner`'s space (see `EditSpace`), so it is per thread like the rest of the
library's memory; the runner reads it after each pair and clears it. The slots below index it; `json`
names them.
"""

from std.sys import is_defined
from std.time import perf_counter_ns


comptime PROBES = is_defined["DINARA_PROBES"]()
"""Whether this build records probes."""

# Phases, in nanoseconds.
comptime ENCODE_NS = 0
comptime DIAGONAL_NS = 1
comptime TWO_ENDED_NS = 2
comptime GROW_NS = 3
comptime SEEDS_NS = 4
comptime BAND_NS = 5
comptime TRACEBACK_NS = 6
comptime CIGAR_NS = 7
comptime FULL_NS = 8
comptime WAVEFRONT_NS = 9
comptime SWEEP_NS = 10
# What settled the pair (see `ENGINES`) and what the diagonal transition projected.
comptime ENGINE = 11
comptime ESTIMATE = 12
comptime TRUSTED = 13
comptime FLOOR = 14
# Diagonal transition: diagonals stepped, and the score each search reached.
comptime DIAGONAL_STEPS = 15
comptime DIAGONAL_SCORE = 16
comptime TWO_ENDED_STEPS = 17
comptime TWO_ENDED_SCORE = 18
comptime GROW_STEPS = 19
# The seed heuristic: none (0), exact (1) or inexact (2) seeds, how many, how many chain, its value at
# the origin, whether exact seeds were rebuilt as inexact ones, matches found and kept by local pruning,
# and whether too many matches made it give up.
comptime SEED_KIND = 20
comptime SEEDS = 21
comptime SEEDS_COUNTED = 22
comptime HEURISTIC_ORIGIN = 23
comptime SEED_REBUILT = 24
comptime SEED_MATCHES = 25
comptime SEED_KEPT = 26
comptime SEED_GAVE_UP = 27
# Band doubling: rounds, the first and last bound, every round's bound summed, tiles swept, words times
# columns swept (cells over 64), checkpoints that gave a round up, bounds lifted to the whole matrix.
comptime ROUNDS = 28
comptime FIRST_BOUND = 29
comptime LAST_BOUND = 30
comptime BOUND_SUM = 31
comptime TILES = 32
comptime WORD_COLUMNS = 33
comptime CHECKPOINT_GIVE_UPS = 34
comptime WIDE_LIFTS = 35
comptime DISTANCE = 36
# The band's traceback: tiles traced, tiles swept again for edges the trail did not keep and their words
# times columns, tiles traced by the wavefront and tiles recomputed, the windows those took and their
# words times columns.
comptime TRACED_TILES = 37
comptime REPLAYED_TILES = 38
comptime REPLAYED_WORD_COLUMNS = 39
comptime FORWARD_SEGMENTS = 40
comptime RECOMPUTED_SEGMENTS = 41
comptime RECOMPUTE_WINDOWS = 42
comptime RECOMPUTE_WORD_COLUMNS = 43
# A full-matrix sweep's words times columns.
comptime FULL_WORD_COLUMNS = 44
# Gap-affine: cost steps of the forward and backward wavefronts, diagonals they computed, meeting
# checks, splits for the memory limit, and the difference sweep's cells.
comptime FORWARD_STEPS = 45
comptime BACKWARD_STEPS = 46
comptime WAVEFRONT_DIAGONALS = 47
comptime MEETS = 48
comptime SPLITS = 49
comptime SWEEP_CELLS = 50
# Band doubling again: the words times columns of the round that succeeded, the rest spent on rounds that
# failed.
comptime LAST_ROUND_WORD_COLUMNS = 51
comptime SLOTS = 52


# The engines `ENGINE` names.
comptime ENGINE_DIAGONAL = 1
comptime ENGINE_TWO_ENDED = 2
comptime ENGINE_SHORT_FULL = 3
comptime ENGINE_BAND = 4
comptime ENGINE_WIDE_FULL = 5
comptime ENGINE_WAVEFRONT = 6
comptime ENGINE_BIDIRECTIONAL = 7
comptime ENGINE_SWEEP = 8
comptime ENGINE_LANES = 9
comptime ENGINE_EMPTY = 10


@always_inline
def clock() -> Int:
    """Nanoseconds now in a build with probes; 0, at no cost, in any other."""
    comptime if PROBES:
        return Int(perf_counter_ns())
    return 0


struct Tally(Movable):
    """One pair's probes, by slot; empty and untouched in a build without probes."""

    var values: List[Int]

    def __init__(out self):
        """Every slot zero, or no slots at all without probes."""
        comptime if PROBES:
            self.values = List[Int](length=SLOTS, fill=0)
        else:
            self.values = List[Int]()

    @always_inline
    def add[slot: Int](mut self, amount: Int):
        """Adds `amount` to `slot`."""
        comptime if PROBES:
            self.values[slot] += amount

    @always_inline
    def get[slot: Int](self) -> Int:
        """`slot`'s value, 0 without probes."""
        comptime if PROBES:
            return self.values[slot]
        return 0

    @always_inline
    def set[slot: Int](mut self, value: Int):
        """Sets `slot` to `value`."""
        comptime if PROBES:
            self.values[slot] = value

    @always_inline
    def since[slot: Int](mut self, started: Int):
        """Adds the nanoseconds since `started`, a `clock()`, to `slot`."""
        comptime if PROBES:
            self.values[slot] += clock() - started

    def clear(mut self):
        """Every slot back to zero, for the next pair."""
        comptime if PROBES:
            for slot in range(SLOTS):
                self.values[slot] = 0

    def json(self) -> String:
        """The slots that are not zero as a JSON object's members, the engine by name."""
        var out = String()
        comptime if PROBES:
            var names: List[StaticString] = [
                "encode_ns",
                "diagonal_ns",
                "two_ended_ns",
                "grow_ns",
                "seeds_ns",
                "band_ns",
                "traceback_ns",
                "cigar_ns",
                "full_ns",
                "wavefront_ns",
                "sweep_ns",
                "engine",
                "estimate",
                "trusted",
                "floor",
                "diagonal_steps",
                "diagonal_score",
                "two_ended_steps",
                "two_ended_score",
                "grow_steps",
                "seed_kind",
                "seeds",
                "seeds_counted",
                "heuristic_origin",
                "seed_rebuilt",
                "seed_matches",
                "seed_kept",
                "seed_gave_up",
                "rounds",
                "first_bound",
                "last_bound",
                "bound_sum",
                "tiles",
                "word_columns",
                "checkpoint_give_ups",
                "wide_lifts",
                "distance",
                "traced_tiles",
                "replayed_tiles",
                "replayed_word_columns",
                "forward_segments",
                "recomputed_segments",
                "recompute_windows",
                "recompute_word_columns",
                "full_word_columns",
                "forward_steps",
                "backward_steps",
                "wavefront_diagonals",
                "meets",
                "splits",
                "sweep_cells",
                "last_round_word_columns",
            ]
            var engines: List[StaticString] = [
                "none",
                "diagonal",
                "two-ended",
                "short full",
                "band",
                "wide full",
                "wavefront",
                "bidirectional",
                "sweep",
                "lanes",
                "empty",
            ]
            for slot in range(SLOTS):
                var value = self.values[slot]
                if value == 0:
                    continue
                if out.byte_length() > 0:
                    out += ", "
                if slot == ENGINE:
                    out.write('"engine": "', engines[value], '"')
                else:
                    out.write('"', names[slot], '": ', value)
        return out^


def front_steps(score: Int, columns: Int, rows: Int) -> Int:
    """The diagonals a front of a `columns` by `rows` matrix steps over scores 0 to `score`: each score
    `s` spans the diagonals from `-min(s, rows)` to `min(s, columns)`."""
    var steps = 0
    for at in range(score + 1):
        steps += min(at, columns) + min(at, rows) + 1
    return steps
