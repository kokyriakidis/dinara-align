# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
Probes for profiling, compiled in only by a build with `-D DINARA_PROBES` (see `benchmarks/profile.py`):
per pair, the engine that settled it, the time each phase took and the work each did. A normal build
defines no `DINARA_PROBES`, every probe is then a constant branch the compiler folds away, and a
`Tally` holds an empty list it never touches.

A `Tally` lives in each `Aligner`'s space (see `EditSpace`), so it is per thread like the rest of the
library's memory; the runner reads it after each pair and clears it. The slots below index it; `json`
names them. Each phase also counts the instructions and cycles it retired, read by `rdpmc` from counters
the tally opens through `perf_event_open` on Linux x86 the first time it is started, so a phase costs a
few dozen cycles more, not a system call.

A build may also force a policy on the engines' choices, `-D DINARA_FORCE=<policy>` (see `FORCE`), so that
`profile.py regret` can time each pair under every choice the library could have made for it.
"""

from std.ffi import external_call
from std.sys import get_defined_string, is_defined, llvm_intrinsic
from std.sys.info import CompilationTarget
from std.time import perf_counter_ns


comptime PROBES = is_defined["DINARA_PROBES"]()
"""Whether this build records probes."""

comptime FORCE = get_defined_string["DINARA_FORCE", ""]()
"""A policy forced on the engines' choices, every one of which keeps the answer exact, or none:

- `band`: no diagonal transition settles a pair; band doubling takes it from the one front's projection.
- `diagonal`: diagonal transition never hands over on a projection, only at its own limit.
- `full`: the whole matrix, swept for a distance, or as one band for an alignment.
- `seeds-none`, `seeds-exact`, `seeds-inexact`: the band's seed heuristic, whatever the gate would choose.
- `sweep`: a gap-affine pair the difference sweep can take takes it from the start.
- `no-sweep`: the wavefront never hands a gap-affine pair to the sweep.
- `bidirectional`: a gap-affine pair skips the one-search phase, both searches growing from the start.
"""
comptime FORCE_BAND = FORCE == "band"
comptime FORCE_DIAGONAL = FORCE == "diagonal"
comptime FORCE_FULL = FORCE == "full"
comptime FORCE_SEEDS = FORCE.startswith("seeds-")
comptime FORCE_SWEEP = FORCE == "sweep"
comptime FORCE_NO_SWEEP = FORCE == "no-sweep"
comptime FORCE_BIDIRECTIONAL = FORCE == "bidirectional"
comptime FORCED_STEP_TENTHS = 1 << 30
"""The step budget a forced diagonal transition takes: past any projection, so only its limit stops it."""

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
# Phases of the paths beside the global ones: the bit-parallel search for an end (prefix, suffix, infix),
# and the scored modes (local, extension, substitution tables).
comptime SEARCH_NS = 52
comptime SCORED_NS = 53
comptime SLOTS = 54


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
comptime ENGINE_SEARCH = 11
comptime ENGINE_SCORED = 12


@always_inline
def clock() -> Int:
    """Nanoseconds now in a build with probes; 0, at no cost, in any other."""
    comptime if PROBES:
        return Int(perf_counter_ns())
    return 0


def perf_event_open(config: Int) -> Int:
    """A hardware counter of `config` (`PERF_TYPE_HARDWARE`) for this thread, user space only, counting from
    now; -1 where none."""
    comptime if not CompilationTarget.is_linux():
        return -1
    # `struct perf_event_attr` as its first version, 64 bytes: type, size, config, the sample period and
    # type, the read format, then the flags, exclude_kernel and exclude_hv.
    var attribute = List[UInt64](length=8, fill=0)
    attribute[0] = UInt64(64) << 32
    attribute[1] = UInt64(config)
    attribute[5] = (1 << 5) | (1 << 6)
    comptime number = 298 if CompilationTarget.is_x86() else 241
    return external_call["syscall", Int](number, attribute.unsafe_ptr(), 0, -1, -1, 0)


def counter_page(config: Int) -> Int:
    """The address of the page through which user space reads a hardware counter of `config` by `rdpmc`, or
    0 where it cannot: off Linux x86, without the event, or without `cap_user_rdpmc`."""
    comptime if not (CompilationTarget.is_linux() and CompilationTarget.is_x86()):
        return 0
    var descriptor = perf_event_open(config)
    if descriptor < 0:
        return 0
    # PROT_READ, MAP_SHARED
    var page = external_call["mmap", Int](0, 4096, 1, 1, descriptor, 0)
    if page == -1 or page == 0:
        return 0
    # `struct perf_event_mmap_page`: capabilities at byte 40, bit 2 `cap_user_rdpmc`.
    var capabilities = MutPointer[UInt64, MutAnyOrigin](unsafe_from_address=page + 40).unsafe_load[volatile=True]()
    if (capabilities >> 2) & 1 == 0:
        return 0
    return page


@always_inline
def read_counter(page: Int) -> Int:
    """The counter behind `page` now, by `rdpmc` under the page's sequence lock; 0 where there is none."""
    comptime if not CompilationTarget.is_x86():
        return 0
    if page == 0:
        return 0
    var words = MutPointer[UInt32, MutAnyOrigin](unsafe_from_address=page)
    while True:
        var sequence = words.unsafe_load[volatile=True](2)
        var index = words.unsafe_load[volatile=True](3)
        var offset = MutPointer[Int64, MutAnyOrigin](unsafe_from_address=page + 16).unsafe_load[volatile=True]()
        var width = Int(MutPointer[UInt16, MutAnyOrigin](unsafe_from_address=page + 48).unsafe_load[volatile=True]())
        if index == 0:
            return 0
        var count = llvm_intrinsic["llvm.x86.rdpmc", Int64](Int32(index - 1))
        count = (count << Int64(64 - width)) >> Int64(64 - width)
        if words.unsafe_load[volatile=True](2) == sequence:
            return Int(count + offset)


struct Tally(Movable):
    """One pair's probes, by slot, and each phase's instructions and cycles; empty and untouched in a build
    without probes."""

    var values: List[Int]
    var instructions: List[Int]
    """Instructions retired in each phase slot."""
    var cycles: List[Int]
    """Cycles in each phase slot."""
    var instruction_page: Int
    """The instruction counter's page (see `counter_page`), 0 where there is none."""
    var cycle_page: Int
    """The cycle counter's page."""
    var opened: Bool
    """Whether the counters have been opened."""
    var marked_instructions: Int
    """The instructions counted when the current phase started."""
    var marked_cycles: Int
    """The cycles counted when it started."""

    def __init__(out self):
        """Every slot zero, or no slots at all without probes; the counters opened when first started."""
        comptime if PROBES:
            self.values = List[Int](length=SLOTS, fill=0)
            self.instructions = List[Int](length=SLOTS, fill=0)
            self.cycles = List[Int](length=SLOTS, fill=0)
        else:
            self.values = List[Int]()
            self.instructions = List[Int]()
            self.cycles = List[Int]()
        self.instruction_page = 0
        self.cycle_page = 0
        self.opened = False
        self.marked_instructions = 0
        self.marked_cycles = 0

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
    def start(mut self) -> Int:
        """Marks the start of a phase: its counters now, and the time returned, for `since`."""
        comptime if PROBES:
            if not self.opened:
                self.opened = True
                # PERF_COUNT_HW_INSTRUCTIONS, PERF_COUNT_HW_CPU_CYCLES
                self.instruction_page = counter_page(1)
                self.cycle_page = counter_page(0)
            self.marked_instructions = read_counter(self.instruction_page)
            self.marked_cycles = read_counter(self.cycle_page)
        return clock()

    @always_inline
    def since[slot: Int](mut self, started: Int):
        """Adds the nanoseconds since `started`, a `start()`, to `slot`, and the instructions and cycles
        since that mark to the phase's counts."""
        comptime if PROBES:
            self.values[slot] += clock() - started
            if self.instruction_page != 0:
                self.instructions[slot] += read_counter(self.instruction_page) - self.marked_instructions
            if self.cycle_page != 0:
                self.cycles[slot] += read_counter(self.cycle_page) - self.marked_cycles

    def clear(mut self):
        """Every slot back to zero, for the next pair."""
        comptime if PROBES:
            for slot in range(SLOTS):
                self.values[slot] = 0
                self.instructions[slot] = 0
                self.cycles[slot] = 0

    def json(self) -> String:
        """The slots that are not zero as a JSON object's members, the engine by name, and each phase's
        instructions and cycles."""
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
                "search_ns",
                "scored_ns",
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
                "search",
                "scored",
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
            # The phase slots, whose instructions and cycles are counted beside their time.
            var phases: List[Int] = [
                ENCODE_NS,
                DIAGONAL_NS,
                TWO_ENDED_NS,
                GROW_NS,
                SEEDS_NS,
                BAND_NS,
                TRACEBACK_NS,
                CIGAR_NS,
                FULL_NS,
                WAVEFRONT_NS,
                SWEEP_NS,
                SEARCH_NS,
                SCORED_NS,
            ]
            for slot in phases:
                var stem = names[slot].removesuffix("_ns")
                if self.instructions[slot] != 0:
                    out.write(', "', stem, '_instructions": ', self.instructions[slot])
                if self.cycles[slot] != 0:
                    out.write(', "', stem, '_cycles": ', self.cycles[slot])
        return out^


def front_steps(score: Int, columns: Int, rows: Int) -> Int:
    """The diagonals a front of a `columns` by `rows` matrix steps over scores 0 to `score`: each score
    `s` spans the diagonals from `-min(s, rows)` to `min(s, columns)`."""
    var steps = 0
    for at in range(score + 1):
        steps += min(at, columns) + min(at, rows) + 1
    return steps
