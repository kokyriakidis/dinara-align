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
- `one-ended`, `two-ended`: an alignment's diagonal transition never gives way to the two-ended search, or
  skips the one front for it from the start.
- `aim-search`: a band aims its first bound at the search's own projection, from every edit it reached, and
  trusts it, rather than at the one front's from its first `PROBE_START`.
- `trace-wavefront`, `trace-recompute`: a band's traceback traces each tile by the wavefront wherever it fits
  its work cap, or always recomputes it.
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
# The search's own projection of the distance, from all the edits it reached, where `ESTIMATE` holds the one
# front's from its first `PROBE_START` that the band aims at: the two `trusted_projection` compares.
comptime SEARCH_ESTIMATE = 54
comptime SLOTS = 55


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


trait Counting(Movable):
    """A struct's own probes, by slot: real counts in a build with probes, nothing at all in any other."""

    def __init__(out self):
        """Every count zero."""
        ...

    def add(mut self, slot: Int, amount: Int):
        """Adds `amount` to `slot`."""
        ...

    def get(self, slot: Int) -> Int:
        """`slot`'s count."""
        ...

    def reset(mut self):
        """Every count back to zero."""
        ...


struct Counts(Counting):
    """Eight counts, the most any struct's probes take; a slot known when the call is inlined picks its field
    with no branch left."""

    var c0: Int
    var c1: Int
    var c2: Int
    var c3: Int
    var c4: Int
    var c5: Int
    var c6: Int
    var c7: Int

    def __init__(out self):
        """Every count zero."""
        self.c0 = 0
        self.c1 = 0
        self.c2 = 0
        self.c3 = 0
        self.c4 = 0
        self.c5 = 0
        self.c6 = 0
        self.c7 = 0

    @always_inline
    def add(mut self, slot: Int, amount: Int):
        """Adds `amount` to `slot`, 0 to 7."""
        if slot == 0:
            self.c0 += amount
        elif slot == 1:
            self.c1 += amount
        elif slot == 2:
            self.c2 += amount
        elif slot == 3:
            self.c3 += amount
        elif slot == 4:
            self.c4 += amount
        elif slot == 5:
            self.c5 += amount
        elif slot == 6:
            self.c6 += amount
        else:
            self.c7 += amount

    @always_inline
    def get(self, slot: Int) -> Int:
        """`slot`'s count, 0 to 7."""
        if slot == 0:
            return self.c0
        if slot == 1:
            return self.c1
        if slot == 2:
            return self.c2
        if slot == 3:
            return self.c3
        if slot == 4:
            return self.c4
        if slot == 5:
            return self.c5
        if slot == 6:
            return self.c6
        return self.c7

    def reset(mut self):
        """Every count back to zero."""
        self = Self()


struct NoCounts(Counting):
    """No counts: a field of this type takes no bytes, so a struct carrying it keeps a normal build's layout,
    and every call on it is gone once inlined."""

    def __init__(out self):
        """Nothing to count."""
        pass

    @always_inline
    def add(mut self, slot: Int, amount: Int):
        """Nothing."""
        pass

    @always_inline
    def get(self, slot: Int) -> Int:
        """Nothing counted: 0."""
        return 0

    def reset(mut self):
        """Nothing."""
        pass


comptime ProbeCounts: Counting = Counts if PROBES else NoCounts
"""The probes a hot struct carries itself (see `Counting`): real in a build with probes, zero bytes otherwise.
Adding plain fields to those structs changed a normal build's layout and code, about 1% more instructions on
long pairs."""


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


comptime UNCOUNTED = -1
"""What `read_counter` gives for a counter the kernel has not got on the processor right now, as when more
events are open than it has counters and it rotates them: no `rdpmc` reading is then meaningful."""


@always_inline
def running_gap(page: Int) -> Int:
    """How much longer the counter behind `page` has been enabled than running, in nanoseconds: a change in it
    across a phase means the kernel rotated the counter off the processor during it."""
    if page == 0:
        return 0
    var enabled = MutPointer[UInt64, MutAnyOrigin](unsafe_from_address=page + 24).unsafe_load[volatile=True]()
    var running = MutPointer[UInt64, MutAnyOrigin](unsafe_from_address=page + 32).unsafe_load[volatile=True]()
    return Int(enabled) - Int(running)


@always_inline
def read_counter(page: Int) -> Int:
    """The counter behind `page` now, by `rdpmc` under the page's sequence lock; 0 where there is none, and
    `UNCOUNTED` while the kernel has it off the processor."""
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
            return UNCOUNTED
        var count = llvm_intrinsic["llvm.x86.rdpmc", Int64](Int32(index - 1))
        count = (count << Int64(64 - width)) >> Int64(64 - width)
        if words.unsafe_load[volatile=True](2) == sequence:
            return Int(count + offset)


trait Tallying(Movable):
    """What a pair's probes offer the library (see `Tally`)."""

    def __init__(out self):
        """Nothing counted yet."""
        ...

    def use_counters(mut self, instruction_page: Int, cycle_page: Int):
        """Counts the phases on counters a caller already opened."""
        ...

    def add(mut self, slot: Int, amount: Int):
        """Adds `amount` to `slot`."""
        ...

    def get(self, slot: Int) -> Int:
        """`slot`'s value."""
        ...

    def set(mut self, slot: Int, value: Int):
        """Sets `slot` to `value`."""
        ...

    def start(mut self) -> Int:
        """Marks the start of a phase, for `since`."""
        ...

    def since(mut self, slot: Int, started: Int):
        """Adds the time, instructions and cycles since `started` to the phase `slot`."""
        ...

    def clear(mut self):
        """Every slot back to zero."""
        ...

    def json(self) -> String:
        """The probes as a JSON object's members."""
        ...


struct Tallied(Tallying):
    """One pair's probes, by slot, and each phase's instructions and cycles: the `Tally` of a build with probes."""

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
    var uncounted: Bool
    """Whether a phase of this pair met a counter off the processor, its counts then left out."""
    var marked_gap: Int
    """The instruction and cycle counters' gaps between enabled and running when the phase started (see
    `running_gap`), summed."""

    def __init__(out self):
        """Every slot zero; the counters opened when first started."""
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
        self.uncounted = False
        self.marked_gap = 0

    def use_counters(mut self, instruction_page: Int, cycle_page: Int):
        """Counts the phases on counters a caller already opened (see `counter_page`), rather than open its
        own: two of each event would need more counters than the processor has, and rotated counters
        give no reading."""
        self.instruction_page = instruction_page
        self.cycle_page = cycle_page
        self.opened = True

    @always_inline
    def add(mut self, slot: Int, amount: Int):
        """Adds `amount` to `slot`."""
        comptime if PROBES:
            self.values[slot] += amount

    @always_inline
    def get(self, slot: Int) -> Int:
        """`slot`'s value."""
        comptime if PROBES:
            return self.values[slot]
        return 0

    @always_inline
    def set(mut self, slot: Int, value: Int):
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
            self.marked_gap = running_gap(self.instruction_page) + running_gap(self.cycle_page)
            self.marked_instructions = read_counter(self.instruction_page)
            self.marked_cycles = read_counter(self.cycle_page)
        return clock()

    @always_inline
    def since(mut self, slot: Int, started: Int):
        """Adds the nanoseconds since `started`, a `start()`, to `slot`, and the instructions and cycles
        since that mark to the phase's counts."""
        comptime if PROBES:
            self.values[slot] += clock() - started
            var instructions = read_counter(self.instruction_page)
            var cycles = read_counter(self.cycle_page)
            # A counter rotated off the processor during the phase, or off it now: no reading means anything.
            var rotated = running_gap(self.instruction_page) + running_gap(self.cycle_page) != self.marked_gap
            if (
                rotated
                or instructions == UNCOUNTED
                or cycles == UNCOUNTED
                or self.marked_instructions == UNCOUNTED
                or self.marked_cycles == UNCOUNTED
            ):
                self.uncounted = True
                return
            self.instructions[slot] += instructions - self.marked_instructions
            self.cycles[slot] += cycles - self.marked_cycles

    def clear(mut self):
        """Every slot back to zero, for the next pair."""
        comptime if PROBES:
            for slot in range(SLOTS):
                self.values[slot] = 0
                self.instructions[slot] = 0
                self.cycles[slot] = 0
            self.uncounted = False

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
                "search_estimate",
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
            if self.uncounted:
                out.write(', "phase_counters": "uncounted"')
                return out^
            for slot in phases:
                var stem = names[slot].removesuffix("_ns")
                if self.instructions[slot] != 0:
                    out.write(', "', stem, '_instructions": ', self.instructions[slot])
                if self.cycles[slot] != 0:
                    out.write(', "', stem, '_cycles": ', self.cycles[slot])
        return out^


struct Untallied(Tallying):
    """No probes: the `Tally` of a normal build, which takes no bytes, so `EditSpace` keeps its size, and every
    call on it is gone once inlined."""

    def __init__(out self):
        """Nothing to count."""
        pass

    def use_counters(mut self, instruction_page: Int, cycle_page: Int):
        """Nothing."""
        pass

    @always_inline
    def add(mut self, slot: Int, amount: Int):
        """Nothing."""
        pass

    @always_inline
    def get(self, slot: Int) -> Int:
        """Nothing counted: 0."""
        return 0

    @always_inline
    def set(mut self, slot: Int, value: Int):
        """Nothing."""
        pass

    @always_inline
    def start(mut self) -> Int:
        """Nothing: 0."""
        return 0

    @always_inline
    def since(mut self, slot: Int, started: Int):
        """Nothing."""
        pass

    def clear(mut self):
        """Nothing."""
        pass

    def json(self) -> String:
        """No probes: empty."""
        return String()


comptime Tally: Tallying = Tallied if PROBES else Untallied
"""One pair's probes as the library carries them, in each `Aligner`'s space (see `EditSpace`): `Tallied` in a
build with probes, `Untallied`, no bytes and no work, in any other."""


def front_steps(score: Int, columns: Int, rows: Int) -> Int:
    """The diagonals a front of a `columns` by `rows` matrix steps over scores 0 to `score`: each score
    `s` spans the diagonals from `-min(s, rows)` to `min(s, columns)`."""
    var steps = 0
    for at in range(score + 1):
        steps += min(at, columns) + min(at, rows) + 1
    return steps
