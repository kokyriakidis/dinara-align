# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
What the runner built with probes counts around each pair, beside the library's own probes: the
processor's counters (instructions, cycles, branch misses, first-level and last-level cache misses), read through Linux's `perf_event_open`, user space only; the allocations and bulk copies that
reach the C library, from `probes/counting.c` when it is preloaded; and the page faults, from
`getrusage`. Each is read before and after a pair, outside its timed span, and the runner prints the
differences. Where a source is missing (another system, `perf_event_paranoid` too strict, no preload),
its counts are left out rather than reported as zero.
"""

from std.ffi import external_call
from std.sys.info import CompilationTarget

from dinara_align.probes import counter_page, read_counter, UNCOUNTED


comptime EVENTS = 5
"""Hardware events counted: see `Counters.json`."""
comptime SHIMS = 10
"""`counting.c`'s counts, in its `dinara_counts` order: see `Counters.json`."""


def perf_event_open(kind: Int, config: Int) -> Int:
    """A counter of `kind` and `config` for this process on any CPU, user space only, counting from now,
    with the times it was enabled and running so a multiplexed count can be told apart; -1 where none."""
    comptime if not CompilationTarget.is_linux():
        return -1
    # `struct perf_event_attr` as its first version, 64 bytes: type, size, config, then the sample
    # period and type, the read format and the flags word.
    var attribute = List[UInt64](length=8, fill=0)
    attribute[0] = UInt64(kind) | (UInt64(64) << 32)
    attribute[1] = UInt64(config)
    # PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING
    attribute[4] = 3
    # exclude_kernel | exclude_hv
    attribute[5] = (1 << 5) | (1 << 6)
    comptime number = 298 if CompilationTarget.is_x86() else 241
    return external_call["syscall", Int](number, attribute.unsafe_ptr(), 0, -1, -1, 0)


struct Counters(Movable):
    """The hardware counters, the shim's counts and the page faults, read around each pair into buffers
    allocated once, so that reading them allocates nothing the shim would count against the pair."""

    var descriptors: List[Int]
    """Each event's file descriptor, -1 where it could not be opened."""
    var before: List[Int]
    """The snapshot a pair's counts start from: the events, the shim's counts, the minor and major faults."""
    var after: List[Int]
    """The snapshot they end at."""
    var buffer: List[UInt64]
    """One event's value and its times enabled and running."""
    var usage: List[Int64]
    """`getrusage`'s answer."""
    var shim: Int
    """The address of `counting.c`'s `dinara_counts`, 0 where it is not preloaded."""
    var instruction_page: Int
    """The instruction counter's page, which an aligner's tally counts its phases on (see `Tally.use_counters`)."""
    var cycle_page: Int
    """The cycle counter's page, likewise."""
    var multiplexed: Bool
    """Whether any event ran for less time than it was enabled, its count then scaled."""

    def __init__(out self):
        """Opens every event and finds the shim."""
        self.descriptors = List[Int]()
        # PERF_TYPE_HARDWARE, then PERF_TYPE_HW_CACHE as cache | read << 8 | miss << 16. Five events, two on the
        # fixed counters and three on general ones, fit the Skylake-X beside the NMI watchdog; a sixth made
        # the kernel rotate them, and a rotated count is an estimate.
        var kinds = [0, 0, 0, 3, 3]
        var configs = [1, 0, 5, 0x10000, 0x10002]
        # Instructions and cycles by `rdpmc` from pages the library opens, shared with each aligner's tally
        # (see `Tally.use_counters`), so no event is open twice and every one fits on the processor at once.
        self.instruction_page = counter_page(1)
        self.cycle_page = counter_page(0)
        for event in range(EVENTS):
            if event < 2:
                self.descriptors.append(-1)
            else:
                self.descriptors.append(perf_event_open(kinds[event], configs[event]))
        var name = String("dinara_counts")
        var handle = -2 if CompilationTarget.is_macos() else 0
        self.shim = external_call["dlsym", Int](handle, name.as_c_string_span().ptr())
        self.before = List[Int](length=EVENTS + SHIMS + 2, fill=-1)
        self.after = List[Int](length=EVENTS + SHIMS + 2, fill=-1)
        self.buffer = List[UInt64](length=3, fill=0)
        self.usage = List[Int64](length=20, fill=0)
        self.multiplexed = False

    def put[end: Bool](mut self, index: Int, value: Int):
        """Sets one count of the snapshot `end` says: the one a pair ends at, or the one it starts from."""
        comptime if end:
            self.after[index] = value
        else:
            self.before[index] = value

    def snapshot[end: Bool](mut self):
        """Every count now into the snapshot `end` says, -1 for an event or source that is missing; allocates
        nothing."""
        for event in range(EVENTS):
            self.put[end](event, -1)
            var descriptor = self.descriptors[event]
            if descriptor < 0:
                continue
            if external_call["read", Int](descriptor, self.buffer.unsafe_ptr(), 24) != 24:
                continue
            if self.buffer[2] < self.buffer[1]:
                self.multiplexed = True
                self.put[end](
                    event, Int(Float64(self.buffer[0]) * Float64(self.buffer[1]) / Float64(max(self.buffer[2], 1)))
                )
            else:
                self.put[end](event, Int(self.buffer[0]))
        for event in range(2):
            var page = self.instruction_page if event == 0 else self.cycle_page
            if page == 0:
                continue
            var value = read_counter(page)
            if value == UNCOUNTED:
                self.multiplexed = True
            else:
                self.put[end](event, value)
        for slot in range(SHIMS):
            self.put[end](EVENTS + slot, -1)
            if self.shim != 0:
                self.put[end](
                    EVENTS + slot,
                    Int(MutPointer[UInt64, MutAnyOrigin](unsafe_from_address=self.shim)[unsafe_offset=slot]),
                )
        # `struct rusage`: two 16-byte times, then maxrss, ixrss, idrss, isrss, minflt, majflt.
        _ = external_call["getrusage", Int32](Int32(0), self.usage.unsafe_ptr())
        self.put[end](EVENTS + SHIMS, Int(self.usage[8]))
        self.put[end](EVENTS + SHIMS + 1, Int(self.usage[9]))

    def start(mut self):
        """Takes the snapshot a pair's counts are measured from."""
        self.snapshot[False]()

    def stop(mut self):
        """Takes the snapshot they end at."""
        self.snapshot[True]()

    def json(self) -> String:
        """The counts from `start` to `stop`, as a JSON object's members, the missing ones left out."""
        # The events, `counting.c`'s counts, then the page faults, in snapshot order.
        var names: List[StaticString] = [
            "instructions",
            "cycles",
            "branch_misses",
            "l1d_misses",
            "llc_misses",
            "allocations",
            "allocated_bytes",
            "frees",
            "reallocations",
            "copies",
            "copied_bytes",
            "moves",
            "moved_bytes",
            "sets",
            "set_bytes",
            "minor_faults",
            "major_faults",
        ]
        var out = String()
        for index in range(len(names)):
            if self.after[index] < 0 or self.before[index] < 0:
                continue
            if out.byte_length() > 0:
                out += ", "
            out.write('"', names[index], '": ', self.after[index] - self.before[index])
        if self.multiplexed:
            out += ', "multiplexed": true'
        return out^
