# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Errors: every call that can fail raises an `AlignmentError`, its `kind` saying which way it failed and
its `detail` what it failed on.

One type for the whole package, since a Mojo function declares a single error type and a caller would
otherwise convert between them at every call.
"""


@fieldwise_init
struct ErrorKind(Equatable, ImplicitlyCopyable, TrivialRegisterPassable, Writable):
    """The way a call failed, which callers branch on: the C API turns each into its own code."""

    var id: UInt8

    comptime UNKNOWN_SYMBOL = Self(1)
    """A sequence holds a letter its `Scoring`'s alphabet lacks, or a byte no sequence may hold."""
    comptime ALPHABET_TOO_LARGE = Self(2)
    """An alphabet with more letters than a substitution table holds."""
    comptime SEQUENCE_TOO_LONG = Self(3)
    """The input needs more memory than it may take: past `max_memory`, or past what the device can
    allocate at once."""
    comptime SCRATCH_TOO_SMALL = Self(4)
    """Device memory sized for a smaller problem than the one launched on it."""
    comptime LENGTH_MISMATCH = Self(5)
    """Two lists, or two gapped rows, of unequal length where each item pairs with one of the other."""
    comptime INVALID_SCORING = Self(6)
    """Costs or scores no search can use: an edit that costs nothing, a reward that costs, a value past
    the range its arithmetic holds."""
    comptime INVALID_ARGUMENT = Self(7)
    """An argument outside what the call takes: a mode, a band, a cap or a count it refuses."""
    comptime OUTSIDE_BAND = Self(8)
    """Every alignment leaves the band of diagonals asked for."""

    def phrase(self) -> StaticString:
        """The failure in a few words, the start of every message of this kind."""
        if self == Self.UNKNOWN_SYMBOL:
            return "a letter the alphabet lacks"
        if self == Self.ALPHABET_TOO_LARGE:
            return "too many letters for a substitution table"
        if self == Self.SEQUENCE_TOO_LONG:
            return "the memory this input needs passes what it may take"
        if self == Self.SCRATCH_TOO_SMALL:
            return "device memory too small for its launch"
        if self == Self.LENGTH_MISMATCH:
            return "inputs of unequal length"
        if self == Self.INVALID_SCORING:
            return "costs or scores no search can use"
        if self == Self.INVALID_ARGUMENT:
            return "an argument outside what the call takes"
        if self == Self.OUTSIDE_BAND:
            return "no alignment stays inside the band"
        return "an unknown failure"

    def write_to(self, mut writer: Some[Writer]):
        """The kind's phrase."""
        writer.write(self.phrase())


@fieldwise_init
struct AlignmentError(Copyable, ImplicitlyCopyable, Writable):
    """A failed call: which way it failed, and what it failed on, the value or the input to look at."""

    var kind: ErrorKind
    var detail: String

    def write_to(self, mut writer: Some[Writer]):
        """`dinara-align: <phrase> [<detail>]`, the form the C API, the Python package and the command line
        pass on."""
        writer.write("dinara-align: ", self.kind.phrase(), " [", self.detail, "]")
