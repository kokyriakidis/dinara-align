"""
The edit distance of a batch of pairs on the GPU: one thread a pair, each running Myers' bit-vector
recurrence over its pair, as Hyyrö's blocks and Edlib's `calculateBlock` write it for patterns longer
than a word.

The shorter sequence of each pair is the pattern, its letters down the bits of up to
`MAX_PATTERN_WORDS` words, and the longer one the text, a column a step. The host builds each pair's
match vectors, one row of words for each symbol the pattern holds and one, empty, for every other,
and codes the text by the same symbols; a thread then reads only its own pair's rows and text. The
pairs go out longest first, so a warp's threads run about as long as each other. A pair whose
pattern is longer stays on the host (see `api.distances`).
"""

from std.math import ceildiv

from max.algorithm import parallelize
from max.gpu import global_idx

from .common import DeviceScope, allocate, hardware_threads, upload, zeroed

comptime MAX_PATTERN_WORDS = 64
"""Words of the longest pattern a thread holds, 4,096 letters: its two vectors stay in its own memory."""

comptime THREADS_PER_BLOCK = 128
"""Threads, and so pairs, a block carries."""


def myers_kernel(
    rows: ImmPointer[UInt64, MutAnyOrigin],
    row_offsets: ImmPointer[Int64, MutAnyOrigin],
    words_of: ImmPointer[Int32, MutAnyOrigin],
    lengths: ImmPointer[Int32, MutAnyOrigin],
    text: ImmPointer[UInt8, MutAnyOrigin],
    text_offsets: ImmPointer[Int64, MutAnyOrigin],
    results: MutPointer[Int32, MutAnyOrigin],
    pairs: Int32,
):
    """Pair `global_idx.x`'s edit distance into `results`: its pattern of `lengths` letters in
    `words_of` words, its match rows from `row_offsets`, a symbol's row of words at the symbol's code
    times the words, and its text's codes between two `text_offsets`. The pattern's whole first column
    costs its letters, and each text letter one more on the top row, a global alignment's borders."""
    var pair = global_idx.x
    if pair >= Int(pairs):
        return
    var words = Int(words_of[unsafe_offset=pair])
    var length = Int(lengths[unsafe_offset=pair])
    var base = Int(row_offsets[unsafe_offset=pair])
    var last_bit = UInt64((length - 1) % 64)
    var positive = Array[UInt64, MAX_PATTERN_WORDS](fill=~UInt64(0))
    var negative = Array[UInt64, MAX_PATTERN_WORDS](fill=UInt64(0))
    var score = length
    for position in range(Int(text_offsets[unsafe_offset=pair]), Int(text_offsets[unsafe_offset=pair + 1])):
        var row = base + Int(text[unsafe_offset=position]) * words
        # The top row's cell grows by one a column: a horizontal step of +1 into the first word.
        var carry = 1
        for word in range(words):
            var equal = rows[unsafe_offset=row + word]
            var vertical_plus = positive[word]
            var vertical_minus = negative[word]
            # The step into this word's top from the word above, -1, 0 or +1.
            var carry_negative = UInt64(1) if carry < 0 else UInt64(0)
            var carry_positive = UInt64(1) if carry > 0 else UInt64(0)
            var crossing = equal | vertical_minus
            equal |= carry_negative
            var horizontal = (((equal & vertical_plus) + vertical_plus) ^ vertical_plus) | equal
            var plus = vertical_minus | ~(horizontal | vertical_plus)
            var minus = vertical_plus & horizontal
            # The step out of its bottom: the last word's at the pattern's last letter.
            var bit = last_bit if word == words - 1 else UInt64(63)
            carry = Int((plus >> bit) & 1) - Int((minus >> bit) & 1)
            plus = (plus << 1) | carry_positive
            minus = (minus << 1) | carry_negative
            positive[word] = minus | ~(crossing | plus)
            negative[word] = plus & crossing
        score += carry
    results[unsafe_offset=pair] = Int32(score)


def device_edit_distances(scope: DeviceScope, patterns: List[String], texts: List[String]) raises -> List[Int]:
    """Every pair's edit distance on the device, `patterns[i]` against `texts[i]`, each pattern of one
    to `64 MAX_PATTERN_WORDS` letters."""
    var pairs = len(patterns)
    # First each pair's symbols, which size its rows, then where its rows and its text go, then both
    # written, every pair on its own, over every thread.
    var symbols = List[Int](length=pairs, fill=0)
    var counts = symbols.unsafe_ptr()

    var workers = hardware_threads()
    # Chunks of pairs, several a thread, so a thread is never a task a pair.
    var chunks = max(min(pairs, workers * 8), 1)

    def count_symbols(chunk: Int) {imm patterns, imm counts, imm pairs, imm chunks}:
        for pair in range(pairs * chunk // chunks, pairs * (chunk + 1) // chunks):
            var seen = SIMD[DType.uint64, 4](0)
            var count = 0
            for letter in patterns[pair].as_bytes():
                var lane = Int(letter) // 64
                var bit = UInt64(1) << UInt64(Int(letter) % 64)
                if seen[lane] & bit == 0:
                    seen[lane] |= bit
                    count += 1
            counts[unsafe_offset=pair] = count

    parallelize(count_symbols, chunks, workers)
    var row_offsets = List[Int64](capacity=pairs + 1)
    var words_of = List[Int32](capacity=pairs)
    var lengths = List[Int32](capacity=pairs)
    var text_offsets = List[Int64](capacity=pairs + 1)
    var row_total = 0
    var text_total = 0
    for pair in range(pairs):
        var words = ceildiv(patterns[pair].byte_length(), 64)
        row_offsets.append(Int64(row_total))
        text_offsets.append(Int64(text_total))
        words_of.append(Int32(words))
        lengths.append(Int32(patterns[pair].byte_length()))
        row_total += (symbols[pair] + 1) * words
        text_total += texts[pair].byte_length()
    row_offsets.append(Int64(row_total))
    text_offsets.append(Int64(text_total))
    var rows = List[UInt64](length=row_total, fill=UInt64(0))
    var text = List[UInt8](length=text_total, fill=0)
    var rows_at = rows.unsafe_ptr()
    var text_at = text.unsafe_ptr()

    def fill(
        chunk: Int,
    ) {
        imm patterns,
        imm texts,
        imm rows_at,
        imm text_at,
        imm row_offsets,
        imm text_offsets,
        imm words_of,
        imm pairs,
        imm chunks,
    }:
        for pair in range(pairs * chunk // chunks, pairs * (chunk + 1) // chunks):
            fill_pair(pair, patterns, texts, rows_at, text_at, row_offsets, text_offsets, words_of)

    parallelize(fill, chunks, workers)

    var rows_buffer = upload(scope, rows)
    var offsets_buffer = upload(scope, row_offsets)
    var words_buffer = upload(scope, words_of)
    var lengths_buffer = upload(scope, lengths)
    var text_buffer = upload(scope, text)
    var text_offsets_buffer = upload(scope, text_offsets)
    var results_buffer = zeroed[DType.int32](scope, pairs)
    scope.context.enqueue_function[myers_kernel](
        rows_buffer.unsafe_ptr(),
        offsets_buffer.unsafe_ptr(),
        words_buffer.unsafe_ptr(),
        lengths_buffer.unsafe_ptr(),
        text_buffer.unsafe_ptr(),
        text_offsets_buffer.unsafe_ptr(),
        results_buffer.unsafe_ptr(),
        Int32(pairs),
        grid_dim=ceildiv(pairs, THREADS_PER_BLOCK),
        block_dim=THREADS_PER_BLOCK,
    )
    scope.context.synchronize()
    var results = List[Int](capacity=pairs)
    with results_buffer.map_to_host() as host:
        for index in range(pairs):
            results.append(Int(host[index]))
    return results^


def fill_pair(
    pair: Int,
    patterns: List[String],
    texts: List[String],
    rows_at: MutPointer[UInt64, _],
    text_at: MutPointer[UInt8, _],
    row_offsets: List[Int64],
    text_offsets: List[Int64],
    words_of: List[Int32],
):
    """One pair's match rows and its text's codes, in their places."""
    var pattern = patterns[pair].as_bytes()
    var words = Int(words_of[pair])
    var start = Int(row_offsets[pair])
    # Each symbol the pattern holds a code, in the order met; every other symbol the last code.
    var codes = Array[UInt8, 256](fill=0xFF)
    var count = 0
    for letter in pattern:
        if codes[Int(letter)] == 0xFF:
            codes[Int(letter)] = UInt8(count)
            count += 1
    for index in range(len(pattern)):
        var at = start + Int(codes[Int(pattern[index])]) * words + index // 64
        rows_at[unsafe_offset=at] = rows_at[unsafe_offset=at] | (UInt64(1) << UInt64(index % 64))
    var position = Int(text_offsets[pair])
    for letter in texts[pair].as_bytes():
        var code = codes[Int(letter)]
        text_at[unsafe_offset=position] = UInt8(count) if code == 0xFF else code
        position += 1
