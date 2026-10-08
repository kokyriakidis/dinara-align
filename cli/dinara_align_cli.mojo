# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""
The command line of dinara-align: pairs of sequences in, their exact alignments out, as a table, SAM or
PAF.

    pixi run build-cli                                   # build/cli: the program and its runtime libraries
    build/cli/dinara-align reference.fa reads.fq      # every read against the one reference
    build/cli/dinara-align -r ACGTACGT -q ACGTTCGT    # one pair given whole
    build/cli/dinara-align --pairs pairs.tsv          # name, reference, query a line

Run with `--help` for every option (see `USAGE`).
"""

from std.sys import argv, exit, stderr

from dinara_align.common import hardware_threads

from dinara_align import (
    Alignment,
    Anchor,
    Band,
    Costs,
    DEFAULT_MAX_MEMORY,
    Mode,
    Ties,
    alignments,
    distances,
)

comptime USAGE = """usage: dinara-align [options] REFERENCE QUERY
       dinara-align [options] -r SEQUENCE -q SEQUENCE
       dinara-align [options] --pairs PAIRS.tsv

Aligns each query to its reference exactly and prints one line per pair.

Inputs:
  REFERENCE QUERY      FASTA or FASTQ files: the i-th query against the i-th reference, or every query
                       against the reference when the reference file holds one sequence
  -r, -q SEQUENCE      one reference and one query, given whole
  --pairs FILE         a file of `name<TAB>reference<TAB>query` lines

Costs (an alignment minimizes them):
  --costs edit                     unit costs, the edit distance (the default)
  --costs linear:X,G               a substitution X, every gapped letter G
  --costs affine:X,O,E             a gap of k letters O + kE (WFA2-lib's default is affine:4,6,2)
  --costs two-piece:X,O,E,O2,E2    a gap the less of O + kE and O2 + kE2 (minimap2's -O, -E)
  --deletions O,E[,O2,E2]          deletions, reference letters alone, priced apart (bwa's -O del,ins)

Modes:
  --mode global                    both sequences end to end (the default)
  --mode infix                     the whole query inside the reference
  --mode prefix | suffix           the whole query against a prefix, or a suffix, of the reference
  --mode reference-in-query        the whole reference inside the query
  --mode ends:RS,RE,QS,QE          so many letters free at each end of each
  --mode local:A                   Smith-Waterman, a match earning A
  --mode overlap:A                 every end gap free, a match earning A
  --mode extension:A[,end]         fixed at the start (or end) of both, a match earning A
  --match-score A                  free ends rewarded: the best score rather than the least cost
  --zdrop Z                        an extension gives up past a fall of Z (KSW2's Z-drop)
  --end-bonus B                    an extension reaches the read's end when within B of its best (KSW2's)

Options:
  --band W | --band LOW,HIGH       every move within diagonals -W..=W, or LOW..=HIGH
  --max-cost C                     pairs costing more than C reported as unaligned
  --ties left | right              equally good alignments: indels left (minimap2) or right (WFA2-lib)
  --cigar-m                        M for matches and substitutions alike, not = and X
  --both-strands                   also align each query's reverse complement, keep the better
  --distance                       the cost alone, no alignment
  --format tsv | sam | paf         the output (tsv by default)
  --threads N                      threads, every one this process may use by default
  --max-memory BYTES               the fronts kept for a traceback, about 80 MB by default
  -h, --help                       this text
"""


@fieldwise_init
struct Record(Copyable, Movable):
    """A FASTA or FASTQ record: its name, up to the first space, and its sequence."""

    var name: String
    var sequence: String


def fail(message: String):
    """Prints `message` to standard error and exits with status 2, as for a bad command line."""
    print("dinara-align:", message, file=stderr)
    exit(2)


def read_records(path: String) raises -> List[Record]:
    """A FASTA or FASTQ file's records, each sequence upper case, whitespace dropped."""
    var text = open(path, "r").read()
    var lines = text.split("\n")
    var records = List[Record]()
    var index = 0
    while index < len(lines):
        var line = String(lines[index].strip())
        index += 1
        if line.byte_length() == 0:
            continue
        if line.startswith(">"):
            var name = String(line.removeprefix(">").split(" ")[0])
            var sequence = String()
            while index < len(lines) and not String(lines[index]).startswith(">"):
                sequence += String(lines[index].strip())
                index += 1
            records.append(Record(name, sequence.upper()))
        elif line.startswith("@"):
            var name = String(line.removeprefix("@").split(" ")[0])
            if index >= len(lines):
                raise Error(String("a FASTQ record without its sequence in ", path))
            var sequence = String(lines[index].strip()).upper()
            # The sequence line, its `+` line and the qualities.
            index += 3
            records.append(Record(name, sequence))
        else:
            raise Error(String(path, ": neither FASTA nor FASTQ at `", line, "`"))
    return records^


def numbers(text: String) raises -> List[Int]:
    """The integers of a comma-separated list, as `--band` and `--costs` take them."""
    var out = List[Int]()
    for part in text.split(","):
        out.append(Int(String(part)))
    return out^


def costs_of(text: String) raises -> Costs:
    """The `Costs` a `--costs` value names: `edit`, `linear:X,G`, `affine:X,O,E` or `two-piece:X,O,E,O2,E2`."""
    var kind = String(text.split(":")[0])
    var values = numbers(String(text.split(":")[1])) if ":" in text else List[Int]()
    if kind == "edit":
        return Costs.edit()
    if kind == "linear" and len(values) == 2:
        return Costs.linear(values[0], values[1])
    if kind == "affine" and len(values) == 3:
        return Costs.affine(values[0], values[1], values[2])
    if kind == "two-piece" and len(values) == 5:
        return Costs.two_piece(values[0], values[1], values[2], values[3], values[4])
    raise Error(String("--costs ", text, ": edit, linear:X,G, affine:X,O,E or two-piece:X,O,E,O2,E2"))


def mode_of(text: String, match_score: Int, zdrop: Int, end_bonus: Int) raises -> Mode:
    """The `Mode` a `--mode` value names; free ends take `--match-score` when above zero, an extension
    `--zdrop` and `--end-bonus` when zero or more, and the other modes their reward from the value itself."""
    var kind = String(text.split(":")[0])
    var values = numbers(String(text.split(":")[1])) if ":" in text and not text.endswith(",end") else List[Int]()
    var mode: Mode
    if kind == "global":
        mode = Mode.GLOBAL
    elif kind == "infix":
        mode = Mode.INFIX
    elif kind == "prefix":
        mode = Mode.PREFIX
    elif kind == "suffix":
        mode = Mode.SUFFIX
    elif kind == "reference-in-query":
        mode = Mode.REFERENCE_IN_QUERY
    elif kind == "ends" and len(values) == 4:
        mode = Mode.ends_free(
            reference_start=values[0], reference_end=values[1], query_start=values[2], query_end=values[3]
        )
    elif kind == "local" and len(values) == 1:
        return Mode.local(values[0])
    elif kind == "overlap" and len(values) == 1:
        return Mode.overlap(values[0])
    elif kind == "extension":
        var reward = Int(String(String(text.split(":")[1]).split(",")[0]))
        var anchor = Anchor.END if text.endswith(",end") else Anchor.START
        var drop = Optional[Int](zdrop) if zdrop >= 0 else None
        var bonus = Optional[Int](end_bonus) if end_bonus >= 0 else None
        return Mode.extension(reward, anchor, zdrop=drop, end_bonus=bonus)
    else:
        raise Error(String("--mode ", text, ": see --help"))
    return mode.with_match_score(match_score) if match_score > 0 else mode


def reverse_complement(sequence: String) -> String:
    """`sequence` read backward with A and T, C and G swapped; any other letter is kept as it is."""
    var bytes = sequence.as_bytes()
    var out = List[UInt8](capacity=len(bytes))
    for index in range(len(bytes) - 1, -1, -1):
        var letter = bytes[index]
        if letter == UInt8(ord("A")):
            letter = UInt8(ord("T"))
        elif letter == UInt8(ord("T")):
            letter = UInt8(ord("A"))
        elif letter == UInt8(ord("C")):
            letter = UInt8(ord("G"))
        elif letter == UInt8(ord("G")):
            letter = UInt8(ord("C"))
        out.append(letter)
    return String(unsafe_from_utf8=out^)


def better(first: Alignment, second: Alignment, scored: Bool) -> Bool:
    """Whether `second` beats `first`: a higher score for a mode with a reward, else a lower cost."""
    return second.score > first.score if scored else second.cost < first.cost


def main():
    """Runs the command line, a failure printed as one line with exit status 2."""
    try:
        run()
    except error:
        # The library's refusals and a bad file both come back as a line, not a stack trace.
        fail(String(String(error).removeprefix("dinara-align: ")))


def run() raises:
    """Reads the options and the pairs, aligns every pair and prints a line for each in the format asked."""
    var arguments = argv()
    var positional = List[String]()
    var literal_reference = Optional[String]()
    var literal_query = Optional[String]()
    var pairs_file = Optional[String]()
    var costs_text = String("edit")
    var deletions = Optional[String]()
    var mode_text = String("global")
    var match_score = 0
    var zdrop = -1
    var end_bonus = -1
    var band = Band()
    var max_cost = -1
    var ties = Ties.LEFT
    var eqx = True
    var both_strands = False
    var cost_only = False
    var format = String("tsv")
    var threads = Optional[Int]()
    var memory = DEFAULT_MAX_MEMORY
    var index = 1
    while index < len(arguments):
        var argument = String(arguments[index])
        var has_value = index + 1 < len(arguments)
        var value = String(arguments[index + 1]) if has_value else String()
        if argument == "-h" or argument == "--help":
            print(USAGE)
            return
        var takes_value = argument in [
            "-r",
            "-q",
            "--pairs",
            "--costs",
            "--deletions",
            "--mode",
            "--match-score",
            "--zdrop",
            "--end-bonus",
            "--band",
            "--max-cost",
            "--ties",
            "--format",
            "--threads",
            "--max-memory",
        ]
        if takes_value and not has_value:
            fail(String(argument, " needs a value"))
        if argument == "-r":
            literal_reference = value.upper()
        elif argument == "-q":
            literal_query = value.upper()
        elif argument == "--pairs":
            pairs_file = value
        elif argument == "--costs":
            costs_text = value
        elif argument == "--deletions":
            deletions = value
        elif argument == "--mode":
            mode_text = value
        elif argument == "--match-score":
            match_score = Int(value)
        elif argument == "--zdrop":
            zdrop = Int(value)
        elif argument == "--end-bonus":
            end_bonus = Int(value)
        elif argument == "--band":
            var edges = numbers(value)
            band = Band.around(edges[0]) if len(edges) == 1 else Band(edges[0], edges[1])
        elif argument == "--max-cost":
            max_cost = Int(value)
        elif argument == "--ties":
            if value != "left" and value != "right":
                fail("--ties: left or right")
            ties = Ties.RIGHT if value == "right" else Ties.LEFT
        elif argument == "--format":
            if value != "tsv" and value != "sam" and value != "paf":
                fail("--format: tsv, sam or paf")
            format = value
        elif argument == "--threads":
            threads = Int(value)
        elif argument == "--max-memory":
            memory = Int(value)
        elif argument == "--cigar-m":
            eqx = False
        elif argument == "--both-strands":
            both_strands = True
        elif argument == "--distance":
            cost_only = True
        elif argument.startswith("-"):
            fail(String("unknown option ", argument, "; see --help"))
        else:
            positional.append(argument)
        index += 2 if takes_value else 1

    var costs = costs_of(costs_text)
    if deletions:
        var values = numbers(deletions.value())
        costs = costs.with_deletions(values[0], values[1], values[2], values[3]) if len(
            values
        ) == 4 else costs.with_deletions(values[0], values[1])
    var mode = mode_of(mode_text, match_score, zdrop, end_bonus)

    # The pairs: names, references and queries.
    var names = List[String]()
    var reference_names = List[String]()
    var references = List[String]()
    var queries = List[String]()
    if literal_reference or literal_query:
        if not literal_reference or not literal_query:
            fail("-r and -q come together")
        names.append("query")
        reference_names.append("reference")
        references.append(literal_reference.value())
        queries.append(literal_query.value())
    elif pairs_file:
        for line in open(pairs_file.value(), "r").read().split("\n"):
            if line.byte_length() == 0:
                continue
            var fields = line.split("\t")
            if len(fields) < 3:
                fail(String(pairs_file.value(), ": a line without name, reference and query"))
            names.append(String(fields[0]))
            reference_names.append(String(fields[0]))
            references.append(String(fields[1]).upper())
            queries.append(String(fields[2]).upper())
    elif len(positional) == 2:
        var targets = read_records(positional[0])
        var reads = read_records(positional[1])
        if len(targets) != 1 and len(targets) != len(reads):
            fail("the reference file holds one sequence, or as many as the query file")
        for item in range(len(reads)):
            ref target = targets[0 if len(targets) == 1 else item]
            names.append(reads[item].name)
            reference_names.append(target.name)
            references.append(target.sequence)
            queries.append(reads[item].sequence)
    else:
        print(USAGE)
        exit(2)

    # The executable, not the library, spreads the work: every thread this process may use, unless asked.
    if not threads:
        threads = hardware_threads()
    if cost_only:
        var found = distances(
            references,
            queries,
            costs,
            mode,
            max_cost=max_cost if max_cost >= 0 else Int.MAX,
            band=band,
            threads=threads,
        )
        print("query\treference\tcost")
        for pair in range(len(found)):
            print(names[pair], reference_names[pair], String(found[pair].value()) if found[pair] else "*", sep="\t")
        return

    var cap = max_cost if max_cost >= 0 else Int.MAX
    var forward = alignments(
        references,
        queries,
        costs,
        mode,
        max_cost=cap,
        band=band,
        ties=ties,
        eqx=eqx,
        threads=threads,
        max_memory=memory,
    ) if not mode.is_scored() else _scored(references, queries, costs, mode, band, ties, eqx, threads, memory)
    var reverse = List[Optional[Alignment]]()
    var flipped = List[String]()
    if both_strands:
        for query in queries:
            flipped.append(reverse_complement(query))
        reverse = alignments(
            references,
            flipped,
            costs,
            mode,
            max_cost=cap,
            band=band,
            ties=ties,
            eqx=eqx,
            threads=threads,
            max_memory=memory,
        ) if not mode.is_scored() else _scored(references, flipped, costs, mode, band, ties, eqx, threads, memory)

    if format == "sam":
        print("@HD\tVN:1.6\tSO:unsorted")
        var seen = List[String]()
        for pair in range(len(references)):
            if reference_names[pair] not in seen:
                seen.append(reference_names[pair])
                print(String("@SQ\tSN:", reference_names[pair], "\tLN:", references[pair].byte_length()))
        print("@PG\tID:dinara-align\tPN:dinara-align")
    elif format == "tsv":
        print("query\treference\tstrand\tcost\tscore\treference_start\treference_end\tquery_start\tquery_end\tcigar")
    for pair in range(len(references)):
        var strand = "+"
        var found = forward[pair].copy()
        var query = queries[pair]
        if (
            both_strands
            and reverse[pair]
            and (not found or better(found.value(), reverse[pair].value(), mode.is_scored()))
        ):
            found = reverse[pair].copy()
            strand = "-"
            query = flipped[pair]
        var length = query.byte_length()
        if not found:
            if format == "sam":
                # Flag 4: the read is unmapped.
                print(names[pair], 4, "*", 0, 0, "*", "*", 0, 0, query, "*", sep="\t")
            elif format == "tsv":
                print(names[pair], reference_names[pair], "*", "*", "*", "*", "*", "*", "*", "*", sep="\t")
            continue
        ref hit = found.value()
        if format == "tsv":
            print(
                names[pair],
                reference_names[pair],
                strand,
                hit.cost,
                hit.score,
                hit.reference_start,
                hit.reference_end,
                hit.query_start,
                hit.query_end,
                hit.cigar,
                sep="\t",
            )
        elif format == "sam":
            # Flag 16 marks the reverse strand; a mapping quality of 255 is none given.
            var flag = 16 if strand == "-" else 0
            print(
                names[pair],
                flag,
                reference_names[pair],
                hit.reference_start + 1,
                255,
                hit.clipped_cigar(length) if hit.cigar.byte_length() > 0 else "*",
                "*",
                0,
                0,
                query,
                "*",
                String("NM:i:", hit.edit_distance(references[pair], query)),
                String("MD:Z:", hit.mismatch_string(references[pair], query)),
                String("AS:i:", hit.score),
                sep="\t",
            )
        else:
            var counted = hit.counts(references[pair], query)
            var columns = counted.matches + counted.mismatches + counted.deleted + counted.inserted
            # PAF counts the query's coordinates on the forward strand.
            var query_start = hit.query_start if strand == "+" else length - hit.query_end
            var query_end = hit.query_end if strand == "+" else length - hit.query_start
            print(
                names[pair],
                length,
                query_start,
                query_end,
                strand,
                reference_names[pair],
                references[pair].byte_length(),
                hit.reference_start,
                hit.reference_end,
                counted.matches,
                columns,
                255,
                String("NM:i:", hit.edit_distance(references[pair], query)),
                String("AS:i:", hit.score),
                String("cg:Z:", hit.cigar),
                sep="\t",
            )


def _scored(
    references: List[String],
    queries: List[String],
    costs: Costs,
    mode: Mode,
    band: Band,
    ties: Ties,
    eqx: Bool,
    threads: Optional[Int],
    memory: Int,
) raises -> List[Optional[Alignment]]:
    """Every pair's alignment under a mode with a match score, which takes no cap, so every pair is aligned."""
    var out = List[Optional[Alignment]]()
    for found in alignments(
        references, queries, costs, mode, band=band, ties=ties, eqx=eqx, threads=threads, max_memory=memory
    ):
        out.append(found.copy())
    return out^
