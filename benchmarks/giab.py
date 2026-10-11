# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0. If a copy of the
# MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Real Illumina and PacBio HiFi pairs, the two technologies QuickEd's evaluation adds to A*PA2's nanopore
reads: reads of HG002 from Genome in a Bottle (Zook et al., Scientific Data 2016), each against the stretch of
HG002's own genome it came from, the T2T HG002 Q100 assembly v1.1 (Hansen et al., 2025).

Against its own genome a read differs only by its sequencing errors; against GRCh38 it would also carry
HG002's variants and whatever the reference's biases add. So:

1. The reads: one region of chromosome 20 streamed from GIAB's indexed GRCh38 alignments over HTTPS, the
   primary alignments of mapping quality 20 or more, a region whose coverage is ordinary (the one first
   tried, chr20:31.0 Mbp, beside the centromere, held 470 reads in 25 kbp and 107,000 in the next 75). That
   alignment only picks reads from a unique stretch of the genome; it decides no pair.
2. The genome: chromosome 20 of both HG002 haplotypes, `chr20_MATERNAL` and `chr20_PATERNAL`, from the
   assembly at `DINARA_HG002` if set, else cut from the published one over HTTPS.
3. The pairs: each read mapped to the two haplotypes by minimap2 (`-x sr` for Illumina, `-x map-hifi` for
   HiFi), its primary alignment kept whatever its mapping quality, which is zero wherever the haplotypes
   agree; its soft-clipped ends trimmed, and the haplotype's stretch its alignment covers cut out. A read
   whose alignment leaves more than a tenth of it clipped, or holds a letter past ACGT, is left out.

Each pair is written as a `.seq` pair, the genome's stretch `>` and the read `<`; the files are built once and
kept.

    python3 benchmarks/giab.py      # what bench.py does the first time: a minute or two
"""

import os
import subprocess
import sys
from pathlib import Path

GIAB = "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son"
HG002 = "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz"
CONTIGS = ["chr20_MATERNAL", "chr20_PATERNAL"]

DATASETS = {
    # 2x250 bp reads, novoalign against GRCh38 (GIAB's index AJtrio_Illumina_2x250bps_novoalign_GRCh37_GRCh38):
    # about 3,000 reads in each 10 kbp of this region.
    "illumina-250": (f"{GIAB}/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam", "chr20:40000000-40020000",
                     "sr"),
    # Sequel II HiFi reads of 15 to 20 kbp, pbmm2 against GRCh38 (index AJtrio_PacBio_CCS_15kb_20kb_chemistry2):
    # about 900 reads in these 200 kbp.
    "pacbio-hifi": (
        f"{GIAB}/PacBio_CCS_15kb_20kb_chemistry2/GRCh38/HG002.SequelII.merged_15kb_20kb.pbmm2.GRCh38.haplotag.10x.bam",
        "chr20:40000000-40200000",
        "map-hifi",
    ),
}

EXCLUDED_FLAGS = 0xF04
"""Unmapped, secondary, failing quality checks, duplicate and supplementary records: the primary alignment alone."""
MIN_MAPQ = 20
"""The GRCh38 mapping quality a read must have to be picked, so it comes from a unique stretch of the genome."""
MAX_CLIPPED = 0.1
"""The most of a read its alignment to HG002 may leave soft-clipped."""


def tool(name: str) -> list[str]:
    """A bioconda tool through pixi, as the benchmark's own environment carries none."""
    return [str(Path.home() / ".pixi/bin/pixi"), "exec", "-c", "conda-forge", "-c", "bioconda", "--spec", name, name]


def genome(directory: Path) -> tuple[Path, dict[str, str]]:
    """HG002's chromosome 20, both haplotypes, as a FASTA file and as upper-case sequences by name."""
    path = directory / "hg002v1.1_chr20.fa"
    if not path.exists():
        source = os.environ.get("DINARA_HG002", HG002)
        print(f"cutting chr20 of both haplotypes from {source} ...", file=sys.stderr, flush=True)
        partial = path.with_suffix(".part")
        with open(partial, "w") as out:
            subprocess.run(tool("samtools") + ["faidx", source, *CONTIGS], check=True, stdout=out)
        partial.rename(path)
    sequences, name, parts = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                sequences[name] = "".join(parts).upper()
            name, parts = line[1:].split()[0], []
        else:
            parts.append(line.strip())
    sequences[name] = "".join(parts).upper()
    return path, sequences


def operations(cigar: str) -> list[tuple[int, str]]:
    """A CIGAR's operations, each its count and letter."""
    found, number = [], ""
    for char in cigar:
        if char.isdigit():
            number += char
        else:
            found.append((int(number), char))
            number = ""
    return found


def reference_span(cigar: str) -> tuple[int, int, int]:
    """The reference letters an alignment covers, and the read letters soft-clipped at its start and end."""
    # Hard clips are gone from the read already; a soft clip sits next to them at either end.
    kept = [(count, op) for count, op in operations(cigar) if op != "H"]
    span = sum(count for count, op in kept if op in "MDN=X")
    head = kept[0][0] if kept and kept[0][1] == "S" else 0
    tail = kept[-1][0] if len(kept) > 1 and kept[-1][1] == "S" else 0
    return span, head, tail


def build(directory: Path) -> list[tuple[str, list[Path]]]:
    """Every GIAB dataset by name and its `.seq` file, built the first time."""
    directory.mkdir(parents=True, exist_ok=True)
    found = []
    for name, (bam, region, preset) in DATASETS.items():
        path = directory / f"{name}.seq"
        if not path.exists():
            reference, haplotypes = genome(directory)
            print(f"building {name}: GIAB's reads of {region} against HG002 ...", file=sys.stderr, flush=True)
            reads = directory / f"{name}.fastq"
            picked = subprocess.run(tool("samtools") + ["view", "-u", "-F", str(EXCLUDED_FLAGS), "-q", str(MIN_MAPQ), bam,
                                                        region], check=True, capture_output=True).stdout
            with open(reads, "wb") as out:
                subprocess.run(tool("samtools") + ["fastq", "-"], input=picked, check=True, stdout=out,
                               stderr=subprocess.DEVNULL)
            mapped = subprocess.run(tool("minimap2") + ["-a", "-x", preset, "--secondary=no", "-t", "4", str(reference),
                                                        str(reads)], check=True, capture_output=True, text=True).stdout
            pairs, skipped = [], 0
            for line in mapped.splitlines():
                if line.startswith("@"):
                    continue
                fields = line.split("\t")
                flags = int(fields[1])
                if flags & 0x904:
                    continue
                contig, start, cigar, read = fields[2], int(fields[3]) - 1, fields[5], fields[9]
                span, head, tail = reference_span(cigar)
                if head + tail > MAX_CLIPPED * len(read):
                    skipped += 1
                    continue
                read = read[head:len(read) - tail]
                segment = haplotypes[contig][start:start + span]
                if not read or len(segment) != span or set(read + segment) - set("ACGT"):
                    skipped += 1
                    continue
                pairs.append(f">{segment}\n<{read}\n")
            partial = path.with_suffix(".part")
            partial.write_text("".join(pairs))
            partial.rename(path)
            reads.unlink()
            print(f"{name}: {len(pairs)} pairs, {skipped} left out", file=sys.stderr, flush=True)
        found.append((name, [path]))
    return found


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).parent))
    from harness import DATA
    build(DATA / "giab")
