#!/usr/bin/env python3
"""A small counting run on a synthetic data set, for checking an installed gbcms.

It builds a reference, two paired-end BAMs and a variant file (an SNV, an MNV, an
insertion and a deletion) with a fixed pseudo-random generator, so every machine writes
the same inputs. It runs ``gbcms dna`` with mFSD on (its statistics call the C maths
library, so they show a difference between C libraries), checks the counts of the
variants the data set was built for, and prints a SHA-256 for each output file.

Two installs that print the same lines counted the same bytes. The lines left out of the
hash are the ones that hold paths (the command, the reference): they differ by where the
run happened, not by what it counted.

Usage::

    python scripts/wheel_fixture_run.py WORKDIR   # WORKDIR must not exist or be empty

Exit status: 0 when the counts are as built and the outputs exist; 1 otherwise.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import subprocess
import sys
from pathlib import Path

import pysam

__all__ = ["main"]

logger = logging.getLogger("wheel_fixture_run")

CONTIG = "1"
LENGTH = 1200
READ_LEN = 100
BASES = "ACGT"

# 0-based anchor position, kind, and the bases involved.
SNV_POS, MNV_POS, INS_POS, DEL_POS = 300, 450, 600, 800
INSERTED = "GT"
DELETED_LEN = 2

# (reference fragments, alternate fragments) built for each variant, per sample.
SAMPLES = {"tumor": (60, 14), "normal": (40, 3)}

# Output lines that hold paths (they differ by where the run happened).
PATH_LINE_PREFIXES = ("#command", "##gbcms_command", "##reference")


class Lcg:
    """A 64-bit linear congruential generator: the same stream on every Python."""

    def __init__(self, seed: int) -> None:
        self.state = seed & (2**64 - 1)

    def below(self, n: int) -> int:
        """Return the next integer in ``[0, n)``."""
        self.state = (self.state * 6364136223846793005 + 1442695040888963407) & (2**64 - 1)
        return (self.state >> 33) % n

    def between(self, low: int, high: int) -> int:
        """Return the next integer in ``[low, high]``."""
        return low + self.below(high - low + 1)


def make_reference(rng: Lcg) -> str:
    """A reference without runs of one base, so no indel can slide when it is normalised."""
    seq: list[str] = []
    while len(seq) < LENGTH:
        base = BASES[rng.below(4)]
        if len(seq) >= 2 and seq[-1] == base and seq[-2] == base:
            continue
        seq.append(base)
    return "".join(seq)


def other_base(base: str, shift: int) -> str:
    """A base that differs from ``base``."""
    return BASES[(BASES.index(base) + shift) % 4]


def variants(ref: str) -> list[dict]:
    """The four variants: 0-based anchor, VCF REF and ALT, and how a read carries ALT."""
    snv_alt = other_base(ref[SNV_POS], 1)
    mnv_ref = ref[MNV_POS : MNV_POS + 2]
    mnv_alt = other_base(mnv_ref[0], 2) + other_base(mnv_ref[1], 3)
    return [
        {"kind": "snv", "pos": SNV_POS, "ref": ref[SNV_POS], "alt": snv_alt},
        {"kind": "mnv", "pos": MNV_POS, "ref": mnv_ref, "alt": mnv_alt},
        {"kind": "ins", "pos": INS_POS, "ref": ref[INS_POS], "alt": ref[INS_POS] + INSERTED},
        {
            "kind": "del",
            "pos": DEL_POS,
            "ref": ref[DEL_POS : DEL_POS + 1 + DELETED_LEN],
            "alt": ref[DEL_POS],
        },
    ]


def read_for(ref: str, var: dict, start: int, alt: bool) -> tuple[str, list[tuple[int, int]]]:
    """Return the bases and CIGAR of a read of READ_LEN bases starting at ``start``.

    The read carries ``var``'s ALT when ``alt`` is true, else the reference.
    """
    pos, kind = var["pos"], var["kind"]
    idx = pos - start  # the read index of the anchor base
    if not alt:
        return ref[start : start + READ_LEN], [(0, READ_LEN)]
    if kind == "snv":
        seq = ref[start : start + READ_LEN]
        return seq[:idx] + var["alt"] + seq[idx + 1 :], [(0, READ_LEN)]
    if kind == "mnv":
        seq = ref[start : start + READ_LEN]
        return seq[:idx] + var["alt"] + seq[idx + 2 :], [(0, READ_LEN)]
    if kind == "ins":
        head = ref[start : pos + 1] + INSERTED
        seq = head + ref[pos + 1 : pos + 1 + READ_LEN - len(head)]
        return seq, [(0, idx + 1), (1, len(INSERTED)), (0, READ_LEN - idx - 1 - len(INSERTED))]
    head = ref[start : pos + 1]
    seq = head + ref[pos + 1 + DELETED_LEN : pos + 1 + DELETED_LEN + READ_LEN - len(head)]
    return seq, [(0, idx + 1), (2, DELETED_LEN), (0, READ_LEN - idx - 1)]


def fragments(ref: str, vars_: list[dict], sample: str, rng: Lcg) -> list[pysam.AlignedSegment]:
    """Paired reads over every variant: reference fragments longer than alternate ones."""
    n_ref, n_alt = SAMPLES[sample]
    reads: list[pysam.AlignedSegment] = []
    for var in vars_:
        for k in range(n_ref + n_alt):
            alt = k >= n_ref
            size = rng.between(120, 200) if alt else rng.between(150, 260)
            start1 = var["pos"] - rng.between(12, 70)
            start2 = start1 + size - READ_LEN
            if start1 < 0 or start2 + READ_LEN > LENGTH:
                continue
            name = f"{sample}_{var['kind']}_{'alt' if alt else 'ref'}_{k}"
            seq1, cigar1 = read_for(ref, var, start1, alt)
            seq2 = ref[start2 : start2 + READ_LEN]
            for mate, (seq, cigar, start, flag, tlen) in enumerate(
                [(seq1, cigar1, start1, 99, size), (seq2, [(0, READ_LEN)], start2, 147, -size)]
            ):
                read = pysam.AlignedSegment()
                read.query_name = name
                read.query_sequence = seq
                read.flag = flag
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = 60
                read.cigartuples = cigar
                read.query_qualities = pysam.qualitystring_to_array(
                    "".join(chr(33 + rng.between(28, 38)) for _ in seq)
                )
                read.next_reference_id = 0
                read.next_reference_start = start2 if mate == 0 else start1
                read.template_length = tlen
                reads.append(read)
    return reads


def write_bam(path: Path, reads: list[pysam.AlignedSegment]) -> None:
    """Write ``reads`` as a coordinate-sorted, indexed BAM."""
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"LN": LENGTH, "SN": CONTIG}]}
    unsorted = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(unsorted), "wb", header=header) as out:
        for read in reads:
            out.write(read)
    pysam.sort("-o", str(path), str(unsorted))
    pysam.index(str(path))
    unsorted.unlink()


def build_inputs(work: Path) -> tuple[Path, Path, list[Path], list[dict]]:
    """Write the reference, the BAMs and the variant file; return their paths."""
    rng = Lcg(20261030)
    ref = make_reference(rng)
    vars_ = variants(ref)
    indir = work / "in"
    indir.mkdir(parents=True)
    fasta = indir / "ref.fa"
    fasta.write_text(f">{CONTIG}\n{ref}\n")
    pysam.faidx(str(fasta))
    bams = []
    for sample in SAMPLES:
        bam = indir / f"{sample}.bam"
        write_bam(bam, fragments(ref, vars_, sample, rng))
        bams.append(bam)
    vcf = indir / "variants.vcf"
    lines = [
        "##fileformat=VCFv4.2",
        f"##contig=<ID={CONTIG},length={LENGTH}>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    lines += [f"{CONTIG}\t{v['pos'] + 1}\t.\t{v['ref']}\t{v['alt']}\t.\t.\t." for v in vars_]
    vcf.write_text("\n".join(lines) + "\n")
    return fasta, vcf, bams, vars_


def run_gbcms(fasta: Path, vcf: Path, bams: list[Path], out: Path) -> None:
    """Run ``gbcms dna`` from the interpreter's own environment."""
    exe = Path(sys.executable).parent / "gbcms"
    cmd = [str(exe), "dna", "-v", str(vcf), "-f", str(fasta), "-o", str(out)]
    cmd += ["--format", "maf", "--mfsd", "--threads", "2"]
    for bam in bams:
        cmd += ["-b", str(bam)]
    logger.info("running: %s", " ".join(cmd))
    done = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if done.returncode != 0:
        sys.stderr.write(done.stdout + done.stderr)
        raise SystemExit(f"gbcms exited with status {done.returncode}")


def check_counts(out: Path, vars_: list[dict]) -> None:
    """Every variant of the tumour sample counts the alternate fragments built, no fewer
    reference fragments than built (reads of the other variants' fragments cover the site
    too), and a depth that holds both."""
    n_ref, n_alt = SAMPLES["tumor"]
    maf = next(out.glob("tumor*.maf"))
    with open(maf, newline="") as fh:
        rows = list(csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"))
    by_pos = {int(r["Start_Position"]): r for r in rows}
    # MAF puts an insertion at its anchor base and a deletion at its first deleted base
    shift = {"snv": 1, "mnv": 1, "ins": 1, "del": 2}
    for var in vars_:
        row = by_pos.get(var["pos"] + shift[var["kind"]])
        if row is None:
            raise SystemExit(f"no output row for the {var['kind']} (rows at {sorted(by_pos)})")
        ref, alt, depth = int(row["ref_count"]), int(row["alt_count"]), int(row["total_count"])
        if alt != n_alt or ref < n_ref or depth < ref + alt:
            raise SystemExit(
                f"{var['kind']}: counted ref {ref}, alt {alt}, depth {depth}; "
                f"built ref {n_ref}, alt {n_alt}"
            )
        logger.info("%s: ref %d alt %d depth %d", var["kind"], ref, alt, depth)


def digests(out: Path) -> list[str]:
    """SHA-256 of each output file, the lines that hold paths left out."""
    lines = []
    for path in sorted(out.iterdir()):
        if not path.is_file():
            continue
        body = [
            ln
            for ln in path.read_text().splitlines(keepends=True)
            if not ln.startswith(PATH_LINE_PREFIXES)
        ]
        lines.append(f"{hashlib.sha256(''.join(body).encode()).hexdigest()}  {path.name}")
    return lines


def main(argv: list[str] | None = None) -> int:
    """Build the inputs, run gbcms, check the counts, print the digests."""
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("workdir", type=Path, help="a folder that does not exist or is empty")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s", stream=sys.stderr)
    work: Path = args.workdir
    if work.exists() and any(work.iterdir()):
        raise SystemExit(f"{work} is not empty")
    fasta, vcf, bams, vars_ = build_inputs(work)
    out = work / "out"
    out.mkdir()
    run_gbcms(fasta, vcf, bams, out)
    check_counts(out, vars_)
    found = digests(out)
    if not found:
        raise SystemExit("gbcms wrote no output file")
    print("\n".join(found))
    return 0


if __name__ == "__main__":
    sys.exit(main())
