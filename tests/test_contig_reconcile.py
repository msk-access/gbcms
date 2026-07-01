"""Contig-naming reconciliation (regression guard for the silent zero-count bug).

A BAM whose contigs are named differently from the variants — UCSC ``chr1`` vs
Ensembl/b37 ``1`` — must still find reads, not silently return zero counts and exit 0.
Before the fix, the binning step looked up the (normalized) variant chrom in the BAM
header verbatim; a ``chr1``-named BAM with a ``1`` lookup found no tid, built 0 bins,
and every count came back 0 while the run exited 0. ``resolve_tid`` now reconciles the
naming via ``normalize_contig``.
"""

import glob

import pysam
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()


def _run_with_bam_contig(tmp_path, bam_contig: str):
    """Reference + VCF use ``chr1``; the BAM header contig is ``bam_contig``.

    Returns (ref_count, alt_count, total_count) from the emitted MAF row.
    """
    # Reference: 500 bp chr1, base 'A' at 0-based pos 100 (the variant site).
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">chr1\n" + "A" * 500 + "\n")
    pysam.faidx(str(fasta))

    # 5 REF ('A' at read index 5) + 3 ALT ('T'), 10M reads starting at 95 → cover pos 100.
    reads = [make_read(f"ref_{i}", "AAAAAAAAAA", start=95, cigar=((0, 10),)) for i in range(5)]
    reads += [make_read(f"alt_{i}", "AAAAATAAAA", start=95, cigar=((0, 10),)) for i in range(3)]
    unsorted = tmp_path / "s.unsorted.bam"
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 500, "SN": bam_contig}]}
    with pysam.AlignmentFile(unsorted, "wb", header=header) as outf:
        for r in reads:
            outf.write(r)
    bam = tmp_path / "s.sorted.bam"
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))

    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=500>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t101\t.\tA\tT\t.\t.\t.\n"
    )

    out = tmp_path / "out"
    out.mkdir()
    res = runner.invoke(
        app,
        [
            "dna",
            "-v",
            str(vcf),
            "-b",
            str(bam),
            "-f",
            str(fasta),
            "-o",
            str(out),
            "--format",
            "maf",
        ],
    )
    assert res.exit_code == 0, res.output
    row = next(read_maf_output(glob.glob(str(out / "*.maf"))[0]))
    return int(row["ref_count"]), int(row["alt_count"]), int(row["total_count"])


def test_ucsc_named_bam_reconciles_to_variants(tmp_path):
    """A ``chr1``-named BAM finds reads for the variant (was 0 counts before the fix)."""
    assert _run_with_bam_contig(tmp_path, "chr1") == (5, 3, 8)


def test_ensembl_named_bam_unaffected(tmp_path):
    """The already-matching case (``1``-named BAM) is unchanged."""
    assert _run_with_bam_contig(tmp_path, "1") == (5, 3, 8)
