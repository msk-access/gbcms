"""End-to-end test that the engine's wrong-length evidence reaches the user.

The wrong-length rule classifies a read carrying a pure indel of a different
length at the variant anchor as ``neither`` + partial evidence (a distinct
slippage allele in the same tract). This test proves the whole reporting
chain works, not just the counting: engine ``partial_alt`` > ``ad`` →
``_compute_diagnostics`` → ``ZERO_ALT;PARTIAL_DOMINANT`` in the MAF's
``gbcms_diagnostic`` column. Without that chain a tract locus whose support
is entirely a different allele would read as an unremarkable zero.

Scenario (the clinical headline shape): an A(10) homopolymer, an annotated
2bp deletion, and a BAM whose only non-reference reads carry a 1bp deletion —
a coexisting allele the annotation does not describe.

Harness conventions follow ``test_e2e_dna_maf.py``: FASTA/VCF use ``chr1``,
the BAM header uses the normalized name ``1``, reads are MAPQ 60 / Q30 so
default gates pass. Everything is synthetic and lives in ``tmp_path``.
"""

import glob
import random

import pysam
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

BAM_CONTIG = "1"  # normalized name; the pipeline strips 'chr' before fetch

ANCHOR_0BASED = 199  # 'C' anchor; the A(10) tract spans 200-209
EXPECTED_DEL_LEN = 2  # annotated CAA>C
READ_DEL_LEN = 1  # what the reads actually carry
N_REF_READS = 8
N_WRONG_LEN_READS = 6
READ_LEN = 100


def _reference_sequence():
    """600bp seeded-random contig with a pinned A(10) tract at 200-209.

    The anchor (199) and right flank (210) are pinned to non-A bases so the
    tract length is exactly 10 regardless of the random draw.
    """
    rng = random.Random(1)
    seq = [rng.choice("ACGT") for _ in range(600)]
    seq[199] = "C"
    seq[200:210] = list("A" * 10)
    seq[210] = "G"
    return "".join(seq)


def _build_reference(tmp_path, ref):
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">chr1\n" + ref + "\n")
    pysam.faidx(str(fasta))
    return fasta


def _build_reads_bam(tmp_path, ref):
    """8 pure-reference reads + 6 reads carrying a 1bp deletion at the tract."""
    reads = []
    for i in range(N_REF_READS):
        s = 150 + (i % 5)
        reads.append(make_read(f"ref_{i}", ref[s : s + READ_LEN], s, ((0, READ_LEN),)))
    p0 = ANCHOR_0BASED + 1  # first deleted base of the read-level D(1)
    for i in range(N_WRONG_LEN_READS):
        s = p0 - 30 - (i % 5)
        left = p0 - s
        seq = ref[s:p0] + ref[p0 + READ_DEL_LEN : p0 + READ_DEL_LEN + (READ_LEN - left)]
        reads.append(
            make_read(f"wl_{i}", seq, s, ((0, left), (2, READ_DEL_LEN), (0, READ_LEN - left)))
        )

    unsorted = tmp_path / "sample.unsorted.bam"
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 600, "SN": BAM_CONTIG}]}
    with pysam.AlignmentFile(unsorted, "wb", header=header) as outf:
        for r in sorted(reads, key=lambda a: a.reference_start):
            outf.write(r)
    sorted_bam = tmp_path / "sample.sorted.bam"
    pysam.sort("-o", str(sorted_bam), str(unsorted))
    pysam.index(str(sorted_bam))
    return str(sorted_bam)


def _build_vcf(tmp_path, ref):
    """One 2bp deletion: 1-based POS = anchor+1, REF = anchor + deleted bases."""
    vcf = tmp_path / "variants.vcf"
    pos_1based = ANCHOR_0BASED + 1
    ref_allele = ref[ANCHOR_0BASED : ANCHOR_0BASED + 1 + EXPECTED_DEL_LEN]
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=600>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t{pos_1based}\t.\t{ref_allele}\t{ref_allele[0]}\t.\t.\t.\n"
    )
    return vcf


def test_wrong_length_locus_reports_partial_dominant(tmp_path):
    ref = _reference_sequence()
    fasta = _build_reference(tmp_path, ref)
    bam = _build_reads_bam(tmp_path, ref)
    vcf = _build_vcf(tmp_path, ref)
    outdir = tmp_path / "out"
    outdir.mkdir()

    result = runner.invoke(
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
            str(outdir),
            "--format",
            "maf",
        ],
    )
    assert result.exit_code == 0, result.output

    maf_files = glob.glob(str(outdir / "*.maf"))
    assert len(maf_files) == 1, f"Expected exactly one MAF, found: {maf_files}"
    rows = list(read_maf_output(maf_files[0]))
    assert len(rows) == 1
    row = rows[0]

    ref_count = int(row["ref_count"])
    alt_count = int(row["alt_count"])
    partial_alt = int(row["partial_alt"])
    any_alt = int(row["any_alt"])
    total_count = int(row["total_count"])

    # The wrong-length reads are a distinct allele: never RD, never AD.
    assert ref_count == N_REF_READS
    assert alt_count == 0
    assert partial_alt == N_WRONG_LEN_READS
    assert any_alt == alt_count + partial_alt
    assert total_count >= ref_count + alt_count

    # House counting invariants (AGENTS.md): strand splits sum, and
    # fragment-level DP includes discarded fragments.
    assert int(row["ref_count_forward"]) + int(row["ref_count_reverse"]) == ref_count
    assert int(row["alt_count_forward"]) + int(row["alt_count_reverse"]) == alt_count
    assert int(row["total_count_fragment"]) >= int(row["ref_count_fragment"]) + int(
        row["alt_count_fragment"]
    )

    # The reporting chain: engine partial dominance surfaces as diagnostics.
    diagnostic = row["gbcms_diagnostic"]
    assert "PARTIAL_DOMINANT" in diagnostic, diagnostic
    assert "ZERO_ALT" in diagnostic, diagnostic
