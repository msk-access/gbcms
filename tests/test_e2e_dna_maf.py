"""First real end-to-end test of the ``gbcms dna`` pipeline producing MAF output.

This is the repository's first *working* end-to-end test that drives the full
``gbcms dna`` CLI command against synthetic inputs and asserts real, non-zero
counts in the primary MSK output format (MAF).

Why it exists: the committed ``tests/testdata/integration_test_*`` fixtures have
a contig-name mismatch and reject every variant, so they never exercise a
successful count. This test builds everything from scratch in ``tmp_path`` — a
tiny reference FASTA, a synthetic BAM with known REF/ALT reads, and a one-SNP
VCF — runs the real CLI via ``CliRunner``, and checks the exact ref/alt/total
counts in the emitted MAF row.

Everything is synthetic and lives in ``tmp_path``; nothing is committed and no
real/patient data is touched.

Coordinate/count gotchas encountered and locked in by this test:
- VCF is 1-based; the engine works 0-based internally. The variant base is at
  0-based reference position 100, which is VCF ``POS=101``. The test passing is
  the confirmation of this convention.
- CONTIG NORMALIZATION (the exact bug that breaks the committed integration
  fixtures): the pipeline strips the ``chr`` prefix from variant chromosome names
  (``CoordinateKernel.normalize_chromosome``), so a ``chr1`` variant is looked up
  in the BAM header as ``1``. If the BAM header contig is ``chr1``, the binning
  step finds no matching tid, builds 0 genomic bins, fetches no reads, and every
  count comes back 0 (the run still exits 0). The FASTA/VCF therefore use ``chr1``
  (so REF validation against the FASTA succeeds) while the BAM header uses the
  normalized name ``1`` (so read fetch succeeds). This is why the test builds the
  BAM by hand instead of using the ``build_bam`` helper, whose header is ``chr1``.
- Reads are single-end, MAPQ 60, all bases Q30, and fully match the reference
  (10M, no soft-clips), so every read passes the default MAPQ/BQ gates and none
  is filtered. 5 REF ('A') + 3 ALT ('T') => ref_count=5, alt_count=3,
  total_count=8 with no 'neither' reads.
- MAF count columns are unprefixed by default: ``ref_count`` (=RD), ``alt_count``
  (=AD), ``total_count`` (=DP). See ``src/gbcms/io/output.py``. The MAF writer
  also emits the normalized contig name, so the Chromosome column reads ``1``.
"""

import glob

import pysam
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

# BAM header contig name. Must be the *normalized* name (no 'chr' prefix) because
# the pipeline strips 'chr' from variant chroms before fetching reads. See the
# module docstring's CONTIG NORMALIZATION note.
BAM_CONTIG = "1"

runner = CliRunner()

# 0-based reference position of the variant base. VCF POS is this + 1.
VARIANT_POS_0BASED = 100
REF_BASE = "A"
ALT_BASE = "T"
N_REF_READS = 5
N_ALT_READS = 3


def _build_reference(tmp_path):
    """Write a 500 bp single-contig (chr1) FASTA with a known REF base at pos 100.

    The contig is all 'A' so the REF base at position 100 is 'A' (matching the
    REF reads). Indexed with ``pysam.faidx`` so a ``.fai`` sits beside it.
    """
    length = 500
    seq = ["A"] * length
    seq[VARIANT_POS_0BASED] = REF_BASE  # explicit/known base at the variant site
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">chr1\n" + "".join(seq) + "\n")
    pysam.faidx(str(fasta))
    assert (tmp_path / "ref.fasta.fai").exists()
    return fasta


def _build_reads_bam(tmp_path):
    """Build a BAM with 5 REF and 3 ALT single-end reads spanning the variant.

    Each read is a 10M alignment starting at 0-based 95, placing the variant base
    at read index 5 (95 + 5 = 100). All reads are forward, MAPQ 60, Q30.

    The header contig is ``BAM_CONTIG`` (= "1", normalized). We do not use the
    ``build_bam`` helper because its header hardcodes ``chr1``, which would not
    match the pipeline's normalized variant chrom and would yield zero counts.
    """
    reads = []
    # Reference base is at query index 5 for a read starting at 95 (95+5 = 100).
    for i in range(N_REF_READS):
        seq = "AAAAA" + REF_BASE + "AAAA"
        reads.append(make_read(f"ref_{i}", seq, start=95, cigar=((0, 10),)))
    for i in range(N_ALT_READS):
        seq = "AAAAA" + ALT_BASE + "AAAA"
        reads.append(make_read(f"alt_{i}", seq, start=95, cigar=((0, 10),)))

    unsorted = tmp_path / "sample.unsorted.bam"
    header = {
        "HD": {"VN": "1.0", "SO": "coordinate"},
        "SQ": [{"LN": 500, "SN": BAM_CONTIG}],
    }
    with pysam.AlignmentFile(unsorted, "wb", header=header) as outf:
        for r in reads:
            outf.write(r)
    sorted_bam = tmp_path / "sample.sorted.bam"
    pysam.sort("-o", str(sorted_bam), str(unsorted))
    pysam.index(str(sorted_bam))
    return str(sorted_bam)


def _build_vcf(tmp_path):
    """Write a minimal one-SNP VCF. VCF is 1-based, so 0-based 100 => POS=101."""
    vcf = tmp_path / "variants.vcf"
    pos_1based = VARIANT_POS_0BASED + 1
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=500>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t{pos_1based}\t.\t{REF_BASE}\t{ALT_BASE}\t.\t.\t.\n"
    )
    return vcf


def test_e2e_dna_maf_counts(tmp_path):
    """Run the full ``gbcms dna`` pipeline and assert exact MAF counts."""
    fasta = _build_reference(tmp_path)
    bam = _build_reads_bam(tmp_path)
    vcf = _build_vcf(tmp_path)
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

    # Output is written as "{bam_stem}{suffix}.maf"; glob to avoid coupling to
    # the exact stem (the helper sorts to "sample.sorted.bam").
    maf_files = glob.glob(str(outdir / "*.maf"))
    assert len(maf_files) == 1, f"Expected exactly one MAF, found: {maf_files}"

    rows = list(read_maf_output(maf_files[0]))
    assert len(rows) == 1, f"Expected exactly one variant row, got {len(rows)}"
    row = rows[0]

    # The variant was accepted at the right coordinate (1-based in MAF).
    # NOTE: the MAF writer normalizes the contig name and drops the 'chr' prefix,
    # so input "chr1" is emitted as "1" in the Chromosome column.
    assert row["Chromosome"] == "1"
    assert int(row["Start_Position"]) == VARIANT_POS_0BASED + 1
    assert row["Reference_Allele"] == REF_BASE
    assert row["Tumor_Seq_Allele2"] == ALT_BASE

    ref_count = int(row["ref_count"])
    alt_count = int(row["alt_count"])
    total_count = int(row["total_count"])

    # Core assertion: exact synthetic counts.
    assert ref_count == N_REF_READS, f"ref_count: expected {N_REF_READS}, got {ref_count}"
    assert alt_count == N_ALT_READS, f"alt_count: expected {N_ALT_READS}, got {alt_count}"
    assert (
        total_count == N_REF_READS + N_ALT_READS
    ), f"total_count: expected {N_REF_READS + N_ALT_READS}, got {total_count}"

    # Counting invariants (per AGENTS.md): DP >= RD + AD, strand splits sum.
    assert total_count >= ref_count + alt_count
    ref_fwd = int(row["ref_count_forward"])
    ref_rev = int(row["ref_count_reverse"])
    alt_fwd = int(row["alt_count_forward"])
    alt_rev = int(row["alt_count_reverse"])
    assert ref_fwd + ref_rev == ref_count
    assert alt_fwd + alt_rev == alt_count
    # All reads are forward-strand in this construction.
    assert ref_rev == 0
    assert alt_rev == 0
