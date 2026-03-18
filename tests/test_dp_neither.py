"""
Tests for Gap 1D: Non-REF Total Depth (DP) Dropout fix.

Verifies that reads classified as neither REF nor ALT (e.g., third alleles,
duplex N bases) are still included in DP and DPF counts, but excluded from
RD/AD counts. Before this fix, such reads hit a `continue` and were
completely excluded from all counting including depth.

Key invariant: dp >= rd + ad (not dp == rd + ad).
"""

import pysam
import pytest
from helpers import count_both

from gbcms._rs import Variant, count_bam


@pytest.fixture
def neither_bam(tmp_path):
    """
    Generates a synthetic BAM for testing DP counting of 'neither' reads.

    Reference: chr1, 1000 bases of 'A'.
    Variant: chr1:100 (0-based), REF=A, ALT=T (SNP)

    Read layout:
        3 reads with 'A' at pos 100 → REF
        2 reads with 'T' at pos 100 → ALT
        2 reads with 'C' at pos 100 → NEITHER (third allele)
        1 read  with 'G' at pos 100 → NEITHER (third allele)
    Total: 8 reads, dp=8, rd=3, ad=2, neither=3
    """
    bam_path = tmp_path / "neither.bam"
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:

        def make_read(qname, base_at_pos, is_reverse=False):
            a = pysam.AlignedSegment()
            a.query_name = qname
            # Place the target base at index 5 (pos 95+5=100)
            a.query_sequence = "AAAAA" + base_at_pos + "AAAA"
            a.flag = 16 if is_reverse else 0
            a.reference_id = 0
            a.reference_start = 95
            a.mapping_quality = 60
            a.cigartuples = [
                (0, 10),
            ]  # 10M
            a.query_qualities = [30] * 10  # type: ignore[assignment]
            outf.write(a)

        # REF reads (A)
        make_read("read_ref_1", "A")
        make_read("read_ref_2", "A")
        make_read("read_ref_3", "A", is_reverse=True)

        # ALT reads (T)
        make_read("read_alt_1", "T")
        make_read("read_alt_2", "T", is_reverse=True)

        # NEITHER reads (C — third allele, not REF or ALT)
        make_read("read_neither_c1", "C")
        make_read("read_neither_c2", "C", is_reverse=True)

        # NEITHER reads (G — fourth allele)
        make_read("read_neither_g1", "G")

    sorted_bam = tmp_path / "neither.sorted.bam"
    pysam.sort("-o", str(sorted_bam), str(bam_path))
    pysam.index(str(sorted_bam))
    return str(sorted_bam)


def test_dp_includes_neither_reads(neither_bam):
    """DP must include all quality-filtered reads, including 'neither' reads."""
    variant = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")

    results = count_bam(
        neither_bam,
        [variant],
        decomposed=[None],
        min_mapq=20,
        min_baseq=20,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )

    counts = results[0]

    # REF: 3, ALT: 2, NEITHER: 3 → DP must be 8, not 5
    assert counts.rd == 3, f"Expected rd=3, got {counts.rd}"
    assert counts.ad == 2, f"Expected ad=2, got {counts.ad}"
    assert counts.dp == 8, f"Expected dp=8 (3+2+3 neither), got {counts.dp}"

    # Key invariant: DP >= RD + AD
    assert (
        counts.dp >= counts.rd + counts.ad
    ), f"dp ({counts.dp}) must be >= rd+ad ({counts.rd + counts.ad})"


def test_dp_strand_includes_neither(neither_bam):
    """Strand-specific DP must include 'neither' reads."""
    variant = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")

    results = count_bam(
        neither_bam,
        [variant],
        decomposed=[None],
        min_mapq=20,
        min_baseq=20,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )

    counts = results[0]

    # Fwd: ref_1, ref_2, alt_1, neither_c1, neither_g1 = 5
    # Rev: ref_3, alt_2, neither_c2 = 3
    assert counts.dp_fwd == 5, f"Expected dp_fwd=5, got {counts.dp_fwd}"
    assert counts.dp_rev == 3, f"Expected dp_rev=3, got {counts.dp_rev}"
    assert counts.dp_fwd + counts.dp_rev == counts.dp


def test_dpf_includes_neither_fragments(neither_bam):
    """DPF must include fragments where reads are 'neither'."""
    variant = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")

    results = count_bam(
        neither_bam,
        [variant],
        decomposed=[None],
        min_mapq=20,
        min_baseq=20,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )

    counts = results[0]

    # All 8 reads have unique QNAMEs → 8 fragments
    # REF fragments: 3, ALT fragments: 2, NEITHER fragments: 3
    assert counts.dpf == 8, f"Expected dpf=8, got {counts.dpf}"
    assert counts.rdf == 3, f"Expected rdf=3, got {counts.rdf}"
    assert counts.adf == 2, f"Expected adf=2, got {counts.adf}"

    # Key invariant: DPF >= RDF + ADF
    assert counts.dpf >= counts.rdf + counts.adf


# ── count_bam_binned parity tests ────────────────────────────────────────


def test_dp_includes_neither_reads_binned(neither_bam):
    """count_bam_binned: DP includes third-allele reads."""
    variant = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")
    counts = count_both(
        neither_bam,
        [variant],
        min_mapq=20,
        min_baseq=20,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
    )[0]
    assert counts.rd == 3
    assert counts.ad == 2
    assert counts.dp == 8
    assert counts.dp >= counts.rd + counts.ad


def test_dp_strand_includes_neither_binned(neither_bam):
    """count_bam_binned: Strand DP includes third-allele reads."""
    variant = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")
    counts = count_both(
        neither_bam,
        [variant],
        min_mapq=20,
        min_baseq=20,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
    )[0]
    assert counts.dp_fwd == 5
    assert counts.dp_rev == 3
    assert counts.dp_fwd + counts.dp_rev == counts.dp


def test_dpf_includes_neither_fragments_binned(neither_bam):
    """count_bam_binned: DPF includes third-allele fragments."""
    variant = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")
    counts = count_both(
        neither_bam,
        [variant],
        min_mapq=20,
        min_baseq=20,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
    )[0]
    assert counts.dpf == 8
    assert counts.rdf == 3
    assert counts.adf == 2
    assert counts.dpf >= counts.rdf + counts.adf
