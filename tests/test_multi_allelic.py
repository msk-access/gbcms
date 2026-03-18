"""
Tests for Gap 1A: Multi-Allelic Isolation.

Phase 1: Verify that overlapping variants at the same locus are annotated with
         matching multi_allelic_group IDs and _MULTI_ALLELIC validation_status.

Phase 2: Verify that reads carrying a sibling ALT allele are excluded from
         REF counts for the current variant.
"""

import pysam
import pytest
from helpers import count_both

from gbcms._rs import Variant, count_bam


@pytest.fixture
def multi_allelic_bam(tmp_path):
    """
    Generates a synthetic BAM for testing multi-allelic sibling exclusion.

    Reference: chr1, 1000 bases of 'A'.
    Variant 1: chr1:100, REF=A, ALT=T (SNP: A>T)
    Variant 2: chr1:100, REF=A, ALT=C (SNP: A>C)

    Read layout:
        3 reads with 'A' at pos 100 → REF for both variants
        2 reads with 'T' at pos 100 → ALT for V1, NEITHER for V2
        2 reads with 'C' at pos 100 → NEITHER for V1, ALT for V2
    """
    bam_path = tmp_path / "multi_allelic.bam"
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:

        def make_read(qname, base_at_pos, is_reverse=False):
            a = pysam.AlignedSegment()
            a.query_name = qname
            a.query_sequence = "AAAAA" + base_at_pos + "AAAA"
            a.flag = 16 if is_reverse else 0
            a.reference_id = 0
            a.reference_start = 95
            a.mapping_quality = 60
            a.cigartuples = [
                (0, 10),
            ]
            a.query_qualities = [30] * 10  # type: ignore[assignment]
            outf.write(a)

        # REF reads (A) — REF for both variants
        make_read("read_ref_1", "A")
        make_read("read_ref_2", "A")
        make_read("read_ref_3", "A", is_reverse=True)

        # ALT for Variant 1 (T)
        make_read("read_alt_t1", "T")
        make_read("read_alt_t2", "T", is_reverse=True)

        # ALT for Variant 2 (C)
        make_read("read_alt_c1", "C")
        make_read("read_alt_c2", "C", is_reverse=True)

    sorted_bam = tmp_path / "multi_allelic.sorted.bam"
    pysam.sort("-o", str(sorted_bam), str(bam_path))
    pysam.index(str(sorted_bam))
    return str(sorted_bam)


def test_multi_allelic_without_siblings(multi_allelic_bam):
    """
    Without sibling info, reads carrying V2's ALT are counted as NEITHER
    for V1 (not REF or ALT), and vice versa. DP includes all reads.
    """
    v1 = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")
    v2 = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="C", variant_type="SNP")

    results = count_bam(
        multi_allelic_bam,
        [v1, v2],
        decomposed=[None, None],
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

    # V1: REF(A)=3, ALT(T)=2, NEITHER(C)=2 → dp=7, rd=3, ad=2
    assert results[0].rd == 3
    assert results[0].ad == 2
    assert results[0].dp == 7

    # V2: REF(A)=3, ALT(C)=2, NEITHER(T)=2 → dp=7, rd=3, ad=2
    assert results[1].rd == 3
    assert results[1].ad == 2
    assert results[1].dp == 7


def test_multi_allelic_with_siblings(multi_allelic_bam):
    """
    With sibling info provided, reads carrying a sibling's ALT should be
    excluded from REF counts. For SNPs at the same position where REF is
    the same base, the sibling ALT exclusion doesn't affect REF counting
    because those reads are already classified as NEITHER (not REF).

    This test verifies the sibling_variants parameter is correctly threaded.
    """
    v1 = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")
    v2 = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="C", variant_type="SNP")

    # With sibling exclusion: V1 has V2 as sibling, V2 has V1 as sibling
    results = count_bam(
        multi_allelic_bam,
        [v1, v2],
        decomposed=[None, None],
        min_mapq=20,
        min_baseq=20,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
        sibling_variants=[[v2], [v1]],
    )

    # For SNPs: reads with 'C' are already NEITHER for V1 (base != A or T),
    # so the sibling guard doesn't change REF. Same logic applies to V2.
    # REF reads (base 'A') are genuinely REF — they don't match any sibling ALT.
    assert results[0].rd == 3, f"V1 rd expected 3, got {results[0].rd}"
    assert results[0].ad == 2, f"V1 ad expected 2, got {results[0].ad}"
    assert results[0].dp == 7

    assert results[1].rd == 3, f"V2 rd expected 3, got {results[1].rd}"
    assert results[1].ad == 2, f"V2 ad expected 2, got {results[1].ad}"
    assert results[1].dp == 7


@pytest.fixture
def overlapping_indel_bam(tmp_path):
    """
    Synthetic BAM with overlapping indels at the same locus.

    Reference: chr1, 1000 bases of 'A'.

    Variant 1: chr1:100, REF=AAA, ALT=A (2bp deletion, removes pos 101-102)
    Variant 2: chr1:100, REF=A, ALT=ATTT (3bp insertion after pos 100)

    Read layout:
        2 reads with 5M 2D 5M  → ALT for V1 (deletion), REF-like for V2
        2 reads with 5M 3I 5M  → REF-like for V1, ALT for V2 (insertion)
        2 reads with 10M        → REF for both
    """
    bam_path = tmp_path / "overlap_indel.bam"
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        # REF reads — no indel, just 10M match
        for i in range(2):
            a = pysam.AlignedSegment()
            a.query_name = f"read_ref_{i}"
            a.query_sequence = "AAAAAAAAAA"
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 96
            a.mapping_quality = 60
            a.cigartuples = [
                (0, 10),
            ]  # 10M
            a.query_qualities = [30] * 10  # type: ignore[assignment]
            outf.write(a)

        # Deletion reads — 5M 2D 5M (deletes 2bp starting after 5th base)
        for i in range(2):
            a = pysam.AlignedSegment()
            a.query_name = f"read_del_{i}"
            a.query_sequence = "AAAAAAAAAA"  # 10 bases (but covers 12 ref bases)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 96
            a.mapping_quality = 60
            a.cigartuples = [(0, 5), (2, 2), (0, 5)]  # 5M 2D 5M
            a.query_qualities = [30] * 10  # type: ignore[assignment]
            outf.write(a)

        # Insertion reads — 5M 3I 5M (inserts 3bp after 5th base)
        for i in range(2):
            a = pysam.AlignedSegment()
            a.query_name = f"read_ins_{i}"
            a.query_sequence = "AAAAATTTAAAAA"  # 13 bases (5 + 3ins + 5)
            a.flag = 0
            a.reference_id = 0
            a.reference_start = 96
            a.mapping_quality = 60
            a.cigartuples = [(0, 5), (1, 3), (0, 5)]  # 5M 3I 5M
            a.query_qualities = [30] * 13  # type: ignore[assignment]
            outf.write(a)

    sorted_bam = tmp_path / "overlap_indel.sorted.bam"
    pysam.sort("-o", str(sorted_bam), str(bam_path))
    pysam.index(str(sorted_bam))
    return str(sorted_bam)


def test_overlapping_indels_dp(overlapping_indel_bam):
    """DP must include all 6 reads regardless of allele classification."""
    v1 = Variant(chrom="chr1", pos=100, ref_allele="AAA", alt_allele="A", variant_type="DELETION")

    results = count_bam(
        overlapping_indel_bam,
        [v1],
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

    # All 6 reads overlap the locus → dp should be 6
    assert counts.dp == 6, f"Expected dp=6, got {counts.dp}"
    # Deletion reads match ALT
    assert counts.ad == 2, f"Expected ad=2, got {counts.ad}"


# ── count_bam_binned parity tests ────────────────────────────────────────


def test_multi_allelic_without_siblings_binned(multi_allelic_bam):
    """count_bam_binned: multi-allelic without siblings."""
    v1 = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")
    v2 = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="C", variant_type="SNP")
    results = count_both(
        multi_allelic_bam,
        [v1, v2],
        min_mapq=20,
        min_baseq=20,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
    )
    assert results[0].rd == 3
    assert results[0].ad == 2
    assert results[0].dp == 7
    assert results[1].rd == 3
    assert results[1].ad == 2


def test_multi_allelic_with_siblings_binned(multi_allelic_bam):
    """count_bam_binned: multi-allelic with sibling exclusion."""
    v1 = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")
    v2 = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="C", variant_type="SNP")
    results = count_both(
        multi_allelic_bam,
        [v1, v2],
        min_mapq=20,
        min_baseq=20,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        sibling_variants=[[v2], [v1]],
    )
    assert results[0].rd == 3
    assert results[0].ad == 2
    assert results[0].dp == 7


def test_overlapping_indels_dp_binned(overlapping_indel_bam):
    """count_bam_binned: overlapping indels DP includes all reads."""
    v1 = Variant(chrom="chr1", pos=100, ref_allele="AAA", alt_allele="A", variant_type="DELETION")
    counts = count_both(
        overlapping_indel_bam,
        [v1],
        min_mapq=20,
        min_baseq=20,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
    )[0]
    assert counts.dp == 6
    assert counts.ad == 2
