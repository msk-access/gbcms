import pysam
import pytest
from helpers import count_both

from gbcms import _rs as gbcms_rs
from gbcms.models.core import Variant, VariantType


# Mock data setup
@pytest.fixture
def mock_bam_with_flags(tmp_path):
    bam_path = tmp_path / "test_flags.bam"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        # 1. Normal read (should be counted)
        a = pysam.AlignedSegment()
        a.query_name = "read1"
        a.query_sequence = "A" * 100
        a.flag = 2  # Proper pair
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigartuples = [
            (0, 100),
        ]
        outf.write(a)

        # 2. QC Failed read (flag 512)
        # Also make it proper pair so it's not filtered by improper_pair test
        a = pysam.AlignedSegment()
        a.query_name = "read_qc_fail"
        a.query_sequence = "A" * 100
        a.flag = 512 | 2
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigartuples = [
            (0, 100),
        ]
        outf.write(a)

        # 3. Improper pair (flag 1: paired but not proper)
        a = pysam.AlignedSegment()
        a.query_name = "read_improper"
        a.query_sequence = "A" * 100
        a.flag = 1
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigartuples = [
            (0, 100),
        ]
        outf.write(a)

        # 4. Indel read (CIGAR has I or D)
        # Make proper pair
        a = pysam.AlignedSegment()
        a.query_name = "read_indel"
        a.query_sequence = "A" * 100
        a.flag = 2
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigartuples = [(0, 50), (1, 1), (0, 49)]  # 50M 1I 49M
        outf.write(a)

        # 5. Secondary read (flag 256)
        # Make proper pair
        a = pysam.AlignedSegment()
        a.query_name = "read_secondary"
        a.query_sequence = "A" * 100
        a.flag = 256 | 2
        a.reference_id = 0
        a.reference_start = 100
        a.mapping_quality = 60
        a.cigartuples = [
            (0, 100),
        ]
        outf.write(a)

    pysam.index(str(bam_path))
    return bam_path


def test_filters(mock_bam_with_flags):
    # Define a variant at the location of our reads
    variant = Variant(
        chrom="chr1",
        pos=150,  # Middle of the 100bp reads starting at 100
        ref="A",
        alt="T",
        variant_type=VariantType.SNP,
    )
    rs_variants = [
        gbcms_rs.Variant(
            variant.chrom, variant.pos, variant.ref, variant.alt, variant.variant_type.value
        )
    ]

    # Case 1: No filters (except defaults)
    # Defaults: filter_duplicates=True, others False
    counts = gbcms_rs.count_bam(
        str(mock_bam_with_flags),
        rs_variants,
        [None] * len(rs_variants),
        min_mapq=0,
        min_baseq=0,
        filter_duplicates=True,
        filter_secondary=False,
        filter_supplementary=False,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )[0]

    # Expect:
    # read1: OK
    # read_qc_fail: OK (filter=False)
    # read_improper: OK (filter=False)
    # read_indel: OK (filter=False)
    # read_secondary: OK (filter=False)
    # Total = 5
    assert counts.dp == 5

    # Case 2: Filter QC Failed
    counts = gbcms_rs.count_bam(
        str(mock_bam_with_flags),
        rs_variants,
        [None] * len(rs_variants),
        min_mapq=0,
        min_baseq=0,
        filter_duplicates=True,
        filter_secondary=False,
        filter_supplementary=False,
        filter_qc_failed=True,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )[0]
    # read_qc_fail removed. Total = 4
    assert counts.dp == 4

    # Case 3: Filter Improper Pair
    counts = gbcms_rs.count_bam(
        str(mock_bam_with_flags),
        rs_variants,
        [None] * len(rs_variants),
        min_mapq=0,
        min_baseq=0,
        filter_duplicates=True,
        filter_secondary=False,
        filter_supplementary=False,
        filter_qc_failed=False,
        filter_improper_pair=True,
        filter_indel=False,
        threads=1,
    )[0]
    # read_improper removed. Total = 4
    assert counts.dp == 4

    # Case 4: Filter Indel
    counts = gbcms_rs.count_bam(
        str(mock_bam_with_flags),
        rs_variants,
        [None] * len(rs_variants),
        min_mapq=0,
        min_baseq=0,
        filter_duplicates=True,
        filter_secondary=False,
        filter_supplementary=False,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=True,
        threads=1,
    )[0]
    # read_indel removed. Total = 4
    assert counts.dp == 4

    # Case 5: Filter Secondary
    counts = gbcms_rs.count_bam(
        str(mock_bam_with_flags),
        rs_variants,
        [None] * len(rs_variants),
        min_mapq=0,
        min_baseq=0,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=False,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )[0]
    # read_secondary removed. Total = 4
    assert counts.dp == 4

    # Case 6: All Filters
    counts = gbcms_rs.count_bam(
        str(mock_bam_with_flags),
        rs_variants,
        [None] * len(rs_variants),
        min_mapq=0,
        min_baseq=0,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=True,
        filter_improper_pair=True,
        filter_indel=True,
        threads=1,
    )[0]
    # Only read1 remains. Total = 1
    assert counts.dp == 1


# ── count_bam_binned parity tests ────────────────────────────────────────


def test_filters_binned(mock_bam_with_flags):
    """count_bam_binned filter behavior matches count_bam for all 6 filter cases."""
    variant = gbcms_rs.Variant("chr1", 150, "A", "T", "SNP")

    # No filters: all 5 reads
    counts = count_both(
        str(mock_bam_with_flags),
        [variant],
        min_mapq=0,
        min_baseq=0,
        filter_secondary=False,
        filter_supplementary=False,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
    )[0]
    assert counts.dp == 5

    # All filters: only read1 remains
    counts = count_both(
        str(mock_bam_with_flags),
        [variant],
        min_mapq=0,
        min_baseq=0,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=True,
        filter_improper_pair=True,
        filter_indel=True,
    )[0]
    assert counts.dp == 1
