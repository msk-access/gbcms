import pysam
import pytest
from helpers import build_bam, count_both, make_read

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
    # read_secondary: SKIPPED — supplementary/secondary never increment read-level
    #   depth, independent of the filter flag (they are not first-class observations).
    # Total = 4
    assert counts.dp == 4

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
    # read_qc_fail removed; read_secondary always skipped. Total = 3
    assert counts.dp == 3

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
    # read_improper removed; read_secondary always skipped. Total = 3
    assert counts.dp == 3

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
    # read_indel removed; read_secondary always skipped. Total = 3
    assert counts.dp == 3

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

    # No filters: 4 reads count (read_secondary is always skipped at read level)
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
    assert counts.dp == 4

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


def test_supplementary_shared_qname_not_double_counted(tmp_path):
    """The core double-count case: a supplementary segment sharing a QNAME with its
    primary must not inflate read-level DP even with --no-filter-supplementary.
    count_both asserts binned↔legacy parity, so both paths must agree on DP=1."""
    bam_path = tmp_path / "supp.bam"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        for flag in (2, 0x800 | 2):  # primary (proper pair), then its supplementary
            a = pysam.AlignedSegment()
            a.query_name = "frag1"  # same QNAME → same fragment
            a.query_sequence = "A" * 100
            a.flag = flag
            a.reference_id = 0
            a.reference_start = 100
            a.mapping_quality = 60
            a.cigartuples = [(0, 100)]
            outf.write(a)
    pysam.index(str(bam_path))

    variant = gbcms_rs.Variant("chr1", 150, "A", "T", "SNP")
    counts = count_both(
        str(bam_path),
        [variant],
        min_mapq=0,
        min_baseq=0,
        filter_secondary=False,
        filter_supplementary=False,
    )[0]
    assert counts.dp == 1, f"supplementary must not double-count DP, got {counts.dp}"


def test_supplementary_only_locus_is_seen_at_fragment_level_when_opted_in(tmp_path):
    """The other half of the contract: opting out must actually admit the evidence.

    Pairs with `test_supplementary_shared_qname_not_double_counted` above. That one pins
    read-level DP against inflation when a primary and its supplementary overlap the *same*
    locus. This covers the opposite arrangement — a locus reached **only** by the
    supplementary segment, which is what a molecule spanning a large deletion looks like.

    Previously both flags were inert: an unconditional skip dropped these records before
    fragment evidence, so such a locus reported `dpf=0` — not a filtered read but a wrong
    answer, since a molecule demonstrably covers it. Now `filter_supplementary=False`
    admits it to fragment evidence (where the QNAME hash makes double-counting impossible)
    while read-level DP keeps its promise and stays 0.
    """
    bam_path = tmp_path / "split.bam"
    header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 2000, "SN": "chr1"}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        # primary covers 100–200; its supplementary segment covers 1000–1100 only
        for flag, start in ((2, 100), (0x800 | 2, 1000)):
            a = pysam.AlignedSegment()
            a.query_name = "frag1"
            a.query_sequence = "A" * 100
            a.flag = flag
            a.reference_id = 0
            a.reference_start = start
            a.mapping_quality = 60
            a.cigartuples = [(0, 100)]
            outf.write(a)
    pysam.index(str(bam_path))

    far = gbcms_rs.Variant("chr1", 1050, "A", "T", "SNP")  # only the supplementary reaches it
    filtered = count_both(str(bam_path), [far], min_mapq=0, min_baseq=0)[0]
    assert filtered.dpf == 0, "default must still exclude it"

    admitted = count_both(
        str(bam_path), [far], min_mapq=0, min_baseq=0, filter_supplementary=False
    )[0]
    assert admitted.dpf == 1, "opting out did not admit the supplementary to fragment evidence"
    assert admitted.dp == 0, (
        "read-level DP must stay 0 — a supplementary is not a first-class read observation, "
        "independent of the filter flag"
    )


# ── qc-failed filter contract ────────────────────────────────────────────


def _qc_contract_reads(flag_alt_and_some_ref):
    """Mixed population at chr1:150 (A>T): 6 REF reads, 4 ALT reads.

    With flag_alt_and_some_ref, the 4 ALT reads and 2 of the REF reads carry
    the QC-fail bit (0x200) — same reads, same sequences, only the flag differs.
    """
    qc = 0x200 if flag_alt_and_some_ref else 0
    reads = []
    for i in range(6):
        flag = qc if flag_alt_and_some_ref and i < 2 else 0
        reads.append(make_read(f"ref{i}", "A" * 100, 100, ((0, 100),), flag=flag))
    alt_seq = "A" * 50 + "T" + "A" * 49
    for i in range(4):
        reads.append(make_read(f"alt{i}", alt_seq, 100, ((0, 100),), flag=qc))
    return reads


def _assert_counting_invariants(c):
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev


def test_qc_failed_filter_contract(tmp_path):
    """Three-way contract for filter_qc_failed, at read AND fragment level.

    1. filter on drops exactly the QC-flagged reads (here: all ALT evidence);
    2. filter off restores counts identical to the same reads unflagged;
    3. the filter has no effect on a BAM without QC-fail flags.
    """
    variant = gbcms_rs.Variant("chr1", 150, "A", "T", "SNP")
    clean_bam = build_bam(tmp_path, _qc_contract_reads(False), "qc_clean.bam")
    flagged_bam = build_bam(tmp_path, _qc_contract_reads(True), "qc_flagged.bam")

    def count(bam, qc_filter):
        c = count_both(
            bam,
            [variant],
            min_mapq=0,
            min_baseq=0,
            filter_secondary=False,
            filter_supplementary=False,
            filter_qc_failed=qc_filter,
            filter_improper_pair=False,
            filter_indel=False,
        )[0]
        _assert_counting_invariants(c)
        return c

    baseline = count(clean_bam, False)
    assert (baseline.dp, baseline.rd, baseline.ad) == (10, 6, 4)

    # 1. Filter on: the 6 flagged reads (2 REF + 4 ALT) vanish — ALT evidence
    #    disappears entirely, so a miscounted no-op here would fake a variant call.
    on = count(flagged_bam, True)
    assert (on.dp, on.rd, on.ad) == (4, 4, 0)
    assert (on.dpf, on.rdf, on.adf) == (4, 4, 0)

    # 2. Filter off: flagged reads are admitted; every read- and fragment-level
    #    field matches the unflagged baseline exactly.
    off = count(flagged_bam, False)
    for field in ("dp", "rd", "ad", "dpf", "rdf", "adf", "rd_fwd", "rd_rev", "ad_fwd", "ad_rev"):
        assert getattr(off, field) == getattr(baseline, field), field

    # 3. No QC flags present: the filter setting must not change anything.
    clean_on = count(clean_bam, True)
    for field in ("dp", "rd", "ad", "dpf", "rdf", "adf"):
        assert getattr(clean_on, field) == getattr(baseline, field), field


# ── secondary/supplementary read-level exclusion contract ────────────────


def test_secondary_admission_is_fragment_only(tmp_path):
    """Admitted secondary alignments must never reach read-level counts.

    With --no-filter-secondary, secondary records are admitted to FRAGMENT
    evidence only. Before the fix they also received read-level ref/alt calls
    while DP correctly excluded them, so shipped output violated DP >= RD+AD
    (issue #89; observed on real RNA BAMs as e.g. rd+ad=182 > dp=124).

    Population at chr1:150 (A>T): 4 primary REF, 2 primary ALT, and 3
    secondary ALT alignments with distinct QNAMEs (multimapper-style, so each
    secondary is also a distinct fragment key).
    """
    reads = []
    for i in range(4):
        reads.append(make_read(f"ref{i}", "A" * 100, 100, ((0, 100),)))
    alt_seq = "A" * 50 + "T" + "A" * 49
    for i in range(2):
        reads.append(make_read(f"alt{i}", alt_seq, 100, ((0, 100),)))
    for i in range(3):
        reads.append(make_read(f"sec{i}", alt_seq, 100, ((0, 100),), flag=0x100))
    bam = build_bam(tmp_path, reads, "sec_leak.bam")
    variant = gbcms_rs.Variant("chr1", 150, "A", "T", "SNP")

    admitted = count_both(
        bam,
        [variant],
        min_mapq=0,
        min_baseq=0,
        filter_secondary=False,
        filter_supplementary=False,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
    )[0]
    # Read level: primaries only — 6 first-class reads, 2 of them ALT.
    assert (
        admitted.dp >= admitted.rd + admitted.ad
    ), f"invariant violated: {admitted.rd}+{admitted.ad} > {admitted.dp}"
    assert (admitted.dp, admitted.rd, admitted.ad) == (6, 4, 2)
    assert admitted.any_alt == admitted.ad + admitted.partial_alt
    # Fragment level: the 3 secondary records are distinct molecules, admitted.
    assert admitted.dpf >= admitted.rdf + admitted.adf
    assert (admitted.dpf, admitted.rdf, admitted.adf) == (9, 4, 5)

    # Default filters: secondaries excluded everywhere.
    default = count_both(
        bam,
        [variant],
        min_mapq=0,
        min_baseq=0,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
    )[0]
    assert (default.dp, default.rd, default.ad) == (6, 4, 2)
    assert (default.dpf, default.rdf, default.adf) == (6, 4, 2)
