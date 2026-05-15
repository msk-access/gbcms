"""
Tests for INDEL-aware fragment consensus (structural ALT priority).

Verifies that when R1 and R2 of a read pair disagree on an insertion or
deletion, the read with direct CIGAR evidence (I/D op) wins the fragment
consensus unconditionally — the anchor base quality comparison is
semantically meaningless for INDEL evidence.

Prior behavior: Conflict fragments where R1=REF(M-block) and R2=ALT(I/D op)
had identical anchor BQ (e.g., both q79 in duplex), fell within the
quality threshold, and were discarded. This systematically under-counted
INDEL evidence at the fragment level.

Fix: ClassifyResult.is_structural propagates through observe() to
FragmentEvidence.has_structural_alt, and resolve() prioritizes structural
ALT over non-structural REF regardless of base quality.

Test coverage:
  1. INS conflict: R1=M-block, R2=I-op → ALT wins (core fix)
  2. INS agreement: both reads have I-op → ALT (no conflict)
  3. DEL conflict: R1=M-block, R2=D-op → ALT wins (core fix)
  4. DEL agreement: both reads have D-op → ALT (no conflict)
  5. SNP conflict tie: quality-weighted consensus unchanged (regression)
  6. INS singleton: single read with I-op → ALT (no conflict)
  7. DEL singleton: single read with D-op → ALT (no conflict)
  8. INS + REF agree: both reads M-block at INS site → REF
  9. DEL + REF agree: both reads M-block at DEL site → REF
"""

from helpers import build_bam, count_both, make_read

from gbcms._rs import Variant

# ── Helpers ──────────────────────────────────────────────────────────────


def _paired_flag(is_read1: bool, is_reverse: bool) -> int:
    """Build SAM flag for a properly paired read."""
    flag = 0x1 | 0x2  # paired + proper pair
    flag |= 0x40 if is_read1 else 0x80  # read1 or read2
    if is_reverse:
        flag |= 0x10  # reverse strand
        flag |= 0x20  # mate reverse (FR orientation: mate is forward)
    else:
        flag |= 0x20  # mate reverse (FR orientation: mate is reverse)
    return flag


def _make_paired(name, seq, start, cigar, is_read1, is_reverse, quals=None):
    """Create a properly paired read with mate info."""
    return make_read(
        name=name,
        seq=seq,
        start=start,
        cigar=cigar,
        flag=_paired_flag(is_read1, is_reverse),
        quals=quals,
    )


def _count_indel_variant(bam_path, ref_allele, alt_allele, variant_type, pos=100):
    """Count a single variant using both APIs with parity assertion."""
    variant = Variant(
        chrom="chr1",
        pos=pos,
        ref_allele=ref_allele,
        alt_allele=alt_allele,
        variant_type=variant_type,
    )
    return count_both(
        bam_path,
        [variant],
        min_mapq=0,
        min_baseq=0,
        filter_qc_failed=True,
        filter_improper_pair=False,
        filter_indel=False,
    )[0]


def _assert_invariants(counts, label=""):
    """Assert the 4 counting invariants from code-quality.md."""
    prefix = f"[{label}] " if label else ""
    assert (
        counts.dp >= counts.rd + counts.ad
    ), f"{prefix}DP invariant: dp={counts.dp} < rd+ad={counts.rd + counts.ad}"
    assert (
        counts.dpf >= counts.rdf + counts.adf
    ), f"{prefix}DPF invariant: dpf={counts.dpf} < rdf+adf={counts.rdf + counts.adf}"
    assert (
        counts.rd == counts.rd_fwd + counts.rd_rev
    ), f"{prefix}RD strand: rd={counts.rd} != {counts.rd_fwd} + {counts.rd_rev}"
    assert (
        counts.ad == counts.ad_fwd + counts.ad_rev
    ), f"{prefix}AD strand: ad={counts.ad} != {counts.ad_fwd} + {counts.ad_rev}"


# ── Test 1: INS conflict — R1=M-block, R2=I-op → ALT wins ──────────────


def test_ins_conflict_structural_alt_wins(tmp_path):
    """
    Core fix: When R1 has M-block (REF) and R2 has matching I-op (ALT),
    the structural ALT wins unconditionally.

    Prior behavior: both reads report anchor BQ (e.g., q79), difference
    is within threshold → fragment discarded. Now: structural ALT wins.

    Variant: chr1:100, REF=A, ALT=AC (1bp insertion of C after pos 100).
    R1: 20M covering pos 90-109, all 'A' → REF (no I-op)
    R2: 11M1I9M covering pos 90-109 with I(C) at pos 100 → ALT
    """
    # R1: plain M-block, no insertion → REF
    r1 = _make_paired(
        "frag1", "A" * 20, 90, ((0, 20),), is_read1=True, is_reverse=False, quals=[79] * 20
    )
    # R2: 11M + 1I + 9M → has insertion of 1bp at position 100 (11 bases from start=90)
    # seq has 21 bases: 11 + 1(inserted) + 9
    r2 = _make_paired(
        "frag1",
        "A" * 11 + "C" + "A" * 9,
        90,
        ((0, 11), (1, 1), (0, 9)),
        is_read1=False,
        is_reverse=True,
        quals=[79] * 21,
    )
    r2.next_reference_start = 90
    r2.template_length = 20
    r1.next_reference_start = 90
    r1.template_length = 20

    bam_path = build_bam(tmp_path, [r1, r2], filename="ins_conflict.bam")
    counts = _count_indel_variant(bam_path, "A", "AC", "INS")

    _assert_invariants(counts, "ins_conflict")
    assert counts.dpf == 1, f"Expected 1 fragment, got {counts.dpf}"
    assert counts.adf == 1, f"Expected adf=1 (structural ALT wins), got {counts.adf}"
    assert counts.rdf == 0, f"Expected rdf=0, got {counts.rdf}"


# ── Test 2: INS agreement — both reads have I-op → ALT ──────────────────


def test_ins_agree_both_alt(tmp_path):
    """
    Both R1 and R2 have matching I-op → ALT. No conflict.
    """
    seq = "A" * 11 + "C" + "A" * 9  # 21 bases
    cigar = ((0, 11), (1, 1), (0, 9))

    r1 = _make_paired("frag2", seq, 90, cigar, is_read1=True, is_reverse=False, quals=[30] * 21)
    r2 = _make_paired("frag2", seq, 90, cigar, is_read1=False, is_reverse=True, quals=[30] * 21)
    r1.next_reference_start = 90
    r1.template_length = 20
    r2.next_reference_start = 90
    r2.template_length = 20

    bam_path = build_bam(tmp_path, [r1, r2], filename="ins_agree.bam")
    counts = _count_indel_variant(bam_path, "A", "AC", "INS")

    _assert_invariants(counts, "ins_agree")
    assert counts.adf == 1, f"Expected adf=1 (both agree ALT), got {counts.adf}"
    assert counts.rdf == 0, f"Expected rdf=0, got {counts.rdf}"


# ── Test 3: DEL conflict — R1=M-block, R2=D-op → ALT wins ──────────────


def test_del_conflict_structural_alt_wins(tmp_path):
    """
    Core fix (DEL): When R1 has M-block (REF) and R2 has matching D-op (ALT),
    the structural ALT wins unconditionally.

    Variant: chr1:100, REF=AC, ALT=A (1bp deletion of C at pos 101).
    R1: 20M covering pos 90-109, all 'A' → REF (no D-op)
    R2: 11M1D9M → has deletion of 1bp at position 101 (11 bases from start=90)
    """
    # R1: plain M-block → REF
    r1 = _make_paired(
        "frag3", "A" * 20, 90, ((0, 20),), is_read1=True, is_reverse=False, quals=[79] * 20
    )
    # R2: 11M + 1D + 9M → deletion of 1bp at ref pos 101
    # seq has 20 bases (deletion consumes ref, not query)
    r2 = _make_paired(
        "frag3",
        "A" * 20,
        90,
        ((0, 11), (2, 1), (0, 9)),
        is_read1=False,
        is_reverse=True,
        quals=[79] * 20,
    )
    r1.next_reference_start = 90
    r1.template_length = 20
    r2.next_reference_start = 90
    r2.template_length = 21  # spans 21 ref bases due to D-op

    bam_path = build_bam(tmp_path, [r1, r2], filename="del_conflict.bam")
    counts = _count_indel_variant(bam_path, "AC", "A", "DEL")

    _assert_invariants(counts, "del_conflict")
    assert counts.dpf == 1, f"Expected 1 fragment, got {counts.dpf}"
    assert counts.adf == 1, f"Expected adf=1 (structural ALT wins), got {counts.adf}"
    assert counts.rdf == 0, f"Expected rdf=0, got {counts.rdf}"


# ── Test 4: DEL agreement — both reads have D-op → ALT ──────────────────


def test_del_agree_both_alt(tmp_path):
    """
    Both R1 and R2 have matching D-op → ALT. No conflict.
    """
    cigar = ((0, 11), (2, 1), (0, 9))

    r1 = _make_paired(
        "frag4", "A" * 20, 90, cigar, is_read1=True, is_reverse=False, quals=[30] * 20
    )
    r2 = _make_paired(
        "frag4", "A" * 20, 90, cigar, is_read1=False, is_reverse=True, quals=[30] * 20
    )
    r1.next_reference_start = 90
    r1.template_length = 21
    r2.next_reference_start = 90
    r2.template_length = 21

    bam_path = build_bam(tmp_path, [r1, r2], filename="del_agree.bam")
    counts = _count_indel_variant(bam_path, "AC", "A", "DEL")

    _assert_invariants(counts, "del_agree")
    assert counts.adf == 1, f"Expected adf=1 (both agree ALT), got {counts.adf}"
    assert counts.rdf == 0, f"Expected rdf=0, got {counts.rdf}"


# ── Test 5: SNP conflict tie — unchanged behavior ───────────────────────


def test_snp_conflict_tie_unchanged(tmp_path):
    """
    Regression test: SNP conflicts still use quality-weighted consensus.
    Equal BQ → fragment discarded. Structural flag is NOT set for SNPs.
    """
    ref_seq = "A" * 20
    alt_seq = "A" * 10 + "T" + "A" * 9  # pos 100 = T

    r1 = _make_paired(
        "frag5", ref_seq, 90, ((0, 20),), is_read1=True, is_reverse=False, quals=[30] * 20
    )
    r2 = _make_paired(
        "frag5", alt_seq, 90, ((0, 20),), is_read1=False, is_reverse=True, quals=[30] * 20
    )
    r1.next_reference_start = 90
    r1.template_length = 20
    r2.next_reference_start = 90
    r2.template_length = 20

    bam_path = build_bam(tmp_path, [r1, r2], filename="snp_tie.bam")
    counts = _count_indel_variant(bam_path, "A", "T", "SNP")

    _assert_invariants(counts, "snp_tie")
    assert counts.dpf == 1, f"Expected 1 fragment, got {counts.dpf}"
    # Equal BQ (30 vs 30), within threshold (10) → discarded
    assert counts.adf == 0, f"Expected adf=0 (tie → discarded), got {counts.adf}"
    assert counts.rdf == 0, f"Expected rdf=0 (tie → discarded), got {counts.rdf}"
    discarded = counts.dpf - (counts.rdf + counts.adf)
    assert discarded == 1, f"Expected 1 discarded fragment, got {discarded}"


# ── Test 6: INS singleton — single read with I-op → ALT ─────────────────


def test_ins_singleton(tmp_path):
    """
    Single read (no mate) with I-op → ALT. No conflict to resolve.
    """
    seq = "A" * 11 + "C" + "A" * 9
    r1 = make_read("frag6", seq, 90, ((0, 11), (1, 1), (0, 9)), quals=[30] * 21)

    bam_path = build_bam(tmp_path, [r1], filename="ins_singleton.bam")
    counts = _count_indel_variant(bam_path, "A", "AC", "INS")

    _assert_invariants(counts, "ins_singleton")
    assert counts.adf == 1, f"Expected adf=1 (singleton ALT), got {counts.adf}"


# ── Test 7: DEL singleton — single read with D-op → ALT ─────────────────


def test_del_singleton(tmp_path):
    """
    Single read (no mate) with D-op → ALT. No conflict to resolve.
    """
    r1 = make_read("frag7", "A" * 20, 90, ((0, 11), (2, 1), (0, 9)), quals=[30] * 20)

    bam_path = build_bam(tmp_path, [r1], filename="del_singleton.bam")
    counts = _count_indel_variant(bam_path, "AC", "A", "DEL")

    _assert_invariants(counts, "del_singleton")
    assert counts.adf == 1, f"Expected adf=1 (singleton ALT), got {counts.adf}"


# ── Test 8: INS REF agreement — both reads M-block → REF ────────────────


def test_ins_ref_agreement(tmp_path):
    """
    Both reads are M-block at INS site → REF. Structural flag not set.
    """
    r1 = _make_paired(
        "frag8", "A" * 20, 90, ((0, 20),), is_read1=True, is_reverse=False, quals=[30] * 20
    )
    r2 = _make_paired(
        "frag8", "A" * 20, 90, ((0, 20),), is_read1=False, is_reverse=True, quals=[30] * 20
    )
    r1.next_reference_start = 90
    r1.template_length = 20
    r2.next_reference_start = 90
    r2.template_length = 20

    bam_path = build_bam(tmp_path, [r1, r2], filename="ins_ref_agree.bam")
    counts = _count_indel_variant(bam_path, "A", "AC", "INS")

    _assert_invariants(counts, "ins_ref_agree")
    assert counts.rdf == 1, f"Expected rdf=1 (both agree REF), got {counts.rdf}"
    assert counts.adf == 0, f"Expected adf=0, got {counts.adf}"


# ── Test 9: DEL REF agreement — both reads M-block → REF ────────────────


def test_del_ref_agreement(tmp_path):
    """
    Both reads are M-block at DEL site → REF. Structural flag not set.
    """
    r1 = _make_paired(
        "frag9", "A" * 20, 90, ((0, 20),), is_read1=True, is_reverse=False, quals=[30] * 20
    )
    r2 = _make_paired(
        "frag9", "A" * 20, 90, ((0, 20),), is_read1=False, is_reverse=True, quals=[30] * 20
    )
    r1.next_reference_start = 90
    r1.template_length = 20
    r2.next_reference_start = 90
    r2.template_length = 20

    bam_path = build_bam(tmp_path, [r1, r2], filename="del_ref_agree.bam")
    counts = _count_indel_variant(bam_path, "AC", "A", "DEL")

    _assert_invariants(counts, "del_ref_agree")
    assert counts.rdf == 1, f"Expected rdf=1 (both agree REF), got {counts.rdf}"
    assert counts.adf == 0, f"Expected adf=0, got {counts.adf}"
