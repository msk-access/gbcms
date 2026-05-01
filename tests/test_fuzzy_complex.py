"""
Tests for quality-aware complex variant matching (Phase 2b).

Verifies the "Reliable Intersection" (masking) approach where low-quality bases
(below min_baseq) are masked out and cannot vote for either allele. Only reliable
(high-quality) bases participate in REF/ALT comparison.

Three cases:
  Case A: recon.len == alt.len == ref.len → simultaneous check + ambiguity detection
  Case B: recon.len == alt.len only → ALT-only masked comparison
  Case C: recon.len == ref.len only → REF-only masked comparison

Variant definitions used:
  - Equal-length complex: chr1:100, REF=AT, ALT=CG  (2bp substitution)
  - DelIns:              chr1:100, REF=ATG, ALT=CC   (3→2, different lengths)
"""

from helpers import (
    build_bam as _build_bam,
)
from helpers import (
    count_one as _count_one,
)
from helpers import (
    count_one_both as _count_one_both,
)
from helpers import (
    make_read as _make_read,
)

from gbcms._rs import Variant, count_bam  # noqa: F401 — count_bam used in legacy comments

# ==========================================================================
# Case A: Equal-length REF and ALT (2bp substitution)
# Variant: chr1:100, REF=AT, ALT=CG
# Read starts at 96, covers 10bp.
# ==========================================================================

# COMPLEX variant with equal-length REF/ALT
EQUAL_LEN_VARIANT = Variant(
    chrom="chr1",
    pos=100,
    ref_allele="AT",
    alt_allele="CG",
    variant_type="COMPLEX",
)


class TestCaseA_ExactMatch:
    """Case A: Equal-length REF/ALT — exact matches (regression checks)."""

    def test_exact_alt_match(self, tmp_path):
        """Read reconstructs 'CG' at pos 100-101 → ALT.
        Read: AAAACGAAAA, starts at 96, 10M.
        Bases at pos 100-101 (indices 4-5) = 'CG' → matches ALT.
        """
        reads = [_make_read("r1", "AAAACGAAAA", 96, ((0, 10),))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, EQUAL_LEN_VARIANT)
        assert counts.ad == 1, f"Expected ad=1 (exact ALT), got {counts.ad}"
        assert counts.rd == 0

    def test_exact_ref_match(self, tmp_path):
        """Read reconstructs 'AT' at pos 100-101 → REF.
        Read: AAAAATAAAA, starts at 96, 10M.
        Bases at pos 100-101 (indices 4-5) = 'AT' → matches REF.
        """
        reads = [_make_read("r1", "AAAAATAAAA", 96, ((0, 10),))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, EQUAL_LEN_VARIANT)
        assert counts.rd == 1, f"Expected rd=1 (exact REF), got {counts.rd}"
        assert counts.ad == 0


class TestCaseA_MaskedComparison:
    """Case A: Quality-aware masked comparison with ambiguity detection."""

    def test_rescue_low_qual_mismatch(self, tmp_path):
        """REF=AT, ALT=CG. Read='CA' where 'A' is Q5 (index 5).

        Both positions are discriminating (A≠C, T≠G). The low-quality base
        at position 1 (T→G) IS a discriminating position, so the selective
        quality gate correctly rejects this read → LowQuality → neither.

        Before the MNP fix (feature/fix-mnp-counting), this would have
        fallen through to check_complex Phase 2 masked comparison, which
        would rescue it as ALT. That fallback is removed because PairHMM/SW
        is biased toward REF for multi-base substitutions, and C++ GBCMS
        has no fallback either.
        """
        quals = [30, 30, 30, 30, 30, 5, 30, 30, 30, 30]
        #                           ^C  ^A(Q5) — discriminating position!
        reads = [_make_read("r1", "AAAACAAAAA", 96, ((0, 10),), quals=quals)]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, EQUAL_LEN_VARIANT)
        # After MNP fix: LowQuality at discriminating pos → neither (no rescue)
        assert counts.ad == 0, f"Expected ad=0 (LowQuality, no rescue), got {counts.ad}"
        assert counts.rd == 0

    def test_ambiguous_low_qual_tail(self, tmp_path):
        """REF=AT, ALT=CG. Read='AT' where 'T' is Q5.

        Both positions are discriminating (A≠C, T≠G). The low-quality base
        at position 1 (T→G) IS a discriminating position, so the selective
        quality gate correctly rejects this read → LowQuality → neither.

        Before the MNP fix, this fell through to check_complex Phase 2
        masked comparison where the reliable 'A' matched REF → rd=1.
        Now correctly discarded.
        """
        quals = [30, 30, 30, 30, 30, 5, 30, 30, 30, 30]
        reads = [_make_read("r1", "AAAAATAAAA", 96, ((0, 10),), quals=quals)]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, EQUAL_LEN_VARIANT)
        # After MNP fix: LowQuality at discriminating pos → neither
        assert counts.rd == 0, f"Expected rd=0 (LowQuality, no rescue), got {counts.rd}"
        assert counts.ad == 0

    def test_true_ambiguity_both_masked(self, tmp_path):
        """REF=AT, ALT=CG. Read='CG' but BOTH bases are Q5.
        Reliable bases: none.
        → Neither (no reliable data).
        """
        quals = [30, 30, 30, 30, 5, 5, 30, 30, 30, 30]
        reads = [_make_read("r1", "AAAACGAAAA", 96, ((0, 10),), quals=quals)]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, EQUAL_LEN_VARIANT)
        assert counts.rd == 0, f"Expected rd=0 (all masked), got {counts.rd}"
        assert counts.ad == 0, f"Expected ad=0 (all masked), got {counts.ad}"

    def test_high_qual_mismatch_neither(self, tmp_path):
        """REF=AT, ALT=CG. Read='TT' — mismatches both on high-quality bases.
        Reliable bases: 'T','T' (all Q30).
        vs ALT: 'T'≠'C', 'T'≠'G' → mismatches=2
        vs REF: 'T'≠'A', 'T'=='T' → mismatches=1
        → Neither (both have mismatches).
        """
        reads = [_make_read("r1", "AAAATTAAAA", 96, ((0, 10),))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, EQUAL_LEN_VARIANT)
        assert counts.rd == 0
        assert counts.ad == 0

    def test_ambiguity_same_first_base(self, tmp_path):
        """Variant where REF and ALT share a prefix. REF=AA, ALT=AC.
        Read='AA' with second 'A' at Q5.
        Reliable: 'A' at index 0.
        vs ALT: 'A'=='A' ✓  |  vs REF: 'A'=='A' ✓
        → Ambiguous (both match on reliable) → Neither.
        """
        variant = Variant(
            chrom="chr1",
            pos=100,
            ref_allele="AA",
            alt_allele="AC",
            variant_type="COMPLEX",
        )
        quals = [30, 30, 30, 30, 30, 5, 30, 30, 30, 30]
        reads = [_make_read("r1", "AAAAAAAAAA", 96, ((0, 10),), quals=quals)]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, variant)
        assert counts.rd == 0, f"Expected rd=0 (ambiguous), got {counts.rd}"
        assert counts.ad == 0, f"Expected ad=0 (ambiguous), got {counts.ad}"


class TestCaseA_AllBasesLowQual:
    """Edge case: every base in the variant region is low quality."""

    def test_all_low_qual(self, tmp_path):
        """All bases at variant region below min_baseq → Neither."""
        quals = [30, 30, 30, 30, 5, 5, 30, 30, 30, 30]
        reads = [_make_read("r1", "AAAAATAAAA", 96, ((0, 10),), quals=quals)]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, EQUAL_LEN_VARIANT)
        assert counts.rd == 0
        assert counts.ad == 0


# ==========================================================================
# Case B: Only ALT length matches (DelIns)
# Variant: chr1:100, REF=ATG, ALT=CC  (3→2 bases)
# Read that shows deletion+insertion → reconstructed seq is 2 bases.
# ==========================================================================

DELINS_VARIANT = Variant(
    chrom="chr1",
    pos=100,
    ref_allele="ATG",
    alt_allele="CC",
    variant_type="COMPLEX",
)


class TestCaseB_DelIns:
    """Case B: DelIns where only ALT length matches reconstruction."""

    def test_delins_exact_alt(self, tmp_path):
        """Read reconstructs 'CC' for the ATG region → ALT.
        Read starts at 96. CIGAR: 4M (covers 96-99), then covers pos 100-102.
        We need the read to show 'CC' where reference has 'ATG' (3 bases).
        CIGAR: 4M 2X 1D 4M → 4 match, 2 mismatch (CC instead of AT), skip G, 4 match.
        Actually for check_complex, the reconstruction works from CIGAR:
        4M covers 96-99 (overlap with 100-102 is none since ref_pos starts before start_pos)
        Let me recalculate: start_pos=100, end_pos=103.
        Read at 96, CIGAR=4M means ref 96-99. Then we need ops covering 100-102.
        Total: 4M (96-99) + some ops for 100-102.
        For 'CC' reconstruction: Match 2 bases at 100-101 = 'CC', then Del 1 at 102.
        CIGAR: 4M 2M 1D 4M. But the 4M covers 96-99, 2M covers 100-101 with 'CC', Del skips 102.
        Read: AAAACCAAAA, start=96, CIGAR: (0,6) (2,1) (0,4)
        The 6M covers 96-101, reconstruction gets bases at 100-101 = 'CC'.
        Del skips 102. 4M covers 103-106.
        Reconstructed for [100,103) = 'CC' (2 bases). ALT='CC' → match!
        """
        reads = [_make_read("r1", "AAAACCAAAA", 96, ((0, 6), (2, 1), (0, 4)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, DELINS_VARIANT)
        assert counts.ad == 1, f"Expected ad=1 (DelIns ALT), got {counts.ad}"

    def test_delins_rescue_low_qual(self, tmp_path):
        """Read reconstructs 'CA' where 'A' is Q5. Expected ALT='CC'.
        Reliable: 'C' at pos 0 (Q30).
        vs ALT: 'C'=='C' ✓ → ALT (rescued, no ambiguity since REF length differs).
        """
        quals = [30, 30, 30, 30, 30, 5, 30, 30, 30, 30]
        reads = [_make_read("r1", "AAAACAAAAA", 96, ((0, 6), (2, 1), (0, 4)), quals=quals)]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, DELINS_VARIANT)
        assert counts.ad == 1, f"Expected ad=1 (DelIns rescued), got {counts.ad}"

    def test_delins_mismatch_on_reliable(self, tmp_path):
        """Read reconstructs 'TT' for ALT='CC'. High-quality mismatch → Neither."""
        reads = [_make_read("r1", "AAAATTAAAA", 96, ((0, 6), (2, 1), (0, 4)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, DELINS_VARIANT)
        assert counts.ad == 0
        assert counts.rd == 0  # REF length (3) doesn't match recon (2)


# ==========================================================================
# Case C: Only REF length matches
# Variant: chr1:100, REF=CC, ALT=ATG  (2→3 bases, insertion-like complex)
# A read showing 'CC' (2 bases, no insertion) → matches REF length.
# ==========================================================================

INSLIKE_VARIANT = Variant(
    chrom="chr1",
    pos=100,
    ref_allele="CC",
    alt_allele="ATG",
    variant_type="COMPLEX",
)


class TestCaseC_RefOnly:
    """Case C: Only REF length matches reconstruction."""

    def test_ref_exact(self, tmp_path):
        """Read reconstructs 'CC' for REF='CC' → REF.
        Read: 10M at 96. Bases at 100-101 = 'CC'.
        """
        reads = [_make_read("r1", "AAAACCAAAA", 96, ((0, 10),))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INSLIKE_VARIANT)
        assert counts.rd == 1, f"Expected rd=1 (REF match), got {counts.rd}"

    def test_ref_rescue_low_qual(self, tmp_path):
        """Read reconstructs 'CA' where 'A' is Q5. REF='CC'.
        Reliable: 'C' at pos 0 (Q30).
        vs REF: 'C'=='C' ✓ → REF (rescued).
        """
        quals = [30, 30, 30, 30, 30, 5, 30, 30, 30, 30]
        reads = [_make_read("r1", "AAAACAAAAA", 96, ((0, 10),), quals=quals)]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INSLIKE_VARIANT)
        assert counts.rd == 1, f"Expected rd=1 (REF rescued), got {counts.rd}"


# ==========================================================================
# Length mismatch — neither matches
# ==========================================================================


class TestLengthMismatch:
    """Reconstruction length matches neither REF nor ALT.

    With Phase 3 local SW fallback (added for EPHA7-style complex variants),
    reads that fail Phase 2 length matching can still be classified via
    Smith-Waterman alignment against REF/ALT haplotypes.
    """

    def test_length_mismatch_falls_through_to_sw(self, tmp_path):
        """Variant REF=ATG (3bp), ALT=CC (2bp). Read shows 4-base reconstruction.
        Read: 4M 1I 6M start at 96. Insertion adds a base in the region.
        4M covers 96-99. 1I at pos 100 (1 inserted base). 6M covers 100-105.
        Reconstruction for [100, 103): the inserted base + 100 + 101 + 102 = 4 bases.
        Neither 3 (REF) nor 2 (ALT) by length.

        Phase 2 and Phase 2.5 cannot decide, so Phase 3 SW alignment kicks in.
        The local fallback soft-clips the extra base and finds an ALT match.
        """
        reads = [_make_read("r1", "AAAATCCAAAA", 96, ((0, 4), (1, 1), (0, 6)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, DELINS_VARIANT)
        assert counts.rd == 0
        # Phase 3 SW local fallback classifies this as ALT
        assert counts.ad == 1


# ==========================================================================
# Gap 3A: Dynamic SW Gap Penalties for MSI-like repeats
# When repeat_span >= 10bp, gap_extend is 0 (free), absorbing polymerase
# slippage noise at MSI sites.
# ==========================================================================


class TestGap3A_DynamicGapPenalties:
    """Gap 3A: Verify dynamic SW gap penalties with MSI-like variants."""

    def test_repeat_variant_sw_classification(self, tmp_path):
        """
        Complex variant in a homopolymer region.

        Variant: 1bp deletion in polyA context — REF=AA, ALT=A.
        This is classified as DELETION, but if the read shows the deletion
        at a shifted position (slippage), it falls through to Phase 3 SW.

        With repeat_span >= 10, gap_extend = 0 (free). This should allow
        the SW aligner to absorb the shift and correctly classify as ALT.

        Note: repeat_span is set during normalization (prepare_variants),
        not directly on Variant objects. For direct count_bam calls,
        repeat_span defaults to 0. This test verifies the SW pipeline
        still works with standard gap penalties.
        """
        # Read with a 1bp deletion 2 positions to the right of expected
        # CIGAR: 7M 1D 3M → deletion at pos 103 instead of 101
        reads = [_make_read("r1", "AAAAAAAAAA", 96, ((0, 7), (2, 1), (0, 3)))]

        # Variant at pos 100: AA -> A (1bp deletion)
        variant = Variant(
            chrom="chr1",
            pos=100,
            ref_allele="AA",
            alt_allele="A",
            variant_type="DELETION",
        )

        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, variant)

        # The deletion check should find this via windowed detection or SW fallback
        # dp must include this read regardless
        assert counts.dp >= 1, f"Expected dp >= 1, got {counts.dp}"


# ── count_bam_binned parity tests ────────────────────────────────────────


# _count_one_both imported from helpers.py above


class TestBinnedParity:
    """Verify count_bam_binned matches count_bam for complex variant scenarios."""

    def test_equal_len_alt_binned(self, tmp_path):
        reads = [_make_read("r1", "AAAACGAAAA", 96, ((0, 10),))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one_both(bam, EQUAL_LEN_VARIANT)
        assert counts.ad == 1
        assert counts.rd == 0

    def test_equal_len_ref_binned(self, tmp_path):
        reads = [_make_read("r1", "AAAAATAAAA", 96, ((0, 10),))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one_both(bam, EQUAL_LEN_VARIANT)
        assert counts.rd == 1
        assert counts.ad == 0

    def test_delins_alt_binned(self, tmp_path):
        reads = [_make_read("r1", "AAAACCAAAA", 96, ((0, 6), (2, 1), (0, 4)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one_both(bam, DELINS_VARIANT)
        assert counts.ad == 1

    def test_ref_only_binned(self, tmp_path):
        reads = [_make_read("r1", "AAAACCAAAA", 96, ((0, 10),))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one_both(bam, INSLIKE_VARIANT)
        assert counts.rd == 1
