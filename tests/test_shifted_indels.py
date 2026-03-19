"""
Tests for windowed indel detection (Phase 2a).

Verifies that insertions and deletions shifted by the aligner within ±5bp of the
expected anchor position are still correctly detected, using three safeguards:
  S1: Inserted/deleted sequence must match expected bases exactly
  S2: Closest match wins when multiple candidates exist
  S3: Reference context at shifted position must be biologically valid

Reference genome: chr1, 500 bases of 'A' (all A's for homopolymer tests).
  - Insertion tests: pos 100, REF=A, ALT=AT (insert T after anchor A)
  - Deletion tests:  pos 200, REF=AT, ALT=A (delete T after anchor A)
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

from gbcms._rs import Variant, count_bam  # noqa: F401 — count_bam referenced in comments

# ==========================================================================
# INSERTION TESTS — Variant: chr1:100, REF=A, ALT=AT
# Anchor base is 'A' at pos 100. Insertion of 'T' after anchor.
# ==========================================================================

# ref_context: reference bases [95, 110) = "AAAAAAAAAAAAAAA" (all A's)
INS_VARIANT = Variant(
    chrom="chr1",
    pos=100,
    ref_allele="A",
    alt_allele="AT",
    variant_type="INSERTION",
    ref_context="AAAAAAAAAAAAAAA",  # [95, 110)
    ref_context_start=95,
)


class TestInsertionStrict:
    """Insertion at the exact expected position (strict fast path)."""

    def test_strict_match(self, tmp_path):
        """INS immediately after anchor at pos 100 → ALT."""
        # Read: 5M 1I 4M, starts at 96. Anchor at index 4 → pos 100.
        reads = [_make_read("r1", "AAAAATAAAA", 96, ((0, 5), (1, 1), (0, 4)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INS_VARIANT)
        assert counts.ad == 1, f"Expected ad=1 (strict ALT), got {counts.ad}"
        assert counts.rd == 0

    def test_ref_no_insertion(self, tmp_path):
        """Read covers anchor with no insertion → REF."""
        reads = [_make_read("r1", "AAAAAAAAAA", 96, ((0, 10),))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INS_VARIANT)
        assert counts.rd == 1, f"Expected rd=1 (REF), got {counts.rd}"
        assert counts.ad == 0


class TestInsertionWindowed:
    """Insertion shifted within ±5bp of the expected anchor position."""

    def test_shifted_right_2bp(self, tmp_path):
        """INS 2bp to the right of anchor → ALT (windowed).
        Anchor is at pos 100. Insertion placed at pos 103 (after M block ending at 102).
        Read starts at 96, 7M 1I 2M. Block end at 103, ins after base at 102.
        Shifted anchor at 102 is 'A' → passes S3.
        """
        reads = [_make_read("r1", "AAAAAAATAA", 96, ((0, 7), (1, 1), (0, 2)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INS_VARIANT)
        assert counts.ad == 1, f"Expected ad=1 (windowed ALT), got {counts.ad}"

    def test_shifted_left_3bp(self, tmp_path):
        """INS 3bp to the left of anchor → ALT (windowed).
        Anchor at 100. Insertion placed at pos 98 (after M block).
        Read starts at 96, 2M 1I 7M. Block end at 98.
        Shifted anchor at 97 is 'A' → passes S3.
        """
        reads = [_make_read("r1", "AATAAAAAAA", 96, ((0, 2), (1, 1), (0, 7)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INS_VARIANT)
        assert counts.ad == 1, f"Expected ad=1 (windowed ALT, shifted left), got {counts.ad}"

    def test_wrong_inserted_sequence(self, tmp_path):
        """INS within window but wrong inserted base → REF (S1 rejects).
        Read has insertion of 'G' instead of 'T' near anchor.
        """
        reads = [_make_read("r1", "AAAAAGAAAA", 96, ((0, 5), (1, 1), (0, 4)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INS_VARIANT)
        # The read covers the anchor so it should be counted as REF
        assert counts.rd == 1, f"Expected rd=1 (wrong seq → REF), got {counts.rd}"
        assert counts.ad == 0

    def test_outside_window(self, tmp_path):
        """INS 8bp away from anchor → not detected (outside ±5bp window).
        Read starts at 96, 12M 1I 2M. Insertion at pos 108, well outside window.
        """
        reads = [_make_read("r1", "AAAAAAAAAAAAAT" + "A", 96, ((0, 12), (1, 1), (0, 2)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INS_VARIANT)
        # Read still covers pos 100 with a match, so it's REF
        assert counts.rd == 1, f"Expected rd=1 (outside window → REF), got {counts.rd}"
        assert counts.ad == 0

    def test_anchor_base_mismatch_s3(self, tmp_path):
        """INS within window but reference base at shifted anchor differs → reject (S3).
        Use a ref_context where the base at the shifted position is NOT 'A'.
        """
        # Custom variant with non-homopolymer ref_context
        # ref_context: "AAAAACGAAAA" covering [95, 106)
        # Pos 100 = 'C' (offset 5), Pos 101 = 'G' (offset 6)
        # Anchor base is 'A' (from REF allele). Shifted anchor at pos 100 would be 'C' → S3 reject.
        variant = Variant(
            chrom="chr1",
            pos=98,
            ref_allele="A",
            alt_allele="AT",
            variant_type="INSERTION",
            ref_context="AAAAACGAAAA",  # [93, 104)
            ref_context_start=93,
        )
        # Read has insertion at pos 101 (shifted +2 from anchor at 98).
        # Shifted anchor at pos 100 → ref_context[100-93]=ref_context[7]='A' → actually passes.
        # Need to set up so the shifted anchor is 'C' at offset 5 (pos 98).
        # Let's make the insertion at pos 99 instead, then shifted anchor is at 98 which is 'C' at offset 5.
        # No wait — let me think more carefully.
        # Variant anchor at pos 98. ref_context starts at 93: "AAAAACGAAAA"
        # ref_context[0..11] = A(93) A(94) A(95) A(96) A(97) C(98) G(99) A(100) A(101) A(102) A(103)
        # Original anchor base from REF: 'A'
        # For an insertion shifted to pos 101 (block_end=101), shifted_anchor = 100 → 'A' → passes
        # For an insertion shifted to pos 100 (block_end=100), shifted_anchor = 99 → 'G' → FAILS S3!
        # Read starts at 93, 7M 1I 3M. Block end = 100. Ins at 100. Shifted anchor = 99 → 'G' ≠ 'A'
        reads = [_make_read("r1", "AAAAAAATAAA", 93, ((0, 7), (1, 1), (0, 3)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, variant)
        # S3 rejects because reference at pos 99 is 'G', not 'A' (the anchor base).
        # The read's CIGAR shows an insertion near the variant → ambiguous classification
        # routes to "neither" (DP only), preserving unbiased VAF.
        assert counts.ad == 0, f"Expected ad=0 (S3 reject), got {counts.ad}"
        assert counts.dp >= 1, f"Expected dp>=1 (read covers locus), got {counts.dp}"


class TestInsertionHomopolymer:
    """Insertions in homopolymer runs — the primary use case for windowed detection."""

    def test_homopolymer_shifted(self, tmp_path):
        """INS in AAAA run shifted by aligner → ALT.
        In an all-A reference, any T insertion within ±5bp should be detected.
        """
        # Insertion shifted +4bp from anchor
        reads = [_make_read("r1", "AAAAAAAAATAA", 96, ((0, 9), (1, 1), (0, 2)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INS_VARIANT)
        assert counts.ad == 1, f"Expected ad=1 (homopolymer shift), got {counts.ad}"


# ==========================================================================
# DELETION TESTS — Variant: chr1:200, REF=AT, ALT=A
# Anchor base 'A' at pos 200. Delete 'T' at pos 201.
# ==========================================================================

# For deletion, ref_context covers [195, 210). In an all-A reference, the deleted
# 'T' won't match at shifted positions. We need a ref_context with T's so S3 passes.
# Let's use a realistic ref_context: "AAAAATATATAAA" covering [195, 208)
# This has T at positions 200, 202, 204 (offsets 5, 7, 9)
# Actually, let's keep it simpler. The expected deleted base is 'T' (REF[1:] = "T").
# For windowed to pass S3, the ref base at the shifted deletion position must also be 'T'.
# ref_context: "AAAAATTTTAAAAAA" covering [195, 210)
# Positions: A(195) A(196) A(197) A(198) A(199) T(200) T(201) T(202) T(203) A(204) ...
# Wait, the ref_allele is "AT", meaning pos 200 = 'A' (anchor) and pos 201 = 'T' (deleted).
# In the all-A reference genome, pos 201 should NOT be 'T', so S3 would reject.
# We need a custom ref_context that matches the biology.

DEL_VARIANT = Variant(
    chrom="chr1",
    pos=200,
    ref_allele="AT",
    alt_allele="A",
    variant_type="DELETION",
    # ref_context covers [195, 210) — reference has T's at 201, 202, 203 so windowed shifts pass S3
    ref_context="AAAAATTTTAAAAAA",  # A(195)..A(199) T(200)? No...
    # Actually: pos 200=A (anchor), 201=T, 202=T, 203=T, 204=A...
    # Re-derive: indices [195, 210) → 15 chars
    # For the variant to have anchor 'A' at 200 and delete 'T' at 201,
    # the reference must have 'A' at 200 and 'T' at 201.
    # For windowed S3, shifted deletions at 202, 203 need 'T' there too.
    # ref_context = "AAAAATTTTTAAAAA" → A(195..199) T(200..204) A(205..209)
    # Wait, pos 200 should be 'A' (anchor)!
    # ref_context = "AAAAAATTTTTAAAA" → A(195..200) T(201..205) A(206..209)
    ref_context_start=195,
)

# Fix: rebuild DEL_VARIANT with correct ref_context
# Positions 195-209: A A A A A A T T T T A A A A A
#                    195         200 201 202 203 204
# pos 200='A' (anchor), 201='T' (deleted), 202='T', 203='T'
DEL_VARIANT = Variant(
    chrom="chr1",
    pos=200,
    ref_allele="AT",
    alt_allele="A",
    variant_type="DELETION",
    ref_context="AAAAAATTTTTAAAA",  # [195, 210): 6 A's + 5 T's + 4 A's = 15 chars
    ref_context_start=195,
)


class TestDeletionStrict:
    """Deletion at the exact expected position (strict fast path)."""

    def test_strict_match(self, tmp_path):
        """DEL immediately after anchor at pos 200 → ALT.
        Read: 5M 1D 5M, starts at 196. Anchor at index 4 → pos 200.
        """
        reads = [_make_read("r1", "AAAAAAAAAA", 196, ((0, 5), (2, 1), (0, 5)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, DEL_VARIANT)
        assert counts.ad == 1, f"Expected ad=1 (strict ALT), got {counts.ad}"
        assert counts.rd == 0

    def test_ref_no_deletion(self, tmp_path):
        """Read covers anchor without deletion → REF."""
        reads = [_make_read("r1", "AAAAAAAAAA", 196, ((0, 10),))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, DEL_VARIANT)
        assert counts.rd == 1, f"Expected rd=1 (REF), got {counts.rd}"
        assert counts.ad == 0


class TestDeletionWindowed:
    """Deletion shifted within ±5bp of the expected anchor position."""

    def test_shifted_right_1bp(self, tmp_path):
        """DEL 1bp to the right → ALT (windowed).
        Anchor at 200. Deletion at pos 203 instead of 201.
        Read: 8M 1D 2M, starts at 196. Block end = 204.
        Wait, let me recalculate. starts at 196, 8M → covers 196-203.
        Block end = 204. Del at 204?  That's +3bp from 201.
        Actually, deletion happens at block_end (204). Expected del at anchor+1 = 201.
        del_ref_pos = 204, distance = |204 - 201| = 3.

        Hmm, let me reconsider. We need deletion within ±5bp of anchor (200).
        window is [195, 205]. del_ref_pos = 204 → within window.
        S1: del_len = 1 = expected. ✓
        S3: ref at 204 = ref_context[204-195] = ref_context[9] = 'T' (5th T). expected_del = 'T'. ✓
        """
        reads = [_make_read("r1", "AAAAAAAAAA", 196, ((0, 8), (2, 1), (0, 2)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, DEL_VARIANT)
        assert counts.ad == 1, f"Expected ad=1 (windowed ALT), got {counts.ad}"

    def test_deletion_wrong_ref_context_s3(self, tmp_path):
        """DEL within window but ref base at shifted pos doesn't match → S3 rejects
        the windowed CIGAR match. The read covers the anchor with a Match op and
        no deletion at the expected position → REF.

        Previously this fell through to Phase 3 SW via is_worth_realignment(),
        but CIGAR is definitive for pure deletions: if the Match op covers the
        anchor and no matching D op was found, the read is REF.
        """
        # ref_context: "AAAAAATGAAAAAAAA" covering [195, 210)
        # pos 201='T', 202='G'. If deletion shifts to 202, ref='G' but expected='T' → S3 reject.
        variant = Variant(
            chrom="chr1",
            pos=200,
            ref_allele="AT",
            alt_allele="A",
            variant_type="DELETION",
            ref_context="AAAAAATGAAAAAAA",  # [195, 210): T at 201, G at 202
            ref_context_start=195,
        )
        # Read with deletion at pos 203: 7M 1D 3M, starts at 196.
        # Block end = 203. del_ref_pos = 203.
        # ref_context[203-195] = ref_context[8] = 'A'. expected_del_seq = 'T'. 'A' ≠ 'T' → S3 reject.
        # The read has a Match op covering the anchor (pos 200) → found_ref_coverage=true.
        # CIGAR is definitive for pure DEL: no matching D at anchor → REF.
        reads = [_make_read("r1", "AAAAAAAAAA", 196, ((0, 7), (2, 1), (0, 3)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, variant)
        # CIGAR definitive: read covers anchor with Match, no deletion → REF
        assert counts.rd == 1, f"Expected rd=1 (CIGAR definitive REF), got {counts.rd}"


class TestNoRefContext:
    """When ref_context is None, S3 cannot validate — shifted indels are rejected (fail-closed)."""

    def test_insertion_no_ref_context(self, tmp_path):
        """Windowed insertion without ref_context → S3 rejects (fail-closed safety)."""
        variant = Variant(
            chrom="chr1",
            pos=100,
            ref_allele="A",
            alt_allele="AT",
            variant_type="INSERTION",
        )
        # Shifted insertion +2bp from anchor
        reads = [_make_read("r1", "AAAAAAATAA", 96, ((0, 7), (1, 1), (0, 2)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, variant)
        # S3 cannot validate without ref_context → reject shifted indel
        assert counts.ad == 0, f"Expected ad=0 (no ref_context → S3 reject), got {counts.ad}"

    def test_deletion_no_ref_context(self, tmp_path):
        """Windowed deletion without ref_context → S3 rejects (fail-closed safety)."""
        variant = Variant(
            chrom="chr1",
            pos=200,
            ref_allele="AT",
            alt_allele="A",
            variant_type="DELETION",
        )
        # Shifted deletion +2bp
        reads = [_make_read("r1", "AAAAAAAAAA", 196, ((0, 7), (2, 1), (0, 3)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, variant)
        # S3 cannot validate without ref_context → reject shifted indel
        assert counts.ad == 0, f"Expected ad=0 (no ref_context → S3 reject), got {counts.ad}"


class TestReadDoesNotCover:
    """Reads that do not cover the variant region → neither REF nor ALT."""

    def test_read_far_away(self, tmp_path):
        """Read at distant position → not counted."""
        reads = [_make_read("r1", "AAAAAAAAAA", 50, ((0, 10),))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INS_VARIANT)
        assert counts.rd == 0
        assert counts.ad == 0
        assert counts.dp == 0


# ── count_bam_binned parity tests ────────────────────────────────────────


# _count_one_both imported from helpers.py above


class TestBinnedParity:
    """Verify count_bam_binned matches count_bam for shifted indel scenarios."""

    def test_insertion_strict_binned(self, tmp_path):
        reads = [_make_read("r1", "AAAAATAAAA", 96, ((0, 5), (1, 1), (0, 4)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one_both(bam, INS_VARIANT)
        assert counts.ad == 1
        assert counts.rd == 0

    def test_insertion_windowed_binned(self, tmp_path):
        reads = [_make_read("r1", "AAAAAAATAA", 96, ((0, 7), (1, 1), (0, 2)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one_both(bam, INS_VARIANT)
        assert counts.ad == 1

    def test_deletion_strict_binned(self, tmp_path):
        reads = [_make_read("r1", "AAAAAAAAAA", 196, ((0, 5), (2, 1), (0, 5)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one_both(bam, DEL_VARIANT)
        assert counts.ad == 1
        assert counts.rd == 0
