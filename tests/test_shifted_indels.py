"""
Tests for windowed indel detection (Phase 2a).

Verifies that insertions and deletions shifted by the aligner within ±5bp of the
expected anchor position are still correctly detected, using three safeguards:
  S1: Inserted/deleted length must match the expected event
  S2: Closest match wins when multiple candidates exist
  S3: The shifted placement gives the variant's haplotype (the same sequence
      elsewhere in a repeat; #189)

The reference is each variant's ref_context (chr1, 500bp), and the reads carry its bases.
  - Insertion tests: pos 100 in an A run, REF=A, ALT=AA (an A inserted into the A
    run: any placement is the same haplotype) and ALT=AT (a T: each placement is its
    own haplotype)
  - Deletion tests:  pos 200, REF=AT, ALT=A (one T of the T run after anchor A)
"""

from helpers import (
    build_bam as _build_bam,
)
from helpers import (
    count_one_checked as _count_one,
)
from helpers import (
    make_read as _make_read,
)

from gbcms._rs import Variant

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


# The same run, an A inserted: every placement in the run gives this haplotype.
INS_A_VARIANT = Variant(
    chrom="chr1",
    pos=100,
    ref_allele="A",
    alt_allele="AA",
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
        """An A inserted 2bp to the right of the anchor, in a read that lies inside
        the A run the variant's context holds (start 96, 7M 1I 2M): its bases read
        the same with or without the A, so it counts toward depth only. (Carriers
        that span the run count ALT at every junction: the shifted-carrier battery.)"""
        reads = [_make_read("r1", "AAAAAAAAAA", 96, ((0, 7), (1, 1), (0, 2)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INS_A_VARIANT)
        assert (counts.ad, counts.rd, counts.dp) == (0, 0, 1)

    def test_shifted_left_3bp(self, tmp_path):
        """An A inserted 3bp to the left of the anchor (start 96, 2M 1I 7M), the read
        inside the run: depth only, as above."""
        reads = [_make_read("r1", "AAAAAAAAAA", 96, ((0, 2), (1, 1), (0, 7)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INS_A_VARIANT)
        assert (counts.ad, counts.rd, counts.dp) == (0, 0, 1)

    def test_shifted_other_haplotype_is_not_alt(self, tmp_path):
        """A T inserted 2bp right of where the variant puts its T: in the A run
        that is another haplotype (a T elsewhere), never ALT (#189: the anchor-base
        test used to accept it). It lies outside the variant's discrimination
        window [100, 102), which the read shows as reference: REF."""
        reads = [_make_read("r1", "AAAAAAATAA", 96, ((0, 7), (1, 1), (0, 2)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INS_VARIANT)
        assert counts.ad == 0, f"Expected ad=0 (another haplotype), got {counts.ad}"
        assert (counts.rd, counts.partial_alt) == (1, 0)

    def test_wrong_inserted_sequence(self, tmp_path):
        """Same-length insertion AT the anchor with a confidently wrong base
        ('G' at Q30 instead of expected 'T') is a third allele: neither REF
        (the read provably carries an insertion — rd must not absorb it) nor
        the queried ALT, and it surfaces as partial evidence. The read reads past
        the A run to the C's after it, so its bases show the insertion: ending
        inside the run, they would fit REF with a substitution as well.
        """
        variant = Variant(
            chrom="chr1",
            pos=100,
            ref_allele="A",
            alt_allele="AT",
            variant_type="INSERTION",
            ref_context="AAAAAAAAAACCCCC",  # [95, 110): the A run ends at 104
            ref_context_start=95,
        )
        reads = [_make_read("r1", "AAAAAGAAAACC", 96, ((0, 5), (1, 1), (0, 6)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, variant)
        assert counts.rd == 0, f"Expected rd=0 (third allele, not REF), got {counts.rd}"
        assert counts.ad == 0
        assert counts.partial_alt == 1, (
            f"Expected partial_alt=1 (wrong-sequence insertion surfaced), "
            f"got {counts.partial_alt}"
        )

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

    def test_shifted_placement_giving_another_haplotype_s3(self, tmp_path):
        """INS within the window whose placement gives another haplotype → not ALT
        (S3). The variant's prepared reference [93, 104) is AAAAACGAAAA (it differs
        from the BAM's all-A reference on purpose): the variant puts its T before 99
        (C T G), the read before 100 (C G T).
        """
        variant = Variant(
            chrom="chr1",
            pos=98,
            ref_allele="A",
            alt_allele="AT",
            variant_type="INSERTION",
            ref_context="AAAAACGAAAA",  # [93, 104)
            ref_context_start=93,
        )
        # Read starts at 93, 7M 1I 3M: its T sits before 100, one junction right of
        # the variant's (before 99). With S = ref[99] = G: X + S = TG, S + Y = GT.
        reads = [_make_read("r1", "AAAAAAATAAA", 93, ((0, 7), (1, 1), (0, 3)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, variant)
        assert counts.ad == 0, f"Expected ad=0 (S3 reject), got {counts.ad}"
        assert counts.dp >= 1, f"Expected dp>=1 (read covers locus), got {counts.dp}"


class TestInsertionHomopolymer:
    """Insertions in homopolymer runs — the primary use case for windowed detection."""

    def test_homopolymer_shifted(self, tmp_path):
        """An A inserted into the A run, shifted +4bp by the aligner, in a read
        inside the run: depth only (its bases fit both alleles)."""
        reads = [_make_read("r1", "AAAAAAAAAAAA", 96, ((0, 9), (1, 1), (0, 2)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, INS_A_VARIANT)
        assert (counts.ad, counts.rd, counts.dp) == (0, 0, 1)


# ==========================================================================
# DELETION TESTS — Variant: chr1:200, REF=AT, ALT=A
# Anchor base 'A' at pos 200. Delete 'T' at pos 201.
# ==========================================================================

# The reference is the variant's ref_context, [195, 210): A at 195-200 (the anchor at
# 200), a T run at 201-205, A at 206-209. The reads carry its bases.
DEL_CONTEXT = "AAAAAATTTTTAAAA"
DEL_VARIANT = Variant(
    chrom="chr1",
    pos=200,
    ref_allele="AT",
    alt_allele="A",
    variant_type="DELETION",
    ref_context=DEL_CONTEXT,
    ref_context_start=195,
)


def _ctx(context: str, lo: int, hi: int, start: int = 195) -> str:
    """The reference bases [lo, hi) of a context starting at `start`."""
    return context[lo - start : hi - start]


class TestDeletionStrict:
    """Deletion at the exact expected position (strict fast path)."""

    def test_strict_match(self, tmp_path):
        """DEL immediately after anchor at pos 200 → ALT.
        Read: 5M 1D 5M from 196: the anchor, the T run one T short, its right flank A.
        """
        seq = _ctx(DEL_CONTEXT, 196, 201) + _ctx(DEL_CONTEXT, 202, 207)
        reads = [_make_read("r1", seq, 196, ((0, 5), (2, 1), (0, 5)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, DEL_VARIANT)
        assert counts.ad == 1, f"Expected ad=1 (strict ALT), got {counts.ad}"
        assert counts.rd == 0

    def test_ref_no_deletion(self, tmp_path):
        """Read spans the T-run (201-205) and its right flank without a deletion → REF.

        Covering the anchor alone is not enough: a read ending inside the run
        reads the same with or without one T, so it is uninformative.
        """
        reads = [_make_read("r1", _ctx(DEL_CONTEXT, 196, 207), 196, ((0, 11),))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, DEL_VARIANT)
        assert counts.rd == 1, f"Expected rd=1 (REF), got {counts.rd}"
        assert counts.ad == 0


class TestDeletionWindowed:
    """Deletion shifted within ±5bp of the expected anchor position."""

    def test_shifted_right_1bp(self, tmp_path):
        """DEL written at 204, the run's fourth T, instead of 201 → ALT (windowed).
        Read: 8M 1D 2M from 196: the same bases as a deletion written at 201, the T run
        one T short and its right flank A.
        """
        seq = _ctx(DEL_CONTEXT, 196, 204) + _ctx(DEL_CONTEXT, 205, 207)
        reads = [_make_read("r1", seq, 196, ((0, 8), (2, 1), (0, 2)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, DEL_VARIANT)
        assert counts.ad == 1, f"Expected ad=1 (windowed ALT), got {counts.ad}"

    def test_deletion_wrong_ref_context_s3(self, tmp_path):
        """DEL within window but ref base at shifted pos doesn't match → S3 rejects
        the windowed CIGAR match. The read covers the anchor with a Match op and
        no deletion at the expected position → REF (CIGAR definitive).

        Deletions of 5bp or more that fail S3 are judged by the read's bases (50bp
        or more at the POS breakpoint, RJ-24); short deletions (1-4bp) failing S3 are
        almost certainly spurious noise, so CIGAR remains definitive for them. A 1bp S3-rejected windowed
        Del does NOT trigger Phase 3; found_ref_coverage=True → REF.
        """
        # ref_context: T at 201, G at 202, an A run at 203-209
        context = "AAAAAATGAAAAAAA"
        variant = Variant(
            chrom="chr1",
            pos=200,
            ref_allele="AT",
            alt_allele="A",
            variant_type="DELETION",
            ref_context=context,
            ref_context_start=195,
        )
        # Read with deletion at pos 203: 7M 1D 3M from 196, one A of the A run deleted.
        # It ends inside the run, so its bases read as the reference's: the T at 201
        # is there → REF.
        seq = _ctx(context, 196, 203) + _ctx(context, 204, 207)
        reads = [_make_read("r1", seq, 196, ((0, 7), (2, 1), (0, 3)))]
        bam = _build_bam(tmp_path, reads)
        counts = _count_one(bam, variant)
        # CIGAR definitive: 1bp del, S3 reject, short del → REF
        assert counts.rd == 1, f"Expected rd=1 (CIGAR definitive REF), got {counts.rd}"
        assert counts.ad == 0


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
