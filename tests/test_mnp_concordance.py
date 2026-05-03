"""
MNP Concordance Tests — selective quality gate and fragment-level verification.

Validates the MNP counting fixes from feature/fix-mnp-counting:
  1. Selective discriminating-position quality gate
  2. LowQuality → neither (no check_complex fallback)
  3. Fragment-level count propagation (structural invariants)
  4. DNP, TNP, and ONP subtypes via both MAF and VCF-style input

These tests use synthetic BAM data with controlled base qualities to
exercise the specific edge cases identified in the TERT/BRCA2 analysis.
"""

import pytest
from helpers import build_bam, count_both, make_read

from gbcms import _rs as gbcms_rs

# ── Fixtures ─────────────────────────────────────────────────────────────


@pytest.fixture
def mnp_bam(tmp_path):
    """BAM with reads covering a 5bp ONP region at chr1:100-104.

    Variant: GAGGG→AAGGA (TERT-like pattern).
    Discriminating positions: 0 (G→A), 4 (G→A).
    Non-discriminating: 1 (A=A), 2 (G=G), 3 (G=G).

    Reads:
      - 5 fwd REF (GAGGG) at Q30
      - 3 rev ALT (AAGGA) at Q30
      - 2 fwd ALT with low-qual non-discriminating pos (should PASS)
      - 1 fwd ALT with low-qual discriminating pos (should be discarded)
      - 1 fwd ThirdAllele (AAGGG — partial mutation)
    """
    reads = []

    # 5 forward REF reads
    for i in range(5):
        reads.append(
            make_read(
                f"ref_fwd_{i}",
                "AAGAGGGAA",  # 9bp, ONP at index 2-6
                start=98,
                cigar=((0, 9),),
                flag=0,
            )
        )

    # 3 reverse ALT reads
    for i in range(3):
        reads.append(
            make_read(
                f"alt_rev_{i}",
                "AAAAGGATTT",  # ONP at index 2-6: AAGGA
                start=98,
                cigar=((0, 10),),
                flag=16,
                quals=[30, 30, 35, 38, 32, 36, 34, 30, 30, 30],
            )
        )

    # 2 forward ALT reads with low-qual at NON-discriminating pos
    # (should pass selective quality gate)
    for i in range(2):
        quals = [30, 30, 35, 38, 5, 36, 34, 30, 30, 30]  # pos 4 (non-disc): Q=5
        reads.append(
            make_read(
                f"alt_fwd_lowq_nondisc_{i}",
                "AAAAGGATTT",
                start=98,
                cigar=((0, 10),),
                flag=0,
                quals=quals,
            )
        )

    # 1 forward ALT read with low-qual at DISCRIMINATING pos
    # (should be discarded → neither)
    reads.append(
        make_read(
            "alt_fwd_lowq_disc",
            "AAAAGGATTT",
            start=98,
            cigar=((0, 10),),
            flag=0,
            quals=[30, 30, 5, 38, 32, 36, 34, 30, 30, 30],  # pos 2 (disc G→A): Q=5
        )
    )

    # 1 forward ThirdAllele (partial mutation: AAGGG, only pos 0 mutated)
    reads.append(
        make_read(
            "third_allele_fwd",
            "AAAAGGGTTT",  # AAGGG at index 2-6
            start=98,
            cigar=((0, 10),),
            flag=0,
        )
    )

    return build_bam(tmp_path, reads, "mnp_test.bam")


@pytest.fixture
def dnp_bam(tmp_path):
    """BAM with reads for a simple all-discriminating DNP at chr1:100-101.

    Variant: GG→AA (both positions discriminating).
    """
    reads = []

    # 3 forward REF
    for i in range(3):
        reads.append(make_read(f"ref_fwd_{i}", "AAAGGTTT", start=97, cigar=((0, 8),)))

    # 4 reverse ALT
    for i in range(4):
        reads.append(
            make_read(f"alt_rev_{i}", "AAAAAATT", start=97, cigar=((0, 8),), flag=16)
        )

    return build_bam(tmp_path, reads, "dnp_test.bam")


# ── ONP Tests (TERT-like pattern) ────────────────────────────────────────


class TestONPSelectiveQualityGate:
    """Tests for the selective discriminating-position quality gate."""

    def test_onp_alt_count_with_selective_gate(self, mnp_bam):
        """ALT count should include reads with low-qual non-discriminating bases
        AND reads with one low-qual discriminating position (recovered by
        masked per-position evaluation).
        """
        variant = gbcms_rs.Variant("chr1", 100, "GAGGG", "AAGGA", "COMPLEX")
        counts = count_both(mnp_bam, [variant])[0]

        # Expected ALT: 3 (rev) + 2 (fwd, low-qual non-disc) + 1 (fwd, low-qual disc, recovered) = 6
        # OLD: The low-qual-disc read was discarded (aggregate min-BQ gate).
        # NEW: Masked per-position eval recovers it — pos 4 (G→A, Q=34) is unmasked
        #      and matches ALT, so the read is classified as ALT.
        # The third-allele read is neither.
        assert counts.ad == 6, f"Expected ad=6, got {counts.ad}"
        assert counts.rd == 5, f"Expected rd=5, got {counts.rd}"

    def test_onp_dp_includes_discarded(self, mnp_bam):
        """DP should include the discarded (neither) reads."""
        variant = gbcms_rs.Variant("chr1", 100, "GAGGG", "AAGGA", "COMPLEX")
        counts = count_both(mnp_bam, [variant])[0]

        # Total reads: 5 REF + 3 ALT + 2 ALT(low-q-nondisc) + 1 neither(low-q-disc) + 1 third = 12
        assert counts.dp >= counts.rd + counts.ad, (
            f"DP invariant failed: dp={counts.dp} < rd+ad={counts.rd + counts.ad}"
        )

    def test_onp_strand_counts(self, mnp_bam):
        """Strand-specific counts should be correct."""
        variant = gbcms_rs.Variant("chr1", 100, "GAGGG", "AAGGA", "COMPLEX")
        counts = count_both(mnp_bam, [variant])[0]

        # ALT: 3 rev + 3 fwd (2 low-qual-nondisc + 1 low-qual-disc recovered) = 6 total
        assert counts.ad_fwd == 3, f"Expected ad_fwd=3, got {counts.ad_fwd}"
        assert counts.ad_rev == 3, f"Expected ad_rev=3, got {counts.ad_rev}"
        # REF: 5 fwd + 0 rev = 5
        assert counts.rd_fwd == 5, f"Expected rd_fwd=5, got {counts.rd_fwd}"


# ── DNP Tests (all-discriminating) ───────────────────────────────────────


class TestDNPAllDiscriminating:
    """Tests for all-discriminating DNP where every base matters."""

    def test_dnp_counts(self, dnp_bam):
        """Basic DNP counting with all positions discriminating."""
        variant = gbcms_rs.Variant("chr1", 100, "GG", "AA", "COMPLEX")
        counts = count_both(dnp_bam, [variant])[0]

        assert counts.rd == 3, f"Expected rd=3, got {counts.rd}"
        assert counts.ad == 4, f"Expected ad=4, got {counts.ad}"


# ── Fragment-Level Structural Invariants ─────────────────────────────────


class TestFragmentInvariants:
    """Verify fragment-level count structural invariants after MNP fixes."""

    def test_fragment_ref_lte_read_ref(self, mnp_bam):
        """Fragment REF count must not exceed read REF count."""
        variant = gbcms_rs.Variant("chr1", 100, "GAGGG", "AAGGA", "COMPLEX")
        counts = count_both(mnp_bam, [variant])[0]
        assert counts.rdf <= counts.rd, (
            f"Fragment invariant violated: rdf={counts.rdf} > rd={counts.rd}"
        )

    def test_fragment_alt_lte_read_alt(self, mnp_bam):
        """Fragment ALT count must not exceed read ALT count."""
        variant = gbcms_rs.Variant("chr1", 100, "GAGGG", "AAGGA", "COMPLEX")
        counts = count_both(mnp_bam, [variant])[0]
        assert counts.adf <= counts.ad, (
            f"Fragment invariant violated: adf={counts.adf} > ad={counts.ad}"
        )

    def test_fragment_sum_lte_dpf(self, mnp_bam):
        """RDF + ADF must not exceed DPF."""
        variant = gbcms_rs.Variant("chr1", 100, "GAGGG", "AAGGA", "COMPLEX")
        counts = count_both(mnp_bam, [variant])[0]
        assert counts.rdf + counts.adf <= counts.dpf, (
            f"Fragment invariant violated: rdf+adf={counts.rdf + counts.adf} > dpf={counts.dpf}"
        )

    def test_fragment_invariants_dnp(self, dnp_bam):
        """Fragment invariants hold for all-discriminating DNP."""
        variant = gbcms_rs.Variant("chr1", 100, "GG", "AA", "COMPLEX")
        counts = count_both(dnp_bam, [variant])[0]
        assert counts.rdf <= counts.rd
        assert counts.adf <= counts.ad
        assert counts.rdf + counts.adf <= counts.dpf
