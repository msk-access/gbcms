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
from helpers import build_bam, count_both, count_one, make_read

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
        reads.append(make_read(f"alt_rev_{i}", "AAAAAATT", start=97, cigar=((0, 8),), flag=16))

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
        assert (
            counts.dp >= counts.rd + counts.ad
        ), f"DP invariant failed: dp={counts.dp} < rd+ad={counts.rd + counts.ad}"

    def test_onp_strand_counts(self, mnp_bam):
        """Strand-specific counts should be correct."""
        variant = gbcms_rs.Variant("chr1", 100, "GAGGG", "AAGGA", "COMPLEX")
        counts = count_both(mnp_bam, [variant])[0]

        # ALT: 3 rev + 3 fwd (2 low-qual-nondisc + 1 low-qual-disc recovered) = 6 total
        assert counts.ad_fwd == 3, f"Expected ad_fwd=3, got {counts.ad_fwd}"
        assert counts.ad_rev == 3, f"Expected ad_rev=3, got {counts.ad_rev}"
        # REF: 5 fwd + 0 rev = 5
        assert counts.rd_fwd == 5, f"Expected rd_fwd=5, got {counts.rd_fwd}"


class TestONPCarrierShapes:
    """Which carrier reads a sparse ONP counts as full ALT vs partial.

    GAGGG>AAGGA discriminates only at positions 0 and 4; the interior GG (and
    the A at position 1) match both alleles. A read carrying the whole
    haplotype is full ALT whatever the interior quality. A read carrying only
    ONE of the two G>A changes is partial — it does not carry the annotated
    allele. A sign-out row reporting many more ALT reads than gbcms AD, with
    the gap sitting in partial_alt, therefore means the carriers hold the two
    changes on different molecules (trans / subclonal / merged-SNV annotation),
    not that the masked comparison drops cis carriers.
    """

    LEFT = "TTACCGTACGCATTCACGTA"  # 80..99
    RIGHT = "CTACGTTGCAACGTTACGAT"  # 105..124

    def _bam(self, tmp_path, carrier_block, carrier_quals=None):
        reads = [
            make_read(
                f"ref_{i}",
                self.LEFT + "GAGGG" + self.RIGHT,
                start=80,
                cigar=((0, 45),),
                flag=16 if i % 2 else 0,
            )
            for i in range(20)
        ]
        reads += [
            make_read(
                f"carrier_{i}",
                self.LEFT + carrier_block + self.RIGHT,
                start=80,
                cigar=((0, 45),),
                flag=16 if i % 2 else 0,
                quals=carrier_quals,
            )
            for i in range(10)
        ]
        return build_bam(tmp_path, reads, "onp_shapes.bam")

    @pytest.mark.parametrize(
        "block, low_bq_offset, expected_ad, expected_partial",
        [
            ("AAGGA", None, 10, 0),  # cis carrier
            ("AAGGA", 2, 10, 0),  # low-BQ interior base cannot vote
            ("AAGGA", 4, 10, 0),  # masked discriminating base, other still votes
            ("AAGGG", None, 0, 10),  # only the first G>A
            ("GAGGA", None, 0, 10),  # only the second G>A
        ],
    )
    def test_carrier_shape(self, tmp_path, block, low_bq_offset, expected_ad, expected_partial):
        quals = None
        if low_bq_offset is not None:
            quals = [30] * 45
            quals[len(self.LEFT) + low_bq_offset] = 5
        bam = self._bam(tmp_path, block, quals)
        variant = gbcms_rs.Variant("chr1", 100, "GAGGG", "AAGGA", "ONP")
        counts = count_both(bam, [variant])[0]

        assert counts.rd == 20
        assert counts.ad == expected_ad
        assert counts.partial_alt == expected_partial
        assert counts.dp >= counts.rd + counts.ad
        assert counts.dpf >= counts.rdf + counts.adf
        assert counts.rd == counts.rd_fwd + counts.rd_rev
        assert counts.ad == counts.ad_fwd + counts.ad_rev

    @pytest.mark.parametrize(
        "block, low_bq_offset, expected_confirmed",
        [
            ("AAGGA", None, 10),  # every discriminating base read as ALT
            ("AAGGA", 2, 10),  # interior base is not discriminating
            ("AAGGA", 4, 0),  # ALT on one read base only: not confirmed
            ("AAGGG", None, 0),  # partial, never ALT
        ],
    )
    def test_confirmed_alt_counts_fully_read_carriers_only(
        self, tmp_path, block, low_bq_offset, expected_confirmed
    ):
        quals = None
        if low_bq_offset is not None:
            quals = [30] * 45
            quals[len(self.LEFT) + low_bq_offset] = 5
        bam = self._bam(tmp_path, block, quals)
        variant = gbcms_rs.Variant("chr1", 100, "GAGGG", "AAGGA", "ONP")
        counts = count_both(bam, [variant])[0]
        legacy = count_one(bam, variant)

        assert counts.mnp_confirmed_alt == expected_confirmed
        assert legacy.mnp_confirmed_alt == expected_confirmed
        assert counts.mnp_confirmed_alt <= counts.ad
        assert counts.dp >= counts.rd + counts.ad
        assert counts.dpf >= counts.rdf + counts.adf
        assert counts.rd == counts.rd_fwd + counts.rd_rev
        assert counts.ad == counts.ad_fwd + counts.ad_rev

    def test_indel_inside_block_is_never_confirmed(self, tmp_path):
        """A read carrying both changed bases plus a 1bp insertion inside the
        block is a different allele, not the annotated haplotype: it may count
        toward ad through the complex path, but never as a read that shows the
        whole MNP. Counting it toward ad also keeps the rescue gate's
        partial_alt > ad from opening at such a locus."""
        reads = [
            make_read(
                f"ref_{i}",
                self.LEFT + "GAGGG" + self.RIGHT,
                start=80,
                cigar=((0, 45),),
                flag=16 if i % 2 else 0,
            )
            for i in range(20)
        ]
        reads += [
            make_read(
                f"ins_{i}",
                self.LEFT + "AAG" + "T" + "GA" + self.RIGHT,
                start=80,
                cigar=((0, 23), (1, 1), (0, 22)),
                flag=16 if i % 2 else 0,
            )
            for i in range(10)
        ]
        bam = build_bam(tmp_path, reads, "onp_ins_in_block.bam")
        variant = gbcms_rs.Variant("chr1", 100, "GAGGG", "AAGGA", "ONP")
        counts = count_both(bam, [variant])[0]
        legacy = count_one(bam, variant)

        assert counts.mnp_confirmed_alt == 0
        assert legacy.mnp_confirmed_alt == 0
        assert counts.ad == 10
        assert counts.partial_alt == 0
        assert counts.dp >= counts.rd + counts.ad
        assert counts.dpf >= counts.rdf + counts.adf
        assert counts.rd == counts.rd_fwd + counts.rd_rev
        assert counts.ad == counts.ad_fwd + counts.ad_rev

    def test_confirmed_alt_is_zero_for_non_mnp_variants(self, tmp_path):
        bam = self._bam(tmp_path, "AAGGA")
        snv = gbcms_rs.Variant("chr1", 100, "G", "A", "SNP")
        counts = count_both(bam, [snv])[0]
        assert counts.ad == 10
        assert counts.mnp_confirmed_alt == 0
        assert counts.dp >= counts.rd + counts.ad
        assert counts.dpf >= counts.rdf + counts.adf
        assert counts.rd == counts.rd_fwd + counts.rd_rev
        assert counts.ad == counts.ad_fwd + counts.ad_rev


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
        assert (
            counts.rdf <= counts.rd
        ), f"Fragment invariant violated: rdf={counts.rdf} > rd={counts.rd}"

    def test_fragment_alt_lte_read_alt(self, mnp_bam):
        """Fragment ALT count must not exceed read ALT count."""
        variant = gbcms_rs.Variant("chr1", 100, "GAGGG", "AAGGA", "COMPLEX")
        counts = count_both(mnp_bam, [variant])[0]
        assert (
            counts.adf <= counts.ad
        ), f"Fragment invariant violated: adf={counts.adf} > ad={counts.ad}"

    def test_fragment_sum_lte_dpf(self, mnp_bam):
        """RDF + ADF must not exceed DPF."""
        variant = gbcms_rs.Variant("chr1", 100, "GAGGG", "AAGGA", "COMPLEX")
        counts = count_both(mnp_bam, [variant])[0]
        assert (
            counts.rdf + counts.adf <= counts.dpf
        ), f"Fragment invariant violated: rdf+adf={counts.rdf + counts.adf} > dpf={counts.dpf}"

    def test_fragment_invariants_dnp(self, dnp_bam):
        """Fragment invariants hold for all-discriminating DNP."""
        variant = gbcms_rs.Variant("chr1", 100, "GG", "AA", "COMPLEX")
        counts = count_both(dnp_bam, [variant])[0]
        assert counts.rdf <= counts.rd
        assert counts.adf <= counts.ad
        assert counts.rdf + counts.adf <= counts.dpf
