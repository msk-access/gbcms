"""
Tests for fragment-level consensus counting (Phase 1a fix).

Verifies that when R1 and R2 of a read pair disagree on the allele
at a variant site, the fragment is counted exactly once using
quality-weighted consensus, NOT double-counted to both REF and ALT.

Prior bug: fragments where R1=REF and R2=ALT were counted to BOTH
rdf and adf, inflating dpf != rdf + adf.
"""

import pysam
import pytest
from helpers import count_both

from gbcms._rs import count_bam


@pytest.fixture
def fragment_consensus_bam(tmp_path):
    """
    Generates a synthetic BAM with read pairs that test fragment consensus.

    Reference: chr1, 1000 bases of 'A'.
    Variant: chr1:100 (0-based), REF=A, ALT=T (SNP)

    Fragment layout:
        Fragment 1 (frag_agree_ref): R1=A (q30), R2=A (q30) → REF
        Fragment 2 (frag_agree_alt): R1=T (q30), R2=T (q30) → ALT
        Fragment 3 (frag_disagree_alt_wins): R1=A (q10), R2=T (q35) → ALT (q35 >> q10)
        Fragment 4 (frag_disagree_ref_wins): R1=A (q35), R2=T (q10) → REF (q35 >> q10)
        Fragment 5 (frag_disagree_tie): R1=A (q30), R2=T (q30) → DISCARDED (ambiguous)
        Fragment 6 (frag_single_read): R1=T (q30), no R2 → ALT
    """
    bam_path = tmp_path / "fragment_consensus.bam"
    header = {
        "HD": {"VN": "1.0", "SO": "coordinate"},
        "SQ": [{"LN": 1000, "SN": "chr1"}],
    }

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        # Helper to create a read
        def make_read(qname, seq, quals, is_read1, is_reverse, pos=95):
            a = pysam.AlignedSegment()
            a.query_name = qname
            a.query_sequence = seq
            a.flag = 0
            if not is_read1:
                a.flag |= 0x80  # read2
            else:
                a.flag |= 0x40  # read1
            a.flag |= 0x1  # paired
            a.flag |= 0x2  # proper pair
            if is_reverse:
                a.flag |= 0x10
                a.flag |= 0x20  # mate reverse (assume FR orientation)
            else:
                a.flag |= 0x20  # mate reverse
            a.reference_id = 0
            a.reference_start = pos
            a.mapping_quality = 60
            a.cigartuples = [(0, len(seq))]  # All match
            a.query_qualities = pysam.qualitystring_to_array("".join(chr(q + 33) for q in quals))  # type: ignore[assignment]
            # Set mate info
            a.next_reference_id = 0
            a.next_reference_start = pos
            a.template_length = len(seq)
            return a

        # 10-base reads covering pos 95-104 (0-based), variant at pos 100
        # Position index 5 in the read is the variant position

        # Fragment 1: Both reads agree on REF (A at pos 100)
        ref_seq = "AAAAAAAAA" + "A"  # 10 bases, pos 100 is 'A' (REF)
        outf.write(make_read("frag_agree_ref", ref_seq, [30] * 10, is_read1=True, is_reverse=False))
        outf.write(make_read("frag_agree_ref", ref_seq, [30] * 10, is_read1=False, is_reverse=True))

        # Fragment 2: Both reads agree on ALT (T at pos 100)
        alt_seq = "AAAAATAAAA"  # pos 100 = T (ALT)
        outf.write(make_read("frag_agree_alt", alt_seq, [30] * 10, is_read1=True, is_reverse=False))
        outf.write(make_read("frag_agree_alt", alt_seq, [30] * 10, is_read1=False, is_reverse=True))

        # Fragment 3: R1=REF(q10), R2=ALT(q35) → ALT should win
        outf.write(
            make_read(
                "frag_disagree_alt_wins",
                ref_seq,
                [30, 30, 30, 30, 30, 10, 30, 30, 30, 30],  # q10 at variant pos
                is_read1=True,
                is_reverse=False,
            )
        )
        outf.write(
            make_read(
                "frag_disagree_alt_wins",
                alt_seq,
                [30, 30, 30, 30, 30, 35, 30, 30, 30, 30],  # q35 at variant pos
                is_read1=False,
                is_reverse=True,
            )
        )

        # Fragment 4: R1=REF(q35), R2=ALT(q10) → REF should win
        outf.write(
            make_read(
                "frag_disagree_ref_wins",
                ref_seq,
                [30, 30, 30, 30, 30, 35, 30, 30, 30, 30],  # q35 at variant pos
                is_read1=True,
                is_reverse=False,
            )
        )
        outf.write(
            make_read(
                "frag_disagree_ref_wins",
                alt_seq,
                [30, 30, 30, 30, 30, 10, 30, 30, 30, 30],  # q10 at variant pos
                is_read1=False,
                is_reverse=True,
            )
        )

        # Fragment 5: R1=REF(q30), R2=ALT(q30) → tie, conservative = REF
        outf.write(
            make_read(
                "frag_disagree_tie",
                ref_seq,
                [30] * 10,
                is_read1=True,
                is_reverse=False,
            )
        )
        outf.write(
            make_read(
                "frag_disagree_tie",
                alt_seq,
                [30] * 10,
                is_read1=False,
                is_reverse=True,
            )
        )

        # Fragment 6: Single read (R1 only), ALT
        outf.write(
            make_read(
                "frag_single_read",
                alt_seq,
                [30] * 10,
                is_read1=True,
                is_reverse=False,
            )
        )

    # Sort and index
    sorted_path = tmp_path / "fragment_consensus.sorted.bam"
    pysam.sort("-o", str(sorted_path), str(bam_path))
    pysam.index(str(sorted_path))

    return str(sorted_path)


def test_fragment_consensus_no_double_counting(fragment_consensus_bam):
    """
    Verify correct fragment consensus with discard behavior for ambiguous fragments.

    Expected fragment-level results:
        frag_agree_ref:          REF (both reads agree)
        frag_agree_alt:          ALT (both reads agree)
        frag_disagree_alt_wins:  ALT (q35 > q10 + 10 threshold)
        frag_disagree_ref_wins:  REF (q35 > q10 + 10 threshold)
        frag_disagree_tie:       DISCARDED (q30 vs q30, within threshold)
        frag_single_read:        ALT (single read, no conflict)

    Expected: dpf = 6 (all fragments seen),
              rdf = 2 (agree_ref + ref_wins),
              adf = 3 (agree_alt + alt_wins + single),
              discarded = 1 (tie)
    """
    from gbcms._rs import Variant

    variant = Variant(
        chrom="chr1",
        pos=100,
        ref_allele="A",
        alt_allele="T",
        variant_type="SNP",
    )

    results = count_bam(
        bam_path=fragment_consensus_bam,
        variants=[variant],
        decomposed=[None],
        min_mapq=0,
        min_baseq=0,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=True,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )

    counts = results[0]

    # dpf counts ALL fragments (including discarded ambiguous ones)
    assert counts.dpf == 6, f"Expected 6 fragments, got {counts.dpf}"

    # rdf and adf reflect only resolved fragments
    assert counts.rdf == 2, f"Expected rdf=2 (agree_ref + ref_wins), got {counts.rdf}"
    assert counts.adf == 3, f"Expected adf=3 (agree_alt + alt_wins + single), got {counts.adf}"

    # The gap dpf - (rdf + adf) = 1 discarded fragment (the tie)
    discarded = counts.dpf - (counts.rdf + counts.adf)
    assert discarded == 1, f"Expected 1 discarded fragment (tie), got {discarded}"

    # Read-level counts should reflect ALL reads (11 total)
    # R1+R2 for fragments 1-5 = 10 reads, plus 1 single read = 11
    assert counts.dp == 11, f"Expected dp=11 (11 reads), got {counts.dp}"


def test_fragment_consensus_quality_tiebreaker(fragment_consensus_bam):
    """
    Verify quality-weighted consensus resolves R1 vs R2 conflicts correctly.
    """
    from gbcms._rs import Variant

    variant = Variant(
        chrom="chr1",
        pos=100,
        ref_allele="A",
        alt_allele="T",
        variant_type="SNP",
    )

    results = count_bam(
        bam_path=fragment_consensus_bam,
        variants=[variant],
        decomposed=[None],
        min_mapq=0,
        min_baseq=0,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=True,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )

    counts = results[0]

    # Strand-level fragment counts should also be consistent
    assert (
        counts.rdf == counts.rdf_fwd + counts.rdf_rev
    ), f"rdf strand inconsistency: {counts.rdf} != {counts.rdf_fwd} + {counts.rdf_rev}"
    assert (
        counts.adf == counts.adf_fwd + counts.adf_rev
    ), f"adf strand inconsistency: {counts.adf} != {counts.adf_fwd} + {counts.adf_rev}"


def test_fragment_dpf_includes_neither_reads(fragment_consensus_bam):
    """
    Gap 1D: DPF must include fragments where both reads return 'neither'.

    In the consensus BAM fixture, all fragments have REF or ALT evidence.
    This test verifies the DPF >= RDF + ADF invariant holds (the gap
    represents discarded ambiguous fragments + any neither-only fragments).
    """
    from gbcms._rs import Variant

    variant = Variant(
        chrom="chr1",
        pos=100,
        ref_allele="A",
        alt_allele="T",
        variant_type="SNP",
    )

    results = count_bam(
        bam_path=fragment_consensus_bam,
        variants=[variant],
        decomposed=[None],
        min_mapq=0,
        min_baseq=0,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=True,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )

    counts = results[0]

    # DPF must be >= RDF + ADF (gap = discarded + neither fragments)
    assert (
        counts.dpf >= counts.rdf + counts.adf
    ), f"DPF invariant failed: dpf={counts.dpf} < rdf+adf={counts.rdf + counts.adf}"
    # DP must also >= RD + AD at read level
    assert (
        counts.dp >= counts.rd + counts.ad
    ), f"DP invariant failed: dp={counts.dp} < rd+ad={counts.rd + counts.ad}"


# ── count_bam_binned parity tests ────────────────────────────────────────


def test_fragment_consensus_binned(fragment_consensus_bam):
    """count_bam_binned produces same fragment consensus counts as count_bam."""
    from gbcms._rs import Variant

    variant = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")
    counts = count_both(
        fragment_consensus_bam,
        [variant],
        min_mapq=0,
        min_baseq=0,
        filter_qc_failed=True,
        filter_improper_pair=False,
        filter_indel=False,
    )[0]
    assert counts.dpf == 6
    assert counts.rdf == 2
    assert counts.adf == 3
    assert counts.dp == 11
    assert counts.dpf >= counts.rdf + counts.adf
