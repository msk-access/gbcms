"""Large deletions count end to end, the same under any bin geometry.

A 50bp deletion's carrier and a read spanning the locus without it classify as ALT
and REF, and a bin anchored by the deletion, with an SNV inside its span, counts
both variants. `count_checked` also counts one variant per bin with the
per-variant fetch and asserts every field agrees, so each test checks that the
bin's single fetch covers the anchor's span. The count assertions prove the path
is exercised, so no test passes vacuously on all-zero counts.

Synthetic fixtures only — no patient data.
"""

from pathlib import Path

from helpers import build_bam, count_checked, make_read

from gbcms._rs import Variant

# A 50bp deletion: REF spans [100, 151) (51bp incl. the anchor at 100); ALT keeps
# only the anchor base, so 50 bases [101, 151) are deleted.
_DEL_POS = 100
_DEL_REF = "A" * 51
_DEL_ALT = "A"


def _large_deletion() -> Variant:
    return Variant(
        chrom="chr1",
        pos=_DEL_POS,
        ref_allele=_DEL_REF,
        alt_allele=_DEL_ALT,
        variant_type="DELETION",
    )


def test_a_large_deletion_counts_alt_and_ref(tmp_path: Path) -> None:
    """A read carrying the full 50bp deletion is ALT; a read spanning the locus
    with no deletion is REF."""
    # ALT: 31M 50D 30M — the 50D covers exactly [101, 151), the deleted span.
    alt_read = make_read("alt_frag", "A" * 61, 70, ((0, 31), (2, 50), (0, 30)), quals=[30] * 61)
    # REF: 111M spanning [70, 181) with no deletion at the locus.
    ref_read = make_read("ref_frag", "A" * 111, 70, ((0, 111),), quals=[30] * 111)

    bam_path = build_bam(tmp_path, [alt_read, ref_read], filename="large_del.bam")
    counts = count_checked(bam_path, [_large_deletion()])[0]

    # Both classes were genuinely exercised.
    assert counts.ad == 1, f"expected 1 ALT read, got {counts.ad}"
    assert counts.rd == 1, f"expected 1 REF read, got {counts.rd}"
    assert counts.dp >= counts.rd + counts.ad
    assert counts.ad == counts.ad_fwd + counts.ad_rev
    assert counts.rd == counts.rd_fwd + counts.rd_rev


def test_a_bin_anchored_by_a_large_deletion_counts_both_variants(tmp_path: Path) -> None:
    """A bin anchored by a large deletion with an SNV inside its REF span: both
    variants count, and the same with one variant per bin. The anchor's wide REF
    span drives the bin footprint the inner SNV also relies on."""
    big_del = _large_deletion()
    # SNP at 130 sits inside the deletion's ref span [100, 151).
    inner_snp = Variant(
        chrom="chr1",
        pos=130,
        ref_allele="A",
        alt_allele="T",
        variant_type="SNP",
    )

    # Fragment 1 carries the deletion (ALT for the del; the SNP locus is deleted,
    # so it contributes no SNP coverage).
    del_read = make_read("del_frag", "A" * 61, 70, ((0, 31), (2, 50), (0, 30)), quals=[30] * 61)
    # Fragment 2 spans the locus with no deletion (REF for the del) and carries a
    # 'T' at pos 130 (ALT for the SNP). seq index of pos 130 = 130 - 70 = 60.
    ref_seq = list("A" * 111)
    ref_seq[60] = "T"
    snp_read = make_read("snp_frag", "".join(ref_seq), 70, ((0, 111),), quals=[30] * 111)

    bam_path = build_bam(tmp_path, [del_read, snp_read], filename="del_anchor.bam")
    results = count_checked(bam_path, [big_del, inner_snp])

    del_counts, snp_counts = results
    # Deletion: one ALT (del_read), one REF (snp_read spans without a deletion).
    assert del_counts.ad == 1
    assert del_counts.rd == 1
    # SNP: only snp_read covers pos 130 (del_read deleted it) → one ALT, no REF.
    assert snp_counts.ad == 1, f"expected 1 SNP ALT, got {snp_counts.ad}"
    assert snp_counts.dp >= snp_counts.ad
