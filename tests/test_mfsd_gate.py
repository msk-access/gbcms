"""PF-1: mFSD compute is gated on the ``mfsd`` flag in the binned production path.

When ``mfsd=False`` (the default for plain ``--dna``), the engine skips the
fragment-size statistics entirely — the per-variant size arrays and every
``mfsd_*`` field stay at their ``BaseCounts`` defaults — *without* perturbing the
actual allele counts. When ``mfsd=True`` the stats are computed exactly as before.

This keeps the engine output-aware (no compute-then-discard, AGENTS invariant #3)
and spares the held ``ref_sizes``/``alt_sizes`` arrays — the dominant mFSD memory
cost under Nextflow fan-out, where N gbcms processes run concurrently.
"""

import pysam
import pytest
from helpers import PARITY_FIELDS

from gbcms._rs import Variant, count_bam_binned


@pytest.fixture
def paired_bam(tmp_path):
    """Proper-pair BAM at chr1:100 (REF=A, ALT=T).

    Three REF pairs and two ALT pairs, each mate spanning pos 100 with
    ``template_length=±200`` so every fragment carries a valid cfDNA-range
    (50–1000 bp) insert size — the precondition for an mFSD size sample.
    """
    bam_path = tmp_path / "paired.bam"
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}

    def make_pair(qname, base):
        # Both mates span pos 100 (index 5 of a 10M read at start 95), overlapping
        # like a short cfDNA fragment; TLEN=200 → physical insert size 200.
        reads = []
        for is_r2 in (False, True):
            a = pysam.AlignedSegment()
            a.query_name = qname
            a.query_sequence = "AAAAA" + base + "AAAA"
            # paired + proper; read1 fwd / read2 rev
            a.flag = (1 | 2 | 128 | 16) if is_r2 else (1 | 2 | 64)
            a.reference_id = 0
            a.reference_start = 95
            a.mapping_quality = 60
            a.cigartuples = [(0, 10)]
            a.query_qualities = [30] * 10  # type: ignore[assignment]
            a.next_reference_id = 0
            a.next_reference_start = 95
            a.template_length = -200 if is_r2 else 200
            reads.append(a)
        return reads

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        for i in range(3):
            for r in make_pair(f"ref_{i}", "A"):
                outf.write(r)
        for i in range(2):
            for r in make_pair(f"alt_{i}", "T"):
                outf.write(r)

    sorted_bam = tmp_path / "paired.sorted.bam"
    pysam.sort("-o", str(sorted_bam), str(bam_path))
    pysam.index(str(sorted_bam))
    return str(sorted_bam)


def _count(bam, *, mfsd):
    variant = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")
    return count_bam_binned(
        bam,
        [variant],
        [None],
        min_mapq=20,
        min_baseq=20,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
        mfsd=mfsd,
    )[0]


def test_mfsd_off_skips_compute(paired_bam):
    """Default (mfsd=False): every mFSD field stays at its BaseCounts default.

    (The raw ``ref_sizes``/``alt_sizes`` arrays are internal — exposed only via the
    parquet writer — so the class-count/mean fields stand in for "nothing computed".)
    """
    off = _count(paired_bam, mfsd=False)
    assert off.mfsd_ref_count == 0
    assert off.mfsd_alt_count == 0
    assert off.mfsd_ref_mean == 0.0
    assert off.mfsd_alt_mean == 0.0


def test_mfsd_on_computes_sizes(paired_bam):
    """mfsd=True: the same proper-pair fragments now yield size samples + stats."""
    on = _count(paired_bam, mfsd=True)
    assert on.mfsd_ref_count == 3
    assert on.mfsd_alt_count == 2
    assert on.mfsd_ref_mean == 200.0
    assert on.mfsd_alt_mean == 200.0


def test_mfsd_gate_is_count_neutral(paired_bam):
    """Gating mFSD must not change any allele/depth count — only the mFSD fields."""
    off = _count(paired_bam, mfsd=False)
    on = _count(paired_bam, mfsd=True)
    for field in PARITY_FIELDS:
        assert getattr(off, field) == getattr(on, field), f"{field} drifted when mFSD toggled"
