"""Group 2 read judgment (#201, #202): a read deleting a pure indel's anchor is
judged by its bases across the discrimination window, the aligner's placement of
its gap a tie-break (RJ-15); the ALT written across several ops counts ALT
(RJ-14). See docs/reference/read-judgment.md.
"""

import random

import pytest
from census import census
from helpers import make_read, write_contig

from gbcms import _rs

A = 400


def _contig(motif, seed=211):
    rng = random.Random(seed)
    seq = [rng.choice("CT") for _ in range(900)]
    seq[A - 3 : A - 3 + len(motif)] = list(motif)
    return "".join(seq)


def _count(bam, variant):
    (c,) = _rs.count_bam_binned(
        bam, [variant], [None], 20, 20, True, True, True, False, False, False, 1
    )
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    return c


def _prepared(fa, ref, alt):
    (pv,) = _rs.prepare_variants([_rs.Variant("1", A, ref, alt, "X")], fa, 5, False, 1, True)
    return pv.variant


def _deleting_anchor(contig, n_del, masked_at=()):
    """Four reads deleting n_del bases from the anchor A; bases at reference
    positions in masked_at are N at quality 2."""
    reads = []
    for i in range(4):
        s = A - 50 + i
        seq = contig[s:A] + contig[A + n_del : s + 100 + n_del]
        quals = [30] * len(seq)
        for p in masked_at:
            q = p - s - n_del
            seq = seq[:q] + "N" + seq[q + 1 :]
            quals[q] = 2
        reads.append(
            make_read(
                f"d{i}", seq, s, ((0, A - s), (2, n_del), (0, len(seq) - (A - s))), quals=quals
            )
        )
    return reads


def test_a_masked_tie_at_the_deleted_anchor_is_alt(tmp_path):
    """GA>G in CCA G AA GAA. The read deletes the anchor G and shows N where the
    ALT keeps it: its bases are the ALT's, with the aligner's tie-break placing
    the gap on the G. Placed one base along, the same bases take the strict path
    and count ALT; the call does not depend on the placement."""
    contig = _contig("CCAGAAGAA")
    reads = _deleting_anchor(contig, 1, masked_at=(A + 1,))
    fa, bam = write_contig(tmp_path, contig, reads, "tie")
    v = _prepared(fa, "GA", "G")
    c = _count(bam, v)
    assert (c.rd, c.ad) == (0, 4), (c.rd, c.ad, c.partial_alt)


def test_a_longer_deletion_over_the_anchor_is_not_the_alt(tmp_path):
    """GA>G; the reads delete four bases from the anchor. Their bases are not the
    1bp deletion's: another allele, partial evidence, never ALT or REF."""
    contig = _contig("CCAGAAGTCTCT")
    fa, bam = write_contig(tmp_path, contig, _deleting_anchor(contig, 4), "long")
    v = _prepared(fa, "GA", "G")
    c = _count(bam, v)
    assert (c.rd, c.ad, c.partial_alt) == (0, 0, 4), (c.rd, c.ad, c.partial_alt)
    assert census(bam, contig, v).ad == 0


@pytest.mark.xfail(
    strict=True, reason="C28 #202: Phase 3 credits a read with no base read at the event"
)
def test_a_deleted_anchor_with_every_window_base_masked_is_not_alt(tmp_path):
    """GA>G; the reads delete the anchor and every base left between their aligned
    flanks is N, the shared flank base too. Nothing at the event is read: not ALT."""
    contig = _contig("CCAGAAGTC")
    reads = _deleting_anchor(contig, 1, masked_at=(A + 1, A + 2, A + 3))
    fa, bam = write_contig(tmp_path, contig, reads, "masked")
    v = _prepared(fa, "GA", "G")
    c = _count(bam, v)
    assert c.ad == 0 and c.rd == 0, (c.rd, c.ad, c.partial_alt)
