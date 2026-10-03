"""Group 2 read judgment (#201, #202, #174): a read deleting a pure indel's anchor
is judged by its bases across the discrimination window, the aligner's placement
of its gap a tie-break (RJ-15); the ALT written across several ops counts ALT
(RJ-14); an exact-carrier ALT call needs quality-weighted evidence (RJ-16). See
docs/reference/read-judgment.md.
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


def test_a_deleted_anchor_with_every_window_base_masked_is_not_alt(tmp_path):
    """GA>G; the reads delete the anchor and every base left between their aligned
    flanks is N, the shared flank base too. Nothing at the event is read: not ALT."""
    contig = _contig("CCAGAAGTC")
    reads = _deleting_anchor(contig, 1, masked_at=(A + 1, A + 2, A + 3))
    fa, bam = write_contig(tmp_path, contig, reads, "masked")
    v = _prepared(fa, "GA", "G")
    c = _count(bam, v)
    assert c.ad == 0 and c.rd == 0, (c.rd, c.ad, c.partial_alt)


# ── C16: an exact-carrier ALT call needs quality-weighted evidence ───────────


def _carrier_reads(contig, seq_at, quals_at):
    """Four reads over the event at A, the reference except `seq_at` (position ->
    base), qualities 30 except `quals_at` (position -> quality)."""
    reads = []
    for i in range(4):
        s = A - 50 + i
        seq = list(contig[s : s + 100])
        quals = [30] * 100
        for p, b in seq_at.items():
            seq[p - s] = b
        for p, q in quals_at.items():
            quals[p - s] = q
        reads.append(make_read(f"c{i}", "".join(seq), s, ((0, 100),), quals=quals))
    return reads


def test_a_ref_read_with_one_error_and_a_low_quality_tail_is_not_alt(tmp_path):
    """TC>GCT in unique sequence. The reads are REF molecules read poorly around
    the event (Q10) except one clear error at its first base, which happens to be
    the ALT's G. Masked bases fit anything, so today that one base decides ALT;
    weighed by quality, the low-quality bases show the reference: not ALT."""
    contig = _contig("GACTCAGTTCGA")
    tail = {p: 10 for p in range(A - 8, A + 12) if p != A}
    fa, bam = write_contig(tmp_path, contig, _carrier_reads(contig, {A: "G"}, tail), "err")
    v = _prepared(fa, "TC", "GCT")
    c = _count(bam, v)
    assert c.ad == 0, (c.rd, c.ad, c.partial_alt)


def test_a_carrier_with_one_low_quality_event_base_is_alt(tmp_path):
    """Guard: ALT carriers of TC>GCT with one event base at Q10 keep their ALT
    call: their clearly read bases carry the evidence."""
    contig = _contig("GACTCAGTTCGA")
    reads = []
    for i in range(4):
        s = A - 50 + i
        seq = contig[s:A] + "GCT" + contig[A + 2 : s + 99]
        quals = [30] * len(seq)
        quals[A + 1 - s] = 10
        reads.append(
            make_read(
                f"a{i}", seq, s, ((0, A - s), (1, 1), (0, len(seq) - (A - s) - 1)), quals=quals
            )
        )
    fa, bam = write_contig(tmp_path, contig, reads, "car")
    v = _prepared(fa, "TC", "GCT")
    c = _count(bam, v)
    assert c.ad == 4, (c.rd, c.ad, c.partial_alt)


# ── Found by the adversarial review of group 2 ───────────────────────────────


@pytest.mark.xfail(
    strict=True, reason="RJ-15: a read whose bases are REF, gap at the anchor, counts partial"
)
@pytest.mark.parametrize("order", ["DI", "ID"])
@pytest.mark.parametrize("ref,alt", [("GT", "G"), ("G", "GA")])
def test_a_ref_read_written_as_a_gap_at_the_anchor_is_ref(tmp_path, order, ref, alt):
    """The reads are exactly REF, but the aligner wrote the anchor base deleted and
    re-inserted (D1 I1, or I1 D1). Judged by its bases, a read holding the REF
    allele is REF."""
    contig = _contig("CCAGTCTCTCC")
    reads = []
    for i in range(4):
        s = A - 50 + i
        k = A - s
        mid = ((2, 1), (1, 1)) if order == "DI" else ((1, 1), (2, 1))
        reads.append(
            make_read(f"r{i}", contig[s : s + 100], s, ((0, k),) + mid + ((0, 100 - k - 1),))
        )
    fa, bam = write_contig(tmp_path, contig, reads, "refgap")
    v = _prepared(fa, ref, alt)
    c = _count(bam, v)
    assert (c.rd, c.ad, c.partial_alt) == (4, 0, 0), (c.rd, c.ad, c.partial_alt)
    assert census(bam, contig, v).rd == 4


def _unique_contig(ctx, offset):
    rng = random.Random(7)
    seq = [rng.choice("CT") for _ in range(900)]
    seq[A - offset : A - offset + len(ctx)] = list(ctx)
    return "".join(seq)


@pytest.mark.xfail(
    strict=True, reason="RJ-16: REF and ALT readings over different spans skew the evidence"
)
def test_one_distinguishing_base_at_min_baseq_is_enough(tmp_path):
    """CTC>G carriers whose one distinguishing base is read at exactly --min-baseq
    (Q20), the event's other bases N: one base at the minimum quality is the
    evidence an ALT call needs (RJ-16)."""
    contig = _unique_contig("GATCCTGACTTCGCATGTCCAGTGACTCATGCGTTACAGGCTTAGCCATGCTTGACGTA", 25)
    hap = contig[:A] + "G" + contig[A + 3 :]
    reads = []
    for i in range(4):
        s, k, keep = A - 20, 20, 22
        seq = list(hap[s : s + 100])
        quals = [37] * 100
        for p in range(k - 3, k + 8):
            if p != keep:
                seq[p] = "N"
        quals[keep] = 20
        reads.append(
            make_read(f"g{i}", "".join(seq), s, ((0, k + 1), (2, 2), (0, 99 - k)), quals=quals)
        )
    fa, bam = write_contig(tmp_path, contig, reads, "oneq20")
    v = _prepared(fa, "CTC", "G")
    c = _count(bam, v)
    assert c.ad == 4, (c.rd, c.ad, c.partial_alt)


@pytest.mark.xfail(
    strict=True, reason="RJ-16: the evidence skips a long event's bases past its junction windows"
)
def test_a_long_event_carrier_with_low_quality_junction_bases_is_alt(tmp_path):
    """A 61-base delins to TG. The carriers' bases at the junction are low quality,
    but their next 70 bases clearly read the right flank, not the deleted body:
    the read reaches on through the event (RJ-9), and those bases are evidence."""
    left = "GATCCTGACTTCGCATGTCCAGTGAC"
    body = "".join("ACGGTTCAGT"[(i * 7 + i // 3) % 10] for i in range(61))
    right = "CATGCGTTACAGGCTTAGCCATGCTTGACGTACCTAGGATTC"
    contig = _unique_contig(left + body + right, len(left))
    hap = contig[:A] + "TG" + contig[A + 61 :]
    reads = []
    for i in range(4):
        s = A - 24
        k = A - s
        quals = [37] * 100
        quals[k - 4 : k + 6] = [12, 0, 23, 0, 0, 0, 0, 2, 37, 0]
        reads.append(
            make_read(
                f"l{i}", hap[s : s + 100], s, ((0, k + 2), (2, 59), (0, 100 - k - 2)), quals=quals
            )
        )
    fa, bam = write_contig(tmp_path, contig, reads, "long")
    v = _prepared(fa, body, "TG")
    c = _count(bam, v)
    assert c.ad == 4, (c.rd, c.ad, c.partial_alt)
