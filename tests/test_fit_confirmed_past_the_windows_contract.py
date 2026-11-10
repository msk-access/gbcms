"""A fit resting on fewer bases than tell the alleles from the others beside them is
confirmed past the windows (#278 review).

One base past a pure indel's tract decides (RJ-3, D1), but another allele's extra base
can repeat that flank base, and a read shifted by another length can fit a window
that ends on the flank but for one mismatch RJ-26 tolerates as an error. Past the
windows its bases show the other length. So the ALT where the read's own bases place
no change or another one, a pure indel's REF where the read writes none, and either
allele's fit on a tolerated mismatch must fit the next 8 bases too, but for one error.
Each shape is one read of another allele beside ten REF reads and ten carriers.
"""

import random

import pytest
from helpers import count_checked, make_read, write_contig

from gbcms import _rs as gbcms_rs

READ = 120


def _flank(seed, n=200):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _count(tmp_path, contig, reads, pos, ref, alt):
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", pos, ref, alt, "X")], fa, 5, False, 1, True
    )
    (c,) = count_checked(bam, [pv.variant])
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    return c.rd, c.ad, c.partial_alt


def _ref_alt_reads(contig, p, ref, alt, n=10):
    reads = [
        make_read(f"r{i}", contig[p - 60 + i : p + 60 + i], p - 60 + i, ((0, READ),))
        for i in range(n)
    ]
    hap = contig[:p] + alt + contig[p + len(ref) :]
    d = len(alt) - len(ref)
    for i in range(n):
        s = p - 60 + i
        ll = p + 1 - s
        cig = ((0, ll), (1, d), (0, READ - ll - d)) if d > 0 else ((0, ll), (2, -d), (0, READ - ll))
        reads.append(make_read(f"a{i}", hap[s : s + READ], s, cig))
    return reads


def _unwritten(other, p, run_end, tail, clip):
    """A read of `other` ending `tail` bases past its changed run, its tail written as an
    aligner writes a short one: soft-clipped after the anchor, or as mismatches."""
    end = run_end + tail
    s = end - READ
    ll = p + 1 - s
    cig = ((0, ll), (4, READ - ll)) if clip else ((0, READ),)
    return make_read(f"o{tail}{'S' if clip else 'M'}", other[s:end], s, cig)


@pytest.mark.parametrize("tail", [3, 6])
@pytest.mark.parametrize("clip", [True, False])
def test_plus_g_read_at_a_deletion_joining_two_runs_is_not_ref(tmp_path, tail, clip):
    """CGG>C between C runs (CCC GG CCC -> C6). A read of the +G allele (CCC GGG CCC),
    its tail unwritten, reads past the run bases shifted by one: not REF, though the REF
    window ends on the one flank base, where the shift is one mismatch RJ-26 tolerates."""
    left = _flank(11) + "ATTCC"
    contig = left + "C" + "GG" + "CCCGATCGAT" + _flank(12)
    p = len(left)
    other = contig[: p + 1] + "G" + contig[p + 1 :]
    reads = _ref_alt_reads(contig, p, "CGG", "C") + [_unwritten(other, p, p + 4, tail, clip)]
    rd, ad, _ = _count(tmp_path, contig, reads, p, "CGG", "C")
    assert (rd, ad) == (10, 10)


@pytest.mark.parametrize("tail", [3, 6])
@pytest.mark.parametrize("clip", [True, False])
def test_minus_g_read_beside_an_insert_after_the_first_g_is_not_ref(tmp_path, tail, clip):
    """The locus of carrier.rs's growth unit test (an insert of other bases after the first
    G of GGG). A read of the -G allele (GG), its tail unwritten: not REF."""
    left = _flank(15) + "CTCATCTC"
    contig = left + "GGG" + "TCTCACTC" + _flank(16)
    p = len(left)
    other = contig[: p + 1] + contig[p + 2 :]
    reads = _ref_alt_reads(contig, p, "G", "GATCCACT") + [_unwritten(other, p, p + 2, tail, clip)]
    rd, ad, _ = _count(tmp_path, contig, reads, p, "G", "GATCCACT")
    assert (rd, ad) == (10, 10)


@pytest.mark.parametrize("tail", [3, 6])
def test_a6_read_at_an_a5_deletion_is_not_ref(tmp_path, tail):
    """GA>G in A5: an A6 read with a clipped tail reads its shifted bases past the run: not REF."""
    left = _flank(13) + "GTC"
    contig = left + "G" + "AAAAA" + "CTGTACG" + _flank(14)
    p = len(left)
    other = contig[: p + 1] + "A" + contig[p + 1 :]
    reads = _ref_alt_reads(contig, p, "GA", "G") + [_unwritten(other, p, p + 7, tail, True)]
    rd, ad, _ = _count(tmp_path, contig, reads, p, "GA", "G")
    assert (rd, ad) == (10, 10)


@pytest.mark.parametrize("past", [10, 60])
def test_another_insert_whose_extra_base_is_the_flank_is_not_alt(tmp_path, past):
    """C>CATT before a CAG run. A read inserting ATTC, written with its own 4-base insert
    (its bases place it inside the window: OwnEvent::Inside), spanning `past` bases beyond:
    another allele (RJ-27), never ALT: the ALT window ends on the flank C the extra base
    repeats, so the bases past it decide."""
    left = _flank(21) + "GACAGGAAATGC"
    contig = left + "CAGCAGCAGCAGCAGCCCTTAATGCG" + _flank(22)
    p = len(left) - 1
    other = contig[: p + 1] + "ATTC" + contig[p + 1 :]
    s = p + 1 + 4 + past - READ
    ll = p + 1 - s
    o = make_read("o", other[s : s + READ], s, ((0, ll), (1, 4), (0, READ - ll - 4)))
    rd, ad, partial = _count(
        tmp_path, contig, _ref_alt_reads(contig, p, "C", "CATT") + [o], p, "C", "CATT"
    )
    assert (rd, ad, partial) == (10, 10, 1)


def test_another_insert_whose_extra_base_is_the_flank_is_not_alt_plain(tmp_path):
    """A>AC, a read inserting CG (its G the flank base), written 2I and spanning 60 bases
    past: another allele, not ALT."""
    left = _flank(23) + "GAGACCGACAAA"
    contig = left + "GCGCCGCGATAGAGTAATC" + _flank(24)
    p = len(left) - 1
    other = contig[: p + 1] + "CG" + contig[p + 1 :]
    s = p + 1 + 2 + 60 - READ
    ll = p + 1 - s
    o = make_read("o", other[s : s + READ], s, ((0, ll), (1, 2), (0, READ - ll - 2)))
    rd, ad, partial = _count(
        tmp_path, contig, _ref_alt_reads(contig, p, "A", "AC") + [o], p, "A", "AC"
    )
    assert (rd, ad, partial) == (10, 10, 1)


def _del15_locus():
    rng = random.Random(7)
    while True:
        mid = "".join(rng.choice("ACGT") for _ in range(40))
        if mid[15] == mid[18] and mid[16] != mid[19]:
            break
    left = _flank(61) + "GTTAAAATTCCCGTCGCTATCAA"
    return left + mid + _flank(62), len(left) - 1


@pytest.mark.parametrize("style", ["written", "clip", "allM"])
@pytest.mark.parametrize("past", [20, 60])
def test_another_deletion_from_the_same_anchor_is_not_alt(tmp_path, style, past):
    """A 15-bp deletion whose next base repeats the anchor (an A run on the ALT, beside
    the change): the ALT window is the anchor and that one base. A read
    deleting 18 bp from the same anchor (the base after it the same A), spanning `past`
    bases beyond, is another allele, however it is written."""
    contig, p = _del15_locus()
    P = p + 1
    ref, alt = contig[p : P + 15], contig[p]
    other = contig[:P] + contig[P + 18 :]
    s = P - 60
    seq = other[s : s + 60 + past]
    cig = {
        "written": ((0, 60), (2, 18), (0, past)),
        "clip": ((0, 60), (4, past)),
        "allM": ((0, 60 + past),),
    }[style]
    reads = _ref_alt_reads(contig, p, ref, alt) + [make_read("o", seq, s, cig)]
    rd, ad, partial = _count(tmp_path, contig, reads, p, ref, alt)
    assert (rd, ad) == (10, 10)


@pytest.mark.parametrize("style", ["clip", "allM"])
def test_another_insert_whose_extra_base_is_the_flank_unwritten_is_not_alt(tmp_path, style):
    """C>CATT before a CAG run, a single +ATTC read with its insert unwritten (32 bases
    past): another allele, not ALT (fitting the ALT better than REF is not enough)."""
    left = _flank(21) + "GACAGGAAATGC"
    contig = left + "CAGCAGCAGCAGCAGCCCTTAATGCG" + _flank(22)
    p = len(left) - 1
    other = contig[: p + 1] + "ATTC" + contig[p + 1 :]
    s, L = p - 13, 50
    cig = ((0, 14), (4, L - 14)) if style == "clip" else ((0, L),)
    reads = _ref_alt_reads(contig, p, "C", "CATT") + [make_read("o", other[s : s + L], s, cig)]
    rd, ad, partial = _count(tmp_path, contig, reads, p, "C", "CATT")
    assert (rd, ad) == (10, 10)


def test_carriers_with_a_phased_neighbour_past_the_window_stay_alt(tmp_path):
    """TC>T in a CC run, every carrier with a DNP 8 bases on (a phased neighbour on the
    ALT haplotype, as at a gate row where 230 carriers carry one). Past the windows
    they mismatch the reference at those two bases, unshifted: no other length explains
    them better, so they are the ALT with a separate event beside it, as before."""
    left = _flank(31) + "CACCAGAACTCCAGCT"
    contig = left + "CCTACTACAGCCACGCGGCAGGCCAGGGCACCGG" + _flank(32)
    p = len(left) - 1  # the anchor T
    alt = contig[: p + 1] + contig[p + 2 :]
    dnp = list(alt)
    dnp[p + 8 : p + 10] = list("TA") if alt[p + 8 : p + 10] != "TA" else list("GC")
    dnp = "".join(dnp)
    reads = [
        make_read(f"r{i}", contig[p - 60 + i : p + 60 + i], p - 60 + i, ((0, READ),))
        for i in range(10)
    ]
    for i in range(10):
        s_ = p - 60 + i
        ll = p + 1 - s_
        reads.append(make_read(f"a{i}", dnp[s_ : s_ + READ], s_, ((0, ll), (2, 1), (0, READ - ll))))
    rd, ad, partial = _count(tmp_path, contig, reads, p, "TC", "T")
    assert (rd, ad, partial) == (10, 10, 0)
