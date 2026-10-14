"""Reads judged beside the other alleles at their locus (RJ-25), and a change inside the
repeat an event is grown through (RJ-23, 2026-10-14).

A read that fits the ALT counts ALT only when, for every other allele of the ALT's
length at the locus, it reads a base telling the two apart: a co-annotated sibling, or
an allele read in full on 3 or more fragments. Where the reads show no other allele, a
read reading only part of the ALT's bases counts as before. A change the REF molecules
carry inside the grown event, off the given allele's own bases, is background: masked,
so REF reads and carriers carrying it keep their call.

Each locus is checked against the read census told the true allele set.
"""

import random

from census import assert_matches, census
from helpers import count_checked, make_read, write_contig

from gbcms import _rs as gbcms_rs

READ = 100
P = 300  # 0-based start of every event


def _contig(n=800, seed=277, alphabet="CGT"):
    rng = random.Random(seed)
    return "".join(rng.choice(alphabet) for _ in range(n))


def _prepared(fa, rows):
    pvs = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", p, r, a, "COMPLEX") for p, r, a in rows], fa, 5, False, 1, True
    )
    assert all(pv.gbcms_status == "PASS" for pv in pvs), [pv.gbcms_status_reason for pv in pvs]
    return [pv.variant for pv in pvs]


def _low(read, at, q=5):
    """The read's bases at query positions `at` below min BQ, letters kept."""
    quals = list(read.query_qualities)
    for i in at:
        quals[i] = q
    read.query_qualities = quals
    return read


def _invariants(c):
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev


# ── another allele of the ALT's length (RJ-25) ───────────────────────────


REF_LEN, ALT, OTHER = 28, "ATATGA", "ATCTGA"  # the other allele: the third inserted base
TELL = 2


def _delins_reads(contig, alt, tag, n, first=0, low=False):
    """`n` carriers of the 28-to-6 delins spelled `alt`, written I then D, each its own
    fragment; with `low`, the telling base below min BQ."""
    hap = contig[:P] + alt + contig[P + REF_LEN :]
    out = []
    for i in range(first, first + n):
        s = P - 60 + i
        left = P - s
        r = make_read(
            f"{tag}{i}",
            hap[s : s + READ],
            s,
            ((0, left), (1, len(alt)), (2, REF_LEN), (0, READ - left - len(alt))),
        )
        out.append(_low(r, [left + TELL]) if low else r)
    return out


def _delins_locus(tmp_path, others_on, co_annotated=False):
    contig = _contig()
    ref = contig[P : P + REF_LEN]
    reads = [
        make_read(f"ref{i}", contig[P - 60 + i : P + 40 + i], P - 60 + i, ((0, READ),))
        for i in range(10)
    ]
    reads += _delins_reads(contig, ALT, "alt", 10)
    reads += _delins_reads(contig, ALT, "low", 6, first=10, low=True)
    reads += _delins_reads(contig, OTHER, "oth", others_on, first=20)
    fa, bam = write_contig(tmp_path, contig, reads)
    rows = [(P, ref, ALT)] + ([(P, ref, OTHER)] if co_annotated else [])
    vs = _prepared(fa, rows)
    siblings = [[w for w in vs if w is not v] for v in vs]
    counts = count_checked(bam, vs, sibling_variants=siblings)
    return contig, bam, vs, counts


def test_a_carrier_that_cannot_tell_the_alt_from_another_allele_is_not_alt(tmp_path):
    """Six fragments spell another allele one inserted base apart: the six carriers
    reading that base below min BQ fit both, so they are not the ALT's."""
    contig, bam, vs, counts = _delins_locus(tmp_path, others_on=6)
    c = counts[0]
    _invariants(c)
    assert (c.rd, c.ad) == (10, 10)
    assert_matches(c, census(bam, contig, vs[0], others=[(P, contig[P : P + REF_LEN], OTHER)]))


def test_a_co_annotated_allele_is_another_allele_at_any_depth(tmp_path):
    """The other allele co-annotated, with no read spelling it: the same."""
    contig, bam, vs, counts = _delins_locus(tmp_path, others_on=0, co_annotated=True)
    c = counts[0]
    _invariants(c)
    assert (c.rd, c.ad) == (10, 10)
    assert_matches(c, census(bam, contig, vs[0], others=[(P, contig[P : P + REF_LEN], OTHER)]))


def test_an_allele_on_two_fragments_is_not_looked_for(tmp_path):
    """Two fragments of another allele are under the 3 the rule needs: the carriers
    reading its telling base below min BQ count ALT, as with no other allele."""
    contig, bam, vs, counts = _delins_locus(tmp_path, others_on=2)
    c = counts[0]
    _invariants(c)
    assert (c.rd, c.ad) == (10, 16)
    assert_matches(c, census(bam, contig, vs[0]))


# ── an MNP joining a somatic change to a germline SNP beside it ──────────


def _mnp_contig():
    contig = list(_contig(seed=31))
    contig[P : P + 2] = "GC"
    return "".join(contig)


def _mnp_reads(contig, block, tag, n, first=0, low_at=None):
    """`n` plain reads carrying `block` at the MNP; `low_at` (0 or 1): that block base
    below min BQ."""
    out = []
    for i in range(first, first + n):
        s = P - 50 + (i % 40)
        seq = contig[s:P] + block + contig[P + 2 : s + READ]
        r = make_read(f"{tag}{i}", seq, s, ((0, READ),))
        out.append(_low(r, [P - s + low_at]) if low_at is not None else r)
    return out


def _mnp_locus(tmp_path, somatic, germline):
    contig = _mnp_contig()
    reads = _mnp_reads(contig, "GC", "ref", 20)
    if germline:
        reads += _mnp_reads(contig, "GT", "germ", 20)
        reads += _mnp_reads(contig, "GT", "germlow", 4, first=20, low_at=0)  # shows ?T
    if somatic:
        reads += _mnp_reads(contig, "AT", "som", 6, first=40)
        reads += _mnp_reads(contig, "AT", "somlowA", 4, first=46, low_at=0)  # shows ?T
        reads += _mnp_reads(contig, "AT", "somlowT", 3, first=50, low_at=1)  # shows A?
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, "GC", "AT")])
    (c,) = count_checked(bam, [v])
    return contig, bam, v, c


def test_a_read_showing_only_the_germline_base_of_a_merged_mnp_is_not_alt(tmp_path):
    """GC>AT written for a somatic G>A beside a germline C>T: a read whose A is below
    min BQ shows only the T, which the germline allele GT (on 24 fragments) shows too.
    Reads showing the A count, the T unread or not."""
    contig, bam, v, c = _mnp_locus(tmp_path, somatic=True, germline=True)
    _invariants(c)
    assert c.ad == 9
    assert_matches(c, census(bam, contig, v, others=[(P, "GC", "GT")]))


def test_beside_a_germline_allele_the_absent_given_allele_counts_no_alt(tmp_path):
    """No AT molecule: the germline reads with their G below min BQ are no ALT."""
    contig, bam, v, c = _mnp_locus(tmp_path, somatic=False, germline=True)
    _invariants(c)
    assert c.ad == 0
    assert_matches(c, census(bam, contig, v, others=[(P, "GC", "GT")]))


def test_with_no_other_allele_an_mnp_read_reading_one_base_counts(tmp_path):
    """No germline allele: a read showing ?T or A? is the MNP's, as before."""
    contig, bam, v, c = _mnp_locus(tmp_path, somatic=True, germline=False)
    _invariants(c)
    assert c.ad == 13
    assert_matches(c, census(bam, contig, v))


# ── a change both haplotypes carry inside the grown event ────────────────


def _g_run_contig():
    """CAC GGGGGG AGCA at P - 3: a GG>A at P + 3 grows through the G run."""
    contig = list(_contig(seed=72))
    contig[P - 3 : P + 10] = "CACGGGGGGAGCA"
    return "".join(contig)


def test_a_germline_snp_in_the_run_leaves_ref_reads_ref(tmp_path):
    """GG>A at the run's fourth G, a germline G>C at its second on half the REF reads
    (the other haplotype): those REF reads hold the given REF allele and count REF."""
    contig = _g_run_contig()
    snp = P + 1
    hap = contig[: P + 3] + "A" + contig[P + 5 :]
    reads = []
    for i in range(20):
        s = P - 50 + i
        seq = list(contig[s : s + READ])
        if i % 2:
            seq[snp - s] = "C"
        reads.append(make_read(f"ref{i}", "".join(seq), s, ((0, READ),)))
    for i in range(10):
        s = P - 50 + i
        left = P + 3 - s
        reads.append(
            make_read(
                f"alt{i}", hap[s : s + READ], s, ((0, left + 1), (2, 1), (0, READ - left - 1))
            )
        )
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P + 3, "GG", "A")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert (c.rd, c.ad) == (20, 10)
    assert_matches(c, census(bam, contig, v))


def _ca_run_contig():
    """G CACACA T at P: a G>TCA grows through the CA repeat."""
    contig = list(_contig(seed=73))
    contig[P : P + 8] = "GCACACAT"
    return "".join(contig)


def _ca_locus(tmp_path, ref_snp, alt_snp):
    """REF reads and G>TCA carriers, half of each with the germline A>G at the
    repeat's fourth base when `ref_snp` / `alt_snp`."""
    contig = _ca_run_contig()
    snp = P + 4
    hap = contig[:P] + "TCA" + contig[P + 1 :]
    reads = []
    for i in range(20):
        s = P - 50 + i
        seq = list(contig[s : s + READ])
        if ref_snp and i % 2:
            seq[snp - s] = "G"
        reads.append(make_read(f"ref{i}", "".join(seq), s, ((0, READ),)))
    for i in range(20):
        s = P - 50 + i
        left = P - s
        seq = list(hap[s : s + READ])
        if alt_snp and i % 2:
            seq[snp + 2 - s] = "G"  # past the ALT's two extra bases
        reads.append(
            make_read(f"alt{i}", "".join(seq), s, ((0, left + 1), (1, 2), (0, READ - left - 3)))
        )
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, "G", "TCA")])
    (c,) = count_checked(bam, [v])
    return contig, bam, v, c


def test_a_germline_snp_in_the_repeat_on_both_haplotypes_keeps_every_call(tmp_path):
    """The SNP on half the REF reads and half the carriers: background, masked in both
    windows (shifted by the insertion's two bases on the ALT's), so every read keeps
    its call."""
    contig, bam, v, c = _ca_locus(tmp_path, ref_snp=True, alt_snp=True)
    _invariants(c)
    assert (c.rd, c.ad) == (20, 20)
    assert_matches(c, census(bam, contig, v))


def test_a_change_in_the_repeat_on_the_carriers_alone_is_another_allele(tmp_path):
    """The change on half the carriers and no REF read: inside the repeat it may be
    the event's own sequence, so it is not masked and those carriers are another
    allele (partial evidence). The census reads only the substituted bases' anchors,
    so it cannot see this change: the engine's counts are checked directly."""
    _, _, _, c = _ca_locus(tmp_path, ref_snp=False, alt_snp=True)
    _invariants(c)
    assert (c.rd, c.ad, c.partial_alt) == (20, 10, 10)


def test_a_change_both_haplotypes_carry_is_neither_another_allele_nor_a_mismatch(tmp_path):
    """A 33-base tandem duplication with its anchor substituted (x>A + the next 33
    bases), the event grown through both copies, and a change in the copy's REF bases
    on a fifth of every allele's molecules (a separate event): carriers spelling the
    ALT with it are no other allele, so carriers ending before it count ALT; and it is
    background inside the grown event, so REF reads and carriers carrying it keep
    their call."""
    contig = _contig(seed=74)
    dup = contig[P + 1 : P + 34]
    alt = "A" + dup
    hap = contig[:P] + alt + contig[P + 1 :]
    at = P + 1 + 19  # a REF base of the copy; on the ALT haplotype 33 bases further
    reads = []
    for i in range(30):
        s = P - 50 + i
        seq = list(contig[s : s + 150])
        if i % 5 == 0:
            seq[at - s] = "A"
        reads.append(make_read(f"ref{i}", "".join(seq), s, ((0, 150),)))
    for i in range(20):
        s = P - 50 + i
        left = P - s
        seq = list(hap[s : s + 150])
        if i % 5 == 0:
            seq[at + 33 - s] = "A"
        cig = ((0, left), (1, 33), (0, 1), (0, 150 - left - 34))
        reads.append(make_read(f"alt{i}", "".join(seq), s, cig))
    # carriers ending ten bases into the second copy: past the REF window, short of the
    # change
    for i in range(6):
        s = P - 60 + i
        left = P - s
        n = left + 34 + 10
        reads.append(
            make_read(f"short{i}", hap[s : s + n], s, ((0, left), (1, 33), (0, n - left - 33)))
        )
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, contig[P], alt)])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert (c.rd, c.ad, c.partial_alt) == (30, 26, 0)


# GRCh38 chr6:76972991-76973290 (public reference), where HG002 carries a T>C and, on
# the same haplotype, a TG deletion twelve bases on, written as one delins; past it, a
# run of 17 T's.
LOCUS_6 = (
    "ATGACTTTTTTGGCCTGTTAGGAGCTGCATCCAACACACCACAAACTTGACAGTATCAGTTTTCTGTTGTGATGTAATGAACT"
    "ACCCCAACTGATCAGCTTTAAACAACTGATTATTTGCTTAAAATTGTTGGGTGGCATCACTGGGTAG"
    "TTCCTATGCTGTTGTGTTTTTTTTTTTTTTTTTCTGTTTTTCTGGGATCACTCACTCATGTGGTGGTAGTCAGTTGAAATCT"
    "CAGTAGGGCCTGGATGGACTAAGATGACCTCACTCACATGTCTGGCAATTGGTGCTCACTGTAGGCTG"
)


def test_a_carrier_shifted_into_a_run_past_the_event_is_no_ref_molecule(tmp_path):
    """Past the event a molecule shifted by the length change fits REF's bases as well,
    so it is a REF molecule only where it reads a base telling REF from the ALT. The
    carriers' own T after the TG deletion is then no change the REF molecules carry,
    and a carrier reading only that T confidently (its C and the run's end below min
    BQ, as one HG002 read) still tells the ALT from REF by it."""
    contig, p = LOCUS_6, 150
    ref, alt = contig[p : p + 14], "CTCCTATGCTGT"
    hap = contig[:p] + alt + contig[p + 14 :]
    reads = [
        make_read(f"ref{i}", contig[p - 60 + i : p + 40 + i], p - 60 + i, ((0, READ),))
        for i in range(10)
    ]
    for i in range(14):
        s = p - 60 + i
        left = p + 12 - s
        r = make_read(f"alt{i}", hap[s : s + READ], s, ((0, left), (2, 2), (0, READ - left)))
        # as HG002's read: its C, and the run's end where the windows next differ, at Q11
        reads.append(_low(r, [p - s, p + 31 - s, p + 33 - s]) if i >= 10 else r)
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(p, ref, alt)])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert (c.rd, c.ad) == (10, 14)
    assert_matches(c, census(bam, contig, v))
