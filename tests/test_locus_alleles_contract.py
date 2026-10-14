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

import pytest
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
    """Past the event a run lets a molecule shifted by the length change fit REF's
    bases as well, so a molecule is REF's only where its bases, the event's included,
    favour REF over the ALT placed the same way. The carriers' own T after the TG
    deletion is then no change the REF molecules carry, and a carrier reading only that
    T confidently (its C and the run's end below min BQ, as one HG002 read) still tells
    the ALT from REF by it."""
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


# ── review findings (2026-10-14) ──────────────────────────────────────────


def _ca_reads(contig, ref_snp=(), alt_snp=(), alt_change=None, n=20):
    """REF reads and G>TCA carriers at the CA repeat: the reads numbered in `ref_snp` /
    `alt_snp` carry the germline A>G at the repeat's fourth base; with `alt_change` (a
    REF offset and a base), every carrier carries that change too."""
    snp = P + 4
    hap = contig[:P] + "TCA" + contig[P + 1 :]
    reads = []
    for i in range(n):
        s = P - 50 + i
        seq = list(contig[s : s + READ])
        if i in ref_snp:
            seq[snp - s] = "G"
        reads.append(make_read(f"ref{i}", "".join(seq), s, ((0, READ),)))
    for i in range(n):
        s = P - 50 + i
        left = P - s
        seq = list(hap[s : s + READ])
        if i in alt_snp:
            seq[snp + 2 - s] = "G"
        if alt_change:
            seq[alt_change[0] + 2 - s] = alt_change[1]
        reads.append(
            make_read(f"alt{i}", "".join(seq), s, ((0, left + 1), (1, 2), (0, READ - left - 3)))
        )
    return reads


def _other_base(contig, at):
    return next(b for b in "ACGT" if b not in contig[at - 1 : at + 2])


def test_an_allele_differing_from_the_alt_only_at_masked_bases_is_not_another(tmp_path):
    """A co-annotated row merging the event with the germline SNP the REF molecules
    carry in the repeat (masked background) differs from the ALT only there: it is the
    ALT with a separate event, so the carriers stay the ALT's."""
    contig = _ca_run_contig()
    reads = _ca_reads(contig, ref_snp=range(1, 20, 2))
    fa, bam = write_contig(tmp_path, contig, reads)
    vs = _prepared(fa, [(P, "G", "TCA"), (P, contig[P : P + 5], "TCACACG")])
    counts = count_checked(bam, vs, sibling_variants=[[vs[1]], [vs[0]]])
    _invariants(counts[0])
    assert (counts[0].rd, counts[0].ad) == (20, 20)


def test_masked_background_and_an_isolated_change_leave_every_carrier_alt(tmp_path):
    """The germline SNP in the repeat on half the REF reads and carriers, and every
    carrier with an isolated change two bases past the event (its own, masked and
    named): an allele differing only at those masked bases is no other allele."""
    contig = _ca_run_contig()
    change = (P + 8, _other_base(contig, P + 8))
    reads = _ca_reads(contig, ref_snp=range(1, 20, 2), alt_snp=range(1, 20, 2), alt_change=change)
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, "G", "TCA")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert c.ad == 20


def test_background_on_every_carrier_beside_a_phased_change_leaves_them_alt(tmp_path):
    """Every carrier carries the germline SNP (on half the REF reads too) and an
    isolated change of its own two bases past the event: the guard reads past the event
    through the masked repeat, so it masks that change and the carriers count ALT."""
    contig = _ca_run_contig()
    change = (P + 8, _other_base(contig, P + 8))
    reads = _ca_reads(contig, ref_snp=range(1, 20, 2), alt_snp=range(20), alt_change=change)
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, "G", "TCA")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert c.ad == 20


def test_a_second_snp_past_the_event_leaves_the_background_found(tmp_path):
    """The REF molecules carrying the repeat's germline SNP carry another in phase three
    bases past the event: one change past the event does not stop them being read as
    REF molecules, so the background is still found and every REF read counts."""
    contig = _ca_run_contig()
    far = P + 10
    reads = _ca_reads(contig, ref_snp=range(1, 20, 2))
    new = _other_base(contig, far)
    for i in range(1, 20, 2):
        r = reads[i]
        seq = list(r.query_sequence)
        seq[far - r.reference_start] = new
        quals = r.query_qualities
        r.query_sequence = "".join(seq)
        r.query_qualities = quals
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, "G", "TCA")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert (c.rd, c.ad) == (20, 20)


def test_a_germline_allele_with_a_second_snp_past_the_mnp_is_still_another(tmp_path):
    """The germline GT carries a second SNP seven bases past the block: it is still
    another allele, so with no AT molecule its reads showing only the T are no ALT."""
    contig = _mnp_contig()
    far = P + 7
    new = _other_base(contig, far)
    reads = _mnp_reads(contig, "GC", "ref", 20)
    germ = _mnp_reads(contig, "GT", "germ", 20) + _mnp_reads(
        contig, "GT", "germlow", 4, first=20, low_at=0
    )
    for r in germ:
        seq, quals = list(r.query_sequence), r.query_qualities
        seq[far - r.reference_start] = new
        r.query_sequence = "".join(seq)
        r.query_qualities = quals
    fa, bam = write_contig(tmp_path, contig, reads + germ)
    (v,) = _prepared(fa, [(P, "GC", "AT")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert c.ad == 0


def test_a_recurring_sequencing_error_is_not_another_allele(tmp_path):
    """CC>TT at 300 REF fragments, 3 of them (1%) with a C>T error spelling CT, under the
    share an allele the reads spell must reach: the 5 carriers whose first T is below
    min BQ count ALT, at any depth."""
    contig = list(_contig(seed=43))
    contig[P - 1 : P + 3] = "GCCA"
    contig = "".join(contig)
    reads = []
    for i in range(300):
        s = P - 60 + (i % 50)
        block = "CT" if i < 3 else "CC"
        reads.append(
            make_read(f"ref{i}", contig[s:P] + block + contig[P + 2 : s + READ], s, ((0, READ),))
        )
    for i in range(20):
        s = P - 60 + (i % 50)
        r = make_read(f"alt{i}", contig[s:P] + "TT" + contig[P + 2 : s + READ], s, ((0, READ),))
        reads.append(_low(r, [P - s]) if i < 5 else r)
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, "CC", "TT")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert c.ad == 20


def test_an_mnp_read_spelling_another_allele_where_the_alt_keeps_the_reference_is_not_alt(
    tmp_path,
):
    """GCA>ACT beside AGT on 6 fragments, read plainly: AGT differs from the ALT only at
    the base the ALT keeps, which the MNP check does not compare; RJ-25 does, so those
    reads are another allele, not the ALT."""
    contig = list(_contig(seed=41))
    contig[P : P + 3] = "GCA"
    contig = "".join(contig)
    reads = []
    for tag, block, n in (("ref", "GCA", 20), ("alt", "ACT", 10), ("oth", "AGT", 6)):
        for i in range(n):
            s = P - 50 + i
            reads.append(
                make_read(
                    f"{tag}{i}", contig[s:P] + block + contig[P + 3 : s + READ], s, ((0, READ),)
                )
            )
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, "GCA", "ACT")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert (c.rd, c.ad) == (20, 10)
    assert_matches(c, census(bam, contig, v, others=[(P, "GCA", "AGT")]))


def test_an_mnp_read_set_aside_for_another_allele_keeps_its_n(tmp_path):
    """GC>AT beside the germline GT: reads with N for the A are another allele's or the
    ALT's, and still count toward the N reads."""
    contig = _mnp_contig()
    reads = _mnp_reads(contig, "GC", "ref", 20) + _mnp_reads(contig, "GT", "germ", 20)
    reads += _mnp_reads(contig, "AT", "som", 6, first=40)
    for r in _mnp_reads(contig, "NT", "somN", 4, first=46):
        reads.append(r)
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, "GC", "AT")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert (c.ad, c.n_count) == (6, 4)


def test_the_alleles_at_a_locus_do_not_depend_on_bin_geometry(tmp_path):
    """G>AT after a 30-A run grows left past the variant's read window. A read aligned
    there with the event in its clip is in a bin's reads for a neighbouring row, not in
    this variant's: the other alleles are found from the variant's own reads, so the
    counts match one variant per bin (count_checked) and AC on two fragments is none."""
    contig = list(_contig(seed=91))
    contig[P - 30 : P] = "A" * 30
    contig[P] = "G"
    contig = "".join(contig)
    hap, oth = (contig[:P] + b + contig[P + 1 :] for b in ("AT", "AC"))
    reads = [
        make_read(f"ref{i}", contig[P - 50 + i : P + 50 + i], P - 50 + i, ((0, READ),))
        for i in range(10)
    ]
    for tag, h, n in (("alt", hap, 16), ("oth", oth, 2)):
        for i in range(n):
            s = P - 50 + i
            left = P - s + 1
            r = make_read(
                f"{tag}{i}", h[s : s + READ], s, ((0, left), (1, 1), (0, READ - left - 1))
            )
            reads.append(_low(r, [left]) if tag == "alt" and i >= 10 else r)
    reads.append(make_read("othclip", oth[P - 40 : P + 20], P - 40, ((0, 30), (4, 30))))
    fa, bam = write_contig(tmp_path, contig, reads)
    snv = contig[P - 100]
    vs = _prepared(fa, [(P - 100, snv, "A" if snv != "A" else "C"), (P, "G", "AT")])
    counts = count_checked(bam, vs)
    _invariants(counts[1])
    assert (counts[1].rd, counts[1].ad, counts[1].partial_alt) == (10, 16, 2)


# ── second review (2026-10-14) ────────────────────────────────────────────


def _locus6_reads(snp_at, new, on_ref, on_alt):
    """LOCUS_6's delins, REF reads and carriers (20 each), a germline change at `snp_at`
    (REF offset) on half the REF reads and / or half the carriers."""
    contig, p = LOCUS_6, 150
    alt = "CTCCTATGCTGT"
    d = len(alt) - 14
    hap = contig[:p] + alt + contig[p + 14 :]
    reads = []
    for i in range(20):
        s = p - 60 + i
        seq = list(contig[s : s + READ])
        if on_ref and i % 2:
            seq[snp_at - s] = new
        reads.append(make_read(f"ref{i}", "".join(seq), s, ((0, READ),)))
    for i in range(20):
        s = p - 60 + i
        left = p + 12 - s
        seq = list(hap[s : s + READ])
        if on_alt and i % 2:
            seq[snp_at + d - s] = new
        reads.append(make_read(f"alt{i}", "".join(seq), s, ((0, left), (2, 2), (0, READ - left))))
    return contig, p, alt, reads


@pytest.mark.parametrize("on_alt", [False, True])
def test_a_germline_snp_in_the_repeat_before_a_run_is_background(tmp_path, on_alt):
    """LOCUS_6: past the event a run tells neither allele from the other, so REF
    molecules are told by the bases inside the event; a germline G>A in the TG repeat
    the event grew through, on half the REF reads (and carriers), is then background."""
    contig, p, alt, reads = _locus6_reads(165, "A", True, on_alt)
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(p, contig[p : p + 14], alt)])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert (c.rd, c.ad) == (20, 20)


def test_a_germline_snp_beside_a_large_deletion_is_background(tmp_path):
    """A 60-base delins grown through a G run, a germline G>C in the run on half the REF
    reads: the ALT's haplotype ends before the bases past the event, so REF molecules
    are told by their own deleted bases."""
    p, n, alt = 500, 60, "ATCTA"
    c = list(_contig(n=1200, seed=76))
    c[p - 6 : p + 1] = "GGGGGGG"
    contig = "".join(c)
    hap = contig[:p] + alt + contig[p + n :]
    reads = []
    for i in range(20):
        s = p - 30 + i
        seq = list(contig[s : s + 150])
        if i % 2:
            seq[p - 4 - s] = "C"
        reads.append(make_read(f"ref{i}", "".join(seq), s, ((0, 150),)))
    for i in range(20):
        s = p - 30 + i
        left = p - s
        cig = ((0, left), (1, len(alt)), (2, n), (0, 150 - left - len(alt)))
        reads.append(make_read(f"alt{i}", hap[s : s + 150], s, cig))
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(p, contig[p : p + n], alt)])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert (c.rd, c.ad) == (20, 20)


@pytest.mark.parametrize("half", [False, True])
def test_an_mnp_with_a_phased_change_beside_its_bases_keeps_its_carriers(tmp_path, half):
    """CC>TT where the carriers carry CCC>TTT (the input under-describes a UV triplet):
    the third T is beside the MNP's own bases, a phased neighbour, not another allele,
    so every carrier counts ALT, as beside any event."""
    contig = list(_contig(seed=34))
    contig[P - 2 : P + 4] = "GACCCA"
    contig = "".join(contig)
    reads = [
        make_read(f"ref{i}", contig[P - 50 + i : P + 50 + i], P - 50 + i, ((0, READ),))
        for i in range(20)
    ]
    for i in range(10):
        s = P - 50 + i
        block = "TTT" if not half or i % 2 else "TTC"
        reads.append(
            make_read(f"alt{i}", contig[s:P] + block + contig[P + 3 : s + READ], s, ((0, READ),))
        )
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, "CC", "TT")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert (c.rd, c.ad) == (20, 10)


def test_a_germline_allele_with_a_phased_dnp_past_the_mnp_is_still_another(tmp_path):
    """The germline GT carries a DNP five bases past the block: still another allele."""
    contig = _mnp_contig()
    reads = _mnp_reads(contig, "GC", "ref", 20)
    germ = _mnp_reads(contig, "GT", "germ", 20) + _mnp_reads(
        contig, "GT", "germlow", 4, first=20, low_at=0
    )
    for r in germ:
        seq, quals = list(r.query_sequence), r.query_qualities
        for at in (P + 5, P + 6):
            seq[at - r.reference_start] = _other_base(contig, at)
        r.query_sequence = "".join(seq)
        r.query_qualities = quals
    fa, bam = write_contig(tmp_path, contig, reads + germ)
    (v,) = _prepared(fa, [(P, "GC", "AT")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert c.ad == 0


def test_a_germline_dnp_past_the_event_leaves_the_background_found(tmp_path):
    """The REF molecules carrying the repeat's germline SNP carry a DNP in phase past the
    event: two changes beside it, still REF molecules, the background still found."""
    contig = _ca_run_contig()
    reads = _ca_reads(contig, ref_snp=range(1, 20, 2))
    for i in range(1, 20, 2):
        r = reads[i]
        seq, quals = list(r.query_sequence), r.query_qualities
        for at in (P + 10, P + 11):
            seq[at - r.reference_start] = _other_base(contig, at)
        r.query_sequence = "".join(seq)
        r.query_qualities = quals
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, "G", "TCA")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert (c.rd, c.ad) == (20, 20)


@pytest.mark.parametrize("n_ref", [50, 300, 1000])
def test_another_allele_as_common_as_a_low_alt_is_found_at_any_depth(tmp_path, n_ref):
    """The 28-to-6 delins and the allele one inserted base apart, each on 16 fragments
    (6 of the carriers with that base below min BQ), beside 50 to 1000 REF reads: the
    other allele recurs far above the ALT's errors, so it is found at any depth."""
    contig = _contig()
    ref = contig[P : P + REF_LEN]
    reads = [
        make_read(
            f"ref{i}",
            contig[P - 60 + (i % 40) : P + 40 + (i % 40)],
            P - 60 + (i % 40),
            ((0, READ),),
        )
        for i in range(n_ref)
    ]
    reads += _delins_reads(contig, ALT, "alt", 10)
    reads += _delins_reads(contig, ALT, "low", 6, first=10, low=True)
    reads += _delins_reads(contig, OTHER, "oth", 16, first=20)
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, ref, ALT)])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert c.ad == 10


# ── third review (2026-10-14) ─────────────────────────────────────────────


def _uv_contig():
    """GACCCA at P - 2: a CC>TT at P beside a third C."""
    contig = list(_contig(seed=34))
    contig[P - 2 : P + 4] = "GACCCA"
    return "".join(contig)


def _delins_errors(tmp_path, flank_change, ref_carry, n_alt, n_err, n_low):
    """The 28-to-6 delins on `n_alt` carriers: `n_err` with a confident error at the
    third inserted base, `n_low` more reading it below min BQ; with `flank_change`, every
    carrier carries a change at P - 1 (in cis; with `ref_carry`, half the REF reads too)."""
    contig = _contig()
    ref = contig[P : P + REF_LEN]
    hap = list(contig[:P] + ALT + contig[P + REF_LEN :])
    if flank_change:
        hap[P - 1] = _other_base(contig, P - 1)
    hap = "".join(hap)
    reads = []
    for i in range(40):
        s = P - 60 + i
        seq = list(contig[s : s + READ])
        if ref_carry and i % 2:
            seq[P - 1 - s] = hap[P - 1]
        reads.append(make_read(f"ref{i}", "".join(seq), s, ((0, READ),)))
    for i in range(n_alt):
        s = P - 60 + (i % 40)
        left = P - s
        seq = list(hap[s : s + READ])
        if i < n_err:
            seq[left + TELL] = "C"
        cig = ((0, left), (1, len(ALT)), (2, REF_LEN), (0, READ - left - len(ALT)))
        r = make_read(f"alt{i}", "".join(seq), s, cig)
        reads.append(_low(r, [left + TELL]) if n_err <= i < n_err + n_low else r)
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, ref, ALT)])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    return c


@pytest.mark.parametrize("flank_change,ref_carry", [(False, False), (True, False), (True, True)])
def test_a_recurring_error_among_carriers_sharing_a_masked_change_is_no_allele(
    tmp_path, flank_change, ref_carry
):
    """200 carriers, 6 with one error at an inserted base: with every carrier carrying a
    masked change in the window (a phased neighbour, a germline SNP in cis), they still
    spell the ALT, so 6 is under 5% of them and the 10 reading that base below min BQ
    count ALT."""
    c = _delins_errors(tmp_path, flank_change, ref_carry, n_alt=200, n_err=6, n_low=10)
    assert c.ad == 194


def test_a_recurring_error_among_deep_carriers_is_no_allele(tmp_path):
    """1000 carriers with a masked change in cis, 3 with one error: no allele."""
    c = _delins_errors(tmp_path, True, False, n_alt=1000, n_err=3, n_low=30)
    assert c.ad == 997


def test_a_recurring_error_among_mnp_carriers_with_a_neighbour_is_no_allele(tmp_path):
    """CC>TT carriers spelling TTT (a phased neighbour beside the MNP), 6 of 200 with an
    error at the first T: under 5% of the carriers, no allele."""
    contig = _uv_contig()
    reads = [
        make_read(f"ref{i}", contig[P - 50 + i : P + 50 + i], P - 50 + i, ((0, READ),))
        for i in range(40)
    ]
    for i in range(200):
        s = P - 50 + (i % 40)
        block = "GTT" if i < 6 else "TTT"
        r = make_read(f"alt{i}", contig[s:P] + block + contig[P + 3 : s + READ], s, ((0, READ),))
        reads.append(_low(r, [P - s]) if 6 <= i < 16 else r)
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, "CC", "TT")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert c.ad == 194


@pytest.mark.parametrize("hom_at", [P - 1, P + 2])
def test_ref_errors_beside_a_homozygous_snp_are_no_allele(tmp_path, hom_at):
    """GC>AT with a homozygous SNP beside it (masked background, on every read): 4 of
    1000 REF reads with an error spelling AC are no allele, so the 10 carriers reading
    only the A count ALT."""
    contig = _mnp_contig()
    hom = _other_base(contig, hom_at)
    reads = []
    for i in range(1000):
        s = P - 50 + (i % 40)
        seq = list(contig[s : s + READ])
        seq[hom_at - s] = hom
        if i < 4:
            seq[P - s] = "A"
        reads.append(make_read(f"ref{i}", "".join(seq), s, ((0, READ),)))
    for i in range(30):
        s = P - 50 + (i % 40)
        seq = list(contig[s : s + READ])
        seq[hom_at - s] = hom
        seq[P - s : P - s + 2] = "AT"
        r = make_read(f"alt{i}", "".join(seq), s, ((0, READ),))
        reads.append(_low(r, [P + 1 - s]) if i < 10 else r)
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, "GC", "AT")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert c.ad == 30


@pytest.mark.parametrize("shape", ["plain", "clip", "indel", "all indel"])
@pytest.mark.parametrize("block", ["TTC", "TTT"])
def test_an_mnp_carrier_with_a_neighbour_counts_however_it_is_aligned(tmp_path, block, shape):
    """CC>TT carriers spelling TTC, or TTT (a phased neighbour beside the MNP's own
    bases): 6 aligned plainly and 4 aligned as `shape` (the block's end clipped, or an
    indel at the block; or all 10 with the indel), all ALT: the exact-carrier route reads
    beside an MNP as beside any event, the neighbour found from the carriers' bases."""
    contig = _uv_contig()
    reads = [
        make_read(f"ref{i}", contig[P - 50 + i : P + 50 + i], P - 50 + i, ((0, READ),))
        for i in range(20)
    ]
    for i in range(10):
        s = P - 60 + i
        seq = contig[s:P] + block + contig[P + 3 : s + READ]
        kind = "indel" if shape == "all indel" else shape if i >= 6 else "plain"
        if kind == "plain":
            reads.append(make_read(f"alt{i}", seq, s, ((0, READ),)))
        elif kind == "clip":
            e = P + 7
            short = contig[s:P] + block + contig[P + 3 : e]
            reads.append(make_read(f"alt{i}", short, s, ((0, P - s + 1), (4, e - P - 1))))
        else:
            cig = ((0, P - s), (1, 1), (0, 2), (2, 1), (0, len(seq) - (P - s) - 3))
            reads.append(make_read(f"alt{i}", seq, s, cig))
    fa, bam = write_contig(tmp_path, contig, reads)
    (v,) = _prepared(fa, [(P, "CC", "TT")])
    (c,) = count_checked(bam, [v])
    _invariants(c)
    assert c.ad == 10
