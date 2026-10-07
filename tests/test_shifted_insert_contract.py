"""C36 #243: a same-length insert of other bases is judged by its bases.

A read whose insertion has the variant's length but other readable bases, written
near the variant (the windowed "same length, other bases" case), went to Phase 3,
which picks the closer of the REF and ALT haplotypes. REF pays for the gap, so
length could win ALT against the read's own bases (Smith-Waterman under bio 3 and
4, PairHMM more often under bio 4), and a read carrying another allele could be
absorbed into REF. On the RC DNA set every Phase-3 read at an insertion row came
this way (50 of 94,530 calls): 2 ALT and 31 partial whose bases read REF.

Contract (operator, 2026-10-06), bases first, never Phase 3:
- the read's bases across the window spell the ALT (an event written one junction
  off with a compensating mismatch): ALT;
- otherwise, inside the discrimination window: another allele, neither with
  partial evidence;
- outside it: a separate event, REF where the window reads REF (RJ-8).
Both backends give the same call.
"""

import random

import pytest
from census import assert_matches, census
from helpers import count_checked, make_read, write_contig

from gbcms import _rs

L = 100
U, R, Q, P = 300, 700, 500, 850  # +8 unique; +CA in (CA)x6; A>AT before C; +A in A10
INS8 = "ATGGACTC"
BACKENDS = ["pairhmm", "sw"]


def _contig():
    rng = random.Random(243)
    c = [rng.choice("ACGT") for _ in range(1000)]
    c[U - 1], c[U], c[U + 1] = "T", "G", "G"  # ...T G | ATGGACTC | G...: no shift
    c[R - 1 : R + 14] = "TG" + "CA" * 6 + "T"  # G>GCA left-aligned at R
    c[Q - 1 : Q + 3] = "GACG"  # A>AT at Q: ...G A | T | C G...
    c[P - 1 : P + 12] = "TG" + "A" * 10 + "C"  # G>GA left-aligned at P
    return "".join(c)


CONTIG = _contig()


def _variant(fa, pos, ref, alt):
    (pv,) = _rs.prepare_variants([_rs.Variant("1", pos, ref, alt, "INS")], fa, 5, False, 1, True)
    v = pv.variant
    assert (pv.gbcms_status, v.pos, v.ref_allele, v.alt_allele) == ("PASS", pos, ref, alt)
    return v


def read(name, i, anchor, ins, junction, *, quals=None, sub=None, start=None, clip=0):
    """A read with `ins` inserted before reference position `junction`. `sub` is one
    (position, base[, quality]) or a list of them, replacing reference bases to write
    compensating mismatches; `start` overrides the read start; `clip` soft-clips the
    read's last bases (aligned bases after the insert become a clip)."""
    s = anchor - 40 - (i % 5) if start is None else start
    left = junction - s
    ref_part = list(CONTIG[s : s + L])
    q_ref = [37] * L
    for one in ([sub] if sub and isinstance(sub[0], int) else sub or []):
        ref_part[one[0] - s] = one[1]
        if len(one) > 2:
            q_ref[one[0] - s] = one[2]
    seq = "".join(ref_part[:left]) + ins + "".join(ref_part[left : L - len(ins)])
    q = q_ref[:left] + [37] * len(ins) + q_ref[left : L - len(ins)]
    for k, qq in enumerate(quals or []):
        q[left + k] = qq
    right = L - left - len(ins)
    cigar = ((0, left), (1, len(ins)), (0, right - clip)) + (((4, clip),) if clip else ())
    if right == clip:
        cigar = ((0, left), (1, len(ins)), (4, clip))
    return make_read(name, seq, s, cigar, flag=16 * (i % 2), quals=q)


def ref_read(name, i, anchor):
    s = anchor - 40 - (i % 5)
    return make_read(name, CONTIG[s : s + L], s, ((0, L),), flag=16 * (i % 2))


def _count(tmp_path, name, reads, pos, ref, alt, backend, with_census=True):
    fa, bam = write_contig(tmp_path, CONTIG, reads, name)
    v = _variant(fa, pos, ref, alt)
    (counts,) = count_checked(bam, [v], alignment_backend=backend)
    if with_census:
        assert_matches(counts, census(bam, CONTIG, v), name)
    return counts.rd, counts.ad, counts.partial_alt


@pytest.mark.parametrize("backend", BACKENDS)
def test_other_bases_outside_the_window_are_a_separate_event(tmp_path, backend):
    """A non-repeat +8: inserts of other bases one junction left of the variant
    (the review's shape, partly masked) and four junctions right (all readable) lie
    outside its discrimination window; the window reads REF: REF, as a readable
    insert there is."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS8
    reads = [
        read(f"l{i}", i, a, "TNNACTCA", a, quals=[30, 2, 2, 12, 10, 30, 21, 2]) for i in range(4)
    ]
    reads += [read(f"r{i}", i, a, "TCCAGTAA", a + 5) for i in range(4)]
    reads += [read(f"c{i}", i, a, INS8, a + 1) for i in range(4)]
    assert _count(tmp_path, f"o{backend}", reads, a, ref, alt, backend) == (8, 4, 0)


@pytest.mark.parametrize("backend", BACKENDS)
def test_other_bases_inside_a_tract_are_another_allele(tmp_path, backend):
    """G>GCA before (CA)x6: a +GT written at a junction inside the tract is another
    allele (neither, partial), never REF or ALT; a +CA there is the variant."""
    a, ref, alt = R, "G", "GCA"
    reads = [read(f"g{i}", i, a, "GT", a + 5) for i in range(4)]
    reads += [read(f"c{i}", i, a, "CA", a + 5) for i in range(4)]
    reads += [ref_read(f"r{i}", i, a) for i in range(4)]
    assert _count(tmp_path, f"t{backend}", reads, a, ref, alt, backend) == (4, 4, 4)


@pytest.mark.parametrize("backend", BACKENDS)
def test_the_alt_written_one_junction_off_is_alt(tmp_path, backend):
    """The ALT molecule written with the gap one junction right and a compensating
    mismatch (the reference G after the anchor shown as the ALT's first inserted A):
    its bases across the window spell the ALT, so it is ALT however the gap was
    placed. No census: it trusts the aligned flank, which here holds the mismatch."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS8
    body = INS8[1:] + CONTIG[a + 1]
    reads = [read(f"w{i}", i, a, body, a + 2, sub=(a + 1, INS8[0])) for i in range(4)]
    assert _count(tmp_path, f"w{backend}", reads, a, ref, alt, backend, with_census=False) == (
        0,
        4,
        0,
    )


# ── From the review: judged where the read's own insert can sit ───────────────
# The aligner's placement is one of several equally good ones: an insert slides one
# junction when the base it moves onto fits as well as the base aligned there did.
# A read whose insert can sit at the variant's junction is judged there, as the
# strict path would; one that can sit inside the window is another allele. The
# census trusts the aligned flank, so it cannot judge these shapes.


@pytest.mark.parametrize("backend", BACKENDS)
def test_an_alt_read_with_an_error_written_off_junction_is_another_allele(tmp_path, backend):
    """The ALT with a sequencing error at its first inserted base (A>G) is written
    one junction right with no mismatch, and with one at its last base (C>G) one
    junction left: both slide onto the junction with a confident mismatch, another
    allele, as the same reads written at the junction are."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS8
    reads = [read(f"f{i}", i, a, "TGGACTCG", a + 2) for i in range(4)]
    reads += [read(f"l{i}", i, a, "GATGGACT", a) for i in range(4)]
    reads += [read(f"j{i}", i, a, "GTGGACTC", a + 1) for i in range(4)]
    got = _count(tmp_path, f"e{backend}", reads, a, ref, alt, backend, with_census=False)
    assert got == (0, 0, 12)


@pytest.mark.parametrize("backend", BACKENDS)
def test_a_masked_compensating_base_slides_back_to_the_alt(tmp_path, backend):
    """The ALT with its first inserted base masked (N at Q2), written one junction
    right with the N aligned on the base after the anchor: slid back, its seven
    readable inserted bases are the ALT's, as at the junction (RJ-20)."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS8
    body = INS8[1:] + CONTIG[a + 1]
    reads = [read(f"m{i}", i, a, body, a + 2, sub=(a + 1, "N", 2)) for i in range(4)]
    reads += [read(f"j{i}", i, a, "N" + INS8[1:], a + 1, quals=[2]) for i in range(4)]
    got = _count(tmp_path, f"m{backend}", reads, a, ref, alt, backend, with_census=False)
    assert got == (0, 8, 0)


@pytest.mark.parametrize("backend", BACKENDS)
def test_carriers_written_further_off_are_alt(tmp_path, backend):
    """ALT carriers the aligner wrote two junctions off with two compensating
    mismatches, one written one junction off that ends in a clip after its insert,
    and one that starts on the anchor: each slides back onto the junction as the ALT."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS8
    two = INS8[2:] + CONTIG[a + 1] + CONTIG[a + 2]
    one = INS8[1:] + CONTIG[a + 1]
    reads = [
        read(f"t{i}", i, a, two, a + 3, sub=[(a + 1, INS8[0]), (a + 2, INS8[1])]) for i in range(4)
    ]
    reads += [
        read(f"c{i}", i, a, one, a + 2, sub=(a + 1, INS8[0]), start=a - 60 - i, clip=30 - i)
        for i in range(4)
    ]
    reads += [read(f"s{i}", i, a, one, a + 2, sub=(a + 1, INS8[0]), start=a) for i in range(4)]
    got = _count(tmp_path, f"f{backend}", reads, a, ref, alt, backend, with_census=False)
    assert got == (0, 12, 0)


@pytest.mark.parametrize("backend", BACKENDS)
def test_another_allele_written_outside_the_window_is_not_ref(tmp_path, backend):
    """A>AT before a C: a read inserting a C after the C (the A>AC allele, written
    one junction right) slides onto the junction: another allele, never REF."""
    a, ref, alt = Q, "A", "AT"
    reads = [read(f"c{i}", i, a, "C", a + 2) for i in range(4)]
    reads += [ref_read(f"r{i}", i, a) for i in range(4)]
    got = _count(tmp_path, f"s{backend}", reads, a, ref, alt, backend, with_census=False)
    assert got == (4, 0, 4)


# ── Second review: readability is the aligner's inserted bases' (RJ-20) ───────
@pytest.mark.parametrize("backend", BACKENDS)
def test_a_readable_insert_does_not_become_unreadable_by_sliding(tmp_path, backend):
    """G>GA before A10: an 11-A read with one run base masked (Q5, or N at Q2) and
    its readable inserted A written right of the mask. Its own inserted base is
    read, so it is ALT (RJ-20), as written at the junction; a slide that ends on the
    masked base must not make it unreadable."""
    a, ref, alt = P, "G", "GA"
    reads = [read(f"q{i}", i, a, "A", a + 8, sub=(a + 5, "A", 5)) for i in range(4)]
    reads += [read(f"n{i}", i, a, "A", a + 11, sub=(a + 5, "N", 2)) for i in range(4)]
    reads += [ref_read(f"r{i}", i, a) for i in range(4)]
    assert _count(tmp_path, f"b1{backend}", reads, a, ref, alt, backend) == (4, 8, 0)


@pytest.mark.parametrize("backend", BACKENDS)
def test_a_carrier_with_a_masked_end_base_slides_back(tmp_path, backend):
    """The ALT written one junction right with a compensating mismatch, its last
    inserted base (where the reference base after the anchor belongs) misread as T
    at Q5: masked, it fits; its readable inserted bases are the ALT's once it slides
    back, as at the junction."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS8
    body = INS8[1:] + "T"
    reads = [
        read(f"m{i}", i, a, body, a + 2, sub=(a + 1, INS8[0]), quals=[37] * 7 + [5])
        for i in range(4)
    ]
    got = _count(tmp_path, f"b2{backend}", reads, a, ref, alt, backend, with_census=False)
    assert got == (0, 4, 0)
