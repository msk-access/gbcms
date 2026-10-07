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
U, R = 300, 700  # a non-repeat +8, and a +CA in a (CA)x6 tract
INS8 = "ATGGACTC"
BACKENDS = ["pairhmm", "sw"]


def _contig():
    rng = random.Random(243)
    c = [rng.choice("ACGT") for _ in range(1000)]
    c[U - 1], c[U], c[U + 1] = "T", "G", "G"  # ...T G | ATGGACTC | G...: no shift
    c[R - 1 : R + 14] = "TG" + "CA" * 6 + "T"  # G>GCA left-aligned at R
    return "".join(c)


CONTIG = _contig()


def _variant(fa, pos, ref, alt):
    (pv,) = _rs.prepare_variants([_rs.Variant("1", pos, ref, alt, "INS")], fa, 5, False, 1, True)
    v = pv.variant
    assert (pv.gbcms_status, v.pos, v.ref_allele, v.alt_allele) == ("PASS", pos, ref, alt)
    return v


def read(name, i, anchor, ins, junction, *, quals=None, sub=None):
    """A read with `ins` inserted before reference position `junction`; `sub`
    replaces the reference base at (position, base) to write a compensating
    mismatch."""
    s = anchor - 40 - (i % 5)
    left = junction - s
    ref_part = list(CONTIG[s : s + L])
    if sub:
        ref_part[sub[0] - s] = sub[1]
    seq = "".join(ref_part[:left]) + ins + "".join(ref_part[left : L - len(ins)])
    q = [37] * L
    for k, qq in enumerate(quals or []):
        q[left + k] = qq
    cigar = ((0, left), (1, len(ins)), (0, L - left - len(ins)))
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


@pytest.mark.xfail(strict=True, reason="Phase 3 picks the closer haplotype")
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


@pytest.mark.xfail(strict=True, reason="Phase 3 picks the closer haplotype")
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
