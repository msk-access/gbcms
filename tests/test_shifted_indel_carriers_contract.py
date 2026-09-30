"""A carrier is ALT wherever its aligner places an equivalent indel (#189).

An indel in a repeat can be written at any junction of its shift region with the
same haplotype. The input is left-aligned; aligners and callers may not be. The
windowed checks accepted a shifted placement only by a proxy:

- an insertion when the reference base before it equalled the input's anchor base,
  never true inside the repeat, so carriers with the input's inserted bases were
  counted **REF** (rotated bases, as `AC` for `CA`, failed the base check and
  reached the haplotype fallback, which counted them);
- a deletion when the reference bases it removes equalled the input's deleted
  bases, so a rotated placement (removing `AC` for `CA`) under 5bp was counted
  REF.

A placement is now accepted when it gives the input's haplotype: the read's own
bases carry the allele. Placements that do not (other bases, outside the repeat,
a read with another gap or a splice between the two placements, an ALT that also
substitutes its anchor base) are not ALT; the variant's bases inserted outside
its discrimination window are a separate event, and the read is REF. Every
placement of the shift region is scanned, however far it slides.

Committed red (xfail-strict) before the implementation; flipped green with it.
"""

import random

import pysam
import pytest
from helpers import count_both, make_read

from gbcms import _rs as gbcms_rs

RED = pytest.mark.xfail(strict=True, reason="read-side and deletion-reach gaps found in review")

READ, L, A = 100, 800, 400  # read length, contig length, the repeat's anchor (0-based)


def _contig(motif):
    rng = random.Random(189)
    seq = [rng.choice("CGT") for _ in range(L)]
    seq[A : A + len(motif)] = list(motif)
    return "".join(seq)


HOMOPOLYMER = _contig("GAAAAAT")  # G, A x5 at A+1..A+5, T
STR = _contig("GCACACACAT")  # G, (CA) x4 at A+1..A+8, T
ANCHOR_SUB = _contig("ACCT")  # A at A, C C at A+1..A+2, T
DUP = "ACGTTGCA"  # no motif of 6 or fewer bases repeats in it
DUPLICATION = _contig("G" + DUP + DUP + "T")  # G, two copies at A+1..A+16, T
SPAN10 = "ATCGGATCTA"  # a unique 10bp stretch (slides neither way between G and C)
SPAN60 = "A" + "".join(random.Random(60).choice("ACGT") for _ in range(58)) + "A"
UNIQUE10 = _contig("G" + SPAN10 + "C")
UNIQUE60 = _contig("G" + SPAN60 + "C")


def _files(tmp_path, ref, reads):
    fa = tmp_path / "ref.fa"
    fa.write_text(f">1\n{ref}\n")
    pysam.faidx(str(fa))
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": L}]}
    raw = tmp_path / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=hdr) as fh:
        for r in reads:
            fh.write(r)
    bam = tmp_path / "s.bam"
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    return str(fa), str(bam)


def _count_reads(tmp_path, ref, ref_allele, alt_allele, carrier):
    """Five REF reads and five carriers, `carrier(s)` giving each carrier's
    (sequence, CIGAR) from its start `s`; counted in both paths (count_both asserts
    parity)."""
    starts = range(A - 50, A - 45)
    reads = [make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),)) for i, s in enumerate(starts)]
    for i, s in enumerate(starts):
        seq, cigar = carrier(s)
        reads.append(make_read(f"a{i}", seq, s, cigar))
    fa, bam = _files(tmp_path, ref, reads)
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", A, ref_allele, alt_allele, "X")], fa, 5, False, 1, True
    )
    assert pv.gbcms_status == "PASS" and pv.variant.pos == A, pv.gbcms_status_reason
    (c,) = count_both(bam, [pv.variant])
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    return c.rd, c.ad, c.partial_alt


def _count(tmp_path, ref, ref_allele, alt_allele, junction, op, bases):
    """Carriers whose CIGAR places the event (`op` 'I' with inserted `bases`, or 'D'
    of `bases` reference bases) before `junction`."""

    def carrier(s):
        left = junction - s
        if op == "I":
            seq = (ref[s:junction] + bases + ref[junction:])[:READ]
            return seq, ((0, left), (1, len(bases)), (0, READ - left - len(bases)))
        seq = (ref[s:junction] + ref[junction + bases :])[:READ]
        return seq, ((0, left), (2, bases), (0, READ - left))

    return _count_reads(tmp_path, ref, ref_allele, alt_allele, carrier)


def _ops(ref, events):
    """A carrier from ordered events, each placed before its reference position:
    (pos, "I", bases), (pos, "D", n) or (pos, "N", n)."""

    def carrier(s):
        seq, cigar, p = "", [], s
        for pos, op, x in events:
            if pos > p:
                seq += ref[p:pos]
                cigar.append((0, pos - p))
                p = pos
            if op == "I":
                seq += x
                cigar.append((1, len(x)))
            else:
                cigar.append((2 if op == "D" else 3, x))
                p += x
        rest = READ - len(seq)
        return seq + ref[p : p + rest], (*cigar, (0, rest))

    return carrier


def _j(n):
    """Junction offsets 1..n after the anchor."""
    return list(range(1, n + 1))


# ── Insertions ────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("junction", _j(6))
def test_homopolymer_insertion_carriers_at_every_junction(tmp_path, junction):
    """G>GA in G AAAAA T: an A inserted anywhere in the run is the same allele."""
    assert _count(tmp_path, HOMOPOLYMER, "G", "GA", A + junction, "I", "A") == (5, 5, 0)


def test_a_two_base_homopolymer_insertion_placed_mid_run(tmp_path):
    assert _count(tmp_path, HOMOPOLYMER, "G", "GAA", A + 3, "I", "AA") == (5, 5, 0)


@pytest.mark.parametrize("junction", _j(9))
def test_str_insertion_carriers_at_every_phase(tmp_path, junction):
    """G>GCA in G (CA)x4 T: placed at an odd offset the read inserts CA (the input's
    bases), at an even one AC (a rotation); both are the same allele."""
    bases = "CA" if junction % 2 else "AC"
    assert _count(tmp_path, STR, "G", "GCA", A + junction, "I", bases) == (5, 5, 0)


def test_str_insertion_of_its_bases_out_of_phase_is_another_allele(tmp_path):
    """Guard: G>GCA, but the carriers insert CA at an even offset (GC CA ACACAT,
    not GCACACACAT). Inside the repeat that is another allele: partial evidence,
    never ALT."""
    assert _count(tmp_path, STR, "G", "GCA", A + 2, "I", "CA") == (5, 0, 5)


def test_an_insertion_of_other_bases_in_the_run_is_not_alt(tmp_path):
    """Guard: a G inserted mid-run is another allele: partial evidence, never ALT."""
    rd, ad, partial = _count(tmp_path, HOMOPOLYMER, "G", "GA", A + 3, "I", "G")
    assert (ad, partial) == (0, 5)


def test_the_bases_inserted_just_past_the_run_are_another_event(tmp_path):
    """G>GA, but the carriers insert an A after the T past the run: they show the
    run and both its flanks as reference, so they are REF."""
    assert _count(tmp_path, HOMOPOLYMER, "G", "GA", A + 7, "I", "A") == (10, 0, 0)


def test_the_bases_inserted_before_the_anchor_are_another_event(tmp_path):
    """G>GA, but the carriers insert the A before the G (AGAAAAAT, not GAAAAAAT):
    they show the run as reference, so they are REF."""
    assert _count(tmp_path, HOMOPOLYMER, "G", "GA", A, "I", "A") == (10, 0, 0)


def test_a_read_with_a_deletion_between_the_placements_is_not_alt(tmp_path):
    """G>GA, but the carriers delete the run's first A and insert an A three bases
    on: their bases are the reference run, not the variant written elsewhere."""

    def carrier(s):
        left = A + 1 - s
        seq = (HOMOPOLYMER[s : A + 1] + HOMOPOLYMER[A + 2 : A + 4] + "A" + HOMOPOLYMER[A + 4 :])[
            :READ
        ]
        return seq, ((0, left), (2, 1), (0, 2), (1, 1), (0, READ - left - 3))

    rd, ad, partial = _count_reads(tmp_path, HOMOPOLYMER, "G", "GA", carrier)
    assert ad == 0


def test_a_read_spliced_over_the_run_start_is_not_alt(tmp_path):
    """G>GA, but the carriers are spliced from before the anchor to mid-run and
    insert an A there: they do not show the run's start, so their insertion is not
    the variant written elsewhere."""

    def carrier(s):
        left = A - 10 - s  # N over [A-10, A+3)
        seq = (HOMOPOLYMER[s : A - 10] + "A" + HOMOPOLYMER[A + 3 :])[:READ]
        return seq, ((0, left), (3, 13), (1, 1), (0, READ - left - 1))

    rd, ad, partial = _count_reads(tmp_path, HOMOPOLYMER, "G", "GA", carrier)
    assert ad == 0


def test_an_anchor_substituting_insertion_is_not_matched_elsewhere(tmp_path):
    """A>CCC in A CC T: the carriers keep the anchor A and insert CC after the two
    C's (ACCCCT, not CCCCCT). No placement elsewhere substitutes the anchor: REF."""
    assert _count(tmp_path, ANCHOR_SUB, "A", "CCC", A + 3, "I", "CC") == (10, 0, 0)


@pytest.mark.parametrize("junction", [9, 13, 17])
def test_a_long_duplication_insertion_at_every_junction(tmp_path, junction):
    """G>G+ACGTTGCA over two copies of it: the insertion slides over 16 junctions,
    past the scan's repeat-span reach (no short motif repeats here); each carrier
    inserts the 8 bases before its junction."""
    j = A + junction
    dup = DUPLICATION[j - 8 : j]
    assert _count(tmp_path, DUPLICATION, "G", "G" + DUP, j, "I", dup) == (5, 5, 0)


# ── Deletions ─────────────────────────────────────────────────────────────────
@pytest.mark.parametrize("junction", [1, 3])
def test_homopolymer_deletion_carriers_in_the_run(tmp_path, junction):
    """Guard: GA>G, the deleted A placed anywhere in the run (counted before too)."""
    assert _count(tmp_path, HOMOPOLYMER, "GA", "G", A + junction, "D", 1) == (5, 5, 0)


@pytest.mark.parametrize("junction", [1, 2, 3, 4])
def test_str_deletion_carriers_at_every_phase(tmp_path, junction):
    """GCA>G in G (CA)x4 T: removing AC at an even offset is the same allele as
    removing CA."""
    assert _count(tmp_path, STR, "GCA", "G", A + junction, "D", 2) == (5, 5, 0)


def test_a_deletion_just_past_the_run_is_not_alt(tmp_path):
    """Guard: GA>G, but the carriers delete the T after the run: another allele."""
    rd, ad, partial = _count(tmp_path, HOMOPOLYMER, "GA", "G", A + 6, "D", 1)
    assert ad == 0


def test_a_long_duplication_deletion_placed_on_its_second_copy(tmp_path):
    """G+ACGTTGCA>G over two copies of it: deleting the second copy is the same
    allele, past the scan's repeat-span reach."""
    assert _count(tmp_path, DUPLICATION, "G" + DUP, "G", A + 9, "D", 8) == (5, 5, 0)


# ── The read's own haplotype across the discrimination window ────────────────
# A placement is the variant written elsewhere only when it is the read's one
# change across the window: another gap, insertion or splice there makes another
# haplotype (or none the read shows).
@RED
def test_an_insertion_with_a_deletion_in_the_window_is_not_alt(tmp_path):
    """G>GA, but the carriers insert an A at A+2 and delete one at A+4: their bases
    are the reference run."""
    events = [(A + 2, "I", "A"), (A + 4, "D", 1)]
    rd, ad, partial = _count_reads(tmp_path, HOMOPOLYMER, "G", "GA", _ops(HOMOPOLYMER, events))
    assert ad == 0


@RED
def test_two_insertions_in_the_run_are_another_allele(tmp_path):
    """G>GA, but the carriers insert an A at A+2 and another at A+4: a +AA allele,
    partial evidence for the +A row, never its ALT."""
    events = [(A + 2, "I", "A"), (A + 4, "I", "A")]
    assert _count_reads(tmp_path, HOMOPOLYMER, "G", "GA", _ops(HOMOPOLYMER, events)) == (5, 0, 5)


@RED
def test_an_insertion_followed_by_a_splice_in_the_run_is_not_alt(tmp_path):
    """G>GA, but the carriers insert an A at A+3 and are spliced right after it:
    they do not show the rest of the run."""
    events = [(A + 3, "I", "A"), (A + 3, "N", 57)]
    rd, ad, partial = _count_reads(tmp_path, HOMOPOLYMER, "G", "GA", _ops(HOMOPOLYMER, events))
    assert ad == 0


@RED
def test_two_deletions_in_the_run_are_another_allele(tmp_path):
    """GA>G, but the carriers delete one A at A+2 and another at A+4: a -AA
    allele, partial evidence for the -A row, never its ALT."""
    events = [(A + 2, "D", 1), (A + 4, "D", 1)]
    assert _count_reads(tmp_path, HOMOPOLYMER, "GA", "G", _ops(HOMOPOLYMER, events)) == (5, 0, 5)


@RED
def test_two_str_deletions_are_another_allele(tmp_path):
    """GCA>G in (CA)x4, but the carriers delete AC at A+2 and again at A+6: a -4
    allele (no -4 row annotated), partial evidence, never the -2 row's ALT."""
    events = [(A + 2, "D", 2), (A + 6, "D", 2)]
    assert _count_reads(tmp_path, STR, "GCA", "G", _ops(STR, events)) == (5, 0, 5)


@RED
def test_a_deletion_cancelled_by_an_insertion_is_not_alt(tmp_path):
    """GCA>G in (CA)x4, but the carriers delete AC at A+2 and insert AC at A+6:
    their bases are the reference repeat."""
    events = [(A + 2, "D", 2), (A + 6, "I", "AC")]
    rd, ad, partial = _count_reads(tmp_path, STR, "GCA", "G", _ops(STR, events))
    assert ad == 0


@RED
@pytest.mark.parametrize(
    "ref, alt_len, events",
    [
        (HOMOPOLYMER, 1, [(A - 10, "N", 12), (A + 2, "D", 1)]),
        (STR, 2, [(A - 10, "N", 12), (A + 2, "D", 2)]),
        (DUPLICATION, 8, [(A - 10, "N", 15), (A + 5, "D", 8)]),
    ],
    ids=["homopolymer", "str", "duplication"],
)
def test_a_deletion_after_a_splice_over_the_anchor_is_not_alt(tmp_path, ref, alt_len, events):
    """The carriers are spliced from before the anchor into the repeat and delete
    there: they do not show the repeat's start, so their deletion is not the
    variant written elsewhere."""
    rd, ad, partial = _count_reads(
        tmp_path, ref, ref[A : A + alt_len + 1], ref[A], _ops(ref, events)
    )
    assert ad == 0


# ── The scan reaches the variant's placements, not its whole deleted span ─────
@RED
@pytest.mark.parametrize(
    "ref, n, start", [(UNIQUE10, 10, 8), (UNIQUE60, 60, 20)], ids=["10bp", "60bp"]
)
def test_another_deletion_inside_a_deletions_span_is_not_alt(tmp_path, ref, n, start):
    """A unique deletion, but the carriers delete a different stretch of its length
    starting inside it: another haplotype, never ALT (a deletion slides over its
    shift region only up to its last placement)."""
    assert _count(tmp_path, ref, ref[A : A + n + 1], ref[A], A + start, "D", n) == (10, 0, 0)
