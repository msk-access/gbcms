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

RED = pytest.mark.xfail(strict=True, reason="placement equivalence gaps found in review")

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


def test_an_insertion_of_other_bases_in_the_run_is_not_alt(tmp_path):
    """Guard: a G inserted mid-run is another allele: partial evidence, never ALT."""
    rd, ad, partial = _count(tmp_path, HOMOPOLYMER, "G", "GA", A + 3, "I", "G")
    assert (ad, partial) == (0, 5)


@RED
def test_the_bases_inserted_just_past_the_run_are_another_event(tmp_path):
    """G>GA, but the carriers insert an A after the T past the run: they show the
    run and both its flanks as reference, so they are REF."""
    assert _count(tmp_path, HOMOPOLYMER, "G", "GA", A + 7, "I", "A") == (10, 0, 0)


@RED
def test_the_bases_inserted_before_the_anchor_are_another_event(tmp_path):
    """G>GA, but the carriers insert the A before the G (AGAAAAAT, not GAAAAAAT):
    they show the run as reference, so they are REF."""
    assert _count(tmp_path, HOMOPOLYMER, "G", "GA", A, "I", "A") == (10, 0, 0)


@RED
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


@RED
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


@RED
def test_an_anchor_substituting_insertion_is_not_matched_elsewhere(tmp_path):
    """A>CCC in A CC T: the carriers keep the anchor A and insert CC after the two
    C's (ACCCCT, not CCCCCT). No placement elsewhere substitutes the anchor: REF."""
    assert _count(tmp_path, ANCHOR_SUB, "A", "CCC", A + 3, "I", "CC") == (10, 0, 0)


@RED
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


@RED
def test_a_long_duplication_deletion_placed_on_its_second_copy(tmp_path):
    """G+ACGTTGCA>G over two copies of it: deleting the second copy is the same
    allele, past the scan's repeat-span reach."""
    assert _count(tmp_path, DUPLICATION, "G" + DUP, "G", A + 9, "D", 8) == (5, 5, 0)
