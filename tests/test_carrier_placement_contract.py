"""C37 #245: a junction insert keeps its ALT; a carrier behind a substitution is not REF.

Measured on RC (ABRA2-realigned), the WES loci (BWA only) and FORTE (STAR), every read
judged by its local haplotype (measured in PR #249):
- the strict path credits ALT to an insert at the variant's junction with the ALT's
  bases even when a slide could absorb a flank substitution into another insert: RC
  28 reads, each the given ALT plus one molecule's substitution, no other allele
  recurring (ABRA2 put them there: a singleton gets no contig of its own);
- a readable insert of the variant's length written just outside the window behind a
  substitution, which no slide can carry back (the base it would place onto the
  reference does not fit), is a separate event, so the read counts REF though its
  fewest-mismatch placement is at the variant's junction: RC 19 reads at 2 rows (one
  recurring other allele), WES 2, FORTE 2.

Contract (operator, 2026-10-07):
- the strict path keeps its junction ALT: every surveyed tool credits it;
- a readable same-length insert whose fewest-mismatch placements, every junction it
  can reach scored by its readable mismatches (not only those a slide reaches base by
  base),
  include one inside the discrimination window is another allele (neither, with
  partial evidence): never REF, and never ALT, since a per-read rule cannot tell a
  singleton error from a recurring other allele (RJ-22).
The census trusts the aligned flank, so it cannot judge the second shape.
"""

import random

import pytest
from helpers import count_checked, make_read, write_contig
from test_shifted_insert_contract import BACKENDS, CONTIG, INS8, L, P, Q, R, U, read, ref_read

from gbcms import _rs


def _count(tmp_path, name, reads, pos, ref, alt, backend, contig=CONTIG):
    """(rd, ad, partial_alt) for one insertion row, with the counting invariants every
    counting test asserts (AGENTS.md). No census: see the module docstring."""
    fa, bam = write_contig(tmp_path, contig, reads, name)
    (pv,) = _rs.prepare_variants([_rs.Variant("1", pos, ref, alt, "INS")], fa, 5, False, 1, True)
    assert (pv.gbcms_status, pv.variant.pos, pv.variant.alt_allele) == ("PASS", pos, alt)
    (c,) = count_checked(bam, [pv.variant], alignment_backend=backend)
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    return c.rd, c.ad, c.partial_alt


@pytest.mark.parametrize("backend", BACKENDS)
def test_the_strict_path_keeps_a_junction_alt_with_a_substitution_in_the_tract(tmp_path, backend):
    """G>GA before A10, and G>GCA before (CA)x6: the ALT written at the junction with a
    substitution inside the run (A>C) or tract (C>G). Sliding the insert along the
    repeat would take the substitution in as another insert with no mismatch, but the
    read is the given ALT plus one substitution, and it stays ALT, as a REF read with a
    substitution in the window stays REF. (The census needs the window to be one
    allele exactly, so it reads both as neither: its first documented limit.)"""
    reads = [read(f"a{i}", i, P, "A", P + 1, sub=(P + 4, "C")) for i in range(4)]
    reads += [ref_read(f"r{i}", i, P) for i in range(4)]
    assert _count(tmp_path, f"ka{backend}", reads, P, "G", "GA", backend) == (4, 4, 0)
    reads = [read(f"c{i}", i, R, "CA", R + 1, sub=(R + 3, "G")) for i in range(4)]
    reads += [ref_read(f"r{i}", i, R) for i in range(4)]
    assert _count(tmp_path, f"kc{backend}", reads, R, "G", "GCA", backend) == (4, 4, 0)


@pytest.mark.parametrize("backend", BACKENDS)
def test_a_carrier_behind_an_anchor_substitution_is_another_allele(tmp_path, backend):
    """The ALT haplotype with its anchor substituted, the aligner writing the insert one
    junction left of the anchor (a tie: one mismatch either way, at the anchor). The
    base it would place onto the anchor does not fit, so no slide carries it back:
    G>GA before A10 (an A inserted, the anchor read as A); the non-repeat +8 (the
    anchor read as T); G>GCA before (CA)x6 (the anchor read as T). At its fewest
    mismatches the insert can sit at the variant's junction: another allele, not REF;
    and not ALT, though that placement shows the ALT's bases."""
    reads = [read(f"a{i}", i, P, "A", P, sub=(P, "A")) for i in range(4)]
    reads += [ref_read(f"r{i}", i, P) for i in range(4)]
    got = _count(tmp_path, f"da{backend}", reads, P, "G", "GA", backend)
    assert got == (4, 0, 4)
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS8
    reads = [read(f"u{i}", i, a, "T" + INS8[:-1], a, sub=(a, INS8[-1])) for i in range(4)]
    reads += [ref_read(f"r{i}", i, a) for i in range(4)]
    got = _count(tmp_path, f"du{backend}", reads, a, ref, alt, backend)
    assert got == (4, 0, 4)
    reads = [read(f"c{i}", i, R, "TC", R, sub=(R, "A")) for i in range(4)]
    reads += [ref_read(f"r{i}", i, R) for i in range(4)]
    got = _count(tmp_path, f"dc{backend}", reads, R, "G", "GCA", backend)
    assert got == (4, 0, 4)


@pytest.mark.parametrize("backend", BACKENDS)
def test_an_insert_whose_window_placement_costs_more_stays_separate(tmp_path, backend):
    """The non-repeat +8: an insert of other bases one junction left of the anchor,
    exact there. Placed at the variant's junction it costs a mismatch the written
    placement does not have, so its fewest-mismatch placement stays outside the
    window: a separate event, REF (RJ-8), as C36 decided."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS8
    reads = [read(f"o{i}", i, a, "CTTAGCCT", a) for i in range(4)]
    reads += [ref_read(f"r{i}", i, a) for i in range(4)]
    got = _count(tmp_path, f"so{backend}", reads, a, ref, alt, backend)
    assert got == (8, 0, 0)


# ── From the review: the scored slide's other branches ─────────────────────────


def _gct_contig():
    """...A G T T T | G C T G C T G C T G C (the anchor C at 300) | T G T G T A T T A A..."""
    rng = random.Random(245)
    c = [rng.choice("ACGT") for _ in range(1000)]
    left, right, a = "AGTTT" + "GCTGCTGCTGC", "TGTGTATTAA", 300
    c[a - len(left) + 1 : a + 1] = left
    c[a + 1 : a + 1 + len(right)] = right
    return "".join(c), a


@pytest.mark.parametrize("off", [7, 10])
@pytest.mark.parametrize("backend", BACKENDS)
def test_a_carrier_written_past_the_scan_window_is_another_allele(tmp_path, backend, off):
    """C>CTGT after (GCT)n, as on FORTE STAR: the ALT with its anchor read as T, written
    as +TCT before a G of the tract, 7 or 10 junctions left (past the scan pad, the tract
    not the variant's own repeat). No slide reaches the window base by base; scored, its
    fewest-mismatch placements include the variant's junction, so the scan-window gate
    admits it and it is another allele, not REF."""
    contig, a = _gct_contig()
    j, reads = a - off, []
    for i in range(4):
        s = a - 40 - i
        ref_part = list(contig[s : s + L])
        ref_part[a - s] = "T"
        left = j - s
        seq = "".join(ref_part[:left]) + "TCT" + "".join(ref_part[left : L - 3])
        cig = ((0, left), (1, 3), (0, L - 3 - left))
        reads.append(make_read(f"c{i}", seq, s, cig, flag=16 * (i % 2)))
    reads += [
        make_read(f"r{i}", contig[a - 40 - i : a + 60 - i], a - 40 - i, ((0, L),)) for i in range(4)
    ]
    got = _count(tmp_path, f"g{backend}{off}", reads, a, "C", "CTGT", backend, contig=contig)
    assert got == (4, 0, 4)


@pytest.mark.parametrize("backend", BACKENDS)
def test_the_left_slide_scores_past_a_misfit(tmp_path, backend):
    """A>AT before a C: a +G written one junction right of the variant's junction, the
    C between read as T. Scored left, the G placed onto that base costs one and the T
    taken in resolves one: the insert can sit at the variant's junction at the written
    cost, so the read (REF with C>TG there, or an insert beside a substitution) is
    another allele, not REF."""
    assert CONTIG[Q + 1] == "C"
    reads = [read(f"q{i}", i, Q, "G", Q + 2, sub=(Q + 1, "T")) for i in range(4)]
    reads += [ref_read(f"r{i}", i, Q) for i in range(4)]
    assert _count(tmp_path, f"l{backend}", reads, Q, "A", "AT", backend) == (4, 0, 4)
