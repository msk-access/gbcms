"""C37 #245: a junction insert keeps its ALT; a carrier behind a substitution is not REF.

Measured on RC (ABRA2-realigned), the WES loci (BWA only) and FORTE (STAR), every read
judged by its local haplotype (CYCLE_6.6.0_PLAN.md, "C37 #245 and C38 #246"):
- the strict path credits ALT to an insert at the variant's junction with the ALT's
  bases even when a slide could absorb a flank substitution into another insert: RC
  28 reads, each the given ALT plus one molecule's substitution, no other allele
  recurring (ABRA2 put them there: a singleton gets no contig of its own);
- a readable insert of the variant's length written just outside the window behind a
  substitution, which no slide can carry back (the base it would place onto the
  reference does not fit), is a separate event, so the read counts REF though its
  fewest-mismatch placement is at the variant's junction: RC 20 reads (one recurring
  other allele), WES 2, FORTE 1.

Contract (operator, 2026-10-07):
- the strict path keeps its junction ALT: every surveyed tool credits it;
- a readable same-length insert whose fewest-mismatch placements, every junction
  scored by its readable mismatches (not only those a slide reaches base by base),
  include one inside the discrimination window is another allele (neither, with
  partial evidence): never REF, and never ALT, since a per-read rule cannot tell a
  singleton error from a recurring other allele (RJ-22).
The census trusts the aligned flank, so it cannot judge the second shape.
"""

import pytest
from test_shifted_insert_contract import BACKENDS, CONTIG, INS8, P, R, U, _count, read, ref_read


@pytest.mark.parametrize("backend", BACKENDS)
def test_the_strict_path_keeps_a_junction_alt_with_a_substitution_in_the_tract(tmp_path, backend):
    """G>GA before A10, and G>GCA before (CA)x6: the ALT written at the junction with a
    substitution inside the run (A>C) or tract (C>G). Sliding the insert along the
    repeat would take the substitution in as another insert with no mismatch, but the
    read is the given ALT plus one substitution, and it stays ALT, as a REF read with a
    substitution in the window stays REF. No census: it needs the window to be one
    allele exactly, so it reads both as neither (its first documented limit)."""
    reads = [read(f"a{i}", i, P, "A", P + 1, sub=(P + 4, "C")) for i in range(4)]
    reads += [ref_read(f"r{i}", i, P) for i in range(4)]
    got = _count(tmp_path, f"ka{backend}", reads, P, "G", "GA", backend, with_census=False)
    assert got == (4, 4, 0)
    reads = [read(f"c{i}", i, R, "CA", R + 1, sub=(R + 3, "G")) for i in range(4)]
    reads += [ref_read(f"r{i}", i, R) for i in range(4)]
    got = _count(tmp_path, f"kc{backend}", reads, R, "G", "GCA", backend, with_census=False)
    assert got == (4, 4, 0)


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
    got = _count(tmp_path, f"da{backend}", reads, P, "G", "GA", backend, with_census=False)
    assert got == (4, 0, 4)
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS8
    reads = [read(f"u{i}", i, a, "T" + INS8[:-1], a, sub=(a, INS8[-1])) for i in range(4)]
    reads += [ref_read(f"r{i}", i, a) for i in range(4)]
    got = _count(tmp_path, f"du{backend}", reads, a, ref, alt, backend, with_census=False)
    assert got == (4, 0, 4)
    reads = [read(f"c{i}", i, R, "TC", R, sub=(R, "A")) for i in range(4)]
    reads += [ref_read(f"r{i}", i, R) for i in range(4)]
    got = _count(tmp_path, f"dc{backend}", reads, R, "G", "GCA", backend, with_census=False)
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
    got = _count(tmp_path, f"so{backend}", reads, a, ref, alt, backend, with_census=False)
    assert got == (8, 0, 0)
