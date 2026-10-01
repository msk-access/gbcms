"""A pure indel's carrier holds the ALT across the variant's window (#188, #191, #192, #121).

The pure-indel checks credit a read from its CIGAR. Measured on the RC set, WES and
FORTE (cluster 1 of the 6.6.0 triage), these gaps remain; each is small in real
data, and each breaks "judge a read's bases, not its placement":

- C20 (#188): a read that is not informative (spans neither of C10's windows) is
  credited ALT from its gap even when its own bases fit both alleles. A read whose
  bases discriminate stays ALT (a long insertion's carriers often cannot span the
  windows yet hold the inserted bases).
- C20 (#188), strict path: an indel at the variant's own junction counts ALT even
  when the read carries another insertion or deletion across the window.
- C22 (#191): a same-length deletion of 5bp or more that gives another haplotype
  goes to Phase 3, which often calls it ALT; it is a distinct allele.
- C23 (#192): a distinct allele keeps REF when the event's repeat has a motif over
  6bp, because "in a repeat" is decided by repeat_span; the shift region decides.
- C8 (#121): a one-base-REF ALT that changes its anchor base (A>CCC) is judged as an
  insertion of its tail; reads that keep the anchor count ALT.

Committed red (xfail-strict) before the implementation; flipped green with it.
"""

import pytest
from test_shifted_indel_carriers_contract import (
    ANCHOR_SUB,
    DUP,
    DUPLICATION,
    HOMOPOLYMER,
    READ,
    A,
    _contig,
    _count,
    _count_reads,
    _ops,
)

RED = pytest.mark.xfail(strict=True, reason="cluster 1: pure-indel read judgment")

SPAN8 = "ATCGGATA"  # a unique 8bp stretch between G and C
UNIQUE8 = _contig("G" + SPAN8 + "C")
INSERT10 = "ACGTTGCATC"  # not low-complexity; differs from the flank at once


# ── C20: an uninformative carrier is not ALT unless its bases discriminate ──────
@RED
def test_a_carrier_ending_inside_the_run_is_uninformative(tmp_path):
    """GA>G in G AAAAA T: the carriers delete an A after the G and end two bases
    later, inside the run. Their bases (G A A) fit both alleles: depth, neither."""

    def carrier(s):
        left = A + 1 - s
        return HOMOPOLYMER[s : A + 1] + HOMOPOLYMER[A + 2 : A + 4], ((0, left), (2, 1), (0, 2))

    c = _count_reads(tmp_path, HOMOPOLYMER, "GA", "G", carrier, full=True)
    assert (c.rd, c.ad, c.partial_alt) == (5, 0, 0)
    assert c.dp == 10, "uninformative carriers still count toward depth"


def test_a_carrier_whose_inserted_bases_discriminate_stays_alt(tmp_path):
    """Guard: G>G+ACGTTGCATC in unique sequence; the carriers end inside the insert
    (a truncated insertion) and span neither window, but their bases (G then A C...
    where the reference continues otherwise) are the ALT's."""
    ref = _contig("G")  # a unique G anchor in random C/G/T flank

    def carrier(s):
        left = A + 1 - s
        return ref[s : A + 1] + INSERT10[:6], ((0, left), (1, 6))

    assert _count_reads(tmp_path, ref, "G", "G" + INSERT10, carrier) == (5, 5, 0)


# ── C20, strict path: the indel must be the read's only change in the window ────
@RED
def test_an_insertion_at_the_junction_with_a_deletion_in_the_window_is_not_alt(tmp_path):
    """G>GA, the carriers insert the A at the junction and delete one 3 bases on:
    their bases are the reference run."""
    events = [(A + 1, "I", "A"), (A + 4, "D", 1)]
    rd, ad, partial = _count_reads(tmp_path, HOMOPOLYMER, "G", "GA", _ops(HOMOPOLYMER, events))
    assert ad == 0


@RED
def test_a_split_plus_two_read_is_another_allele_at_the_junction(tmp_path):
    """G>GA, the carriers insert an A at the junction and another 3 bases on: a +AA
    allele, partial evidence for the +A row."""
    events = [(A + 1, "I", "A"), (A + 4, "I", "A")]
    assert _count_reads(tmp_path, HOMOPOLYMER, "G", "GA", _ops(HOMOPOLYMER, events)) == (5, 0, 5)


@RED
def test_a_split_minus_two_read_is_another_allele_at_the_junction(tmp_path):
    """GA>G, the carriers delete an A at the junction and another 3 bases on: a -AA
    allele, partial evidence for the -A row."""
    events = [(A + 1, "D", 1), (A + 4, "D", 1)]
    assert _count_reads(tmp_path, HOMOPOLYMER, "GA", "G", _ops(HOMOPOLYMER, events)) == (5, 0, 5)


# ── C22: a non-equivalent same-length deletion is a distinct allele ─────────────
@RED
def test_a_non_equivalent_same_length_deletion_is_a_distinct_allele(tmp_path):
    """G+ATCGGATA>G in unique sequence, the carriers delete 8 other bases starting
    2 in: another haplotype, never ALT; in unique sequence the anchor stays REF and
    the read is partial evidence."""
    assert _count(tmp_path, UNIQUE8, "G" + SPAN8, "G", A + 3, "D", 8) == (10, 0, 5)


# ── C23: "in a repeat" is decided by the shift region ──────────────────────────
@RED
def test_a_distinct_allele_in_a_long_period_duplication_is_not_ref(tmp_path):
    """G>G+ACGTTGCA over two copies of it (no motif of 6bp or fewer repeats): the
    carriers insert the 8 bases unrotated 13 bases in, another haplotype inside the
    window. The event slides, so it is a repeat: neither + partial, not REF."""
    assert _count(tmp_path, DUPLICATION, "G", "G" + DUP, A + 13, "I", DUP) == (5, 0, 5)


@RED
def test_a_wrong_length_deletion_in_a_long_period_duplication_is_not_ref(tmp_path):
    """G+ACGTTGCA>G over two copies of it: the carriers delete 5 bases starting 5 in
    (a wrong-length deletion inside the region): neither + partial, not REF."""
    assert _count(tmp_path, DUPLICATION, "G" + DUP, "G", A + 5, "D", 5) == (5, 0, 5)


# ── C8: an anchor-changing one-base-REF ALT is judged by its whole allele ──────
@RED
def test_a_read_keeping_the_anchor_is_not_an_anchor_changing_alt(tmp_path):
    """A>CCC in A CC T: the carriers keep the anchor A and insert CC right after it
    (ACCCCT, not CCCCCT): not the given allele."""
    rd, ad, partial = _count(tmp_path, ANCHOR_SUB, "A", "CCC", A + 1, "I", "CC")
    assert ad == 0


def test_true_anchor_changing_carriers_count_alt(tmp_path):
    """Guard: A>CCC, the carriers read C at the anchor and insert CC after it (the
    aligner writes a mismatch and an insertion): the given allele."""

    def carrier(s):
        left = A + 1 - s
        seq = (ANCHOR_SUB[s:A] + "C" + "CC" + ANCHOR_SUB[A + 1 :])[:READ]
        return seq, ((0, left), (1, 2), (0, READ - left - 2))

    assert _count_reads(tmp_path, ANCHOR_SUB, "A", "CCC", carrier) == (5, 5, 0)
