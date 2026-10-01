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

import random

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

SPAN8 = "ATCGGATA"  # a unique 8bp stretch between G and C
UNIQUE8 = _contig("G" + SPAN8 + "C")
INSERT10 = "ACGTTGCATC"  # not low-complexity; differs from the flank at once
_RNG66, _RNG57 = random.Random(66), random.Random(57)
COPY66 = "A" + "".join(_RNG66.choice("ACGT") for _ in range(64)) + "C"
ITD66 = _contig("G" + COPY66 + "T")  # the ALT duplicates this 66bp copy in tandem
SPAN57 = "A" + "".join(_RNG57.choice("ACGT") for _ in range(55)) + "A"
UNIQUE57 = _contig("G" + SPAN57 + "C")
HOMOPOLYMER6 = _contig("GAAAAAAT")  # G, A x6 at A+1..A+6, T
_RNG30 = random.Random(30)
UNIT30 = "".join(_RNG30.choice("ACGT") for _ in range(30))
while UNIT30[0] == "G" or UNIT30[-1] == "G":
    UNIT30 = "".join(_RNG30.choice("ACGT") for _ in range(30))
TANDEM30 = _contig("G" + UNIT30 * 5 + "T")  # G, five copies of a 30bp unit, T
GCCCT = _contig("GCCCT")  # G at A, C C C T after it
_RNG60 = random.Random(60)
INSERT60 = INSERT10 + "".join(_RNG60.choice("ACGT") for _ in range(49)) + "C"  # not G: stays put
UNIQUE_GT = _contig("GT")  # a unique G anchor, then T
_CG7 = _contig("C" + "G" * 7 + "C")
CCG7 = _CG7[: A - 1] + "C" + _CG7[A:]  # C C GGGGGGG C, the second C (the anchor) at A


# ── C20: an uninformative carrier is not ALT unless its bases discriminate ──────
def test_a_carrier_ending_inside_the_run_is_uninformative(tmp_path):
    """GA>G in G AAAAA T: the carriers delete an A after the G and end two bases
    later, inside the run. Their bases (G A A) fit both alleles: depth, neither."""

    def carrier(s):
        left = A + 1 - s
        return HOMOPOLYMER[s : A + 1] + HOMOPOLYMER[A + 2 : A + 4], ((0, left), (2, 1), (0, 2))

    c = _count_reads(tmp_path, HOMOPOLYMER, "GA", "G", carrier, full=True)
    assert (c.rd, c.ad, c.partial_alt) == (5, 0, 0)
    assert c.dp == 10, "uninformative carriers still count toward depth"


def test_a_carrier_ending_on_the_discriminating_base_stays_alt(tmp_path):
    """G>GA in G AAAAA T: the carriers insert an A at the junction and end on the
    last A of the run, where REF has the T. They read the deciding base (an A, not
    the T) unmasked, and their CIGAR already says ALT: ALT. (C10's extra margin
    base protects CIGAR-only REF calls from a hidden terminal mismatch; here the
    base itself is read.)"""

    def carrier(s):
        left = A + 1 - s
        return HOMOPOLYMER[s : A + 1] + "A" + HOMOPOLYMER[A + 1 : A + 6], (
            (0, left),
            (1, 1),
            (0, 5),
        )

    assert _count_reads(tmp_path, HOMOPOLYMER, "G", "GA", carrier) == (5, 5, 0)


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
def test_an_insertion_at_the_junction_with_a_deletion_in_the_window_is_not_alt(tmp_path):
    """G>GA, the carriers insert the A at the junction and delete one 3 bases on:
    their bases are the reference run."""
    events = [(A + 1, "I", "A"), (A + 4, "D", 1)]
    rd, ad, partial = _count_reads(tmp_path, HOMOPOLYMER, "G", "GA", _ops(HOMOPOLYMER, events))
    assert ad == 0


def test_a_split_plus_two_read_is_another_allele_at_the_junction(tmp_path):
    """G>GA, the carriers insert an A at the junction and another 3 bases on: a +AA
    allele, partial evidence for the +A row."""
    events = [(A + 1, "I", "A"), (A + 4, "I", "A")]
    assert _count_reads(tmp_path, HOMOPOLYMER, "G", "GA", _ops(HOMOPOLYMER, events)) == (5, 0, 5)


def test_a_split_minus_two_read_is_another_allele_at_the_junction(tmp_path):
    """GA>G, the carriers delete an A at the junction and another 3 bases on: a -AA
    allele, partial evidence for the -A row."""
    events = [(A + 1, "D", 1), (A + 4, "D", 1)]
    assert _count_reads(tmp_path, HOMOPOLYMER, "GA", "G", _ops(HOMOPOLYMER, events)) == (5, 0, 5)


# ── C22: a non-equivalent same-length deletion is a distinct allele ─────────────
def test_a_non_equivalent_same_length_deletion_is_a_distinct_allele(tmp_path):
    """G+ATCGGATA>G in unique sequence, the carriers delete 8 other bases starting
    2 in: another haplotype inside the discrimination window, so neither REF nor
    ALT: partial evidence."""
    assert _count(tmp_path, UNIQUE8, "G" + SPAN8, "G", A + 3, "D", 8) == (5, 0, 5)


@pytest.mark.parametrize("shift", [1, 2, 3])
def test_a_misplaced_deletion_whose_bases_are_the_alt_counts_alt(tmp_path, shift):
    """G+ATCGGATA>G, the carriers' bases are exactly the ALT haplotype but the
    aligner wrote the D(8) `shift` bases to the right with compensating mismatches:
    their bases carry the allele, so they are ALT."""

    def carrier(s):
        left = A + 1 + shift - s
        seq = (UNIQUE8[s : A + 1] + UNIQUE8[A + 9 :])[:READ]
        return seq, ((0, left), (2, 8), (0, READ - left))

    assert _count_reads(tmp_path, UNIQUE8, "G" + SPAN8, "G", carrier) == (5, 5, 0)


# ── C23: "in a repeat" is decided by the shift region ──────────────────────────
def test_a_distinct_allele_in_a_long_period_duplication_is_not_ref(tmp_path):
    """G>G+ACGTTGCA over two copies of it (no motif of 6bp or fewer repeats): the
    carriers insert the 8 bases unrotated 13 bases in, another haplotype inside the
    window. The event slides, so it is a repeat: neither + partial, not REF."""
    assert _count(tmp_path, DUPLICATION, "G", "G" + DUP, A + 13, "I", DUP) == (5, 0, 5)


def test_a_wrong_length_deletion_in_a_long_period_duplication_is_not_ref(tmp_path):
    """G+ACGTTGCA>G over two copies of it: the carriers delete 5 bases starting 5 in
    (a wrong-length deletion inside the region): neither + partial, not REF."""
    assert _count(tmp_path, DUPLICATION, "G" + DUP, "G", A + 5, "D", 5) == (5, 0, 5)


# ── C8: an anchor-changing one-base-REF ALT is judged by its whole allele ──────
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


# ── Review findings ──────────────────────────────────────────────────────────
def test_a_long_duplication_carrier_ending_with_the_insert_is_uninformative(tmp_path):
    """G>G+COPY66, a tandem duplication of the 66 bases after the G: the carriers
    write the whole insert at the junction and end there. Their bases (G then one
    copy) are also the reference, so they fit both alleles: depth only. The rule
    must reach long events (it judged nothing from about 58bp up)."""

    def carrier(s):
        left = A + 1 - s
        return ITD66[s : A + 1] + COPY66, ((0, left), (1, 66))

    c = _count_reads(tmp_path, ITD66, "G", "G" + COPY66, carrier, full=True)
    assert c.ad == 0
    assert c.dp == 10


def test_a_masked_first_inserted_base_does_not_withdraw_a_carrier(tmp_path):
    """G>G+ACGTTGCATC in unique sequence, carriers ending inside the insert whose
    first inserted base is below --min-baseq: the next nine discriminate, so ALT."""
    ref = _contig("G")

    def carrier(s):
        left = A + 1 - s
        quals = [30] * left + [5] + [30] * 9
        return ref[s : A + 1] + INSERT10, ((0, left), (1, 10)), quals

    assert _count_reads(tmp_path, ref, "G", "G" + INSERT10, carrier) == (5, 5, 0)


def test_a_large_deletion_with_a_small_insertion_at_its_far_junction_stays_alt(tmp_path):
    """A unique 57bp deletion, the carriers write D(57) then I(TT): within the
    large-deletion band's tolerance (≤3 changed bases), as D(56)+I(TT) already is."""
    events = [(A + 1, "D", 57), (A + 58, "I", "TT")]
    ref = UNIQUE57
    assert _count_reads(tmp_path, ref, ref[A : A + 58], ref[A], _ops(ref, events)) == (5, 5, 0)


def test_soft_clipped_bases_behind_a_hard_clip_are_not_read(tmp_path):
    """GA>G in G AAAAA T, carriers ending two bases into the run, the rest of the
    read soft-clipped and then hard-clipped: the clipped bases are not read, so the
    carriers fit both alleles and are uninformative."""

    def carrier(s):
        left = A + 1 - s
        seq = HOMOPOLYMER[s : A + 1] + HOMOPOLYMER[A + 2 : A + 4] + HOMOPOLYMER[A + 4 : A + 10]
        return seq, ((0, left), (2, 1), (0, 2), (4, 6), (5, 3))

    c = _count_reads(tmp_path, HOMOPOLYMER, "GA", "G", carrier, full=True)
    assert (c.rd, c.ad, c.partial_alt) == (5, 0, 0)


def test_a_misplaced_deletion_ambiguous_at_a_masked_base_is_not_alt(tmp_path):
    """CG>C deletes a G of the run after C C. The carriers' CIGAR deletes the first
    C instead, and the base where "deleted a C" (C GGGGGGG) and "deleted a G"
    (CC GGGGGG) differ is an N: their bases cannot tell the two alleles apart, so
    not ALT (a real RC case)."""

    def carrier(s):
        left = A - 1 - s
        hap = CCG7[s : A - 1] + CCG7[A:]  # the first C deleted
        seq = (hap[: left + 1] + "N" + hap[left + 2 :])[:READ]
        return seq, ((0, left), (2, 1), (0, READ - left))

    rd, ad, partial = _count_reads(tmp_path, CCG7, "CG", "C", carrier)
    assert ad == 0


def test_a_misplaced_deletion_whose_read_base_settles_it_is_alt(tmp_path):
    """Guard: the same carriers with a C read (unmasked) at that base spell the ALT
    (CC GGGGGG): their bases carry the allele wherever the CIGAR put the gap."""

    def carrier(s):
        left = A - 1 - s
        hap = CCG7[s : A - 1] + CCG7[A:]
        seq = (hap[: left + 1] + "C" + hap[left + 2 :])[:READ]
        return seq, ((0, left), (2, 1), (0, READ - left))

    assert _count_reads(tmp_path, CCG7, "CG", "C", carrier) == (5, 5, 0)


# ── Second review findings ────────────────────────────────────────────────────
@pytest.mark.parametrize(
    "ref, ref_allele, kept, reads_after",
    [
        (HOMOPOLYMER6, "GAA", 2, 2),  # -AA in A6: G AA D(2) AA, the read ends in the run
        (DUPLICATION, "G" + DUP, 0, 7),  # one copy of a tandem 8bp unit, ending 1 short
    ],
    ids=["homopolymer", "duplication"],
)
def test_a_deletion_carrier_ending_inside_the_region_is_uninformative(
    tmp_path, ref, ref_allele, kept, reads_after
):
    """The carriers delete the event and end before the region does. Their reference
    extent counts their own gap, so it spans C10's reference windows; their bases
    (the anchor, then fewer bases of the repeat than either allele holds) fit both
    alleles: depth only."""
    n = len(ref_allele) - 1

    def carrier(s):
        left = A + 1 + kept - s
        seq = ref[s : A + 1 + kept] + ref[A + 1 + kept + n : A + 1 + kept + n + reads_after]
        return seq, ((0, left), (2, n), (0, reads_after))

    c = _count_reads(tmp_path, ref, ref_allele, "G", carrier, full=True)
    assert (c.ad, c.dp) == (0, 10)


def test_a_deletion_carrier_reading_past_the_run_stays_alt(tmp_path):
    """Guard: -AA in A6, the carriers read G AAAA and then the T, where REF has a
    fifth A: their bases discriminate, so ALT."""

    def carrier(s):
        left = A + 3 - s
        return HOMOPOLYMER6[s : A + 3] + HOMOPOLYMER6[A + 5 : A + 8], ((0, left), (2, 2), (0, 3))

    assert _count_reads(tmp_path, HOMOPOLYMER6, "GAA", "G", carrier)[1] == 5


@pytest.mark.parametrize(
    "ref_allele, events, expected",
    [
        # The 57bp deletion written after a 2bp insertion: within the large-deletion
        # band, as the same pair written D-first is.
        (UNIQUE57[A : A + 58], [(A + 1, "I", "TT"), (A + 1, "D", 57)], (5, 5, 0)),
        # An 8bp deletion written after a 1bp insertion: another allele.
        ("G" + SPAN8, [(A + 1, "I", "T"), (A + 1, "D", 8)], (5, 0, 5)),
    ],
    ids=["57bp-band", "8bp"],
)
def test_a_deletion_written_after_an_insertion_is_judged(tmp_path, ref_allele, events, expected):
    """The aligner may write an I/D pair at the junction in either order; the
    deletion after the insertion must be seen, not read as reference."""
    ref = UNIQUE57 if len(ref_allele) == 58 else UNIQUE8
    assert _count_reads(tmp_path, ref, ref_allele, "G", _ops(ref, events)) == expected


def test_an_insertion_written_after_a_deletion_is_judged(tmp_path):
    """G>G+ACGTTGCATC, the carriers delete the base after the anchor and then insert
    the ten bases (M D(1) I(10)): another allele, as M I(10) D(1) is, not REF."""
    ref = _contig("GT")
    events = [(A + 1, "D", 1), (A + 2, "I", INSERT10)]
    assert _count_reads(tmp_path, ref, "G", "G" + INSERT10, _ops(ref, events)) == (5, 0, 5)


def test_a_sliding_large_deletion_with_another_indel_far_in_its_window_is_not_alt(tmp_path):
    """Two copies of a 30bp unit deleted from five: the carriers write the exact
    D(60) at the anchor and insert 10 bases 80 bases on, inside the discrimination
    window (the region is 150bp) but past the band's old reach. The read nets -50:
    not this allele."""
    events = [(A + 1, "D", 60), (A + 80, "I", "ACGTACGTAC")]
    ref_allele = "G" + UNIT30 * 2
    rd, ad, partial = _count_reads(tmp_path, TANDEM30, ref_allele, "G", _ops(TANDEM30, events))
    assert ad == 0


def test_a_base_past_the_reference_stretch_decides_nothing(tmp_path):
    """G>G+ACGTTGCATC before C C C T, carriers ending four bases into the insert with
    the first three masked: they read G ? ? ? T, and the reference also has a T
    there, so they fit both alleles: depth only."""

    def carrier(s):
        left = A + 1 - s
        quals = [30] * left + [5, 5, 5, 30]
        return GCCCT[s : A + 1] + INSERT10[:4], ((0, left), (1, 4)), quals

    c = _count_reads(tmp_path, GCCCT, "G", "G" + INSERT10, carrier, full=True)
    assert (c.ad, c.dp) == (0, 10)


# ── Third review findings ─────────────────────────────────────────────────────
def test_a_deletion_then_a_same_length_insertion_of_other_bases_is_not_alt(tmp_path):
    """G>G+ACGTA, the carriers delete the base after the anchor and insert ACGTC
    (M D(1) I(5) M): the read changes length by 4, never the ALT's 5. Another
    allele, as M I(5) D(1) is."""
    events = [(A + 1, "D", 1), (A + 2, "I", "ACGTC")]
    ref = UNIQUE_GT
    assert _count_reads(tmp_path, ref, "G", "GACGTA", _ops(ref, events)) == (5, 0, 5)


def test_a_long_insertion_carrier_masked_past_the_first_difference_stays_alt(tmp_path):
    """G>G+INSERT60 in unique sequence, carriers ending ten bases into the insert
    with the first three masked: the next seven discriminate, as they do for a short
    insertion. The REF stretch must reach as far as the ALT's for a 60bp insert."""

    def carrier(s):
        left = A + 1 - s
        quals = [30] * left + [5, 5, 5] + [30] * 7
        return UNIQUE_GT[s : A + 1] + INSERT60[:10], ((0, left), (1, 10)), quals

    assert _count_reads(tmp_path, UNIQUE_GT, "G", "G" + INSERT60, carrier) == (5, 5, 0)


def test_a_deleted_anchor_reinserted_with_the_insert_is_alt(tmp_path):
    """G>G+ACGTTGCATC, the carriers delete the anchor G and insert G+ACGTTGCATC
    after it (M D(1) I(11) M): their bases are exactly the ALT, so ALT."""
    events = [(A, "D", 1), (A + 1, "I", "G" + INSERT10)]
    ref = UNIQUE_GT
    assert _count_reads(tmp_path, ref, "G", "G" + INSERT10, _ops(ref, events)) == (5, 5, 0)


def test_a_deleted_anchor_replaced_before_the_insert_is_not_alt(tmp_path):
    """Guard: the carriers delete the anchor G and insert T+ACGTTGCATC (an anchor
    substitution with the insert): not the given allele."""
    events = [(A, "D", 1), (A + 1, "I", "T" + INSERT10)]
    ref = UNIQUE_GT
    assert _count_reads(tmp_path, ref, "G", "G" + INSERT10, _ops(ref, events))[1] == 0


# ── Reads that start on the flank base ─────────────────────────────────────────
@pytest.mark.xfail(strict=True, reason="H3 decision 5: red until fixed")
def test_a_carrier_starting_on_the_flank_whose_bases_discriminate_is_alt(tmp_path):
    """G>GA in G AAAAA T: carriers whose first base is the G (the run's left flank,
    here the anchor), then the inserted A and the run's five A's: six A's after the
    G where REF has five and then the T, so their bases carry the ALT. The ALT
    side's windows already start at the flank; reading the bases starts there too
    (on real data a 66bp duplication's carriers starting on the anchor, holding the
    whole insert, were withdrawn)."""
    from helpers import count_bam_checked, make_read, write_contig

    from gbcms import _rs as gbcms_rs

    seq = "G" + "A" + HOMOPOLYMER[A + 1 : A + 6]
    reads = [make_read(f"a{i}", seq, A, ((0, 1), (1, 1), (0, 5))) for i in range(5)]
    fa, bam = write_contig(tmp_path, HOMOPOLYMER, reads, "flank")
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", A, "G", "GA", "X")], fa, 5, False, 1, True
    )
    (c,) = count_bam_checked(
        bam, [pv.variant], [None], 20, 20, True, True, True, False, False, False, 1
    )
    assert (c.ad, c.dp) == (5, 5)


def test_a_flank_start_whose_bases_fit_both_stays_uninformative(tmp_path):
    """Guard: the same carriers ending one A earlier (five A's after the G) fit both
    alleles: depth only."""
    from helpers import count_bam_checked, make_read, write_contig

    from gbcms import _rs as gbcms_rs

    seq = "G" + "A" + HOMOPOLYMER[A + 1 : A + 5]
    reads = [make_read(f"a{i}", seq, A, ((0, 1), (1, 1), (0, 4))) for i in range(5)]
    fa, bam = write_contig(tmp_path, HOMOPOLYMER, reads, "flank2")
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", A, "G", "GA", "X")], fa, 5, False, 1, True
    )
    (c,) = count_bam_checked(
        bam, [pv.variant], [None], 20, 20, True, True, True, False, False, False, 1
    )
    assert (c.ad, c.dp) == (0, 5)
