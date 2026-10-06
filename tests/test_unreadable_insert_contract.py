"""C35 #240: an insertion's ALT needs one of the read's own inserted bases read.

A read whose inserted bases are all N (fgbio masks duplex disagreements to N at
Q2) or all below min BQ carries the insertion's length but not its sequence. The
engine credited it ALT on length alone: the strict path found the bases
unverifiable and handed the read to Phase 3, where REF pays for a gap, under both
backends; shifted placements went the same way; the deleted-anchor and split-op
reading (`read_bases_fit_alt`) needed one readable base anywhere between the
flanks; truncations ignored BQ. On the RC DNA set that was 127 reads, mostly 1bp
homopolymer insertions with a duplex-masked N inside the run.

Contract (operator, 2026-10-06): such a read is partial (depth and partial_alt,
never REF or ALT), on every path, as the post-splice path already does. A read
whose readable inserted bases match the ALT stays ALT. "Readable" is one gate for
every base: not N and at or above min BQ. The read census holds the same rule.
"""

import random

import pytest
from census import Verdict, assert_matches, census, judge, window_for
from helpers import count_checked, make_read, write_contig
from read_judgment_cases import _build, _mask_inserted

from gbcms import _rs

L = 100  # read length
U, H, T = 300, 600, 800  # non-repeat 5bp insertion, homopolymer +A, 8bp insertion
INS5, INS8 = "ATGGA", "ATGGACTC"
BACKENDS = ["pairhmm", "sw"]


def _contig():
    rng = random.Random(240)
    c = [rng.choice("ACGT") for _ in range(1000)]
    c[U], c[U + 1] = "C", "C"  # neither end of ATGGA repeats a neighbour: no shift
    c[H : H + 8] = "G" + "A" * 6 + "C"  # G>GA left-aligned at H
    c[T], c[T + 1] = "G", "G"
    return "".join(c)


CONTIG = _contig()


def _variant(fa, pos, ref, alt):
    (pv,) = _rs.prepare_variants([_rs.Variant("1", pos, ref, alt, "INS")], fa, 5, False, 1, True)
    v = pv.variant
    assert (pv.gbcms_status, v.pos, v.ref_allele, v.alt_allele) == ("PASS", pos, ref, alt)
    return v


def _quals(n, masked, q):
    out = [37] * n
    for i in masked:
        out[i] = q
    return out


def ins_read(name, i, anchor, ins, *, shift=0, low=(), q=2):
    """A read with `ins` inserted `shift` bases after the anchor (a later junction of
    a run); the inserted bases at offsets `low` get quality `q`."""
    s = anchor - 40 - (i % 5)
    j = anchor + 1 + shift
    left = j - s
    seq = CONTIG[s:j] + ins + CONTIG[j : j + L - left - len(ins)]
    quals = _quals(L, [left + k for k in low], q)
    cigar = ((0, left), (1, len(ins)), (0, L - left - len(ins)))
    return make_read(name, seq, s, cigar, flag=16 * (i % 2), quals=quals)


def ref_read(name, i, anchor):
    s = anchor - 40 - (i % 5)
    return make_read(name, CONTIG[s : s + L], s, ((0, L),), flag=16 * (i % 2))


def deleted_anchor_read(name, i, anchor, ins, *, low=(), q=2):
    """The aligner deletes the anchor and re-inserts it with the insert: M D(1) I M."""
    s = anchor - 40 - (i % 5)
    left = anchor - s
    body = CONTIG[anchor] + ins
    seq = CONTIG[s:anchor] + body + CONTIG[anchor + 1 : anchor + 1 + L - left - len(body)]
    quals = _quals(L, [left + 1 + k for k in low], q)
    cigar = ((0, left), (2, 1), (1, len(body)), (0, L - left - len(body)))
    return make_read(name, seq, s, cigar, flag=16 * (i % 2), quals=quals)


def _count(tmp_path, name, reads, pos, ref, alt, backend="pairhmm", with_census=True):
    fa, bam = write_contig(tmp_path, CONTIG, reads, name)
    v = _variant(fa, pos, ref, alt)
    (counts,) = count_checked(bam, [v], alignment_backend=backend)
    if with_census:
        assert_matches(counts, census(bam, CONTIG, v), name)
    return counts.rd, counts.ad, counts.partial_alt


@pytest.mark.parametrize("backend", BACKENDS)
def test_an_unreadable_insert_at_the_junction_is_partial(tmp_path, backend):
    """Four each: clean ALT, NNNNN at Q2, the right bases at Q2, REF."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS5
    reads = [ins_read(f"c{i}", i, a, INS5) for i in range(4)]
    reads += [ins_read(f"n{i}", i, a, "N" * 5, low=range(5)) for i in range(4)]
    reads += [ins_read(f"q{i}", i, a, INS5, low=range(5)) for i in range(4)]
    reads += [ref_read(f"r{i}", i, a) for i in range(4)]
    assert _count(tmp_path, f"j{backend}", reads, a, ref, alt, backend) == (4, 4, 8)


def test_a_partly_masked_insert_stays_alt(tmp_path):
    """Guard: readable inserted bases that match the ALT are evidence (two masked of
    five, as N at Q2 or as letters at Q5)."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS5
    reads = [ins_read(f"n{i}", i, a, "NN" + INS5[2:], low=(0, 1)) for i in range(4)]
    reads += [ins_read(f"q{i}", i, a, INS5, low=(0, 1), q=5) for i in range(4)]
    reads += [ref_read(f"r{i}", i, a) for i in range(4)]
    assert _count(tmp_path, "pm", reads, a, ref, alt) == (4, 8, 0)


@pytest.mark.parametrize("backend", BACKENDS)
def test_an_unreadable_insert_shifted_in_a_run_is_partial(tmp_path, backend):
    """G>GA before A6: four each of a clean +A at the left junction, a readable +A
    written two bases into the run, the run's extra base masked there (N at Q2, A at
    Q5; the aligner writes the masked base as the insertion), and REF."""
    a, ref, alt = H, "G", "GA"
    reads = [ins_read(f"c{i}", i, a, "A") for i in range(4)]
    reads += [ins_read(f"s{i}", i, a, "A", shift=2) for i in range(4)]
    reads += [ins_read(f"n{i}", i, a, "N", shift=2, low=(0,)) for i in range(4)]
    reads += [ins_read(f"q{i}", i, a, "A", shift=2, low=(0,), q=5) for i in range(4)]
    reads += [ref_read(f"r{i}", i, a) for i in range(4)]
    assert _count(tmp_path, f"hp{backend}", reads, a, ref, alt, backend) == (4, 8, 8)


def test_an_unreadable_insert_after_a_deleted_anchor_is_partial(tmp_path):
    """M D(1) I M: the read deletes the anchor and re-inserts it with the insert. Its
    readable re-inserted anchor is not the insert: four readable carriers are ALT,
    four with the insert masked are partial."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS5
    reads = [deleted_anchor_read(f"d{i}", i, a, INS5) for i in range(4)]
    reads += [deleted_anchor_read(f"n{i}", i, a, "N" * 5, low=range(5)) for i in range(4)]
    reads += [ref_read(f"r{i}", i, a) for i in range(4)]
    assert _count(tmp_path, "da", reads, a, ref, alt) == (4, 4, 4)


def test_an_unreadable_truncation_is_partial(tmp_path):
    """A truncated copy of an 8bp insert (its first five bases at the junction):
    readable it is the same event (ALT, the identity band unchanged), unreadable
    (letters at Q5) it is partial."""
    a, ref, alt = T, CONTIG[T], CONTIG[T] + INS8
    reads = [ins_read(f"c{i}", i, a, INS8) for i in range(4)]
    reads += [ins_read(f"t{i}", i, a, INS8[:5]) for i in range(4)]
    reads += [ins_read(f"q{i}", i, a, INS8[:5], low=range(5), q=5) for i in range(4)]
    reads += [ref_read(f"r{i}", i, a) for i in range(4)]
    # No census: it knows two haplotypes of exact length, and the truncation band is
    # engine policy (a readable truncation is ALT there, neither in the census).
    assert _count(tmp_path, "tr", reads, a, ref, alt, with_census=False) == (4, 8, 4)


def test_the_census_holds_the_same_rule(tmp_path):
    """The census judges by bases, two haplotypes, so length alone made an unreadable
    insert fit the ALT. Its policy now: an insertion's ALT verdict needs the read's own
    inserted bases (its I ops between the anchors) to hold more readable bases than
    the bases it deletes there."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS5
    fa, _ = write_contig(tmp_path, CONTIG, [], "cj")
    w = window_for(_variant(fa, a, ref, alt), CONTIG)
    assert judge(ins_read("c", 0, a, INS5), w) == Verdict.ALT
    assert judge(ins_read("n", 0, a, "N" * 5, low=range(5)), w) == Verdict.UNREADABLE
    assert judge(ins_read("q", 0, a, INS5, low=range(5), q=5), w) == Verdict.UNREADABLE
    assert judge(ins_read("p", 0, a, INS5, low=(0, 1), q=5), w) == Verdict.ALT
    assert judge(deleted_anchor_read("d", 0, a, INS5), w) == Verdict.ALT
    assert judge(deleted_anchor_read("e", 0, a, "N" * 5, low=range(5)), w) == Verdict.UNREADABLE


# ── From the review ───────────────────────────────────────────────────────────
def _shaped(name, i, anchor, events, mask=None):
    """A read from `anchor - 40 - i % 5` with read-judgment events (pos, "I", bases)
    or (pos, "D", n), its inserted bases optionally masked (kind, first, count)."""
    s = anchor - 40 - (i % 5)
    seq, cigar = _build(CONTIG, s, L, events)
    read = make_read(name, seq, s, cigar, flag=16 * (i % 2))
    if mask:
        _mask_inserted(read, mask)
    return read


@pytest.mark.xfail(strict=True, reason="a far deletion offsets the anchor")
def test_an_unrelated_deletion_does_not_cancel_the_rule(tmp_path):
    """The deleted-anchor shape plus an unrelated D2 twenty bases on: a deletion
    outside the flanks re-inserts nothing, so it must not offset the anchor."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS5
    events = [(a, "D", 1), (a + 1, "I", CONTIG[a] + INS5), (a + 21, "D", 2)]
    reads = [_shaped(f"n{i}", i, a, events, ("N", 1, 5)) for i in range(4)]
    reads += [_shaped(f"c{i}", i, a, events) for i in range(4)]
    assert _count(tmp_path, "ud", reads, a, ref, alt) == (0, 4, 4)


@pytest.mark.xfail(strict=True, reason="flagged anywhere in the scan window")
@pytest.mark.parametrize("backend", BACKENDS)
def test_an_unreadable_insert_outside_the_window_is_a_separate_event(tmp_path, backend):
    """An unreadable same-length insert outside the discrimination window cannot be
    the variant at another placement: a separate event, REF where the window is
    reference (RJ-8), as a readable one is. A +5 four bases past the junction; a +1
    after the run's far flank."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS5
    reads = [ins_read(f"n{i}", i, a, "N" * 5, shift=4, low=range(5)) for i in range(4)]
    reads += [ins_read(f"c{i}", i, a, INS5, shift=4) for i in range(4)]
    assert _count(tmp_path, f"ow{backend}", reads, a, ref, alt, backend) == (8, 0, 0)
    reads = [ins_read(f"h{i}", i, H, "N", shift=7, low=(0,)) for i in range(4)]
    assert _count(tmp_path, f"oh{backend}", reads, H, "G", "GA", backend) == (4, 0, 0)


@pytest.mark.xfail(strict=True, reason="the census nets only deletions between its anchors")
def test_the_census_nets_deletions_across_the_bases_it_reads(tmp_path):
    """M D3 I M: the read deletes the anchor and the two bases before it and
    re-inserts them with a masked insert. Its readable inserted bases are the three
    it deleted, so none is the insert: partial, engine and census alike."""
    a, ref, alt = U, CONTIG[U], CONTIG[U] + INS5
    events = [(a - 2, "D", 3), (a + 1, "I", CONTIG[a - 2 : a + 1] + INS5)]
    reads = [_shaped(f"n{i}", i, a, events, ("N", 3, 5)) for i in range(4)]
    assert _count(tmp_path, "d3", reads, a, ref, alt) == (0, 0, 4)
