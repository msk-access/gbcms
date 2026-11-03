"""RNA: a spliced read's exon piece that holds only part of a deletion's span is depth only.

Decided (operator, 2026-10-30, #292): a read whose exon piece covers only part of a
deletion's span must not count REF for it, with BAQ on or off. Mechanism: the
piece starts at a junction placed at (or inside) the event's edge and ends inside the
event, so no flank base of its own places it (RJ-3); the aligner put it where its few
bases fit, so they cannot tell REF from ALT spliced at another acceptor. It is the REF
mirror of RJ-17's decided ALT rule (a deletion written right beside the read's splice
has no aligned flank: depth only). Until a splice stopped lowering base quality (RJ-31,
2026-11-03) BAQ's junction penalty hid such a piece's bases; with BAQ off, or once the
penalty was gone, the exact-carrier rule counted it REF. The CIGAR checks' span-aligned
REF testimony already required the whole span.

Out of scope (pinned as they are): a piece that holds the whole span, a piece with a
flank base of its own before the event (judged by its bases, as an unspliced read ending
inside a unique deletion is: decided spec cases `u-8`, `u+10`), and DNA.
"""

import random

import pytest
from helpers import count_checked, make_read, write_contig

from gbcms import _rs

RL = 100


def _bases(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


# exon 1 [0, I0), intron [I0, E), exon 2 from E; exon 2 starts with three bases that
# differ from the three after them, and exon 1 ends likewise
EX1 = _bases(81, 257) + "AGATCC"
INTRON = "GT" + _bases(82, 76) + "AG"
EX2 = "TCAGGCTCC" + _bases(83, 290)
CONTIG = EX1 + INTRON + EX2
I0, E = len(EX1), len(EX1) + len(INTRON)
assert all(CONTIG[E + k] != CONTIG[E + 3 + k] for k in range(3))
assert all(CONTIG[I0 - 3 + k] != CONTIG[I0 - 6 + k] for k in range(3))
# neither deletion shifts: not left (anchor against last deleted base), not right (first
# deleted base against the base after the span)
assert CONTIG[E - 1] != CONTIG[E + 2] and CONTIG[E] != CONTIG[E + 3]
assert CONTIG[I0 - 4] != CONTIG[I0 - 1] and CONTIG[I0 - 3] != CONTIG[I0]


def _count(tmp_path, reads, row, apply_baq, prepared=True, name="s"):
    """(RD, AD, DP) of the row: prepared, the exact-carrier rule judges it; without a
    prepared reference, the CIGAR checks do."""
    tag = f"{name}{int(apply_baq)}{int(prepared)}"
    fa, bam = write_contig(tmp_path, CONTIG, reads, tag)
    pos, ref, alt = row
    assert CONTIG[pos : pos + len(ref)] == ref
    v = _rs.Variant("1", pos, ref, alt, "X")
    if prepared:
        (pv,) = _rs.prepare_variants([v], fa, 5, False, 1, True)
        v = pv.variant
    (c,) = count_checked(bam, [v], mode="rna", reference_fasta=fa, apply_baq=apply_baq, min_mapq=0)
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    return c.rd, c.ad, c.dp


def _unspliced(lo: int, hi: int, n: int, prefix: str) -> list:
    """Reads of the genome covering [lo, hi), each RL long."""
    out = []
    for i in range(n):
        s = lo - 30 + 3 * i
        out.append(make_read(f"{prefix}{i}", CONTIG[s : s + RL], s, ((0, RL),)))
    return out


# ── a deletion of exon 2's first bases, its anchor the intron's last base ──────


def _acceptor_pieces(k: int, n: int, prefix: str) -> list:
    """Reads spliced from exon 1 into exon 2 at the annotated junction, holding the
    first k bases of exon 2 and ending there."""
    out = []
    for i in range(n):
        left = RL - k - i
        s = I0 - left
        seq = CONTIG[s:I0] + CONTIG[E : E + k]
        out.append(make_read(f"{prefix}{i}", seq, s, ((0, left), (3, E - I0), (0, k))))
    return out


DEL3 = (E - 1, CONTIG[E - 1 : E + 3], CONTIG[E - 1])  # deletes exon 2's first three bases


@pytest.mark.parametrize("prepared", [True, False], ids=["exact-carrier", "cigar-checks"])
@pytest.mark.parametrize("apply_baq", [True, False], ids=["baq", "no-baq"])
@pytest.mark.parametrize("k", [1, 2])
def test_a_piece_holding_part_of_the_span_after_its_junction_is_depth_only(
    k, apply_baq, prepared, tmp_path
):
    """Four reads hold k of the three deleted bases after their junction and end there:
    in DP, never REF. Four unspliced reads through the event are REF."""
    reads = _acceptor_pieces(k, 4, "p") + _unspliced(E - 5, E + 8, 4, "u")
    got = _count(tmp_path, reads, DEL3, apply_baq, prepared)
    assert got == (4, 0, 8)


@pytest.mark.parametrize("apply_baq", [True, False], ids=["baq", "no-baq"])
def test_a_piece_holding_the_whole_span_is_not_this_rule(apply_baq, tmp_path):
    """Out of scope: a piece holding all three deleted bases and five bases past them is
    read as before: REF by its bases, with BAQ on or off."""
    reads = _acceptor_pieces(8, 4, "w") + _unspliced(E - 5, E + 8, 4, "u")
    rd, ad, dp = _count(tmp_path, reads, DEL3, apply_baq)
    assert (rd, ad, dp) == (8, 0, 8)


# ── the mirror geometry: a deletion of exon 1's last bases ──────────────────────


def _donor_pieces(k: int, n: int, prefix: str) -> list:
    """Reads that start k bases before exon 1's end, splice there, and run on in exon 2."""
    out = []
    for i in range(n):
        right = RL - k
        seq = CONTIG[I0 - k : I0] + CONTIG[E : E + right]
        out.append(make_read(f"{prefix}{i}", seq, I0 - k, ((0, k), (3, E - I0), (0, right))))
    return out


DEL3_END = (I0 - 4, CONTIG[I0 - 4 : I0], CONTIG[I0 - 4])  # deletes exon 1's last three bases


@pytest.mark.parametrize("apply_baq", [True, False], ids=["baq", "no-baq"])
@pytest.mark.parametrize("k", [1, 2])
def test_a_piece_holding_part_of_the_span_before_its_junction_never_counts_ref(
    k, apply_baq, tmp_path
):
    """The mirror: reads starting inside exon 1's last three (deleted) bases and
    splicing at exon 1's end hold k of them and no flank base. They start past the
    deletion's anchor, so they do not reach the locus: never REF, not in DP (as today)."""
    reads = _donor_pieces(k, 4, "q") + _unspliced(I0 - 8, I0 + 5, 4, "u")
    assert _count(tmp_path, reads, DEL3_END, apply_baq) == (4, 0, 4)


# ── the ALT mirror (RJ-17) ─────────────────────────────────────────────────────


@pytest.mark.parametrize("apply_baq", [True, False], ids=["baq", "no-baq"])
def test_a_deletion_written_beside_its_splice_stays_depth_only(apply_baq, tmp_path):
    """RJ-17: carriers of the deletion of exon 2's first three bases, written as the
    junction then a D, have no aligned flank: depth only, with BAQ on or off."""
    carriers = []
    for i in range(4):
        left = 40 + i
        s = I0 - left
        right = RL - left
        seq = CONTIG[s:I0] + CONTIG[E + 3 : E + 3 + right]
        carriers.append(make_read(f"c{i}", seq, s, ((0, left), (3, E - I0), (2, 3), (0, right))))
    reads = carriers + _unspliced(E - 5, E + 8, 4, "u")
    assert _count(tmp_path, reads, DEL3, apply_baq) == (4, 0, 8)


# ── not this rule: a piece with a flank base of its own before the event ────────


def test_a_piece_with_its_own_flank_base_is_judged_by_its_bases(tmp_path):
    """Out of scope: reads spliced into exon 2 that hold exon 2's first bases as the
    flank of a deletion further in, and end inside its span, are judged by their bases
    as an unspliced read ending inside a unique deletion is (decided spec `u-8`): REF,
    BAQ off."""
    row = (E + 2, CONTIG[E + 2 : E + 6], CONTIG[E + 2])  # deletes exon 2's bases 4-6
    reads = _acceptor_pieces(4, 4, "f") + _unspliced(E - 2, E + 12, 4, "u")
    rd, ad, dp = _count(tmp_path, reads, row, False)
    assert (rd, ad, dp) == (8, 0, 8)
