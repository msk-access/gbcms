"""RNA: a splice never lowers a base's quality (RJ-31, #292).

Heuristic BAQ used to lower the five read bases each side of every CIGAR `N` by 20,
so a base near a junction passed the run's `--min-baseq` only where its quality was 20
above the bar: Q40 data passed at bar 20 and nothing passed at bar 30. A splice moves
no base, so a base is judged by its own quality; the read's insertions and deletions
keep the indel penalty. These tests pin the shapes that decision touches, on spliced
reads without a GTF (with a GTF, BAQ is skipped within five bases of an annotated exon
edge, so the shapes there never differed).

Out of scope, pinned as they are: the penalty of a read's own indel, and the short
piece an aligner writes at an annotated junction where the molecule spliced at a
novel one a few bases on, and a short piece between two junctions (strict xfails, #409).
"""

import random

import pytest
from helpers import count_checked, make_read, write_contig

from gbcms import _rs as gbcms_rs

RL = 100


def _bases(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


# exon 1, a 100-base intron (GTA...AG), exon 2 from E; the intron's first three bases
# differ from exon 2's first three at every base
EX1 = _bases(71, 260) + "CAGTC"
INTRON = "GTA" + _bases(72, 95) + "AG"
EX2 = "CCTGGAACTG" + _bases(73, 290)
CONTIG = EX1 + INTRON + EX2
I0, E = len(EX1), len(EX1) + len(INTRON)  # the intron, [I0, E); exon 2 from E
MRNA = list(range(I0)) + list(range(E, len(CONTIG)))
assert all(CONTIG[I0 + k] != CONTIG[E + k] for k in range(3))


def _cigar(genome: list) -> tuple:
    """The CIGAR of a read whose bases sit at `genome` (None: inserted)."""
    ops: list = []
    prev = None
    for g in genome:
        if g is not None and prev is not None and g != prev + 1:
            gap = g - prev - 1
            ops.append([3 if gap > 20 else 2, gap])
        op = 1 if g is None else 0
        if ops and ops[-1][0] == op:
            ops[-1][1] += 1
        else:
            ops.append([op, 1])
        prev = g if g is not None else prev
    return tuple(map(tuple, ops))


def _spliced(name: str, edits: dict, offset: int, qual: int = 30):
    """A spliced read of the mRNA from exon 1 into exon 2, its first base 40 +
    `offset` bases before the junction; `edits` maps a genome position to its bases
    ("" deletes it; more than one inserts the rest after it)."""
    j = MRNA.index(E)
    seg = MRNA[j - 40 - offset : j - 40 - offset + RL + 1]
    bases, genome = [], []
    for g in seg:
        b = edits.get(g, CONTIG[g])
        if b == "":
            continue
        bases.append(b[0])
        genome.append(g)
        bases += list(b[1:])  # an insertion after g
        genome += [None] * (len(b) - 1)
    bases, genome = bases[:RL], genome[:RL]
    return make_read(name, "".join(bases), genome[0], _cigar(genome), quals=[qual] * len(bases))


def _count(tmp_path, reads, row, apply_baq, min_baseq=20, gtf=False, contig=CONTIG):
    """The row's counts, the four counting invariants asserted."""
    tag = f"s{int(apply_baq)}{min_baseq}{int(gtf)}"
    fa, bam = write_contig(tmp_path, contig, reads, tag)
    pos, ref, alt = row
    assert contig[pos : pos + len(ref)] == ref
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", pos, ref, alt, "X")], fa, 5, False, 1, True
    )
    extra = {}
    if gtf:
        path = tmp_path / f"{tag}.gtf"
        path.write_text(
            "".join(
                f'1\tT\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "G"; transcript_id "T";\n'
                for s, e in ((0, I0), (E, len(CONTIG)))
            )
        )
        extra["gtf_path"] = str(path)
    (c,) = count_checked(
        bam,
        [pv.variant],
        mode="rna",
        reference_fasta=fa,
        apply_baq=apply_baq,
        min_baseq=min_baseq,
        **extra,
    )
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    return c


# ── a delins at exon 1's last bases, read by spliced reads ──────────────────────

DELINS = (I0 - 2, CONTIG[I0 - 2 : I0], "G")
assert DELINS[1] == "TC"


@pytest.mark.parametrize("apply_baq", [True, False], ids=["baq", "no-baq"])
def test_a_delins_at_an_exon_end_on_spliced_reads(apply_baq, tmp_path):
    """Eight spliced REF and eight spliced ALT reads at a 2>1 delins of exon 1's last
    two bases: the ALT's base and the junction's next bases are read at their own
    quality, so every ALT read counts, with BAQ on or off."""
    reads = [
        _spliced(f"{k}{i}", {I0 - 2: "G", I0 - 1: ""} if k == "A" else {}, 3 * i)
        for k in "RA"
        for i in range(8)
    ]
    c = _count(tmp_path, reads, DELINS, apply_baq)
    assert (c.rd, c.ad) == (8, 8)


# ── an SNV two bases into exon 2, read by spliced reads ─────────────────────────

SNV = (E + 2, CONTIG[E + 2], "A" if CONTIG[E + 2] != "A" else "C")


def _snv_reads(extra: dict | None = None) -> list:
    alt = {SNV[0]: SNV[2]}
    return [
        _spliced(f"{k}{i}", {**(extra or {}), **(alt if k == "A" else {})}, 3 * i)
        for k in "RA"
        for i in range(8)
    ]


@pytest.mark.parametrize("gtf", [False, True], ids=["no-gtf", "gtf"])
@pytest.mark.parametrize("min_baseq", [10, 20, 30])
@pytest.mark.parametrize("apply_baq", [True, False], ids=["baq", "no-baq"])
def test_an_snv_beside_a_junction_is_read_at_its_own_quality(apply_baq, min_baseq, gtf, tmp_path):
    """The reads are Q30. The SNV's base lies two bases past their junction: it counts
    ALT at every bar up to 30, as it does with BAQ off. (The old junction penalty put
    it at Q10, under the bar at 20 and 30.)"""
    c = _count(tmp_path, _snv_reads(), SNV, apply_baq, min_baseq, gtf)
    assert (c.rd, c.ad) == (8, 8)


def test_an_snv_beside_a_junction_and_a_deletion_keeps_the_indel_penalty(tmp_path):
    """A 1-base deletion every read carries three bases past the SNV: its own penalty
    reaches the SNV's base, so BAQ on leaves the SNV depth only; BAQ off counts it. The
    junction adds nothing to either."""
    reads = _snv_reads({E + 5: ""})
    c = _count(tmp_path, reads, SNV, True)
    assert (c.rd, c.ad) == (0, 0)
    c = _count(tmp_path, reads, SNV, False)
    assert (c.rd, c.ad) == (8, 8)


# ── the exon-edge exception: what it still lifts ────────────────────────────────


def _tx_counts(c) -> dict:
    """{transcript: (AD, RD, DP)} from `T:AD,RD,DP|...`."""
    out = {}
    for entry in filter(None, c.transcript_read_counts.split("|")):
        tx, nums = entry.rsplit(":", 1)
        out[tx] = tuple(int(x) for x in nums.split(","))
    return out


def test_the_exon_edge_exception_lifts_a_reads_own_indel_penalty(tmp_path):
    """The SNV sits two bases inside exon 2's first base; every read has a 1-base deletion
    three bases on, whose own penalty reaches the SNV's base. Without a GTF BAQ applies
    and the SNV is depth only (0, 0). With a GTF the SNV lies within five bases of an
    annotated exon edge, BAQ is skipped, and every read counts (8, 8): in the main counts,
    the per-transcript counts and ASJD alike. (Reads without an indel of their own show
    nothing of the exception: a splice lowers no quality, with a GTF or without.)"""
    reads = _snv_reads({E + 5: ""})
    plain = _count(tmp_path, reads, SNV, True)
    assert (plain.rd, plain.ad) == (0, 0)
    edge = _count(tmp_path, reads, SNV, True, gtf=True)
    assert (edge.rd, edge.ad) == (8, 8)
    assert _tx_counts(edge) == {"T": (8, 8, 16)}
    assert (edge.asjd_n_ref_total, edge.asjd_n_alt_total) == (8, 8)


# ── a short piece at an annotated junction, judged at a delins ──────────────────
# Its own contig. A hom-REF molecule splices at a novel acceptor a few bases from the
# annotated one and holds only 6 or 8 bases past it; the aligner writes it at the
# annotated acceptor, so its bases there are the intron's last bases and the exon's
# first (mismatches). Delins GC>A-like rows in exon 2's first bases: with the old
# penalty the junction's ten read bases were masked and a masked base fits any allele,
# so the few bases left fitted the ALT by chance and the read counted ALT. Read at their
# own quality they fit neither allele: partial evidence, as with BAQ off.

X1 = _bases(7, 300)
XI = "GT" + _bases(8, 96) + "AG"
X2 = _bases(9, 300)
XCONTIG = X1 + XI + X2
XI0, XE = len(X1), len(X1) + len(XI)
XMRNA = list(range(XI0)) + list(range(XE, len(XCONTIG)))
Q40 = 40


def _x_spliced(name: str, j_before: int):
    """A spliced REF read, `j_before` bases before the junction."""
    j = XMRNA.index(XE)
    seg = XMRNA[j - j_before : j - j_before + RL]
    return make_read(name, "".join(XCONTIG[g] for g in seg), seg[0], _cigar(seg), quals=[Q40] * RL)


def _short_piece(name: str, past: int, shift: int):
    """A read of a molecule spliced at exon 2's base `shift` from the annotated acceptor,
    holding `past` bases past it, written at the annotated acceptor."""
    start = XI0 - (RL - past)
    seq = XCONTIG[start:XI0] + XCONTIG[XE + shift : XE + shift + past]
    cigar = ((0, RL - past), (3, XE - XI0), (0, past))
    return make_read(name, seq, start, cigar, quals=[Q40] * RL)


def _x_delins(side: str, k: int):
    p = XI0 - 2 - k if side == "don" else XE + k
    ref = XCONTIG[p : p + 2]
    return p, ref, next(b for b in "ACGT" if b not in ref)


SHORT_PIECES = {
    "acceptor 0, 6 past, 6 early": ("acc", 0, 6, -6),
    "acceptor 1, 6 past, 3 early": ("acc", 1, 6, -3),
    "acceptor 2, 8 past, 4 early": ("acc", 2, 8, -4),
    "acceptor 2, 8 past, 6 early": ("acc", 2, 8, -6),
    "acceptor 3, 8 past, 6 early": ("acc", 3, 8, -6),
    "donor 0, 6 past, 6 early": ("don", 0, 6, -6),
}


@pytest.mark.parametrize("apply_baq", [True, False], ids=["baq", "no-baq"])
@pytest.mark.parametrize("shape", list(SHORT_PIECES))
def test_a_short_piece_fitting_neither_allele_of_a_delins_is_partial(shape, apply_baq, tmp_path):
    """Q40 reads at bar 30: ten spliced REF reads and one short piece. The piece's bases
    fit neither allele of the delins; it is partial evidence, never ALT, with BAQ on or
    off (invariant 7: the read's own bases must carry the ALT)."""
    side, k, past, shift = SHORT_PIECES[shape]
    row = _x_delins(side, k)
    reads = [_x_spliced(f"b{i}", 10 + 8 * i) for i in range(10)] + [_short_piece("t", past, shift)]
    c = _count(tmp_path, reads, row, apply_baq, 30, contig=XCONTIG)
    assert (c.rd, c.ad, c.partial_alt) == (10, 0, 1)


# ── the limit: a short piece at an annotated junction, judged at an SNV (#409) ──
# Pieces too short to place their own bases are RJ-3's class. The old penalty hid their
# bases (and half the depth beside the junction); now they count at their quality. The
# rule that keeps them depth only is ticketed: this flips when it lands.

SHORT = (E + 1, CONTIG[E + 1], CONTIG[E + 4])
assert SHORT[1] != SHORT[2]


@pytest.mark.xfail(
    strict=True,
    reason="#409: a piece with 3 bases past its junction, written at the annotated acceptor "
    "where its molecule spliced three bases on, counts ALT at a hom-REF row",
)
def test_short_pieces_at_an_annotated_acceptor_are_depth_only(tmp_path):
    """Eight spliced REF reads; four reads whose molecules splice three bases past the
    annotated acceptor and hold three bases past it, written at the annotated one: their
    bases there are exon 2's from three bases on, mismatches. Any ALT is false."""
    reads = [_spliced(f"S{i}", {}, 3 * i) for i in range(8)]
    for i in range(4):
        reads.append(_spliced(f"N{i}", {E + k: CONTIG[E + 3 + k] for k in range(3)}, 57))
    c = _count(tmp_path, reads, SHORT, True)
    assert (c.rd, c.ad) == (8, 0)


# ── the same limit for a piece between two junctions (#409) ─────────────────────
# A piece too short to place its own bases is RJ-3's class wherever it lies in a read: at
# its end (above) or between two junctions. Here an aligner writes a 6-base piece at
# 500-506 between two introns; the bases it holds are misplaced, and one of them
# (502) reads as the ALT at a hom-REF row. Ten unspliced reads show REF there. With the
# old penalty both junctions' radii stacked on the piece and masked it whole.

M1 = _bases(1, 300)
MI = "GT" + _bases(2, 396) + "AG"  # the introns' span, [300, 700)
M3 = _bases(3, 300)
MCONTIG = M1 + MI + M3
MPIECE = 500  # the written piece [500, 506)
MPOS = MPIECE + 2


@pytest.mark.xfail(
    strict=True,
    reason="#409: a 6-base piece between two junctions, its bases misplaced, counts ALT at a "
    "hom-REF row (develop: depth only, the stacked junction penalties masked it)",
)
@pytest.mark.parametrize("qual,min_baseq", [(30, 20), (40, 30), (40, 20)])
def test_a_short_piece_between_two_junctions_is_depth_only(qual, min_baseq, tmp_path):
    """Ten unspliced REF reads and four reads `47M 200N 6M 194N 47M` whose 6-base middle
    piece carries a base that is not the reference's. Any ALT is false."""
    alt = next(b for b in "ACGT" if b != MCONTIG[MPOS])
    reads = []
    for i in range(10):
        st = MPOS - 50 + 3 * i
        reads.append(make_read(f"u{i}", MCONTIG[st : st + RL], st, ((0, RL),), quals=[qual] * RL))
    for i in range(4):
        piece = list(MCONTIG[MPIECE : MPIECE + 6])
        piece[MPOS - MPIECE] = alt
        seq = MCONTIG[253:300] + "".join(piece) + MCONTIG[700:747]
        cigar = ((0, 47), (3, MPIECE - 300), (0, 6), (3, 700 - (MPIECE + 6)), (0, 47))
        reads.append(make_read(f"m{i}", seq, 253, cigar, quals=[qual] * RL))
    c = _count(tmp_path, reads, (MPOS, MCONTIG[MPOS], alt), True, min_baseq, contig=MCONTIG)
    assert (c.rd, c.ad) == (10, 0)
