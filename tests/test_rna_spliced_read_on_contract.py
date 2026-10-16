"""RNA pure indels beside an exon edge: a spliced read reads on along its own mRNA.

A read whose windows fit both alleles, or whose fit rests on too few bases, is read
on past the windows (RJ-27). For a read spliced just past them, the bases there are
the next exon's: judged against the genome's continuation (the intron), a REF read
misses REF and fits a shifted allele by chance, and a read whose bases at the
junction BAQ masks has nothing left to decide on. A read spliced within the reach of
its windows' reading is judged on the reference spliced along its junctions, as a
read spliced inside the windows is, so it counts as the same read without a junction
does. No GTF: BAQ is suppressed only within five bases of an annotated exon edge.
"""

import random

import pytest
from helpers import count_checked, make_read, write_contig

from gbcms import _rs as gbcms_rs

RL = 100


def _flank(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


EX1 = _flank(41, 260) + "CAGTC"
INTRON = "GT" + _flank(42, 36) + "AG"
EX2 = "ACCTGGA" + _flank(43, 300)
CONTIG = EX1 + INTRON + EX2
I0, E = len(EX1), len(EX1) + len(INTRON)  # the intron, [I0, E); exon 2 from E
MRNA = list(range(I0)) + list(range(E, len(CONTIG)))


def _hap(p: int, ref: str, alt: str) -> list:
    """The mRNA with `ref`>`alt` at genome `p`: (base, genome position, None inserted)."""
    out = []
    for g in MRNA:
        if p < g <= p + len(ref) - len(alt):
            continue
        out.append((CONTIG[g], g))
        if g == p:
            out += [(b, None) for b in alt[1 : len(alt) - len(ref) + 1]]
    return out


def _read(name: str, seg: list, ins_after_n: bool = False):
    """A read of the haplotype slice `seg`, its CIGAR from the genome positions; an
    insertion right before a splice is written after it with `ins_after_n`."""
    ops: list = []
    prev = None
    for _, g in seg:
        if g is not None and prev is not None and g != prev + 1:
            gap = g - prev - 1
            ops.append([3 if gap > 20 else 2, gap])
        op = 1 if g is None else 0
        if ops and ops[-1][0] == op:
            ops[-1][1] += 1
        else:
            ops.append([op, 1])
        prev = g if g is not None else prev
    for i in range(len(ops) - 1):
        if ins_after_n and ops[i][0] == 1 and ops[i + 1][0] == 3:
            ops[i], ops[i + 1] = ops[i + 1], ops[i]
    start = next(g for _, g in seg if g is not None)
    return make_read(name, "".join(b for b, _ in seg), start, tuple(map(tuple, ops)))


def _count(p: int, ref: str, alt: str, apply_baq: bool, tmp_path, ins_after_n=False):
    """Eight spliced REF and eight spliced ALT reads from exon 1 into exon 2."""
    assert CONTIG[p : p + len(ref)] == ref
    haps = {"R": _hap(p, ref, ref), "A": _hap(p, ref, alt)}
    j = next(i for i, (_, g) in enumerate(haps["R"]) if g == E)
    reads = [
        _read(f"{k}{i}", h[j - 40 + 3 * i : j - 40 + 3 * i + RL], ins_after_n)
        for k, h in haps.items()
        for i in range(8)
    ]
    fa, bam = write_contig(tmp_path, CONTIG, reads, f"rna{int(apply_baq)}")
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", p, ref, alt, "X")], fa, 5, False, 1, True
    )
    (c,) = count_checked(bam, [pv.variant], mode="rna", reference_fasta=fa, apply_baq=apply_baq)
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    return c.rd, c.ad, c.partial_alt


# (anchor, REF, ALT) beside exon 1's end: an insertion of TG after the anchor two to
# four bases before the edge, and a 2-base deletion ending one base before it
NEAR_EXON_END = {
    "ins two before the end": (I0 - 2, "T", "TTG"),
    "ins three before the end": (I0 - 3, "G", "GTG"),
    "ins four before the end": (I0 - 4, "A", "ATG"),
    "del ending one before the end": (I0 - 4, "AGT", "A"),
}


@pytest.mark.parametrize("apply_baq", [False, True], ids=["no-baq", "baq"])
@pytest.mark.parametrize("shape", list(NEAR_EXON_END))
def test_spliced_reads_near_an_exon_end_read_on_into_the_next_exon(shape, apply_baq, tmp_path):
    """Spliced just past the windows: read on along the mRNA, REF reads are REF and
    ALT reads ALT, with BAQ masking the bases at the junction or not."""
    assert _count(*NEAR_EXON_END[shape], apply_baq, tmp_path) == (8, 8, 0)


# beside exon 2's start: the anchor is the exon's first base (#291)
NEAR_EXON_START = {
    "ins after the first base": (E, "A", "ATG"),
    "del after the first base": (E, "ACC", "A"),
}


@pytest.mark.parametrize("shape", list(NEAR_EXON_START))
def test_spliced_reads_at_an_exon_first_base_anchor_read_on_past_the_masked_junction(
    shape, tmp_path
):
    """BAQ masks the five bases each side of the junction, the whole windows: the
    spliced reads read on along exon 2, as they do with BAQ off."""
    assert _count(*NEAR_EXON_START[shape], True, tmp_path) == (8, 8, 0)
    assert _count(*NEAR_EXON_START[shape], False, tmp_path) == (8, 8, 0)


# at exon 1's very end: TG inserted after its last base, and its last two bases deleted
AT_EXON_END = {
    "ins after the last base": ((I0 - 1, "C", "CTG"), (8, 8, 0)),
    # the ALT reads write the deletion as a splice two bases early: REF spliced at a
    # donor there reads the same, so they show neither allele (RJ-17)
    "del of the last two bases": ((I0 - 3, "GTC", "G"), (8, 0, 0)),
}


@pytest.mark.parametrize("apply_baq", [False, True], ids=["no-baq", "baq"])
@pytest.mark.parametrize("shape", list(AT_EXON_END))
def test_an_exon_end_at_the_event_leaves_its_left_flank_whole(shape, apply_baq, tmp_path):
    """An exon end where the event starts cuts its right flank, not its left: the
    windows read from the left and on along the mRNA, past the bases BAQ masks at the
    junction, instead of from the right, where the alleles differ only there."""
    variant, counts = AT_EXON_END[shape]
    assert _count(*variant, apply_baq, tmp_path) == counts


@pytest.mark.parametrize("apply_baq", [False, True], ids=["no-baq", "baq"])
@pytest.mark.parametrize("written", ["before the N", "after the N"])
def test_an_insertion_at_an_exon_end_counts_however_it_is_written_beside_the_splice(
    written, apply_baq, tmp_path
):
    """TG after exon 1's last base sits between exon 1 and exon 2 on the mRNA whether the
    aligner writes it before the splice or after: the same bases, the same counts. An
    insertion written after the splice is the event on the spliced reference, not a
    separate indel past the windows that would end their reading on."""
    after = written == "after the N"
    assert _count(I0 - 1, "C", "CTG", apply_baq, tmp_path, ins_after_n=after) == (8, 8, 0)
