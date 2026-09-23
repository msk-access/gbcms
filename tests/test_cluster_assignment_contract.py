"""Target contract for exclusive assignment at co-annotated tract clusters.

The rule these tests pin: when several same-type indels are annotated within
one repeat-tract / scan-window neighborhood, a molecule belongs to exactly
ONE of them — the row whose canonical (left-aligned position + bases) form
its CIGAR op matches. For every other co-annotated row the molecule is a
DISTINCT allele (neither + partial evidence), never full ALT and never REF.
The windowed S3 shift-tolerance exists for aligner-vs-annotation
representation shifts of the SAME allele; it must not let a molecule hop
between distinct annotated alleles.

Adjudicated on real data (issue #92, 2026-09-22): canonical exclusive
assignment reproduces sign-out exactly at a five-deletion cluster where
per-row windowed counting over-attributed 2-5x.

Runs through the production CLI (pipeline sibling assembly): sibling paths
are exempt from the binned<->legacy parity oracle, so the oracle is NOT used
here. Committed red (xfail-strict) before the fix.
"""

import glob
import random

import pysam
import pytest
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

XFAIL = pytest.mark.xfail(strict=True, reason="cluster exclusive assignment: fix pending")

BAM_CONTIG = "1"
READ_LEN = 100


def _mk_ref(n=1000, seed=17, plants=()):
    rng = random.Random(seed)
    ref = [rng.choice("ACGT") for _ in range(n)]
    for pos, motif in plants:
        ref[pos : pos + len(motif)] = list(motif)
    return "".join(ref)


def _bam(tmp_path, ref, reads, name="cluster.bam"):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": BAM_CONTIG, "LN": len(ref)}]}
    p = tmp_path / name
    with pysam.AlignmentFile(p, "wb", header=header) as fh:
        for a in sorted(reads, key=lambda a: a.reference_start):
            fh.write(a)
    sp = tmp_path / name.replace(".bam", ".s.bam")
    pysam.sort("-o", str(sp), str(p))
    pysam.index(str(sp))
    return sp


def _fasta(tmp_path, ref):
    fa = tmp_path / "ref.fasta"
    fa.write_text(">chr1\n" + ref + "\n")
    pysam.faidx(str(fa))
    return fa


def _vcf(tmp_path, rows):
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(
            f"chr1\t{pos1}\t.\t{ref_al}\t{alt_al}\t.\t.\t.\n" for pos1, ref_al, alt_al in rows
        )
    )
    return vcf


def _run(tmp_path, vcf, bam, fasta, outname="out"):
    outdir = tmp_path / outname
    outdir.mkdir(exist_ok=True)
    result = runner.invoke(
        app,
        [
            "dna",
            "-v",
            str(vcf),
            "-b",
            f"S:{bam}",
            "-f",
            str(fasta),
            "-o",
            str(outdir),
            "--format",
            "maf",
        ],
    )
    assert result.exit_code == 0, result.output
    rows = list(read_maf_output(glob.glob(str(outdir / "*.maf"))[0]))
    for r in rows:
        assert int(r["total_count"]) >= int(r["ref_count"]) + int(r["alt_count"])
        assert int(r["any_alt"]) == int(r["alt_count"]) + int(r["partial_alt"])
    return {(r["Chromosome"], r["vcf_pos"]): r for r in rows}


def _del_reads(ref, p0, length, n, prefix):
    """Reads carrying an exact deletion of `length` starting at 0-based p0."""
    out = []
    for i in range(n):
        s = p0 - 35 - (i % 5)
        left = p0 - s
        seq = ref[s:p0] + ref[p0 + length : p0 + length + (READ_LEN - left)]
        out.append(
            make_read(f"{prefix}{i}", seq, s, ((0, left), (2, length), (0, READ_LEN - left)))
        )
    return out


def _ref_reads(ref, center, n, prefix="wt"):
    out = []
    for i in range(n):
        s = center - 50 + (i % 5)
        out.append(make_read(f"{prefix}{i}", ref[s : s + READ_LEN], s, ((0, READ_LEN),)))
    return out


# ── Cluster geometry ─────────────────────────────────────────────────────
# Four co-annotated 2bp deletions in one scan-window neighborhood, each
# already-canonical (no annotation left-aligns away) and canonically
# DISTINCT. Planted at 0-based 298:  G T A G A C A G T T A C G
#                                    298     302     306     310
#   A: del AG @[300,302)   B: del AC @[302,304)   C: del AG @[304,306)
#   D: del AC @[308,310)
# A/B/C span-chain (strict multi-allelic group today); D sits 3bp past C's
# span — inside the scan window, outside span-touch — so D is the probe for
# the widened (window-overlap) grouping.
# Pre-fix leaks by windowed S3: A<->C (same bases AG, both directions) and
# D<-B (same bases AC, one direction) — carriers count ALT for tract-mates.
PLANT_POS = 298
PLANT = "GTAGACAGTTACG"
D_A = 300  # deletes "AG"
D_B = 302  # deletes "AC"
D_C = 304  # deletes "AG"
D_D = 308  # deletes "AC"


def _cluster_setup(tmp_path, n_a=8, n_b=6, n_c=4, n_d=3, n_wt=10):
    ref = _mk_ref(plants=((PLANT_POS, PLANT),))
    rows = [
        (p0, ref[p0 - 1 : p0 + 2], ref[p0 - 1])  # 1-based POS = p0 (0-based anchor p0-1)
        for p0 in (D_A, D_B, D_C, D_D)
    ]
    reads = (
        _del_reads(ref, D_A, 2, n_a, "a")
        + _del_reads(ref, D_B, 2, n_b, "b")
        + _del_reads(ref, D_C, 2, n_c, "c")
        + _del_reads(ref, D_D, 2, n_d, "d")
        + _ref_reads(ref, D_B, n_wt)
    )
    return ref, rows, reads


def _canon(ref, p0, ln):
    while p0 > 0 and ref[p0 - 1] == ref[p0 + ln - 1]:
        p0 -= 1
    return p0, ref[p0 : p0 + ln]


def test_cluster_alleles_are_canonically_distinct(tmp_path):
    """Meta-guard: the synthetic cluster must mirror the adjudicated real
    case — every annotation already canonical (it must not left-align away
    during prep) and all four canonical alleles DISTINCT. If this fails the
    geometry is wrong, not the engine."""
    ref, _, _ = _cluster_setup(tmp_path)
    forms = {}
    for p in (D_A, D_B, D_C, D_D):
        cp, bases = _canon(ref, p, 2)
        assert cp == p, f"annotation at {p} left-aligns to {cp} — geometry broken"
        forms[(cp, bases)] = p
    assert len(forms) == 4, forms


@XFAIL
def test_cluster_rows_count_only_their_own_molecules(tmp_path):
    """Each co-annotated row counts exactly its own carriers as ALT; every
    tract-mate carrier appears as partial evidence at most."""
    ref, rows, reads = _cluster_setup(tmp_path)
    res = _run(tmp_path, _vcf(tmp_path, rows), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    ads = {pos: int(res[("1", str(pos))]["alt_count"]) for pos, _, _ in rows}
    assert ads == {D_A: 8, D_B: 6, D_C: 4, D_D: 3}, f"per-row ad must be own carriers only: {ads}"


@XFAIL
def test_cluster_sum_bounded_by_distinct_molecules(tmp_path):
    """Invariant: sum of per-row ad across a co-annotated cluster never
    exceeds the number of distinct ALT-carrying molecules (18 here)."""
    ref, rows, reads = _cluster_setup(tmp_path)
    res = _run(tmp_path, _vcf(tmp_path, rows), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    total = sum(int(res[("1", str(pos))]["alt_count"]) for pos, _, _ in rows)
    assert total <= 21, f"sum(ad)={total} exceeds 21 distinct ALT molecules"


@XFAIL
def test_tract_mate_carriers_surface_as_partial(tmp_path):
    """A tract-mate's carriers are structural indel evidence of a DIFFERENT
    allele: they must surface in partial_alt (not vanish, not count REF)."""
    ref, rows, reads = _cluster_setup(tmp_path)
    res = _run(tmp_path, _vcf(tmp_path, rows), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    r_a = res[("1", str(D_A))]
    # carriers of B and C (6+4) are distinct-allele evidence for row A
    assert (
        int(r_a["partial_alt"]) >= 10
    ), f"tract-mate carriers must appear as partial for row A, got {r_a['partial_alt']}"
    # and they never count REF for row A
    assert int(r_a["ref_count"]) == 10, f"only true WT reads count REF, got {r_a['ref_count']}"


def test_distant_same_alleles_unaffected(tmp_path):
    """Guard (green now, green after): two identical-sequence deletions far
    apart (no window overlap, no grouping) keep counting independently."""
    ref = _mk_ref(plants=((PLANT_POS, PLANT), (498, PLANT)))
    d_far = 504  # AG at the second plant, same local geometry as D_C
    rows = [
        (D_C, ref[D_C - 1 : D_C + 2], ref[D_C - 1]),
        (d_far, ref[d_far - 1 : d_far + 2], ref[d_far - 1]),
    ]
    reads = (
        _del_reads(ref, D_C, 2, 5, "n")
        + _del_reads(ref, d_far, 2, 7, "f")
        + _ref_reads(ref, D_C, 4, "w1")
        + _ref_reads(ref, d_far, 4, "w2")
    )
    res = _run(tmp_path, _vcf(tmp_path, rows), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    assert int(res[("1", str(D_C))]["alt_count"]) == 5
    assert int(res[("1", str(d_far))]["alt_count"]) == 7


def test_shifted_self_representation_still_rescued(tmp_path):
    """Guard (green now, green after): a read whose D op is a shifted
    representation of the SAME canonical allele still counts ALT — claiming
    scopes S3 tolerance to self, it must not destroy it."""
    # AT-tract: annotation left-aligns to tract start; reads carry the D
    # right-shifted within the tract — same canonical allele.
    ref = _mk_ref(plants=((320, "CATATATATATG"),))
    anchor = 320  # 'C', 0-based; deletion of AT from the tract
    rows = [(anchor + 1, ref[anchor : anchor + 3], ref[anchor])]
    reads = []
    for i in range(5):
        p0 = 325  # shifted 4bp right of annotated span start (321), same tract
        s = p0 - 35 - i
        left = p0 - s
        seq = ref[s:p0] + ref[p0 + 2 : p0 + 2 + (READ_LEN - left)]
        reads.append(make_read(f"sh{i}", seq, s, ((0, left), (2, 2), (0, READ_LEN - left))))
    reads += _ref_reads(ref, anchor, 4)
    res = _run(tmp_path, _vcf(tmp_path, rows), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    assert (
        int(res[("1", str(anchor + 1))]["alt_count"]) == 5
    ), f"shifted SELF representations must stay ALT: {res}"
