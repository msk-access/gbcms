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

from gbcms import _rs as gbcms_rs
from gbcms.cli import app

runner = CliRunner()


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
    # Key on the normalized OUTPUT Start_Position so the MAF and VCF input
    # paths index identically (vcf_pos exists only for VCF input, and the
    # chromosome label is passed through differently per input format).
    return {int(r["Start_Position"]): r for r in rows}


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


def test_cluster_rows_count_only_their_own_molecules(tmp_path):
    """Each co-annotated row counts exactly its own carriers as ALT; every
    tract-mate carrier appears as partial evidence at most."""
    ref, rows, reads = _cluster_setup(tmp_path)
    res = _run(tmp_path, _vcf(tmp_path, rows), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    ads = {pos: int(res[pos + 1]["alt_count"]) for pos, _, _ in rows}
    assert ads == {D_A: 8, D_B: 6, D_C: 4, D_D: 3}, f"per-row ad must be own carriers only: {ads}"


def test_cluster_sum_bounded_by_distinct_molecules(tmp_path):
    """Invariant: sum of per-row ad across a co-annotated cluster never
    exceeds the number of distinct ALT-carrying molecules (18 here)."""
    ref, rows, reads = _cluster_setup(tmp_path)
    res = _run(tmp_path, _vcf(tmp_path, rows), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    total = sum(int(res[pos + 1]["alt_count"]) for pos, _, _ in rows)
    assert total <= 21, f"sum(ad)={total} exceeds 21 distinct ALT molecules"


def test_tract_mate_carriers_surface_as_partial(tmp_path):
    """A tract-mate's carriers are structural indel evidence of a DIFFERENT
    allele: they must surface in partial_alt (not vanish, not count REF)."""
    ref, rows, reads = _cluster_setup(tmp_path)
    res = _run(tmp_path, _vcf(tmp_path, rows), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    r_a = res[D_A + 1]
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
    assert int(res[D_C + 1]["alt_count"]) == 5
    assert int(res[d_far + 1]["alt_count"]) == 7


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
        int(res[anchor + 2]["alt_count"]) == 5
    ), f"shifted SELF representations must stay ALT: {res}"


# ── Insertion tract cluster ──────────────────────────────────────────────
# Two co-annotated insertions of the SAME bases at nearby, canonically
# distinct positions. A cross carrier's I op S3-matches the other row's
# windowed scan pre-fix; post-fix the true row (span-explanation cost 0)
# claims it and the other row records partial evidence.
INS_PLANT_POS = 298
INS_PLANT = "GTCGACTGTC"  # ensures ins "CA" is canonically stable at both sites
INS_A = 301  # 0-based: insertion between ref[300] and ref[301]
INS_B = 307


def _ins_reads(ref, p0, ins, n, prefix):
    out = []
    for i in range(n):
        s = p0 - 35 - (i % 5)
        left = p0 - s
        seq = ref[s:p0] + ins + ref[p0 : p0 + (READ_LEN - left - len(ins))]
        out.append(
            make_read(
                f"{prefix}{i}", seq, s, ((0, left), (1, len(ins)), (0, READ_LEN - left - len(ins)))
            )
        )
    return out


def test_insertion_cluster_exclusive_assignment(tmp_path):
    """INS twin of the deletion cluster: each row keeps exactly its own
    carriers; cross carriers surface as partial, never full ALT."""
    ref = _mk_ref(plants=((INS_PLANT_POS, INS_PLANT),))
    ins = "CA"
    for p in (INS_A, INS_B):  # meta-guard: canonically stable annotations
        assert ref[p - 1] != ins[-1], f"ins at {p} would left-align — geometry broken"
    rows = [
        (p, ref[p - 1], ref[p - 1] + ins) for p in (INS_A, INS_B)
    ]  # 1-based POS p, anchor ref[p-1]
    reads = (
        _ins_reads(ref, INS_A, ins, 7, "ia")
        + _ins_reads(ref, INS_B, ins, 4, "ib")
        + _ref_reads(ref, INS_A, 6)
    )
    res = _run(tmp_path, _vcf(tmp_path, rows), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    r_a, r_b = res[INS_A], res[INS_B]
    assert int(r_a["alt_count"]) == 7, f"row A must keep only its own carriers: {r_a['alt_count']}"
    assert int(r_b["alt_count"]) == 4, f"row B must keep only its own carriers: {r_b['alt_count']}"
    assert int(r_a["alt_count"]) + int(r_b["alt_count"]) <= 11


# ── Equivalent double-annotation guard ───────────────────────────────────
def test_equivalent_representations_both_keep_carriers(tmp_path):
    """Guard (adjudicated on real double-annotated data): when two
    co-annotated rows are equivalent representations of the SAME event —
    a delins TT>GTC and (T>G mismatch + ins C) describing one haplotype —
    the contest's equal-cost tie must NOT zero either row. Carriers whose
    reads literally show the shared haplotype count on both rows."""
    ref = _mk_ref(plants=((298, "ACTTGAC"),))
    p0 = 300  # 0-based: ref TT at [300, 302)
    assert ref[300:302] == "TT"
    rows = [
        (p0 + 1, "TT", "GTC"),  # complex representation, 1-based POS 301
        (p0 + 2, ref[p0 + 1], ref[p0 + 1] + "C"),  # ins C after 302 (anchor T)
    ]
    reads = []
    for i in range(6):  # carriers show G T C: mismatch at 300, ins C after 301
        s0 = p0 - 35 - i
        left = p0 - s0
        seq = ref[s0:p0] + "GT" + "C" + ref[p0 + 2 : p0 + 2 + (READ_LEN - left - 3)]
        reads.append(
            make_read(f"eq{i}", seq, s0, ((0, left + 2), (1, 1), (0, READ_LEN - left - 3)))
        )
    reads += _ref_reads(ref, p0, 5)
    res = _run(tmp_path, _vcf(tmp_path, rows), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    ads = sorted(int(r["alt_count"]) for r in res.values())
    # Neither row may be zeroed; the exact split (both-full for true
    # equivalence, or exclusive when representations differ) is engine
    # policy, but carriers must be visible as ad or partial on every row.
    for r in res.values():
        assert (
            int(r["alt_count"]) + int(r["partial_alt"]) >= 6
        ), f"carriers vanished from a row: ad={r['alt_count']} partial={r['partial_alt']}"
    assert max(ads) >= 6, f"no row kept the carriers as full AD: {ads}"


# ── Complex (delins) sibling — the key cross-type case ───────────────────
# Plant at 298:  G T A G T C A G T T
#                298     302     306
#   A: pure del AG @[300,302)      (anchor 299 'T' — already canonical)
#   X: delins CAG@[303,306) -> T   (anchor-substituted -> check_complex)
# A delins carrier aligns as (mismatch C->T)(D2 deleting AG@[304,306)):
# its D2's deleted bases 'AG' S3-match row A's expected 'AG' at a windowed
# position -> pre-fix the pure-del row counts delins carriers as ALT.
CPLX_PLANT_POS = 298
CPLX_PLANT = "GTAGTCAGTT"
CPLX_DEL = 300  # 0-based first deleted base of the pure del
CPLX_X = 303  # 0-based first ref base of the delins


def _complex_setup(tmp_path, n_del=6, n_x=5, n_wt=8):
    ref = _mk_ref(plants=((CPLX_PLANT_POS, CPLX_PLANT),))
    reads = _del_reads(ref, CPLX_DEL, 2, n_del, "pd")
    for i in range(n_x):  # delins carriers: M(..mismatch at 303)-D2-M
        s = CPLX_X - 35 - (i % 5)
        left = CPLX_X - s
        seq = ref[s:CPLX_X] + "T" + ref[CPLX_X + 3 : CPLX_X + 3 + (READ_LEN - left - 1)]
        reads.append(make_read(f"xc{i}", seq, s, ((0, left + 1), (2, 2), (0, READ_LEN - left - 1))))
    reads += _ref_reads(ref, CPLX_DEL, n_wt)
    return ref, reads


def _complex_vcf(tmp_path, ref):
    return _vcf(
        tmp_path,
        [
            (CPLX_DEL, ref[CPLX_DEL - 1 : CPLX_DEL + 2], ref[CPLX_DEL - 1]),
            (CPLX_X + 1, ref[CPLX_X : CPLX_X + 3], "T"),  # CAG>T, 1-based POS 304
        ],
    )


def test_complex_sibling_claims_windowed_carrier(tmp_path):
    """The key cross-type case: a delins tract-mate's carriers S3-match the
    pure-del row's windowed scan pre-fix. Post-fix the delins sibling claims
    them (full-ALT for the complex row beats a windowed match), and the
    pure-del row records them as partial evidence."""
    ref, reads = _complex_setup(tmp_path)
    res = _run(
        tmp_path, _complex_vcf(tmp_path, ref), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref)
    )
    r_del = res[CPLX_DEL + 1]
    # The writer re-synthesizes VCF-sourced delins coordinates (DEL-style
    # representation), so the delins row is found by exclusion, not by POS.
    r_x = next(r for pos, r in res.items() if pos != CPLX_DEL + 1)
    assert int(r_x["alt_count"]) == 5, f"delins row must keep its own carriers: {r_x['alt_count']}"
    assert (
        int(r_del["alt_count"]) == 6
    ), f"pure-del row must not absorb delins carriers, got ad={r_del['alt_count']}"
    assert (
        int(r_del["partial_alt"]) >= 5
    ), f"claimed delins carriers must surface as partial, got {r_del['partial_alt']}"


def _cluster_maf(tmp_path, ref, rows):
    """MAF twin of _vcf for deletion rows: Start = first deleted base
    (1-based), REF = deleted bases, ALT = '-'."""
    maf = tmp_path / "variants.maf"
    lines = [
        "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
        "Reference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode"
    ]
    for pos1, ref_al, alt_al in rows:
        if len(alt_al) == 1 and len(ref_al) > 1 and ref_al[0] == alt_al:
            # pure deletion (VCF anchor form: anchor base preserved)
            deleted = ref_al[1:]
            start = pos1 + 1
            lines.append(f"GENE\tchr1\t{start}\t{start + len(deleted) - 1}\t{deleted}\t-\tS")
        else:  # complex delins: identical representation in MAF (no '-')
            lines.append(f"GENE\tchr1\t{pos1}\t{pos1 + len(ref_al) - 1}\t{ref_al}\t{alt_al}\tS")
    maf.write_text("\n".join(lines) + "\n")
    return maf


def test_cluster_maf_and_vcf_paths_agree(tmp_path):
    """Guard (green now, green after): the cluster counts identically through
    the MAF input path and the VCF input path — whatever the assignment
    semantics, the two front doors must agree."""
    ref, rows, reads = _cluster_setup(tmp_path)
    bam = _bam(tmp_path, ref, reads)
    fasta = _fasta(tmp_path, ref)
    res_v = _run(tmp_path, _vcf(tmp_path, rows), bam, fasta, outname="out_v")
    res_m = _run(tmp_path, _cluster_maf(tmp_path, ref, rows), bam, fasta, outname="out_m")
    for pos, _, _ in rows:
        rv, rm = res_v[pos + 1], res_m[pos + 1]
        for col in ("total_count", "ref_count", "alt_count", "partial_alt", "any_alt"):
            assert rv[col] == rm[col], f"{col} differs at {pos}: vcf={rv[col]} maf={rm[col]}"


def test_complex_cluster_maf_and_vcf_paths_agree(tmp_path):
    """Guard: same MAF/VCF agreement for the delins-bearing cluster (the MAF
    complex representation carries no '-' alleles and its own coordinate
    conventions)."""
    ref, reads = _complex_setup(tmp_path)
    bam = _bam(tmp_path, ref, reads)
    fasta = _fasta(tmp_path, ref)
    rows = [
        (CPLX_DEL, ref[CPLX_DEL - 1 : CPLX_DEL + 2], ref[CPLX_DEL - 1]),
        (CPLX_X + 1, ref[CPLX_X : CPLX_X + 3], "T"),
    ]
    res_v = _run(tmp_path, _complex_vcf(tmp_path, ref), bam, fasta, outname="out_v")
    res_m = _run(tmp_path, _cluster_maf(tmp_path, ref, rows), bam, fasta, outname="out_m")
    # The writer renders VCF-sourced complex rows differently from
    # MAF-sourced ones (re-synthesized vs original coordinates), so pair
    # rows by sorted position and compare COUNTS only.
    assert len(res_v) == len(res_m) == 2
    for kv, km in zip(sorted(res_v), sorted(res_m), strict=True):
        for col in ("total_count", "ref_count", "alt_count", "partial_alt", "any_alt"):
            assert (
                res_v[kv][col] == res_m[km][col]
            ), f"{col} differs at {kv}/{km}: vcf={res_v[kv][col]} maf={res_m[km][col]}"


# ── An SNV between window-joined deletions stays out of the cluster ─────
# Two 3bp deletions whose scan windows overlap while their spans do not, with
# an annotated SNV in the gap between them. The SNV must not join their group
# (it would lose every deletion carrier from its REF count): a candidate joins
# on a member's own span or window, not on the box from the group's first
# member to its last.
GAP_D1, GAP_SNV, GAP_D2 = 151, 157, 163  # 0-based: first deleted base / SNV base


def _gap_setup(n_d1=6, n_d2=5, n_snv=4, n_wt=8):
    ref = _mk_ref(seed=29, plants=((150, "TCGG"), (157, "A"), (162, "GACA")))
    rows = [(151, "TCGG", "T"), (158, "A", "C"), (163, "GACA", "G")]
    snv_reads = []
    for i in range(n_snv):
        s = GAP_SNV - 50 + (i % 5)
        seq = ref[s:GAP_SNV] + "C" + ref[GAP_SNV + 1 : s + READ_LEN]
        snv_reads.append(make_read(f"s{i}", seq, s, ((0, READ_LEN),)))
    reads = (
        _del_reads(ref, GAP_D1, 3, n_d1, "p")
        + _del_reads(ref, GAP_D2, 3, n_d2, "q")
        + snv_reads
        + _ref_reads(ref, GAP_SNV, n_wt)
    )
    return ref, rows, reads


@pytest.mark.xfail(
    strict=True, reason="the group's span box pulls in the SNV between the deletions"
)
def test_snv_between_window_joined_deletions_stays_ungrouped(tmp_path):
    ref, rows, _ = _gap_setup()
    fasta = _fasta(tmp_path, ref)
    prepared = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", p - 1, r, a, "SNP") for p, r, a in rows],
        str(fasta),
        5,
        False,
        1,
        True,
    )
    got = [(pv.multi_allelic_group, pv.gbcms_status_reason) for pv in prepared]
    assert got == [(1, "TRACT_CLUSTER"), (None, ""), (1, "TRACT_CLUSTER")]


@pytest.mark.xfail(strict=True, reason="deletion carriers leave the SNV's REF count")
def test_snv_between_window_joined_deletions_keeps_its_ref_reads(tmp_path):
    """Deletion carriers show the reference base at the SNV: they are REF for it."""
    ref, rows, reads = _gap_setup()
    res = _run(tmp_path, _vcf(tmp_path, rows), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    snv = res[GAP_SNV + 1]
    assert (int(snv["ref_count"]), int(snv["alt_count"]), int(snv["partial_alt"])) == (19, 4, 0)
    assert int(snv["ref_count"]) == int(snv["ref_count_forward"]) + int(snv["ref_count_reverse"])
    assert int(snv["alt_count"]) == int(snv["alt_count_forward"]) + int(snv["alt_count_reverse"])
    assert int(snv["total_count_fragment"]) >= int(snv["ref_count_fragment"]) + int(
        snv["alt_count_fragment"]
    )
