"""Target contract for ASJD-2: splice-disruption signals the junction tallies
cannot see (issue #97).

ASJD compares junction usage between REF- and ALT-classified reads. At a
splice-site variant the informative reads are often exactly the ones the
splice-aware evidence rule excludes: spliced reads whose CIGAR N spans the
variant observe nothing there. Two markers surface what that excluded
population shows, through the existing ``asjd_diagnostic`` string (no new
columns; RNA + GTF only):

- ``RETENTION_DOMINANT(n)``: the variant touches an annotated splice site,
  spliced-over fragments (n) outnumber the allele-classified ones, and the
  classified reads are predominantly junction-free — the counts at this
  locus describe the intron-retaining population only, so ``vaf`` is the VAF
  within that population, not allelic balance.
- ``NOVEL_JUNC_AT_SPLICE_LOSS(n@start-end)``: the variant touches an
  annotated splice site and the excluded population's top unannotated
  junction — anchored to an annotated site, and not the deletion itself
  written as a splice — is carried by more fragments (n) than confirm ALT.

Geometries mirror three signed-out RNA cases measured locally: a donor SNV
with allele-specific retention, an acceptor SNV expressed partly as exon
skipping, and an exon-removing deletion expressed as an exon skip. Both
variant input paths (VCF and MAF) are exercised. Committed red
(xfail-strict) before the implementation.
"""

import glob
import random

import pysam
import pytest
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

XFAIL = pytest.mark.xfail(strict=True, reason="ASJD-2 markers: implementation pending")

BAM_CONTIG = "1"
READ_LEN = 100

# ── Gene model (0-based, half-open) ─────────────────────────────────────
# E1 [100,300)  intron1 [300,500)  E2 [500,600)  intron2 [600,800)  E3 [800,1000)
# Canonical motifs planted: GT at 300 / AG at 498; GT at 600 / AG at 798.
E1, E2, E3 = (100, 300), (500, 600), (800, 1000)
PLANTS = ((300, "GT"), (498, "AG"), (600, "GT"), (798, "AG"))


def _mk_ref(n=1200, seed=23):
    rng = random.Random(seed)
    ref = [rng.choice("ACGT") for _ in range(n)]
    for pos, motif in PLANTS:
        ref[pos : pos + len(motif)] = list(motif)
    return "".join(ref)


def _gtf(tmp_path):
    gtf = tmp_path / "gene.gtf"
    gtf.write_text(
        "".join(
            f'chr1\tTEST\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
            for s, e in (E1, E2, E3)
        )
    )
    return gtf


def _bam(tmp_path, ref, reads, name="rna.bam"):
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
        + "".join(f"chr1\t{p}\t.\t{r}\t{a}\t.\t.\t.\n" for p, r, a in rows)
    )
    return vcf


def _maf(tmp_path, rows):
    """MAF twin of _vcf: anchor-preserved deletions become Start = first
    deleted base, REF = deleted bases, ALT = '-'; SNVs are identical."""
    maf = tmp_path / "variants.maf"
    lines = [
        "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
        "Reference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode"
    ]
    for p, r, a in rows:
        if len(a) == 1 and len(r) > 1 and r[0] == a:
            d = r[1:]
            lines.append(f"G1\tchr1\t{p + 1}\t{p + len(d)}\t{d}\t-\tS")
        else:
            lines.append(f"G1\tchr1\t{p}\t{p + len(r) - 1}\t{r}\t{a}\tS")
    maf.write_text("\n".join(lines) + "\n")
    return maf


def _run(tmp_path, variants, bam, fasta, gtf, outname="out"):
    outdir = tmp_path / outname
    outdir.mkdir(exist_ok=True)
    result = runner.invoke(
        app,
        [
            "rna",
            "-v",
            str(variants),
            "-b",
            f"S:{bam}",
            "-f",
            str(fasta),
            "-o",
            str(outdir),
            "--format",
            "maf",
            "--gtf",
            str(gtf),
        ],
    )
    assert result.exit_code == 0, result.output
    rows = list(read_maf_output(glob.glob(str(outdir / "*.maf"))[0]))
    for r in rows:
        assert int(r["total_count"]) >= int(r["ref_count"]) + int(r["alt_count"])
        assert int(r["total_count_fragment"]) >= int(r["ref_count_fragment"]) + int(
            r["alt_count_fragment"]
        )
        assert int(r["ref_count"]) == int(r["ref_count_forward"]) + int(r["ref_count_reverse"])
        assert int(r["alt_count"]) == int(r["alt_count_forward"]) + int(r["alt_count_reverse"])
        assert int(r["any_alt"]) == int(r["alt_count"]) + int(r["partial_alt"])
    return rows


def _flags(row):
    return [f for f in row["asjd_diagnostic"].split(";") if f]


# ── Read builders ────────────────────────────────────────────────────────
# Reads are reverse-strand single-end (flag 16): sense for a '+' gene under
# the default dUTP ('reverse') protocol, so the RNA defaults
# (--enforce-strandedness, BAQ with GTF boundary suppression) apply as in
# production rather than filtering the synthetic reads as antisense.
SENSE = 16


def _spliced(ref, donor, acceptor, n, prefix):
    """Reads whose left M block ends at `donor` (exclusive), then N over
    [donor, acceptor), then M from `acceptor`."""
    out = []
    for i in range(n):
        s = donor - 40 - (i % 5)
        left = donor - s
        right = READ_LEN - left
        seq = ref[s:donor] + ref[acceptor : acceptor + right]
        out.append(
            make_read(
                f"{prefix}{i}", seq, s, ((0, left), (3, acceptor - donor), (0, right)), flag=SENSE
            )
        )
    return out


def _through(ref, pos, base, n, prefix):
    """Unspliced (M-only) reads covering `pos`, carrying `base` there."""
    out = []
    for i in range(n):
        s = pos - 50 + (i % 5)
        seq = ref[s:pos] + base + ref[pos + 1 : s + READ_LEN]
        out.append(make_read(f"{prefix}{i}", seq, s, ((0, READ_LEN),), flag=SENSE))
    return out


# ═════════════ Donor SNV with allele-specific intron retention ══════════
DONOR_SNV = 300  # 0-based, the G of intron-1's GT donor


def _retention_setup(ref):
    rows = [(DONOR_SNV + 1, "G", "A")]
    reads = (
        _through(ref, DONOR_SNV, "A", 20, "ret_alt")  # mutant allele retains intron 1
        + _through(ref, DONOR_SNV, "G", 2, "ret_ref")
        + _spliced(ref, E1[1], E2[0], 60, "wt")  # normal E1->E2 splicing
        + _spliced(ref, E1[1], E3[0], 5, "skip")  # minority exon-2 skip
    )
    return rows, reads


@XFAIL
def test_donor_snv_retention_marker(tmp_path):
    """The allele-classified reads are all intron-retaining (junction-free)
    while 65 spliced fragments skip the locus: RETENTION_DOMINANT(65)."""
    ref = _mk_ref()
    rows, reads = _retention_setup(ref)
    res = _run(
        tmp_path,
        _vcf(tmp_path, rows),
        _bam(tmp_path, ref, reads),
        _fasta(tmp_path, ref),
        _gtf(tmp_path),
    )
    (r,) = res
    assert int(r["alt_count"]) == 20
    assert "RETENTION_DOMINANT(65)" in _flags(r), r["asjd_diagnostic"]


def test_minority_novel_junction_below_alt_is_silent(tmp_path):
    """Guard (green now, green after): 5 fragments on an anchored novel skip
    junction do not outnumber the 20 ALT (retained) fragments — the mutant
    allele's dominant outcome is retention, so no novel-junction marker."""
    ref = _mk_ref()
    rows, reads = _retention_setup(ref)
    res = _run(
        tmp_path,
        _vcf(tmp_path, rows),
        _bam(tmp_path, ref, reads),
        _fasta(tmp_path, ref),
        _gtf(tmp_path),
    )
    (r,) = res
    assert not any(f.startswith("NOVEL_JUNC_AT_SPLICE_LOSS") for f in _flags(r)), r[
        "asjd_diagnostic"
    ]


# ═════════════ Acceptor SNV expressed partly as exon skipping ═══════════
ACCEPTOR_SNV = 499  # 0-based, the G of intron-1's AG acceptor


@XFAIL
def test_acceptor_snv_skip_and_retention_markers(tmp_path):
    """Acceptor SNV: 5 ALT + 5 REF retained reads, 40 normally spliced and
    12 exon-2-skip fragments spanning the locus. The anchored novel skip
    junction (12) outnumbers ALT (5) -> NOVEL marker; spliced fragments
    (52) dominate the junction-free classified population -> RETENTION."""
    ref = _mk_ref()
    rows = [(ACCEPTOR_SNV + 1, "G", "A")]
    reads = (
        _through(ref, ACCEPTOR_SNV, "A", 5, "ra")
        + _through(ref, ACCEPTOR_SNV, "G", 5, "rr")
        + _spliced(ref, E1[1], E2[0], 40, "wt")
        + _spliced(ref, E1[1], E3[0], 12, "skip")
    )
    res = _run(
        tmp_path,
        _vcf(tmp_path, rows),
        _bam(tmp_path, ref, reads),
        _fasta(tmp_path, ref),
        _gtf(tmp_path),
    )
    (r,) = res
    flags = _flags(r)
    assert f"NOVEL_JUNC_AT_SPLICE_LOSS(12@{E1[1]}-{E3[0]})" in flags, r["asjd_diagnostic"]
    assert "RETENTION_DOMINANT(52)" in flags, r["asjd_diagnostic"]


# ═════════════ Exon-removing deletion expressed as exon skipping ════════
DEL_START, DEL_END = 480, 620  # deletes E2 entirely (both of its splice sites)


def _exon_del_rows(ref):
    return [(DEL_START, ref[DEL_START - 1 : DEL_END], ref[DEL_START - 1])]


def _exon_del_reads(ref):
    return (
        _spliced(ref, E1[1], E3[0], 15, "skip")  # mutant transcripts: E1->E3
        + _spliced(ref, E1[1], E2[0], 30, "wt")  # WT E1->E2 (does not span the deletion)
        + _through(ref, 550, ref[550], 30, "ex2")  # WT exon-2 coverage
    )


@XFAIL
def test_exon_deletion_expressed_as_skip(tmp_path):
    """The deletion leaves ad=0 honestly (its carriers splice over the
    locus), but the excluded population's anchored novel E1->E3 junction
    (15 fragments) is the splicing consequence: surface it."""
    ref = _mk_ref()
    res = _run(
        tmp_path,
        _vcf(tmp_path, _exon_del_rows(ref)),
        _bam(tmp_path, ref, _exon_del_reads(ref)),
        _fasta(tmp_path, ref),
        _gtf(tmp_path),
    )
    (r,) = res
    assert int(r["alt_count"]) == 0
    assert f"NOVEL_JUNC_AT_SPLICE_LOSS(15@{E1[1]}-{E3[0]})" in _flags(r), r["asjd_diagnostic"]


def test_exon_deletion_markers_agree_across_input_paths(tmp_path):
    """Guard (green now, green after): the VCF and MAF front doors produce
    identical counts and ASJD diagnostics for the exon deletion."""
    ref = _mk_ref()
    bam = _bam(tmp_path, ref, _exon_del_reads(ref))
    fasta, gtf = _fasta(tmp_path, ref), _gtf(tmp_path)
    rows = _exon_del_rows(ref)
    (rv,) = _run(tmp_path, _vcf(tmp_path, rows), bam, fasta, gtf, outname="out_v")
    (rm,) = _run(tmp_path, _maf(tmp_path, rows), bam, fasta, gtf, outname="out_m")
    for col in ("ref_count", "alt_count", "total_count", "partial_alt", "asjd_diagnostic"):
        assert rv[col] == rm[col], (col, rv[col], rm[col])


def test_deletion_written_as_splice_is_not_a_novel_junction(tmp_path):
    """Guard (green now, green after): a 30bp deletion starting at the donor,
    written by the aligner as N(30) on its carriers, produces a junction
    anchored at an annotated site — but it is the deletion itself, not a
    splicing consequence. SPLICE_SKIP_DOMINANT already explains it; the
    novel-junction marker must stay silent."""
    ref = _mk_ref()
    rows = [(E1[1], ref[E1[1] - 1 : E1[1] + 30], ref[E1[1] - 1])]
    reads = _spliced(ref, E1[1], E1[1] + 30, 12, "delN") + _spliced(ref, E1[1], E2[0], 20, "wt")
    res = _run(
        tmp_path,
        _vcf(tmp_path, rows),
        _bam(tmp_path, ref, reads),
        _fasta(tmp_path, ref),
        _gtf(tmp_path),
    )
    (r,) = res
    assert not any(f.startswith("NOVEL_JUNC_AT_SPLICE_LOSS") for f in _flags(r)), r[
        "asjd_diagnostic"
    ]
    assert "SPLICE_SKIP_DOMINANT(" in r["gbcms_diagnostic"], r["gbcms_diagnostic"]


# ═════════════ Guards away from splice loss ═════════════════════════════
def test_exon_side_boundary_snv_with_normal_splicing_is_silent(tmp_path):
    """Guard: an SNV on the last exon base is aligned in every spliced read,
    so nothing is excluded — neither marker may fire."""
    ref = _mk_ref()
    pos = E1[1] - 1
    alt = "A" if ref[pos] != "A" else "C"
    wt = _spliced(ref, E1[1], E2[0], 30, "wt")
    carriers = []
    for a in _spliced(ref, E1[1], E2[0], 10, "car"):
        seq = list(a.query_sequence)
        seq[pos - a.reference_start] = alt
        a.query_sequence = "".join(seq)
        a.query_qualities = [30] * len(seq)
        carriers.append(a)
    rows = [(pos + 1, ref[pos], alt)]
    res = _run(
        tmp_path,
        _vcf(tmp_path, rows),
        _bam(tmp_path, ref, wt + carriers),
        _fasta(tmp_path, ref),
        _gtf(tmp_path),
    )
    (r,) = res
    assert int(r["alt_count"]) == 10
    flags = _flags(r)
    assert not any(
        f.startswith(("RETENTION_DOMINANT", "NOVEL_JUNC_AT_SPLICE_LOSS")) for f in flags
    ), flags


def test_mid_exon_snv_is_silent(tmp_path):
    """Guard: a mid-exon SNV (50bp from any splice site) is out of scope
    for both markers regardless of nearby junction traffic."""
    ref = _mk_ref()
    pos = 550
    alt = "A" if ref[pos] != "A" else "C"
    reads = _through(ref, pos, alt, 8, "alt") + _through(ref, pos, ref[pos], 8, "ref")
    reads += _spliced(ref, E1[1], E2[0], 20, "wt") + _spliced(ref, E1[1], E3[0], 20, "skip")
    rows = [(pos + 1, ref[pos], alt)]
    res = _run(
        tmp_path,
        _vcf(tmp_path, rows),
        _bam(tmp_path, ref, reads),
        _fasta(tmp_path, ref),
        _gtf(tmp_path),
    )
    (r,) = res
    flags = _flags(r)
    assert not any(
        f.startswith(("RETENTION_DOMINANT", "NOVEL_JUNC_AT_SPLICE_LOSS")) for f in flags
    ), flags


def test_markers_absent_without_gtf(tmp_path):
    """Guard: the markers are GTF-gated like every ASJD field — without
    --gtf the ASJD columns are absent altogether (column gating rule)."""
    ref = _mk_ref()
    rows, reads = _retention_setup(ref)
    outdir = tmp_path / "nogtf"
    outdir.mkdir()
    result = runner.invoke(
        app,
        [
            "rna",
            "-v",
            str(_vcf(tmp_path, rows)),
            "-b",
            f"S:{_bam(tmp_path, ref, reads)}",
            "-f",
            str(_fasta(tmp_path, ref)),
            "-o",
            str(outdir),
            "--format",
            "maf",
        ],
    )
    assert result.exit_code == 0, result.output
    (r,) = list(read_maf_output(glob.glob(str(outdir / "*.maf"))[0]))
    assert "asjd_diagnostic" not in r
