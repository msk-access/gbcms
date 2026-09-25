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
(xfail-strict) before the implementation; flipped green with it.
"""

import glob

from helpers import read_maf_output
from rna_fixtures import (
    E1,
    E2,
    E3,
    mk_ref,
    run_rna,
    spliced,
    through,
    write_bam,
    write_fasta,
    write_gtf,
    write_maf,
    write_vcf,
)
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()


def _flags(row):
    return [f for f in row["asjd_diagnostic"].split(";") if f]


# ═════════════ Donor SNV with allele-specific intron retention ══════════
DONOR_SNV = 300  # 0-based, the G of intron-1's GT donor


def _retention_setup(ref):
    rows = [(DONOR_SNV + 1, "G", "A")]
    reads = (
        through(ref, DONOR_SNV, "A", 20, "ret_alt")  # mutant allele retains intron 1
        + through(ref, DONOR_SNV, "G", 2, "ret_ref")
        + spliced(ref, E1[1], E2[0], 60, "wt")  # normal E1->E2 splicing
        + spliced(ref, E1[1], E3[0], 5, "skip")  # minority exon-2 skip
    )
    return rows, reads


def test_donor_snv_retention_marker(tmp_path):
    """The allele-classified reads are all intron-retaining (junction-free)
    while 65 spliced fragments skip the locus: RETENTION_DOMINANT(65)."""
    ref = mk_ref()
    rows, reads = _retention_setup(ref)
    res = run_rna(
        tmp_path,
        write_vcf(tmp_path, rows),
        write_bam(tmp_path, ref, reads),
        write_fasta(tmp_path, ref),
        write_gtf(tmp_path),
    )
    (r,) = res
    assert int(r["alt_count"]) == 20
    assert "RETENTION_DOMINANT(65)" in _flags(r), r["asjd_diagnostic"]


def test_minority_novel_junction_below_alt_is_silent(tmp_path):
    """Guard (green now, green after): 5 fragments on an anchored novel skip
    junction do not outnumber the 20 ALT (retained) fragments — the mutant
    allele's dominant outcome is retention, so no novel-junction marker."""
    ref = mk_ref()
    rows, reads = _retention_setup(ref)
    res = run_rna(
        tmp_path,
        write_vcf(tmp_path, rows),
        write_bam(tmp_path, ref, reads),
        write_fasta(tmp_path, ref),
        write_gtf(tmp_path),
    )
    (r,) = res
    assert not any(f.startswith("NOVEL_JUNC_AT_SPLICE_LOSS") for f in _flags(r)), r[
        "asjd_diagnostic"
    ]


# ═════════════ Acceptor SNV expressed partly as exon skipping ═══════════
ACCEPTOR_SNV = 499  # 0-based, the G of intron-1's AG acceptor


def test_acceptor_snv_skip_and_retention_markers(tmp_path):
    """Acceptor SNV: 5 ALT + 5 REF retained reads, 40 normally spliced and
    12 exon-2-skip fragments spanning the locus. The anchored novel skip
    junction (12) outnumbers ALT (5) -> NOVEL marker; spliced fragments
    (52) dominate the junction-free classified population -> RETENTION."""
    ref = mk_ref()
    rows = [(ACCEPTOR_SNV + 1, "G", "A")]
    reads = (
        through(ref, ACCEPTOR_SNV, "A", 5, "ra")
        + through(ref, ACCEPTOR_SNV, "G", 5, "rr")
        + spliced(ref, E1[1], E2[0], 40, "wt")
        + spliced(ref, E1[1], E3[0], 12, "skip")
    )
    res = run_rna(
        tmp_path,
        write_vcf(tmp_path, rows),
        write_bam(tmp_path, ref, reads),
        write_fasta(tmp_path, ref),
        write_gtf(tmp_path),
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
        spliced(ref, E1[1], E3[0], 15, "skip")  # mutant transcripts: E1->E3
        + spliced(ref, E1[1], E2[0], 30, "wt")  # WT E1->E2 (does not span the deletion)
        + through(ref, 550, ref[550], 30, "ex2")  # WT exon-2 coverage
    )


def test_exon_deletion_expressed_as_skip(tmp_path):
    """The deletion leaves ad=0 honestly (its carriers splice over the
    locus), but the excluded population's anchored novel E1->E3 junction
    (15 fragments) is the splicing consequence: surface it."""
    ref = mk_ref()
    res = run_rna(
        tmp_path,
        write_vcf(tmp_path, _exon_del_rows(ref)),
        write_bam(tmp_path, ref, _exon_del_reads(ref)),
        write_fasta(tmp_path, ref),
        write_gtf(tmp_path),
    )
    (r,) = res
    assert int(r["alt_count"]) == 0
    assert f"NOVEL_JUNC_AT_SPLICE_LOSS(15@{E1[1]}-{E3[0]})" in _flags(r), r["asjd_diagnostic"]


def test_exon_deletion_markers_agree_across_input_paths(tmp_path):
    """Guard (green now, green after): the VCF and MAF front doors produce
    identical counts and ASJD diagnostics for the exon deletion."""
    ref = mk_ref()
    bam = write_bam(tmp_path, ref, _exon_del_reads(ref))
    fasta, gtf = write_fasta(tmp_path, ref), write_gtf(tmp_path)
    rows = _exon_del_rows(ref)
    (rv,) = run_rna(tmp_path, write_vcf(tmp_path, rows), bam, fasta, gtf, outname="out_v")
    (rm,) = run_rna(tmp_path, write_maf(tmp_path, rows), bam, fasta, gtf, outname="out_m")
    for col in ("ref_count", "alt_count", "total_count", "partial_alt", "asjd_diagnostic"):
        assert rv[col] == rm[col], (col, rv[col], rm[col])


def test_deletion_written_as_splice_is_not_a_novel_junction(tmp_path):
    """Guard (green now, green after): a 30bp deletion starting at the donor,
    written by the aligner as N(30) on its carriers, produces a junction
    anchored at an annotated site — but it is the deletion itself, not a
    splicing consequence. SPLICE_SKIP_DOMINANT already explains it; the
    novel-junction marker must stay silent."""
    ref = mk_ref()
    rows = [(E1[1], ref[E1[1] - 1 : E1[1] + 30], ref[E1[1] - 1])]
    reads = spliced(ref, E1[1], E1[1] + 30, 12, "delN") + spliced(ref, E1[1], E2[0], 20, "wt")
    res = run_rna(
        tmp_path,
        write_vcf(tmp_path, rows),
        write_bam(tmp_path, ref, reads),
        write_fasta(tmp_path, ref),
        write_gtf(tmp_path),
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
    ref = mk_ref()
    pos = E1[1] - 1
    alt = "A" if ref[pos] != "A" else "C"
    wt = spliced(ref, E1[1], E2[0], 30, "wt")
    carriers = []
    for a in spliced(ref, E1[1], E2[0], 10, "car"):
        seq = list(a.query_sequence)
        seq[pos - a.reference_start] = alt
        a.query_sequence = "".join(seq)
        a.query_qualities = [30] * len(seq)
        carriers.append(a)
    rows = [(pos + 1, ref[pos], alt)]
    res = run_rna(
        tmp_path,
        write_vcf(tmp_path, rows),
        write_bam(tmp_path, ref, wt + carriers),
        write_fasta(tmp_path, ref),
        write_gtf(tmp_path),
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
    ref = mk_ref()
    pos = 550
    alt = "A" if ref[pos] != "A" else "C"
    reads = through(ref, pos, alt, 8, "alt") + through(ref, pos, ref[pos], 8, "ref")
    reads += spliced(ref, E1[1], E2[0], 20, "wt") + spliced(ref, E1[1], E3[0], 20, "skip")
    rows = [(pos + 1, ref[pos], alt)]
    res = run_rna(
        tmp_path,
        write_vcf(tmp_path, rows),
        write_bam(tmp_path, ref, reads),
        write_fasta(tmp_path, ref),
        write_gtf(tmp_path),
    )
    (r,) = res
    flags = _flags(r)
    assert not any(
        f.startswith(("RETENTION_DOMINANT", "NOVEL_JUNC_AT_SPLICE_LOSS")) for f in flags
    ), flags


def test_markers_absent_without_gtf(tmp_path):
    """Guard: the markers are GTF-gated like every ASJD field — without
    --gtf the ASJD columns are absent altogether (column gating rule)."""
    ref = mk_ref()
    rows, reads = _retention_setup(ref)
    outdir = tmp_path / "nogtf"
    outdir.mkdir()
    result = runner.invoke(
        app,
        [
            "rna",
            "-v",
            str(write_vcf(tmp_path, rows)),
            "-b",
            f"S:{write_bam(tmp_path, ref, reads)}",
            "-f",
            str(write_fasta(tmp_path, ref)),
            "-o",
            str(outdir),
            "--format",
            "maf",
        ],
    )
    assert result.exit_code == 0, result.output
    (r,) = list(read_maf_output(glob.glob(str(outdir / "*.maf"))[0]))
    assert "asjd_diagnostic" not in r


def test_markers_need_junction_evidence_floors(tmp_path):
    """Guard: 3 spliced-over fragments vs 1 classified, and a 2-fragment
    anchored novel junction vs 0 ALT, are too little evidence to speak on.
    Both markers stay silent below ASJD's own junction floors (10 REF-side
    fragments for the spliced population, 5 ALT-side for a novel junction)."""
    ref = mk_ref()
    rows = [(DONOR_SNV + 1, "G", "A")]
    reads = (
        through(ref, DONOR_SNV, "G", 1, "ret")
        + spliced(ref, E1[1], E2[0], 3, "wt")
        + spliced(ref, E1[1], E3[0], 2, "skip")
    )
    res = run_rna(
        tmp_path,
        write_vcf(tmp_path, rows),
        write_bam(tmp_path, ref, reads),
        write_fasta(tmp_path, ref),
        write_gtf(tmp_path),
    )
    (r,) = res
    flags = _flags(r)
    assert not any(
        f.startswith(("RETENTION_DOMINANT", "NOVEL_JUNC_AT_SPLICE_LOSS")) for f in flags
    ), flags


def _pair(reads_r1, reads_r2):
    """Give mate pairs a shared QNAME and proper-pair flags. R1 is
    reverse-strand and R2 forward: both sense for a '+' gene under dUTP."""
    for i, (a, b) in enumerate(zip(reads_r1, reads_r2, strict=True)):
        name = f"{a.query_name}_frag{i}"
        a.query_name, b.query_name = name, name
        a.flag = 1 | 2 | 16 | 64  # paired, proper, reverse, read1
        b.flag = 1 | 2 | 32 | 128  # paired, proper, mate reverse, read2
    return reads_r1 + reads_r2


def test_marker_counts_are_per_fragment(tmp_path):
    """Both mates of every fragment carry the same evidence: 12 fragments
    spliced over the donor (24 reads) and 8 retained ALT fragments (16
    reads). The marker count is fragments, not reads: RETENTION_DOMINANT(12)."""
    ref = mk_ref()
    rows = [(DONOR_SNV + 1, "G", "A")]
    wt = _pair(spliced(ref, E1[1], E2[0], 12, "wt1"), spliced(ref, E1[1], E2[0], 12, "wt2"))
    ret = _pair(through(ref, DONOR_SNV, "A", 8, "ra1"), through(ref, DONOR_SNV, "A", 8, "ra2"))
    res = run_rna(
        tmp_path,
        write_vcf(tmp_path, rows),
        write_bam(tmp_path, ref, wt + ret),
        write_fasta(tmp_path, ref),
        write_gtf(tmp_path),
    )
    (r,) = res
    assert int(r["alt_count"]) == 16, "read-level ALT counts both mates"
    assert int(r["alt_count_fragment"]) == 8
    assert "RETENTION_DOMINANT(12)" in _flags(r), r["asjd_diagnostic"]
