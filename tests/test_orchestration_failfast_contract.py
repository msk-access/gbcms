"""Target contract for orchestration fail-fast and silent-failure fixes.

Each xfail(strict=True) test documents a verified defect in the Python
orchestration layer and states the behavior the fix must produce; they fail
on the current code for the reasons given and must flip green with the fix:

  - a rejected (FAIL) variant must not crash ``--mfsd`` or RNA ``--gtf``
    runs (the zero-count stub must satisfy every column the writers read)
  - ``--mfsd-parquet`` must work when any variant was rejected (stubs are
    not BaseCounts and must not reach the Rust parquet writer)
  - ``--bam-list`` honors fail-fast: a missing entry or an unreadable list
    file exits non-zero unless ``--lenient-bam`` was given
  - duplicate sample names are an error, never a silent overwrite
  - input-MAF columns colliding with gbcms output names are warned about
  - ``merge`` must not coerce float-formatted count strings to 0
  - ragged MAF rows fail loudly in the batch readers
  - a failed sample's error report names the exception type

Guards pin adjacent behavior that must NOT change.
"""

import glob
import gzip  # noqa: F401  (documents that GTFs here are plain text on purpose)

import pysam
import pytest
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

XFAIL = pytest.mark.xfail(
    strict=True, reason="orchestration silent-failure: fix pending (issue #92 class 1)"
)

BAM_CONTIG = "1"  # normalized name; the pipeline strips 'chr' before fetch
SNV_POS0 = 100  # 0-based PASS SNV position
FAIL_POS0 = 300  # 0-based position of the deliberately REF-mismatched variant


# ── shared fixtures ──────────────────────────────────────────────────────
def _build_reference(tmp_path):
    """600bp all-A chr1 FASTA with pinned non-A bases at both variant sites."""
    seq = ["A"] * 600
    seq[SNV_POS0] = "C"
    seq[FAIL_POS0] = "G"
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">chr1\n" + "".join(seq) + "\n")
    pysam.faidx(str(fasta))
    return fasta, "".join(seq)


def _build_bam(tmp_path, ref, name="sample.bam"):
    """Reads over both loci: 4 REF + 2 ALT at the SNV, 4 REF at the FAIL site."""
    reads = []
    for i in range(4):
        s = SNV_POS0 - 50 + i
        reads.append(make_read(f"r{i}", ref[s : s + 100], s, ((0, 100),)))
    for i in range(2):
        s = SNV_POS0 - 40 + i
        seq = ref[s:SNV_POS0] + "T" + ref[SNV_POS0 + 1 : s + 100]
        reads.append(make_read(f"a{i}", seq, s, ((0, 100),)))
    for i in range(4):
        s = FAIL_POS0 - 50 + i
        reads.append(make_read(f"f{i}", ref[s : s + 100], s, ((0, 100),)))
    unsorted = tmp_path / f"unsorted_{name}"
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 600, "SN": BAM_CONTIG}]}
    with pysam.AlignmentFile(unsorted, "wb", header=header) as fh:
        for r in sorted(reads, key=lambda a: a.reference_start):
            fh.write(r)
    bam = tmp_path / name
    pysam.sort("-o", str(bam), str(unsorted))
    pysam.index(str(bam))
    return bam


def _build_vcf(tmp_path, with_fail=True):
    """One PASS SNV; optionally one variant whose REF mismatches the FASTA
    (rejected at preparation with gbcms_status=FAIL)."""
    rows = [f"chr1\t{SNV_POS0 + 1}\t.\tC\tT\t.\t.\t."]
    if with_fail:
        rows.append(f"chr1\t{FAIL_POS0 + 1}\t.\tTTT\tT\t.\t.\t.")  # FASTA has G here
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=600>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" + "\n".join(rows) + "\n"
    )
    return vcf


def _build_gtf(tmp_path):
    """Minimal one-transcript, two-exon GTF covering both loci."""
    gtf = tmp_path / "tiny.gtf"
    gtf.write_text(
        'chr1\tTEST\texon\t50\t250\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        'chr1\tTEST\texon\t280\t550\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
    )
    return gtf


def _invoke(args):
    return runner.invoke(app, args)


def _dna_args(vcf, bam, fasta, outdir, *extra):
    return ["dna", "-v", str(vcf), "-b", str(bam), "-f", str(fasta), "-o", str(outdir), *extra]


# ═════════════ rejected-variant stub completeness ════════════════════════
def test_fail_variant_without_gated_modes_is_fine(tmp_path):
    """Guard: a FAIL variant in a plain DNA run already works — the row is
    emitted with zero counts and gbcms_status=FAIL. Must not change."""
    fasta, ref = _build_reference(tmp_path)
    bam = _build_bam(tmp_path, ref)
    outdir = tmp_path / "out"
    outdir.mkdir()
    result = _invoke(_dna_args(_build_vcf(tmp_path), bam, fasta, outdir, "--format", "maf"))
    assert result.exit_code == 0, result.output
    rows = list(read_maf_output(glob.glob(str(outdir / "*.maf"))[0]))
    assert len(rows) == 2
    statuses = sorted(r["gbcms_status"] for r in rows)
    assert statuses[0].startswith("FAIL") and statuses[1] == "PASS"


@XFAIL
def test_mfsd_run_survives_a_rejected_variant(tmp_path):
    """--mfsd + one FAIL variant: the zero-count stub lacks the mFSD q-value
    and sub/mono-nucleosomal fields the MAF writer reads, so the sample
    crashes and the run exits non-zero. The fix completes the stub: exit 0,
    both rows written, FAIL row present."""
    fasta, ref = _build_reference(tmp_path)
    bam = _build_bam(tmp_path, ref)
    outdir = tmp_path / "out"
    outdir.mkdir()
    result = _invoke(
        _dna_args(_build_vcf(tmp_path), bam, fasta, outdir, "--format", "maf", "--mfsd")
    )
    assert result.exit_code == 0, result.output
    rows = list(read_maf_output(glob.glob(str(outdir / "*.maf"))[0]))
    assert len(rows) == 2


@XFAIL
def test_rna_gtf_run_survives_a_rejected_variant(tmp_path):
    """RNA + --gtf + one FAIL variant: the stub lacks exon_boundary_dist,
    transcript_* and asjd_* fields read by the RNA/GTF columns → crash.
    Fixed: exit 0 with both rows written."""
    fasta, ref = _build_reference(tmp_path)
    bam = _build_bam(tmp_path, ref)
    outdir = tmp_path / "out"
    outdir.mkdir()
    result = _invoke(
        [
            "rna",
            "-v",
            str(_build_vcf(tmp_path)),
            "-b",
            str(bam),
            "-f",
            str(fasta),
            "-o",
            str(outdir),
            "--format",
            "maf",
            "--gtf",
            str(_build_gtf(tmp_path)),
        ]
    )
    assert result.exit_code == 0, result.output
    rows = list(read_maf_output(glob.glob(str(outdir / "*.maf"))[0]))
    assert len(rows) == 2


@XFAIL
def test_mfsd_parquet_survives_a_rejected_variant(tmp_path):
    """--mfsd-parquet + one FAIL variant: the merged counts list contains a
    Python stub that PyO3 cannot cast as BaseCounts → the parquet write (or
    the row write before it) fails. Fixed: exit 0, parquet written with the
    PASS rows only, exclusion logged."""
    fasta, ref = _build_reference(tmp_path)
    bam = _build_bam(tmp_path, ref)
    outdir = tmp_path / "out"
    outdir.mkdir()
    result = _invoke(
        _dna_args(
            _build_vcf(tmp_path),
            bam,
            fasta,
            outdir,
            "--format",
            "maf",
            "--mfsd",
            "--mfsd-parquet",
        )
    )
    assert result.exit_code == 0, result.output
    parquets = glob.glob(str(outdir / "*.fsd.parquet"))
    assert len(parquets) == 1, "mFSD parquet must be written despite the rejected variant"
    import polars as pl

    assert pl.read_parquet(parquets[0]).height == 1  # PASS variant only


# ═════════════ --bam-list fail-fast contract ═════════════════════════════
@XFAIL
def test_bam_list_missing_entry_fails_fast(tmp_path):
    """A missing BAM in --bam-list without --lenient-bam must exit non-zero
    (the --lenient-bam help text and the parser docstring both promise
    fail-fast). Today the entry is skipped and the run exits 0 with fewer
    samples."""
    fasta, ref = _build_reference(tmp_path)
    good = _build_bam(tmp_path, ref, name="good.bam")
    listfile = tmp_path / "bams.list"
    listfile.write_text(f"good {good}\nmissing {tmp_path / 'nope.bam'}\n")
    outdir = tmp_path / "out"
    outdir.mkdir()
    result = _invoke(
        [
            "dna",
            "-v",
            str(_build_vcf(tmp_path, with_fail=False)),
            "-L",
            str(listfile),
            "-f",
            str(fasta),
            "-o",
            str(outdir),
            "--format",
            "maf",
        ]
    )
    assert result.exit_code != 0, "missing list entry must fail fast without --lenient-bam"


def test_bam_list_missing_entry_lenient_skips(tmp_path):
    """Guard: with --lenient-bam the missing entry is skipped (ERROR logged)
    and the good sample still processes."""
    fasta, ref = _build_reference(tmp_path)
    good = _build_bam(tmp_path, ref, name="good.bam")
    listfile = tmp_path / "bams.list"
    listfile.write_text(f"good {good}\nmissing {tmp_path / 'nope.bam'}\n")
    outdir = tmp_path / "out"
    outdir.mkdir()
    result = _invoke(
        [
            "dna",
            "-v",
            str(_build_vcf(tmp_path, with_fail=False)),
            "-L",
            str(listfile),
            "-f",
            str(fasta),
            "-o",
            str(outdir),
            "--format",
            "maf",
            "--lenient-bam",
        ]
    )
    assert result.exit_code == 0, result.output
    assert glob.glob(str(outdir / "good*.maf")), "good sample must still be processed"


@XFAIL
def test_bam_list_read_error_fails_fast(tmp_path):
    """An unreadable --bam-list (here: a directory) raises OSError, which is
    currently logged and swallowed — the run continues with whatever --bam
    supplied and exits 0. A half-read (or unreadable) list is as fatal as a
    missing one: must exit non-zero."""
    fasta, ref = _build_reference(tmp_path)
    good = _build_bam(tmp_path, ref, name="good.bam")
    bad_list = tmp_path / "list_is_a_dir"
    bad_list.mkdir()
    outdir = tmp_path / "out"
    outdir.mkdir()
    result = _invoke(
        [
            "dna",
            "-v",
            str(_build_vcf(tmp_path, with_fail=False)),
            "-b",
            str(good),
            "-L",
            str(bad_list),
            "-f",
            str(fasta),
            "-o",
            str(outdir),
            "--format",
            "maf",
        ]
    )
    assert result.exit_code != 0, "unreadable bam-list must fail fast"


@XFAIL
def test_duplicate_sample_names_are_an_error(tmp_path):
    """Two --bam arguments whose sample names collide (same file stem in
    different directories) silently overwrite each other today — one BAM is
    never processed. Must exit non-zero naming the collision."""
    fasta, ref = _build_reference(tmp_path)
    d1 = tmp_path / "runA"
    d2 = tmp_path / "runB"
    d1.mkdir()
    d2.mkdir()
    bam1 = _build_bam(d1, ref, name="tumor.bam")
    bam2 = _build_bam(d2, ref, name="tumor.bam")
    outdir = tmp_path / "out"
    outdir.mkdir()
    result = _invoke(
        [
            "dna",
            "-v",
            str(_build_vcf(tmp_path, with_fail=False)),
            "-b",
            str(bam1),
            "-b",
            str(bam2),
            "-f",
            str(fasta),
            "-o",
            str(outdir),
            "--format",
            "maf",
        ]
    )
    assert result.exit_code != 0, "duplicate sample name must be a hard error"


# ═════════════ writer and merge silent data loss ═════════════════════════
@XFAIL
def test_maf_writer_warns_on_colliding_input_columns(tmp_path, caplog):
    """An input MAF that already carries a gbcms output column name (e.g.
    ref_count from a previous genotyping run) has its values replaced in
    place with no signal, while the writer's docstring claims originals are
    never overwritten. The fix warns, listing the collisions."""
    import logging

    from gbcms.core.kernel import CoordinateKernel
    from gbcms.io.output import MafWriter
    from gbcms.pipeline import _zero_counts

    v = CoordinateKernel.vcf_to_internal(chrom="chr1", pos=101, ref="C", alt="T")
    v = v.model_copy(
        update={
            "metadata": {
                "Hugo_Symbol": "G1",
                "Chromosome": "chr1",
                "Start_Position": "101",
                "End_Position": "101",
                "Reference_Allele": "C",
                "Tumor_Seq_Allele2": "T",
                "Tumor_Sample_Barcode": "S1",
                "ref_count": "999",  # collides with a gbcms output column
            }
        }
    )
    out = tmp_path / "out.maf"
    with caplog.at_level(logging.WARNING, logger="gbcms.io.output"):
        w = MafWriter(out)
        w.write(v, _zero_counts(), sample_name="S1")
        w.close()
    assert any(
        "ref_count" in rec.message and "replace" in rec.message.lower() for rec in caplog.records
    ), "collision with input-MAF columns must be warned about"


@XFAIL
def test_merge_does_not_coerce_float_count_strings_to_zero(tmp_path):
    """A count value of '12.0' (pandas/R round-trip formatting) is cast with
    strict=False and silently becomes null → 0 in the combined sums. The fix
    tolerates float-formatted integers."""
    from gbcms.merge import COMBINED_ADDITIVE_ALL

    metric = COMBINED_ADDITIVE_ALL[0]
    header = (
        "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
        f"Reference_Allele\tTumor_Seq_Allele2\t{metric}\n"
    )
    duplex = tmp_path / "d.maf"
    duplex.write_text(header + "G1\t1\t101\t101\tC\tT\t12.0\n")
    simplex = tmp_path / "s.maf"
    simplex.write_text(header + "G1\t1\t101\t101\tC\tT\t5\n")
    out = tmp_path / "merged.maf"
    result = _invoke(
        [
            "merge",
            "-i",
            f"duplex:{duplex}",
            "-i",
            f"simplex:{simplex}",
            "-o",
            str(out),
        ]
    )
    assert result.exit_code == 0, result.output
    rows = list(read_maf_output(out))
    combined = rows[0][f"simplex_duplex_{metric}"]
    assert combined == "17", f"12.0 + 5 must merge to 17, got {combined!r}"


@XFAIL
def test_batch_reader_rejects_ragged_rows(tmp_path):
    """A row with more fields than the header currently has its overflow
    silently truncated (and gbcms count columns are the trailing columns).
    The batch readers must fail loudly on ragged rows."""
    import polars as pl

    from gbcms.io.batch import read_maf

    p = tmp_path / "ragged.maf"
    p.write_text("A\tB\tC\nx\ty\tz\textra_field\n")
    with pytest.raises(pl.exceptions.ComputeError):
        read_maf(p)


# ═════════════ failure diagnostics ═══════════════════════════════════════
@XFAIL
def test_failed_sample_report_names_the_exception_type(tmp_path, monkeypatch):
    """A sample failure is reported via str(e) only — for a KeyError that is
    just the key, with no exception type and no traceback. The report must
    name the type."""
    from gbcms import pipeline as pipeline_mod

    def _boom(self, *a, **k):
        raise KeyError("some_missing_column")

    monkeypatch.setattr(pipeline_mod.Pipeline, "_write_output", _boom)
    fasta, ref = _build_reference(tmp_path)
    bam = _build_bam(tmp_path, ref)
    outdir = tmp_path / "out"
    outdir.mkdir()
    result = _invoke(
        _dna_args(_build_vcf(tmp_path, with_fail=False), bam, fasta, outdir, "--format", "maf")
    )
    assert result.exit_code != 0
    assert "KeyError" in result.output, (
        "failure report must include the exception type, got: " + result.output[-500:]
    )
