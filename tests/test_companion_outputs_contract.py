"""Companion files are named from the main output's stem, every dot-suffix kept.

The mFSD Parquet, the mFSD report and the observations Parquet sit beside the
MAF/VCF as ``{sample}{suffix}.<kind>``. A dotted ``--suffix`` (``.genotyped``) or a
dotted sample name (``S.1``) must survive, so the companions pair with their output
and two samples never write the same companion file. ``--mfsd-report`` reads the
main output as a MAF, so it is refused with ``--format vcf``.
"""

import pysam
import pytest
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

POS0 = 100  # 0-based SNV position (REF A, ALT T)


def _reference(tmp_path):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "A" * 1000 + "\n")
    pysam.faidx(str(fasta))
    return fasta


def _pair(qname, base):
    """A proper pair whose mates both cover the SNV (insert size 200)."""
    reads = []
    for is_r2 in (False, True):
        a = pysam.AlignedSegment()
        a.query_name = qname
        a.query_sequence = "A" * 5 + base + "A" * 4
        a.flag = (1 | 2 | 128 | 16) if is_r2 else (1 | 2 | 64)
        a.reference_id = 0
        a.reference_start = POS0 - 5
        a.mapping_quality = 60
        a.cigartuples = [(0, 10)]
        a.query_qualities = [30] * 10  # type: ignore[assignment]
        a.next_reference_id = 0
        a.next_reference_start = POS0 - 5
        a.template_length = -200 if is_r2 else 200
        reads.append(a)
    return reads


def _bam(tmp_path, name):
    raw = tmp_path / f"{name}.raw.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"LN": 1000, "SN": "chr1"}]}
    with pysam.AlignmentFile(str(raw), "wb", header=header) as fh:
        for i in range(3):
            for r in _pair(f"ref{i}", "A"):
                fh.write(r)
        for i in range(3):
            for r in _pair(f"alt{i}", "T"):
                fh.write(r)
    bam = tmp_path / f"{name}.bam"
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    return bam


def _vcf(tmp_path):
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t{POS0 + 1}\t.\tA\tT\t.\t.\t.\n"
    )
    return vcf


def _run(tmp_path, samples, *extra):
    out = tmp_path / "out"
    args = ["dna", "-v", str(_vcf(tmp_path)), "-f", str(_reference(tmp_path)), "-o", str(out)]
    for s in samples:
        args += ["-b", f"{s}:{_bam(tmp_path, s.replace('.', '_'))}"]
    return runner.invoke(app, [*args, *extra]), out


def test_mfsd_companions_keep_a_dotted_suffix(tmp_path):
    pytest.importorskip("polars")
    res, out = _run(
        tmp_path,
        ["S1"],
        "--format",
        "maf",
        "--suffix",
        ".genotyped",
        "--mfsd-report",
        "--mfsd-report-min-alt",
        "1",
    )
    assert res.exit_code == 0, res.output
    names = sorted(p.name for p in out.iterdir())
    assert "S1.genotyped.maf" in names, names
    assert "S1.genotyped.fsd.parquet" in names, names
    assert "S1.genotyped.mfsd_report.html" in names, names


def test_dotted_sample_names_do_not_share_a_companion(tmp_path):
    res, out = _run(tmp_path, ["S.1", "S.2"], "--format", "maf", "--mfsd", "--mfsd-parquet")
    assert res.exit_code == 0, res.output
    parquets = sorted(p.name for p in out.glob("*.fsd.parquet"))
    assert parquets == ["S.1.fsd.parquet", "S.2.fsd.parquet"], parquets


def test_observations_parquet_honours_the_suffix(tmp_path):
    pytest.importorskip("pyarrow.parquet")
    res, out = _run(
        tmp_path, ["S.1"], "--format", "maf", "--suffix", ".genotyped", "--observations-parquet"
    )
    assert res.exit_code == 0, res.output
    names = sorted(p.name for p in out.iterdir())
    assert "S.1.genotyped.maf" in names, names
    assert "S.1.genotyped.observations.parquet" in names, names


def test_rna_observations_parquet_keeps_a_dotted_suffix_and_sample_name(tmp_path):
    """``gbcms rna`` names its companion as ``gbcms dna`` does."""
    pytest.importorskip("pyarrow.parquet")
    out = tmp_path / "out"
    res = runner.invoke(
        app,
        ["rna", "-v", str(_vcf(tmp_path)), "-f", str(_reference(tmp_path)), "-o", str(out)]
        + ["-b", f"S.1:{_bam(tmp_path, 'S_1')}"]
        + ["--format", "maf", "--suffix", ".genotyped", "--observations-parquet"],
    )
    assert res.exit_code == 0, res.output
    names = sorted(p.name for p in out.iterdir())
    assert names == ["S.1.genotyped.maf", "S.1.genotyped.observations.parquet"], names


def test_mfsd_report_with_vcf_output_is_refused_at_the_cli(tmp_path):
    res, out = _run(tmp_path, ["S1"], "--format", "vcf", "--mfsd-report")
    assert res.exit_code != 0, res.output
    assert not out.exists() or not list(out.iterdir()), "nothing is written"


def test_mfsd_report_with_vcf_output_is_refused_by_the_model(tmp_path):
    from pydantic import ValidationError

    from gbcms.models.core import OutputConfig, OutputFormat

    with pytest.raises(ValidationError, match="mfsd_report"):
        OutputConfig(
            directory=tmp_path,
            format=OutputFormat.VCF,
            mfsd=True,
            mfsd_parquet=True,
            mfsd_report=True,
        )
    # the MAF report is unchanged
    OutputConfig(
        directory=tmp_path,
        format=OutputFormat.MAF,
        mfsd=True,
        mfsd_parquet=True,
        mfsd_report=True,
    )
