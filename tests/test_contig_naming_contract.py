"""Target contract for issue #103: output keeps the input's contig naming.

Readers normalize contig names internally (strip ``chr``) so the variant file,
FASTA and BAM reconcile whatever their naming. The output must still use the
naming of the input it came from:

- VCF input -> MAF output: ``Chromosome`` and ``vcf_region`` use the input's name.
- VCF or MAF input -> VCF output: ``CHROM`` uses the input's name, and every
  record's contig is declared in a ``##contig`` line (reference lengths kept) —
  otherwise the VCF is malformed (htslib: "Contig ... is not defined in the
  header").
- MAF input -> MAF output already passes ``Chromosome`` through; guarded.
- Unprefixed (b37-style) input is unchanged everywhere.
- One INFO line per run names the two conventions when the input's naming
  differs from the reference's.

Committed red (xfail-strict) before the implementation.
"""

import glob
import random

import pysam
import pytest
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

XFAIL = pytest.mark.xfail(strict=True, reason="contig naming (#103): implementation pending")

READ_LEN = 100
POS1 = 201  # 1-based SNV position


def _ref():
    rng = random.Random(103)
    return "".join(rng.choice("ACGT") for _ in range(600))


def _alt(ref):
    base = ref[POS1 - 1]
    return "A" if base != "A" else "C"


def _fasta(tmp_path, ref, contig):
    fa = tmp_path / f"ref_{contig}.fasta"
    fa.write_text(f">{contig}\n{ref}\n")
    pysam.faidx(str(fa))
    return fa


def _bam(tmp_path, ref, contig):
    reads = []
    alt = _alt(ref)
    for i in range(8):
        s = POS1 - 60 + i
        seq = list(ref[s : s + READ_LEN])
        if i % 2:
            seq[POS1 - 1 - s] = alt
        reads.append(make_read(f"r{i}", "".join(seq), s, ((0, READ_LEN),)))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": contig, "LN": len(ref)}]}
    p = tmp_path / f"reads_{contig}.bam"
    with pysam.AlignmentFile(p, "wb", header=header) as fh:
        for a in sorted(reads, key=lambda a: a.reference_start):
            fh.write(a)
    sp = tmp_path / f"reads_{contig}.s.bam"
    pysam.sort("-o", str(sp), str(p))
    pysam.index(str(sp))
    return sp


def _vcf_input(tmp_path, ref, contig):
    vcf = tmp_path / "in.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"{contig}\t{POS1}\t.\t{ref[POS1 - 1]}\t{_alt(ref)}\t.\t.\t.\n"
    )
    return vcf


def _maf_input(tmp_path, ref, contig):
    maf = tmp_path / "in.maf"
    maf.write_text(
        "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\t"
        "Tumor_Seq_Allele2\tTumor_Sample_Barcode\n"
        f"G\t{contig}\t{POS1}\t{POS1}\t{ref[POS1 - 1]}\t{_alt(ref)}\tS\n"
    )
    return maf


def _run(tmp_path, variants, bam, fasta, fmt, outname):
    outdir = tmp_path / outname
    outdir.mkdir()
    result = runner.invoke(
        app,
        [
            "dna",
            "-v",
            str(variants),
            "-b",
            f"S:{bam}",
            "-f",
            str(fasta),
            "-o",
            str(outdir),
            "--format",
            fmt,
        ],
    )
    assert result.exit_code == 0, result.output
    out = glob.glob(str(outdir / f"*.{fmt}"))
    assert len(out) == 1, out
    return out[0], " ".join(result.output.split())


def _maf_row(path):
    (row,) = list(read_maf_output(path))
    assert int(row["alt_count"]) == 4, "counting must be unaffected by naming"
    return row


def _vcf_records(path, capfd):
    capfd.readouterr()
    with pysam.VariantFile(path) as vf:
        declared = set(vf.header.contigs)
        recs = list(vf)
    err = capfd.readouterr().err
    assert "not defined in the header" not in err, err
    assert len(recs) == 1
    return recs[0], declared


# ── #103 reproduction: chr-named input, chr FASTA, unprefixed BAM ─────────
@XFAIL
def test_vcf_input_maf_output_keeps_chr_naming(tmp_path):
    ref = _ref()
    path, _ = _run(
        tmp_path,
        _vcf_input(tmp_path, ref, "chr1"),
        _bam(tmp_path, ref, "1"),
        _fasta(tmp_path, ref, "chr1"),
        "maf",
        "o",
    )
    row = _maf_row(path)
    assert row["Chromosome"] == "chr1"
    assert row["vcf_region"] == f"chr1:{POS1}"


@XFAIL
def test_vcf_input_vcf_output_is_well_formed_with_chr_naming(tmp_path, capfd):
    ref = _ref()
    path, _ = _run(
        tmp_path,
        _vcf_input(tmp_path, ref, "chr1"),
        _bam(tmp_path, ref, "1"),
        _fasta(tmp_path, ref, "chr1"),
        "vcf",
        "o",
    )
    rec, declared = _vcf_records(path, capfd)
    assert rec.chrom == "chr1"
    assert "chr1" in declared


@XFAIL
def test_maf_input_vcf_output_keeps_chr_naming(tmp_path, capfd):
    ref = _ref()
    path, _ = _run(
        tmp_path,
        _maf_input(tmp_path, ref, "chr1"),
        _bam(tmp_path, ref, "1"),
        _fasta(tmp_path, ref, "chr1"),
        "vcf",
        "o",
    )
    rec, declared = _vcf_records(path, capfd)
    assert rec.chrom == "chr1"
    assert "chr1" in declared


def test_maf_input_maf_output_keeps_chr_naming(tmp_path):
    """Guard (green now, green after): MAF rows pass Chromosome through."""
    ref = _ref()
    path, _ = _run(
        tmp_path,
        _maf_input(tmp_path, ref, "chr1"),
        _bam(tmp_path, ref, "1"),
        _fasta(tmp_path, ref, "chr1"),
        "maf",
        "o",
    )
    assert _maf_row(path)["Chromosome"] == "chr1"


@XFAIL
def test_vcf_header_declares_input_naming_when_reference_differs(tmp_path, capfd):
    """chr-named input against an unprefixed reference: the ##contig line
    for the reference's '1' is written under the input's name 'chr1', with
    the reference's length, so the record's CHROM is declared."""
    ref = _ref()
    path, _ = _run(
        tmp_path,
        _vcf_input(tmp_path, ref, "chr1"),
        _bam(tmp_path, ref, "1"),
        _fasta(tmp_path, ref, "1"),
        "vcf",
        "o",
    )
    rec, declared = _vcf_records(path, capfd)
    assert rec.chrom == "chr1"
    assert declared == {"chr1"}
    with pysam.VariantFile(path) as vf:
        assert vf.header.contigs["chr1"].length == len(ref)


@XFAIL
def test_naming_difference_is_logged_once(tmp_path):
    ref = _ref()
    _, log = _run(
        tmp_path,
        _vcf_input(tmp_path, ref, "chr1"),
        _bam(tmp_path, ref, "1"),
        _fasta(tmp_path, ref, "1"),
        "maf",
        "o",
    )
    assert log.count("contig naming") == 1, log


# ── b37-style naming is untouched ──────────────────────────────────────────
@pytest.mark.parametrize("fmt", ["maf", "vcf"])
def test_unprefixed_naming_unchanged(tmp_path, capfd, fmt):
    """Guard (green now, green after): '1' everywhere stays '1', and nothing
    about naming is logged."""
    ref = _ref()
    path, log = _run(
        tmp_path,
        _vcf_input(tmp_path, ref, "1"),
        _bam(tmp_path, ref, "1"),
        _fasta(tmp_path, ref, "1"),
        fmt,
        "o",
    )
    if fmt == "maf":
        row = _maf_row(path)
        assert row["Chromosome"] == "1" and row["vcf_region"] == f"1:{POS1}"
    else:
        rec, declared = _vcf_records(path, capfd)
        assert rec.chrom == "1" and declared == {"1"}
    assert "contig naming" not in log
