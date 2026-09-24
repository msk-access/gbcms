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

Committed red (xfail-strict) before the implementation; flipped green with it.
"""

import glob
import logging
import random

import polars as pl
import pysam
import pytest
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app
from gbcms.core.kernel import CoordinateKernel
from gbcms.merge import merge_mafs
from gbcms.models.core import MergeConfig, Variant, VariantType
from gbcms.pipeline import _declared_contigs

runner = CliRunner()

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


def test_maf_input_vcf_output_declares_input_naming_when_reference_differs(tmp_path, capfd):
    """MAF input takes the same route as VCF input: its 'chr1' is declared with
    the unprefixed reference's length."""
    ref = _ref()
    path, _ = _run(
        tmp_path,
        _maf_input(tmp_path, ref, "chr1"),
        _bam(tmp_path, ref, "1"),
        _fasta(tmp_path, ref, "1"),
        "vcf",
        "o",
    )
    rec, declared = _vcf_records(path, capfd)
    assert rec.chrom == "chr1"
    assert declared == {"chr1"}


# ── Mitochondrion: chrM and MT are the same contig ─────────────────────────
def test_chrM_input_against_MT_reference(tmp_path, capfd):
    """The engine reconciles chrM with MT (not just a 'chr' strip); the header
    and the naming log must pair them the same way, so 'chrM' is declared with
    MT's length and the difference is logged."""
    ref = _ref()
    path, log = _run(
        tmp_path,
        _vcf_input(tmp_path, ref, "chrM"),
        _bam(tmp_path, ref, "MT"),
        _fasta(tmp_path, ref, "MT"),
        "vcf",
        "o",
    )
    rec, declared = _vcf_records(path, capfd)
    assert rec.chrom == "chrM"
    assert declared == {"chrM"}
    with pysam.VariantFile(path) as vf:
        assert vf.header.contigs["chrM"].length == len(ref)
    assert log.count("contig naming") == 1, log


def test_reference_aliases_declare_each_name_once():
    """A reference listing one contig under two aliases (chr1 and 1) must not
    produce a duplicate ##contig line for the input's name."""
    v = Variant(chrom="1", pos=0, ref="A", alt="C", variant_type=VariantType.SNP)
    v.original_chrom = "chr1"
    assert _declared_contigs([("chr1", 600), ("1", 600), ("2", 900)], [v]) == [
        ("chr1", 600),
        ("2", 900),
    ]


# ── gbcms merge across inputs that name contigs differently ────────────────
def test_merge_joins_across_contig_naming(tmp_path, caplog):
    """Runs from a 'chr'-named and an unprefixed variant file describe the same
    variant: merge joins them into one row in the first input's naming, keeps
    a row only the second input has in that input's naming, and logs it."""
    cols = [
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "ref_count",
        "alt_count",
    ]

    def maf(name, rows):
        p = tmp_path / name
        p.write_text("\n".join(["\t".join(cols), *("\t".join(r) for r in rows)]) + "\n")
        return p

    duplex = maf("duplex.maf", [["chr1", "100", "100", "A", "T", "20", "10"]])
    simplex = maf(
        "simplex.maf",
        [["1", "100", "100", "A", "T", "5", "2"], ["MT", "50", "50", "G", "C", "7", "1"]],
    )
    out = tmp_path / "merged.maf"
    with caplog.at_level(logging.INFO, logger="gbcms.merge"):
        merge_mafs(
            MergeConfig(
                inputs={"duplex": duplex, "simplex": simplex}, output=out, add_combined=False
            )
        )
    result = pl.read_csv(out, separator="\t", infer_schema_length=0).sort("Chromosome")
    assert result["Chromosome"].to_list() == ["MT", "chr1"]
    joined = result.filter(pl.col("Chromosome") == "chr1").row(0, named=True)
    assert (joined["duplex_ref_count"], joined["simplex_ref_count"]) == ("20", "5")
    assert not any(c.startswith("_") for c in result.columns), result.columns
    naming = [r.message for r in caplog.records if "name contigs differently" in r.message]
    assert len(naming) == 1, [r.message for r in caplog.records]


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


# ── Mitochondrion: every source reconciles chrM ~ M ~ MT ~ chrMT ───────────
@pytest.mark.xfail(strict=True, reason="FASTA fetch does not fold the mitochondrial aliases")
@pytest.mark.parametrize(
    "variant_name,fasta_name,bam_name",
    [("chrM", "MT", "MT"), ("MT", "chrM", "chrM"), ("chrMT", "chrM", "MT"), ("M", "MT", "chrM")],
)
def test_mitochondrial_names_count_across_sources(tmp_path, variant_name, fasta_name, bam_name):
    """Counting reconciles the mitochondrion's spellings on every side (BAM and
    FASTA alike), so the naming log's 'counting reconciles them' is true and
    the variant counts; the output keeps the input's name."""
    ref = _ref()
    path, _ = _run(
        tmp_path,
        _vcf_input(tmp_path, ref, variant_name),
        _bam(tmp_path, ref, bam_name),
        _fasta(tmp_path, ref, fasta_name),
        "maf",
        "o",
    )
    row = _maf_row(path)
    assert row["Chromosome"] == variant_name


# ── contig_key mirrors the engine's normalize_contig ───────────────────────
@pytest.mark.parametrize(
    "name,key",
    [
        # rust/src/shared/contig.rs tests, mirrored
        ("chr1", "1"),
        ("CHR7", "7"),
        ("chrX", "X"),
        ("1", "1"),
        ("MT", "MT"),
        ("M", "MT"),
        ("chrM", "MT"),
        ("chrMT", "MT"),
        ("chrm", "MT"),
        ("mt", "MT"),
        ("MT1", "MT1"),
        ("GL000220", "GL000220"),
        ("chrUn_gl000220", "Un_gl000220"),
    ],
)
def test_contig_key_mirrors_engine_rule(name, key):
    assert CoordinateKernel.contig_key(name) == key


# ── gbcms merge names each contig one way ──────────────────────────────────
_MERGE_COLS = [
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "ref_count",
    "alt_count",
]


def _merge_maf(tmp_path, name, rows, cols=_MERGE_COLS):
    p = tmp_path / name
    p.write_text("\n".join(["\t".join(cols), *("\t".join(r) for r in rows)]) + "\n")
    return p


@pytest.mark.xfail(strict=True, reason="rows only a later input has keep that input's naming")
def test_merge_writes_each_contig_one_way(tmp_path, caplog):
    """A row only a later input has takes the name the first input uses for
    that contig, so one file never names a contig two ways; each later input's
    naming difference is logged against the name actually written."""
    duplex = _merge_maf(tmp_path, "d.maf", [["chr1", "100", "100", "A", "T", "20", "10"]])
    simplex = _merge_maf(
        tmp_path,
        "s.maf",
        [["1", "100", "100", "A", "T", "5", "2"], ["1", "200", "200", "G", "C", "7", "1"]],
    )
    standard = _merge_maf(
        tmp_path,
        "t.maf",
        [["CHR1", "200", "200", "G", "C", "3", "1"], ["CHR1", "300", "300", "C", "A", "4", "0"]],
    )
    out = tmp_path / "merged.maf"
    with caplog.at_level(logging.INFO, logger="gbcms.merge"):
        merge_mafs(
            MergeConfig(
                inputs={"duplex": duplex, "simplex": simplex, "standard": standard},
                output=out,
                add_combined=False,
            )
        )
    result = pl.read_csv(out, separator="\t", infer_schema_length=0)
    assert result["Chromosome"].to_list() == ["chr1", "chr1", "chr1"]
    naming = [r.message for r in caplog.records if "name contigs differently" in r.message]
    assert len(naming) == 2, naming
    assert "'simplex'" in naming[0] and "2 row(s)" in naming[0], naming
    assert "'standard'" in naming[1] and "2 row(s)" in naming[1], naming


@pytest.mark.xfail(strict=True, reason="helper column names are not reserved")
@pytest.mark.parametrize("col", ["_contig_key", "_row_duplex", "_chrom_simplex"])
def test_merge_rejects_reserved_helper_columns(tmp_path, col):
    rows = [["chr1", "100", "100", "A", "T", "5", "1", "x"]]
    duplex = _merge_maf(tmp_path, "d.maf", rows, [*_MERGE_COLS, col])
    simplex = _merge_maf(tmp_path, "s.maf", [r[:-1] for r in rows])
    with pytest.raises(ValueError, match=col):
        merge_mafs(
            MergeConfig(
                inputs={"duplex": duplex, "simplex": simplex},
                output=tmp_path / "m.maf",
                add_combined=False,
            )
        )


@pytest.mark.xfail(strict=True, reason="one input naming a contig two ways is silent")
def test_merge_warns_when_one_input_names_a_contig_two_ways(tmp_path, caplog):
    """chrM and MT rows in one MAF share a contig key, so the same variant under
    both names joins twice (duplicate merged rows): say so."""
    duplex = _merge_maf(tmp_path, "d.maf", [["chrM", "100", "100", "A", "T", "20", "10"]])
    simplex = _merge_maf(
        tmp_path,
        "s.maf",
        [["chrM", "100", "100", "A", "T", "11", "1"], ["MT", "100", "100", "A", "T", "12", "1"]],
    )
    with caplog.at_level(logging.WARNING, logger="gbcms.merge"):
        merge_mafs(
            MergeConfig(
                inputs={"duplex": duplex, "simplex": simplex},
                output=tmp_path / "m.maf",
                add_combined=False,
            )
        )
    warns = [r.message for r in caplog.records if r.levelno == logging.WARNING]
    assert any("'simplex'" in w and "chrM" in w and "MT" in w for w in warns), warns
