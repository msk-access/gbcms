"""`--rescue-homopolymer` is deprecated in 6.7.0 and removed in 6.8.0.

Every entry point accepts the flag, warns once and ignores it: the twin is never
built, `WARN_HOMOPOLYMER_DECOMP` is never set, and the output is byte-identical to
a run without the flag. (hp1: the twin won at 1 of 20 candidate rows, AD +2.)
"""

import logging
from pathlib import Path

import pytest
from helpers import make_read
from rna_fixtures import (
    SENSE,
    run_rna,
    write_bam,
    write_fasta,
    write_gtf,
    write_vcf,
)
from test_decomposed_twin_contract import RNA_RUN, RUN_LEN, _ref_read, _rna_ref, _twin_carrier
from test_observed_allele_contract import (
    SOX2_ROW,
    _ref,
    _ref_reads,
    _run,
    _twin_carriers,
)
from typer.main import get_command

import gbcms
from gbcms.cli import app
from gbcms.models.core import GbcmsDnaConfig, OutputConfig
from gbcms.pipeline import Pipeline


class _Notices(logging.Handler):
    """Collects the deprecation notices (the CLI reconfigures the `gbcms` logger, so
    caplog on the root logger does not see them)."""

    def __init__(self):
        super().__init__(logging.WARNING)
        self.messages = []

    def emit(self, record):
        self.messages.append(record.getMessage())

    @property
    def count(self):
        return sum(NOTICE in m for m in self.messages)


@pytest.fixture
def notices():
    h = _Notices()
    log = logging.getLogger("gbcms.deprecations")
    log.addHandler(h)
    yield h
    log.removeHandler(h)


ROOT = Path(__file__).resolve().parents[1]
NOTICE = "deprecated in 6.7.0 and removed in 6.8.0"


def _data_lines(path):
    """The output bytes without the `#` provenance lines (they hold the command line)."""
    return [ln for ln in Path(path).read_bytes().splitlines() if not ln.startswith(b"#")]


def _maf(out):
    (path,) = sorted(Path(out).glob("*.maf"))
    return _data_lines(path)


def test_dna_warns_and_output_is_identical(tmp_path, notices):
    ref = _ref()
    reads = _twin_carriers(ref, 20) + _ref_reads(ref, 10)
    (flagged,) = _run(tmp_path, ref, reads, [SOX2_ROW], "--rescue-homopolymer")
    assert notices.count == 1, notices.messages
    (tmp_path / "plain").mkdir()
    (plain,) = _run(tmp_path / "plain", ref, reads, [SOX2_ROW])
    assert notices.count == 1, "the run without the flag does not warn"
    assert flagged == plain
    assert _maf(tmp_path / "out") == _maf(tmp_path / "plain" / "out")


def test_rna_warns_and_output_is_identical(tmp_path, notices):
    ref = _rna_ref()
    reads = [
        _twin_carrier(f"s{i}", ref, RNA_RUN, RNA_RUN - 60 + i, flag=SENSE) for i in range(8)
    ] + [_ref_read(f"r{i}", ref, RNA_RUN - 50 + i, flag=SENSE) for i in range(10)]
    args = (
        write_vcf(tmp_path, [(RNA_RUN + 1, "C" * RUN_LEN, "T")]),
        write_bam(tmp_path, ref, reads),
        write_fasta(tmp_path, ref),
        write_gtf(tmp_path),
    )
    run_rna(tmp_path, *args, outname="flag", extra=("--rescue-homopolymer",))
    assert notices.count == 1, notices.messages
    run_rna(tmp_path, *args, outname="plain")
    assert notices.count == 1, "the run without the flag does not warn"
    assert _maf(tmp_path / "flag") == _maf(tmp_path / "plain")


def test_pipeline_warns(tmp_path, notices):
    (tmp_path / "in.txt").write_text("")

    def config(**over):
        return GbcmsDnaConfig(
            variant_file=tmp_path / "in.txt",
            bam_files={"s": tmp_path / "in.txt"},
            reference_fasta=tmp_path / "in.txt",
            output=OutputConfig(directory=tmp_path / "o"),
            **over,
        )

    Pipeline(config())
    assert notices.count == 0
    Pipeline(config(rescue_homopolymer=True))
    assert notices.count == 1, notices.messages


def test_observe_molecules_warns_and_is_identical(tmp_path):
    import pysam

    from gbcms.models.core import Variant, VariantType

    ref = "A" * 100 + "CCCCCC" + "T" + "A" * 100
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + ref + "\n")
    pysam.faidx(str(fasta))
    from helpers import build_bam

    reads = []
    for i in range(6):
        seq = list(ref[90:120])
        seq[15] = "T"
        reads.append(make_read(f"alt{i}", "".join(seq), 90, ((0, 30),)))
    reads += [make_read(f"ref{i}", ref[90:120], 90, ((0, 30),)) for i in range(4)]
    bam = build_bam(tmp_path, reads, filename="d.bam")
    variants = [
        Variant(chrom="chr1", pos=100, ref="CCCCCC", alt="T", variant_type=VariantType.DELETION)
    ]
    plain = gbcms.observe_molecules(bam, variants, reference_fasta=str(fasta))
    with pytest.warns(DeprecationWarning, match="rescue-homopolymer"):
        flagged = gbcms.observe_molecules(
            bam, variants, reference_fasta=str(fasta), rescue_homopolymer=True
        )

    def key(r):
        return [(o.variant_index, o.molecule_hash, o.allele, o.best_qual) for o in r.observations]

    assert key(flagged) == key(plain)
    assert flagged.variant_status == plain.variant_status


def test_help_and_nextflow_docs_say_deprecated():
    for cmd in ("dna", "rna"):
        (opt,) = [
            p for p in get_command(app).commands[cmd].params if "--rescue-homopolymer" in p.opts
        ]
        assert "Deprecated in 6.7.0" in opt.help and "removed in 6.8.0" in opt.help
        assert "OBSERVED_ALLELE" in opt.help and "COEXISTING_ALLELE" in opt.help
    for rel in ("nextflow/README.md", "docs/nextflow/parameters.md"):
        (line,) = [
            ln for ln in (ROOT / rel).read_text().splitlines() if "`--rescue_homopolymer`" in ln
        ]
        assert "eprecated" in line and "6.8.0" in line and "OBSERVED_ALLELE" in line, rel
    config = (ROOT / "nextflow/nextflow.config").read_text()
    assert "rescue_homopolymer" in config
    assert "deprecated in 6.7.0" in config
