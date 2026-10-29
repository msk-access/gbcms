"""Preparation records each row's rows within a read's reach only for the sequence judge (#304).

Mechanism. `Variant.reach_rows` (the input's other rows within a read's reach of a row,
whole-input so they never depend on the bin) is read by the sequence judge alone. Under the
default `--read-judge rules` nothing reads it, so preparation that is told the judge skips it
(the engine is output-aware: the intent is plumbed, nothing is computed to be discarded).
Told nothing, preparation records it as before, so a caller that prepares and then counts
under the judge gets the same rows.
"""

import random

import pytest
from helpers import make_read, write_contig
from typer.testing import CliRunner

from gbcms import _rs
from gbcms.cli import app


def _contig():
    rng = random.Random(356)
    return "".join(rng.choice("ACGT") for _ in range(800))


def _rows(contig):
    """An SNV, a deletion 30 bases on (whose REF reaches it), and an SNV far away."""
    alt = lambda p: "A" if contig[p] != "A" else "C"  # noqa: E731
    return [
        _rs.Variant("1", 300, contig[300], alt(300), "SNP"),
        _rs.Variant("1", 329, contig[329:340], contig[329], "DEL"),
        _rs.Variant("1", 700, contig[700], alt(700), "SNP"),
    ]


def _prepare(tmp_path, **kw):
    contig = _contig()
    fa, _ = write_contig(tmp_path, contig, [make_read("r", contig[:100], 0, ((0, 100),))])
    return _rs.prepare_variants(_rows(contig), fa, 5, False, 1, True, **kw)


def test_the_rules_path_prepares_no_reach_rows(tmp_path):
    prepared = _prepare(tmp_path, read_judge="rules")
    assert all(p.gbcms_status == "PASS" for p in prepared)
    assert [p.variant.reach_rows for p in prepared] == [[], [], []]


@pytest.mark.parametrize("judge", ["sequence", "haplotype"])
def test_the_judge_path_prepares_the_reach_rows_as_before(tmp_path, judge):
    """Told the sequence judge (or its former name), the rows are those preparation records
    when told nothing: the near SNV and deletion each list the other; the far SNV none."""
    told = _prepare(tmp_path, read_judge=judge)
    untold = _prepare(tmp_path)
    assert [p.variant.reach_rows for p in told] == [p.variant.reach_rows for p in untold]
    snv, deletion, far = (p.variant.reach_rows for p in told)
    assert [r[0] for r in snv] == [329] and [r[0] for r in deletion] == [300] and far == []


def test_an_unknown_read_judge_fails_loudly(tmp_path):
    with pytest.raises(ValueError, match="read_judge"):
        _prepare(tmp_path, read_judge="sequnce")


@pytest.mark.parametrize("judge", ["rules", "sequence"])
def test_the_pipeline_prepares_for_the_judge_it_counts_with(tmp_path, monkeypatch, judge):
    """A `gbcms dna` run hands its read judge to preparation."""
    contig = _contig()
    reads = [make_read(f"r{i}", contig[250 + i : 350 + i], 250 + i, ((0, 100),)) for i in range(5)]
    fa, bam = write_contig(tmp_path, contig, reads)
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=800>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"1\t301\t.\t{contig[300]}\t{'A' if contig[300] != 'A' else 'C'}\t.\t.\t.\n"
    )
    seen = []
    real = _rs.prepare_variants

    def spy(*args, **kwargs):
        seen.append(kwargs.get("read_judge"))
        return real(*args, **kwargs)

    monkeypatch.setattr(_rs, "prepare_variants", spy)
    res = CliRunner().invoke(
        app,
        ["dna", "-v", str(vcf), "-b", f"S:{bam}", "-f", fa, "-o", str(tmp_path / "out")]
        + ["--format", "maf", "--read-judge", judge],
    )
    assert res.exit_code == 0, res.output
    assert seen and set(seen) == {judge}
