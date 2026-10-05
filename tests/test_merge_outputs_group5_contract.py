"""Group 5 (merge, outputs, observability) contracts.

Operator decisions (M4 #194 2026-09-30; the rest 2026-10-05):
- M4 #194: a missing or non-numeric count cell in either flavour makes the
  combined cell NA, and merge warns once per column with the number of rows; the
  combined VAFs have the writers' four decimals.
- M5 #221: a row only a later input has takes the first input's columns from the
  earliest later input that has the row; the column set stays the first input's.
- M6 #223: every column gbcms writes is prefixed per input; merge's set is the
  writer's, in every mode.
- M2 #129: the merged MAF carries its own provenance and each input's; inputs
  from different gbcms versions warn; the pre-6.5.0 VCF-input shape mixed with
  another is refused; the log counts the rows each input lacks; a build names
  its commit.
- H1 #148: every output is written to a temp file and renamed: a failed run
  leaves nothing at the final path, and no temp file.
- O1 #130: the --rescue-mnp recount does not repeat the per-BAM warnings;
  records without bases are counted once each.
- O2 #131: the run start logs every resolved option, what count-changing options
  imply, and per-BAM facts; it warns when --min-baseq removes more than 10% of
  sampled bases, the header has ALT contigs, or more than 1% of primaries are
  hard-clipped, and never for unmarked duplicates.
- D6 #156: one QC-flags page lists every flag the code emits.
"""

import csv
import logging
import pathlib
import re

import pysam
import pytest
from helpers import make_read, read_maf_output
from test_input_group4_contract import READ, REF, P, _files, _maf, _ref_reads, _row
from typer.testing import CliRunner

from gbcms import __version__, _rs
from gbcms.cli import app

runner = CliRunner()
ENGINE_LOGGER = "_rs.counting.engine"
ROOT = pathlib.Path(__file__).resolve().parents[1]

_KEY = ["Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"]
_COUNTS = ["ref_count", "alt_count", "ref_count_fragment", "alt_count_fragment"]
_MCOLS = [*_KEY, *_COUNTS]


def _write(path, cols, rows, head=()):
    lines = [*head, "\t".join(cols), *("\t".join(r) for r in rows)]
    path.write_text("\n".join(lines) + "\n")
    return path


def _merge(tmp_path, inputs, caplog=None, add_combined=True):
    """Merge hand-written MAFs ({type: (columns, rows[, header lines])}) with the
    library call; the merged rows."""
    from gbcms.merge import merge_mafs
    from gbcms.models.core import MergeConfig

    paths = {}
    for name, spec in inputs.items():
        cols, rows = spec[0], spec[1]
        head = spec[2] if len(spec) > 2 else ()
        paths[name] = _write(tmp_path / f"{name}.maf", cols, rows, head)
    out = tmp_path / "merged.maf"
    config = MergeConfig(inputs=paths, output=out, add_combined=add_combined)
    if caplog is None:
        merge_mafs(config)
    else:
        with caplog.at_level(logging.INFO, logger="gbcms.merge"):
            merge_mafs(config)
    return list(read_maf_output(out))


def _cli_merge(tmp_path, inputs):
    args = ["merge"]
    for name, spec in inputs.items():
        cols, rows = spec[0], spec[1]
        head = spec[2] if len(spec) > 2 else ()
        args += ["--input", f"{name}:{_write(tmp_path / f'{name}.maf', cols, rows, head)}"]
    out = tmp_path / "merged.maf"
    res = runner.invoke(app, [*args, "--output", str(out)], env={"COLUMNS": "3000"})
    return res, out


def _dna(tmp_path, variants, bam, fa, *extra, fmt="maf", name="o"):
    res = runner.invoke(
        app,
        ["dna", "-v", str(variants), "-b", f"S:{bam}", "-f", str(fa), "-o", str(tmp_path / name)]
        + ["--format", fmt, *extra],
        env={"COLUMNS": "3000"},
    )
    return res, " ".join(res.output.split())


# ── M4 #194: a missing count is not a zero ───────────────────────────────────


@pytest.mark.parametrize("flavor", ["duplex", "simplex"])
@pytest.mark.parametrize("cell", ["NA", "nan", "inf", "", "text"])
def test_a_missing_count_cell_makes_the_combined_cell_na(tmp_path, caplog, cell, flavor):
    good = ["1", "100", "100", "A", "T", "5", "2", "4", "2"]
    bad = ["1", "100", "100", "A", "T", cell, "10", "18", "9"]
    rows = {"duplex": bad if flavor == "duplex" else good, "simplex": good}
    if flavor == "simplex":
        rows = {"duplex": good, "simplex": bad}
    merged = _merge(tmp_path, {t: (_MCOLS, [r]) for t, r in rows.items()}, caplog)
    row = merged[0]
    assert (
        row["simplex_duplex_ref_count"],
        row["simplex_duplex_total_count"],
        row["simplex_duplex_vaf"],
    ) == ("NA", "NA", "NA")
    assert row["simplex_duplex_alt_count"] == "12"
    warned = [
        r.message
        for r in caplog.records
        if r.levelno == logging.WARNING and "ref_count" in r.message and "1 row" in r.message
    ]
    assert len(warned) == 1, [r.message for r in caplog.records]


def test_combined_vafs_have_four_decimals(tmp_path):
    """The combined VAFs are written as the writers write theirs (4 decimals)."""
    merged = _merge(
        tmp_path,
        {
            "duplex": (_MCOLS, [["1", "100", "100", "A", "T", "20", "10", "20", "10"]]),
            "simplex": (_MCOLS, [["1", "100", "100", "A", "T", "5", "2", "5", "2"]]),
        },
    )
    assert merged[0]["simplex_duplex_vaf"] == "0.3243"
    assert merged[0]["simplex_duplex_vaf_fragment"] == "0.3243"


def test_numeric_and_absent_counts_still_combine(tmp_path):
    """Guard: '12.0' cells sum as integers, and a row one input lacks counts 0
    for that input (its counts are absent, not missing)."""
    merged = _merge(
        tmp_path,
        {
            "duplex": (_MCOLS, [["1", "100", "100", "A", "T", "12.0", "3", "12", "3"]]),
            "simplex": (
                _MCOLS,
                [["1", "100", "100", "A", "T", "5", "2", "5", "2"]]
                + [["1", "200", "200", "G", "C", "7", "1", "7", "1"]],
            ),
        },
    )
    assert [(r["simplex_duplex_ref_count"], r["simplex_duplex_alt_count"]) for r in merged] == [
        ("17", "5"),
        ("7", "1"),
    ]


# ── M5 #221: a later-only row keeps its annotations ──────────────────────────

_ACOLS = ["Hugo_Symbol", "Tumor_Sample_Barcode", *_MCOLS]


def test_a_row_only_a_later_input_has_keeps_its_annotations(tmp_path):
    """Its first-input columns come from the later input; a column only the
    later input has is not added; a row the first input has keeps its own."""
    merged = _merge(
        tmp_path,
        {
            "duplex": (_ACOLS, [["G1", "S", "1", "100", "100", "A", "T", "20", "10", "20", "10"]]),
            "simplex": (
                [*_ACOLS, "Only_Simplex"],
                [
                    ["X", "S", "1", "100", "100", "A", "T", "5", "2", "5", "2", "a"],
                    ["G2", "S", "1", "200", "200", "G", "C", "7", "1", "7", "1", "b"],
                ],
            ),
        },
    )
    assert [(r["Hugo_Symbol"], r["Tumor_Sample_Barcode"]) for r in merged] == [
        ("G1", "S"),
        ("G2", "S"),
    ]
    assert "Only_Simplex" not in merged[0]


def test_the_earliest_later_input_fills_the_row(tmp_path):
    rows = {
        "duplex": [["G1", "S", "1", "100", "100", "A", "T", "1", "1", "1", "1"]],
        "simplex": [["G2", "S", "1", "200", "200", "G", "C", "2", "1", "2", "1"]],
        "standard": [["G3", "S", "1", "200", "200", "G", "C", "3", "1", "3", "1"]],
    }
    merged = _merge(tmp_path, {t: (_ACOLS, r) for t, r in rows.items()}, add_combined=False)
    assert [r["Hugo_Symbol"] for r in merged] == ["G1", "G2"]


# ── M6 #223: every gbcms column is per input ─────────────────────────────────


def test_mfsd_and_rna_columns_are_prefixed_per_input(tmp_path):
    cols = [*_MCOLS, "mfsd_ref_mean", "mfsd_ks_valid", "rna_sense_depth", "exon_boundary_dist"]
    merged = _merge(
        tmp_path,
        {
            "duplex": (
                cols,
                [
                    ["1", "100", "100", "A", "T", "20", "10", "20", "10"]
                    + ["167.0", "True", "30", "4"]
                ],
            ),
            "simplex": (
                cols,
                [["1", "100", "100", "A", "T", "5", "2", "5", "2"] + ["171.0", "False", "7", "4"]],
            ),
        },
    )
    row = merged[0]
    assert (row["duplex_mfsd_ref_mean"], row["simplex_mfsd_ref_mean"]) == ("167.0", "171.0")
    assert (row["duplex_rna_sense_depth"], row["simplex_rna_sense_depth"]) == ("30", "7")
    assert not {"mfsd_ref_mean", "mfsd_ks_valid", "rna_sense_depth"} & set(row)


def test_merge_knows_every_column_the_writers_emit(tmp_path):
    from gbcms.io.output import MafWriter
    from gbcms.merge import ALL_GBCMS_BASENAMES

    modes = [
        {},
        {"mfsd": True},
        {"mode": "rna"},
        {"mode": "rna", "has_gtf": True},
        {"rescue_mnp": True},
        {"show_normalization": True},
    ]
    emitted = set()
    for i, kwargs in enumerate(modes):
        w = MafWriter(tmp_path / f"w{i}.maf", **kwargs)
        emitted |= set(w._gbcms_column_names())
        w.close()
    assert emitted - ALL_GBCMS_BASENAMES == set()


# ── M2 #129: inputs from different gbcms versions ────────────────────────────

_ROW_D = ["1", "100", "100", "A", "T", "20", "10", "20", "10"]
_ROW_S = ["1", "100", "100", "A", "T", "5", "2", "5", "2"]


def test_the_merged_maf_carries_its_own_and_each_inputs_provenance(tmp_path, monkeypatch):
    monkeypatch.setattr("sys.argv", ["gbcms", "merge", "--input", "duplex:d.maf"])
    res, out = _cli_merge(
        tmp_path,
        {
            "duplex": (_MCOLS, [_ROW_D], ["#gbcms v6.5.0", "#command gbcms dna"]),
            "simplex": (_MCOLS, [_ROW_S], ["#gbcms v6.5.0", "#command gbcms dna"]),
        },
    )
    assert res.exit_code == 0, res.output
    head = [x.rstrip("\n") for x in open(out) if x.startswith("#")]
    assert head[0].startswith(f"#gbcms v{__version__}")
    assert any(x.startswith("#command ") and " merge " in x for x in head)
    for t in ("duplex", "simplex"):
        assert any(x.startswith(f"#input {t}:") and "gbcms v6.5.0" in x for x in head), head
    assert len(list(read_maf_output(out))) == 1


def test_merge_warns_when_the_inputs_come_from_different_versions(tmp_path, caplog):
    _merge(
        tmp_path,
        {
            "duplex": (_MCOLS, [_ROW_D], ["#gbcms v6.4.0"]),
            "simplex": (_MCOLS, [_ROW_S], ["#gbcms v6.5.0"]),
        },
        caplog,
    )
    warned = [r.message for r in caplog.records if r.levelno == logging.WARNING]
    assert any("v6.4.0" in m and "v6.5.0" in m for m in warned), warned


def test_merge_of_one_version_does_not_warn_about_versions(tmp_path, caplog):
    """Guard (green now, green after)."""
    _merge(
        tmp_path,
        {
            "duplex": (_MCOLS, [_ROW_D], ["#gbcms v6.5.0"]),
            "simplex": (_MCOLS, [_ROW_S], ["#gbcms v6.5.0"]),
        },
        caplog,
    )
    assert not [
        r for r in caplog.records if r.levelno == logging.WARNING and "version" in r.message
    ]


def test_merge_refuses_the_pre_650_vcf_input_shape_mixed_with_another(tmp_path):
    """A VCF-input MAF from before 6.5.0 has vcf_pos but no vcf_ref/vcf_alt; its
    rows do not join a later one's, so merge stops instead of writing split rows."""
    old = [*_MCOLS, "vcf_pos"]
    new = [*_MCOLS, "vcf_pos", "vcf_ref", "vcf_alt"]
    res, out = _cli_merge(
        tmp_path,
        {
            "duplex": (old, [[*_ROW_D, "100"]]),
            "simplex": (new, [[*_ROW_S, "100", "A", "T"]]),
        },
    )
    assert res.exit_code != 0
    assert "6.5.0" in " ".join(res.output.split())
    assert not out.exists()


def test_merge_of_two_pre_650_vcf_input_mafs_still_joins(tmp_path):
    """Guard: one representation on both sides joins."""
    old = [*_MCOLS, "vcf_pos"]
    merged = _merge(
        tmp_path, {"duplex": (old, [[*_ROW_D, "100"]]), "simplex": (old, [[*_ROW_S, "100"]])}
    )
    assert len(merged) == 1


def test_merge_log_counts_the_rows_an_input_lacks(tmp_path, caplog):
    """Not rows whose REF count is 0: a present row with ref_count 0 is not
    missing from its input."""
    _merge(
        tmp_path,
        {
            "duplex": (
                _MCOLS,
                [
                    ["1", "100", "100", "A", "T", "0", "10", "0", "10"],
                    _ROW_D[:1] + ["200", "200", "G", "C", "3", "1", "3", "1"],
                ],
            ),
            "simplex": (_MCOLS, [["1", "100", "100", "A", "T", "0", "2", "0", "2"]]),
        },
        caplog,
    )
    text = " | ".join(r.message for r in caplog.records)
    assert re.search(r"'simplex' lacks 1 of 2", text), text
    assert "'duplex' lacks" not in text


def test_the_provenance_line_names_the_build_commit():
    """A build carries its commit (from git, or GBCMS_BUILD_COMMIT where the
    source has no git, as in Docker), so two builds of one version differ."""
    from gbcms.io.output import provenance_line

    commit = _rs.build_commit()
    assert re.fullmatch(r"[0-9a-f]{7,40}|", commit)
    assert provenance_line() == f"#gbcms v{__version__}" + (f" ({commit})" if commit else "")


# ── H1 #148: no partial output files ─────────────────────────────────────────


def _five_snvs(tmp_path):
    rows = [
        _row(P + 1 + 10 * k, REF[P + 10 * k], "A" if REF[P + 10 * k] != "A" else "C")
        for k in range(5)
    ]
    fa, bam = _files(tmp_path, _ref_reads())
    return _maf(tmp_path, rows), fa, bam


@pytest.mark.parametrize("fmt", ["maf", "vcf"])
def test_a_failed_write_leaves_no_output_file(tmp_path, monkeypatch, fmt):
    from gbcms.io import output as out_mod

    cls = out_mod.VcfWriter if fmt == "vcf" else out_mod.MafWriter
    orig, calls = cls.write, {"n": 0}

    def fail_third(self, *a, **k):
        calls["n"] += 1
        if calls["n"] == 3:
            raise OSError(28, "No space left on device")
        return orig(self, *a, **k)

    monkeypatch.setattr(cls, "write", fail_third)
    maf, fa, bam = _five_snvs(tmp_path)
    res, _ = _dna(tmp_path, maf, bam, fa, fmt=fmt)
    assert res.exit_code != 0
    assert sorted(p.name for p in (tmp_path / "o").iterdir()) == []


@pytest.mark.parametrize("fmt", ["maf", "vcf"])
def test_a_run_leaves_only_its_outputs(tmp_path, fmt):
    """Guard: success writes the output at its name, and no temp file."""
    maf, fa, bam = _five_snvs(tmp_path)
    res, _ = _dna(tmp_path, maf, bam, fa, fmt=fmt)
    assert res.exit_code == 0, res.output
    assert sorted(p.name for p in (tmp_path / "o").iterdir()) == [f"S.{fmt}"]


def test_a_failed_merge_write_leaves_no_output_file(tmp_path, monkeypatch):
    import polars as pl

    orig = pl.DataFrame.write_csv

    def partial(self, file, *a, **k):
        orig(self.head(0), file, *a, **k)
        raise OSError(28, "No space left on device")

    monkeypatch.setattr(pl.DataFrame, "write_csv", partial)
    with pytest.raises(OSError):
        _merge(tmp_path, {"duplex": (_MCOLS, [_ROW_D]), "simplex": (_MCOLS, [_ROW_S])})
    assert not (tmp_path / "merged.maf").exists()
    assert not [p for p in tmp_path.iterdir() if "merged" in p.name]


@pytest.mark.parametrize("command", ["normalize", "convert"])
def test_a_failed_normalize_or_convert_write_leaves_no_output_file(tmp_path, monkeypatch, command):
    orig, calls = csv.DictWriter.writerow, {"n": 0}

    def fail_second(self, row):
        calls["n"] += 1
        if calls["n"] == 2:
            raise OSError(28, "No space left on device")
        return orig(self, row)

    maf, fa, _ = _five_snvs(tmp_path)
    out = tmp_path / ("n.tsv" if command == "normalize" else "c.vcf")
    if command == "convert":
        out = tmp_path / "c.maf"
        maf = tmp_path / "v.vcf"
        maf.write_text(
            f"##fileformat=VCFv4.2\n##contig=<ID=1,length={len(REF)}>\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "".join(
                f"1\t{P + 1 + 10 * k}\t.\t{REF[P + 10 * k]}\tA\t.\t.\t.\n"
                for k in range(5)
                if REF[P + 10 * k] != "A"
            )
        )
    monkeypatch.setattr(csv.DictWriter, "writerow", fail_second)
    res = runner.invoke(app, [command, "-v", str(maf), "-f", str(fa), "-o", str(out)])
    assert res.exit_code != 0
    assert not out.exists()
    assert not [
        p
        for p in tmp_path.iterdir()
        if p.name.startswith(".") or "partial" in p.name or "tmp" in p.name
    ]


def test_a_failed_parquet_write_leaves_no_output_file(tmp_path, monkeypatch):
    """The parquet writers write to a temp path the run renames: a failure
    inside one leaves nothing at the final path."""
    from gbcms import pipeline as pl_mod

    real = pl_mod._get_rs()

    class Failing:
        def __getattr__(self, name):
            return getattr(real, name)

        @staticmethod
        def write_fsd_parquet(path, *a, **k):
            pathlib.Path(path).write_bytes(b"PAR1partial")
            raise OSError(28, "No space left on device")

    monkeypatch.setattr(pl_mod, "_get_rs", lambda: Failing())
    maf, fa, bam = _five_snvs(tmp_path)
    res, _ = _dna(tmp_path, maf, bam, fa, "--mfsd", "--mfsd-parquet")
    assert res.exit_code != 0
    assert not list((tmp_path / "o").glob("*.parquet*"))
    assert not [p for p in (tmp_path / "o").iterdir() if p.name.startswith(".")]


# ── O1 #130: per-BAM warnings once ───────────────────────────────────────────


def _with_record_without_bases(tmp_path, reads, at):
    fa, bam = _files(tmp_path, reads)
    raw = tmp_path / "raw2.bam"
    with (
        pysam.AlignmentFile(str(bam)) as src,
        pysam.AlignmentFile(str(raw), "wb", template=src) as dst,
    ):
        for r in src:
            dst.write(r)
        x = pysam.AlignedSegment(dst.header)
        x.query_name, x.flag, x.reference_id, x.reference_start = "nobases", 0, 0, at
        x.mapping_quality, x.cigartuples = 60, ((0, READ),)
        dst.write(x)
    out = tmp_path / "s2.bam"
    pysam.sort("-o", str(out), str(raw))
    pysam.index(str(out))
    return fa, out


def test_the_rescue_recount_does_not_repeat_the_per_bam_warnings(tmp_path):
    a0 = "A" if REF[P] != "A" else "C"
    a1 = "A" if REF[P + 1] != "A" else "C"
    hap = REF[:P] + a0 + REF[P + 1 :]
    reads = _ref_reads() + [
        make_read(f"c{i}", hap[P - 40 + i : P - 40 + i + READ], P - 40 + i, ((0, READ),))
        for i in range(5)
    ]
    fa, bam = _with_record_without_bases(tmp_path, reads, P - 30)
    maf = _maf(tmp_path, [_row(P + 1, REF[P : P + 2], a0 + a1)])
    res, text = _dna(tmp_path, maf, bam, fa, "--rescue-mnp", "--umi-tag", "RX")
    assert res.exit_code == 0, res.output
    assert "now reports component" in text  # the recount ran
    assert text.count("stored without bases") == 1, text
    assert text.count("carries the RX tag") == 1, text


def test_a_record_without_bases_is_counted_once_across_bins(tmp_path, caplog):
    """Two one-variant bins both fetch the record: it is one record."""
    fa, bam = _with_record_without_bases(tmp_path, _ref_reads(), P - 30)
    vs = [
        _rs.Variant("1", P + k, REF[P + k], "A" if REF[P + k] != "A" else "C", "SNP")
        for k in (0, 10)
    ]
    args = {
        "min_mapq": 20,
        "min_baseq": 20,
        "filter_duplicates": True,
        "filter_secondary": True,
        "filter_supplementary": True,
        "filter_qc_failed": False,
        "filter_improper_pair": False,
        "filter_indel": False,
        "threads": 1,
    }
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _rs.count_bam_binned(str(bam), vs, [None, None], bin_max_variants=1, **args)
    msgs = [r.message for r in caplog.records if "without bases" in r.message]
    assert len(msgs) == 1 and "skipped 1 record(s)" in msgs[0], msgs
    assert "upper bound" not in msgs[0]


# ── O2 #131: the run start says what the run does ────────────────────────────


def _leaf_fields(model, prefix=""):
    from pydantic import BaseModel

    out = []
    for name, field in model.model_fields.items():
        ann = field.annotation
        if isinstance(ann, type) and issubclass(ann, BaseModel):
            out += _leaf_fields(ann, f"{name}.")
        else:
            out.append(name)
    return out


def test_the_run_start_logs_every_resolved_option(tmp_path):
    from gbcms.models.core import GbcmsDnaConfig

    maf, fa, bam = _five_snvs(tmp_path)
    res, text = _dna(tmp_path, maf, bam, fa)
    assert res.exit_code == 0, res.output
    block = text[text.index("Run settings") :]
    missing = [
        f
        for f in _leaf_fields(GbcmsDnaConfig)
        if f not in ("command_line",) and f"{f}=" not in block
    ]
    assert missing == [], missing


def test_the_run_start_says_what_count_changing_options_imply(tmp_path):
    maf, fa, bam = _five_snvs(tmp_path)
    res, text = _dna(tmp_path, maf, bam, fa, "--no-filter-duplicates", "--min-mapq", "0", "--mfsd")
    assert res.exit_code == 0, res.output
    assert "duplicate reads are counted" in text
    assert "--min-mapq 0" in text and "multi-mapped" in text
    assert "--mfsd adds" in text


def _bam_with(tmp_path, reads, extra_sq=()):
    fa = tmp_path / "ref.fa"
    fa.write_text(f">1\n{REF}\n")
    pysam.faidx(str(fa))
    sq = [{"SN": "1", "LN": len(REF)}] + [{"SN": n, "LN": 1000} for n in extra_sq]
    raw = tmp_path / "raw.bam"
    with pysam.AlignmentFile(
        str(raw), "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": sq}
    ) as fh:
        for r in reads:
            fh.write(r)
    bam = tmp_path / "s.bam"
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    return fa, bam


def test_the_run_start_logs_each_bams_properties(tmp_path):
    fa, bam = _bam_with(tmp_path, _ref_reads())
    maf = _maf(tmp_path, [_row(P + 1, REF[P], "A" if REF[P] != "A" else "C")])
    res, text = _dna(tmp_path, maf, bam, fa)
    assert res.exit_code == 0, res.output
    for fact in (
        "duplicates flagged",
        "base qualities",
        "below --min-baseq",
        "ALT contigs",
        "hard-clipped",
    ):
        assert fact in text, fact
    assert "duplicate" not in " ".join(x for x in res.output.splitlines() if "WARN" in x)


@pytest.mark.parametrize("case", ["min-baseq", "alt contigs", "hard clips"])
def test_the_run_start_warns_on_bam_properties_that_change_counts(tmp_path, case):
    reads = _ref_reads()
    extra = ()
    if case == "min-baseq":
        reads = [
            make_read(
                f"q{i}",
                REF[P - 50 + i : P - 50 + i + READ],
                P - 50 + i,
                ((0, READ),),
                quals=[15] * READ,
            )
            for i in range(6)
        ]
    elif case == "alt contigs":
        extra = ("chr1_KI270706v1_alt",)
    else:
        reads += [
            make_read(
                f"h{i}",
                REF[P - 45 + i : P - 45 + i + READ - 5],
                P - 45 + i,
                ((5, 5), (0, READ - 5)),
            )
            for i in range(3)
        ]
    fa, bam = _bam_with(tmp_path, reads, extra)
    maf = _maf(tmp_path, [_row(P + 1, REF[P], "A" if REF[P] != "A" else "C")])
    res, text = _dna(tmp_path, maf, bam, fa)
    assert res.exit_code == 0, res.output
    want = {
        "min-baseq": "--min-baseq 20 removes",
        "alt contigs": "ALT contig",
        "hard clips": "hard-clipped",
    }[case]
    warns = " ".join(x for x in res.output.splitlines() if "WARN" in x)
    assert want in " ".join(warns.split()), res.output


def test_no_bam_property_warning_at_ordinary_settings(tmp_path):
    """Guard: Q30 reads, no duplicates flagged, no ALT contigs, no hard clips."""
    fa, bam = _bam_with(tmp_path, _ref_reads())
    maf = _maf(tmp_path, [_row(P + 1, REF[P], "A" if REF[P] != "A" else "C")])
    res, _ = _dna(tmp_path, maf, bam, fa)
    warns = " ".join(x for x in res.output.splitlines() if "WARN" in x)
    for word in ("--min-baseq", "ALT contig", "hard-clipped", "duplicate"):
        assert word not in warns, res.output


# ── D6 #156: one page for every QC flag ──────────────────────────────────────

_NOT_FLAGS = {"UTF-8", "VARIANT_KEY", "BH-FDR"}


def _emitted_flags():
    lit = re.compile(r'[rbf]?"(?:[^"\\\n]|\\.)*"')
    noise = re.compile(r"^(GBCMS_|RUST_|PYO3|DEFAULT_|MAX_|MIN_|LOG_|ALLELE_|PHASE_|SB_|FSB_)")
    flags = set()
    files = list((ROOT / "rust" / "src").rglob("*.rs")) + list(
        (ROOT / "src" / "gbcms").rglob("*.py")
    )
    for p in files:
        s = p.read_text()
        if p.suffix == ".rs":
            # The unit-test module, not the first #[cfg(test)] (some guard an import).
            m = re.search(r"#\[cfg\(test\)\]\s*(?:pub\s+)?mod\s+\w+\s*\{", s)
            s = s[: m.start()] if m else s
        code = "\n".join(x for x in s.splitlines() if not x.lstrip().startswith(("//", "#")))
        for literal in lit.findall(code):
            for t in re.findall(r"\b([A-Z][A-Z0-9]+(?:[_-][A-Z0-9]+)+)\b", literal):
                if not noise.match(t) and t not in _NOT_FLAGS:
                    flags.add(t)
    return flags


def test_the_qc_flags_page_lists_every_flag_the_code_emits():
    page = ROOT / "docs" / "reference" / "qc-flags.md"
    text = page.read_text()
    flags = _emitted_flags()
    assert len(flags) > 30
    assert sorted(f for f in flags if f not in text) == []


def test_the_output_reference_and_glossary_link_the_qc_flags_page():
    for page in ("output-formats.md", "glossary.md"):
        assert "qc-flags.md" in (ROOT / "docs" / "reference" / page).read_text(), page
