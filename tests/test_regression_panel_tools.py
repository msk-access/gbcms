"""The regression panel's comparison and attribution tools (scripts/regression_panel/,
docs/development/regression-panel.md), on synthetic tier and output files: no BAMs.
"""

import csv
import importlib.util
import os
from pathlib import Path

import pytest

TOOLS = Path(__file__).resolve().parents[1] / "scripts" / "regression_panel"
KEY = ("Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2")
PANEL_COLS = KEY + ("signout_t_ref_count", "signout_t_alt_count", "panel_strata")
OUT_COLS = PANEL_COLS + (
    "ref_count",
    "alt_count",
    "partial_alt",
    "total_count",
    "ref_count_fragment",
    "alt_count_fragment",
    "gbcms_status",
)


def _load(name):
    spec = importlib.util.spec_from_file_location(name, TOOLS / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _write(path, cols, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        w.writerows(rows)


def _read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# Panel variants of one sample: an SNV, an insertion 30 bp away, a deletion 500 bp away.
VARIANTS = [
    ("1", "100", "100", "A", "T", "90", "10", "shape:SNV;vaf:VAF10-40%"),
    ("1", "130", "131", "-", "G", "80", "20", "shape:INS_1;context:INS:HOMOPOLYMER"),
    ("1", "630", "631", "AC", "-", "95", "5", "shape:DEL_2-5;context:DEL:unique"),
]


def _tier(tmp_path):
    tier = tmp_path / "tier"
    _write(tier / "mafs" / "r001.maf", PANEL_COLS, VARIANTS)
    _write(tier / "mafs" / "r002.maf", PANEL_COLS, VARIANTS)
    _write(
        tier / "runs.tsv",
        ("run_id", "tag", "arm", "mode", "root", "bam_relpath", "variants", "extra_args"),
        [
            ("r001_dna_tumor", "r001", "dna", "dna", "dmp", "a.bam", "mafs/r001.maf", ""),
            ("r002_dna_duplex", "r002", "dna", "dna", "dmp", "b.bam", "mafs/r002.maf", ""),
            ("r002_dna_simplex", "r002", "dna", "dna", "dmp", "c.bam", "mafs/r002.maf", ""),
            ("r001_normal", "r001", "normal", "dna", "dmp", "n.bam", "mafs/r001.maf", ""),
        ],
    )
    return tier


def _out(root, build, rid, alts):
    """One run's MAF output: the panel rows with the given ALT counts (read = fragment)."""
    rows = [
        v + ("50", str(a), "0", str(50 + a), "50", str(a), "PASS")
        for v, a in zip(VARIANTS, alts, strict=True)
    ]
    _write(root / build / rid / "S.maf", OUT_COLS, rows)


def test_compare_sums_access_flavors_and_reports_concordance_by_stratum(tmp_path):
    compare = _load("compare_panel")
    tier, out = _tier(tmp_path), tmp_path / "out"
    for build, alts in (("base", (10, 12, 5)), ("new", (10, 20, 5))):
        _out(out, build, "r001_dna_tumor", alts)
        _out(out, build, "r002_dna_duplex", (6, 12, 3))
        _out(out, build, "r002_dna_simplex", (4, 8, 2))
        _out(out, build, "r001_normal", (0, 0, 0))
        (out / build / "times.0.tsv").write_text(
            "run_id\tarm\texit\twall_s\tmax_rss_kb\nr001_dna_tumor\tdna\t0\t2.0\t1000\n"
        )
    (out / "new" / "failed.0.txt").write_text("r001_normal\tnormal\n")
    rep = tmp_path / "rep"
    import sys

    argv = sys.argv
    sys.argv = ["compare_panel.py", str(tier), str(out), "base", "new", str(rep)]
    try:
        compare.main()
    finally:
        sys.argv = argv
    summ = {r["arm"]: r for r in _read(rep / "version_summary.tsv")}
    assert summ["dna"]["rows"] == "9" and summ["dna"]["rows count-changed"] == "1"
    cells = _read(rep / "version_rows.tsv")
    assert {(c["run_id"], c["column"], c["base"], c["new"]) for c in cells} >= {
        ("r001_dna_tumor", "alt_count", "12", "20")
    }
    conc = {(r["stratum"], r["level"]): r for r in _read(rep / "concordance_by_stratum.tsv")}
    # The insertion: r001 base 12 / new 20 against 20; r002 duplex + simplex 12 + 8 = 20.
    ins = conc[("shape:INS_1", "fragment")]
    assert ins["n"] == "2" and ins["within_base"] == "0.500" and ins["within_new"] == "1.000"
    assert conc[("all", "read")]["n"] == "6"
    assert len(_read(rep / "normals.tsv")) == 3
    gate = (rep / "gate.txt").read_text()
    assert "failed runs: base 0, new 1" in gate


def test_attribution_prepare_check_and_trail(tmp_path):
    attribute = _load("attribute")
    tier, out, attrib = _tier(tmp_path), tmp_path / "out", tmp_path / "attrib"
    # BASE -> NEW changes the insertion's alt_count (12 -> 20) in r001 only.
    _out(out, "base", "r001_dna_tumor", (10, 12, 5))
    _out(out, "new", "r001_dna_tumor", (10, 20, 5))
    for rid in ("r002_dna_duplex", "r002_dna_simplex", "r001_normal"):
        for build in ("base", "new"):
            _out(out, build, rid, (1, 1, 1))
    attribute.prepare(str(tier), str(out), "base", "new", str(attrib), 100)
    runs = _read(attrib / "runs.tsv")
    assert [r["run_id"] for r in runs] == ["r001_dna_tumor"]
    reduced = _read(tier / runs[0]["variants"])
    # The changed insertion and the SNV 30 bp away; not the deletion 500 bp away.
    assert [v["Start_Position"] for v in reduced] == ["100", "130"]

    def reduced_out(build, alts):
        rows = [
            v + ("50", str(a), "0", str(50 + a), "50", str(a), "PASS")
            for v, a in zip(VARIANTS[:2], alts, strict=True)
        ]
        _write(attrib / "out" / build / "r001_dna_tumor" / "S.maf", OUT_COLS, rows)

    reduced_out("base", (10, 12))
    reduced_out("new", (10, 20))
    attribute.check(str(tier), str(out), "base", "new", str(attrib))
    # Checkpoints: the insertion moves 12 -> 15 at cp1 and 15 -> 20 at cp3.
    cps = tmp_path / "checkpoints.tsv"
    _write(
        cps,
        ("order", "label", "sha", "interval_merges"),
        [
            ("1", "#1 a", "aaa1111", "#1"),
            ("2", "#2 b", "bbb2222", "#2"),
            ("3", "#3 c", "ccc3333", "#3"),
        ],
    )
    reduced_out("cp_aaa1111", (10, 15))
    reduced_out("cp_bbb2222", (10, 15))
    reduced_out("cp_ccc3333", (10, 20))
    attribute.attribute(str(tier), str(out), "base", "new", str(attrib), str(cps))
    cells = _read(attrib / "attribution_cells.tsv")
    alt = [c for c in cells if c["column"] == "alt_count"]
    assert len(alt) == 1
    assert alt[0]["changed_at"] == "#1 a [#1]: 12 -> 15 ; #3 c [#3]: 15 -> 20"
    assert _read(attrib / "unattributed.tsv") == []


def test_attribution_flags_a_last_checkpoint_that_is_not_new(tmp_path):
    attribute = _load("attribute")
    tier, out, attrib = _tier(tmp_path), tmp_path / "out", tmp_path / "attrib"
    _out(out, "base", "r001_dna_tumor", (10, 12, 5))
    _out(out, "new", "r001_dna_tumor", (10, 20, 5))
    for rid in ("r002_dna_duplex", "r002_dna_simplex", "r001_normal"):
        for build in ("base", "new"):
            _out(out, build, rid, (1, 1, 1))
    attribute.prepare(str(tier), str(out), "base", "new", str(attrib), 100)
    cps = tmp_path / "checkpoints.tsv"
    _write(cps, ("order", "label", "sha", "interval_merges"), [("1", "#1 a", "aaa1111", "#1")])
    rows = [
        v + ("50", str(a), "0", str(50 + a), "50", str(a), "PASS")
        for v, a in zip(VARIANTS[:2], (10, 18), strict=True)
    ]
    _write(attrib / "out" / "cp_aaa1111" / "r001_dna_tumor" / "S.maf", OUT_COLS, rows)
    attribute.attribute(str(tier), str(out), "base", "new", str(attrib), str(cps))
    un = _read(attrib / "unattributed.tsv")
    assert [u["column"] for u in un if u["column"] == "alt_count"] == ["alt_count"]


@pytest.mark.skipif(os.name != "posix", reason="bash runner")
def test_the_runner_shards_and_skips_finished_runs(tmp_path):
    """run_panel.sh with a stub gbcms: each shard runs its share once; a second pass
    runs nothing (finished runs are skipped)."""
    tier = _tier(tmp_path)
    stub = tmp_path / "gbcms"
    stub.write_text(
        '#!/bin/sh\nwhile [ $# -gt 0 ]; do [ "$1" = "-o" ] && mkdir -p "$2" && touch "$2/S.maf"; shift; done\n'
    )
    stub.chmod(0o755)
    out = tmp_path / "out"
    import subprocess

    def run(shard):
        env = {**os.environ, "SHARD": str(shard), "NSHARDS": "2"}
        subprocess.run(
            [
                "bash",
                str(TOOLS / "run_panel.sh"),
                str(tier),
                str(tier / "runs.tsv"),
                "b",
                str(stub),
                "/dmp",
                "/forte",
                "b37.fa",
                "hg38.fa",
                "x.gtf",
                str(out),
            ],
            env=env,
            check=True,
            capture_output=True,
        )

    run(0)
    run(1)
    done = sorted(p.name for p in (out / "b").iterdir() if p.is_dir())
    assert done == ["r001_dna_tumor", "r001_normal", "r002_dna_duplex", "r002_dna_simplex"]
    timed = sum(len(_read(p)) for p in (out / "b").glob("times.*.tsv"))
    assert timed == 4
    run(0)
    run(1)
    assert sum(len(_read(p)) for p in (out / "b").glob("times.*.tsv")) == 4
