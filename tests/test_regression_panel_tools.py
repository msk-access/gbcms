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
    # Matched: reads for the IMPACT run, duplex + simplex fragments for ACCESS.
    assert conc[("shape:INS_1", "matched")]["within_new"] == "1.000"
    assert len(_read(rep / "normals.tsv")) == 3
    gate = (rep / "gate.txt").read_text()
    assert "failed runs: base 0, new 1" in gate


def test_compare_reads_vcf_output_field_by_field(tmp_path):
    """A VCF run: a new INFO field is a header difference, not a changed record; a
    changed FORMAT value is one changed cell."""
    compare = _load("compare_panel")
    tier = tmp_path / "tier"
    _write(tier / "mafs" / "r001.maf", PANEL_COLS, VARIANTS)
    _write(
        tier / "runs.tsv",
        ("run_id", "tag", "arm", "mode", "root", "bam_relpath", "variants", "extra_args"),
        [("r001_vcf", "r001", "vcf", "dna", "dmp", "a.bam", "mafs/r001.maf", "--format vcf")],
    )
    head = '##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description="d">\n'
    cols = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    for build, extra_head, info, ad in (
        ("base", "", "DP=60", "50,10"),
        (
            "new",
            '##INFO=<ID=MAF_START,Number=1,Type=Integer,Description="m">\n',
            "DP=60;MAF_START=100",
            "50,12",
        ),
    ):
        d = tmp_path / "out" / build / "r001_vcf"
        d.mkdir(parents=True)
        (d / "S.vcf").write_text(
            head + extra_head + cols + f"1\t100\t.\tA\tT\t.\t.\t{info}\tAD\t{ad}\n"
        )
    rep = tmp_path / "rep"
    import sys

    argv = sys.argv
    sys.argv = ["compare_panel.py", str(tier), str(tmp_path / "out"), "base", "new", str(rep)]
    try:
        compare.main()
    finally:
        sys.argv = argv
    cells = _read(rep / "version_rows.tsv")
    assert [(c["column"], c["base"], c["new"]) for c in cells] == [("FORMAT:AD", "50,10", "50,12")]
    hdr = {(r["arm"], r["format"]): r for r in _read(rep / "header_diff.tsv")}
    assert hdr[("vcf", "vcf")]["only_in_new"] == "INFO:MAF_START"


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


# ── From the review: each asserts the correct behaviour ─────────────────────────
RUN_COLS = ("run_id", "tag", "arm", "mode", "root", "bam_relpath", "variants", "extra_args")
V3 = [("1", "100", "100", "A", "T"), ("1", "130", "131", "-", "G"), ("1", "160", "160", "C", "G")]


def _compare(tier, out, rep):
    import sys

    compare = _load("compare_panel")
    argv = sys.argv
    sys.argv = ["compare_panel.py", str(tier), str(out), "base", "new", str(rep)]
    try:
        compare.main()
    except SystemExit:
        pass
    finally:
        sys.argv = argv


def _tier_one(tmp_path, arm="dna", rid="r001_dna_tumor", extra=""):
    tier = tmp_path / "tier"
    _write(tier / "mafs" / "r001.maf", KEY, V3)
    _write(
        tier / "runs.tsv",
        RUN_COLS,
        [(rid, "r001", arm, "dna", "dmp", "a.bam", "mafs/r001.maf", extra)],
    )
    return tier


@pytest.mark.xfail(
    strict=True, reason="the mfsd arm is reduced; its BH q-values depend on the row set"
)
def test_arms_with_bh_qvalues_are_attributed_on_their_full_variant_files(tmp_path, capsys):
    attribute = _load("attribute")
    tier = _tier_one(tmp_path, "mfsd", "r001_mfsd_duplex", "--mfsd")
    out, attrib = tmp_path / "out", tmp_path / "attrib"
    cols = KEY + ("alt_count", "mfsd_pval_alt_ref", "mfsd_qval_alt_ref")
    full_b = [
        V3[0] + ("5", "0.01", "0.03"),
        V3[1] + ("7", "0.02", "0.03"),
        V3[2] + ("1", "0.04", "0.04"),
    ]
    full_n = [
        V3[0] + ("5", "0.01", "0.03"),
        V3[1] + ("9", "0.02", "0.03"),
        V3[2] + ("1", "0.04", "0.04"),
    ]
    for b, rows in (("base", full_b), ("new", full_n)):
        _write(out / b / "r001_mfsd_duplex" / "S.maf", cols, rows)
        _write(attrib / "out" / b / "r001_mfsd_duplex" / "S.maf", cols, rows)
    attribute.prepare(str(tier), str(out), "base", "new", str(attrib), 20)  # pad 20 would drop 160
    runs = _read(attrib / "runs.tsv")
    assert len(_read(tier / runs[0]["variants"])) == 3
    attribute.check(str(tier), str(out), "base", "new", str(attrib))
    assert "reduction exact" in capsys.readouterr().out


@pytest.mark.xfail(strict=True, reason="GNU time's status line reaches the time records")
@pytest.mark.skipif(os.name != "posix", reason="bash runner")
def test_a_failed_run_under_gnu_time_does_not_crash_compare(tmp_path):
    import subprocess

    fake_time = tmp_path / "gnu_time"  # GNU time: a status line on failure, then the format
    fake_time.write_text(
        "#!/bin/sh\n"
        'if [ "$1" = "-f" ] && [ "$#" -eq 3 ]; then exit 0; fi\n'
        'fmt=$2; o=$4; shift 4; "$@"; rc=$?\n'
        '[ $rc -ne 0 ] && echo "Command exited with non-zero status $rc" > "$o" || : > "$o"\n'
        'echo "$fmt" | sed "s/%e/0.01/; s/%M/1234/" >> "$o"; exit $rc\n'
    )
    fake_time.chmod(0o755)
    runner = tmp_path / "run_panel.sh"
    runner.write_text((TOOLS / "run_panel.sh").read_text().replace("/usr/bin/time", str(fake_time)))
    stub = tmp_path / "gbcms"
    stub.write_text("#!/bin/sh\nexit 1\n")
    stub.chmod(0o755)
    tier, out = _tier_one(tmp_path), tmp_path / "out"
    for b in ("base", "new"):
        subprocess.run(
            [
                "bash",
                str(runner),
                str(tier),
                str(tier / "runs.tsv"),
                b,
                str(stub),
                "/d",
                "/f",
                "a",
                "b",
                "g",
                str(out),
            ],
            check=True,
            capture_output=True,
        )
    _compare(tier, out, tmp_path / "rep")
    assert "failed runs: base 1, new 1" in (tmp_path / "rep" / "gate.txt").read_text()


@pytest.mark.xfail(strict=True, reason="a MAF in place counts as finished")
@pytest.mark.skipif(os.name != "posix", reason="bash runner")
def test_a_run_that_failed_after_writing_its_maf_is_retried_and_counted(tmp_path):
    import subprocess

    stub = tmp_path / "gbcms"
    stub.write_text(
        '#!/bin/sh\nwhile [ $# -gt 0 ]; do [ "$1" = "-o" ] && o=$2; shift; done\n'
        'mkdir -p "$o"; touch "$o/S.maf"; exit 1\n'
    )
    stub.chmod(0o755)
    tier = _tier_one(tmp_path, "mfsd", "r001_mfsd_duplex", "--mfsd --mfsd-parquet")
    out = tmp_path / "out"
    for b in ("base", "new", "new"):  # the operator re-submits the candidate
        subprocess.run(
            [
                "bash",
                str(TOOLS / "run_panel.sh"),
                str(tier),
                str(tier / "runs.tsv"),
                b,
                str(stub),
                "/d",
                "/f",
                "a",
                "b",
                "g",
                str(out),
            ],
            check=True,
            capture_output=True,
        )
    assert (
        sum(len(_read(p)) for p in (out / "new").glob("times.*.tsv")) == 2
    ), "the failed run was not retried"
    _compare(tier, out, tmp_path / "rep")
    assert "failed runs: base 1, new 1" in (tmp_path / "rep" / "gate.txt").read_text()


@pytest.mark.xfail(strict=True, reason="records only in base are not reported")
def test_compare_reports_vcf_records_only_in_base(tmp_path):
    tier, out = _tier_one(tmp_path, "vcf", "r001_vcf", "--format vcf"), tmp_path / "out"
    head = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    recs = {
        "base": "1\t100\t.\tA\tT\t.\tPASS\tDP=60\tAD\t50,10\n1\t160\t.\tC\tG\t.\tPASS\tDP=60\tAD\t50,10\n",
        "new": "1\t100\t.\tA\tT\t.\tPASS\tDP=60\tAD\t50,10\n",
    }
    for b, body in recs.items():
        (out / b / "r001_vcf").mkdir(parents=True)
        (out / b / "r001_vcf" / "S.vcf").write_text(head + body)
    _compare(tier, out, tmp_path / "rep")
    summ = _read(tmp_path / "rep" / "version_summary.tsv")[0]
    assert summ.get("records only in base") == "1", summ


@pytest.mark.xfail(strict=True, reason="rows are matched by position")
def test_compare_aligns_rows_by_locus_when_row_counts_differ(tmp_path):
    tier, out = _tier_one(tmp_path), tmp_path / "out"
    cols = KEY + ("alt_count",)
    _write(out / "base" / "r001_dna_tumor" / "S.maf", cols, [V3[0] + ("5",), V3[2] + ("3",)])
    _write(
        out / "new" / "r001_dna_tumor" / "S.maf",
        cols,
        [V3[0] + ("5",), V3[1] + ("4",), V3[2] + ("3",)],
    )
    _compare(tier, out, tmp_path / "rep")
    changed = {r["column"] for r in _read(tmp_path / "rep" / "version_rows.tsv")}
    assert not changed & set(KEY), f"key columns reported as changed cells: {changed}"
    summ = _read(tmp_path / "rep" / "version_summary.tsv")[0]
    assert summ.get("rows only in new") == "1", summ


@pytest.mark.xfail(strict=True, reason="prepare assumes the same rows in both versions")
def test_attribute_prepare_survives_a_row_count_difference(tmp_path):
    attribute = _load("attribute")
    tier, out = _tier_one(tmp_path), tmp_path / "out"
    cols = KEY + ("alt_count",)
    _write(out / "base" / "r001_dna_tumor" / "S.maf", cols, [V3[0] + ("5",), V3[2] + ("3",)])
    _write(
        out / "new" / "r001_dna_tumor" / "S.maf",
        cols,
        [V3[0] + ("6",), V3[1] + ("4",), V3[2] + ("3",)],
    )
    attribute.prepare(str(tier), str(out), "base", "new", str(tmp_path / "attrib"), 100)
    assert len(_read(tmp_path / "attrib" / "runs.tsv")) == 1


@pytest.mark.xfail(strict=True, reason="a row only in one version is skipped silently")
def test_attribute_traces_a_row_present_in_only_one_version(tmp_path):
    attribute = _load("attribute")
    tier, out, attrib = _tier_one(tmp_path), tmp_path / "out", tmp_path / "attrib"
    cols = KEY + ("alt_count",)
    moved = ("1", "130", "130", "-", "G")  # NEW writes End_Position differently
    _write(out / "base" / "r001_dna_tumor" / "S.maf", cols, [V3[0] + ("5",), V3[1] + ("4",)])
    _write(out / "new" / "r001_dna_tumor" / "S.maf", cols, [V3[0] + ("5",), moved + ("9",)])
    _write(
        attrib / "runs.tsv",
        RUN_COLS,
        [("r001_dna_tumor", "r001", "dna", "dna", "dmp", "a.bam", "x", "")],
    )
    cps = tmp_path / "cps.tsv"
    _write(cps, ("order", "label", "sha", "interval_merges"), [("1", "#1 a", "aaa1111", "#1")])
    _write(
        attrib / "out" / "base" / "r001_dna_tumor" / "S.maf", cols, [V3[0] + ("5",), V3[1] + ("4",)]
    )
    _write(
        attrib / "out" / "cp_aaa1111" / "r001_dna_tumor" / "S.maf",
        cols,
        [V3[0] + ("5",), moved + ("9",)],
    )
    attribute.attribute(str(tier), str(out), "base", "new", str(attrib), str(cps))
    cells = _read(attrib / "attribution_cells.tsv")
    assert {c["column"] for c in cells} == {"(row)"} and len(cells) == 2, cells


@pytest.mark.xfail(strict=True, reason="a missing checkpoint output hides its step")
def test_attribute_does_not_skip_a_missing_checkpoint(tmp_path):
    attribute = _load("attribute")
    tier, out, attrib = _tier_one(tmp_path), tmp_path / "out", tmp_path / "attrib"
    cols = KEY + ("alt_count",)
    _write(out / "base" / "r001_dna_tumor" / "S.maf", cols, [V3[0] + ("12",)])
    _write(out / "new" / "r001_dna_tumor" / "S.maf", cols, [V3[0] + ("20",)])
    _write(
        attrib / "runs.tsv",
        RUN_COLS,
        [("r001_dna_tumor", "r001", "dna", "dna", "dmp", "a.bam", "x", "")],
    )
    cps = tmp_path / "cps.tsv"
    _write(
        cps,
        ("order", "label", "sha", "interval_merges"),
        [
            ("1", "#1 a", "aaa1111", "#1"),
            ("2", "#2 b", "bbb2222", "#2"),
            ("3", "#3 c", "ccc3333", "#3"),
        ],
    )
    _write(attrib / "out" / "base" / "r001_dna_tumor" / "S.maf", cols, [V3[0] + ("12",)])
    _write(attrib / "out" / "cp_aaa1111" / "r001_dna_tumor" / "S.maf", cols, [V3[0] + ("15",)])
    _write(
        attrib / "out" / "cp_ccc3333" / "r001_dna_tumor" / "S.maf", cols, [V3[0] + ("20",)]
    )  # cp2 failed
    attribute.attribute(str(tier), str(out), "base", "new", str(attrib), str(cps))
    assert _read(attrib / "attribution_cells.tsv") == []
    assert [u["column"] for u in _read(attrib / "unattributed.tsv")] == ["alt_count"]


@pytest.mark.xfail(strict=True, reason="gbcms inherits the run list as its stdin")
@pytest.mark.skipif(os.name != "posix", reason="bash runner")
def test_a_command_that_reads_stdin_does_not_eat_the_run_list(tmp_path):
    import subprocess

    stub = tmp_path / "gbcms"
    stub.write_text(
        "#!/bin/sh\ncat > /dev/null\n"
        'while [ $# -gt 0 ]; do [ "$1" = "-o" ] && mkdir -p "$2" && touch "$2/S.maf"; shift; done\n'
    )
    stub.chmod(0o755)
    tier, out = tmp_path / "t i–er", tmp_path / "o u–t"  # spaces and an en-dash
    _write(tier / "mafs" / "r001.maf", KEY, V3)
    _write(
        tier / "runs.tsv",
        RUN_COLS,
        [
            ("r001_dna_tumor", "r001", "dna", "dna", "dmp", "a.bam", "mafs/r001.maf", ""),
            ("r002_dna_tumor", "r002", "dna", "dna", "dmp", "b.bam", "mafs/r001.maf", ""),
            ("f001_rna", "f001", "rna", "rna", "forte", "c.bam", "mafs/r001.maf", ""),
        ],
    )
    subprocess.run(
        [
            "bash",
            str(TOOLS / "run_panel.sh"),
            str(tier),
            str(tier / "runs.tsv"),
            "b",
            str(stub),
            "/d m p",
            "/f",
            "a",
            "b",
            "g",
            str(out),
        ],
        check=True,
        capture_output=True,
    )
    done = sorted(p.name for p in (out / "b").iterdir() if p.is_dir())
    assert done == ["f001_rna", "r001_dna_tumor", "r002_dna_tumor"], done
