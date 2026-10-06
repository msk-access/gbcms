"""Group 7 (D4 #139): GTF loading without noodles-gtf or a cache.

Measured on GRCh38.111 (CYCLE_6.6.0_PLAN.md, "D4 GTF loading"): noodles-gtf parsed
every line before the feature check, 5.9 s of a 6.3 s load; a byte-level loop that
checks the feature column first loads the same exons in 0.7 s, so the bincode cache
and its build step save about 0.5 s a sample. The config accepts ``.gtf.gz``, but
the parser read it as text and failed ("stream did not contain valid UTF-8").

Contracts:
- a gzip or BGZF ``.gtf.gz`` loads exactly as the plain file;
- noodles' attribute grammar, except that whitespace runs may separate key and
  value (noodles read ``transcript_id  "T1"`` as the ID `` "T1"``);
- malformed exon lines (start > end, a coordinate past i32, an empty
  ``transcript_id``) are rejected and warned about, never loaded: noodles kept
  them, so their bogus boundaries moved ``exon_boundary_dist``;
- the cache is deprecated: ``--gtf-cache-dir`` and ``build-gtf-cache`` are
  accepted, warn, and write nothing; noodles-gtf, bincode, the cache module, its
  binding and the Nextflow step are gone.
"""

import gzip
import logging
import re
from pathlib import Path

import pysam
import pytest
from helpers import make_read
from rna_fixtures import E1, READ_LEN, SENSE, mk_ref, run_rna, write_bam, write_fasta, write_vcf
from typer.testing import CliRunner

from gbcms import _rs
from gbcms.cli import app

ROOT = Path(__file__).resolve().parents[1]
POS = E1[0] + 100  # 0-based SNV, mid-exon: 100 bases from either E1 edge
GTF_LOGGER = "_rs.annotation.gtf"
runner = CliRunner()

_EXONS = ((101, 300), (501, 600), (801, 1000))  # 1-based closed, as rna_fixtures


def _gtf_text(attrs='gene_id "G1"; transcript_id "T1";', extra=""):
    return "".join(f"chr1\tTEST\texon\t{s}\t{e}\t.\t+\t.\t{attrs}\n" for s, e in _EXONS) + extra


def _reads(ref):
    alt = next(b for b in "ACGT" if b != ref[POS])
    out = []
    for i in range(16):
        s = POS - 50 + (i % 9)
        base = alt if i < 6 else ref[POS]
        seq = ref[s:POS] + base + ref[POS + 1 : s + READ_LEN]
        out.append(make_read(f"r{i}", seq, s, ((0, READ_LEN),), flag=SENSE))
    return out, alt


def _run(tmp_path, gtf_text, name="gene.gtf", compress=None, extra=()):
    """One SNV row from `gbcms rna` with a GTF of `gtf_text` (in its own dir)."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    ref = mk_ref()
    reads, alt = _reads(ref)
    fa, bam = write_fasta(tmp_path, ref), write_bam(tmp_path, ref, reads)
    plain = tmp_path / "plain.gtf"
    plain.write_text(gtf_text)
    gtf = tmp_path / name
    if compress == "gzip":
        with gzip.open(gtf, "wt") as fh:
            fh.write(gtf_text)
    elif compress == "bgzf":
        pysam.tabix_compress(str(plain), str(gtf), force=True)
    else:
        gtf = plain
    (row,) = run_rna(
        tmp_path, write_vcf(tmp_path, [(POS + 1, ref[POS], alt)]), bam, fa, gtf, extra=extra
    )
    return row


# ── Input: .gtf.gz ───────────────────────────────────────────────────────────


@pytest.mark.xfail(strict=True, reason="the parser reads .gtf.gz as text and fails")
@pytest.mark.parametrize("compress", ["gzip", "bgzf"])
def test_gzipped_gtf_loads_like_the_plain_file(tmp_path, compress):
    """Every output field of a .gtf.gz run equals the plain run's."""
    plain = _run(tmp_path / "plain", _gtf_text())
    packed = _run(tmp_path / "packed", _gtf_text(), name="gene.gtf.gz", compress=compress)
    assert plain["transcript_read_counts"].startswith("T1:")
    assert packed == plain


# ── Grammar ──────────────────────────────────────────────────────────────────


@pytest.mark.xfail(strict=True, reason='noodles reads `transcript_id  "T1"` as the ID ` "T1"`')
def test_whitespace_runs_between_key_and_value_read_the_id(tmp_path):
    row = _run(tmp_path, _gtf_text(attrs='gene_id  "G1";  transcript_id  "T1";'))
    assert row["transcript_read_counts"].startswith("T1:"), row["transcript_read_counts"]


def test_unquoted_values_and_quoted_semicolons_load(tmp_path):
    """Guard: noodles' grammar is kept — an unquoted value runs to the next ';',
    and a quoted value may hold one."""
    row = _run(tmp_path, _gtf_text(attrs='note "a; b"; gene_id G1; transcript_id T1;'))
    assert row["transcript_read_counts"].startswith("T1:"), row["transcript_read_counts"]


# ── Malformed exon lines ─────────────────────────────────────────────────────
# Each junk exon puts a boundary 3 bases from the SNV (whose true distance is 100)
# if it is loaded.
_JUNK = {
    "start-after-end": f"chr1\tTEST\texon\t{POS + 4}\t{POS - 10}\t.\t+\t.\t"
    'gene_id "G9"; transcript_id "T9";\n',
    "end-past-i32": f"chr1\tTEST\texon\t{POS + 4}\t3000000000\t.\t+\t.\t"
    'gene_id "G9"; transcript_id "T9";\n',
    "empty-transcript-id": f"chr1\tTEST\texon\t{POS + 4}\t{POS + 30}\t.\t+\t.\t"
    'gene_id "G9"; transcript_id "";\n',
}


@pytest.mark.xfail(
    strict=True, reason="noodles loads these lines; their boundaries move the distance"
)
@pytest.mark.parametrize("kind", sorted(_JUNK))
def test_malformed_exon_lines_are_not_loaded(tmp_path, kind):
    clean = _run(tmp_path / "clean", _gtf_text())
    junk = _run(tmp_path / "junk", _gtf_text(extra=_JUNK[kind]))
    assert clean["exon_boundary_dist"] == "100"
    assert junk == clean


@pytest.mark.xfail(strict=True, reason="malformed exon lines are skipped at debug level only")
def test_rejected_exon_lines_are_warned_with_a_count(tmp_path, caplog):
    _rs.reset_log_caching()
    with caplog.at_level(logging.WARNING, logger=GTF_LOGGER):
        _run(tmp_path, _gtf_text(extra="".join(_JUNK[k] for k in sorted(_JUNK))))
    msgs = [r.getMessage() for r in caplog.records if r.name == GTF_LOGGER]
    hit = [m for m in msgs if "malformed" in m]
    assert hit, msgs
    assert re.search(r"\b3 exon lines?\b", hit[0]), hit[0]
    assert re.search(r"line \d+", hit[0]), hit[0]


# ── The cache is deprecated ──────────────────────────────────────────────────


@pytest.mark.xfail(strict=True, reason="--gtf-cache-dir still writes a bincode cache")
def test_gtf_cache_dir_is_accepted_ignored_and_warned(tmp_path):
    cache = tmp_path / "cache"
    clean = _run(tmp_path / "clean", _gtf_text())
    d = tmp_path / "flagged"
    flagged = _run(d, _gtf_text(), extra=("--gtf-cache-dir", str(cache)))
    assert flagged == clean
    assert not cache.exists() or not any(cache.iterdir())
    # the warning reaches the user in the CLI's output
    (tmp_path / "again").mkdir()
    args = ["rna", "-v", str(d / "variants.vcf"), "-b", f"S:{d / 'rna.s.bam'}"]
    args += ["-f", str(d / "ref.fasta"), "-o", str(tmp_path / "again")]
    args += ["--gtf", str(d / "plain.gtf"), "--gtf-cache-dir", str(cache)]
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    assert "deprecated" in result.output.lower(), result.output


@pytest.mark.xfail(strict=True, reason="noodles-gtf and bincode are still dependencies")
def test_no_noodles_gtf_or_bincode_dependency():
    cargo = (ROOT / "rust/Cargo.toml").read_text()
    deps = re.search(r"^\[dependencies\]\n(.*?)(?=^\[|\Z)", cargo, re.M | re.S).group(1)
    names = set(re.findall(r"^([A-Za-z0-9_-]+)\s*=", deps, re.M))
    assert "noodles-gtf" not in names
    assert "bincode" not in names
    assert not (ROOT / "rust/src/annotation/cache.rs").exists()


@pytest.mark.xfail(strict=True, reason="the binding still exposes the cache")
def test_the_binding_has_no_cache():
    stub = (ROOT / "src/gbcms/_rs.pyi").read_text()
    assert "gtf_cache_dir" not in stub
    assert "def build_gtf_cache" not in stub
    assert not hasattr(_rs, "build_gtf_cache")


@pytest.mark.xfail(strict=True, reason="the Nextflow pipeline still builds the cache")
def test_nextflow_has_no_cache_step_and_warns_on_the_old_param():
    nf = ROOT / "nextflow"
    assert not (nf / "modules/local/gbcms/build_gtf_cache").exists()
    main = (nf / "main.nf").read_text()
    assert "GBCMS_BUILD_GTF_CACHE" not in main
    assert re.search(
        r"log\.warn[^\n]*--gtf_cache", main
    ), "main.nf warns that --gtf_cache is deprecated"
    assert "gtf_cache" not in (nf / "modules/local/gbcms/rna/main.nf").read_text()
    assert "gtf_cache" not in (nf / "workflows/rna.nf").read_text()
    config = (nf / "nextflow.config").read_text()
    assert re.search(r"^\s*gtf_cache\s*=\s*null\b", config, re.M), "the old param defaults to null"
