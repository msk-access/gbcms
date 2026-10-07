"""Group 7 (D4 #139): GTF loading without noodles-gtf or a cache.

Measured on GRCh38.111 (PR #230): noodles-gtf parsed
every line before the feature check; the byte-level loader reads each line's
feature and chromosome first and loads a whole Ensembl GTF in 2 s instead of 9 s, so
the bincode cache and its build step saved 1-1.4 s a sample. The config accepted
``.gtf.gz``, but the parser read it as text and failed ("stream did not contain
valid UTF-8").

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
runner = CliRunner()

_EXONS = ((101, 300), (501, 600), (801, 1000))  # 1-based closed, as rna_fixtures


def _gtf_text(attrs='gene_id "G1"; transcript_id "T1";', extra="", filler=0):
    """The gene model's exon lines (after `filler` gene lines on chr2, which push
    them past the first compressed block) plus `extra`."""
    head = "".join(
        f'chr2\tTEST\tgene\t{100 * i + 1}\t{100 * i + 50}\t.\t+\t.\tgene_id "F{i}";\n'
        for i in range(filler)
    )
    exons = "".join(f"chr1\tTEST\texon\t{s}\t{e}\t.\t+\t.\t{attrs}\n" for s, e in _EXONS)
    return head + exons + extra


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
    data = gtf_text if isinstance(gtf_text, bytes) else gtf_text.encode()
    plain = tmp_path / "plain.gtf"
    plain.write_bytes(data)
    gtf = tmp_path / name
    if compress == "gzip":
        gtf.write_bytes(gzip.compress(data))
    elif compress == "members":  # concatenated gzip members, split mid-file
        half = data.index(b"\n", len(data) // 2) + 1
        gtf.write_bytes(gzip.compress(data[:half]) + gzip.compress(data[half:]))
    elif compress == "bgzf":
        pysam.tabix_compress(str(plain), str(gtf), force=True)
    else:
        gtf = plain
    (row,) = run_rna(
        tmp_path, write_vcf(tmp_path, [(POS + 1, ref[POS], alt)]), bam, fa, gtf, extra=extra
    )
    return row


def _invoke(d, out, *extra):
    """Re-run `gbcms rna` on the files `_run` wrote in `d`, for its output text
    (unwrapped: Rich wraps log lines at the terminal width)."""
    out.mkdir()
    args = ["rna", "-v", str(d / "variants.vcf"), "-b", f"S:{d / 'rna.s.bam'}"]
    args += ["-f", str(d / "ref.fasta"), "-o", str(out), "--gtf", str(d / "plain.gtf"), *extra]
    result = runner.invoke(app, args, env={"COLUMNS": "3000"})
    assert result.exit_code == 0, result.output
    return result.output


# ── Input: .gtf.gz ───────────────────────────────────────────────────────────


@pytest.mark.parametrize("compress", ["gzip", "bgzf", "members"])
def test_gzipped_gtf_loads_like_the_plain_file(tmp_path, compress):
    """Every output field of a .gtf.gz run equals the plain run's. The exons sit
    past ~250 KB of other lines, so a decoder that stopped after the first BGZF
    block or gzip member would lose them."""
    text = _gtf_text(filler=5000)
    assert len(text) > 3 * 65536
    plain = _run(tmp_path / "plain", text)
    packed = _run(tmp_path / "packed", text, name="gene.gtf.gz", compress=compress)
    assert plain["transcript_read_counts"].startswith("T1:")
    assert packed == plain


# ── Grammar ──────────────────────────────────────────────────────────────────


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
    # noodles rejected an unquoted last value with no ';' (GTF2.2 ends every
    # attribute with one); "." strips to an empty ID like "" does.
    "unquoted-final-value": f"chr1\tTEST\texon\t{POS + 4}\t{POS + 30}\t.\t+\t.\t"
    'gene_id "G9"; transcript_id "T9"; exon_number 1\n',
    "dot-transcript-id": f"chr1\tTEST\texon\t{POS + 4}\t{POS + 30}\t.\t+\t.\t"
    'gene_id "G9"; transcript_id ".";\n',
}


@pytest.mark.parametrize("kind", sorted(_JUNK))
def test_malformed_exon_lines_are_not_loaded(tmp_path, kind):
    clean = _run(tmp_path / "clean", _gtf_text())
    junk = _run(tmp_path / "junk", _gtf_text(extra=_JUNK[kind]))
    assert clean["exon_boundary_dist"] == "100"
    assert junk == clean


def test_rejected_exon_lines_are_warned_with_a_count(tmp_path):
    d = tmp_path / "junk"
    _run(d, _gtf_text(extra="".join(_JUNK[k] for k in sorted(_JUNK))))
    hit = [m for m in _invoke(d, tmp_path / "again").splitlines() if "malformed" in m]
    assert hit, "a warning names the rejected lines"
    assert re.search(rf"\b{len(_JUNK)} exon lines\b", hit[0]), hit[0]
    assert re.search(r"line \d+", hit[0]), hit[0]


def test_a_non_utf8_exon_line_is_rejected_and_warned(tmp_path):
    """The old parser stopped the run on a non-UTF-8 line; dropping it silently
    would derive an intron through the missing exon. It is a rejected exon line."""
    lines = _gtf_text().encode().splitlines(keepends=True)
    lines[1] = lines[1].replace(b'transcript_id "T1";', b'transcript_id "T1"; note "caf\xe9";')
    d = tmp_path / "latin1"
    _run(d, b"".join(lines))
    hit = [m for m in _invoke(d, tmp_path / "again").splitlines() if "malformed" in m]
    assert hit and re.search(r"\b1 exon line\b", hit[0]) and "UTF-8" in hit[0], hit


def test_an_index_emptied_by_rejections_says_so(tmp_path):
    d = tmp_path / "all-junk"
    _run(d, "".join(_JUNK[k] for k in ("start-after-end", "end-past-i32")))
    out = _invoke(d, tmp_path / "again")
    empty = [m for m in out.splitlines() if "inert" in m]
    assert empty and "rejected" in empty[0], empty
    assert "wrong file" not in out and "lacks these contigs" not in out


# ── The cache is deprecated ──────────────────────────────────────────────────


def test_gtf_cache_dir_is_accepted_ignored_and_warned(tmp_path):
    cache = tmp_path / "cache"
    clean = _run(tmp_path / "clean", _gtf_text())
    d = tmp_path / "flagged"
    flagged = _run(d, _gtf_text(), extra=("--gtf-cache-dir", str(cache)))
    assert flagged == clean
    assert not cache.exists() or not any(cache.iterdir())
    output = _invoke(d, tmp_path / "again", "--gtf-cache-dir", str(cache))
    assert "deprecated" in output.lower(), output


def test_no_noodles_gtf_or_bincode_dependency():
    cargo = (ROOT / "rust/Cargo.toml").read_text()
    deps = re.search(r"^\[dependencies\]\n(.*?)(?=^\[|\Z)", cargo, re.M | re.S).group(1)
    names = set(re.findall(r"^([A-Za-z0-9_-]+)\s*=", deps, re.M))
    assert "noodles-gtf" not in names
    assert "bincode" not in names
    assert not (ROOT / "rust/src/annotation/cache.rs").exists()


def test_the_binding_has_no_cache():
    stub = (ROOT / "src/gbcms/_rs.pyi").read_text()
    assert "gtf_cache_dir" not in stub
    assert "def build_gtf_cache" not in stub
    assert not hasattr(_rs, "build_gtf_cache")


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
