"""Exact speed-ups: the same outputs, less work.

Contracts:
- the BAM properties line (sampled base qualities, hard clips, duplicates) is the line the
  per-record Counter loop gave, on a BAM with binned and unbinned qualities, hard clips,
  secondary records, duplicates and a record without qualities;
- a run parses the GTF and loads the RNA editing database once: the first BAM's call parses
  and loads them, each later call reuses them and says so at INFO, the parse's warnings are said
  again, and every BAM's row equals a one-BAM run's row;
- a GTF that changes between two counting calls of one process is read again, never served
  from the once-per-run index (the same path, another length), and the same relative GTF name read
  from another working folder is that folder's file (the same length and modification time);
- importing ``gbcms`` (and its CLI) does not load polars, and the lazily read names still work.
"""

import glob
import logging
import os
import subprocess
import sys
from collections import Counter

import pysam
import pytest
from helpers import read_maf_output
from test_gtf_loading_group7_contract import _JUNK, POS, _gtf_text, _run
from typer.testing import CliRunner

from gbcms import pipeline as pl_mod
from gbcms.cli import app


def _old_properties_line(bam, min_baseq, sample=20_000):
    """The previous algorithm, kept here as the reference: a Counter over each primary
    record's qualities, merged record by record."""
    with pysam.AlignmentFile(str(bam)) as fh:
        records = dups = primaries = hard = 0
        quals: dict[int, int] = {}
        for i, r in enumerate(fh):
            if i >= sample:
                break
            records += 1
            dups += bool(r.flag & 0x400)
            if r.flag & 0x904:
                continue
            primaries += 1
            if r.cigartuples and any(op == 5 for op, _ in r.cigartuples):
                hard += 1
            q = r.query_qualities
            if q is not None:
                for value, n in Counter(q).items():
                    quals[value] = quals.get(value, 0) + n
    return records, dups, primaries, hard, dict(sorted(quals.items()))


def _write_mixed_bam(tmp_path, qual_values):
    path = tmp_path / "props.bam"
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 5000, "SN": "chr1"}]}
    raw = tmp_path / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as out:
        for i in range(60):
            a = pysam.AlignedSegment()
            a.query_name = f"r{i}"
            a.reference_id = 0
            a.reference_start = 10 + i * 5
            a.mapping_quality = 60
            n = 50
            a.query_sequence = "ACGT" * 12 + "AC"
            if i % 7 == 0:
                a.cigartuples = [(5, 10), (0, n)]  # hard-clipped
            else:
                a.cigartuples = [(0, n)]
            a.flag = 0
            if i % 5 == 0:
                a.flag |= 0x400  # duplicate
            if i % 11 == 0:
                a.flag |= 0x100  # secondary: not a primary
            if i == 13:
                a.query_qualities = None  # no base qualities
            else:
                a.query_qualities = [qual_values[(i + j) % len(qual_values)] for j in range(n)]  # type: ignore[assignment]
            out.write(a)
    pysam.sort("-o", str(path), str(raw))
    return path


@pytest.mark.parametrize("chunk", [None, 7], ids=["one-piece", "7-byte-pieces"])
@pytest.mark.parametrize(
    "qual_values",
    [[37], [2, 12, 23, 37], [0, 1, 7, 20, 25, 30, 33, 37, 40, 41, 60, 93]],
    ids=["one-value", "binned", "unbinned"],
)
def test_bam_properties_line_equals_the_counter_loop(
    tmp_path, caplog, monkeypatch, qual_values, chunk
):
    """The line is the same whether the sampled qualities are tallied in one piece or in small
    pieces (the bound that keeps a long-read BAM's sample from holding every quality)."""
    if chunk:
        monkeypatch.setattr(pl_mod, "BAM_QUAL_CHUNK", chunk)
    bam = _write_mixed_bam(tmp_path, qual_values)
    records, dups, primaries, hard, quals = _old_properties_line(bam, 20)
    bases = sum(quals.values())
    below = sum(n for v, n in quals.items() if v < 20) / bases
    alphabet = sorted(quals)
    bq = (
        f"{len(alphabet)} value{'' if len(alphabet) == 1 else 's'}"
        f" ({alphabet[0]}-{alphabet[-1]})" + (" (binned)" if 0 < len(alphabet) <= 8 else "")
    )
    want = (
        f"BAM S: first {records} records: duplicates flagged {100 * dups / records:.1f}%; "
        f"base qualities {bq}; bases below --min-baseq 20: {100 * below:.1f}%; ALT contigs 0; "
        f"hard-clipped primaries {100 * hard / primaries:.1f}%"
    )
    with caplog.at_level(logging.INFO, logger=pl_mod.logger.name):
        pl_mod.log_bam_properties("S", bam, min_baseq=20)
    lines = [r.getMessage() for r in caplog.records if r.getMessage().startswith("BAM S: first")]
    assert lines == [want]


def _cli_rna(d, out, bams):
    """`gbcms rna` in a fresh process (so no parse is kept from an earlier call) on the files
    `_run` wrote in `d`, with the GTF and an editing database: its log and its rows by sample."""
    out.mkdir()
    cmd = [
        sys.executable,
        "-c",
        "from gbcms.cli import app; app()",
        "rna",
        "-v",
        str(d / "variants.vcf"),
    ]
    for tag in bams:
        cmd += ["-b", f"{tag}:{d / 'rna.s.bam'}"]
    cmd += [
        "-f",
        str(d / "ref.fasta"),
        "-o",
        str(out),
        "--format",
        "maf",
        "--gtf",
        str(d / "plain.gtf"),
    ]
    cmd += ["--rna-editing-db", str(d / "edits.txt")]
    done = subprocess.run(
        cmd, capture_output=True, text=True, env={**os.environ, "COLUMNS": "3000"}
    )
    assert done.returncode == 0, done.stderr[-2000:]
    rows = {}
    for f in glob.glob(str(out / "*.maf")):
        for r in read_maf_output(f):
            rows[r["Tumor_Sample_Barcode"]] = r
    return done.stdout + done.stderr, rows


def test_a_run_parses_the_gtf_once_and_says_so(tmp_path):
    d = tmp_path / "files"
    one = _run(d, _gtf_text(extra=_JUNK["start-after-end"]))
    alt = one["Tumor_Seq_Allele2"]
    ref_base = one["Reference_Allele"]
    (d / "edits.txt").write_text(
        f"Accession\tRegion\tPosition\tRef\tEd\tStrand\nE1\tchr1\t{POS + 1}\t{ref_base}\t{alt}\t+\n"
    )
    log1, rows1 = _cli_rna(d, tmp_path / "one", ["A"])
    log2, rows2 = _cli_rna(d, tmp_path / "two", ["A", "B"])

    # One BAM: the parse and the load, as before, and no reuse.
    assert log1.count("Loading GTF annotation from") == 1
    assert log1.count("GTF annotation: built index from") == 1
    assert "reusing the index" not in log1 and "Reusing " not in log1
    # Two BAMs: one parse and one load, then one reuse of each, said at INFO.
    assert log2.count("Loading GTF annotation from") == 1
    assert log2.count("GTF annotation loaded in") == 1
    assert log2.count("GTF annotation: built index from") == 1
    assert log2.count("GTF annotation: reusing the index parsed earlier in this run from") == 1
    assert log2.count("Reusing 1 RNA editing sites loaded earlier in this run from") == 1
    # The parse's warning (the rejected exon line) is said for each BAM.
    assert sum("rejected as malformed" in ln for ln in log1.splitlines()) == 1
    assert sum("rejected as malformed" in ln for ln in log2.splitlines()) == 2

    # Every BAM's row is the one-BAM row: the reused index and sites give the same answers.
    assert rows1["A"]["rna_editing_site"] in ("True", "true", "1"), rows1["A"]["rna_editing_site"]
    want = {k: v for k, v in rows1["A"].items() if k != "Tumor_Sample_Barcode"}
    for tag in ("A", "B"):
        got = {k: v for k, v in rows2[tag].items() if k != "Tumor_Sample_Barcode"}
        assert got == want, tag
        r = rows2[tag]
        assert int(r["total_count"]) >= int(r["ref_count"]) + int(r["alt_count"])
        assert int(r["total_count_fragment"]) >= int(r["ref_count_fragment"]) + int(
            r["alt_count_fragment"]
        )
        assert int(r["ref_count"]) == int(r["ref_count_forward"]) + int(r["ref_count_reverse"])
        assert int(r["alt_count"]) == int(r["alt_count_forward"]) + int(r["alt_count_reverse"])


def test_a_gtf_changed_between_two_calls_is_read_again(tmp_path):
    """The once-per-run index is keyed by the file: the same path with other content (so
    another length) gives the other file's answers."""
    d = tmp_path / "same_path"
    before = _run(d, _gtf_text())
    near = f'chr1\tTEST\texon\t{POS + 5}\t{POS + 40}\t.\t+\t.\tgene_id "G9"; transcript_id "T9";\n'
    after = _run(d, _gtf_text(extra=near))
    fresh = _run(tmp_path / "fresh_path", _gtf_text(extra=near))
    assert before["exon_boundary_dist"] == "100"
    assert after["exon_boundary_dist"] != "100"
    assert after == fresh
    again = _run(d, _gtf_text())
    assert again == before


def test_the_same_relative_gtf_name_in_another_folder_is_that_folders_file(tmp_path, monkeypatch):
    """The once-per-run index is keyed by the file a name resolves to. A Python caller that
    moves to another working folder and passes the same relative GTF name, a file of the same
    length and modification time, gets that folder's answers."""
    a, b = tmp_path / "a", tmp_path / "b"
    near_a = (
        f'chr1\tTEST\texon\t{POS + 5}\t{POS + 40}\t.\t+\t.\tgene_id "G9"; transcript_id "T9";\n'
    )
    near_b = (
        f'chr1\tTEST\texon\t{POS + 6}\t{POS + 40}\t.\t+\t.\tgene_id "G9"; transcript_id "T9";\n'
    )
    assert len(near_a) == len(near_b)
    want_a = _run(a, _gtf_text(extra=near_a))
    want_b = _run(b, _gtf_text(extra=near_b))
    assert want_a["exon_boundary_dist"] != want_b["exon_boundary_dist"]
    st = (a / "plain.gtf").stat()
    os.utime(b / "plain.gtf", ns=(st.st_atime_ns, st.st_mtime_ns))
    assert (b / "plain.gtf").stat().st_size == st.st_size

    def relative_run(folder, out):
        monkeypatch.chdir(folder)
        args = ["rna", "-v", "variants.vcf", "-b", "S:rna.s.bam", "-f", "ref.fasta", "-o", out]
        result = CliRunner().invoke(app, args + ["--format", "maf", "--gtf", "plain.gtf"])
        assert result.exit_code == 0, result.output
        (row,) = read_maf_output(glob.glob(str(folder / out / "*.maf"))[0])
        return row

    assert relative_run(a, "rel") == want_a
    assert relative_run(b, "rel") == want_b, "the first folder's index was reused"


def test_importing_gbcms_does_not_load_polars():
    code = (
        "import sys, gbcms, gbcms.cli, gbcms.io\n"
        "assert 'polars' not in sys.modules, 'polars was loaded at import'\n"
        "from gbcms import merge_mafs\n"
        "from gbcms.io import read_maf, read_parquet, scan_maf, write_maf\n"
        "assert callable(merge_mafs) and all(callable(f) for f in (read_maf, read_parquet, scan_maf, write_maf))\n"
        "assert 'polars' in sys.modules\n"
        "import gbcms.io\n"
        "try:\n"
        "    gbcms.io.no_such_name\n"
        "except AttributeError:\n"
        "    pass\n"
        "else:\n"
        "    raise SystemExit('a missing name must raise AttributeError')\n"
    )
    done = subprocess.run([sys.executable, "-I", "-c", code], capture_output=True, text=True)
    assert done.returncode == 0, done.stderr
