"""
Tests for centralized Polars batch I/O helpers (``gbcms.io.batch``).

These tests verify the batch read/write functions used by the merge
engine and mFSD report generator. They do NOT test the streaming
``csv``-based I/O in ``io/input.py`` and ``io/output.py``.
"""

import polars as pl
import pytest

from gbcms.io.batch import read_maf, read_parquet, scan_maf, write_maf

# ── Helpers ──────────────────────────────────────────────────────────────────


def _write_test_maf(path, rows, columns=None):
    """Write a minimal test MAF file.

    Args:
        path: Output file path.
        rows: List of dicts with column values.
        columns: Optional explicit column order. Defaults to keys of first row.
    """
    if not rows:
        columns = columns or ["Chromosome", "Start_Position"]
        path.write_text("\t".join(columns) + "\n")
        return

    columns = columns or list(rows[0].keys())
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(str(row.get(c, "")) for c in columns))
    path.write_text("\n".join(lines) + "\n")


# ── Tests ────────────────────────────────────────────────────────────────────


def test_read_maf(tmp_path):
    """read_maf loads a tab-separated MAF file with all columns as strings."""
    maf = tmp_path / "test.maf"
    _write_test_maf(
        maf,
        [
            {"Chromosome": "chr1", "Start_Position": "100", "ref_count": "10"},
            {"Chromosome": "chr2", "Start_Position": "200", "ref_count": "20"},
        ],
    )

    df = read_maf(maf)
    assert df.height == 2
    assert df.width == 3
    # All columns should be strings (infer_schema_length=0)
    assert all(dtype == pl.Utf8 for dtype in df.dtypes), f"Expected all Utf8, got {df.dtypes}"
    assert df["Chromosome"].to_list() == ["chr1", "chr2"]


def test_read_maf_skips_comments(tmp_path):
    """read_maf skips lines starting with '#'."""
    maf = tmp_path / "commented.maf"
    maf.write_text(
        "# version 2.4\n" "# Hugo_Symbol filter\n" "Chromosome\tStart_Position\n" "chr1\t100\n"
    )

    df = read_maf(maf)
    assert df.height == 1
    assert df["Chromosome"].to_list() == ["chr1"]


def test_read_maf_missing_file(tmp_path):
    """read_maf raises FileNotFoundError with the path in the message."""
    missing = tmp_path / "nonexistent.maf"

    with pytest.raises(FileNotFoundError, match=str(missing)):
        read_maf(missing)


def test_scan_maf_lazy(tmp_path):
    """scan_maf returns a LazyFrame that is not materialized until .collect()."""
    maf = tmp_path / "lazy.maf"
    _write_test_maf(
        maf,
        [{"Chromosome": "chr1", "Start_Position": "100"}],
    )

    lf = scan_maf(maf)
    assert isinstance(lf, pl.LazyFrame)

    # Collecting materializes the data
    df = lf.collect()
    assert df.height == 1
    assert df["Chromosome"].to_list() == ["chr1"]


def test_scan_maf_missing_file(tmp_path):
    """scan_maf raises FileNotFoundError for missing files."""
    missing = tmp_path / "nonexistent.maf"

    with pytest.raises(FileNotFoundError, match=str(missing)):
        scan_maf(missing)


def test_read_parquet(tmp_path):
    """read_parquet loads a Parquet file with native types."""
    parquet = tmp_path / "test.parquet"
    df_write = pl.DataFrame({"chrom": ["chr1", "chr2"], "pos": [100, 200]})
    df_write.write_parquet(parquet)

    df = read_parquet(parquet)
    assert df.height == 2
    assert df["chrom"].to_list() == ["chr1", "chr2"]
    assert df["pos"].to_list() == [100, 200]


def test_read_parquet_missing_file(tmp_path):
    """read_parquet raises FileNotFoundError with the path in the message."""
    missing = tmp_path / "nonexistent.parquet"

    with pytest.raises(FileNotFoundError, match=str(missing)):
        read_parquet(missing)


def test_write_maf_roundtrip(tmp_path):
    """write_maf writes TSV that read_maf can round-trip."""
    original = pl.DataFrame(
        {
            "Chromosome": ["chr1", "chr2"],
            "Start_Position": ["100", "200"],
            "ref_count": ["10", "20"],
        }
    )

    out = tmp_path / "output.maf"
    write_maf(original, out)

    assert out.exists()
    roundtrip = read_maf(out)
    assert roundtrip.height == 2
    assert roundtrip["Chromosome"].to_list() == ["chr1", "chr2"]
    assert roundtrip["ref_count"].to_list() == ["10", "20"]
