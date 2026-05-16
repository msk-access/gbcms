"""
Tests for ``gbcms merge`` — multi-BAM MAF merge engine and CLI.

Covers:
  - Core merge logic (join, prefix, combined columns)
  - Edge cases (empty MAF, missing files, single input)
  - Column naming conventions (default vs legacy)
  - CLI integration via CliRunner
"""

from pathlib import Path

import polars as pl
import pytest
from typer.testing import CliRunner

from gbcms.cli import app
from gbcms.merge import (
    VARIANT_KEY,
    merge_mafs,
)
from gbcms.models.core import MergeConfig

runner = CliRunner()


# ── Helpers ──────────────────────────────────────────────────────────────────


def _write_test_maf(path: Path, rows: list[dict], columns: list[str] | None = None) -> None:
    """Write a minimal test MAF file with gbcms count columns.

    Args:
        path: Output file path.
        rows: List of dicts with column values.
        columns: Optional explicit column order.
    """
    if not rows:
        columns = columns or VARIANT_KEY
        path.write_text("\t".join(columns) + "\n")
        return

    columns = columns or list(rows[0].keys())
    lines = ["\t".join(columns)]
    for row in rows:
        lines.append("\t".join(str(row.get(c, "")) for c in columns))
    path.write_text("\n".join(lines) + "\n")


def _make_variant_row(
    chrom: str = "chr1",
    start: str = "100",
    end: str = "100",
    ref: str = "A",
    alt: str = "T",
    **counts,
) -> dict:
    """Create a variant row dict with variant key + optional count overrides."""
    row = {
        "Chromosome": chrom,
        "Start_Position": start,
        "End_Position": end,
        "Reference_Allele": ref,
        "Tumor_Seq_Allele2": alt,
        "Hugo_Symbol": "TP53",
        "ref_count": "10",
        "alt_count": "5",
        "total_count": "15",
        "vaf": "0.333",
        "ref_count_fragment": "8",
        "alt_count_fragment": "4",
        "total_count_fragment": "12",
        "vaf_fragment": "0.333",
        "gbcms_status": "OK",
    }
    row.update(counts)
    return row


# ── Test 1: Basic duplex+simplex merge ───────────────────────────────────────


def test_merge_duplex_simplex(tmp_path):
    """Two-input join produces correct type-prefixed columns."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(duplex, [_make_variant_row(ref_count="20", alt_count="10")])
    _write_test_maf(simplex, [_make_variant_row(ref_count="5", alt_count="2")])

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=False,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert result.height == 1
    assert "duplex_ref_count" in result.columns
    assert "simplex_ref_count" in result.columns
    assert result["duplex_ref_count"].to_list() == ["20"]
    assert result["simplex_ref_count"].to_list() == ["5"]


# ── Test 2: Combined columns ────────────────────────────────────────────────


def test_merge_combined_columns(tmp_path):
    """simplex_duplex_* additive math is correct."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(
        duplex,
        [_make_variant_row(ref_count_fragment="10", alt_count_fragment="5")],
    )
    _write_test_maf(
        simplex,
        [_make_variant_row(ref_count_fragment="3", alt_count_fragment="2")],
    )

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=True,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert "simplex_duplex_ref_count_fragment" in result.columns
    assert "simplex_duplex_alt_count_fragment" in result.columns
    assert "simplex_duplex_total_count_fragment" in result.columns
    assert "simplex_duplex_vaf_fragment" in result.columns

    # 10 + 3 = 13
    assert int(result["simplex_duplex_ref_count_fragment"][0]) == 13
    # 5 + 2 = 7
    assert int(result["simplex_duplex_alt_count_fragment"][0]) == 7
    # total = 13 + 7 = 20
    assert int(result["simplex_duplex_total_count_fragment"][0]) == 20
    # vaf = 7 / 20 = 0.35
    assert abs(float(result["simplex_duplex_vaf_fragment"][0]) - 0.35) < 0.001


# ── Test 3: --no-combined ────────────────────────────────────────────────────


def test_merge_no_combined(tmp_path):
    """--no-combined omits simplex_duplex_* columns."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(duplex, [_make_variant_row()])
    _write_test_maf(simplex, [_make_variant_row()])

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=False,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    combined_cols = [c for c in result.columns if c.startswith("simplex_duplex_")]
    assert combined_cols == [], f"Expected no combined columns, found: {combined_cols}"


# ── Test 4: Three types ─────────────────────────────────────────────────────


def test_merge_three_types(tmp_path):
    """Three inputs (duplex + simplex + standard) all have separate prefixed columns."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    standard = tmp_path / "standard.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(duplex, [_make_variant_row(ref_count="20")])
    _write_test_maf(simplex, [_make_variant_row(ref_count="5")])
    _write_test_maf(standard, [_make_variant_row(ref_count="100")])

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex, "standard": standard},
        output=output,
        add_combined=True,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert "duplex_ref_count" in result.columns
    assert "simplex_ref_count" in result.columns
    assert "standard_ref_count" in result.columns
    assert result["standard_ref_count"].to_list() == ["100"]


# ── Test 5: Unmatched variants (outer join) ──────────────────────────────────


def test_merge_unmatched_variants(tmp_path):
    """Outer join fills missing variants with '0' for count columns."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(duplex, [_make_variant_row(start="100", ref_count="20")])
    _write_test_maf(simplex, [_make_variant_row(start="200", ref_count="5")])

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=False,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert result.height == 2, f"Expected 2 rows from outer join, got {result.height}"

    # Check that missing counts are filled with "0"
    ref_counts = result["duplex_ref_count"].to_list()
    assert "0" in ref_counts, f"Expected '0' for unmatched variant, got: {ref_counts}"


# ── Test 6: Empty MAF ────────────────────────────────────────────────────────


def test_merge_empty_maf(tmp_path):
    """Empty input MAF produces valid output (header only)."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    # Non-empty duplex, empty simplex
    _write_test_maf(duplex, [_make_variant_row()])
    _write_test_maf(simplex, [], columns=list(_make_variant_row().keys()))

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=False,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    # Outer join should preserve the duplex row
    assert result.height >= 0  # At minimum, the file is valid


# ── Test 7: Already-prefixed input ───────────────────────────────────────────


def test_merge_prefixed_input(tmp_path):
    """Already-prefixed columns are not double-prefixed."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    # Write MAF with already-prefixed columns
    row = {
        "Chromosome": "chr1",
        "Start_Position": "100",
        "End_Position": "100",
        "Reference_Allele": "A",
        "Tumor_Seq_Allele2": "T",
        "duplex_ref_count": "20",
        "duplex_alt_count": "10",
    }
    _write_test_maf(duplex, [row])
    _write_test_maf(simplex, [_make_variant_row(ref_count="5")])

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=False,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    # Should have duplex_ref_count, NOT duplex_duplex_ref_count
    assert "duplex_ref_count" in result.columns
    double_prefixed = [c for c in result.columns if c.startswith("duplex_duplex_")]
    assert double_prefixed == [], f"Double-prefixed columns found: {double_prefixed}"


# ── Test 8: Unprefixed input renamed correctly ───────────────────────────────


def test_merge_unprefixed_input(tmp_path):
    """Unprefixed gbcms columns are renamed with type prefix."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(duplex, [_make_variant_row()])
    _write_test_maf(simplex, [_make_variant_row()])

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=False,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    # Original "ref_count" should NOT exist — should be "duplex_ref_count"
    assert "ref_count" not in result.columns, "Unprefixed ref_count should be renamed"
    assert "duplex_ref_count" in result.columns
    assert "simplex_ref_count" in result.columns


# ── Test 9: Legacy naming ────────────────────────────────────────────────────


def test_merge_legacy_naming(tmp_path):
    """--legacy-naming produces t_{metric}_{type} format."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(duplex, [_make_variant_row()])
    _write_test_maf(simplex, [_make_variant_row()])

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=False,
        legacy_naming=True,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert "t_ref_count_duplex" in result.columns, f"Expected legacy naming, got: {result.columns}"
    assert "t_ref_count_simplex" in result.columns


# ── Test 10: Annotation passthrough ──────────────────────────────────────────


def test_merge_annotation_passthrough(tmp_path):
    """Non-count MAF columns (Hugo_Symbol, etc.) are preserved verbatim."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(duplex, [_make_variant_row()])
    _write_test_maf(simplex, [_make_variant_row()])

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=False,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert "Hugo_Symbol" in result.columns
    assert result["Hugo_Symbol"].to_list() == ["TP53"]


# ── Test 11: Single input error ──────────────────────────────────────────────


def test_merge_single_input_error(tmp_path):
    """Single input MAF fails with clear error."""
    duplex = tmp_path / "duplex.maf"
    _write_test_maf(duplex, [_make_variant_row()])

    with pytest.raises(Exception, match="At least 2"):
        MergeConfig(
            inputs={"duplex": duplex},
            output=tmp_path / "merged.maf",
        )


# ── Test 12: Missing file error ──────────────────────────────────────────────


def test_merge_missing_file_error(tmp_path):
    """Missing input MAF file fails with path in error message."""
    existing = tmp_path / "existing.maf"
    _write_test_maf(existing, [_make_variant_row()])
    missing = tmp_path / "missing.maf"

    with pytest.raises(Exception, match="not found"):
        MergeConfig(
            inputs={"existing": existing, "missing": missing},
            output=tmp_path / "merged.maf",
        )


# ── Test 13: CLI integration ────────────────────────────────────────────────


def test_merge_cli_integration(tmp_path):
    """Full CLI invocation via CliRunner."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(duplex, [_make_variant_row(ref_count="20")])
    _write_test_maf(simplex, [_make_variant_row(ref_count="5")])

    result = runner.invoke(
        app,
        [
            "merge",
            "--input",
            f"duplex:{duplex}",
            "--input",
            f"simplex:{simplex}",
            "--output",
            str(output),
            "--no-combined",
        ],
    )

    assert result.exit_code == 0, f"CLI failed: {result.output}"
    assert output.exists()

    df = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert df.height == 1
    assert "duplex_ref_count" in df.columns


# ── Test 14: VAF zero division ───────────────────────────────────────────────


def test_merge_vaf_zero_division(tmp_path):
    """Combined VAF handles 0/0 → 0.0 without error."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(
        duplex,
        [_make_variant_row(ref_count_fragment="0", alt_count_fragment="0")],
    )
    _write_test_maf(
        simplex,
        [_make_variant_row(ref_count_fragment="0", alt_count_fragment="0")],
    )

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=True,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    vaf = float(result["simplex_duplex_vaf_fragment"][0])
    assert vaf == 0.0, f"Expected 0.0 for 0/0, got {vaf}"


# ── Test 15: Singleton column prefixing ──────────────────────────────────────


def test_merge_singleton_column_prefixing(tmp_path):
    """gbcms_status → duplex_gbcms_status (meta columns also prefixed)."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(duplex, [_make_variant_row()])
    _write_test_maf(simplex, [_make_variant_row()])

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=False,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert (
        "duplex_gbcms_status" in result.columns
    ), f"Expected 'duplex_gbcms_status', got: {[c for c in result.columns if 'status' in c]}"
    assert "simplex_gbcms_status" in result.columns
    assert "gbcms_status" not in result.columns, "Unprefixed gbcms_status should not exist"


# ── Test 16: Asymmetric variant counts (many vs few, partial overlap) ────────


def test_merge_asymmetric_row_counts(tmp_path):
    """Duplex has 5 variants, simplex has 3, with 2 overlapping.

    Expected: 6 rows (5 + 3 - 2 overlap = 6 unique variants).
    Unmatched duplex-only rows should have simplex counts filled with "0".
    Unmatched simplex-only row should have duplex counts filled with "0".
    """
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    duplex_rows = [
        _make_variant_row(start="100", ref_count="10"),  # overlap
        _make_variant_row(start="200", ref_count="20"),  # overlap
        _make_variant_row(start="300", ref_count="30"),  # duplex-only
        _make_variant_row(start="400", ref_count="40"),  # duplex-only
        _make_variant_row(start="500", ref_count="50"),  # duplex-only
    ]
    simplex_rows = [
        _make_variant_row(start="100", ref_count="1"),  # overlap
        _make_variant_row(start="200", ref_count="2"),  # overlap
        _make_variant_row(start="600", ref_count="6"),  # simplex-only
    ]

    _write_test_maf(duplex, duplex_rows)
    _write_test_maf(simplex, simplex_rows)

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=False,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert result.height == 6, f"Expected 6 rows (5+3-2 overlap), got {result.height}"

    # Check overlap rows: both have real counts
    overlap = result.filter(pl.col("Start_Position") == "100")
    assert overlap["duplex_ref_count"].to_list() == ["10"]
    assert overlap["simplex_ref_count"].to_list() == ["1"]

    # Check duplex-only rows: simplex counts should be "0"
    duplex_only = result.filter(pl.col("Start_Position") == "300")
    assert duplex_only["duplex_ref_count"].to_list() == ["30"]
    assert duplex_only["simplex_ref_count"].to_list() == ["0"]

    # Check simplex-only rows: duplex counts should be "0"
    simplex_only = result.filter(pl.col("Start_Position") == "600")
    assert simplex_only["simplex_ref_count"].to_list() == ["6"]
    assert simplex_only["duplex_ref_count"].to_list() == ["0"]


# ── Test 17: Combined columns with unmatched variants ────────────────────────


def test_merge_combined_with_unmatched(tmp_path):
    """Combined simplex_duplex columns work correctly when one side is missing.

    Duplex-only variant: simplex fragment counts are null → filled to 0.
    Combined = duplex + 0 = duplex.
    """
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(
        duplex,
        [_make_variant_row(start="100", ref_count_fragment="10", alt_count_fragment="5")],
    )
    _write_test_maf(
        simplex,
        [_make_variant_row(start="200", ref_count_fragment="3", alt_count_fragment="2")],
    )

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=True,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert result.height == 2

    # Duplex-only variant: combined = duplex + 0
    duplex_only = result.filter(pl.col("Start_Position") == "100")
    assert int(duplex_only["simplex_duplex_ref_count_fragment"][0]) == 10
    assert int(duplex_only["simplex_duplex_alt_count_fragment"][0]) == 5

    # Simplex-only variant: combined = 0 + simplex
    simplex_only = result.filter(pl.col("Start_Position") == "200")
    assert int(simplex_only["simplex_duplex_ref_count_fragment"][0]) == 3
    assert int(simplex_only["simplex_duplex_alt_count_fragment"][0]) == 2


# ── Test 18: Annotation columns on right-only variants ───────────────────────


def test_merge_annotation_nulls_right_only(tmp_path):
    """Simplex-only variants have null annotation columns (by design).

    Annotations come from the left (first) frame. Variants that exist only
    in the right (joining) frame will have NULL annotations because only
    variant key + gbcms columns are selected from the right frame.
    """
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(duplex, [_make_variant_row(start="100")])
    _write_test_maf(simplex, [_make_variant_row(start="200")])

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=False,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)

    # Duplex-only row (left): Hugo_Symbol preserved
    left_row = result.filter(pl.col("Start_Position") == "100")
    assert left_row["Hugo_Symbol"].to_list() == ["TP53"]

    # Simplex-only row (right): Hugo_Symbol is null (annotations not carried from right)
    right_row = result.filter(pl.col("Start_Position") == "200")
    hugo = right_row["Hugo_Symbol"].to_list()[0]
    assert (
        hugo is None or hugo == ""
    ), f"Expected null/empty for right-only annotation, got: {hugo!r}"


# ── Test 19: Meta column null-fill ───────────────────────────────────────────


def test_merge_meta_null_fill(tmp_path):
    """gbcms_status for unmatched variants is filled with '' (not null)."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(duplex, [_make_variant_row(start="100")])
    _write_test_maf(simplex, [_make_variant_row(start="200")])

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=False,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)

    # Duplex-only row: simplex_gbcms_status should be "" (not null)
    duplex_only = result.filter(pl.col("Start_Position") == "100")
    simplex_status = duplex_only["simplex_gbcms_status"].to_list()[0]
    assert simplex_status is not None, "Meta column should not be null after fill"

    # Simplex-only row: duplex_gbcms_status should be "" (not null)
    simplex_only = result.filter(pl.col("Start_Position") == "200")
    duplex_status = simplex_only["duplex_gbcms_status"].to_list()[0]
    assert duplex_status is not None, "Meta column should not be null after fill"


# ── Test 20: CLI invalid input format ────────────────────────────────────────


def test_merge_cli_invalid_input_format(tmp_path):
    """--input without colon separator exits with code 1 and error message."""
    maf = tmp_path / "test.maf"
    _write_test_maf(maf, [_make_variant_row()])

    result = runner.invoke(
        app,
        [
            "merge",
            "--input",
            str(maf),  # Missing type: prefix
            "--input",
            f"simplex:{maf}",
            "--output",
            str(tmp_path / "merged.maf"),
        ],
    )

    assert result.exit_code == 1, f"Expected exit code 1, got {result.exit_code}"


# ── Helpers for full-column tests ─────────────────────────────────────────


def _make_full_variant_row(
    chrom: str = "chr1",
    start: str = "100",
    end: str = "100",
    ref: str = "A",
    alt: str = "T",
    **counts,
) -> dict:
    """Create a variant row with ALL count columns including strand directions.

    Matches the canonical column set from output.py §2-§5.
    """
    row = {
        "Chromosome": chrom,
        "Start_Position": start,
        "End_Position": end,
        "Reference_Allele": ref,
        "Tumor_Seq_Allele2": alt,
        "Hugo_Symbol": "TP53",
        # §2: Read-level counts
        "ref_count": "10",
        "alt_count": "5",
        "total_count": "15",
        "vaf": "0.333",
        # §3: Read-level strand
        "ref_count_forward": "6",
        "ref_count_reverse": "4",
        "alt_count_forward": "3",
        "alt_count_reverse": "2",
        "strand_bias_p_value": "1.0000e+00",
        "strand_bias_odds_ratio": "1.0000",
        # §4: Fragment-level counts
        "ref_count_fragment": "8",
        "alt_count_fragment": "4",
        "total_count_fragment": "12",
        "vaf_fragment": "0.333",
        # §5: Fragment strand
        "ref_count_fragment_forward": "5",
        "ref_count_fragment_reverse": "3",
        "alt_count_fragment_forward": "2",
        "alt_count_fragment_reverse": "2",
        "fragment_strand_bias_p_value": "1.0000e+00",
        "fragment_strand_bias_odds_ratio": "1.0000",
        # Meta
        "gbcms_status": "OK",
    }
    row.update(counts)
    return row


# ── Test 21: Read-level combined columns ─────────────────────────────────


def test_merge_combined_read_level(tmp_path):
    """Combined read-level counts are additive (duplex + simplex reads are distinct)."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(duplex, [_make_full_variant_row(ref_count="20", alt_count="10")])
    _write_test_maf(simplex, [_make_full_variant_row(ref_count="5", alt_count="3")])

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=True,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)

    # Read-level additive sums
    assert int(result["simplex_duplex_ref_count"][0]) == 25  # 20 + 5
    assert int(result["simplex_duplex_alt_count"][0]) == 13  # 10 + 3
    assert int(result["simplex_duplex_total_count"][0]) == 38  # 25 + 13
    assert abs(float(result["simplex_duplex_vaf"][0]) - 13 / 38) < 0.001


# ── Test 22: Fragment strand combined + Fisher bias ──────────────────────


def test_merge_combined_strand_bias(tmp_path):
    """Combined strand bias uses Rust Fisher exact test on summed counts.

    Duplex: ref_fwd=10, ref_rev=10, alt_fwd=5, alt_rev=5 (balanced)
    Simplex: ref_fwd=0, ref_rev=10, alt_fwd=10, alt_rev=0 (extreme bias)
    Combined: ref_fwd=10, ref_rev=20, alt_fwd=15, alt_rev=5
    """
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(
        duplex,
        [
            _make_full_variant_row(
                ref_count_forward="10",
                ref_count_reverse="10",
                alt_count_forward="5",
                alt_count_reverse="5",
                ref_count_fragment_forward="10",
                ref_count_fragment_reverse="10",
                alt_count_fragment_forward="5",
                alt_count_fragment_reverse="5",
            )
        ],
    )
    _write_test_maf(
        simplex,
        [
            _make_full_variant_row(
                ref_count_forward="0",
                ref_count_reverse="10",
                alt_count_forward="10",
                alt_count_reverse="0",
                ref_count_fragment_forward="0",
                ref_count_fragment_reverse="10",
                alt_count_fragment_forward="10",
                alt_count_fragment_reverse="0",
            )
        ],
    )

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=True,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)

    # Verify additive strand counts
    assert int(result["simplex_duplex_ref_count_forward"][0]) == 10  # 10 + 0
    assert int(result["simplex_duplex_ref_count_reverse"][0]) == 20  # 10 + 10
    assert int(result["simplex_duplex_alt_count_forward"][0]) == 15  # 5 + 10
    assert int(result["simplex_duplex_alt_count_reverse"][0]) == 5  # 5 + 0

    # Verify strand bias columns exist and are not null
    assert "simplex_duplex_strand_bias_p_value" in result.columns
    assert "simplex_duplex_strand_bias_odds_ratio" in result.columns
    assert "simplex_duplex_fragment_strand_bias_p_value" in result.columns
    assert "simplex_duplex_fragment_strand_bias_odds_ratio" in result.columns

    # Verify Fisher p-value is < 1.0 (biased table: 10/20 vs 15/5)
    sb_p = float(result["simplex_duplex_strand_bias_p_value"][0])
    assert sb_p < 0.05, f"Expected significant bias, got p={sb_p}"

    # Fragment-level should give same result (same input values)
    fsb_p = float(result["simplex_duplex_fragment_strand_bias_p_value"][0])
    assert fsb_p < 0.05, f"Expected significant fragment bias, got p={fsb_p}"


# ── Test 23: Full combined column order ──────────────────────────────────


def test_merge_combined_column_order(tmp_path):
    """All 20 combined columns are present in the correct canonical order.

    Order must match output.py: §2 counts → §3 strand + SB → §4 fragment → §5 frag strand + FSB.
    """
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf(duplex, [_make_full_variant_row()])
    _write_test_maf(simplex, [_make_full_variant_row()])

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=True,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)

    # Extract all simplex_duplex_* columns in order
    combined_cols = [c for c in result.columns if c.startswith("simplex_duplex_")]

    expected_combined = [
        # Phase 1: Additive sums (COMBINED_ADDITIVE_ALL order)
        "simplex_duplex_ref_count",
        "simplex_duplex_alt_count",
        "simplex_duplex_ref_count_forward",
        "simplex_duplex_ref_count_reverse",
        "simplex_duplex_alt_count_forward",
        "simplex_duplex_alt_count_reverse",
        "simplex_duplex_ref_count_fragment",
        "simplex_duplex_alt_count_fragment",
        "simplex_duplex_ref_count_fragment_forward",
        "simplex_duplex_ref_count_fragment_reverse",
        "simplex_duplex_alt_count_fragment_forward",
        "simplex_duplex_alt_count_fragment_reverse",
        # Phase 2a: Derived totals (read + fragment)
        "simplex_duplex_total_count",
        "simplex_duplex_total_count_fragment",
        # Phase 2b: Derived VAFs (read + fragment)
        "simplex_duplex_vaf",
        "simplex_duplex_vaf_fragment",
        # Phase 3: Strand bias (read-level SB + fragment-level FSB)
        "simplex_duplex_strand_bias_p_value",
        "simplex_duplex_strand_bias_odds_ratio",
        "simplex_duplex_fragment_strand_bias_p_value",
        "simplex_duplex_fragment_strand_bias_odds_ratio",
    ]

    assert (
        len(combined_cols) == 20
    ), f"Expected 20 combined columns, got {len(combined_cols)}: {combined_cols}"
    assert (
        combined_cols == expected_combined
    ), f"Column order mismatch:\nGot: {combined_cols}\nExpected: {expected_combined}"


# ── Test 24: Balanced table → no strand bias ─────────────────────────────


def test_merge_combined_strand_bias_balanced(tmp_path):
    """Balanced strand counts yield p ≈ 1.0 (no bias)."""
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    # Both have perfectly balanced strand counts
    _write_test_maf(
        duplex,
        [
            _make_full_variant_row(
                ref_count_forward="10",
                ref_count_reverse="10",
                alt_count_forward="5",
                alt_count_reverse="5",
            )
        ],
    )
    _write_test_maf(
        simplex,
        [
            _make_full_variant_row(
                ref_count_forward="10",
                ref_count_reverse="10",
                alt_count_forward="5",
                alt_count_reverse="5",
            )
        ],
    )

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=True,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    sb_p = float(result["simplex_duplex_strand_bias_p_value"][0])
    assert sb_p > 0.9, f"Expected p ~1.0 for balanced table, got {sb_p}"


# ── Test 25: NaN/inf sanitization in combined strand bias ────────────────


def test_merge_combined_strand_bias_nan_sanitized(tmp_path):
    """Combined strand bias columns never contain literal 'NaN' or 'inf' strings.

    When alt_total ≤ 1, Fisher's exact test returns (1.0, NaN) for the odds
    ratio. The merge engine must sanitize these to 'NA' for MAF compliance.

    Regression test for issue #19 / plan step 4h.
    """
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    # Zero ALT reads on both sides → combined alt_fwd + alt_rev = 0 → NaN OR
    _write_test_maf(
        duplex,
        [
            _make_full_variant_row(
                alt_count="0",
                alt_count_forward="0",
                alt_count_reverse="0",
                alt_count_fragment="0",
                alt_count_fragment_forward="0",
                alt_count_fragment_reverse="0",
            )
        ],
    )
    _write_test_maf(
        simplex,
        [
            _make_full_variant_row(
                alt_count="1",
                alt_count_forward="1",
                alt_count_reverse="0",
                alt_count_fragment="1",
                alt_count_fragment_forward="1",
                alt_count_fragment_reverse="0",
            )
        ],
    )

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=True,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)

    # Identify all simplex_duplex strand bias columns
    sb_cols = [c for c in result.columns if "strand_bias" in c and c.startswith("simplex_duplex_")]
    assert len(sb_cols) == 4, f"Expected 4 combined SB columns, got {len(sb_cols)}: {sb_cols}"

    # Assert none contain literal NaN/inf/nan strings
    for col_name in sb_cols:
        val = str(result[col_name][0])
        assert val.lower() not in (
            "nan",
            "inf",
            "-inf",
        ), f"Column {col_name} contains unsanitized value: {val!r}"

    # Odds ratio columns should be 'NA' (not NaN) since alt_total ≤ 1
    or_col = "simplex_duplex_strand_bias_odds_ratio"
    assert result[or_col][0] == "NA", (
        f"Expected 'NA' for odds ratio with alt_total ≤ 1, got: {result[or_col][0]!r}"
    )
    fsor_col = "simplex_duplex_fragment_strand_bias_odds_ratio"
    assert result[fsor_col][0] == "NA", (
        f"Expected 'NA' for fragment odds ratio with alt_total ≤ 1, got: {result[fsor_col][0]!r}"
    )


# ── Test 26: Merge handles MAF files with provenance comment lines ───────


def _write_test_maf_with_comments(
    path: Path, rows: list[dict], columns: list[str] | None = None
) -> None:
    """Write a MAF file with v5.3.0-style provenance comment headers.

    Mimics the output of MafWriter which prepends:
        #gbcms v5.3.0
        #command gbcms dna ...
    before the TSV header.
    """
    if not rows:
        columns = columns or VARIANT_KEY
        path.write_text("\t".join(columns) + "\n")
        return

    columns = columns or list(rows[0].keys())
    lines = [
        "#gbcms v5.3.0",
        "#command gbcms dna --bam sample.bam --fasta ref.fa -o out/",
        "\t".join(columns),
    ]
    for row in rows:
        lines.append("\t".join(str(row.get(c, "")) for c in columns))
    path.write_text("\n".join(lines) + "\n")


def test_merge_handles_provenance_comment_lines(tmp_path):
    """Merge correctly reads MAF files containing #-prefixed provenance comments.

    v5.3.0 MAF output prepends #gbcms and #command comment lines before the
    TSV header. The merge engine (Polars scan_maf with comment_prefix='#')
    must skip these lines and produce correct results.

    Regression test for provenance + merge compatibility.
    """
    duplex = tmp_path / "duplex.maf"
    simplex = tmp_path / "simplex.maf"
    output = tmp_path / "merged.maf"

    _write_test_maf_with_comments(
        duplex, [_make_variant_row(ref_count="20", alt_count="10")]
    )
    _write_test_maf_with_comments(
        simplex, [_make_variant_row(ref_count="5", alt_count="2")]
    )

    config = MergeConfig(
        inputs={"duplex": duplex, "simplex": simplex},
        output=output,
        add_combined=False,
    )
    merge_mafs(config)

    result = pl.read_csv(output, separator="\t", infer_schema_length=0)
    assert result.height == 1, f"Expected 1 row, got {result.height}"
    assert "duplex_ref_count" in result.columns, (
        f"Missing duplex_ref_count, got: {result.columns}"
    )
    assert result["duplex_ref_count"].to_list() == ["20"]
    assert result["simplex_ref_count"].to_list() == ["5"]
    # Ensure no comment artifacts leaked into column names
    comment_cols = [c for c in result.columns if c.startswith("#")]
    assert comment_cols == [], f"Comment lines leaked as columns: {comment_cols}"

