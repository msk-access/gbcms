"""Tests for mFSD interactive HTML report generation.

Uses synthetic test fixtures with generic sample identifiers (no patient data).
"""

from pathlib import Path

import pytest

from gbcms.report import generate_mfsd_report

# ── Fixture paths ──────────────────────────────────────────────────────────────
TESTDATA = Path(__file__).parent / "testdata"
MULTI_PARQUET = TESTDATA / "mfsd_multi.parquet"
MULTI_MAF = TESTDATA / "mfsd_multi.maf"
SINGLE_PARQUET = TESTDATA / "mfsd_single.parquet"
SINGLE_MAF = TESTDATA / "mfsd_single.maf"


@pytest.fixture
def multi_report(tmp_path: Path) -> Path:
    """Generate a multi-variant (3 variant) mFSD report."""
    out = tmp_path / "multi.mfsd_report.html"
    return generate_mfsd_report(
        parquet_path=MULTI_PARQUET,
        maf_path=MULTI_MAF,
        output_path=out,
        min_alt=3,
        max_variants=20,
        sample_name="TEST-SAMPLE",
    )


@pytest.fixture
def single_report(tmp_path: Path) -> Path:
    """Generate a single-variant mFSD report."""
    out = tmp_path / "single.mfsd_report.html"
    return generate_mfsd_report(
        parquet_path=SINGLE_PARQUET,
        maf_path=SINGLE_MAF,
        output_path=out,
        min_alt=3,
        max_variants=20,
        sample_name="TEST-SAMPLE",
    )


# ── Tests ──────────────────────────────────────────────────────────────────────


class TestReportCreation:
    """Verify report files are created successfully."""

    def test_output_file_created(self, multi_report: Path) -> None:
        """Report should exist on disk as a non-empty HTML file."""
        assert multi_report.exists()
        assert multi_report.stat().st_size > 0
        assert multi_report.suffix == ".html"

    def test_single_variant_creates_file(self, single_report: Path) -> None:
        """Single-variant report should also be created."""
        assert single_report.exists()
        assert single_report.stat().st_size > 0


class TestNavigator:
    """Verify the variant navigator UI."""

    def test_multi_variant_has_navigator(self, multi_report: Path) -> None:
        """Multi-variant reports should include the sticky navigator bar."""
        html = multi_report.read_text()
        assert 'id="variant-nav"' in html
        assert 'id="nav-dropdown"' in html
        assert 'id="nav-prev"' in html
        assert 'id="nav-next"' in html
        assert 'id="nav-counter"' in html
        assert 'id="nav-mode-toggle"' in html

    def test_single_variant_no_navigator(self, single_report: Path) -> None:
        """Single-variant reports should NOT include the navigator bar."""
        html = single_report.read_text()
        assert 'id="variant-nav"' not in html


class TestPlotly:
    """Verify Plotly chart integration."""

    def test_report_contains_plotly(self, multi_report: Path) -> None:
        """Report should include the Plotly.js CDN script tag."""
        html = multi_report.read_text()
        assert "plotly" in html.lower()
        assert "<script" in html

    def test_variant_card_data_attributes(self, multi_report: Path) -> None:
        """Each variant card should have data-variant-index attributes."""
        html = multi_report.read_text()
        assert 'data-variant-index="0"' in html
        assert "data-variant-label=" in html


class TestThemeAndBranding:
    """Verify theme toggle and branding elements."""

    def test_report_theme_toggle(self, multi_report: Path) -> None:
        """Report should include a theme toggle button."""
        html = multi_report.read_text()
        assert "toggleTheme()" in html
        assert "data-theme=" in html

    def test_report_branded_footer(self, multi_report: Path) -> None:
        """Report should include the branded footer."""
        html = multi_report.read_text()
        assert "Ronak Shah" in html or "MSK" in html


class TestSummaryCards:
    """Verify the summary dashboard section."""

    def test_report_summary_cards(self, multi_report: Path) -> None:
        """Report should contain Fragment Origin Signal classification cards."""
        html = multi_report.read_text()
        assert "TUMOR-LIKE" in html
        assert "Fragment Origin Signal" in html or "Per-Variant Analysis" in html


class TestFiltering:
    """Verify min_alt filtering logic."""

    def test_min_alt_filtering(self, tmp_path: Path) -> None:
        """High min_alt threshold should exclude variants with few ALT fragments."""
        out = tmp_path / "filtered.html"
        report = generate_mfsd_report(
            parquet_path=MULTI_PARQUET,
            maf_path=MULTI_MAF,
            output_path=out,
            min_alt=500,  # Only APC/EGFR have enough ALT fragments
            max_variants=20,
            sample_name="TEST-SAMPLE",
        )
        html = report.read_text()
        # TP53 (1 ALT fragment) should be excluded
        # The number of variant cards should be less than 3
        card_count = html.count('class="card"')
        assert card_count < 3

    def test_max_variants_limit(self, tmp_path: Path) -> None:
        """max_variants should cap the number of variant cards."""
        out = tmp_path / "limited.html"
        report = generate_mfsd_report(
            parquet_path=MULTI_PARQUET,
            maf_path=MULTI_MAF,
            output_path=out,
            min_alt=3,
            max_variants=1,
            sample_name="TEST-SAMPLE",
        )
        html = report.read_text()
        card_count = html.count('class="card"')
        assert card_count == 1


class TestErrorHandling:
    """Verify graceful error handling."""

    def test_missing_parquet_raises(self, tmp_path: Path) -> None:
        """Should raise FileNotFoundError for missing parquet."""
        with pytest.raises(FileNotFoundError):
            generate_mfsd_report(
                parquet_path=tmp_path / "nonexistent.parquet",
                maf_path=MULTI_MAF,
                output_path=tmp_path / "out.html",
            )

    def test_missing_maf_raises(self, tmp_path: Path) -> None:
        """Should raise FileNotFoundError for missing MAF."""
        with pytest.raises(FileNotFoundError):
            generate_mfsd_report(
                parquet_path=MULTI_PARQUET,
                maf_path=tmp_path / "nonexistent.maf",
                output_path=tmp_path / "out.html",
            )
