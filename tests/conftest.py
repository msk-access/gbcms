"""
Shared test fixtures for gbcms test suite.

Provides path fixtures for testdata/ directory assets and ensures
tests/ is on sys.path so helpers.py can be imported.
"""

import sys
from pathlib import Path

import pytest

# Make tests/ importable so helpers.py can be imported from test modules
sys.path.insert(0, str(Path(__file__).parent))

from helpers import SUITE_READ_JUDGE, UNJUDGED, install_suite_read_judge  # noqa: E402

# GBCMS_TEST_READ_JUDGE=haplotype: the whole suite under the haplotype judge (helpers.py).
if SUITE_READ_JUDGE == "haplotype":
    install_suite_read_judge()
elif SUITE_READ_JUDGE != "rules":
    raise pytest.UsageError(
        f"GBCMS_TEST_READ_JUDGE={SUITE_READ_JUDGE!r}: expected 'rules' or 'haplotype'"
    )


def pytest_report_header(config):
    return f"read judge: {SUITE_READ_JUDGE} (GBCMS_TEST_READ_JUDGE)"


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    """Under the haplotype judge, the engine calls it could not take (no FASTA beside the
    BAM, of its contigs, to read haplotypes from): those tests ran under the rules, not the
    judge."""
    if SUITE_READ_JUDGE != "haplotype" or not UNJUDGED:
        return
    tests = sorted({t for t, _ in UNJUDGED})
    terminalreporter.section("not run under the haplotype judge")
    terminalreporter.write_line(
        f"{len(tests)} tests made {len(UNJUDGED)} engine calls under the rules "
        "(no single FASTA of the BAM's contigs beside it):"
    )
    for t in tests:
        terminalreporter.write_line(f"  {t}")


# ── Path Fixtures ─────────────────────────────────────────────────────────


@pytest.fixture
def testdata_dir() -> Path:
    """Path to the testdata/ directory containing real BAM/VCF/FASTA files."""
    return Path(__file__).parent / "testdata"


@pytest.fixture
def sample_bam(testdata_dir: Path) -> str:
    """Path to sample1_integration_test.bam."""
    return str(testdata_dir / "sample1_integration_test.bam")


@pytest.fixture
def sample_fasta(testdata_dir: Path) -> str:
    """Path to integration_test_reference.fa."""
    return str(testdata_dir / "integration_test_reference.fa")


@pytest.fixture
def sample_vcf(testdata_dir: Path) -> str:
    """Path to integration_test_variants.vcf."""
    return str(testdata_dir / "integration_test_variants.vcf")
