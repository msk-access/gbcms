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
