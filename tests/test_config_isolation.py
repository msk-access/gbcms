"""
Tests for DNA vs RNA configuration isolation.

Verifies that GbcmsDnaConfig and GbcmsRnaConfig:
- Have correct mode identifiers
- Have appropriate default values
- Don't leak fields between modes
- Apply correct validation rules
"""

import pytest
from pydantic import ValidationError

from gbcms.models.core import (
    GbcmsBaseConfig,
    GbcmsConfig,
    GbcmsDnaConfig,
    GbcmsRnaConfig,
    OutputConfig,
)

# ── Mode Identity ─────────────────────────────────────────────────────────


def test_dna_config_mode():
    """GbcmsDnaConfig.mode must be 'dna'."""
    assert GbcmsDnaConfig.model_fields["mode"].default == "dna"


def test_rna_config_mode():
    """GbcmsRnaConfig.mode must be 'rna'."""
    assert GbcmsRnaConfig.model_fields["mode"].default == "rna"


def test_deprecated_alias():
    """GbcmsConfig is a deprecated alias for GbcmsDnaConfig."""
    assert GbcmsConfig is GbcmsDnaConfig


# ── RNA-Specific Fields ──────────────────────────────────────────────────


def test_rna_has_strandedness_field():
    """GbcmsRnaConfig has enforce_strandedness (default True)."""
    assert "enforce_strandedness" in GbcmsRnaConfig.model_fields
    assert GbcmsRnaConfig.model_fields["enforce_strandedness"].default is True


def test_rna_has_editing_db_field():
    """GbcmsRnaConfig has rna_editing_db (default None)."""
    assert "rna_editing_db" in GbcmsRnaConfig.model_fields
    assert GbcmsRnaConfig.model_fields["rna_editing_db"].default is None


def test_dna_lacks_rna_fields():
    """GbcmsDnaConfig must NOT have RNA-specific fields."""
    assert "enforce_strandedness" not in GbcmsDnaConfig.model_fields
    assert "rna_editing_db" not in GbcmsDnaConfig.model_fields


# ── Default Values ───────────────────────────────────────────────────────


def test_rna_default_mapq():
    """RNA mode defaults to MAPQ=1 (STAR assigns low MAPQ to novel junctions)."""
    factory = GbcmsRnaConfig.model_fields["quality"].default_factory
    assert factory is not None
    quality = factory()  # type: ignore[call-arg]
    assert quality.min_mapping_quality == 1


def test_dna_default_mapq():
    """DNA mode defaults to MAPQ=20."""
    factory = GbcmsDnaConfig.model_fields["quality"].default_factory
    assert factory is not None
    quality = factory()  # type: ignore[call-arg]
    assert quality.min_mapping_quality == 20


def test_rna_default_backend():
    """RNA mode defaults to PairHMM alignment backend (normalized to 'hmm')."""
    factory = GbcmsRnaConfig.model_fields["alignment"].default_factory
    assert factory is not None
    alignment = factory()  # type: ignore[call-arg]
    # The model stores "pairhmm" but the validator normalizes it to "hmm"
    assert alignment.backend == "hmm"


def test_dna_default_backend():
    """DNA mode now defaults to PairHMM alignment backend."""
    factory = GbcmsDnaConfig.model_fields["alignment"].default_factory
    assert factory is not None
    alignment = factory()  # type: ignore[call-arg]
    # Pydantic default_factory doesn't trigger field_validator, stores 'pairhmm'
    assert alignment.backend == "pairhmm"


# ── Inheritance ──────────────────────────────────────────────────────────


def test_both_inherit_from_base():
    """Both configs must inherit from GbcmsBaseConfig."""
    assert issubclass(GbcmsDnaConfig, GbcmsBaseConfig)
    assert issubclass(GbcmsRnaConfig, GbcmsBaseConfig)


def test_shared_fields_exist():
    """Both configs must have shared fields (apply_baq, umi_tag, threads)."""
    for field in ("apply_baq", "umi_tag", "threads"):
        assert field in GbcmsDnaConfig.model_fields, f"DNA config missing '{field}'"
        assert field in GbcmsRnaConfig.model_fields, f"RNA config missing '{field}'"


# ── Validation ───────────────────────────────────────────────────────────


def test_rna_editing_db_validates_path(tmp_path):
    """rna_editing_db rejects non-existent paths."""
    with pytest.raises(ValidationError, match="RNA editing database not found"):
        GbcmsRnaConfig(
            variant_file=tmp_path / "dummy.vcf",  # won't reach this validator
            bam_files={},
            reference_fasta=tmp_path / "dummy.fa",
            output=OutputConfig(directory=tmp_path),
            rna_editing_db=tmp_path / "nonexistent.vcf",
        )
