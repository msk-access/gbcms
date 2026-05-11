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


def test_rna_default_baq():
    """RNA mode defaults to apply_baq=True (no upstream BQSR/consensus)."""
    assert GbcmsRnaConfig.model_fields["apply_baq"].default is True


def test_dna_default_baq():
    """DNA mode defaults to apply_baq=False (upstream BQSR/consensus)."""
    assert GbcmsDnaConfig.model_fields["apply_baq"].default is False


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


# ── v5.0.0: Library Type + GTF Fields ───────────────────────────────────


def test_rna_has_library_type_field():
    """GbcmsRnaConfig has library_type with default 'capture'."""
    assert "library_type" in GbcmsRnaConfig.model_fields
    assert GbcmsRnaConfig.model_fields["library_type"].default == "capture"


def test_rna_library_type_validator_rejects_invalid():
    """library_type rejects unsupported values at field validation."""
    with pytest.raises(ValueError, match="Invalid library_type"):
        GbcmsRnaConfig.validate_library_type("pcr")


def test_rna_library_type_validator_normalizes_case():
    """library_type normalizes to lowercase."""
    assert GbcmsRnaConfig.validate_library_type("Amplicon") == "amplicon"
    assert GbcmsRnaConfig.validate_library_type("CAPTURE") == "capture"


def test_dna_lacks_library_type_field():
    """GbcmsDnaConfig must NOT have library_type."""
    assert "library_type" not in GbcmsDnaConfig.model_fields


def test_rna_has_gtf_field():
    """GbcmsRnaConfig has gtf field with default None."""
    assert "gtf" in GbcmsRnaConfig.model_fields
    assert GbcmsRnaConfig.model_fields["gtf"].default is None


def test_dna_lacks_gtf_field():
    """GbcmsDnaConfig must NOT have gtf."""
    assert "gtf" not in GbcmsDnaConfig.model_fields


def test_amplicon_auto_disables_strandedness(tmp_path):
    """model_validator auto-disables enforce_strandedness when library_type='amplicon'."""
    # Create real dummy files so Pydantic path validators pass
    vcf = tmp_path / "dummy.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")
    fasta = tmp_path / "dummy.fa"
    fasta.write_text(">chr1\nACGT\n")
    fai = tmp_path / "dummy.fa.fai"
    fai.write_text("chr1\t4\t6\t4\t5\n")

    config = GbcmsRnaConfig(
        variant_file=vcf,
        bam_files={},
        reference_fasta=fasta,
        output=OutputConfig(directory=tmp_path),
        library_type="amplicon",
        enforce_strandedness=True,  # should be auto-overridden
    )
    assert config.library_type == "amplicon"
    assert config.enforce_strandedness is False, (
        "amplicon mode should auto-disable enforce_strandedness via model_validator"
    )


# ── v5.0.0: rescue_mnp_threshold Config Validation ─────────────────────


def test_rescue_mnp_threshold_default():
    """rescue_mnp_threshold defaults to 1.0 (permissive, C++ compatible)."""
    assert GbcmsBaseConfig.model_fields["rescue_mnp_threshold"].default == 1.0


def test_rescue_mnp_threshold_shared():
    """Both DNA and RNA configs inherit rescue_mnp_threshold from base."""
    assert "rescue_mnp_threshold" in GbcmsDnaConfig.model_fields
    assert "rescue_mnp_threshold" in GbcmsRnaConfig.model_fields


def test_rescue_mnp_threshold_rejects_above_1():
    """rescue_mnp_threshold > 1.0 is rejected by Pydantic ge/le constraint."""
    with pytest.raises(ValidationError, match="rescue_mnp_threshold"):
        GbcmsBaseConfig.model_validate(
            {"rescue_mnp_threshold": 1.5},
            context={"_skip_file_validation": True},
        )


def test_rescue_mnp_threshold_rejects_below_0():
    """rescue_mnp_threshold < 0.0 is rejected by Pydantic ge/le constraint."""
    with pytest.raises(ValidationError, match="rescue_mnp_threshold"):
        GbcmsBaseConfig.model_validate(
            {"rescue_mnp_threshold": -0.1},
            context={"_skip_file_validation": True},
        )

