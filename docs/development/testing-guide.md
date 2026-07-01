# Testing Guide

This guide covers running tests, adding new tests, and accuracy validation for gbcms.

## Running Tests

### Quick Test

```bash
# Run all tests
pytest -v

# Run with coverage
pytest --cov=gbcms --cov-report=html

# Run specific test file
pytest tests/test_accuracy.py -v
```

### Test Categories

| Category | Files | Purpose |
|:---------|:------|:--------|
| Accuracy | `test_accuracy.py` | SNP, indel, complex variant counting, DP invariant |
| Shifted Indels | `test_shifted_indels.py` | Windowed indel detection (±5bp), 3-layer safeguards |
| Complex Masking | `test_fuzzy_complex.py` | Quality-aware masked comparison, ambiguity detection, MSI gap penalties |
| Fragment Consensus | `test_fragment_consensus.py` | Quality-weighted R1/R2 conflict resolution, DPF invariant |
| INDEL Fragment Consensus | `test_indel_fragment_consensus.py` | Structural ALT priority for INDEL conflicts, Phase 3 dispatch, SNP regression |
| Normalization | `test_normalization.py` | Left-alignment, REF validation, homopolymer detection, dynamic window expansion |
| DP Neither | `test_dp_neither.py` | Gap 1D: DP includes third-allele/neither reads |
| Multi-Allelic | `test_multi_allelic.py` | Gap 1A: Sibling ALT exclusion, overlapping indel DP |
| CLI | `test_cli_sample_id.py` | Argument parsing, input validation, error paths, `--lenient-bam`, subcommands |
| CLI DNA/RNA | `test_cli_dna_rna.py` | Command existence, option isolation, error handling |
| Filters | `test_filters.py` | Read filtering logic |
| MAF | `test_maf_*.py` | MAF column preservation, reader |
| Pipeline | `test_pipeline_v2.py` | End-to-end workflow |
| Strand | `test_strand_counts.py` | Strand-specific counts |
| Alignment Backend | `test_alignment_backend.py` | PairHMM default, backend integration |
| Config Isolation | `test_config_isolation.py` | DNA vs RNA mode, defaults, field validation |
| RNA Output | `test_rna_output.py` | MAF/VCF RNA column presence/absence |
| BAQ | `test_baq.py` | BAQ quality downgrade feature |
| Phase 2 Output | `test_phase2_output.py` | Diagnostic columns (any_alt, partial_alt, n_count) MAF/VCF |
| Diagnostic Flags | `test_diagnostic_flags.py` | gbcms_diagnostic flag computation (ZERO_ALT, PARTIAL_DOMINANT, etc.) |
| MNP Rescue | `test_rescue_mnp.py` | --rescue-mnp flag, conditional columns, candidate filtering, audit trail |

### Rust-Level Tests

```bash
# Run Rust unit tests (normalize + counting inline tests)
cd rust && cargo test

# Run a specific Rust test
cargo test test_window_expansion_long_homopolymer
```

Rust tests live inside `#[cfg(test)]` modules and cover:

| Area | Tests | Purpose |
|:-----|:------|:--------|
| Left-alignment | 10+ | SNP passthrough, homopolymer shifts, offset handling |
| Repeat detection | 3 | `find_tandem_repeat()` edge cases |
| Adaptive padding | 3 | Context padding from repeat spans |
| Window expansion | 1 | Gap 1B: >100bp repeat normalization |
| MNP classification | 17 | Masked quality, partial counting, N-base handling |
| N-base coverage | 8 | SNP/MNP/PairHMM N guards, n_count propagation |
| Fragment consensus | 11 | `resolve()` structural ALT priority, `observe()` sticky flag, quality paths |
| Invariants | 1 | `any_alt = ad + partial_alt`, depth decomposition |

---

## Test Structure

```
tests/
├── conftest.py                  # Shared pytest fixtures (paths, RNA BAM, editing DB)
├── helpers.py                   # Shared helpers (build_bam, make_read, count_one, count_both)
├── test_accuracy.py             # Variant type accuracy + DP invariant
├── test_alignment_backend.py    # PairHMM default, backend integration
├── test_baq.py                  # BAQ quality downgrade
├── test_cli_dna_rna.py          # dna/rna command existence, option isolation
├── test_cli_sample_id.py        # CLI argument parsing
├── test_config_isolation.py      # DNA vs RNA config isolation + validation
├── test_dp_neither.py           # Gap 1D: DP includes third-allele reads
├── test_filters.py              # Read filtering
├── test_fragment_consensus.py   # Fragment-level quality consensus + DPF invariant
├── test_indel_fragment_consensus.py  # INDEL structural ALT priority + Phase 3 dispatch
├── test_fuzzy_complex.py        # Quality-aware masked complex matching + MSI penalties
├── test_maf_preservation.py     # MAF column preservation
├── test_maf_reader.py           # MAF input parsing
├── test_mfsd_flag.py            # mFSD flag behavior
├── test_multi_allelic.py        # Gap 1A: Sibling ALT exclusion
├── test_normalization.py        # Left-alignment, REF validation, window expansion
├── test_pipeline_v2.py          # End-to-end pipeline
├── test_phase2_output.py        # Diagnostic columns (any_alt, partial_alt, n_count)
├── test_diagnostic_flags.py     # v4.2.0 gbcms_diagnostic flag computation
├── test_rescue_mnp.py           # v4.3.0 --rescue-mnp MNP decomposition rescue
├── test_rna_output.py           # RNA MAF/VCF column presence/absence
├── test_shifted_indels.py       # Windowed indel detection (±5bp)
└── test_strand_counts.py        # Strand-specific counts
```

---

## Writing Tests

### Basic Test Template

```python
import pytest
from pathlib import Path

def test_my_feature(tmp_path):
    """Test description."""
    # Arrange
    input_file = tmp_path / "input.txt"
    input_file.write_text("test data")
    
    # Act
    result = my_function(input_file)
    
    # Assert
    assert result.success
    assert result.count == 42
```

### Accuracy Test Template

```python
def test_snp_accuracy():
    """Verify SNP counting against known BAM."""
    # Create variant
    variant = Variant("chr1", 100, "A", "T", "SNP")
    
    # Run counting
    results = count_bam(bam_path, [variant], decomposed=[None], ...)
    
    # Validate allele counts
    assert results[0].rd == 50
    assert results[0].ad == 10
    # Gap 1D invariant: DP includes ALL reads (including 'neither')
    assert results[0].dp >= results[0].rd + results[0].ad
```

### Multi-Allelic Test Template

```python
def test_with_siblings():
    """Verify sibling ALT exclusion at multi-allelic sites."""
    v1 = Variant("chr1", 100, "A", "T", "SNP")
    v2 = Variant("chr1", 100, "A", "C", "SNP")
    
    results = count_bam(
        bam_path, [v1, v2], decomposed=[None, None],
        sibling_variants=[[v2], [v1]],  # Gap 1A: sibling info
        ...
    )
```

### Key Invariants to Assert

All counting tests should verify:

| Invariant | Description |
|:----------|:------------|
| `dp >= rd + ad` | DP includes 'neither' reads (Gap 1D) |
| `dpf >= rdf + adf` | DPF includes discarded ambiguous fragments |
| `rd == rd_fwd + rd_rev` | Strand consistency |
| `ad == ad_fwd + ad_rev` | Strand consistency |
| `any_alt = ad + partial_alt` | Phase 2 decomposition invariant |
| `any_alt >= ad` | Partial count is non-negative |
| `dp >= rd + ad + partial_alt + n_count` | Depth decomposition with N-base diagnostic |

### No Silent Failures Matrix

Every N/masked/partial path must produce a deterministic, traceable outcome:

| Scenario | Classification | n_count | partial_alt | Log Level |
|:---------|:---------------|:--------|:------------|:----------|
| SNP with N base | `neither_n` (uninformative) | +1 | — | `trace` |
| MNP with N at 1 position | ALT/REF via unmasked | +1 | depends on match | `trace` |
| MNP with N at ALL positions | `LowQuality` (DP only) | +1 | 0 | `trace` |
| MNP with N + low-BQ | `LowQuality` (DP only) | +1 | 0 | `trace` |
| Complex with N in haplotype | via masked compare | +1 | depends on match | `trace` |
| ALT = "N" in input VCF/MAF | `FAIL` + reason `ALT_CONTAINS_N` | — | — | `warn` (validation) |
| ThirdAllele with partial match | neither + partial | — | +1 | `trace` |

---

## Manual Validation

### Using samtools for Spot-Check

```bash
# Check counts at specific position
samtools mpileup -r chr1:100-100 -q 20 \
    -f ref.fa sample.bam 2>/dev/null | \
    awk '{print "DP="$4}'
```

### Comparing with gbcms Output

```bash
# Run gbcms
gbcms dna -v variants.maf -b sample.bam -f ref.fa -o output/

# Check output
awk -F'\t' 'NR==2 {print "REF="$41, "ALT="$42}' output/*.maf
```

---

## Accuracy Validation

### Variant Types Tested

| Type | Test | Status |
|:-----|:-----|:------:|
| SNP | `test_snp_accuracy` | ✅ |
| Insertion | `test_insertion_accuracy` | ✅ |
| Deletion | `test_deletion_accuracy` | ✅ |
| Complex | `test_complex_accuracy` | ✅ |
| MNP | `test_mnp_accuracy` | ✅ |
| Shifted Indels | `test_shifted_indels.py` (15 cases) | ✅ |
| Complex Masking | `test_fuzzy_complex.py` (15 cases) | ✅ |
| DP Neither | `test_dp_neither.py` (3 cases) | ✅ |
| Multi-Allelic | `test_multi_allelic.py` (4 cases) | ✅ |
| Fragment Consensus | `test_fragment_consensus.py` (3 cases) | ✅ |
| INDEL Fragment Consensus | `test_indel_fragment_consensus.py` (11 cases) | ✅ |
| Window Expansion | `test_normalization.py` (9 cases) | ✅ |
| MSI Gap Penalties | `test_fuzzy_complex.py::TestGap3A` | ✅ |

### Real-World Validation

```bash
# Compare gbcms vs samtools for a SNP
# Position: chr1:11168293 G>A

# gbcms output
awk -F'\t' '$5=="1" && $6=="11168293" {print "REF="$41, "ALT="$42}' output.maf

# samtools output
samtools mpileup -r 1:11168293-11168293 -q 20 -f ref.fa sample.bam | \
    awk '{gsub(/\^.|\$/,"",$5); print "DP="$4, "Pileup="$5}'
```

---

## Coverage Targets

| Module | Target | Current |
|:-------|:------:|:-------:|
| cli.py | 90% | ~76% |
| pipeline.py | 70% | 16% |
| io/input.py | 85% | 80% |
| io/output.py | 90% | 84% |
| models/core.py | 90% | 96% |

**Current totals**: 266 Python + 161 Rust tests.

Run coverage report:

```bash
pytest --cov=gbcms --cov-report=html
open htmlcov/index.html
```

---

## Related

- [Developer Guide](developer-guide.md) — Setup, build commands, and project layout
- [Contributing](contributing.md) — Code standards and pull request process
- [Architecture](../reference/architecture.md) — Module structure and data flow
- [Allele Classification](../reference/allele-classification.md) — Engine logic being tested
