# Integration Test Data

This directory contains real genomic data files for integration testing of gbcms.

## Files

- **integration_test_reference.fa**: 20kb region from chr1:11,180,000-11,200,000
- **integration_test_variants.vcf**: 155 real variants from the same region
- **sample1_integration_test.bam**: 41,556 alignments from Sample1 in the region (0.96 MB)
- **sample2_integration_test.bam**: 2,762 alignments from Sample2 in the region (0.11 MB)
- **masked_repeat_reads.json**: repeat shapes for `test_masked_repeat_reads_contract.py`.
  Each is a reference and two haplotypes (synthetic VNTRs and tandem duplications, and two
  GRCh38 reference segments), with paired reads simulated from each haplotype and placed
  by BWA-MEM 0.7.17 as production runs it (`-K 1000000 -Y`). A read is stored as where BWA
  wrote it (flag, position, CIGAR, MAPQ, mate) and its haplotype offset; its bases are the
  haplotype's from there. No sample data.

## Usage

These files are used by the pytest fixtures in `conftest.py` for integration testing:

```python
@pytest.fixture
def sample_fasta(temp_dir: Path) -> Path:
    test_data_dir = Path(__file__).parent / "testdata"
    return test_data_dir / "integration_test_reference.fa"

@pytest.fixture  
def sample_bam(temp_dir: Path, sample_fasta: Path) -> Path:
    test_data_dir = Path(__file__).parent / "testdata"
    return test_data_dir / "sample1_integration_test.bam"

@pytest.fixture
def multi_sample_bams(temp_dir: Path) -> dict:
    test_data_dir = Path(__file__).parent / "testdata"
    return {
        "sample1": test_data_dir / "sample1_integration_test.bam",
        "sample2": test_data_dir / "sample2_integration_test.bam",
    }
```

## Testing Commands

```bash
# Single sample test
gbcms count run --fasta tests/testdata/integration_test_reference.fa \
                --bam sample1:tests/testdata/sample1_integration_test.bam \
                --vcf tests/testdata/integration_test_variants.vcf \
                --output /tmp/test_output.vcf

# Multi-sample test  
gbcms count run --fasta tests/testdata/integration_test_reference.fa \
                --bam sample1:tests/testdata/sample1_integration_test.bam \
                --bam sample2:tests/testdata/sample2_integration_test.bam \
                --vcf tests/testdata/integration_test_variants.vcf \
                --output /tmp/test_multi_output.vcf

# Fillout format
gbcms count run --fasta tests/testdata/integration_test_reference.fa \
                --bam sample1:tests/testdata/sample1_integration_test.bam \
                --vcf tests/testdata/integration_test_variants.vcf \
                --fillout --output /tmp/test_fillout.maf
```

## Data Source

These files were extracted from real genomic data:
- **FASTA**: Human genome reference (hg19) chr1 region
- **BAMs**: Real sequencing data from tumor/normal samples  
- **VCF**: Real variant calls from the same samples
- **Region**: chr1:11,180,000-11,200,000 (contains multiple variants for testing)

