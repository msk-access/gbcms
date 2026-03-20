"""
Test strand count output in VCF and MAF formats.
"""

import csv

import pytest

from gbcms.io.output import MafWriter, VcfWriter
from gbcms.models.core import Variant, VariantType


# Mock objects for testing
class MockCounts:
    def __init__(self):
        # Basic counts
        self.dp = 20
        self.rd = 15
        self.ad = 5

        # Strand counts
        self.dp_fwd = 12
        self.dp_rev = 8
        self.rd_fwd = 10
        self.rd_rev = 5
        self.ad_fwd = 2
        self.ad_rev = 3

        # Fragment counts
        self.dpf = 10
        self.rdf = 7
        self.adf = 3

        # Fragment strand counts
        self.rdf_fwd = 4
        self.rdf_rev = 3
        self.adf_fwd = 1
        self.adf_rev = 2

        # Stats
        self.sb_pval = 0.05
        self.sb_or = 1.5
        self.fsb_pval = 0.1
        self.fsb_or = 1.2

        # mFSD fields (default to zero counts, NaN stats)
        _nan = float("nan")
        self.mfsd_ref_count = 0
        self.mfsd_alt_count = 0
        self.mfsd_nonref_count = 0
        self.mfsd_n_count = 0
        self.mfsd_ref_mean = _nan
        self.mfsd_alt_mean = _nan
        self.mfsd_nonref_mean = _nan
        self.mfsd_n_mean = _nan
        self.mfsd_alt_llr = _nan
        self.mfsd_ref_llr = _nan
        self.mfsd_delta_alt_ref = _nan
        self.mfsd_ks_alt_ref = _nan
        self.mfsd_pval_alt_ref = _nan
        self.mfsd_delta_alt_nonref = _nan
        self.mfsd_ks_alt_nonref = _nan
        self.mfsd_pval_alt_nonref = _nan
        self.mfsd_delta_ref_nonref = _nan
        self.mfsd_ks_ref_nonref = _nan
        self.mfsd_pval_ref_nonref = _nan
        self.mfsd_delta_alt_n = _nan
        self.mfsd_ks_alt_n = _nan
        self.mfsd_pval_alt_n = _nan
        self.mfsd_delta_ref_n = _nan
        self.mfsd_ks_ref_n = _nan
        self.mfsd_pval_ref_n = _nan
        self.mfsd_delta_nonref_n = _nan
        self.mfsd_ks_nonref_n = _nan
        self.mfsd_pval_nonref_n = _nan

        # ── RNA-specific fields ───────────────────────────────────────────────
        # Required by both VcfWriter(mode="rna") and MafWriter(mode="rna").
        # Use realistic non-zero values so formatting code paths are exercised.
        self.sense_depth = 12
        self.antisense_depth = 3
        self.sense_strand_alt_count = 4
        self.rna_editing_site_overlap = False
        self.splice_spanning_count = 2


@pytest.fixture
def mock_variant():
    return Variant(
        chrom="chr1", pos=100, ref="A", alt="T", variant_type=VariantType.SNP, original_id="rs123"
    )


@pytest.fixture
def mock_counts():
    return MockCounts()


def test_vcf_output(tmp_path, mock_variant, mock_counts):
    output_path = tmp_path / "test.vcf"
    writer = VcfWriter(output_path)
    writer.write(mock_variant, mock_counts)
    writer.close()

    assert output_path.exists()

    with open(output_path) as f:
        lines = f.readlines()

    # Check header
    header_lines = [line for line in lines if line.startswith("##")]
    format_lines = [line for line in header_lines if line.startswith("##FORMAT")]

    # Verify Number=2 for count fields
    assert any("ID=DP,Number=2" in line for line in format_lines)
    assert any("ID=RD,Number=2" in line for line in format_lines)
    assert any("ID=AD,Number=2" in line for line in format_lines)
    assert any("ID=RDF,Number=2" in line for line in format_lines)
    assert any("ID=ADF,Number=2" in line for line in format_lines)

    # Verify Number=1 for VAF/FAF
    assert any("ID=VAF,Number=1" in line for line in format_lines)
    assert any("ID=FAF,Number=1" in line for line in format_lines)

    # Check data line
    data_line = [line for line in lines if not line.startswith("#")][0]
    fields = data_line.strip().split("\t")

    # FORMAT is column 8 (0-based), Sample is column 9
    format_str = fields[8]
    sample_data = fields[9]

    format_keys = format_str.split(":")
    sample_values = sample_data.split(":")

    data_dict = dict(zip(format_keys, sample_values, strict=True))

    # Verify values
    assert data_dict["DP"] == "15,5"  # ref_total, alt_total
    assert data_dict["RD"] == "10,5"  # ref_fwd, ref_rev
    assert data_dict["AD"] == "2,3"  # alt_fwd, alt_rev
    assert data_dict["RDF"] == "4,3"  # ref_frag_fwd, ref_frag_rev
    assert data_dict["ADF"] == "1,2"  # alt_frag_fwd, alt_frag_rev

    # Verify VAF calculations
    # VAF = AD / (RD + AD) = 5 / 20 = 0.25
    assert float(data_dict["VAF"]) == 0.2500

    # FAF = ADF / (RDF + ADF) = 3 / 10 = 0.3
    assert float(data_dict["FAF"]) == 0.3000


def test_maf_output(tmp_path, mock_variant, mock_counts):
    """Test VCF→MAF output with strand counts and GDC coordinate conversion."""
    output_path = tmp_path / "test.maf"
    writer = MafWriter(output_path)
    writer.write(mock_variant, mock_counts)
    writer.close()

    assert output_path.exists()

    with open(output_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        row = next(reader)

    # VCF→MAF coordinate conversion (SNP at 0-based pos=100 → 1-based pos=101)
    assert row["Start_Position"] == "101"
    assert row["End_Position"] == "101"
    assert row["Variant_Type"] == "SNP"
    assert row["Reference_Allele"] == "A"
    assert row["Tumor_Seq_Allele2"] == "T"

    # Check standard count fields (default: no prefix)
    assert row["ref_count"] == "15"
    assert row["alt_count"] == "5"
    assert row["total_count"] == "20"

    # Check VAF
    assert float(row["vaf"]) == 0.2500
    assert float(row["vaf_fragment"]) == 0.3000

    # Check strand count fields (default: no prefix)
    assert row["ref_count_forward"] == "10"
    assert row["ref_count_reverse"] == "5"
    assert row["alt_count_forward"] == "2"
    assert row["alt_count_reverse"] == "3"

    assert row["ref_count_fragment_forward"] == "4"
    assert row["ref_count_fragment_reverse"] == "3"
    assert row["alt_count_fragment_forward"] == "1"
    assert row["alt_count_fragment_reverse"] == "2"
