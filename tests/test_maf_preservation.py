"""
Test preservation of MAF input columns in output.

Verifies that MafWriter in MAF→MAF mode:
1. Preserves ALL original input columns (positional, allele, type, extras).
2. Appends gbcms count columns with configurable prefix.
3. Never overwrites original values.
"""

import csv

from gbcms.io.input import MafReader
from gbcms.io.output import MafWriter


# Mock counts object matching Rust BaseCounts interface
class MockCounts:
    def __init__(self):
        self.dp = 100
        self.rd = 80
        self.ad = 20
        self.dp_fwd = 50
        self.dp_rev = 50
        self.rd_fwd = 40
        self.rd_rev = 40
        self.ad_fwd = 10
        self.ad_rev = 10
        self.dpf = 50
        self.rdf = 40
        self.adf = 10
        self.rdf_fwd = 20
        self.rdf_rev = 20
        self.adf_fwd = 5
        self.adf_rev = 5
        self.sb_pval = 1.0
        self.sb_or = 1.0
        self.fsb_pval = 1.0
        self.fsb_or = 1.0
        # Decomposed ALT counting
        self.any_alt = 20  # same as ad for SNP (no partial)
        self.partial_alt = 0
        self.n_count = 0

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
        # Needed so this MockCounts works with MafWriter(mode='rna') without AttributeError.
        self.sense_depth = 0
        self.antisense_depth = 0
        self.sense_strand_alt_count = 0
        self.rna_editing_site_overlap = False
        self.splice_spanning_count = 0


def test_maf_column_preservation(tmp_path):
    """Original MAF columns must be preserved in output, including extras."""
    # 1. Create a dummy MAF input with standard + extra columns
    input_maf = tmp_path / "input.maf"
    with open(input_maf, "w") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
                "Hugo_Symbol",
                "Entrez_Gene_Id",
                "Variant_Type",
                "Extra_Column",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "Chromosome": "chr1",
                "Start_Position": "100",
                "End_Position": "100",
                "Reference_Allele": "A",
                "Tumor_Seq_Allele2": "T",
                "Hugo_Symbol": "TP53",
                "Entrez_Gene_Id": "7157",
                "Variant_Type": "SNP",
                "Extra_Column": "custom_value",
            }
        )

    # 2. Read the MAF
    reader = MafReader(input_maf)
    variants = list(reader)
    assert len(variants) == 1
    variant = variants[0]

    # Verify metadata was captured
    assert variant.metadata["Hugo_Symbol"] == "TP53"
    assert variant.metadata["Entrez_Gene_Id"] == "7157"
    assert variant.metadata["Extra_Column"] == "custom_value"

    # 3. Write the MAF (default: no prefix)
    output_maf = tmp_path / "output.maf"
    maf_writer = MafWriter(output_maf)
    maf_writer.write(variant, MockCounts())
    maf_writer.close()

    # 4. Verify output
    with open(output_maf) as f:
        csv_reader = csv.DictReader(f, delimiter="\t")
        row = next(csv_reader)

        # Original columns preserved exactly
        assert row["Hugo_Symbol"] == "TP53"
        assert row["Entrez_Gene_Id"] == "7157"
        assert row["Start_Position"] == "100", "Original Start_Position must be preserved"
        assert row["End_Position"] == "100", "Original End_Position must be preserved"
        assert row["Reference_Allele"] == "A", "Original Reference_Allele must be preserved"
        assert row["Tumor_Seq_Allele2"] == "T", "Original Tumor_Seq_Allele2 must be preserved"
        assert row["Variant_Type"] == "SNP", "Original Variant_Type must be preserved"

        # Extra input columns must be preserved
        assert row["Extra_Column"] == "custom_value", "Extra input columns must survive"

        # gbcms count columns appended (no prefix by default)
        assert row["total_count"] == "100"
        assert row["ref_count"] == "80"
        assert row["alt_count"] == "20"


def test_maf_column_prefix(tmp_path):
    """Verify --column-prefix controls count column names."""
    input_maf = tmp_path / "input.maf"
    with open(input_maf, "w") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "Chromosome": "chr1",
                "Start_Position": "100",
                "End_Position": "100",
                "Reference_Allele": "A",
                "Tumor_Seq_Allele2": "T",
            }
        )

    reader = MafReader(input_maf)
    variants = list(reader)
    variant = variants[0]

    # Test with 't_' prefix (legacy compat)
    output_maf = tmp_path / "output_legacy.maf"
    maf_writer = MafWriter(output_maf, column_prefix="t_")
    maf_writer.write(variant, MockCounts())
    maf_writer.close()

    with open(output_maf) as f:
        csv_reader = csv.DictReader(f, delimiter="\t")
        row = next(csv_reader)
        assert row["t_ref_count"] == "80"
        assert row["t_alt_count"] == "20"
        assert row["t_total_count"] == "100"
        assert "ref_count" not in row, "Unprefixed columns should not appear"


def test_maf_preserve_barcode(tmp_path):
    """Verify --preserve-barcode controls Tumor_Sample_Barcode handling."""
    input_maf = tmp_path / "input.maf"
    with open(input_maf, "w") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "Chromosome",
                "Start_Position",
                "End_Position",
                "Reference_Allele",
                "Tumor_Seq_Allele2",
                "Tumor_Sample_Barcode",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "Chromosome": "chr1",
                "Start_Position": "100",
                "End_Position": "100",
                "Reference_Allele": "A",
                "Tumor_Seq_Allele2": "T",
                "Tumor_Sample_Barcode": "C-XULUC7-L001-d01",
            }
        )

    reader = MafReader(input_maf)
    variant = list(reader)[0]

    # Default: BAM sample name overrides original barcode
    out_default = tmp_path / "out_default.maf"
    w1 = MafWriter(out_default)
    w1.write(variant, MockCounts(), sample_name="Plasma_CtDNA")
    w1.close()

    with open(out_default) as f:
        row = next(csv.DictReader(f, delimiter="\t"))
        assert (
            row["Tumor_Sample_Barcode"] == "Plasma_CtDNA"
        ), "Default should override barcode with BAM name"

    # --preserve-barcode: original barcode preserved
    out_preserve = tmp_path / "out_preserve.maf"
    w2 = MafWriter(out_preserve, preserve_barcode=True)
    w2.write(variant, MockCounts(), sample_name="Plasma_CtDNA")
    w2.close()

    with open(out_preserve) as f:
        row = next(csv.DictReader(f, delimiter="\t"))
        assert (
            row["Tumor_Sample_Barcode"] == "C-XULUC7-L001-d01"
        ), "preserve_barcode=True should keep original barcode"


def test_vcf_to_maf_always_uses_sample_name(tmp_path):
    """VCF→MAF path always uses BAM sample name regardless of preserve_barcode.

    VCF records have no Tumor_Sample_Barcode field, so MafWriter must always
    populate it from the sample_name argument. The preserve_barcode flag only
    applies to MAF→MAF (variant.metadata is populated from MAF input rows).
    """
    from gbcms.models.core import Variant, VariantType

    # VCF-origin variant: metadata is empty (no MAF row dict)
    vcf_variant = Variant(
        chrom="chr1",
        pos=99,
        ref="C",
        alt="G",
        variant_type=VariantType.SNP,
        original_id="rs42",
    )
    assert not vcf_variant.metadata, "VCF-origin variant must have no metadata"

    out_vcf_maf = tmp_path / "vcf_origin.maf"
    # preserve_barcode=True should NOT affect VCF→MAF — no metadata to preserve
    w = MafWriter(out_vcf_maf, preserve_barcode=True)
    w.write(vcf_variant, MockCounts(), sample_name="RNASample_01")
    w.close()

    with open(out_vcf_maf) as f:
        row = next(csv.DictReader(f, delimiter="\t"))
    assert row["Tumor_Sample_Barcode"] == "RNASample_01", (
        "VCF→MAF must always use BAM sample_name for Tumor_Sample_Barcode, "
        f"got: {row['Tumor_Sample_Barcode']!r}"
    )
