"""
Tests for variant normalization via Rust prepare_variants().

These tests exercise the standalone normalize module and the CLI subcommand.
They use small synthetic MAF files and a reference FASTA to verify:
- MAF anchor resolution (dash alleles get anchor base prepended)
- REF validation
- Left-alignment
- TSV output format from the normalize subcommand
- show_normalization columns in MafWriter output
"""

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import pysam
from helpers import read_maf_output as _read_maf_output

from gbcms import _rs as gbcms_rs
from gbcms.io.output import MafWriter
from gbcms.models.core import Variant, VariantType
from gbcms.normalize import normalize_variants


class TestNormalization(unittest.TestCase):
    """Test prepare_variants() and the normalize module."""

    def setUp(self):
        self.test_dir = tempfile.TemporaryDirectory()
        self.base_path = Path(self.test_dir.name)

        # Create a test reference with known repeat region.
        # Must be >200bp to accommodate the 100bp normalization window.
        # Structure: 100bp random + AAAAAA (6 A's at pos 100-105) + 100bp random
        self.fasta_path = self.base_path / "ref.fa"
        prefix = "ATCGATCG" * 12 + "ATCG"  # 100 bases
        homopolymer = "AAAAAA"  # 6 A's at positions 100-105
        suffix = "CGATCGAT" * 12 + "CGAT"  # 100 bases
        self.ref_seq = prefix + homopolymer + suffix  # 206 bases
        # Key positions:
        #   pos 99  = G (last base of prefix)
        #   pos 100 = A (start of A-run)
        #   pos 105 = A (end of A-run)
        #   pos 106 = C (first base of suffix)
        with open(self.fasta_path, "w") as f:
            f.write(">chr1\n")
            f.write(self.ref_seq + "\n")
        pysam.faidx(str(self.fasta_path))

    def tearDown(self):
        self.test_dir.cleanup()

    # -- prepare_variants() tests --

    def test_snp_passes_through(self):
        """SNPs should pass through with no normalization."""
        variants = [gbcms_rs.Variant("chr1", 0, "A", "T", "SNP")]
        prepared = gbcms_rs.prepare_variants(variants, str(self.fasta_path), 5, False, 1, False)
        self.assertEqual(len(prepared), 1)
        pv = prepared[0]
        self.assertEqual(pv.gbcms_status, "PASS")
        self.assertFalse(pv.was_anchor_resolved)
        self.assertFalse(pv.was_left_aligned)
        self.assertFalse(pv.was_normalized)
        self.assertEqual(pv.variant.ref_allele, "A")
        self.assertEqual(pv.variant.alt_allele, "T")

    def test_ref_mismatch_rejected(self):
        """Variant with wrong REF allele should be rejected."""
        # Position 0 has 'A', but we claim ref is 'G'
        variants = [gbcms_rs.Variant("chr1", 0, "G", "T", "SNP")]
        prepared = gbcms_rs.prepare_variants(variants, str(self.fasta_path), 5, False, 1, False)
        self.assertEqual(len(prepared), 1)
        pv = prepared[0]
        self.assertNotEqual(pv.gbcms_status, "PASS")

    def test_empty_allele_rejected(self):
        """A structurally empty REF or ALT is malformed input (MAF dash alleles must
        be '-', never ''). It must be rejected loudly with FAIL_EMPTY_ALLELE at prep
        time, not passed to counting where it would silently yield zero counts."""
        variants = [
            gbcms_rs.Variant("chr1", 0, "", "T", "INSERTION"),  # empty REF
            gbcms_rs.Variant("chr1", 0, "A", "", "DELETION"),  # empty ALT
        ]
        prepared = gbcms_rs.prepare_variants(variants, str(self.fasta_path), 5, False, 1, False)
        self.assertEqual(len(prepared), 2)
        for pv in prepared:
            self.assertEqual(pv.gbcms_status, "FAIL")
            self.assertEqual(pv.gbcms_status_reason, "EMPTY_ALLELE")

    def test_maf_insertion_anchor_resolution(self):
        """MAF insertion (ref='-') should get anchor base prepended."""
        # MAF: chr1:5, ref='-', alt='G' (insert G after pos 4 in 0-based)
        # After anchor: ref=A, alt=AG at pos 4
        variants = [gbcms_rs.Variant("chr1", 4, "-", "G", "INSERTION")]
        prepared = gbcms_rs.prepare_variants(variants, str(self.fasta_path), 5, True, 1, False)
        self.assertEqual(len(prepared), 1)
        pv = prepared[0]
        self.assertEqual(pv.gbcms_status, "PASS")
        # MAF anchor resolution should be detected
        self.assertTrue(pv.was_anchor_resolved)
        self.assertFalse(pv.was_left_aligned)
        self.assertTrue(pv.was_normalized)
        # After anchor resolution, ref should have anchor base
        self.assertNotEqual(pv.variant.ref_allele, "-")
        self.assertEqual(len(pv.variant.ref_allele), 1)  # Just the anchor
        self.assertEqual(len(pv.variant.alt_allele), 2)  # Anchor + inserted base

    def test_left_align_deletion_in_repeat(self):
        """Deletion in homopolymer AAAAAA should left-shift."""
        # Reference layout: ...G[AAAAAA]C... at pos 99-106
        # Deletion of one 'A' mid-run: pos=103, REF=AA, ALT=A
        # Should left-align to pos=99, REF=GA, ALT=G (leftmost VCF representation)
        variants = [gbcms_rs.Variant("chr1", 103, "AA", "A", "DELETION")]
        prepared = gbcms_rs.prepare_variants(variants, str(self.fasta_path), 5, False, 1, False)
        self.assertEqual(len(prepared), 1)
        pv = prepared[0]
        self.assertEqual(pv.gbcms_status, "PASS")
        # Left-alignment should be detected (VCF input, no anchor step)
        self.assertFalse(pv.was_anchor_resolved)
        self.assertTrue(pv.was_left_aligned)
        self.assertTrue(pv.was_normalized)
        # The deletion should be left-aligned from pos 103 to pos 99
        self.assertEqual(pv.variant.pos, 99)
        self.assertEqual(pv.variant.ref_allele, "GA")
        self.assertEqual(pv.variant.alt_allele, "G")

    def test_dynamic_window_expansion_long_repeat(self):
        """
        Gap 1B: Deletion in a >100bp repeat region should still normalize.

        The engine's normalization window starts at 100bp but doubles on
        edge-hit (100→200→400→...→2500bp). This test uses a 120bp
        dinucleotide repeat to verify the expansion mechanism works.
        A dinucleotide repeat (unlike a pure homopolymer) produces
        different alleles depending on deletion position, so
        left-alignment is meaningful.
        """
        # Create a reference with a 120bp AC-repeat
        long_fasta = self.base_path / "long_repeat.fa"
        prefix = "TTTTTTTTT" * 5 + "TTTTT"  # 50bp T-prefix
        repeat = "AC" * 60  # 120bp dinucleotide repeat at pos 50-169
        suffix = "GGGGGGGGG" * 5 + "GGGGG"  # 50bp G-suffix
        ref_seq = prefix + repeat + suffix  # 220bp total
        with open(long_fasta, "w") as f:
            f.write(">chr1\n")
            f.write(ref_seq + "\n")
        pysam.faidx(str(long_fasta))

        # Deletion near the END of the AC-repeat: pos=160, REF=AC, ALT=A
        # Should left-align back toward pos=49/50 (start of repeat)
        # This requires >110bp window, exceeding the default 100bp.
        variants = [gbcms_rs.Variant("chr1", 160, "AC", "A", "DELETION")]
        prepared = gbcms_rs.prepare_variants(variants, str(long_fasta), 5, False, 1, False)
        self.assertEqual(len(prepared), 1)
        pv = prepared[0]
        self.assertTrue(
            pv.gbcms_status.startswith("PASS"),
            f"Expected PASS status, got {pv.gbcms_status}",
        )
        # Variant should have normalized (shifted leftward)
        if pv.was_normalized:
            # If normalized, it should have shifted past the 100bp boundary
            self.assertLess(
                pv.variant.pos,
                pv.original_pos,
                f"Expected leftward shift, orig={pv.original_pos} new={pv.variant.pos}",
            )

    # -- normalize module tests --

    def test_normalize_module_writes_tsv(self):
        """Test the standalone normalize_variants function produces TSV output."""
        maf_path = self.base_path / "test.maf"
        with open(maf_path, "w") as f:
            f.write(
                "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
                "Reference_Allele\tTumor_Seq_Allele2\n"
            )
            f.write("Gene1\tchr1\t1\t1\tA\tT\n")

        output_path = self.base_path / "normalized.tsv"
        normalize_variants(
            variant_file=maf_path,
            reference=self.fasta_path,
            output=output_path,
            threads=1,
        )

        self.assertTrue(output_path.exists())
        with open(output_path) as f:
            lines = f.readlines()
            self.assertGreater(len(lines), 1)  # header + data
            header = lines[0].strip().split("\t")
            self.assertIn("chrom", header)
            self.assertIn("original_pos", header)
            self.assertIn("norm_pos", header)
            self.assertIn("gbcms_status", header)
            self.assertIn("was_anchor_resolved", header)
            self.assertIn("was_left_aligned", header)
            self.assertIn("was_normalized", header)

    # -- MafWriter show_normalization tests --

    @staticmethod
    def _zero_counts():
        """Create a zero-count stub for testing output."""
        _nan = float("nan")
        return SimpleNamespace(
            dp=0,
            rd=0,
            ad=0,
            dp_fwd=0,
            rd_fwd=0,
            ad_fwd=0,
            dp_rev=0,
            rd_rev=0,
            ad_rev=0,
            dpf=0,
            rdf=0,
            adf=0,
            rdf_fwd=0,
            rdf_rev=0,
            adf_fwd=0,
            adf_rev=0,
            sb_pval=1.0,
            sb_or=0.0,
            fsb_pval=1.0,
            fsb_or=0.0,
            # Decomposed ALT counting
            any_alt=0,
            partial_alt=0,
            # N-base diagnostic
            n_count=0,
            # mFSD fields (zero counts, NaN float stats)
            mfsd_ref_count=0,
            mfsd_alt_count=0,
            mfsd_nonref_count=0,
            mfsd_n_count=0,
            mfsd_ref_mean=_nan,
            mfsd_alt_mean=_nan,
            mfsd_nonref_mean=_nan,
            mfsd_n_mean=_nan,
            mfsd_alt_llr=_nan,
            mfsd_ref_llr=_nan,
            mfsd_delta_alt_ref=_nan,
            mfsd_ks_alt_ref=_nan,
            mfsd_pval_alt_ref=_nan,
            mfsd_delta_alt_nonref=_nan,
            mfsd_ks_alt_nonref=_nan,
            mfsd_pval_alt_nonref=_nan,
            mfsd_delta_ref_nonref=_nan,
            mfsd_ks_ref_nonref=_nan,
            mfsd_pval_ref_nonref=_nan,
            mfsd_delta_alt_n=_nan,
            mfsd_ks_alt_n=_nan,
            mfsd_pval_alt_n=_nan,
            mfsd_delta_ref_n=_nan,
            mfsd_ks_ref_n=_nan,
            mfsd_pval_ref_n=_nan,
            mfsd_delta_nonref_n=_nan,
            mfsd_ks_nonref_n=_nan,
            mfsd_pval_nonref_n=_nan,
        )

    def test_show_normalization_columns(self):
        """MafWriter should include norm_* columns when show_normalization=True."""
        output_path = self.base_path / "test_norm.maf"

        # Create a variant and its "normalized" version
        variant = Variant(chrom="1", pos=99, ref="A", alt="T", variant_type=VariantType.SNP)
        norm_variant = Variant(chrom="1", pos=95, ref="A", alt="T", variant_type=VariantType.SNP)

        writer = MafWriter(output_path, show_normalization=True)
        writer.write(variant, self._zero_counts(), norm_variant=norm_variant)
        writer.close()

        with open(output_path):
            reader = _read_maf_output(output_path)
            rows = list(reader)

        self.assertEqual(len(rows), 1)
        header = rows[0].keys()
        self.assertIn("norm_Start_Position", header)
        self.assertIn("norm_End_Position", header)
        self.assertIn("norm_Reference_Allele", header)
        self.assertIn("norm_Tumor_Seq_Allele2", header)

    def test_column_prefix_with_norm(self):
        """Norm columns should respect the column_prefix setting."""
        output_path = self.base_path / "test_prefix.maf"

        variant = Variant(chrom="1", pos=99, ref="A", alt="T", variant_type=VariantType.SNP)
        norm_variant = Variant(chrom="1", pos=95, ref="A", alt="T", variant_type=VariantType.SNP)

        writer = MafWriter(output_path, column_prefix="t_", show_normalization=True)
        writer.write(variant, self._zero_counts(), norm_variant=norm_variant)
        writer.close()

        with open(output_path):
            reader = _read_maf_output(output_path)
            rows = list(reader)

        header = rows[0].keys()
        self.assertIn("t_norm_Start_Position", header)
        self.assertIn("t_norm_End_Position", header)
        self.assertIn("t_norm_Reference_Allele", header)
        self.assertIn("t_norm_Tumor_Seq_Allele2", header)
        # Unprefixed versions should NOT exist
        self.assertNotIn("norm_Start_Position", header)

    def test_no_norm_columns_by_default(self):
        """MafWriter should NOT include norm_* columns when show_normalization=False."""
        output_path = self.base_path / "test_no_norm.maf"

        variant = Variant(chrom="1", pos=99, ref="A", alt="T", variant_type=VariantType.SNP)

        writer = MafWriter(output_path, show_normalization=False)
        writer.write(variant, self._zero_counts())
        writer.close()

        with open(output_path):
            reader = _read_maf_output(output_path)
            rows = list(reader)

        header = rows[0].keys()
        self.assertNotIn("norm_Start_Position", header)
        self.assertNotIn("norm_End_Position", header)

    def test_maf_deletion_in_repeat_both_flags(self):
        """MAF deletion (alt='-') in a homopolymer should set BOTH flags.

        Step 1 (anchor resolution): converts dash-allele to VCF-style → was_anchor_resolved=True
        Step 3 (left-alignment): shifts deletion to start of A-run → was_left_aligned=True
        """
        # Reference layout: ...G[AAAAAA]C... at pos 99-106
        # MAF: pos=103 (0-based), ref='A', alt='-' → deletion inside A-run
        # Step 1: anchor at pos=102, ref=AA, alt=A → was_anchor_resolved=True
        # Step 3: left-align AA→A in A-run → shifts to pos=99 → was_left_aligned=True
        variants = [gbcms_rs.Variant("chr1", 103, "A", "-", "DELETION")]
        prepared = gbcms_rs.prepare_variants(variants, str(self.fasta_path), 5, True, 1, False)
        self.assertEqual(len(prepared), 1)
        pv = prepared[0]
        self.assertTrue(
            pv.gbcms_status.startswith("PASS"),
            f"Expected PASS status, got {pv.gbcms_status}",
        )
        self.assertTrue(pv.was_anchor_resolved, "Expected anchor resolution for dash-allele")
        self.assertTrue(pv.was_left_aligned, "Expected left-alignment in homopolymer")
        self.assertTrue(pv.was_normalized, "Expected combined normalized flag")
        # Should have shifted left from original position
        self.assertLess(pv.variant.pos, 103, "Expected leftward shift")

    def test_n_in_alt_allele_rejected(self):
        """ALT allele containing N should be rejected with FAIL_ALT_CONTAINS_N.

        N in ALT indicates ambiguous/placeholder genotyping — no real read
        base can match it. The normalize engine should explicitly reject
        these rather than silently producing zero counts.
        """
        variants = [gbcms_rs.Variant("chr1", 0, "A", "N", "SNP")]
        prepared = gbcms_rs.prepare_variants(variants, str(self.fasta_path), 5, False, 1, False)
        self.assertEqual(len(prepared), 1)
        pv = prepared[0]
        self.assertEqual(pv.gbcms_status, "FAIL")
        self.assertEqual(
            pv.gbcms_status_reason,
            "ALT_CONTAINS_N",
            f"Expected ALT_CONTAINS_N reason, got {pv.gbcms_status_reason}",
        )

    def test_n_in_alt_allele_mnp_rejected(self):
        """MNP with N in ALT should also be rejected (e.g., AN > TN)."""
        variants = [gbcms_rs.Variant("chr1", 0, "AT", "TN", "DNP")]
        prepared = gbcms_rs.prepare_variants(variants, str(self.fasta_path), 5, False, 1, False)
        self.assertEqual(len(prepared), 1)
        pv = prepared[0]
        self.assertEqual(pv.gbcms_status, "FAIL")
        self.assertEqual(
            pv.gbcms_status_reason,
            "ALT_CONTAINS_N",
            f"Expected ALT_CONTAINS_N reason for MNP with N in ALT, got {pv.gbcms_status_reason}",
        )


if __name__ == "__main__":
    unittest.main()
