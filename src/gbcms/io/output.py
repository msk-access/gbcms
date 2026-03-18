"""
Output Writers: Formatting results for VCF and MAF.

This module provides classes to write processed variants and their counts
to output files, handling format-specific columns and headers.
"""

import csv
import logging
from pathlib import Path
from typing import Any

from ..core.kernel import CoordinateKernel
from ..models.core import Variant

__all__ = ["OutputWriter", "MafWriter", "VcfWriter"]

logger = logging.getLogger(__name__)


def _fmt(v: float) -> str:
    """Format a float for MAF output. NaN → 'NA' (standard missing value for tabular formats)."""
    import math

    return "NA" if (isinstance(v, float) and math.isnan(v)) else f"{v:.4f}"


def _fmt_vcf(v: float) -> str:
    """Format a float for VCF INFO fields. NaN → '.' (VCF spec missing value sentinel)."""
    import math

    return "." if (isinstance(v, float) and math.isnan(v)) else f"{v:.4f}"


class OutputWriter:
    """Abstract base class for output writers."""

    def write(self, variant: Variant, counts: Any):
        raise NotImplementedError

    def close(self):
        pass


class MafWriter(OutputWriter):
    """
    Writes results to a MAF-like file (Fillout format).

    Supports two output strategies based on input format:
    - MAF→MAF: Preserves all original MAF columns, appends gbcms count columns.
      Original positional/allele/type columns are NEVER overwritten.
    - VCF→MAF: Generates GDC-compliant MAF coordinates from internal VCF-style
      representation using CoordinateKernel.internal_to_maf().

    Count column names are controlled by the column_prefix parameter:
    - Default (empty): 'ref_count', 'alt_count', 'total_count', etc.
    - Legacy ('t_'):   't_ref_count', 't_alt_count', 't_total_count', etc.
    """

    # Default MAF columns for VCF→MAF output (minimal GDC-compatible set)
    _DEFAULT_MAF_HEADERS = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Strand",
        "Variant_Classification",
        "Variant_Type",
        "Reference_Allele",
        "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2",
        "Tumor_Sample_Barcode",
        "Matched_Norm_Sample_Barcode",
    ]

    def __init__(
        self,
        path: Path,
        column_prefix: str = "",
        preserve_barcode: bool = False,
        show_normalization: bool = False,
        mfsd: bool = False,
        mode: str = "dna",
    ):
        """
        Initialize MafWriter.

        Args:
            path: Output file path.
            column_prefix: Prefix for gbcms count columns (e.g., '', 't_', 'gbcms_').
            preserve_barcode: If True, keep original Tumor_Sample_Barcode from
                input MAF. If False (default), override with BAM sample name.
                Only applies to MAF→MAF; VCF→MAF always uses BAM name.
            show_normalization: If True, append norm_* columns showing
                left-aligned coordinates in the output.
            mfsd: If True, append 31 mFSD columns (KS test, LLR, mean sizes,
                pairwise comparisons, derived metrics). Controlled by --mfsd flag.
            mode: Counting mode ('dna' or 'rna'). In RNA mode, 5 RNA-specific
                columns are appended instead of mFSD columns.
        """
        self.path = path
        self.column_prefix = column_prefix
        self.preserve_barcode = preserve_barcode
        self.show_normalization = show_normalization
        self.mfsd = mfsd
        self.mode = mode
        self.file = open(path, "w")
        self.writer: csv.DictWriter | None = None
        self._headers_written = False
        logger.debug(
            "MafWriter initialized: path=%s, column_prefix='%s', "
            "preserve_barcode=%s, show_normalization=%s, mfsd=%s, mode=%s",
            path,
            column_prefix,
            preserve_barcode,
            show_normalization,
            mfsd,
            mode,
        )

    def _gbcms_column_names(self) -> list[str]:
        """
        Build the list of gbcms-generated count column names with the configured prefix.

        Returns:
            Ordered list of gbcms column names.
        """
        p = self.column_prefix
        cols = [
            # Validation status
            "validation_status",
            # Core counts
            f"{p}ref_count",
            f"{p}alt_count",
            f"{p}total_count",
            f"{p}vaf",
            # Fragment counts
            f"{p}ref_count_fragment",
            f"{p}alt_count_fragment",
            f"{p}total_count_fragment",
            f"{p}vaf_fragment",
            # Strand bias (unprefixed — always unique)
            "strand_bias_p_value",
            "strand_bias_odds_ratio",
            "fragment_strand_bias_p_value",
            "fragment_strand_bias_odds_ratio",
            # Strand counts
            f"{p}ref_count_forward",
            f"{p}ref_count_reverse",
            f"{p}alt_count_forward",
            f"{p}alt_count_reverse",
            f"{p}ref_count_fragment_forward",
            f"{p}ref_count_fragment_reverse",
            f"{p}alt_count_fragment_forward",
            f"{p}alt_count_fragment_reverse",
        ]
        if self.mfsd:
            # ── mFSD: Mutant Fragment Size Distribution (31 columns) ──────────
            # Only appended when --mfsd is set. Without the flag these columns
            # are completely absent from output (not NA-filled or zero-filled).
            cols += [
                # Raw counts — fragments in each class with valid insert size (50–1000 bp)
                "mfsd_ref_count",
                "mfsd_alt_count",
                "mfsd_nonref_count",
                "mfsd_n_count",
                # LLR relative to healthy/tumor cfDNA Gaussian model
                "mfsd_alt_llr",
                "mfsd_ref_llr",
                # Mean fragment size per class (bp)
                "mfsd_ref_mean",
                "mfsd_alt_mean",
                "mfsd_nonref_mean",
                "mfsd_n_mean",
                # Pairwise KS comparisons: 6 pairs × (delta, D-stat, p-value)
                "mfsd_delta_alt_ref",
                "mfsd_ks_alt_ref",
                "mfsd_pval_alt_ref",
                "mfsd_delta_alt_nonref",
                "mfsd_ks_alt_nonref",
                "mfsd_pval_alt_nonref",
                "mfsd_delta_ref_nonref",
                "mfsd_ks_ref_nonref",
                "mfsd_pval_ref_nonref",
                "mfsd_delta_alt_n",
                "mfsd_ks_alt_n",
                "mfsd_pval_alt_n",
                "mfsd_delta_ref_n",
                "mfsd_ks_ref_n",
                "mfsd_pval_ref_n",
                "mfsd_delta_nonref_n",
                "mfsd_ks_nonref_n",
                "mfsd_pval_nonref_n",
                # Derived quality metrics (computed in Python from Rust exports)
                "mfsd_error_rate",
                "mfsd_n_rate",
                "mfsd_size_ratio",
                "mfsd_quality_score",
                "mfsd_alt_confidence",
                "mfsd_ks_valid",
            ]
        if self.show_normalization:
            cols.extend(self._norm_column_names())

        # ── RNA-specific columns (5) ──────────────────────────────────────────
        # Only appended in RNA mode. These replace mFSD as the RNA-specific
        # diagnostic columns. When mode="dna" these are completely absent.
        if self.mode == "rna":
            cols += [
                "rna_sense_depth",
                "rna_antisense_depth",
                "rna_alt_sense_count",
                "rna_editing_site",
                "rna_splice_spanning",
            ]

        return cols

    def _norm_column_names(self) -> list[str]:
        """Normalization columns (only appended when --show-normalization)."""
        p = self.column_prefix
        return [
            f"{p}norm_Start_Position",
            f"{p}norm_End_Position",
            f"{p}norm_Reference_Allele",
            f"{p}norm_Tumor_Seq_Allele2",
        ]

    def _init_writer(self, original_headers: list[str]) -> None:
        """
        Initialize the CSV writer with dynamically constructed headers.

        Header order: original MAF columns first, then gbcms columns appended.
        Duplicate column names (already present in original) are skipped.

        Args:
            original_headers: Column names from the input MAF (or defaults for VCF→MAF).
        """
        gbcms_cols = self._gbcms_column_names()
        existing = set(original_headers)

        # Only append gbcms columns not already in the original headers
        new_cols = [c for c in gbcms_cols if c not in existing]
        self.fieldnames = list(original_headers) + new_cols

        self.writer = csv.DictWriter(
            self.file,
            fieldnames=self.fieldnames,
            delimiter="\t",
            extrasaction="ignore",
        )
        self.writer.writeheader()
        self._headers_written = True

        logger.debug(
            "MafWriter headers: %d original + %d gbcms = %d total columns",
            len(original_headers),
            len(new_cols),
            len(self.fieldnames),
        )

    def _populate_gbcms_counts(self, counts: Any) -> dict[str, str]:
        """
        Build the gbcms count columns dictionary with the configured prefix.

        Calculates VAF values and formats all required count data as strings.
        mFSD derived metrics are only computed when self.mfsd is True — skipping
        the NaN arithmetic and attribute accesses when mFSD is disabled.

        Args:
            counts: BaseCounts object from the Rust engine.

        Returns:
            Dictionary mapping prefixed column names to string values.
        """
        p = self.column_prefix

        # Calculate VAFs with zero-division protection
        total_reads = counts.rd + counts.ad
        vaf = counts.ad / total_reads if total_reads > 0 else 0.0

        total_frags = counts.rdf + counts.adf
        vaf_frag = counts.adf / total_frags if total_frags > 0 else 0.0

        result: dict[str, str] = {
            # Core counts
            f"{p}ref_count": str(counts.rd),
            f"{p}alt_count": str(counts.ad),
            f"{p}total_count": str(counts.dp),
            f"{p}vaf": f"{vaf:.4f}",
            # Fragment counts
            f"{p}ref_count_fragment": str(counts.rdf),
            f"{p}alt_count_fragment": str(counts.adf),
            f"{p}total_count_fragment": str(counts.dpf),
            f"{p}vaf_fragment": f"{vaf_frag:.4f}",
            # Strand bias (unprefixed)
            "strand_bias_p_value": f"{counts.sb_pval:.4e}",
            "strand_bias_odds_ratio": f"{counts.sb_or:.4f}",
            "fragment_strand_bias_p_value": f"{counts.fsb_pval:.4e}",
            "fragment_strand_bias_odds_ratio": f"{counts.fsb_or:.4f}",
            # Strand counts
            f"{p}ref_count_forward": str(counts.rd_fwd),
            f"{p}ref_count_reverse": str(counts.rd_rev),
            f"{p}alt_count_forward": str(counts.ad_fwd),
            f"{p}alt_count_reverse": str(counts.ad_rev),
            f"{p}ref_count_fragment_forward": str(counts.rdf_fwd),
            f"{p}ref_count_fragment_reverse": str(counts.rdf_rev),
            f"{p}alt_count_fragment_forward": str(counts.adf_fwd),
            f"{p}alt_count_fragment_reverse": str(counts.adf_rev),
        }

        # ── RNA-specific count columns ─────────────────────────────────────────
        if self.mode == "rna":
            result.update(
                {
                    "rna_sense_depth": str(counts.sense_depth),
                    "rna_antisense_depth": str(counts.antisense_depth),
                    "rna_alt_sense_count": str(counts.sense_strand_alt_count),
                    "rna_editing_site": str(counts.rna_editing_site_overlap),
                    "rna_splice_spanning": str(counts.splice_spanning_count),
                }
            )

        if self.mfsd:
            # ── mFSD derived metrics ───────────────────────────────────────────
            # Computed here in Python from already-exported Rust counts.
            # Only computed when --mfsd is set — avoids NaN arithmetic overhead
            # on every variant when mFSD analysis is not requested.
            _nan = float("nan")
            total_mfsd = (
                counts.mfsd_ref_count
                + counts.mfsd_alt_count
                + counts.mfsd_nonref_count
                + counts.mfsd_n_count
            )
            mfsd_error_rate = counts.mfsd_nonref_count / total_mfsd if total_mfsd > 0 else _nan
            mfsd_n_rate = counts.mfsd_n_count / total_mfsd if total_mfsd > 0 else _nan
            # Size ratio: mean(ALT) / mean(REF); NaN if either is 0/missing
            mfsd_size_ratio = (
                counts.mfsd_alt_mean / counts.mfsd_ref_mean
                if counts.mfsd_ref_mean > 0 and counts.mfsd_alt_count > 0
                else _nan
            )
            # Quality score: 1 - error_rate - n_rate; NaN if either is NaN
            mfsd_quality_score = (
                1.0 - mfsd_n_rate - mfsd_error_rate
                if not (mfsd_n_rate != mfsd_n_rate or mfsd_error_rate != mfsd_error_rate)
                else _nan
            )
            # Categorical confidence based on ALT fragment count
            if counts.mfsd_alt_count >= 5:
                mfsd_alt_confidence = "HIGH"
            elif counts.mfsd_alt_count >= 1:
                mfsd_alt_confidence = "LOW"
            else:
                mfsd_alt_confidence = "NONE"
            # KS test validity: both ALT and REF need >= 5 fragments
            mfsd_ks_valid = counts.mfsd_alt_count >= 5 and counts.mfsd_ref_count >= 5

            result.update(
                {
                    # Raw counts
                    "mfsd_ref_count": str(counts.mfsd_ref_count),
                    "mfsd_alt_count": str(counts.mfsd_alt_count),
                    "mfsd_nonref_count": str(counts.mfsd_nonref_count),
                    "mfsd_n_count": str(counts.mfsd_n_count),
                    # LLR
                    "mfsd_alt_llr": _fmt(counts.mfsd_alt_llr),
                    "mfsd_ref_llr": _fmt(counts.mfsd_ref_llr),
                    # Mean sizes
                    "mfsd_ref_mean": _fmt(counts.mfsd_ref_mean),
                    "mfsd_alt_mean": _fmt(counts.mfsd_alt_mean),
                    "mfsd_nonref_mean": _fmt(counts.mfsd_nonref_mean),
                    "mfsd_n_mean": _fmt(counts.mfsd_n_mean),
                    # Pairwise KS comparisons: 6 pairs × 3 values
                    "mfsd_delta_alt_ref": _fmt(counts.mfsd_delta_alt_ref),
                    "mfsd_ks_alt_ref": _fmt(counts.mfsd_ks_alt_ref),
                    "mfsd_pval_alt_ref": _fmt(counts.mfsd_pval_alt_ref),
                    "mfsd_delta_alt_nonref": _fmt(counts.mfsd_delta_alt_nonref),
                    "mfsd_ks_alt_nonref": _fmt(counts.mfsd_ks_alt_nonref),
                    "mfsd_pval_alt_nonref": _fmt(counts.mfsd_pval_alt_nonref),
                    "mfsd_delta_ref_nonref": _fmt(counts.mfsd_delta_ref_nonref),
                    "mfsd_ks_ref_nonref": _fmt(counts.mfsd_ks_ref_nonref),
                    "mfsd_pval_ref_nonref": _fmt(counts.mfsd_pval_ref_nonref),
                    "mfsd_delta_alt_n": _fmt(counts.mfsd_delta_alt_n),
                    "mfsd_ks_alt_n": _fmt(counts.mfsd_ks_alt_n),
                    "mfsd_pval_alt_n": _fmt(counts.mfsd_pval_alt_n),
                    "mfsd_delta_ref_n": _fmt(counts.mfsd_delta_ref_n),
                    "mfsd_ks_ref_n": _fmt(counts.mfsd_ks_ref_n),
                    "mfsd_pval_ref_n": _fmt(counts.mfsd_pval_ref_n),
                    "mfsd_delta_nonref_n": _fmt(counts.mfsd_delta_nonref_n),
                    "mfsd_ks_nonref_n": _fmt(counts.mfsd_ks_nonref_n),
                    "mfsd_pval_nonref_n": _fmt(counts.mfsd_pval_nonref_n),
                    # Derived quality metrics
                    "mfsd_error_rate": _fmt(mfsd_error_rate),
                    "mfsd_n_rate": _fmt(mfsd_n_rate),
                    "mfsd_size_ratio": _fmt(mfsd_size_ratio),
                    "mfsd_quality_score": _fmt(mfsd_quality_score),
                    "mfsd_alt_confidence": mfsd_alt_confidence,
                    "mfsd_ks_valid": str(mfsd_ks_valid),
                }
            )

        return result

    def write(
        self,
        variant: Variant,
        counts: Any,
        sample_name: str = "TUMOR",
        validation_status: str = "PASS",
        norm_variant: Variant | None = None,
    ) -> None:
        """
        Write a single variant row to the MAF output.

        Two output strategies:
        - MAF→MAF (variant.metadata populated): Pass through all original columns,
          append gbcms count columns. Original values are NEVER overwritten.
        - VCF→MAF (no metadata): Generate GDC-compliant MAF coordinates from
          internal representation using CoordinateKernel.internal_to_maf().

        Args:
            variant: Normalized Variant with optional metadata from input MAF.
            counts: BaseCounts object from the Rust engine.
            sample_name: Sample name for Tumor_Sample_Barcode column.
            validation_status: Validation status from prepare_variants().
            norm_variant: Optional left-aligned Variant (for --show-normalization).
        """
        # Initialize writer on first variant (headers depend on input format)
        if not self._headers_written:
            if variant.metadata:
                # MAF→MAF: use original input headers
                self._init_writer(list(variant.metadata.keys()))
            else:
                # VCF→MAF: use default GDC MAF headers + VCF-origin fields
                vcf_headers = self._DEFAULT_MAF_HEADERS + [
                    "vcf_id",
                    "vcf_pos",
                    "vcf_region",
                ]
                self._init_writer(vcf_headers)

        assert self.writer is not None

        # Build the output row based on input format
        if variant.metadata:
            # MAF→MAF: start with ALL original metadata (preserves every column)
            row = dict(variant.metadata)
        else:
            # VCF→MAF: build row from internal representation
            row = dict.fromkeys(self.fieldnames, "")

            # Convert internal coordinates to GDC MAF format
            maf_coords = CoordinateKernel.internal_to_maf(variant)
            row.update(maf_coords)
            row["Chromosome"] = variant.chrom

            # VCF-origin tracking fields
            vcf_pos = variant.pos + 1
            row["vcf_pos"] = str(vcf_pos)
            row["vcf_region"] = f"{variant.chrom}:{vcf_pos}"
            if variant.original_id:
                row["vcf_id"] = variant.original_id

        # Set sample barcode:
        # - VCF→MAF: always use BAM sample name (no barcode in VCF)
        # - MAF→MAF + preserve_barcode: keep original from input metadata
        # - MAF→MAF + no preserve_barcode: override with BAM sample name
        if not (variant.metadata and self.preserve_barcode):
            row["Tumor_Sample_Barcode"] = sample_name

        # Append gbcms count columns (both paths, never overwrites originals)
        row["validation_status"] = validation_status
        row.update(self._populate_gbcms_counts(counts))

        # Normalization columns (only when --show-normalization is enabled)
        if self.show_normalization and norm_variant:
            maf_norm = CoordinateKernel.internal_to_maf(norm_variant)
            p = self.column_prefix
            row[f"{p}norm_Start_Position"] = maf_norm["Start_Position"]
            row[f"{p}norm_End_Position"] = maf_norm["End_Position"]
            row[f"{p}norm_Reference_Allele"] = maf_norm["Reference_Allele"]
            row[f"{p}norm_Tumor_Seq_Allele2"] = maf_norm["Tumor_Seq_Allele2"]

        self.writer.writerow(row)

    def close(self) -> None:
        """Close the output file."""
        self.file.close()
        logger.debug("MafWriter closed: %s", self.path)


class VcfWriter(OutputWriter):
    """Writes results to a VCF file."""

    def __init__(
        self,
        path: Path,
        sample_name: str = "SAMPLE",
        show_normalization: bool = False,
        mfsd: bool = False,
        mode: str = "dna",
    ):
        self.path = path
        self.sample_name = sample_name
        self.show_normalization = show_normalization
        self.mfsd = mfsd
        self.mode = mode
        self.file = open(path, "w")
        self._headers_written = False
        logger.debug(
            "VcfWriter initialized: path=%s, sample=%s, show_normalization=%s, mfsd=%s, mode=%s",
            path,
            sample_name,
            show_normalization,
            mfsd,
            mode,
        )

    def _write_header(self):
        """Write VCF header lines.

        mFSD ##INFO fields (7 lines) are only included when self.mfsd is True.
        """
        headers = [
            "##fileformat=VCFv4.2",
            "##source=gbcms_v2",
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
            '##INFO=<ID=VS,Number=1,Type=String,Description="Validation status from prepare_variants">',
            '##INFO=<ID=SB_PVAL,Number=1,Type=Float,Description="Fisher strand bias p-value">',
            '##INFO=<ID=SB_OR,Number=1,Type=Float,Description="Fisher strand bias odds ratio">',
            '##INFO=<ID=FSB_PVAL,Number=1,Type=Float,Description="Fisher fragment strand bias p-value">',
            '##INFO=<ID=FSB_OR,Number=1,Type=Float,Description="Fisher fragment strand bias odds ratio">',
        ]
        if self.mfsd:
            # mFSD INFO fields (7 primary diagnostics). VCF key = MAF column name uppercased.
            # Only added when --mfsd is set — keeps VCF header minimal for standard runs.
            headers.extend(
                [
                    '##INFO=<ID=MFSD_DELTA_ALT_REF,Number=1,Type=Float,Description="mFSD mean(ALT) − mean(REF) fragment size delta (bp)">',
                    '##INFO=<ID=MFSD_KS_ALT_REF,Number=1,Type=Float,Description="mFSD 2-sample KS D-statistic (ALT vs REF)">',
                    '##INFO=<ID=MFSD_PVAL_ALT_REF,Number=1,Type=Float,Description="mFSD KS p-value (ALT vs REF)">',
                    '##INFO=<ID=MFSD_ALT_LLR,Number=1,Type=Float,Description="mFSD LLR for ALT fragments: Σ log(P_tumor/P_healthy); positive=tumor-like">',
                    '##INFO=<ID=MFSD_REF_LLR,Number=1,Type=Float,Description="mFSD LLR for REF fragments">',
                    '##INFO=<ID=MFSD_ALT_COUNT,Number=1,Type=Integer,Description="ALT-classified fragments in mFSD window (50–1000 bp)">',
                    '##INFO=<ID=MFSD_REF_COUNT,Number=1,Type=Integer,Description="REF-classified fragments in mFSD window (50–1000 bp)">',
                ]
            )
        if self.show_normalization:
            headers.extend(
                [
                    '##INFO=<ID=NORM_POS,Number=1,Type=Integer,Description="Left-aligned VCF position (1-based)">',
                    '##INFO=<ID=NORM_REF,Number=1,Type=String,Description="Left-aligned REF allele">',
                    '##INFO=<ID=NORM_ALT,Number=1,Type=String,Description="Left-aligned ALT allele">',
                ]
            )
        # ── RNA-specific INFO/FORMAT headers ────────────────────────────────
        if self.mode == "rna":
            headers.extend(
                [
                    '##INFO=<ID=SEN,Number=1,Type=Integer,Description="Reads on the transcript sense strand">',
                    '##INFO=<ID=ANT,Number=1,Type=Integer,Description="Reads on the antisense strand">',
                    '##INFO=<ID=ASEN,Number=1,Type=Integer,Description="ALT reads on the transcript sense strand">',
                    '##INFO=<ID=RED,Number=0,Type=Flag,Description="Locus is a candidate A-to-I RNA editing site (A>G on + strand or T>C on - strand)">',
                    '##INFO=<ID=SPL,Number=1,Type=Integer,Description="ALT reads spanning a splice junction (CIGAR N)">',
                ]
            )
        headers.extend(
            [
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allelic depths for the ref and alt alleles (fwd,rev)">',
                '##FORMAT=<ID=DP,Number=2,Type=Integer,Description="Approximate read depth (ref_total,alt_total)">',
                '##FORMAT=<ID=RD,Number=2,Type=Integer,Description="Reference read depth (fwd,rev)">',
                '##FORMAT=<ID=RDF,Number=2,Type=Integer,Description="Ref Fragment Count (fwd,rev)">',
                '##FORMAT=<ID=ADF,Number=2,Type=Integer,Description="Alt Fragment Count (fwd,rev)">',
                '##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant Allele Fraction (read level)">',
                '##FORMAT=<ID=FAF,Number=1,Type=Float,Description="Variant Allele Fraction (fragment level)">',
            ]
        )
        if self.mode == "rna":
            headers.extend(
                [
                    '##FORMAT=<ID=SEN,Number=1,Type=Integer,Description="Sense strand depth">',
                    '##FORMAT=<ID=ANT,Number=1,Type=Integer,Description="Antisense strand depth">',
                    '##FORMAT=<ID=ASEN,Number=1,Type=Integer,Description="ALT sense strand count">',
                    '##FORMAT=<ID=SPL,Number=1,Type=Integer,Description="Splice-spanning ALT count">',
                ]
            )
        headers.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{self.sample_name}")
        self.file.write("\n".join(headers) + "\n")
        self._headers_written = True

    def write(
        self,
        variant: Variant,
        counts: Any,
        sample_name: str = "SAMPLE",
        validation_status: str = "PASS",
        norm_variant: Variant | None = None,
    ):
        if not self._headers_written:
            self._write_header()

        # VCF POS is 1-based
        pos = variant.pos + 1

        # INFO fields (VCF spec: missing values use '.' not 'NA')
        info_parts = [
            f"DP={counts.dp}",
            f"VS={validation_status}",
            f"SB_PVAL={counts.sb_pval:.4e}",
            f"SB_OR={counts.sb_or:.4f}",
            f"FSB_PVAL={counts.fsb_pval:.4e}",
            f"FSB_OR={counts.fsb_or:.4f}",
        ]
        if self.mfsd:
            # mFSD primary diagnostic INFO fields (7 values).
            # Only populated when --mfsd is set; '.' for NaN per VCF spec.
            info_parts.extend(
                [
                    f"MFSD_DELTA_ALT_REF={_fmt_vcf(counts.mfsd_delta_alt_ref)}",
                    f"MFSD_KS_ALT_REF={_fmt_vcf(counts.mfsd_ks_alt_ref)}",
                    f"MFSD_PVAL_ALT_REF={_fmt_vcf(counts.mfsd_pval_alt_ref)}",
                    f"MFSD_ALT_LLR={_fmt_vcf(counts.mfsd_alt_llr)}",
                    f"MFSD_REF_LLR={_fmt_vcf(counts.mfsd_ref_llr)}",
                    f"MFSD_ALT_COUNT={counts.mfsd_alt_count}",
                    f"MFSD_REF_COUNT={counts.mfsd_ref_count}",
                ]
            )
        if self.show_normalization and norm_variant:
            info_parts.extend(
                [
                    f"NORM_POS={norm_variant.pos + 1}",
                    f"NORM_REF={norm_variant.ref}",
                    f"NORM_ALT={norm_variant.alt}",
                ]
            )
        # ── RNA-specific INFO fields ──────────────────────────────────────
        if self.mode == "rna":
            info_parts.extend(
                [
                    f"SEN={counts.sense_depth}",
                    f"ANT={counts.antisense_depth}",
                    f"ASEN={counts.sense_strand_alt_count}",
                    f"SPL={counts.splice_spanning_count}",
                ]
            )
            if counts.rna_editing_site_overlap:
                info_parts.append("RED")
        info = ";".join(info_parts)

        # FORMAT fields
        gt = "0/1" if counts.ad > 0 else "0/0"
        dp = f"{counts.rd},{counts.ad}"
        rd = f"{counts.rd_fwd},{counts.rd_rev}"
        ad = f"{counts.ad_fwd},{counts.ad_rev}"
        rdf = f"{counts.rdf_fwd},{counts.rdf_rev}"
        adf = f"{counts.adf_fwd},{counts.adf_rev}"

        total_reads = counts.rd + counts.ad
        vaf = counts.ad / total_reads if total_reads > 0 else 0.0

        total_frags = counts.rdf + counts.adf
        faf = counts.adf / total_frags if total_frags > 0 else 0.0

        if self.mode == "rna":
            format_str = "GT:DP:RD:AD:RDF:ADF:VAF:FAF:SEN:ANT:ASEN:SPL"
            sample_data = (
                f"{gt}:{dp}:{rd}:{ad}:{rdf}:{adf}:{vaf:.4f}:{faf:.4f}"
                f":{counts.sense_depth}:{counts.antisense_depth}"
                f":{counts.sense_strand_alt_count}:{counts.splice_spanning_count}"
            )
        else:
            format_str = "GT:DP:RD:AD:RDF:ADF:VAF:FAF"
            sample_data = f"{gt}:{dp}:{rd}:{ad}:{rdf}:{adf}:{vaf:.4f}:{faf:.4f}"

        row = [
            variant.chrom,
            str(pos),
            variant.original_id or ".",
            variant.ref,
            variant.alt,
            ".",  # QUAL
            ".",  # FILTER
            info,
            format_str,
            sample_data,
        ]

        self.file.write("\t".join(row) + "\n")

    def close(self):
        self.file.close()
