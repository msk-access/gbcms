"""
Output Writers: Formatting results for VCF and MAF.

This module provides classes to write processed variants and their counts
to output files, handling format-specific columns and headers.
"""

import csv
import logging
import math
from pathlib import Path
from typing import Any

from ..core.kernel import CoordinateKernel
from ..models.core import Variant

__all__ = ["OutputWriter", "MafWriter", "VcfWriter"]

logger = logging.getLogger(__name__)

# ── Clonal Hematopoiesis (CH) gene set ──────────────────────────────────────
# Well-established CH driver genes. Variants in these genes flagged as
# "CH-associated" in mFSD output for CH-vs-ctDNA fragment size interpretation.
# Source: Steensma et al. (2015), Jaiswal et al. (2014), Bolton et al. (2020).
CH_GENES: frozenset[str] = frozenset(
    {
        "DNMT3A",
        "TET2",
        "ASXL1",
        "PPM1D",
        "TP53",
        "JAK2",
        "SF3B1",
        "SRSF2",
        "U2AF1",
        "CBL",
        "ATM",
        "BCOR",
        "EZH2",
        "IDH1",
        "IDH2",
        "GNAS",
        "GNB1",
        "BCORL1",
        "SETD2",
        "STAG2",
    }
)


def _fmt(v: float) -> str:
    """Format a float for MAF output.

    NaN/Inf → 'NA' (standard missing value for tabular formats).
    Guards against both NaN and Inf which can arise from Fisher strand
    bias when ALT total ≤ 1 (OR undefined, see issue #19).
    """
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return "NA"
    return f"{v:.4f}"


def _fmt_vcf(v: float) -> str:
    """Format a float for VCF INFO fields.

    NaN/Inf → '.' (VCF 4.2 spec missing value sentinel).
    The VCF spec does not support 'inf' or 'nan' in Float fields —
    downstream tools (bcftools, GATK, VariantAnnotation) will reject them.
    """
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return "."
    return f"{v:.4f}"


def _fmt_sci(v: float) -> str:
    """Format a float in scientific notation for MAF output.

    NaN/Inf → 'NA'. Used for strand bias p-values which need scientific
    notation for very small values (e.g., 2.4000e-01).
    """
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return "NA"
    return f"{v:.4e}"


def _fmt_vcf_sci(v: float) -> str:
    """Format a float in scientific notation for VCF INFO fields.

    NaN/Inf → '.' (VCF 4.2 spec). Used for strand bias p-values.
    """
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return "."
    return f"{v:.4e}"


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
        rescue_mnp: bool = False,
        has_gtf: bool = False,
        command_line: str = "",
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
            mfsd: If True, append 34 mFSD columns (KS test, LLR, mean sizes,
                pairwise comparisons, derived metrics). Controlled by --mfsd flag.
            mode: Counting mode ('dna' or 'rna'). In RNA mode, 5 RNA-specific
                columns are appended instead of mFSD columns.
            command_line: Full CLI command for #command provenance header.
        """
        self.path = path
        self.column_prefix = column_prefix
        self.preserve_barcode = preserve_barcode
        self.show_normalization = show_normalization
        self.mfsd = mfsd
        self.mode = mode
        self.rescue_mnp = rescue_mnp
        self.has_gtf = has_gtf
        self.command_line = command_line
        self.file = open(path, "w")

        # Write provenance comment headers before TSV data (issue #19).
        # These are #-prefixed lines that downstream readers skip via
        # comment_prefix="#" (e.g., Polars read_maf in batch.py).
        from .. import __version__

        self.file.write(f"#gbcms v{__version__}\n")
        if self.command_line:
            self.file.write(f"#command {self.command_line}\n")

        self.writer: csv.DictWriter | None = None
        self._headers_written = False
        logger.debug(
            "MafWriter initialized: path=%s, column_prefix='%s', "
            "preserve_barcode=%s, show_normalization=%s, mfsd=%s, mode=%s, rescue_mnp=%s, has_gtf=%s",
            path,
            column_prefix,
            preserve_barcode,
            show_normalization,
            mfsd,
            mode,
            rescue_mnp,
            has_gtf,
        )

    def _gbcms_column_names(self) -> list[str]:
        """
        Build the list of gbcms-generated count column names with the configured prefix.

        Column order (v5.1):
          1. Status & diagnostic flags
          2. Core read counts + ALT decomposition (any_alt, partial_alt, n_count)
          3. Read-level strand counts → derived strand bias
          4. Core fragment counts
          5. Fragment-level strand counts → derived fragment strand bias
          6. Optional: mFSD, RNA, normalization

        Design rationale:
          - any_alt/partial_alt/n_count immediately after alt_count for discoverability
          - Strand counts before strand bias (source data before derived metric)
          - Read and fragment layers fully separated (no interleaving)

        Returns:
            Ordered list of gbcms column names.
        """
        p = self.column_prefix

        # ── 1. Status & diagnostic flags ──────────────────────────────────────
        cols = [
            "gbcms_status",
            "gbcms_diagnostic",
        ]
        # gbcms_rescue column is only present when --rescue-mnp is enabled (design §5)
        if self.rescue_mnp:
            cols.append("gbcms_rescue")

        cols.extend(
            [
                # ── 2. Core read counts + ALT decomposition ──────────────────
                f"{p}ref_count",
                f"{p}alt_count",
                # Decomposed ALT counting: any_alt = alt_count + partial_alt
                # (see types.rs invariant). Placed here for immediate visibility
                # alongside alt_count — critical for Phase 3 INDEL diagnostics.
                f"{p}any_alt",
                f"{p}partial_alt",
                # N-base diagnostic: reads with N at discriminating position (duplex masking QC)
                f"{p}n_count",
                f"{p}total_count",
                f"{p}vaf",
                # ── 3. Read-level strand counts → strand bias ─────────────────
                # Source counts first, then the derived Fisher test statistics.
                f"{p}ref_count_forward",
                f"{p}ref_count_reverse",
                f"{p}alt_count_forward",
                f"{p}alt_count_reverse",
                # Strand bias (unprefixed — always unique)
                "strand_bias_p_value",
                "strand_bias_odds_ratio",
                # ── 4. Core fragment counts ───────────────────────────────────
                f"{p}ref_count_fragment",
                f"{p}alt_count_fragment",
                f"{p}total_count_fragment",
                f"{p}vaf_fragment",
                # ── 5. Fragment strand counts → fragment strand bias ──────────
                f"{p}ref_count_fragment_forward",
                f"{p}ref_count_fragment_reverse",
                f"{p}alt_count_fragment_forward",
                f"{p}alt_count_fragment_reverse",
                "fragment_strand_bias_p_value",
                "fragment_strand_bias_odds_ratio",
            ]
        )
        if self.mfsd:
            # ── mFSD: Mutant Fragment Size Distribution (40 columns) ──────────
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
                "mfsd_qval_alt_ref",
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
                # Sub-nucleosomal / mono-nucleosomal fractions
                # for CH-vs-ctDNA differentiation (computed in Rust)
                "mfsd_sub_nuc_ref_frac",
                "mfsd_sub_nuc_alt_frac",
                "mfsd_sub_nuc_enrichment",
                "mfsd_mono_nuc_ref_frac",
                "mfsd_mono_nuc_alt_frac",
                # CH gene flag (computed in Python from Hugo_Symbol)
                "mfsd_ch_flag",
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
            # GTF annotation columns — only present when --gtf is provided
            if self.has_gtf:
                cols += [
                    "exon_boundary_dist",
                    "transcript_read_counts",
                    "transcript_fragment_counts",
                    # P4c: ASJD columns
                    "asjd_flag",
                    "asjd_pval",
                    "asjd_qval",
                    "asjd_ref_junction",
                    "asjd_alt_junction",
                    "asjd_ref_motif",
                    "asjd_alt_motif",
                    "asjd_ref_known",
                    "asjd_alt_known",
                    "asjd_n_ref_junc",
                    "asjd_n_alt_junc",
                    "asjd_n_ref_total",
                    "asjd_n_alt_total",
                    "asjd_diagnostic",
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

    def _populate_gbcms_counts(
        self,
        counts: Any,
        hugo_symbol: str = "",
    ) -> dict[str, str]:
        """
        Build the gbcms count columns dictionary with the configured prefix.

        Calculates VAF values and formats all required count data as strings.
        mFSD derived metrics are only computed when self.mfsd is True — skipping
        the NaN arithmetic and attribute accesses when mFSD is disabled.

        Args:
            counts: BaseCounts object from the Rust engine.
            hugo_symbol: Gene symbol for CH gene flagging (from MAF Hugo_Symbol column).

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
            # Strand bias (unprefixed) — use _fmt/_fmt_sci guards for NaN/Inf (#19)
            "strand_bias_p_value": _fmt_sci(counts.sb_pval),
            "strand_bias_odds_ratio": _fmt(counts.sb_or),
            "fragment_strand_bias_p_value": _fmt_sci(counts.fsb_pval),
            "fragment_strand_bias_odds_ratio": _fmt(counts.fsb_or),
            # Strand counts
            f"{p}ref_count_forward": str(counts.rd_fwd),
            f"{p}ref_count_reverse": str(counts.rd_rev),
            f"{p}alt_count_forward": str(counts.ad_fwd),
            f"{p}alt_count_reverse": str(counts.ad_rev),
            f"{p}ref_count_fragment_forward": str(counts.rdf_fwd),
            f"{p}ref_count_fragment_reverse": str(counts.rdf_rev),
            f"{p}alt_count_fragment_forward": str(counts.adf_fwd),
            f"{p}alt_count_fragment_reverse": str(counts.adf_rev),
            # Decomposed ALT counting
            f"{p}any_alt": str(counts.any_alt),
            f"{p}partial_alt": str(counts.partial_alt),
            # N-base diagnostic
            f"{p}n_count": str(counts.n_count),
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
            # GTF annotation count columns — only populated when --gtf is provided
            if self.has_gtf:
                result.update(
                    {
                        "exon_boundary_dist": (
                            str(counts.exon_boundary_dist)
                            if counts.exon_boundary_dist is not None
                            else ""
                        ),
                        "transcript_read_counts": counts.transcript_read_counts,
                        "transcript_fragment_counts": counts.transcript_fragment_counts,
                        # P4c: ASJD
                        "asjd_flag": str(counts.asjd_flag),
                        "asjd_pval": f"{counts.asjd_pval:.4e}" if counts.asjd_pval < 1.0 else "",
                        "asjd_qval": f"{counts.asjd_qval:.4e}" if counts.asjd_qval < 1.0 else "",
                        "asjd_ref_junction": counts.asjd_ref_junction,
                        "asjd_alt_junction": counts.asjd_alt_junction,
                        "asjd_ref_motif": counts.asjd_ref_motif,
                        "asjd_alt_motif": counts.asjd_alt_motif,
                        "asjd_ref_known": str(counts.asjd_ref_known),
                        "asjd_alt_known": str(counts.asjd_alt_known),
                        "asjd_n_ref_junc": str(counts.asjd_n_ref_junc),
                        "asjd_n_alt_junc": str(counts.asjd_n_alt_junc),
                        "asjd_n_ref_total": str(counts.asjd_n_ref_total),
                        "asjd_n_alt_total": str(counts.asjd_n_alt_total),
                        "asjd_diagnostic": counts.asjd_diagnostic,
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
            # KS test validity: the Rust D-statistic (mfsd_ks_alt_ref) is NaN exactly
            # when a fragment class fell below MIN_FOR_KS (ks_test returns NaN), so
            # derive validity from it rather than re-hardcoding the threshold here. A
            # literal `>= 5` would silently drift from Rust's MIN_FOR_KS if that changed,
            # re-opening the bug where insufficient KS tests are reported as valid.
            mfsd_ks_valid = not math.isnan(counts.mfsd_ks_alt_ref)

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
                    "mfsd_qval_alt_ref": _fmt(counts.mfsd_qval_alt_ref),
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
                    # Sub-nucleosomal / mono-nucleosomal fractions (from Rust)
                    "mfsd_sub_nuc_ref_frac": _fmt(counts.mfsd_sub_nuc_ref_frac),
                    "mfsd_sub_nuc_alt_frac": _fmt(counts.mfsd_sub_nuc_alt_frac),
                    "mfsd_sub_nuc_enrichment": _fmt(counts.mfsd_sub_nuc_enrichment),
                    "mfsd_mono_nuc_ref_frac": _fmt(counts.mfsd_mono_nuc_ref_frac),
                    "mfsd_mono_nuc_alt_frac": _fmt(counts.mfsd_mono_nuc_alt_frac),
                    # CH gene flag (Python-side; True if Hugo_Symbol in CH_GENES)
                    "mfsd_ch_flag": str(hugo_symbol.upper() in CH_GENES if hugo_symbol else False),
                }
            )

        return result

    def write(
        self,
        variant: Variant,
        counts: Any,
        sample_name: str = "TUMOR",
        gbcms_status: str = "PASS",
        gbcms_diagnostic: str = "",
        gbcms_rescue: str = "",
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
            gbcms_status: Normalization/counting status from prepare_variants().
            gbcms_diagnostic: Post-counting diagnostic flags.
            gbcms_rescue: Rescue audit trail.
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
        row["gbcms_status"] = gbcms_status
        row["gbcms_diagnostic"] = gbcms_diagnostic
        # gbcms_rescue only present when --rescue-mnp is enabled (design §5)
        if self.rescue_mnp:
            row["gbcms_rescue"] = gbcms_rescue
        # Extract Hugo_Symbol from MAF metadata for CH gene flagging (empty for VCF input)
        hugo = variant.metadata.get("Hugo_Symbol", "") if variant.metadata else ""
        row.update(self._populate_gbcms_counts(counts, hugo_symbol=hugo))

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
        rescue_mnp: bool = False,
        has_gtf: bool = False,
        command_line: str = "",
        reference_fasta: str = "",
        contigs: list[tuple[str, int]] | None = None,
    ):
        self.path = path
        self.sample_name = sample_name
        self.show_normalization = show_normalization
        self.mfsd = mfsd
        self.mode = mode
        self.rescue_mnp = rescue_mnp
        self.has_gtf = has_gtf
        self.command_line = command_line
        self.reference_fasta = reference_fasta
        self.contigs = contigs or []
        self.file = open(path, "w")
        self._headers_written = False
        logger.debug(
            "VcfWriter initialized: path=%s, sample=%s, show_normalization=%s, "
            "mfsd=%s, mode=%s, rescue_mnp=%s, has_gtf=%s",
            path,
            sample_name,
            show_normalization,
            mfsd,
            mode,
            rescue_mnp,
            has_gtf,
        )

    def _write_header(self):
        """Write VCF header lines.

        Includes provenance (version, command), reference, contig, and FILTER
        headers per VCF 4.2 spec. mFSD ##INFO fields (7 lines) are only
        included when self.mfsd is True.
        """
        from .. import __version__

        headers = [
            "##fileformat=VCFv4.2",
            f"##source=gbcms v{__version__}",
        ]
        # Provenance headers (issue #19)
        if self.command_line:
            headers.append(f"##gbcms_command={self.command_line}")
        if self.reference_fasta:
            headers.append(f"##reference=file://{self.reference_fasta}")
        # Contig headers — recommended by VCF 4.2 spec, required by some tools
        for name, length in self.contigs:
            headers.append(f"##contig=<ID={name},length={length}>")
        # FILTER header — required by VCF 4.2 spec even when only PASS is used
        headers.append('##FILTER=<ID=PASS,Description="All filters passed">')
        # INFO fields
        headers.extend(
            [
                '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
                '##INFO=<ID=GS,Number=1,Type=String,Description="gbcms normalization/counting status">',
                '##INFO=<ID=GD,Number=1,Type=String,Description="gbcms post-counting diagnostic flags">',
            ]
        )
        # GR INFO header only included when --rescue-mnp is enabled (design §5)
        if self.rescue_mnp:
            headers.append(
                '##INFO=<ID=GR,Number=1,Type=String,Description="gbcms rescue audit trail">'
            )
        # INFO field order: DP → status → ALT decomposition → strand bias → mFSD
        headers.extend(
            [
                '##INFO=<ID=AAD,Number=1,Type=Integer,Description="Any ALT Depth: reads with evidence of ALT at >=1 discriminating position (any_alt = ad + partial_alt)">',
                '##INFO=<ID=PAD,Number=1,Type=Integer,Description="Partial ALT Depth: reads matching ALT at some but not all discriminating positions">',
                '##INFO=<ID=NAD,Number=1,Type=Integer,Description="N-base Depth: reads with N base at discriminating position (duplex masking QC)">',
                '##INFO=<ID=SB_PVAL,Number=1,Type=Float,Description="Fisher strand bias p-value">',
                '##INFO=<ID=SB_OR,Number=1,Type=Float,Description="Fisher strand bias odds ratio">',
                '##INFO=<ID=FSB_PVAL,Number=1,Type=Float,Description="Fisher fragment strand bias p-value">',
                '##INFO=<ID=FSB_OR,Number=1,Type=Float,Description="Fisher fragment strand bias odds ratio">',
            ]
        )
        if self.mfsd:
            # mFSD INFO fields (7 primary diagnostics). VCF key = MAF column name uppercased.
            # Only added when --mfsd is set — keeps VCF header minimal for standard runs.
            headers.extend(
                [
                    '##INFO=<ID=MFSD_DELTA_ALT_REF,Number=1,Type=Float,Description="mFSD mean(ALT) − mean(REF) fragment size delta (bp)">',
                    '##INFO=<ID=MFSD_KS_ALT_REF,Number=1,Type=Float,Description="mFSD 2-sample KS D-statistic (ALT vs REF)">',
                    '##INFO=<ID=MFSD_PVAL_ALT_REF,Number=1,Type=Float,Description="mFSD KS p-value (ALT vs REF)">',
                    '##INFO=<ID=MFSD_QVAL_ALT_REF,Number=1,Type=Float,Description="mFSD KS q-value (Benjamini-Hochberg FDR across variants, ALT vs REF; drives TUMOR-LIKE/CH-LIKE)">',
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
            # GTF annotation INFO headers — only when --gtf is provided
            if self.has_gtf:
                headers.extend(
                    [
                        '##INFO=<ID=EBD,Number=1,Type=Integer,Description="Distance (bp) to nearest annotated exon boundary. Missing (.) when no GTF provided.">',
                        '##INFO=<ID=TXRC,Number=1,Type=String,Description="Per-transcript read counts. Format: ENST:AD,RD,DP|ENST:AD,RD,DP. Empty when no GTF or no overlap.">',
                        '##INFO=<ID=TXFC,Number=1,Type=String,Description="Per-transcript fragment counts. Format: ENST:ADF,RDF,DPF|ENST:ADF,RDF,DPF. Empty when no GTF or no overlap.">',
                        # P4c: ASJD INFO headers
                        '##INFO=<ID=ASJD,Number=0,Type=Flag,Description="Allele-specific junction divergence detected (Fisher p<0.05)">',
                        '##INFO=<ID=ASJDP,Number=1,Type=Float,Description="ASJD raw Fisher exact p-value">',
                        '##INFO=<ID=ASJDQ,Number=1,Type=Float,Description="ASJD BH-corrected q-value">',
                        '##INFO=<ID=ASJDRJ,Number=1,Type=String,Description="ASJD dominant REF junction (start-end)">',
                        '##INFO=<ID=ASJDAJ,Number=1,Type=String,Description="ASJD dominant ALT junction (start-end)">',
                        '##INFO=<ID=ASJDRM,Number=1,Type=String,Description="ASJD REF splice motif">',
                        '##INFO=<ID=ASJDAM,Number=1,Type=String,Description="ASJD ALT splice motif">',
                        '##INFO=<ID=ASJDRK,Number=1,Type=Integer,Description="ASJD REF junction known in GTF (1/0)">',
                        '##INFO=<ID=ASJDAK,Number=1,Type=Integer,Description="ASJD ALT junction known in GTF (1/0)">',
                        '##INFO=<ID=ASJDNR,Number=1,Type=Integer,Description="ASJD REF reads on dominant junction">',
                        '##INFO=<ID=ASJDNA,Number=1,Type=Integer,Description="ASJD ALT reads on dominant junction">',
                        '##INFO=<ID=ASJDD,Number=1,Type=String,Description="ASJD diagnostic flags (pipe-separated)">',
                    ]
                )
        # FORMAT fields: VCF 4.2 spec-compliant layout.
        # DP = single int (total depth), AD = Number=R (ref,alt totals),
        # ADF/ADR = strand-by-allele (bcftools convention),
        # FAD/FADF/FADR = fragment-level equivalents.
        headers.extend(
            [
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total read depth">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (ref,alt)">',
                '##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Allelic depths on forward strand (ref_fwd,alt_fwd)">',
                '##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Allelic depths on reverse strand (ref_rev,alt_rev)">',
                '##FORMAT=<ID=VAF,Number=1,Type=Float,Description="Variant allele fraction (read level)">',
                '##FORMAT=<ID=FAD,Number=R,Type=Integer,Description="Fragment allelic depths (ref_frag,alt_frag)">',
                '##FORMAT=<ID=FADF,Number=R,Type=Integer,Description="Fragment depths on forward strand (ref_frag_fwd,alt_frag_fwd)">',
                '##FORMAT=<ID=FADR,Number=R,Type=Integer,Description="Fragment depths on reverse strand (ref_frag_rev,alt_frag_rev)">',
                '##FORMAT=<ID=FAF,Number=1,Type=Float,Description="Variant allele fraction (fragment level)">',
                '##FORMAT=<ID=AAD,Number=1,Type=Integer,Description="Any ALT depth (alt + partial_alt)">',
                '##FORMAT=<ID=PAD,Number=1,Type=Integer,Description="Partial ALT depth">',
                '##FORMAT=<ID=NAD,Number=1,Type=Integer,Description="N-base depth (reads with N at discriminating position)">',
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
        gbcms_status: str = "PASS",
        gbcms_diagnostic: str = "",
        gbcms_rescue: str = "",
        norm_variant: Variant | None = None,
    ):
        if not self._headers_written:
            self._write_header()

        # VCF POS is 1-based
        pos = variant.pos + 1

        # INFO fields (VCF spec: missing values use '.' not 'NA')
        # VCF uses ';' as the INFO field delimiter, so multi-value fields
        # (GS, GD) convert ';' → '|' to avoid parser mis-splitting.
        # GR is handled conditionally below (only when --rescue-mnp).
        gs_vcf = gbcms_status.replace(";", "|")
        gd_vcf = gbcms_diagnostic.replace(";", "|") if gbcms_diagnostic else "."
        info_parts = [
            f"DP={counts.dp}",
            f"GS={gs_vcf}",
            f"GD={gd_vcf}",
        ]
        # GR INFO value only included when --rescue-mnp is enabled (design §5)
        if self.rescue_mnp:
            gr_vcf = gbcms_rescue.replace(";", "|") if gbcms_rescue else "."
            info_parts.append(f"GR={gr_vcf}")
        # INFO order: ALT decomposition (near DP), then strand bias
        info_parts.extend(
            [
                f"AAD={counts.any_alt}",
                f"PAD={counts.partial_alt}",
                f"NAD={counts.n_count}",
                f"SB_PVAL={_fmt_vcf_sci(counts.sb_pval)}",
                f"SB_OR={_fmt_vcf(counts.sb_or)}",
                f"FSB_PVAL={_fmt_vcf_sci(counts.fsb_pval)}",
                f"FSB_OR={_fmt_vcf(counts.fsb_or)}",
            ]
        )
        if self.mfsd:
            # mFSD primary diagnostic INFO fields (7 values).
            # Only populated when --mfsd is set; '.' for NaN per VCF spec.
            info_parts.extend(
                [
                    f"MFSD_DELTA_ALT_REF={_fmt_vcf(counts.mfsd_delta_alt_ref)}",
                    f"MFSD_KS_ALT_REF={_fmt_vcf(counts.mfsd_ks_alt_ref)}",
                    f"MFSD_PVAL_ALT_REF={_fmt_vcf(counts.mfsd_pval_alt_ref)}",
                    f"MFSD_QVAL_ALT_REF={_fmt_vcf(counts.mfsd_qval_alt_ref)}",
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
            # GTF annotation INFO values — only when --gtf is provided
            if self.has_gtf:
                # EBD: exon boundary distance (GTF-informed, '.' when no GTF)
                ebd = counts.exon_boundary_dist
                info_parts.append(f"EBD={ebd if ebd is not None else '.'}")
                # TXRC/TXFC: per-transcript counts (empty → '.', ';' → '|' for VCF)
                txrc = counts.transcript_read_counts
                txfc = counts.transcript_fragment_counts
                info_parts.append(f"TXRC={txrc.replace(';', '|') if txrc else '.'}")
                info_parts.append(f"TXFC={txfc.replace(';', '|') if txfc else '.'}")
                # P4c: ASJD VCF INFO values
                if counts.asjd_flag:
                    info_parts.append("ASJD")
                if counts.asjd_pval < 1.0:
                    info_parts.append(f"ASJDP={counts.asjd_pval:.4e}")
                    info_parts.append(f"ASJDQ={counts.asjd_qval:.4e}")
                if counts.asjd_ref_junction:
                    info_parts.append(f"ASJDRJ={counts.asjd_ref_junction}")
                if counts.asjd_alt_junction:
                    info_parts.append(f"ASJDAJ={counts.asjd_alt_junction}")
                if counts.asjd_ref_motif:
                    info_parts.append(f"ASJDRM={counts.asjd_ref_motif}")
                if counts.asjd_alt_motif:
                    info_parts.append(f"ASJDAM={counts.asjd_alt_motif}")
                info_parts.append(f"ASJDRK={int(counts.asjd_ref_known)}")
                info_parts.append(f"ASJDAK={int(counts.asjd_alt_known)}")
                if counts.asjd_n_ref_junc > 0 or counts.asjd_n_alt_junc > 0:
                    info_parts.append(f"ASJDNR={counts.asjd_n_ref_junc}")
                    info_parts.append(f"ASJDNA={counts.asjd_n_alt_junc}")
                if counts.asjd_diagnostic:
                    info_parts.append(f"ASJDD={counts.asjd_diagnostic.replace(';', '|')}")
        info = ";".join(info_parts)

        # FORMAT fields — VCF 4.2 spec-compliant layout (v5.1).
        # DP = single int (total depth), AD = ref,alt (Number=R),
        # ADF/ADR = strand-by-allele (bcftools convention),
        # FAD/FADF/FADR = fragment-level equivalents.
        gt = "0/1" if counts.ad > 0 else "0/0"
        dp = str(counts.dp)  # single int (VCF spec)
        ad = f"{counts.rd},{counts.ad}"  # Number=R: ref,alt
        adf = f"{counts.rd_fwd},{counts.ad_fwd}"  # fwd strand: ref_fwd,alt_fwd
        adr = f"{counts.rd_rev},{counts.ad_rev}"  # rev strand: ref_rev,alt_rev
        fad = f"{counts.rdf},{counts.adf}"  # fragment: ref_frag,alt_frag
        fadf = f"{counts.rdf_fwd},{counts.adf_fwd}"  # frag fwd: ref_frag_fwd,alt_frag_fwd
        fadr = f"{counts.rdf_rev},{counts.adf_rev}"  # frag rev: ref_frag_rev,alt_frag_rev

        total_reads = counts.rd + counts.ad
        vaf = counts.ad / total_reads if total_reads > 0 else 0.0

        total_frags = counts.rdf + counts.adf
        faf = counts.adf / total_frags if total_frags > 0 else 0.0

        if self.mode == "rna":
            format_str = "GT:DP:AD:ADF:ADR:VAF:FAD:FADF:FADR:FAF:AAD:PAD:NAD:SEN:ANT:ASEN:SPL"
            sample_data = (
                f"{gt}:{dp}:{ad}:{adf}:{adr}:{vaf:.4f}"
                f":{fad}:{fadf}:{fadr}:{faf:.4f}"
                f":{counts.any_alt}:{counts.partial_alt}:{counts.n_count}"
                f":{counts.sense_depth}:{counts.antisense_depth}"
                f":{counts.sense_strand_alt_count}:{counts.splice_spanning_count}"
            )
        else:
            format_str = "GT:DP:AD:ADF:ADR:VAF:FAD:FADF:FADR:FAF:AAD:PAD:NAD"
            sample_data = (
                f"{gt}:{dp}:{ad}:{adf}:{adr}:{vaf:.4f}"
                f":{fad}:{fadf}:{fadr}:{faf:.4f}"
                f":{counts.any_alt}:{counts.partial_alt}:{counts.n_count}"
            )

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
