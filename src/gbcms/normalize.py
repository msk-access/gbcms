"""
Standalone Variant Normalization.

This module provides a standalone ``normalize`` workflow that reads variants
from a VCF or MAF file, applies MAF anchor resolution, REF validation,
and bcftools-style left-alignment via the Rust ``prepare_variants()``
function, and writes the results to a TSV file showing both original and
normalized coordinates.

Usage via CLI::

    gbcms normalize -v variants.maf -f ref.fa -o normalized.tsv
"""

import csv
import logging
from pathlib import Path

from gbcms import _rs as gbcms_rs

from .io.input import MafReader, VcfReader

logger = logging.getLogger(__name__)

__all__ = ["normalize_variants"]


def normalize_variants(
    variant_file: Path,
    reference: Path,
    output: Path,
    threads: int = 1,
    context_padding: int = 5,
) -> None:
    """
    Normalize variants and write results to a TSV file.

    Args:
        variant_file: Path to VCF or MAF file.
        reference: Path to reference FASTA.
        output: Output TSV path.
        threads: Number of threads for parallel processing.
        context_padding: Flanking bases for ref_context.
    """
    # 1. Load raw variants
    suffix = variant_file.suffix.lower()
    is_maf = suffix == ".maf"

    reader: VcfReader | MafReader
    if suffix in [".vcf", ".gz"]:
        reader = VcfReader(variant_file)
    elif is_maf:
        reader = MafReader(variant_file)
    else:
        raise ValueError(f"Unsupported variant file format: {suffix}")

    variants = list(reader)
    if hasattr(reader, "close"):
        reader.close()

    logger.info("Loaded %d variants from %s", len(variants), variant_file)

    if not variants:
        logger.error("No variants found. Exiting.")
        return

    # 2. Prepare via Rust
    rs_input = [
        gbcms_rs.Variant(v.chrom, v.pos, v.ref, v.alt, v.variant_type.value) for v in variants
    ]
    prepared = gbcms_rs.prepare_variants(
        rs_input,
        str(reference),
        context_padding,
        is_maf,
        threads,
        True,  # adaptive_context — always enabled for standalone normalization
    )

    # 3. Write TSV
    fieldnames = [
        "chrom",
        "original_pos",
        "original_ref",
        "original_alt",
        "norm_pos",
        "norm_ref",
        "norm_alt",
        "variant_type",
        "validation_status",
        "was_normalized",
    ]

    with open(output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for pv in prepared:
            writer.writerow(
                {
                    "chrom": pv.variant.chrom,
                    "original_pos": pv.original_pos + 1,  # 1-based for display
                    "original_ref": pv.original_ref,
                    "original_alt": pv.original_alt,
                    "norm_pos": pv.variant.pos + 1,  # 1-based for display
                    "norm_ref": pv.variant.ref_allele,
                    "norm_alt": pv.variant.alt_allele,
                    "variant_type": pv.variant.variant_type,
                    "validation_status": pv.validation_status,
                    "was_normalized": pv.was_normalized,
                }
            )

    # Log summary
    n_pass = sum(1 for p in prepared if p.validation_status.startswith("PASS"))
    n_norm = sum(1 for p in prepared if p.was_normalized)
    logger.info(
        "Normalization complete: %d variants, %d PASS, %d left-aligned → %s",
        len(prepared),
        n_pass,
        n_norm,
        output,
    )
