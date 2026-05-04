"""
Input Adapters: Handling VCF and MAF inputs.

This module provides classes to read variants from VCF and MAF files,
converting them into the internal normalized representation using CoordinateKernel.

Note:
    MAF indel normalization (anchor base fetch, REF validation, left-alignment)
    is now handled by the Rust ``prepare_variants()`` function.  MafReader
    yields raw MAF coordinates; the pipeline calls ``prepare_variants`` to
    normalise them before counting.
"""

import csv
import logging
import sys
from collections.abc import Iterator
from pathlib import Path

import pysam
from pydantic import ValidationError

from ..core.kernel import CoordinateKernel
from ..models.core import Variant

# Production MAF files (e.g. data_mutations_extended.txt) can contain fields
# exceeding Python's default 128 KB CSV field limit (long COMMENTS columns).
csv.field_size_limit(sys.maxsize)

logger = logging.getLogger(__name__)

__all__ = ["VariantReader", "VcfReader", "MafReader"]


class VariantReader:
    """Abstract base class for variant readers."""

    def __iter__(self) -> Iterator[Variant]:
        raise NotImplementedError


class VcfReader(VariantReader):
    """Reads variants from a VCF file."""

    def __init__(self, path: Path):
        self.path = path
        self._vcf = pysam.VariantFile(str(path))

    def __iter__(self) -> Iterator[Variant]:
        skipped = 0
        for record in self._vcf:
            # pysam record.pos is 0-based. VCF POS is 1-based.
            # CoordinateKernel.vcf_to_internal expects 1-based VCF POS.
            for alt in record.alts or []:
                if not record.ref:
                    skipped += 1
                    if skipped <= 5:
                        logger.warning(
                            "Skipped VCF record with empty REF: %s:%d",
                            record.chrom,
                            record.pos + 1,
                        )
                    continue

                yield CoordinateKernel.vcf_to_internal(
                    chrom=record.chrom,
                    pos=record.pos,
                    ref=record.ref,
                    alt=alt,
                    original_id=record.id,
                )
        if skipped > 5:
            logger.warning("... and %d more VCF records with empty REF", skipped - 5)
        if skipped:
            logger.info("VcfReader: skipped %d records with empty REF allele", skipped)

    def close(self):
        self._vcf.close()


class MafReader(VariantReader):
    """Reads variants from a MAF file.

    Yields raw MAF coordinates as internal ``Variant`` objects using
    ``maf_to_internal()``.  Anchor base resolution, REF validation,
    and left-alignment are performed downstream by the Rust
    ``prepare_variants()`` function.

    Args:
        path: Path to the MAF file.
    """

    def __init__(self, path: Path):
        self.path = path

    def __iter__(self) -> Iterator[Variant]:
        skipped = 0
        with open(self.path) as f:
            # Skip comment lines
            while True:
                pos = f.tell()
                line = f.readline()
                if not line.startswith("#"):
                    f.seek(pos)
                    break

            reader = csv.DictReader(f, delimiter="\t")

            for row_num, row in enumerate(reader, start=1):
                try:
                    chrom = row["Chromosome"]
                    start_pos = int(row["Start_Position"])
                    end_pos = int(row["End_Position"])
                    ref = row["Reference_Allele"]
                    alt = row["Tumor_Seq_Allele2"]

                    yield CoordinateKernel.maf_to_internal(
                        chrom=chrom,
                        start_pos=start_pos,
                        end_pos=end_pos,
                        ref=ref,
                        alt=alt,
                    ).model_copy(update={"metadata": row})

                except (KeyError, ValueError, TypeError, ValidationError) as exc:
                    skipped += 1
                    if skipped <= 5:
                        logger.warning(
                            "Skipped MAF row %d: %s — %s",
                            row_num,
                            type(exc).__name__,
                            exc,
                        )
                    continue

        if skipped > 5:
            logger.warning("... and %d more malformed MAF rows", skipped - 5)
        if skipped:
            logger.info("MafReader: skipped %d malformed rows out of file %s", skipped, self.path)

    def close(self):
        pass
