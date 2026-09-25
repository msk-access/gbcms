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
import re
import sys
from collections import Counter
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

# A countable allele: bases only (N included; preparation rejects it visibly).
_SEQUENCE_ALLELE = re.compile(r"[ACGTNacgtn]+")

__all__ = ["VariantReader", "VcfReader", "MafReader"]


class VariantReader:
    """Abstract base class for variant readers."""

    def __iter__(self) -> Iterator[Variant]:
        raise NotImplementedError


class VcfReader(VariantReader):
    """Reads variants from a VCF file, one Variant per ALT allele.

    ALT alleles that name no sequence to count are skipped with a WARNING, not
    yielded: ``*`` (an overlapping deletion), symbolic ``<...>`` alleles,
    breakends, a missing ALT (``.``), and any other allele with a base outside
    A/C/G/T/N. ``N`` is kept so preparation reports it as a visible FAIL row
    (``ALT_CONTAINS_N``). A record whose REF is not a base sequence (pysam
    reports an empty REF as ``.``) is skipped the same way. Each skip is
    counted by reason and the totals are logged once per file.
    """

    def __init__(self, path: Path):
        self.path = path
        self._vcf = pysam.VariantFile(str(path))

    @staticmethod
    def _uncountable(alt: str) -> str | None:
        """Why an ALT allele cannot be counted, or None when it can."""
        if alt == ".":
            return "missing ALT '.'"
        if alt == "*":
            return "overlapping deletion '*'"
        if alt.startswith("<"):
            return "symbolic allele"
        if "[" in alt or "]" in alt or alt.startswith(".") or alt.endswith("."):
            return "breakend"
        if not _SEQUENCE_ALLELE.fullmatch(alt):
            return "non-sequence allele"
        return None

    def __iter__(self) -> Iterator[Variant]:
        skipped: Counter[str] = Counter()

        def skip(record: pysam.VariantRecord, alt: str, reason: str) -> None:
            skipped[reason] += 1
            if sum(skipped.values()) <= 5:
                logger.warning(
                    "Skipped VCF ALT %s:%d %s>%s: %s — not countable",
                    record.chrom,
                    record.pos,
                    record.ref,
                    alt,
                    reason,
                )

        for record in self._vcf:
            # pysam record.pos is the 1-based VCF POS (record.start is 0-based),
            # which is what CoordinateKernel.vcf_to_internal expects.
            if not record.alts:
                skip(record, ".", "missing ALT '.'")
                continue
            ref = record.ref or ""
            for alt in record.alts:
                if not _SEQUENCE_ALLELE.fullmatch(ref):
                    skip(record, alt, "REF not a base sequence")
                    continue
                reason = self._uncountable(alt)
                if reason:
                    skip(record, alt, reason)
                    continue
                yield CoordinateKernel.vcf_to_internal(
                    chrom=record.chrom,
                    pos=record.pos,
                    ref=ref,
                    alt=alt,
                    original_id=record.id,
                )
        if skipped:
            logger.warning(
                "VcfReader: skipped %d ALT alleles that cannot be counted (%s) in %s",
                sum(skipped.values()),
                ", ".join(f"{reason}: {n}" for reason, n in skipped.most_common()),
                self.path,
            )

    def close(self):
        self._vcf.close()


class MafReader(VariantReader):
    """Reads variants from a MAF file.

    Yields raw MAF coordinates as internal ``Variant`` objects using
    ``maf_to_internal()``.  Anchor base resolution, REF validation,
    and left-alignment are performed downstream by the Rust
    ``prepare_variants()`` function. Alleles are read as maf2vcf reads them
    (:meth:`CoordinateKernel.maf_alleles`): ``Tumor_Seq_Allele1`` is the
    variant allele when ``Tumor_Seq_Allele2`` is empty or the reference, and
    the number of rows read that way is logged.

    Args:
        path: Path to the MAF file.
    """

    def __init__(self, path: Path):
        self.path = path

    def __iter__(self) -> Iterator[Variant]:
        skipped = allele1_rows = 0
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
                    allele2 = row["Tumor_Seq_Allele2"]
                    ref, alt = CoordinateKernel.maf_alleles(
                        row["Reference_Allele"], row.get("Tumor_Seq_Allele1") or "", allele2
                    )
                    # Counted when the fallback took Allele1: the result differs
                    # from reading Allele2 alone.
                    if alt != CoordinateKernel.maf_alleles(ref, "", allele2)[1]:
                        allele1_rows += 1

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

        if allele1_rows:
            logger.warning(
                "MafReader: Tumor_Seq_Allele1 is the variant allele for %d row(s) whose "
                "Tumor_Seq_Allele2 is empty or the reference (maf2vcf's reading) in %s",
                allele1_rows,
                self.path,
            )
        if skipped > 5:
            logger.warning("... and %d more malformed MAF rows", skipped - 5)
        if skipped:
            logger.info("MafReader: skipped %d malformed rows out of file %s", skipped, self.path)

    def close(self):
        pass
