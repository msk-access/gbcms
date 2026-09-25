"""
Reference base lookups for variant representation.

A MAF row written as VCF needs the reference base maf2vcf prepends (see
:meth:`gbcms.core.kernel.CoordinateKernel.maf_to_vcf`). Counting never reads
the reference from Python — the engine prepares variants against the FASTA
itself — so this is the only Python-side FASTA access.
"""

import logging
from pathlib import Path

import pysam

from ..core.kernel import CoordinateKernel

logger = logging.getLogger(__name__)

__all__ = ["ReferenceBases"]


class ReferenceBases:
    """Single reference bases by 1-based position.

    Contigs are matched with the engine's naming rule
    (:meth:`CoordinateKernel.contig_key`: ``chr1``~``1``, ``chrM``~``MT``). A
    base that cannot be fetched (contig absent from the FASTA, position off the
    contig) is returned as ``N``, which keeps the record valid VCF — a counted
    row at such a locus has already FAILed preparation against the same FASTA.
    Each one is counted; :meth:`close` logs the total.
    """

    def __init__(self, fasta: str | Path):
        self._fasta = pysam.FastaFile(str(fasta))
        self._names = {CoordinateKernel.contig_key(n): n for n in self._fasta.references}
        self.unfetched = 0

    @property
    def contigs(self) -> list[tuple[str, int]]:
        """(name, length) of every FASTA contig, in FASTA order."""
        return list(zip(self._fasta.references, self._fasta.lengths, strict=True))

    def base(self, chrom: str, pos: int) -> str:
        """The upper-case base at 1-based ``pos`` of ``chrom``, or ``N``."""
        name = self._names.get(CoordinateKernel.contig_key(chrom))
        base = ""
        if name is not None and pos >= 1:
            try:
                base = self._fasta.fetch(name, pos - 1, pos)
            except ValueError:
                base = ""
        if base:
            return base.upper()
        self.unfetched += 1
        if self.unfetched <= 5:
            logger.warning(
                "No reference base at %s:%d for a MAF-to-VCF anchor — written as N",
                chrom,
                pos,
            )
        return "N"

    def close(self) -> None:
        if self.unfetched:
            logger.warning(
                "%d MAF-to-VCF anchor base(s) could not be fetched from the reference "
                "and were written as N",
                self.unfetched,
            )
        self._fasta.close()
