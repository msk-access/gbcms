"""
Standalone VCF <-> MAF conversion.

Writes the representation gbcms's own writers use, without counting:

- VCF input -> MAF: each ALT allele as vcf2maf writes it (leading bases REF and
  ALT share are trimmed; Start/End and Variant_Type follow the trimmed
  alleles), in the columns a VCF-input ``gbcms dna --format maf`` row has
  before its count columns, including the record itself (``vcf_pos`` /
  ``vcf_ref`` / ``vcf_alt``). ALT alleles that cannot be counted are skipped
  the same way (``VcfReader``).
- MAF input -> VCF: each row as maf2vcf writes it (the reference base is
  prepended to ``-`` alleles, and to unequal alleles with different first
  bases), in input order.

Representation only: REF is not checked against the reference and nothing is
left-aligned (``gbcms normalize`` does both).

Usage via CLI::

    gbcms convert -v variants.vcf -o variants.maf
    gbcms convert -v variants.maf -f ref.fa -o variants.vcf
"""

import csv
import logging
from functools import partial
from pathlib import Path

from . import __version__
from .core.kernel import CoordinateKernel
from .io.output import MafWriter, declared_contigs, vcf_contig_lines
from .io.reference import ReferenceBases
from .pipeline import read_variant_file

logger = logging.getLogger(__name__)

__all__ = ["vcf_to_maf_file", "maf_to_vcf_file"]


def vcf_to_maf_file(variant_file: Path, output: Path, command_line: str = "") -> int:
    """Write every countable ALT allele of a VCF as a MAF row; returns the row count."""
    variants = read_variant_file(variant_file)
    with open(output, "w", newline="") as fh:
        fh.write(f"#gbcms v{__version__}\n")
        if command_line:
            fh.write(f"#command {command_line}\n")
        writer = csv.DictWriter(fh, fieldnames=MafWriter.vcf_input_headers(), delimiter="\t")
        writer.writeheader()
        for variant in variants:
            row = dict.fromkeys(writer.fieldnames, "")
            row.update(MafWriter.vcf_input_fields(variant))
            writer.writerow(row)
    logger.info("Converted %s to %d MAF rows: %s", variant_file, len(variants), output)
    return len(variants)


def maf_to_vcf_file(
    variant_file: Path, reference: Path, output: Path, command_line: str = ""
) -> int:
    """Write every MAF row as a VCF record; returns the record count."""
    variants = read_variant_file(variant_file)
    bases = ReferenceBases(reference)
    try:
        header = ["##fileformat=VCFv4.2", f"##source=gbcms v{__version__}"]
        if command_line:
            header.append(f"##gbcms_command={command_line}")
        header.append(f"##reference=file://{reference}")
        # Contigs in the MAF's own naming, as gbcms's VCF output declares them.
        header.extend(vcf_contig_lines(declared_contigs(bases.contigs, variants)))
        header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
        with open(output, "w") as fh:
            fh.write("\n".join(header) + "\n")
            for v in variants:
                pos, ref, alt = CoordinateKernel.maf_to_vcf(
                    v.pos + 1, v.ref, v.alt, partial(bases.base, v.chrom)
                )
                fh.write("\t".join([v.output_chrom, str(pos), ".", ref, alt, ".", ".", "."]))
                fh.write("\n")
    finally:
        bases.close()
    logger.info("Converted %s to %d VCF records: %s", variant_file, len(variants), output)
    return len(variants)
