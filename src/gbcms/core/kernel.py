"""
Coordinate Kernel: The source of truth for genomic coordinate systems.

Handles conversion between:
- VCF (1-based)
- MAF (1-based)
- Internal (0-based, half-open [start, end))

Ensures consistent representation of variants:
- SNPs: 0-based index of the base.
- Insertions: 0-based index of the ANCHOR base (preceding the insertion).
- Deletions: 0-based index of the ANCHOR base (preceding the deletion).
"""

from gbcms.models.core import Variant, VariantType

__all__ = ["CoordinateKernel"]


class CoordinateKernel:
    """
    Stateless utility for coordinate transformations and normalization.
    """

    @staticmethod
    def vcf_to_internal(
        chrom: str, pos: int, ref: str, alt: str, original_id: str | None = None
    ) -> Variant:
        """
        Convert VCF coordinates (1-based) to internal normalized Variant.

        Args:
            chrom: Chromosome name
            pos: 1-based position from VCF
            ref: Reference allele
            alt: Alternate allele
            original_id: Optional VCF ID

        Returns:
            Normalized Variant object
        """
        norm_chrom = CoordinateKernel.normalize_chromosome(chrom)

        # Determine variant type and internal position
        if len(ref) == 1 and len(alt) == 1:
            vtype = VariantType.SNP
            # SNP: VCF POS is the base itself.
            # 1-based 10 -> 0-based 9
            internal_pos = pos - 1

        elif len(ref) == 1 and len(alt) > 1:
            vtype = VariantType.INSERTION
            # Insertion: VCF POS is the base BEFORE the insertion (the anchor).
            # VCF: POS=10, REF=A, ALT=AT (Insertion of T after A at 10)
            # Internal: 0-based index of the ANCHOR base.
            # 1-based 10 -> 0-based 9
            internal_pos = pos - 1

        elif len(ref) > 1 and len(alt) == 1:
            vtype = VariantType.DELETION
            # Deletion: VCF POS is the base BEFORE the deletion (the anchor).
            # VCF: POS=10, REF=AT, ALT=A (Deletion of T after A at 10)
            # Internal POS: 0-based index of the ANCHOR base.
            # 1-based 10 -> 0-based 9
            internal_pos = pos - 1

        else:
            vtype = VariantType.COMPLEX
            # Complex: Treat start as 0-based index of first ref base
            internal_pos = pos - 1

        return Variant(
            chrom=norm_chrom,
            pos=internal_pos,
            ref=ref,
            alt=alt,
            variant_type=vtype,
            original_id=original_id,
        )

    @staticmethod
    def maf_to_internal(chrom: str, start_pos: int, end_pos: int, ref: str, alt: str) -> Variant:
        """
        Convert MAF coordinates (1-based inclusive) to internal normalized Variant.

        MAF coordinates are generally 1-based inclusive [start, end].
        """
        norm_chrom = CoordinateKernel.normalize_chromosome(chrom)

        # Handle MAF indels which often use '-'
        if ref == "-" or alt == "-":
            # MAF Insertion: Start_Position is the base BEFORE the insertion (anchor).
            # ref='-', alt='T' → insertion of T after the anchor.
            # Without FASTA lookup, we cannot determine the anchor base,
            # so ref remains '-'. Use --fasta for proper VCF-style normalization.
            if ref == "-":  # Insertion
                vtype = VariantType.INSERTION
                # MAF Start_Position is the anchor base (1-based).
                # Convert to 0-based: internal_pos = Start_Position - 1.
                internal_pos = start_pos - 1
            else:  # Deletion (alt == '-')
                vtype = VariantType.DELETION
                # MAF Start_Position is the FIRST DELETED base (1-based).
                # For VCF-style anchor-based representation, the anchor is
                # at Start_Position - 1. However, without FASTA we cannot
                # fetch the anchor base, so we store the deletion start.
                # WARNING: This produces non-VCF coordinates. Use --fasta
                # for proper anchor-based normalization via MafReader.
                internal_pos = start_pos - 1

        elif len(ref) == len(alt) == 1:
            vtype = VariantType.SNP
            internal_pos = start_pos - 1

        else:
            # MNPs (same-length multi-base substitutions: DNP, TNP, ONP) and
            # true DelIns variants (e.g., EPHA7: GC→T) are both classified as
            # COMPLEX here. The Rust counting engine dispatches by allele
            # length (ref_len == alt_len → check_mnp), not by this type label.
            # See commit 5f2d4f6 for the original COMPLEX introduction.
            vtype = VariantType.COMPLEX
            internal_pos = start_pos - 1

        return Variant(chrom=norm_chrom, pos=internal_pos, ref=ref, alt=alt, variant_type=vtype)

    @staticmethod
    def _gdc_variant_type(variant: Variant) -> str:
        """
        Map internal variant type to GDC MAF TCGAv2 Variant_Type.

        GDC standard defines: SNP, DNP, TNP, ONP, INS, DEL.
        COMPLEX variants are mapped based on allele length comparison:
        - ref > alt → DEL (net deletion)
        - alt > ref → INS (net insertion)
        - equal length: DNP (len=2), TNP (len=3), ONP (len>3)

        Reference: https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format_TCGAv2/

        Args:
            variant: Internal normalized Variant.

        Returns:
            GDC-compliant Variant_Type string.
        """
        if variant.variant_type == VariantType.SNP:
            return "SNP"
        elif variant.variant_type == VariantType.INSERTION:
            return "INS"
        elif variant.variant_type == VariantType.DELETION:
            return "DEL"
        else:  # COMPLEX — infer from allele lengths
            ref_len = len(variant.ref)
            alt_len = len(variant.alt)
            if ref_len > alt_len:
                return "DEL"  # net deletion
            elif alt_len > ref_len:
                return "INS"  # net insertion
            elif ref_len == 2:
                return "DNP"  # di-nucleotide polymorphism
            elif ref_len == 3:
                return "TNP"  # tri-nucleotide polymorphism
            else:
                return "ONP"  # oligo-nucleotide polymorphism

    @staticmethod
    def internal_to_maf(variant: Variant) -> dict[str, str]:
        """
        Convert internal 0-based VCF-style variant to GDC MAF coordinates.

        This is the reverse of maf_to_internal() / vcf_to_internal().
        Produces MAF-compliant fields following GDC MAF TCGAv2 specification:
        - Start_Position / End_Position are 1-based inclusive
        - Deletions: REF = deleted bases (no anchor), ALT = '-'
        - Insertions: REF = '-', ALT = inserted bases (no anchor)
        - SNP/DNP/TNP/ONP: alleles as-is

        GDC validation rules this satisfies:
        - DEL: End - Start + 1 == len(Reference_Allele), len(ref) >= len(alt)
        - INS: End - Start == 1, len(ref) <= len(alt)
        - SNP/DNP/TNP/ONP: alleles don't contain '-'

        Args:
            variant: Internal normalized Variant (0-based, VCF-style alleles).

        Returns:
            Dictionary with MAF coordinate fields:
            Start_Position, End_Position, Reference_Allele,
            Tumor_Seq_Allele2, Variant_Type (all as strings).
        """
        pos_1based = variant.pos + 1
        vtype = CoordinateKernel._gdc_variant_type(variant)

        if variant.variant_type == VariantType.DELETION:
            # VCF: POS=10, REF=AT, ALT=A → strip anchor base
            # MAF: Start=11 (first deleted base), End=11, REF=T, ALT=-
            maf_ref = variant.ref[1:]  # strip anchor base
            return {
                "Start_Position": str(pos_1based + 1),
                "End_Position": str(pos_1based + len(variant.ref) - 1),
                "Reference_Allele": maf_ref,
                "Tumor_Seq_Allele2": "-",
                "Variant_Type": vtype,
            }
        elif variant.variant_type == VariantType.INSERTION:
            # VCF: POS=10, REF=A, ALT=AT → strip anchor base
            # MAF: Start=10 (anchor), End=11, REF=-, ALT=T
            maf_alt = variant.alt[1:]  # strip anchor base
            return {
                "Start_Position": str(pos_1based),
                "End_Position": str(pos_1based + 1),
                "Reference_Allele": "-",
                "Tumor_Seq_Allele2": maf_alt,
                "Variant_Type": vtype,
            }
        else:
            # SNP, COMPLEX (DNP/TNP/ONP/DEL/INS by GDC label)
            # Coordinates span the reference allele, alleles written as-is
            return {
                "Start_Position": str(pos_1based),
                "End_Position": str(pos_1based + len(variant.ref) - 1),
                "Reference_Allele": variant.ref,
                "Tumor_Seq_Allele2": variant.alt,
                "Variant_Type": vtype,
            }

    @staticmethod
    def normalize_chromosome(chrom: str) -> str:
        """
        Normalize chromosome name (remove 'chr' prefix).
        """
        if chrom.lower().startswith("chr"):
            return chrom[3:]
        return chrom
