"""
Coordinate Kernel: The source of truth for genomic coordinate systems and
variant representation.

Handles conversion between:
- VCF (1-based)
- MAF (1-based)
- Internal (0-based, half-open [start, end))

Ensures consistent representation of variants:
- SNPs: 0-based index of the base.
- Insertions: 0-based index of the ANCHOR base (preceding the insertion).
- Deletions: 0-based index of the ANCHOR base (preceding the deletion).

VCF <-> MAF conversion follows vcf2maf / maf2vcf (github.com/mskcc/vcf2maf):
:meth:`CoordinateKernel.vcf_to_maf` and :meth:`CoordinateKernel.maf_to_vcf`
are the one place either direction is decided, for the writers and for
``gbcms convert``. Type labels come from the alleles
(:meth:`CoordinateKernel.allele_type`); the engine applies the same rule to
its prepared variants.

Note: This module is pure Python and performs coordinate transformations
only. Allele counting is dispatched via :mod:`gbcms.pipeline` →
``gbcms_rs.count_bam_*`` (Rust FFI). The GIL is released before rayon
parallel iteration in the Rust layer (see ``py.allow_threads``).
"""

from collections.abc import Callable

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

        # VCF POS is the first REF base: the anchor for an indel, the first
        # substituted base otherwise. 1-based 10 -> 0-based 9.
        internal_pos = pos - 1
        vtype = CoordinateKernel.allele_type(ref, alt)

        return Variant(
            chrom=norm_chrom,
            pos=internal_pos,
            ref=ref,
            alt=alt,
            variant_type=vtype,
            original_id=original_id,
            original_chrom=chrom,
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

        else:
            # Sequence alleles on both sides are used as written (the engine
            # anchors only '-' alleles), so they take the VCF-style label.
            vtype = CoordinateKernel.allele_type(ref, alt)
            internal_pos = start_pos - 1

        return Variant(
            chrom=norm_chrom,
            pos=internal_pos,
            ref=ref,
            alt=alt,
            variant_type=vtype,
            original_chrom=chrom,
        )

    @staticmethod
    def allele_type(ref: str, alt: str) -> VariantType:
        """Type of a VCF-style allele pair (bases compared case-insensitively).

        INSERTION / DELETION only when the one-base allele is the other's first
        base — the shared anchor — so the rest of the longer allele is exactly
        what was inserted or deleted. Anything else of unequal length (a delins,
        with or without a shared first base) and every multi-base substitution
        is COMPLEX. The engine's ``variant_type_for`` (rust/src/normalize) is the
        same rule for prepared variants; the counting engine itself dispatches
        on the alleles, never on this label.
        """
        r, a = ref.upper(), alt.upper()
        if len(r) == 1 and len(a) == 1:
            return VariantType.SNP
        if len(r) == 1 and len(a) > 1 and a[0] == r[0]:
            return VariantType.INSERTION
        if len(a) == 1 and len(r) > 1 and r[0] == a[0]:
            return VariantType.DELETION
        return VariantType.COMPLEX

    @staticmethod
    def vcf_to_maf(pos: int, ref: str, alt: str) -> dict[str, str]:
        """MAF coordinates for a VCF record, exactly as vcf2maf writes them.

        The leading bases REF and ALT share are trimmed (never trailing ones),
        advancing the position; an allele trimmed to nothing becomes ``-``.
        Equal trimmed lengths are SNP / DNP / TNP / ONP by length and span the
        allele; otherwise the type is INS (ALT longer) or DEL. A DEL, or an INS
        that still has REF bases, spans its REF bases; an INS whose REF trimmed
        to ``-`` spans the two bases around the insertion point. The trim is
        decided from the alleles alone, so a delins with no shared first base
        keeps every base (``TTAC>A`` is a 4bp DEL, ``C>TA`` a 1bp INS).

        Args:
            pos: 1-based VCF POS.
            ref: VCF REF allele.
            alt: VCF ALT allele (one allele).

        Returns:
            Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2
            and Variant_Type, as strings.
        """
        ref_len, alt_len = len(ref), len(alt)
        while ref and alt and ref[0].upper() == alt[0].upper() and ref.upper() != alt.upper():
            ref, alt = ref[1:] or "-", alt[1:] or "-"
            ref_len, alt_len, pos = ref_len - 1, alt_len - 1, pos + 1
        if ref_len == alt_len:
            start, end = pos, pos + alt_len - 1
            vtype = {1: "SNP", 2: "DNP", 3: "TNP"}.get(alt_len, "ONP")
        elif ref_len < alt_len:
            start, end = (pos - 1, pos) if ref == "-" else (pos, pos + ref_len - 1)
            vtype = "INS"
        else:
            start, end = pos, pos + ref_len - 1
            vtype = "DEL"
        return {
            "Start_Position": str(start),
            "End_Position": str(end),
            "Reference_Allele": ref,
            "Tumor_Seq_Allele2": alt,
            "Variant_Type": vtype,
        }

    @staticmethod
    def maf_to_vcf(
        start: int, ref: str, alt: str, base_at: Callable[[int], str]
    ) -> tuple[int, str, str]:
        """VCF POS / REF / ALT for a MAF row, exactly as maf2vcf writes them.

        The reference base is prepended to both alleles when an allele is
        ``-``, or when the lengths differ and the first bases differ: the base
        AT Start for a ``-`` insertion (MAF Start is the base before the
        insertion), else the base before Start, which becomes POS. Anything
        else (an SNP or MNP, or unequal alleles that already share their first
        base) is written as-is at Start, and no base is fetched.

        Args:
            start: MAF Start_Position (1-based).
            ref: Reference_Allele (``-`` for an insertion).
            alt: Tumor_Seq_Allele2 (``-`` for a deletion).
            base_at: Returns the reference base at a 1-based position of the
                row's contig.

        Returns:
            (POS, REF, ALT), POS 1-based.
        """
        ref, alt = ("" if ref == "-" else ref), ("" if alt == "-" else alt)
        if ref and alt and (len(ref) == len(alt) or ref[0].upper() == alt[0].upper()):
            return start, ref, alt
        pos = start if not ref else start - 1
        anchor = base_at(pos)
        return pos, anchor + ref, anchor + alt

    @staticmethod
    def contig_key(chrom: str) -> str:
        """Naming-independent key for comparing contig names across sources.

        Mirrors the engine's ``normalize_contig`` (rust/src/shared/contig.rs):
        strips a ``chr`` prefix in any case and folds the mitochondrial aliases
        ``M`` / ``MT`` to ``MT``. Used to pair names written differently by the
        variant file, reference and BAMs (``chrM`` vs ``MT``); it is never a
        name written to output.
        """
        bare = chrom[3:] if chrom[:3].lower() == "chr" else chrom
        return "MT" if bare.upper() in ("M", "MT") else bare

    @staticmethod
    def normalize_chromosome(chrom: str) -> str:
        """
        Normalize chromosome name (remove 'chr' prefix).
        """
        if chrom.lower().startswith("chr"):
            return chrom[3:]
        return chrom
