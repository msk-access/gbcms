"""A position-aware read census: what each read's own bases say about a variant,
independent of the engine (#171).

The validation behind this cycle's counting fixes judged every read by its own
bases across the variant's whole repeat tract, and found every real
classification defect that way; comparing two engines that share a classifier
found none. This module is that judge, for synthetic BAMs in tests.

The census reads the bases strictly between two anchors.
- **Pure indel:** the anchors are the two flank bases of its tract (its shift
  region), which both alleles share. The tract is read whole, so a carrier of the
  indel written anywhere in it reads the same.
- **Anything else:** the anchors are the bases either side of the substituted
  bases.

How a read is judged:
- **Both anchors aligned:** the stretch between them must be one allele exactly.
- **One anchor aligned:** the read is read from that anchor as far as the longer
  allele reaches, through soft-clipped bases, and compared with the same length
  of each allele.
- **Masking:** a base below `min_baseq`, or N, fits anything.
- **Splices:** reading stops at a splice N, as if the read ended there.

Two of the engine's decided rules are encoded here, because they are policy, not
bases:
- a REF verdict needs one base past the first base where the alleles differ (the
  margin that guards CIGAR-only REF calls against a hidden terminal mismatch);
- an ALT verdict needs only that base, read unmasked.

Engine REF/ALT counts must equal the census REF/ALT counts. Every other verdict is
neither, whether partial or uninformative. The census is exact for clean bases,
with two limits:
- A pure indel's REF call stands on its CIGAR and extent, so a sequencing error or
  a masked base inside the window can make the two differ.
- The census knows two haplotypes. A read carrying a third allele that ends before
  that allele differs from the ALT fits the ALT, while the engine reads its CIGAR.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from enum import Enum

import pysam

SKIP_FLAGS = (
    0x4 | 0x100 | 0x200 | 0x400 | 0x800
)  # unmapped, secondary, QC-fail, duplicate, supplementary


class Verdict(Enum):
    REF = "REF"
    ALT = "ALT"
    FITS_BOTH = "fits both"
    CONTRADICTS_BOTH = "contradicts both"
    NOT_ANCHORED = "not anchored"
    SPLICED = "spliced"


@dataclass(frozen=True)
class Window:
    """Reference positions of the two anchors (0-based), each allele's bases
    strictly between them, and each haplotype read on past one anchor (for a read
    holding only the other): rightwards from the left anchor (`*_right`), and
    leftwards towards the right anchor (`*_left`, in reference order)."""

    left: int
    right: int
    ref: str
    alt: str
    ref_right: str
    alt_right: str
    ref_left: str
    alt_left: str


@dataclass
class Census:
    window: Window
    verdicts: dict[tuple[str, int], Verdict] = field(default_factory=dict)

    @property
    def counts(self) -> Counter:
        return Counter(self.verdicts.values())

    @property
    def rd(self) -> int:
        return self.counts[Verdict.REF]

    @property
    def ad(self) -> int:
        return self.counts[Verdict.ALT]


def is_pure_indel(ref: str, alt: str) -> bool:
    """One allele is the other plus bases after their shared first base."""
    short, long_ = sorted((ref.upper(), alt.upper()), key=len)
    return len(short) != len(long_) and len(short) >= 1 and long_.startswith(short)


def window_for(variant, contig: str) -> Window:
    """The census window for a prepared `Variant` on `contig` (the reference
    sequence its positions index). A pure indel's anchors are the two flank bases
    of its tract, which both alleles share."""
    pos, ref, alt = variant.pos, variant.ref_allele.upper(), variant.alt_allele.upper()
    contig = contig.upper()
    if is_pure_indel(ref, alt):
        lo, hi = variant.shift_region or (pos + 1, pos + len(ref))  # [lo, hi), the tract
        left, right = lo - 1, hi
    else:
        left, right = pos - 1, pos + len(ref)
    left, right = max(left, -1), min(right, len(contig))
    # The ALT haplotype, sliced at the same anchors: reference positions past the
    # REF allele sit `delta` further along it.
    hap, delta = contig[:pos] + alt + contig[pos + len(ref) :], len(alt) - len(ref)
    ext = abs(delta) + 2  # a one-sided reading may run past the far anchor
    lo_ext, hi_ext = max(left + 1 - ext, 0), min(right + ext, len(contig))
    return Window(
        left,
        right,
        contig[left + 1 : right],
        hap[left + 1 : right + delta],
        contig[left + 1 : hi_ext],
        hap[left + 1 : hi_ext + delta],
        contig[lo_ext:right],
        hap[lo_ext : right + delta],
    )


def _fits(stretch: str, allele: str) -> bool:
    return len(stretch) == len(allele) and all(
        s in ("N", a) for s, a in zip(stretch, allele, strict=True)
    )


def _first_difference(a: str, b: str) -> int | None:
    return next((i for i, (x, y) in enumerate(zip(a, b, strict=False)) if x != y), None)


def _one_side(stretch: str, ref: str, alt: str) -> tuple[bool, bool]:
    """(fits ALT, fits REF) for a stretch read from one anchor, against the same
    length of each haplotype read from that anchor. REF needs one base past the
    first difference."""
    k = len(stretch)
    fits_alt, fits_ref = _fits(stretch, alt[:k]), _fits(stretch, ref[:k])
    if fits_ref and not fits_alt:
        d = _first_difference(ref[:k], alt[:k])
        if d is not None and k < d + 2:
            return True, True  # ends on the deciding base: not decisive for REF
    return fits_alt, fits_ref


def _verdict(fits_alt: bool, fits_ref: bool) -> Verdict:
    if fits_alt and fits_ref:
        return Verdict.FITS_BOTH
    if fits_alt:
        return Verdict.ALT
    if fits_ref:
        return Verdict.REF
    return Verdict.CONTRADICTS_BOTH


def _splices(read) -> list[tuple[int, int, int]]:
    """(query offset, reference start, reference end) of each splice N."""
    out, rp, qp = [], read.reference_start, 0
    for op, n in read.cigartuples:
        if op == 3:
            out.append((qp, rp, rp + n))
        if op in (0, 2, 3, 7, 8):
            rp += n
        if op in (0, 1, 4, 7, 8):
            qp += n
    return out


def judge(read, window: Window, min_baseq: int = 20) -> Verdict:
    seq = read.query_sequence
    if not seq:
        return Verdict.NOT_ANCHORED
    quals = read.query_qualities
    bases = "".join(
        "N" if quals is not None and quals[i] < min_baseq else b for i, b in enumerate(seq.upper())
    )
    pairs = {r: q for q, r in read.get_aligned_pairs(matches_only=True)}
    splices = _splices(read)
    stops = [q for q, _, _ in splices]
    ql, qr = pairs.get(window.left), pairs.get(window.right)
    if window.left < 0:  # the window starts at the contig start: the read's first base anchors
        ql = -1 if read.reference_start == 0 else None
    if ql is not None and qr is not None and not any(ql < q <= qr for q in stops):
        stretch = bases[ql + 1 : qr]
        return _verdict(_fits(stretch, window.alt), _fits(stretch, window.ref))
    n = max(len(window.ref), len(window.alt))
    sides = []
    if ql is not None:
        end = min([q for q in stops if q > ql] + [len(bases)])
        stretch = bases[ql + 1 : min(ql + 1 + n, end)]
        sides.append(_one_side(stretch, window.ref_right, window.alt_right))
    if qr is not None:
        start = max([q for q in stops if q <= qr] + [0])
        stretch = bases[max(qr - n, start) : qr]
        sides.append(_one_side(stretch[::-1], window.ref_left[::-1], window.alt_left[::-1]))
    if not sides:
        covered = any(a <= window.left + 1 and b >= window.right for _, a, b in splices)
        return Verdict.SPLICED if covered else Verdict.NOT_ANCHORED
    return _verdict(all(a for a, _ in sides), all(r for _, r in sides))


def census(
    bam_path: str,
    contig: str,
    variant,
    *,
    min_mapq: int = 20,
    min_baseq: int = 20,
    skip_flags: int = SKIP_FLAGS,
) -> Census:
    """Judge every admitted read overlapping the variant's anchor base (the reads
    whose depth the engine reports), keyed by (name, mate)."""
    window = window_for(variant, contig)
    result = Census(window)
    with pysam.AlignmentFile(bam_path) as fh:
        for read in fh.fetch(variant.chrom, variant.pos, variant.pos + 1):
            if read.flag & skip_flags or read.mapping_quality < min_mapq or not read.query_sequence:
                continue
            if not (read.reference_start <= variant.pos < read.reference_end):
                continue
            mate = 2 if read.is_read2 else 1
            result.verdicts[(read.query_name, mate)] = judge(read, window, min_baseq)
    return result


def assert_matches(counts, result: Census, what: str = "") -> None:
    """Engine REF and ALT equal the census's; everything else is neither."""
    assert (counts.rd, counts.ad) == (result.rd, result.ad), (
        f"{what}: engine rd/ad {counts.rd}/{counts.ad}, census {result.rd}/{result.ad} "
        f"({dict(result.counts)})"
    )
    assert counts.dp >= counts.rd + counts.ad
    assert counts.dpf >= counts.rdf + counts.adf
    assert counts.rd == counts.rd_fwd + counts.rd_rev
    assert counts.ad == counts.ad_fwd + counts.ad_rev
