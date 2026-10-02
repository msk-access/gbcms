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
- **Soft clips:** for a pure indel, reading stops at a clip (the rule reads aligned
  bases only); for other variants clipped bases are read where they sit.
- **Fragment end:** a read ends at its fragment end. Bases past it (read-through
  into adapter, when the insert is shorter than the read) are not read, and the
  read's extent stops at its last aligned base inside the fragment.

The tract is the census's own slide of the indel along the reference, not prep's
`shift_region`, so a prep error cannot move the census with it.

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
    indel: bool
    contig: str = field(repr=False)
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


def tract(pos: int, ref: str, alt: str, contig: str) -> tuple[int, int]:
    """A pure indel's tract `[lo, hi)`, found by sliding it along `contig` (no use
    of prep): a deletion's removed bases over every equivalent start, an
    insertion's junctions `lo..=hi` (empty, `lo == hi`, in unique sequence)."""
    ref, alt, contig = ref.upper(), alt.upper(), contig.upper()
    if len(alt) > len(ref):
        bases, j = alt[len(ref) :], pos + len(ref)  # inserted before reference position j
        lo, hi, x = j, j, bases
        while hi < len(contig) and contig[hi] == x[0]:
            x, hi = x[1:] + x[0], hi + 1
        x = bases
        while lo > 0 and contig[lo - 1] == x[-1]:
            x, lo = x[-1] + x[:-1], lo - 1
        return lo, hi
    n, start = len(ref) - len(alt), pos + len(alt)  # deleted [start, start + n)
    lo, hi = start, start
    while hi + n < len(contig) and contig[hi] == contig[hi + n]:
        hi += 1
    while lo > 0 and contig[lo - 1] == contig[lo - 1 + n]:
        lo -= 1
    return lo, hi + n


def window_for(variant, contig: str) -> Window:
    """The census window for a `Variant` on `contig` (the reference sequence its
    positions index). A pure indel's anchors are the two flank bases of its tract,
    which both alleles share; the tract is the census's own slide."""
    pos, ref, alt = variant.pos, variant.ref_allele.upper(), variant.alt_allele.upper()
    contig = contig.upper()
    if is_pure_indel(ref, alt):
        lo, hi = tract(pos, ref, alt, contig)
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
        is_pure_indel(ref, alt),
        contig,
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


def _alt_reading(stretch: str, ref: str, alt: str) -> bool | None:
    """The ALT side's reading from one anchor: True once a base where the alleles
    differ is read unmasked and fits the ALT, every unmasked base before it fitting
    the ALT too (bases after it are not read); False at the first unmasked base
    that does not fit the ALT; None when the stretch ends undecided."""
    for i, b in enumerate(stretch):
        if i >= min(len(ref), len(alt)):
            return None
        if b == "N":
            continue
        if b != alt[i]:
            return False
        if ref[i] != alt[i]:
            return True
    return None


def _one_side(stretch: str, ref: str, alt: str, margin: bool) -> tuple[bool, bool]:
    """(fits ALT, fits REF) for a stretch read from one anchor, against the same
    length of each haplotype read from that anchor. For a pure indel, ALT is settled
    at its first deciding base (the bases after it are not read), and REF needs one
    base past the first difference."""
    k = len(stretch)
    fits_alt, fits_ref = _fits(stretch, alt[:k]), _fits(stretch, ref[:k])
    if margin:
        reading = _alt_reading(stretch, ref, alt)
        if reading is not None:
            fits_alt = reading
            if reading:
                return True, False
    if margin and fits_ref and not fits_alt:
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


def _clips(read) -> tuple[int, int]:
    """Query offsets bounding the aligned bases: soft clips at either end lie
    outside `[lead, trail)`."""
    ops = read.cigartuples
    lead = sum(n for op, n in _leading(ops) if op == 4)
    trail = len(read.query_sequence) - sum(n for op, n in _leading(ops[::-1]) if op == 4)
    return lead, trail


def _leading(ops):
    for op, n in ops:
        if op not in (4, 5):
            return
        yield op, n


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


def fragment_span(read) -> tuple[int, int] | None:
    """Query offsets `[lo, hi)` inside the read's fragment, when the fragment is
    well defined: paired, mate mapped on the same contig, the pair facing inward.
    TLEN runs from the forward read's 5' end to the reverse read's 5' end,
    positive on the forward read (an outward pair defines no fragment), so a
    forward read keeps its bases up to the last aligned one
    before the fragment end and a reverse read from the first aligned one at or
    after the fragment start (an insertion beside the boundary lies outside). Past
    an aligned end the boundary falls in the soft clip, base for base."""
    if (
        not read.is_paired
        or read.mate_is_unmapped
        or read.template_length == 0
        or read.next_reference_id != read.reference_id
        or read.is_reverse == read.mate_is_reverse
        or (read.template_length > 0) == read.is_reverse
    ):
        return None
    aligned = read.get_aligned_pairs(matches_only=True)
    if not aligned:
        return None
    n, tlen = len(read.query_sequence), abs(read.template_length)
    start, end = read.reference_start, read.reference_end
    first_q, after_q = aligned[0][0], aligned[-1][0] + 1
    if not read.is_reverse:
        frag_end = start + tlen
        if frag_end >= end:
            return 0, min(after_q + frag_end - end, n)
        return 0, max((q for q, r in aligned if r < frag_end), default=-1) + 1
    frag_start = end - tlen
    if frag_start <= start:
        return max(first_q - (start - frag_start), 0), n
    return min((q for q, r in aligned if r >= frag_start), default=after_q), n


def _molecule(read) -> tuple[int, int, dict[int, int]]:
    """The read's query offsets inside its fragment `[lo, hi)` and its aligned
    pairs there, reference to query."""
    lo, hi = fragment_span(read) or (0, len(read.query_sequence))
    pairs = {r: q for q, r in read.get_aligned_pairs(matches_only=True) if lo <= q < hi}
    return lo, hi, pairs


def judge(read, window: Window, min_baseq: int = 20) -> Verdict:
    seq = read.query_sequence
    if not seq:
        return Verdict.NOT_ANCHORED
    quals = read.query_qualities
    bases = "".join(
        "N" if quals is not None and quals[i] < min_baseq else b for i, b in enumerate(seq.upper())
    )
    frag_lo, frag_hi, pairs = _molecule(read)
    splices = _splices(read)
    stops = [q for q, _, _ in splices]
    ql, qr = pairs.get(window.left), pairs.get(window.right)
    if window.indel:
        # A read that starts (or ends) on a flank anchors on it only when it reads that
        # base: masked it fits anything, and for a deletion sliding through a repeat the
        # flank is what tells the read from REF placed further along the run.

        def reads_flank(q: int, pos: int) -> bool:
            return bases[q] != "N" and pos < len(window.contig) and bases[q] == window.contig[pos]

        if ql is not None and window.left - 1 not in pairs and not reads_flank(ql, window.left):
            ql = None
        if qr is not None and window.right + 1 not in pairs and not reads_flank(qr, window.right):
            qr = None
    n = max(len(window.ref), len(window.alt))
    # A pure indel's rules do not read soft-clipped bases (an aligner's clip of a
    # pure indel's carrier is the clip-borne-carrier work, not this rule's).
    lead, trail = _clips(read) if window.indel else (0, len(bases))
    lead, trail = max(lead, frag_lo), min(trail, frag_hi)

    def from_left():
        end = min([q for q in stops if q > ql] + [trail])
        reach = len(window.ref_right) if window.indel else n  # through the far flank and on
        stretch = bases[ql + 1 : min(ql + 1 + reach, end)]
        return _one_side(stretch, window.ref_right, window.alt_right, window.indel)

    def from_right():
        start = max([q for q in stops if q <= qr] + [lead])
        reach = len(window.ref_left) if window.indel else n
        stretch = bases[max(qr - reach, start) : qr]
        return _one_side(stretch[::-1], window.ref_left[::-1], window.alt_left[::-1], window.indel)

    if ql is not None and qr is not None and not any(ql < q <= qr for q in stops):
        stretch = bases[ql + 1 : qr]
        fits_alt, fits_ref = _fits(stretch, window.alt), _fits(stretch, window.ref)
        # REF needs one base past where the alleles first differ, read from either
        # side: inside the tract for a deletion, beyond the far flank for an insertion.
        if window.indel and fits_ref and not fits_alt:
            if not any(r and not a for a, r in (from_left(), from_right())):
                return Verdict.FITS_BOTH
        return _verdict(fits_alt, fits_ref)
    sides = []
    if ql is not None:
        sides.append(from_left())
    if qr is not None:
        sides.append(from_right())
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
    """Judge every admitted read whose extent inside its fragment overlaps the
    variant's anchor base (the reads whose depth the engine reports), keyed by
    (name, mate)."""
    window = window_for(variant, contig)
    result = Census(window)
    with pysam.AlignmentFile(bam_path) as fh:
        for read in fh.fetch(variant.chrom, variant.pos, variant.pos + 1):
            if read.flag & skip_flags or read.mapping_quality < min_mapq or not read.query_sequence:
                continue
            _, _, pairs = _molecule(read)
            if not pairs or not (min(pairs) <= variant.pos <= max(pairs)):
                continue
            mate = 2 if read.is_read2 else 1
            result.verdicts[(read.query_name, mate)] = judge(read, window, min_baseq)
    return result


def assert_matches(counts, result: Census, what: str = "", depth: bool = True) -> None:
    """Engine REF and ALT equal the census's; everything else is neither. With
    `depth`, DP equals the census's read count (unpaired DNA reads over the anchor,
    none admitted by clipped bases or spliced over it)."""
    assert (counts.rd, counts.ad) == (result.rd, result.ad), (
        f"{what}: engine rd/ad {counts.rd}/{counts.ad}, census {result.rd}/{result.ad} "
        f"({dict(result.counts)})"
    )
    if depth:
        assert counts.dp == len(
            result.verdicts
        ), f"{what}: engine dp {counts.dp}, census reads {len(result.verdicts)}"
    assert counts.dp >= counts.rd + counts.ad
    assert counts.dpf >= counts.rdf + counts.adf
    assert counts.rd == counts.rd_fwd + counts.rd_rev
    assert counts.ad == counts.ad_fwd + counts.ad_rev
