"""Partial-haplotype contract for multi-nucleotide variants (DNP/ONP).

Pins the cis-phasing behavior that distinguishes gbcms from legacy per-position
counters: a multi-nt ALT is counted only when the FULL haplotype co-occurs on
one read. Reads matching the ALT at some-but-not-all changed positions are
"partial" — surfaced transparently via partial_alt/any_alt, never credited to ad.

Three read-level patterns, each verified against real data before being
synthesized here (all reads below are constructed; no patient data):

1. Spurious annotation — no read carries ANY changed position. The annotated
   multimer has zero support; ad, partial_alt, and any_alt must all be 0.
2. Partial haplotype — every ALT-ish read mutates only ONE position of the
   multimer (an over-specified annotation wrapping a real SNV). Per-position
   counters credit these reads to the multimer; gbcms must not (ad=0), while
   still reporting the partial signal (partial_alt > 0).
3. Mixed cis + partial — ad counts exactly the full-haplotype reads.

Engine invariant asserted throughout: any_alt == ad + partial_alt.
"""

from helpers import build_bam, count_both, make_read

from gbcms import _rs as gbcms_rs

# chr1:100 on the synthetic contig; reads are 10bp starting at 98,
# so the variant's first base sits at read index 2.
ONP = gbcms_rs.Variant("chr1", 100, "GAGGG", "AAGGA", "COMPLEX")  # changed: pos 0 and 4
DNP = gbcms_rs.Variant("chr1", 100, "GG", "AA", "COMPLEX")  # both positions changed

ONP_REF = "AAGAGGGAAT"  # GAGGG at index 2-6
ONP_CIS = "AAAAGGATTT"  # AAGGA — full ALT haplotype
ONP_PARTIAL = "AAAAGGGTTT"  # AAGGG — pos 0 mutated only
DNP_REF = "AAGGAATTTT"  # GG at index 2-3
DNP_CIS = "AAAAAATTTT"  # AA — full ALT haplotype
DNP_PARTIAL = "AAGAAATTTT"  # GA — pos 1 mutated only


def _count(bam_path, variant):
    c = count_both(bam_path, [variant], min_mapq=0, min_baseq=0)[0]
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    assert c.any_alt == c.ad + c.partial_alt
    return c


def _reads(seq, n, prefix):
    return [make_read(f"{prefix}{i}", seq, start=98, cigar=((0, 10),)) for i in range(n)]


def test_spurious_multimer_has_no_signal_at_all(tmp_path):
    """All reads are full REF: the annotated DNP is simply absent.

    ad=0 alone is ambiguous (could be partial support); the contract is that
    partial_alt and any_alt are ALSO 0, so downstream can tell "variant absent"
    from "variant mis-specified".
    """
    bam = build_bam(tmp_path, _reads(DNP_REF, 8, "ref"), "spurious.bam")
    c = _count(bam, DNP)
    assert (c.dp, c.rd, c.ad) == (8, 8, 0)
    assert c.partial_alt == 0
    assert c.any_alt == 0


def test_partial_onp_not_credited_to_ad(tmp_path):
    """5-mer ONP where reads mutate only position 0 (real SNV, over-specified
    annotation). Per-position counters credit all 5 reads to the ONP; the
    full-haplotype rule must not, while keeping the partial signal visible."""
    reads = _reads(ONP_REF, 4, "ref") + _reads(ONP_PARTIAL, 5, "partial")
    bam = build_bam(tmp_path, reads, "partial_onp.bam")
    c = _count(bam, ONP)
    assert (c.dp, c.rd, c.ad) == (9, 4, 0)
    assert c.partial_alt == 5
    assert c.any_alt == 5


def test_partial_dnp_not_credited_to_ad(tmp_path):
    """DNP where reads mutate only position 1 — same contract as the ONP case."""
    reads = _reads(DNP_REF, 4, "ref") + _reads(DNP_PARTIAL, 5, "partial")
    bam = build_bam(tmp_path, reads, "partial_dnp.bam")
    c = _count(bam, DNP)
    assert (c.dp, c.rd, c.ad) == (9, 4, 0)
    assert c.partial_alt == 5
    assert c.any_alt == 5


def test_mixed_cis_and_partial_split_correctly(tmp_path):
    """ad counts exactly the full-haplotype reads; partials stay partial."""
    reads = (
        _reads(DNP_REF, 3, "ref") + _reads(DNP_CIS, 2, "cis") + _reads(DNP_PARTIAL, 4, "partial")
    )
    bam = build_bam(tmp_path, reads, "mixed_dnp.bam")
    c = _count(bam, DNP)
    assert (c.dp, c.rd, c.ad) == (9, 3, 2)
    assert c.partial_alt == 4
    assert c.any_alt == 6


def test_cis_onp_positive_control(tmp_path):
    """The counter DOES fire when the full ONP haplotype is genuinely cis —
    guards against the trivial way to pass the partial tests (never counting)."""
    reads = _reads(ONP_REF, 3, "ref") + _reads(ONP_CIS, 6, "cis")
    bam = build_bam(tmp_path, reads, "cis_onp.bam")
    c = _count(bam, ONP)
    assert (c.dp, c.rd, c.ad) == (9, 3, 6)
    assert c.partial_alt == 0
    assert c.any_alt == 6
