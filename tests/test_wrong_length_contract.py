"""Target contract for wrong-length pure-indel evidence (issue #91).

A read whose CIGAR proves a pure insertion/deletion of a DIFFERENT length at
the variant anchor is evidence of a *different allele* (or a truncated
representation, for complex insertion sequences) — never a definitive REF or
full ALT call:

  - pure DEL, single wrong-length D op (outside the >=50bp +/-3bp band)
      -> neither + partial_alt
  - pure INS, wrong length, observed insert NOT a high-identity substring of a
    non-low-complexity expected insert -> neither + partial_alt
  - pure INS truncation of a complex expected insert -> ALT (same event)
  - split/multi-op reads whose net matches the expected event -> Phase-3
    (unchanged; regime-2 protection)
  - delins/complex variants -> Phase-3 (unchanged)

Until the fix lands, the tests documenting the bug are xfail(strict=True):
they fail on develop for the reasons stated, and must flip green with the fix.
Every counting assertion runs through count_both (binned<->legacy parity) and
asserts the counting invariants.
"""

import random

import pysam
import pytest
from helpers import build_bam as _repo_build_bam  # noqa: F401  (500bp variant)
from helpers import count_both, make_read

from gbcms._rs import Variant

XFAIL = pytest.mark.xfail(
    strict=True, reason="issue #91: wrong-length pure indel gets a definitive call"
)


# ── local fixtures: contigs longer than helpers' 500bp default ───────────
def _mk_ref(n=600, seed=1, plants=()):
    rng = random.Random(seed)
    ref = [rng.choice("ACGT") for _ in range(n)]
    for pos, motif in plants:
        ref[pos : pos + len(motif)] = list(motif)
    return "".join(ref)


def _bam(tmp_path, ref, reads, name="wl.bam"):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": len(ref)}]}
    p = tmp_path / name
    reads = sorted(reads, key=lambda a: a.reference_start)
    with pysam.AlignmentFile(p, "wb", header=header) as fh:
        for a in reads:
            fh.write(a)
    sp = tmp_path / name.replace(".bam", ".s.bam")
    pysam.sort("-o", str(sp), str(p))
    pysam.index(str(sp))
    return str(sp)


def _del_variant(ref, p0, length, pad=5):
    """VCF-style pure deletion: anchor at p0-1, REF = anchor+deleted, ALT = anchor."""
    anchor = p0 - 1
    return Variant(
        chrom="chr1",
        pos=anchor,
        ref_allele=ref[anchor : p0 + length],
        alt_allele=ref[anchor],
        variant_type="DELETION",
        ref_context=ref[anchor - pad : p0 + length + pad],
        ref_context_start=anchor - pad,
    )


def _ins_variant(ref, anchor, ins_seq, pad=5):
    """VCF-style pure insertion after `anchor` (0-based)."""
    return Variant(
        chrom="chr1",
        pos=anchor,
        ref_allele=ref[anchor],
        alt_allele=ref[anchor] + ins_seq,
        variant_type="INSERTION",
        ref_context=ref[anchor - pad : anchor + 1 + pad],
        ref_context_start=anchor - pad,
    )


def _ref_reads(ref, center, n=8, rl=100):
    return [
        make_read(
            f"ref{i}",
            ref[center - 50 + (i % 5) : center - 50 + (i % 5) + rl],
            center - 50 + (i % 5),
            ((0, rl),),
        )
        for i in range(n)
    ]


def _del_reads(ref, p0, length, n=6, rl=100):
    out = []
    for i in range(n):
        s = p0 - 30 - (i % 5)
        left = p0 - s
        seq = ref[s:p0] + ref[p0 + length : p0 + length + (rl - left)]
        out.append(make_read(f"del{i}", seq, s, ((0, left), (2, length), (0, rl - left))))
    return out


def _split_del_reads(ref, p0, l1, mid, l2, n=6, rl=100):
    """Two D ops separated by `mid` matched bases: net = -(l1+l2)."""
    out = []
    for i in range(n):
        s = p0 - 30 - (i % 5)
        left = p0 - s
        right = rl - left - mid
        seq = (
            ref[s:p0]
            + ref[p0 + l1 : p0 + l1 + mid]
            + ref[p0 + l1 + mid + l2 : p0 + l1 + mid + l2 + right]
        )
        out.append(
            make_read(f"split{i}", seq, s, ((0, left), (2, l1), (0, mid), (2, l2), (0, right)))
        )
    return out


def _ins_reads(ref, anchor, ins_seq, n=6, rl=100):
    out = []
    for i in range(n):
        s = anchor - 30 - (i % 5)
        left = anchor + 1 - s
        L = len(ins_seq)
        seq = ref[s : anchor + 1] + ins_seq + ref[anchor + 1 : anchor + 1 + (rl - left - L)]
        out.append(make_read(f"ins{i}", seq, s, ((0, left), (1, L), (0, rl - left - L))))
    return out


def _count(bam, variant):
    c = count_both(bam, [variant], min_mapq=0, min_baseq=0)[0]
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    assert c.any_alt == c.ad + c.partial_alt
    return c


HOMOPOLY = [(199, "C"), (200, "A" * 10), (210, "G")]  # pinned non-A anchor/flank
GGC = [(199, "T"), (200, "GGC" * 8), (224, "A")]  # pinned non-repeat anchor/flank

# a complex (non-low-complexity) 19bp insert for the containment cases
COMPLEX_INS = "CTTAGTCACCTTCGTGGCA"


# ═════════════════════════ pure deletions ════════════════════════════════
def test_del_exact_control(tmp_path):
    ref = _mk_ref(plants=HOMOPOLY)
    bam = _bam(tmp_path, ref, _ref_reads(ref, 200) + _del_reads(ref, 200, 2))
    c = _count(bam, _del_variant(ref, 200, 2))
    assert (c.rd, c.ad, c.partial_alt) == (8, 6, 0)


@XFAIL
def test_del_homopolymer_1bp_vs_2bp(tmp_path):
    """The clinical headline case: 1bp slippage allele vs annotated 2bp del."""
    ref = _mk_ref(plants=HOMOPOLY)
    bam = _bam(tmp_path, ref, _ref_reads(ref, 200) + _del_reads(ref, 200, 1))
    c = _count(bam, _del_variant(ref, 200, 2))
    assert (c.rd, c.ad, c.partial_alt) == (8, 0, 6)


@XFAIL
def test_del_homopolymer_3bp_vs_2bp(tmp_path):
    ref = _mk_ref(plants=HOMOPOLY)
    bam = _bam(tmp_path, ref, _ref_reads(ref, 200) + _del_reads(ref, 200, 3))
    c = _count(bam, _del_variant(ref, 200, 2))
    assert (c.rd, c.ad, c.partial_alt) == (8, 0, 6)


@XFAIL
def test_del_ggc_tract_6bp_vs_3bp(tmp_path):
    ref = _mk_ref(plants=GGC)
    bam = _bam(tmp_path, ref, _ref_reads(ref, 200) + _del_reads(ref, 200, 6))
    c = _count(bam, _del_variant(ref, 200, 3))
    assert (c.rd, c.ad, c.partial_alt) == (8, 0, 6)


def test_del_unique_context_1bp_vs_2bp(tmp_path):
    """Guard: unique context already resolves to neither+partial today."""
    ref = _mk_ref()
    bam = _bam(tmp_path, ref, _ref_reads(ref, 200) + _del_reads(ref, 200, 1))
    c = _count(bam, _del_variant(ref, 200, 2))
    assert (c.rd, c.ad, c.partial_alt) == (8, 0, 6)


def test_del_large_exact_control(tmp_path):
    ref = _mk_ref(n=800)
    bam = _bam(tmp_path, ref, _ref_reads(ref, 300) + _del_reads(ref, 300, 100))
    c = _count(bam, _del_variant(ref, 300, 100))
    assert (c.rd, c.ad, c.partial_alt) == (8, 6, 0)


@XFAIL
def test_del_large_60_vs_100(tmp_path):
    """Same-anchor D(60) for an annotated 100bp del: real large dels show zero
    length wobble (9 events, 12-539bp), so 60 != 100 is a different event."""
    ref = _mk_ref(n=800)
    bam = _bam(tmp_path, ref, _ref_reads(ref, 300) + _del_reads(ref, 300, 60))
    c = _count(bam, _del_variant(ref, 300, 100))
    assert (c.rd, c.ad, c.partial_alt) == (8, 0, 6)


@XFAIL
def test_del_large_40_vs_100(tmp_path):
    ref = _mk_ref(n=800)
    bam = _bam(tmp_path, ref, _ref_reads(ref, 300) + _del_reads(ref, 300, 40))
    c = _count(bam, _del_variant(ref, 300, 100))
    assert (c.rd, c.ad, c.partial_alt) == (8, 0, 6)


def test_del_large_split_representation_stays_alt(tmp_path):
    """Regime-2 protection: D(50)+2M+D(48) net=-98~-100 with wobble... net
    -(50+48)=-98 is outside the +/-3 band of -100 on purpose? No: keep the
    protected case exact — net -100 via D(50)+2M+D(50) skipping 2 matched
    bases is not expressible; use D(60)+2M+D(40): net -100. Must stay ALT."""
    ref = _mk_ref(n=800)
    bam = _bam(tmp_path, ref, _ref_reads(ref, 300) + _split_del_reads(ref, 300, 60, 2, 40))
    c = _count(bam, _del_variant(ref, 300, 102))
    # deleted span = 60 + 40 with 2 matched bases between; annotate the exact
    # composite: expected REF spans 102 bases of which 100 are deleted. The
    # contract: this read pattern must keep counting as ALT (Phase-3 arbitration),
    # not be demoted by the wrong-length rule.
    assert c.ad == 6


# ═════════════════════════ pure insertions ═══════════════════════════════
def test_ins_exact_control(tmp_path):
    ref = _mk_ref(plants=HOMOPOLY)
    bam = _bam(tmp_path, ref, _ref_reads(ref, 199) + _ins_reads(ref, 199, "AA"))
    c = _count(bam, _ins_variant(ref, 199, "AA"))
    assert (c.rd, c.ad, c.partial_alt) == (8, 6, 0)


@XFAIL
def test_ins_homopolymer_1bp_vs_2bp(tmp_path):
    """+A vs annotated +AA in the homopolymer: distinct slippage alleles.
    The read carrying +A is NOT reference — rd must not absorb it."""
    ref = _mk_ref(plants=HOMOPOLY)
    bam = _bam(tmp_path, ref, _ref_reads(ref, 199) + _ins_reads(ref, 199, "A"))
    c = _count(bam, _ins_variant(ref, 199, "AA"))
    assert (c.rd, c.ad, c.partial_alt) == (8, 0, 6)


def test_ins_unique_context_1bp_vs_2bp(tmp_path):
    ref = _mk_ref()
    bam = _bam(tmp_path, ref, _ref_reads(ref, 200) + _ins_reads(ref, 200, "T"))
    c = _count(bam, _ins_variant(ref, 200, "TT"))
    assert (c.rd, c.ad, c.partial_alt) == (8, 0, 6)


@XFAIL
def test_ins_truncation_of_complex_insert_counts_alt(tmp_path):
    """Regime 3: sequencing errors truncate long inserted sequences, producing
    shorter I ops whose sequence is a prefix of the expected insert. Sign-out
    counts these as the same event; so must we (containment rule)."""
    ref = _mk_ref()
    bam = _bam(tmp_path, ref, _ref_reads(ref, 200) + _ins_reads(ref, 200, COMPLEX_INS[:11]))
    c = _count(bam, _ins_variant(ref, 200, COMPLEX_INS))
    assert (c.rd, c.ad, c.partial_alt) == (8, 6, 0)


@XFAIL
def test_ins_unrelated_insert_stays_partial(tmp_path):
    """Wrong-length insert whose sequence is NOT contained in the expected
    insert: different event -> partial."""
    ref = _mk_ref()
    bam = _bam(tmp_path, ref, _ref_reads(ref, 200) + _ins_reads(ref, 200, "GGGGGGGGGGG"))
    c = _count(bam, _ins_variant(ref, 200, COMPLEX_INS))
    assert (c.rd, c.ad, c.partial_alt) == (8, 0, 6)


# ═════════════════════════ delins guard (regime 2) ═══════════════════════
def test_delins_reads_keep_phase3_routing(tmp_path):
    """Complex AAA>G vs plain 1bp-del reads: the wrong-length rule must NOT
    apply to delins — Phase-3 owns this class (validated on the curated
    internal set). This test pins the routing by asserting the wrong-length
    reads are not force-converted to partial via the pure-indel rule:
    whatever Phase-3 decides, ad stays 0 and the classification is stable
    before and after the fix."""
    ref = _mk_ref(plants=HOMOPOLY)
    bam = _bam(tmp_path, ref, _ref_reads(ref, 200) + _del_reads(ref, 200, 1))
    v = Variant(
        chrom="chr1",
        pos=199,
        ref_allele=ref[199:203],
        alt_allele=ref[199] + "G",
        variant_type="COMPLEX",
        ref_context=ref[194:208],
        ref_context_start=194,
    )
    c = _count(bam, v)
    assert c.ad == 0
