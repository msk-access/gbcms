"""Complex variants count only reads that carry the whole allele (C1 #141).

A complex variant (delins, a deletion whose anchor also changes, or an MNP read
with an indel) is REF or ALT for a read only when the read's own bases carry
that allele across the whole event, with two reference bases of flank on each
side. Soft-clipped bases count; bases below min BQ are the only tolerance.
Reads that end inside the event, or carry a near-match, are neither (a read
closer to ALT is partial evidence). REF and ALT windows have equal length, so
neither allele is favoured by read placement.

This replaces tolerant scoring (local alignment, likelihoods, an edit-distance
margin, anchor-only REF) that credited reads carrying other alleles.
"""

import random

import pysam
import pytest
from helpers import count_both, make_read

from gbcms import _rs as gbcms_rs

READ = 100
POS = 300  # 0-based start of the event


def _ref(n=900, seed=141):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _files(tmp_path, ref, reads):
    fa = tmp_path / "ref.fa"
    fa.write_text(f">1\n{ref}\n")
    pysam.faidx(str(fa))
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": len(ref)}]}
    raw = tmp_path / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=hdr) as fh:
        for r in reads:
            fh.write(r)
    bam = tmp_path / "s.bam"
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    return str(fa), str(bam)


def _prepared(fa, ref_allele, alt_allele, pos=POS):
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", pos, ref_allele, alt_allele, "COMPLEX")], fa, 5, False, 1, True
    )
    assert pv.gbcms_status == "PASS", pv.gbcms_status_reason
    return pv.variant


def _distinct_alt(ref, n, seed=7):
    """An ALT of length n whose first and last bases differ from the REF's."""
    rng = random.Random(seed)
    while True:
        alt = "".join(rng.choice("ACGT") for _ in range(n))
        if alt[0] != ref[POS] and alt[-1] != ref[POS + 1]:
            return alt


def _invariants(c):
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev


def _alt_read(name, hap, s, ref_len, alt_len, quals=None):
    """An ALT read starting at or before the event, aligned as matches up to it,
    then an I or D for the length change, then matches (mismatches where the ALT
    differs). A read that ends before the event, or inside an insert, is
    aligned as far as it goes."""
    left = POS - s
    d = alt_len - ref_len
    if left >= READ:
        cig = ((0, READ),)
    elif d > 0:
        cig = (
            ((0, left), (1, d), (0, READ - left - d))
            if READ - left - d > 0
            else ((0, left), (1, READ - left))
        )
    elif d < 0:
        cig = ((0, left), (2, -d), (0, READ - left))
    else:
        cig = ((0, READ),)
    return make_read(name, hap[s : s + READ], s, cig, quals=quals)


# ── counting rules ───────────────────────────────────────────────────────


def _delins_case(ref, alt, extra_alt=()):
    hap = ref[:POS] + alt + ref[POS + 2 :]
    reads = [
        make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),))
        for i, s in enumerate(range(240, 250))
    ]
    reads += [_alt_read(f"a{i}", hap, s, 2, len(alt)) for i, s in enumerate(range(240, 250))]
    return hap, reads + list(extra_alt)


def test_reads_ending_inside_the_event_are_not_ref_or_alt(tmp_path):
    ref = _ref()
    alt = _distinct_alt(ref, 3)
    hap, reads = _delins_case(ref, alt)
    for i, s in enumerate(range(POS - READ + 1, POS - READ + 5)):  # last base POS..POS+3
        reads.append(make_read(f"rp{i}", ref[s : s + READ], s, ((0, READ),)))
        reads.append(make_read(f"ap{i}", hap[s : s + READ], s, ((0, READ),)))
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 2], alt)])[0]
    _invariants(c)
    assert (c.rd, c.ad) == (10, 10)


def test_a_read_one_confident_base_off_is_not_alt(tmp_path):
    ref = _ref()
    alt = _distinct_alt(ref, 3)
    off = alt[0] + ("A" if alt[1] != "A" else "C") + alt[2]
    off_hap = ref[:POS] + off + ref[POS + 2 :]
    near = [_alt_read(f"n{i}", off_hap, s, 2, 3) for i, s in enumerate(range(250, 258))]
    hap, reads = _delins_case(ref, alt, near)
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 2], alt)])[0]
    _invariants(c)
    assert c.ad == 10  # the 8 near-matches are not the given allele
    assert c.partial_alt >= 8  # but they are ALT-like evidence


def test_a_low_quality_differing_base_is_masked(tmp_path):
    ref = _ref()
    alt = _distinct_alt(ref, 3)
    off = alt[0] + ("A" if alt[1] != "A" else "C") + alt[2]
    off_hap = ref[:POS] + off + ref[POS + 2 :]
    near = []
    for i, s in enumerate(range(250, 258)):
        quals = [30] * READ
        quals[POS - s + 1] = 5  # the differing base is below min BQ
        near.append(_alt_read(f"n{i}", off_hap, s, 2, 3, quals=quals))
    hap, reads = _delins_case(ref, alt, near)
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 2], alt)])[0]
    assert c.ad == 18


def test_a_soft_clipped_carrier_counts(tmp_path):
    """The aligner reads the event's first base as a mismatch and clips the rest:
    the clipped bases carry the ALT and its flank, so the read is a carrier."""
    ref = _ref()
    alt = _distinct_alt(ref, 3)
    hap, reads = _delins_case(ref, alt)
    for i, s in enumerate(range(236, 240)):
        left = POS - s + 1
        reads.append(make_read(f"c{i}", hap[s : s + READ], s, ((0, left), (4, READ - left))))
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 2], alt)])[0]
    _invariants(c)
    assert c.ad == 14


def test_a_carrier_clipped_before_the_event_is_outside_depth(tmp_path):
    """Aligned bases stop just before the event and the whole ALT is clipped: the
    read does not overlap the variant position, so it is outside DP (depth at the
    position, as in a pileup) and cannot count as REF or ALT either."""
    ref = _ref()
    alt = _distinct_alt(ref, 3)
    hap, reads = _delins_case(ref, alt)
    for i, s in enumerate(range(236, 240)):
        left = POS - s
        reads.append(make_read(f"c{i}", hap[s : s + READ], s, ((0, left), (4, READ - left))))
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 2], alt)])[0]
    _invariants(c)
    assert (c.ad, c.dp) == (10, 20)


def test_one_haplotype_aligned_two_ways_counts_the_same(tmp_path):
    ref = _ref()
    alt = _distinct_alt(ref, 3)
    hap = ref[:POS] + alt + ref[POS + 2 :]
    reads = [_alt_read(f"a{i}", hap, s, 2, 3) for i, s in enumerate(range(240, 250))]
    for i, s in enumerate(range(250, 260)):  # same bases, the I placed after the mismatches
        left = POS + 2 - s
        reads.append(
            make_read(f"b{i}", hap[s : s + READ], s, ((0, left), (1, 1), (0, READ - left - 1)))
        )
    reads += [
        make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),))
        for i, s in enumerate(range(240, 250))
    ]
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 2], alt)])[0]
    assert (c.rd, c.ad) == (10, 20)


# ── no placement bias ────────────────────────────────────────────────────

SHAPES = {
    "deletion-type 27>6": (27, 6),
    "insertion-type 2>9": (2, 9),
    "MNP-like 3>3 with an indel read": (3, 4),
    "long 60>5": (60, 5),
}


@pytest.mark.parametrize("shape", list(SHAPES))
def test_vaf_is_unbiased_by_read_placement(tmp_path, shape):
    """50% VAF, one read per molecule, starts uniform over the locus: the counted
    VAF stays near 50% whatever the event's shape."""
    ref_len, alt_len = SHAPES[shape]
    ref = _ref()
    rng = random.Random(ref_len * 31 + alt_len)
    alt = ""
    while not alt or alt[0] == ref[POS] or alt[-1] == ref[POS + ref_len - 1]:
        alt = "".join(rng.choice("ACGT") for _ in range(alt_len))
    hap = ref[:POS] + alt + ref[POS + ref_len :]
    reads = []
    for i, s in enumerate(range(POS - 140, POS + 1)):  # one start per base, both alleles
        reads.append(make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),)))
        reads.append(_alt_read(f"a{i}", hap, s, ref_len, alt_len))
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + ref_len], alt)])[0]
    _invariants(c)
    assert c.rd + c.ad >= 40, "too few informative reads to judge"
    assert abs(c.ad / (c.rd + c.ad) - 0.5) <= 0.06, (c.rd, c.ad)


# ── one rule across backends ─────────────────────────────────────────────


def test_both_backends_count_complex_reads_alike(tmp_path):
    ref = _ref()
    alt = _distinct_alt(ref, 3)
    off = alt[0] + ("A" if alt[1] != "A" else "C") + alt[2]
    off_hap = ref[:POS] + off + ref[POS + 2 :]
    near = [_alt_read(f"n{i}", off_hap, s, 2, 3) for i, s in enumerate(range(250, 258))]
    hap, reads = _delins_case(ref, alt, near)
    fa, bam = _files(tmp_path, ref, reads)
    v = _prepared(fa, ref[POS : POS + 2], alt)
    counts = {}
    for backend in ("pairhmm", "sw"):
        (c,) = gbcms_rs.count_bam_binned(
            bam,
            [v],
            [None],
            20,
            20,
            True,
            True,
            True,
            False,
            False,
            False,
            1,
            alignment_backend=backend,
        )
        counts[backend] = (c.rd, c.ad, c.partial_alt)
    assert counts["pairhmm"] == counts["sw"] == (10, 10, 8)


# ── the other complex branches ───────────────────────────────────────────


def test_an_mnp_read_with_an_indel_counts_only_as_an_exact_carrier(tmp_path):
    """A 3bp MNP. Reads carrying it written as a 1bp deletion plus a 1bp
    insertion inside the block (same bases) are carriers; reads with an extra
    inserted base inside the block carry a different allele."""
    ref = _ref()
    block = ref[POS : POS + 3]
    alt = "".join("A" if b != "A" else "C" for b in block)
    hap = ref[:POS] + alt + ref[POS + 3 :]
    reads = [
        make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),))
        for i, s in enumerate(range(240, 250))
    ]
    reads += [
        make_read(f"a{i}", hap[s : s + READ], s, ((0, READ),))
        for i, s in enumerate(range(240, 250))
    ]
    for i, s in enumerate(range(250, 256)):  # same haplotype, D1 then I1 inside the block
        left = POS + 1 - s
        reads.append(
            make_read(
                f"w{i}", hap[s : s + READ], s, ((0, left), (2, 1), (1, 1), (0, READ - left - 1))
            )
        )
    other = ref[:POS] + alt[:2] + "T" + alt[2] + ref[POS + 3 :]  # an extra T inside the block
    for i, s in enumerate(range(256, 264)):
        left = POS + 2 - s
        reads.append(
            make_read(f"x{i}", other[s : s + READ], s, ((0, left), (1, 1), (0, READ - left - 1)))
        )
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, block, alt)])[0]
    _invariants(c)
    assert c.ad == 16
    assert c.partial_alt == 8
    assert c.mnp_confirmed_alt <= c.ad


def test_a_deletion_with_a_substituted_anchor_counts_exact_carriers(tmp_path):
    """REF of 3 bases -> one base unlike REF's first and last (so no trimming
    makes it a plain deletion): 10 exact carriers; 8 reads with a different
    substituted base are partial evidence."""
    ref = _ref()
    ref3 = ref[POS : POS + 3]
    alt = next(b for b in "ACGT" if b not in (ref3[0], ref3[2]))
    other = next(b for b in "ACGT" if b not in (ref3[0], ref3[2], alt))
    reads = [
        make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),))
        for i, s in enumerate(range(240, 250))
    ]
    for tag, base, starts in (("a", alt, range(240, 250)), ("x", other, range(250, 258))):
        hap = ref[:POS] + base + ref[POS + 3 :]
        for i, s in enumerate(starts):
            left = POS + 1 - s
            reads.append(
                make_read(f"{tag}{i}", hap[s : s + READ], s, ((0, left), (2, 2), (0, READ - left)))
            )
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref3, alt)])[0]
    _invariants(c)
    assert (c.rd, c.ad) == (10, 10)
    assert c.partial_alt == 8
