"""Complex variants count carriers whose allele lies in soft-clipped bases (C12, #167).

Aligners soft-clip ALT reads whose allele sits near a read end, while REF reads
at the same positions align fully. A read whose aligned span stops before the
variant position is admitted to DP and to REF/ALT when the exact-carrier rule
decides it from its own bases, clipped ones included, with an aligned window
anchor, inside a well-defined fragment (paired, mate mapped in the opposite
orientation, TLEN set). Clipped bases past the fragment end are adapter and
cannot admit a read; a clipped read the rule cannot decide stays out, as before.
DNA only: an RNA read's clip may hold the next exon's bases.
"""

import random

from helpers import count_both, make_read
from test_complex_exact_contract import (
    POS,
    READ,
    _distinct_alt,
    _files,
    _invariants,
    _prepared,
    _ref,
)

from gbcms import _rs as gbcms_rs


def _paired(read, tlen=300, reverse=False):
    """Mark a read as one mate of a proper pair spanning `tlen` bases."""
    read.flag = 0x1 | 0x2 | 0x40 | (0x10 if reverse else 0x20)
    read.next_reference_id = 0
    read.next_reference_start = read.reference_start + (0 if reverse else max(0, tlen - READ))
    read.template_length = -tlen if reverse else tlen
    return read


def _right_clipped(name, seq_hap, s, tlen=300):
    """Aligned up to the base before POS; the rest (the event onward) soft-clipped."""
    left = POS - s
    return _paired(make_read(name, seq_hap[s : s + READ], s, ((0, left), (4, READ - left))), tlen)


def _case(ref_len=2, alt_len=3):
    ref = _ref()
    alt = _distinct_alt(ref, alt_len) if ref_len == 2 else None
    if alt is None:
        rng = random.Random(ref_len * 31 + alt_len)
        alt = ""
        while not alt or alt[0] == ref[POS] or alt[-1] == ref[POS + ref_len - 1]:
            alt = "".join(rng.choice("ACGT") for _ in range(alt_len))
    hap = ref[:POS] + alt + ref[POS + ref_len :]
    return ref, alt, hap


def _overlapping_ref(ref, n=10):
    return [
        _paired(make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),)))
        for i, s in enumerate(range(240, 240 + n))
    ]


def test_right_clipped_alt_carriers_count(tmp_path):
    """ALT reads aligned up to the event and clipped from it: they carry the whole
    ALT window in their clipped bases, so they are ALT and in DP."""
    ref, alt, hap = _case()
    reads = _overlapping_ref(ref) + [
        _right_clipped(f"a{i}", hap, s) for i, s in enumerate(range(POS - 90, POS - 70))
    ]
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 2], alt)])[0]
    _invariants(c)
    assert (c.rd, c.ad, c.dp) == (10, 20, 30)


def test_ref_reads_clipped_the_same_way_count_alike(tmp_path):
    """The same admission for REF: REF reads clipped at the event are REF."""
    ref, alt, _ = _case()
    reads = [_right_clipped(f"c{i}", ref, s) for i, s in enumerate(range(POS - 90, POS - 70))]
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 2], alt)])[0]
    _invariants(c)
    assert (c.rd, c.ad, c.dp) == (20, 0, 20)


def test_left_clipped_alt_carriers_count(tmp_path):
    """ALT reads whose alignment starts after the event, the event clipped at
    their start: judged back from the right anchor."""
    ref, alt, hap = _case()
    reads = _overlapping_ref(ref)
    after = POS + 2  # first REF base after the event; hap offset POS + 3
    for i, hs in enumerate(range(POS - 40, POS - 20)):  # read start on the haplotype
        clip = POS + 3 - hs
        reads.append(
            _paired(
                make_read(f"l{i}", hap[hs : hs + READ], after, ((4, clip), (0, READ - clip))),
                reverse=True,
            )
        )
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 2], alt)])[0]
    _invariants(c)
    assert (c.rd, c.ad) == (10, 20)


def test_clipped_bases_past_the_fragment_are_adapter(tmp_path):
    """Read-through: the fragment ends before the clipped window does, so those
    bases are adapter and cannot admit the read."""
    ref, alt, hap = _case()
    reads = _overlapping_ref(ref)
    reads += [
        _right_clipped(f"a{i}", hap, s, tlen=POS - s + 2)
        for i, s in enumerate(range(POS - 90, POS - 70))
    ]
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 2], alt)])[0]
    _invariants(c)
    assert (c.rd, c.ad, c.dp) == (10, 0, 10)


def test_reads_without_a_fragment_are_not_admitted_by_clips(tmp_path):
    """Unpaired reads have no defined fragment end: their clips admit nothing."""
    ref, alt, hap = _case()
    reads = _overlapping_ref(ref)
    for i, s in enumerate(range(POS - 90, POS - 70)):
        left = POS - s
        reads.append(make_read(f"u{i}", hap[s : s + READ], s, ((0, left), (4, READ - left))))
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 2], alt)])[0]
    _invariants(c)
    assert (c.rd, c.ad, c.dp) == (10, 0, 10)


def test_undecided_clipped_reads_stay_out_of_depth(tmp_path):
    """Clipped reads carrying a third allele: the rule cannot decide them, so DP
    is unchanged and they are not partial evidence."""
    ref, alt, _ = _case()
    other = next(
        a
        for a in ("TTT", "GGG", "CCC", "AAA")
        if a != alt and a[0] != ref[POS] and a[-1] != ref[POS + 1]
    )
    third = ref[:POS] + other + ref[POS + 2 :]
    reads = _overlapping_ref(ref) + [
        _right_clipped(f"x{i}", third, s) for i, s in enumerate(range(POS - 90, POS - 70))
    ]
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 2], alt)])[0]
    _invariants(c)
    assert (c.rd, c.ad, c.dp, c.partial_alt) == (10, 0, 10, 0)


def test_rna_reads_are_not_admitted_by_clips(tmp_path):
    """RNA: a clip may hold the next exon's bases, so clips admit nothing."""
    ref, alt, hap = _case()
    reads = _overlapping_ref(ref) + [
        _right_clipped(f"a{i}", hap, s) for i, s in enumerate(range(POS - 90, POS - 70))
    ]
    fa, bam = _files(tmp_path, ref, reads)
    v = _prepared(fa, ref[POS : POS + 2], alt)
    (c,) = gbcms_rs.count_bam_binned(
        bam, [v], [None], 20, 20, True, True, True, False, False, False, 1, mode="rna"
    )
    assert (c.rd, c.ad) == (10, 0)


def test_a_long_deletions_split_carriers_count(tmp_path):
    """A 60bp REF to a 10bp ALT, carriers aligned up to the event with the ALT and
    the far flank clipped (a split read's primary alignment): judged at the
    junction they show, so they are ALT."""
    ref, alt, hap = _case(60, 10)
    reads = _overlapping_ref(ref) + [
        _right_clipped(f"a{i}", hap, s) for i, s in enumerate(range(POS - 80, POS - 60))
    ]
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 60], alt)])[0]
    _invariants(c)
    assert (c.rd, c.ad) == (10, 20)


def test_vaf_is_unbiased_when_aligners_clip_alt_near_read_ends(tmp_path):
    """50% sample, one read per start on each allele; ALT reads with fewer than 20
    bases past the event are clipped from it, as aligners do. Admission restores
    the balance the clipping took away."""
    ref, alt, hap = _case()
    reads = []
    for i, s in enumerate(range(POS - READ + 5, POS + 1)):
        reads.append(_paired(make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),))))
        left = POS - s
        right = READ - left - 3
        if right < 20:
            reads.append(_right_clipped(f"a{i}", hap, s))
        else:
            reads.append(
                _paired(
                    make_read(
                        f"a{i}", hap[s : s + READ], s, ((0, left), (1, 1), (0, READ - left - 1))
                    )
                )
            )
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 2], alt)])[0]
    _invariants(c)
    assert c.rd + c.ad >= 80
    assert abs(c.ad / (c.rd + c.ad) - 0.5) <= 0.03, (c.rd, c.ad)


def test_unclipped_reads_starting_inside_a_long_deletion_stay_out(tmp_path):
    """REF reads starting inside a 60bp deletion hold its right junction with
    aligned bases only: nothing of theirs is clipped, and no ALT read can start
    there, so they stay outside depth as before (admission is for clipped bases)."""
    ref, alt, _ = _case(60, 10)
    reads = _overlapping_ref(ref)
    reads += [
        _paired(make_read(f"i{i}", ref[s : s + READ], s, ((0, READ),)))
        for i, s in enumerate(range(POS + 10, POS + 30))
    ]
    fa, bam = _files(tmp_path, ref, reads)
    c = count_both(bam, [_prepared(fa, ref[POS : POS + 60], alt)])[0]
    _invariants(c)
    assert (c.rd, c.ad, c.dp) == (10, 0, 10)
