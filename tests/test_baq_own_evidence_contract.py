"""BAQ spares the variant's own indel evidence (#166).

Heuristic BAQ lowers base quality by 20 within 5 bases of any insertion,
deletion or splice junction in a read, to protect SNV calls near indels. When the
read's indel IS the variant being counted, that penalty lands on the bases that
show the allele: at Q37 they fall below min BQ and every ALT read of a small indel
or complex variant was lost, while REF reads (no indel) kept full quality. BAQ now
leaves a read's insertions and deletions within the variant's own event alone; it
still penalizes splice junctions and indels elsewhere in the read.
"""

import pytest
from helpers import make_read
from test_complex_exact_contract import POS, READ, _alt_read, _files, _invariants, _prepared, _ref

from gbcms import _rs as gbcms_rs

SHAPES = {
    "insertion 1>2": (1, 2),
    "deletion 2>1": (2, 1),
    "delins 2>3": (2, 3),
    "delins 3>2": (3, 2),
}


def _case(ref_len, alt_len, q):
    ref = _ref()
    if ref_len == 1:  # pure insertion after POS's base
        ref_allele = ref[POS]
        alt = ref_allele + next(b for b in "ACGT" if b not in (ref[POS], ref[POS + 1]))
    elif alt_len == 1 and ref_len == 2:  # pure deletion of the base after POS
        ref_allele, alt = ref[POS : POS + 2], ref[POS]
    else:
        ref_allele = ref[POS : POS + ref_len]
        alt = next(
            a
            for a in ("GCC", "TGG", "CAA", "ATT", "GC", "TG", "CA", "AT")
            if len(a) == alt_len and a[0] != ref_allele[0] and a[-1] != ref_allele[-1]
        )
    hap = ref[:POS] + alt + ref[POS + ref_len :]
    starts = range(POS - 70, POS - 40)
    reads = [
        make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),), quals=[q] * READ)
        for i, s in enumerate(starts)
    ]
    for i, s in enumerate(starts):
        if ref_len == 1 or (
            alt_len == 1 and ref_len == 2
        ):  # anchor-preserved indel: I/D right after POS's base
            k = POS - s + 1
            d = alt_len - ref_len
            cig = ((0, k), (1, d), (0, READ - k - d)) if d > 0 else ((0, k), (2, -d), (0, READ - k))
            reads.append(make_read(f"a{i}", hap[s : s + READ], s, cig, quals=[q] * READ))
        else:
            reads.append(_alt_read(f"a{i}", hap, s, ref_len, alt_len, quals=[q] * READ))
    return ref, ref_allele, alt, reads


@pytest.mark.parametrize("shape", list(SHAPES))
def test_rna_baq_keeps_the_variants_own_indel_reads_at_q37(tmp_path, shape):
    ref_len, alt_len = SHAPES[shape]
    ref, ref_allele, alt, reads = _case(ref_len, alt_len, 37)
    fa, bam = _files(tmp_path, ref, reads)
    v = _prepared(fa, ref_allele, alt)
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
        mode="rna",
        apply_baq=True,
    )
    _invariants(c)
    assert (c.rd, c.ad) == (30, 30)


@pytest.mark.parametrize("shape", list(SHAPES))
def test_dna_baq_keeps_the_variants_own_indel_reads_at_q37(tmp_path, shape):
    """The same with BAQ opted into in DNA mode (the legacy parity path never runs
    BAQ, so the binned counter is called directly)."""
    ref_len, alt_len = SHAPES[shape]
    ref, ref_allele, alt, reads = _case(ref_len, alt_len, 37)
    fa, bam = _files(tmp_path, ref, reads)
    (c,) = gbcms_rs.count_bam_binned(
        bam,
        [_prepared(fa, ref_allele, alt)],
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
        apply_baq=True,
    )
    _invariants(c)
    assert (c.rd, c.ad) == (30, 30)


def test_baq_still_penalizes_an_snv_beside_another_indel(tmp_path):
    """An SNV three bases from a read's unrelated deletion: BAQ still lowers the
    SNV base (Q37 -> 17 < 20), so those ALT reads are not counted."""
    ref = _ref()
    snv_alt = next(b for b in "ACGT" if b != ref[POS])
    hap = ref[:POS] + snv_alt + ref[POS + 1 :]
    reads = []
    for i, s in enumerate(range(POS - 70, POS - 40)):
        k = POS - s - 3  # a 1-base deletion three bases before the SNV
        seq = hap[s : s + k] + hap[s + k + 1 : s + READ + 1]
        reads.append(make_read(f"a{i}", seq, s, ((0, k), (2, 1), (0, READ - k)), quals=[37] * READ))
    fa, bam = _files(tmp_path, ref, reads)
    (c,) = gbcms_rs.count_bam_binned(
        bam,
        [_prepared(fa, ref[POS], snv_alt)],
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
        apply_baq=True,
    )
    _invariants(c)
    assert c.ad == 0


@pytest.mark.parametrize("shape", list(SHAPES))
def test_baq_at_q40_is_unchanged(tmp_path, shape):
    """Q40 - 20 = 20 passes min BQ 20 either way: FORTE-like data count as before."""
    ref_len, alt_len = SHAPES[shape]
    ref, ref_allele, alt, reads = _case(ref_len, alt_len, 40)
    fa, bam = _files(tmp_path, ref, reads)
    v = _prepared(fa, ref_allele, alt)
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
        mode="rna",
        apply_baq=True,
    )
    assert (c.rd, c.ad) == (30, 30)
