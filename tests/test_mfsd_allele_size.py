"""mFSD sizes an ALT fragment as the molecule that carries the ALT (6.7.0).

A carrier of a deletion written as a soft clip (or judged by its bases) has no D op
in either mate, so TLEN, a reference span, holds the deleted bases: before 6.7.0 its
fragment read longer than the molecule by the deletion's length, and a fragment of a
deletion longer than about 830 bases fell out of mFSD altogether. Every molecule here
is 167 bases; REF and ALT fragments must size alike. Synthetic and PHI-free.
"""

import random

import pytest
from helpers import make_read, write_contig

from gbcms import _rs

A, L, MOLECULE, READ, N = 600, 3000, 167, 100, 12


def _contig():
    rng = random.Random(254)
    seq = [rng.choice("CGT") for _ in range(L)]
    return "".join(seq)


def _pair(name, seq1, pos1, cigar1, seq2, pos2, cigar2, tlen):
    """A forward first read and its reverse mate."""
    r1 = make_read(name, seq1, pos1, cigar1, flag=0x1 | 0x2 | 0x20 | 0x40)
    r2 = make_read(name, seq2, pos2, cigar2, flag=0x1 | 0x2 | 0x10 | 0x80)
    r1.next_reference_id = r2.next_reference_id = 0
    r1.next_reference_start, r2.next_reference_start = pos2, pos1
    r1.template_length, r2.template_length = tlen, -tlen
    return [r1, r2]


@pytest.mark.parametrize("n", [60, 300, 1000])
def test_a_clip_written_carrier_is_sized_as_its_molecule(tmp_path, n):
    """Each carrier's first read is clipped at the junction (its aligned bases end at
    the deletion), its mate aligned past the deletion: TLEN is 167 + n."""
    contig = _contig()
    if contig[A - 1] == contig[A + n - 1]:  # no slide: the deletion sits at A
        contig = contig[: A - 1] + next(b for b in "CGT" if b != contig[A + n - 1]) + contig[A:]
    hap = contig[:A] + contig[A + n :]
    reads = []
    for i in range(N):
        s = A - 60 + 2 * i
        # REF molecule: both reads on the reference
        reads += _pair(
            f"ref{i}",
            contig[s : s + READ],
            s,
            ((0, READ),),
            contig[s + MOLECULE - READ : s + MOLECULE],
            s + MOLECULE - READ,
            ((0, READ),),
            MOLECULE,
        )
        # ALT molecule: the first read clipped at the junction, its mate past the deletion
        left = A - s
        reads += _pair(
            f"alt{i}",
            hap[s : s + READ],
            s,
            ((0, left), (4, READ - left)),
            hap[s + MOLECULE - READ : s + MOLECULE],
            s + MOLECULE - READ + n,
            ((0, READ),),
            MOLECULE + n,
        )
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A - 1, contig[A - 1 : A + n], contig[A - 1], "DEL")],
        fa,
        5,
        False,
        1,
        True,
    )
    (c,) = _rs.count_bam_binned(
        bam, [pv.variant], [None], 20, 20, True, True, True, False, False, False, 1, mfsd=True
    )
    assert (c.rd, c.ad) == (N, N)
    assert (c.mfsd_ref_count, c.mfsd_alt_count) == (N, N)
    assert c.mfsd_ref_mean == pytest.approx(MOLECULE)
    assert c.mfsd_alt_mean == pytest.approx(MOLECULE)
