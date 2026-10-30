"""The sequence judge (`--read-judge sequence`, #304): an indel beside X's site.

Mechanism. An indel a read carries outside X's site is a separate event (RJ-8), wherever
its aligner wrote it: one placed by its bases beside the site is beside it, even where
the read writes it inside the site with mismatches its slide resolves (RJ-27). The read
is judged at X by its bases there, REF where they are REF's. Before this the judge took
no indel the reads write beside X's site as a context, so such a read carried a gap on
every candidate sequence alike, or, written inside the site, sat between REF and X and
was depth only.

Shapes: X a +A in an A run followed by a C run and a T; beside it, past the C run, the T
doubled or deleted (an indel whose every placement lies clear of X's discrimination
window: one that could sit on the run's first C would be inside it, another allele by
RJ-7); half the reads carrying it written there, half written inside X's run as an
indel after the G with the mismatches the shift leaves; both strands. A known partial
stands as a strict xfail: where no read writes the indel beside the site, the judge has
no candidate holding it (placing a read's own indel by its bases is the rules' RJ-27
machinery, not yet the judge's).
"""

import random

import pytest
from helpers import count_checked, make_read, write_contig

from gbcms import _rs

L = 100
Z = 420  # G>GA at Z: G, A10, C5, T


def _contig():
    rng = random.Random(8)
    c = [rng.choice("ACGT") for _ in range(900)]
    c[Z - 1 : Z + 18] = "TG" + "A" * 10 + "C" * 5 + "TG"
    return "".join(c)


CONTIG = _contig()


def _beside(kind):
    """The REF molecule with the indel beside X's site: the T after the C run doubled or
    deleted. (its bases, the genome position where it is written)."""
    t = Z + 16  # the T after the C run (0-based)
    if kind == "insertion":
        return CONTIG[:t] + "T" + CONTIG[t:], t
    return CONTIG[:t] + CONTIG[t + 1 :], t


def _reads(kind, n, written_beside, first=0):
    hap, j = _beside(kind)
    out = []
    for i in range(first, first + n):
        s = Z - 40 - (i % 5)
        seq = hap[s : s + L]
        flag = 16 if i % 2 else 0
        if written_beside:
            left = j - s
            cig = (
                ((0, left), (1, 1), (0, L - left - 1))
                if kind == "insertion"
                else ((0, left), (2, 1), (0, L - left))
            )
            if kind == "deletion":
                seq = hap[s : s + L]
        else:
            # written inside X's run: as X's +A (insertion) or a −A (deletion) after the G,
            # the shifted bases read as mismatches
            left = Z + 1 - s
            cig = (
                ((0, left), (1, 1), (0, L - left - 1))
                if kind == "insertion"
                else ((0, left), (2, 1), (0, L - left))
            )
        out.append(
            make_read(f"{kind[0]}{i}{'b' if written_beside else 'w'}", seq, s, cig, flag=flag)
        )
    return out


def _count(tmp_path, reads):
    fa, bam = write_contig(tmp_path, CONTIG, reads)
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", Z, "G", "GA", "INS")], fa, 5, False, 1, True, read_judge="sequence"
    )
    (c,) = count_checked(bam, [pv.variant], read_judge="sequence", reference_fasta=fa)
    assert c.dp >= c.rd + c.ad and c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    return c


@pytest.mark.parametrize("kind", ["insertion", "deletion"])
def test_reads_carrying_an_indel_beside_the_site_are_ref_however_written(tmp_path, kind):
    """4 reads write the indel where it lies, 4 inside X's run: all 8 are REF."""
    c = _count(tmp_path, _reads(kind, 4, True) + _reads(kind, 4, False, first=4))
    assert (c.rd, c.ad, c.partial_alt) == (8, 0, 0)


def test_reads_writing_a_beside_deletion_only_inside_the_site_are_ref(tmp_path):
    """No read writes the deletion where it lies: written as a deletion after the G, each
    still fits REF with one gap better than X by its bases: REF."""
    c = _count(tmp_path, _reads("deletion", 8, False))
    assert (c.rd, c.ad, c.partial_alt) == (8, 0, 0)


@pytest.mark.xfail(
    strict=True,
    reason="known partial: the judge does not place a read's own indel by its bases (RJ-27)",
)
def test_reads_writing_a_beside_insertion_only_inside_the_site_are_ref(tmp_path):
    """No read writes the insertion where it lies: written as X's +A, the shift's two
    mismatches cost X about what REF's gap costs, so with no candidate holding the
    insertion beside the site each read is depth only, not REF."""
    c = _count(tmp_path, _reads("insertion", 8, False))
    assert (c.rd, c.ad, c.partial_alt) == (8, 0, 0)
