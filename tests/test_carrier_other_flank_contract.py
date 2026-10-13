"""How the exact-carrier rule places a read whose reading flank is not aligned
(RJ-23, amended 2026-10-12; docs/reference/read-judgment.md).

A read whose reading flank lies in its soft clip, under its own deletion, or at the
event's other copy is placed from the other flank, each allele's window at its own
length. A placement counts only where the read's confident bases favour that allele
over the other one placed from the same anchor by at least one base: a molecule whose
bases there are below --min-baseq (or N) shows nothing either way. Every case counts
through `count_checked` (binning invariance). Synthetic and PHI-free.
"""

import itertools
import random

from helpers import count_checked, make_read, write_contig
from read_judgment_cases import _grows

from gbcms import _rs

A, L = 400, 900


def _contig(seed=253):
    rng = random.Random(seed)
    return "".join(rng.choice("CGT") for _ in range(L))


def _pair(read, forward=True):
    read.flag = 0x1 | 0x2 | (0x20 | 0x40 if forward else 0x10 | 0x80)
    read.next_reference_id = 0
    read.next_reference_start = read.reference_start + (150 if forward else -200)
    read.template_length = 300 if forward else -300
    return read


def _count(tmp_path, contig, reads, ref, alt, pos=A):
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants([_rs.Variant("1", pos, ref, alt, "X")], fa, 5, False, 1, True)
    (c,) = count_checked(bam, [pv.variant])
    return c


def test_ref_molecules_behind_low_quality_bases_are_never_alt(tmp_path):
    """REF molecules of a 28-to-6 delins starting inside the event (from twelve bases
    before its end to six past it), their first twelve bases random at Q12 (a NovaSeq
    bin, below --min-baseq): soft-clipped, or aligned in full as an aligner writes a
    reverse read's poor 3' tail. Whatever placement fits them, their confident bases
    are REF's: none is ALT."""
    contig, n = _contig(), 28
    rng = random.Random(1)
    reads = []
    for i in range(240):
        g0 = A + n - 12 + i % 19  # aligned start of the confident bases
        junk = "".join(rng.choice("ACGT") for _ in range(12))
        seq = junk + contig[g0 : g0 + 88]
        quals = [12] * 12 + [30] * 88
        if i % 2:
            r = make_read(f"c{i}", seq, g0, ((4, 12), (0, 88)), quals=quals)
            reads.append(_pair(r))
        else:
            r = make_read(f"a{i}", seq, g0 - 12, ((0, 100),), quals=quals)
            reads.append(_pair(r, forward=False))
    c = _count(tmp_path, contig, reads, contig[A : A + n], "ATATGA")
    assert c.ad == 0, (c.rd, c.ad, c.dp)


def test_a_carrier_written_as_two_deletions_around_retained_alt_bases_is_alt(tmp_path):
    """A 28-to-6 delins whose ALT is REF's bases 19-24: an aligner writes a carrier
    as nineteen bases deleted, the six kept, three deleted (as BWA did on a panel
    delins), with no mismatch. Its bases are the ALT's, so it is ALT, as REF reads
    from the same starts are REF."""
    contig, n = list(_contig()), 28
    ref = "".join(contig[A : A + n])
    alt = ref[19:25]
    contig = "".join(contig)
    hap = contig[:A] + alt + contig[A + n :]
    reads = []
    for i in range(6):
        s = A - 60 + 3 * i
        seq = hap[s : s + 120]
        cig = ((0, A - s), (2, 19), (0, 6), (2, 3), (0, 120 - (A - s) - 6))
        reads.append(_pair(make_read(f"alt{i}", seq, s, cig)))
        reads.append(_pair(make_read(f"ref{i}", contig[s : s + 120], s, ((0, 120),))))
    c = _count(tmp_path, contig, reads, ref, alt)
    assert (c.rd, c.ad) == (6, 6)


def test_a_germline_change_of_two_bases_in_the_flank_keeps_both_alleles(tmp_path):
    """Every read, REF and ALT, shows two adjacent changes in the left flank (a
    germline MNP): background on both haplotypes, masked as a germline SNP is, not
    another allele of the event. Both alleles still count. The contig is the first
    seed where the event grows through no repeat on its left on either allele, so the
    flank is the two bases before POS."""
    n, alt = 28, "ATATGA"
    for seed in itertools.count(253):
        contig = _contig(seed)
        hap = contig[:A] + alt + contig[A + n :]
        if not any(_grows(x, A, max(6, n - len(alt)), True) for x in (contig, hap)):
            break

    def mnp(seq, s):
        b = list(seq)
        for g in (A - 2, A - 1):
            b[g - s] = "A"  # the contig has no A
        return "".join(b)

    reads = []
    for i in range(6):
        s = A - 60 + 3 * i
        cig = ((0, A - s), (1, len(alt)), (2, n), (0, 120 - (A - s) - len(alt)))
        reads.append(_pair(make_read(f"alt{i}", mnp(hap[s : s + 120], s), s, cig)))
        reads.append(_pair(make_read(f"ref{i}", mnp(contig[s : s + 120], s), s, ((0, 120),))))
    c = _count(tmp_path, contig, reads, contig[A : A + n], alt)
    assert (c.rd, c.ad) == (6, 6)


def test_a_long_event_read_from_the_right_counts_clip_written_carriers(tmp_path):
    """GTT>CA ending a 70-G run is read from the right flank. Carriers starting before
    the run and ending two to nine bases past the ALT, written with the gap or with
    their ALT and right flank soft-clipped (placed from the left flank), count alike:
    their bases are the ALT's either way. REF reads end where the carriers do."""
    c = list(_contig())
    c[A - 69 : A + 3] = "G" * 70 + "TT"
    contig = "".join(c)
    hap = contig[:A] + "CA" + contig[A + 3 :]
    counts = {}
    for form in ("gap", "clip"):
        reads = []
        for i, past in enumerate(range(2, 10)):
            s = A - 85 - i
            if form == "gap":
                cig = ((0, A - s), (1, 2), (2, 3), (0, past))
            else:
                cig = ((0, A - s), (4, 2 + past))
            reads.append(_pair(make_read(f"alt{i}", hap[s : A + 2 + past], s, cig)))
            ref = contig[s : A + 3 + past]
            reads.append(_pair(make_read(f"ref{i}", ref, s, ((0, len(ref)),))))
        d = tmp_path / form
        d.mkdir()
        got = _count(d, contig, reads, "GTT", "CA")
        counts[form] = (got.rd, got.ad)
    assert counts["gap"][1] > 0 and counts["gap"][0] == counts["gap"][1], counts
    assert counts["clip"] == counts["gap"], counts
