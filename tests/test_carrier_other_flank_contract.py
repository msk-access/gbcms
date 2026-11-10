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


def _flank_like_insert(contig, k=20, subs=(3, 9, 15)):
    """A 40-base insert after base A whose last `k` bases are the left flank's last `k`
    (ending at the anchor) with substitutions: an insert whose end resembles the flank
    (a repeat element beside its own kind, a near target-site duplication), as aligners
    misplace."""
    rng = random.Random(279)
    head = "".join(rng.choice("ACG") for _ in range(40 - k))
    tail = list(contig[A - k + 1 : A + 1])
    for i in subs:
        tail[i] = {"C": "A", "G": "A", "T": "A", "A": "C"}[tail[i]]
    return head + "".join(tail)


def _over_flank(name, ins, contig, k, flank, clip, right):
    """A carrier as an aligner writes it when the insert's end resembles the left flank:
    `flank` bases of the left flank and the insert's first bases soft-clipped (`clip`
    of the insert's head), its last `k` aligned onto the flank's last `k` (as mismatches),
    then the right flank."""
    head = ins[40 - k - clip : 40 - k]
    seq = contig[A + 1 - flank : A + 1] + head + ins[40 - k :] + contig[A + 1 : A + 1 + right]
    return _pair(make_read(name, seq, A - k + 1, ((4, flank + clip), (0, k + right))))


def _flank_like_row(contig, ins, carriers):
    """Six REF reads and six carriers written with the insertion, besides `carriers`."""
    hap = contig[: A + 1] + ins + contig[A + 1 :]
    reads = list(carriers)
    for i in range(6):
        s = A - 60 - 3 * i
        reads.append(_pair(make_read(f"ref{i}", contig[s : s + 150], s, ((0, 150),))))
        n = A + 1 - s
        cig = ((0, n), (1, 40), (0, 150 - n - 40))
        reads.append(_pair(make_read(f"alt{i}", hap[s : s + 150], s, cig)))
    return reads


def test_a_carrier_reaching_past_an_insert_aligned_over_the_flank_is_alt(tmp_path):
    """A carrier whose clip holds the left flank and the insert's head, its aligned bases
    over the left flank the insert's end (it resembles the flank), fits neither allele
    read from that flank. Read whole from the right flank, at the ALT's length, its
    bases are the ALT's, and the read reaches past the insert's start, as a carrier read
    from the left flank does: ALT (#279, operator 2026-10-16: spanning carriers only)."""
    contig = _contig()
    ins = _flank_like_insert(contig)
    # their clips reach past the window's left flank (the tract grows four bases left,
    # through the bases the insert's end shares with the flank)
    carriers = [_over_flank(f"sp{i}", ins, contig, 20, 8 + i, 20, 60) for i in range(6)]
    c = _count(tmp_path, contig, _flank_like_row(contig, ins, carriers), contig[A], contig[A] + ins)
    assert (c.rd, c.ad, c.partial_alt) == (6, 12, 0)


def test_a_carrier_starting_inside_the_insert_is_depth_only_never_ref(tmp_path):
    """A carrier starting inside the insert, its end aligned over the left flank: its two
    bases at the REF window happen to fit REF (the insert ends in the anchor's base), but
    its bases before them are the insert's, which the ALT read from the right flank fits
    and REF does not. It shows only the right junction, which a REF molecule starting
    there has no counterpart of, so it is depth only (RJ-23's one reading flank), not
    REF."""
    contig = _contig()
    ins = _flank_like_insert(contig)
    carriers = [_over_flank(f"in{i}", ins, contig, 20, 0, 8 + 2 * i, 60) for i in range(6)]
    c = _count(tmp_path, contig, _flank_like_row(contig, ins, carriers), contig[A], contig[A] + ins)
    assert (c.rd, c.ad, c.partial_alt) == (6, 6, 0)


def test_reads_of_another_insert_with_the_same_flank_like_end_are_not_alt(tmp_path):
    """Guard: written the same way, but their clipped bases are another insert's: the
    ALT placed from the right flank does not fit them, so none is ALT."""
    contig = _contig()
    ins = _flank_like_insert(contig)
    other = "".join("T" if b != "T" else "C" for b in ins[:20]) + ins[20:]
    reads = [_over_flank(f"o{i}", other, contig, 20, 8 + i, 20, 60) for i in range(6)]
    c = _count(tmp_path, contig, reads, contig[A], contig[A] + ins)
    assert c.ad == 0


def test_ref_molecules_behind_a_poor_flank_are_not_placed_from_the_other_one(tmp_path):
    """Guard: REF molecules whose bases over the left flank's end are random at Q12 and
    whose clip is junk at Q12: their confident bases show nothing either way there, so
    no placement from the right flank makes them ALT."""
    contig = _contig()
    ins = _flank_like_insert(contig)
    rng = random.Random(7)
    reads = []
    for i in range(6):
        clip = 12 + 2 * i
        junk = "".join(rng.choice("ACGT") for _ in range(clip + 20))
        seq = junk + contig[A + 1 : A + 61]
        quals = [12] * (clip + 20) + [30] * 60
        reads.append(_pair(make_read(f"p{i}", seq, A - 19, ((4, clip), (0, 80)), quals=quals)))
    c = _count(tmp_path, contig, reads, contig[A], contig[A] + ins)
    assert c.ad == 0


def test_a_ref_molecule_starting_inside_a_deletion_whose_end_resembles_the_flank_is_not_alt(
    tmp_path,
):
    """The mirror for deletions: 40 deleted bases whose end resembles the left flank. A
    REF molecule starting inside them, their end aligned over the left flank as
    mismatches, fits the ALT's two bases at the anchor; its bases before them are the
    deleted sequence, which REF read from the right flank fits and the ALT does not. It
    shows only REF's right junction: depth only, never ALT."""
    base = _contig()
    dele = _flank_like_insert(base)
    contig = base[: A + 1] + dele + base[A + 1 :]  # REF holds the 40 bases
    reads = []
    for i in range(6):
        s = A - 60 - 3 * i
        reads.append(_pair(make_read(f"ref{i}", contig[s : s + 150], s, ((0, 150),))))
        n = A + 1 - s
        alt = contig[s : A + 1] + contig[A + 41 : A + 41 + 150 - n]
        reads.append(_pair(make_read(f"alt{i}", alt, s, ((0, n), (2, 40), (0, 150 - n)))))
        clip = 8 + 2 * i
        seq = dele[20 - clip : 20] + dele[20:] + contig[A + 41 : A + 101]
        reads.append(_pair(make_read(f"in{i}", seq, A - 19, ((4, clip), (0, 80)))))
    c = _count(tmp_path, contig, reads, contig[A : A + 41], contig[A])
    assert (c.rd, c.ad) == (6, 6)
