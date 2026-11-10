"""Shapes the 6.7.0 gate found at pure indels judged by their bases (#278).

- A germline change that joins the run an insertion slides through to the next run
  (CAA -> CAG between two CAG runs), on the REF molecules and every carrier: REF
  molecules and carriers fit as far as the guard compares, so it reads on until they
  differ; the change is then background, masked, and both keep their call.
- A REF-length read whose substitution spells a deletion's shorter ALT window (AT>CA
  next to a CAG run turns the run's last ATG into CAG): its bases past the windows are
  unshifted REF, so it is no carrier.
- The same REF-length reads written as the deletion, their tail clipped where the run
  shows one unit more: their bases place no deletion, so they are read on as unwritten
  reads are, and are no carriers.
- A REF read with two indels upstream of the windows: each is a separate event;
  summed into one length change they would fit no placement and make it another allele.
"""

import random

from helpers import count_checked, make_read, write_contig

from gbcms import _rs as gbcms_rs

READ = 120


def _flank(seed: int, n: int = 300) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _count(tmp_path, contig, reads, pos, ref, alt):
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", pos, ref, alt, "X")], fa, 5, False, 1, True
    )
    (c,) = count_checked(bam, [pv.variant])
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    return c.rd, c.ad, c.partial_alt


def test_a_germline_change_joining_an_insertions_run_to_the_next_is_background(tmp_path):
    """C>CCAG at CAG x7, CAA, CAG x5: a third of the REF molecules and every carrier read
    the CAA as CAG. Both alleles keep their reads: REF 30, ALT 20."""
    left = _flank(4) + "CAGCTAAACCAGTTC"
    contig = left + "CAG" * 7 + "CAA" + "CAG" * 5 + "C" + _flank(5)
    p, germ = len(left) - 1, len(left) + 23  # the anchor C; the CAA's last A

    def hap(carrier: bool, germline: bool) -> str:
        h = list(contig)
        if germline:
            h[germ] = "G"
        if carrier:
            h[p + 1 : p + 1] = list("CAG")
        return "".join(h)

    reads = []
    for i in range(30):
        s = p - 60 + i
        reads.append(make_read(f"r{i}", hap(False, i < 10)[s : s + READ], s, ((0, READ),)))
    for i in range(20):
        s = p - 60 + i
        left_len = p + 1 - s
        cigar = ((0, left_len), (1, 3), (0, READ - left_len - 3))
        reads.append(make_read(f"a{i}", hap(True, True)[s : s + READ], s, cigar))
    assert _count(tmp_path, contig, reads, p, "C", "CCAG") == (30, 20, 0)


def test_clean_ref_reads_ending_inside_the_joined_run_are_ref(tmp_path):
    """The same locus, plus REF molecules without the germline change whose reads end
    inside the run it joins (a few bases past the CAA). The masks take either base at
    the change, so no masked base tells the alleles apart before the joined run's end;
    but these reads show the reference's A where the ALT reads its inserted CAG's G:
    their own bases decide, REF (they were depth only)."""
    left = _flank(4) + "CAGCTAAACCAGTTC"
    contig = left + "CAG" * 7 + "CAA" + "CAG" * 5 + "C" + _flank(5)
    p, germ = len(left) - 1, len(left) + 23  # the anchor C; the CAA's last A

    def hap(carrier: bool, germline: bool) -> str:
        h = list(contig)
        if germline:
            h[germ] = "G"
        if carrier:
            h[p + 1 : p + 1] = list("CAG")
        return "".join(h)

    reads = []
    for i in range(30):
        s = p - 60 + i
        reads.append(make_read(f"r{i}", hap(False, i < 10)[s : s + READ], s, ((0, READ),)))
    for i in range(20):
        s = p - 60 + i
        left_len = p + 1 - s
        cigar = ((0, left_len), (1, 3), (0, READ - left_len - 3))
        reads.append(make_read(f"a{i}", hap(True, True)[s : s + READ], s, cigar))
    for i in range(10):  # clean REF reads ending 3-12 bases past the CAA
        end = germ + 4 + i
        s = end - READ
        reads.append(make_read(f"e{i}", contig[s:end], s, ((0, READ),)))
    assert _count(tmp_path, contig, reads, p, "C", "CCAG") == (40, 20, 0)


def test_a_ref_length_substitution_spelling_a_deletions_window_is_not_alt(tmp_path):
    """TATG>T (ATG x3 then CAG x6): REF reads reading the last ATG as CAG fit the short
    ALT window, but past it their bases are unshifted REF: not carriers."""
    left = _flank(7) + "GGGTGGCT"
    run = "ATG" * 3 + "CAG" * 6 + "CAAC"
    contig = left + run + _flank(8)
    p = len(left) - 1  # the anchor T
    mnp = list(contig)
    mnp[p + 7 : p + 9] = list("CA")  # the third ATG's AT read as CA
    mnp = "".join(mnp)
    alt = contig[: p + 1] + contig[p + 4 :]
    reads = [
        make_read(f"r{i}", contig[p - 60 + i : p - 60 + i + READ], p - 60 + i, ((0, READ),))
        for i in range(10)
    ]
    for i in range(8):
        s = p - 60 + i
        reads.append(make_read(f"m{i}", mnp[s : s + READ], s, ((0, READ),)))
    for i in range(6):
        s = p - 60 + i
        left_len = p + 1 - s
        reads.append(
            make_read(f"a{i}", alt[s : s + READ], s, ((0, left_len), (2, 3), (0, READ - left_len)))
        )
    rd, ad, _ = _count(tmp_path, contig, reads, p, contig[p : p + 4], contig[p])
    assert (rd, ad) == (10, 6)


def test_a_ref_length_substitution_written_as_the_deletion_is_not_alt(tmp_path):
    """TATG>T (ATG x3 then CAG x6): the REF-length reads that read the last ATG as CAG,
    written by the aligner as the deletion (M, D3 at the anchor, M through the run) with
    the tail clipped where they read one CAG more than the ALT. Their bases place no
    3-base deletion (the clip shows the run's REF length): read on past the windows,
    they fit REF's continuation, not the ALT's, so they are no carriers."""
    left = _flank(7) + "GGGTGGCT"
    run = "ATG" * 3 + "CAG" * 6 + "CAAC"
    contig = left + run + _flank(8)
    p = len(left) - 1  # the anchor T
    mnp = list(contig)
    mnp[p + 7 : p + 9] = list("CA")
    mnp = "".join(mnp)
    alt = contig[: p + 1] + contig[p + 4 :]
    reads = [
        make_read(f"r{i}", contig[p - 60 + i : p - 60 + i + READ], p - 60 + i, ((0, READ),))
        for i in range(10)
    ]
    for i in range(8):
        s = p - 60 + i
        left_len = p + 1 - s
        aligned = 6 + 18  # ATG ATG and six CAG fit REF past the deletion
        cigar = ((0, left_len), (2, 3), (0, aligned), (4, READ - left_len - aligned))
        reads.append(make_read(f"m{i}", mnp[s : s + READ], s, cigar))
    for i in range(6):
        s = p - 60 + i
        left_len = p + 1 - s
        reads.append(
            make_read(f"a{i}", alt[s : s + READ], s, ((0, left_len), (2, 3), (0, READ - left_len)))
        )
    rd, ad, _ = _count(tmp_path, contig, reads, p, contig[p : p + 4], contig[p])
    assert (rd, ad) == (10, 6)


def test_two_indels_upstream_are_separate_events(tmp_path):
    """A 6-base deletion in unique sequence: REF reads that delete 11 bases 40 upstream
    and 1 base 15 upstream show REF across the windows: REF, as their aligner wrote
    them (summed into one 12-base change, no placement fits and they were partial)."""
    left = _flank(9, 400)
    contig = left + "GACTTGACC" + _flank(10)
    p = len(left)  # the anchor G
    alt = contig[: p + 1] + contig[p + 7 :]
    reads = [
        make_read(f"r{i}", contig[p - 80 + i : p - 80 + i + READ], p - 80 + i, ((0, READ),))
        for i in range(10)
    ]
    up = contig[: p - 40] + contig[p - 29 : p - 15] + contig[p - 14 :]  # 11 then 1 deleted
    for i in range(6):
        s = p - 80 + i
        a = (p - 40) - s
        cigar = ((0, a), (2, 11), (0, 14), (2, 1), (0, READ - a - 14))
        reads.append(make_read(f"u{i}", up[s : s + READ], s, cigar))
    for i in range(6):
        s = p - 80 + i
        left_len = p + 1 - s
        reads.append(
            make_read(f"a{i}", alt[s : s + READ], s, ((0, left_len), (2, 6), (0, READ - left_len)))
        )
    assert _count(tmp_path, contig, reads, p, contig[p : p + 7], contig[p]) == (16, 6, 0)


def test_a_germline_change_inside_a_deletions_bases_is_background_on_ref(tmp_path):
    """A 6-base deletion in unique sequence: 12 of the 30 REF molecules carry a germline
    change inside the deleted bases, which the ALT has none of. They are REF molecules:
    background on REF, masked there, so REF keeps them (VAF 10/40, not 10/28)."""
    left = _flank(11, 400)
    contig = left + "GACTTGACC" + _flank(12)
    p = len(left)  # the anchor G
    germ = p + 3  # the third deleted base (T)
    alt = contig[: p + 1] + contig[p + 7 :]
    reads = []
    for i in range(30):
        s = p - 60 + i
        h = list(contig)
        if i < 12:
            h[germ] = "C"
        reads.append(make_read(f"r{i}", "".join(h)[s : s + READ], s, ((0, READ),)))
    for i in range(10):
        s = p - 60 + i
        left_len = p + 1 - s
        cigar = ((0, left_len), (2, 6), (0, READ - left_len))
        reads.append(make_read(f"a{i}", alt[s : s + READ], s, cigar))
    assert _count(tmp_path, contig, reads, p, contig[p : p + 7], contig[p]) == (30, 10, 0)
