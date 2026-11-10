"""An RNA pure indel at an exon start whose run continues into the previous exon.

A>AA where exon 2 starts with an A run and exon 1 ends in A's: on the mRNA the run
spans the junction (eight A's on REF, nine on the ALT). An exon edge cuts the left
flank, so the windows read from the right; but a read starting at the exon start, or
one whose short exon-1 overhang the aligner writes into the intron as mismatches,
reads from inside that run: its A's cannot tell the alleles apart. The windows of a
read starting at the exon start, or within the five bases of overhang an aligner
writes into the intron, end at the exon edge, as a spliced read's do, so REF and ALT
molecules starting there are judged alike, and neither counts. A spliced read reads
the run whole on its spliced reference and counts; a pre-mRNA read through the intron
reads the intron's base before the run.
"""

import random

from helpers import count_checked, make_read, write_contig

from gbcms import _rs as gbcms_rs

RL = 100


def _flank(seed: int, n: int) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("CGT") for _ in range(n))


EX1 = _flank(41, 260) + "CAAAA"
INTRON = "GT" + _flank(42, 36) + "AG"
EX2 = "AAAACTGGA" + _flank(43, 300)
CONTIG = EX1 + INTRON + EX2
I0, E = len(EX1), len(EX1) + len(INTRON)  # the intron, [I0, E); exon 2 from E
GAP = E - I0
MRNA_REF = EX1 + EX2
MRNA_ALT = EX1 + EX2[:1] + "A" + EX2[1:]


def _count(tmp_path, reads):
    fa, bam = write_contig(tmp_path, CONTIG, reads, "rna")
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", E, "A", "AA", "X")], fa, 5, False, 1, True
    )
    (c,) = count_checked(bam, [pv.variant], mode="rna", reference_fasta=fa)
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    return c.rd, c.ad, c.partial_alt


def _spliced(name, mrna, s):
    """An mRNA read from exon 1 at `s`, spliced into exon 2 (the ALT's extra A in it)."""
    left = I0 - s
    seq = mrna[s : s + RL]
    if mrna is MRNA_REF:
        return make_read(name, seq, s, ((0, left), (3, GAP), (0, RL - left)))
    return make_read(name, seq, s, ((0, left), (3, GAP), (0, 1), (1, 1), (0, RL - left - 2)))


def test_spliced_reads_count_and_exon_start_reads_do_not(tmp_path):
    """Spliced REF and ALT reads count; reads starting at the exon start (REF with its
    four A's, the ALT with five) read from inside the run and count neither."""
    reads = [_spliced(f"sR{i}", MRNA_REF, I0 - 40 + 3 * i) for i in range(6)]
    reads += [_spliced(f"sA{i}", MRNA_ALT, I0 - 40 + 3 * i) for i in range(6)]
    for i in range(4):
        reads.append(make_read(f"uR{i}", MRNA_REF[I0 : I0 + RL], E, ((0, RL),)))
        reads.append(make_read(f"uA{i}", MRNA_ALT[I0 : I0 + RL], E, ((0, 1), (1, 1), (0, RL - 2))))
    rd, ad, _ = _count(tmp_path, reads)
    assert (rd, ad) == (6, 6)


def test_a_short_exon_one_overhang_written_into_the_intron_is_not_alt(tmp_path):
    """REF mRNA reads starting one to three bases before the junction, the overhang
    written into the intron as mismatches: the A's they read past the exon edge are the
    run's, from exon 1. Not ALT (one base of overhang reads five A's, the ALT's count on
    exon 2 alone)."""
    reads = [_spliced(f"sR{i}", MRNA_REF, I0 - 40 + 5 * i) for i in range(6)]
    for ov in (1, 2, 3):
        for i in range(2):
            reads.append(
                make_read(f"o{ov}_{i}", MRNA_REF[I0 - ov : I0 - ov + RL], E - ov, ((0, RL),))
            )
    rd, ad, _ = _count(tmp_path, reads)
    assert (rd, ad) == (6, 0)


def test_a_pre_mrna_read_through_the_intron_still_reads_ref(tmp_path):
    """Guard: pre-mRNA REF reads through the intron read its AG before the run, which
    the longer allele does not have there: REF."""
    reads = [_spliced(f"sR{i}", MRNA_REF, I0 - 40 + 3 * i) for i in range(6)]
    for i in range(4):
        s = E - 30 - i
        reads.append(make_read(f"p{i}", CONTIG[s : s + RL], s, ((0, RL),)))
    rd, ad, _ = _count(tmp_path, reads)
    assert (rd, ad) == (10, 0)
