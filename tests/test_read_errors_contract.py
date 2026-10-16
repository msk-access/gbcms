"""One substitution in a read is an error, on REF and ALT reads alike (operator,
2026-10-15, #278).

A read whose bases fit an allele but for one confident mismatch, where the other
allele mismatches the same read base too, decides nothing by that base: it counts for
the allele its other bases show. Inside the event the mismatch is no error when it
recurs: carried by two or more fragments and 5% of those fitting the allele, it is
another allele at the locus (RJ-25), so the reads carrying it are neither (partial
evidence when closer to the ALT). In the flank a shared substitution is a phased
neighbour (RJ-23) and the reads keep their call. Past the shorter allele's window the
other allele's base is the reference's next, which both alleles continue with, so a
deletion row's REF reads and an insertion row's ALT reads tolerate an error there as
the others do.
"""

import random

from helpers import count_checked, make_read, write_contig

from gbcms import _rs as gbcms_rs

READ = 100
P = 300  # the anchor


def _contig() -> str:
    """Random C/G/T sequence: an A anywhere is a substitution no flank repeats."""
    rng = random.Random(4)
    seq = "".join(rng.choice("CGT") for _ in range(800))
    # the deletion row's six bases neither slide left nor right
    assert seq[P + 1] != seq[P + 7] and seq[P] != seq[P + 6]
    return seq


CONTIG = _contig()
INSERT = "ACGTTGCA"  # its A ends cannot slide onto C/G/T flanks
DELETED = 6


def _prepared(fa, ref, alt):
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", P, ref, alt, "X")], fa, 5, False, 1, True
    )
    assert pv.gbcms_status == "PASS", pv.gbcms_status_reason
    return pv.variant


def _reads(hap: str, tag: str, n: int, first: int = 0, sub: int | None = None, cigar=None):
    """`n` reads of `hap` from P - 40 - i, each its own fragment; `sub` (an offset into
    the read's haplotype) read as A; `cigar(left)` writes its ops."""
    out = []
    for i in range(first, first + n):
        s = P - 40 - i
        seq = list(hap[s : s + READ])
        if sub is not None:
            seq[sub - s] = "A"
        left = P + 1 - s
        out.append(make_read(f"{tag}{i}", "".join(seq), s, cigar(left) if cigar else ((0, READ),)))
    return out


def _insertion(tmp_path, reads):
    fa, bam = write_contig(tmp_path, CONTIG, reads, "ins")
    (c,) = count_checked(bam, [_prepared(fa, CONTIG[P], CONTIG[P] + INSERT)])
    _invariants(c)
    return c.rd, c.ad, c.partial_alt


ALT_HAP = CONTIG[: P + 1] + INSERT + CONTIG[P + 1 :]


def _carriers(n, first=0, sub=None):
    return _reads(
        ALT_HAP,
        "a",
        n,
        first,
        sub,
        lambda left: ((0, left), (1, len(INSERT)), (0, READ - left - len(INSERT))),
    )


WRONG = P + 1 + 3  # the insert's fourth base (T), on the ALT's haplotype


def _invariants(c):
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev


def test_one_substitution_on_a_carrier_is_an_error(tmp_path):
    """Beside six carriers, one more reads the base before the anchor as A, where
    REF's window mismatches it too: an error, ALT."""
    reads = _reads(CONTIG, "r", 6) + _carriers(6) + _carriers(1, first=6, sub=P - 1)
    assert _insertion(tmp_path, reads) == (6, 7, 0)


def test_a_flank_substitution_two_fragments_share_keeps_their_call(tmp_path):
    """Two carriers read the same flank base as A: a phased neighbour beside the event
    (too few to mask), not another allele of the event: ALT, as GATK, bcftools,
    freebayes and bam-readcount count them."""
    reads = _reads(CONTIG, "r", 6) + _carriers(6) + _carriers(2, first=6, sub=P - 1)
    assert _insertion(tmp_path, reads) == (6, 8, 0)


def test_one_inserted_base_read_wrong_is_an_error(tmp_path):
    """One carrier reads an inserted base as A (RJ-21's error case): ALT."""
    reads = _reads(CONTIG, "r", 6) + _carriers(6) + _carriers(1, first=6, sub=WRONG)
    assert _insertion(tmp_path, reads) == (6, 7, 0)


def test_an_inserted_base_two_fragments_read_alike_is_another_allele(tmp_path):
    """Two carriers read the same inserted base as A: inside the event that is another
    allele at the locus, neither, closer to the ALT (partial evidence)."""
    reads = _reads(CONTIG, "r", 6) + _carriers(6) + _carriers(2, first=6, sub=WRONG)
    assert _insertion(tmp_path, reads) == (6, 6, 2)


def test_one_substitution_on_a_ref_read_is_an_error(tmp_path):
    """One REF read reads the base before the anchor as A: an error, REF."""
    reads = _reads(CONTIG, "r", 6) + _reads(CONTIG, "e", 1, first=6, sub=P - 1) + _carriers(6)
    assert _insertion(tmp_path, reads) == (7, 6, 0)


def test_an_error_past_the_shorter_window_is_tolerated_on_ref_reads(tmp_path):
    """A six-base deletion: REF's window runs six bases past the ALT's, whose reading
    ends first. A REF read with an error on REF's window's last base is an error as one
    on a base both windows read is: the ALT continues with the reference there too."""
    alt_hap = CONTIG[: P + 1] + CONTIG[P + 1 + DELETED :]
    carriers = _reads(
        alt_hap, "a", 6, cigar=lambda left: ((0, left), (2, DELETED), (0, READ - left))
    )
    erred = _reads(CONTIG, "e", 1, first=6, sub=P + DELETED + 2)
    reads = _reads(CONTIG, "r", 6) + erred + carriers
    fa, bam = write_contig(tmp_path, CONTIG, reads, "del")
    (c,) = count_checked(bam, [_prepared(fa, CONTIG[P : P + 1 + DELETED], CONTIG[P])])
    _invariants(c)
    assert (c.rd, c.ad, c.partial_alt) == (7, 6, 0)
