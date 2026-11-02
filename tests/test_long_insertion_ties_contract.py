"""The sequence judge at an insertion of 50 bp or more: a read that ties the ALT and another
allele holds the ALT's equal windows too (RJ-23, amended 2026-11-02).

Mechanism. A REF or ALT call holds its allele's equal window from the anchor (the base before
the insert and as many bases on as REF's window holds). A read inside the insert, or holding
only its far junction, holds no window and is depth only, whichever allele it fits. The first
build checked the allele the call names. A read that ties the ALT and another allele (REF
beaten) is partial evidence of that allele, and that allele's haplotype, which has no long
insert, looks like REF's: the read held its window by its alignment there, though its alignment
to the ALT lay inside the insert. HG002's chr18:74617096 (a 1,997 bp insertion that begins with
a C and then the reference's own bases, beside a found 1 bp insertion A>AC) and chr4:660975
each gave partial evidence from such reads.

Shape (synthetic, 8-base repeat unit TGACTGTC, the tract 17 units long after the anchor A at
0-based 500): the ALT inserts C and the tract's first 60 bases, so it begins as the 1 bp
insertion A>AC does, and the tract's bases 58 and 59 (A, C) are the anchor and the inserted
base: a read starting at the ALT's inserted bases 59 and 60 and running on into the tract
reads "AC" and the tract from its start, as the 1 bp insertion's carrier does. It fits both
alleles exactly; REF is beaten by the insertion. The windows are the anchor and the inserted
C. The 1 bp insertion is listed as its own input row.
"""

import collections
import random
import re

from helpers import count_bam_checked, make_read, write_contig

from gbcms import _rs

ENGINE_LOGGER = "_rs.counting.engine"
ARGS = (20, 20, True, True, True, False, False, False, 1)
CALL = re.compile(
    r"read call \S+ \S+>(\S+) read=([a-z]+)\d+ mate=(\d) ref=(\w+) alt=(\w+) .*?partial=(\d+)"
)
A = 500  # the anchor (0-based), the contig's base 'A'
UNIT = "TGACTGTC"  # the tract's repeat unit; its bases 2 and 3 are 'A', 'C'
UNITS = 17
M = 60  # bases of the tract the ALT's insert copies after its first base, C
N = 100  # read length


def _shape():
    """The contig, the long insertion's ALT, the 1 bp insertion's ALT, and the reads."""
    rng = random.Random(2641)
    left = "".join(rng.choice("ACGT") for _ in range(A)) + "A"
    right = "".join(rng.choice("ACGT") for _ in range(500))
    tract = UNIT * UNITS
    contig = left + tract + right
    insert = "C" + tract[:M]
    assert tract[M - 2 : M] == "AC" and tract[0] != "C" and tract[M - 1] != "A"
    alt = "A" + insert
    one = "AC"
    hap_alt = contig[: A + 1] + insert + contig[A + 1 :]
    hap_one = contig[: A + 1] + "C" + contig[A + 1 :]
    reads = []

    def add(name, seq, start, cigar):
        reads.append(make_read(name, seq, start, cigar))

    for i in range(6):  # REF molecules over the anchor
        s = 440 + 4 * i
        add(f"ref{i}", contig[s : s + N], s, ((0, N),))
    for i in range(6):  # carriers of the long insertion written with it (flank, the whole insert)
        s = 470 + 3 * i
        m1 = A + 1 - s
        add(
            f"thr{i}", hap_alt[s : s + N], s, ((0, m1), (1, len(insert)), (0, N - m1 - len(insert)))
        )
    # fragments of the insert with its far junction: the ALT's inserted bases 59 and 60 and the
    # tract from its start, which is also what a carrier of the 1 bp insertion reads from the
    # anchor on; aligned as that insertion
    for i in range(4):
        far = hap_alt[A + 1 + 59 : A + 1 + 59 + N]  # the anchor-like "AC", then the tract
        assert far[:2] == "AC" and far[2 : 2 + 40] == tract[:40]
        add(f"far{i}", far, A, ((0, 1), (1, 1), (0, N - 2)))
    for i in range(
        4
    ):  # carriers of the 1 bp insertion from the flank on (they fit the ALT too, to its window)
        s = 445 + 3 * i
        m1 = A + 1 - s
        add(f"one{i}", hap_one[s : s + N], s, ((0, m1), (1, 1), (0, N - m1 - 1)))
    return contig, alt, one, reads


def _run(tmp_path, caplog):
    """Counts the long insertion and the 1 bp insertion (two rows of one call) under the
    sequence judge, asserting the counting invariants on every row. Returns the counts and the
    tally of read calls by (kind of read, call) per row ALT."""
    contig, alt, one, reads = _shape()
    fa, bam = write_contig(tmp_path, contig, reads, "tie")
    alts = [alt, one]
    pvs = _rs.prepare_variants(
        [_rs.Variant("1", A, contig[A], a, "INDEL") for a in alts],
        fa,
        5,
        False,
        1,
        True,
        read_judge="sequence",
    )
    caplog.clear()
    with caplog.at_level(5, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        counts = count_bam_checked(
            bam,
            [pv.variant for pv in pvs],
            [None] * len(pvs),
            *ARGS,
            read_judge="sequence",
            reference_fasta=fa,
        )
    for c in counts:
        assert c.dp >= c.rd + c.ad and c.dpf >= c.rdf + c.adf
        assert c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    tallies = {a: collections.Counter() for a in alts}
    seen = {}
    for record in caplog.records:
        if m := CALL.search(record.getMessage()):
            row, kind, mate, ref, alt_, partial = m.groups()
            call = (
                "ALT"
                if alt_ == "true"
                else "REF" if ref == "true" else "partial" if partial != "0" else "depth"
            )
            seen[(row, record.getMessage().split(" read=")[1].split()[0], mate)] = (kind, call)
    for (row, _, _), kind_call in seen.items():
        tallies[row][kind_call] += 1
    return counts, tallies, alt


def test_the_shape_is_a_tie_of_the_long_insertion_and_the_one_base_insertion():
    """Guard on the shape itself: the far-junction read is the 1 bp insertion's carrier from the
    anchor on, and the long insertion's inserted bases 59 and 60 read the anchor and C."""
    contig, alt, one, reads = _shape()
    tract = UNIT * UNITS
    insert = alt[1:]
    assert alt[0] == contig[A] == "A" and insert == "C" + tract[:M] and len(insert) == M + 1
    assert insert[M - 1 : M + 1] == "AC"
    assert contig[A + 1 : A + 1 + len(tract)] == tract


def test_a_read_tying_the_long_insertion_with_another_allele_holds_the_windows_of_both(
    tmp_path, caplog
):
    """The far-junction reads fit the long insertion and the 1 bp insertion exactly, REF by a
    gap. Their alignment to the ALT lies inside the insert, which holds no window from the
    anchor, so they are depth only on the long insertion's row, not partial evidence of the 1 bp
    insertion; no read adds an ALT count. Carriers written with the insert are ALT; REF molecules
    are REF; the 1 bp insertion's carriers from the flank fit both alleles over the windows and
    stay partial evidence."""
    counts, tallies, alt = _run(tmp_path, caplog)
    long_row = tallies[alt]
    assert long_row[("thr", "ALT")] == 6, long_row
    assert long_row[("ref", "REF")] == 6, long_row
    assert long_row[("far", "depth")] == 4 and long_row[("far", "partial")] == 0, long_row
    assert long_row[("one", "partial")] == 4, long_row
    long_counts = counts[0]
    assert (long_counts.rd, long_counts.ad, long_counts.partial_alt) == (6, 6, 4)
