"""The sequence judge at an insertion of 50 bp or more: a REF call that the reach rule made holds
REF's equal window alone (RJ-23, amended 2026-11-02; the class plan decides it for every shape).

Mechanism. A read that ties REF with another allele it never reaches (the two differ only past a
repeat's end and the read ends inside it, or the read writes an indel of its own across the
other allele) is REF: it shows REF where REF and X differ and nothing of the other allele. At a
pure insertion of 50 bp or more a REF call holds REF's equal window from the anchor. The tie
check of the previous fix asked the windows of every allele a call ties, and so also the other
allele's, measured on its own haplotype: where the other allele is a longer insertion, its
haplotype carries the insertion before the called insertion's anchor, and a read that lies
inside the repeat aligns to its copy, so it holds no window there. Those REF calls became depth
only. They stay as on develop: REF's window alone.

Shape (synthetic): a tandem repeat of eight copies of a 43-base unit, each copy with up to two
bases changed, so a read inside it fits the reference and a copy shifted by a unit alike (one
base apart). The input lists a 172 bp insertion (four units) and an 86 bp insertion (two units)
at the repeat's start. Reads of 100 and 151 bases sweep the start, some with a base changed.
"""

import collections
import random
import re

from helpers import count_bam_checked, make_read, write_contig

from gbcms import _rs

ARGS = (20, 20, True, True, True, False, False, False, 1)
CALL = re.compile(
    r"read call \S+ (\S+)>(\S+) read=(\S+) mate=(\d) ref=(\w+) alt=(\w+) .*?partial=(\d+)"
)
A = 300  # the repeat's first base (0-based)


def _shape():
    """The contig, the two insertions' (position, REF, ALT), and the reads."""
    rng = random.Random(17)
    unit_len = rng.choice([12, 18, 24, 30, 43])
    copies_n = rng.choice([6, 8, 10])
    unit = "".join(rng.choice("ACGT") for _ in range(unit_len))
    copies = []
    for _ in range(copies_n):
        c = list(unit)
        for _ in range(rng.choice([0, 1, 1, 2])):
            c[rng.randrange(unit_len)] = rng.choice("ACGT")
        copies.append("".join(c))

    def flank(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    contig = flank(A) + "".join(copies) + flank(900)
    anchor = A - 1 + rng.choice([0, 1])
    x_units, z_units = rng.choice([2, 3, 4]), rng.choice([1, 2])
    assert (unit_len, copies_n, anchor, x_units, z_units) == (43, 8, 300, 4, 2)
    long_ins = contig[anchor + 1 : anchor + 1 + x_units * unit_len]
    short_ins = contig[anchor + 1 : anchor + 1 + z_units * unit_len]
    rows = [
        (anchor, contig[anchor], contig[anchor] + long_ins),
        (anchor, contig[anchor], contig[anchor] + short_ins),
    ]
    if rng.random() < 0.5:  # the short one written a base to the left
        left = anchor - 1
        rows[1] = (
            left,
            contig[left],
            contig[left] + contig[left + 1 : left + 1 + z_units * unit_len],
        )
    reads = []
    for start in range(anchor - 120, anchor + 60, 7):
        for length in (100, 151):
            seq = list(contig[start : start + length])
            if rng.random() < 0.5:
                seq[rng.randrange(length)] = rng.choice("ACGT")
            reads.append(make_read(f"r{start}x{length}", "".join(seq), start, ((0, length),)))
    return contig, rows, reads


def test_a_ref_call_made_by_the_reach_rule_holds_refs_window_alone(tmp_path, caplog):
    """On the 86 bp insertion's row seven reads lie inside the repeat and tie REF with the 172 bp
    insertion, which they never reach: they are REF, as on develop (the previous fix made them
    depth only). No read is ALT or partial evidence, and the other row has no REF call. Counts
    are binning invariant."""
    contig, rows, reads = _shape()
    fa, bam = write_contig(tmp_path, contig, reads, "reach")
    variants = [_rs.Variant("1", pos, ref, alt, "INDEL") for pos, ref, alt in rows]
    pvs = _rs.prepare_variants(variants, fa, 5, False, 1, True, read_judge="sequence")
    caplog.clear()
    with caplog.at_level(5, logger="_rs"):
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
    seen = {}  # the read calls are logged by both counts (production bins, one variant per bin)
    for record in caplog.records:
        m = record.getMessage()
        if g := CALL.search(m):
            _, alt, _, mate, is_ref, is_alt, partial = g.groups()
            call = (
                "ALT"
                if is_alt == "true"
                else "REF" if is_ref == "true" else "partial" if partial != "0" else "depth"
            )
            seen[(alt, m.split(" read=")[1].split()[0], mate)] = call
    by_row = collections.defaultdict(dict)
    for (alt, name, mate), call in seen.items():
        by_row[len(alt)][(name, mate)] = call
    long_row, short_row = by_row[len(rows[0][2])], by_row[len(rows[1][2])]
    assert collections.Counter(long_row.values()) == {"depth": 34}, long_row
    assert collections.Counter(short_row.values()) == {"depth": 27, "REF": 7}, short_row
    assert sorted(k for k, v in short_row.items() if v == "REF") == EXPECTED_REF, short_row
    assert (counts[1].rd, counts[1].ad, counts[1].partial_alt) == (7, 0, 0)
    assert (counts[0].rd, counts[0].ad, counts[0].partial_alt) == (0, 0, 0)


# the reads that are REF on develop and stay so (name, mate): reads of 151 bases that start 43
# to 1 bases before the repeat
EXPECTED_REF = [(f"r{start}x151", "2") for start in range(257, 300, 7)]
