"""The sequence judge at an insertion of 50 bp or more: a change the read cannot place is neither
carried nor barred (RJ-23, amended 2026-11-03).

Mechanism. A long insertion is read over equal windows from its anchor, and the guard lays the
recurring changes beside it on every candidate haplotype as copies of a zone. A read carries a
listed change only where its own bases span the change's place. A change the read cannot place
(the read ends inside its repeat, or before it) is a plain alternative: no mask moves for it, for
any state, in either pass. Only changes the read does span can bar another copy of the zone. A
read that spans none is scored on plain graphs, as on develop.

The shapes are made up and small. Each is checked against the counts of develop, which the tests
assert: the previous fix (a read carries a change only where it spans its place) put the molecules
below on the wrong side of a bar and changed their calls.

Every count is checked for the four invariants and for binning invariance.
"""

import collections
import re

from helpers import count_bam_checked, make_read, write_contig

from gbcms import _rs

ARGS = (20, 20, True, True, True, False, False, False, 1)
CALL = re.compile(
    r"read call \d+:(\d+) (\S+)>(\S+) read=(\S+) mate=(\d) ref=(\w+) alt=(\w+) .*?partial=(\d+)"
)


def _molecule(contig, inserted, snps=()):
    """[(base, reference position | None)] of a molecule: `inserted` maps a reference position to
    the bases inserted after it; `snps` are positions that carry another base."""
    out = []
    for pos, base in enumerate(contig):
        if pos in snps:
            base = next(b for b in "ACGT" if b != base)
        out.append((base, pos))
        out.extend((b, None) for b in inserted.get(pos, ""))
    return out


def _read(name, molecule, start, length):
    """A read of `length` bases of `molecule` from reference position `start`, with the CIGAR the
    molecule gives it."""
    first = next(k for k, (_, pos) in enumerate(molecule) if pos == start)
    window = molecule[first : first + length]
    assert len(window) == length
    cigar, prev = [], None
    for _, pos in window:
        if pos is None:
            op = 1
        else:
            if prev is not None and pos > prev + 1:
                cigar.append([2, pos - prev - 1])
            op, prev = 0, pos
        if cigar and cigar[-1][0] == op:
            cigar[-1][1] += 1
        else:
            cigar.append([op, 1])
    seq = "".join(b for b, _ in window)
    return make_read(name, seq, start, tuple((op, n) for op, n in cigar))


def _row(by_row, alt_length):
    """The calls at the one row whose ALT is `alt_length` long."""
    (calls,) = (v for (_, alt), v in by_row.items() if len(alt) == alt_length)
    return calls


def _calls(tmp_path, caplog, contig, rows, reads, label):
    """Count `rows` (position, REF, ALT) on `reads` with the sequence judge. Returns the counts
    and, for each row (keyed by its position and ALT), the call of every read: {(name, mate): call}.
    """
    fa, bam = write_contig(tmp_path, contig, reads, label)
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
            pos, _, alt, _, mate, is_ref, is_alt, partial = g.groups()
            call = (
                "ALT"
                if is_alt == "true"
                else "REF" if is_ref == "true" else "partial" if partial != "0" else "depth"
            )
            seen[(int(pos) - 1, alt, m.split(" read=")[1].split()[0], mate)] = call
    by_row = collections.defaultdict(dict)
    for (pos, alt, name, mate), call in seen.items():
        by_row[(pos, alt)][(name, mate)] = call
    return counts, by_row


# ── A 60 bp insertion beside a repeat of four 12-base units, and molecules that carry the repeat's
# own unit insertions: they never carry the listed insertion, and cannot place it ───────────────

LEFT = (
    "CAGCTGAGGCTTTTTTACTAATTGTGAACACACTGTTGAGTATCGATAGCTTTTAAGGATACTAGATCTCAAAGGCTGGATTGCACTCC"
    "GGTAAAATCGACATTCTTAGGAGTGTTCCCTAATCGATGACAAGGTGGTCCGAAACAGTGTTGACGGATTATTACCGGGTGGTGTAGAA"
    "TCTAAAGATGTTCCCACTCTAT"
)
UNITS = ["CTATAGTGGTCG", "CTATAGTGGTCG", "CTATAGTGATCG", "CTATAGTGGTCA"]
RIGHT = (
    "CAACATAATGATTATCGTCCACCAATACAAATTAGCCCTACCTGCTTGCGGAAATTTGGAACGGAGATCCGCGTCCCGATACATGCATG"
    "ACGATGGGTTAACTGTCTACTCGCGGACTTGGGTGAGTCATAGTCAACGCTATATGCCCGAATAAGGAAAGATTGGTTTGGACGAGGAA"
    "TTTACCGACTGGAATTGGCAGGAGGAAGACCTGGTCGCAGGTGATCAAGCAAGCCTGTACCCGTGTTAGAGTATGGATAAGGGGTCCAC"
    "GGGTGTTATACAGGCTACATAAGATGCTGACCGAATCACTGGCCTGAAGGAGCCGGGATAAGGAACCATCGGGTTGGGTGAGACCGGGA"
    "CGGCCCTAGGATAAGGTAAACCGTTCATACTGCGGGTTATCTGT"
)
X_INSERT = "CCTTTCCATACATCTCCTACTGAAAGCTACCGTTCCCAGGCCAGATGCTTACCATCTAGT"
ANCHOR = len(LEFT)  # the base before the repeat (0-based)
UNIT = 12
ANCHOR_Z = 200  # the 90 bp insertion's anchor in CONTIG_Z


def _unit_shape():
    """The contig and the rows: the 60 bp insertion at the repeat's start, then its listed
    neighbours: three insertions of a unit (at the units' first bases), a deletion of 3 bases
    and an insertion of 6 beyond the repeat."""
    contig = LEFT + "G" + "".join(UNITS) + RIGHT
    assert contig[ANCHOR] == "G" and len(X_INSERT) == 60
    rows = [(ANCHOR, "G", "G" + X_INSERT)]
    for k in (1, 2, 3):
        pos = ANCHOR + k * UNIT + 1
        rows.append((pos, contig[pos], contig[pos] + contig[pos + 1 : pos + 1 + UNIT]))
    rows.append((ANCHOR + 55, contig[ANCHOR + 55 : ANCHOR + 59], contig[ANCHOR + 55]))
    rows.append(
        (ANCHOR + 95, contig[ANCHOR + 95], contig[ANCHOR + 95] + contig[ANCHOR + 96 : ANCHOR + 102])
    )
    return contig, rows


def _unit_reads(contig):
    """Molecules of the reference with two bases changed in the repeat (refs), of the repeat with
    its last unit written twice (nbb, refna), and of the repeat with two units written twice more
    (refnb). The insertions are the repeat's own units, so a read that ends inside the repeat
    cannot say which copy of the repeat it lies on."""
    last = ANCHOR + 1 + 3 * UNIT  # the last unit's first base
    first_twice = ANCHOR + 1 + UNIT
    ref = _molecule(contig, {}, (ANCHOR + 14, ANCHOR + 21))
    one = _molecule(contig, {last: "TATAGTGGTCAC"})
    two = _molecule(contig, {first_twice: "TATAGTGGTCGCTATAGTGATCGC", last: "TATAGTGGTCAC"})
    return [
        _read("refs11x151", ref, ANCHOR - 29, 151),
        _read("refs12x100", ref, ANCHOR - 18, 100),
        _read("refs12x151", ref, ANCHOR - 18, 151),
        _read("nbb0x151", one, ANCHOR - 100, 151),
        _read("nbb3x100", one, ANCHOR - 49, 100),
        _read("refna4x151", one, ANCHOR - 98, 151),
        _read("refna8x100", one, ANCHOR - 46, 100),
        _read("refnb7x151", two, ANCHOR - 59, 151),
        _read("refnb8x151", two, ANCHOR - 46, 151),
        _read("refnb9x151", two, ANCHOR - 33, 151),
    ]


def test_a_read_that_cannot_place_a_neighbour_stays_ref_where_develop_calls_it_ref(
    tmp_path, caplog
):
    """Seven reads of molecules that never carry the 60 bp insertion are REF at its row on develop,
    among them every read that writes the repeat's last unit twice: that change is one the read
    cannot place (the read ends inside the repeat or a unit short of it), so it is a plain
    alternative. Four of the seven went to depth when such a change was treated as carried for
    every read (a REF-reading read must stay REF). The three reads that write units in cis twice
    are partial on develop and stay so; none is ALT."""
    contig, rows = _unit_shape()
    counts, by_row = _calls(tmp_path, caplog, contig, rows, _unit_reads(contig), "unit")
    row = _row(by_row, len(rows[0][2]))
    assert collections.Counter(row.values()) == {"REF": 7, "partial": 3}, row
    assert sorted(k for k, v in row.items() if v == "partial") == [
        ("refnb7x151", "2"),
        ("refnb8x151", "2"),
        ("refnb9x151", "2"),
    ]
    for name in ("nbb0x151", "nbb3x100", "refna4x151", "refna8x100"):
        assert row[(name, "2")] == "REF", (name, row)
    assert (counts[0].rd, counts[0].ad) == (7, 0)


# ── A 90 bp insertion and a 150 bp insertion seven bases after it: the shorter one's inserted bases
# are the first 90 of the longer's. A read of a molecule with the longer insertion is not a carrier
# of the shorter ──────────────────────────────────────────────────────────────────────────────────

CONTIG_Z = (
    "TGGCATAGGGGGGTAAGTGTTTGTGTGTCTTTTGGACCCCCGGGCTAGAGTGGGACGGGGAACAACCTAAGGCCGCGATAATGTCCAT"
    "ATTGTGAGCTCCAGCAGACGGCATTCGGTATGTGCACGTCCAACGATTTTCCGGCGAAGTGTCACGGTAATTATGCTATGCCGCCCAG"
    "CGCCCCATGAGAGCCCTGGATCAAGGTGGTGTCCTAGATTTCATTTTAGGATCTGGTGGTGTGCTAGATTTCATTTTATGATCTGGTG"
    "GTGTGCTAGATTTCATTTTATGATCTGGTGGTGTGCTAGATTTCATTTTCTGATCTGGTGGTGTGCTAGATTTCATTTTCTGATCTGG"
    "TGGTGTGCTAGATTTCATTTTATGATCTGAACCTGGGTAGAATGACGACGTTAGTTAGGTAAAAGTCTCCCTGACCCCGATCACGCGG"
    "AACCAGCAAAAAATTAGCTGCTCGGAATCGTATCAACCTGGACTATTCGATGCACCCTCACGGTCAGGGTTTAGGATCCTCCAGCGCA"
    "TTTCAATAGTTCTATTACGCGTAGGCTAGTTCATAACAGACAGCAATTTATTCGCCACTGCGTGGTAGTACTCTAGATTAGCACGGGT"
    "AATGCCCGGATTGGAACGCACAGAAGGCTCGGCG"
)
Z_INSERT = (
    "AATGACGCCATCACGCTTTACGGGGCCAGGTCCACCGCATTAACCTAAATTCGAACCTCCCCAGTCCGGAGCCCGAAGCGGGTTGCTG"
    "GGTACATGTGGGGCAAGACTGTCCGGTATTGTGAGCTACTAACCGTTTGCGCACACATACCT"
)
X_LENGTH = 90  # the shorter insertion's inserted bases: Z_INSERT[:X_LENGTH]
Z_OFFSET = 7  # the longer insertion's anchor is this far after the shorter one's


def _two_insertions_shape():
    """The contig, the two rows (the 90 bp insertion, the 150 bp insertion) and the reads: three of
    the shorter insertion's molecule, which has two bases changed beside its anchor (the changes
    recur, so the guard lays them as masks), and one of the longer insertion's molecule, which
    ends 12 bases into the insertion."""
    contig = CONTIG_Z
    assert contig[ANCHOR_Z] == "G" and contig[ANCHOR_Z + Z_OFFSET] == "T"
    rows = [
        (ANCHOR_Z, "G", "G" + Z_INSERT[:X_LENGTH]),
        (ANCHOR_Z + Z_OFFSET, "T", "T" + Z_INSERT),
    ]
    shorter = _molecule(contig, {ANCHOR_Z: Z_INSERT[:X_LENGTH]}, (ANCHOR_Z + 6, ANCHOR_Z + 13))
    longer = _molecule(contig, {ANCHOR_Z + Z_OFFSET: Z_INSERT})
    reads = [
        _read(f"short{k}", shorter, ANCHOR_Z - back, 151) for k, back in enumerate((40, 29, 18))
    ]
    reads.append(_read("longer", longer, ANCHOR_Z - 80, 100))
    return contig, rows, reads


def test_a_read_of_a_neighbouring_longer_insertion_is_not_alt_at_the_shorter_row(tmp_path, caplog):
    """The read of the 150 bp insertion's molecule ends 12 bases into the insertion. It is REF at
    the 90 bp insertion's row on develop: it does not carry that insertion, and a read is never
    credited to a row's ALT unless its own bases carry the ALT. It was ALT there when a change was
    treated as carried for every read. At its own row it is ALT, and the three reads of the shorter
    insertion are ALT at theirs and REF here."""
    contig, rows, reads = _two_insertions_shape()
    counts, by_row = _calls(tmp_path, caplog, contig, rows, reads, "two")
    shorter_row, longer_row = _row(by_row, len(rows[0][2])), _row(by_row, len(rows[1][2]))
    assert shorter_row == {
        ("short0", "2"): "ALT",
        ("short1", "2"): "ALT",
        ("short2", "2"): "ALT",
        ("longer", "2"): "REF",
    }, shorter_row
    assert longer_row == {
        ("short0", "2"): "REF",
        ("short1", "2"): "REF",
        ("short2", "2"): "REF",
        ("longer", "2"): "ALT",
    }, longer_row
    assert (counts[0].rd, counts[0].ad, counts[1].rd, counts[1].ad) == (1, 3, 3, 1)
