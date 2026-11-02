"""The sequence judge at an insertion of 50 bp or more: the windows read on to the base that tells
them apart are counted in the columns of the path a read lies on (RJ-23, RJ-27, amended
2026-11-02).

Mechanism. A pure insertion of 50 bp or more is read over equal windows from its anchor. Where
the guard's masks hide the base that tells REF's window from X's, both windows read on to the next
base that does, and a read short of it is depth only. The length they read on to was measured on
REF's and X's haplotypes without the listed changes beside X, while the masks are laid by the
columns of the path a read lies on. A deletion listed between the anchor and that base takes
bases out of the haplotype of a read that carries it, so the base is nearer by as many columns
there: a REF molecule carrying the deletion that reached the base was short of a window measured
without it (depth only), while one without the deletion reaching it was REF.

Shape: a 60 bp insertion of new sequence after the anchor at 0-based 500 (a run of none, so the
windows hold the anchor and one base). A germline change on half the molecules at each of 501,
502, 506 and 507, which the guard masks, and X's inserted bases at the other columns up to 508 are
REF's, so no base tells the windows apart until 509: both read on to it, 10 columns from the anchor
with the anchor. A 2 bp deletion listed after 503 (504 and 505) takes two bases out of the
haplotype of a read carrying it: 509 is two columns nearer there, and the windows hold 8.
"""

import collections
import random
import re

from helpers import count_bam_checked, make_read, write_contig

from gbcms import _rs

ARGS = (20, 20, True, True, True, False, False, False, 1)
CALL = re.compile(
    r"read call \S+ (\S+)>(\S+) read=([a-z_]+)\d+ mate=(\d) ref=(\w+) alt=(\w+) .*?partial=(\d+)"
)
A = 500  # the anchor (0-based)
N = 60  # the insertion's length
DEL = A + 3  # the listed deletion removes the two bases after this one (504 and 505)
MASKED = (A + 1, A + 2, A + 6, A + 7)  # germline changes the guard masks
STOP = A + 9  # the base that tells the windows apart


def _shape():
    """The contig, X's ALT, and the reads."""
    rng = random.Random(2647)
    contig = "".join(rng.choice("ACGT") for _ in range(1200))

    def differs(b):
        return next(t for t in "ACGT" if t != b)

    ins = [rng.choice("ACGT") for _ in range(N)]
    # X's inserted base i sits at column i + 1 of its haplotype: REF's base there where the
    # guard does not mask it (so it tells nothing), and not REF's at 509
    for i, at in ((2, A + 3), (3, A + 4), (4, A + 5), (5, A + 8), (7, A + 8)):
        ins[i] = contig[at]
    ins[6] = differs(contig[STOP])
    ins[8] = differs(contig[STOP])
    ins[0] = differs(contig[A + 1])  # no run: the insertion does not slide
    ins[N - 1] = differs(contig[A])
    ins = "".join(ins)
    snp = {p: differs(contig[p]) for p in MASKED}

    def haplotype(x: bool, deleted: bool, snps: bool) -> str:
        seq = list(contig)
        if snps:
            for p, b in snp.items():
                seq[p] = b
        seq = "".join(seq)
        if deleted:
            seq = seq[: DEL + 1] + seq[DEL + 3 :]
        if x:
            seq = seq[: A + 1] + ins + seq[A + 1 :]
        return seq

    reads = []
    # REF molecules over the anchor and the masked bases, with the germline changes or without
    for i in range(6):
        s = 440 + 3 * i
        reads.append(
            make_read(f"ref{i}", haplotype(False, False, False)[s : s + 100], s, ((0, 100),))
        )
        reads.append(
            make_read(f"refsnp{i}", haplotype(False, False, True)[s : s + 100], s, ((0, 100),))
        )
    # carriers of the insertion, written with it
    for i in range(4):
        s = 470 + 3 * i
        m = A + 1 - s
        reads.append(
            make_read(
                f"thr{i}",
                haplotype(True, False, False)[s : s + 100],
                s,
                ((0, m), (1, N), (0, 100 - m - N)),
            )
        )
    # molecules ending at 509 (the base that tells) or one base short of it: REF molecules with
    # the deletion (written with it) and without
    deleted, plain = haplotype(False, True, False), haplotype(False, False, False)
    for i in range(3):
        s = 490 - i
        m = DEL + 1 - s  # the bases up to the deletion's anchor
        reads.append(make_read(f"delto{i}", deleted[s : STOP - 1], s, ((0, m), (2, 2), (0, 4))))
        reads.append(make_read(f"delshort{i}", deleted[s : STOP - 2], s, ((0, m), (2, 2), (0, 3))))
        reads.append(make_read(f"plainto{i}", plain[s : STOP + 1], s, ((0, STOP + 1 - s),)))
        reads.append(make_read(f"plainshort{i}", plain[s:STOP], s, ((0, STOP - s),)))
    return contig, contig[A] + ins, reads


def test_the_windows_read_on_are_counted_in_the_columns_of_the_path_a_read_lies_on(
    tmp_path, caplog
):
    """The REF molecules carrying the deletion that end at 509 hold the windows on their own
    haplotype and are REF, as the plain molecules ending there are; one base short of it, both are
    depth only. Carriers of the insertion are ALT; REF molecules with the germline changes or
    without them are REF. Counts are binning invariant."""
    contig, alt, reads = _shape()
    fa, bam = write_contig(tmp_path, contig, reads, "columns")
    variants = [
        _rs.Variant("1", A, contig[A], alt, "INDEL"),
        _rs.Variant("1", DEL, contig[DEL : DEL + 3], contig[DEL], "INDEL"),
    ]
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
    messages = [r.getMessage() for r in caplog.records]
    assert any("both read on, 2 -> 10 bases from the anchor" in m for m in messages)
    seen = {}  # the read calls are logged by both counts (production bins, one variant per bin)
    for m in messages:
        if g := CALL.search(m):
            _, alt_, kind, mate, is_ref, is_alt, partial = g.groups()
            call = (
                "ALT"
                if is_alt == "true"
                else "REF" if is_ref == "true" else "partial" if partial != "0" else "depth"
            )
            seen[(alt_, m.split(" read=")[1].split()[0], mate)] = (kind, call)
    calls = collections.defaultdict(collections.Counter)
    for (alt_, _, _), kind_call in seen.items():
        calls[alt_][kind_call] += 1
    row = calls[alt]
    assert row[("thr", "ALT")] == 4, row
    for kind in ("ref", "refsnp"):
        assert row[(kind, "REF")] == 6, (kind, row)
    assert row[("plainto", "REF")] == 3 and row[("plainshort", "depth")] == 3, row
    assert row[("delto", "REF")] == 3 and row[("delshort", "depth")] == 3, row
    x = counts[0]
    assert (x.rd, x.ad) == (18, 4), (x.rd, x.ad, x.partial_alt)
