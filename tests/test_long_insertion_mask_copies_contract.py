"""The sequence judge at an insertion of 50 bp or more: a read lies on one copy of the masks'
zone for every candidate (RJ-23, RJ-27, amended 2026-11-02).

Mechanism. The guard's masks are laid by offset from the anchor on every candidate haplotype,
by the columns of the path a read lies on (a listed insertion or deletion between the anchor and
a masked base moves the offset on a read that carries it). The graph holds a copy of the span
the masks lie over per set of such listed changes, and a state's best path took whichever copy
masked its own mismatches: an X read ending in the insert, long before a deletion listed past
it, took the deletion's copy for X (masking a base where the read has a third base) and the other
copy for a listed sibling of X's length that differs from X at a masked offset. X scored 0, the
sibling -4, and the read was ALT from a deletion it cannot carry. A read now lies on the copy of
the changes it carries or spans, the same for every state.

Shape: a 60 bp insertion of new sequence after the anchor at 0-based 500 (a run of none), a
listed sibling that differs from it at inserted base 19, a 2 bp deletion listed after 505 (a
context), and a germline change at 520 (20 bases from the anchor, so inserted base 19) on half
the molecules, which the guard masks (it reads 25 bases beside the event). Through the deletion
the masked base is two nearer, at inserted base 17. Reads of X's haplotype end inside the insert
with a third base at inserted base 17 (written three ways: clipped, inserted, written through),
or carry X and the deletion.
"""

import collections
import random
import re

import pytest
from helpers import count_bam_checked, make_read, write_contig

from gbcms import _rs

ENGINE_LOGGER = "_rs.counting.engine"
ARGS = (20, 20, True, True, True, False, False, False, 1)
CALL = re.compile(
    r"read call \S+ (\S+)>(\S+) read=([a-z_]+)\d+ mate=(\d) ref=(\w+) alt=(\w+) .*?partial=(\d+)"
)
A = 500  # the anchor (0-based)
N = 60  # the insertion's length
DEL = A + 5  # the listed deletion removes the two bases after this one
SNP = A + 20  # the germline change the guard masks


def _shape(repeat: bool = False):
    """The contig, the rows' (REF, ALT) alleles (X, the sibling, the deletion), and the reads.
    With `repeat`, the six bases after the deletion's padding base are one base, a tract the
    deletion slides over."""
    rng = random.Random(2645)
    bases = [rng.choice("ACGT") for _ in range(1200)]
    if repeat:
        tract = next(b for b in "ACGT" if b != bases[DEL] and b != bases[DEL + 7])
        bases[DEL + 1 : DEL + 7] = [tract] * 6
    contig = "".join(bases)
    ins = "".join(rng.choice("ACGT") for _ in range(N))
    while ins[0] == contig[A + 1] or ins[-1] == contig[A]:
        ins = ins[1:] + ins[0]
    other = next(b for b in "ACGT" if b != ins[19])
    sibling = ins[:19] + other + ins[20:]
    snp_alt = next(b for b in "ACGT" if b != contig[SNP])
    reads = []

    def haplotype(x: bool, deleted: bool, snp: bool) -> str:
        seq = list(contig)
        if snp:
            seq[SNP] = snp_alt
        seq = "".join(seq)
        if deleted:
            seq = seq[: DEL + 1] + seq[DEL + 3 :]
        if x:
            seq = seq[: A + 1] + ins + seq[A + 1 :]
        return seq

    def add(kind, i, hap, start, cigar, length):
        reads.append(make_read(f"{kind}{i}", hap[start : start + length], start, cigar))

    # REF molecules over the anchor and past the germline change, with it or without
    for i in range(4):
        start = 466 + 3 * i
        add("ref", i, haplotype(False, False, False), start, ((0, 100),), 100)
        add("refsnp", i, haplotype(False, False, True), start, ((0, 100),), 100)
    hx = haplotype(True, False, False)
    hxd = haplotype(True, True, False)
    third = next(b for b in "ACGT" if b != hx[A + 1 + 17])
    # X reads ending inside the insert, a third base at inserted base 17; the aligner clips the
    # insert, writes it, or writes the bases through
    for i in range(2):
        start = 490 - i
        m = A + 1 - start  # the bases up to the anchor
        seq = hx[start : start + 45 + i]
        seq = seq[: A + 1 + 17 - start] + third + seq[A + 1 + 17 - start + 1 :]
        for kind, cigar in (
            ("clipped", ((0, m), (4, len(seq) - m))),
            ("inserted", ((0, m), (1, len(seq) - m))),
            ("through", ((0, len(seq)),)),
        ):
            reads.append(make_read(f"{kind}_x{i}", seq, start, cigar))
    # carriers of X and the deletion, written with both
    for i in range(4):
        start = 486 + i
        m = A + 1 - start
        before = DEL + 1 - (A + 1)  # the bases between the anchor and the deletion's
        rest = 100 - m - N - before
        add("xdel", i, hxd, start, ((0, m), (1, N), (0, before), (2, 2), (0, rest)), 100)
    return contig, ins, sibling, reads


def test_a_read_lies_on_the_masks_of_the_changes_it_carries_or_spans(tmp_path, caplog):
    """X reads ending in the insert, with a third base at the base the deletion's copy masks, are
    partial evidence of X and the sibling, whichever way the aligner wrote them (they were ALT:
    X took the deletion's copy, the sibling the other). Carriers of X and the deletion are told X
    from the sibling at inserted base 19, which their copy leaves unmasked: ALT. REF molecules
    are REF, with the germline change or without it. Counts are binning invariant."""
    contig, ins, sibling, reads = _shape()
    fa, bam = write_contig(tmp_path, contig, reads, "copies")
    x_alt = contig[A] + ins
    variants = [
        _rs.Variant("1", A, contig[A], x_alt, "INDEL"),
        _rs.Variant("1", A, contig[A], contig[A] + sibling, "INDEL"),
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
    laid = [
        m for m in messages if "the guard's masks laid on every candidate haplotype by offset" in m
    ]
    assert laid and all("through 1 context changes" in m for m in laid), laid
    seen = {}  # the read calls are logged by both counts (production bins, one variant per bin)
    for m in messages:
        if g := CALL.search(m):
            ref_, alt_, kind, mate, is_ref, is_alt, partial = g.groups()
            call = (
                "ALT"
                if is_alt == "true"
                else "REF" if is_ref == "true" else "partial" if partial != "0" else "depth"
            )
            seen[(alt_, m.split(" read=")[1].split()[0], mate)] = (kind, call)
    calls = collections.defaultdict(collections.Counter)
    for (alt_, _, _), kind_call in seen.items():
        calls[alt_][kind_call] += 1
    row = calls[x_alt]
    for kind in ("clipped_x", "inserted_x", "through_x"):
        assert row[(kind, "ALT")] == 0, (kind, row)
    assert row[("xdel", "ALT")] == 4, row
    for kind in ("ref", "refsnp"):
        assert row[(kind, "REF")] == 4, (kind, row)
    x = counts[0]
    assert (x.rd, x.ad) == (8, 4), (x.rd, x.ad, x.partial_alt)


@pytest.mark.parametrize("repeat", [False, True], ids=["unique", "repeat"])
def test_a_read_touching_a_listed_deletions_place_carries_nothing(tmp_path, caplog, repeat):
    """The copy laid for the deletion cuts the reference under it and holds its padding base, so
    a read of X ending on that base (or, with a tract after it, two bases into the tract) went
    through the deletion's segment at no cost and was held to its copy, which masks inserted base
    17: a sequencing error there was hidden, X scored 0 and the sibling -4, and the read was ALT
    (3 of 3 in both shapes). A read carries a change only where its bases span the change's
    place (the anchor to the base past the repeat it slides over), so these reads lie on the
    masks without it, and the error weighs on X and the sibling alike: partial evidence, with
    the error or without it. A sequencing error never raises a read's class."""
    contig, ins, sibling, shape_reads = _shape(repeat)
    reads = [r for r in shape_reads if r.query_name.startswith("ref")]
    hx = contig[: A + 1] + ins + contig[A + 1 :]
    third = next(b for b in "ACGT" if b != hx[A + 1 + 17])
    end = DEL + 2 if repeat else DEL  # the last genome base of the reads
    for i in range(3):
        start = 488 - i
        seq = hx[start : A + 1 + N + end - A]
        cigar = ((0, A + 1 - start), (1, N), (0, end - A))
        at = A + 1 + 17 - start
        reads.append(make_read(f"padplain{i}", seq, start, cigar))
        reads.append(make_read(f"padthird{i}", seq[:at] + third + seq[at + 1 :], start, cigar))
    fa, bam = write_contig(tmp_path, contig, reads, "pad")
    x_alt = contig[A] + ins
    variants = [
        _rs.Variant("1", A, contig[A], x_alt, "INDEL"),
        _rs.Variant("1", A, contig[A], contig[A] + sibling, "INDEL"),
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
    assert any(
        "the guard's masks laid on every candidate haplotype by offset" in m for m in messages
    )
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
    row = calls[x_alt]
    assert row[("padthird", "ALT")] == 0 and row[("padplain", "ALT")] == 0, row
    assert row[("padthird", "partial")] == 3 and row[("padplain", "partial")] == 3, row
    for kind in ("ref", "refsnp"):
        assert row[(kind, "REF")] == 4, (kind, row)
    x = counts[0]
    assert (x.rd, x.ad) == (8, 0), (x.rd, x.ad, x.partial_alt)
