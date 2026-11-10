"""The sequence judge at an insertion of 50 bp or more: the guard's masks are laid by offset
from the anchor on every candidate haplotype along the same path (RJ-23, RJ-27, amended
2026-11-02).

Mechanism. A germline change both haplotypes carry, past the insertion's run, is masked by the
guard in both windows by offset from the anchor: on X's haplotype its place is at the offset it
has on REF's, plus the insert's length, and the base at the other offset is masked on each so
that it weighs on none. A listed insertion or deletion between the anchor and the masked base
moves that offset on every haplotype a read carrying it lies on. The first build laid the masks
at the offsets the haplotypes have without it. They are laid by the columns of each path now,
the graph holding a copy of the span the masks lie over per set of such changes that can lie
together (the unit tests of `rust/src/counting/judge/mod.rs` and `haplotypes.rs` pin the
offsets and the scores; no read's call changed on the shapes they were swept over).

Shape: a 60 bp duplication after the anchor at 0-based 500 (a copy of the 60 bases after it but
its last, so the base past its run tells the alleles apart), a germline change at 0-based 564 on
half the molecules, and a 2 bp deletion listed after 561 (two bases before the change), which
some molecules carry. The counts are the carriers of the duplication and the REF molecules, with
and without the deletion and the change.
"""

import collections
import random
import re

from helpers import count_bam_checked, make_read, write_contig

from gbcms import _rs

ARGS = (20, 20, True, True, True, False, False, False, 1)
CALL = re.compile(
    r"read call \S+ (\S+)>(\S+) read=([a-z]+)\d+ mate=(\d) ref=(\w+) alt=(\w+) .*?partial=(\d+)"
)
A = 500  # the anchor (0-based)
N = 60  # the duplication's length
SNP = 564  # the germline change
DEL = 561  # the listed deletion removes the two bases after this one


def _shape():
    """The contig, the duplication's ALT, the deletion's (REF, ALT), the read names by kind,
    and the reads."""
    rng = random.Random(2644)
    contig = "".join(rng.choice("ACGT") for _ in range(1200))
    last = next(b for b in "ACGT" if b != contig[A + N] and b != contig[A])
    ins = contig[A + 1 : A + N] + last
    snp_alt = next(b for b in "ACGT" if b != contig[SNP])
    reads = []

    def haplotype(dup: bool, deleted: bool, snp: bool) -> str:
        seq = list(contig)
        if snp:
            seq[SNP] = snp_alt
        seq = "".join(seq)
        if deleted:
            seq = seq[: DEL + 1] + seq[DEL + 3 :]
        if dup:
            seq = seq[: A + 1] + ins + seq[A + 1 :]
        return seq

    def add(kind, i, hap, start, cigar):
        reads.append(make_read(f"{kind}{i}", hap[start : start + 100], start, cigar))

    # REF molecules over the anchor and past the change: with it or without, with the listed
    # deletion (written as a 2 bp deletion) or without
    for i in range(4):
        start = 466 + 3 * i
        add("ref", i, haplotype(False, False, False), start, ((0, 100),))
        add("refsnp", i, haplotype(False, False, True), start, ((0, 100),))
        m1 = DEL + 1 - start
        cigar = ((0, m1), (2, 2), (0, 100 - m1))
        add("refdel", i, haplotype(False, True, False), start, cigar)
        add("refdelsnp", i, haplotype(False, True, True), start, cigar)
    # carriers of the duplication written with the insert
    for i in range(8):
        start = 463 + 3 * i
        m1 = A + 1 - start
        add("thr", i, haplotype(True, False, False), start, ((0, m1), (1, N), (0, 100 - m1 - N)))
    return contig, contig[A] + ins, reads


def test_a_listed_change_between_the_anchor_and_a_masked_base_leaves_the_counts(tmp_path, caplog):
    """Carriers written with the insert are ALT; REF molecules are REF with the germline change
    or without it, with the listed deletion or without; the masks are laid by the columns of
    each path (the judge says so at debug, with the deletion among the changes it counts); the
    counts are the carriers and the REF molecules, one variant per bin or the production bins
    alike."""
    contig, alt, reads = _shape()
    fa, bam = write_contig(tmp_path, contig, reads, "dupdel")
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
    laid = [
        m for m in messages if "the guard's masks laid on every candidate haplotype by offset" in m
    ]
    assert laid and all("through 1 context changes" in m for m in laid), laid
    seen = {}  # the read calls are logged by both counts (production bins, one variant per bin)
    for m in messages:
        if g := CALL.search(m):
            _, alt_, kind, mate, ref, alt_hit, partial = g.groups()
            call = (
                "ALT"
                if alt_hit == "true"
                else "REF" if ref == "true" else "partial" if partial != "0" else "depth"
            )
            seen[(alt_, m.split(" read=")[1].split()[0], mate)] = (kind, call)
    calls = collections.defaultdict(collections.Counter)
    for (alt_, _, _), kind_call in seen.items():
        calls[alt_][kind_call] += 1
    row = calls[alt]
    assert row[("thr", "ALT")] == 8, row
    for kind in ("ref", "refsnp", "refdel", "refdelsnp"):
        assert row[(kind, "REF")] == 4, (kind, row)
    assert sum(n for (_, call), n in row.items() if call == "ALT") == 8, row
    dup = counts[0]
    assert (dup.rd, dup.ad) == (16, 8), (dup.rd, dup.ad, dup.partial_alt)
