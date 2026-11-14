"""The sequence judge (`--read-judge sequence`): whether the reads need a found indel.

Mechanism. A found indel beside X's site stays a context only where the reads reaching it
need it for their calls: it takes a read from REF or ALT to another call when it is left
out, on at least two reads and a twentieth of all the reads reaching it. The reads weighed
are every read reaching the item, in no order: the count must not depend on how deep the
site is, or on which reads a coordinate sort puts first.

Shapes (made-up sequence only): two SNVs 9 bases apart (1501 C>T-like and 1510 G>A-like,
by the contig's own bases), a found deletion starting 10 bases after the second SNV and
not listed as a row; 600 REF reads, 200 ALT reads without the deletion, 184 carriers of
both SNVs and the deletion (their left piece 1 to 149 bases), 150-base reads, both strands.
A 100-base deletion is needed by many carriers, so the item is a context and the judge
gives the rules path's AD at both rows. A known partial stands as a strict xfail: a
30-base deletion whose needing carriers are under a twentieth of the reads reaching it.
"""

import random

import pytest
from helpers import count_checked, make_read, write_contig

from gbcms import _rs

L = 150
P = 1500  # SNV1 (0-based)
S2 = P + 9  # SNV2
GAP = 10  # bases from SNV2 to the first deleted base


def _shape(dlen, ref_reads=600, alt_nodel=200, carriers=184, seed=7):
    """(contig, reads, rows): the shape's made-up molecule, reads and the two SNV rows."""
    rng = random.Random(seed)
    contig = [rng.choice("ACGT") for _ in range(6000)]
    d0 = S2 + GAP  # first deleted base
    nxt = {"A": "C", "C": "G", "G": "T", "T": "A"}
    # no slide: the anchor base differs from the last deleted base, the first deleted
    # base from the base after it
    if contig[d0 - 1] == contig[d0 + dlen - 1]:
        contig[d0 + dlen - 1] = nxt[contig[d0 - 1]]
    if contig[d0] == contig[d0 + dlen]:
        contig[d0 + dlen] = nxt[contig[d0]]
    contig = "".join(contig)
    tr = {"A": "G", "C": "T", "G": "A", "T": "C"}
    alt1, alt2 = tr[contig[P]], tr[contig[S2]]
    hap_a = list(contig)
    hap_a[P], hap_a[S2] = alt1, alt2
    hap_a = "".join(hap_a)
    hap_d = hap_a[:d0] + hap_a[d0 + dlen :]
    reads = []
    lo, hi = S2 - L + 1, P
    for i in range(ref_reads):
        s = rng.randint(lo, hi)
        reads.append(make_read(f"ref{i}", contig[s : s + L], s, ((0, L),), flag=16 if i % 2 else 0))
    for i in range(alt_nodel):
        s = rng.randint(lo, hi)
        reads.append(make_read(f"alt{i}", hap_a[s : s + L], s, ((0, L),), flag=16 if i % 2 else 0))
    for i in range(carriers):
        k = rng.randint(1, L - 1)  # bases before the deletion
        s = d0 - k
        cig = ((0, k), (2, dlen), (0, L - k))
        reads.append(make_read(f"del{i}_k{k}", hap_d[s : s + L], s, cig, flag=16 if i % 2 else 0))
    rows = [(P, contig[P], alt1, "SNP"), (S2, contig[S2], alt2, "SNP")]
    return contig, reads, rows


class AdDiffers(AssertionError):
    """The judge's AD differs from the rules path's (the only failure the xfail expects)."""


def _alt_by_bases(reads, rows):
    """Per row, the reads whose own aligned bases carry the row's ALT."""
    out = []
    for pos, _, alt, _ in rows:
        n = 0
        for r in reads:
            at = {rp: qp for qp, rp in r.get_aligned_pairs(matches_only=True)}
            base = at.get(pos)
            n += base is not None and r.query_sequence[base] == alt
        out.append(n)
    return out


def _ad(tmp_path, judge, dlen):
    contig, reads, rows = _shape(dlen)
    fa, bam = write_contig(tmp_path, contig, reads, name=judge)
    vs = [_rs.Variant("1", p, r, a, t) for p, r, a, t in rows]
    pvs = _rs.prepare_variants(vs, fa, 5, False, 1, True, read_judge=judge)
    cs = count_checked(bam, [p.variant for p in pvs], read_judge=judge, reference_fasta=fa)
    for c in cs:
        assert c.dp >= c.rd + c.ad and c.dpf >= c.rdf + c.adf
        assert c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    return [c.ad for c in cs]


def test_a_deletion_every_needing_carrier_reads_is_context_at_any_depth(tmp_path):
    """The 100-base deletion is needed by the carriers over the SNVs: weighing only the
    first reads reaching it left it out as no context and lost reads of each SNV's ALT."""
    rules = _ad(tmp_path, "rules", 100)
    judge = _ad(tmp_path, "sequence", 100)
    assert judge == rules
    # the shape's truth by bases: 200 ALT reads without the deletion plus the carriers
    # whose left piece reaches each SNV
    _, reads, rows = _shape(100)
    by_bases = _alt_by_bases(reads, rows)
    assert by_bases == [369, 378]
    assert judge == by_bases


@pytest.mark.xfail(
    strict=True,
    raises=AdDiffers,
    reason=(
        "known partial: a deletion whose needing carriers are under a twentieth of the "
        "reads that reach it is not a context, so the reads carrying it lose their ALT"
    ),
)
def test_a_deletion_few_carriers_need_is_context(tmp_path):
    rules = _ad(tmp_path, "rules", 30)
    judge = _ad(tmp_path, "sequence", 30)
    if judge != rules:
        raise AdDiffers(f"judge {judge} against rules {rules}")
