"""The sequence judge (`--read-judge sequence`, #304): a co-annotated sibling at X's site.

Mechanism. A sibling the input lists over X's bases (an MNP's component, an SNV inside a
delins) is another allele only where the reads show it: on at least one fragment, and
on as many as an error's share (5%) of those spelling the allele it is nearest (REF or
X, the larger when as near; RJ-25, operator 2026-10-14). A component the reads never
show apart from the event is part of it, and at depth one fragment spelling it is an
error's share. Before this the judge made every listed sibling a state, so reads that
read only the bases it shares with X (the rest below min BQ) tied it with X and lost
their ALT call.

Shapes: MNPs of 2 and 3 bases, a component at each base they change; the component
spelled apart on 0, 1 or 5 fragments at 20 and 100 REF fragments; the masked carriers on
both strands. Counts are under the judge, explicitly.
"""

import random

import pytest
from helpers import count_checked, make_read, write_contig

from gbcms import _rs

READ = 100
P = 300

MNPS = {"2 bases": ("CC", "GT"), "3 bases": ("GCA", "ACT")}
CASES = [(name, j) for name, (r, a) in MNPS.items() for j in range(len(r)) if r[j] != a[j]]


def _contig(ref_block, seed=31, n=800):
    rng = random.Random(seed)
    c = [rng.choice("ACGT") for _ in range(n)]
    c[P - 1] = next(b for b in "ACGT" if b != ref_block[0])
    c[P : P + len(ref_block)] = ref_block
    c[P + len(ref_block)] = next(b for b in "ACGT" if b != ref_block[-1])
    return "".join(c)


def _reads(contig, ref_len, block, tag, n, first=0, low=()):
    out = []
    for i in range(first, first + n):
        s = P - 50 + (i % 40)
        seq = contig[s:P] + block + contig[P + ref_len : s + READ]
        r = make_read(f"{tag}{i}", seq, s, ((0, READ),), flag=16 if i % 2 else 0)
        if low:
            q = list(r.query_qualities)
            for k in low:
                q[P - s + k] = 5
            r.query_qualities = q
        out.append(r)
    return out


def _alt_count(tmp_path, name, j, alone, n_ref):
    """The MNP with its component at offset `j` co-annotated: `n_ref` REF fragments, 10
    carriers read in full, 4 reading only the component's base (X's other changes below
    min BQ), and `alone` fragments spelling the component apart."""
    ref_b, alt_b = MNPS[name]
    contig = _contig(ref_b)
    others = [k for k in range(len(ref_b)) if k != j and ref_b[k] != alt_b[k]]
    comp = ref_b[:j] + alt_b[j] + ref_b[j + 1 :]
    reads = _reads(contig, len(ref_b), ref_b, "ref", n_ref)
    reads += _reads(contig, len(ref_b), alt_b, "alt", 10, first=n_ref)
    reads += _reads(contig, len(ref_b), alt_b, "altlow", 4, first=n_ref + 10, low=others)
    reads += _reads(contig, len(ref_b), comp, "comp", alone, first=n_ref + 14)
    fa, bam = write_contig(tmp_path, contig, reads)
    pvs = _rs.prepare_variants(
        [
            _rs.Variant("1", P, ref_b, alt_b, "COMPLEX"),
            _rs.Variant("1", P + j, ref_b[j], alt_b[j], "SNP"),
        ],
        fa,
        5,
        False,
        1,
        True,
    )
    vs = [pv.variant for pv in pvs]
    c = count_checked(
        bam, vs, sibling_variants=[[vs[1]], [vs[0]]], read_judge="sequence", reference_fasta=fa
    )[0]
    assert c.dp >= c.rd + c.ad and c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    return c.ad


@pytest.mark.parametrize("name, j", CASES, ids=[f"{n} component {j}" for n, j in CASES])
def test_a_component_never_spelled_apart_is_part_of_the_mnp(tmp_path, name, j):
    """No fragment spells the component apart from the MNP: the 4 carriers reading only
    its base count ALT."""
    assert _alt_count(tmp_path, name, j, 0, 20) == 14


@pytest.mark.parametrize("name, j", CASES, ids=[f"{n} component {j}" for n, j in CASES])
@pytest.mark.parametrize(
    "alone, n_ref, ad",
    [(1, 20, 10), (1, 100, 14), (5, 100, 10)],
    ids=["1 of 20", "1 of 100", "5 of 100"],
)
def test_a_component_is_another_allele_only_above_an_errors_share(
    tmp_path, name, j, alone, n_ref, ad
):
    """Spelled apart on one fragment beside 20 REF fragments (5%) it is another allele
    (the 4 carriers reading only its base tie it with the MNP); beside 100 one fragment is
    an error's share and five are not."""
    assert _alt_count(tmp_path, name, j, alone, n_ref) == ad
