"""The general own-gap rule (`GBCMS_JUDGE_OWNGAP=all`, measurement only; the default `tie`
is the decided tie rule): a state the read reaches only through an indel of its own at that
state's change is not the read's, whatever it ties with or beats.

X a T inserted in a TT run beside a listed 33 bp deletion whose right placement ends at X's
anchor (another allele at X's site). A read with an unlisted 11 bp deletion inside the
listed deletion's bases reaches that allele only through its own gap: REF at X it is REF,
carrying X's T it is ALT, on either strand and under either rule; the deletion's carriers
are another allele. The judge reads bases, not strands, so the call is the strand's alike;
the counts split by strand as the reads do.
"""

import pytest
from helpers import count_checked, make_read, write_contig

from gbcms import _rs as gbcms_rs

MOTIF = "CAGAGAACAAAAGAAAGAAAGATTTTCTTACTTCAGAGAATTCTTTGCCACGTATTTCTAGC"


def _contig():
    import random

    rng = random.Random(59)
    c = [rng.choice("ACGT") for _ in range(600)]
    c[261 : 261 + len(MOTIF)] = MOTIF
    return "".join(c)


def _prepared(fa, rows):
    pvs = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", p, r, a, "COMPLEX") for p, r, a in rows], fa, 5, False, 1, True
    )
    assert all(pv.gbcms_status == "PASS" for pv in pvs), [pv.gbcms_status_reason for pv in pvs]
    return [pv.variant for pv in pvs]


@pytest.mark.parametrize("rule", ["tie", "all"])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("carries_x", [False, True])
def test_a_read_beside_another_allele_by_its_own_indel_keeps_its_call(
    tmp_path, monkeypatch, rule, reverse, carries_x
):
    monkeypatch.setenv("GBCMS_JUDGE_OWNGAP", rule)
    c = _contig()
    hx = c[:301] + "T" + c[301:]
    h_del = c[:261] + c[294:]
    reads = [make_read(f"ref{i}", c[240 + i : 340 + i], 240 + i, ((0, 100),)) for i in range(6)]
    reads += [
        make_read(f"alt{i}", hx[250 + i : 350 + i], 250 + i, ((0, 51 - i), (1, 1), (0, 48 + i)))
        for i in range(4)
    ]
    reads += [
        make_read(f"del{i}", h_del[200 + i : 300 + i], 200 + i, ((0, 61 - i), (2, 33), (0, 39 + i)))
        for i in range(3)
    ]
    h = hx if carries_x else c
    seq = h[275:281] + h[292:386]
    cigar = ((0, 6), (2, 11), (0, 9), (1, 1), (0, 84)) if carries_x else ((0, 6), (2, 11), (0, 94))
    reads.append(make_read("own", seq, 275, cigar, flag=16 if reverse else 0))
    fa, bam = write_contig(tmp_path, c, reads)
    x, _ = _prepared(fa, [(300, "A", "AT"), (260, c[260:294], c[260])])
    got = count_checked(bam, [x], read_judge="sequence", reference_fasta=fa)[0]
    assert got.dp >= got.rd + got.ad
    assert got.dpf >= got.rdf + got.adf
    assert got.rd == got.rd_fwd + got.rd_rev
    assert got.ad == got.ad_fwd + got.ad_rev
    assert (got.rd, got.ad) == (6 + (not carries_x), 4 + carries_x)
    own_rev = (got.ad_rev if carries_x else got.rd_rev) == (1 if reverse else 0)
    assert own_rev
