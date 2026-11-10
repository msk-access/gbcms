"""The sequence judge's graph scorer.

It scores a read, per state of X's site, by one dynamic programme over a local sequence
graph (reference segments; the other listed rows, the guard's background and recurring
indels as bubbles a path may take or skip; isolated masks as wildcards), so every
combination of the rows is weighed at once. It replaced enumerating the contexts and
scoring each haplotype (`GBCMS_JUDGE_SCORER=enumerated`, since removed).

Each shape below is one the enumeration got wrong before its context machinery (ranking,
caps, molecule contexts) was built, or a cost it paid for that machinery. Counts
invariants hold in every call.
"""

import random
import time

import pytest
import test_haplotype_judge_contract as hj
from helpers import make_read, write_contig

from gbcms import _rs


def _invariants(c):
    assert c.dp >= c.rd + c.ad and c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev


# ── The shapes the enumeration got wrong ──


@pytest.mark.parametrize("x_pos, carrier", [(32702397, "h0"), (32702395, "h1")])
def test_the_hla_two_insertion_site_reads_each_haplotype_in_its_own_context(
    tmp_path, caplog, x_pos, carrier
):
    """HLA-DQB1 (chr6:32702395/97): two SNVs in trans among 40 rows, each haplotype with its
    own insertion in one CA repeat. The enumeration tied X's carriers across every state
    until each molecule's rows were a context; the graph holds every combination."""
    tally = hj._two_haplotype_calls(tmp_path, caplog, hj.HLA, hj.HLA_G0, hj.HLA_ROWS, x_pos)
    other = "h1" if carrier == "h0" else "h0"
    assert tally[(other, "ALT")] == 0, dict(tally)
    assert tally[(carrier, "ALT")] >= 3 * tally[(carrier, "neither")], dict(tally)
    assert tally[(carrier, "REF")] == 0, dict(tally)


def test_the_chr19_insertion_haplotype_is_never_alt_without_its_molecules_context(tmp_path, caplog):
    """chr19:48569347 A>G (#246): the insertion's haplotype (an insertion, an SNV and a
    26-bp deletion) read as X beside the other haplotype's rows wherever its own rows were
    not a context together (the enumeration's missing context)."""
    anchor = 48569346 - hj.CHR19_G0 + 1

    def bwa(i):
        if i < anchor:
            return i
        if i == anchor:
            return None
        return i - 1 if i < anchor + 17 else i + 20

    tally = hj._two_haplotype_calls(
        tmp_path, caplog, hj.CHR19, hj.CHR19_G0, hj.CHR19_ROWS, 48569347, {1: bwa}
    )
    assert tally[("h1", "ALT")] == 0, dict(tally)
    assert tally[("h0", "ALT")] > 0 and tally[("h1", "REF")] > 0, dict(tally)


def test_the_chr3_snv_beside_an_insertion_is_read_whichever_allele_a_read_holds(tmp_path, caplog):
    """chr3:112233052 T>C one base after a 21-bp insertion in a CA/TA repeat (#246)."""
    tally = hj._two_haplotype_calls(tmp_path, caplog, hj.CHR3, hj.CHR3_G0, hj.CHR3_ROWS, 112233052)
    assert tally[("h0", "REF")] == 0 and tally[("h1", "ALT")] == 0, dict(tally)
    assert tally[("h0", "ALT")] > 0, dict(tally)


@pytest.mark.parametrize(
    "shape",
    [
        "test_an_insertion_carrier_written_through_counts_for_the_insertion_only",
        "test_a_molecule_carrying_both_in_cis_counts_at_both_rows",
        "test_two_insertions_at_one_anchor_count_their_own_carriers",
    ],
)
def test_the_246_pairs_in_trans_and_in_cis(tmp_path, shape):
    """#246: an SNV one base after an insertion's anchor, its ALT the first inserted base,
    in trans and in cis, and two insertions at one anchor: the judge contract's shapes."""
    getattr(hj, shape)(tmp_path)


# ── A dense cluster, in bounded time ──


def _dense_cluster(tmp_path, n_rows, n_reads, seed, indels=False):
    """`n_rows` rows (SNVs six bases apart; with `indels`, every fifth a 2-bp deletion or
    a 1-bp insertion) on one haplotype, `n_reads` reads of 150 bases from it and from REF,
    a third of the bases at the rows read at quality 2 (each read shows its own set)."""
    rng = random.Random(seed)
    contig = "".join(rng.choice("ACGT") for _ in range(1200 + 6 * n_rows))
    rows = []
    for k in range(n_rows):
        p = 500 + 6 * k
        if indels and k % 5 == 2:
            rows.append((p, contig[p : p + 3], contig[p]))
        elif indels and k % 5 == 4:
            rows.append((p, contig[p], contig[p] + rng.choice("ACGT")))
        else:
            rows.append((p, contig[p], "ACGT"[("ACGT".index(contig[p]) + 1 + k % 3) % 4]))
    hap, origin, at = [], [], 0
    for p, r, a in rows:
        hap += contig[at:p]
        origin += list(range(at, p))
        hap += a
        origin += [p + i if i < min(len(r), len(a)) else None for i in range(len(a))]
        at = p + len(r)
    hap += contig[at:]
    origin += list(range(at, len(contig)))
    hap = "".join(hap)
    positions = {p for p, _, _ in rows}
    reads = []
    for i in range(n_reads):
        quals = [30] * 150
        if i % 2:
            st = 360 + rng.randrange(0, 6 * n_rows + 100)
            seq = contig[st : st + 150]
            for j in range(150):
                if st + j in positions and rng.random() < 0.33:
                    quals[j] = 2
            reads.append(make_read(f"r{i}", seq, st, ((0, 150),), quals=quals))
        else:
            s = 360 + rng.randrange(0, 6 * n_rows + 100)
            while origin[s] is None:
                s += 1
            seq = hap[s : s + 150]
            for j in range(150):
                if origin[s + j] in positions and rng.random() < 0.33:
                    quals[j] = 2
            reads.append(make_read(f"h{i}", seq, origin[s], ((0, 150),), quals=quals))
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    pvs = _rs.prepare_variants(
        [_rs.Variant("1", p, r, a, "X") for p, r, a in rows],
        fa,
        5,
        False,
        1,
        True,
        read_judge="sequence",
    )
    return fa, bam, [pv.variant for pv in pvs]


def _timed(bam, fa, variants):
    t0 = time.monotonic()
    counts = _rs.count_bam_binned(
        bam,
        variants,
        [None] * len(variants),
        20,
        20,
        True,
        True,
        True,
        False,
        False,
        False,
        1,
        reference_fasta=fa,
        read_judge="sequence",
    )
    return counts, time.monotonic() - t0


def test_a_60_row_cluster_row_is_judged_in_a_few_seconds(tmp_path):
    """The judge contract's 60-row cluster (4,000 reads; the enumeration's bound was 40 s
    for one row): one row judged within 5 s on one thread, its carriers ALT, REF reads REF."""
    fa, bam, vs = _dense_cluster(tmp_path, 60, 4000, 93470262)
    (c,), elapsed = _timed(bam, fa, [vs[30]])
    _invariants(c)
    assert c.ad > 0 and c.rd > 0, (c.rd, c.ad)
    assert elapsed < 5, f"{elapsed:.1f}s"


def test_a_105_row_cluster_with_indels_is_judged_whole_in_bounded_time(tmp_path):
    """The chr3:93.47 Mb cluster's size (105 rows, SNVs with deletions and insertions among
    them; 2,000 reads): every row judged within two minutes on one thread, each row's
    carriers counted. The enumerated scorer ran past ten minutes here (c46d took 2 h 16 min
    on HG002's cluster)."""
    fa, bam, vs = _dense_cluster(tmp_path, 105, 2000, 9347, indels=True)
    counts, elapsed = _timed(bam, fa, vs)
    for c in counts:
        _invariants(c)
    assert sum(c.ad > 0 for c in counts) >= 0.9 * len(counts), [c.ad for c in counts]
    assert elapsed < 120, f"{elapsed:.1f}s"
