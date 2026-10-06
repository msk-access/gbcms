"""S3 #239: Fisher's exact test at any depth.

`fisher_exact_2x2` (strand bias, fragment strand bias, the ASJD junction test, and
merge's combined strand bias) returned p = 0 whenever C(n, a + c) exceeded f64's
range: statrs's binomial overflowed and the sum stayed 0. That needs both columns
(the forward and reverse totals) wide, so on strand-balanced tables it began near
1,030 reads, while a table with one thin strand stayed correct at any depth. On the RC
cfDNA rows 316 of 1,060 read-level strand-bias p-values were 0 where the exact p is
above 0 (298 of them >= 0.05). It also counted every table less likely than 1e-10
regardless of the observed one, so strongly biased tables were floored near 1e-10.

Contract (operator, 2026-10-06): the exact two-sided p on the raw counts, R's
`fisher.test` definition, i.e. the sum over tables whose probability is at most the
observed one's times (1 + 1e-7). The oracle here is exact: integer binomials.
"""

import math
import random
from fractions import Fraction

import pytest
from helpers import build_bam, count_one_checked, make_read

from gbcms._rs import Variant, fisher_exact_2x2

REL = 1e-9
# Deep tables: a log recurrence started at one end of the range drifts as its sum
# grows (to ~n ln 2 on balanced tables); started at the mode it does not.
REL_DEEP = 1e-12


def oracle(a, b, c, d):
    """Exact two-sided p (R's relative tie tolerance 1 + 1e-7), with gbcms's guards:
    1.0 when row 2 holds <= 1 observation or a marginal is empty or full."""
    n = a + b + c + d
    if n == 0 or c + d <= 1:
        return 1.0
    r1, c1 = a + b, a + c
    if r1 in (0, n) or c1 in (0, n):
        return 1.0
    lo, hi = max(0, r1 + c1 - n), min(r1, c1)
    # Hypergeometric weights C(r1, k) C(n - r1, c1 - k) by the exact recurrence
    # (each division is exact, since every weight is an integer).
    w = [math.comb(r1, lo) * math.comb(n - r1, c1 - lo)]
    for k in range(lo, hi):
        w.append(w[-1] * (r1 - k) * (c1 - k) // ((k + 1) * (n - r1 - c1 + k + 1)))
    wo = w[a - lo]
    num = sum(x for x in w if x * 10**7 <= wo * (10**7 + 1))
    return float(Fraction(num, math.comb(n, c1)))


def close(got, want):
    return abs(got - want) <= REL * max(want, 1e-300)


DEEP = [
    (258, 258, 258, 258),  # balanced, n = 1032: exact p = 1.0
    (300, 280, 260, 270),  # n = 1110
    (30, 25, 500, 520),  # n = 1075
    (1500, 1400, 40, 35),  # a deep REF row, a thin ALT row
    (5000, 4800, 30, 45),  # n = 9875
    (12000, 11500, 400, 360),  # RNA depth
]
BIASED = [
    (250, 50, 20, 280),  # reported 2.1e-10; exact 7.0e-90
    (300, 10, 5, 290),
    (200, 100, 10, 200),
    (150, 150, 2, 200),
]


@pytest.mark.parametrize("table", DEEP, ids=[f"n{sum(t)}" for t in DEEP])
def test_deep_tables_give_the_exact_p(table):
    got, want = fisher_exact_2x2(*table)[0], oracle(*table)
    assert close(got, want), (table, got, want)


@pytest.mark.parametrize("table", BIASED, ids=[str(t) for t in BIASED])
def test_strongly_biased_tables_are_not_floored(table):
    got, want = fisher_exact_2x2(*table)[0], oracle(*table)
    assert close(got, want), (table, got, want)


BALANCED_DEEP = [
    (5100, 4900, 4950, 5050),  # n = 20,000
    (10200, 9800, 9900, 10100),  # n = 40,000
]


@pytest.mark.parametrize("table", BALANCED_DEEP, ids=[f"n{sum(t)}" for t in BALANCED_DEEP])
def test_precision_holds_at_depth(table):
    got, want = fisher_exact_2x2(*table)[0], oracle(*table)
    assert abs(got - want) <= REL_DEEP * want, (table, got, want)


# Tables with another table as likely as the observed one to within R's 1 + 1e-7 but
# not scipy's 1 + 1e-14: R counts it, so R's p (and gbcms's) differs from scipy's.
NEAR_TIES = [
    # table, R's fisher.test p, scipy 1.15's p
    ((531, 913, 484, 623), 4.461780029e-4, 3.86606e-4),
    ((202, 287, 360, 498), 0.8631643753, 0.818544),
    ((397, 404, 799, 804), 0.9310402918, 0.896817),
]


@pytest.mark.parametrize("table,r_p,scipy_p", NEAR_TIES, ids=[str(t[0]) for t in NEAR_TIES])
def test_near_ties_follow_rs_rule(table, r_p, scipy_p):
    got = fisher_exact_2x2(*table)[0]
    assert close(got, oracle(*table)), (table, got)
    assert got == pytest.approx(r_p, rel=1e-9)
    assert got != pytest.approx(scipy_p, rel=1e-3)


def test_extreme_margins_stay_bounded():
    """Cells at u32 max: terms more than e^-750 below the mode add exactly 0, so the
    sum stops there and time and memory stay bounded (the whole range is 4.3e9 tables)."""
    m = 2**32 - 1
    assert fisher_exact_2x2(m, m, m, m)[0] == 1.0
    assert fisher_exact_2x2(m, 0, 0, m)[0] == 0.0
    assert 0.0 < fisher_exact_2x2(m, m - 200_000, m - 200_000, m)[0] < 1.0


def test_a_grid_of_tables_gives_the_exact_p():
    """Every table with cells 0-9, plus 300 seeded random tables up to n = 3,000."""
    tables = [
        (a, b, c, d) for a in range(10) for b in range(10) for c in range(10) for d in range(10)
    ]
    rng = random.Random(239)
    tables += [
        tuple(rng.randint(0, 1500) for _ in range(2)) + tuple(rng.randint(0, 60) for _ in range(2))
        for _ in range(300)
    ]
    pairs = ((t, fisher_exact_2x2(*t)[0], oracle(*t)) for t in tables)
    bad = [x for x in pairs if not close(x[1], x[2])]
    assert not bad, f"{len(bad)} of {len(tables)} tables differ, e.g. {bad[:3]}"


def test_the_guards_and_odds_ratio_are_unchanged():
    """Guard: row 2 with <= 1 observation and empty or full marginals give p = 1;
    the odds ratio is (a*d)/(b*c), NaN when b*c = 0."""
    assert fisher_exact_2x2(0, 0, 0, 0)[0] == 1.0
    assert fisher_exact_2x2(50, 40, 1, 0)[0] == 1.0  # row 2 total 1
    assert fisher_exact_2x2(0, 0, 5, 7)[0] == 1.0  # row 1 empty
    assert fisher_exact_2x2(5, 0, 7, 0)[0] == 1.0  # column 2 empty
    assert fisher_exact_2x2(10, 20, 30, 40)[1] == pytest.approx(10 * 40 / (20 * 30))
    assert math.isnan(fisher_exact_2x2(10, 0, 30, 40)[1])


def test_strand_bias_columns_are_exact_at_depth(tmp_path):
    """1,240 single-end reads over an SNV, strand-balanced on both alleles (600/600
    REF, 20/20 ALT): both strand-bias tests give p = 1, not 0."""
    contig = "ACGT" * 50
    pos = 100
    alt = "T" if contig[pos] != "T" else "G"
    reads = []
    for i in range(1240):
        is_alt = i >= 1200
        rev = (i % 2) == 1
        start = pos - 40 + (i % 20)
        seq = contig[start : start + 80]
        if is_alt:
            seq = seq[: pos - start] + alt + seq[pos - start + 1 :]
        reads.append(make_read(f"r{i}", seq, start, ((0, 80),), flag=16 if rev else 0))
    bam = build_bam(tmp_path, reads)
    counts = count_one_checked(bam, Variant("chr1", pos, contig[pos], alt, "SNP"))
    assert counts.dp >= counts.rd + counts.ad
    assert counts.dpf >= counts.rdf + counts.adf
    assert counts.rd == counts.rd_fwd + counts.rd_rev
    assert counts.ad == counts.ad_fwd + counts.ad_rev
    assert (counts.rd_fwd, counts.rd_rev, counts.ad_fwd, counts.ad_rev) == (600, 600, 20, 20)
    assert counts.sb_pval == pytest.approx(oracle(600, 600, 20, 20)) == 1.0
    assert counts.fsb_pval == pytest.approx(
        oracle(counts.rdf_fwd, counts.rdf_rev, counts.adf_fwd, counts.adf_rev)
    )
