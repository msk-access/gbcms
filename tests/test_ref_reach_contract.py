"""Tell what the read reaches (RJ-25, amended 2026-11-02; RJ-24 for a delins).

Mechanism. A listed allele whose span meets X's site is always another allele at that site
(the input is taken as correct, invariant 7): the 1-fragment-and-5% floor of 2026-10-14 is
gone for listed alleles, on the rules and under `--read-judge sequence`. A read counts REF
where its confident bases tell REF from X and from every REF-length allele whose difference
its span reaches; an allele past the read's end does not contest it (RJ-9). Reached and not
told, the read is neither. A read counts ALT only where it tells X from every other allele
of X's length (a tie is never credited). An SNV row where the input lists an indel or delins
that meets its site is counted by the sequence judge also under the rules (RJ-25, 2026-11-03);
an MNV row stays on the rules. At a delins (or deletion) whose changed REF bases number 50 or more, a substitution inside
its core (4 bases in from both ends) is a context: a read that reaches it and spells it is
REF.

Counts are (RD, AD, partial_alt) of X's row. Every row of every call is checked against the
four counting invariants, on both paths, which must agree.
"""

import random

import pytest
from helpers import count_checked, make_read, write_contig

from gbcms import _rs

READ, P = 100, 400
PATHS = ["rules", "sequence"]
# develop's rules counts of X (RD, AD, partial) in the two MNV shapes below, measured on 3b63e6d0
MNV_DELETION_ON_DEVELOP = (10, 18, 0)
MNV_INSERTION_ON_DEVELOP = (10, 6, 18)


def _contig(seed=7, n=1100):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _other(b, avoid=""):
    return next(x for x in "ACGT" if x != b and x not in avoid)


def _flag(i):
    return 16 if i % 2 else 0


def _counts(tmp_path, judge, contig, rows, reads, siblings=True):
    """The counts of every row (pos, ref, alt), each the others' sibling, under `judge`;
    the four counting invariants asserted on each."""
    fa, bam = write_contig(tmp_path, contig, reads)
    pvs = _rs.prepare_variants(
        [_rs.Variant("1", p, r, a, "COMPLEX") for p, r, a in rows],
        fa,
        5,
        False,
        1,
        True,
        read_judge=judge,
    )
    assert all(pv.gbcms_status == "PASS" for pv in pvs), [pv.gbcms_status_reason for pv in pvs]
    vs = [pv.variant for pv in pvs]
    sibs = [[w for w in vs if w is not v] if siblings else [] for v in vs]
    got = count_checked(bam, vs, sibling_variants=sibs, read_judge=judge, reference_fasta=fa)
    for c in got:
        assert c.dp >= c.rd + c.ad
        assert c.dpf >= c.rdf + c.adf
        assert c.rd == c.rd_fwd + c.rd_rev
        assert c.ad == c.ad_fwd + c.ad_rev
    return got


def _call(c):
    return (c.rd, c.ad, c.partial_alt)


def _made(root, names):
    """A folder under `root` for each name."""
    for n in names:
        (root / n).mkdir()
    return names


# ── a delins X with a listed SNV S among the bases it replaces ───────────


def _delins(contig, ref_len, alt_len):
    """X: `ref_len` REF bases replaced by `alt_len` others, neither end shared with REF."""
    rng = random.Random(ref_len * 1000 + alt_len)
    ref = contig[P : P + ref_len]
    alt = _other(ref[0]) + "".join(rng.choice("ACGT") for _ in range(alt_len - 2)) + _other(ref[-1])
    return (P, ref, alt)


def _carriers(contig, x, n, tag="car"):
    """`n` reads carrying X, written as an insertion then a deletion."""
    _, ref, alt = x
    out = []
    left_max = min(60, READ - len(alt) - 12)
    for i in range(n):
        left = left_max - (i % 20)
        s = P - left
        seq = contig[s:P] + alt + contig[P + len(ref) : P + len(ref) + READ - left - len(alt)]
        cigar = ((0, left), (1, len(alt)), (2, len(ref)), (0, READ - left - len(alt)))
        out.append(make_read(f"{tag}{i}", seq, s, cigar, flag=_flag(i)))
    return out


def _ref_reads(contig, tag, n, start, span=1, at=None, base=None):
    """`n` REF reads starting at P + start + (i % span); with `at`, the base at P + at is `base`
    (a read spelling S)."""
    out = []
    for i in range(n):
        s = P + start + (i % span if span > 1 else i)
        seq = list(contig[s : s + READ])
        if at is not None:
            seq[P + at - s] = base
        out.append(make_read(f"{tag}{i}", "".join(seq), s, ((0, READ),), flag=_flag(i)))
    return out


def _example(tmp_path, judge, ref_len, spellers, alt_len=2, at=19):
    """The plan's worked example: X replacing `ref_len` bases, a listed SNV S at the base
    `at` of the replaced ones; 20 carriers; 30 REF reads that reach S's base and read REF;
    50 REF reads that end before it; `spellers` REF-molecule reads carrying S's base, each
    holding the whole REF window. Returns X's counts."""
    contig = _contig()
    x = _delins(contig, ref_len, alt_len)
    s = (P + at, contig[P + at], _other(contig[P + at]))
    reads = _carriers(contig, x, 20)
    reads += _ref_reads(contig, "reach", 30, -60)
    reads += _ref_reads(contig, "short", 50, -90, span=10)
    reads += _ref_reads(contig, "spell", spellers, -37, at=at, base=s[2])
    return _counts(tmp_path, judge, contig, [x, s], reads)[0]


# (REF bases replaced, the replaced base S changes, what S is): the plan's worked example at 30>2,
# its 60>2 twin with S three bases from the end (a junction: a state), and the same S at the
# 20th replaced base (inside the core of an event of 50 bases or more: a context).
SHAPES = {
    "30>2": (30, 19, False),
    "60>2 at a junction": (60, 57, False),
    "60>2 in the core": (60, 19, True),
}


def _want(core, spellers):
    """(RD, AD, partial) with `spellers` S-spelling reads: they are S's (partial) where S is
    another allele, REF's where S is a context."""
    return (80 + spellers, 20, 0) if core else (80, 20, spellers)


@pytest.mark.parametrize("judge", PATHS)
@pytest.mark.parametrize("shape", SHAPES)
def test_a_listed_allele_no_read_spells_leaves_the_short_ref_reads_ref(tmp_path, judge, shape):
    """S is listed and no read spells it. The 50 REF reads that end before S's base cannot
    tell REF from S, but they do not reach it: they count REF, whatever the length of X
    (a rule that made S contest every REF read would give RD 30 here)."""
    ref_len, at, _ = SHAPES[shape]
    assert _call(_example(tmp_path, judge, ref_len, 0, at=at)) == (80, 20, 0)


@pytest.mark.parametrize("judge", PATHS)
@pytest.mark.parametrize("shape", SHAPES)
def test_the_worked_example_with_s_spelled_on_six_fragments(tmp_path, judge, shape):
    """S spelled on 6 fragments, each holding the whole REF window: S is another allele on
    both paths, so its 6 readers are partial evidence, the 30 REF reads that reach S's base
    read REF there, and the 50 that end before it do not reach it: RD 80, AD 20. (Inside the
    core of a long event S is a context and the 6 readers are REF: RD 86.)"""
    ref_len, at, core = SHAPES[shape]
    assert _call(_example(tmp_path, judge, ref_len, 6, at=at)) == _want(core, 6)


@pytest.mark.parametrize("shape", SHAPES)
@pytest.mark.parametrize("spellers", [1, 6])
def test_both_paths_agree_on_the_worked_example(tmp_path, shape, spellers):
    """One fragment spelling S already makes S another allele, the same on both paths."""
    ref_len, at, core = SHAPES[shape]
    got = {}
    for j in PATHS:
        (tmp_path / j).mkdir()
        got[j] = _call(_example(tmp_path / j, j, ref_len, spellers, at=at))
    assert got["rules"] == got["sequence"] == _want(core, spellers), got


# ── a found allele of REF's length, past where short REF reads end ───────


@pytest.mark.parametrize("judge", PATHS)
def test_a_found_ref_length_allele_past_a_ref_reads_end_does_not_contest_it(tmp_path, judge):
    """A 28>6 delins; 6 fragments (above the floor of an unlisted allele: 3 and 5%) carry a
    change at the 21st replaced base. 8 REF reads end at the 13th replaced base, before it,
    and 20 read it as REF's: all 28 count REF. The readers of the change are neither."""
    contig = _contig(seed=277)
    x = _delins(contig, 28, 6)
    at = 20
    new = _other(contig[P + at])
    reads = _carriers(contig, x, 10)
    reads += _ref_reads(contig, "full", 20, -50)
    reads += _ref_reads(contig, "end", 8, -87, span=4)
    reads += _ref_reads(contig, "found", 6, -30, at=at, base=new)
    (c,) = _counts(tmp_path, judge, contig, [x], reads, siblings=False)[:1]
    assert (c.rd, c.ad) == (28, 10)


# ── An SNV or MNV row beside a listed allele that changes the length ──


def _substitution_row(tmp_path, judge, shape):
    """X a substitution; S a listed allele that spells X's bases and changes the length (or,
    for 'holds', keeps it). 10 reads carry X, 10 carry S, 10 are REF."""
    contig = _contig(seed=1)
    c = contig
    if shape.startswith("mnv"):
        x = (P, c[P : P + 2], _other(c[P]) + _other(c[P + 1]))
    else:
        x = (P, c[P], _other(c[P]))
    end = P + len(x[1])
    if shape == "snv + insertion":
        s = (P, c[P], x[2] + _other(c[P + 1]) + _other(c[P + 2]))
    elif shape == "snv + deletion":
        s = (P, c[P : P + 3], x[2] + c[P + 2])
    elif shape == "mnv + deletion":
        tail = c[end : end + 5]
        s = (P, x[1] + tail, x[2] + tail[:3] + tail[4:])
    else:  # holds: S keeps the length and carries X's change with another beside it
        tail = c[end : end + 4]
        s = (P, x[1] + tail, x[2] + tail[:3] + _other(tail[3]))

    def reads_of(tag, v, n):
        _, r, a = v
        out = []
        for i in range(n):
            st = P - 50 + (i % 40)
            left = P - st
            seq = (c[st:P] + a + c[P + len(r) :])[:READ]
            cigar = [(0, left + min(len(a), len(r)))]
            if len(a) > len(r):
                cigar.append((1, len(a) - len(r)))
            if len(a) < len(r):
                cigar.append((2, len(r) - len(a)))
            cigar.append((0, READ - left - len(a)))
            out.append(make_read(f"{tag}{i}", seq, st, tuple(cigar), flag=_flag(i)))
        return out

    reads = reads_of("x", x, 10) + reads_of("s", s, 10) + _ref_reads(c, "ref", 10, -50, span=40)
    return _counts(tmp_path, judge, contig, [x, s], reads)[0]


@pytest.mark.parametrize("judge", PATHS)
@pytest.mark.parametrize("shape", ["snv + insertion", "snv + deletion"])
def test_a_read_spelling_a_listed_length_changing_allele_is_not_alt_at_an_snv_row(
    tmp_path, judge, shape
):
    """X is an SNV and a listed allele that changes the length meets its site. Its 10 readers
    show X's base but spell that allele whole: partial evidence, not AD. On the rules path
    prep routes such a row to the sequence judge (RJ-25, 2026-11-03): the single-base check
    reads one base and weighs no listed allele that changes the length."""
    assert _call(_substitution_row(tmp_path, judge, shape)) == (10, 10, 10)


@pytest.mark.parametrize("judge", PATHS)
def test_a_listed_allele_of_the_same_length_that_holds_x_still_counts_at_both_rows(tmp_path, judge):
    """RJ-28: S keeps the length and holds X's change, so its readers are X's ALT too."""
    assert _call(_substitution_row(tmp_path, judge, "holds")) == (10, 20, 0)


def _padded_row(tmp_path, judge, shape):
    """X an SNV given with a base of padding after it (prep keeps a same-length row as given:
    its span is two bases, one differing). 'delins': S a delins over the SNV that deletes the
    padding base. 'insertion': S a pure insertion anchored on X's last base. 10 readers each."""
    c = _contig(seed=1)
    x1 = _other(c[P])
    x = (P, c[P : P + 2], x1 + c[P + 1])
    s = (
        (P, c[P : P + 3], x1 + c[P + 2])
        if shape == "delins"
        else (P + 1, c[P + 1], c[P + 1] + "TT")
    )

    def alt_reads(tag, v, n):
        _, r, a = v
        out = []
        for i in range(n):
            st = P - 50 + (i % 40)
            left = P - st
            seq = (c[st:P] + a + c[P + len(r) :])[:READ] if v[0] == P else None
            if seq is None:  # an insertion anchored one base on
                seq = (c[st : P + 1] + a + c[P + 2 :])[:READ]
                left += 1
            cigar = [(0, left + min(len(a), len(r)))]
            if len(a) > len(r):
                cigar.append((1, len(a) - len(r)))
            if len(a) < len(r):
                cigar.append((2, len(r) - len(a)))
            cigar.append((0, READ - left - len(a)))
            out.append(make_read(f"{tag}{i}", seq, st, tuple(cigar), flag=_flag(i)))
        return out

    reads = alt_reads("x", x, 10) + alt_reads("s", s, 10) + _ref_reads(c, "ref", 10, -50, span=40)
    return _counts(tmp_path, judge, c, [x, s], reads)[0]


def test_an_snv_given_with_padding_is_routed_like_a_bare_one(tmp_path):
    """The row counts as an SNV by the one base that differs, not by the length of its alleles:
    a listed delins over it routes it to the judge (develop's rules gave (10, 10, 10), the judge
    gives (10, 20, 0): the delins holds the row's change at its first base)."""
    got = {j: _call(_padded_row(tmp_path / j, j, "delins")) for j in _made(tmp_path, PATHS)}
    assert got["rules"] == got["sequence"] == (10, 20, 0), got


def test_a_pure_insertion_anchored_on_the_last_base_of_a_padded_snv_does_not_route_it(tmp_path):
    """The judge re-anchors it on the next base, where it no longer meets the row: the row keeps
    the rules' count, as on develop (the judge's own count is (10, 10, 10))."""
    assert _call(_padded_row(tmp_path, "rules", "insertion")) == (10, 10, 0)


def _mnv_row(tmp_path, judge, shape):
    """X a 2-base MNV; S a listed allele that changes the length. 'mnv + deletion': S a delins
    from X's first base. 'mnv + insertion': S an 18-base insertion anchored on X's first base,
    and 12 reads carry both (the first base X's, then the insertion, then X's second base)."""
    contig = _contig(seed=1)
    c = contig
    x = (P, c[P : P + 2], _other(c[P]) + _other(c[P + 1]))
    ins = "".join(_other(c[P + 1 + k]) for k in range(18))
    if shape == "mnv + deletion":
        end = P + 2
        tail = c[end : end + 5]
        s = (P, x[1] + tail, x[2] + tail[:3] + tail[4:])
        both = None
    else:
        s = (P, c[P], c[P] + ins)
        both = (c[:P] + x[2][0] + ins + x[2][1] + c[P + 2 :], P + 1, 18)

    def alt_reads(tag, v, n):
        _, r, a = v
        out = []
        for i in range(n):
            st = P - 50 + (i % 40)
            left = P - st
            seq = (c[st:P] + a + c[P + len(r) :])[:READ]
            cigar = [(0, left + min(len(a), len(r)))]
            if len(a) > len(r):
                cigar.append((1, len(a) - len(r)))
            if len(a) < len(r):
                cigar.append((2, len(r) - len(a)))
            cigar.append((0, READ - left - len(a)))
            out.append(make_read(f"{tag}{i}", seq, st, tuple(cigar), flag=_flag(i)))
        return out

    reads = alt_reads("x", x, 6) + alt_reads("s", s, 6) + _ref_reads(c, "ref", 10, -50, span=40)
    if both:
        hap, cut, n_ins = both
        for i in range(12):
            st = P - 50 + (i % 40)
            seq = hap[st : st + READ]
            reads.append(
                make_read(
                    f"b{i}",
                    seq,
                    st,
                    ((0, cut - st), (1, n_ins), (0, READ - (cut - st) - n_ins)),
                    flag=_flag(i),
                )
            )
    else:
        reads += alt_reads("y", x, 6)
    return _counts(tmp_path, judge, contig, [x, s], reads)[0]


@pytest.mark.parametrize("shape", ["mnv + deletion", "mnv + insertion"])
def test_an_mnv_row_beside_a_listed_length_changing_allele_stays_on_the_rules_path(tmp_path, shape):
    """Only an SNV row is routed to the judge: an MNV row keeps the rules' count whatever is
    listed at its site, as on develop (beside a listed delins the readers of both count ALT; the
    sequence judge would call 6 of them partial evidence)."""
    want = {"mnv + deletion": MNV_DELETION_ON_DEVELOP, "mnv + insertion": MNV_INSERTION_ON_DEVELOP}
    assert _call(_mnv_row(tmp_path, "rules", shape)) == want[shape]


# The reference of a short tandem tract (public GRCh38, HG002 truth haplotypes): 100 bases.
TRACT = (
    "TAACTAAAATGCTTGTTGTGCCCTGGGTTGCGCGCGCGCGCGCACACACACACACACACACACACACACACACACACGTGTT"
    "TTGAAAGAAGAGTGGGAA"
)
LEFT = 300


def _tract_locus(tmp_path, judge):
    """An SNV X (G>A, the 42nd base of the tract) beside two listed insertions: a 6-base one
    3 bases before it and a 10-base one at its anchor. The X haplotype carries X, the SNV two
    bases before it and the 6-base insertion; the other haplotype carries the 10-base insertion.
    12 reads of each haplotype, 6 REF reads. Returns X's counts."""
    rng = random.Random(15)
    flank = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
    contig = flank(LEFT) + TRACT + flank(600)
    i6, i2, ix = (
        LEFT + 37,
        LEFT + 39,
        LEFT + 41,
    )  # the G the 6-base insertion follows; two SNV bases
    assert contig[i6] == contig[i2] == contig[ix] == "G"
    rows = [
        (i6, "G", "GCACACA"),
        (i2, "G", "A"),
        (ix, "G", "A"),
        (ix, "G", "GCACACACACA"),
    ]

    def read(tag, n, edits, ins_at, ins):
        out = []
        for k in range(n):
            s = LEFT + 37 - 50 + (k % 30)
            hap = list(contig[s : s + READ + 20])
            for g, base in edits:
                hap[g - s] = base
            left = ins_at + 1 - s
            seq = "".join(hap[:left]) + ins + "".join(hap[left:])
            seq = seq[:READ]
            cigar = ((0, left), (1, len(ins)), (0, READ - left - len(ins))) if ins else ((0, READ),)
            out.append(make_read(f"{tag}{k}", seq, s, cigar, flag=_flag(k)))
        return out

    reads = read("x", 12, [(i2, "A"), (ix, "A")], i6, "CACACA")
    reads += read("o", 12, [], ix, "CACACACACA")
    reads += read("ref", 6, [], 0, "")
    return _counts(tmp_path, judge, contig, rows, reads)[2]


@pytest.mark.parametrize("judge", PATHS)
def test_reads_carrying_x_and_a_listed_insertion_nearby_stay_alt_on_both_paths(tmp_path, judge):
    """The 12 reads of the haplotype that holds X and a 6-base insertion nearby carry X's base
    and are X's ALT; the 12 reads of the other haplotype (a 10-base insertion at X's anchor)
    show REF's base at X and are not. A listed insertion that the reads do not carry makes
    none of X's carriers partial evidence."""
    assert _tract_locus(tmp_path, judge).ad == 12


# ── RJ-24 for a delins: the core of a long event is a context ────────────


def _core(tmp_path, judge, ref_len, alt_len, at, spellers=6):
    """X replaces `ref_len` bases with `alt_len`; a listed SNV S at the replaced base `at`;
    20 carriers, 30 REF reads, `spellers` REF-molecule reads carrying S's base that reach it."""
    contig = _contig()
    x = _delins(contig, ref_len, alt_len)
    s = (P + at, contig[P + at], _other(contig[P + at]))
    reads = _carriers(contig, x, 20)
    reads += _ref_reads(contig, "ref", 30, -60)
    reads += _ref_reads(contig, "spell", spellers, -45, at=at, base=s[2])
    return _counts(tmp_path, judge, contig, [x, s], reads)[0]


@pytest.mark.parametrize("judge", PATHS)
@pytest.mark.parametrize(
    "ref_len, alt_len, at",
    [(60, 30, 40), (50, 2, 25)],
    ids=["60>30", "50>2 (the bar)"],
)
def test_a_substitution_in_the_core_of_a_long_delins_is_a_context(
    tmp_path, judge, ref_len, alt_len, at
):
    """The REF bases replaced number 50 or more and S lies 4 or more bases in from both ends:
    S changes bases only REF's molecules hold, so the 6 reads that reach and spell it are
    REF, as at a deletion of 50 bp or more (RJ-24)."""
    assert _call(_core(tmp_path, judge, ref_len, alt_len, at)) == (36, 20, 0)


@pytest.mark.parametrize("judge", PATHS)
@pytest.mark.parametrize(
    "ref_len, alt_len, at, why",
    [
        (49, 2, 24, "under the bar"),
        (5, 60, 2, "5>60: five REF bases replaced"),
        (60, 30, 2, "S within 4 bases of the end: at a junction"),
    ],
    ids=["49>2", "5>60", "60>30 at a junction"],
)
def test_outside_the_core_a_listed_substitution_still_makes_its_readers_partial(
    tmp_path, judge, ref_len, alt_len, at, why
):
    """No core (under 50 REF bases replaced, or S at a junction): the 6 readers of S are S's,
    partial evidence, as before."""
    assert _call(_core(tmp_path, judge, ref_len, alt_len, at)) == (30, 20, 6), why
