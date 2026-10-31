"""A long deletion beside a run of its own unit: another unit count present at the site
holds the run (#297, RJ-24 amended 2026-10-30, the narrow rule).

Where the bases past a deletion of 50 bp or more repeat with the period of another allele
at its breakpoint (one more or one fewer unit of a tandem run whose deleted copy differs
from it, so the windows stop at the copy's first differing base), that allele reads like
the given one over the junction windows and on until the haplotypes first differ. A read
ending sooner fits either. The rule holds such a read to that base only where the other
allele is present at the site: listed in the input (always, whatever the reads show: rule X,
operator 2026-10-31; the input is taken as correct), or, unlisted, spelled in full by its own
bases on 3 or more fragments and 5% of those spelling the ALT there (RJ-25's floor). A read fitting the
ALT that does not tell it from that allele, by a confident base where the two differ,
is the ALT or that allele: neither, with partial evidence (RJ-25). A site with no such
allele keeps its reads as before. REF is told from both at the junction, so its reads
are unchanged, and the sibling claim keeps REF's window (RJ-24). Where a read of the
library's length cannot reach the first difference, nothing is held: the counts stand
and the row is named, UNIT_COUNT_UNSPANNED(change:bases/longest read).

Shapes (synthetic, PHI-free, paired inward, 6 reads per class unless named):
REF = flank, u1, the run of `m` copies of u, flank; u1 is u but for its third base.
The row deletes u1 (one unit); "two" deletes u1 and the run's first u. "short" reads hold
50 bases past the junction (100 bp from 50 before it); "holding" reads 251 bp from 20
before it. Counts per row: (ref_count, alt_count, partial_alt), with total_count where it
matters. Every case runs on the rules path and under the sequence judge. The last section
pins rule X on #297's S1 example (no molecule carries the listed count), through the engine
and through the CLI.
"""

import glob
import random

import pytest
from helpers import count_checked, make_read, read_maf_output, write_contig
from typer.testing import CliRunner

from gbcms import _rs
from gbcms.cli import app

JUDGES = ["rules", "sequence"]
SHORT, HOLDING = (100, 50), (251, 20)


def _other(b):
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[b]


def _shape(n, m, own=None, seed=297):
    """REF and the haplotypes (REF, "one": the row's deletion, "two": one more unit) for a
    unit of `n` bases; the deleted copy (`own` units of it, default one) differs from the
    unit at its third base. Returns (REF, the anchor's next base `a`, haps, rows)."""
    rng = random.Random(seed * 1000 + n * 10 + m)
    while True:
        u = "".join(rng.choice("ACGT") for _ in range(n))
        left = "".join(rng.choice("ACGT") for _ in range(400))
        right = "".join(rng.choice("ACGT") for _ in range(400))
        if (n < 3 or len(set(u)) > 2) and right[0] != u[0] and left[-1] != u[-1]:
            break
    own = own or 1
    deleted = list(u * own)
    deleted[2] = _other(deleted[2])
    deleted = "".join(deleted)
    a = len(left)
    ref = left + deleted + u * m + right
    k = len(deleted)
    haps = {"REF": ref, "one": ref[:a] + ref[a + k :], "two": ref[:a] + ref[a + k + n :]}
    rows = {
        "one": (a - 1, ref[a - 1 : a + k], ref[a - 1]),
        "two": (a - 1, ref[a - 1 : a + k + n], ref[a - 1]),
    }
    return ref, a, haps, rows


def _reads(a, name, hap, shape, *, size, k=6, tag="", clip=False):
    """`k` reads of `hap` (`shape`: length, bases before the junction `a`), written as an
    aligner would: REF ungapped; a carrier with its deletion of `size` bases at the junction
    (or soft-clipped there, `clip`) where it holds it, placed past it otherwise."""
    rl, before = shape
    out = []
    for i in range(k):
        s = a - before - i
        seq = hap[s : s + rl]
        if name == "REF" or s + rl <= a:
            pos, cig = s, ((0, rl),)
        elif s < a and clip:
            pos, cig = s, ((0, a - s), (4, rl - (a - s)))
        elif s < a:
            pos, cig = s, ((0, a - s), (2, size), (0, rl - (a - s)))
        else:
            pos, cig = s + size, ((0, rl),)
        r = make_read(f"{name}{tag}{i}", seq, pos, cig)
        r.flag = 0x1 | 0x2 | 0x20 | 0x40
        r.next_reference_id, r.next_reference_start, r.template_length = 0, pos + 250, 450
        out.append(r)
    return out


def _kind(ref_allele, alt_allele):
    """The variant kind of a row, from its alleles."""
    if len(ref_allele) == len(alt_allele):
        return "SNP"
    return "DEL" if len(ref_allele) > len(alt_allele) else "INS"


def _invariants(c):
    """The four counting invariants (AGENTS.md), asserted on one row's counts."""
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev


def _count(tmp_path, ref, rows, reads, judge, listed=False):
    """Counts per row, each with its UNIT_COUNT_UNSPANNED field; `listed`: the rows are
    co-annotated siblings of one another."""
    fa, bam = write_contig(tmp_path, ref, reads, "c")
    vs = [_rs.Variant("1", p, r, al, _kind(r, al)) for (p, r, al) in rows]
    pvs = _rs.prepare_variants(vs, fa, 5, False, 1, True, read_judge=judge)
    variants = [pv.variant for pv in pvs]
    kw = {"reference_fasta": fa} if judge == "sequence" else {}
    if listed:
        kw["sibling_variants"] = [[w for w in variants if w is not v] for v in variants]
    got = count_checked(bam, variants, read_judge=judge, **kw)
    for c in got:
        _invariants(c)
    return [((c.rd, c.ad, c.partial_alt, c.dp), c.unit_count_unspanned) for c in got]


def _ticket(n=60, m=2, clip=False):
    """The ticket's reads: 6 REF holding, 6 of two units and 6 of the given one, short
    (`clip`: the two-unit carriers soft-clipped at the junction, as BWA writes many)."""
    ref, a, h, rows = _shape(n, m)
    reads = _reads(a, "REF", h["REF"], HOLDING, size=0)
    reads += _reads(a, "two", h["two"], SHORT, size=2 * n, clip=clip)
    reads += _reads(a, "one", h["one"], SHORT, size=n)
    return ref, a, h, rows, reads


# ── #297's table, with and without the other unit count ────────────────────────────


@pytest.mark.parametrize("judge", JUDGES)
def test_no_other_count_keeps_the_reads(tmp_path, judge):
    """The given count alone, short of the run and holding it: every carrier is ALT, as
    before (no hold where no other count is at the site)."""
    ref, a, h, rows = _shape(60, 2)
    reads = _reads(a, "REF", h["REF"], HOLDING, size=0)
    reads += _reads(a, "one", h["one"], SHORT, size=60)
    reads += _reads(a, "one", h["one"], HOLDING, size=60, k=4, tag="h")
    assert _count(tmp_path, ref, [rows["one"]], reads, judge) == [((6, 10, 0, 16), "")]


@pytest.mark.parametrize("judge", JUDGES)
def test_short_carriers_of_another_count_alone_do_not_hold_the_run(tmp_path, judge):
    """#297's table as filed: no read spells the two-unit deletion in full (all end inside
    the run), so it is not found at the site and the reads count as before: 6/12/0."""
    ref, a, h, rows, reads = _ticket()
    assert _count(tmp_path, ref, [rows["one"]], reads, judge) == [((6, 12, 0, 18), "")]


@pytest.mark.parametrize("judge", JUDGES)
@pytest.mark.parametrize("clip", [False, True], ids=["written", "clipped"])
@pytest.mark.parametrize("given", [0, 4], ids=["alone", "beside-holding-carriers"])
def test_another_count_present_holds_the_run(tmp_path, judge, clip, given):
    """#297's table with the two-unit deletion present: 4 of its carriers hold the run (an
    aligner writes every two-unit carrier as the deletion, or soft-clips each at the
    junction: then no read writes it, and only their bases show it). The short carriers of
    either count fit both: partial evidence, never ALT (RJ-25). Carriers of the given count
    holding the run count ALT; the two-unit ones holding it carry another allele. REF reads
    are unchanged."""
    ref, a, h, rows, reads = _ticket(clip=clip)
    reads += _reads(a, "two", h["two"], HOLDING, size=120, k=4, tag="h", clip=clip)
    reads += _reads(a, "one", h["one"], HOLDING, size=60, k=given, tag="h")
    got = _count(tmp_path, ref, [rows["one"]], reads, judge)
    assert got == [((6, given, 16, 22 + given), "")]


@pytest.mark.parametrize("judge", JUDGES)
def test_one_fewer_unit_present_holds_the_run(tmp_path, judge):
    """The mirror: the row deletes two units (the differing copy and the run's first), and
    the one-unit deletion is present (4 carriers holding the run of three copies). The two
    read alike for two copies; the short carriers of either fit both: partial; carriers
    holding the run count for the count they show."""
    ref, a, h, rows = _shape(60, 3)
    reads = _reads(a, "REF", h["REF"], HOLDING, size=0)
    reads += _reads(a, "two", h["two"], SHORT, size=120)
    reads += _reads(a, "one", h["one"], SHORT, size=60)
    reads += _reads(a, "one", h["one"], HOLDING, size=60, k=4, tag="h")
    reads += _reads(a, "two", h["two"], HOLDING, size=120, k=4, tag="h")
    got = _count(tmp_path, ref, [rows["two"]], reads, judge)
    assert got == [((6, 4, 16, 26), "")]


@pytest.mark.parametrize("judge", JUDGES)
def test_an_inserted_unit_is_told_at_the_junction(tmp_path, judge):
    """An allele inserting one unit keeps the differing copy, so the junction tells it from
    the one-unit deletion: with 4 of its carriers holding the run, nothing is held and the
    short carriers of the given count stay ALT."""
    ref, a, h, rows = _shape(60, 2)
    plus = ref[: a + 60] + ref[a + 60 : a + 120] + ref[a + 60 :]
    reads = _reads(a, "REF", h["REF"], HOLDING, size=0)
    reads += _reads(a, "one", h["one"], SHORT, size=60)
    for i in range(4):
        s = a - 20 - i
        r = make_read(f"plus{i}", plus[s : s + 251], s, ((0, 80 + i), (1, 60), (0, 111 - i)))
        r.flag = 0x1 | 0x2 | 0x20 | 0x40
        r.next_reference_id, r.next_reference_start, r.template_length = 0, s + 250, 450
        reads.append(r)
    (((rd, ad, _, _), diag),) = _count(tmp_path, ref, [rows["one"]], reads, judge)
    assert (rd, ad, diag) == (6, 6, "")


@pytest.mark.parametrize("judge", JUDGES)
@pytest.mark.parametrize("holding", [1, 2, 3])
def test_the_other_count_is_found_on_three_fragments(tmp_path, judge, holding):
    """Found as RJ-25 finds another allele: on 3 or more fragments. With three, the short
    carriers of either count (they cannot tell the counts apart) and the three that hold the
    run are partial evidence: (6, 0, 15). With fewer the count is not at the site, and the
    floor decides only the reads that cannot tell: the 12 short carriers stay ALT. A read
    whose own bases spell the other count is no ALT whatever the floor (F3, RJ-6): partial
    evidence on the rules path (the judge, where an allele under the floor is no state,
    leaves it depth only). Before F3 the rules path counted the holders ALT: (6, 14) for two."""
    ref, a, h, rows, reads = _ticket()
    reads += _reads(a, "two", h["two"], HOLDING, size=120, k=holding, tag="h")
    (((rd, ad, partial, _), _),) = _count(tmp_path, ref, [rows["one"]], reads, judge)
    if holding >= 3:
        assert (rd, ad, partial) == (6, 0, 12 + holding)
    else:
        assert (rd, ad) == (6, 12)
        assert partial == (holding if judge == "rules" else 0)


@pytest.mark.parametrize("judge", JUDGES)
def test_the_other_count_is_found_on_five_percent_of_the_alt_fragments(tmp_path, judge):
    """And on 5% of the fragments spelling the ALT there: 3 two-unit carriers holding the
    run beside 80 of the given count holding it (3.75%) leave the site as before; beside
    40 (7.5%) they hold it."""
    ref, a, h, rows = _shape(60, 2)
    base = _reads(a, "REF", h["REF"], HOLDING, size=0)
    base += _reads(a, "one", h["one"], SHORT, size=60)
    base += _reads(a, "two", h["two"], HOLDING, size=120, k=3, tag="h")
    for many, short_alt in ((80, True), (40, False)):
        d = tmp_path / str(many)
        d.mkdir()
        reads = base + _reads(a, "one", h["one"], HOLDING, size=60, k=many, tag="m")
        (((rd, ad, partial, _), _),) = _count(d, ref, [rows["one"]], reads, judge)
        assert rd == 6
        if short_alt:
            # the 3 two-unit carriers spell the count under the floor: no ALT (F3), partial
            # on the rules path, depth only under the judge; the short carriers stay ALT
            assert ad == many + 6, (many, ad)
            assert partial == (3 if judge == "rules" else 0), (many, partial)
        else:
            assert ad == many, (many, ad)


@pytest.mark.parametrize("judge", JUDGES)
def test_a_listed_other_count_holds_the_run(tmp_path, judge):
    """The two-unit deletion listed beside the row: the short carriers fit both rows, so
    neither row's ALT (partial evidence at each); carriers holding the run count for the
    row whose count they show (RJ-5, RJ-25)."""
    ref, a, h, rows, reads = _ticket()
    reads += _reads(a, "two", h["two"], HOLDING, size=120, k=4, tag="h")
    reads += _reads(a, "one", h["one"], HOLDING, size=60, k=4, tag="h")
    got = _count(tmp_path, ref, [rows["one"], rows["two"]], reads, judge, listed=True)
    assert [g[0][:3] for g in got] == [(6, 4, 16), (6, 4, 16)]


# ── A run longer than the read ─────────────────────────────────────────────────────

LONG = (150, 20)  # every read 150 bp: none reaches the end of a four-unit run


def _long_run(listed):
    ref, a, h, rows = _shape(60, 4)
    reads = _reads(a, "REF", h["REF"], LONG, size=0)
    reads += _reads(a, "one", h["one"], LONG, size=60)
    reads += _reads(a, "two", h["two"], LONG, size=120)
    return ref, ([rows["one"], rows["two"]] if listed else [rows["one"]]), reads


@pytest.mark.parametrize("judge", JUDGES)
def test_a_run_longer_than_the_read_keeps_the_counts_and_names_the_row(tmp_path, judge):
    """Four units after the deleted copy: one more unit deleted reads like the given
    deletion for 184 bases from the windows' first base, and a 150-bp read cannot hold
    the 208 it takes to spell the count in full (184 and the 24 past the first difference,
    F7; 188 before). Nothing is held (every carrier ALT, as before) and the row is named."""
    ref, rows, reads = _long_run(False)
    assert _count(tmp_path, ref, rows, reads, judge) == [((6, 12, 0, 18), "-120:208/150")]


@pytest.mark.parametrize("judge", JUDGES)
def test_a_listed_count_past_the_read_holds_the_reads_and_names_the_rows(tmp_path, judge):
    """The two-unit deletion listed too, past what a read holds: no read can tell the given
    count from the listed one, so each carrier of either is partial evidence at each row, on
    both paths, and both rows are named, each with the other's count (F2: the input is taken
    as correct, whatever the reads reach). Before, the rules path counted every carrier ALT
    at the one-unit row, (6, 12, 0), and at the two-unit row left the one-unit carriers
    partial, (6, 6, 6); the sequence judge held both rows, (6, 0, 12)."""
    ref, rows, reads = _long_run(True)
    got = _count(tmp_path, ref, rows, reads, judge, listed=True)
    assert got == [((6, 0, 12, 18), "-120:208/150"), ((6, 0, 12, 18), "-60:208/150")]


# ── The sibling claim keeps REF's window ───────────────────────────────────────────


@pytest.mark.parametrize("judge", JUDGES)
def test_the_sibling_claim_keeps_refs_window_beside_another_count(tmp_path, judge):
    """A co-annotated SNV 20 bases into the deleted copy, on every REF read, with the
    two-unit deletion present (4 carriers holding the run): the run held for the ALT
    leaves REF's window where it was, so the SNV, past it, claims no REF read (RJ-24);
    REF stays 6 (the first build grew the window and the claim with it)."""
    ref, a, h, rows, reads = _ticket()
    reads += _reads(a, "two", h["two"], HOLDING, size=120, k=4, tag="h")
    g = a + 20
    snv_base = _other(ref[g])
    with_snv = ref[:g] + snv_base + ref[g + 1 :]
    reads = [r for r in reads if not r.query_name.startswith("REF")]
    reads += _reads(a, "REF", with_snv, HOLDING, size=0)
    fa, bam = write_contig(tmp_path, ref, reads, "c")
    vs = [_rs.Variant("1", *rows["one"], "DEL"), _rs.Variant("1", g, ref[g], snv_base, "SNP")]
    deletion, snv = (
        pv.variant for pv in _rs.prepare_variants(vs, fa, 5, False, 1, True, read_judge=judge)
    )
    kw = {"reference_fasta": fa} if judge == "sequence" else {}
    got, _ = count_checked(
        bam, [deletion, snv], read_judge=judge, sibling_variants=[[snv], [deletion]], **kw
    )
    assert (got.rd, got.ad, got.partial_alt) == (6, 0, 16)


# ── Short units: the first difference, not a run end, sets the hold ────────────────


def _short_unit(n, run, own, seed=11):
    """A deletion of `own` copies of an `n`-base unit (the copies' third base changed),
    before a run of `run` exact copies: one more unit deleted reads like the given
    deletion until the run's last unit."""
    return _shape(n, run, own=own, seed=seed)


SMALL = [
    pytest.param(1, 30, 50, id="unit-1"),
    pytest.param(2, 20, 25, id="unit-2"),
    pytest.param(4, 10, 13, id="unit-4"),
    pytest.param(5, 8, 10, id="unit-5"),
]


@pytest.mark.parametrize("judge", JUDGES)
@pytest.mark.parametrize("n, run, own", SMALL)
@pytest.mark.parametrize("present", [False, True], ids=["alone", "present"])
def test_short_units(tmp_path, judge, n, run, own, present):
    """A 50-to-52-bp deletion of copies of a 1-, 2-, 4- or 5-base unit before a run of
    it. Short carriers hold 15 bases past the junction, past the windows but short of
    where one more unit deleted first differs. With that allele present (4 carriers holding
    the run), they fit both: partial; alone, ALT as before."""
    ref, a, h, rows = _short_unit(n, run, own)
    size = n * own
    short = (60, 45)  # 15 bases past the junction
    reads = _reads(a, "REF", h["REF"], HOLDING, size=0)
    reads += _reads(a, "one", h["one"], short, size=size)
    reads += _reads(a, "one", h["one"], HOLDING, size=size, k=4, tag="h")
    if present:
        reads += _reads(a, "two", h["two"], HOLDING, size=size + n, k=4, tag="h")
    (((rd, ad, partial, _), _),) = _count(tmp_path, ref, [rows["one"]], reads, judge)
    assert (rd, ad, partial) == ((6, 4, 10) if present else (6, 10, 0))


# ── Errors over a long stretch (RJ-26) ─────────────────────────────────────────────


@pytest.mark.parametrize("judge", JUDGES)
@pytest.mark.parametrize(
    "error, alt",
    [
        ("none", 4),
        ("in the run", 4),  # where both counts read alike: costs both, decides nothing
        ("a third base past the difference", 4),  # neither count's: an error
        ("the difference unread", 0),  # every telling base N (unread): a tie
    ],
)
def test_errors_over_the_held_run(tmp_path, judge, error, alt):
    """The two-unit deletion present (4 carriers holding the run); 4 carriers of the given
    count hold it with one error. Only the bases where the two counts differ tell them:
    an error in the run both read alike, or a third base where they differ, leaves the
    call ALT; with those bases unread (N, as a duplex disagreement is masked) the read fits
    both."""
    ref, a, h, rows, reads = _ticket()
    reads += _reads(a, "two", h["two"], HOLDING, size=120, k=4, tag="h")
    one = h["one"]
    first = next(i for i in range(a, len(one)) if one[i] != h["two"][i])

    hap = one
    if error == "the difference unread":
        hap = one[:first] + "N" * (len(one) - first)
    elif error == "in the run":
        hap = one[: a + 30] + _other(one[a + 30]) + one[a + 31 :]
    elif error == "a third base past the difference":
        third = next(b for b in "ACGT" if b not in (one[first + 1], h["two"][first + 1]))
        hap = one[: first + 1] + third + one[first + 2 :]
    reads += _reads(a, "one", hap, HOLDING, size=60, k=4, tag="h")
    (((rd, ad, _, _), _),) = _count(tmp_path, ref, [rows["one"]], reads, judge)
    assert (rd, ad) == (6, alt)


# ── A listed neighbour is no other count ───────────────────────────────────────────


def _neighbour_shape(seed=5):
    """The ALT continues T: 12-base units for 52 bases, an A, then 13-base units. One more
    unit deleted reads T from its 12th base, which from T's 40th base on is the ALT with an A
    inserted after its 39th: the carriers of the deletion with that phased insertion read as
    one more unit deleted (HG002 4:3610329's shape)."""
    rng = random.Random(seed)
    cgt = lambda n: "".join(rng.choice("CGT") for _ in range(n))  # noqa: E731
    left, right = cgt(400), cgt(400)
    deleted = "A" + cgt(59)
    if left[-1] == deleted[-1]:
        left = left[:-1] + ("C" if deleted[-1] != "C" else "G")
    t = list(cgt(12))
    while len(t) < 240:
        i = len(t)
        t.append(t[i - 12] if i < 52 else "A" if i == 52 else t[i - 13])
    t = "".join(t)
    a = len(left)
    ref = left + deleted + t + right
    carrier = left + t[:40] + "A" + t[40:] + right
    deletion = (a - 1, ref[a - 1 : a + 60], ref[a - 1])
    g = a + 60 + 39  # the insertion's anchor, on REF past the deleted bases
    insertion = (g, ref[g], ref[g] + "A")
    return ref, a, carrier, deletion, insertion


@pytest.mark.parametrize("judge", JUDGES)
@pytest.mark.parametrize("listed", [True, False], ids=["insertion-listed", "not-listed"])
def test_a_listed_neighbour_makes_no_other_count(tmp_path, judge, listed):
    """Carriers of the deletion with a phased 1-bp insertion past the junction, which the
    input lists: they are the deletion with a listed neighbour (RJ-24), never one more unit
    deleted, so they count ALT. Not listed, their bases are one more unit deleted's over most
    of what they hold: that allele is found on them and held (on the rules path every carrier
    is partial evidence; the sequence judge, weighing each read whole with its own insertion,
    counts ALT the 3 that read on past where the two part, near the run's end)."""
    ref, a, carrier, deletion, insertion = _neighbour_shape()
    reads = []
    for i in range(6):
        s = a - 20 - i
        r = make_read(f"REF{i}", ref[s : s + 251], s, ((0, 251),))
        c = make_read(
            f"alt{i}",
            carrier[s : s + 251],
            s,
            ((0, 20 + i), (2, 60), (0, 40), (1, 1), (0, 190 - i)),
        )
        for x in (r, c):
            x.flag = 0x1 | 0x2 | 0x20 | 0x40
            x.next_reference_id, x.next_reference_start, x.template_length = (
                0,
                x.reference_start + 250,
                450,
            )
        reads += [r, c]
    fa, bam = write_contig(tmp_path, ref, reads, "c")
    kinds = [("DEL", deletion)] + ([("INS", insertion)] if listed else [])
    vs = [_rs.Variant("1", p, r, al, kind) for kind, (p, r, al) in kinds]
    pvs = _rs.prepare_variants(vs, fa, 5, False, 1, True, read_judge=judge)
    kw = {"reference_fasta": fa} if judge == "sequence" else {}
    got = count_checked(bam, [pv.variant for pv in pvs], read_judge=judge, **kw)[0]
    assert got.dp >= got.rd + got.ad and got.rd == got.rd_fwd + got.rd_rev
    assert (got.rd, got.ad) == ((6, 6) if listed else (6, 0) if judge == "rules" else (6, 3))


# ── Edge cases: a reference N in the run; the named AD-low limit ───────────────────


@pytest.mark.parametrize("judge", JUDGES)
def test_a_reference_n_inside_the_run_tells_nothing(tmp_path, judge):
    """The reference holds an N 30 bases into the run (the reads read a base there): an N
    tells no allele from another, so one more unit deleted still first differs at the run's
    copy end, and the counts are those with no N: (6, 0, 16)."""
    ref, a, h, rows, reads = _ticket()
    reads += _reads(a, "two", h["two"], HOLDING, size=120, k=4, tag="h")
    g = a + 60 + 30
    with_n = ref[:g] + "N" + ref[g + 1 :]
    assert _count(tmp_path, with_n, [rows["one"]], reads, judge) == [((6, 0, 16, 22), "")]


@pytest.mark.parametrize("judge", JUDGES)
def test_ad_reads_low_beside_a_present_count(tmp_path, judge):
    """Named limit (RJ-24, measured on the spike-in arm, D12 phase 1b): REF keeps its
    windows while an ALT read must reach the first difference. One REF read and one carrier
    of the given count start at every base from 144 to 4 before the junction (150 bp each;
    REF holds its windows from all 141), and the two-unit deletion is listed: only the 86
    carriers starting 89 bases before the junction or later reach where the counts differ,
    so AD/(RD+AD) reads 86/227 = 0.379 for a true 0.5."""
    ref, a, h, rows = _shape(60, 2)
    reads = []
    for i, s in enumerate(range(a - 144, a - 3)):
        for name, hap in (("REF", h["REF"]), ("one", h["one"])):
            cig = ((0, 150),) if name == "REF" else ((0, a - s), (2, 60), (0, 150 - (a - s)))
            r = make_read(f"{name}{i}", hap[s : s + 150], s, cig)
            r.flag = 0x1 | 0x2 | 0x20 | 0x40
            r.next_reference_id, r.next_reference_start, r.template_length = 0, s + 250, 450
            reads.append(r)
    got = _count(tmp_path, ref, [rows["one"], rows["two"]], reads, judge, listed=True)
    rd, ad = got[0][0][:2]
    assert (rd, ad) == (141, 86)


@pytest.mark.parametrize("judge", JUDGES)
def test_reads_the_aligner_places_elsewhere_spell_no_count(tmp_path, judge):
    """RJ-29's placement test: the 4 two-unit carriers holding the run (soft-clipped, so no
    read writes that deletion) are placed elsewhere as well by their aligner (XS equal to
    AS): they may come from another copy, so they do not make the two-unit deletion present,
    and the site keeps its reads as before (the 12 short carriers ALT; the holding
    two-unit carriers depth only, on both paths)."""
    ref, a, h, rows, reads = _ticket(clip=True)
    doubted = _reads(a, "two", h["two"], HOLDING, size=120, k=4, tag="h", clip=True)
    for r in doubted:
        r.set_tag("AS", 230)
        r.set_tag("XS", 230)
    got = _count(tmp_path, ref, [rows["one"]], reads + doubted, judge)
    # the doubted holders spell the two-unit count, so they are no ALT (F3), and as the aligner
    # doubts their placement they are depth only on both paths (before F3 the rules counted
    # them ALT: AD 16 against the judge's 12)
    assert got == [((6, 12, 0, 22), "")]


@pytest.mark.parametrize("judge", JUDGES)
@pytest.mark.parametrize("given", [0, 4], ids=["alone", "beside-holding-carriers"])
def test_the_count_that_holds_the_calls_is_named(tmp_path, judge, given):
    """The two-unit deletion found on its 4 holding carriers takes the ALT calls of the 16
    carriers that do not tell the given count from it, and is named for the row (the
    OBSERVED_ALLELE form, at the deletion's anchor: 16 fragments whose ALT call it takes, and
    0 or 4 that tell the given count from it), on both paths; a site where it is absent names
    nothing (the no-other-count test)."""
    ref, a, h, rows, reads = _ticket()
    reads += _reads(a, "two", h["two"], HOLDING, size=120, k=4, tag="h")
    reads += _reads(a, "one", h["one"], HOLDING, size=60, k=given, tag="h")
    fa, bam = write_contig(tmp_path, ref, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", *rows["one"], "DEL")], fa, 5, False, 1, True, read_judge=judge
    )
    kw = {"reference_fasta": fa} if judge == "sequence" else {}
    (c,) = count_checked(bam, [pv.variant], read_judge=judge, **kw)
    p, r, al = rows["two"]
    assert c.unit_count_allele == f"{p + 1}:{r}>{al}:16/{given}"


# ── Rule X (operator, 2026-10-31): a listed other count is always a candidate ──────


def _counts(tmp_path, ref, rows, reads, judge, listed=False):
    """The `BaseCounts` of each row, the four counting invariants asserted on every one."""
    fa, bam = write_contig(tmp_path, ref, reads, "c")
    vs = [_rs.Variant("1", p, r, al, _kind(r, al)) for (p, r, al) in rows]
    pvs = _rs.prepare_variants(vs, fa, 5, False, 1, True, read_judge=judge)
    variants = [pv.variant for pv in pvs]
    kw = {"reference_fasta": fa} if judge == "sequence" else {}
    if listed:
        kw["sibling_variants"] = [[w for w in variants if w is not v] for v in variants]
    got = count_checked(bam, variants, read_judge=judge, **kw)
    for c in got:
        _invariants(c)
    return got


def _s1_example(listed):
    """#297's S1 worked example: the given deletion (one unit) alone in the molecules. 6 REF
    reads, 6 short carriers of the given count (they end inside the run), and 4 carriers of
    it that read through the run. No molecule carries the two-unit deletion, which the input
    lists (`listed`) or leaves out."""
    ref, a, h, rows = _shape(60, 2)
    reads = _reads(a, "REF", h["REF"], HOLDING, size=0)
    reads += _reads(a, "one", h["one"], SHORT, size=60)
    reads += _reads(a, "one", h["one"], HOLDING, size=60, k=4, tag="h")
    return ref, ([rows["one"], rows["two"]] if listed else [rows["one"]]), reads, rows


@pytest.mark.parametrize("judge", JUDGES)
def test_a_listed_count_no_molecule_carries_still_holds_the_short_carriers(tmp_path, judge):
    """Rule X: the input lists the two-unit deletion, so it is an allele that may be at the
    site (invariant 7: the input is taken as correct), whether or not a read shows it. The 6
    short carriers of the given count end inside the run, before the base where the two
    deletions first differ: each fits either, so each is partial evidence, never ALT. The 4
    carriers that read through the run tell the given count: ALT. (RD, AD, partial) =
    (6, 4, 6), on the rules path and under the sequence judge. The listed count is named for
    the row: COEXISTING_ALLELE form, 6 fragments whose ALT call it takes, 4 that tell the given
    count from it."""
    ref, rows_in, reads, rows = _s1_example(listed=True)
    got = _counts(tmp_path, ref, rows_in, reads, judge, listed=True)
    c = got[0]
    assert (c.rd, c.ad, c.partial_alt) == (6, 4, 6)
    p, r, al = rows["two"]
    assert c.unit_count_allele == f"{p + 1}:{r}>{al}:6/4"


@pytest.mark.parametrize("judge", JUDGES)
def test_the_same_count_unlisted_and_unshown_holds_nothing(tmp_path, judge):
    """The same reads with the two-unit deletion left out of the input: no row lists it and
    no read shows it (RJ-25's floor: spelled in full on 3 or more fragments), so it is not
    at the site and every carrier of the given count is ALT: (6, 10, 0), nothing named."""
    ref, rows_in, reads, _ = _s1_example(listed=False)
    (c,) = _counts(tmp_path, ref, rows_in, reads, judge)
    assert (c.rd, c.ad, c.partial_alt) == (6, 10, 0)
    assert c.unit_count_allele == ""


@pytest.mark.parametrize("judge", JUDGES)
@pytest.mark.parametrize("listed", [True, False], ids=["listed", "unlisted"])
def test_the_s1_example_through_the_cli_names_the_listed_count(tmp_path, judge, listed):
    """The same reads through `gbcms dna` (a VCF of the deletion, and of the two-unit deletion
    when listed): the MAF row of the given deletion has ref_count, alt_count and partial_alt
    (6, 4, 6) and `gbcms_diagnostic` names the listed count, COEXISTING_ALLELE at the
    deletion's anchor (6 fragments whose ALT call it takes, 4 that tell the given count from
    it). Unlisted: (6, 10, 0) and no allele named."""
    ref, rows_in, reads, rows = _s1_example(listed)
    fa, bam = write_contig(tmp_path, ref, reads, "c")
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        f"##fileformat=VCFv4.2\n##contig=<ID=1,length={len(ref)}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(f"1\t{p + 1}\t.\t{r}\t{al}\t.\t.\t.\n" for p, r, al in rows_in)
    )
    out = tmp_path / "out"
    res = CliRunner().invoke(
        app,
        ["dna", "-v", str(vcf), "-b", f"S:{bam}", "-f", fa, "-o", str(out)]
        + ["--format", "maf", "--read-judge", judge],
    )
    assert res.exit_code == 0, res.output
    p, r, al = rows["one"]
    (row,) = [
        x
        for x in read_maf_output(glob.glob(str(out / "*.maf"))[0])
        if x["vcf_pos"] == str(p + 1) and x["vcf_ref"] == r
    ]
    counts = (row["ref_count"], row["alt_count"], row["partial_alt"])
    flags = row["gbcms_diagnostic"].split(";")
    if listed:
        p2, r2, al2 = rows["two"]
        assert counts == ("6", "4", "6")
        assert f"COEXISTING_ALLELE(1:{p2 + 1}:{r2}>{al2}:6/4)" in flags
    else:
        assert counts == ("6", "10", "0")
        assert not any(f.startswith(("OBSERVED_ALLELE", "COEXISTING_ALLELE")) for f in flags)


# ── Closing the rule on one principle (#297 close): a read is ALT only where its own ──
# ── confident bases tell the given ALT from every candidate ──────────────────────────


def _tie_carriers(a, h, mode, k=6):
    """`k` carriers of the given count that cannot tell it from the listed two-unit count:
    251 bp from 20 bases before the junction, written with the deletion. Past the junction
    their bases are read at Q5 (unread) but for two bases. Mode "A": the bases where the two
    counts agree are read confidently and every base where they differ is unread. Mode "B":
    one telling base shows the given count's base and one shows the other count's, so one
    base against one: a tie."""
    out = []
    for i in range(k):
        s = a - 20 - i
        seq = list(h["one"][s : s + 251])
        qual = [30] * 251
        start = (a - s) + 60  # the first base past the junction, in read coordinates
        diffs = [x for x in range(start, 251) if h["one"][s + x] != h["two"][s + x]]
        for x in range(start, 251):
            agree = h["one"][s + x] == h["two"][s + x]
            qual[x] = 30 if (mode == "A" and agree) else 5
        if mode == "B":
            qual[diffs[0]] = 30  # one telling base, the given count's
            seq[diffs[1]] = h["two"][s + diffs[1]]
            qual[diffs[1]] = 30  # one telling base, the other count's
        cig = ((0, a - s), (2, 60), (0, 251 - (a - s)))
        r = make_read(f"one{mode}{i}", "".join(seq), s, cig, quals=qual)
        r.flag = 0x1 | 0x2 | 0x20 | 0x40
        r.next_reference_id, r.next_reference_start, r.template_length = 0, s + 250, 450
        out.append(r)
    return out


def _unrelated_row(kind, ref, a):
    """One input row that no read carries and that has nothing to do with the two counts:
    a 1-bp deletion 270 bases past the junction (past what a read holds), or an SNV 10 bases
    into the run."""
    if kind == "deletion":
        far = a + 60 + 120 + 150 - 1
        return (far, ref[far : far + 2], ref[far])
    return (a + 70, ref[a + 70], _other(ref[a + 70]))


@pytest.mark.parametrize("judge", JUDGES)
@pytest.mark.parametrize("mode", ["A", "B"], ids=["bases-agree", "one-against-one"])
@pytest.mark.parametrize("extra", ["none", "deletion", "snv"])
def test_an_unrelated_listed_row_makes_no_tie_read_the_deletion_with_a_neighbour(
    tmp_path, judge, mode, extra
):
    """F1. 6 REF reads and 6 carriers of the given count that cannot tell it from the listed
    two-unit count (S1's shape): each is partial evidence, (RD, AD, partial) = (6, 0, 6), on
    both paths. One more input row that no read carries (a 1-bp deletion past every read's
    reach, or an SNV 10 bases into the run) changes nothing: a read is "the deletion with its
    listed neighbour" only where its confident bases show the neighbour's change (RJ-24,
    RJ-28: a row's count does not depend on an unrelated input row). The rules path gave
    (6, 6, 0) with such a row: the neighbour's bases equal the ALT's over everything a read
    holds, so every tie read fitted it."""
    ref, a, h, rows = _shape(60, 2)
    reads = _reads(a, "REF", h["REF"], HOLDING, size=0) + _tie_carriers(a, h, mode)
    rws = [rows["one"], rows["two"]] + ([] if extra == "none" else [_unrelated_row(extra, ref, a)])
    fa, bam = write_contig(tmp_path, ref, reads, "c")
    vs = [_rs.Variant("1", p, r, al, _kind(r, al)) for (p, r, al) in rws]
    variants = [
        pv.variant for pv in _rs.prepare_variants(vs, fa, 5, False, 1, True, read_judge=judge)
    ]
    siblings = [[variants[1]], [variants[0]]] + [[] for _ in rws[2:]]
    kw = {"reference_fasta": fa} if judge == "sequence" else {}
    got = count_checked(bam, variants, read_judge=judge, sibling_variants=siblings, **kw)
    for c in got:
        _invariants(c)
    c = got[0]
    assert (c.rd, c.ad, c.partial_alt) == (6, 0, 6)


def _bases_next_to_other_indels(reads, first=12, period=8):
    """Give each read a 1-bp deletion in its alignment every `period` bases from `first` bases
    past its junction deletion on (its sequence is untouched): heuristic BAQ downgrades the
    bases within 5 of each, so past the junction none is readable at --min-baseq 20."""
    for r in reads:
        cig = list(r.cigartuples)
        k = max(i for i, (op, _) in enumerate(cig) if op == 2)
        left = cig[k + 1][1] - first
        new = cig[: k + 1] + [(0, first)]
        while left > 0:
            take = min(period, left)
            new.append((0, take))
            left -= take
            if left > 0:
                new.append((2, 1))
        r.cigartuples = new
    return reads


@pytest.mark.parametrize("judge", JUDGES)
@pytest.mark.parametrize("baq", [False, True], ids=["stored", "apply-baq"])
def test_the_presence_test_reads_the_qualities_the_hold_reads(tmp_path, judge, baq):
    """F6 (invariant 2: one quality contract). The two-unit deletion's 4 carriers holding the
    run each carry an alignment indel every 8 bases past the junction, so under `--apply-baq`
    no base they hold past it is readable (the classifiers read the downgraded qualities). The
    presence test read the stored qualities: 4 fragments spelled the count, it was present, and
    every ALT-fitting read was held, (RD, AD, partial) = (6, 0, 16). Read as the hold reads them,
    no carrier spells anything: the count is not at the site and the 16 reads that cannot tell
    stay ALT, (6, 16, 0). Stored qualities: (6, 0, 16), as before. Both paths."""
    ref, a, h, rows, reads = _ticket()
    reads += _bases_next_to_other_indels(
        _reads(a, "two", h["two"], HOLDING, size=120, k=4, tag="h")
    )
    fa, bam = write_contig(tmp_path, ref, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", *rows["one"], "DEL")], fa, 5, False, 1, True, read_judge=judge
    )
    kw = {"reference_fasta": fa} if judge == "sequence" else {}
    (c,) = count_checked(bam, [pv.variant], read_judge=judge, apply_baq=baq, **kw)
    _invariants(c)
    assert (c.rd, c.ad, c.partial_alt, c.dp) == ((6, 16, 0, 22) if baq else (6, 0, 16, 22))


@pytest.mark.parametrize("judge", JUDGES)
@pytest.mark.parametrize(
    "copies, named",
    [(3, ""), (4, "-35:168/150"), (5, "-105:168/150")],
    ids=["in-reach", "4-to-23-bases-short", "far"],
)
def test_an_unlisted_count_no_read_can_spell_in_full_is_named(tmp_path, judge, copies, named):
    """F7. A 70-base deletion (two copies of a 35-base unit) before a run of `copies` of the
    unit; every read 150 bp, none of the other count. The nearest other count first differs from
    the ALT 109, 144 or 179 bases from the windows' first base (the one fewer unit deleted, or,
    in the last run, the one more: the nearest on either side is named). A read spells an
    allele in full only 24 bases past its first difference, so with 144 (4 to 23 bases before
    the longest read's end) no read can: the reach test read the spelling's 4 bases, not its
    24, and the count was neither held nor named. It is named now, with the bases a read must
    hold (168). Counts unchanged; with 109 a read can spell it, so it is not named."""
    ref, a, h, rows = _shape(35, copies, own=2)
    reads = _reads(a, "REF", h["REF"], LONG, size=0) + _reads(a, "one", h["one"], LONG, size=70)
    assert _count(tmp_path, ref, [rows["one"]], reads, judge) == [((6, 6, 0, 12), named)]
