"""The sequence judge (`--read-judge sequence`, #304) at a change at a base the event keeps.

Mechanism. A substitution or delins X changes some bases of its span and keeps others (REF
and X agree there). A change at a kept base that the REF molecules carry too (a germline
SNP, the guard's background) is a separate event on both haplotypes, not another allele
(RJ-25, operator 2026-10-14): it explains a read's base there on REF and on X alike, so a
carrier holding it counts ALT and a REF molecule holding it counts REF, and reads it
cannot explain (a third base there that no REF molecule carries) are another allele.
Before this, the judge took X's span from its first to its last change as X's own bases:
a background change inside it was dropped as meeting X, a mask there was taken as X's
base, and an ALT call needed X's base there while a REF call needed none.

Shapes: MNPs of 2, 3 and 4 bases (the 2-base one keeps no base: its change beside the
block is the control), the change at each kept base; the change on the carriers (cis) or
not (trans); read on both strands, below min BQ, and by reads ending at the block's last
base or starting at its first; written plainly or with an indel at the block (the found
alleles' route); a delins whose carriers write the change (a found allele), its masked
reads judged by the deltas they add. Counts are under the judge, explicitly.
"""

import random

import pytest
from helpers import count_checked, make_read, write_contig

from gbcms import _rs

READ = 100
P = 300  # 0-based start of every event


def _contig(seed, n=800):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(n))


def _other(b, *avoid):
    return next(x for x in "ACGT" if x != b and x not in avoid)


def _invariants(c):
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev


def _reads(contig, ref_len, block, tag, n, *, low=None, ends=False, indel=False, first=0):
    """`n` reads carrying `block` in place of the `ref_len` REF bases at P, alternating
    strands. `low`: block offsets read below min BQ. `ends`: every third read ends on
    the block's last base and every third starts on its first. `indel`: the reads write
    the block as an insertion and a deletion (a same-length block only)."""
    out = []
    for i in range(first, first + n):
        if ends and i % 3 == 1:
            s = P + len(block) - READ  # ends on the block's last base
        elif ends and i % 3 == 2:
            s = P  # starts on its first
        else:
            s = P - 50 + (i * 7) % 40
        seq = (contig[s:P] + block + contig[P + ref_len : s + READ + ref_len])[:READ]
        left = P - s
        if len(block) != ref_len:
            k = len(block)
            cig = ((0, left), (1, k), (2, ref_len), (0, READ - left - k))
        elif indel and 0 < left and left + len(block) < READ:
            cig = ((0, left), (1, 1), (0, len(block) - 1), (2, 1), (0, READ - left - len(block)))
        else:
            cig = ((0, READ),)
        r = make_read(f"{tag}{i}", seq, s, cig, flag=16 if i % 2 else 0)
        if low is not None:
            q = list(r.query_qualities)
            for k in low:
                if 0 <= left + k < READ:
                    q[left + k] = 5
            r.query_qualities = q
        out.append(r)
    return out


def _count(tmp_path, contig, row, groups, name="s"):
    tmp_path.mkdir(parents=True, exist_ok=True)
    reads = [r for g in groups for r in g]
    fa, bam = write_contig(tmp_path, contig, reads, name)
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", P, row[0], row[1], "COMPLEX")],
        fa,
        5,
        False,
        1,
        True,
        read_judge="sequence",
    )
    assert pv.gbcms_status == "PASS", pv.gbcms_status_reason
    (c,) = count_checked(bam, [pv.variant], read_judge="sequence", reference_fasta=fa)
    _invariants(c)
    return c


def _with(block, k, base):
    return block[:k] + base + block[k + 1 :]


# REF, ALT, the offsets the MNP keeps
MNPS = {
    "3 bases": ("GCA", "ACT", (1,)),
    "4 bases": ("GCAT", "ACAA", (1, 2)),
}
CASES = [(name, k) for name, (_, _, kept) in MNPS.items() for k in kept]


def _mnp(seed, ref_block):
    c = list(_contig(seed))
    c[P - 1] = _other(ref_block[0])  # nothing at the block extends it
    c[P : P + len(ref_block)] = ref_block
    c[P + len(ref_block)] = _other(ref_block[-1])
    return "".join(c)


@pytest.mark.parametrize("ends", [False, True], ids=["mid-read", "read ends"])
@pytest.mark.parametrize("name, k", CASES, ids=[f"{n} kept {k}" for n, k in CASES])
def test_carriers_of_a_change_both_haplotypes_carry_at_a_kept_base_count_alt(
    tmp_path, name, k, ends
):
    """The change at the kept base on half the REF molecules and on every carrier (cis):
    a separate event, so the carriers are ALT and the REF molecules REF, either way."""
    ref_b, alt_b, _ = MNPS[name]
    contig = _mnp(11, ref_b)
    g = _other(ref_b[k])
    c = _count(
        tmp_path,
        contig,
        (ref_b, alt_b),
        [
            _reads(contig, len(ref_b), ref_b, "ref", 20, ends=ends),
            _reads(contig, len(ref_b), _with(ref_b, k, g), "germ", 20, ends=ends, first=20),
            _reads(contig, len(ref_b), _with(alt_b, k, g), "alt", 10, ends=ends, first=40),
        ],
    )
    assert (c.rd, c.ad, c.partial_alt) == (40, 10, 0)


@pytest.mark.parametrize("name, k", CASES, ids=[f"{n} kept {k}" for n, k in CASES])
def test_a_third_base_at_the_kept_base_does_not_take_the_ref_molecules_reads(tmp_path, name, k):
    """As above with 5 fragments of a third base at the kept base (no REF molecule's
    change): they are another allele, and the REF molecules carrying the change stay REF."""
    ref_b, alt_b, _ = MNPS[name]
    contig = _mnp(12, ref_b)
    g = _other(ref_b[k])
    t = _other(ref_b[k], g, alt_b[k])
    c = _count(
        tmp_path,
        contig,
        (ref_b, alt_b),
        [
            _reads(contig, len(ref_b), ref_b, "ref", 20),
            _reads(contig, len(ref_b), _with(ref_b, k, g), "germ", 20, first=20),
            _reads(contig, len(ref_b), _with(alt_b, k, g), "alt", 10, first=40),
            _reads(contig, len(ref_b), _with(ref_b, k, t), "third", 5, first=50),
        ],
    )
    assert (c.rd, c.ad) == (40, 10)
    assert c.partial_alt == 5


@pytest.mark.parametrize("name, k", CASES, ids=[f"{n} kept {k}" for n, k in CASES])
def test_the_change_read_below_min_bq_counts_on_either_allele(tmp_path, name, k):
    """Carriers and REF molecules reading the kept base below min BQ count as those
    reading the change: 4 more of each."""
    ref_b, alt_b, _ = MNPS[name]
    contig = _mnp(13, ref_b)
    g = _other(ref_b[k])
    c = _count(
        tmp_path,
        contig,
        (ref_b, alt_b),
        [
            _reads(contig, len(ref_b), ref_b, "ref", 20),
            _reads(contig, len(ref_b), _with(ref_b, k, g), "germ", 20, first=20),
            _reads(contig, len(ref_b), _with(alt_b, k, g), "alt", 10, first=40),
            _reads(contig, len(ref_b), _with(ref_b, k, g), "germlow", 4, low=[k], first=50),
            _reads(contig, len(ref_b), _with(alt_b, k, g), "altlow", 4, low=[k], first=54),
        ],
    )
    assert (c.rd, c.ad, c.partial_alt) == (44, 14, 0)


@pytest.mark.parametrize("name, k", CASES, ids=[f"{n} kept {k}" for n, k in CASES])
def test_carriers_written_with_an_indel_at_the_block_count_alike(tmp_path, name, k):
    """Half the carriers written with an insertion and a deletion at the block (their
    bases the same): ALT as the plain ones."""
    ref_b, alt_b, _ = MNPS[name]
    contig = _mnp(14, ref_b)
    g = _other(ref_b[k])
    c = _count(
        tmp_path,
        contig,
        (ref_b, alt_b),
        [
            _reads(contig, len(ref_b), ref_b, "ref", 20),
            _reads(contig, len(ref_b), _with(ref_b, k, g), "germ", 20, first=20),
            _reads(contig, len(ref_b), _with(alt_b, k, g), "alt", 5, first=40),
            _reads(contig, len(ref_b), _with(alt_b, k, g), "altgap", 5, indel=True, first=45),
        ],
    )
    assert (c.rd, c.ad) == (40, 10)


@pytest.mark.parametrize("name, k", CASES, ids=[f"{n} kept {k}" for n, k in CASES])
def test_a_change_only_the_ref_molecules_carry_is_background(tmp_path, name, k):
    """Control (trans): the carriers hold X's bases, half the REF molecules the change:
    REF 40, ALT 10."""
    ref_b, alt_b, _ = MNPS[name]
    contig = _mnp(15, ref_b)
    g = _other(ref_b[k])
    c = _count(
        tmp_path,
        contig,
        (ref_b, alt_b),
        [
            _reads(contig, len(ref_b), ref_b, "ref", 20),
            _reads(contig, len(ref_b), _with(ref_b, k, g), "germ", 20, first=20),
            _reads(contig, len(ref_b), alt_b, "alt", 10, first=40),
        ],
    )
    assert (c.rd, c.ad, c.partial_alt) == (40, 10, 0)


def test_a_two_base_mnp_keeps_no_base_and_its_neighbour_is_a_context(tmp_path):
    """Control: GC>AT changes both its bases; a germline change on the base after it, on
    the carriers and half the REF molecules, is beside it: REF 40, ALT 10."""
    ref_b, alt_b = "GC", "AT"
    contig = _mnp(16, ref_b)
    after = contig[P + 2]
    g = _other(after)

    def block(b):
        return b + g

    c = _count(
        tmp_path,
        contig,
        (ref_b, alt_b),
        [
            _reads(contig, 3, ref_b + after, "ref", 20),
            _reads(contig, 3, block(ref_b), "germ", 20, first=20),
            _reads(contig, 3, block(alt_b), "alt", 10, first=40),
        ],
    )
    assert (c.rd, c.ad, c.partial_alt) == (40, 10, 0)


def _delins_groups(contig, extra=()):
    """CAGTC>GAGC at P, a germline A>T at the base the ALT keeps (offset 1): 20 REF reads,
    20 REF molecules with the T, 10 carriers with it (a found allele, written as an
    insertion and a deletion), then `extra` groups."""
    groups = [
        _reads(contig, 5, "CAGTC", "ref", 20),
        _reads(contig, 5, "CTGTC", "germ", 20, first=20),
        _reads(contig, 5, "GTGC", "alt", 10, first=40),
    ]
    return groups + list(extra)


def _delins_contig():
    c = list(_contig(17))
    c[P - 1] = "T"
    c[P : P + 5] = "CAGTC"
    c[P + 5] = "A"
    return "".join(c)


def test_delins_carriers_reading_the_kept_base_below_min_bq_count_alt(tmp_path):
    """6 carriers reading the kept base below min BQ add 6 ALT and no partial evidence:
    the change is a separate event, so the carriers' found allele is no other allele."""
    contig = _delins_contig()
    base = _count(tmp_path / "a", contig, ("CAGTC", "GAGC"), _delins_groups(contig), "a")
    low = _reads(contig, 5, "GTGC", "altlow", 6, low=[1], first=50)
    more = _count(tmp_path / "b", contig, ("CAGTC", "GAGC"), _delins_groups(contig, [low]), "b")
    assert (more.ad - base.ad, more.partial_alt - base.partial_alt) == (6, 0)


def test_delins_ref_molecules_reading_the_kept_base_below_min_bq_count_ref(tmp_path):
    """6 REF molecules reading the kept base below min BQ add 6 REF."""
    contig = _delins_contig()
    base = _count(tmp_path / "a", contig, ("CAGTC", "GAGC"), _delins_groups(contig), "a")
    low = _reads(contig, 5, "CAGTC", "reflow", 6, low=[1], first=50)
    more = _count(tmp_path / "b", contig, ("CAGTC", "GAGC"), _delins_groups(contig, [low]), "b")
    assert (more.rd - base.rd, more.partial_alt - base.partial_alt) == (6, 0)
