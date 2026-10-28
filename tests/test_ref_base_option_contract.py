"""The sequence judge's REF-base option (`GBCMS_JUDGE_REFBASE=partial`, measurement only;
the default `margin` is the judge as decided), #304.

Mechanism. At a base an MNP or a delins changes, a read's confident base that is not the
called allele's is one of two things. The other allele's base there (under the read's own
placement: at a delins, the other allele's base at the same telling index counted from
either end, where the read's own gap may put it) is that allele's evidence, not an error:
a read showing one allele's bases elsewhere and the other's there carries neither, so it
is another allele, partial evidence (RJ-25: an allele lacking part of the event; RJ-28 at
an MNP). A third base, neither allele's, is a one-read error that costs both alleles alike,
and the margin forgives it (RJ-26), at an MNP as at a delins.

The default differs at two shapes, both decided: at a delins the margin alone decides
(2026-10-24), so a carrier with REF's base at a telling base stays ALT (and a REF read with
the ALT's base stays REF); at an MNP an ALT call needs every base X changes read as X's
(proposal §1, phase C), so a carrier with a third base there is depth only. Everything
else is the same under both: a base below `--min-baseq` weighs nothing, a base X keeps
tells nothing (a germline change there is a separate event, RJ-25), and at an MNP a read
mixing REF's and X's bases is another allele already (RJ-28).

Shapes: delins under 20 bases (the 3-to-1 shape of a patient normal's row), 20-49 and 50
or more; MNPs of 2, 3 and 4 bases; the other allele's base against a third base, on
carriers and on REF reads; both strands, a base below quality, a read ending a few bases
past the event, and the odd base in a soft clip (a paired read, the clip inside its
fragment).
"""

import random

import pytest
from helpers import count_bam_checked, make_read, write_contig

from gbcms import _rs

ARGS = (20, 20, True, True, True, False, False, False, 1)
READ = 150
P = 300
FLANK = "A"  # the bases either side of the event, unlike any allele's first or last base


def _contig(ref_al, seed=304):
    rng = random.Random(seed)
    bases = [rng.choice("ACGT") for _ in range(700)]
    bases[P - 1] = FLANK
    bases[P : P + len(ref_al)] = list(ref_al)
    bases[P + len(ref_al)] = FLANK
    return "".join(bases)


def _alleles(ref_n, alt_n, seed):
    """A REF of `ref_n` and an ALT of `alt_n` bases from C, G and T, sharing no first or
    last base (so neither end is trimmed)."""
    rng = random.Random(seed)
    while True:
        r = "".join(rng.choice("CGT") for _ in range(ref_n))
        a = "".join(rng.choice("CGT") for _ in range(alt_n))
        if r[0] != a[0] and r[-1] != a[-1]:
            return r, a


DELINS = {
    "under 20 (3>1)": ("CGC", "T"),
    "20-49 (25>10)": _alleles(25, 10, 2049),
    "50+ (60>20)": _alleles(60, 20, 5099),
}
MNPS = {"2": ("AC", "GT"), "3": ("ACG", "GTA"), "4": ("ACGT", "GTAC")}


def _paired(read, tlen=400):
    """One mate of a proper pair spanning `tlen` bases, its mate downstream."""
    read.flag = 0x1 | 0x2 | 0x40 | 0x20
    read.next_reference_id = 0
    read.next_reference_start = read.reference_start + tlen - READ
    read.template_length = tlen
    return read


class Locus:
    """REF > ALT at P: 20 REF reads and 10 carriers (half of each reverse), and the reads
    a case adds, each an allele read with one base changed."""

    def __init__(self, ref_al, alt_al):
        self.ref_al, self.alt_al = ref_al, alt_al
        self.contig = _contig(ref_al)
        self.hap = self.contig[:P] + alt_al + self.contig[P + len(ref_al) :]
        self.reads = []
        for i in range(20):
            self.reads.append(self.read(f"ref{i}", False, i))
        for i in range(10):
            self.reads.append(self.read(f"alt{i}", True, i))

    def read(self, name, carrier, i, k=None, base=None, rev=False, low=False, end=None, clip=False):
        """An allele read starting 30 to 39 bases before P; with `k`, its base at
        offset `k` into the allele set to `base` (quality 10 where `low`); `end`, bases
        past the allele it reads; `clip`, aligned up to P and soft-clipped from there."""
        s = P - 30 - (i % 10)
        seq_from = self.hap if carrier else self.contig
        allele = self.alt_al if carrier else self.ref_al
        n = READ if end is None else (P - s) + len(allele) + end
        seq = list(seq_from[s : s + n])
        quals = [30] * n
        if k is not None:
            seq[P - s + k] = base
            if low:
                quals[P - s + k] = 10
        left = P - s
        if clip:
            cig = ((0, left), (4, n - left))
        elif carrier and len(self.alt_al) != len(self.ref_al):
            cig = (
                (0, left),
                (1, len(self.alt_al)),
                (2, len(self.ref_al)),
                (0, n - left - len(self.alt_al)),
            )
        else:
            cig = ((0, n),)
        r = make_read(name, "".join(seq), s, cig, flag=16 if rev else 0, quals=quals)
        return _paired(r) if clip else r

    def counts(self, tmp_path, monkeypatch, rule):
        if rule is None:
            monkeypatch.delenv("GBCMS_JUDGE_REFBASE", raising=False)
        else:
            monkeypatch.setenv("GBCMS_JUDGE_REFBASE", rule)
        d = tmp_path / (rule or "default")
        d.mkdir(exist_ok=True)
        fa, bam = write_contig(d, self.contig, self.reads)
        (pv,) = _rs.prepare_variants(
            [_rs.Variant("1", P, self.ref_al, self.alt_al, "X")], fa, 5, False, 1, True
        )
        (c,) = count_bam_checked(
            bam, [pv.variant], [None], *ARGS, read_judge="sequence", reference_fasta=fa
        )
        assert c.dp >= c.rd + c.ad and c.dpf >= c.rdf + c.adf
        assert c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
        return (c.rd, c.ad, c.partial_alt, c.dp)


def _third(*taken):
    return next(b for b in "ACGT" if b not in taken)


def _other_at(own, other, k):
    """The other allele's base at telling index `k` of `own`, counted from the left."""
    return other[k] if k < len(other) else other[len(other) - (len(own) - k)]


def _both(locus, tmp_path, monkeypatch):
    return locus.counts(tmp_path, monkeypatch, None), locus.counts(tmp_path, monkeypatch, "partial")


# ── delins ───────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("rev", [False, True], ids=["forward", "reverse"])
@pytest.mark.parametrize("size", list(DELINS))
def test_a_delins_carrier_with_refs_base_at_a_telling_base_is_another_allele(
    tmp_path, monkeypatch, size, rev
):
    """The read carries the length change and REF's base where the ALT has another (the
    normal's 3>1 row: the 2-base deletion with REF's C where the ALT's T is): it lacks
    part of the event. Default: ALT (the margin forgives it). `partial`: another allele."""
    ref_al, alt_al = DELINS[size]
    loc = Locus(ref_al, alt_al)
    k = 0
    loc.reads.append(loc.read("odd", True, 3, k=k, base=ref_al[k], rev=rev))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == (20, 11, 0), default
    assert partial[:3] == (20, 10, 1), partial
    assert default[3] == partial[3] == 31


@pytest.mark.parametrize("size", list(DELINS))
def test_a_delins_carrier_with_a_third_base_stays_alt(tmp_path, monkeypatch, size):
    """A base neither allele holds at that index from either end is an error the margin
    forgives, under both."""
    ref_al, alt_al = DELINS[size]
    loc = Locus(ref_al, alt_al)
    k = len(alt_al) // 2
    third = _third(alt_al[k], _other_at(alt_al, ref_al, k), ref_al[len(ref_al) - (len(alt_al) - k)])
    loc.reads.append(loc.read("odd", True, 3, k=k, base=third))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == partial[:3] == (20, 11, 0), (default, partial)


@pytest.mark.parametrize("size", list(DELINS))
def test_a_delins_ref_read_with_the_alts_base_at_a_telling_base_is_another_allele(
    tmp_path, monkeypatch, size
):
    """The mirror: REF's length with the ALT's base at a REF telling base. Default: REF.
    `partial`: another allele."""
    ref_al, alt_al = DELINS[size]
    loc = Locus(ref_al, alt_al)
    loc.reads.append(loc.read("odd", False, 3, k=0, base=alt_al[0]))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == (21, 10, 0), default
    assert partial[:3] == (20, 10, 1), partial


@pytest.mark.parametrize("size", list(DELINS))
def test_a_delins_ref_read_with_a_third_base_stays_ref(tmp_path, monkeypatch, size):
    ref_al, alt_al = DELINS[size]
    loc = Locus(ref_al, alt_al)
    loc.reads.append(loc.read("odd", False, 3, k=0, base=_third(ref_al[0], alt_al[0], alt_al[-1])))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == partial[:3] == (21, 10, 0), (default, partial)


@pytest.mark.parametrize("size", list(DELINS))
def test_a_delins_odd_base_below_quality_weighs_nothing(tmp_path, monkeypatch, size):
    ref_al, alt_al = DELINS[size]
    loc = Locus(ref_al, alt_al)
    loc.reads.append(loc.read("odd", True, 3, k=0, base=ref_al[0], low=True))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == partial[:3] == (20, 11, 0), (default, partial)


def test_a_delins_carrier_ending_past_the_event_is_read_alike(tmp_path, monkeypatch):
    """The odd read ends 8 bases past the event (the 25>10 shape)."""
    ref_al, alt_al = DELINS["20-49 (25>10)"]
    loc = Locus(ref_al, alt_al)
    loc.reads.append(loc.read("odd", True, 3, k=0, base=ref_al[0], end=8))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == (20, 11, 0), default
    assert partial[:3] == (20, 10, 1), partial


def test_a_delins_odd_base_in_a_soft_clip_is_read_alike(tmp_path, monkeypatch):
    """The odd read is aligned up to the event and soft-clipped from it, inside its
    fragment, so its clipped bases are read (RJ-13); the 25>10 shape."""
    ref_al, alt_al = DELINS["20-49 (25>10)"]
    loc = Locus(ref_al, alt_al)
    loc.reads.append(loc.read("odd", True, 3, k=0, base=ref_al[0], clip=True))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == (20, 11, 0), default
    assert partial[:3] == (20, 10, 1), partial


# ── MNPs ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("rev", [False, True], ids=["forward", "reverse"])
@pytest.mark.parametrize("n", list(MNPS))
def test_an_mnp_carrier_with_a_third_base_at_a_changed_base_counts_alt(
    tmp_path, monkeypatch, n, rev
):
    """Default: depth only (an ALT call needs every changed base read as X's).
    `partial`: ALT, the error forgiven as at a delins."""
    ref_al, alt_al = MNPS[n]
    loc = Locus(ref_al, alt_al)
    loc.reads.append(loc.read("odd", True, 3, k=0, base=_third(ref_al[0], alt_al[0]), rev=rev))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == (20, 10, 0), default
    assert partial[:3] == (20, 11, 0), partial
    assert default[3] == partial[3] == 31


@pytest.mark.parametrize("n", list(MNPS))
def test_an_mnp_carrier_with_refs_base_at_a_changed_base_is_another_allele(
    tmp_path, monkeypatch, n
):
    """Guard: one component short of the MNV is another allele under both (RJ-28)."""
    ref_al, alt_al = MNPS[n]
    loc = Locus(ref_al, alt_al)
    loc.reads.append(loc.read("odd", True, 3, k=0, base=ref_al[0]))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == partial[:3] == (20, 10, 1), (default, partial)


@pytest.mark.parametrize("n", list(MNPS))
def test_an_mnp_ref_read_with_a_third_base_stays_ref(tmp_path, monkeypatch, n):
    ref_al, alt_al = MNPS[n]
    loc = Locus(ref_al, alt_al)
    loc.reads.append(loc.read("odd", False, 3, k=0, base=_third(ref_al[0], alt_al[0])))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == partial[:3] == (21, 10, 0), (default, partial)


@pytest.mark.parametrize("n", list(MNPS))
def test_an_mnp_ref_read_with_one_of_the_alts_bases_is_another_allele(tmp_path, monkeypatch, n):
    """Guard: one component of the MNV is another allele under both (RJ-28)."""
    ref_al, alt_al = MNPS[n]
    loc = Locus(ref_al, alt_al)
    loc.reads.append(loc.read("odd", False, 3, k=0, base=alt_al[0]))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == partial[:3] == (20, 10, 1), (default, partial)


@pytest.mark.parametrize("n", list(MNPS))
def test_an_mnp_third_base_below_quality_weighs_nothing(tmp_path, monkeypatch, n):
    ref_al, alt_al = MNPS[n]
    loc = Locus(ref_al, alt_al)
    loc.reads.append(loc.read("odd", True, 3, k=0, base=_third(ref_al[0], alt_al[0]), low=True))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == partial[:3] == (20, 11, 0), (default, partial)


def test_an_mnp_carrier_ending_past_the_event_is_read_alike(tmp_path, monkeypatch):
    ref_al, alt_al = MNPS["3"]
    loc = Locus(ref_al, alt_al)
    loc.reads.append(loc.read("odd", True, 3, k=1, base=_third(ref_al[1], alt_al[1]), end=8))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == (20, 10, 0), default
    assert partial[:3] == (20, 11, 0), partial


def test_an_mnp_third_base_in_a_soft_clip_is_read_alike(tmp_path, monkeypatch):
    ref_al, alt_al = MNPS["3"]
    loc = Locus(ref_al, alt_al)
    loc.reads.append(loc.read("odd", True, 3, k=1, base=_third(ref_al[1], alt_al[1]), clip=True))
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == (20, 10, 0), default
    assert partial[:3] == (20, 11, 0), partial


# ── a base the event keeps ───────────────────────────────────────────────────


@pytest.mark.parametrize("ref_al, alt_al", [("CAG", "TAC"), ("CAGT", "GACA")], ids=["mnp", "mnp-4"])
def test_a_germline_change_at_a_kept_base_stays_a_separate_event(
    tmp_path, monkeypatch, ref_al, alt_al
):
    """A germline change at the base the MNP keeps, on half the REF reads and half the
    carriers: a separate event on both haplotypes (RJ-25). Its carriers stay ALT and its
    REF reads REF, under both."""
    loc = Locus(ref_al, alt_al)
    kept = next(i for i in range(len(ref_al)) if ref_al[i] == alt_al[i])
    germ = _third(ref_al[kept], "A")
    loc.reads = [
        loc.read(r.query_name, r.query_name.startswith("alt"), i, k=kept, base=germ) if i % 2 else r
        for i, r in enumerate(loc.reads)
    ]
    default, partial = _both(loc, tmp_path, monkeypatch)
    assert default[:3] == partial[:3] == (20, 10, 0), (default, partial)
