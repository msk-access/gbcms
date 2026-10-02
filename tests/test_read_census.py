"""The engine's pure-indel counts agree with the read census (#171).

`census.py` judges each read by its own bases across the indel's whole tract. This
module checks it on hand-built reads, then generates reads around ten pure indels
(homopolymer, STR, unique and duplication contexts, insertions and deletions). The
reads start before the anchor and end anywhere, from just past it to past the
tract. For each class of read, the engine's REF and ALT counts must equal the
census's:
- REF reads;
- carriers with the indel written at any equivalent placement;
- wrong-length indels of the same kind where the rules are settled;
- the indel's bases written well outside the tract.

Read shapes whose rules are still open decisions are strict xfails. They turn
green only when the decision lands:
- C27 #201: the ALT spelled across several ops;
- C28 #202: a read that deletes the anchor falls back to Phase 3.

Clean bases only (Q30): a pure indel's REF call stands on its CIGAR and extent, so
sequencing errors inside the window are not what this compares. Wrong-length reads
run past where their allele differs from the ALT: the census knows two
haplotypes, so a read that ends before then fits the ALT, while the engine reads
its CIGAR's length and calls it another allele.
"""

import random

import pytest
from census import Verdict, assert_matches, census, judge, tract, window_for
from helpers import make_read, write_contig

from gbcms import _rs as gbcms_rs

A, L = 400, 900  # anchor position, contig length
ROWS = {
    "hp-A": ("GAAAAAT", "GA", "G"),
    "hp-AA": ("GAAAAAAT", "GAA", "G"),
    "hp+A": ("GAAAAAT", "G", "GA"),
    "hp+AA": ("GAAAAAT", "G", "GAA"),
    "u-8": ("GATCGGATAC", "GATCGGATA", "G"),
    "u+10": ("GT", "G", "GACGTTGCATC"),
    "ca-CA": ("G" + "CA" * 5 + "T", "GCA", "G"),
    "ca+CA": ("G" + "CA" * 5 + "T", "G", "GCA"),
    "dup-8": ("G" + "ACGTTGCA" * 2 + "T", "GACGTTGCA", "G"),
    "dup+8": ("G" + "ACGTTGCA" * 2 + "T", "G", "GACGTTGCA"),
}
ARGS = {
    "min_mapq": 20,
    "min_baseq": 20,
    "filter_duplicates": True,
    "filter_secondary": True,
    "filter_supplementary": True,
    "filter_qc_failed": False,
    "filter_improper_pair": False,
    "filter_indel": False,
    "threads": 1,
}


def _contig(motif, seed=189):
    rng = random.Random(seed)
    seq = [rng.choice("CGT") for _ in range(L)]
    seq[A : A + len(motif)] = list(motif)
    return "".join(seq)


def _read(contig, start, length, events, aligned=None):
    """A read of `length` bases from `start` with `events` placed before their
    reference positions: (pos, "I", bases) or (pos, "D", n), soft-clipped after
    `aligned` bases when given. None when the aligned part would end inside an
    insertion or in a deletion."""
    seq, cigar, p = "", [], start
    for pos, op, x in events:
        if pos > p:
            seq += contig[p:pos]
            cigar.append([0, pos - p])
            p = pos
        if op == "I":
            seq += x
            cigar.append([1, len(x)])
        else:
            cigar.append([2, x])
            p += x
    seq += contig[p : p + length]
    cigar.append([0, length])
    keep = length if aligned is None else min(aligned, length)
    out, q = [], 0
    for op, n in cigar:
        if q >= keep:
            break
        if op == 2:
            out.append([2, n])
            continue
        take = min(n, keep - q)
        if op == 1 and take < n:
            return None
        out.append([op, take])
        q += take
    if not out or out[-1][0] == 2:
        return None
    merged = []  # adjacent aligned blocks join; split indel ops stay split
    for op, n in out:
        if merged and merged[-1][0] == op == 0:
            merged[-1][1] += n
        else:
            merged.append([op, n])
    if keep < length:
        merged.append([4, length - keep])
    return seq[:length], tuple(map(tuple, merged))


def _engine_and_census(tmp_path, name, contig, ref, alt, shapes):
    """Reads from `shapes` [(start, length, events)] in one BAM; the engine's
    counts and the census for the row."""
    reads = []
    for i, (s, n, events, *aligned) in enumerate(shapes):
        built = _read(contig, s, n, events, *aligned)
        if built is not None:
            reads.append(make_read(f"{name}{i}", built[0], s, built[1]))
    fa, bam = write_contig(tmp_path, contig, reads, name)
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", A, ref, alt, "X")], fa, 5, False, 1, True
    )
    assert pv.gbcms_status == "PASS" and pv.variant.pos == A, pv.gbcms_status_reason
    (counts,) = gbcms_rs.count_bam_binned(bam, [pv.variant], [None], **ARGS)
    return counts, census(bam, contig, pv.variant), pv.variant, len(reads)


# ── The census on hand-built reads ────────────────────────────────────────────
def _verdict(tmp_path, row, events, start, length):
    motif, ref, alt = ROWS[row]
    contig = _contig(motif)
    built = _read(contig, start, length, events)
    fa, _ = write_contig(tmp_path, contig, [], f"v{start}{length}")
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", A, ref, alt, "X")], fa, 5, False, 1, True
    )
    return judge(make_read("x", built[0], start, built[1]), window_for(pv.variant, contig))


@pytest.mark.parametrize(
    "row, events, length, expected",
    [
        ("hp-A", [], 100, Verdict.REF),  # spans the run and the flank
        ("hp-A", [(A + 1, "D", 1)], 100, Verdict.ALT),
        ("hp-A", [(A + 4, "D", 1)], 100, Verdict.ALT),  # written elsewhere in the run
        ("hp-A", [(A + 1, "D", 1)], 54, Verdict.FITS_BOTH),  # ends inside the run
        ("hp-A", [(A + 1, "D", 2)], 100, Verdict.CONTRADICTS_BOTH),  # -AA on a -A row
        ("hp+A", [(A + 1, "I", "A")], 57, Verdict.ALT),  # ends on the deciding base (the 6th A)
        ("hp+A", [], 56, Verdict.FITS_BOTH),  # REF ending inside the run
        ("hp+A", [], 57, Verdict.REF),  # REF ending on the flank: the whole run read
        ("u+10", [(A + 1, "I", "ACGTTGCATC")], 100, Verdict.ALT),
    ],
)
def test_the_census_judges_hand_built_reads(tmp_path, row, events, length, expected):
    assert _verdict(tmp_path, row, events, A - 50, length) == expected


def test_a_masked_deciding_base_is_read_past(tmp_path):
    motif, ref, alt = ROWS["u+10"]
    contig = _contig(motif)
    seq, cigar = _read(contig, A - 50, 61, [(A + 1, "I", "ACGTTGCATC")])  # ends with the insert
    fa, _ = write_contig(tmp_path, contig, [], "m")
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", A, ref, alt, "X")], fa, 5, False, 1, True
    )
    w = window_for(pv.variant, contig)
    quals = [30] * len(seq)
    quals[51] = 5  # the first inserted base
    assert judge(make_read("m", seq, A - 50, cigar, quals=quals), w) == Verdict.ALT
    quals[51:61] = [5] * 10  # every inserted base: nothing discriminates
    assert judge(make_read("m", seq, A - 50, cigar, quals=quals), w) == Verdict.FITS_BOTH


def test_reading_stops_at_a_splice(tmp_path):
    """A read spliced out across the whole window says nothing about it."""
    motif, ref, alt = ROWS["hp-A"]
    contig = _contig(motif)
    fa, _ = write_contig(tmp_path, contig, [], "n")
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", A, ref, alt, "X")], fa, 5, False, 1, True
    )
    w = window_for(pv.variant, contig)
    s = A - 50
    seq = contig[s : A - 2] + contig[A + 30 : A + 30 + 48]
    spliced = make_read("n", seq, s, ((0, A - 2 - s), (3, 32), (0, 48)))
    assert judge(spliced, w) == Verdict.SPLICED


# ── The engine against the census ─────────────────────────────────────────────
def _placements(contig, ref, alt):
    """Equivalent junctions (insertions) or starts (deletions) of the indel, from
    the census's own slide."""
    lo, hi = tract(A, ref, alt, contig)
    n = abs(len(alt) - len(ref))
    return list(range(lo, hi + 1)) if len(alt) > len(ref) else list(range(lo, hi - n + 1))


def _settled_shapes(contig, ref, alt, rng):
    """(class, [(start, length, events[, aligned])]) for the settled classes. Reads
    start on the tract's left flank (the anchor) or anywhere before it."""
    lo, hi = tract(A, ref, alt, contig)
    ins = len(alt) > len(ref)
    n = abs(len(alt) - len(ref))
    shapes = {"ref": [], "alt": [], "wrong": [], "elsewhere": [], "clipped": []}
    for k in range(40):
        # Every read overlaps the anchor; the first three start on it (the flank).
        s = lo - 1 if k < 3 else lo - 1 - rng.randint(1, 60)
        length = rng.randint(A + 1 - s, 100)  # and ends anywhere past it
        shapes["ref"].append((s, length, []))
        j = rng.choice(_placements(contig, ref, alt))
        if ins:
            # The inserted bases are the haplotype's at this junction (a rotation).
            k = j - (A + 1)
            hap = contig[: A + 1] + alt[1:] + contig[A + 1 :]
            shapes["alt"].append((s, length, [(j, "I", hap[j : j + n] if k >= 0 else alt[1:])]))
            wrong = (j, "I", hap[j : j + n] + hap[j : j + 1])  # one base longer
            shapes["wrong"].append((s, 100, [wrong]))  # read on past where it differs from the ALT
        else:
            shapes["alt"].append((s, length, [(j, "D", n)]))
            if n >= 5:
                shapes["wrong"].append((s, 100, [(j, "D", n - 1)]))
            else:
                shapes["wrong"].append((s, 100, [(A + 1, "D", n + 1)]))  # at the junction
        far = hi + 6 + rng.randint(0, 10)  # well past the tract and its flank
        far_event = [(far, "I", alt[1:])] if ins else [(far, "D", n)]
        shapes["elsewhere"].append((s, rng.randint(A + 1 - s, far - s + 20), far_event))
        # Soft-clipped from a point past the anchor: a pure indel's rules read only
        # aligned bases, so the clipped bases decide nothing.
        events = [] if rng.random() < 0.5 else shapes["alt"][-1][2]
        shapes["clipped"].append((s, 100, events, rng.randint(A + 1 - s, 100)))
    return shapes


@pytest.mark.parametrize("row", sorted(ROWS))
def test_prep_finds_the_census_tract(tmp_path, row):
    """Prep's shift region is the tract the census slides independently."""
    motif, ref, alt = ROWS[row]
    contig = _contig(motif)
    fa, _ = write_contig(tmp_path, contig, [], "prep")
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", A, ref, alt, "X")], fa, 5, False, 1, True
    )
    assert tuple(pv.variant.shift_region) == tract(A, ref, alt, contig)


@pytest.mark.parametrize("row", sorted(ROWS))
def test_engine_counts_equal_the_census(tmp_path, row):
    motif, ref, alt = ROWS[row]
    contig = _contig(motif)
    rng = random.Random(row)
    for cls, shapes in _settled_shapes(contig, ref, alt, rng).items():
        counts, result, _, n_reads = _engine_and_census(
            tmp_path, f"{row}_{cls}", contig, ref, alt, shapes
        )
        assert n_reads > 0
        assert_matches(counts, result, f"{row} {cls}")


# ── Open decisions: strict xfails until they land ─────────────────────────────
def _open_case(tmp_path, row, events, name):
    motif, ref, alt = ROWS[row]
    contig = _contig(motif)
    shapes = [(A - 50 + i, 100, events) for i in range(5)]
    counts, result, _, _ = _engine_and_census(tmp_path, name, contig, ref, alt, shapes)
    assert_matches(counts, result, name)


@pytest.mark.parametrize(
    "row, events",
    [
        ("hp-AA", [(A + 3, "D", 1)]),  # -A written in the run, off the junction, on a -AA row
        ("hp+A", [(A + 1, "D", 1)]),  # -A carriers on a +A row
    ],
    ids=["short-wrong-length-off-junction", "opposite-kind"],
)
def test_c26_short_indels_in_the_window(tmp_path, row, events):
    _open_case(tmp_path, row, events, "c26")


@pytest.mark.xfail(
    strict=True, reason="C27 #201: the ALT spelled across several ops counts partial"
)
def test_c27_the_alt_split_across_ops(tmp_path):
    _open_case(tmp_path, "hp-AA", [(A + 1, "D", 1), (A + 2, "D", 1)], "c27")


@pytest.mark.xfail(strict=True, reason="C28 #202: a read deleting the anchor falls back to Phase 3")
@pytest.mark.parametrize(
    "events",
    [
        [(A, "D", 3), (A + 5, "I", "T")],  # holds neither allele; the closer is the ALT
        [(A, "D", 3)],  # holds neither allele; the closer is the REF
    ],
    ids=["credited-alt", "credited-ref"],
)
def test_c28_a_read_deleting_the_anchor(tmp_path, events):
    """Phase 3 credits whichever haplotype is closer to a read that deletes the
    anchor, so reads holding neither allele count ALT or REF. (A read deleting the
    anchor whose bases spell the ALT one base along is ALT, and is: judge bases.)"""
    _open_case(tmp_path, "u-8", events, "c28")


def test_a_pure_indels_clipped_bases_decide_nothing(tmp_path):
    """The rule reads aligned bases only: a REF-haplotype read aligned up to the run
    and soft-clipped from it fits both alleles, as the engine counts it."""
    motif, ref, alt = ROWS["hp+A"]
    contig = _contig(motif)
    s = A - 50
    seq, cigar = _read(contig, s, 100, [], aligned=A + 1 - s)
    fa, _ = write_contig(tmp_path, contig, [], "clip")
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", A, ref, alt, "X")], fa, 5, False, 1, True
    )
    assert judge(make_read("c", seq, s, cigar), window_for(pv.variant, contig)) == Verdict.FITS_BOTH


def test_a_substitution_needs_no_margin(tmp_path):
    """An SNV's REF read ending on the SNV's base reads it: REF (the margin is a
    pure indel's rule)."""
    contig = _contig("GCT")
    fa, _ = write_contig(tmp_path, contig, [], "snv")
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", A + 1, "C", "A", "SNP")], fa, 5, False, 1, True
    )
    s = A - 50
    read = make_read("s", contig[s : A + 2], s, ((0, A + 2 - s),))  # ends on the C
    assert judge(read, window_for(pv.variant, contig)) == Verdict.REF


@pytest.mark.parametrize("row", sorted(ROWS))
def test_reads_starting_on_the_flank_equal_the_census(tmp_path, row):
    """Reads whose first base is the tract's left flank (the anchor): REF reads and
    carriers ending anywhere, and carriers whose flank base is masked or wrong. The
    read-by-bases rule reads from a flank only when the read reads it (unmasked, the
    reference's), so the engine and the census agree on each kind."""
    motif, ref, alt = ROWS[row]
    contig = _contig(motif)
    lo, hi = tract(A, ref, alt, contig)
    s = lo - 1
    rng = random.Random(f"flank-{row}")
    placements = _placements(contig, ref, alt)
    n = abs(len(alt) - len(ref))
    ins = len(alt) > len(ref)
    kinds = {}
    for kind in ("ref", "alt", "masked", "wrong"):
        reads = []
        for i in range(12):
            length = rng.randint(2, 100)
            if kind == "ref":
                events = []
            else:
                j = rng.choice(placements)
                if ins:
                    hap = contig[: A + 1] + alt[1:] + contig[A + 1 :]
                    events = [(j, "I", hap[j : j + n])]
                else:
                    events = [(j, "D", n)]
            built = _read(contig, s, length, events)
            if built is None:
                continue
            seq, cigar = built
            quals = [30] * len(seq)
            if kind == "masked":
                quals[0] = 2
            elif kind == "wrong":
                seq = ("T" if seq[0] != "T" else "C") + seq[1:]
            reads.append(make_read(f"{kind}{i}", seq, s, cigar, quals=quals))
        fa, bam = write_contig(tmp_path, contig, reads, f"{row}_{kind}")
        (pv,) = gbcms_rs.prepare_variants(
            [gbcms_rs.Variant("1", A, ref, alt, "X")], fa, 5, False, 1, True
        )
        (counts,) = gbcms_rs.count_bam_binned(bam, [pv.variant], [None], **ARGS)
        result = census(bam, contig, pv.variant)
        assert_matches(counts, result, f"{row} flank {kind}")
        kinds[kind] = result
    if not ins:
        # A deletion's flank is all that tells such a read from REF one base along
        # (an insertion's longer run decides it before the flank).
        assert kinds["masked"].ad == 0 and kinds["wrong"].ad == 0
