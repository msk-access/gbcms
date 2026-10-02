"""What a read contributes (group 1: #176, #182, #183).

- RJ-10 (C17): a read ends at its fragment end. Bases past it (read-through into
  adapter when the insert is shorter than the read) are neither bases nor reach.
- RJ-11 (C19): a record with absent base qualities (QUAL '*', stored as 0xFF) is
  dropped by the read filter and warned once per BAM.
- RJ-12 (O7): an unmapped record is not an alignment: dropped, so it reaches no
  count (mq0_count included), also at --min-mapq 0. A mapped MAPQ-0 alignment
  still counts there (pseudogene loci such as PMS2 are run at --min-mapq 0).

Committed red (xfail-strict) before the fixes; see docs/reference/read-judgment.md.
"""

import logging
import random

import pytest
from census import census
from helpers import make_read, write_contig

from gbcms import _rs

A = 400
ENGINE_LOGGER = "_rs.counting.engine"


def _contig(motif, seed=189):
    rng = random.Random(seed)
    seq = [rng.choice("CGT") for _ in range(900)]
    seq[A : A + len(motif)] = list(motif)
    return "".join(seq)


def _count(bam, variant, min_mapq=20):
    (c,) = _rs.count_bam_binned(
        bam, [variant], [None], min_mapq, 20, True, True, True, False, False, False, 1
    )
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    return c


def _pair(name, contig, s, frag_end, r1_seq=None):
    """R1 forward from s, 100 bases, its fragment [s, frag_end): bases past
    frag_end are adapter. R2 reverse ends at frag_end."""
    a = make_read(name, r1_seq or contig[s : s + 100], s, ((0, 100),), flag=99)
    b = make_read(name, contig[frag_end - 100 : frag_end], frag_end - 100, ((0, 100),), flag=147)
    t = frag_end - s
    a.template_length, b.template_length = t, -t
    a.next_reference_id = b.next_reference_id = 0
    a.next_reference_start, b.next_reference_start = frag_end - 100, s
    return [a, b]


# ── RJ-10: a read ends at its fragment end ─────────────────────────────────


def test_an_adapter_base_on_an_snv_is_not_the_reads(tmp_path):
    """R1 reads 50 bases past its fragment's end, onto the SNV: those bases are
    adapter, so R1 neither shows an allele there nor reaches the SNV (depth)."""
    contig = _contig("GCT")
    reads = []
    for i in range(4):
        s = A - 60 + i
        seq = list(contig[s : s + 100])
        seq[A + 1 - s] = "A"  # the ALT, on an adapter base
        reads += _pair(f"f{i}", contig, s, s + 50, "".join(seq))
    _, bam = write_contig(tmp_path, contig, reads, "snv")
    c = _count(bam, _rs.Variant("1", A + 1, contig[A + 1], "A", "SNP"))
    assert (c.dp, c.rd, c.ad) == (0, 0, 0)


def test_a_molecule_ending_inside_the_run_is_not_ref_on_its_adapter(tmp_path):
    """GA>G in G A*5 T. The molecule ends inside the run (A+3); R1's adapter bases
    align on past it. Its molecule cannot tell the alleles apart: depth only, as
    the census says (the census reads only the molecule's bases)."""
    contig = _contig("GAAAAAT")
    reads = []
    for i in range(4):
        reads += _pair(f"m{i}", contig, A - 60 + i, A + 3)
    fa, bam = write_contig(tmp_path, contig, reads, "run")
    (pv,) = _rs.prepare_variants([_rs.Variant("1", A, "GA", "G", "X")], fa, 5, False, 1, True)
    c = _count(bam, pv.variant)
    assert c.rd == 0, (c.rd, c.ad)
    result = census(bam, contig, pv.variant)
    assert result.rd == 0, dict(result.counts)


def test_a_molecule_spanning_the_run_is_ref(tmp_path):
    """Guard: the same pairs with the molecule reaching past the run are REF."""
    contig = _contig("GAAAAAT")
    reads = []
    for i in range(4):
        reads += _pair(f"m{i}", contig, A - 60 + i, A + 20)
    fa, bam = write_contig(tmp_path, contig, reads, "span")
    (pv,) = _rs.prepare_variants([_rs.Variant("1", A, "GA", "G", "X")], fa, 5, False, 1, True)
    assert _count(bam, pv.variant).rd == 8


def _counts(c):
    return (c.dp, c.rd, c.ad, c.partial_alt, c.rdf, c.adf)


@pytest.mark.parametrize(
    "ref,alt", [("TA", "GG"), ("TAC", "GGT"), ("TAC", "GG"), ("T", "GA")], ids=str
)
def test_a_read_past_its_fragment_end_counts_as_the_read_trimmed_there(tmp_path, ref, alt):
    """A forward molecule ends on the event's first base, which carries the ALT's
    first base; R1 reads on into adapter. Its call must be the call of the same
    read trimmed at its fragment end: the adapter neither matches nor reaches."""
    contig = _contig("TACGT")
    counts = []
    for trimmed in (False, True):
        reads = []
        for i in range(4):
            s, end = A - 60 + i, A + 1
            seq = list(contig[s : s + 100])
            seq[A - s] = alt[0]
            r1 = "".join(seq)
            pair = _pair(f"t{i}", contig, s, end, r1)
            if trimmed:
                n = end - s
                pair[0].query_sequence = r1[:n]
                pair[0].cigartuples = ((0, n),)
                pair[0].query_qualities = [30] * n
            reads += pair
        fa, bam = write_contig(tmp_path, contig, reads, f"cx{int(trimmed)}")
        (pv,) = _rs.prepare_variants([_rs.Variant("1", A, ref, alt, "X")], fa, 5, False, 1, True)
        counts.append(_counts(_count(bam, pv.variant)))
    assert counts[0] == counts[1], counts
    assert counts[0][2] == 0, counts


def test_an_outward_pair_keeps_its_bases(tmp_path):
    """R1 forward at s, R2 reverse ending before s: an outward-facing (RF) pair,
    as at a tandem-duplication junction. BWA writes TLEN 5' to 5' with R1's
    negative; no fragment is defined, so nothing past it is adapter."""
    contig = _contig("GCT")
    reads = []
    for i in range(6):
        s = A - 20 + i
        seq = list(contig[s : s + 100])
        seq[A + 1 - s] = "A"
        a = make_read(f"o{i}", "".join(seq), s, ((0, 100),), flag=0x1 | 0x20 | 0x40)
        b = make_read(
            f"o{i}", contig[s - 111 : s - 11], s - 111, ((0, 100),), flag=0x1 | 0x10 | 0x80
        )
        a.template_length, b.template_length = -13, 13
        a.next_reference_id = b.next_reference_id = 0
        a.next_reference_start, b.next_reference_start = s - 111, s
        reads += [a, b]
    _, bam = write_contig(tmp_path, contig, reads, "rf")
    c = _count(bam, _rs.Variant("1", A + 1, contig[A + 1], "A", "SNP"))
    assert (c.dp, c.ad) == (6, 6), (c.dp, c.rd, c.ad)


# ── RJ-11: absent base qualities ───────────────────────────────────────────


def test_a_record_with_absent_qualities_is_dropped_and_warned(tmp_path, caplog):
    contig = _contig("GCT")
    reads = []
    for i in range(4):
        s = A - 50 + i
        seq = list(contig[s : s + 100])
        seq[A + 1 - s] = "A"
        a = make_read(f"q{i}", "".join(seq), s, ((0, 100),))
        a.query_qualities = None  # QUAL '*': the BAM stores 0xFF
        reads.append(a)
    _, bam = write_contig(tmp_path, contig, reads, "qual")
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        c = _count(bam, _rs.Variant("1", A + 1, contig[A + 1], "A", "SNP"))
    assert (c.dp, c.rd, c.ad) == (0, 0, 0)
    warns = [r.getMessage() for r in caplog.records if "base qualities" in r.getMessage()]
    assert len(warns) == 1, warns


# ── RJ-12: unmapped records ────────────────────────────────────────────────


def _snv_with_unmapped_mate(tmp_path, name, mapq0_read=False):
    contig = _contig("GCT")
    reads = [
        make_read(f"r{i}", contig[A - 50 + i : A + 50 + i], A - 50 + i, ((0, 100),))
        for i in range(10)
    ]
    seq = list(contig[A + 1 : A + 101])
    seq[0] = "A"
    # Placed at the variant, MAPQ 0, carrying a CIGAR and the ALT base.
    reads.append(make_read("u", "".join(seq), A + 1, ((0, 100),), flag=0x1 | 0x4 | 0x80, mapq=0))
    if mapq0_read:  # a mapped multi-mapper (pseudogene-style), REF
        reads.append(make_read("m", contig[A - 20 : A + 80], A - 20, ((0, 100),), mapq=0))
    _, bam = write_contig(tmp_path, contig, reads, name)
    return bam, _rs.Variant("1", A + 1, contig[A + 1], "A", "SNP")


def test_an_unmapped_mate_counts_nowhere(tmp_path):
    bam, v = _snv_with_unmapped_mate(tmp_path, "um")
    c = _count(bam, v)
    assert (c.dp, c.rd, c.ad, c.mq0_count) == (10, 10, 0, 0)


def test_an_unmapped_mate_counts_nowhere_at_mapq0(tmp_path):
    bam, v = _snv_with_unmapped_mate(tmp_path, "um0")
    c = _count(bam, v, min_mapq=0)
    assert (c.dp, c.rd, c.ad, c.mq0_count) == (10, 10, 0, 0)


def test_a_mapped_mapq0_read_counts_at_mapq0(tmp_path):
    """Guard: a mapped MAPQ-0 alignment counts at --min-mapq 0, in depth, REF and
    mq0_count (pseudogene loci such as PMS2 rely on it)."""
    bam, v = _snv_with_unmapped_mate(tmp_path, "mq0", mapq0_read=True)
    c = _count(bam, v, min_mapq=0)
    assert c.rd >= 11 and c.mq0_count >= 1, (c.dp, c.rd, c.ad, c.mq0_count)
