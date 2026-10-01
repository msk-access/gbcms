"""Missing data never changes a count silently (#204).

The 6.6.0 code-quality sweep found decision paths where an unavailable reference
changed a call with nothing logged:

- An ALT call on a pure indel from a read that spans neither window stood on the
  CIGAR's placement alone when the reference around the event was unavailable:
  near a contig end the ALT side's reading stretch failed to fetch, so carriers
  whose bases fit both alleles counted ALT while the same REF reads were
  withdrawn. Without any prepared reference, the same, with no warning.
- An insertion whose ALT also changes the anchor base (C>TA) goes to the
  exact-carrier rule, but prep never widened its reference for that rule, so in a
  long run the rule fell back to the previous classifier with only a trace: reads
  that keep the anchor counted REF.

Committed red (xfail-strict) before the fix.
"""

import logging
import random

import pytest
from helpers import count_bam_checked, make_read, write_contig

from gbcms import _rs

A = 400
ENGINE_LOGGER = "_rs.counting.engine"
ARGS = (20, 20, True, True, True, False, False, False, 1)


def _flank(seed, n=400):
    rng = random.Random(seed)
    return "".join(rng.choice("CGT") for _ in range(n))


def _prepared(fa, ref, alt):
    (pv,) = _rs.prepare_variants([_rs.Variant("1", A, ref, alt, "X")], fa, 5, False, 1, True)
    assert pv.gbcms_status == "PASS", pv.gbcms_status_reason
    return pv.variant


def _invariants(c):
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev


def _ending_in_the_run(ref):
    """G>GA in G AAAAA T: carriers inserting the A and REF reads, all ending on
    the fourth A, where both alleles read the same."""
    reads = []
    for i in range(6):
        s = A - 60 + i
        seq = ref[s : A + 1] + "A" + ref[A + 1 : A + 5]
        reads.append(make_read(f"alt{i}", seq, s, ((0, A + 1 - s), (1, 1), (0, 4))))
        reads.append(make_read(f"ref{i}", ref[s : A + 5], s, ((0, A + 5 - s),)))
    return reads


@pytest.mark.xfail(strict=True, reason="H3 S1: red until fixed")
def test_carriers_ending_in_the_run_at_a_contig_end_are_uninformative(tmp_path):
    """The same run and reads mid-contig and two bases from the contig end: both
    alleles' readers are uninformative in both places, so neither REF nor ALT."""
    pre, post = _flank(7), _flank(8)
    for name, ref in (("mid", pre + "GAAAAATC" + post), ("end", pre + "GAAAAATC")):
        fa, bam = write_contig(tmp_path, ref, _ending_in_the_run(ref), name)
        (c,) = count_bam_checked(bam, [_prepared(fa, "G", "GA")], [None], *ARGS)
        _invariants(c)
        assert (c.rd, c.ad, c.dp) == (0, 0, 12), name


@pytest.mark.xfail(strict=True, reason="H3 S1: red until fixed")
def test_an_alt_read_without_a_reference_to_read_is_depth_only_and_warned(tmp_path, caplog):
    """An unprepared +A row (no reference around it): carriers ending with the
    inserted base span neither window, and their bases cannot be read against any
    reference. Depth only, and one warning for the variant."""
    ref = _flank(9) + "GAAAAATC" + _flank(10)
    reads = []
    for i in range(5):
        s = A - 60 + i
        reads.append(make_read(f"alt{i}", ref[s : A + 1] + "A", s, ((0, A + 1 - s), (1, 1))))
    _, bam = write_contig(tmp_path, ref, reads, "bare")
    v = _rs.Variant("1", A, "G", "GA", "INSERTION")  # not prepared
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        (c,) = count_bam_checked(bam, [v], [None], *ARGS)
    _invariants(c)
    assert (c.ad, c.dp) == (0, 5)
    warns = [r.getMessage() for r in caplog.records if "could not be judged" in r.getMessage()]
    assert len(warns) == 2 and "1:401" in warns[0], warns  # once per variant, per count


def _ins_snv(tmp_path, run_len, name):
    """C>TA before an A run: exact carriers (T, inserted A), reads keeping the C
    with the inserted A (another allele), and REF reads."""
    ref = _flank(11) + "C" + "A" * run_len + "G" + _flank(12)
    reads = []
    for i in range(6):
        s = A - 40 + i
        tail = 100 - 2 - (A - s)
        seq = ref[s:A] + "TA" + ref[A + 1 : A + 1 + tail]
        reads.append(make_read(f"alt{i}", seq, s, ((0, A + 1 - s), (1, 1), (0, tail))))
        seq = ref[s : A + 1] + "A" + ref[A + 1 : A + 1 + tail]
        reads.append(make_read(f"keep{i}", seq, s, ((0, A + 1 - s), (1, 1), (0, tail))))
        reads.append(make_read(f"ref{i}", ref[s : s + 100], s, ((0, 100),)))
    fa, bam = write_contig(tmp_path, ref, reads, name)
    (c,) = count_bam_checked(bam, [_prepared(fa, "C", "TA")], [None], *ARGS)
    _invariants(c)
    return c.rd, c.ad, c.partial_alt


@pytest.mark.xfail(strict=True, reason="H3 S2: red until fixed")
def test_an_anchor_changing_insertion_in_a_long_run_is_judged_by_its_whole_allele(tmp_path):
    """The same reads before a 30bp and an 80bp A run count the same: the
    exact-carrier rule judges both (reads keeping the anchor are not REF)."""
    short = _ins_snv(tmp_path, 30, "run30")
    assert short[0] == 6, short  # only the REF reads
    assert _ins_snv(tmp_path, 80, "run80") == short
