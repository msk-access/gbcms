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
  long run the rule fell back to the previous classifier with only a trace. Prep
  now widens it; a variant the rule still cannot judge is warned once.

Committed red (xfail-strict) before the fix. Long complex events' junction
windows read on through the read (C25 #199).
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


def test_carriers_ending_in_the_run_at_a_contig_end_are_uninformative(tmp_path):
    """The same run and reads mid-contig and two bases from the contig end: both
    alleles' readers are uninformative in both places, so neither REF nor ALT."""
    pre, post = _flank(7), _flank(8)
    for name, ref in (("mid", pre + "GAAAAATC" + post), ("end", pre + "GAAAAATC")):
        fa, bam = write_contig(tmp_path, ref, _ending_in_the_run(ref), name)
        (c,) = count_bam_checked(bam, [_prepared(fa, "G", "GA")], [None], *ARGS)
        _invariants(c)
        assert (c.rd, c.ad, c.dp) == (0, 0, 12), name


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


def _traces(caplog, call):
    caplog.clear()
    with caplog.at_level(5, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        out = call()
    return out, [r.getMessage() for r in caplog.records]


def test_an_anchor_changing_insertion_in_a_long_run_is_judged_by_the_exact_carrier_rule(
    tmp_path, caplog
):
    """Prep widens the reference for every variant the exact-carrier rule judges,
    anchor-changing insertions included, so a long run does not send the rule to
    the previous classifier."""
    _, logs = _traces(caplog, lambda: _ins_snv(tmp_path, 80, "run80"))
    assert not [m for m in logs if "previous complex classifier" in m], "no fallback"


def test_an_unjudgeable_complex_variant_is_warned_once(tmp_path, caplog):
    """Without a prepared reference the exact-carrier rule cannot judge C>TA: the
    previous classifier counts its reads, and one warning per count says so."""
    ref = _flank(13) + "C" + "A" * 10 + "G" + _flank(14)
    reads = [
        make_read(f"r{i}", ref[A - 40 + i : A + 60 + i], A - 40 + i, ((0, 100),)) for i in range(4)
    ]
    _, bam = write_contig(tmp_path, ref, reads, "bare")
    v = _rs.Variant("1", A, "C", "TA", "COMPLEX")  # not prepared
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        (c,) = count_bam_checked(bam, [v], [None], *ARGS)
    _invariants(c)
    warns = [r.getMessage() for r in caplog.records if "exact-carrier rule" in r.getMessage()]
    assert len(warns) == 2 and "1:401 C>TA: 4 read(s)" in warns[0], warns
    assert "no prepared reference" in warns[0], "the warning names why the rule could not judge"


def test_an_anchor_changing_insertion_in_a_long_run_judges_reads_by_their_bases(tmp_path):
    """C>TA before an A run: reads keeping the C with the inserted A show it only
    once they reach the run's end. Before a 30bp run they do, and are neither;
    before an 80bp run they end inside it, base for base the REF reads, and count
    REF like them: bases decide, not the CIGAR's insertion (operator, 2026-10-01)."""
    assert _ins_snv(tmp_path, 30, "run30") == (6, 6, 0)
    assert _ins_snv(tmp_path, 80, "run80") == (12, 6, 0)


def test_a_read_spliced_over_an_anchor_changing_insertions_anchor_is_out_of_depth(tmp_path):
    """A>CCC changes its anchor base, so the anchor is what a read must show: an RNA
    read spliced over it, resuming at the next base, observes nothing of the allele
    and leaves depth, as for any delins (operator decision, 2026-10-01)."""
    rng = random.Random(3)
    ref = "".join(rng.choice("GT") for _ in range(900))
    ref = ref[:A] + "A" + ref[A + 1 :]
    reads = []
    for i in range(5):
        s = A - 340 + i
        seq = ref[s : s + 40] + ref[A + 1 : A + 61]
        reads.append(make_read(f"r{i}", seq, s, ((0, 40), (3, A + 1 - (s + 40)), (0, 60))))
    fa, bam = write_contig(tmp_path, ref, reads, "spliced")
    v = _prepared(fa, "A", "CCC")
    (c,) = count_bam_checked(
        bam, [v], [None], 0, 0, True, True, True, False, False, False, 1, mode="rna"
    )
    assert (c.dp, c.splice_skip_excluded) == (0, 5)


@pytest.mark.parametrize("backend", ["sw", "pairhmm"])
def test_a_soft_masked_reference_counts_the_same(tmp_path, backend):
    """A FASTA soft-masked (lower case) over the locus is the same reference: the
    reads, including those judged by alignment against the reference context (a
    deletion that removes the anchor), count the same as against upper case."""
    from helpers import assert_same_counts

    rng = random.Random(5)
    up = "".join(rng.choice("ACGT") for _ in range(900))
    up = up[:A] + "GAAAAAAT" + up[A + 8 :]
    low = up[:300] + up[300:520].lower() + up[520:]
    reads = []
    for i in range(6):
        s = A - 50 + i
        tail = 100 - (A - s)
        reads.append(
            make_read(
                f"d{i}", up[s:A] + up[A + 3 : A + 3 + tail], s, ((0, A - s), (2, 3), (0, tail))
            )
        )
        reads.append(
            make_read(
                f"c{i}",
                up[s : A + 1] + up[A + 2 : A + 1 + tail],
                s,
                ((0, A + 1 - s), (2, 1), (0, tail - 1)),
            )
        )
        reads.append(make_read(f"r{i}", up[s : s + 100], s, ((0, 100),)))
    counts = []
    for name, ref in (("upper", up), ("lower", low)):
        fa, bam = write_contig(tmp_path, ref, reads, name)
        counts.append(
            count_bam_checked(
                bam, [_prepared(fa, "GA", "G")], [None], *ARGS, alignment_backend=backend
            )
        )
    assert_same_counts(counts[1], counts[0], f"soft-masked vs upper case, {backend}")


def test_identical_runs_give_identical_mfsd_bits(tmp_path):
    """The fragment-size statistics are summed in a fixed order: twenty identical
    runs give bit-identical values (they varied in the last digits with the
    fragments' hash order)."""
    import struct

    rng = random.Random(9)
    ref = "".join(rng.choice("ACGT") for _ in range(1200))
    alt = "T" if ref[A] != "T" else "G"
    reads = []
    for n in range(120):
        frag, s = rng.randint(120, 320), A - rng.randint(5, 95)
        seq = ref[s : s + 100]
        if n % 3 == 0:
            seq = seq[: A - s] + alt + seq[A - s + 1 :]
        m = s + frag - 100
        r1 = make_read(f"f{n}", seq, s, ((0, 100),), flag=0x1 | 0x2 | 0x40 | 0x20)
        r2 = make_read(f"f{n}", ref[m : m + 100], m, ((0, 100),), flag=0x1 | 0x2 | 0x80 | 0x10)
        for r, mate, tlen in ((r1, m, frag), (r2, s, -frag)):
            r.next_reference_id, r.next_reference_start, r.template_length = 0, mate, tlen
        reads += [r1, r2]
    fa, bam = write_contig(tmp_path, ref, reads, "mfsd")
    v = _prepared(fa, ref[A], alt)
    fields = ("mfsd_ref_llr", "mfsd_alt_llr", "mfsd_ref_mean", "mfsd_alt_mean", "mfsd_pval_alt_ref")
    seen = set()
    for _ in range(20):
        (c,) = _rs.count_bam_binned(bam, [v], [None], *ARGS, mfsd=True)
        seen.add(tuple(struct.pack("d", getattr(c, f)) for f in fields))
    assert len(seen) == 1, f"{len(seen)} distinct results over 20 runs"


def test_a_decomposed_twin_in_a_long_run_is_judged_by_the_exact_carrier_rule(tmp_path, caplog):
    """Prep gives the homopolymer twin (AAA>C read as AAA>AAC) its own widened
    reference, so in a long run the exact-carrier rule judges the twin's reads too
    rather than falling back, with a warning naming an allele the user never gave."""
    run = 80
    ref = _flank(41) + "G" + "A" * run + "C" + _flank(42)
    end = A + 1 + run  # the C
    pos = end - 3  # REF: the run's last three A's
    reads = []
    for i in range(6):
        s = pos - 60 + i
        seq = ref[s : end - 2] + ref[end : s + 102]
        reads.append(make_read(f"d{i}", seq, s, ((0, end - 2 - s), (2, 2), (0, s + 102 - end))))
        reads.append(make_read(f"r{i}", ref[s : s + 100], s, ((0, 100),)))
    fa, bam = write_contig(tmp_path, ref, reads, "twin")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", pos, "AAA", "C", "X")], fa, 5, False, 1, True, rescue_homopolymer=True
    )
    assert pv.decomposed_variant is not None and pv.decomposed_variant.event_ref is not None
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        count_bam_checked(bam, [pv.variant], [pv.decomposed_variant], *ARGS)
    assert not [r for r in caplog.records if "exact-carrier rule" in r.getMessage()]


@pytest.mark.parametrize(
    "ref_allele, alt_allele, read_hap",
    [
        # C>TA before a 60-A run: the reads carry only C>T and span the run and the G
        # (T + 60 A + G; the ALT is T + 61 A + G).
        ("C", "TA", lambda ref: ref[:A] + "T" + ref[A + 1 :]),
        # CA>T before the same run: the reads carry only C>T, not the deleted A.
        ("CA", "T", lambda ref: ref[:A] + "T" + ref[A + 1 :]),
    ],
    ids=["ins-snv", "del-snv"],
)
def test_a_long_events_left_junction_alone_does_not_make_alt(
    tmp_path, ref_allele, alt_allele, read_hap
):
    """Reads whose bases carry only the anchor substitution, spanning the whole
    run and its far flank, contradict the ALT there: not ALT. The exact-carrier
    rule's long-event junction windows let the left junction alone decide."""
    ref = _flank(61) + "C" + "A" * 60 + "G" + _flank(62, 500)
    s, e = A - 50, A + 62  # the last base read is the G after the run
    hap = read_hap(ref)
    reads = [make_read(f"snv{i}", hap[s:e], s, ((0, e - s),)) for i in range(3)]
    fa, bam = write_contig(tmp_path, ref, reads, "junction")
    (c,) = count_bam_checked(bam, [_prepared(fa, ref_allele, alt_allele)], [None], *ARGS)
    _invariants(c)
    assert c.ad == 0


def test_a_deciding_base_past_the_contig_end_is_depth_only_and_warned(tmp_path, caplog):
    """+A in a run of six A's that ends the contig: carriers show seven A's. REF has
    no base past the end, but a circular contig's continues at its start, so the
    reads stay depth only (operator decision, 2026-10-01) — counted and warned once
    per variant, never silent."""
    ref = _flank(15) + "GAAAAAA"
    reads = []
    for i in range(5):
        s = A - 60 + i
        reads.append(
            make_read(
                f"e{i}",
                ref[s : A + 1] + "A" + ref[A + 1 : A + 7],
                s,
                ((0, A + 1 - s), (1, 1), (0, 6)),
            )
        )
    fa, bam = write_contig(tmp_path, ref, reads, "edge")
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        (c,) = count_bam_checked(bam, [_prepared(fa, "G", "GA")], [None], *ARGS)
    _invariants(c)
    assert (c.ad, c.dp) == (0, 5)
    warns = [r.getMessage() for r in caplog.records if "a contig edge" in r.getMessage()]
    assert len(warns) == 2 and "5 ALT read(s)" in warns[0], warns


def test_the_pass_warns_once_about_rows_it_can_only_judge_degraded(tmp_path, caplog):
    """Rows passed to the engine unprepared are said once per counting pass, with
    their count: indels and complex variants without a prepared reference context,
    and MNPs whose REF equals ALT."""
    ref = _flank(21) + "C" + "A" * 10 + "G" + _flank(22)
    reads = [
        make_read(f"r{i}", ref[A - 40 + i : A + 60 + i], A - 40 + i, ((0, 100),)) for i in range(4)
    ]
    _, bam = write_contig(tmp_path, ref, reads, "degraded")
    rows = [
        _rs.Variant("1", A, "C", "TA", "COMPLEX"),
        _rs.Variant("1", A, "CA", "C", "DEL"),
        _rs.Variant("1", A + 1, "AA", "aa", "MNP"),
    ]
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        for c in _rs.count_bam_binned(bam, rows, [None] * 3, *ARGS):
            _invariants(c)
    msgs = [r.getMessage() for r in caplog.records]
    unprepared = [m for m in msgs if "carry no prepared reference context" in m]
    degenerate = [m for m in msgs if "have REF equal to ALT" in m]
    assert len(unprepared) == 1 and unprepared[0].startswith("2 "), unprepared
    assert len(degenerate) == 1 and degenerate[0].startswith("1 "), degenerate


def test_each_route_to_the_exact_carrier_rule_is_named_per_read(tmp_path, caplog):
    """Every read the exact-carrier rule judges has a trace naming its read, the
    route that sent it there, and the outcome."""
    _, logs = _traces(caplog, lambda: _ins_snv(tmp_path, 10, "route"))
    routed = [m for m in logs if "exact-carrier rule (an insertion whose anchor changes)" in m]
    assert routed and all(" read=" in m and " alt=" in m for m in routed), routed[:3]


def test_why_an_alt_call_is_kept_is_traced(tmp_path, caplog):
    """An ALT call from a read spanning neither ALT-side window stands on its bases,
    and the trace says so: carriers of GA>G in G A*40 T reading from the G to the T."""
    ref = _flank(31) + "G" + "A" * 40 + "T" + _flank(32)
    seq = "G" + "A" * 39 + "T"
    reads = [make_read(f"k{i}", seq, A, ((0, 20), (2, 1), (0, 21))) for i in range(3)]
    fa, bam = write_contig(tmp_path, ref, reads, "kept")
    v = _prepared(fa, "GA", "G")
    (c,), logs = _traces(caplog, lambda: count_bam_checked(bam, [v], [None], *ARGS))
    _invariants(c)
    assert c.ad == 3
    kept = [m for m in logs if "ALT kept — its own bases tell the alleles apart" in m]
    assert kept and all(" read=k" in m for m in kept), kept
    totals = [m for m in logs if "depth reads by deciding rule" in m]
    assert totals and all("ALT kept by its bases 3;" in m for m in totals), totals


def test_a_row_with_an_empty_allele_counts_no_allele_and_is_warned(tmp_path, caplog):
    """Prep rejects empty alleles; a row passed to the engine directly with one shows
    no allele (neither REF nor ALT, no panic), and the pass says so once."""
    ref = _flank(41) + "C" + "A" * 10 + "G" + _flank(42)
    reads = [
        make_read(f"r{i}", ref[A - 40 + i : A + 60 + i], A - 40 + i, ((0, 100),)) for i in range(4)
    ]
    _, bam = write_contig(tmp_path, ref, reads, "empty")
    rows = [_rs.Variant("1", A, "", "T", "X"), _rs.Variant("1", A, "", "", "X")]
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        counts = _rs.count_bam_binned(bam, rows, [None] * 2, *ARGS)
    for c in counts:
        _invariants(c)
        assert (c.rd, c.ad, c.partial_alt) == (0, 0, 0)
    empty = [r.getMessage() for r in caplog.records if "have an empty allele" in r.getMessage()]
    assert len(empty) == 1 and empty[0].startswith("2 "), empty
