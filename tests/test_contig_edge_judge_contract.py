"""The sequence judge (`--read-judge sequence`, #304) at a contig edge (RJ-3).

Mechanism. Past a contig's edge the reference is not known: a circular contig runs on at
its start, so a base a read places past the edge tells no allele from another, and a
read whose deciding base lies there is depth only, counted and warned (H3 #204, operator
2026-10-01). Before this the judge's candidate sequences simply stopped at the edge, so
a read running past it paid a gap on the allele the contig ended in and counted for the
other.

A read reaching the edge can only be called for the longer allele (the shorter one fits
it as a read ending early), and that call rests on the shorter allele holding no base
past the edge. At an insertion that call is ALT: the judge withdraws it to depth and
warns, as the rules do. At a deletion it is REF, and two decided tests count such REF
reads REF (`test_contig_end_n_read_contract`, the MAF deletion at Start 1 in
`test_input_group4_contract`), so the judge keeps them; the rules' exact-carrier route
makes them depth. That difference is open for the operator, not pinned here.

Shapes: an insertion and a deletion in a run that ends the contig, read on both
strands.
"""

import logging
import random

import pytest
from helpers import count_bam_checked, make_read, write_contig

from gbcms import _rs

ENGINE_LOGGER = "_rs.counting.engine"
ARGS = (20, 20, True, True, True, False, False, False, 1)


def _flank(n, seed=204):
    rng = random.Random(seed)
    return "".join(rng.choice("CGT") for _ in range(n))


def _edge_counts(tmp_path, caplog, ref_al, alt_al, carriers_ins):
    """`ref_al` > `alt_al` at a G before A6 that ends a contig of 400 bases; 5 carriers
    (their indel written right after the G) and 5 REF reads, all reading to the edge."""
    contig = _flank(393) + "GAAAAAA"
    a = 393  # the G
    reads = []
    for i in range(5):
        s = a - 60 + i
        flag = 16 if i % 2 else 0
        if carriers_ins:
            seq = contig[s : a + 1] + "A" + contig[a + 1 :]
            cig = ((0, a + 1 - s), (1, 1), (0, 6))
        else:
            seq = contig[s : a + 1] + contig[a + 2 :]
            cig = ((0, a + 1 - s), (2, 1), (0, 5))
        reads.append(make_read(f"c{i}", seq, s, cig, flag=flag))
        reads.append(make_read(f"r{i}", contig[s:], s, ((0, len(contig) - s),), flag=flag))
    fa, bam = write_contig(tmp_path, contig, reads, "edge")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", a, ref_al, alt_al, "INDEL")], fa, 5, False, 1, True
    )
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        (c,) = count_bam_checked(
            bam, [pv.variant], [None], *ARGS, read_judge="sequence", reference_fasta=fa
        )
    assert c.dp >= c.rd + c.ad and c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    warns = [r.getMessage() for r in caplog.records if "a contig edge" in r.getMessage()]
    return c, warns


def test_an_insertions_reads_running_to_a_contig_end_are_depth_only_and_warned(tmp_path, caplog):
    """Every read reaches the edge, where the base telling the run's length would lie:
    none is REF or ALT, all are depth, and the carriers are warned about."""
    c, warns = _edge_counts(tmp_path, caplog, "G", "GA", True)
    assert (c.rd, c.ad, c.dp) == (0, 0, 10)
    assert len(warns) == 2 and "5 ALT read(s)" in warns[0], warns


def test_a_deletions_carriers_running_to_a_contig_end_are_never_alt(tmp_path, caplog):
    """Guard: the carriers fit REF as reads ending a base early, so none is ALT; every
    read is in depth."""
    c, _ = _edge_counts(tmp_path, caplog, "GA", "G", False)
    assert (c.ad, c.dp) == (0, 10)


def test_an_insertion_in_a_run_starting_the_contig_is_depth_only_and_warned(tmp_path, caplog):
    """The left edge: +A in a run of six A's that starts the contig. Carriers (the
    insert first) and REF reads start on the contig's first base, where the base telling
    the run's length would lie before it: all are depth, and the carriers are warned
    about."""
    contig = "AAAAAAG" + _flank(393)
    reads = []
    for i in range(5):
        flag = 16 if i % 2 else 0
        reads.append(
            make_read(f"c{i}", "A" + contig[: 60 + i], 0, ((1, 1), (0, 60 + i)), flag=flag)
        )
        reads.append(make_read(f"r{i}", contig[: 60 + i], 0, ((0, 60 + i),), flag=flag))
    fa, bam = write_contig(tmp_path, contig, reads, "edge")
    (pv,) = _rs.prepare_variants([_rs.Variant("1", 6, "G", "AG", "INDEL")], fa, 5, False, 1, True)
    assert (pv.variant.pos, pv.variant.ref_allele, pv.variant.alt_allele) == (0, "A", "AA")
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        (c,) = count_bam_checked(
            bam, [pv.variant], [None], *ARGS, read_judge="sequence", reference_fasta=fa
        )
    assert c.dp >= c.rd + c.ad and c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    assert (c.rd, c.ad, c.dp) == (0, 0, 10)
    warns = [r.getMessage() for r in caplog.records if "a contig edge" in r.getMessage()]
    assert len(warns) == 2 and "5 ALT read(s)" in warns[0], warns


def test_alike_carriers_at_a_contig_end_are_each_depth_only_and_warned(tmp_path, caplog):
    """Reads alike in bases, qualities and placement (one start: duplex copies) share one
    judged verdict per row; the edge is weighed on each of them, so every carrier is depth
    only and warned, as when each is judged afresh."""
    contig = _flank(393) + "GAAAAAA"
    a = 393
    s = a - 60
    seq = contig[s : a + 1] + "A" + contig[a + 1 :]
    reads = [make_read(f"c{i}", seq, s, ((0, a + 1 - s), (1, 1), (0, 6))) for i in range(5)]
    reads += [make_read(f"r{i}", contig[s:], s, ((0, len(contig) - s),)) for i in range(5)]
    fa, bam = write_contig(tmp_path, contig, reads, "edge")
    (pv,) = _rs.prepare_variants([_rs.Variant("1", a, "G", "GA", "INDEL")], fa, 5, False, 1, True)
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        (c,) = count_bam_checked(
            bam, [pv.variant], [None], *ARGS, read_judge="sequence", reference_fasta=fa
        )
    assert c.dp >= c.rd + c.ad and c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    assert (c.rd, c.ad, c.dp) == (0, 0, 10)
    warns = [r.getMessage() for r in caplog.records if "a contig edge" in r.getMessage()]
    assert len(warns) == 2 and "5 ALT read(s)" in warns[0], warns


# ── an indel sitting at the contig's first or last base (the edge class, operator 2026-10-29) ──
#
# The class is a pure indel whose placements touch a contig edge, its repeat running to it
# or the indel itself sitting there: an insertion before the first base (a MAF insertion at
# Start 0, VCF's base-after form at position 1) or after the last one, a deletion of the
# first or last bases. A read reaching that edge can be called only for the longer allele,
# resting on the shorter one holding no base past the edge. An insertion's carriers there are
# depth only, counted and warned (RJ-3); a deletion's REF reads stand (decided tests).

INS2 = "GT"


def _edge_contig():
    return "C" + _flank(398, seed=219) + "C"


def _count_at(tmp_path, caplog, contig, reads, pos, ref_al, alt_al):
    fa, bam = write_contig(tmp_path, contig, reads, "edge")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", pos, ref_al, alt_al, "INDEL")], fa, 5, False, 1, True
    )
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        (c,) = count_bam_checked(
            bam, [pv.variant], [None], *ARGS, read_judge="sequence", reference_fasta=fa
        )
    assert c.dp >= c.rd + c.ad and c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    warns = [r.getMessage() for r in caplog.records if "a contig edge" in r.getMessage()]
    return c, warns


@pytest.mark.parametrize("alike", [False, True], ids=["distinct", "alike"])
def test_an_insertion_before_the_contigs_first_base_is_depth_only_and_warned(
    tmp_path, caplog, alike
):
    """The insertion at Start 0 (VCF: position 1, REF the first base, ALT the inserted bases
    then it): carriers open with the inserted bases, written as an insertion or a soft clip;
    REF reads start on the first base. Each read's deciding bases lie before the contig:
    every read is depth only, and the carriers are warned about. Alike carriers (one length,
    one set of qualities) share the judge's verdict and are each withdrawn."""
    contig = _edge_contig()
    n = 100
    reads = []
    for i in range(3):
        k = n if alike else n - i
        seq = INS2 + contig[: k - len(INS2)]
        reads.append(make_read(f"i{i}", seq, 0, ((1, len(INS2)), (0, k - len(INS2)))))
        reads.append(make_read(f"s{i}", seq, 0, ((4, len(INS2)), (0, k - len(INS2)))))
        reads.append(make_read(f"r{i}", contig[:k], 0, ((0, k),)))
        reads.append(make_read(f"q{i}", contig[:k], 0, ((0, k),)))
    c, warns = _count_at(tmp_path, caplog, contig, reads, 0, contig[0], INS2 + contig[0])
    assert (c.rd, c.ad, c.dp) == (0, 0, 12)
    assert len(warns) == 2 and "6 ALT read(s)" in warns[0], warns


def test_an_insertion_after_the_contigs_last_base_is_depth_only_and_warned(tmp_path, caplog):
    """The mirror: the inserted bases after the last base; carriers end with them."""
    contig = _edge_contig()
    L = len(contig)
    reads = []
    for i in range(4):
        s = L - 98 + i
        k = L - s
        reads.append(make_read(f"i{i}", contig[s:] + INS2, s, ((0, k), (1, len(INS2)))))
        reads.append(make_read(f"r{i}", contig[s:], s, ((0, k),)))
    c, warns = _count_at(tmp_path, caplog, contig, reads, L - 1, contig[-1], contig[-1] + INS2)
    assert (c.rd, c.ad, c.dp) == (0, 0, 8)
    assert len(warns) == 2 and "4 ALT read(s)" in warns[0], warns


def test_a_deletion_of_the_contigs_first_bases_keeps_its_ref_reads(tmp_path, caplog):
    """Guard: the first two bases deleted (VCF: position 1, REF them and the base after,
    ALT that base). REF reads starting on the first base stand REF (the MAF deletion at
    Start 1, `test_input_group4_contract.py`); none is warned about."""
    contig = _edge_contig()
    reads = [make_read(f"r{i}", contig[:100], 0, ((0, 100),)) for i in range(6)]
    c, warns = _count_at(tmp_path, caplog, contig, reads, 0, contig[:3], contig[2])
    assert (c.rd, c.ad, c.dp) == (6, 0, 6)
    assert warns == []
