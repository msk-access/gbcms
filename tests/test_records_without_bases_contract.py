"""A record stored without its bases (SEQ '*') is not an observation (#172).

BAM allows a record with no sequence (and no qualities): aligners write
secondary alignments that way, and a stripped primary can look the same. The
SNP check indexed the empty sequence and panicked, heuristic BAQ sliced its
empty qualities and panicked, and the indel checks counted it REF (and depth)
from its CIGAR alone. Such a record cannot show any allele, so it is dropped
where every read enters counting (the shared read filter), in every path:
counts equal those of the same BAM without it.

Committed red (xfail-strict) before the implementation; flipped green with it.
"""

import logging

import pysam
import pytest
from helpers import count_checked, make_read
from test_complex_exact_contract import POS, READ, _files, _prepared, _ref

from gbcms import _rs as gbcms_rs

ENGINE_LOGGER = "_rs.counting.engine"

_REF = _ref()


def _alt_base(b):
    return "A" if b != "A" else "C"


# (REF, ALT) at POS for every variant type the engine dispatches on.
SHAPES = {
    "SNV": (_REF[POS], _alt_base(_REF[POS])),
    "MNP": (_REF[POS : POS + 2], _alt_base(_REF[POS]) + _alt_base(_REF[POS + 1])),
    "insertion": (_REF[POS], _REF[POS] + "GT"),
    "deletion": (_REF[POS : POS + 3], _REF[POS]),
    "delins": (_REF[POS : POS + 2], "GCC" if _REF[POS] != "G" else "TCC"),
}


def _ref_reads(n=10):
    return [
        make_read(f"r{i}", _REF[s : s + READ], s, ((0, READ),))
        for i, s in enumerate(range(POS - 60, POS - 60 + n))
    ]


def _no_bases(flag=0, cigar=((0, READ),), name="noseq"):
    """A record over POS with no SEQ and no QUAL."""
    x = pysam.AlignedSegment()
    x.query_name, x.flag, x.reference_id, x.reference_start = name, flag, 0, POS - 40
    x.mapping_quality, x.cigartuples = 60, cigar
    return x


def _fields(c):
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    return (c.dp, c.rd, c.ad, c.dpf, c.rdf, c.adf, c.partial_alt, c.mq0_count)


def _pair(tmp_path, shape, extra, **kw):
    """(counts without, counts with) the extra records, each checked for binning
    invariance (count_checked)."""
    ref, alt = SHAPES[shape]
    out = []
    for tag, reads in (("without", _ref_reads()), ("with", _ref_reads() + extra)):
        d = tmp_path / tag
        d.mkdir()
        fa, bam = _files(d, _REF, reads)
        (c,) = count_checked(bam, [_prepared(fa, ref, alt)], **kw)
        out.append(c)
    return out


# SNV: panicked. Insertion, deletion: counted REF and depth from the CIGAR alone.
@pytest.mark.parametrize("shape", list(SHAPES))
def test_a_primary_record_without_bases_is_not_counted(tmp_path, shape):
    without, with_ = _pair(tmp_path, shape, [_no_bases()])
    assert _fields(with_) == _fields(without)
    assert (with_.rd, with_.ad) == (10, 0)


@pytest.mark.parametrize("shape", list(SHAPES))
def test_a_secondary_record_without_bases_is_not_counted_when_secondaries_are_kept(tmp_path, shape):
    """With --no-filter-secondary a secondary alignment reaches fragment evidence;
    aligners store secondaries without SEQ. Its primary maps elsewhere (a QNAME
    no read here shares), so a call on it would be a fragment of its own. (The
    insertion and deletion checks counted it a REF fragment from its CIGAR.)"""
    without, with_ = _pair(tmp_path, shape, [_no_bases(flag=0x100)], filter_secondary=False)
    assert _fields(with_) == _fields(without)


@pytest.mark.parametrize("shape", list(SHAPES))
def test_rna_with_baq_skips_a_spliced_record_without_bases(tmp_path, shape):
    """RNA mode with heuristic BAQ (the binned path only): BAQ walks the record's
    qualities at its splice junction."""
    ref, alt = SHAPES[shape]
    spliced = _no_bases(cigar=((0, 60), (3, 200), (0, READ - 60)))
    counts = []
    for tag, reads in (("without", _ref_reads()), ("with", _ref_reads() + [spliced])):
        d = tmp_path / tag
        d.mkdir()
        fa, bam = _files(d, _REF, reads)
        (c,) = gbcms_rs.count_bam_binned(
            bam,
            [_prepared(fa, ref, alt)],
            [None],
            20,
            20,
            True,
            True,
            True,
            False,
            False,
            False,
            1,
            mode="rna",
            apply_baq=True,
            mfsd=True,
        )
        counts.append(c)
    assert _fields(counts[1]) == _fields(counts[0])


def test_skipped_records_are_warned_once_per_counting_pass(tmp_path, caplog):
    """Dropping them is said once, above DEBUG: a BAM stripped of its sequences
    would otherwise count nothing silently."""
    ref, alt = SHAPES["delins"]
    fa, bam = _files(tmp_path, _REF, _ref_reads() + [_no_bases(name="a"), _no_bases(name="b")])
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        gbcms_rs.count_bam_binned(
            bam, [_prepared(fa, ref, alt)], [None], 20, 20, True, True, True, False, False, False, 1
        )
    warns = [r.message for r in caplog.records if "without bases" in r.message]
    assert len(warns) == 1 and "in 2 bin fetch" in warns[0], warns
