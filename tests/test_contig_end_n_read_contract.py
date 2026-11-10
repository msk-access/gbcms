"""A read that calls no base over a row's window is depth only, at a contig end too.

At a deletion whose repeat runs to the end of its contig, the prepared reference
cannot hold the event with its flank, so the exact-carrier rule cannot judge the
row's reads and the CIGAR checks classify them. Their reference-coverage rules
(a pure deletion's anchor covered with no gap; a deletion-direction delins' anchor
covered by a clean match) rest on placement, so a read whose bases there are all N
was counted REF. An N is no base call: such a read counts toward depth (and the N
class) only, as the exact-carrier rule counts it away from the contig end.
"""

import random

import pytest
from helpers import count_checked, make_read, write_contig

from gbcms import _rs as gbcms_rs

READ_LEN = 100
# 692 seeded bases, then a GC repeat that runs to the contig's last base: deleting
# GC after the T at 0-based 694 is TGC>T, 5 bp from the end of the 700-bp contig.
CONTIG = "".join(random.Random(4).choice("ACGT") for _ in range(692)) + "CGTGCGCG"
POS0 = 694
START = len(CONTIG) - READ_LEN  # reads end on the contig's last base


def _counts(tmp_path, reads, ref_allele, alt_allele):
    fasta, bam = write_contig(tmp_path, CONTIG, reads)
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", POS0, ref_allele, alt_allele, "X")], fasta, 5, False, 1, True
    )
    assert pv.gbcms_status == "PASS", pv.gbcms_status_reason
    (counts,) = count_checked(bam, [pv.variant], decomposed=[pv.decomposed_variant])
    assert counts.dp >= counts.rd + counts.ad
    assert counts.dpf >= counts.rdf + counts.adf
    assert counts.rd == counts.rd_fwd + counts.rd_rev
    assert counts.ad == counts.ad_fwd + counts.ad_rev
    return counts


def _all_n():
    return make_read("all_n", "N" * READ_LEN, START, ((0, READ_LEN),))


def _n_over_the_window():
    # Its own bases up to 690, N over the last ten (the whole window and the anchor).
    seq = CONTIG[START:690] + "N" * (len(CONTIG) - 690)
    return make_read("window_n", seq, START, ((0, READ_LEN),))


def _reference():
    return make_read("ref", CONTIG[START:], START, ((0, READ_LEN),))


ROWS = [
    pytest.param("TGC", "T", id="deletion"),
    pytest.param("TG", "C", id="deletion-direction-delins"),
]


def test_the_row_sits_where_the_rule_cannot_judge():
    assert CONTIG[POS0 : POS0 + 3] == "TGC" and len(CONTIG) - POS0 == 6


@pytest.mark.parametrize("ref_allele,alt_allele", ROWS)
@pytest.mark.parametrize("read", [_all_n, _n_over_the_window], ids=["all-N", "N-over-window"])
def test_a_read_calling_no_base_over_the_window_is_depth_only(
    tmp_path, ref_allele, alt_allele, read
):
    c = _counts(tmp_path, [read()], ref_allele, alt_allele)
    assert (c.dp, c.rd, c.ad) == (1, 0, 0), (c.dp, c.rd, c.ad)
    assert (c.dpf, c.rdf, c.adf) == (1, 0, 0), (c.dpf, c.rdf, c.adf)
    assert c.n_count == 1  # in the N class, as the exact-carrier rule puts it


@pytest.mark.parametrize("ref_allele,alt_allele", ROWS)
def test_a_reference_read_there_is_still_ref(tmp_path, ref_allele, alt_allele):
    c = _counts(tmp_path, [_reference()], ref_allele, alt_allele)
    assert (c.dp, c.rd, c.ad) == (1, 1, 0), (c.dp, c.rd, c.ad)
    assert c.n_count == 0
