"""`MNV_SHARED` counts only reads that spell the longer MNV whole, confidently (#304, RJ-28).

Mechanism. Under the sequence judge an MNV listed beside a shorter MNP or SNV that holds its
change is that row's ALT, and the row's diagnostic names how many of its ALT reads spell the
MNV whole. An ALT call does not guarantee that the row's own changed bases were read
confidently: an MNP read with a changed base below `--min-baseq` counts ALT (RJ-25). So a
read counts as spelling the MNV only where every base the MNV changes, the row's own
included, is read at or above `--min-baseq` as the MNV's. Before this the row's own bases were
taken as read, and such a read was named as sharing the MNV.

Shape: CC>GT with a co-listed CCA>GTC; six GTC carriers read confidently, three with the T
at Q15, five REF reads. All nine carriers are the row's ALT reads; how many of them spell
CCA>GTC whole depends on `--min-baseq` (a user setting, swept).
"""

import random

import pytest
from helpers import count_bam_checked, make_read, write_contig

from gbcms import _rs

READ = 100
P = 300  # 0-based first base of the MNP


def _contig():
    rng = random.Random(363)
    bases = [rng.choice("ACGT") for _ in range(700)]
    bases[P - 1] = "T"
    bases[P : P + 3] = list("CCA")
    bases[P + 3] = "T"
    return "".join(bases)


def _counts(tmp_path, min_baseq):
    contig = _contig()
    hap = contig[:P] + "GTC" + contig[P + 3 :]
    reads = []
    for i in range(6):
        s = P - 40 + i
        reads.append(
            make_read(f"a{i}", hap[s : s + READ], s, ((0, READ),), flag=16 if i % 2 else 0)
        )
    for i in range(3):
        s = P - 45 + i
        quals = [30] * READ
        quals[P + 1 - s] = 15
        reads.append(make_read(f"q{i}", hap[s : s + READ], s, ((0, READ),), quals=quals))
    for i in range(5):
        s = P - 35 + i
        reads.append(make_read(f"r{i}", contig[s : s + READ], s, ((0, READ),)))
    fa, bam = write_contig(tmp_path, contig, reads, f"q{min_baseq}")
    rows = [_rs.Variant("1", P, "CC", "GT", "MNP"), _rs.Variant("1", P, "CCA", "GTC", "MNP")]
    prepared = [pv.variant for pv in _rs.prepare_variants(rows, fa, 5, False, 1, True)]
    counts = count_bam_checked(
        bam,
        prepared,
        [None, None],
        20,
        min_baseq,
        True,
        True,
        True,
        False,
        False,
        False,
        1,
        read_judge="sequence",
        reference_fasta=fa,
    )
    for c in counts:
        assert c.dp >= c.rd + c.ad and c.dpf >= c.rdf + c.adf
        assert c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    return counts[0]


@pytest.mark.parametrize(
    "min_baseq, spelled", [(10, 9), (15, 9), (16, 6), (20, 6)], ids=lambda v: str(v)
)
def test_a_read_spells_the_listed_mnv_only_where_every_base_it_changes_is_confident(
    tmp_path, min_baseq, spelled
):
    """The row counts all nine carriers ALT; the three whose T is at Q15 spell CCA>GTC
    whole only where Q15 is at or above `--min-baseq`."""
    c = _counts(tmp_path, min_baseq)
    assert (c.rd, c.ad) == (5, 9)
    assert c.mnv_shared == f"{P + 1}:CCA>GTC:{spelled}/9"
