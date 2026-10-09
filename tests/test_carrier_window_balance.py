"""Complex variants: windows of equal length at one anchor judge REF and ALT molecules
alike (docs/reference/read-judgment.md, RJ-23).

A heterozygous locus read at every start: one REF read and one ALT read starting at
each position along their own haplotypes, all the same length. Whatever reads the
rule leaves as depth, it must leave as many of each allele, so AD equals RD exactly:
a window of one allele that a read can hold from more (or fewer) starts than the
other's biases the VAF toward that allele. Each read is the forward read of a
proper pair, so its fragment admits the ALT bases it holds in a soft clip, as a
real pair's does. Synthetic and PHI-free.
"""

import random
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
from helpers import make_read, write_contig  # noqa: E402

from gbcms import _rs  # noqa: E402

A, L, READ = 400, 900, 100
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
# (reference bases replaced, ALT): shrinking, growing, a long shrinking one (windows
# past 50 bases) and a small one.
SHAPES = {
    "28>6": (28, "ATATGA"),
    "6>28": (6, "ATATGAATTATAGGATAATTGTAAATAT"),
    "35>1": (35, "A"),
    "96>2": (96, "AT"),
    "3>5": (3, "ATTAA"),
}


def _contig() -> str:
    rng = random.Random(253)
    return "".join(rng.choice("CGT") for _ in range(L))


def _reads(contig: str, n: int, alt: str):
    """One REF and one ALT read at every start along each haplotype that reaches the
    event's neighbourhood; a read starting inside the ALT event is written with the
    event bases it holds soft-clipped, as an aligner without the event would."""
    reads, k = [], 0
    for s in range(A - READ - 10, A + n + 20):
        seq = contig[s : s + READ]
        reads.append(make_read(f"ref{k}", seq, s, ((0, READ),)))
        k += 1
    hap = contig[:A] + alt + contig[A + n :]
    for h in range(A - READ - 10, A + len(alt) + 20):
        seq = hap[h : h + READ]
        if h + READ <= A:
            cig, pos = ((0, READ),), h
        elif h < A and A - h + len(alt) >= READ:  # ends inside the ALT bases: clipped
            before = A - h
            cig, pos = ((0, before), (4, READ - before)), h
        elif h < A:  # crosses the event: the insertion and the deletion as written
            before = A - h
            cig, pos = ((0, before), (1, len(alt)), (2, n), (0, READ - before - len(alt))), h
        elif h < A + len(alt):  # starts inside the ALT bases: they are clipped
            clip = A + len(alt) - h
            cig, pos = ((4, clip), (0, READ - clip)), A + n
        else:
            cig, pos = ((0, READ),), h - len(alt) + n
        cig = tuple((op, ln) for op, ln in cig if ln > 0)
        reads.append(make_read(f"alt{k}", seq, pos, cig))
        k += 1
    for r in reads:
        _pair(r)
    return reads


def _pair(read) -> None:
    """The forward, first read of a proper pair whose fragment runs well past it."""
    read.flag = 0x1 | 0x2 | 0x20 | 0x40
    read.next_reference_id = 0
    read.next_reference_start = read.reference_start + 250
    read.template_length = 400


@pytest.mark.parametrize("shape", sorted(SHAPES))
def test_equal_windows_judge_ref_and_alt_molecules_alike(tmp_path, shape):
    n, alt = SHAPES[shape]
    contig = _contig()
    fa, bam = write_contig(tmp_path, contig, _reads(contig, n, alt), "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A, contig[A : A + n], alt, "X")], fa, 5, False, 1, True
    )
    (c,) = _rs.count_bam_binned(bam, [pv.variant], [None], **ARGS)
    assert c.dp >= c.rd + c.ad
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    assert c.ad > 0, f"{shape}: no ALT read judged"
    assert c.ad == c.rd, f"{shape}: AD {c.ad} vs RD {c.rd} from equal molecules"


def test_a_carrier_whose_deletion_covers_the_reading_flank_is_read_from_the_other(tmp_path):
    """A 30-to-1 delins whose left flank `CA` recurs at the end of the REF allele: an
    aligner may delete the flank and the first 28 REF bases, align the carrier's `CA`
    to the REF's last two, and write the ALT `G` as a mismatch (as BWA did on an 82-to-1
    panel delins). The reading flank is under the read's own deletion, so the read is
    placed from the other flank at each allele's length, as a read whose clip holds
    the flank is: its bases hold the whole ALT window, so it is ALT."""
    n = 30
    c = list(_contig())
    c[A - 2 : A] = "CA"
    c[A] = "T"
    c[A + n - 3 : A + n + 3] = "CAAGTT"
    contig = "".join(c)
    hap = contig[:A] + "G" + contig[A + n :]
    reads = []
    for i in range(6):
        s = A - 90 + i
        seq = hap[s : A + 4]  # ends three bases past the ALT: CA G GTT
        reads.append(make_read(f"alt{i}", seq, s, ((0, A - 2 - s), (2, n - 1), (0, 6))))
    reads += [
        make_read(f"ref{i}", contig[s : s + READ], s, ((0, READ),))
        for i, s in enumerate(range(A - 60, A - 54))
    ]
    for r in reads:
        _pair(r)
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A, contig[A : A + n], "G", "X")], fa, 5, False, 1, True
    )
    (got,) = _rs.count_bam_binned(bam, [pv.variant], [None], **ARGS)
    assert got.rd == got.rd_fwd + got.rd_rev and got.ad == got.ad_fwd + got.ad_rev
    assert got.dp >= got.rd + got.ad
    assert (got.rd, got.ad) == (6, 6)
