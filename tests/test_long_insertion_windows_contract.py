"""The sequence judge at an insertion of 50 bp or more: RJ-23's equal windows from the anchor.

Mechanism (operator, 2026-10-30). A read is not expected to hold a long insertion whole, so
both alleles are read over windows of one length from the anchor: REF's is the anchor, the
run the insert slides over and the base past it; X's is the anchor and as many of X's bases.
A REF or ALT call holds its allele's window. A read lying inside the insert, or holding only
its far junction, holds none and is depth only: counting it would give the longer allele read
starts the shorter one cannot have. Before this the judge counted such reads ALT.

A change the guard masks past the run, where the reads written straight through the insertion
show X's own base (a haplotype it found carries it), stays masked (RJ-27, #294). On X's
haplotype the mask lands past the insert, so read there the base would tell REF's molecules
apart and never X's: both windows read on to the next base telling them apart, and a read short
of it is depth only, a REF read as a carrier (RJ-27's read-on, decided 2026-10-30 for both
alleles alike). The `cis`
option (`GBCMS_JUDGE_LONG_INS=cis`, measurement only, not adopted: it credits such a carrier
on other fragments' evidence, the cis rescue of REJ-20261023-001) reads that base on REF's
haplotypes as the reference's where three or more REF molecules over the anchor and past it,
told REF by their bases, read the reference's base there and none another; with a REF
molecule carrying the change (a germline change on REF's haplotype) the mask stays.

Shape: an 89 bp duplication after the anchor at 0-based 500, a copy of the 89 bases after it
but its last, so the run is 501-588 and 589 tells the alleles apart. The duplication's
haplotype carries a change at 589 in cis (its second copy is the first's twin), as HG002's
duplications at chr21:41865095 and chr22:39653046 do.
"""

import collections
import random
import re

import pysam
from helpers import count_bam_checked, make_read, write_contig

from gbcms import _rs

ENGINE_LOGGER = "_rs.counting.engine"
ARGS = (20, 20, True, True, True, False, False, False, 1)
CALL = re.compile(
    r"read call \S+ \S+ read=([a-z]+)\d+ mate=(\d) ref=(\w+) alt=(\w+) .*?partial=(\d+)"
)
A = 500  # the anchor (0-based)
N = 89  # the insert


def _shape(germline: bool):
    """The contig, the row's ALT, and the reads: REF molecules over the anchor and past 589
    (`ref`), REF reads elsewhere (`refb`), carriers written with the insert (`thr`), carriers
    over the anchor ending at 589 written straight (`end`), REF reads over the anchor ending at
    589 (`rend`), the second copy's reads placed
    apart with the change at 589 (`two`), and fragments starting inside the insert with its
    tail clipped, placed after the anchor (`clip`, read 1 of a proper pair). With `germline`,
    three REF molecules carrying the change at 589 too (`germ`)."""
    rng = random.Random(264)
    contig = "".join(rng.choice("ACGT") for _ in range(1200))
    unit = contig[A + 1 : A + 1 + N]
    last = next(b for b in "ACGT" if b != contig[A + N] and b != contig[A])
    ins = unit[:-1] + last
    dup = contig[: A + 1] + ins + ins + contig[A + 1 + N :]
    reads = []

    def add(name, seq, start, cigar, **kw):
        reads.append(make_read(name, seq, start, cigar, **kw))
        return reads[-1]

    for i in range(8):
        s = 468 + 3 * i
        add(f"ref{i}", contig[s : s + 151], s, ((0, 151),))
    for i in range(8):
        s = 380 + 25 * i
        add(f"refb{i}", contig[s : s + 151], s, ((0, 151),))
    for i in range(4):
        s = 455 + 5 * i
        m1 = A + 1 - s
        add(f"thr{i}", dup[s : s + 151], s, ((0, m1), (1, N), (0, 151 - m1 - N)))
    for i in range(4):
        s = 439 + 2 * i
        seq = dup[s : A + 1 + N]
        add(f"end{i}", seq, s, ((0, len(seq)),))
    for i in range(3):
        s = 440 + 2 * i
        add(f"rend{i}", contig[s : A + 1 + N], s, ((0, A + 1 + N - s),))
    for i in range(6):
        o = 10 + 10 * i
        add(f"two{i}", dup[A + 1 + N + o : A + 1 + N + o + 151], A + 1 + o, ((0, 151),))
    for i in range(2):
        o = 78 + 2 * i
        r1 = add(
            f"clip{i}",
            dup[A + 1 + o : A + 1 + o + 151],
            A + 1,
            ((4, N - o), (0, 151 - N + o)),
            flag=99,
        )
        mate = 820 + 5 * i
        r2 = add(f"clip{i}", contig[mate : mate + 151], mate, ((0, 151),), flag=147)
        span = mate + 151 - (A + 1 - (N - o))
        r1.next_reference_id, r1.next_reference_start, r1.template_length = 0, mate, span
        r2.next_reference_id, r2.next_reference_start, r2.template_length = 0, A + 1, -span
    if germline:
        g = contig[: A + N] + last + contig[A + N + 1 :]
        for i in range(3):
            s = 470 + 4 * i
            add(f"germ{i}", g[s : s + 151], s, ((0, 151),))
    return contig, contig[A] + ins, reads


def _run(tmp_path, caplog, germline: bool):
    contig, alt, reads = _shape(germline)
    fa, bam = write_contig(tmp_path, contig, reads, "dup")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A, contig[A], alt, "INDEL")], fa, 5, False, 1, True, read_judge="sequence"
    )
    caplog.clear()
    with caplog.at_level(5, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        (c,) = count_bam_checked(
            bam, [pv.variant], [None], *ARGS, read_judge="sequence", reference_fasta=fa
        )
    assert c.dp >= c.rd + c.ad and c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    calls = {}
    for record in caplog.records:
        if m := CALL.search(record.getMessage()):
            kind, mate, ref, alt_, partial = m.groups()
            call = (
                "ALT"
                if alt_ == "true"
                else "REF" if ref == "true" else "partial" if partial != "0" else "depth"
            )
            calls[(record.getMessage().split(" read=")[1].split()[0], mate)] = (kind, call)
    tally = collections.Counter(calls.values())
    return c, tally


def test_a_long_insertions_reads_are_judged_over_equal_windows_from_its_anchor(tmp_path, caplog):
    """Carriers written with the insert are ALT; fragments starting inside the insert are
    depth only, though their pair admits them; reads ending at the masked base past the run,
    carriers and REF reads alike, are short of the windows read on and are depth only; REF
    molecules over the anchor and past 589 are REF."""
    c, tally = _run(tmp_path, caplog, germline=False)
    assert tally[("thr", "ALT")] == 4, tally
    assert tally[("clip", "depth")] == 2 and tally[("clip", "ALT")] == 0, tally
    assert tally[("end", "depth")] == 4 and tally[("rend", "depth")] == 3, tally
    assert tally[("ref", "REF")] == 8, tally
    assert not any(call == "ALT" for (kind, call) in tally if kind in ("ref", "refb", "two")), tally
    assert (c.rd, c.ad) == (10, 4)


def test_the_cis_option_reads_a_change_the_ref_molecules_lack_as_the_references_on_ref(
    tmp_path, caplog, monkeypatch
):
    """`GBCMS_JUDGE_LONG_INS=cis` (measurement only): the REF molecules over the anchor and
    past 589 read the reference's base there, so REF's haplotypes do, and the carriers ending
    at 589 are ALT and the REF reads ending there REF; the rest stand."""
    monkeypatch.setenv("GBCMS_JUDGE_LONG_INS", "cis")
    c, tally = _run(tmp_path, caplog, germline=False)
    assert tally[("thr", "ALT")] == 4 and tally[("end", "ALT")] == 4, tally
    assert tally[("rend", "REF")] == 3, tally
    assert tally[("clip", "depth")] == 2 and tally[("ref", "REF")] == 8, tally
    assert (c.rd, c.ad) == (13, 8)


def test_the_cis_option_keeps_a_change_a_ref_molecule_carries_masked(tmp_path, caplog, monkeypatch):
    """Under the `cis` option, three REF molecules carrying the change at 589 too (a germline
    change on REF's haplotype) keep the mask on REF: both windows read on, and the reads ending
    at 589 stay depth only."""
    monkeypatch.setenv("GBCMS_JUDGE_LONG_INS", "cis")
    c, tally = _run(tmp_path, caplog, germline=True)
    assert tally[("end", "depth")] == 4 and tally[("end", "ALT")] == 0, tally
    assert tally[("rend", "depth")] == 3, tally
    assert tally[("thr", "ALT")] == 4 and tally[("germ", "REF")] == 3, tally
    assert (c.rd, c.ad) == (13, 4)


def test_the_shape_is_a_tandem_copy_with_one_telling_base():
    """Guard on the shape itself: the insert slides over 501-588 and 589 tells it apart."""
    contig, alt, _ = _shape(False)
    ins = alt[1:]
    assert contig[A + 1 : A + N] == ins[:-1] and contig[A + N] != ins[-1] and contig[A] != ins[-1]
    assert pysam.__version__  # pysam builds the BAM
