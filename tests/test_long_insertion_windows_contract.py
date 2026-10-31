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
alleles alike). Reading that base on REF's haplotypes as the reference's, where REF's
molecules lack the change, would credit such a carrier on other fragments' evidence (the cis
rescue, vetoed as REJ-20261023-001); no option does.

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
    r"read call \S+ \S+>(\S+) read=([a-z]+)\d+ mate=(\d) ref=(\w+) alt=(\w+) .*?partial=(\d+)"
)
A = 500  # the anchor (0-based)
N = 89  # the insert


def _shape(sibling: bool = False):
    """The contig, the row's ALT, and the reads: REF molecules over the anchor and past 589
    (`ref`), REF reads elsewhere (`refb`), carriers written with the insert (`thr`), carriers
    over the anchor ending at 589 written straight (`end`), REF reads over the anchor ending at
    589 (`rend`), the second copy's reads placed
    apart with the change at 589 (`two`), and fragments starting inside the insert with its
    tail clipped, placed after the anchor (`clip`, read 1 of a proper pair). With `sibling`,
    carriers of a second allele of the same length that differs from the first at the insert's
    last base alone (`sib`, its ALT is `_sibling_alt`), and carriers of either allele with a
    third base there (`thx`)."""
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
    if sibling:
        other, third = _sibling_bases(contig, last)
        for kind, base in (("sib", other), ("thx", third)):
            for i in range(4):
                s = 455 + 5 * i
                m1 = A + 1 - s
                seq = list(dup[s : s + 151])
                seq[m1 + N - 1] = base  # the insert's last base, as written in the read
                add(f"{kind}{i}", "".join(seq), s, ((0, m1), (1, N), (0, 151 - m1 - N)))
    return contig, contig[A] + ins, reads


def _sibling_bases(contig, last):
    """The sibling's base at the insert's last place, and the third base there: neither is the
    first allele's `last`, REF's base at 589, or the other."""
    other = next(b for b in "ACGT" if b not in (last, contig[A + N]))
    third = next(b for b in "ACGT" if b not in (last, contig[A + N], other))
    return other, third


def _sibling_alt(contig, alt):
    """The second allele's ALT: the first's insert with `_sibling_bases`' other base last."""
    return alt[:-1] + _sibling_bases(contig, alt[-1])[0]


def _run(tmp_path, caplog, sibling: bool = False):
    """Counts the row (and with `sibling` its second allele, listed in the same call) under the
    sequence judge, asserting the counting invariants on every row. Returns the first row's
    counts and the tally of read calls by (kind of read, call) per row ALT."""
    contig, alt, reads = _shape(sibling)
    fa, bam = write_contig(tmp_path, contig, reads, "dup")
    alts = [alt, _sibling_alt(contig, alt)] if sibling else [alt]
    pvs = _rs.prepare_variants(
        [_rs.Variant("1", A, contig[A], a, "INDEL") for a in alts],
        fa,
        5,
        False,
        1,
        True,
        read_judge="sequence",
    )
    caplog.clear()
    with caplog.at_level(5, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        counts = count_bam_checked(
            bam,
            [pv.variant for pv in pvs],
            [None] * len(pvs),
            *ARGS,
            read_judge="sequence",
            reference_fasta=fa,
        )
    for c in counts:
        assert c.dp >= c.rd + c.ad and c.dpf >= c.rdf + c.adf
        assert c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    calls = {}
    for record in caplog.records:
        if m := CALL.search(record.getMessage()):
            row, kind, mate, ref, alt_, partial = m.groups()
            call = (
                "ALT"
                if alt_ == "true"
                else "REF" if ref == "true" else "partial" if partial != "0" else "depth"
            )
            calls[(row, record.getMessage().split(" read=")[1].split()[0], mate)] = (kind, call)
    tallies = {a: collections.Counter() for a in alts}
    for (row, _, _), kind_call in calls.items():
        tallies[row][kind_call] += 1
    if sibling:
        return counts, tallies
    return counts[0], tallies[alt]


def test_a_long_insertions_reads_are_judged_over_equal_windows_from_its_anchor(tmp_path, caplog):
    """Carriers written with the insert are ALT; fragments starting inside the insert are
    depth only, though their pair admits them; reads ending at the masked base past the run,
    carriers and REF reads alike, are short of the windows read on and are depth only; REF
    molecules over the anchor and past 589 are REF."""
    c, tally = _run(tmp_path, caplog)
    assert tally[("thr", "ALT")] == 4, tally
    assert tally[("clip", "depth")] == 2 and tally[("clip", "ALT")] == 0, tally
    assert tally[("end", "depth")] == 4 and tally[("rend", "depth")] == 3, tally
    assert tally[("ref", "REF")] == 8, tally
    assert not any(call == "ALT" for (kind, call) in tally if kind in ("ref", "refb", "two")), tally
    assert (c.rd, c.ad) == (10, 4)


def test_a_sibling_allele_takes_the_masks_laid_on_both_windows_too(tmp_path, caplog):
    """A second allele of the insertion's length, listed with it and differing from it at the
    insert's last base alone, which the guard masks at 589 (89 bases from the anchor): it takes
    the mask at that offset as REF and the first allele do (RJ-27: the masks are read on every
    haplotype alike). No base tells the two alleles apart, so a carrier of either, or a read
    with a third base there, is partial evidence on both rows and ALT on neither (before, a
    mask on the sibling's haplotype was missing: each row counted its own carriers ALT, and
    the reads with a third base, no allele's, ALT on both). REF reads keep their calls."""
    counts, tallies = _run(tmp_path, caplog, sibling=True)
    for tally in tallies.values():
        for kind in ("thr", "sib", "thx"):
            assert tally[(kind, "ALT")] == 0 and tally[(kind, "partial")] == 4, (kind, tally)
        assert tally[("ref", "REF")] == 8, tally
    assert [(c.rd, c.ad, c.partial_alt) for c in counts] == [(10, 0, 12)] * 2


def test_the_shape_is_a_tandem_copy_with_one_telling_base():
    """Guard on the shape itself: the insert slides over 501-588 and 589 tells it apart."""
    contig, alt, _ = _shape()
    ins = alt[1:]
    assert contig[A + 1 : A + N] == ins[:-1] and contig[A + N] != ins[-1] and contig[A] != ins[-1]
    assert pysam.__version__  # pysam builds the BAM


def test_a_long_insertion_no_base_tells_apart_is_depth_only_and_named(tmp_path):
    """A 60 bp insertion of ten units of a 6-base repeat that runs on past the judge's reach:
    no base within the haplotypes the judge reads tells REF's window from the ALT's, so every
    read is depth only, and the row names it, NO_TELLING_BASE (no column added)."""
    import glob

    from helpers import read_maf_output
    from typer.testing import CliRunner

    from gbcms.cli import app

    rng = random.Random(2643)
    unit = "ACGTTG"
    left = "".join(rng.choice("ACGT") for _ in range(499)) + "C"
    contig = left + unit * 200 + "".join(rng.choice("ACGT") for _ in range(400))
    a = 499  # the anchor (0-based), before the repeat
    ins = unit * 10
    hx = contig[: a + 1] + ins + contig[a + 1 :]
    reads = []
    for i in range(6):
        s = 440 + 7 * i
        reads.append(make_read(f"ref{i}", contig[s : s + 150], s, ((0, 150),)))
        reads.append(
            make_read(f"alt{i}", hx[s : s + 150], s, ((0, a + 1 - s), (1, 60), (0, 89 - (a - s))))
        )
    fa, bam = write_contig(tmp_path, contig, reads, "rep")
    vcf = tmp_path / "rows.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"1\t{a + 1}\t.\t{contig[a]}\t{contig[a] + ins}\t.\t.\t.\n"
    )
    out = tmp_path / "out"
    args = ["dna", "-v", str(vcf), "-b", f"S:{bam}", "-f", fa, "-o", str(out), "--format", "maf"]
    result = CliRunner().invoke(app, args + ["--read-judge", "sequence"], env={"COLUMNS": "400"})
    assert result.exit_code == 0, result.output
    (row,) = list(read_maf_output(glob.glob(str(out / "*.maf"))[0]))
    rd, ad, dp = int(row["ref_count"]), int(row["alt_count"]), int(row["total_count"])
    assert dp >= rd + ad and (rd, ad) == (0, 0) and dp == 12, row
    assert int(row["ref_count_forward"]) + int(row["ref_count_reverse"]) == rd
    assert "NO_TELLING_BASE" in row["gbcms_diagnostic"].split(";"), row["gbcms_diagnostic"]
