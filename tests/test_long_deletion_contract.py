"""Pure deletions of 50 bp or more are judged by their bases at the POS breakpoint
(RJ-24, C40 #254; docs/reference/read-judgment.md).

The engine's calls on every long-deletion read shape of the read-judgment spec
must equal the read census's, which judges the same windows from the bases alone
(clips and the hypothesis placement included), and must not depend on how the
variants are binned (`count_checked`). Two kinds of shape the census does not model are
left to the spec: changes recurring across reads (masked in the engine) and a
read's own 5+ bp indel inside the deleted span (C26, from its alignment).
A recurring change on the carriers alone names the larger allele they carry.
Synthetic and PHI-free.
"""

import glob

import pytest
import read_judgment_cases as rj
from census import assert_matches, census
from helpers import count_checked, make_read, read_maf_output, write_contig
from typer.testing import CliRunner

from gbcms import _rs
from gbcms.cli import app

A = rj.A
CENSUS_LIMITS = {
    "carriers with a recurring change two bases past the microhomology",
    "carriers with a recurring change in the flank",
    "right-clipped carriers with a recurring change in the flank",
    "carriers with a recurring change on the window's last base",
    "carriers with a run of changes from the window's last base on",
    "carriers with a run of three from the window's last base on",
    "carriers with a phased DNP inside the window",
    "carriers with a DNP at the junction, its first base the deleted one's",
    "REF reads with a recurring change inside the window",
    "REF reads with a run of changes from the window's last base on",
    "REF reads with a phased DNP inside the window",
    "REF reads with a run of changes in the flank",
    "carriers and REF reads with a run of changes in the flank",
    "REF read with its own 20-bp deletion inside the deleted span",
}
SHAPES = [
    (v, shape)
    for v, shapes in rj.LONG_DEL_SHAPES.items()
    for shape in shapes
    if shape not in CENSUS_LIMITS
]


def _prepared(tmp_path, v, shape):
    n, mh = rj.LONG_DEL[v]
    contig = rj._long_contig(n, mh)
    fa, bam = write_contig(tmp_path, contig, rj._long_del_reads(n, mh, shape), "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A - 1, contig[A - 1 : A + n], contig[A - 1], "DEL")],
        fa,
        5,
        False,
        1,
        True,
    )
    return contig, bam, pv.variant


@pytest.mark.parametrize("v, shape", SHAPES, ids=[f"{v} | {s}" for v, s in SHAPES])
def test_the_engine_matches_the_census(tmp_path, v, shape):
    contig, bam, variant = _prepared(tmp_path, v, shape)
    (c,) = count_checked(bam, [variant])
    assert_matches(c, census(bam, contig, variant), f"{v} | {shape}")


def _canonical(pos: int, ref: str, alt: str, contig: str) -> tuple[int, str, str]:
    """observed::canonical: trim shared last bases (extending left when an allele
    empties), then shared leading bases keeping one; 0-based position."""
    while True:
        if ref and alt and ref[-1] == alt[-1]:
            ref, alt = ref[:-1], alt[:-1]
        elif not ref or not alt:
            pos, ref, alt = pos - 1, contig[pos - 1] + ref, contig[pos - 1] + alt
        else:
            break
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        pos, ref, alt = pos + 1, ref[1:], alt[1:]
    return pos, ref, alt


@pytest.mark.parametrize(
    "shape, at",
    [
        ("carriers with a recurring change two bases past the microhomology", "past"),
        ("right-clipped carriers with a recurring change in the flank", "flank"),
        ("carriers with a run of changes from the window's last base on", "run"),
    ],
)
def test_a_recurring_change_on_the_carriers_names_the_larger_allele(tmp_path, shape, at):
    """Every carrier shows a change no REF read can (two bases past the junction,
    where REF reads read deleted bases; or in the flank, on carriers alone; or a run
    from the window's last base on that rejoins the flank). The change is masked, so
    the carriers count ALT, and the diagnostic names the deletion with the change:
    from the windows' first base through the change (a run's last base, past the
    window included), and at least through the junction."""
    n, mh = rj.LONG_DEL["del 60"]
    contig = rj._long_contig(n, mh)
    reads = rj._long_del_reads(n, mh, shape) + rj._long_del_reads(n, mh, "REF read, long flanks")
    for i, r in enumerate(reads):
        r.query_name = f"m{i}"
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        f"##fileformat=VCFv4.2\n##contig=<ID=1,length={len(contig)}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"1\t{A}\t.\t{contig[A - 1 : A + n]}\t{contig[A - 1]}\t.\t.\t.\n"
    )
    out = tmp_path / "out"
    res = CliRunner().invoke(
        app, ["dna", "-v", str(vcf), "-b", f"S:{bam}", "-f", fa, "-o", str(out), "--format", "maf"]
    )
    assert res.exit_code == 0, res.output
    (row,) = read_maf_output(glob.glob(str(out / "*.maf"))[0])
    assert (int(row["ref_count"]), int(row["alt_count"])) == (4, 4)
    first = A - 4  # the windows' first base, 0-based
    if at == "past":
        changed = A + n + mh + rj.RECUR_PAST
        ref_al = contig[first : changed + 1]
        alt_al = contig[first:A] + contig[A + n : changed] + rj._other(contig[changed])
    elif at == "run":
        changed = A + n + mh + 3  # the window's last base, and the base past it
        ref_al = contig[first : changed + 2]
        alt_al = (
            contig[first:A]
            + contig[A + n : changed]
            + rj._other(contig[changed])
            + rj._other(contig[changed + 1])
        )
    else:
        changed = A - 2  # through the junction: the first base past it
        ref_al = contig[first : A + n + 1]
        alt_al = (
            contig[first:changed]
            + rj._other(contig[changed])
            + contig[changed + 1 : A]
            + contig[A + n]
        )
    pos, ref_al, alt_al = _canonical(first, ref_al, alt_al, contig)
    assert len(ref_al) - len(alt_al) == n
    larger = f"OBSERVED_ALLELE(1:{pos + 1}:{ref_al}>{alt_al}:4/0)"
    assert larger in row["gbcms_diagnostic"].split(";"), row["gbcms_diagnostic"]


def test_a_germline_snp_at_the_anchor_keeps_clipped_carriers(tmp_path):
    """Every read, REF and ALT, shows a germline SNP at the base before the deletion,
    where the flank's other three bases equal the deleted bases before the right flank:
    the guard masks the SNP, and a carrier clipped at the junction (placed from the
    right flank) still shows the junction by that base against the deleted tail.
    D-written and clipped carriers count alike."""
    n, mh = rj.LONG_DEL["del 60"]
    c = list(rj._long_contig(n, mh))
    c[A - 4 : A - 1] = c[A + n - 4 : A + n - 1]
    contig = "".join(c)
    assert contig[A - 1] != contig[A + n - 1]  # still placed leftmost

    def snp(seq, s):
        return seq[: A - 1 - s] + "A" + seq[A - s :]  # the contig has no A

    hap = contig[:A] + contig[A + n :]
    reads = []
    for i in range(4):
        s = A - 50 - i
        reads.append(
            make_read(
                f"d{i}", snp(hap[s : s + 120], s), s, ((0, A - s), (2, n), (0, 120 - (A - s)))
            )
        )
        reads.append(make_read(f"r{i}", snp(contig[s : s + 120], s), s, ((0, 120),)))
        k = 20 + i  # bases before the junction, in the clip
        seq = snp(hap[A - k : A - k + 120], A - k)
        reads.append(make_read(f"c{i}", seq, A + n, ((4, k), (0, 120 - k))))
    for r in reads:
        r.flag = 0x1 | 0x2 | 0x20 | 0x40
        r.next_reference_id, r.next_reference_start, r.template_length = (
            0,
            r.reference_start + 250,
            400,
        )
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A - 1, contig[A - 1 : A + n], contig[A - 1], "DEL")],
        fa,
        5,
        False,
        1,
        True,
    )
    (got,) = count_checked(bam, [pv.variant])
    assert (got.rd, got.ad) == (4, 8)


@pytest.mark.parametrize("min_baseq", [0, 1])
def test_an_n_clip_shows_no_flank_at_any_minimum_quality(tmp_path, min_baseq):
    """REF reads starting inside the deleted span behind an N clip show nothing of the
    left flank: at --min-baseq 0 or 1 one base's evidence weighs nothing, and a tie
    places no read."""
    n, mh = rj.LONG_DEL["del 60"]
    contig = rj._long_contig(n, mh)
    reads = rj._long_del_reads(n, mh, "REF read starting inside the deletion behind an N clip")
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A - 1, contig[A - 1 : A + n], contig[A - 1], "DEL")],
        fa,
        5,
        False,
        1,
        True,
    )
    (got,) = count_checked(bam, [pv.variant], min_baseq=min_baseq)
    assert got.ad == 0


@pytest.mark.parametrize(
    "offset, both, rd",
    [
        (20, False, 4),  # deep in the deleted span
        (4, False, 4),  # the first deleted base past REF's window
        (3, False, 0),  # REF's window's last base
        (-3, True, 4),  # the left flank both windows read, on both alleles
    ],
)
def test_a_sibling_claims_ref_reads_only_where_ref_is_read(tmp_path, offset, both, rd):
    """Every REF read (and, in the flank, every carrier) carries a co-annotated SNV at
    `offset` from the first deleted base. Inside REF's window, past the shared flank,
    the SNV's carriers are its molecules there (RJ-5) and leave the deletion's REF.
    Further into the deleted bases the windows read nothing, for either allele, and
    no carrier holds them; in the flank both alleles read it: there the REF reads
    count REF, as the carriers count ALT."""
    n, mh = rj.LONG_DEL["del 60"]
    contig = rj._long_contig(n, mh)
    hap = contig[:A] + contig[A + n :]
    g = A + offset  # the SNV, 0-based; the contig has no A

    def snv(seq, s, at):
        return seq[: at - s] + "A" + seq[at - s + 1 :]

    reads = []
    for i in range(4):
        s = A - 50 - i
        carrier = hap[s : s + 100]
        if both:
            carrier = snv(carrier, s, g)
        reads.append(make_read(f"d{i}", carrier, s, ((0, A - s), (2, n), (0, 100 - (A - s)))))
        reads.append(make_read(f"r{i}", snv(contig[s : s + 100], s, g), s, ((0, 100),)))
    for r in reads:
        r.flag = 0x1 | 0x2 | 0x20 | 0x40
        r.next_reference_id, r.next_reference_start, r.template_length = (
            0,
            r.reference_start + 250,
            400,
        )
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    deletion, sibling = (
        pv.variant
        for pv in _rs.prepare_variants(
            [
                _rs.Variant("1", A - 1, contig[A - 1 : A + n], contig[A - 1], "DEL"),
                _rs.Variant("1", g, contig[g], "A", "SNP"),
            ],
            fa,
            5,
            False,
            1,
            True,
        )
    )
    got, _ = count_checked(bam, [deletion, sibling], sibling_variants=[[sibling], [deletion]])
    assert got.dp >= got.rd + got.ad and got.dpf >= got.rdf + got.adf
    assert got.rd == got.rd_fwd + got.rd_rev and got.ad == got.ad_fwd + got.ad_rev
    assert (got.rd, got.ad) == (rd, 4)


@pytest.mark.parametrize(
    "carriers, counts",
    [
        # a deletion one base longer at the same start: its bases differ from the given
        # deletion's at the window's last base and the base past it, then the A run
        # matches them whether shifted or not, but the deletion one base longer gives
        # the whole stretch: another breakpoint, partial
        ("one base longer", (4, 4, 4)),
        # the given deletion with a phased TT there: no other junction gives it, so it
        # rejoins the flank, a phased neighbour masked as a recurring change is
        ("a phased DNP", (4, 8, 0)),
    ],
)
def test_a_run_rejoins_the_flank_only_at_the_given_junction(tmp_path, carriers, counts):
    """The right flank is GGGG C, then a run of 20 A. Four exact carriers, four REF
    reads, and four carriers of either a deletion one base longer (GGG C A A...) or the
    given deletion with TT on the window's last base and the base past it. Both
    differ from the given ALT on those two bases alone and match the A run after
    them; only the phased TT is the given deletion with a neighbour (#285)."""
    n, mh = rj.LONG_DEL["del 60"]
    c = list(rj._long_contig(n, mh))
    c[A + n : A + n + 25] = "GGGGC" + "A" * 20
    c[A] = "T"  # the first deleted base differs from the flank's: no microhomology
    contig = "".join(c)
    assert contig[A - 1] != contig[A + n - 1]
    hap = contig[:A] + contig[A + n :]
    longer = contig[:A] + contig[A + n + 1 :]
    reads = []
    for i in range(4):
        s = A - 50 - i
        cig = ((0, A - s), (2, n), (0, 100 - (A - s)))
        reads.append(make_read(f"d{i}", hap[s : s + 100], s, cig))
        reads.append(make_read(f"r{i}", contig[s : s + 100], s, ((0, 100),)))
        if carriers == "one base longer":
            other = make_read(
                f"x{i}", longer[s : s + 100], s, ((0, A - s), (2, n + 1), (0, 100 - (A - s)))
            )
        else:
            seq = hap[s : s + 100]
            k = A + 3 - s  # the window's last base on the read
            other = make_read(f"x{i}", seq[:k] + "TT" + seq[k + 2 :], s, cig)
        reads.append(other)
    for r in reads:
        r.flag = 0x1 | 0x2 | 0x20 | 0x40
        r.next_reference_id, r.next_reference_start, r.template_length = (
            0,
            r.reference_start + 250,
            400,
        )
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A - 1, contig[A - 1 : A + n], contig[A - 1], "DEL")],
        fa,
        5,
        False,
        1,
        True,
    )
    (got,) = count_checked(bam, [pv.variant])
    assert (got.rd, got.ad, got.partial_alt) == counts


@pytest.mark.parametrize(
    "with_given, counts",
    [
        (False, (4, 0, 4)),
        (True, (4, 4, 4)),
    ],
)
def test_another_unit_count_of_a_vntr_is_another_breakpoint(tmp_path, with_given, counts):
    """Three 60-bp units: the second is the first with its third base changed, the
    third the second with its sixth and seventh changed. The given deletion removes one
    unit; carriers removing two show the third unit where the given ALT has the second:
    a DNP on the window's last base and the base past it, then the same unit on. The
    deletion of two units gives their bases whole, past the 60 bases held beside an
    event: another breakpoint, partial, beside carriers of the given one or alone."""
    n, mh = rj.LONG_DEL["del 60"]
    c = list(rj._long_contig(n, mh))
    u1 = c[A : A + n]
    u2 = list(u1)
    u2[2] = rj._other(u1[2])
    u3 = list(u2)
    u3[5], u3[6] = rj._other(u2[5]), rj._other(u2[6])
    c[A + n : A + 2 * n], c[A + 2 * n : A + 3 * n] = u2, u3
    contig = "".join(c)
    assert contig[A - 1] != contig[A + n - 1]
    one, two = contig[:A] + contig[A + n :], contig[:A] + contig[A + 2 * n :]
    reads = []
    for i in range(4):
        s = A - 50 - i
        reads.append(make_read(f"r{i}", contig[s : s + 100], s, ((0, 100),)))
        reads.append(
            make_read(f"x{i}", two[s : s + 100], s, ((0, A - s), (2, 2 * n), (0, 100 - (A - s))))
        )
        if with_given:
            reads.append(
                make_read(f"d{i}", one[s : s + 100], s, ((0, A - s), (2, n), (0, 100 - (A - s))))
            )
    for r in reads:
        r.flag = 0x1 | 0x2 | 0x20 | 0x40
        r.next_reference_id, r.next_reference_start, r.template_length = (
            0,
            r.reference_start + 250,
            400,
        )
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A - 1, contig[A - 1 : A + n], contig[A - 1], "DEL")],
        fa,
        5,
        False,
        1,
        True,
    )
    (got,) = count_checked(bam, [pv.variant])
    assert (got.rd, got.ad, got.partial_alt) == counts


def test_a_fragment_reads_on_past_a_run_from_either_mate(tmp_path):
    """Overlapping mates of four carriers with a phased DNP on the window's last base and
    the base past it: the left mate ends 14 bases past the junction, inside the bases
    the run must rejoin over; the right mate reads well past them. The fragment's bases
    are read from whichever mate reads them, so the run rejoins and every mate counts
    ALT, as the right mates alone do."""
    n, mh = rj.LONG_DEL["del 60"]
    contig = rj._long_contig(n, mh)
    hap = list(contig[:A] + contig[A + n :])
    hap[A + 3], hap[A + 4] = rj._other(hap[A + 3]), rj._other(hap[A + 4])
    hap = "".join(hap)
    reads = []
    for i in range(4):
        s = A - 50 - i
        reads.append(make_read(f"r{i}", contig[s : s + 100], s, ((0, 100),)))
        s1, s2 = A + 14 - 100, A - 10 - i
        left = make_read(
            f"f{i}", hap[s1 : s1 + 100], s1, ((0, A - s1), (2, n), (0, 100 - (A - s1)))
        )
        right = make_read(
            f"f{i}", hap[s2 : s2 + 100], s2, ((0, A - s2), (2, n), (0, 100 - (A - s2)))
        )
        left.flag, right.flag = 0x1 | 0x2 | 0x20 | 0x40, 0x1 | 0x2 | 0x10 | 0x80
        size = s2 + 100 + n - s1
        left.next_reference_id, left.next_reference_start, left.template_length = 0, s2, size
        right.next_reference_id, right.next_reference_start, right.template_length = 0, s1, -size
        reads += [left, right]
    for r in reads[0::3]:
        r.flag = 0x1 | 0x2 | 0x20 | 0x40
        r.next_reference_id, r.next_reference_start, r.template_length = (
            0,
            r.reference_start + 250,
            400,
        )
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A - 1, contig[A - 1 : A + n], contig[A - 1], "DEL")],
        fa,
        5,
        False,
        1,
        True,
    )
    (got,) = count_checked(bam, [pv.variant])
    assert (got.rd, got.ad, got.partial_alt) == (4, 8, 0)


def test_another_allele_reading_like_alt_past_the_windows_is_not_alt(tmp_path):
    """REF's and ALT's continuations agree past their first two bases. Four reads of
    another allele hold ALT's first continuation base, REF's second, a change two bases
    before the junction, ALT's twelve bases past the windows, then one base deleted: their
    four continuation bases miss each allele by one, their twelve bases past the windows
    read as ALT's, and they go on to differ from it. They are not placed as carriers, so
    their changes are not masked as neighbours: neither allele, beside four carriers and
    four REF reads."""
    n, mh = rj.LONG_DEL["del 60"]
    c = list(rj._long_contig(n, mh))
    c[A + 2], c[A + 3] = c[A + n + 2], c[A + n + 3]
    if c[A + 1] == c[A + n + 1]:
        c[A + 1] = next(b for b in "CGT" if b != c[A + n + 1])
    contig = "".join(c)
    assert contig[A] != contig[A + n] and contig[A - 1] != contig[A + n - 1]
    hap = contig[:A] + contig[A + n :]
    other = list(hap)
    other[A - 2] = rj._other(other[A - 2])
    other[A + 1] = contig[A + 1]
    del other[A + 4 + 12 + 3]  # one base deleted three past the twelve
    other = "".join(other)
    reads = []
    for i in range(4):
        s = A - 50 - i
        cig = ((0, A - s), (2, n), (0, 100 - (A - s)))
        reads.append(make_read(f"d{i}", hap[s : s + 100], s, cig))
        reads.append(make_read(f"x{i}", other[s : s + 100], s, cig))  # aligned as written
        reads.append(make_read(f"r{i}", contig[s : s + 100], s, ((0, 100),)))
    for r in reads:
        r.flag = 0x1 | 0x2 | 0x20 | 0x40
        r.next_reference_id, r.next_reference_start, r.template_length = (
            0,
            r.reference_start + 250,
            400,
        )
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", A - 1, contig[A - 1 : A + n], contig[A - 1], "DEL")],
        fa,
        5,
        False,
        1,
        True,
    )
    (got,) = count_checked(bam, [pv.variant])
    assert (got.rd, got.ad) == (4, 4)
