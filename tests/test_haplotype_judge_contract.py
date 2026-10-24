"""Package 2 (#304): the haplotype judge, `--read-judge haplotype`.

For a row X and each admitted read, the read is scored whole against the local haplotypes:
every context of the other listed rows within read reach, crossed with every state of X's
site (REF, X, each other allele overlapping X). It counts ALT where the best haplotype
carrying X beats the best one without it by one confident base, REF where the best with
REF at the site beats every other site state by as much; otherwise it is depth only (or
partial evidence where another site allele wins). Because every context appears in every
state, only the read's bases at X decide: a linked variant elsewhere never credits it.

Each shape is one the size-specific rules got wrong or that the decision record (#304)
names; counts are (ref_count, alt_count) per row in input order. Runs through the CLI.
"""

import glob
import random

import pysam
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()
READ_LEN = 100
COUNT_COLUMNS = (
    "total_count",
    "ref_count",
    "alt_count",
    "total_count_fragment",
    "ref_count_fragment",
    "alt_count_fragment",
    "ref_count_forward",
    "ref_count_reverse",
    "alt_count_forward",
    "alt_count_reverse",
)


def _contig(seed, length=1200, motifs=()):
    """A random contig with `motifs` ((0-based start, sequence), ...) written in."""
    rng = random.Random(seed)
    c = [rng.choice("ACGT") for _ in range(length)]
    for start, seq in motifs:
        c[start : start + len(seq)] = seq
    return "".join(c)


def _apply(ref, rows):
    """`ref` with VCF rows ((1-based pos, REF, ALT), ...) applied, right to left."""
    h = ref
    for pos, r, a in sorted(rows, reverse=True):
        assert h[pos - 1 : pos - 1 + len(r)] == r, (pos, r)
        h = h[: pos - 1] + a + h[pos - 1 + len(r) :]
    return h


def _plain(name, hap, start, quals=None, length=READ_LEN):
    """A read of `hap` from `start`, written gapless at `start` (on REF, or through an
    event as an aligner without the gap would)."""
    seq = hap[start : start + length]
    return make_read(name, seq, start, ((0, len(seq)),), quals=quals)


def _gapped(name, ref, pos, alt_ins, start, length=READ_LEN):
    """A carrier of an insertion after the 1-based anchor `pos`, its gap written there."""
    left = pos - start
    seq = ref[start:pos] + alt_ins + ref[pos : start + length - len(alt_ins)]
    return make_read(
        name, seq, start, ((0, left), (1, len(alt_ins)), (0, length - left - len(alt_ins)))
    )


def _run(tmp_path, ref, rows, reads):
    fa = tmp_path / "ref.fasta"
    fa.write_text(">1\n" + ref + "\n")
    pysam.faidx(str(fa))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": len(ref)}]}
    raw, bam = tmp_path / "raw.bam", tmp_path / "reads.bam"
    with pysam.AlignmentFile(raw, "wb", header=header) as fh:
        for a in reads:
            fh.write(a)
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    vcf = tmp_path / "rows.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(f"1\t{p}\t.\t{r}\t{a}\t.\t.\t.\n" for p, r, a in rows)
    )
    out = tmp_path / "out"
    result = runner.invoke(
        app,
        [
            "dna",
            "-v",
            str(vcf),
            "-b",
            f"S:{bam}",
            "-f",
            str(fa),
            "-o",
            str(out),
            "--format",
            "maf",
            "--read-judge",
            "haplotype",
        ],
        env={"COLUMNS": "400"},
    )
    assert result.exit_code == 0, result.output
    got = list(read_maf_output(glob.glob(str(out / "*.maf"))[0]))
    for r in got:
        n = {k: int(r[k]) for k in COUNT_COLUMNS}
        assert n["total_count"] >= n["ref_count"] + n["alt_count"]
        assert n["total_count_fragment"] >= n["ref_count_fragment"] + n["alt_count_fragment"]
        assert n["ref_count"] == n["ref_count_forward"] + n["ref_count_reverse"]
        assert n["alt_count"] == n["alt_count_forward"] + n["alt_count_reverse"]
    return [(int(r["ref_count"]), int(r["alt_count"]), int(r["total_count"])) for r in got]


# ── An SNV one base after an insertion's anchor, its ALT the first inserted base (#246) ──

T0 = 300  # 1-based anchor C at 0-based 299; the SNV is the A after it
REF246 = _contig(246, motifs=((T0 - 1, "CA"),))
INS246 = (T0, "C", "CTGGAC")
SNV246 = (T0 + 1, "A", "T")


def test_an_insertion_carrier_written_through_counts_for_the_insertion_only(tmp_path):
    """In trans: 4 insertion carriers with the gap at the anchor, 6 written through (their
    first inserted T on the SNV's base), 5 SNV carriers, 4 REF reads. The insertion row
    counts its 10 carriers, REF for the SNV and REF reads; the SNV row counts its 5."""
    ins_hap, snv_hap = _apply(REF246, [INS246]), _apply(REF246, [SNV246])
    reads = [_gapped(f"j{i}", REF246, T0, "TGGAC", T0 - 40 - i) for i in range(4)]
    reads += [_plain(f"t{i}", ins_hap, T0 - 40 - i) for i in range(6)]
    reads += [_plain(f"s{i}", snv_hap, T0 - 50 + i) for i in range(5)]
    reads += [_plain(f"r{i}", REF246, T0 - 50 + i) for i in range(4)]
    (ins_rd, ins_ad, _), (snv_rd, snv_ad, _) = _run(tmp_path, REF246, [INS246, SNV246], reads)
    assert (ins_rd, ins_ad) == (9, 10)
    assert (snv_rd, snv_ad) == (14, 5)


def test_a_molecule_carrying_both_in_cis_counts_at_both_rows(tmp_path):
    """In cis: 6 molecules carry the insertion and, after it, the SNV's T; 4 REF reads."""
    both = _apply(REF246, [INS246, SNV246])
    reads = [_plain(f"b{i}", both, T0 - 40 - i) for i in range(6)]
    reads += [_plain(f"r{i}", REF246, T0 - 50 + i) for i in range(4)]
    (ins_rd, ins_ad, _), (snv_rd, snv_ad, _) = _run(tmp_path, REF246, [INS246, SNV246], reads)
    assert (ins_rd, ins_ad, snv_rd, snv_ad) == (4, 6, 4, 6)


# ── Two alleles at one site ──


def test_two_insertions_at_one_anchor_count_their_own_carriers(tmp_path):
    """A compound het: +TGGAC and +AATTG after one anchor, 5 carriers each, 3 REF reads.
    Each row counts its own carriers; the other's are partial evidence, never REF."""
    a, b = (T0, "C", "CTGGAC"), (T0, "C", "CAATTG")
    reads = [_gapped(f"a{i}", REF246, T0, "TGGAC", T0 - 40 - i) for i in range(5)]
    reads += [_gapped(f"b{i}", REF246, T0, "AATTG", T0 - 40 - i) for i in range(5)]
    reads += [_plain(f"r{i}", REF246, T0 - 50 + i) for i in range(3)]
    (a_rd, a_ad, _), (b_rd, b_ad, _) = _run(tmp_path, REF246, [a, b], reads)
    assert (a_rd, a_ad, b_rd, b_ad) == (3, 5, 3, 5)


# ── A repeat: reads that cannot tell the unit count ──

R0 = 400  # 1-based anchor G, then (CA)x8, then a unique flank
TRACT = "CA" * 8
REF_REP = _contig(290, motifs=((R0 - 1, "G" + TRACT + "TTGACG"),))
INS_REP = (R0, "G", "GCA")


def test_reads_that_cannot_tell_the_unit_count_are_depth_only(tmp_path):
    """+CA in a (CA)x8 tract: 5 reads spanning it with 9 units, 4 with 8, 3 ending inside
    the tract (either length), 2 spanning it with 10 units (another allele, unlisted)."""
    alt = _apply(REF_REP, [INS_REP])
    other = _apply(REF_REP, [(R0, "G", "GCACA")])
    reads = [_plain(f"a{i}", alt, R0 - 60 + i) for i in range(5)]
    reads += [_plain(f"r{i}", REF_REP, R0 - 60 + i) for i in range(4)]
    reads += [_plain(f"x{i}", REF_REP, R0 + 6 - READ_LEN + i) for i in range(3)]  # end in the tract
    reads += [_plain(f"o{i}", other, R0 - 60 + i) for i in range(2)]
    [(rd, ad, dp)] = _run(tmp_path, REF_REP, [INS_REP], reads)
    assert (rd, ad) == (4, 5)
    assert dp == 14


# ── Base quality decides nothing below the gate ──

S0 = 500
REF_SNV = _contig(500)
SNV_Q = (S0, REF_SNV[S0 - 1], "G" if REF_SNV[S0 - 1] != "G" else "T")


def test_a_low_quality_base_at_the_only_telling_base_is_depth_only(tmp_path):
    """5 carriers read the SNV at Q30, 4 at Q2, 3 REF reads."""
    alt = _apply(REF_SNV, [SNV_Q])
    reads = [_plain(f"a{i}", alt, S0 - 50 + i) for i in range(5)]
    for i in range(4):
        q = [30] * READ_LEN
        q[S0 - 1 - (S0 - 50 + i)] = 2
        reads.append(_plain(f"q{i}", alt, S0 - 50 + i, quals=q))
    reads += [_plain(f"r{i}", REF_SNV, S0 - 50 + i) for i in range(3)]
    [(rd, ad, dp)] = _run(tmp_path, REF_SNV, [SNV_Q], reads)
    assert (rd, ad, dp) == (3, 5, 12)


def test_a_linked_variant_never_credits_a_read_whose_bases_at_x_say_nothing(tmp_path):
    """Two SNVs 20 bp apart in cis: 5 molecules read both; 4 read the second but have N at
    the first; 3 REF reads. The first row counts 5, the second 9: phase is no evidence."""
    y = (S0 + 20, REF_SNV[S0 + 19], "C" if REF_SNV[S0 + 19] != "C" else "A")
    both = _apply(REF_SNV, [SNV_Q, y])
    reads = [_plain(f"b{i}", both, S0 - 50 + i) for i in range(5)]
    for i in range(4):
        start = S0 - 50 + i
        seq = list(both[start : start + READ_LEN])
        seq[S0 - 1 - start] = "N"
        reads.append(make_read(f"n{i}", "".join(seq), start, ((0, READ_LEN),)))
    reads += [_plain(f"r{i}", REF_SNV, S0 - 50 + i) for i in range(3)]
    (x_rd, x_ad, _), (y_rd, y_ad, _) = _run(tmp_path, REF_SNV, [SNV_Q, y], reads)
    assert (x_rd, x_ad) == (3, 5)
    assert (y_rd, y_ad) == (3, 9)


# ── Long events ──

L0 = 600
REF_LONG = _contig(144, length=2400)
INS60 = "".join(random.Random(60).choice("ACGT") for _ in range(60))


def test_a_long_insertion_carrier_written_as_a_soft_clip_counts(tmp_path):
    """A 60-bp insertion: 4 carriers with the gap at the anchor, 4 aligned up to the anchor
    with the insert and the right flank soft-clipped, 4 REF reads."""
    row = (L0, REF_LONG[L0 - 1], REF_LONG[L0 - 1] + INS60)
    reads = [_gapped(f"g{i}", REF_LONG, L0, INS60, L0 - 30 - i) for i in range(4)]
    for i in range(4):
        start = L0 - 30 - i
        seq = REF_LONG[start:L0] + INS60 + REF_LONG[L0 : start + READ_LEN - 60]
        left = L0 - start
        reads.append(make_read(f"c{i}", seq, start, ((0, left), (4, READ_LEN - left))))
    reads += [_plain(f"r{i}", REF_LONG, L0 - 50 + i) for i in range(4)]
    [(rd, ad, _)] = _run(tmp_path, REF_LONG, [row], reads)
    assert (rd, ad) == (4, 8)


def test_a_long_deletion_counts_ref_on_reads_over_its_anchor_only(tmp_path):
    """A 400-bp deletion: 5 carriers read across its junction, 4 REF reads over the anchor,
    6 reads wholly inside the deleted bases (REF's molecules, but no read of the ALT could
    lie there: they are not admitted, as the size-specific rules did not admit them)."""
    row = (L0, REF_LONG[L0 - 1 : L0 + 400], REF_LONG[L0 - 1])
    alt = _apply(REF_LONG, [row])
    reads = [_plain(f"a{i}", alt, L0 - 50 + i) for i in range(5)]
    reads += [_plain(f"r{i}", REF_LONG, L0 - 50 + i) for i in range(4)]
    reads += [_plain(f"x{i}", REF_LONG, L0 + 100 + 20 * i) for i in range(6)]
    [(rd, ad, dp)] = _run(tmp_path, REF_LONG, [row], reads)
    assert (rd, ad, dp) == (4, 5, 9)


# ── Prep: the rows a read is weighed against come from the whole input ──


def test_prep_lists_every_other_row_within_a_reads_reach(tmp_path):
    """A row lists the input's other rows whose REF span lies within 300 bases of its own
    (judge::READ_REACH), whatever the bins; never itself, a farther row or another
    chromosome's."""
    from helpers import write_contig

    from gbcms import _rs

    contig = _contig(77, length=3000)
    fa, _ = write_contig(tmp_path, contig, [], "c")
    x = 1000

    def snv(pos):
        return _rs.Variant("1", pos, contig[pos], "A" if contig[pos] != "A" else "C", "SNP")

    deletion = _rs.Variant("1", x - 305, contig[x - 305 : x - 295], contig[x - 305], "DEL")
    rows = [snv(x), snv(x + 100), snv(x + 299), snv(x + 302), snv(x - 299), deletion, snv(x + 900)]
    prepared = [pv.variant for pv in _rs.prepare_variants(rows, fa, 5, False, 1, True)]
    near = {(p, r, a) for p, r, a in prepared[0].reach_rows}
    want = {
        (v.pos, v.ref_allele.upper(), v.alt_allele.upper())
        for v in (prepared[1], prepared[2], prepared[4], prepared[5])
    }
    assert near == want, "within reach: +100, +299, -299 and the deletion ending 295 bases before X"
    assert all(
        (p, r, a) != (prepared[0].pos, prepared[0].ref_allele, prepared[0].alt_allele)
        for p, r, a in prepared[0].reach_rows
    )
    assert (
        prepared[6].reach_rows == []
    ), "a row 900 bases away lists nothing it cannot share a read with"


# ── Admission: a read that misses the anchor counts only where its alignment to the
# haplotype it is judged on holds the anchor ──

A0 = 700
REF_ADM = _contig(700, length=1600)
INS30 = "".join(random.Random(30).choice("ACGT") for _ in range(30))


def test_reads_past_a_deletions_anchor_are_not_admitted(tmp_path):
    """A 30-bp deletion: 5 carriers over the junction, 4 REF reads over the anchor; 6 REF
    reads starting inside the deleted bases and 6 carriers aligned after the junction hold
    the anchor on neither haplotype, so they are not depth (as the size-specific rules)."""
    row = (A0, REF_ADM[A0 - 1 : A0 + 30], REF_ADM[A0 - 1])
    alt = _apply(REF_ADM, [row])
    reads = [_plain(f"a{i}", alt, A0 - 50 + i) for i in range(5)]
    reads += [_plain(f"r{i}", REF_ADM, A0 - 50 + i) for i in range(4)]
    reads += [_plain(f"x{i}", REF_ADM, A0 + 5 + i) for i in range(6)]
    for i in range(6):
        hs = A0 + 5 + i  # haplotype offset past the junction; on the genome 30 further
        reads.append(make_read(f"p{i}", alt[hs : hs + READ_LEN], hs + 30, ((0, READ_LEN),)))
    [(rd, ad, dp)] = _run(tmp_path, REF_ADM, [row], reads)
    assert (rd, ad, dp) == (4, 5, 9)


def test_reads_past_an_insertions_anchor_are_not_admitted(tmp_path):
    """A 30-bp insertion: 4 carriers with the gap at the anchor, 4 REF reads over it; 5
    carriers aligned from the base after the anchor with inserted bases soft-clipped before
    it, and 5 REF reads starting there, hold the anchor on neither haplotype."""
    row = (A0, REF_ADM[A0 - 1], REF_ADM[A0 - 1] + INS30)
    alt = _apply(REF_ADM, [row])
    reads = [_gapped(f"g{i}", REF_ADM, A0, INS30, A0 - 30 - i) for i in range(4)]
    reads += [_plain(f"r{i}", REF_ADM, A0 - 50 + i) for i in range(4)]
    for i in range(5):
        k = 10 + i  # inserted bases in the clip
        hs = A0 + 30 - k
        reads.append(make_read(f"c{i}", alt[hs : hs + READ_LEN], A0, ((4, k), (0, READ_LEN - k))))
    reads += [_plain(f"x{i}", REF_ADM, A0) for i in range(5)]
    [(rd, ad, dp)] = _run(tmp_path, REF_ADM, [row], reads)
    assert (rd, ad, dp) == (4, 4, 8)


# ── An MNP's ALT needs every changed base read ──


def test_an_mnp_read_showing_one_of_its_changes_is_not_alt(tmp_path):
    """A 3-base MNP: 10 carriers over it, 5 carriers ending on its first base (depth), 5
    starting on its last (not over the anchor), 10 REF reads."""
    m0 = 1000
    ref = _contig(5, length=2000)
    alt3 = "".join("A" if b != "A" else "C" for b in ref[m0 - 1 : m0 + 2])
    row = (m0, ref[m0 - 1 : m0 + 2], alt3)
    alt = _apply(ref, [row])
    reads = [_plain(f"f{i}", alt, m0 - 60 + i) for i in range(10)]
    reads += [_plain(f"e{i}", alt, m0 - READ_LEN) for i in range(5)]
    reads += [_plain(f"s{i}", alt, m0 + 1) for i in range(5)]
    reads += [_plain(f"r{i}", ref, m0 - 30 + i) for i in range(10)]
    [(rd, ad, dp)] = _run(tmp_path, ref, [row], reads)
    assert (rd, ad, dp) == (10, 10, 25)


# ── An SNV in a repeat beside a listed expansion of it ──


def test_an_snv_in_a_repeat_counts_its_carriers_whichever_side_the_expansion_lies(tmp_path):
    """(CA)x8 with an SNV in unit 5 and a listed +CA left-aligned at the anchor: 10
    carriers of the SNV with the expansion after it, 10 with it before. Each read's bases
    say G at the SNV."""
    p0 = 1001  # 1-based anchor G, (CA)x8 after it
    ref = _contig(7, length=2400, motifs=((p0 - 1, "G" + "CA" * 8 + "T"),))
    snv = (p0 + 10, "A", "G")  # unit 5's A
    ins = (p0, "G", "GCA")
    right = _apply(ref, [snv, (p0 + 16, "A", "ACA")])
    left = _apply(ref, [snv, ins])
    reads = []
    for i in range(10):
        s = p0 - 61 + i
        k = p0 + 17 - 1 - s
        reads.append(make_read(f"r{i}", right[s : s + 120], s, ((0, k), (1, 2), (0, 120 - k - 2))))
        k = p0 - s
        reads.append(make_read(f"l{i}", left[s : s + 120], s, ((0, k), (1, 2), (0, 120 - k - 2))))
    (snv_rd, snv_ad, snv_dp), _ = _run(tmp_path, ref, [snv, ins], reads)
    assert (snv_rd, snv_ad, snv_dp) == (0, 20, 20)


# ── #246 among many rows: the insertion's haplotype must be one the judge builds ──


def test_an_snv_on_an_insertions_first_base_among_many_rows_counts_its_own_carriers(tmp_path):
    """#246 in a CTTT repeat as HG002 carries it (chr6:36681273/4): a 12-bp insertion at
    the run's start whose first base is the SNV's ALT, the two in trans, each haplotype
    with its own rows around them (shared, and its own) within a read's reach, so the
    contexts are capped. The insertion's carriers, written 151M through the repeat, are
    REF for the SNV: the haplotype carrying the insertion with its neighbours is among
    those kept, so no state stands in for a missing row."""
    a0 = 700  # 1-based anchor T, the SNV the T after it
    tract = "GCTTTTCTTTC" + "T" + "TTCTTTTTTTC" + "TTTC" * 12 + "TTTTCTTTCTTTCTTTT"
    ref = _contig(36681273, length=1600, motifs=((a0 - 12, tract),))
    alt_of = lambda b: "C" if b != "C" else "A"  # noqa: E731
    snv_at = lambda p: (p, ref[p - 1], alt_of(ref[p - 1]))  # noqa: E731
    ins = (a0, "T", "TCTCTTTCTTTCC")
    snv = (a0 + 1, "T", "C")
    shared = [snv_at(a0 - 130), snv_at(a0 - 76), snv_at(a0 + 7), snv_at(a0 + 135)]
    h0_only = [snv_at(a0 - 40), snv_at(a0 + 61), snv_at(a0 + 138)]
    h1_only = [snv_at(a0 - 33), snv_at(a0 + 66), snv_at(a0 + 134)]
    rows = [ins, snv] + shared + h0_only + h1_only
    h0 = _apply(ref, [ins] + shared + h0_only)
    h1 = _apply(ref, [snv] + shared + h1_only)
    reads = []
    for i in range(12):
        s = a0 - 120 + 9 * i
        reads.append(make_read(f"i{i}", h0[s : s + 151], s, ((0, 151),)))
        reads.append(make_read(f"s{i}", h1[s : s + 151], s, ((0, 151),)))
    got = _run(tmp_path, ref, rows, reads)
    snv_rd, snv_ad, _ = got[1]
    assert (snv_rd, snv_ad) == (12, 12), got
