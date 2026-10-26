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
import pytest
from helpers import count_checked, make_read, read_maf_output, write_contig
from typer.testing import CliRunner

from gbcms import _rs
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


def _run(tmp_path, ref, rows, reads, partial=False, maf_rows=False):
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
    if maf_rows:
        return got
    if partial:
        return [
            (int(r["ref_count"]), int(r["alt_count"]), int(r["total_count"]), int(r["partial_alt"]))
            for r in got
        ]
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
_r60 = random.Random(60)  # one generator: a fresh Random(60) per base gives a G run
INS60 = "".join(_r60.choice("ACGT") for _ in range(60))


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


def test_a_long_deletions_carriers_ending_inside_it_are_alt_not_another_allele(tmp_path):
    """A 60-bp unit and its copy one base longer side by side; the row deletes the unit.
    6 carriers span the copy (the deletion written at the anchor), 6 end inside the copy
    (aligned through the unit with its extra base written as a 1-bp insertion: an allele
    their bases spell, on reads that cannot tell it from the deletion), 6 REF reads. The
    spelled insertion is the deletion's own spelling, never another allele: (6, 12, 0)."""
    rng = random.Random(85)
    left = "".join(rng.choice("ACGT") for _ in range(400))
    unit = "".join(rng.choice("ACGT") for _ in range(60))
    extra = next(b for b in "ACGT" if b not in (unit[9], unit[10]))
    copy = unit[:10] + extra + unit[10:]
    ref = left + unit + copy + "".join(rng.choice("ACGT") for _ in range(600))
    row = (400, ref[399:460], ref[399])
    alt = _apply(ref, [row])
    reads = []
    for i in range(6):
        s = 380 - i
        reads.append(
            make_read(
                f"a{i}", alt[s : s + READ_LEN], s, ((0, 400 - s), (2, 60), (0, READ_LEN - 400 + s))
            )
        )
        s = 330 + 2 * i  # ends inside the copy, past its extra base
        k = 410 - s
        reads.append(
            make_read(f"p{i}", alt[s : s + READ_LEN], s, ((0, k), (1, 1), (0, READ_LEN - k - 1)))
        )
        reads.append(_plain(f"r{i}", ref, 360 + i))
    [(rd, ad, _, partial)] = _run(tmp_path, ref, [row], reads, partial=True)
    assert (rd, ad, partial) == (6, 12, 0)


@pytest.mark.parametrize("kind", ["insertion carriers", "deletion REF reads"])
def test_a_sporadic_indel_inside_a_long_events_own_bases_is_no_other_allele(tmp_path, kind):
    """A 120-bp insertion (its carriers soft-clipped at the junction) or a 400-bp deletion;
    6 reads of each allele over the junction, and 6 more of the allele named carrying a
    1-bp deletion each at its own place 30-40 bases inside the event's own bases (the
    inserted bases, or the deleted ones REF holds): a change on that molecule, away from
    the junction that tells the alleles apart, never another allele of the row. Every read
    counts for the allele its junction shows: 12 of that allele, 6 of the other."""
    if kind == "insertion carriers":
        rng = random.Random(120)
        ins = "".join(rng.choice("ACGT") for _ in range(120))
        row = (L0, REF_LONG[L0 - 1], REF_LONG[L0 - 1] + ins)
    else:
        row = (L0, REF_LONG[L0 - 1 : L0 + 400], REF_LONG[L0 - 1])
    alt = _apply(REF_LONG, [row])
    noisy_hap = alt if kind == "insertion carriers" else REF_LONG
    reads = []
    for i in range(6):
        s = L0 - 30 - i
        if kind == "insertion carriers":
            reads.append(
                make_read(f"a{i}", alt[s : s + READ_LEN], s, ((0, L0 - s), (4, READ_LEN - L0 + s)))
            )
        else:
            reads.append(_plain(f"a{i}", alt, s))
        reads.append(_plain(f"r{i}", REF_LONG, s - 10))
        cut = L0 + 30 + 2 * i  # the base each noisy read lacks, inside the event
        seq = noisy_hap[s:cut] + noisy_hap[cut + 1 : s + READ_LEN + 1]
        if kind == "insertion carriers":
            reads.append(make_read(f"n{i}", seq, s, ((0, L0 - s), (4, READ_LEN - L0 + s))))
        else:
            reads.append(
                make_read(f"n{i}", seq, s, ((0, cut - s), (2, 1), (0, READ_LEN - cut + s)))
            )
    [(rd, ad, _, partial)] = _run(tmp_path, REF_LONG, [row], reads, partial=True)
    assert (rd, ad, partial) == ((6, 12, 0) if kind == "insertion carriers" else (12, 6, 0))


def test_a_long_insertions_partial_spelling_stays_x_though_a_noisy_carrier_prefers_it(tmp_path):
    """A 120-bp insertion: the downstream 119 bases with one base (e) after their first 60,
    the next 10 bases repeating its first 10 but for one. 4 carriers ending inside it past
    e, written as a 1-bp insertion of e (an allele their bases spell, X's over their span);
    3 carriers reaching 10 bases past it (151 bp); 1 more with two errors, one turning it to that
    spelling's haplotype there, so it fits the spelling one base better than X yet not
    whole; 6 REF reads. The spelling is X's: one noisy read is no other allele. The 7
    carriers count ALT, the noisy one neither: (6, 7, 0)."""
    rng = random.Random(119)
    left = "".join(rng.choice("ACGT") for _ in range(500))
    d = [rng.choice("ACGT") for _ in range(400)]
    d[119:129] = d[0:10]
    q = 5  # the one base where the 10 after the insertion differ from its first 10
    d[119 + q] = next(b for b in "ACGT" if b != d[q])
    d = "".join(d)
    e = next(b for b in "ACGT" if b not in (d[59], d[60]))
    ins = d[:60] + e + d[60:119]
    ref = left + d + "".join(rng.choice("ACGT") for _ in range(400))
    a = len(left)  # 0-based start of d; the anchor is the base before it
    row = (a, ref[a - 1], ref[a - 1] + ins)
    alt = _apply(ref, [row])
    reads = []
    for i in range(4):  # end inside the insertion past e: REF + e as written
        s = a - 30 - 2 * i
        k = a + 60 - s
        reads.append(
            make_read(f"w{i}", alt[s : s + READ_LEN], s, ((0, k), (1, 1), (0, READ_LEN - k - 1)))
        )
    end = a + 120 + 10  # 10 bases past the insertion, on the alternate haplotype: 151-bp
    for i in range(3):  # reads from the anchor's left
        reads.append(_plain(f"c{i}", alt, end - 151 - i, length=151))
    s = end - 151
    noisy = list(alt[s : s + 151])
    noisy[a + 120 + q - s] = d[119 + q]  # X's base after the insertion read as the spelling's
    noisy[60] = next(b for b in "ACGT" if b != noisy[60])  # and one error elsewhere
    reads.append(make_read("n0", "".join(noisy), s, ((0, 151),)))
    reads += [_plain(f"r{i}", ref, a - 20 + 3 * i) for i in range(6)]  # past e's place
    [(rd, ad, _, partial)] = _run(tmp_path, ref, [row], reads, partial=True)
    assert (rd, ad, partial) == (6, 7, 0)


# ── A sibling allele two of its reads write whole is another allele ──


def test_a_sibling_two_reads_write_whole_is_another_allele(tmp_path):
    """(AT)10 T23 on REF; the row makes it (AT)11 T28, a sibling (unlisted) (AT)12 T28 (the
    chr4:71385191 shape). 5 row carriers and 4 REF reads span it. Of 6 sibling carriers, 2
    span it, their alignment writing the sibling whole; 3 soft-clip the T run's last bases
    and the flank, writing +ATAT (a short-T allele that fits none of them); 1 ends inside
    the T run, written gapless, missing the row by one base and REF by two. Two reads
    writing an allele make it a candidate; the reads telling it from the row and REF make
    it a state: every sibling read is another allele, never the row's ALT: (4, 5, 6)."""
    rng = random.Random(301)
    left = "".join(rng.choice("ACGT") for _ in range(398)) + "GC"
    right = "GAG" + "".join(rng.choice("ACGT") for _ in range(397))
    ref = left + "AT" * 10 + "T" * 23 + right
    a = len(left) + 20  # 1-based anchor: the run's last T
    row = (a, "T", "T" + "AT" + "T" * 5)
    sib_ins = "ATAT" + "T" * 5
    sib = _apply(ref, [(a, "T", "T" + sib_ins)])
    alt = _apply(ref, [row])
    reads = [_plain(f"x{i}", alt, len(left) - 30 - 2 * i) for i in range(5)]
    reads += [_plain(f"r{i}", ref, len(left) - 31 - 2 * i) for i in range(4)]
    for i in range(2):  # the sibling written whole
        reads.append(_gapped(f"w{i}", ref, a, sib_ins, len(left) - 35 - 3 * i))
    for i in range(3):  # aligned through the run, +ATAT, then 23 T's, the rest clipped
        s = len(left) - 40 + 2 * i
        m = a - s
        clip = READ_LEN - m - 4 - 23
        reads.append(
            make_read(f"s{i}", sib[s : s + READ_LEN], s, ((0, m), (1, 4), (0, 23), (4, clip)))
        )
    reads.append(_plain("t0", sib, a + 4 + 19 - READ_LEN))  # ends after 19 of its T's
    [(rd, ad, _, partial)] = _run(tmp_path, ref, [row], reads, partial=True)
    assert (rd, ad, partial) == (4, 5, 6)


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
_r30 = random.Random(30)
INS30 = "".join(_r30.choice("ACGT") for _ in range(30))


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


def test_an_mnp_carrier_with_a_changed_base_below_the_quality_floor_counts_alt(tmp_path):
    """A 2-base MNP: 6 carriers read both changed bases, 4 read the second at Q2 (the
    first decides, as the rules count it, RJ-25), 5 REF reads."""
    m0 = 1000
    ref = _contig(7, length=2000)
    row = (
        m0,
        ref[m0 - 1 : m0 + 1],
        "".join("A" if b != "A" else "C" for b in ref[m0 - 1 : m0 + 1]),
    )
    alt = _apply(ref, [row])
    reads = [_plain(f"f{i}", alt, m0 - 60 + i) for i in range(6)]
    for i in range(4):
        q = [30] * READ_LEN
        q[m0 - (m0 - 50 + i)] = 2
        reads.append(_plain(f"q{i}", alt, m0 - 50 + i, quals=q))
    reads += [_plain(f"r{i}", ref, m0 - 30 + i) for i in range(5)]
    [(rd, ad, dp)] = _run(tmp_path, ref, [row], reads)
    assert (rd, ad, dp) == (5, 10, 15)


# ── An SNV listed beside an MNV holding its change (RJ-28) ──


def test_an_snv_beside_an_mnv_holding_its_change_counts_the_mnvs_carriers(tmp_path):
    """An SNV and a listed DNP whose first base is the SNV's: 6 DNP carriers, 4 carrying
    the SNV alone, 5 REF reads. The SNV's row counts all 10 carriers (RJ-28) and names the
    6 it shares with the DNP; the DNP's row counts its 6, the SNV-only reads partial."""
    m0 = 1000
    ref = _contig(11, length=2000)
    a1 = "A" if ref[m0 - 1] != "A" else "C"
    a2 = "G" if ref[m0] != "G" else "T"
    snv, dnp = (m0, ref[m0 - 1], a1), (m0, ref[m0 - 1 : m0 + 1], a1 + a2)
    h_dnp, h_snv = _apply(ref, [dnp]), _apply(ref, [snv])
    reads = [_plain(f"d{i}", h_dnp, m0 - 60 + i) for i in range(6)]
    reads += [_plain(f"s{i}", h_snv, m0 - 55 + i) for i in range(4)]
    reads += [_plain(f"r{i}", ref, m0 - 30 + i) for i in range(5)]
    got = {
        (r["Reference_Allele"], r["Tumor_Seq_Allele2"]): r
        for r in _run(tmp_path, ref, [snv, dnp], reads, maf_rows=True)
    }
    s, d = got[(snv[1], snv[2])], got[(dnp[1], dnp[2])]
    assert (int(s["ref_count"]), int(s["alt_count"])) == (5, 10)
    assert (int(d["ref_count"]), int(d["alt_count"]), int(d["partial_alt"])) == (5, 6, 4)
    assert f"MNV_SHARED(1:{m0}:{dnp[1]}>{dnp[2]}:6/10)" in s["gbcms_diagnostic"].split(";")
    assert "MNV_SHARED" not in d["gbcms_diagnostic"]


# ── The fit rule: a read no candidate sequence explains ──


def test_a_read_fitting_no_candidate_sequence_counts_only_where_its_aligner_is_sure(tmp_path):
    """The fit rule (xs, decided 2026-10-25): an SNV's carriers carrying twelve confident
    mismatches away from it fit none of the row's candidate sequences within ten. The 3
    the aligner places elsewhere as well (XS = AS) are depth only; the 4 it places here
    alone (XS < AS) are judged by their bases at the row, ALT. 6 clean carriers, 5 REF
    reads. UNEXPLAINED_READS names the 7 of the row's 18 depth reads."""
    m0 = 600
    ref = _contig(53, length=1200)
    alt = "A" if ref[m0 - 1] != "A" else "C"
    snv = (m0, ref[m0 - 1], alt)
    h = _apply(ref, [snv])
    flip = {"A": "C", "C": "G", "G": "T", "T": "A"}

    def noisy(name, start, best, elsewhere):
        # mismatches every 8 read bases, X's base kept: starts one apart never share one
        x = m0 - 1 - start
        seq = list(h[start : start + READ_LEN])
        for i in [k for k in range(0, READ_LEN, 8) if k != x][:12]:
            seq[i] = flip[seq[i]]
        a = make_read(name, "".join(seq), start, ((0, READ_LEN),))
        a.set_tag("AS", best)
        a.set_tag("XS", elsewhere)
        return a

    reads = [_plain(f"c{i}", h, m0 - 70 + 3 * i) for i in range(6)]
    reads += [noisy(f"sure{i}", m0 - 60 + i, 52, 20) for i in range(4)]
    reads += [noisy(f"doubt{i}", m0 - 56 + i, 52, 52) for i in range(3)]
    reads += [_plain(f"r{i}", ref, m0 - 40 + 3 * i) for i in range(5)]
    (row,) = _run(tmp_path, ref, [snv], reads, maf_rows=True)
    assert (int(row["ref_count"]), int(row["alt_count"]), int(row["total_count"])) == (5, 10, 18)
    assert "UNEXPLAINED_READS(7/18)" in row["gbcms_diagnostic"].split(";")


# ── A long deletion with a listed row inside its deleted bases ──


def test_a_listed_snv_inside_a_long_deletions_bases_leaves_its_ref_reads_ref(tmp_path):
    """A 113-bp deletion and a listed SNV 86 bases into its deleted bases: 8 REF reads
    over the anchor short of the SNV, 4 over it carrying the SNV, 6 deletion carriers. The
    deletion's row counts the 12 REF reads REF, the SNV's a change on REF's molecules."""
    d0 = 1000
    ref = _contig(13, length=2000)
    dele = (d0, ref[d0 - 1 : d0 + 113], ref[d0 - 1])
    snv = (d0 + 87, ref[d0 + 86], "A" if ref[d0 + 86] != "A" else "C")
    h_snv = _apply(ref, [snv])
    reads = [_plain(f"r{i}", ref, d0 - 70 + i) for i in range(8)]
    reads += [_plain(f"s{i}", h_snv, d0 - 10 + i) for i in range(4)]
    for i in range(6):
        s = d0 - 60 + i
        seq = ref[s:d0] + ref[d0 + 113 : d0 + 113 + READ_LEN - (d0 - s)]
        reads.append(make_read(f"a{i}", seq, s, ((0, d0 - s), (2, 113), (0, READ_LEN - (d0 - s)))))
    (del_rd, del_ad, _), _snv = _run(tmp_path, ref, [dele, snv], reads)
    assert (del_rd, del_ad) == (12, 6)


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


# ── The guard's recurring changes beside the event explain reads (contexts) ──
# A change the guard masks may take any base: every haplotype takes any base there, so a
# germline haplotype carrying it never reads as X; a recurring indel beside the event is a
# context, crossed with every state.


def _judged_fixture(tmp_path, caplog, monkeypatch, name, **kw):
    """A masked-repeat fixture (tests/testdata/masked_repeat_reads.json, BWA-placed reads)
    counted under the haplotype judge: its counts and per-haplotype read calls."""
    from pathlib import Path

    import test_masked_repeat_reads_contract as masked

    def judged(bam, variants, **k):
        k.setdefault("read_judge", "haplotype")
        k.setdefault("reference_fasta", str(Path(bam).with_suffix(".fa")))
        return count_checked(bam, variants, **k)

    monkeypatch.setattr(masked, "count_checked", judged)
    counts, calls = masked._counted(tmp_path, caplog, name, **kw)
    assert counts.dp >= counts.rd + counts.ad
    assert counts.dpf >= counts.rdf + counts.adf
    assert counts.rd == counts.rd_fwd + counts.rd_rev
    assert counts.ad == counts.ad_fwd + counts.ad_rev
    return counts, calls


def test_a_masked_germline_snp_reading_as_the_insertion_never_counts_alt(
    tmp_path, caplog, monkeypatch
):
    """HG002's 28-bp VNTR row (chr10:129773386) in miniature: no insertion, a germline SNP
    whose haplotype reads as the one-unit insertion over a short read. The guard masks it
    (any base); masked in every haplotype, its reads never count ALT."""
    counts, calls = _judged_fixture(tmp_path, caplog, monkeypatch, "vntr28b_snp_alone")
    assert calls[("h1", "ALT")] == 0, dict(calls)
    assert counts.ad == 0, (counts.rd, counts.ad)


@pytest.mark.parametrize("seed", range(3))
def test_a_minority_masked_snp_haplotype_in_trans_never_counts_alt(
    tmp_path, caplog, monkeypatch, seed
):
    """The insertion on half the molecules, the SNP on a quarter of the others: the
    SNP's reads never count ALT."""
    _, calls = _judged_fixture(
        tmp_path, caplog, monkeypatch, "vntr28b_ins_and_snp", thin=("h2", 0.25, True), seed=seed
    )
    assert calls[("h2", "ALT")] == 0, dict(calls)


def _mnp_beside_indel(tmp_path, reads_of, extra=()):
    """GC>TA in unique sequence (test_complex_exact_contract's shape), counted under the
    judge: (rd, ad, partial)."""
    from test_complex_exact_contract import _files, _flanked, _prepared

    left, ref = _flanked("ACGTCAGT" + "GC" + "TGCAGGTCATTGACCGTAGCTTAG", 7)
    p = len(left) + 8
    hap = ref[:p] + "TA" + ref[p + 2 :]
    fa, bam = _files(tmp_path, ref, reads_of(ref, hap, p) + list(extra))
    (c,) = count_checked(
        bam, [_prepared(fa, "GC", "TA", pos=p)], read_judge="haplotype", reference_fasta=fa
    )
    assert c.dp >= c.rd + c.ad and c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    return c.rd, c.ad, c.partial_alt


@pytest.mark.parametrize("on", ["REF", "ALT"])
def test_a_germline_deletion_beside_an_mnp_explains_its_reads(tmp_path, on):
    """A 1-bp deletion two bases before the block on one haplotype, written by the
    aligner (the guard's recurring indel): its reads keep their call, 40 and 40."""

    def reads_of(ref, hap, p):
        out = []
        for tag, h in (("r", ref), ("a", hap)):
            deleted = (tag == "r") == (on == "REF")
            for i, s in enumerate(range(p - 70, p - 30)):
                if deleted:
                    g = p - 2
                    seq = h[s:g] + h[g + 1 : s + READ_LEN + 1]
                    cig = ((0, g - s), (2, 1), (0, READ_LEN - (g - s)))
                else:
                    seq, cig = h[s : s + READ_LEN], ((0, READ_LEN),)
                out.append(make_read(f"{tag}{i}", seq, s, cig))
        return out

    assert _mnp_beside_indel(tmp_path, reads_of)[:2] == (40, 40)


def test_a_germline_insertion_written_right_after_an_mnp_is_beside_it(tmp_path):
    """A 24-bp duplication of the right flank written as an insertion right after the
    block on 20 of 46 REF reads (the guard's recurring indel, anchored on the block's
    last base): beside the block, it explains those reads, which count REF."""

    def reads_of(ref, hap, p):
        dup = ref[p + 2 : p + 26]
        out = [
            make_read(f"r{i}", ref[s : s + READ_LEN], s, ((0, READ_LEN),))
            for i, s in enumerate(range(p - 70, p - 50))
        ]
        out += [
            make_read(f"a{i}", hap[s : s + READ_LEN], s, ((0, READ_LEN),))
            for i, s in enumerate(range(p - 70, p - 64))
        ]
        for i in range(20):
            s = p - 70 + i
            seq = (ref[s : p + 2] + dup + ref[p + 2 : s + READ_LEN])[:READ_LEN]
            m = p + 2 - s
            out.append(make_read(f"i{i}", seq, s, ((0, m), (1, 24), (0, READ_LEN - m - 24))))
        return out

    assert _mnp_beside_indel(tmp_path, reads_of) == (40, 6, 0)


# ── Another allele at X's site: never REF, never X (#290, #297, #301) ──


def test_another_length_of_the_repeat_unlisted_is_partial_evidence(tmp_path):
    """+CA in (CA)x8 listed; 5 carriers, 4 REF reads, and 4 reads of +CACA (not listed),
    all spanning the tract: the +CACA reads are another allele (#290), never REF or ALT."""
    alt = _apply(REF_REP, [INS_REP])
    other = _apply(REF_REP, [(R0, "G", "GCACA")])
    reads = [_plain(f"a{i}", alt, R0 - 60 + i) for i in range(5)]
    reads += [_plain(f"r{i}", REF_REP, R0 - 60 + i) for i in range(4)]
    reads += [_plain(f"o{i}", other, R0 - 60 + i) for i in range(4)]
    [(rd, ad, dp, partial)] = _run(tmp_path, REF_REP, [INS_REP], reads, partial=True)
    assert (rd, ad, partial) == (4, 5, 4)
    assert dp == 13


def _imperfect_repeat_counts(tmp_path, where):
    """#290's synthetic: an imperfect (GAAGAGG)n repeat; the row inserts its last 33
    bases at its 3' end; 6 REF, 6 ALT and 6 reads of a +15 expansion of the repeat written
    `where` bases upstream of the junction. (rd, ad, partial)."""
    rng = random.Random(290)
    units = ["GAAGAGG"] * 10
    units[3], units[6] = "GAAGGGG", "GAAGAGA"
    left = "".join(rng.choice("ACGT") for _ in range(400))
    right = "".join(rng.choice("ACGT") for _ in range(400))
    contig = left + "".join(units) + right
    r1 = len(left) + 70
    dup = contig[r1 - 33 : r1]
    alt = contig[:r1] + dup + contig[r1:]
    at = r1 - where
    exp = contig[:at] + contig[at - 15 : at] + contig[at:]

    def ins_read(name, hap, ins_at, n, s):
        k = ins_at - s
        return make_read(name, hap[s : s + 151], s, ((0, k), (1, n), (0, 151 - k - n)))

    reads = []
    for i in range(6):
        s = len(left) - 30 + 2 * i
        reads.append(make_read(f"r{i}", contig[s : s + 151], s, ((0, 151),)))
        reads.append(ins_read(f"a{i}", alt, r1, 33, s))
        reads.append(ins_read(f"x{i}", exp, at, 15, s))
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", r1 - 1, contig[r1 - 1], contig[r1 - 1] + dup, "INS")],
        fa,
        5,
        False,
        1,
        True,
    )
    (c,) = count_checked(bam, [pv.variant], read_judge="haplotype", reference_fasta=fa)
    assert c.dp >= c.rd + c.ad and c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    return c.rd, c.ad, c.partial_alt


@pytest.mark.parametrize("where", [0, 30, 60])
def test_another_length_of_an_imperfect_repeat_is_partial_wherever_written(tmp_path, where):
    """Decided on #290: the +15 reads are another allele wherever the aligner wrote the
    expansion in the repeat: (6, 6, 6)."""
    assert _imperfect_repeat_counts(tmp_path, where) == (6, 6, 6)


@pytest.mark.parametrize("differ", [(5, 6), (5,)], ids=["DNP", "one base"])
def test_another_unit_count_of_a_vntr_is_partial_evidence(tmp_path, differ):
    """#297: three 60-bp units after a long deletion's anchor whose units 2 and 3 differ;
    the row deletes one unit. 4 REF reads, 4 carriers of the two-unit deletion (another
    allele, never ALT), 4 of the row's: (4, 4, 4)."""
    import read_judgment_cases as rj

    a = rj.A
    n, mh = rj.LONG_DEL["del 60"]
    c = list(rj._long_contig(n, mh))
    u1 = c[a : a + n]
    u2 = list(u1)
    u2[2] = rj._other(u1[2])
    u3 = list(u2)
    for k in differ:
        u3[k] = rj._other(u2[k])
    c[a + n : a + 2 * n], c[a + 2 * n : a + 3 * n] = u2, u3
    contig = "".join(c)
    one, two = contig[:a] + contig[a + n :], contig[:a] + contig[a + 2 * n :]
    reads = []
    for i in range(4):
        s = a - 50 - i
        reads.append(make_read(f"r{i}", contig[s : s + 100], s, ((0, 100),)))
        reads.append(
            make_read(f"x{i}", two[s : s + 100], s, ((0, a - s), (2, 2 * n), (0, 100 - (a - s))))
        )
        reads.append(
            make_read(f"d{i}", one[s : s + 100], s, ((0, a - s), (2, n), (0, 100 - (a - s))))
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
        [_rs.Variant("1", a - 1, contig[a - 1 : a + n], contig[a - 1], "DEL")],
        fa,
        5,
        False,
        1,
        True,
    )
    (got,) = count_checked(bam, [pv.variant], read_judge="haplotype", reference_fasta=fa)
    assert got.dp >= got.rd + got.ad and got.rd == got.rd_fwd + got.rd_rev
    assert got.ad == got.ad_fwd + got.ad_rev
    assert (got.rd, got.ad, got.partial_alt) == (4, 4, 4)


# ── A delins: reads ending inside the event are judged alike on either allele ──


def test_reads_ending_inside_a_delins_count_alike_on_either_allele(tmp_path):
    """A 24>6 delins in unique sequence: 10 REF and 10 ALT reads spanning it, and 10 of
    each ending three bases into its event (its first bases differ). Each read is judged
    by the bases it holds, on REF and ALT alike: (20, 20)."""
    ref = _contig(274, length=1600)
    p0 = 800  # 0-based start of the 24 REF bases
    rng = random.Random(6)
    alt6 = ""
    while not alt6 or alt6[0] == ref[p0] or alt6[-1] == ref[p0 + 23]:
        alt6 = "".join(rng.choice("ACGT") for _ in range(6))
    row = (p0 + 1, ref[p0 : p0 + 24], alt6)
    alt = _apply(ref, [row])
    reads = [_plain(f"r{i}", ref, p0 - 50 + i) for i in range(10)]
    reads += [_plain(f"a{i}", alt, p0 - 50 + i) for i in range(10)]
    reads += [_plain(f"e{i}", ref, p0 + 3 - READ_LEN) for i in range(10)]
    reads += [_plain(f"f{i}", alt, p0 + 3 - READ_LEN) for i in range(10)]
    [(rd, ad, _)] = _run(tmp_path, ref, [row], reads)
    assert (rd, ad) == (20, 20)


# ── Clip admission is bounded by the read's fragment (RJ-13) ──


def test_a_clip_admits_a_read_only_inside_its_fragment(tmp_path):
    """An SNV; 6 REF reads over it; 8 carriers aligned from the base after it with it in
    a 20-base leading soft clip: paired inward, their clip is inside the fragment and
    admits them; unpaired, there is no fragment to bound it, and they are not admitted."""
    ref = _contig(13, length=1600)
    x = (S0, ref[S0 - 1], "G" if ref[S0 - 1] != "G" else "T")
    alt = _apply(ref, [x])

    def carriers(paired):
        out = []
        for i in range(8):
            hs = S0 - 20 + i  # the clip holds X and the bases before it
            k = S0 - hs
            r = make_read(f"c{i}", alt[hs : hs + READ_LEN], S0, ((4, k), (0, READ_LEN - k)))
            if paired:
                r.flag = 0x1 | 0x2 | 0x20 | 0x40
                r.next_reference_id, r.next_reference_start, r.template_length = (
                    0,
                    S0 + 150,
                    250 + k,
                )
            out.append(r)
        return out

    base = [_plain(f"r{i}", ref, S0 - 50 + i) for i in range(6)]
    for name, paired, want in (("paired", True, (6, 8, 14)), ("unpaired", False, (6, 0, 6))):
        d = tmp_path / name
        d.mkdir()
        [(rd, ad, dp)] = _run(d, ref, [x], base + carriers(paired))
        assert (rd, ad, dp) == want, name


# ── Contexts: the haplotypes a read's molecule carries are formed whole ──
# HG002 (public) at two #246 rows, its GRCh38 reference over +-420 bases and its phased
# truth rows there; reads tiled from each truth haplotype, written with its indels.

HLA_G0 = 32701977
HLA = (
    "TGTCAGCCTATGATTGATGAAAGATTCTGCTTAAACAATTGAGTCAGTAAGGCTTCCACACTTTGCCATTGAACTGAGTGGGATGTGG"
    "GAAATGCTTTCAAAGTTCCGGGACTTTACAGGTATGTCCCAATTTTTACTTTCTGCCTTCACATGTCCTTATATTCAGCCAGGTCTCT"
    "CCTGAGCAGGCTTAACCAAGCTGATGCAGACAGCTTTCCACACCACTGGGGATAAATGAGATTTTAGCAAGGCTATTTTTGACTATCT"
    "CATTCTCTGGATCTCCCTGGTAAAATTCCTGACTGGCCTGCCATCATGACTTCAAGCTAGTTGAGATGTTAGCCTTCCATAATTATTT"
    "GCCCTGAGGTCTATGCTGTTTTCTACAAAATGCCAAGACGTTATTTTGAGGGTGTGTGTGTGTGTGTGTGTGTGTGTGTGCATGCATG"
    "TGCGTGAGCTCTTTTCCAAAAAATGTGTGTGAGCTCCCTCTGGCAGGTTAAGAAGAGCTGTCTAATCCCATGGCCTGTCCTCCCACAC"
    "AGGTAGAAATTCTGCACCATAAGTTGAGGTGAAGGAAGTTGGTGAAAGCATACTCTGGCCAAAATAATACTGATTTTGCTGTTTTTAT"
    "TAGGATTCAGTTTAGATTCTCTTAACTAAATGCTTCTCAATGTTCGTATCCTTTGGTCAATTTCCAGAGTCTTTAAATGATTGTTTTT"
    "GGCAACGTTGTCCAATGTTCTCACTGCTTTTTAGATTCGTTGAAATGTTTTGATGTCAATATATGTGGTATATTTTTATGTATTTTCT"
    "ACGTGAATTTGAAAATAAGCTGTATATAGCATTTGGTACATACACACT"
)
HLA_ROWS = [
    (32701984, "C", "G", 1, 0),
    (32702013, "A", "AC", 1, 1),
    (32702040, "T", "C", 0, 1),
    (32702050, "C", "G", 1, 1),
    (32702056, "G", "A", 0, 1),
    (32702060, "T", "C", 1, 1),
    (32702061, "G", "T", 1, 1),
    (32702101, "T", "C", 1, 1),
    (32702106, "A", "G", 1, 1),
    (32702178, "G", "T", 1, 1),
    (32702260, "G", "T", 1, 1),
    (32702333, "T", "C", 1, 1),
    (32702359, "T", "C", 1, 1),
    (32702368, "G", "A", 1, 1),
    (32702374, "T", "C", 0, 1),
    (32702380, "G", "T", 1, 1),
    (32702384, "G", "A", 0, 1),
    (32702394, "G", "GCGCA", 1, 0),
    (32702395, "T", "C", 0, 1),
    (32702396, "G", "GCA", 0, 1),
    (32702397, "T", "C", 1, 0),
    (32702398, "G", "A", 1, 0),
    (32702399, "T", "C", 1, 1),
    (32702400, "G", "A", 0, 1),
    (32702401, "T", "C", 0, 1),
    (32702407, "T", "C", 1, 0),
    (32702408, "G", "A", 1, 0),
    (32702410, "A", "ACG", 0, 1),
    (32702419, "C", "T", 1, 1),
    (32702437, "A", "C", 1, 1),
    (32702467, "A", "G", 1, 1),
    (32702531, "A", "G", 1, 1),
    (32702591, "A", "G", 0, 1),
    (32702637, "C", "T", 0, 1),
    (32702687, "G", "T", 1, 1),
    (32702700, "C", "A", 1, 1),
    (32702718, "C", "T", 0, 1),
    (32702771, "G", "A", 1, 0),
    (32702787, "G", "T", 1, 1),
    (32702799, "A", "G", 1, 1),
]
CHR19_G0 = 48568927  # 1-based genome position of the window's first base
CHR19 = (
    "GCTCCCTCTCTCGCTACCTTCAGTCTTTGCTCAGATGTCGAGTTCCCAGAGGGGCCTCCCTGACGCCACCGTTCTAGCAGCATTTAGCA"
    "TTTAGATAAATGACAAATTTTAGATTAAATGTTAGATAAATTTTAGGCTGGGCGTGGTGGCTCACGCCTGTAATCCCAGCACTTTAGGA"
    "GGCTGAGGCAGGCAGATCACAAGGTCAGGAGATCGAGACCATCCTGGGTAACATGGTGAAACCCTGTCTCTACTAAAAATACAAAAAAT"
    "TAGCAGGGCGTGGTGGCGGGCGCCTGTAGTCCCAGCTACTCAGGAGGCTGAGGCAGCAGAATGGCGTGAATCCAGGAGGCGGAGCTTGC"
    "AGTGAGCCAAGATCGCGCCACTGCACTCCAGCCTGGGCGACAGAGAGAGACTCCGTCTCAAAAAAAAAAAAAAAAAAAAACATATATAT"
    "ATATATATATATATATATATATATATATGATAAATTTTAGATAAATTTTTATCTAAAATTGTAACTCACTGGATTCAGTGACATCTATC"
    "TGTAATCCCAGCACTTTAGGAGGCTGAGGTGGGAGGATCGCTTGAGCCCAGGAGTTTTCTTTCTTTTTTTTCTTATTTTTTAAAGAGAT"
    "TGGGTTCTCACTCTGTGTAGCCCAGGCTGGTCTTGAACTCCTGGCCTCAAGCAGTCCTCCCATCTAGGCCTCCCAAAGTGTTGGGATTA"
    "CAGGCATGAGCCACCATGCCTGGCCTCAAGATTGGTTTTTTGTTCTGAGACACAGACTCACTCTGTCACCCAGGCTGGAGTGCAGTGAC"
    "ATGATCCCAGCTCACCATGACCTCCGCCTCCCAGGTTAA"
)
# (1-based POS, REF, ALT, on haplotype 0, on haplotype 1)
CHR19_ROWS = [
    (48569203, "G", "A", 1, 1),
    (48569210, "C", "T", 1, 1),
    (48569299, "G", "C", 1, 1),
    (48569300, "C", "A", 1, 1),
    (48569328, "A", "C", 1, 1),
    (48569346, "A", "AGAAAAG", 0, 1),
    (48569347, "A", "G", 1, 0),
    (48569351, "A", "G", 0, 1),
    (48569352, "A", "G", 1, 0),
    (48569357, "A", "G", 1, 0),
    (48569359, "AAAACATATATATATATATATATATAT", "A", 0, 1),
    (48569361, "AACATATATATATATATATAT", "A", 1, 0),
    (48569383, "T", "A", 1, 0),
    (48569398, "A", "G", 1, 1),
    (48569399, "T", "A", 1, 1),
    (48569433, "T", "C", 1, 1),
    (48569453, "CATCT", "C", 1, 1),
    (48569459, "T", "C", 1, 1),
    (48569537, "T", "C", 1, 1),
]
CHR3_G0 = 112232632
CHR3 = (
    "TATTGTGTTTTGCAAAAAAAACAAAAACAAAAACATTGTAAAAACCATTGAAAACTGGACATGATGTAGGACACGGAAAAGTGCAGCAG"
    "CTGTGAAGACCCCAGGACAACTTTTATTATTTCCCTTTAGGGAAACCATGACTGATTAGTGTGCTGAGTGTTCTTAAAAGAATGAGTCT"
    "GACTGCAGATATCTTAAAGCATATCTCAGAAGCAGAAGGTGACATCTCCCACAGCTAAACTTTGGCTTTGAATTCTATTTTATATTGA"
    "AAATATTTTTAACCATAAATTAAGGTCTAGCAAATGCTTAGATACACAATAGCCTAAATAAATAAACAAAGACCCCTGACCCCTCCCAA"
    "TTTGGTCAAAAACTAGCCATTTTCTAAGGTTCACTTTCATACACACACACACACACACACACACATATATATATATATATATTATATTT"
    "TTTTTTTTTTTAAGAAAGGGTATCTCTTTGCTGCCCAGGCTGGAGTGCAGTGGCACAATCACGGTTCACCACAGCCTCAACCTCCCAGG"
    "CTTGAGTGATCTTCCTACCTCAGTCTCCCGAGTAGCGGGGACTACAGATGCACACCACCACCCAACTAATTTTTTGTATATTTTATAGA"
    "CACATGGTTTTGGCATGTTGCCAGGTTCGTCTCAAACTGCTTGGCTCAAGCAATTCACCTGCCTTGGCTTCCCAAAGTGTTAAGATTAC"
    "AGGCATGAGCCACTGTGCTGACCTATGATGTCTTTTAAAGGGAGCTCCAATACACAGCCTCCTCTTTTGAATTTAAACATATATATTTT"
    "AAACATACACATAAAACATATATTTTATTTTCAAATTTTT"
)
CHR3_ROWS = [
    (112233026, "T", "TACAC", 1, 0),
    (112233051, "A", "ACACACACATATATATATATAT", 0, 1),
    (112233052, "T", "C", 1, 0),
    (112233054, "T", "C", 1, 0),
    (112233067, "ATTATATTT", "A", 0, 1),
    (112233073, "T", "TATATATATA", 1, 0),
    (112233180, "T", "C", 1, 1),
]


def _haplotype_with_origins(ref, rows):
    """`ref` with `rows` ((1-based local POS, REF, ALT), ...) applied: the haplotype, and
    per base its reference offset (None for an inserted base)."""
    hap, origin, at = [], [], 0
    for pos, r, a in sorted(rows):
        p = pos - 1
        assert ref[p : p + len(r)] == r, (pos, r)
        hap += ref[at:p]
        origin += list(range(at, p))
        pre = 0
        while pre < min(len(r), len(a)) and r[pre] == a[pre]:
            pre += 1
        hap += a
        shared = min(len(r), len(a))
        origin += [p + k if k < shared else None for k in range(len(a))]
        at = p + len(r)
    hap += ref[at:]
    origin += list(range(at, len(ref)))
    return "".join(hap), origin


def _read_on(name, hap, origin, s, length=151):
    """A read of `hap` from offset `s` (its first base a reference one), with the CIGAR
    its origins give (deletions where the reference offsets jump, insertions where none)."""
    while origin[s] is None:
        s += 1
    seq = hap[s : s + length]
    cig, last = [], None
    for o in origin[s : s + len(seq)]:
        if o is None:
            op = 1
        else:
            if last is not None and o > last + 1:
                cig.append([2, o - last - 1])
            op, last = 0, o
        if cig and cig[-1][0] == op:
            cig[-1][1] += 1
        else:
            cig.append([op, 1])
    while cig and cig[-1][0] == 1:  # a read ends on an aligned base
        seq = seq[: -cig[-1][1]]
        cig.pop()
    return make_read(name, seq, origin[s], tuple(tuple(c) for c in cig))


def _two_haplotype_calls(tmp_path, caplog, window, g0, rows, x_pos, written=None):
    """Every row listed; reads tiled every 4 bases from each truth haplotype, written
    with its own indels or, for a haplotype in `written`, as that path gives each of its
    bases a reference offset (an aligner's other representation). Per haplotype, X's
    read calls (ALT / REF / neither)."""
    import collections
    import re

    local = [(p - g0 + 1, r, a) for p, r, a, _, _ in rows]
    reads = []
    for k in (0, 1):
        on = [(p - g0 + 1, r, a) for p, r, a, h0, h1 in rows if (h0, h1)[k]]
        hap, origin = _haplotype_with_origins(window, on)
        if written and k in written:
            origin = [written[k](i) for i in range(len(hap))]
        reads += [_read_on(f"h{k}_{s}", hap, origin, s) for s in range(0, len(hap) - 151, 4)]
    fa, bam = write_contig(tmp_path, window, reads, "w")
    pvs = _rs.prepare_variants(
        [_rs.Variant("1", p - 1, r, a, "X") for p, r, a in local], fa, 5, False, 1, True
    )
    vs = [pv.variant for pv in pvs]
    xi = next(i for i, (p, _, _) in enumerate(local) if p == x_pos - g0 + 1)
    caplog.clear()
    with caplog.at_level(5, logger="_rs.counting.engine"):
        _rs.reset_log_caching()
        (c,) = count_checked(bam, [vs[xi]], read_judge="haplotype", reference_fasta=fa)
    assert c.dp >= c.rd + c.ad and c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    key = f"1:{vs[xi].pos + 1} {vs[xi].ref_allele}>{vs[xi].alt_allele} "
    call = re.compile(r"read=(h[01])_\d+ mate=\d ref=(\w+) alt=(\w+)")
    tally = collections.Counter()
    for record in caplog.records:
        msg = record.getMessage()
        if "read call " in msg and key in msg and (m := call.search(msg)):
            hap, is_ref, is_alt = m.groups()
            tally[
                (hap, "ALT" if is_alt == "true" else "REF" if is_ref == "true" else "neither")
            ] += 1
    return tally


@pytest.mark.parametrize("x_pos, carrier", [(32702397, "h0"), (32702395, "h1")])
def test_each_haplotypes_rows_are_a_context_at_a_dense_hla_site(tmp_path, caplog, x_pos, carrier):
    """chr6:32702395/97 (HLA-DQB1), two SNVs in trans among 40 rows, each haplotype with
    its own insertions in one CA repeat (+CGCA, +CA). Each read's insertion is chosen
    with the rows its molecule carries, never both insertions together (each explains
    part of the other's length), so each haplotype's rows are a context: X's carriers
    count ALT (they tied across every state, depth only), the other haplotype's never."""
    tally = _two_haplotype_calls(tmp_path, caplog, HLA, HLA_G0, HLA_ROWS, x_pos)
    other = "h1" if carrier == "h0" else "h0"
    assert tally[(other, "ALT")] == 0, dict(tally)
    assert tally[(carrier, "ALT")] >= 3 * tally[(carrier, "neither")], dict(tally)
    assert tally[(carrier, "REF")] == 0, dict(tally)


def test_a_molecules_whole_context_is_formed_from_reads_that_each_show_part_of_it(tmp_path, caplog):
    """chr19:48569347 A>G (#246: one base after a +GAAAAG whose first base is G): on the
    other haplotype the insertion, an SNV and a 26-bp deletion; on X's, SNVs and a
    21-bp deletion at the same run. Reads of the insertion's haplotype read as X beside
    the other haplotype's rows unless their own molecule's rows are a context together:
    they are never ALT."""
    # The insertion's haplotype as BWA-MEM wrote HG002's reads there (76M1I16M21D58M):
    # one base inserted after the anchor, sixteen aligned, then 21 deleted.
    anchor = 48569346 - CHR19_G0 + 1  # its first base after the anchor, 0-based

    def bwa(i):
        if i < anchor:
            return i
        if i == anchor:
            return None
        return i - 1 if i < anchor + 17 else i + 20

    tally = _two_haplotype_calls(tmp_path, caplog, CHR19, CHR19_G0, CHR19_ROWS, 48569347, {1: bwa})
    assert tally[("h1", "ALT")] == 0, dict(tally)
    assert tally[("h0", "ALT")] > 0 and tally[("h1", "REF")] > 0, dict(tally)


def test_a_reads_context_is_chosen_whichever_allele_it_holds_at_x(tmp_path, caplog):
    """chr3:112233052 T>C beside a 21-bp insertion one base before it (#246) in a CA/TA
    repeat, each haplotype with its own indels: X's carriers fit the other haplotype's
    rows as well as their own with REF at X, so the context is found with X too; they
    are never REF, and the other haplotype's reads never ALT."""
    tally = _two_haplotype_calls(tmp_path, caplog, CHR3, CHR3_G0, CHR3_ROWS, 112233052)
    assert tally[("h0", "REF")] == 0 and tally[("h1", "ALT")] == 0, dict(tally)
    assert tally[("h0", "ALT")] > 0, dict(tally)


# ── A dense cluster: the context enumeration is bounded ──


def test_a_dense_cluster_of_rows_is_judged_in_bounded_time(tmp_path):
    """The HG002 pericentromeric cluster (chr3:93.46-93.49 Mb, a listed SNV on nearly
    every base, 5,000+ reads) in miniature: 60 SNV rows six bases apart on one
    haplotype, 4,000 reads of 150 bases from it and from REF, a third of the bases at
    the rows read at quality 2, so nearly every read shows its own set of rows. The sets
    reads show together are enumerated within bounds: one row is judged in well under a
    minute (the unbounded enumeration took about a minute here), its carriers ALT."""
    import time

    rng = random.Random(93470262)
    contig = "".join(rng.choice("ACGT") for _ in range(1800))
    pos = [500 + 6 * k for k in range(60)]
    rows = [
        (p, contig[p], "ACGT"[("ACGT".index(contig[p]) + 1 + k % 3) % 4]) for k, p in enumerate(pos)
    ]
    hap = list(contig)
    for p, _, a in rows:
        hap[p] = a
    hap = "".join(hap)
    reads = []
    for i in range(4000):
        h = (hap, contig)[i % 2]
        st = 360 + rng.randrange(0, 500)
        quals = [30] * 150
        for p in pos:
            if st <= p < st + 150 and rng.random() < 0.33:
                quals[p - st] = 2
        reads.append(make_read(f"r{i}", h[st : st + 150], st, ((0, 150),), quals=quals))
    fa, bam = write_contig(tmp_path, contig, reads, "c")
    pvs = _rs.prepare_variants(
        [_rs.Variant("1", p, r, a, "SNP") for p, r, a in rows], fa, 5, False, 1, True
    )
    x = pvs[30].variant
    t0 = time.monotonic()
    (c,) = _rs.count_bam_binned(
        bam,
        [x],
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
        reference_fasta=fa,
        read_judge="haplotype",
    )
    elapsed = time.monotonic() - t0
    assert c.dp >= c.rd + c.ad and c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev and c.ad == c.ad_fwd + c.ad_rev
    assert elapsed < 40, f"{elapsed:.1f}s"
    assert c.ad > 0 and c.rd > 0, (c.rd, c.ad)
