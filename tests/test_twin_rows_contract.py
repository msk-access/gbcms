"""C38 #246 follow-up: one allele given twice is not its own sibling.

An input can give one allele as two rows: an insertion of a repeat unit written at
either end of the repeat (FORTE's indel probes list +CTG at both ends of a (CTG)x11
tract). Preparation normalizes both to one variant, and the pipeline made each row
the other's sibling. A carrier written at the junction is an exact op the AD guard
never contests, so it counted at both rows; a carrier the aligner wrote anywhere
else in the tract fit both rows' haplotypes exactly, so the guard called it "either
allele, neither row's AD" and it counted at neither (measured: +CTG AD 7 with the
twin, 62 without it, as when counted alone; not binning).

Contract (operator, 2026-10-07): a sibling identical to the row after preparation
is left out of its sibling list, so each twin row counts every carrier, and the run
warns, naming the rows. A genuinely different allele at the site stays a sibling.

Runs through the production CLI (pipeline sibling assembly).
"""

import glob
import random

import pysam
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

T0 = 300  # 0-based anchor A, then (CTG)x8, then A
UNITS = 8
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


def _contig():
    rng = random.Random(246)
    c = [rng.choice("ACGT") for _ in range(800)]
    c[T0 - 1 : T0 + 2 + 3 * UNITS] = "GA" + "CTG" * UNITS + "A"
    return "".join(c)


REF = _contig()
END = T0 + 3 * UNITS  # 0-based: the tract's last G


def _ins_read(name, i, junction, ins):
    s = T0 - 40 - (i % 5)
    left = junction - s
    seq = REF[s:junction] + ins + REF[junction : s + READ_LEN - len(ins)]
    return make_read(name, seq, s, ((0, left), (1, len(ins)), (0, READ_LEN - left - len(ins))))


def _run(tmp_path, rows, reads, log=False):
    fa = tmp_path / "ref.fasta"
    fa.write_text(">1\n" + REF + "\n")
    pysam.faidx(str(fa))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": len(REF)}]}
    raw, bam = tmp_path / "raw.bam", tmp_path / "twins.bam"
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
        ["dna", "-v", str(vcf), "-b", f"S:{bam}", "-f", str(fa), "-o", str(out), "--format", "maf"],
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
    counts = [(int(r["ref_count"]), int(r["alt_count"]), int(r["partial_alt"])) for r in got]
    return (counts, result.output) if log else counts


def test_one_allele_given_at_both_ends_of_a_repeat_counts_every_carrier(tmp_path):
    """+CTG given after the anchor A and after the tract's last G (VCF rows in input
    order), and +CTGCTG, a genuinely different allele. Carriers of +CTG: 4 written at
    the junction, 6 mid-tract; 3 carriers of +CTGCTG; 6 REF reads. Both +CTG rows count
    all 10 carriers; the +CTGCTG carriers stay another allele there."""
    rows = [(T0 + 1, "A", "ACTG"), (END + 1, "G", "GCTG"), (T0 + 1, "A", "ACTGCTG")]
    reads = [_ins_read(f"j{i}", i, T0 + 1, "CTG") for i in range(4)]
    reads += [_ins_read(f"m{i}", i, T0 + 13, "CTG") for i in range(6)]
    reads += [_ins_read(f"d{i}", i, T0 + 1, "CTGCTG") for i in range(3)]
    reads += [
        make_read(f"r{i}", REF[T0 - 50 + i : T0 + 50 + i], T0 - 50 + i, ((0, READ_LEN),))
        for i in range(6)
    ]
    assert _run(tmp_path, rows, reads) == [(6, 10, 3), (6, 10, 3), (6, 3, 10)]


def test_twin_rows_are_found_after_normalization():
    """The pipeline's twin finder groups valid rows by their prepared allele."""
    from gbcms.pipeline import duplicate_alleles

    class V:
        def __init__(self, chrom, pos, ref, alt):
            self.chrom, self.pos, self.ref_allele, self.alt_allele = chrom, pos, ref, alt

    class P:
        def __init__(self, v):
            self.variant = v

    prepared = [
        P(V("1", 300, "A", "ACTG")),
        P(V("1", 300, "A", "ACTGCTG")),
        P(V("1", 300, "a", "actg")),
    ]
    assert duplicate_alleles(prepared, [0, 1, 2]) == [[0, 2]]
    assert duplicate_alleles(prepared, [0, 1]) == []


# ── From the review ──────────────────────────────────────────────────────────


def _masked_mid_tract(name, i):
    """+CTG written mid-tract with its base on the tract's first C (where +CTG and
    +TTG differ) at BQ 5: it fits +CTG and +TTG exactly."""
    a = _ins_read(name, i, T0 + 13, "CTG")
    q = list(a.query_qualities)
    q[T0 + 1 - a.reference_start] = 5
    a.query_qualities = q
    return a


def test_a_genuine_sibling_still_claims_what_it_explains_as_well(tmp_path):
    """+CTG given twice and +TTG, a same-length allele one base apart. Carriers of +CTG
    masked at that base, written mid-tract, are either allele: neither row's AD at
    the +CTG twins or at +TTG, as when +CTG is given once (the twin is left out, the
    genuine sibling is not)."""
    rows = [(T0 + 1, "A", "ACTG"), (END + 1, "G", "GCTG"), (T0 + 1, "A", "ATTG")]
    reads = [_ins_read(f"j{i}", i, T0 + 1, "CTG") for i in range(4)]
    reads += [_masked_mid_tract(f"x{i}", i) for i in range(5)]
    reads += [_ins_read(f"d{i}", i, T0 + 1, "TTG") for i in range(3)]
    reads += [
        make_read(f"r{i}", REF[T0 - 50 + i : T0 + 50 + i], T0 - 50 + i, ((0, READ_LEN),))
        for i in range(6)
    ]
    assert _run(tmp_path, rows, reads) == [(6, 4, 8), (6, 4, 8), (6, 3, 9)]


def test_an_mnp_given_twice_is_rescued_as_when_given_once(tmp_path):
    """With --rescue-mnp, an MNP row given twice is rescued at both rows, as when given
    once: its group holds no other allele, so no co-annotated row owns its reads."""
    import test_rescue_mnp as m

    def counts(rows, tag):
        d = tmp_path / tag
        d.mkdir()
        outdir, _ = m._invoke(d, {"S": m._component_carrier_reads()}, rows, True, "maf")
        out = list(read_maf_output(glob.glob(str(outdir / "S.maf"))[0]))
        return [(r["ref_count"], r["alt_count"], r["partial_alt"], r["gbcms_rescue"]) for r in out]

    (alone,) = counts([m.MNP_ROW], "once")
    assert "outcome=rescued" in alone[3]
    assert counts([m.MNP_ROW, m.MNP_ROW], "twice") == [alone, alone]


def test_a_row_given_verbatim_is_noted_and_two_ways_is_warned(tmp_path):
    """A row given verbatim twice (a cohort MAF listing a recurrent indel once per
    sample) is noted at INFO; one allele given two ways is a WARNING naming both."""
    reads = [_ins_read(f"j{i}", i, T0 + 1, "CTG") for i in range(4)]

    def noted(rows, tag):
        (tmp_path / tag).mkdir()
        _, out = _run(tmp_path / tag, rows, reads, log=True)
        return [line for line in out.splitlines() if " given " in line]

    verbatim = noted([(T0 + 1, "A", "ACTG"), (T0 + 1, "A", "ACTG")], "verbatim")
    assert verbatim and all(" INFO " in line for line in verbatim)
    two_ways = noted([(T0 + 1, "A", "ACTG"), (END + 1, "G", "GCTG")], "two")
    assert two_ways and all(" WARNING " in line for line in two_ways)
    assert f"1:{END + 1} G>GCTG" in two_ways[0]
