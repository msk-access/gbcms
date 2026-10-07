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
import pytest
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

T0 = 300  # 0-based anchor A, then (CTG)x8, then A
UNITS = 8
READ_LEN = 100


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


def _run(tmp_path, rows, reads):
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
    )
    assert result.exit_code == 0, result.output
    got = list(read_maf_output(glob.glob(str(out / "*.maf"))[0]))
    for r in got:
        assert int(r["total_count"]) >= int(r["ref_count"]) + int(r["alt_count"])
    return [(int(r["ref_count"]), int(r["alt_count"]), int(r["partial_alt"])) for r in got]


@pytest.mark.xfail(
    strict=True, reason="each twin is the other's sibling: mid-tract carriers count at neither"
)
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


@pytest.mark.xfail(strict=True, reason="no twin finder yet")
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
