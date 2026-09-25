"""Count the given allele; name the allele the reads carry (O4 #163, C3 #111).

gbcms takes the input allele as correct and counts it. When the reads carry a
different allele (a curator's edit, a caller merging events), the counts stay
the given allele's, and ``gbcms_diagnostic`` names what the reads show:
``OBSERVED_ALLELE(chrom:pos:REF>ALT:n/m)``, with n reads carrying it exactly and
m reads carrying the given ALT exactly.

The homopolymer twin (a delins ``CCCCCC>T`` also counted as ``CCCCCT``, the
winner reported) is opt-in: ``--rescue-homopolymer``. The case it was built
for, SOX2, is a 1bp deletion plus C>T (``CCCCT``) that the twin only matched
by tolerance; the diagnostic names it exactly.
"""

import glob
import random

import pysam
import pytest
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

RUN = 300  # 0-based start of the C run
RUN_LEN = 6
READ = 100


def _ref():
    # A/G-only flanks keep the C run and its following T unambiguous.
    rng = random.Random(163)
    left = "".join(rng.choice("AG") for _ in range(RUN))
    right = "".join(rng.choice("AG") for _ in range(300))
    return left + "C" * RUN_LEN + "T" + right


def _files(tmp_path, ref, reads):
    fa = tmp_path / "ref.fa"
    fa.write_text(f">1\n{ref}\n")
    pysam.faidx(str(fa))
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": len(ref)}]}
    raw = tmp_path / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=hdr) as fh:
        for r in reads:
            fh.write(r)
    bam = tmp_path / "s.bam"
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    return str(fa), str(bam)


def _run(tmp_path, ref, reads, rows, *extra):
    fa, bam = _files(tmp_path, ref, reads)
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=607>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(f"1\t{p}\t.\t{r}\t{a}\t.\t.\t.\n" for p, r, a in rows)
    )
    out = tmp_path / "out"
    res = CliRunner().invoke(
        app,
        ["dna", "-v", str(vcf), "-b", f"S:{bam}", "-f", fa, "-o", str(out), "--format", "maf"]
        + list(extra),
    )
    assert res.exit_code == 0, res.output
    return list(read_maf_output(glob.glob(str(out / "*.maf"))[0]))


def _ref_reads(ref, n, start=250):
    return [
        make_read(f"r{i}", ref[start + i : start + i + READ], start + i, ((0, READ),))
        for i in range(n)
    ]


def _sox2_carriers(ref, n):
    """A 1bp deletion in the run plus its last C read as T: CCCCCC -> CCCCT."""
    hap = ref[:RUN] + "CCCCT" + ref[RUN + RUN_LEN :]
    reads = []
    for i in range(n):
        s = 240 + i
        left = RUN + 4 - s  # bases before the deleted C (0-based 304)
        reads.append(
            make_read(f"c{i}", hap[s : s + READ], s, ((0, left), (2, 1), (0, READ - left)))
        )
    return reads


SOX2_ROW = (RUN + 1, "C" * RUN_LEN, "T")


@pytest.mark.xfail(strict=True, reason="the twin is dual-counted by default")
def test_the_given_allele_is_counted_by_default(tmp_path):
    ref = _ref()
    (row,) = _run(tmp_path, ref, _sox2_carriers(ref, 20) + _ref_reads(ref, 10), [SOX2_ROW])
    assert "WARN_HOMOPOLYMER_DECOMP" not in row["gbcms_status_reason"]


@pytest.mark.xfail(strict=True, reason="--rescue-homopolymer does not exist yet")
def test_rescue_homopolymer_keeps_the_twin(tmp_path):
    ref = _ref()
    reads = _sox2_carriers(ref, 20) + _ref_reads(ref, 10)
    (row,) = _run(tmp_path, ref, reads, [SOX2_ROW], "--rescue-homopolymer")
    assert "WARN_HOMOPOLYMER_DECOMP" in row["gbcms_status_reason"]


@pytest.mark.xfail(strict=True, reason="no OBSERVED_ALLELE diagnostic yet")
def test_the_diagnostic_names_what_the_reads_carry(tmp_path):
    ref = _ref()
    (row,) = _run(tmp_path, ref, _sox2_carriers(ref, 20) + _ref_reads(ref, 10), [SOX2_ROW])
    # CCCCCC(T) -> CCCCT(T): trimmed, the change is CC>T at 1-based 305.
    assert "OBSERVED_ALLELE(1:305:CC>T:20/0)" in row["gbcms_diagnostic"].split(";")


def test_a_correct_allele_gets_no_flag(tmp_path):
    ref = _ref()
    hap = ref[:RUN] + "T" + ref[RUN + RUN_LEN :]
    carriers = [
        make_read(
            f"a{i}", hap[s : s + READ], s, ((0, RUN - s + 1), (2, 5), (0, READ - (RUN - s + 1)))
        )
        for i, s in enumerate(range(240, 260))
    ]
    (row,) = _run(tmp_path, ref, carriers + _ref_reads(ref, 10), [SOX2_ROW])
    assert "OBSERVED_ALLELE" not in row["gbcms_diagnostic"]


def test_scattered_errors_do_not_name_an_allele(tmp_path):
    """Twelve reads, each with a different high-quality error inside the run
    (6 positions x 2 bases): no sequence recurs, so there is no allele to name."""
    ref = _ref()
    reads = []
    for i in range(12):
        s = 250 + i
        seq = list(ref[s : s + READ])
        seq[RUN - s + i % RUN_LEN] = "AG"[i // RUN_LEN]
        reads.append(make_read(f"e{i}", "".join(seq), s, ((0, READ),)))
    (row,) = _run(tmp_path, ref, reads + _ref_reads(ref, 10, start=230), [SOX2_ROW])
    assert "OBSERVED_ALLELE" not in row["gbcms_diagnostic"]


@pytest.mark.xfail(strict=True, reason="no OBSERVED_ALLELE diagnostic yet")
def test_a_mis_described_snv_is_named(tmp_path):
    """Given C>T at the run's last base; the reads carry C>G."""
    ref = _ref()
    hap = ref[: RUN + 5] + "G" + ref[RUN + 6 :]
    carriers = [
        make_read(f"g{i}", hap[s : s + READ], s, ((0, READ),))
        for i, s in enumerate(range(240, 255))
    ]
    (row,) = _run(tmp_path, ref, carriers + _ref_reads(ref, 10), [(RUN + 6, "C", "T")])
    assert "OBSERVED_ALLELE(1:306:C>G:15/0)" in row["gbcms_diagnostic"].split(";")


def test_a_germline_snp_beside_the_event_is_not_named(tmp_path):
    """REF reads carrying a SNP three bases past the event: outside the event's
    core, so it is not an allele of this row."""
    ref = _ref()
    p = RUN + RUN_LEN + 3
    alt_base = "A" if ref[p] != "A" else "G"
    hap = ref[:p] + alt_base + ref[p + 1 :]
    reads = [
        make_read(f"s{i}", hap[s : s + READ], s, ((0, READ),))
        for i, s in enumerate(range(240, 260))
    ]
    (row,) = _run(tmp_path, ref, reads + _ref_reads(ref, 10), [SOX2_ROW])
    assert "OBSERVED_ALLELE" not in row["gbcms_diagnostic"]
