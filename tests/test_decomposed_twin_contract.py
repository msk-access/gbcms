"""Target contract: the homopolymer-decomposed twin in RNA and across samples.

Where a called delins looks like a miscollapsed event in a homopolymer run
(REF a run of one base, ALT the single base that follows the run), the engine
also counts a corrected allele — the run with its last base replaced by that
following base (``CCCCCC>T`` is dual-counted as ``CCCCCC>CCCCCT``) — and
reports whichever form has more ALT support, flagging ``WARN_HOMOPOLYMER_DECOMP``.

- RNA + GTF: the twin takes its original's gene strand, so
  ``--enforce-strandedness`` applies to it. It was counted without one, so
  where it won, antisense carriers were counted.
- ``WARN_HOMOPOLYMER_DECOMP`` is per sample. The prepared variants are shared
  by every sample of a run, and the flag set for one sample stayed on every
  later sample's row.
- Guard: the legacy parity oracle and the binned engine count the same twin.

Committed red (xfail-strict) before the implementation; flipped green with it.
"""

import glob
import random

import pysam
from helpers import count_both, make_read, read_maf_output
from rna_fixtures import E1, SENSE, mk_ref, run_rna, write_bam, write_fasta, write_gtf, write_vcf
from typer.testing import CliRunner

from gbcms import _rs
from gbcms.cli import app

RUN_LEN = 6  # homopolymer run: C * RUN_LEN, then T
READ_LEN = 100
DNA_RUN = 300  # 0-based run start in the DNA fixture


def _dna_ref():
    # A/G-only flanks keep the C run and its following T unambiguous.
    rng = random.Random(107)
    left = "".join(rng.choice("AG") for _ in range(DNA_RUN))
    right = "".join(rng.choice("AG") for _ in range(300))
    return left + "C" * RUN_LEN + "T" + right


def _twin_carrier(name, ref, run_at, start, flag=0):
    """A read of the corrected allele: the run's last base reads T."""
    hap = ref[: run_at + RUN_LEN - 1] + "T" + ref[run_at + RUN_LEN :]
    return make_read(name, hap[start : start + READ_LEN], start, ((0, READ_LEN),), flag=flag)


def _ref_read(name, ref, start, flag=0):
    return make_read(name, ref[start : start + READ_LEN], start, ((0, READ_LEN),), flag=flag)


def _bam(path, ref, reads):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": len(ref)}]}
    raw = path.with_suffix(".raw.bam")
    with pysam.AlignmentFile(raw, "wb", header=header) as fh:
        for a in sorted(reads, key=lambda a: a.reference_start):
            fh.write(a)
    pysam.sort("-o", str(path), str(raw))
    pysam.index(str(path))
    return str(path)


def _fasta(tmp_path, ref):
    fa = tmp_path / "ref.fa"
    fa.write_text(f">1\n{ref}\n")
    pysam.faidx(str(fa))
    return fa


def _invariants(c):
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev


# ── parity guard ──────────────────────────────────────────────────────────
def test_twin_parity_between_engines(tmp_path):
    """Guard: where the twin wins, the legacy oracle and the binned engine
    report the same dual-count (count_both asserts field-by-field parity)."""
    ref = _dna_ref()
    fa = _fasta(tmp_path, ref)
    variant = _rs.Variant("1", DNA_RUN, "C" * RUN_LEN, "T", "COMPLEX")
    (prepared,) = _rs.prepare_variants([variant], str(fa), 5, False, 1, True)
    assert prepared.decomposed_variant.alt_allele == "C" * (RUN_LEN - 1) + "T"
    reads = [_twin_carrier(f"t{i}", ref, DNA_RUN, DNA_RUN - 60 + i) for i in range(8)]
    reads += [_ref_read(f"r{i}", ref, DNA_RUN - 50 + i) for i in range(6)]
    (counts,) = count_both(
        _bam(tmp_path / "r.bam", ref, reads),
        [prepared.variant],
        decomposed=[prepared.decomposed_variant],
    )
    _invariants(counts)
    assert counts.used_decomposed and counts.ad == 8


# ── WARN_HOMOPOLYMER_DECOMP is per sample ─────────────────────────────────
def test_decomp_flag_is_per_sample(tmp_path):
    """Sample A's reads carry the corrected allele (its twin wins); sample B's
    are reference (its original stands). A is counted first, so a flag left on
    the shared prepared variant would land on B's row."""
    ref = _dna_ref()
    fa = _fasta(tmp_path, ref)
    bam_a = _bam(
        tmp_path / "a.bam",
        ref,
        [_twin_carrier(f"t{i}", ref, DNA_RUN, DNA_RUN - 60 + i) for i in range(8)]
        + [_ref_read(f"r{i}", ref, DNA_RUN - 50 + i) for i in range(6)],
    )
    bam_b = _bam(
        tmp_path / "b.bam", ref, [_ref_read(f"r{i}", ref, DNA_RUN - 50 + i) for i in range(10)]
    )
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"1\t{DNA_RUN + 1}\t.\t{'C' * RUN_LEN}\tT\t.\t.\t.\n"
    )
    out = tmp_path / "out"
    out.mkdir()
    result = CliRunner().invoke(
        app,
        [
            "dna",
            "-v",
            str(vcf),
            "-b",
            f"A:{bam_a}",
            "-b",
            f"B:{bam_b}",
            "-f",
            str(fa),
            "-o",
            str(out),
            "--format",
            "maf",
        ],
    )
    assert result.exit_code == 0, result.output
    reasons = {}
    for sample in ("A", "B"):
        (path,) = [
            p for p in glob.glob(str(out / "*.maf")) if p.rsplit("/", 1)[1].startswith(sample)
        ]
        (row,) = list(read_maf_output(path))
        reasons[sample] = row["gbcms_status_reason"]
    assert "WARN_HOMOPOLYMER_DECOMP" in reasons["A"], reasons
    assert "WARN_HOMOPOLYMER_DECOMP" not in reasons["B"], reasons


# ── gene strand (RNA + GTF, --enforce-strandedness default) ────────────────
RNA_RUN = 200  # 0-based, inside E1 [100, 300), far from its edges


def _rna_ref():
    ref = mk_ref()
    return ref[: RNA_RUN - 1] + "A" + "C" * RUN_LEN + "T" + ref[RNA_RUN + RUN_LEN + 1 :]


def test_twin_respects_strandedness(tmp_path):
    """8 sense and 5 antisense reads carry the corrected allele and the twin
    wins. Under --enforce-strandedness (RNA default) only the sense carriers
    count, as for any variant."""
    ref = _rna_ref()
    assert E1[0] + 25 < RNA_RUN < E1[1] - RUN_LEN - 25
    sense = [_twin_carrier(f"s{i}", ref, RNA_RUN, RNA_RUN - 60 + i, flag=SENSE) for i in range(8)]
    anti = [_twin_carrier(f"a{i}", ref, RNA_RUN, RNA_RUN - 55 + i, flag=0) for i in range(5)]
    refs = [_ref_read(f"r{i}", ref, RNA_RUN - 50 + i, flag=SENSE) for i in range(10)]
    (row,) = run_rna(
        tmp_path,
        write_vcf(tmp_path, [(RNA_RUN + 1, "C" * RUN_LEN, "T")]),
        write_bam(tmp_path, ref, sense + anti + refs),
        write_fasta(tmp_path, ref),
        write_gtf(tmp_path),
    )
    assert "WARN_HOMOPOLYMER_DECOMP" in row["gbcms_status_reason"], "twin must win"
    assert int(row["alt_count"]) == 8
    assert int(row["total_count"]) == 18, "antisense reads are not depth either"
