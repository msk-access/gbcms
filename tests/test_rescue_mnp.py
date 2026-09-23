"""Tests for the MNP rescue pass (--rescue-mnp).

Rescue exists for sign-out MNPs whose carriers hold only a component of the
annotated haplotype: the engine correctly reports the haplotype as absent
(carriers land in ``partial_alt``), and rescue reports the best-supported
component instead. The end-to-end battery drives the real CLI on a synthetic
TERT-shaped geometry (GAGGG>AAGGA, discriminating at block positions 0 and 4).

Covers:
1.  Default config: rescue_mnp is False
2.  MAF column / VCF GR header present only with the flag
3.  Column count with rescue_mnp=True
4.  Component carriers + one masked stray full read → rescued; the row is the
    winning component's own counts (counting invariants hold), the audit keeps
    the MNP forensics, diagnostics describe the written row
5.  Cis carriers (haplotype dominates) → not a candidate
6.  MNP in a co-annotated group → skipped, exclusive assignment untouched
7.  Two BAMs in one run → a sample's audit never leaks into the next sample
8.  Flag off → MNP counts as the engine produced them, no rescue column
9.  Indel-disrupted partial evidence (no component carriers) → no_improvement,
    counts and diagnostics untouched
10. Outcome resolution: tie-break, no_improvement, and ref_validation_failed
    (a component SNV failing preparation — not reachable through a BAM)
11. Audit format
12. Confirmed haplotype (reads carrying every change, all read) → not rescued;
    error-level confirmed reads within the base-quality allowance → rescued
"""

import csv
import glob
import io
import random
import re
import types

import pysam
from helpers import make_read
from typer.testing import CliRunner

from gbcms.cli import app
from gbcms.io.output import MafWriter, VcfWriter
from gbcms.models.core import GbcmsBaseConfig
from gbcms.pipeline import _format_rescue_audit, _resolve_mnp_rescue

runner = CliRunner()

# ── Test 1: Default config ──────────────────────────────────────────────────


def test_rescue_flag_absent_by_default():
    """GbcmsBaseConfig.rescue_mnp defaults to False."""
    # Test via direct field default — no file I/O needed
    assert GbcmsBaseConfig.model_fields["rescue_mnp"].default is False


# ── Test 2: MAF column absent without flag ───────────────────────────────────


def test_rescue_column_absent_without_flag():
    """MafWriter omits gbcms_rescue column when rescue_mnp=False."""
    writer = MafWriter.__new__(MafWriter)
    writer.column_prefix = ""
    writer.mfsd = False
    writer.show_normalization = False
    writer.mode = "dna"
    writer.rescue_mnp = False

    cols = writer._gbcms_column_names()
    assert (
        "gbcms_rescue" not in cols
    ), f"gbcms_rescue should not be in columns without --rescue-mnp: {cols}"


# ── Test 3: MAF column present with flag ─────────────────────────────────────


def test_rescue_column_present_with_flag():
    """MafWriter includes gbcms_rescue column when rescue_mnp=True."""
    writer = MafWriter.__new__(MafWriter)
    writer.column_prefix = ""
    writer.mfsd = False
    writer.show_normalization = False
    writer.mode = "dna"
    writer.rescue_mnp = True

    cols = writer._gbcms_column_names()
    assert "gbcms_rescue" in cols, f"gbcms_rescue should be in columns with --rescue-mnp: {cols}"


# ── Test 4: VCF GR absent without flag ───────────────────────────────────────


def test_vcf_gr_absent_without_flag(tmp_path):
    """VcfWriter omits GR INFO header and value when rescue_mnp=False."""
    vcf_path = tmp_path / "no_rescue.vcf"
    writer = VcfWriter(vcf_path, sample_name="TUMOR", rescue_mnp=False)
    writer._write_header()
    writer.close()

    header_text = vcf_path.read_text()
    assert "ID=GR," not in header_text, "GR INFO header should not appear without --rescue-mnp"


# ── Test 5: VCF GR present with flag ─────────────────────────────────────────


def test_vcf_gr_present_with_flag(tmp_path):
    """VcfWriter includes GR INFO header when rescue_mnp=True."""
    vcf_path = tmp_path / "with_rescue.vcf"
    writer = VcfWriter(vcf_path, sample_name="TUMOR", rescue_mnp=True)
    writer._write_header()
    writer.close()

    header_text = vcf_path.read_text()
    assert "ID=GR," in header_text, "GR INFO header should appear with --rescue-mnp"


# ── Test 6: Column count with rescue ────────────────────────────────────────


def test_column_count_with_rescue():
    """With rescue_mnp=True, there should be 27 columns (26 base + gbcms_rescue)."""
    writer = MafWriter.__new__(MafWriter)
    writer.column_prefix = ""
    writer.mfsd = False
    writer.show_normalization = False
    writer.mode = "dna"
    writer.rescue_mnp = True

    cols = writer._gbcms_column_names()
    assert len(cols) == 27, f"Expected 27 gbcms MAF columns with rescue, got {len(cols)}: {cols}"
    assert cols[3] == "gbcms_rescue", f"gbcms_rescue should be the 4th column, got {cols[3]}"


# ── End-to-end battery: TERT-shaped sparse ONP ──────────────────────────────
#
# Harness conventions follow test_e2e_partial_dominant.py: FASTA/VCF use
# ``chr1``, the BAM header uses the normalized name ``1``, reads are MAPQ 60 /
# Q30 so default gates pass. Reads are unpaired, so each read is its own
# fragment. The block GAGGG sits at 0-based 200..204; REF and ALT differ only
# at block offsets 0 and 4 (1-based 201 and 205).

BAM_CONTIG = "1"
BLOCK_START = 200
REF_BLOCK = "GAGGG"
ALT_BLOCK = "AAGGA"
READ_LEN = 100
N_REF = 20


def _reference_sequence():
    """600bp seeded-random contig with the pinned GAGGG block at 200..204.

    The flanking bases are pinned to C so no flank base can extend or mimic a
    block allele.
    """
    rng = random.Random(7)
    seq = [rng.choice("ACGT") for _ in range(600)]
    seq[BLOCK_START - 1] = "C"
    seq[BLOCK_START : BLOCK_START + 5] = list(REF_BLOCK)
    seq[BLOCK_START + 5] = "C"
    return "".join(seq)


REF = _reference_sequence()


def _read(name, block, i, low_bq_block_offset=None, block_cigar=((0, 5),)):
    """A forward/reverse-alternating read carrying ``block`` (aligned as
    ``block_cigar``) in place of the reference block, with 40-44bp of left
    flank and flank to 100 reference bases on the right."""
    start = BLOCK_START - 40 - (i % 5)
    left = REF[start:BLOCK_START]
    right = REF[BLOCK_START + 5 : start + READ_LEN]
    quals = [30] * (len(left) + len(block) + len(right))
    if low_bq_block_offset is not None:
        quals[len(left) + low_bq_block_offset] = 5
    cigar = ((0, len(left)),) + tuple(block_cigar) + ((0, len(right)),)
    return make_read(name, left + block + right, start, cigar, flag=16 if i % 2 else 0, quals=quals)


def _component_carrier_reads():
    """20 REF + 9 first-change-only carriers + 1 carrier whose second
    discriminating base is BQ 5 (masked, so the engine counts it full ALT —
    the stray ``ad = 1`` that the old ``ad == 0`` gate tripped on)."""
    reads = [_read(f"ref_{i}", REF_BLOCK, i) for i in range(N_REF)]
    reads += [_read(f"comp_{i}", "AAGGG", i) for i in range(9)]
    reads.append(_read("stray_masked", "AAGGG", 9, low_bq_block_offset=4))
    return reads


def _cis_carrier_reads():
    """20 REF + 10 full-haplotype carriers."""
    reads = [_read(f"ref_{i}", REF_BLOCK, i) for i in range(N_REF)]
    reads += [_read(f"cis_{i}", ALT_BLOCK, i) for i in range(10)]
    return reads


def _germline_component_reads():
    """20 REF + 30 reads carrying only the first change (a heterozygous
    germline SNP merged into the annotation) + 10 carrying both changes (the
    real somatic MNP, every discriminating base read)."""
    reads = [_read(f"ref_{i}", REF_BLOCK, i) for i in range(N_REF)]
    reads += [_read(f"germ_{i}", "AAGGG", i) for i in range(30)]
    reads += [_read(f"cis_{i}", ALT_BLOCK, i) for i in range(10)]
    return reads


def _component_reads_with_error_level_cis():
    """20 REF + 200 first-change-only carriers + 1 fully read full-haplotype
    read — the level a sequencing error at the second position produces."""
    reads = [_read(f"ref_{i}", REF_BLOCK, i) for i in range(N_REF)]
    reads += [_read(f"comp_{i}", "AAGGG", i) for i in range(200)]
    reads.append(_read("err_cis", ALT_BLOCK, 0))
    return reads


def _indel_disrupted_reads():
    """20 REF + 10 reads carrying a 1bp insertion inside the block and REF
    bases otherwise: no component of the MNP, but the complex path counts them
    REF with nearby-indel evidence, so partial_alt dominates ad."""
    reads = [_read(f"ref_{i}", REF_BLOCK, i) for i in range(N_REF)]
    reads += [
        _read(f"ins_{i}", "GAGTGG", i, block_cigar=((0, 3), (1, 1), (0, 2))) for i in range(10)
    ]
    return reads


def _write_bam(tmp_path, name, reads):
    unsorted = tmp_path / f"{name}.unsorted.bam"
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 600, "SN": BAM_CONTIG}]}
    with pysam.AlignmentFile(unsorted, "wb", header=header) as outf:
        for r in sorted(reads, key=lambda a: a.reference_start):
            outf.write(r)
    sorted_bam = tmp_path / f"{name}.bam"
    pysam.sort("-o", str(sorted_bam), str(unsorted))
    pysam.index(str(sorted_bam))
    return sorted_bam


MNP_ROW = ("chr1", BLOCK_START + 1, REF_BLOCK, ALT_BLOCK)
SNV_AT_BLOCK_START = ("chr1", BLOCK_START + 1, "G", "A")


def _run(tmp_path, bams, vcf_rows, rescue=True):
    """Run ``gbcms dna`` on named BAMs; return {sample: [MAF rows]} in VCF order."""
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">chr1\n" + REF + "\n")
    pysam.faidx(str(fasta))
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=600>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(f"{c}\t{p}\t.\t{r}\t{a}\t.\t.\t.\n" for c, p, r, a in vcf_rows)
    )
    outdir = tmp_path / ("out_rescue" if rescue else "out_plain")
    args = ["dna", "-v", str(vcf), "-f", str(fasta), "-o", str(outdir), "--format", "maf"]
    for name, reads in bams.items():
        args += ["-b", f"{name}:{_write_bam(tmp_path, name, reads)}"]
    if rescue:
        args.append("--rescue-mnp")
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output

    rows = {}
    for name in bams:
        (path,) = glob.glob(str(outdir / f"{name}.maf"))
        with open(path) as f:
            body = "".join(line for line in f if not line.startswith("#"))
        rows[name] = list(csv.DictReader(io.StringIO(body), delimiter="\t"))
    return rows


def _assert_counting_invariants(row):
    """The four AGENTS.md counting invariants, on the row as written."""
    rd, ad = int(row["ref_count"]), int(row["alt_count"])
    assert int(row["total_count"]) >= rd + ad
    assert int(row["total_count_fragment"]) >= int(row["ref_count_fragment"]) + int(
        row["alt_count_fragment"]
    )
    assert rd == int(row["ref_count_forward"]) + int(row["ref_count_reverse"])
    assert ad == int(row["alt_count_forward"]) + int(row["alt_count_reverse"])


def _audit(row):
    return dict(part.split("=", 1) for part in row["gbcms_rescue"].split(";"))


def test_component_carriers_rescued_as_one_coherent_genotype(tmp_path):
    rows = _run(tmp_path, {"S": _component_carrier_reads()}, [MNP_ROW])
    (row,) = rows["S"]

    # The row is the winning component (block offset 0, G>A at 201) counted
    # as an SNV: every ALT-derived column moves together.
    assert int(row["ref_count"]) == N_REF
    assert int(row["alt_count"]) == 10
    assert int(row["alt_count_fragment"]) == 10
    assert int(row["partial_alt"]) == 0
    assert int(row["any_alt"]) == 10
    _assert_counting_invariants(row)

    audit = _audit(row)
    assert audit["method"] == "decomposed"
    assert audit["outcome"] == "rescued"
    # The MNP's own evaluation survives only here: 20 REF, the one masked
    # stray counted full ALT, the nine single-change carriers as partial.
    assert (audit["original_ref"], audit["original_alt"], audit["original_partial"]) == (
        "20",
        "1",
        "9",
    )
    # The stray's other discriminating base was masked: nothing showed the haplotype.
    assert audit["original_confirmed"] == "0"
    assert re.fullmatch(r"(chr)?1:201\(G>A\)", audit["adopted"]), audit["adopted"]
    assert re.fullmatch(r"(chr)?1:201\(G>A\):10,(chr)?1:205\(G>A\):0", audit["positions"]), audit[
        "positions"
    ]

    # Diagnostics describe the written row, not the pre-rescue MNP counts.
    assert row["gbcms_diagnostic"] == "MNP_DISC_RATIO(2/5);MNP_RESCUE_ELIGIBLE"


def test_flag_off_reports_engine_counts_without_rescue_column(tmp_path):
    rows = _run(tmp_path, {"S": _component_carrier_reads()}, [MNP_ROW], rescue=False)
    (row,) = rows["S"]
    assert "gbcms_rescue" not in row
    assert (int(row["ref_count"]), int(row["alt_count"]), int(row["partial_alt"])) == (20, 1, 9)
    assert "PARTIAL_DOMINANT" in row["gbcms_diagnostic"].split(";")
    _assert_counting_invariants(row)


def test_cis_carriers_are_not_candidates(tmp_path):
    rows = _run(tmp_path, {"S": _cis_carrier_reads()}, [MNP_ROW])
    (row,) = rows["S"]
    assert row["gbcms_rescue"] == ""
    assert (int(row["ref_count"]), int(row["alt_count"]), int(row["partial_alt"])) == (20, 10, 0)
    _assert_counting_invariants(row)


def test_grouped_mnp_is_skipped_and_keeps_exclusive_assignment(tmp_path):
    variants = [MNP_ROW, SNV_AT_BLOCK_START]
    plain = _run(tmp_path, {"S": _component_carrier_reads()}, variants, rescue=False)["S"]
    rescued = _run(tmp_path, {"S": _component_carrier_reads()}, variants)["S"]

    mnp_plain, mnp = plain[0], rescued[0]
    count_cols = ("ref_count", "alt_count", "partial_alt", "alt_count_fragment")
    assert {c: mnp[c] for c in count_cols} == {c: mnp_plain[c] for c in count_cols}
    audit = _audit(mnp)
    assert audit["outcome"] == "skipped_grouped"
    assert "positions" not in audit
    assert audit["original_alt"] == mnp_plain["alt_count"]

    snv = rescued[1]
    assert snv["gbcms_rescue"] == ""
    assert {c: snv[c] for c in count_cols} == {c: plain[1][c] for c in count_cols}


def test_indel_partial_evidence_is_declined_not_rescued(tmp_path):
    plain = _run(tmp_path, {"S": _indel_disrupted_reads()}, [MNP_ROW], rescue=False)["S"][0]
    (row,) = _run(tmp_path, {"S": _indel_disrupted_reads()}, [MNP_ROW])["S"]

    assert int(row["partial_alt"]) > int(row["alt_count"])  # a candidate
    count_cols = ("ref_count", "alt_count", "partial_alt", "alt_count_fragment", "total_count")
    assert {c: row[c] for c in count_cols} == {c: plain[c] for c in count_cols}
    assert row["gbcms_diagnostic"] == plain["gbcms_diagnostic"]
    audit = _audit(row)
    assert audit["outcome"] == "no_improvement"
    assert "adopted" not in audit
    assert re.fullmatch(r"(chr)?1:201\(G>A\):0,(chr)?1:205\(G>A\):0", audit["positions"])
    _assert_counting_invariants(row)


def test_confirmed_haplotype_is_not_rescued(tmp_path):
    plain = _run(tmp_path, {"S": _germline_component_reads()}, [MNP_ROW], rescue=False)["S"][0]
    (row,) = _run(tmp_path, {"S": _germline_component_reads()}, [MNP_ROW])["S"]

    assert (int(row["alt_count"]), int(row["partial_alt"])) == (10, 30)  # a partial-dominant row
    count_cols = ("ref_count", "alt_count", "partial_alt", "alt_count_fragment", "total_count")
    assert {c: row[c] for c in count_cols} == {c: plain[c] for c in count_cols}
    assert row["gbcms_diagnostic"] == plain["gbcms_diagnostic"]
    audit = _audit(row)
    assert audit["outcome"] == "haplotype_confirmed"
    assert audit["original_confirmed"] == "10"
    assert "positions" not in audit  # no component re-count was run
    _assert_counting_invariants(row)


def test_error_level_confirmed_reads_do_not_block_rescue(tmp_path):
    (row,) = _run(tmp_path, {"S": _component_reads_with_error_level_cis()}, [MNP_ROW])["S"]
    audit = _audit(row)
    assert audit["outcome"] == "rescued"
    assert audit["original_confirmed"] == "1"  # within ceil(200 x 1%) = 2
    assert int(row["alt_count"]) == 201
    _assert_counting_invariants(row)


def test_rescue_audit_does_not_leak_into_later_samples(tmp_path):
    rows = _run(tmp_path, {"A": _component_carrier_reads(), "B": _cis_carrier_reads()}, [MNP_ROW])
    assert _audit(rows["A"][0])["outcome"] == "rescued"
    assert rows["B"][0]["gbcms_rescue"] == ""
    assert int(rows["B"][0]["alt_count"]) == 10


# ── Outcome resolution and audit format (pure helpers) ──────────────────────


def test_resolve_adopts_best_component_leftmost_on_ties():
    assert _resolve_mnp_rescue(1, [10, 0]) == ("rescued", 0)
    assert _resolve_mnp_rescue(1, [3, 7]) == ("rescued", 1)
    assert _resolve_mnp_rescue(0, [5, 5]) == ("rescued", 0)
    assert _resolve_mnp_rescue(0, [None, 4]) == ("rescued", 1)


def test_resolve_never_adopts_a_component_that_does_not_beat_the_haplotype():
    assert _resolve_mnp_rescue(8, [4, 8]) == ("no_improvement", None)
    assert _resolve_mnp_rescue(0, [0, None]) == ("no_improvement", None)


def test_resolve_reports_when_no_component_could_be_counted():
    assert _resolve_mnp_rescue(0, [None, None]) == ("ref_validation_failed", None)


def test_audit_format():
    original = types.SimpleNamespace(rd=486, ad=1, partial_alt=88, mnp_confirmed_alt=0)
    assert _format_rescue_audit("skipped_grouped", original) == (
        "method=decomposed;outcome=skipped_grouped;"
        "original_ref=486;original_alt=1;original_partial=88;original_confirmed=0"
    )
    assert _format_rescue_audit(
        "rescued",
        original,
        ["5:1295250(G>A):87", "5:1295254(G>A):ref_fail"],
        "5:1295250(G>A)",
    ) == (
        "method=decomposed;outcome=rescued;original_ref=486;original_alt=1;"
        "original_partial=88;original_confirmed=0;adopted=5:1295250(G>A);"
        "positions=5:1295250(G>A):87,5:1295254(G>A):ref_fail"
    )


def test_confirmed_error_allowance_follows_the_base_quality_threshold():
    from gbcms.pipeline import _confirmed_error_allowance

    assert _confirmed_error_allowance(850, 20) == 9  # 8.5 at Q20 (1% error)
    assert _confirmed_error_allowance(238, 20) == 3
    assert _confirmed_error_allowance(200, 20) == 2
    assert _confirmed_error_allowance(100, 30) == 1  # 0.1 at Q30, rounded up
    assert _confirmed_error_allowance(0, 20) == 0
    assert _confirmed_error_allowance(100, 0) == 100  # no quality gate: any base may be wrong
