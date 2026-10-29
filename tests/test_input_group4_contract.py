"""Group 4 (input and representation) contracts.

Operator decisions (2026-09-25 for I1, I3, I4; 2026-10-05 for the rest):
- I1 #123: a MAF row with an allele that is not a base sequence (an IUPAC code,
  a stray character) is a visible FAIL row, NON_SEQUENCE_ALLELE, kept in MAF
  output; lowercase bases are bases. VCF output, which cannot carry such an
  allele, writes a symbolic record: REF the reference base at POS, ALT
  <NON_SEQUENCE> (declared in the header).
- MAF input, VCF output: every record carries the MAF row it came from
  (MAF_START, MAF_REF, MAF_ALT), as VCF input's MAF output carries vcf_pos,
  vcf_ref and vcf_alt, so a result can be looked up by its input.
- C9 #122: a MAF deletion at Start_Position 1 (no base before it) is resolved to
  the VCF spec's form at position 1, the base after the event, and counted, as
  the same event given as VCF already is.
- I7 #219: a MAF dash insertion at Start_Position 0 (before the contig's first base) is
  kept: resolved to the VCF spec's form at position 1, the base after the event, and
  counted, as the C9 deletion at Start 1 is. It was dropped before it reached
  preparation, so the output had no row for it. Any other row at Start 0 (a
  substitution, a deletion) names a base that does not exist: a visible FAIL row,
  START_BEFORE_CONTIG, not a dropped one.
- I2 #124: End_Position is optional. gbcms merge joins on contig, Start and
  alleles, not End_Position, and fills it from the inputs that have each row.
- I3 #125: VCF input's MAF output fills Tumor_Seq_Allele1 with the reference
  allele (MSK's convention on every sign-out row; maf2vcf reads an empty one as
  the reference).
- I4 #126: a Tumor_Seq_Allele1 differing from both REF and Allele2 is not a
  second allele: one row, one allele (vcf2maf's reading).
- #147: a decomposed or sibling list shorter than the variants is padded; a
  longer one is a caller error.
"""

import csv
import glob
import logging
import random

import pysam
import pytest
from helpers import make_read, read_maf_output
from pydantic import ValidationError
from typer.testing import CliRunner

from gbcms import _rs
from gbcms.cli import app
from gbcms.core.kernel import CoordinateKernel
from gbcms.io.input import MafReader
from gbcms.io.output import MafWriter
from gbcms.models.core import Variant, VariantType

runner = CliRunner()
_RNG = random.Random(404)
REF = "".join(_RNG.choice("ACGT") for _ in range(600))
P = 300  # 0-based site used by most rows
READ = 100


def _files(tmp_path, reads):
    fa = tmp_path / "ref.fa"
    fa.write_text(f">1\n{REF}\n")
    pysam.faidx(str(fa))
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": len(REF)}]}
    raw = tmp_path / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=hdr) as fh:
        for r in reads:
            fh.write(r)
    bam = tmp_path / "s.bam"
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    return fa, bam


def _ref_reads(start=P - 50, n=6):
    return [
        make_read(f"r{i}", REF[start + i : start + i + READ], start + i, ((0, READ),))
        for i in range(n)
    ]


def _maf(tmp_path, rows, header=None):
    head = header or [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Tumor_Sample_Barcode",
    ]
    lines = ["\t".join(head)] + ["\t".join(str(row.get(c, "")) for c in head) for row in rows]
    path = tmp_path / "v.maf"
    path.write_text("\n".join(lines) + "\n")
    return path


def _row(start, ref, alt, end=None, **extra):
    end = end if end is not None else start + max(len(ref), 1) - 1
    return {
        "Hugo_Symbol": "G",
        "Chromosome": "1",
        "Start_Position": start,
        "End_Position": end,
        "Reference_Allele": ref,
        "Tumor_Seq_Allele2": alt,
        "Tumor_Sample_Barcode": "S",
        **extra,
    }


def _run(tmp_path, variants, reads, fmt="maf", name="out"):
    fa, bam = _files(tmp_path, reads)
    out = tmp_path / name
    res = runner.invoke(
        app,
        [
            "dna",
            "-v",
            str(variants),
            "-b",
            f"S:{bam}",
            "-f",
            str(fa),
            "-o",
            str(out),
            "--format",
            fmt,
        ],
    )
    assert res.exit_code == 0, res.output
    return glob.glob(str(out / f"*.{fmt}"))[0]


def _maf_rows(path):
    return list(read_maf_output(path))


def _vcf_records(path):
    with pysam.VariantFile(path) as fh:
        return list(fh)


# ── I1: alleles that are not base sequences ──────────────────────────────────

_NON_SEQUENCE = [
    pytest.param(_row(P + 1, REF[P], "R"), id="IUPAC ALT"),
    pytest.param(_row(P + 1, REF[P], "LU"), id="letters ALT"),
    pytest.param(_row(P + 1, REF[P], "."), id="dot ALT"),
    pytest.param(_row(P + 1, REF[P], "T[N]{?}"), id="brackets ALT"),
    pytest.param(_row(P + 1, REF[P : P + 2] + "V", REF[P]), id="IUPAC REF"),
]


@pytest.mark.parametrize("row", _NON_SEQUENCE)
def test_a_non_sequence_maf_allele_is_a_fail_row(tmp_path, row):
    """The row is kept in MAF output, FAIL with NON_SEQUENCE_ALLELE, not counted."""
    good = _row(P + 21, REF[P + 20], next(b for b in "ACGT" if b != REF[P + 20]))
    out = _maf_rows(_run(tmp_path, _maf(tmp_path, [row, good]), _ref_reads()))
    assert len(out) == 2
    bad = out[0]
    assert (bad["gbcms_status"], bad["gbcms_status_reason"]) == ("FAIL", "NON_SEQUENCE_ALLELE")
    assert (bad["Reference_Allele"], bad["Tumor_Seq_Allele2"]) == (
        row["Reference_Allele"],
        row["Tumor_Seq_Allele2"],
    )
    assert out[1]["gbcms_status"] == "PASS"


def test_lowercase_maf_alleles_are_bases(tmp_path):
    """Guard: lowercase alleles are bases: the row counts as its uppercase twin."""
    alt = next(b for b in "ACGT" if b != REF[P])
    reads = _ref_reads()
    hap = REF[:P] + alt + REF[P + 1 :]
    reads += [
        make_read(f"a{i}", hap[P - 40 + i : P - 40 + i + READ], P - 40 + i, ((0, READ),))
        for i in range(5)
    ]
    upper = _maf_rows(_run(tmp_path, _maf(tmp_path, [_row(P + 1, REF[P], alt)]), reads, name="u"))[
        0
    ]
    (tmp_path / "l").mkdir()
    lower = _maf_rows(
        _run(
            tmp_path / "l", _maf(tmp_path / "l", [_row(P + 1, REF[P].lower(), alt.lower())]), reads
        )
    )[0]
    assert lower["gbcms_status"] == "PASS"
    assert (
        (lower["ref_count"], lower["alt_count"])
        == (upper["ref_count"], upper["alt_count"])
        == ("6", "5")
    )


def test_vcf_output_writes_a_non_sequence_row_as_a_symbolic_record(tmp_path):
    """VCF output stays valid: REF is the reference base at POS, ALT is the
    declared symbolic <NON_SEQUENCE>, GS=FAIL, GSR=NON_SEQUENCE_ALLELE, and the
    MAF row it came from is in MAF_START/MAF_REF/MAF_ALT."""
    path = _run(tmp_path, _maf(tmp_path, [_row(P + 1, REF[P], "LU")]), _ref_reads(), fmt="vcf")
    with pysam.VariantFile(path) as fh:
        assert "NON_SEQUENCE" in fh.header.alts
        (rec,) = list(fh)
    assert (rec.pos, rec.ref, rec.alts) == (P + 1, REF[P], ("<NON_SEQUENCE>",))
    assert (rec.info["GS"], rec.info["GSR"]) == ("FAIL", "NON_SEQUENCE_ALLELE")
    assert (rec.info["MAF_START"], rec.info["MAF_REF"], rec.info["MAF_ALT"]) == (
        P + 1,
        REF[P],
        "LU",
    )


# ── MAF origin in VCF output ─────────────────────────────────────────────────


def test_vcf_output_of_maf_input_carries_the_maf_row(tmp_path):
    """A MAF deletion (Start at its first deleted base) is written as maf2vcf's
    record, anchored one base before; the MAF row it came from rides along."""
    rows = [
        _row(P + 1, REF[P : P + 2], "-"),
        _row(P + 21, REF[P + 20], "A" if REF[P + 20] != "A" else "C"),
    ]
    recs = _vcf_records(_run(tmp_path, _maf(tmp_path, rows), _ref_reads(), fmt="vcf"))
    assert (recs[0].pos, recs[0].ref, recs[0].alts[0]) == (P, REF[P - 1 : P + 2], REF[P - 1])
    got = [(r.info["MAF_START"], r.info["MAF_REF"], r.info["MAF_ALT"]) for r in recs]
    assert got == [
        (row["Start_Position"], row["Reference_Allele"], row["Tumor_Seq_Allele2"]) for row in rows
    ]


def test_vcf_output_of_vcf_input_has_no_maf_row(tmp_path):
    """Guard: a VCF-input record echoes the input record and gains no MAF fields."""
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        f"##fileformat=VCFv4.2\n##contig=<ID=1,length={len(REF)}>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"1\t{P + 1}\t.\t{REF[P]}\t{'A' if REF[P] != 'A' else 'C'}\t.\t.\t.\n"
    )
    (rec,) = _vcf_records(_run(tmp_path, vcf, _ref_reads(), fmt="vcf"))
    assert "MAF_START" not in rec.info


# ── C9: a MAF deletion at Start_Position 1 ───────────────────────────────────


def test_a_maf_deletion_at_start_1_counts_in_the_base_after_form(tmp_path):
    """The contig's first two bases deleted: the MAF row keeps its input columns
    and counts the reads that show those bases as REF (an ALT molecule starts at
    base 3, so no read can show the deletion); VCF output is the base-after
    record with the MAF row in MAF_START/MAF_REF/MAF_ALT."""
    reads = [make_read(f"r{i}", REF[0:READ], 0, ((0, READ),)) for i in range(6)]
    maf = _maf(tmp_path, [_row(1, REF[0:2], "-")])
    (row,) = _maf_rows(_run(tmp_path, maf, reads))
    assert (row["Start_Position"], row["Reference_Allele"], row["Tumor_Seq_Allele2"]) == (
        "1",
        REF[0:2],
        "-",
    )
    assert row["gbcms_status"] == "PASS"
    assert (row["ref_count"], row["alt_count"], row["total_count"]) == ("6", "0", "6")
    (rec,) = _vcf_records(_run(tmp_path, maf, reads, fmt="vcf", name="vcf"))
    assert (rec.pos, rec.ref, rec.alts[0]) == (1, REF[0:3], REF[2])
    assert rec.info["GS"] == "PASS"
    assert (rec.info["MAF_START"], rec.info["MAF_REF"], rec.info["MAF_ALT"]) == (1, REF[0:2], "-")


# ── I7: a MAF dash insertion at Start_Position 0 ─────────────────────────────

INS = "GT"  # the inserted bases


def _start_0_reads():
    """Six REF reads from the contig's first base; three reads that open with the inserted
    bases as an insertion and three that open with them as a soft clip (no base precedes
    them); and REF reads over the site the valid row uses."""
    rest = READ - len(INS)
    return (
        [make_read(f"r{i}", REF[:READ], 0, ((0, READ),)) for i in range(6)]
        + [make_read(f"i{i}", INS + REF[:rest], 0, ((1, len(INS)), (0, rest))) for i in range(3)]
        + [make_read(f"s{i}", INS + REF[:rest], 0, ((4, len(INS)), (0, rest))) for i in range(3)]
        + _ref_reads()
    )


def _good_row():
    return _row(P + 1, REF[P], "A" if REF[P] != "A" else "C")


def _vcf_form_of_the_insertion(tmp_path):
    """The same event given as VCF: the base after the event at POS 1 (REF the first base,
    ALT the inserted bases then that base), and the valid row."""
    g = _good_row()
    vcf = tmp_path / "twin.vcf"
    vcf.write_text(
        f"##fileformat=VCFv4.2\n##contig=<ID=1,length={len(REF)}>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"1\t1\t.\t{REF[0]}\t{INS + REF[0]}\t.\t.\t.\n"
        f"1\t{P + 1}\t.\t{g['Reference_Allele']}\t{g['Tumor_Seq_Allele2']}\t.\t.\t.\n"
    )
    return vcf


def _counts(row):
    return tuple(row[c] for c in ("ref_count", "alt_count", "total_count", "partial_alt"))


def test_a_maf_insertion_at_start_0_is_kept_and_counts_as_its_vcf_form(tmp_path):
    """A MAF insertion's Start is the base before it, so Start 0 puts it before the contig's
    first base. The row was dropped before preparation (a two-row MAF gave one row); it is now
    kept with its input columns and counts as the same event given as VCF does (the base
    after it, at position 1): the REF reads, and the reads that open with the inserted bases."""
    reads = _start_0_reads()
    out = _maf_rows(_run(tmp_path, _maf(tmp_path, [_row(0, "-", INS), _good_row()]), reads))
    assert len(out) == 2
    first, second = out
    assert (first["Start_Position"], first["Reference_Allele"], first["Tumor_Seq_Allele2"]) == (
        "0",
        "-",
        INS,
    )
    assert (first["gbcms_status"], first["gbcms_status_reason"]) == ("PASS", "")
    assert second["Start_Position"] == str(P + 1) and second["gbcms_status"] == "PASS"

    (tmp_path / "twin").mkdir()
    twin = _maf_rows(_run(tmp_path / "twin", _vcf_form_of_the_insertion(tmp_path), reads))
    assert [_counts(r) for r in out] == [_counts(r) for r in twin]
    assert first["ref_count"] == "6" and int(first["alt_count"]) > 0


def test_vcf_output_of_a_maf_insertion_at_start_0_is_the_base_after_form(tmp_path):
    """POS 1, REF the first base, ALT the inserted bases then the first base; the MAF row it
    came from (Start 0, '-') rides along in MAF_START/MAF_REF/MAF_ALT."""
    maf = _maf(tmp_path, [_row(0, "-", INS), _good_row()])
    recs = _vcf_records(_run(tmp_path, maf, _start_0_reads(), fmt="vcf"))
    assert len(recs) == 2
    first = recs[0]
    assert (first.pos, first.ref, first.alts[0]) == (1, REF[0], INS + REF[0])
    assert first.info["GS"] == "PASS"
    assert (first.info["MAF_START"], first.info["MAF_REF"], first.info["MAF_ALT"]) == (0, "-", INS)


def test_convert_writes_a_maf_insertion_at_start_0_in_the_base_after_form(tmp_path):
    """gbcms convert writes the same record, and keeps the row (it dropped it)."""
    fa, _ = _files(tmp_path, [])
    out = tmp_path / "o.vcf"
    maf = _maf(tmp_path, [_row(0, "-", INS), _good_row()])
    res = runner.invoke(app, ["convert", "-v", str(maf), "-f", str(fa), "-o", str(out)])
    assert res.exit_code == 0, res.output
    recs = _vcf_records(str(out))
    assert [(r.pos, r.ref, r.alts[0]) for r in recs][0] == (1, REF[0], INS + REF[0])
    assert len(recs) == 2


def test_maf_to_vcf_for_an_insertion_at_start_0_uses_the_base_after():
    """The base after an insertion at Start 0 is position 1's (an insertion at Start s lies
    between s and s + 1); POS is 1, as for an event at position 1."""
    asked = []

    def base_at(pos):
        asked.append(pos)
        return {1: "C", 2: "G"}[pos]

    assert CoordinateKernel.maf_to_vcf(0, "-", "TT", base_at) == (1, "C", "TTC")
    assert asked == [1]
    # An insertion after position 1 is anchored on it, as before.
    assert CoordinateKernel.maf_to_vcf(1, "-", "TT", base_at) == (1, "C", "CTT")


def test_the_variant_model_allows_position_minus_1_and_nothing_below():
    """Internally a MAF row is placed by Start - 1. Start 0 is -1: before the contig's first
    base, where preparation resolves an insertion and fails any other row. Lower is no
    position at all."""
    for ref, alt in [("-", INS), ("A", "C"), ("AC", "-"), ("A", INS)]:
        v = Variant(chrom="1", pos=-1, ref=ref, alt=alt, variant_type=VariantType.COMPLEX)
        assert v.pos == -1
        with pytest.raises(ValidationError):
            Variant(chrom="1", pos=-2, ref=ref, alt=alt, variant_type=VariantType.COMPLEX)


def test_a_vcf_pos_of_0_is_still_not_a_position():
    """Guard: only a MAF row can start before the first base. A VCF record at POS 0 is an
    error, as before (the readers do not turn it into a row)."""
    with pytest.raises(ValueError):
        CoordinateKernel.vcf_to_internal("1", 0, "A", "C")


def test_maf_reader_reads_every_row_at_start_0_and_skips_a_negative_start(tmp_path, caplog):
    """Start 0 is a position the reader keeps for preparation to judge (an insertion is
    resolved, any other row fails visibly). A negative Start is no position: skipped with a
    warning, as before."""
    rows = [_row(0, "-", INS), _row(0, REF[0], "T"), _row(0, REF[0], "-"), _row(-3, REF[0], "T")]
    with caplog.at_level(logging.WARNING):
        got = [(v.pos, v.ref, v.alt) for v in MafReader(_maf(tmp_path, rows))]
    assert got == [(-1, "-", INS), (-1, REF[0], "T"), (-1, REF[0], "-")]
    assert sum("Skipped MAF row" in m for m in caplog.messages) == 1


def _plain_fasta(tmp_path, seq):
    fa = tmp_path / "p.fa"
    fa.write_text(f">1\n{seq}\n")
    pysam.faidx(str(fa))
    return str(fa)


def test_prepare_resolves_an_insertion_before_the_first_base(tmp_path):
    """The prepared variant is the base-after form at position 1 (0-based 0)."""
    fa = _plain_fasta(tmp_path, REF[:60])
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", -1, "-", INS, "INSERTION")], fa, 5, True, 1, True
    )
    assert (pv.gbcms_status, pv.gbcms_status_reason) == ("PASS", "")
    assert (pv.variant.pos, pv.variant.ref_allele, pv.variant.alt_allele) == (
        0,
        REF[0],
        INS + REF[0],
    )
    assert pv.was_anchor_resolved and pv.original_pos == -1


def test_an_insertion_before_a_first_base_that_is_n_is_a_visible_fail_row(tmp_path):
    """Human contigs open with N: the base-after form has an N in its ALT, which no read
    can carry. The row is kept as FAIL ALT_CONTAINS_N, not dropped."""
    fa = _plain_fasta(tmp_path, "N" + REF[1:60])
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", -1, "-", INS, "INSERTION")], fa, 5, True, 1, True
    )
    assert (pv.gbcms_status, pv.gbcms_status_reason) == ("FAIL", "ALT_CONTAINS_N")


_START_0_FAILS = [
    pytest.param(lambda: _row(0, REF[0], "T" if REF[0] != "T" else "C"), id="SNV"),
    pytest.param(lambda: _row(0, REF[0], "-"), id="deletion"),
    pytest.param(lambda: _row(0, REF[0:2], REF[5]), id="delins"),
    pytest.param(lambda: _row(0, REF[0], INS + REF[0]), id="insertion with sequence alleles"),
]


@pytest.mark.parametrize("make_row", _START_0_FAILS)
def test_a_maf_row_at_start_0_that_is_not_an_insertion_is_a_visible_fail_row(tmp_path, make_row):
    """Position 0 is no base, so only an insertion can start there. Any other row was
    dropped with a warning (a two-row MAF gave one row); it is now kept as FAIL with
    START_BEFORE_CONTIG, its input columns as written, not counted. The valid row beside
    it counts as usual."""
    row = make_row()
    out = _maf_rows(_run(tmp_path, _maf(tmp_path, [row, _good_row()]), _start_0_reads()))
    assert len(out) == 2
    bad, good = out
    assert (bad["gbcms_status"], bad["gbcms_status_reason"]) == ("FAIL", "START_BEFORE_CONTIG")
    assert (bad["Start_Position"], bad["Reference_Allele"], bad["Tumor_Seq_Allele2"]) == (
        "0",
        row["Reference_Allele"],
        row["Tumor_Seq_Allele2"],
    )
    assert (bad["ref_count"], bad["alt_count"], bad["total_count"]) == ("0", "0", "0")
    assert (good["gbcms_status"], good["Start_Position"]) == ("PASS", str(P + 1))
    assert good["ref_count"] == "6"


def test_vcf_output_writes_a_start_0_fail_row_as_a_symbolic_record(tmp_path):
    """VCF cannot place a row at position 0: the record is the declared symbolic
    <NON_SEQUENCE> at POS 1 (REF the first base), GS=FAIL, GSR=START_BEFORE_CONTIG, and
    the MAF row it came from in MAF_START/MAF_REF/MAF_ALT."""
    row = _row(0, REF[0], "-")
    recs = _vcf_records(
        _run(tmp_path, _maf(tmp_path, [row, _good_row()]), _start_0_reads(), fmt="vcf")
    )
    assert len(recs) == 2
    first = recs[0]
    assert (first.pos, first.ref, first.alts) == (1, REF[0], ("<NON_SEQUENCE>",))
    assert (first.info["GS"], first.info["GSR"]) == ("FAIL", "START_BEFORE_CONTIG")
    assert (first.info["MAF_START"], first.info["MAF_REF"], first.info["MAF_ALT"]) == (
        0,
        REF[0],
        "-",
    )


def test_convert_writes_a_start_0_row_as_a_symbolic_record_at_position_1(tmp_path):
    fa, _ = _files(tmp_path, [])
    out = tmp_path / "o.vcf"
    maf = _maf(tmp_path, [_row(0, REF[0], "T" if REF[0] != "T" else "C"), _good_row()])
    res = runner.invoke(app, ["convert", "-v", str(maf), "-f", str(fa), "-o", str(out)])
    assert res.exit_code == 0, res.output
    recs = _vcf_records(str(out))
    assert [(r.pos, r.ref, r.alts[0]) for r in recs][0] == (1, REF[0], "<NON_SEQUENCE>")
    assert len(recs) == 2


@pytest.mark.parametrize("is_maf", [True, False], ids=["maf", "vcf-style"])
def test_prepare_fails_a_row_before_the_first_base_that_is_not_an_insertion(tmp_path, is_maf):
    """Whatever the input format, a variant at position -1 that is not a MAF dash insertion
    is FAIL START_BEFORE_CONTIG (it must never reach a reference fetch at -1). A dash
    allele exists only in a MAF."""
    fa = _plain_fasta(tmp_path, REF[:60])
    variants = [_rs.Variant("1", -1, "A", "C", "SNP"), _rs.Variant("1", -1, "AC", "G", "COMPLEX")]
    if is_maf:
        variants.append(_rs.Variant("1", -1, "A", "-", "DELETION"))
    got = _rs.prepare_variants(variants, fa, 5, is_maf, 1, True)
    assert [(p.gbcms_status, p.gbcms_status_reason) for p in got] == [
        ("FAIL", "START_BEFORE_CONTIG")
    ] * len(variants)
    assert [p.original_pos for p in got] == [-1] * len(variants)


def test_a_row_with_a_malformed_allele_keeps_its_own_reason_at_start_0(tmp_path):
    """Guard: the allele checks run first, so a Start-0 row with an IUPAC ALT is still
    NON_SEQUENCE_ALLELE, and ALT == REF is still ALT_EQUALS_REF."""
    fa = _plain_fasta(tmp_path, REF[:60])
    variants = [_rs.Variant("1", -1, "A", "R", "SNP"), _rs.Variant("1", -1, "-", "-", "INSERTION")]
    got = _rs.prepare_variants(variants, fa, 5, True, 1, True)
    assert [p.gbcms_status_reason for p in got] == ["NON_SEQUENCE_ALLELE", "ALT_EQUALS_REF"]


# ── I2: End_Position optional ────────────────────────────────────────────────


@pytest.mark.parametrize("header", ["no column", "empty values"])
def test_end_position_is_optional(tmp_path, header):
    alt = "A" if REF[P] != "A" else "C"
    if header == "no column":
        head = [
            "Hugo_Symbol",
            "Chromosome",
            "Start_Position",
            "Reference_Allele",
            "Tumor_Seq_Allele2",
            "Tumor_Sample_Barcode",
        ]
        maf = _maf(tmp_path, [_row(P + 1, REF[P], alt)], header=head)
    else:
        maf = _maf(tmp_path, [_row(P + 1, REF[P], alt, end="")])
    (row,) = _maf_rows(_run(tmp_path, maf, _ref_reads()))
    assert (row["gbcms_status"], row["ref_count"]) == ("PASS", "6")


# ── I3: VCF input's MAF output fills Tumor_Seq_Allele1 ───────────────────────


@pytest.mark.parametrize(
    "ref,alt,allele1",
    [("C", "T", "C"), ("C", "CTT", "-"), ("CTT", "C", "TT"), ("CAG", "TTA", "CAG")],
    ids=["SNP", "insertion", "deletion", "MNP"],
)
def test_vcf_input_maf_output_fills_allele1_with_the_reference_allele(ref, alt, allele1):
    """Tumor_Seq_Allele1 is the MAF reference allele (after vcf2maf's trim)."""
    v = Variant(chrom="1", pos=99, ref=ref, alt=alt, variant_type=VariantType.SNP)
    fields = MafWriter.vcf_input_fields(v)
    assert fields["Tumor_Seq_Allele1"] == fields["Reference_Allele"] == allele1


# ── I4: one allele per row ───────────────────────────────────────────────────


def test_a_differing_allele1_is_not_a_second_allele(tmp_path):
    """Guard: REF, Allele1 and Allele2 all differ: one row, counted for Allele2."""
    alt2 = "A" if REF[P] != "A" else "C"
    alt1 = next(b for b in "ACGT" if b not in (REF[P], alt2))
    hap = REF[:P] + alt2 + REF[P + 1 :]
    reads = _ref_reads() + [
        make_read(f"a{i}", hap[P - 40 + i : P - 40 + i + READ], P - 40 + i, ((0, READ),))
        for i in range(4)
    ]
    head = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "End_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele1",
        "Tumor_Seq_Allele2",
        "Tumor_Sample_Barcode",
    ]
    maf = _maf(tmp_path, [_row(P + 1, REF[P], alt2, Tumor_Seq_Allele1=alt1)], header=head)
    out = _maf_rows(_run(tmp_path, maf, reads))
    assert len(out) == 1
    assert (out[0]["Tumor_Seq_Allele2"], out[0]["alt_count"]) == (alt2, "4")


# ── #147: decomposed and sibling lists ───────────────────────────────────────

_ARGS = {
    "min_mapq": 20,
    "min_baseq": 20,
    "filter_duplicates": True,
    "filter_secondary": True,
    "filter_supplementary": True,
    "filter_qc_failed": False,
    "filter_improper_pair": False,
    "filter_indel": False,
    "threads": 1,
}


def _two_snvs(tmp_path):
    fa, bam = _files(tmp_path, _ref_reads())
    vs = [
        _rs.Variant("1", P + k, REF[P + k], "A" if REF[P + k] != "A" else "C", "SNP")
        for k in (0, 5)
    ]
    return str(bam), vs


def test_a_short_decomposed_list_is_padded(tmp_path):
    bam, vs = _two_snvs(tmp_path)
    full = _rs.count_bam_binned(bam, vs, [None, None], **_ARGS)
    short = _rs.count_bam_binned(bam, vs, [None], **_ARGS)
    assert [c.rd for c in short] == [c.rd for c in full] == [6, 6]


@pytest.mark.parametrize("which", ["decomposed", "sibling_variants"])
def test_a_list_longer_than_the_variants_is_an_error(tmp_path, which):
    bam, vs = _two_snvs(tmp_path)
    kwargs = dict(_ARGS)
    decomposed = [None, None]
    if which == "decomposed":
        decomposed = [None, None, None]
    else:
        kwargs["sibling_variants"] = [[], [], []]
    with pytest.raises(ValueError, match="3"):
        _rs.count_bam_binned(bam, vs, decomposed, **kwargs)


def test_merge_joins_outputs_of_mafs_without_end_position(tmp_path):
    """I2's reach: outputs of a MAF without End_Position merge, joined on the
    rest of the variant key (End_Position follows from Start and REF)."""
    from gbcms.merge import merge_mafs
    from gbcms.models.core import MergeConfig

    head = [
        "Hugo_Symbol",
        "Chromosome",
        "Start_Position",
        "Reference_Allele",
        "Tumor_Seq_Allele2",
        "Tumor_Sample_Barcode",
    ]
    rows = [
        _row(P + 1, REF[P], "A" if REF[P] != "A" else "C"),
        _row(P + 11, REF[P + 10 : P + 12], "-"),
    ]
    maf = _maf(tmp_path, rows, header=head)
    d = _run(tmp_path, maf, _ref_reads(), name="d")
    s = _run(tmp_path, maf, _ref_reads(), name="s")
    out = tmp_path / "m.maf"
    merge_mafs(MergeConfig(inputs={"duplex": d, "simplex": s}, output=out))
    merged = list(read_maf_output(out))
    assert [(r["Start_Position"], r["Reference_Allele"]) for r in merged] == [
        (str(r["Start_Position"]), r["Reference_Allele"]) for r in rows
    ]
    assert merged[0]["simplex_duplex_ref_count"] == "12"


_MERGE_COLS = [
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "ref_count",
    "alt_count",
]


def _merge(tmp_path, inputs, caplog=None):
    """Merge hand-written gbcms MAFs ({type: (columns, rows)}); the merged rows."""
    import logging

    from gbcms.merge import merge_mafs
    from gbcms.models.core import MergeConfig

    paths = {}
    for name, (cols, rows) in inputs.items():
        paths[name] = tmp_path / f"{name}.maf"
        paths[name].write_text("\n".join(["\t".join(cols), *("\t".join(r) for r in rows)]) + "\n")
    out = tmp_path / "merged.maf"
    config = MergeConfig(inputs=paths, output=out, add_combined=False)
    if caplog is None:
        merge_mafs(config)
    else:
        with caplog.at_level(logging.INFO, logger="gbcms.merge"):
            merge_mafs(config)
    return list(read_maf_output(out))


def test_merge_joins_inputs_that_disagree_on_end_position(tmp_path, caplog):
    """End_Position follows from Start and REF and is not part of a variant's
    identity: inputs that write it differently for one variant join into one
    row, which keeps the first input's End_Position; the difference is logged."""
    merged = _merge(
        tmp_path,
        {
            "duplex": (_MERGE_COLS, [["1", "100", "100", "A", "T", "20", "10"]]),
            "simplex": (_MERGE_COLS, [["1", "100", "101", "A", "T", "5", "2"]]),
        },
        caplog,
    )
    assert len(merged) == 1, merged
    row = merged[0]
    assert row["End_Position"] == "100"
    assert (row["duplex_ref_count"], row["simplex_ref_count"]) == ("20", "5")
    assert not any(c.startswith("_") for c in row), list(row)
    logged = [r.message for r in caplog.records if "End_Position" in r.message]
    assert any("simplex" in m for m in logged), [r.message for r in caplog.records]


@pytest.mark.parametrize("second_has_end", [True, False])
def test_merge_keeps_end_position_of_a_row_only_a_later_input_has(tmp_path, second_has_end):
    """Guard (green now, green after): a row only a later input has keeps the
    End_Position that input gives; an input without the column gives none."""
    cols = _MERGE_COLS if second_has_end else [c for c in _MERGE_COLS if c != "End_Position"]
    later = [["1", "100", "100", "A", "T", "5", "2"], ["1", "200", "202", "GCA", "-", "7", "1"]]
    if not second_has_end:
        later = [[v for c, v in zip(_MERGE_COLS, r, strict=True) if c in cols] for r in later]
    merged = _merge(
        tmp_path,
        {
            "duplex": (_MERGE_COLS, [["1", "100", "100", "A", "T", "20", "10"]]),
            "simplex": (cols, later),
        },
    )
    assert [(r["Start_Position"], r["End_Position"]) for r in merged] == [
        ("100", "100"),
        ("200", "202" if second_has_end else ""),
    ]
    assert [r["simplex_ref_count"] for r in merged] == ["5", "7"]


def test_convert_maf_to_vcf_carries_the_maf_row(tmp_path):
    """gbcms convert writes VCF output of MAF input too: the MAF row in INFO, and
    a non-sequence allele as the symbolic record."""
    fa, _ = _files(tmp_path, [])
    rows = [_row(P + 1, REF[P : P + 2], "-"), _row(P + 21, REF[P + 20], "LU")]
    out = tmp_path / "o.vcf"
    res = runner.invoke(
        app, ["convert", "-v", str(_maf(tmp_path, rows)), "-f", str(fa), "-o", str(out)]
    )
    assert res.exit_code == 0, res.output
    recs = _vcf_records(str(out))
    assert [(r.pos, r.ref, r.alts[0]) for r in recs] == [
        (P, REF[P - 1 : P + 2], REF[P - 1]),
        (P + 21, REF[P + 20], "<NON_SEQUENCE>"),
    ]
    assert [(r.info["MAF_START"], r.info["MAF_REF"], r.info["MAF_ALT"]) for r in recs] == [
        (P + 1, REF[P : P + 2], "-"),
        (P + 21, REF[P + 20], "LU"),
    ]


def test_maf_origin_values_are_percent_encoded(tmp_path):
    """A hand-edited allele with characters an INFO value cannot hold (';', '=',
    ',', a space) is percent-encoded, so the VCF stays parseable and the value
    decodes to the input."""
    fa, _ = _files(tmp_path, [])
    rows = [_row(P + 1, REF[P], "A;B=C, D")]
    out = tmp_path / "o.vcf"
    res = runner.invoke(
        app, ["convert", "-v", str(_maf(tmp_path, rows)), "-f", str(fa), "-o", str(out)]
    )
    assert res.exit_code == 0, res.output
    line = [x for x in open(out) if not x.startswith("#")][0]
    assert "MAF_ALT=A%3BB%3DC%2C%20D" in line
    (rec,) = _vcf_records(str(out))
    assert rec.alts == ("<NON_SEQUENCE>",)


def _convert_line(tmp_path, rows, header=None):
    fa, _ = _files(tmp_path, [])
    out = tmp_path / "o.vcf"
    res = runner.invoke(
        app, ["convert", "-v", str(_maf(tmp_path, rows, header)), "-f", str(fa), "-o", str(out)]
    )
    assert res.exit_code == 0, res.output
    return [x for x in open(out) if not x.startswith("#")], out


def test_maf_origin_writes_the_rows_own_placeholder_alleles(tmp_path):
    """MAF_REF / MAF_ALT are the row's values as written, so the record can be
    looked up by its input: a placeholder REF ('0', '--') is not rewritten as
    the '-' preparation reads it as."""
    lines, _ = _convert_line(tmp_path, [_row(P + 1, "0", "T"), _row(P + 11, "--", "GG")])
    assert "MAF_REF=0;" in lines[0], lines[0]
    assert "MAF_REF=--;" in lines[1], lines[1]


def test_an_empty_maf_allele_is_missing_in_the_maf_origin(tmp_path):
    """An empty allele (FAIL EMPTY_ALLELE) is written as the VCF missing value
    '.', and a literal '.' allele is encoded, so '.' only ever means missing.
    The symbolic record's header names no single reason (GSR has it)."""
    lines, out = _convert_line(tmp_path, [_row(P + 1, "", "A"), _row(P + 21, REF[P + 20], ".")])
    assert "MAF_REF=.;" in lines[0], lines[0]
    assert "MAF_ALT=%2E" in lines[1], lines[1]
    with pysam.VariantFile(str(out)) as fh:
        assert "GSR" in fh.header.alts["NON_SEQUENCE"].get("Description")


@pytest.mark.parametrize("ref_len", [1, 2])
def test_a_dash_allele_is_not_a_base_outside_maf_input(tmp_path, ref_len):
    """'-' is a MAF dash allele only: given as non-MAF input (the observations
    API) it is FAIL NON_SEQUENCE_ALLELE, not a PASS row that counts 0."""
    fa, _ = _files(tmp_path, [])
    v = _rs.Variant("1", P, REF[P : P + ref_len], "-", "SNP")
    (pv,) = _rs.prepare_variants([v], str(fa), 5, False)
    assert (pv.gbcms_status, pv.gbcms_status_reason) == ("FAIL", "NON_SEQUENCE_ALLELE")


def test_merge_pairs_rows_of_one_input_that_differ_only_in_end_position(tmp_path, caplog):
    """When an input has rows that share contig, Start and alleles but not
    End_Position, merge joins on End_Position too (every input has it), so the
    rows pair with their own counterparts instead of every combination."""
    merged = _merge(
        tmp_path,
        {
            "duplex": (
                _MERGE_COLS,
                [
                    ["1", "100", "100", "A", "T", "20", "10"],
                    ["1", "100", "101", "A", "T", "3", "1"],
                ],
            ),
            "simplex": (
                _MERGE_COLS,
                [["1", "100", "100", "A", "T", "5", "2"], ["1", "100", "101", "A", "T", "7", "1"]],
            ),
        },
        caplog,
    )
    assert [(r["End_Position"], r["duplex_ref_count"], r["simplex_ref_count"]) for r in merged] == [
        ("100", "20", "5"),
        ("101", "3", "7"),
    ]
    assert any("End_Position" in r.message for r in caplog.records)


# ── REF_MISMATCH: where the given REF does sit ───────────────────────────────


def _shifted(offset, n=6):
    """(1-based Start, REF): a REF of n bases that matches the reference
    exactly at Start+offset and nowhere else within 3 bases, and under 90% at
    Start (a REF_MISMATCH)."""
    for q in range(P + 60, len(REF) - 60):
        ref = REF[q : q + n]
        s = q - offset
        at = REF[s : s + n]
        exact = [o for o in range(-3, 4) if REF[s + o : s + o + n] == ref]
        if exact == [offset] and sum(a == b for a, b in zip(at, ref, strict=True)) / n < 0.9:
            return s + 1, ref
    raise AssertionError("no such locus in the test reference")


@pytest.mark.parametrize("offset", [-1, 2, -3])
def test_ref_mismatch_names_where_the_given_ref_sits(tmp_path, offset):
    """The row stays FAIL REF_MISMATCH with zero counts; gbcms_diagnostic says
    the given REF matches the reference exactly `offset` bases from Start."""
    start, ref = _shifted(offset)
    out = _maf_rows(_run(tmp_path, _maf(tmp_path, [_row(start, ref, ref[0])]), _ref_reads()))
    row = out[0]
    assert (row["gbcms_status"], row["gbcms_status_reason"]) == ("FAIL", "REF_MISMATCH")
    assert row["gbcms_diagnostic"] == f"REF_AT_OFFSET({offset:+d})"
    assert (row["ref_count"], row["alt_count"]) == ("0", "0")


def test_ref_mismatch_offset_of_a_maf_dash_deletion_is_the_given_bases(tmp_path):
    """A MAF '-' deletion: the deleted bases as given, from Start (the anchor is
    the reference's own base and is not part of what the row gives)."""
    start, ref = _shifted(1)
    (row,) = _maf_rows(_run(tmp_path, _maf(tmp_path, [_row(start, ref, "-")]), _ref_reads()))
    assert (row["gbcms_status_reason"], row["gbcms_diagnostic"]) == (
        "REF_MISMATCH",
        "REF_AT_OFFSET(+1)",
    )


def test_ref_mismatch_offset_in_vcf_output_and_normalize(tmp_path):
    """The diagnostic reaches VCF output (GD) and gbcms normalize's TSV."""
    start, ref = _shifted(-1)
    maf = _maf(tmp_path, [_row(start, ref, ref[0])])
    (rec,) = _vcf_records(_run(tmp_path, maf, _ref_reads(), fmt="vcf"))
    assert (rec.info["GSR"], rec.info["GD"]) == ("REF_MISMATCH", "REF_AT_OFFSET(-1)")
    fa = tmp_path / "ref.fa"
    tsv = tmp_path / "n.tsv"
    res = runner.invoke(app, ["normalize", "-v", str(maf), "-f", str(fa), "-o", str(tsv)])
    assert res.exit_code == 0, res.output
    (norm,) = list(csv.DictReader(open(tsv), delimiter="\t"))
    assert norm["gbcms_diagnostic"] == "REF_AT_OFFSET(-1)"


def test_ref_mismatch_lists_every_offset_in_a_repeat(tmp_path):
    """In a repeat the given REF can sit at several offsets: all within 3 bases
    are listed, nearest first (the left one first at equal distance)."""
    fa = tmp_path / "rep.fa"
    seq = REF[:300] + "AC" * 20 + REF[340:]
    fa.write_text(f">1\n{seq}\n")
    pysam.faidx(str(fa))
    # The repeat starts at 0-based 300 (an A); 311 is inside it, a C, so the REF
    # "ACACAC" mismatches there and sits at 310, 312, 308 and 314.
    v = _rs.Variant("1", 311, "ACACAC", "A", "COMPLEX")
    (pv,) = _rs.prepare_variants([v], str(fa), 5, False)
    assert pv.gbcms_status_reason == "REF_MISMATCH"
    assert pv.gbcms_diagnostic == "REF_AT_OFFSET(-1/+1/-3/+3)"


@pytest.mark.parametrize(
    "case",
    ["nowhere near", "two bases"],
)
def test_ref_mismatch_without_a_nearby_exact_ref_names_nothing(tmp_path, case):
    """Guard: no offset when the REF matches nowhere within 3 bases, nor for a
    REF under 3 bases (a nearby match is then often chance)."""
    fa, _ = _files(tmp_path, [])
    if case == "two bases":
        start, ref = next(
            (q + 2, REF[q : q + 2])
            for q in range(P + 60, len(REF) - 60)
            if REF[q + 1 : q + 3] != REF[q : q + 2] and REF[q + 1] != REF[q]
        )
    else:
        ref = "".join({"A": "C", "C": "G", "G": "T", "T": "A"}[b] for b in REF[P : P + 8])
        start = P + 1
        assert all(REF[P + o : P + o + 8] != ref for o in range(-3, 4))
    (pv,) = _rs.prepare_variants(
        [_rs.Variant("1", start - 1, ref, ref[0], "COMPLEX")], str(fa), 5, False
    )
    assert pv.gbcms_status_reason == "REF_MISMATCH"
    assert pv.gbcms_diagnostic == ""
