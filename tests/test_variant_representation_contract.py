"""Variant representation contract (#110).

VCF <-> MAF conversion follows vcf2maf / maf2vcf, internal type labels follow
the alleles, and the VCF reader skips ALT alleles it cannot count — loudly.

Oracles: ``vcf2maf.pl`` and ``maf2vcf.pl`` (github.com/mskcc/vcf2maf, commit
589406f), run on the synthetic 17-shape battery regenerated below (no patient
data). The expected tables are their literal outputs:

- VCF -> MAF: vcf2maf trims the leading bases REF and ALT share (never the
  trailing ones); an allele trimmed to nothing becomes ``-``. Equal trimmed
  lengths are SNP/DNP/TNP/ONP by length, otherwise INS or DEL; an insertion
  whose REF trimmed to ``-`` spans the two bases around it.
- MAF -> VCF: maf2vcf prepends the reference base before Start (the base AT
  Start for a ``-`` insertion) when an allele is ``-``, or when the lengths
  differ and the first bases differ; anything else is written as-is at Start.

Counting is not part of this change: the counts guard pins that every shape
still counts 10 ALT / 10 REF from both inputs.
"""

import glob
import logging
import random

import pysam
import pytest
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms import _rs as gbcms_rs
from gbcms.cli import app
from gbcms.core.kernel import CoordinateKernel
from gbcms.io.input import MafReader, VcfReader
from gbcms.io.output import MafWriter, VcfWriter
from gbcms.merge import merge_mafs
from gbcms.models.core import MergeConfig, Variant, VariantType
from gbcms.pipeline import _zero_counts

runner = CliRunner()

# The battery's VCF records: shape -> (POS, REF, ALT).
VCF_INPUT = {
    "S01_snv": (101, "T", "A"),
    "S02_dnp": (221, "AG", "CT"),
    "S03_onp4": (341, "TACG", "ATTA"),
    "S04_del_anchored": (461, "TAA", "T"),
    "S05_ins_anchored": (581, "T", "TGT"),
    "S06_delins_alt1_noshare": (701, "TTAC", "A"),
    "S07_delins_ref1_noshare": (821, "C", "TA"),
    "S08_delins_shared1_netdel": (941, "GAAA", "GT"),
    "S09_cplx_ins_shared1": (1061, "CG", "CCAG"),
    "S10_shared2_to_snv": (1181, "TCT", "TCG"),
    "S11_shared2_to_dnp": (1301, "TATG", "TAAC"),
    "S12_suffix_only_eq": (1421, "TTA", "GTA"),
    "S13_suffix_only_netdel": (1541, "ACG", "CG"),
    "S14_delins_netins_noshare": (1661, "AT", "CGTC"),
    "S15_homopolymer_delins": (1781, "CCCCCC", "T"),
    "S16_shared2_del": (1901, "AGAGC", "AG"),
    "S17_shared2_ins": (2021, "TC", "TCGG"),
}

# vcf2maf.pl output: shape -> (Start, End, Reference_Allele, Tumor_Seq_Allele2, Variant_Type).
VCF2MAF = {
    "S01_snv": (101, 101, "T", "A", "SNP"),
    "S02_dnp": (221, 222, "AG", "CT", "DNP"),
    "S03_onp4": (341, 344, "TACG", "ATTA", "ONP"),
    "S04_del_anchored": (462, 463, "AA", "-", "DEL"),
    "S05_ins_anchored": (581, 582, "-", "GT", "INS"),
    "S06_delins_alt1_noshare": (701, 704, "TTAC", "A", "DEL"),
    "S07_delins_ref1_noshare": (821, 821, "C", "TA", "INS"),
    "S08_delins_shared1_netdel": (942, 944, "AAA", "T", "DEL"),
    "S09_cplx_ins_shared1": (1062, 1062, "G", "CAG", "INS"),
    "S10_shared2_to_snv": (1183, 1183, "T", "G", "SNP"),
    "S11_shared2_to_dnp": (1303, 1304, "TG", "AC", "DNP"),
    "S12_suffix_only_eq": (1421, 1423, "TTA", "GTA", "TNP"),
    "S13_suffix_only_netdel": (1541, 1543, "ACG", "CG", "DEL"),
    "S14_delins_netins_noshare": (1661, 1662, "AT", "CGTC", "INS"),
    "S15_homopolymer_delins": (1781, 1786, "CCCCCC", "T", "DEL"),
    "S16_shared2_del": (1903, 1905, "AGC", "-", "DEL"),
    "S17_shared2_ins": (2022, 2023, "-", "GG", "INS"),
}

# maf2vcf.pl output for the vcf2maf MAF above: shape -> (POS, REF, ALT).
MAF2VCF = {
    "S01_snv": (101, "T", "A"),
    "S02_dnp": (221, "AG", "CT"),
    "S03_onp4": (341, "TACG", "ATTA"),
    "S04_del_anchored": (461, "TAA", "T"),
    "S05_ins_anchored": (581, "T", "TGT"),
    "S06_delins_alt1_noshare": (700, "TTTAC", "TA"),
    "S07_delins_ref1_noshare": (820, "AC", "ATA"),
    "S08_delins_shared1_netdel": (941, "GAAA", "GT"),
    "S09_cplx_ins_shared1": (1061, "CG", "CCAG"),
    "S10_shared2_to_snv": (1183, "T", "G"),
    "S11_shared2_to_dnp": (1303, "TG", "AC"),
    "S12_suffix_only_eq": (1421, "TTA", "GTA"),
    "S13_suffix_only_netdel": (1540, "CACG", "CCG"),
    "S14_delins_netins_noshare": (1660, "GAT", "GCGTC"),
    "S15_homopolymer_delins": (1780, "ACCCCCC", "AT"),
    "S16_shared2_del": (1902, "GAGC", "G"),
    "S17_shared2_ins": (2022, "C", "CGG"),
}


def _battery():
    """The synthetic reference and shape table the oracles ran on (seeded)."""
    rng = random.Random(110)
    ref = [rng.choice("ACGT") for _ in range(2400)]

    def other(b):
        return rng.choice([x for x in "ACGT" if x != b])

    def seq(p, n):
        return "".join(ref[p : p + n])

    rows = {}

    def add(name, p, r, a):
        rows[name] = (p + 1, r, a)

    P = iter(range(100, 2300, 120))
    p = next(P)
    add("S01_snv", p, seq(p, 1), other(ref[p]))
    p = next(P)
    add("S02_dnp", p, seq(p, 2), other(ref[p]) + other(ref[p + 1]))
    p = next(P)
    add("S03_onp4", p, seq(p, 4), "".join(other(ref[p + i]) for i in range(4)))
    p = next(P)
    add("S04_del_anchored", p, seq(p, 3), seq(p, 1))
    p = next(P)
    add("S05_ins_anchored", p, seq(p, 1), seq(p, 1) + "GT")
    p = next(P)
    add("S06_delins_alt1_noshare", p, seq(p, 4), other(ref[p]))
    p = next(P)
    add("S07_delins_ref1_noshare", p, seq(p, 1), other(ref[p]) + "A")
    p = next(P)
    add("S08_delins_shared1_netdel", p, seq(p, 4), seq(p, 1) + other(ref[p + 1]))
    p = next(P)
    add("S09_cplx_ins_shared1", p, seq(p, 2), seq(p, 1) + "CA" + ref[p + 1])
    p = next(P)
    add("S10_shared2_to_snv", p, seq(p, 3), seq(p, 2) + other(ref[p + 2]))
    p = next(P)
    add("S11_shared2_to_dnp", p, seq(p, 4), seq(p, 2) + other(ref[p + 2]) + other(ref[p + 3]))
    p = next(P)
    add("S12_suffix_only_eq", p, seq(p, 3), other(ref[p]) + seq(p + 1, 2))
    p = next(P)
    add("S13_suffix_only_netdel", p, seq(p, 3), other(ref[p]) + ref[p + 2])
    p = next(P)
    add("S14_delins_netins_noshare", p, seq(p, 2), other(ref[p]) + "GTC")
    p = next(P)
    ref[p : p + 7] = list("CCCCCCT")
    ref[p - 1] = "A"
    add("S15_homopolymer_delins", p, "CCCCCC", "T")
    p = next(P)
    add("S16_shared2_del", p, seq(p, 5), seq(p, 2))
    p = next(P)
    add("S17_shared2_ins", p, seq(p, 2), seq(p, 2) + "GG")
    return "".join(ref), rows


REF_SEQ, BATTERY = _battery()


def _fasta(path):
    fasta = path / "ref.fa"
    fasta.write_text(">1\n" + REF_SEQ + "\n")
    pysam.faidx(str(fasta))
    return fasta


def _maf_text(rows=VCF2MAF):
    head = "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode\n"
    body = "".join(
        f"G\t1\t{s}\t{e}\t{r}\t{a}\tT1\n" for s, e, r, a, _ in (rows[k] for k in VCF_INPUT)
    )
    return head + body


def _vcf_text():
    head = (
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=2400>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    return head + "".join(
        f"1\t{p}\t{k}\t{r}\t{a}\t.\tPASS\t.\n" for k, (p, r, a) in VCF_INPUT.items()
    )


def _maf_variant(shape):
    s, e, r, a, _ = VCF2MAF[shape]
    row = {
        "Chromosome": "1",
        "Start_Position": str(s),
        "End_Position": str(e),
        "Reference_Allele": r,
        "Tumor_Seq_Allele2": a,
    }
    return CoordinateKernel.maf_to_internal("1", s, e, r, a).model_copy(update={"metadata": row})


def _vcf_records(path):
    return [line.rstrip("\n").split("\t") for line in open(path) if not line.startswith("#")]


# ── The battery is the one the oracles ran on ────────────────────────────


def test_battery_matches_the_oracle_inputs():
    """The regenerated shapes are the literal VCF the oracles converted."""
    assert BATTERY == VCF_INPUT


# ── Kernel: one representation module ────────────────────────────────────


def test_vcf_to_maf_is_vcf2maf():
    got = {}
    for k, (p, r, a) in VCF_INPUT.items():
        m = CoordinateKernel.vcf_to_maf(p, r, a)
        got[k] = (
            int(m["Start_Position"]),
            int(m["End_Position"]),
            m["Reference_Allele"],
            m["Tumor_Seq_Allele2"],
            m["Variant_Type"],
        )
    assert got == VCF2MAF


def test_maf_to_vcf_is_maf2vcf():
    """The anchor base is fetched only when maf2vcf prepends one."""
    got, fetched = {}, {}
    for k, (s, _e, r, a, _t) in VCF2MAF.items():
        calls = []

        def base_at(pos1, calls=calls):
            calls.append(pos1)
            return REF_SEQ[pos1 - 1]

        got[k] = CoordinateKernel.maf_to_vcf(s, r, a, base_at)
        fetched[k] = bool(calls)
    assert got == MAF2VCF
    padded = {k: len(MAF2VCF[k][1]) != len(VCF2MAF[k][2].replace("-", "")) for k in VCF2MAF}
    assert fetched == padded


def test_allele_type_rule():
    """INSERTION/DELETION only when the single-base allele is the other's
    first base (the shared anchor, compared case-insensitively); a multi-base
    substitution or a delins without a shared anchor is COMPLEX."""
    table = {
        ("T", "A"): VariantType.SNP,
        ("TAA", "T"): VariantType.DELETION,
        ("T", "TGT"): VariantType.INSERTION,
        ("tAA", "T"): VariantType.DELETION,
        ("t", "TG"): VariantType.INSERTION,
        ("TTAC", "A"): VariantType.COMPLEX,
        ("C", "TA"): VariantType.COMPLEX,
        ("GAAA", "GT"): VariantType.COMPLEX,
        ("AG", "CT"): VariantType.COMPLEX,
    }
    assert {k: CoordinateKernel.allele_type(*k) for k in table} == table


def test_readers_label_by_allele_type():
    """Both readers use the one rule; MAF '-' alleles stay INSERTION/DELETION."""
    for p, r, a in VCF_INPUT.values():
        v = CoordinateKernel.vcf_to_internal("1", p, r, a)
        assert v.variant_type == CoordinateKernel.allele_type(r, a), (p, r, a)
    assert CoordinateKernel.vcf_to_internal("1", 701, "TTAC", "A").variant_type == "COMPLEX"
    assert CoordinateKernel.vcf_to_internal("1", 821, "C", "TA").variant_type == "COMPLEX"
    maf = CoordinateKernel.maf_to_internal
    assert maf("1", 462, 463, "AA", "-").variant_type == VariantType.DELETION
    assert maf("1", 581, 582, "-", "GT").variant_type == VariantType.INSERTION
    assert maf("1", 461, 462, "TA", "T").variant_type == VariantType.DELETION
    assert maf("1", 701, 704, "TTAC", "A").variant_type == VariantType.COMPLEX


def test_engine_labels_follow_the_alleles(tmp_path):
    """prepare_variants derives every PASS label from its final alleles with the
    same rule as the kernel, whatever label it was given."""
    fasta = _fasta(tmp_path)
    rs = [gbcms_rs.Variant("1", p - 1, r, a, "SNP") for p, r, a in VCF_INPUT.values()]
    prepared = gbcms_rs.prepare_variants(rs, str(fasta), 5, False, 1, True)
    maf_rs = [
        gbcms_rs.Variant("1", 461, "AA", "-", "SNP"),  # MAF deletion (Start 462)
        gbcms_rs.Variant("1", 580, "-", "GT", "SNP"),  # MAF insertion (Start 581)
    ]
    prepared += gbcms_rs.prepare_variants(maf_rs, str(fasta), 5, True, 1, True)
    assert [pv.gbcms_status for pv in prepared] == ["PASS"] * len(prepared)
    got = [pv.variant.variant_type for pv in prepared]
    want = [
        CoordinateKernel.allele_type(pv.variant.ref_allele, pv.variant.alt_allele).value
        for pv in prepared
    ]
    assert got == want
    assert got[-2:] == ["DELETION", "INSERTION"]


# ── Writers ──────────────────────────────────────────────────────────────


def test_maf_writer_vcf_rows_are_vcf2maf(tmp_path):
    out = tmp_path / "o.maf"
    w = MafWriter(out)
    for k, (p, r, a) in VCF_INPUT.items():
        w.write(CoordinateKernel.vcf_to_internal("1", p, r, a, original_id=k), _zero_counts())
    w.close()
    rows = {row["vcf_id"]: row for row in read_maf_output(out)}
    got = {
        k: (
            int(x["Start_Position"]),
            int(x["End_Position"]),
            x["Reference_Allele"],
            x["Tumor_Seq_Allele2"],
            x["Variant_Type"],
        )
        for k, x in rows.items()
    }
    assert got == VCF2MAF
    # The VCF record itself is kept (vcf2maf's vcf_pos / vcf_ref / vcf_alt).
    assert {k: (int(x["vcf_pos"]), x["vcf_ref"], x["vcf_alt"]) for k, x in rows.items()} == (
        VCF_INPUT
    )


def test_maf_writer_norm_columns_follow_the_alleles(tmp_path):
    """A left-aligned delins handed over with a DELETION label is still written
    by its alleles (TTAC>A keeps its first base: vcf2maf's DEL 701-704)."""
    out = tmp_path / "o.maf"
    w = MafWriter(out, show_normalization=True)
    v = CoordinateKernel.vcf_to_internal("1", 701, "TTAC", "A")
    norm = Variant(chrom="1", pos=700, ref="TTAC", alt="A", variant_type=VariantType.DELETION)
    w.write(v, _zero_counts(), norm_variant=norm)
    w.close()
    row = next(iter(read_maf_output(out)))
    got = tuple(
        row[f"norm_{c}"]
        for c in ("Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2")
    )
    assert got == ("701", "704", "TTAC", "A")


def test_vcf_writer_maf_rows_are_maf2vcf(tmp_path):
    fasta = _fasta(tmp_path)
    out = tmp_path / "o.vcf"
    w = VcfWriter(out, reference_fasta=str(fasta))
    for k in VCF_INPUT:
        w.write(_maf_variant(k), _zero_counts())
    w.close()
    got = {k: (int(r[1]), r[3], r[4]) for k, r in zip(VCF_INPUT, _vcf_records(out), strict=True)}
    assert got == MAF2VCF


def test_vcf_writer_vcf_rows_echo_the_input(tmp_path):
    out = tmp_path / "o.vcf"
    w = VcfWriter(out)
    for k, (p, r, a) in VCF_INPUT.items():
        w.write(CoordinateKernel.vcf_to_internal("1", p, r, a, original_id=k), _zero_counts())
    w.close()
    got = {r[2]: (int(r[1]), r[3], r[4]) for r in _vcf_records(out)}
    assert got == VCF_INPUT


def test_vcf_writer_maf_row_needs_the_reference(tmp_path):
    w = VcfWriter(tmp_path / "o.vcf")
    with pytest.raises(ValueError, match="reference"):
        w.write(_maf_variant("S04_del_anchored"), _zero_counts())
    w.close()


def test_vcf_writer_unfetchable_anchor_is_n_and_warned(tmp_path, caplog):
    """A contig the FASTA lacks: the record stays valid VCF (anchor N) and the
    run says how many anchors it could not fetch."""
    fasta = _fasta(tmp_path)
    out = tmp_path / "o.vcf"
    w = VcfWriter(out, reference_fasta=str(fasta))
    v = CoordinateKernel.maf_to_internal("7", 462, 463, "AA", "-")
    v = v.model_copy(update={"metadata": {"Chromosome": "7"}})
    with caplog.at_level(logging.WARNING):
        w.write(v, _zero_counts())
        w.close()
    assert [(r[0], int(r[1]), r[3], r[4]) for r in _vcf_records(out)] == [("7", 461, "NAA", "N")]
    assert any("anchor" in m and "1" in m for m in caplog.messages)


# ── VCF reader: ALT alleles that cannot be counted ───────────────────────


def test_vcf_reader_skips_uncountable_alts(tmp_path, caplog):
    vcf = tmp_path / "s.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=2400>\n"
        '##ALT=<ID=DEL,Description="Deletion">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t101\t.\tT\t*\t.\t.\t.\n"
        "1\t221\t.\tA\t<DEL>\t.\t.\t.\n"
        "1\t341\t.\tT\tT]1:500]\t.\t.\t.\n"
        "1\t461\t.\tT\t.\t.\t.\t.\n"
        "1\t581\t.\tT\tA,*\t.\t.\t.\n"
        "1\t701\t.\tT\tR\t.\t.\t.\n"
        "1\t821\t.\tC\tN\t.\t.\t.\n"
    )
    with caplog.at_level(logging.WARNING):
        got = [(v.pos + 1, v.ref, v.alt) for v in VcfReader(vcf)]
    # ALT N is kept: preparation rejects it as a visible FAIL row (ALT_CONTAINS_N).
    assert got == [(581, "T", "A"), (821, "C", "N")]
    summary = " ".join(caplog.messages)
    for reason in ("'*'", "symbolic", "breakend", "missing", "non-sequence"):
        assert reason in summary, reason
    assert "6 ALT" in summary


# ── End to end: coordinates change, counts do not ────────────────────────


def _minimal(p0, r, a):
    while r and a and r[0] == a[0]:
        r, a, p0 = r[1:], a[1:], p0 + 1
    while r and a and r[-1] == a[-1]:
        r, a = r[:-1], a[:-1]
    return p0, r, a


def _bam(path):
    """Per shape: 10 reads carrying the ALT haplotype (the trimmed event as
    substitutions plus one D or I) and 10 reference reads, all forward."""
    reads = []
    for k, (pos, r, a) in VCF_INPUT.items():
        q, R, A = _minimal(pos - 1, r, a)
        hap = REF_SEQ[:q] + A + REF_SEQ[q + len(R) :]
        for i in range(10):
            s = q - 50 + i
            common = min(len(R), len(A))
            cig = [(0, q - s + common)]
            if len(R) > len(A):
                cig.append((2, len(R) - len(A)))
            elif len(A) > len(R):
                cig.append((1, len(A) - len(R)))
            cig.append((0, 100 - sum(n for op, n in cig if op in (0, 1))))
            reads.append(make_read(f"{k}a{i}", hap[s : s + 100], s, tuple(cig)))
            reads.append(make_read(f"{k}r{i}", REF_SEQ[s : s + 100], s, ((0, 100),)))
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": len(REF_SEQ)}]}
    raw = path / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=header) as fh:
        for rd in reads:
            fh.write(rd)
    bam = path / "s.bam"
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    return bam


@pytest.fixture(scope="module")
def e2e(tmp_path_factory):
    d = tmp_path_factory.mktemp("repr")
    fasta, bam = _fasta(d), _bam(d)
    (d / "shapes.vcf").write_text(_vcf_text())
    (d / "shapes.maf").write_text(_maf_text())
    out = {}
    for name, inp, fmt in (("vcf_in", "shapes.vcf", "maf"), ("maf_in", "shapes.maf", "vcf")):
        o = d / name
        res = runner.invoke(
            app,
            ["dna", "-v", str(d / inp), "-b", f"S:{bam}", "-f", str(fasta), "-o", str(o)]
            + ["--format", fmt],
        )
        assert res.exit_code == 0, res.output
        out[name] = glob.glob(str(o / f"*.{fmt}"))[0]
    out["dir"] = d
    return out


def _maf_rows(path):
    return {r["vcf_id"]: r for r in read_maf_output(path)}


def _vcf_samples(path):
    out = []
    for rec in _vcf_records(path):
        fmt = dict(zip(rec[8].split(":"), rec[9].split(":"), strict=True))
        out.append((rec, fmt))
    return out


def test_e2e_output_coordinates(e2e):
    rows = _maf_rows(e2e["vcf_in"])
    got = {
        k: (
            int(x["Start_Position"]),
            int(x["End_Position"]),
            x["Reference_Allele"],
            x["Tumor_Seq_Allele2"],
            x["Variant_Type"],
        )
        for k, x in rows.items()
    }
    assert got == VCF2MAF
    recs = _vcf_samples(e2e["maf_in"])
    got_vcf = {k: (int(r[1]), r[3], r[4]) for k, (r, _f) in zip(VCF_INPUT, recs, strict=True)}
    assert got_vcf == MAF2VCF
    # Every record is readable VCF.
    with pysam.VariantFile(e2e["maf_in"]) as fh:
        assert len(list(fh)) == len(VCF_INPUT)


def test_e2e_counts_unchanged(e2e):
    """Every shape counts 10 ALT / 10 REF from VCF and from MAF input."""
    for k, r in _maf_rows(e2e["vcf_in"]).items():
        rd, ad = int(r["ref_count"]), int(r["alt_count"])
        assert (rd, ad) == (10, 10), k
        assert int(r["total_count"]) >= rd + ad
        assert int(r["total_count_fragment"]) >= int(r["ref_count_fragment"]) + int(
            r["alt_count_fragment"]
        )
        assert rd == int(r["ref_count_forward"]) + int(r["ref_count_reverse"])
        assert ad == int(r["alt_count_forward"]) + int(r["alt_count_reverse"])
    for k, (_rec, f) in zip(VCF_INPUT, _vcf_samples(e2e["maf_in"]), strict=True):
        rd, ad = map(int, f["AD"].split(","))
        (rdf, adf), (rdr, adr) = map(int, f["ADF"].split(",")), map(int, f["ADR"].split(","))
        assert (rd, ad) == (10, 10), k
        assert int(f["DP"]) >= rd + ad
        assert rd == rdf + rdr and ad == adf + adr


# ── gbcms convert: the representation on its own ────────────────────────


def test_convert_vcf_to_maf(tmp_path):
    (tmp_path / "s.vcf").write_text(_vcf_text())
    res = runner.invoke(
        app, ["convert", "-v", str(tmp_path / "s.vcf"), "-o", str(tmp_path / "o.maf")]
    )
    assert res.exit_code == 0, res.output
    rows = _maf_rows(tmp_path / "o.maf")
    got = {
        k: (
            int(x["Start_Position"]),
            int(x["End_Position"]),
            x["Reference_Allele"],
            x["Tumor_Seq_Allele2"],
            x["Variant_Type"],
        )
        for k, x in rows.items()
    }
    assert got == VCF2MAF
    assert {k: (int(x["vcf_pos"]), x["vcf_ref"], x["vcf_alt"]) for k, x in rows.items()} == (
        VCF_INPUT
    )


def test_convert_maf_to_vcf(tmp_path):
    fasta = _fasta(tmp_path)
    (tmp_path / "s.maf").write_text(_maf_text())
    out = tmp_path / "o.vcf"
    res = runner.invoke(
        app, ["convert", "-v", str(tmp_path / "s.maf"), "-f", str(fasta), "-o", str(out)]
    )
    assert res.exit_code == 0, res.output
    with pysam.VariantFile(str(out)) as fh:
        got = [(rec.pos, rec.ref, rec.alts[0]) for rec in fh]
    assert got == [MAF2VCF[k] for k in VCF_INPUT]


def test_convert_maf_needs_the_reference(tmp_path):
    (tmp_path / "s.maf").write_text(_maf_text())
    res = runner.invoke(
        app, ["convert", "-v", str(tmp_path / "s.maf"), "-o", str(tmp_path / "o.vcf")]
    )
    assert res.exit_code == 1
    assert not (tmp_path / "o.vcf").exists()


def test_convert_output_must_be_the_other_format(tmp_path):
    (tmp_path / "s.vcf").write_text(_vcf_text())
    res = runner.invoke(
        app, ["convert", "-v", str(tmp_path / "s.vcf"), "-o", str(tmp_path / "o.vcf")]
    )
    assert res.exit_code == 1
    assert not (tmp_path / "o.vcf").exists()


# ── MAF alleles as maf2vcf reads them ────────────────────────────────────


def test_maf_alleles_follow_maf2vcf():
    """The variant allele is Tumor_Seq_Allele2, or Tumor_Seq_Allele1 when
    Allele2 is empty or the reference (older MAFs); alleles made only of
    '-', '?' or '0' are placeholders for an empty allele ('-')."""
    maf = CoordinateKernel.maf_alleles
    assert maf("A", "A", "G") == ("A", "G")
    assert maf("A", "G", "A") == ("A", "G")
    assert maf("A", "G", "") == ("A", "G")
    assert maf("C", "C", "--") == ("C", "-")
    assert maf("C", "C", "0") == ("C", "-")
    assert maf("--", "--", "TT") == ("-", "TT")
    assert maf("-", "GT", "-") == ("-", "GT")
    # Nothing else to use: REF == ALT stays, and preparation rejects it.
    assert maf("A", "A", "A") == ("A", "A")


def test_maf_reader_uses_allele1_when_allele2_is_the_reference(tmp_path, caplog):
    maf = tmp_path / "s.maf"
    maf.write_text(
        "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele1\t"
        "Tumor_Seq_Allele2\n"
        "1\t101\t101\tT\tT\tA\n"
        "1\t221\t222\tAG\tCT\tAG\n"
    )
    with caplog.at_level(logging.WARNING):
        got = [(v.pos + 1, v.ref, v.alt) for v in MafReader(maf)]
    assert got == [(101, "T", "A"), (221, "AG", "CT")]
    assert any("Tumor_Seq_Allele1" in m and "1 row" in m for m in caplog.messages)


def test_maf_to_vcf_at_position_1_uses_the_base_after():
    """VCF spec: an event at position 1 carries the base after it."""
    bases = {1: "C", 2: "G", 3: "T"}
    assert CoordinateKernel.maf_to_vcf(1, "C", "-", bases.__getitem__) == (1, "CG", "G")
    assert CoordinateKernel.maf_to_vcf(1, "CG", "A", bases.__getitem__) == (1, "CGT", "AT")


def test_vcf_reader_skips_non_sequence_ref_and_names_missing_alts(tmp_path, caplog):
    vcf = tmp_path / "s.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=2400>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t101\t.\t.\tA\t.\t.\t.\n"
        "1\t221\t.\tA\tG,.\t.\t.\t.\n"
    )
    with caplog.at_level(logging.WARNING):
        got = [(v.pos + 1, v.ref, v.alt) for v in VcfReader(vcf)]
    assert got == [(221, "A", "G")]
    summary = caplog.messages[-1]
    assert "REF" in summary and "missing ALT" in summary and "breakend" not in summary


def test_prepare_rejects_alt_equal_to_ref(tmp_path):
    """A 'variant' whose ALT is its REF (any case; or '-' for both in a MAF)
    is a visible FAIL row, not a silent PASS with meaningless counts."""
    fasta = _fasta(tmp_path)
    vcf_style = [
        gbcms_rs.Variant("1", 220, "AG", "AG", "COMPLEX"),
        gbcms_rs.Variant("1", 100, "t", "T", "SNP"),
    ]
    maf_style = [gbcms_rs.Variant("1", 580, "-", "-", "INSERTION")]
    prepared = gbcms_rs.prepare_variants(vcf_style, str(fasta), 5, False, 1, True)
    prepared += gbcms_rs.prepare_variants(maf_style, str(fasta), 5, True, 1, True)
    assert [(pv.gbcms_status, pv.gbcms_status_reason) for pv in prepared] == [
        ("FAIL", "ALT_EQUALS_REF")
    ] * 3


# ── Consumers keyed on the VCF record ────────────────────────────────────


def _vcf_input_maf(path, records):
    w = MafWriter(path)
    for p, r, a in records:
        w.write(CoordinateKernel.vcf_to_internal("1", p, r, a), _zero_counts())
    w.close()


def test_merge_joins_vcf_input_rows_by_their_record(tmp_path):
    """TCT>TCG and T>G trim to the same MAF record; each VCF record stays one
    merged row, joined to its own counterpart."""
    records = [(101, "T", "A"), (1181, "TCT", "TCG"), (1183, "T", "G")]
    _vcf_input_maf(tmp_path / "d.maf", records)
    _vcf_input_maf(tmp_path / "s.maf", records)
    out = tmp_path / "m.maf"
    merge_mafs(
        MergeConfig(
            inputs={"duplex": tmp_path / "d.maf", "simplex": tmp_path / "s.maf"}, output=out
        )
    )
    rows = list(read_maf_output(out))
    assert [(r["vcf_pos"], r["vcf_ref"], r["vcf_alt"]) for r in rows] == [
        (str(p), r, a) for p, r, a in records
    ]


@pytest.mark.xfail(strict=True, reason="a later-only row loses its MAF key columns")
def test_merge_by_record_keeps_later_only_rows_whole(tmp_path):
    """A variant only a later input has keeps its MAF coordinates and alleles."""
    _vcf_input_maf(tmp_path / "d.maf", [(101, "T", "A")])
    _vcf_input_maf(tmp_path / "s.maf", [(101, "T", "A"), (1183, "T", "G")])
    out = tmp_path / "m.maf"
    merge_mafs(
        MergeConfig(
            inputs={"duplex": tmp_path / "d.maf", "simplex": tmp_path / "s.maf"}, output=out
        )
    )
    got = [
        tuple(
            r[c]
            for c in ("Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2")
        )
        for r in read_maf_output(out)
    ]
    assert got == [("101", "101", "T", "A"), ("1183", "1183", "T", "G")]


def test_merge_warns_on_duplicate_join_keys(tmp_path, caplog):
    maf = (
        "Chromosome\tStart_Position\tEnd_Position\tReference_Allele\tTumor_Seq_Allele2\t"
        "ref_count\talt_count\n"
        "1\t101\t101\tT\tA\t5\t1\n"
        "1\t101\t101\tT\tA\t5\t1\n"
    )
    (tmp_path / "d.maf").write_text(maf)
    (tmp_path / "s.maf").write_text(maf)
    with caplog.at_level(logging.WARNING):
        merge_mafs(
            MergeConfig(
                inputs={"duplex": tmp_path / "d.maf", "simplex": tmp_path / "s.maf"},
                output=tmp_path / "m.maf",
            )
        )
    assert any("duplicate" in m and "duplex" in m for m in caplog.messages)


def test_mfsd_report_keys_rows_by_the_parquet_record():
    """The fragment-size Parquet is keyed by the variant as genotyped (the VCF
    record for VCF input; the MAF alleles, Allele1 fallback included, for MAF
    input); the report's MAF lookup must use the same key."""
    from gbcms.report.mfsd_report import _maf_row_key

    for k, (p, r, a) in VCF_INPUT.items():
        row = MafWriter.vcf_input_fields(CoordinateKernel.vcf_to_internal("1", p, r, a))
        assert _maf_row_key(row) == f"1:{p}:{r}:{a}", k
    maf_row = {
        "Chromosome": "1",
        "Start_Position": "221",
        "Reference_Allele": "AG",
        "Tumor_Seq_Allele1": "CT",
        "Tumor_Seq_Allele2": "AG",
    }
    assert _maf_row_key(maf_row) == "1:221:AG:CT"


# ── gbcms convert: input checks ──────────────────────────────────────────


def test_convert_checks_its_inputs(tmp_path):
    (tmp_path / "s.maf").write_text(_maf_text())
    (tmp_path / "s.vcf").write_text(_vcf_text())
    missing = runner.invoke(
        app, ["convert", "-v", str(tmp_path / "none.maf"), "-o", str(tmp_path / "o.vcf")]
    )
    assert missing.exit_code == 2 and missing.exception is not None
    unindexed = tmp_path / "raw.fa"
    unindexed.write_text(">1\n" + REF_SEQ + "\n")
    res = runner.invoke(
        app,
        [
            "convert",
            "-v",
            str(tmp_path / "s.maf"),
            "-f",
            str(unindexed),
            "-o",
            str(tmp_path / "o.vcf"),
        ],
    )
    assert res.exit_code == 1 and not (tmp_path / "o.vcf").exists()
    assert not (tmp_path / "raw.fa.fai").exists()
    fasta = _fasta(tmp_path)
    ok = runner.invoke(
        app,
        ["convert", "-v", str(tmp_path / "s.vcf"), "-f", str(fasta), "-o", str(tmp_path / "o.maf")],
    )
    assert ok.exit_code == 0
    assert "--fasta is not used for VCF input" in ok.output
