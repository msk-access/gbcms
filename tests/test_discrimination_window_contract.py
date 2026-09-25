"""Discrimination-window contract for indels (6.6.0 C10 #157, C2 #119).

An indel's discrimination window is the repeat tract holding its first changed
base (its own span in unique sequence), plus one base on each side. It is the
stretch a read must cover to show which allele it carries, as IGV shows it.

- C10: a read that starts or ends inside the window cannot tell the alleles
  apart (its bases match both, so the aligner places no gap). It counts toward
  depth but is neither REF nor ALT. GATK's AD likewise counts only informative
  reads.
- C2: at a co-annotated group, a read is excluded from a row's REF only when it
  carries a sibling's change inside that row's window, for read and fragment
  counts alike. A read carrying a sibling elsewhere still shows REF there.

"""

import glob
import random

import pysam
import pytest
from helpers import count_both, make_read, read_maf_output
from typer.testing import CliRunner

from gbcms import _rs as gbcms_rs
from gbcms.cli import app

runner = CliRunner()
READ = 100


def _ref(plants, n=700, seed=31):
    rng = random.Random(seed)
    ref = [rng.choice("ACGT") for _ in range(n)]
    for pos, motif in plants:
        ref[pos : pos + len(motif)] = list(motif)
    return "".join(ref)


def _files(tmp_path, ref, reads, contig="1"):
    fa = tmp_path / "ref.fa"
    fa.write_text(f">{contig}\n{ref}\n")
    pysam.faidx(str(fa))
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": contig, "LN": len(ref)}]}
    raw = tmp_path / "raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=hdr) as fh:
        for r in reads:
            fh.write(r)
    bam = tmp_path / "s.bam"
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    return str(fa), str(bam)


def _prepared(fa, variants):
    return [pv.variant for pv in gbcms_rs.prepare_variants(variants, fa, 5, False, 1, True)]


def _invariants(c):
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev


# ── C10: reads ending inside the repeat tract ────────────────────────────
# C | AAAAAAAAAA | G   (tract 0-based 200-209; anchor C at 199)
TRACT = ((199, "C" + "A" * 10 + "G"),)


def _deletion_reads(ref):
    alt = ref[:200] + ref[201:]  # one A deleted
    reads = []
    for i in range(10):  # REF, spanning the tract
        s = 120 + i
        reads.append(make_read(f"ref_span{i}", ref[s : s + READ], s, ((0, READ),)))
    for i in range(10):  # ALT, spanning; D placed at the tract start (left-aligned)
        s = 120 + i
        left = 200 - s
        reads.append(
            make_read(f"alt_span{i}", alt[s : s + READ], s, ((0, left), (2, 1), (0, READ - left)))
        )
    for i in range(10):  # REF reads ending inside the tract (last base 200..209)
        s = 101 + i
        reads.append(make_read(f"ref_part{i}", ref[s : s + READ], s, ((0, READ),)))
    for i in range(10):  # ALT carriers ending inside the tract: identical bases, no gap
        s = 101 + i
        reads.append(make_read(f"alt_part{i}", alt[s : s + READ], s, ((0, READ),)))
    return reads


def test_reads_ending_inside_a_homopolymer_deletion_are_uninformative(tmp_path):
    ref = _ref(TRACT)
    fa, bam = _files(tmp_path, ref, _deletion_reads(ref))
    (v,) = _prepared(fa, [gbcms_rs.Variant("1", 199, "CA", "C", "DELETION")])
    c = count_both(bam, [v])[0]
    _invariants(c)
    assert (c.rd, c.ad) == (10, 10)
    assert c.dp == 40  # the 20 uninformative reads stay in depth


def test_reads_ending_inside_a_homopolymer_insertion_are_uninformative(tmp_path):
    ref = _ref(TRACT)
    alt = ref[:200] + "A" + ref[200:]  # one A inserted
    reads = []
    for i in range(10):
        s = 120 + i
        reads.append(make_read(f"ref_span{i}", ref[s : s + READ], s, ((0, READ),)))
    for i in range(10):
        s = 120 + i
        left = 200 - s
        reads.append(
            make_read(
                f"alt_span{i}", alt[s : s + READ], s, ((0, left), (1, 1), (0, READ - left - 1))
            )
        )
    for i in range(10):
        s = 101 + i
        reads.append(make_read(f"ref_part{i}", ref[s : s + READ], s, ((0, READ),)))
        reads.append(make_read(f"alt_part{i}", alt[s : s + READ], s, ((0, READ),)))
    fa, bam = _files(tmp_path, ref, reads)
    (v,) = _prepared(fa, [gbcms_rs.Variant("1", 199, "C", "CA", "INSERTION")])
    c = count_both(bam, [v])[0]
    _invariants(c)
    assert (c.rd, c.ad) == (10, 10)
    assert c.dp == 40


def test_spanning_reads_in_unique_sequence_stay_ref(tmp_path):
    """Guard: a unique-sequence deletion; reads covering its span +1 are REF."""
    ref = _ref(())
    reads = [
        make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),))
        for i, s in enumerate(range(150, 160))
    ]
    fa, bam = _files(tmp_path, ref, reads)
    (v,) = _prepared(fa, [gbcms_rs.Variant("1", 199, ref[199:202], ref[199], "DELETION")])
    c = count_both(bam, [v])[0]
    _invariants(c)
    assert (c.rd, c.ad) == (10, 0)


def test_snv_counts_unchanged_by_the_window(tmp_path):
    """Guard: SNVs have no window; partial overlap is not a concept for them."""
    ref = _ref(TRACT)
    reads = [
        make_read(f"r{i}", ref[s : s + READ], s, ((0, READ),))
        for i, s in enumerate(range(106, 116))  # all cover 205 (inside the tract)
    ]
    fa, bam = _files(tmp_path, ref, reads)
    (v,) = _prepared(fa, [gbcms_rs.Variant("1", 205, "A", "C", "SNP")])
    c = count_both(bam, [v])[0]
    assert (c.rd, c.ad) == (10, 0)


# ── C2: sibling ALT inside vs outside the row's window ───────────────────
# Tract 1: C AAAAAAAA G at 199-208 (A-run 200-207). Tract 2: C TTTTTTTT G at
# 219-228 (T-run 220-227). The two 1bp deletions' scan windows overlap, so
# they form one tract-cluster group; their discrimination windows do not.
TWO_TRACTS = ((199, "C" + "A" * 8 + "G"), (219, "C" + "T" * 8 + "G"))


def _run_cli(tmp_path, ref, reads, rows):
    fa, bam = _files(tmp_path, ref, reads)
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=1,length=700>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(f"1\t{p}\t.\t{r}\t{a}\t.\t.\t.\n" for p, r, a in rows)
    )
    out = tmp_path / "out"
    res = runner.invoke(
        app, ["dna", "-v", str(vcf), "-b", f"S:{bam}", "-f", fa, "-o", str(out), "--format", "maf"]
    )
    assert res.exit_code == 0, res.output
    got = {int(r["vcf_pos"]): r for r in read_maf_output(glob.glob(str(out / "*.maf"))[0])}
    for r in got.values():
        assert int(r["total_count"]) >= int(r["ref_count"]) + int(r["alt_count"])
        assert int(r["total_count_fragment"]) >= int(r["ref_count_fragment"]) + int(
            r["alt_count_fragment"]
        )
    return got


def _span_read(name, hap, s, d_at=None):
    """A 100bp read from haplotype `hap` at `s`; a 1bp D at `d_at` when given."""
    if d_at is None:
        return make_read(name, hap[s : s + READ], s, ((0, READ),))
    left = d_at - s
    return make_read(name, hap[s : s + READ], s, ((0, left), (2, 1), (0, READ - left)))


def test_sibling_alt_outside_the_window_is_ref_here(tmp_path):
    ref = _ref(TWO_TRACTS)
    b_alt = ref[:220] + ref[221:]  # B: one T deleted from tract 2
    reads = [_span_read(f"wt{i}", ref, 150 + i) for i in range(6)]
    reads += [_span_read(f"b{i}", b_alt, 150 + i, d_at=220) for i in range(8)]  # show REF at A
    rows = [(200, "CA", "C"), (220, "CT", "C")]
    got = _run_cli(tmp_path, ref, reads, rows)
    a = got[200]
    assert "TRACT_CLUSTER" in a["gbcms_status_reason"]
    assert (int(a["ref_count"]), int(a["ref_count_fragment"])) == (14, 14)


def _twin_reads(ref, kind, snv_at=200):
    """6 wild-type reads plus 8 carriers of a sibling whose change sits in A's tract.

    - snv: an A>G at `snv_at` (no gap, so A's classifier sees REF).
    - ins: one extra A in the same homopolymer (the del/ins twin geometry).
    """
    reads = [_span_read(f"wt{i}", ref, 150 + i) for i in range(6)]
    for i in range(8):
        s = 150 + i
        if kind == "snv":
            hap = ref[:snv_at] + "G" + ref[snv_at + 1 :]
            reads.append(make_read(f"sib{i}", hap[s : s + READ], s, ((0, READ),)))
        else:
            hap = ref[:200] + "A" + ref[200:]
            left = 200 - s
            reads.append(
                make_read(
                    f"sib{i}", hap[s : s + READ], s, ((0, left), (1, 1), (0, READ - left - 1))
                )
            )
    return reads


# snv: A>G on the deleted base itself, so the two rows share a site (spans overlap).
SIBLING_ROWS = {"snv": (201, "A", "G"), "ins": (200, "C", "CA")}


@pytest.mark.parametrize("kind", ["snv", "ins"])
def test_sibling_alt_inside_the_window_is_not_ref_for_reads_or_fragments(tmp_path, kind):
    ref = _ref(TWO_TRACTS)
    rows = [(200, "CA", "C"), SIBLING_ROWS[kind]]
    _run_cli(tmp_path, ref, _twin_reads(ref, kind), rows)
    out = glob.glob(str(tmp_path / "out" / "*.maf"))[0]
    a = next(r for r in read_maf_output(out) if r["vcf_ref"] == "CA")
    assert "MULTI_ALLELIC" in a["gbcms_status_reason"]
    assert (int(a["ref_count"]), int(a["ref_count_fragment"])) == (6, 6)


def test_snv_outside_the_span_is_not_a_sibling_and_its_carriers_are_ref(tmp_path):
    """Guard: an SNV inside the tract but outside the deletion's VCF span is its
    own site. Its carriers show the reference tract length, so they are REF for
    the deletion (as GATK counts AD at the deletion's site)."""
    ref = _ref(TWO_TRACTS)
    rows = [(200, "CA", "C"), (204, "A", "G")]
    _run_cli(tmp_path, ref, _twin_reads(ref, "snv", snv_at=203), rows)
    out = glob.glob(str(tmp_path / "out" / "*.maf"))[0]
    a = next(r for r in read_maf_output(out) if r["vcf_ref"] == "CA")
    assert "MULTI_ALLELIC" not in a["gbcms_status_reason"]
    assert (int(a["ref_count"]), int(a["ref_count_fragment"])) == (14, 14)
