"""Group 6 statistics contracts (S1 #153, S2 #154).

Operator decisions 2026-09-25 (S1: the LLR as a mean per fragment; S2: CH-LIKE must
not rest on a non-significant KS), refined 2026-10-05: mFSD is graded evidence
(leans somatic / no size evidence / insufficient), never a hard origin call; gene
membership is a note, not a gate; the KS p-value is exact; an empty fragment class
has no mean; ``mfsd_alt_confidence`` names how much data there is.
Measured on labeled ACCESS plasma (``CYCLE_6.6.0_PLAN.md`` § Group 6).
"""

import math
import re
from pathlib import Path

import pysam
import pytest

from gbcms._rs import Variant, count_bam_binned
from gbcms.io.output import MafWriter, VcfWriter
from gbcms.report.mfsd_report import _classify_origin

ROOT = Path(__file__).resolve().parent.parent

# ── a synthetic locus with chosen fragment sizes ────────────────────────────────


def _bam(tmp_path, ref_sizes, alt_sizes):
    """chr1:100 (0-based; REF A, ALT T): one proper pair per fragment, both mates
    covering the site, TLEN = the fragment size (no indels, so the physical size)."""
    path = tmp_path / "frag.bam"
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 5000, "SN": "chr1"}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for kind, sizes, base in (("r", ref_sizes, "A"), ("a", alt_sizes, "T")):
            for i, size in enumerate(sizes):
                for is_r2 in (False, True):
                    a = pysam.AlignedSegment()
                    a.query_name = f"{kind}{i}"
                    a.query_sequence = "AAAAA" + base + "AAAA"
                    a.flag = (1 | 2 | 128 | 16) if is_r2 else (1 | 2 | 64)
                    a.reference_id = 0
                    a.reference_start = 95
                    a.mapping_quality = 60
                    a.cigartuples = [(0, 10)]
                    a.query_qualities = [30] * 10  # type: ignore[assignment]
                    a.next_reference_id = 0
                    a.next_reference_start = 95
                    a.template_length = -size if is_r2 else size
                    out.write(a)
    sorted_path = tmp_path / "frag.sorted.bam"
    pysam.sort("-o", str(sorted_path), str(path))
    pysam.index(str(sorted_path))
    return str(sorted_path)


def _count(bam):
    variant = Variant(chrom="chr1", pos=100, ref_allele="A", alt_allele="T", variant_type="SNP")
    return count_bam_binned(
        bam, [variant], [None], min_mapq=20, min_baseq=20, filter_duplicates=True,
        filter_secondary=True, filter_supplementary=True, filter_qc_failed=False,
        filter_improper_pair=False, filter_indel=False, threads=1, mfsd=True,
    )[0]  # fmt: skip


def _llr_term(x):
    """One fragment's log(P_tumor/P_healthy) under gbcms's two Gaussians."""
    zt, zh = (x - 145.0) / 35.0, (x - 167.0) / 30.0
    return 0.5 * (zh * zh - zt * zt) + math.log(30.0 / 35.0)


def _ks_exact(a, b):
    """Reference two-sample KS: D from the tie-aware ECDF walk, p from the share of
    monotone lattice paths that stay strictly inside the band (exact at any n, m)."""
    a, b = sorted(a), sorted(b)
    n, m = len(a), len(b)
    i = j = 0
    d = 0.0
    while i < n and j < m:
        v = min(a[i], b[j])
        while i < n and a[i] <= v:
            i += 1
        while j < m and b[j] <= v:
            j += 1
        d = max(d, abs(i / n - j / m))
    band = d * n * m - 1e-9
    prev = [0.0] * (m + 1)
    prev[0] = 1.0
    for jj in range(1, m + 1):
        prev[jj] = prev[jj - 1] if jj * n < band else 0.0
    for ii in range(1, n + 1):
        cur = [0.0] * (m + 1)
        cur[0] = prev[0] if ii * m < band else 0.0
        for jj in range(1, m + 1):
            if abs(ii * m - jj * n) < band:
                cur[jj] = prev[jj] * ii / (ii + jj) + cur[jj - 1] * jj / (ii + jj)
        prev = cur
    return d, min(max(1.0 - prev[m], 0.0), 1.0)


_REF = [150 + (k * 7) % 60 for k in range(40)]  # 150-209 bp, a mono-nucleosome-like spread


# ── S1: the LLR is a mean per fragment; an empty class has no mean ──────────────


def test_the_alt_and_ref_llr_are_the_mean_per_fragment(tmp_path):
    alt = [118, 126, 133, 141, 152]
    c = _count(_bam(tmp_path, _REF, alt))
    assert c.mfsd_alt_count == len(alt) and c.mfsd_ref_count == len(_REF)
    assert c.mfsd_alt_llr == pytest.approx(sum(map(_llr_term, alt)) / len(alt), rel=1e-12)
    assert c.mfsd_ref_llr == pytest.approx(sum(map(_llr_term, _REF)) / len(_REF), rel=1e-12)


def test_an_empty_class_has_no_mean_size_or_llr(tmp_path):
    c = _count(_bam(tmp_path, _REF, []))
    assert c.mfsd_alt_count == 0
    assert math.isnan(c.mfsd_alt_mean) and math.isnan(c.mfsd_alt_llr)
    assert math.isnan(c.mfsd_n_mean) and math.isnan(
        c.mfsd_nonref_mean
    )  # no N / other fragments here
    assert c.mfsd_ref_mean == pytest.approx(sum(_REF) / len(_REF))


def test_the_vcf_header_describes_the_llr_as_a_mean(tmp_path):
    path = tmp_path / "out.vcf"
    writer = VcfWriter(path, sample_name="TUMOR", mfsd=True)
    writer._write_header()
    writer.close()
    lines = [
        ln
        for ln in path.read_text().splitlines()
        if re.match(r"##INFO=<ID=MFSD_(ALT|REF)_LLR,", ln)
    ]
    assert len(lines) == 2
    assert all("mean per fragment" in ln and "Σ" not in ln for ln in lines), lines


# ── S2: the KS p-value is exact where it decides ───────────────────────────────


def test_the_ks_p_value_is_exact_where_it_already_was(tmp_path):
    alt = [118, 126, 133, 141, 152]
    c = _count(_bam(tmp_path, _REF, alt))  # 5 x 40 = 200 lattice cells
    d, p = _ks_exact(alt, _REF)
    assert c.mfsd_ks_alt_ref == pytest.approx(d)
    assert c.mfsd_pval_alt_ref == pytest.approx(p, rel=1e-9)


def test_the_ks_p_value_is_exact_for_few_alt_fragments_against_a_deep_ref(tmp_path):
    # 5 ALT vs 2,400 REF fragments: 12,000 lattice cells, past the old exact limit,
    # where the asymptotic series overstated p about 1.7-2.3x (measured).
    ref = [140 + (k * 37) % 80 for k in range(2400)]
    # exact p 0.0172 here; the asymptotic series gives 0.0346.
    alt = [143, 150, 156, 161, 168]
    c = _count(_bam(tmp_path, ref, alt))
    d, p = _ks_exact(alt, ref)
    assert c.mfsd_ks_alt_ref == pytest.approx(d)
    assert c.mfsd_pval_alt_ref == pytest.approx(p, rel=1e-6)


# ── S2: mfsd_alt_confidence names how much data there is ───────────────────────


@pytest.mark.parametrize(
    "alt_count,tier",
    [
        (0, "NONE"),
        (1, "SPARSE"),
        (4, "SPARSE"),
        (5, "TESTABLE"),
    ],
)
def test_alt_confidence_names_the_data_not_a_verdict(tmp_path, alt_count, tier):
    from test_mfsd_flag import _MockCounts, _MockVariant

    counts = _MockCounts(mfsd=True)
    counts.mfsd_alt_count = alt_count
    path = tmp_path / "out.maf"
    writer = MafWriter(path, mfsd=True)
    writer.write(_MockVariant(), counts)
    writer.close()
    header, row = [
        ln.split("\t") for ln in path.read_text().splitlines() if not ln.startswith("#")
    ][:2]
    assert dict(zip(header, row, strict=True))["mfsd_alt_confidence"] == tier


# ── S2: graded classes; the CH gene is a note, not a gate ──────────────────────

_NAN = float("nan")


def test_a_tumor_tp53_variant_with_shorter_fragments_leans_somatic():
    # TP53 is in the CH gene set; that no longer blocks the somatic direction.
    signal, reason = _classify_origin("TP53", 2.0, 0.2, 0.4, 0.001, True, 20, 3)
    assert signal == "LEANS-SOMATIC"
    assert "CH-associated gene" in reason


def test_a_ch_gene_without_a_size_shift_has_no_size_evidence():
    # Today's CH-LIKE: a CH gene, low enrichment, a non-significant KS. Absence of a
    # shift is not evidence for CH (about 30% of tumor variants look REF-like below
    # 50 fragments).
    signal, reason = _classify_origin("DNMT3A", 1.0, 0.2, 0.2, 0.6, True, 10, 3)
    assert signal == "NO-SIZE-EVIDENCE"
    assert "CH-associated gene" in reason


def test_a_significant_shortening_leans_somatic_without_a_further_enrichment_gate():
    signal, _ = _classify_origin("NOTACHGENE", 1.15, 0.20, 0.23, 0.01, True, 30, 3)
    assert signal == "LEANS-SOMATIC"


def test_a_significantly_longer_alt_is_not_somatic_evidence():
    signal, reason = _classify_origin("NOTACHGENE", 0.5, 0.30, 0.15, 0.001, True, 30, 3)
    assert signal == "NO-SIZE-EVIDENCE"
    assert "longer" in reason


def test_no_report_class_is_named_ch_like(tmp_path):
    from gbcms.report import generate_mfsd_report

    testdata = Path(__file__).parent / "testdata"
    out = generate_mfsd_report(
        parquet_path=testdata / "mfsd_multi.parquet", maf_path=testdata / "mfsd_multi.maf",
        output_path=tmp_path / "r.html", min_alt=3, max_variants=20, sample_name="TEST-SAMPLE",
    )  # fmt: skip
    html = out.read_text()
    assert "CH-LIKE" not in html and "CH-like" not in html
    # the summary cards name the graded classes (the badge itself is checked on a
    # report built from known sizes, below)
    assert '<div class="label">Leans somatic</div>' in html
    assert '<div class="label">No size evidence</div>' in html


# ── docs follow the code ───────────────────────────────────────────────────────


def test_the_docs_describe_the_graded_classes_and_the_mean_llr():
    qc = (ROOT / "docs" / "reference" / "qc-flags.md").read_text()
    for term in ("LEANS-SOMATIC", "NO-SIZE-EVIDENCE", "INSUFFICIENT", "TESTABLE", "SPARSE"):
        assert term in qc, term
    stale = [
        p.relative_to(ROOT) for p in (ROOT / "docs").rglob("*.md") if "CH-LIKE" in p.read_text()
    ]
    assert not stale, stale
    metrics = (ROOT / "docs" / "reference" / "counting-metrics.md").read_text()
    row = next(ln for ln in metrics.splitlines() if ln.startswith("| `mfsd_alt_llr`"))
    assert "mean" in row and "Σ" not in row, row


# ── review round: direction, skipped contigs, writers ──────────────────────────


def test_the_direction_comes_from_the_ks_gap_not_the_short_fragment_share():
    # Shares can be equal (both 0) while ALT is clearly shorter, or tilt the other way
    # while a long ALT tail drives the test.
    signal, _ = _classify_origin("NOTACHGENE", _NAN, 0.0, 0.0, 0.0, True, 30, 3, alt_shorter=True)
    assert signal == "LEANS-SOMATIC"
    signal, reason = _classify_origin(
        "NOTACHGENE", 1.1, 0.10, 0.11, 1e-4, True, 30, 3, alt_shorter=False
    )
    assert signal == "NO-SIZE-EVIDENCE"
    assert "longer" in reason


def test_the_report_reads_the_direction_from_the_fragment_sizes(tmp_path):
    import polars as pl

    from gbcms.report import generate_mfsd_report

    # REF 165-204 bp and ALT 150-157 bp: ALT about 31 bp shorter, yet no fragment on
    # either side is under 150 bp, so the sub-nucleosomal shares are both 0.
    ref = [165 + k % 40 for k in range(200)]
    alt = [150 + k % 8 for k in range(30)]
    pl.DataFrame(
        {
            "chrom": ["1"],
            "pos": [1000],
            "ref": ["C"],
            "alt": ["T"],
            "ref_sizes": [ref],
            "alt_sizes": [alt],
        }
    ).write_parquet(tmp_path / "s.fsd.parquet")
    cols = {
        "Hugo_Symbol": "NOTACHGENE", "Chromosome": "1", "Start_Position": "1000",
        "Reference_Allele": "C", "Tumor_Seq_Allele1": "C", "Tumor_Seq_Allele2": "T",
        "mfsd_sub_nuc_enrichment": "NA", "mfsd_pval_alt_ref": "0.0000", "mfsd_qval_alt_ref": "0.0000",
        "mfsd_ks_valid": "True", "mfsd_alt_mean": "153.5", "mfsd_ref_mean": "184.5",
        "mfsd_delta_alt_ref": "-31.0", "mfsd_alt_llr": "0.4", "mfsd_sub_nuc_ref_frac": "0.0000",
        "mfsd_sub_nuc_alt_frac": "0.0000",
    }  # fmt: skip
    (tmp_path / "s.maf").write_text("\t".join(cols) + "\n" + "\t".join(cols.values()) + "\n")
    out = generate_mfsd_report(
        parquet_path=tmp_path / "s.fsd.parquet", maf_path=tmp_path / "s.maf",
        output_path=tmp_path / "r.html", min_alt=3, max_variants=20, sample_name="T",
    )  # fmt: skip
    html = out.read_text()
    assert ">LEANS-SOMATIC</span>" in html
    assert "ALT longer than REF" not in html


def _bam_loci(tmp_path, loci):
    """chr1 loci at the given 0-based sites (REF A, ALT T), each with its fragments."""
    path = tmp_path / "loci.bam"
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 50000, "SN": "chr1"}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for site, ref_sizes, alt_sizes in loci:
            for kind, sizes, base in (("r", ref_sizes, "A"), ("a", alt_sizes, "T")):
                for i, size in enumerate(sizes):
                    for is_r2 in (False, True):
                        a = pysam.AlignedSegment()
                        a.query_name = f"s{site}{kind}{i}"
                        a.query_sequence = "AAAAA" + base + "AAAA"
                        a.flag = (1 | 2 | 128 | 16) if is_r2 else (1 | 2 | 64)
                        a.reference_id = 0
                        a.reference_start = site - 5
                        a.mapping_quality = 60
                        a.cigartuples = [(0, 10)]
                        a.query_qualities = [30] * 10  # type: ignore[assignment]
                        a.next_reference_id = 0
                        a.next_reference_start = site - 5
                        a.template_length = -size if is_r2 else size
                        out.write(a)
    sorted_path = tmp_path / "loci.sorted.bam"
    pysam.sort("-o", str(sorted_path), str(path))
    pysam.index(str(sorted_path))
    return str(sorted_path)


def _count_all(bam, variants):
    return count_bam_binned(
        bam, variants, [None] * len(variants), min_mapq=20, min_baseq=20, filter_duplicates=True,
        filter_secondary=True, filter_supplementary=True, filter_qc_failed=False,
        filter_improper_pair=False, filter_indel=False, threads=1, mfsd=True,
    )  # fmt: skip


def test_rows_on_a_contig_absent_from_the_bam_stay_out_of_the_mfsd_family(tmp_path):
    bam = _bam_loci(
        tmp_path,
        [(100, _REF, [118, 126, 133, 141, 152, 120]), (2000, _REF, [150, 160, 170, 180, 190, 155])],
    )

    def snv(chrom, pos):
        return Variant(chrom=chrom, pos=pos, ref_allele="A", alt_allele="T", variant_type="SNP")

    real = [snv("chr1", 100), snv("chr1", 2000)]
    alone = _count_all(bam, real)
    mixed = _count_all(bam, real + [snv("chr9", 100 + k) for k in range(8)])
    for skipped in mixed[2:]:
        assert math.isnan(skipped.mfsd_ks_alt_ref)  # no test ran, so not in the BH family
        assert math.isnan(skipped.mfsd_alt_mean) and math.isnan(skipped.mfsd_alt_llr)
    for a, b in zip(alone, mixed[:2], strict=True):
        assert b.mfsd_qval_alt_ref == pytest.approx(a.mfsd_qval_alt_ref)


def test_the_writers_write_an_empty_class_as_missing(tmp_path):
    from test_mfsd_flag import _MockCounts, _MockVariant

    counts = _MockCounts(mfsd=True)
    counts.mfsd_alt_count = 0
    counts.mfsd_alt_mean = counts.mfsd_alt_llr = float("nan")
    maf, vcf = tmp_path / "o.maf", tmp_path / "o.vcf"
    for writer in (MafWriter(maf, mfsd=True), VcfWriter(vcf, sample_name="T", mfsd=True)):
        writer.write(_MockVariant(), counts)
        writer.close()
    header, row = [ln.split("\t") for ln in maf.read_text().splitlines() if not ln.startswith("#")][
        :2
    ]
    cells = dict(zip(header, row, strict=True))
    assert cells["mfsd_alt_mean"] == "NA" and cells["mfsd_alt_llr"] == "NA"
    info = [ln for ln in vcf.read_text().splitlines() if not ln.startswith("#")][0].split("\t")[7]
    assert "MFSD_ALT_LLR=." in info.split(";")
