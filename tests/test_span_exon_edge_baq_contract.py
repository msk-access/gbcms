"""The exon-edge BAQ rule measures from the variant's whole REF span (#106).

Heuristic BAQ lowers base quality by 20 within 5 read bases of a CIGAR N, so for
an exon ending at E, bases E-5..E-1 of every read spliced at E are penalized. At a
variant within 5bp of an annotated exon boundary BAQ is skipped: the penalty
would land on exactly the reads that splice there. The distance was measured from
the variant's first base, so a multi-base variant starting more than 5bp from a
right exon edge but reaching into those last 5 bases kept BAQ:

- its edge-side bases were masked in every spliced read, so no spliced carrier
  could show the whole MNP (`mnp_confirmed_alt` 0), and `--rescue-mnp` adopted a
  component where the reads show the annotated haplotype;
- a spliced read carrying only the edge-side change voted on its unmasked bases,
  which match REF, and was counted REF;
- the MNP row and its component re-count at the edge base ran under two BAQ
  rules.

The distance (and the `exon_boundary_dist` column) is now the distance from the
nearest annotated boundary to any base of the REF span, 0 when a boundary lies
inside it (Ensembl VEP's overlap view). SNVs and insertions (one-base REF) are
unchanged, and BAQ still applies to multi-base variants away from edges.

Committed red (xfail-strict) before the implementation; flipped green with it.
"""

import pytest
from helpers import make_read
from rna_fixtures import (
    E1,
    E2,
    READ_LEN,
    SENSE,
    mk_ref,
    run_rna,
    spliced,
    write_bam,
    write_fasta,
    write_gtf,
    write_maf,
    write_vcf,
)

from gbcms import _rs as gbcms_rs

DONOR = E1[1]  # E1's right edge (exclusive end): reads spliced here lose DONOR-5..DONOR-1
MNP4 = DONOR - 7  # 4bp MNP over DONOR-7..DONOR-4: first base 7bp out, last base 4bp
DNP = DONOR - 6  # 2bp MNP over DONOR-6..DONOR-5: first base 6bp out, last base 5bp
TX = "T1"
_REF = mk_ref()


def _mutate(ref, changes):
    """`ref` with each {pos: base} applied."""
    seq = list(ref)
    for p, b in changes.items():
        seq[p] = b
    return "".join(seq)


def _other(base):
    return "A" if base != "A" else "C"


def _mnp(pos, n):
    """(REF, ALT, {pos: alt base}) for an n-base MNP at `pos` changing every base."""
    ref = _REF[pos : pos + n]
    alt = "".join(_other(b) for b in ref)
    return ref, alt, {pos + i: alt[i] for i in range(n)}


def _unspliced(hap, pos, n, prefix):
    """M-only reads over `pos` (pre-mRNA / retained-intron reads): no N, so no
    splice-junction penalty."""
    out = []
    for i in range(n):
        s = pos - 50 + (i % 5)
        out.append(
            make_read(f"{prefix}{i}", hap[s : s + READ_LEN], s, ((0, READ_LEN),), flag=SENSE)
        )
    return out


def _with_quals(reads, q):
    for r in reads:
        r.query_qualities = [q] * r.query_length
    return reads


def _tx(row, col="transcript_read_counts"):
    """{transcript: (AD, RD, DP)} from `ENST:AD,RD,DP|...`."""
    out = {}
    for entry in filter(None, row[col].split("|")):
        tx, nums = entry.rsplit(":", 1)
        out[tx] = tuple(int(x) for x in nums.split(","))
    return out


def _audit(row):
    return dict(part.split("=", 1) for part in row["gbcms_rescue"].split(";") if "=" in part)


def _run(tmp_path, variants, reads, extra=()):
    """CLI RNA run (RNA defaults: BAQ, strandedness, --gtf) of VCF rows
    (0-based pos, REF, ALT); returns the MAF rows."""
    return run_rna(
        tmp_path,
        write_vcf(tmp_path, [(p + 1, r, a) for p, r, a in variants]),
        write_bam(tmp_path, _REF, reads),
        write_fasta(tmp_path, _REF),
        write_gtf(tmp_path),
        extra=extra,
    )


def _engine(tmp_path, pos, ref, alt, reads):
    """Binned engine, RNA mode, BAQ on, GTF given: the BaseCounts (for the
    internal `mnp_confirmed_alt`)."""
    fa = write_fasta(tmp_path, _REF, contig="1")
    bam = write_bam(tmp_path, _REF, reads)
    gtf = write_gtf(tmp_path, contig="1")
    (pv,) = gbcms_rs.prepare_variants(
        [gbcms_rs.Variant("1", pos, ref, alt, "MNP")], str(fa), 5, False, 1, True
    )
    assert pv.gbcms_status == "PASS", pv.gbcms_status_reason
    (c,) = gbcms_rs.count_bam_binned(
        str(bam),
        [pv.variant],
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
        mode="rna",
        apply_baq=True,
        gtf_path=str(gtf),
    )
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    return c


# ── Spliced carriers of an MNP reaching a right exon edge show the haplotype ──
@pytest.mark.parametrize("q", [30, 37, 40])  # Q40 - 20 passes min BQ 20: never masked
def test_spliced_carriers_of_an_edge_mnp_show_the_whole_haplotype(tmp_path, q):
    ref, alt, change = _mnp(MNP4, 4)
    reads = spliced(_REF, DONOR, 500, 30, "ref") + spliced(
        _mutate(_REF, change), DONOR, 500, 10, "alt"
    )
    c = _engine(tmp_path, MNP4, ref, alt, _with_quals(reads, q))
    assert (c.rd, c.ad) == (30, 10)
    assert c.mnp_confirmed_alt == 10, "every changed base read in every spliced carrier"


# ── A spliced read carrying only the edge-side change is not REF ──────────────
def test_spliced_reads_with_only_the_edge_change_are_not_counted_ref(tmp_path):
    """DNP at DONOR-6..DONOR-5; ten spliced reads carry only DONOR-5's change.
    With that base masked they voted REF on DONOR-6; read in full they are
    partial ALT, in the main counts, per-transcript and ASJD alike."""
    ref, alt, change = _mnp(DNP, 2)
    edge_only = _mutate(_REF, {DNP + 1: change[DNP + 1]})
    reads = spliced(_REF, DONOR, 500, 30, "ref") + spliced(edge_only, DONOR, 500, 10, "edge")
    (row,) = _run(tmp_path, [(DNP, ref, alt)], reads)
    assert (int(row["ref_count"]), int(row["alt_count"]), int(row["partial_alt"])) == (30, 0, 10)
    assert _tx(row)[TX][:2] == (0, 30)
    assert int(row["asjd_n_ref_total"]) == 30
    assert row["exon_boundary_dist"] == str(DONOR - (DNP + 1))


# ── --rescue-mnp: one rule for the MNP row and its component re-count ─────────
def test_rescue_keeps_an_edge_mnp_its_spliced_reads_show(tmp_path):
    """Six spliced reads carry the whole DNP; ten unspliced reads carry only
    DONOR-6's change. The BAM shows the annotated haplotype, so rescue keeps it
    (it adopted the DONOR-6 component: the spliced carriers' edge base was
    masked, so none showed the whole MNP)."""
    ref, alt, change = _mnp(DNP, 2)
    hap = _mutate(_REF, change)
    far_only = _mutate(_REF, {DNP: change[DNP]})
    reads = (
        spliced(_REF, DONOR, 500, 30, "ref")
        + spliced(hap, DONOR, 500, 6, "full")
        + _unspliced(far_only, DNP, 10, "far")
    )
    (row,) = _run(tmp_path, [(DNP, ref, alt)], reads, extra=("--rescue-mnp",))
    audit = _audit(row)
    assert audit["outcome"] == "haplotype_confirmed", row["gbcms_rescue"]
    assert audit["original_confirmed"] == "6"
    assert (int(row["alt_count"]), int(row["partial_alt"])) == (6, 10)


def test_rescue_counts_an_edge_mnp_and_its_component_under_one_rule(tmp_path):
    """Only DONOR-5's change is carried, by ten spliced reads. The MNP row sees
    them as partial ALT (it counted them REF), so rescue re-counts the
    components; the DONOR-5 component, counted under the same rule, holds the
    same ten reads."""
    ref, alt, change = _mnp(DNP, 2)
    edge_only = _mutate(_REF, {DNP + 1: change[DNP + 1]})
    reads = spliced(_REF, DONOR, 500, 30, "ref") + spliced(edge_only, DONOR, 500, 10, "edge")
    (row,) = _run(tmp_path, [(DNP, ref, alt)], reads, extra=("--rescue-mnp",))
    audit = _audit(row)
    assert audit.get("outcome") == "rescued", row["gbcms_rescue"]
    assert audit["adopted"] == f"chr1:{DNP + 2}({ref[1]}>{alt[1]})"
    assert (audit["original_ref"], audit["original_partial"]) == ("30", "10")
    assert (int(row["ref_count"]), int(row["alt_count"])) == (30, 10)


def test_rescue_counts_a_far_component_under_the_mnps_rule(tmp_path):
    """Ten reads splice at a cryptic donor 3bp before the annotated one and carry
    only DONOR-6's change. The MNP row (span distance 5: no BAQ) holds them as
    partial ALT. The DONOR-6 component is re-counted under the MNP row's rule, so
    it holds the same ten reads (on its own distance, 6, BAQ masked that base in
    every one of them, and rescue found nothing to adopt)."""
    ref, alt, change = _mnp(DNP, 2)
    far_only = _mutate(_REF, {DNP: change[DNP]})
    reads = spliced(_REF, DONOR, 500, 30, "ref") + spliced(far_only, DONOR - 3, 500, 10, "cryptic")
    (row,) = _run(tmp_path, [(DNP, ref, alt)], reads, extra=("--rescue-mnp",))
    audit = _audit(row)
    assert audit.get("outcome") == "rescued", row["gbcms_rescue"]
    assert audit["adopted"] == f"chr1:{DNP + 1}({ref[0]}>{alt[0]})"
    assert audit["original_partial"] == "10"
    assert int(row["alt_count"]) == 10


def test_a_rescued_row_keeps_the_mnps_distance(tmp_path):
    """Unspliced reads carry only DONOR-6's change, and rescue adopts that
    component's counts. The row keeps the MNP's coordinates, so it reports the
    MNP's distance (5), not the component's own (6)."""
    ref, alt, change = _mnp(DNP, 2)
    far_only = _mutate(_REF, {DNP: change[DNP]})
    reads = spliced(_REF, DONOR, 500, 30, "ref") + _unspliced(far_only, DNP, 10, "far")
    (row,) = _run(tmp_path, [(DNP, ref, alt)], reads, extra=("--rescue-mnp",))
    assert _audit(row).get("outcome") == "rescued", row["gbcms_rescue"]
    assert row["exon_boundary_dist"] == "5"


# ── exon_boundary_dist is the span distance ──────────────────────────────────
# (label, first REF base, end of REF, expected distance). A deletion's REF
# starts at its VCF anchor base; the rest are MNPs changing every base.
SPANS = [
    # A 4bp MNP over DONOR-7..DONOR-4: its last base is 4bp from the edge.
    ("MNP into a right edge", MNP4, MNP4 + 4, DONOR - (MNP4 + 3)),
    # A deletion from 10bp inside E1 to 5bp into the intron.
    ("deletion across a right edge", DONOR - 10, DONOR + 6, 0),
    # A deletion from 8bp inside intron 1 to 3bp into E2.
    ("deletion across a left edge", E2[0] - 8, E2[0] + 3, 0),
    # Guard: a span starting 6bp into an exon is 6bp from its left edge.
    ("MNP inside a left edge", E2[0] + 6, E2[0] + 8, 6),
]


def _span_variant(label, pos, end):
    ref = _REF[pos:end]
    alt = ref[0] if label.startswith("deletion") else "".join(_other(b) for b in ref)
    return pos, ref, alt


@pytest.mark.parametrize("label,pos,end,expected", SPANS, ids=[s[0] for s in SPANS])
def test_the_distance_is_measured_from_the_ref_span(tmp_path, label, pos, end, expected):
    (row,) = _run(tmp_path, [_span_variant(label, pos, end)], spliced(_REF, DONOR, 500, 20, "ref"))
    assert row["exon_boundary_dist"] == str(expected)


def test_maf_input_gives_the_vcf_distance(tmp_path):
    """The same variants as MAF rows (a deletion's Start is its first deleted
    base) are prepared to the same REF span, so they report the same distance."""
    variants = [_span_variant(label, pos, end) for label, pos, end, _ in SPANS]
    rows = [(p + 1, r, a) for p, r, a in variants]
    reads = spliced(_REF, DONOR, 500, 20, "ref")
    common = (write_bam(tmp_path, _REF, reads), write_fasta(tmp_path, _REF), write_gtf(tmp_path))
    got = [
        [r["exon_boundary_dist"] for r in run_rna(tmp_path, v, *common, outname=name)]
        for name, v in (("vcf", write_vcf(tmp_path, rows)), ("maf", write_maf(tmp_path, rows)))
    ]
    assert got[0] == got[1] == [str(e) for *_, e in SPANS]


def test_one_base_ref_variants_keep_the_first_base_distance(tmp_path):
    """Guard: an SNV and an insertion have a one-base REF, so their distance is
    unchanged."""
    snv, ins = DONOR - 7, DONOR - 9
    variants = [
        (snv, _REF[snv], _other(_REF[snv])),
        (ins, _REF[ins], _REF[ins] + "GT"),
    ]
    rows = _run(tmp_path, variants, spliced(_REF, DONOR, 500, 20, "ref"))
    assert [r["exon_boundary_dist"] for r in rows] == ["7", "9"]


# ── Guard: BAQ still applies one base beyond the window ──────────────────────
def test_baq_still_applies_to_an_mnp_ending_six_bases_from_the_edge(tmp_path):
    """DNP at DONOR-7..DONOR-6 (span distance 6). Its carriers splice at a novel
    donor 3bp after it, whose N penalty masks both MNP bases (a junction's
    penalty is never spared), so they are not counted ALT, as before."""
    pos = DONOR - 7
    ref, alt, change = _mnp(pos, 2)
    carriers = spliced(_mutate(_REF, change), pos + 5, 500, 10, "novel")
    (row,) = _run(tmp_path, [(pos, ref, alt)], _unspliced(_REF, pos, 20, "ref") + carriers)
    assert (int(row["ref_count"]), int(row["alt_count"])) == (20, 0)
