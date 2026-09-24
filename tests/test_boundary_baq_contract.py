"""Target contract: one base-quality rule across the RNA counting views, a
deterministic ASJD junction choice, and an exon distance found in any contig
naming.

Heuristic BAQ lowers base qualities within 5bp of a CIGAR N. At a variant
within 5bp of an annotated exon boundary that penalizes exactly the reads that
splice there, so the main counts skip BAQ at such sites. Per-transcript counts
and ASJD read the same reads and must apply the same rule:

- at an exon-edge SNV, the per-transcript counts of the transcript every read
  is compatible with equal the main counts (they undercounted by ~3% on real
  RNA, up to all of a transcript's REF reads);
- ASJD sees the junction evidence there (it saw none: every spliced read was
  masked, so a splice-disrupting variant at the edge had no junction at all);
- away from a boundary BAQ still applies, in every view.

ASJD reports each partition's dominant junction. When junctions tie for the
most fragments the choice must not depend on hash order (it did: the same
input could report a divergence on one run and none on the next). A tie is not
a divergence: when the REF and ALT partitions' top junctions overlap, both
report the shared junction and no test is run.

Committed red (xfail-strict) before the implementation; flipped green with it.
"""

import pytest
from helpers import make_read
from rna_fixtures import (
    E1,
    E2,
    E3,
    READ_LEN,
    SENSE,
    mk_ref,
    run_rna,
    spliced,
    through,
    write_bam,
    write_fasta,
    write_gtf,
    write_vcf,
)

TX = "T1"  # the GTF's one transcript: E1-E2-E3
EDGE_SNV = E1[1] - 2  # 0-based; 2bp inside E1's donor edge
MID_SNV = E1[1] - 20  # 0-based; 20bp from any edge, still covered by E1->x reads
INTERIOR_SNV = 200  # 0-based; mid-E1, ~100bp from any edge
REPEATS = 12  # fresh engine runs per tie case: hash order varies per run


def _with(ref, pos, base):
    return ref[:pos] + base + ref[pos + 1 :]


def _alt(ref, pos):
    return "A" if ref[pos] != "A" else "C"


def _tx(row, col="transcript_read_counts"):
    """{transcript: (AD, RD, DP)} from `ENST:AD,RD,DP|...`; DP includes neither."""
    out = {}
    for entry in filter(None, row[col].split("|")):
        tx, nums = entry.rsplit(":", 1)
        ad, rd, dp = (int(x) for x in nums.split(","))
        assert dp >= ad + rd, entry
        out[tx] = (ad, rd, dp)
    return out


def _one(tmp_path, pos, reads, outname="out", names=None):
    """Run one SNV at `pos`; `names` overrides the contig name per source
    (keys: vcf, fasta, bam, gtf)."""
    names = {"vcf": "chr1", "fasta": "chr1", "bam": "1", "gtf": "chr1", **(names or {})}
    ref = _REF
    (row,) = run_rna(
        tmp_path,
        write_vcf(tmp_path, [(pos + 1, ref[pos], _alt(ref, pos))], contig=names["vcf"]),
        write_bam(tmp_path, ref, reads, name=f"{outname}.bam", contig=names["bam"]),
        write_fasta(tmp_path, ref, contig=names["fasta"]),
        write_gtf(tmp_path, contig=names["gtf"]),
        outname=outname,
    )
    return row


_REF = mk_ref()


# ── Per-transcript counts follow the main counts' BAQ rule ─────────────────
def test_per_transcript_matches_main_at_exon_edge(tmp_path):
    ref, alt = _REF, _alt(_REF, EDGE_SNV)
    reads = spliced(ref, E1[1], E2[0], 30, "ref") + spliced(
        _with(ref, EDGE_SNV, alt), E1[1], E2[0], 10, "alt"
    )
    row = _one(tmp_path, EDGE_SNV, reads)
    assert (int(row["ref_count"]), int(row["alt_count"])) == (30, 10), "main counts"
    ad, rd, _ = _tx(row)[TX]
    assert (ad, rd) == (10, 30)
    adf, rdf, _ = _tx(row, "transcript_fragment_counts")[TX]
    assert (adf, rdf) == (int(row["alt_count_fragment"]), int(row["ref_count_fragment"]))


def test_baq_still_applies_away_from_exon_edges(tmp_path):
    """Guard: mid-exon, a 1bp deletion 2bp from the SNV is inside BAQ's
    radius, so those reads' base at the SNV drops below --min-baseq in every
    view — the boundary rule must not switch BAQ off wholesale."""
    ref, pos = _REF, INTERIOR_SNV
    near_del = []
    for i in range(10):
        s = pos - 50 + (i % 5)
        left = pos + 2 - s
        seq = ref[s : pos + 2] + ref[pos + 3 : s + READ_LEN + 1]
        near_del.append(
            make_read(f"del{i}", seq, s, ((0, left), (2, 1), (0, READ_LEN - left)), flag=SENSE)
        )
    row = _one(tmp_path, pos, near_del + through(ref, pos, ref[pos], 20, "ok"))
    assert int(row["ref_count"]) == 20, "BAQ masks the deletion-adjacent reads in the main counts"
    ad, rd, _ = _tx(row)[TX]
    assert (ad, rd) == (0, 20)


def test_asjd_still_applies_baq_away_from_exon_edges(tmp_path):
    """Guard for ASJD: spliced reads with a 1bp deletion 2bp from a mid-exon
    SNV (20bp from the donor) lose that base to BAQ, so ASJD's REF tally
    counts only the clean spliced reads, as the main counts do."""
    ref, pos = _REF, MID_SNV
    near_del = []
    for i in range(10):
        s = E1[1] - 40 - (i % 5)
        left, mid = pos + 2 - s, E1[1] - (pos + 3)
        right = READ_LEN - left - mid
        seq = ref[s : pos + 2] + ref[pos + 3 : E1[1]] + ref[E2[0] : E2[0] + right]
        cigar = ((0, left), (2, 1), (0, mid), (3, E2[0] - E1[1]), (0, right))
        near_del.append(make_read(f"del{i}", seq, s, cigar, flag=SENSE))
    row = _one(tmp_path, pos, near_del + spliced(ref, E1[1], E2[0], 20, "ok"))
    assert int(row["ref_count"]) == 20
    assert int(row["asjd_n_ref_total"]) == 20


# ── ASJD sees junction evidence at exon edges ──────────────────────────────
def test_asjd_sees_junctions_at_exon_edge(tmp_path):
    """REF splices E1->E2, ALT skips E2 (E1->E3): the classic splice-disrupting
    geometry, at a variant 2bp from the donor."""
    ref, alt = _REF, _alt(_REF, EDGE_SNV)
    reads = spliced(ref, E1[1], E2[0], 30, "ref") + spliced(
        _with(ref, EDGE_SNV, alt), E1[1], E3[0], 10, "skip"
    )
    row = _one(tmp_path, EDGE_SNV, reads)
    assert row["asjd_ref_junction"] == f"{E1[1]}-{E2[0]}"
    assert row["asjd_alt_junction"] == f"{E1[1]}-{E3[0]}"
    assert (int(row["asjd_n_ref_junc"]), int(row["asjd_n_alt_junc"])) == (30, 10)
    assert row["asjd_pval"] != ""


# ── Dominant-junction ties: deterministic, and a tie is not a divergence ───
def _j(acceptor):
    return f"{E1[1]}-{acceptor}"


NEAR_E2 = E2[0] + 2  # an aligner's placement of E2's acceptor 2bp off (same junction)
IN_INTRON1 = 450  # an unrelated acceptor left of E2's (sorts first)
IN_INTRON2 = 700  # an unrelated acceptor between E2 and E3

# (REF groups, ALT groups) as (acceptor, fragments, prefix), then the expected
# (REF junction, ALT junction, divergence tested).
TIE_CASES = {
    # ALT tied between REF's junction and an exon skip: shared, no test.
    "alt_tie": (
        [(E2[0], 30, "ref")],
        [(E2[0], 6, "alt_inc"), (E3[0], 6, "alt_skip")],
        (_j(E2[0]), _j(E2[0]), False),
    ),
    # REF tied; ALT uses one of REF's tied junctions: shared, no test.
    "ref_tie": (
        [(E2[0], 15, "ref_inc"), (E3[0], 15, "ref_skip")],
        [(E3[0], 8, "alt_skip")],
        (_j(E3[0]), _j(E3[0]), False),
    ),
    # ALT tied between REF's junction placed 2bp off and an unrelated one that
    # sorts first: within tolerance it is REF's junction, so no test.
    "near_tie": (
        [(E2[0], 30, "ref")],
        [(IN_INTRON1, 6, "alt_far"), (NEAR_E2, 6, "alt_near")],
        (_j(E2[0]), _j(NEAR_E2), False),
    ),
    # ALT tied between two junctions REF does not use: the leftmost, tested.
    "alt_only_tie": (
        [(E2[0], 30, "ref")],
        [(IN_INTRON2, 6, "alt_a"), (E3[0], 6, "alt_b")],
        (_j(E2[0]), _j(IN_INTRON2), True),
    ),
}


@pytest.mark.parametrize("case", sorted(TIE_CASES))
def test_asjd_dominant_junction_ties(tmp_path, case):
    ref_groups, alt_groups, expected = TIE_CASES[case]
    ref, alt = _REF, _alt(_REF, MID_SNV)
    mutant = _with(ref, MID_SNV, alt)
    reads = [r for acc, n, p in ref_groups for r in spliced(ref, E1[1], acc, n, p)] + [
        r for acc, n, p in alt_groups for r in spliced(mutant, E1[1], acc, n, p)
    ]
    seen = set()
    for k in range(REPEATS):
        row = _one(tmp_path, MID_SNV, reads, outname=f"run{k}")
        seen.add((row["asjd_ref_junction"], row["asjd_alt_junction"], row["asjd_pval"] != ""))
    assert seen == {expected}, seen


# ── exon_boundary_dist: found in any contig naming, never a sentinel ──────
def test_exon_distance_found_for_mitochondrial_naming(tmp_path):
    """A chrM variant against an Ensembl-named (MT) GTF: the distance lookup
    reconciles the name as every other annotation lookup does, so the edge
    rule applies and the per-transcript counts match the main counts."""
    ref, alt = _REF, _alt(_REF, EDGE_SNV)
    reads = spliced(ref, E1[1], E2[0], 30, "ref") + spliced(
        _with(ref, EDGE_SNV, alt), E1[1], E2[0], 10, "alt"
    )
    row = _one(
        tmp_path, EDGE_SNV, reads, names={"vcf": "chrM", "fasta": "chrM", "bam": "MT", "gtf": "MT"}
    )
    assert row["exon_boundary_dist"] == str(E1[1] - EDGE_SNV)
    assert (int(row["ref_count"]), int(row["alt_count"])) == (30, 10)
    ad, rd, _ = _tx(row)[TX]
    assert (ad, rd) == (10, 30)


def test_exon_distance_empty_without_annotation(tmp_path):
    """A variant whose contig the GTF does not annotate has no distance: the
    column is empty, not a sentinel integer."""
    ref = _REF
    row = _one(
        tmp_path,
        INTERIOR_SNV,
        through(ref, INTERIOR_SNV, ref[INTERIOR_SNV], 20, "ok"),
        names={"gtf": "chr5"},
    )
    assert row["exon_boundary_dist"] == ""
