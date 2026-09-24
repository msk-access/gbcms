"""Target contract: one base-quality rule across the RNA counting views, and a
deterministic ASJD junction choice.

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
from test_asjd2_contract import (
    E1,
    E2,
    E3,
    READ_LEN,
    SENSE,
    _bam,
    _fasta,
    _gtf,
    _mk_ref,
    _run,
    _spliced,
    _through,
    _vcf,
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
    """{transcript: (AD, RD, DP)} from `ENST:AD,RD,DP|...`."""
    out = {}
    for entry in filter(None, row[col].split("|")):
        tx, nums = entry.rsplit(":", 1)
        out[tx] = tuple(int(x) for x in nums.split(","))
    return out


def _one(tmp_path, pos, reads, outname="out"):
    ref = _REF
    (row,) = _run(
        tmp_path,
        _vcf(tmp_path, [(pos + 1, ref[pos], _alt(ref, pos))]),
        _bam(tmp_path, ref, reads, name=f"{outname}.bam"),
        _fasta(tmp_path, ref),
        _gtf(tmp_path),
        outname=outname,
    )
    return row


_REF = _mk_ref()


# ── Per-transcript counts follow the main counts' BAQ rule ─────────────────
def test_per_transcript_matches_main_at_exon_edge(tmp_path):
    ref, alt = _REF, _alt(_REF, EDGE_SNV)
    reads = _spliced(ref, E1[1], E2[0], 30, "ref") + _spliced(
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
    row = _one(tmp_path, pos, near_del + _through(ref, pos, ref[pos], 20, "ok"))
    assert int(row["ref_count"]) == 20, "BAQ masks the deletion-adjacent reads in the main counts"
    ad, rd, _ = _tx(row)[TX]
    assert (ad, rd) == (0, 20)


# ── ASJD sees junction evidence at exon edges ──────────────────────────────
def test_asjd_sees_junctions_at_exon_edge(tmp_path):
    """REF splices E1->E2, ALT skips E2 (E1->E3): the classic splice-disrupting
    geometry, at a variant 2bp from the donor."""
    ref, alt = _REF, _alt(_REF, EDGE_SNV)
    reads = _spliced(ref, E1[1], E2[0], 30, "ref") + _spliced(
        _with(ref, EDGE_SNV, alt), E1[1], E3[0], 10, "skip"
    )
    row = _one(tmp_path, EDGE_SNV, reads)
    assert row["asjd_ref_junction"] == f"{E1[1]}-{E2[0]}"
    assert row["asjd_alt_junction"] == f"{E1[1]}-{E3[0]}"
    assert (int(row["asjd_n_ref_junc"]), int(row["asjd_n_alt_junc"])) == (30, 10)
    assert row["asjd_pval"] != ""


# ── A tie for the dominant junction is not a divergence ───────────────────
TIE_CASES = {
    # ALT fragments split evenly between REF's junction and an exon skip.
    "alt_tie": (
        [(E2[0], 30, "ref")],
        [(E2[0], 6, "alt_inc"), (E3[0], 6, "alt_skip")],
        f"{E1[1]}-{E2[0]}",
    ),
    # REF fragments split evenly; ALT uses one of REF's tied junctions.
    "ref_tie": (
        [(E2[0], 15, "ref_inc"), (E3[0], 15, "ref_skip")],
        [(E3[0], 8, "alt_skip")],
        f"{E1[1]}-{E3[0]}",
    ),
}


@pytest.mark.parametrize("case", sorted(TIE_CASES))
def test_asjd_tie_is_shared_junction_not_divergence(tmp_path, case):
    ref_groups, alt_groups, shared = TIE_CASES[case]
    ref, alt = _REF, _alt(_REF, MID_SNV)
    mutant = _with(ref, MID_SNV, alt)
    reads = [r for acc, n, p in ref_groups for r in _spliced(ref, E1[1], acc, n, p)] + [
        r for acc, n, p in alt_groups for r in _spliced(mutant, E1[1], acc, n, p)
    ]
    seen = set()
    for k in range(REPEATS):
        row = _one(tmp_path, MID_SNV, reads, outname=f"run{k}")
        seen.add((row["asjd_ref_junction"], row["asjd_alt_junction"], row["asjd_pval"]))
    assert seen == {(shared, shared, "")}, seen
