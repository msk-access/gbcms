"""RNA pins: decided behaviour where reads carry another allele, and an edit inside an insertion.

Tests only; no count changes with them. Each pins a decision already made (invariant 7:
count the given allele, judge bases, never infer what a wrong input should be).

- Per-transcript counts and ASJD: a read carrying another allele (counted in the main row's
  `partial_alt`) is in a transcript's DP and DPF, in no transcript AD or RD, and in no ASJD
  partition.
- Strandedness: under enforcement an antisense read carrying another allele appears in no
  column; with `--no-strandedness` it counts in `partial_alt`; with no GTF there is no gene
  strand, so enforcement is a no-op and the same read counts.
- Editing database: an entry cannot reach the bases of an insertion (SNV rows only), so the
  reads decide inside the windows (RJ-26): one edited carrier stays ALT, two or more at 5% of
  the carriers are another allele, an edit on the one telling base leaves the read without
  evidence.

Geometry: the three-exon gene of `rna_fixtures`, an 8-A homopolymer at [260, 268) in exon 1,
reads spliced exon 1 to exon 2 (or exon 3), sense under the default dUTP protocol. Every run
goes through `run_rna`, which asserts the four counting invariants on every row.
"""

import glob
import random

import pytest
from helpers import make_read, read_maf_output
from rna_fixtures import (
    E1,
    E2,
    E3,
    PLANTS,
    READ_LEN,
    SENSE,
    run_rna,
    runner,
    write_bam,
    write_fasta,
    write_gtf,
    write_vcf,
)

from gbcms.cli import app

ANTISENSE = 0  # forward single-end: the opposite transcript strand under dUTP ('reverse')
RUN, RUN_LEN = 260, 8  # an 8-A homopolymer at [260, 268) inside E1 [100, 300)
DONOR = E1[1]


def _ref():
    rng = random.Random(23)
    ref = [rng.choice("CGT") for _ in range(1200)]  # no stray A runs
    for pos, motif in PLANTS:
        ref[pos : pos + len(motif)] = list(motif)
    ref[RUN : RUN + RUN_LEN] = list("A" * RUN_LEN)
    return "".join(ref)


REF = _ref()
ANCHOR = RUN - 1  # 0-based base before the run


def _spliced(donor_to, kind, k, n, prefix, flag=SENSE):
    """`n` reads spliced exon 1 -> the exon starting at `donor_to`, carrying a deletion
    (`kind` "del") or an insertion ("ins") of `k` A's in the run, or REF (k = 0)."""
    out = []
    for i in range(n):
        s = 230 + (i % 5)
        if k == 0:
            left = DONOR - s
            cig = ((0, left), (3, donor_to - DONOR), (0, READ_LEN - left))
            seq = REF[s:DONOR] + REF[donor_to : donor_to + READ_LEN - left]
        elif kind == "del":
            q_left = DONOR - s - k
            right = READ_LEN - q_left
            seq = REF[s:RUN] + REF[RUN + k : DONOR] + REF[donor_to : donor_to + right]
            cig = ((0, RUN - s), (2, k), (0, DONOR - RUN - k), (3, donor_to - DONOR), (0, right))
        else:
            q_left = DONOR - s + k
            right = READ_LEN - q_left
            seq = REF[s:RUN] + "A" * k + REF[RUN:DONOR] + REF[donor_to : donor_to + right]
            cig = ((0, RUN - s), (1, k), (0, DONOR - RUN), (3, donor_to - DONOR), (0, right))
        out.append(make_read(f"{prefix}{i}", seq, s, cig, flag=flag))
    return out


# Rows: the 1-A deletion and the 1-A insertion in the run (VCF, anchor kept)
DEL_ROW = (ANCHOR + 1, REF[ANCHOR] + "A", REF[ANCHOR])
INS_ROW = (ANCHOR + 1, REF[ANCHOR], REF[ANCHOR] + "A")


def _gtf_two_transcripts(tmp_path):
    """T1 holds exons 1, 2, 3; T2 skips exon 2 (exons 1 and 3)."""
    gtf = tmp_path / "two.gtf"
    gtf.write_text(
        "".join(
            f'chr1\tTEST\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "G1"; transcript_id "{t}";\n'
            for t, exons in (("T1", (E1, E2, E3)), ("T2", (E1, E3)))
            for s, e in exons
        )
    )
    return gtf


def _run(tmp_path, row, reads, gtf=None, extra=(), name="out"):
    (res,) = run_rna(
        tmp_path,
        write_vcf(tmp_path, [row]),
        write_bam(tmp_path, REF, reads, name=f"{name}.bam"),
        write_fasta(tmp_path, REF),
        gtf if gtf is not None else write_gtf(tmp_path),
        outname=name,
        extra=extra,
    )
    return res


def _main(r):
    return (
        int(r["ref_count"]),
        int(r["alt_count"]),
        int(r["partial_alt"]),
        int(r["total_count"]),
    )


def _flags(r, col):
    return [f for f in r[col].split(";") if f]


# ═══════════════ another allele in the transcript views and ASJD ═══════════════


@pytest.mark.parametrize("kind", ["del", "ins"])
def test_another_length_reads_are_transcript_depth_only(tmp_path, kind):
    """20 REF, 4 ALT (1 A) and 12 reads of another length (2 A), all spliced exon 1 to
    exon 2. Main: RD 20, AD 4, partial 12, DP 36. The one transcript holds every read:
    AD 4, RD 20, DP = RD + AD + partial = 36, in the read and the fragment views. ASJD
    partitions REF and ALT only: 20 and 4 fragments, ALT below the junction floor."""
    row = DEL_ROW if kind == "del" else INS_ROW
    reads = (
        _spliced(E2[0], kind, 0, 20, "ref")
        + _spliced(E2[0], kind, 1, 4, "alt")
        + _spliced(E2[0], kind, 2, 12, "two")
    )
    r = _run(tmp_path, row, reads)
    assert _main(r) == (20, 4, 12, 36)
    assert "PARTIAL_DOMINANT" in _flags(r, "gbcms_diagnostic")
    assert r["transcript_read_counts"] == "T1:4,20,36"
    assert r["transcript_fragment_counts"] == "T1:4,20,36"
    assert 36 == 20 + 4 + 12
    assert (int(r["asjd_n_ref_total"]), int(r["asjd_n_alt_total"])) == (20, 4)
    assert "LOW_ALT_JUNC" in _flags(r, "asjd_diagnostic")


def test_an_exon_skipping_transcript_holds_only_the_reads_on_its_junctions(tmp_path):
    """T1 (exons 1, 2, 3) is compatible only with the reads spliced exon 1 to exon 2, T2
    (exons 1, 3) only with the reads spliced exon 1 to exon 3. Each transcript's DP holds
    its reads, other-length reads included, and no other. The main counts and ASJD hold
    all of them. (Junction-free reads, compatible with both, are in the next test.)"""
    reads = (
        _spliced(E2[0], "del", 0, 20, "ar")
        + _spliced(E2[0], "del", 1, 4, "aa")
        + _spliced(E2[0], "del", 2, 12, "ao")
        + _spliced(E3[0], "del", 0, 6, "br")
        + _spliced(E3[0], "del", 1, 3, "ba")
        + _spliced(E3[0], "del", 2, 2, "bo")
    )
    r = _run(tmp_path, DEL_ROW, reads, gtf=_gtf_two_transcripts(tmp_path))
    assert _main(r) == (26, 7, 14, 47)
    assert r["transcript_read_counts"] == "T1:4,20,36|T2:3,6,11"
    assert r["transcript_fragment_counts"] == "T1:4,20,36|T2:3,6,11"
    assert (int(r["asjd_n_ref_total"]), int(r["asjd_n_alt_total"])) == (26, 7)


def test_a_junction_free_read_fits_every_transcript_and_counts_once_in_the_row(tmp_path):
    """RJ-32: a read with no junction fits every transcript it overlaps. Five unspliced REF
    reads are added to the two-transcript row: each transcript's DP and RD hold them
    (T1 4,20,36 -> 4,25,41; T2 3,6,11 -> 3,11,16), the row's counts hold each once."""
    reads = (
        _spliced(E2[0], "del", 0, 20, "ar")
        + _spliced(E2[0], "del", 1, 4, "aa")
        + _spliced(E2[0], "del", 2, 12, "ao")
        + _spliced(E3[0], "del", 0, 6, "br")
        + _spliced(E3[0], "del", 1, 3, "ba")
        + _spliced(E3[0], "del", 2, 2, "bo")
    )
    free = [
        make_read(f"free{i}", REF[230 + i : 330 + i], 230 + i, ((0, READ_LEN),), flag=SENSE)
        for i in range(5)
    ]
    r = _run(tmp_path, DEL_ROW, reads + free, gtf=_gtf_two_transcripts(tmp_path))
    assert _main(r) == (31, 7, 14, 52)
    assert r["transcript_read_counts"] == "T1:4,25,41|T2:3,11,16"
    assert r["transcript_fragment_counts"] == "T1:4,25,41|T2:3,11,16"


def test_mnp_other_allele_reads_are_transcript_depth_only(tmp_path):
    """MNP twin: 12 reads matching the ALT at one of its two bases (partial evidence) are in
    the transcript's DP, in no transcript AD or RD, and in no ASJD partition."""
    pos = 250  # 0-based, inside exon 1, in every read's first block
    ref2 = REF[pos : pos + 2]
    alt = "".join(next(b for b in "ACGT" if b != c) for c in ref2)
    other = alt[0] + ref2[1]  # matches the ALT at one of its two bases: partial evidence

    def reads(base2, n, prefix):
        out = []
        for i in range(n):
            s = 230 + (i % 5)
            left = DONOR - s
            seq = REF[s:pos] + base2 + REF[pos + 2 : DONOR] + REF[E2[0] : E2[0] + READ_LEN - left]
            cig = ((0, left), (3, E2[0] - DONOR), (0, READ_LEN - left))
            out.append(make_read(f"{prefix}{i}", seq, s, cig, flag=SENSE))
        return out

    allr = reads(ref2, 20, "r") + reads(alt, 4, "a") + reads(other, 12, "o")
    r = _run(tmp_path, (pos + 1, ref2, alt), allr)
    rd, ad, partial, total = _main(r)
    assert (rd, ad, partial, total) == (20, 4, 12, 36)
    assert r["transcript_read_counts"] == "T1:4,20,36"
    assert r["transcript_fragment_counts"] == "T1:4,20,36"
    assert (int(r["asjd_n_ref_total"]), int(r["asjd_n_alt_total"])) == (20, 4)


# ═══════════════ antisense reads carrying another allele ═══════════════


def _antisense_reads():
    return (
        _spliced(E2[0], "del", 0, 20, "ref")
        + _spliced(E2[0], "del", 1, 6, "alt")
        + _spliced(E2[0], "del", 2, 10, "anti", flag=ANTISENSE)
    )


def test_antisense_other_length_reads_appear_in_no_column_when_enforced(tmp_path):
    """20 REF and 6 ALT sense reads plus 10 antisense reads of another length (2 A).
    Enforced: the antisense reads are in no count, no depth column, no flag; the row is
    the row of the sense reads alone."""
    r = _run(tmp_path, DEL_ROW, _antisense_reads(), name="enforced")
    assert _main(r) == (20, 6, 0, 26)
    assert int(r["rna_antisense_depth"]) == 0
    assert "PARTIAL_DOMINANT" not in _flags(r, "gbcms_diagnostic")
    (tmp_path / "sense").mkdir()
    alone = _run(tmp_path / "sense", DEL_ROW, _antisense_reads()[:26], name="sense")
    changed = [c for c in r if r[c] != alone[c]]
    assert changed == [], {c: (alone[c], r[c]) for c in changed}


def test_antisense_other_length_reads_count_as_partial_without_strandedness(tmp_path):
    """The same reads with `--no-strandedness`: 10 reads in `partial_alt`, and the
    flag `PARTIAL_DOMINANT` depends on the setting."""
    r = _run(tmp_path, DEL_ROW, _antisense_reads(), extra=("--no-strandedness",), name="free")
    assert _main(r) == (20, 6, 10, 36)
    assert "PARTIAL_DOMINANT" in _flags(r, "gbcms_diagnostic")
    assert int(r["rna_antisense_depth"]) == 0


def test_with_no_gtf_there_is_no_gene_strand_and_enforcement_is_a_no_op(tmp_path):
    """No GTF: no gene strand, so every read is sense and the antisense reads count as the
    sense ones do (the counts are the same as `--no-strandedness` with a GTF)."""
    outdir = tmp_path / "nogtf"
    outdir.mkdir()
    res = runner.invoke(
        app,
        [
            "rna",
            "-v",
            str(write_vcf(tmp_path, [DEL_ROW])),
            "-b",
            f"S:{write_bam(tmp_path, REF, _antisense_reads())}",
            "-f",
            str(write_fasta(tmp_path, REF)),
            "-o",
            str(outdir),
            "--format",
            "maf",
        ],
    )
    assert res.exit_code == 0, res.output
    (r,) = list(read_maf_output(glob.glob(str(outdir / "*.maf"))[0]))
    assert int(r["total_count"]) >= int(r["ref_count"]) + int(r["alt_count"])
    assert int(r["total_count_fragment"]) >= int(r["ref_count_fragment"]) + int(
        r["alt_count_fragment"]
    )
    assert int(r["ref_count"]) == int(r["ref_count_forward"]) + int(r["ref_count_reverse"])
    assert int(r["alt_count"]) == int(r["alt_count_forward"]) + int(r["alt_count_reverse"])
    assert _main(r) == (20, 6, 10, 36)
    assert "PARTIAL_DOMINANT" in _flags(r, "gbcms_diagnostic")


# ═══════════════ an edit inside an insertion's inserted bases ═══════════════


def _unique_anchor(ins):
    """An exon 1 anchor near 200 where `ins` cannot shift (neither end repeats a flank base)."""
    for pos in range(190, 240):
        if REF[pos] != ins[-1] and REF[pos + 1] != ins[0]:
            return pos
    raise AssertionError("no anchor")


def _insert_reads(pos, alt_ins, edited, k):
    """30 REF reads, 10 carriers of `alt_ins` and `k` carriers of `edited`, unspliced, sense."""
    reads = []
    for i in range(30):
        s = pos - 50 + (i % 5)
        reads.append(make_read(f"ref{i}", REF[s : s + READ_LEN], s, ((0, READ_LEN),), flag=SENSE))
    for tag, ins, n in (("alt", alt_ins, 10), ("edit", edited, k)):
        for i in range(n):
            s = pos - 50 + (i % 5)
            left = pos + 1 - s
            right = READ_LEN - left - len(ins)
            seq = REF[s : pos + 1] + ins + REF[pos + 1 : pos + 1 + right]
            reads.append(
                make_read(f"{tag}{i}", seq, s, ((0, left), (1, len(ins)), (0, right)), flag=SENSE)
            )
    return reads


def _editing_db(tmp_path, pos0):
    """A REDIportal TABLE1 listing the anchor (1-based position, A>G)."""
    db = tmp_path / "red.txt"
    db.write_text(
        "Accession\tRegion\tPosition\tRef\tEd\tStrand\n" f"ED1\tchr1\t{pos0 + 1}\tA\tG\t+\n"
    )
    return db


# (inserted bases, their edited form, anchor, edited carriers) -> RD, AD, partial
_INSERT_CASES = [
    ("TAC", "TGC", None, 1, (30, 11, 0)),
    ("TAC", "TGC", None, 2, (30, 10, 2)),
    ("TAC", "TGC", None, 6, (30, 10, 6)),
    ("CAG", "CGG", 200, 1, (30, 10, 1)),
]


@pytest.mark.parametrize("with_db", [False, True], ids=["no-db", "db-lists-anchor"])
@pytest.mark.parametrize("ins,edited,anchor,k,expected", _INSERT_CASES)
def test_an_edit_inside_inserted_bases_is_decided_by_the_reads(
    tmp_path, ins, edited, anchor, k, expected, with_db
):
    """One edited carrier stays ALT (AD 11); two or more at 5% are another allele (AD 10,
    partial 2 or 6); an edit on the one telling base leaves the carrier without evidence
    (AD 10, partial 1). An editing database listing the anchor changes nothing and the
    row's editing flag stays absent: the catalogue flags SNV rows only."""
    pos = anchor if anchor is not None else _unique_anchor(ins)
    reads = _insert_reads(pos, ins, edited, k)
    extra = ("--rna-editing-db", str(_editing_db(tmp_path, pos))) if with_db else ()
    r = _run(tmp_path, (pos + 1, REF[pos], REF[pos] + ins), reads, extra=extra)
    rd, ad, partial, _ = _main(r)
    assert (rd, ad, partial) == expected
    assert r["rna_editing_site"] == "False"
