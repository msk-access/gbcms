"""Target contract for splice-aware indel evidence in RNA mode.

The rule these tests pin: a read testifies about an indel locus only through
ALIGNED BASES at the discriminating junction. A RefSkip (CIGAR ``N``) over
the span is **no coverage** — samtools-pileup semantics, which the engine's
own DP comment already claims — so such reads classify as neither and are
excluded from DP. A ``D`` op remains deletion evidence (the aligner asserts
deletion); an ``N`` op never is (the aligner asserts splicing). An indel
candidate directly after an ``N`` gets the same windowed inspection as one
after an ``M``.

Each test pins one face of the rule (they were committed red as
xfail(strict) before the fix landed):

  - a read spliced AROUND a deletion (span inside its N gap) carries no
    information: neither, excluded from DP
  - a variant whose anchor sits inside the intron gets depth only from
    reads with aligned bases there
  - the same junction gap counts ALT as D(100) and nothing as N(100) —
    the call must not flip on the aligner's representation choice
  - an indel op directly after a splice N gets the same inspection as one
    after an M block (M-N-D-M carriers count ALT)
  - splicing elsewhere in a read does NOT affect classification at an
    exonic locus the read covers with aligned bases

All through the public ``rna`` CLI with synthetic spliced reads.
"""

import glob
import random

import pysam
import pytest
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

BAM_CONTIG = "1"
READ_LEN = 100


def _mk_ref(n=1200, seed=11, plants=()):
    rng = random.Random(seed)
    ref = [rng.choice("ACGT") for _ in range(n)]
    for pos, motif in plants:
        ref[pos : pos + len(motif)] = list(motif)
    return "".join(ref)


def _bam(tmp_path, ref, reads, name="rna.bam"):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": BAM_CONTIG, "LN": len(ref)}]}
    p = tmp_path / name
    with pysam.AlignmentFile(p, "wb", header=header) as fh:
        for a in sorted(reads, key=lambda a: a.reference_start):
            fh.write(a)
    sp = tmp_path / name.replace(".bam", ".s.bam")
    pysam.sort("-o", str(sp), str(p))
    pysam.index(str(sp))
    return sp


def _fasta(tmp_path, ref):
    fa = tmp_path / "ref.fasta"
    fa.write_text(">chr1\n" + ref + "\n")
    pysam.faidx(str(fa))
    return fa


def _vcf(tmp_path, rows):
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(
            f"chr1\t{pos1}\t.\t{ref_al}\t{alt_al}\t.\t.\t.\n" for pos1, ref_al, alt_al in rows
        )
    )
    return vcf


def _run(tmp_path, vcf, bam, fasta, mode="rna", extra=(), outname="out"):
    outdir = tmp_path / outname
    outdir.mkdir(exist_ok=True)
    result = runner.invoke(
        app,
        [
            mode,
            "-v",
            str(vcf),
            "-b",
            f"S:{bam}",
            "-f",
            str(fasta),
            "-o",
            str(outdir),
            "--format",
            "maf",
            *extra,
        ],
    )
    assert result.exit_code == 0, result.output
    rows = list(read_maf_output(glob.glob(str(outdir / "*.maf"))[0]))
    for r in rows:
        assert int(r["total_count"]) >= int(r["ref_count"]) + int(r["alt_count"])
        assert int(r["total_count_fragment"]) >= int(r["ref_count_fragment"]) + int(
            r["alt_count_fragment"]
        )
        assert int(r["ref_count"]) == int(r["ref_count_forward"]) + int(r["ref_count_reverse"])
        assert int(r["alt_count"]) == int(r["alt_count_forward"]) + int(r["alt_count_reverse"])
        assert int(r["any_alt"]) == int(r["alt_count"]) + int(r["partial_alt"])
    return rows


def _exonic_reads(ref, center, n, prefix="ref"):
    """Plain M-only reads covering `center`."""
    out = []
    for i in range(n):
        s = center - 50 + (i % 5)
        out.append(make_read(f"{prefix}{i}", ref[s : s + READ_LEN], s, ((0, READ_LEN),)))
    return out


def _del_reads(ref, p0, length, n, prefix="alt"):
    """Reads carrying an exact deletion of `length` starting at 0-based p0."""
    out = []
    for i in range(n):
        s = p0 - 30 - (i % 5)
        left = p0 - s
        seq = ref[s:p0] + ref[p0 + length : p0 + length + (READ_LEN - left)]
        out.append(
            make_read(f"{prefix}{i}", seq, s, ((0, left), (2, length), (0, READ_LEN - left)))
        )
    return out


def _spliced_reads(ref, junction, gap, n, prefix="spl"):
    """Reads whose exon1 M block ends at `junction` (exclusive), then N(gap),
    then exon2 M — no aligned bases over [junction, junction+gap)."""
    out = []
    for i in range(n):
        s = junction - 40 - (i % 5)
        left = junction - s
        right = READ_LEN - left
        e2 = junction + gap
        seq = ref[s:junction] + ref[e2 : e2 + right]
        out.append(make_read(f"{prefix}{i}", seq, s, ((0, left), (3, gap), (0, right))))
    return out


# ── anchor/deletion geometry used by most tests ──────────────────────────
# unique context: anchor at 0-based 300, annotated 2bp deletion of 301-302
ANCHOR = 300
DEL_LEN = 2


def _del_variant_rows(ref):
    return [(ANCHOR + 1, ref[ANCHOR : ANCHOR + 1 + DEL_LEN], ref[ANCHOR])]


def test_spliced_around_deletion_carries_no_information(tmp_path):
    """Reads whose N gap covers the whole deleted span have zero aligned
    bases over it: they cannot distinguish REF from ALT. Pre-fix they count
    as definitive REF (the strict lookahead sees N, not D, and falls through
    to ref-coverage). Fixed: neither, excluded from DP."""
    ref = _mk_ref()
    reads = (
        _exonic_reads(ref, ANCHOR, 6)
        + _del_reads(ref, ANCHOR + 1, DEL_LEN, 4)
        # exon1 ends right after the anchor; N covers the deleted span
        + _spliced_reads(ref, ANCHOR + 1, 120, 5)
    )
    rows = _run(
        tmp_path,
        _vcf(tmp_path, _del_variant_rows(ref)),
        _bam(tmp_path, ref, reads),
        _fasta(tmp_path, ref),
    )
    r = rows[0]
    assert (int(r["ref_count"]), int(r["alt_count"])) == (6, 4)
    assert (
        int(r["total_count"]) == 10
    ), f"spliced-around reads must not count DP, got dp={r['total_count']}"


def test_anchor_inside_intron_carries_no_information(tmp_path):
    """A variant whose anchor lies INSIDE the N gap: spliced reads have no
    base at the anchor and no bases over the span. Pre-fix the N-inclusive
    span gate sends them to Phase 3 and counts DP. Fixed: excluded entirely;
    only reads with aligned bases there count."""
    ref = _mk_ref()
    # intron [280, 400); anchor 320 is intronic for the spliced population
    intron_start, gap = 280, 120
    anchor = 320
    rows_v = [(anchor + 1, ref[anchor : anchor + 3], ref[anchor])]
    reads = (
        _spliced_reads(ref, intron_start, gap, 5)
        # 3 reads with genuine aligned bases across the anchor (pre-mRNA-like)
        + _exonic_reads(ref, anchor, 3, prefix="pre")
    )
    rows = _run(tmp_path, _vcf(tmp_path, rows_v), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    r = rows[0]
    assert (
        int(r["total_count"]) == 3
    ), f"only base-covering reads may count DP, got dp={r['total_count']}"
    assert int(r["ref_count"]) == 3


def test_gap_representation_does_not_flip_the_call(tmp_path):
    """The same 100bp gap at the anchor: as D(100) it is deletion evidence
    (in-band → ALT); as N(100) it is asserted splicing — NO evidence either
    way. Pre-fix the N form counts definitive REF, so the REF/ALT ledger
    flips on the aligner's arbitrary D-vs-N representation choice. Fixed:
    D reads → ALT, N reads → excluded."""
    ref = _mk_ref()
    anchor, glen = 300, 100
    # Pin the anchor: if ref[anchor] == ref[anchor+glen], prep left-aligns the
    # deletion away from the junction the reads carry, and the test would
    # exercise the shifted-representation path (S3 reject → Phase 3) instead
    # of its stated intent — the D-vs-N flip at a FIXED junction. (The
    # shifted+spliced combination is a real, separate defect: D6 consensus
    # splicing corrupts ref_context coordinates and Phase 3 then drops these
    # carriers — issue #94 cluster B; it gets its own test with the D6 fix.)
    if ref[anchor] == ref[anchor + glen]:
        swap = "A" if ref[anchor] != "A" else "G"
        ref = _mk_ref(plants=((anchor + glen, swap),))
    assert ref[anchor] != ref[anchor + glen]
    rows_v = [(anchor + 1, ref[anchor : anchor + 1 + glen], ref[anchor])]
    reads = (
        _exonic_reads(ref, anchor, 6)
        + _del_reads(ref, anchor + 1, glen, 4)  # D(100) form
        + _spliced_reads(ref, anchor + 1, glen, 5)  # N(100) form, same gap
    )
    rows = _run(tmp_path, _vcf(tmp_path, rows_v), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    r = rows[0]
    assert int(r["alt_count"]) == 4
    assert (
        int(r["ref_count"]) == 6
    ), f"N-represented gaps must not count REF, got rd={r['ref_count']}"
    assert int(r["total_count"]) == 10
    # The exclusion is not silent: 5 N-represented reads exceed ad=4 at a
    # deletion locus, so the diagnostic warns that carriers may exist as
    # junction reads (STAR writes deletions >= alignIntronMin as N).
    assert "SPLICE_SKIP_DOMINANT(5)" in r["gbcms_diagnostic"]


def test_indel_across_junction_is_examined(tmp_path):
    """M-N-D-M: the deletion op sits directly after the splice N, at the
    genomic position the variant expects (anchor = last intronic base).
    Pre-fix the scans inspect only ops following an M block, so the D is
    structurally invisible. Fixed: a candidate after an N gets the same
    windowed inspection, and these carriers count ALT."""
    ref = _mk_ref()
    intron_start, gap = 280, 40
    anchor = intron_start + gap - 1  # last intronic base, 0-based 319
    p0 = anchor + 1  # deleted bases start at the first exon2 base
    rows_v = [(anchor + 1, ref[anchor : anchor + 1 + DEL_LEN], ref[anchor])]
    carriers = []
    for i in range(4):
        s = intron_start - 40 - (i % 4)
        left = intron_start - s
        right = READ_LEN - left
        seq = ref[s:intron_start] + ref[p0 + DEL_LEN : p0 + DEL_LEN + right]
        carriers.append(make_read(f"c{i}", seq, s, ((0, left), (3, gap), (2, DEL_LEN), (0, right))))
    # REF-side spliced reads: same junction, no deletion in exon2
    refs = _spliced_reads(ref, intron_start, gap, 5, prefix="r")
    rows = _run(
        tmp_path,
        _vcf(tmp_path, rows_v),
        _bam(tmp_path, ref, carriers + refs),
        _fasta(tmp_path, ref),
    )
    r = rows[0]
    assert (
        int(r["alt_count"]) == 4
    ), f"junction-adjacent D carriers must count ALT, got ad={r['alt_count']}"
    # REF-side spliced reads observe the deleted span (aligned exon2 bases)
    # but not the anchor: neither, still in DP. Nothing counts REF.
    assert (int(r["ref_count"]), int(r["total_count"])) == (0, 9)
    # Fragment layer must agree: the carriers' structural evidence survives
    # consensus even though BAQ zeroes the first exon base after the N+D
    # (ad>0 with adf=0 was a live divergence before FragmentEvidence::resolve
    # recognized structural ALT independent of base quality).
    assert (int(r["alt_count_fragment"]), int(r["total_count_fragment"])) == (4, 9)


def test_splicing_elsewhere_does_not_affect_exonic_calls(tmp_path):
    """Guard: reads spliced UPSTREAM of an exonic deletion, whose aligned
    exon2 bases fully cover the locus, classify exactly like unspliced
    reads — REF with M across the junction, ALT with the exact D."""
    ref = _mk_ref()
    intron_start, gap = 180, 60  # intron ends at 240, locus at 300 is exonic
    reads = []
    for i in range(5):  # spliced REF reads covering the locus with M
        s = intron_start - 20 - (i % 4)
        left = intron_start - s
        right = READ_LEN - left
        seq = ref[s:intron_start] + ref[intron_start + gap : intron_start + gap + right]
        reads.append(make_read(f"sr{i}", seq, s, ((0, left), (3, gap), (0, right))))
    for i in range(3):  # spliced ALT reads: M to anchor, D(2), M
        s = intron_start - 10 - (i % 3)
        left1 = intron_start - s
        e2 = intron_start + gap
        left2 = ANCHOR + 1 - e2
        right = READ_LEN - left1 - left2
        seq = (
            ref[s:intron_start]
            + ref[e2 : ANCHOR + 1]
            + ref[ANCHOR + 1 + DEL_LEN : ANCHOR + 1 + DEL_LEN + right]
        )
        reads.append(
            make_read(
                f"sa{i}", seq, s, ((0, left1), (3, gap), (0, left2), (2, DEL_LEN), (0, right))
            )
        )
    rows = _run(
        tmp_path,
        _vcf(tmp_path, _del_variant_rows(ref)),
        _bam(tmp_path, ref, reads),
        _fasta(tmp_path, ref),
    )
    r = rows[0]
    assert (int(r["ref_count"]), int(r["alt_count"])) == (5, 3)
    assert int(r["total_count"]) == 8


def test_intronic_snv_counts_no_spliced_depth(tmp_path):
    """An SNV inside the intron: spliced reads have no base there. They
    already classify neither (find_read_pos is N-aware) but pre-fix they
    still count DP through the N-inclusive overlap gate. Fixed: DP counts
    only reads with an aligned base at the position."""
    ref = _mk_ref()
    intron_start, gap = 280, 120
    snv_pos = 340  # intronic for the spliced population
    ref_base = ref[snv_pos]
    alt_base = "A" if ref_base != "A" else "G"
    rows_v = [(snv_pos + 1, ref_base, alt_base)]
    reads = _spliced_reads(ref, intron_start, gap, 6)
    rows = _run(tmp_path, _vcf(tmp_path, rows_v), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    r = rows[0]
    assert (
        int(r["total_count"]) == 0
    ), f"pure-spliced pileup at an intronic SNV must have DP 0, got {r['total_count']}"


def test_insertion_after_junction_is_examined(tmp_path):
    """M-N-I-M: an insertion at the exon boundary reached through the splice.
    The I op's only left neighbor is the N, so pre-fix no scan ever saw it.
    Fixed: the post-N inspection resolves it at the exact expected junction,
    attributing the evidence to the first inserted base."""
    ref = _mk_ref()
    intron_start, gap = 280, 40
    anchor = intron_start + gap - 1  # last intronic base, 0-based 319
    # Pin against left-alignment: the insert must not end with the anchor base.
    last = "A" if ref[anchor] != "A" else "C"
    ins_seq = "TG" + last
    rows_v = [(anchor + 1, ref[anchor], ref[anchor] + ins_seq)]
    e2 = intron_start + gap  # first exon2 base, 320
    carriers = []
    for i in range(4):
        s = intron_start - 40 - (i % 4)
        left = intron_start - s
        right = READ_LEN - left - len(ins_seq)
        seq = ref[s:intron_start] + ins_seq + ref[e2 : e2 + right]
        carriers.append(
            make_read(f"ic{i}", seq, s, ((0, left), (3, gap), (1, len(ins_seq)), (0, right)))
        )
    refs = _spliced_reads(ref, intron_start, gap, 5, prefix="ir")
    vcf = _vcf(tmp_path, rows_v)
    bam = _bam(tmp_path, ref, carriers + refs)
    fasta = _fasta(tmp_path, ref)
    # Default RNA BAQ stacks splice + indel penalties on the inserted bases
    # (they sit within 5bp of both the N and the I), so the sequence cannot
    # be verified at >= min_baseq: honest outcome is partial evidence, never
    # REF absorption. (With --gtf, boundary suppression restores the quals
    # at annotated junctions.)
    r = _run(tmp_path, vcf, bam, fasta)[0]
    assert (int(r["ref_count"]), int(r["alt_count"])) == (0, 0)
    assert (
        int(r["partial_alt"]) == 4
    ), f"BAQ-unverifiable junction inserts must carry partial evidence, got {r['partial_alt']}"
    # Without BAQ the inserted bases keep their real quality and the post-N
    # inspection confirms the exact-sequence match: full ALT, read and
    # fragment level.
    r = _run(tmp_path, vcf, bam, fasta, extra=("--no-baq",), outname="out_nobaq")[0]
    assert (
        int(r["alt_count"]) == 4
    ), f"junction-adjacent I carriers must count ALT, got ad={r['alt_count']}"
    assert int(r["alt_count_fragment"]) == 4


def test_spliced_over_insertion_carries_no_information(tmp_path):
    """An insertion annotated inside the spliced-out intron: reads whose N
    covers both junction flanks (and carry no indel op near the locus)
    observe nothing there — excluded from DP."""
    ref = _mk_ref()
    intron_start, gap = 280, 120
    anchor = 330  # intronic for the spliced population
    last = "A" if ref[anchor] != "A" else "C"
    rows_v = [(anchor + 1, ref[anchor], ref[anchor] + "TG" + last)]
    reads = _spliced_reads(ref, intron_start, gap, 5) + _exonic_reads(ref, anchor, 3, prefix="pre")
    rows = _run(tmp_path, _vcf(tmp_path, rows_v), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    r = rows[0]
    assert (
        int(r["total_count"]) == 3
    ), f"only flank-observing reads may count DP, got dp={r['total_count']}"
    assert int(r["ref_count"]) == 3


def test_windowed_deletion_after_junction_in_repeat(tmp_path):
    """Shifted representation across a junction: the aligner extends the N
    over the annotated span and writes the D two bases downstream inside the
    same AT tract — the same event, shifted. The triage must NOT exclude the
    read (it carries an indel op inside the scan window) and the post-N
    windowed scan verifies the deleted bases at the shifted position."""
    # ref[339]='C' pins the anchor (no left-shift), tract AT×5 at [340, 350),
    # ref[350]='G' ends the tract.
    ref = _mk_ref(plants=((339, "CATATATATATG"),))
    anchor = 339
    rows_v = [(anchor + 1, ref[anchor : anchor + 3], ref[anchor])]  # del [340,342)
    exon1_end = 300
    carriers = []
    for i in range(4):
        s = exon1_end - 40 - (i % 4)
        left = exon1_end - s
        right = READ_LEN - left
        # N covers [300, 344) — including the annotated span — then D(2) at
        # [344, 346) inside the tract, then exon2 M.
        seq = ref[s:exon1_end] + ref[346 : 346 + right]
        carriers.append(
            make_read(f"wc{i}", seq, s, ((0, left), (3, 344 - exon1_end), (2, 2), (0, right)))
        )
    refs = _exonic_reads(ref, anchor, 5)
    vcf = _vcf(tmp_path, rows_v)
    bam = _bam(tmp_path, ref, carriers + refs)
    fasta = _fasta(tmp_path, ref)
    # DNA mode isolates the cluster-A machinery (no D6 consensus splicing,
    # no BAQ): the triage defers, the post-N windowed scan verifies the
    # deleted bases at the shifted position, carriers count ALT.
    r = _run(tmp_path, vcf, bam, fasta, mode="dna")[0]
    assert (
        int(r["alt_count"]) == 4
    ), f"shifted post-N carriers must count ALT, got ad={r['alt_count']}"
    assert int(r["ref_count"]) == 5
    assert int(r["total_count"]) == 9


@pytest.mark.xfail(
    strict=True,
    reason="D6 consensus splicing corrupts ref_context offsets: the S3 "
    "sequence check reads the spliced context at genomic coordinates and "
    "rejects the shifted candidate (issue #94 cluster B)",
)
def test_windowed_deletion_after_junction_in_repeat_rna(tmp_path):
    """RNA-mode twin of the DNA-mode case above. The 4 junction carriers
    reach the post-N windowed scan (splice_skip_excluded=0 proves the triage
    deferred), but D6 splices the intron out of ref_context without a
    coordinate map, so verify_deleted_bases reads garbage at the shifted
    position and the same-event carriers fall to neither."""
    ref = _mk_ref(plants=((339, "CATATATATATG"),))
    anchor = 339
    rows_v = [(anchor + 1, ref[anchor : anchor + 3], ref[anchor])]
    exon1_end = 300
    carriers = []
    for i in range(4):
        s = exon1_end - 40 - (i % 4)
        left = exon1_end - s
        right = READ_LEN - left
        seq = ref[s:exon1_end] + ref[346 : 346 + right]
        carriers.append(
            make_read(f"wc{i}", seq, s, ((0, left), (3, 344 - exon1_end), (2, 2), (0, right)))
        )
    refs = _exonic_reads(ref, anchor, 5)
    rows = _run(
        tmp_path,
        _vcf(tmp_path, rows_v),
        _bam(tmp_path, ref, carriers + refs),
        _fasta(tmp_path, ref),
    )
    r = rows[0]
    assert (
        int(r["alt_count"]) == 4
    ), f"shifted post-N carriers must count ALT, got ad={r['alt_count']}"


def test_intronic_mnp_counts_no_spliced_depth(tmp_path):
    """A 2bp MNP inside the intron: spliced reads observe neither position —
    excluded from DP entirely."""
    ref = _mk_ref()
    intron_start, gap = 280, 120
    p = 335  # intronic
    alt = "".join("A" if b != "A" else "G" for b in ref[p : p + 2])
    rows_v = [(p + 1, ref[p : p + 2], alt)]
    reads = _spliced_reads(ref, intron_start, gap, 6)
    rows = _run(tmp_path, _vcf(tmp_path, rows_v), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    r = rows[0]
    assert (
        int(r["total_count"]) == 0
    ), f"pure-spliced pileup at an intronic MNP must have DP 0, got {r['total_count']}"


def test_spliced_over_delins_carries_no_information(tmp_path):
    """A delins (check_complex path) whose whole REF span sits inside the N
    gap: spliced reads are excluded before any reconstruction or alignment
    can stitch across the splice."""
    ref = _mk_ref()
    intron_start, gap = 280, 120
    p = 330  # intronic REF span [330, 333)
    alt = "".join("A" if b != "A" else "G" for b in ref[p : p + 2])  # 3bp -> 2bp
    rows_v = [(p + 1, ref[p : p + 3], alt)]
    reads = _spliced_reads(ref, intron_start, gap, 5) + _exonic_reads(ref, p, 3, prefix="pre")
    rows = _run(tmp_path, _vcf(tmp_path, rows_v), _bam(tmp_path, ref, reads), _fasta(tmp_path, ref))
    r = rows[0]
    assert (
        int(r["total_count"]) == 3
    ), f"spliced-over delins reads must not count DP, got dp={r['total_count']}"


def test_legacy_parity_with_spliced_reads(tmp_path):
    """The binned↔legacy parity oracle holds for N-CIGAR reads: the
    splice-skip exclusion and the post-N deletion evidence are mirrored in
    the legacy count_bam path (count_both asserts every parity field)."""
    from helpers import build_bam, count_both

    from gbcms._rs import Variant

    rng = random.Random(7)
    ref = "".join(rng.choice("ACGT") for _ in range(500))
    intron_start, gap = 200, 40
    anchor = intron_start + gap - 1  # last intronic base, 239
    p0 = anchor + 1
    rl = 100
    reads = []
    for i in range(4):  # M-N-D-M carriers
        s = intron_start - 40 - (i % 4)
        left = intron_start - s
        right = rl - left
        seq = ref[s:intron_start] + ref[p0 + 2 : p0 + 2 + right]
        reads.append(make_read(f"pc{i}", seq, s, ((0, left), (3, gap), (2, 2), (0, right))))
    for i in range(5):  # spliced-over-anchor reads, aligned over the span
        s = intron_start - 30 - (i % 5)
        left = intron_start - s
        right = rl - left
        seq = ref[s:intron_start] + ref[intron_start + gap : intron_start + gap + right]
        reads.append(make_read(f"pr{i}", seq, s, ((0, left), (3, gap), (0, right))))
    for i in range(3):  # pre-mRNA reads with aligned bases across the locus
        s = anchor - 50 + i
        reads.append(make_read(f"pm{i}", ref[s : s + rl], s, ((0, rl),)))
    bam = build_bam(tmp_path, reads)
    pad = 8
    v = Variant(
        chrom="chr1",
        pos=anchor,
        ref_allele=ref[anchor : p0 + 2],
        alt_allele=ref[anchor],
        variant_type="DELETION",
        ref_context=ref[anchor - pad : p0 + 2 + pad],
        ref_context_start=anchor - pad,
    )
    c = count_both(bam, [v], min_mapq=0, min_baseq=0)[0]
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev
    assert c.ad == 4, f"M-N-D-M carriers must count ALT in both paths, got ad={c.ad}"
    assert c.rd == 3, f"only base-covering pre-mRNA reads count REF, got rd={c.rd}"
