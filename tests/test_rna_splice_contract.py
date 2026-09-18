"""Target contract for splice-aware indel evidence in RNA mode.

The rule these tests pin: a read testifies about an indel locus only through
ALIGNED BASES at the discriminating junction. A RefSkip (CIGAR ``N``) over
the span is **no coverage** — samtools-pileup semantics, which the engine's
own DP comment already claims — so such reads classify as neither and are
excluded from DP. A ``D`` op remains deletion evidence (the aligner asserts
deletion); an ``N`` op never is (the aligner asserts splicing). An indel
candidate directly after an ``N`` gets the same windowed inspection as one
after an ``M``.

Each xfail(strict=True) test documents a verified defect and states the
fixed behavior; guards pin behavior that must not change:

  - a read spliced AROUND a deletion (span inside its N gap) currently
    counts as definitive REF and inflates DP
  - a variant whose anchor sits inside the intron currently reaches Phase 3
    through the N-inclusive span gate and counts DP with zero bases there
  - the same junction gap counts ALT as D(100) but REF as N(100) — the
    call flips on the aligner's representation choice
  - an indel op separated from the anchor by a splice N is structurally
    invisible (M-N-D-M carriers count as REF)
  - splicing elsewhere in a read must NOT affect classification at an
    exonic locus the read covers with aligned bases (guards)

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

XFAIL = pytest.mark.xfail(strict=True, reason="splice-aware evidence rule: fix pending")

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


def _run(tmp_path, vcf, bam, fasta):
    outdir = tmp_path / "out"
    outdir.mkdir(exist_ok=True)
    result = runner.invoke(
        app,
        [
            "rna",
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
        ],
    )
    assert result.exit_code == 0, result.output
    rows = list(read_maf_output(glob.glob(str(outdir / "*.maf"))[0]))
    for r in rows:
        assert int(r["total_count"]) >= int(r["ref_count"]) + int(r["alt_count"])
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


@XFAIL
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


@XFAIL
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


@XFAIL
def test_gap_representation_does_not_flip_the_call(tmp_path):
    """The same 100bp gap at the anchor: as D(100) it is deletion evidence
    (in-band → ALT); as N(100) it is asserted splicing — NO evidence either
    way. Pre-fix the N form counts definitive REF, so the REF/ALT ledger
    flips on the aligner's arbitrary D-vs-N representation choice. Fixed:
    D reads → ALT, N reads → excluded."""
    ref = _mk_ref()
    anchor, glen = 300, 100
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


@XFAIL
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


@XFAIL
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
