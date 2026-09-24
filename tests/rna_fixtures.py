"""Shared RNA fixtures: a three-exon gene model with planted canonical splice
motifs, file writers for its reference / GTF / BAM / variants, a CLI runner
that asserts the counting invariants, and read builders. Contig names are
parameters (defaults: ``chr1`` variants / FASTA / GTF, ``1`` BAM), so naming
tests can reuse the same geometry.
"""

import glob
import random

import pysam
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

BAM_CONTIG = "1"
READ_LEN = 100

# ── Gene model (0-based, half-open) ─────────────────────────────────────
# E1 [100,300)  intron1 [300,500)  E2 [500,600)  intron2 [600,800)  E3 [800,1000)
# Canonical motifs planted: GT at 300 / AG at 498; GT at 600 / AG at 798.
E1, E2, E3 = (100, 300), (500, 600), (800, 1000)
PLANTS = ((300, "GT"), (498, "AG"), (600, "GT"), (798, "AG"))


def mk_ref(n=1200, seed=23):
    rng = random.Random(seed)
    ref = [rng.choice("ACGT") for _ in range(n)]
    for pos, motif in PLANTS:
        ref[pos : pos + len(motif)] = list(motif)
    return "".join(ref)


def write_gtf(tmp_path, contig="chr1"):
    gtf = tmp_path / "gene.gtf"
    gtf.write_text(
        "".join(
            f'{contig}\tTEST\texon\t{s + 1}\t{e}\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
            for s, e in (E1, E2, E3)
        )
    )
    return gtf


def write_bam(tmp_path, ref, reads, name="rna.bam", contig=BAM_CONTIG):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": contig, "LN": len(ref)}]}
    p = tmp_path / name
    with pysam.AlignmentFile(p, "wb", header=header) as fh:
        for a in sorted(reads, key=lambda a: a.reference_start):
            fh.write(a)
    sp = tmp_path / name.replace(".bam", ".s.bam")
    pysam.sort("-o", str(sp), str(p))
    pysam.index(str(sp))
    return sp


def write_fasta(tmp_path, ref, contig="chr1"):
    fa = tmp_path / "ref.fasta"
    fa.write_text(f">{contig}\n{ref}\n")
    pysam.faidx(str(fa))
    return fa


def write_vcf(tmp_path, rows, contig="chr1"):
    vcf = tmp_path / "variants.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(f"{contig}\t{p}\t.\t{r}\t{a}\t.\t.\t.\n" for p, r, a in rows)
    )
    return vcf


def write_maf(tmp_path, rows):
    """MAF twin of write_vcf: anchor-preserved deletions become Start = first
    deleted base, REF = deleted bases, ALT = '-'; SNVs are identical."""
    maf = tmp_path / "variants.maf"
    lines = [
        "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
        "Reference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode"
    ]
    for p, r, a in rows:
        if len(a) == 1 and len(r) > 1 and r[0] == a:
            d = r[1:]
            lines.append(f"G1\tchr1\t{p + 1}\t{p + len(d)}\t{d}\t-\tS")
        else:
            lines.append(f"G1\tchr1\t{p}\t{p + len(r) - 1}\t{r}\t{a}\tS")
    maf.write_text("\n".join(lines) + "\n")
    return maf


def run_rna(tmp_path, variants, bam, fasta, gtf, outname="out"):
    outdir = tmp_path / outname
    outdir.mkdir(exist_ok=True)
    result = runner.invoke(
        app,
        [
            "rna",
            "-v",
            str(variants),
            "-b",
            f"S:{bam}",
            "-f",
            str(fasta),
            "-o",
            str(outdir),
            "--format",
            "maf",
            "--gtf",
            str(gtf),
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


# ── Read builders ────────────────────────────────────────────────────────
# Reads are reverse-strand single-end (flag 16): sense for a '+' gene under
# the default dUTP ('reverse') protocol, so the RNA defaults
# (--enforce-strandedness, BAQ with GTF boundary suppression) apply as in
# production rather than filtering the synthetic reads as antisense.
SENSE = 16


def spliced(ref, donor, acceptor, n, prefix):
    """Reads whose left M block ends at `donor` (exclusive), then N over
    [donor, acceptor), then M from `acceptor`."""
    out = []
    for i in range(n):
        s = donor - 40 - (i % 5)
        left = donor - s
        right = READ_LEN - left
        seq = ref[s:donor] + ref[acceptor : acceptor + right]
        out.append(
            make_read(
                f"{prefix}{i}", seq, s, ((0, left), (3, acceptor - donor), (0, right)), flag=SENSE
            )
        )
    return out


def through(ref, pos, base, n, prefix):
    """Unspliced (M-only) reads covering `pos`, carrying `base` there."""
    out = []
    for i in range(n):
        s = pos - 50 + (i % 5)
        seq = ref[s:pos] + base + ref[pos + 1 : s + READ_LEN]
        out.append(make_read(f"{prefix}{i}", seq, s, ((0, READ_LEN),), flag=SENSE))
    return out
