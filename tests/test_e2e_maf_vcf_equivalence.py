"""MAF and VCF inputs describing the same variants must count identically.

MAF and VCF encode indels differently (MAF: ``-`` alleles, Start_Position at
the first deleted base / at the anchor; VCF: shared anchor base, POS at the
anchor). Both converge in ``prepare_variants`` (anchor resolution against the
FASTA + left-alignment), so every downstream count must be byte-identical
regardless of which representation the caller used. This test commits that
contract for the two indel classes the wrong-length work touched most: a
homopolymer deletion locus that also carries wrong-length (distinct-allele)
reads, and a unique-context insertion locus.

Harness conventions follow ``test_e2e_dna_maf.py``: FASTA and variant files
use ``chr1``, the BAM header uses the normalized name ``1``, reads are
MAPQ 60 / Q30. Everything is synthetic and lives in ``tmp_path``.
"""

import glob
import random

import pysam
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

BAM_CONTIG = "1"

# Deletion locus: A(10) tract at 200-209, anchor 'C' at 199. Annotated 2bp del.
DEL_ANCHOR = 199
DEL_LEN = 2
# Insertion locus: unique context, +TG after 0-based 400.
INS_ANCHOR = 400
INS_SEQ = "TG"

READ_LEN = 100

# The gbcms count columns whose equality IS the contract.
COUNT_COLS = [
    "total_count",
    "ref_count",
    "alt_count",
    "any_alt",
    "partial_alt",
    "n_count",
    "ref_count_forward",
    "ref_count_reverse",
    "alt_count_forward",
    "alt_count_reverse",
    "total_count_fragment",
    "ref_count_fragment",
    "alt_count_fragment",
    "gbcms_diagnostic",
]


def _reference_sequence():
    rng = random.Random(7)
    seq = [rng.choice("ACGT") for _ in range(600)]
    seq[199] = "C"
    seq[200:210] = list("A" * 10)
    seq[210] = "G"
    return "".join(seq)


def _build_reference(tmp_path, ref):
    fasta = tmp_path / "ref.fasta"
    fasta.write_text(">chr1\n" + ref + "\n")
    pysam.faidx(str(fasta))
    return fasta


def _build_reads_bam(tmp_path, ref):
    """Deletion locus: 8 REF + 4 exact D(2) + 6 wrong-length D(1).
    Insertion locus: 6 REF + 5 exact I(2)."""
    reads = []
    p0 = DEL_ANCHOR + 1  # first deleted base
    for i in range(8):
        s = 150 + (i % 5)
        reads.append(make_read(f"dref_{i}", ref[s : s + READ_LEN], s, ((0, READ_LEN),)))
    for i in range(4):
        s = p0 - 30 - (i % 3)
        left = p0 - s
        seq = ref[s:p0] + ref[p0 + DEL_LEN : p0 + DEL_LEN + (READ_LEN - left)]
        reads.append(
            make_read(f"dalt_{i}", seq, s, ((0, left), (2, DEL_LEN), (0, READ_LEN - left)))
        )
    for i in range(6):
        s = p0 - 30 - (i % 5)
        left = p0 - s
        seq = ref[s:p0] + ref[p0 + 1 : p0 + 1 + (READ_LEN - left)]
        reads.append(make_read(f"dwl_{i}", seq, s, ((0, left), (2, 1), (0, READ_LEN - left))))

    a1 = INS_ANCHOR + 1
    for i in range(6):
        s = 350 + (i % 5)
        reads.append(make_read(f"iref_{i}", ref[s : s + READ_LEN], s, ((0, READ_LEN),)))
    for i in range(5):
        s = a1 - 30 - (i % 4)
        left = a1 - s
        ins_len = len(INS_SEQ)
        seq = ref[s:a1] + INS_SEQ + ref[a1 : a1 + (READ_LEN - left - ins_len)]
        reads.append(
            make_read(
                f"ialt_{i}", seq, s, ((0, left), (1, ins_len), (0, READ_LEN - left - ins_len))
            )
        )

    unsorted = tmp_path / "sample.unsorted.bam"
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 600, "SN": BAM_CONTIG}]}
    with pysam.AlignmentFile(unsorted, "wb", header=header) as outf:
        for r in sorted(reads, key=lambda a: a.reference_start):
            outf.write(r)
    sorted_bam = tmp_path / "sample.sorted.bam"
    pysam.sort("-o", str(sorted_bam), str(unsorted))
    pysam.index(str(sorted_bam))
    return str(sorted_bam)


def _build_vcf(tmp_path, ref):
    """VCF: anchor-based. Deletion POS = anchor+1 (1-based anchor), REF = anchor+deleted.
    Insertion POS = anchor+1, ALT = anchor + inserted."""
    vcf = tmp_path / "variants.vcf"
    del_ref = ref[DEL_ANCHOR : DEL_ANCHOR + 1 + DEL_LEN]
    ins_ref = ref[INS_ANCHOR]
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=600>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"chr1\t{DEL_ANCHOR + 1}\t.\t{del_ref}\t{del_ref[0]}\t.\t.\t.\n"
        f"chr1\t{INS_ANCHOR + 1}\t.\t{ins_ref}\t{ins_ref}{INS_SEQ}\t.\t.\t.\n"
    )
    return vcf


def _build_maf(tmp_path, ref):
    """MAF: '-' alleles. Deletion Start = first deleted base (1-based),
    REF = deleted bases, ALT = '-'. Insertion Start = anchor (1-based),
    REF = '-', ALT = inserted bases."""
    maf = tmp_path / "variants.maf"
    deleted = ref[DEL_ANCHOR + 1 : DEL_ANCHOR + 1 + DEL_LEN]
    del_start = DEL_ANCHOR + 2  # 1-based first deleted base
    del_end = del_start + DEL_LEN - 1
    ins_start = INS_ANCHOR + 1  # 1-based anchor
    maf.write_text(
        "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
        "Reference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode\n"
        f"GENE1\tchr1\t{del_start}\t{del_end}\t{deleted}\t-\tSAMPLE\n"
        f"GENE2\tchr1\t{ins_start}\t{ins_start + 1}\t-\t{INS_SEQ}\tSAMPLE\n"
    )
    return maf


def _run_and_read(tmp_path, variants_file, tag, bam, fasta):
    outdir = tmp_path / f"out_{tag}"
    outdir.mkdir()
    result = runner.invoke(
        app,
        [
            "dna",
            "-v",
            str(variants_file),
            "-b",
            str(bam),
            "-f",
            str(fasta),
            "-o",
            str(outdir),
            "--format",
            "maf",
        ],
    )
    assert result.exit_code == 0, result.output
    maf_files = glob.glob(str(outdir / "*.maf"))
    assert len(maf_files) == 1
    rows = list(read_maf_output(maf_files[0]))
    assert len(rows) == 2, f"expected 2 variant rows from {tag} input, got {len(rows)}"
    return rows


def test_maf_and_vcf_inputs_count_identically(tmp_path):
    ref = _reference_sequence()
    fasta = _build_reference(tmp_path, ref)
    bam = _build_reads_bam(tmp_path, ref)

    vcf_rows = _run_and_read(tmp_path, _build_vcf(tmp_path, ref), "vcf", bam, fasta)
    maf_rows = _run_and_read(tmp_path, _build_maf(tmp_path, ref), "maf", bam, fasta)

    for i, (v, m) in enumerate(zip(vcf_rows, maf_rows, strict=True)):
        for col in COUNT_COLS:
            assert v[col] == m[col], (
                f"row {i}: {col} differs between input representations: "
                f"VCF={v[col]!r} MAF={m[col]!r}"
            )

    # Sanity: the counts themselves are the expected synthetic truth, so the
    # equivalence is not two identically-wrong results.
    del_row, ins_row = vcf_rows
    assert (int(del_row["ref_count"]), int(del_row["alt_count"]), int(del_row["partial_alt"])) == (
        8,
        4,
        6,
    )
    assert (int(ins_row["ref_count"]), int(ins_row["alt_count"]), int(ins_row["partial_alt"])) == (
        6,
        5,
        0,
    )
    # Counting invariants (house rules).
    for row in vcf_rows:
        assert int(row["total_count"]) >= int(row["ref_count"]) + int(row["alt_count"])
        assert int(row["any_alt"]) == int(row["alt_count"]) + int(row["partial_alt"])
        assert int(row["ref_count_forward"]) + int(row["ref_count_reverse"]) == int(
            row["ref_count"]
        )
        assert int(row["alt_count_forward"]) + int(row["alt_count_reverse"]) == int(
            row["alt_count"]
        )
