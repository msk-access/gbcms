"""One reading of a variant file's extension, shared by every command.

``.vcf``, ``.vcf.gz`` and ``.vcf.bgz`` are VCFs, ``.maf`` is a MAF, case aside;
anything else is refused. The CLI accepts ``.vcf.bgz``, so ``gbcms normalize`` (which
had its own check) must read it too.
"""

import shutil
from pathlib import Path

import pysam
import pytest

from gbcms.io.input import variant_format
from gbcms.normalize import normalize_variants


@pytest.mark.parametrize(
    "name,fmt",
    [
        ("v.vcf", "vcf"),
        ("v.vcf.gz", "vcf"),
        ("v.vcf.bgz", "vcf"),
        ("V.VCF.BGZ", "vcf"),
        ("s.1.maf", "maf"),
        ("S.MAF", "maf"),
    ],
)
def test_variant_format_by_extension(name, fmt):
    assert variant_format(Path(name)) == fmt


@pytest.mark.parametrize("name", ["v.txt", "v.gz", "v.maf.gz", "v.bgz", "vcf"])
def test_variant_format_refuses_other_extensions(name):
    with pytest.raises(ValueError, match="Unsupported variant file"):
        variant_format(Path(name))


@pytest.mark.parametrize("ext", [".vcf", ".vcf.gz", ".vcf.bgz"])
def test_normalize_reads_every_vcf_extension(tmp_path, ext):
    fasta = tmp_path / "ref.fa"
    fasta.write_text(">chr1\n" + "ACGT" * 50 + "\n")
    pysam.faidx(str(fasta))
    plain = tmp_path / "plain.vcf"
    plain.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=200>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t2\t.\tC\tT\t.\t.\t.\n"
    )
    variants = tmp_path / f"variants{ext}"
    if ext == ".vcf":
        shutil.copy(plain, variants)
    else:
        pysam.tabix_compress(str(plain), str(variants))

    out = tmp_path / "normalized.tsv"
    normalize_variants(variant_file=variants, reference=fasta, output=out)

    lines = out.read_text().splitlines()
    assert len(lines) == 2, lines  # header + the one variant
    row = dict(zip(lines[0].split("\t"), lines[1].split("\t"), strict=True))
    assert (row["original_pos"], row["original_ref"], row["original_alt"]) == ("2", "C", "T")
    assert row["gbcms_status"] == "PASS"
