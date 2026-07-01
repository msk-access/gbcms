"""M5a: the ``build-gtf-cache`` pre-warm command.

Verifies that the command parses a GTF for the variant chromosomes and writes a
single reusable cache file into ``--gtf-cache-dir`` (which a later run hits to skip
the parse), and that it rejects an unsupported variant-file extension.

The real cold→warm speedup and cross-run cache-hit fidelity are validated on the
full Ensembl GTF outside the unit suite; here we use a tiny inline fixture.
"""

from pathlib import Path

from typer.testing import CliRunner

from gbcms.cli import app

runner = CliRunner()

_GTF = (
    '1\ttest\texon\t100\t200\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
    '1\ttest\texon\t300\t400\t.\t+\t.\tgene_id "ENSG1"; transcript_id "ENST1";\n'
)
_VCF = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=1>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "1\t150\t.\tA\tT\t.\t.\t.\n"
)


def _write(tmp_path: Path, name: str, content: str) -> Path:
    p = tmp_path / name
    p.write_text(content)
    return p


def test_build_gtf_cache_writes_single_reusable_cache(tmp_path):
    gtf = _write(tmp_path, "tiny.gtf", _GTF)
    vcf = _write(tmp_path, "v.vcf", _VCF)
    cache = tmp_path / "cache"

    result = runner.invoke(
        app,
        [
            "build-gtf-cache",
            "--gtf",
            str(gtf),
            "--variants",
            str(vcf),
            "--gtf-cache-dir",
            str(cache),
        ],
    )
    assert result.exit_code == 0, result.output
    files = list(cache.glob("gbcms-gtf-*.idx"))
    assert len(files) == 1, f"expected exactly one cache file, found {files}"

    # Re-running with the same inputs reuses the existing entry (no proliferation).
    result2 = runner.invoke(
        app,
        [
            "build-gtf-cache",
            "--gtf",
            str(gtf),
            "--variants",
            str(vcf),
            "--gtf-cache-dir",
            str(cache),
        ],
    )
    assert result2.exit_code == 0, result2.output
    assert len(list(cache.glob("gbcms-gtf-*.idx"))) == 1


def test_build_gtf_cache_rejects_bad_variant_extension(tmp_path):
    gtf = _write(tmp_path, "tiny.gtf", _GTF)
    bad = _write(tmp_path, "variants.txt", "not a variant file")
    cache = tmp_path / "cache"
    result = runner.invoke(
        app,
        [
            "build-gtf-cache",
            "--gtf",
            str(gtf),
            "--variants",
            str(bad),
            "--gtf-cache-dir",
            str(cache),
        ],
    )
    assert result.exit_code != 0
