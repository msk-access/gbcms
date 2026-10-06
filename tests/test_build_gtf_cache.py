"""The ``build-gtf-cache`` command, deprecated in 6.6.0.

It pre-built a bincode cache of the parsed GTF (M5a). A byte-level parser now loads
the GTF in about a second, so the cache saved ~0.5 s a sample (CYCLE_6.6.0_PLAN.md,
"D4 GTF loading"). The command stays for one release so pipelines that call it keep
working: it checks its options as before, warns that it is deprecated, and writes
nothing. It goes in 6.7.0.
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


def test_build_gtf_cache_is_deprecated_and_writes_nothing(tmp_path):
    gtf = _write(tmp_path, "tiny.gtf", _GTF)
    vcf = _write(tmp_path, "v.vcf", _VCF)
    cache = tmp_path / "cache"

    args = ["build-gtf-cache", "--gtf", str(gtf), "--variants", str(vcf)]
    result = runner.invoke(app, [*args, "--gtf-cache-dir", str(cache)])
    assert result.exit_code == 0, result.output
    assert "deprecated" in result.output.lower(), result.output
    assert not cache.exists() or not any(cache.iterdir())


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
