"""The Rust log cap during a counting call (a build-speed change).

For the length of a counting call the Rust `log` crate's maximum level follows the Python
loggers, so a `trace!` or `debug!` that Python would drop is not built on every read. The
contract: a Rust record that a Python logger would pass is never dropped, inside a counting
call or outside one.

- The CLI loads the Rust module before it sets the log level. The variant preparation logs at
  INFO before any counting call, and those lines stay in the log.
- After a counting call the maximum is open again, so a level lowered later is followed by
  code outside counting with no reset of the cap.
- Inside a counting call a child logger set more verbose than the root is followed.
"""

import logging
import os
import random
import subprocess
import sys
from pathlib import Path

from helpers import count_bam_checked, make_read, write_contig

from gbcms import _rs

DATA = Path(__file__).parent / "testdata"
ENGINE_LOGGER = "_rs.counting.engine"
NORMALIZE_LOGGER = "_rs.normalize.engine"
ARGS = (20, 20, True, True, True, False, False, False, 1)
A = 200


def _contig_and_snv(tmp_path):
    rng = random.Random(7)
    ref = "".join(rng.choice("ACGT") for _ in range(400))
    alt = "C" if ref[A] != "C" else "T"
    reads = []
    for i in range(6):
        s = A - 50 + i
        seq = ref[s : s + 100]
        if i % 2:
            seq = seq[: A - s] + alt + seq[A - s + 1 :]
        reads.append(make_read(f"r{i}", seq, s, ((0, 100),)))
    fa, bam = write_contig(tmp_path, ref, reads, "cap")
    return fa, bam, ref[A], alt


def _invariants(c):
    assert c.dp >= c.rd + c.ad
    assert c.dpf >= c.rdf + c.adf
    assert c.rd == c.rd_fwd + c.rd_rev
    assert c.ad == c.ad_fwd + c.ad_rev


def test_cli_keeps_the_rust_info_lines_logged_before_counting(tmp_path):
    cmd = [
        sys.executable,
        "-c",
        "from gbcms.cli import app; app()",
        "dna",
        "-v",
        str(DATA / "integration_test_variants.vcf"),
        "-b",
        f"S1:{DATA / 'sample1_integration_test.bam'}",
        "-f",
        str(DATA / "integration_test_reference.fa"),
        "-o",
        str(tmp_path),
    ]
    env = {**os.environ, "COLUMNS": "200"}
    r = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    log = r.stdout + r.stderr
    assert "prepare_variants: 3 variants" in log
    assert "prepare_variants complete" in log


def test_the_cap_ends_with_the_counting_call(tmp_path, caplog):
    fa, bam, ref_base, alt = _contig_and_snv(tmp_path)
    with caplog.at_level(logging.WARNING):
        _rs.reset_log_caching()
        (c,) = count_bam_checked(bam, [_rs.Variant("1", A, ref_base, alt, "SNP")], [None], *ARGS)
    _invariants(c)
    assert c.ad == 3
    caplog.clear()
    with caplog.at_level(logging.INFO, logger=NORMALIZE_LOGGER):
        _rs.reset_log_caching()
        _rs.prepare_variants([_rs.Variant("1", A, ref_base, alt, "SNP")], fa, 5, False)
    assert any(
        r.name == NORMALIZE_LOGGER and r.getMessage().startswith("prepare_variants:")
        for r in caplog.records
    )


def test_a_child_logger_more_verbose_than_the_root_is_followed(tmp_path, caplog):
    _, bam, ref_base, alt = _contig_and_snv(tmp_path)
    logging.getLogger().setLevel(logging.WARNING)
    with caplog.at_level(logging.DEBUG, logger=ENGINE_LOGGER):
        _rs.reset_log_caching()
        (c,) = count_bam_checked(bam, [_rs.Variant("1", A, ref_base, alt, "SNP")], [None], *ARGS)
    _invariants(c)
    assert any(
        r.name == ENGINE_LOGGER
        and r.levelno == logging.DEBUG
        and r.getMessage().startswith("Bin tid=")
        for r in caplog.records
    )
