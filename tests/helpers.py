"""
Shared test helpers for gbcms.

Provides:
- BAM construction helpers (build_bam, make_read, write_contig) for synthetic data
- count_bam_checked / count_checked / count_one_checked: production counting,
  checked for binning invariance (every field equal with one variant per bin)
- assert_same_counts: every field of two lists of BaseCounts equal
- install_suite_read_judge: the suite under the haplotype judge
  (GBCMS_TEST_READ_JUDGE=haplotype; installed by conftest.py)
"""

import functools
import glob
import math
import os

import pysam

from gbcms import _rs as gbcms_rs

# The core read and fragment counts, for tests that compare two runs on them.
COUNT_FIELDS = [
    "dp",
    "rd",
    "ad",
    "dp_fwd",
    "rd_fwd",
    "ad_fwd",
    "dp_rev",
    "rd_rev",
    "ad_rev",
    "dpf",
    "rdf",
    "adf",
]


# ── BAM Construction Helpers ─────────────────────────────────────────────


def build_bam(tmp_path, reads, filename="test.bam"):
    """Write reads to a sorted, indexed BAM. Returns path string.

    Creates a single-contig BAM (chr1, 500bp) from the given AlignedSegments,
    sorts by coordinate, and indexes for random access.
    """
    bam_path = tmp_path / filename
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 500, "SN": "chr1"}]}
    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        for r in reads:
            outf.write(r)
    sorted_bam = tmp_path / filename.replace(".bam", ".sorted.bam")
    pysam.sort("-o", str(sorted_bam), str(bam_path))
    pysam.index(str(sorted_bam))
    return str(sorted_bam)


def make_read(name, seq, start, cigar, flag=0, mapq=60, quals=None):
    """Create an AlignedSegment with sensible defaults.

    Args:
        name: Query name.
        seq: Query sequence string.
        start: 0-based reference start position.
        cigar: CIGAR tuples, e.g. ((0, 10),) for 10M.
        flag: SAM flag (default: 0 = forward, unpaired).
        mapq: Mapping quality (default: 60).
        quals: Base quality array (default: all Q30).
    """
    a = pysam.AlignedSegment()
    a.query_name = name
    a.query_sequence = seq
    a.flag = flag
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = mapq
    a.cigartuples = cigar
    a.query_qualities = quals if quals else [30] * len(seq)  # type: ignore[assignment]
    return a


# ── Counting, checked for binning invariance ─────────────────────────────


def count_bam_checked(*args, **kwargs) -> list:
    """`count_bam_binned` with production bin geometry, checked against one
    variant per bin with the per-variant fetch (window 1, cap 1) through the same
    loop: every field of every row must agree. Takes `count_bam_binned`'s
    arguments; returns the production counts."""
    production = gbcms_rs.count_bam_binned(*args, **kwargs)
    per_variant = gbcms_rs.count_bam_binned(*args, **kwargs, bin_window=1, bin_max_variants=1)
    assert_same_counts(per_variant, production, "one variant per bin vs production bins")
    return production


def count_checked(
    bam_path: str,
    variants: list,
    *,
    decomposed: list | None = None,
    min_mapq: int = 20,
    min_baseq: int = 20,
    filter_duplicates: bool = True,
    filter_secondary: bool = True,
    filter_supplementary: bool = True,
    filter_qc_failed: bool = False,
    filter_improper_pair: bool = False,
    filter_indel: bool = False,
    threads: int = 1,
    **kwargs,
) -> list:
    """`count_bam_checked` with the suite's default filters (duplicates,
    secondary and supplementary dropped; MAPQ and BQ 20). Extra keyword arguments
    go to `count_bam_binned`."""
    return count_bam_checked(
        bam_path,
        variants,
        decomposed if decomposed is not None else [None] * len(variants),
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        filter_duplicates=filter_duplicates,
        filter_secondary=filter_secondary,
        filter_supplementary=filter_supplementary,
        filter_qc_failed=filter_qc_failed,
        filter_improper_pair=filter_improper_pair,
        filter_indel=filter_indel,
        threads=threads,
        **kwargs,
    )


def count_one_checked(bam_path, variant, **kwargs):
    """One variant through `count_checked`."""
    return count_checked(bam_path, [variant], **kwargs)[0]


def write_contig(tmp_path, contig, reads, name="s", chrom="1"):
    """A one-contig FASTA (indexed) and a sorted, indexed BAM of `reads` on it.
    Returns (fasta path, bam path) as strings."""
    fa = tmp_path / f"{name}.fa"
    fa.write_text(f">{chrom}\n{contig}\n")
    pysam.faidx(str(fa))
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": chrom, "LN": len(contig)}]}
    raw = tmp_path / f"{name}.raw.bam"
    with pysam.AlignmentFile(str(raw), "wb", header=hdr) as fh:
        for r in reads:
            fh.write(r)
    bam = tmp_path / f"{name}.bam"
    pysam.sort("-o", str(bam), str(raw))
    pysam.index(str(bam))
    return str(fa), str(bam)


# ── Comparing counts ─────────────────────────────────────────────────────

# Benjamini-Hochberg q-values are computed over the rows in one call, so they
# legitimately change when a call holds other rows.
ROW_SET_FIELDS = frozenset({"mfsd_qval_alt_ref", "asjd_qval"})


def _fields(obj) -> list[str]:
    """The public data attributes of a pyo3 counts object."""
    return sorted(n for n in dir(obj) if not n.startswith("_") and not callable(getattr(obj, n)))


def _same(a, b) -> bool:
    if isinstance(a, float) and isinstance(b, float):
        if math.isnan(a) or math.isnan(b):
            return math.isnan(a) and math.isnan(b)
        # Summation order over hash maps moves the last ulps of some statistics;
        # the outputs round them far coarser. Relative only: p-values can be tiny.
        return math.isclose(a, b, rel_tol=1e-9, abs_tol=0.0)
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_same(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b, strict=True))
    return type(a) is type(b) and a == b


def assert_same_counts(got: list, want: list, what: str, skip=frozenset()) -> None:
    """Every field of every BaseCounts equal (floats to 1e-9 relative, NaN equal
    to NaN), except `skip`."""
    assert len(got) == len(want), f"{what}: {len(got)} rows vs {len(want)}"
    for i, (g, w) in enumerate(zip(got, want, strict=True)):
        for f in _fields(w):
            if f in skip:
                continue
            assert _same(
                getattr(g, f), getattr(w, f)
            ), f"{what}: row {i} field {f!r}: {getattr(g, f)!r} vs {getattr(w, f)!r}"


# ── The suite under the haplotype judge ──────────────────────────────────
#
# With GBCMS_TEST_READ_JUDGE=haplotype in the environment, every engine call and every
# CLI `dna` run that names no read judgment counts under the haplotype judge
# (`--read-judge haplotype`, #304), so the whole suite runs under it with no test edited
# (docs/proposals/haplotype-judge.md §8). The judge reads each row's haplotypes from the
# reference: an engine call given none takes the one FASTA beside its BAM; a call with no
# single FASTA there runs under the rules and is reported as not judged (in the session
# summary, and appended to GBCMS_TEST_READ_JUDGE_LOG when set), never passed silently.
# A call that names a read judgment keeps it.

SUITE_READ_JUDGE = os.environ.get("GBCMS_TEST_READ_JUDGE", "rules").strip().lower()

# Positional slots of `count_bam_binned(_observations)` (see `_rs.pyi`).
_BAM_AT, _REFERENCE_AT = 0, 28
_READ_JUDGE_AT = {"count_bam_binned": 33, "count_bam_binned_observations": 34}

# (test, reason) for every engine call the switch left under the rules.
UNJUDGED: list[tuple[str, str]] = []


def _fasta_beside(bam_path) -> str | None:
    """The one FASTA (.fa / .fasta) in the BAM's directory whose contigs the BAM's header
    names with the same lengths (an extract of a region is not the BAM's reference), or
    None."""
    if not bam_path:
        return None
    d = os.path.dirname(os.path.abspath(str(bam_path)))
    found = sorted(set(glob.glob(os.path.join(d, "*.fa")) + glob.glob(os.path.join(d, "*.fasta"))))
    if len(found) != 1:
        return None
    try:
        with pysam.FastaFile(found[0]) as fa, pysam.AlignmentFile(str(bam_path)) as bam:
            lengths = dict(zip(fa.references, fa.lengths, strict=True))
            shared = [
                (n, ln) for n, ln in zip(bam.references, bam.lengths, strict=True) if n in lengths
            ]
            if not shared or any(lengths[n] != ln for n, ln in shared):
                return None
    except (OSError, ValueError):
        return None
    return found[0]


def _note_unjudged(reason: str) -> None:
    test = os.environ.get("PYTEST_CURRENT_TEST", "?").rsplit(" (", 1)[0]
    UNJUDGED.append((test, reason))
    log = os.environ.get("GBCMS_TEST_READ_JUDGE_LOG")
    if log:
        with open(log, "a") as fh:
            fh.write(f"{test}\t{reason}\n")


def _under_suite_judge(engine):
    """`engine` with calls that name no read judgment run under the haplotype judge."""
    judge_at = _READ_JUDGE_AT[engine.__name__]

    @functools.wraps(engine)
    def call(*args, **kwargs):
        if "read_judge" in kwargs or len(args) > judge_at:
            return engine(*args, **kwargs)
        given = kwargs.get("reference_fasta") or (len(args) > _REFERENCE_AT and args[_REFERENCE_AT])
        if not given:
            bam = args[_BAM_AT] if args else kwargs.get("bam_path")
            fasta = _fasta_beside(bam)
            if fasta is None:
                _note_unjudged("no single FASTA of the BAM's contigs beside it")
                return engine(*args, **kwargs)
            kwargs["reference_fasta"] = fasta
        return engine(*args, **kwargs, read_judge="haplotype")

    # pyo3 publishes the signature the stub-parity test reads.
    call.__text_signature__ = engine.__text_signature__  # type: ignore[attr-defined]
    return call


def install_suite_read_judge() -> None:
    """Run the suite under the haplotype judge: the engine entry points, the pipeline's
    configuration default, and CLI `dna` runs through `CliRunner` that name no read
    judgment."""
    from typer.testing import CliRunner

    import gbcms.observations
    from gbcms.models.core import AlignmentConfig

    for name in _READ_JUDGE_AT:
        wrapped = _under_suite_judge(getattr(gbcms_rs, name))
        setattr(gbcms_rs, name, wrapped)
        if hasattr(gbcms.observations, name):
            setattr(gbcms.observations, name, wrapped)
    AlignmentConfig.model_fields["read_judge"].default = "haplotype"
    AlignmentConfig.model_rebuild(force=True)

    invoke = CliRunner.invoke

    def judged_invoke(self, app, args=None, *rest, **kwargs):
        named = args and any(str(a).startswith("--read-judge") for a in args)
        if args and list(args)[0] == "dna" and not named:
            args = [*args, "--read-judge", "haplotype"]
        return invoke(self, app, args, *rest, **kwargs)

    CliRunner.invoke = judged_invoke  # type: ignore[method-assign]


# ── MAF Output Reading ───────────────────────────────────────────────────


def read_maf_output(path):
    """Read MAF output, skipping #-prefixed provenance comment lines.

    MafWriter emits provenance headers (e.g., ``#gbcms v5.3.0``,
    ``#command ...``) before the TSV header row.  These must be
    skipped for ``csv.DictReader`` to parse the file correctly.

    Returns:
        csv.DictReader positioned at the first data row.
    """
    import csv
    import io

    with open(path) as f:
        lines = [line for line in f if not line.startswith("#")]
    return csv.DictReader(io.StringIO("".join(lines)), delimiter="\t")
