"""
Shared test helpers for gbcms.

Provides:
- BAM construction helpers (_build_bam, _make_read) for synthetic test data
- count_one: single-variant counting via legacy API
- count_one_both: single-variant counting via BOTH APIs with parity assertion
- count_both: multi-variant counting via BOTH APIs with parity assertion
"""

import pysam

from gbcms import _rs as gbcms_rs

# Key BaseCounts fields to compare for count_bam vs count_bam_binned parity.
PARITY_FIELDS = [
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


# ── Single-Variant Counting ──────────────────────────────────────────────


def count_one(bam_path, variant):
    """Count a single variant via the legacy count_bam API.

    Applies default filters (filter_duplicates, filter_secondary, filter_supplementary).
    Returns a single BaseCounts object.
    """
    results = gbcms_rs.count_bam(
        bam_path,
        [variant],
        decomposed=[None],
        min_mapq=20,
        min_baseq=20,
        filter_duplicates=True,
        filter_secondary=True,
        filter_supplementary=True,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
        threads=1,
    )
    return results[0]


def count_one_both(bam_path, variant):
    """Count a single variant via BOTH APIs, assert parity, return production result.

    Calls count_both with a single variant and returns the single BaseCounts object.
    """
    return count_both(
        bam_path,
        [variant],
        min_mapq=20,
        min_baseq=20,
        filter_qc_failed=False,
        filter_improper_pair=False,
        filter_indel=False,
    )[0]


# ── Multi-Variant Parity Counting ────────────────────────────────────────


def count_both(
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
    fragment_qual_threshold: int = 10,
    sibling_variants: list | None = None,
) -> list:
    """Call both count_bam and count_bam_binned, assert parity, return results.

    Returns the count_bam_binned results (production path).
    Raises AssertionError if any parity field differs between the two APIs.
    """
    if decomposed is None:
        decomposed = [None] * len(variants)
    if sibling_variants is None:
        sibling_variants = [[] for _ in variants]

    # Legacy path
    results_legacy = gbcms_rs.count_bam(
        bam_path,
        variants,
        decomposed,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        filter_duplicates=filter_duplicates,
        filter_secondary=filter_secondary,
        filter_supplementary=filter_supplementary,
        filter_qc_failed=filter_qc_failed,
        filter_improper_pair=filter_improper_pair,
        filter_indel=filter_indel,
        threads=threads,
    )

    # Production path
    results_binned = gbcms_rs.count_bam_binned(
        bam_path,
        variants,
        decomposed,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
        filter_duplicates=filter_duplicates,
        filter_secondary=filter_secondary,
        filter_supplementary=filter_supplementary,
        filter_qc_failed=filter_qc_failed,
        filter_improper_pair=filter_improper_pair,
        filter_indel=filter_indel,
        threads=threads,
        fragment_qual_threshold=fragment_qual_threshold,
        sibling_variants=sibling_variants,
    )

    # Assert parity on key fields
    assert len(results_legacy) == len(results_binned), (
        f"Result count mismatch: count_bam={len(results_legacy)}, "
        f"count_bam_binned={len(results_binned)}"
    )

    for i, (leg, bn) in enumerate(zip(results_legacy, results_binned, strict=True)):
        for field in PARITY_FIELDS:
            v_leg = getattr(leg, field)
            v_bn = getattr(bn, field)
            assert v_leg == v_bn, (
                f"Variant {i} field '{field}' mismatch: "
                f"count_bam={v_leg}, count_bam_binned={v_bn}"
            )

    return results_binned


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
