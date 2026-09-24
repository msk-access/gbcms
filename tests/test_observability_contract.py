"""Target contract for the 6.5.0 observability bundle (T3/T4/T5).

Three signals that were silent become visible, without new output columns:

- T5 ``SW_FALLBACK(n)`` (gbcms_diagnostic) + a per-variant WARN: under the
  PairHMM backend, reads that could not be scored by the pangenomic
  haplotype matrix fall back to Smith-Waterman. Measured to never fire on
  well-formed input; when it does, upstream input was malformed and the
  counts partly came from a different scorer (or, if SW cannot run either,
  the reads were lost as NEITHER) — the row must say so.
- T3 ``CLIP_CANDIDATES(n)`` (gbcms_diagnostic): at an insertion locus with
  no confirmed ALT, soft clips (>= 8bp) whose clip boundary lies within the
  insert's duplication reach point at carriers the aligner represented as
  clips rather than I ops.
- T4: ``--umi-tag TAG`` that no processed read carries logs one WARN per
  BAM (fragment grouping silently fell back to QNAME).

Both variant input paths (VCF and MAF) are exercised for the CLI geometry.
Committed red (xfail-strict) before the implementation; flipped green with it.
"""

import glob
import logging
import random
import types

import pysam
from helpers import make_read, read_maf_output
from typer.testing import CliRunner

from gbcms import _rs
from gbcms.cli import app
from gbcms.pipeline import Pipeline

runner = CliRunner()

READ_LEN = 100
ENGINE_LOGGER = "_rs.counting.engine"


def _bam(tmp_path, ref, reads, name="obs.bam"):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "1", "LN": len(ref)}]}
    p = tmp_path / name
    with pysam.AlignmentFile(p, "wb", header=header) as fh:
        for a in sorted(reads, key=lambda a: a.reference_start):
            fh.write(a)
    sp = tmp_path / name.replace(".bam", ".s.bam")
    pysam.sort("-o", str(sp), str(p))
    pysam.index(str(sp))
    return sp


def _binned(bam, variants, backend="pairhmm", umi_tag=None):
    return _rs.count_bam_binned(
        str(bam),
        variants,
        [None] * len(variants),
        20,
        20,
        True,
        False,
        False,
        False,
        False,
        False,
        1,
        alignment_backend=backend,
        umi_tag=umi_tag,
    )


# ═════════════ T5: SW fallback under the PairHMM backend ════════════════
SW_POS = 300


def _sw_setup(tmp_path):
    """A delins whose ref_context window ends exactly at the variant, so the
    pangenomic matrix cannot place it; 6 reads carry a D2 inside the window
    (routed to Phase-0 realignment), 4 are clean."""
    rng = random.Random(5)
    ref = "".join(rng.choice("ACGT") for _ in range(600))
    reads = []
    for i in range(6):
        s = 240 + i
        left = 285 - s
        seq = ref[s:285] + ref[287 : 287 + READ_LEN - left]
        reads.append(make_read(f"d{i}", seq, s, ((0, left), (2, 2), (0, READ_LEN - left))))
    for i in range(4):
        reads.append(make_read(f"w{i}", ref[250 + i : 350 + i], 250 + i, ((0, READ_LEN),)))
    bam = _bam(tmp_path, ref, reads)

    def variant(ctx_start, ctx_end):
        return _rs.Variant(
            "1", SW_POS, ref[SW_POS : SW_POS + 3], "T", "COMPLEX", ref[ctx_start:ctx_end], ctx_start
        )

    return bam, variant(270, 300), variant(270, 320)


def test_sw_fallback_is_counted_and_warned(tmp_path, caplog):
    bam, misplaced, _ = _sw_setup(tmp_path)
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        (c,) = _binned(bam, [misplaced])
    assert c.sw_fallback_reads == 6
    warns = [r for r in caplog.records if r.levelno == logging.WARNING and "fallback" in r.message]
    assert len(warns) == 1, [r.message for r in caplog.records]
    assert "1:301" in warns[0].message and "6" in warns[0].message


def test_explicit_sw_backend_is_never_a_fallback(tmp_path):
    """SW chosen via --alignment-backend sw is the primary scorer, not a fallback."""
    bam, misplaced, _ = _sw_setup(tmp_path)
    (c,) = _binned(bam, [misplaced], backend="sw")
    assert c.sw_fallback_reads == 0


def test_well_formed_context_never_falls_back(tmp_path):
    bam, _, good = _sw_setup(tmp_path)
    (c,) = _binned(bam, [good])
    assert c.sw_fallback_reads == 0


def test_sw_fallback_counter_parity_with_legacy(tmp_path):
    bam, misplaced, _ = _sw_setup(tmp_path)
    (b,) = _binned(bam, [misplaced])
    (legacy,) = _rs.count_bam(
        str(bam), [misplaced], [None], 20, 20, True, False, False, False, False, False, 1
    )
    assert b.sw_fallback_reads == legacy.sw_fallback_reads == 6


def test_sw_fallback_counts_only_depth_reads(tmp_path):
    """Two more fallback-routed reads whose alignment ends before the anchor:
    they are classified (they overlap the scan window) but never contribute
    to DP/RD/AD, so they must not count toward SW_FALLBACK — the flag claims
    the row's counts came partly from a different scorer."""
    bam, misplaced, _ = _sw_setup(tmp_path)
    rng = random.Random(5)
    ref = "".join(rng.choice("ACGT") for _ in range(600))
    reads = [
        make_read(a.query_name, a.query_sequence, a.reference_start, a.cigartuples)
        for a in pysam.AlignmentFile(str(bam))
    ]
    for i in range(2):
        s = 196 + i
        left = 285 - s
        right = READ_LEN - left
        seq = ref[s:285] + ref[287 : 287 + right]
        reads.append(make_read(f"na{i}", seq, s, ((0, left), (2, 2), (0, right))))
    assert all(a.reference_end < SW_POS + 1 for a in reads[-2:]), "must end before the anchor"
    bam2 = _bam(tmp_path, ref, reads, name="obs_na.bam")
    (c,) = _binned(bam2, [misplaced])
    assert c.dp == 10
    assert c.sw_fallback_reads == 6


def _no_context_setup(tmp_path):
    """A delins with NO reference context (prep's context fetch failed):
    5 reads carry a 1bp insertion inside the span (a 4-base reconstruction,
    neither allele's length) and must escalate to Phase 3, where no scorer
    can run; 4 reads are clean."""
    rng = random.Random(5)
    ref = "".join(rng.choice("ACGT") for _ in range(600))
    reads = []
    for i in range(5):
        s = 240 + i
        left = SW_POS + 1 - s
        right = READ_LEN - left - 1
        seq = ref[s : SW_POS + 1] + "A" + ref[SW_POS + 1 : SW_POS + 1 + right]
        reads.append(make_read(f"i{i}", seq, s, ((0, left), (1, 1), (0, right))))
    for i in range(4):
        reads.append(make_read(f"w{i}", ref[250 + i : 350 + i], 250 + i, ((0, READ_LEN),)))
    variant = _rs.Variant("1", SW_POS, ref[SW_POS : SW_POS + 3], "T", "COMPLEX")
    return _bam(tmp_path, ref, reads, name="obs_nc.bam"), variant


def test_missing_context_reads_are_flagged_not_silent(tmp_path, caplog):
    """Without a reference context the PairHMM matrix, and SW, cannot run:
    reads needing Phase 3 end NEITHER. That loss must be counted and warned
    (with the real reason), not silent — in both the binned and legacy paths."""
    bam, variant = _no_context_setup(tmp_path)
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        (c,) = _binned(bam, [variant])
    assert c.sw_fallback_reads == 5
    warns = [r.message for r in caplog.records if "SW_FALLBACK(5)" in r.message]
    assert len(warns) == 1 and "no reference context" in warns[0], warns
    (legacy,) = _rs.count_bam(
        str(bam), [variant], [None], 20, 20, True, False, False, False, False, False, 1
    )
    assert legacy.sw_fallback_reads == 5
    (sw,) = _binned(bam, [variant], backend="sw")
    assert sw.sw_fallback_reads == 0, "explicit SW backend: SW is chosen, not fallen into"


def _diag(counts_kwargs, ref_allele="GGG", alt_allele="T"):
    pv = types.SimpleNamespace(
        variant=types.SimpleNamespace(ref_allele=ref_allele, alt_allele=alt_allele),
        gbcms_status="PASS",
        gbcms_status_reason="",
        gbcms_diagnostic="",
    )
    base = {"dp": 10, "ad": 2, "partial_alt": 0, "any_alt": 2, "n_count": 0}
    base.update(counts_kwargs)
    pipe = object.__new__(Pipeline)
    pipe.config = types.SimpleNamespace(rescue_mnp_threshold=1.0)
    pipe._compute_diagnostics([pv], [types.SimpleNamespace(**base)])
    return pv.gbcms_diagnostic.split(";") if pv.gbcms_diagnostic else []


def test_sw_fallback_flag():
    assert "SW_FALLBACK(6)" in _diag({"sw_fallback_reads": 6})


def test_no_sw_fallback_flag_when_zero():
    assert not any(f.startswith("SW_FALLBACK") for f in _diag({"sw_fallback_reads": 0}))


# ═════════════ T3: clip-borne insertion carriers ════════════════════════
ANCHOR = 300  # 0-based anchor; a 30bp tandem duplication inserted after it
DUP = 30


def _itd_ref():
    rng = random.Random(41)
    ref = "".join(rng.choice("ACGT") for _ in range(900))
    # already left-aligned: the base before the insert differs from the
    # dup's last base, so prep cannot shift the annotation
    assert ref[ANCHOR] != ref[ANCHOR + DUP]
    return ref


def _itd_rows(ref):
    dup = ref[ANCHOR + 1 : ANCHOR + 1 + DUP]
    return [(ANCHOR + 1, ref[ANCHOR], ref[ANCHOR] + dup)]


def _clip_carriers(ref, n_right=5, n_left=3):
    """Carrier molecules (ref[..A] + dup + ref[A+1..]) aligned clip-only:
    right-clipped reads align through the first copy and clip at the
    duplication boundary (A+1+DUP); left-clipped reads align from A+1 and
    clip the prefix + inserted copy."""
    dup = ref[ANCHOR + 1 : ANCHOR + 1 + DUP]
    out = []
    for i in range(n_right):
        s = ANCHOR - 40 + i
        m = ANCHOR + 1 + DUP - s
        seq = ref[s : ANCHOR + 1] + dup + ref[ANCHOR + 1 : ANCHOR + 1 + (READ_LEN - m)]
        out.append(make_read(f"cr{i}", seq, s, ((0, m), (4, READ_LEN - m))))
    for i in range(n_left):
        s = ANCHOR - 10 - i
        clip = (ANCHOR + 1 - s) + DUP
        seq = ref[s : ANCHOR + 1] + dup + ref[ANCHOR + 1 : ANCHOR + 1 + (READ_LEN - clip)]
        out.append(make_read(f"cl{i}", seq, ANCHOR + 1, ((4, clip), (0, READ_LEN - clip))))
    return out


def _wt(ref, center, n, prefix="wt"):
    return [
        make_read(
            f"{prefix}{i}",
            ref[center - 50 + i : center + 50 + i],
            center - 50 + i,
            ((0, READ_LEN),),
        )
        for i in range(n)
    ]


def _fasta(tmp_path, ref):
    fa = tmp_path / "ref.fasta"
    fa.write_text(">chr1\n" + ref + "\n")
    pysam.faidx(str(fa))
    return fa


def _vcf(tmp_path, rows):
    vcf = tmp_path / "v.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "".join(f"chr1\t{p}\t.\t{r}\t{a}\t.\t.\t.\n" for p, r, a in rows)
    )
    return vcf


def _maf(tmp_path, rows):
    """MAF twin: an anchor-preserved insertion becomes Start = anchor,
    End = anchor + 1, REF = '-', ALT = inserted bases."""
    maf = tmp_path / "v.maf"
    lines = [
        "Hugo_Symbol\tChromosome\tStart_Position\tEnd_Position\t"
        "Reference_Allele\tTumor_Seq_Allele2\tTumor_Sample_Barcode"
    ]
    for p, r, a in rows:
        if len(r) == 1 and len(a) > 1 and a[0] == r:
            lines.append(f"G\tchr1\t{p}\t{p + 1}\t-\t{a[1:]}\tS")
        else:
            lines.append(f"G\tchr1\t{p}\t{p + len(r) - 1}\t{r}\t{a}\tS")
    maf.write_text("\n".join(lines) + "\n")
    return maf


def _run(tmp_path, variants, bam, fasta, outname="out"):
    outdir = tmp_path / outname
    outdir.mkdir(exist_ok=True)
    result = runner.invoke(
        app,
        [
            "dna",
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
        ],
    )
    assert result.exit_code == 0, result.output
    rows = list(read_maf_output(glob.glob(str(outdir / "*.maf"))[0]))
    for r in rows:
        assert int(r["total_count"]) >= int(r["ref_count"]) + int(r["alt_count"])
        assert int(r["any_alt"]) == int(r["alt_count"]) + int(r["partial_alt"])
    return rows


def _diag_flags(row):
    return [f for f in row["gbcms_diagnostic"].split(";") if f]


def test_clip_only_insertion_carriers_are_flagged(tmp_path):
    ref = _itd_ref()
    bam = _bam(tmp_path, ref, _clip_carriers(ref) + _wt(ref, ANCHOR, 10))
    (r,) = _run(tmp_path, _vcf(tmp_path, _itd_rows(ref)), bam, _fasta(tmp_path, ref))
    assert int(r["alt_count"]) == 0, "precondition: clip-only carriers are not ALT"
    assert "CLIP_CANDIDATES(8)" in _diag_flags(r), r["gbcms_diagnostic"]


def test_clip_candidates_agree_across_input_paths(tmp_path):
    """Guard: VCF and MAF front doors agree on counts and diagnostics."""
    ref = _itd_ref()
    bam = _bam(tmp_path, ref, _clip_carriers(ref) + _wt(ref, ANCHOR, 10))
    fasta, rows = _fasta(tmp_path, ref), _itd_rows(ref)
    (rv,) = _run(tmp_path, _vcf(tmp_path, rows), bam, fasta, "out_v")
    (rm,) = _run(tmp_path, _maf(tmp_path, rows), bam, fasta, "out_m")
    for col in ("ref_count", "alt_count", "total_count", "partial_alt", "gbcms_diagnostic"):
        assert rv[col] == rm[col], (col, rv[col], rm[col])


def test_insertion_with_confirmed_alt_is_never_flagged(tmp_path):
    """Guard: I-op carriers present (ad > 0) — stray clips do not flag."""
    ref = _itd_ref()
    dup = ref[ANCHOR + 1 : ANCHOR + 1 + DUP]
    i_carriers = []
    for i in range(6):
        s = ANCHOR - 40 + i
        left = ANCHOR + 1 - s
        seq = ref[s : ANCHOR + 1] + dup + ref[ANCHOR + 1 : ANCHOR + 1 + (READ_LEN - left - DUP)]
        i_carriers.append(
            make_read(f"ic{i}", seq, s, ((0, left), (1, DUP), (0, READ_LEN - left - DUP)))
        )
    bam = _bam(tmp_path, ref, i_carriers + _clip_carriers(ref, 2, 1) + _wt(ref, ANCHOR, 10))
    (r,) = _run(tmp_path, _vcf(tmp_path, _itd_rows(ref)), bam, _fasta(tmp_path, ref))
    assert int(r["alt_count"]) > 0
    assert not any(f.startswith("CLIP_CANDIDATES") for f in _diag_flags(r)), r["gbcms_diagnostic"]


def test_snv_locus_is_never_flagged(tmp_path):
    """Guard: clip candidates are an insertion concept only."""
    ref = _itd_ref()
    alt = "A" if ref[ANCHOR] != "A" else "C"
    bam = _bam(tmp_path, ref, _clip_carriers(ref) + _wt(ref, ANCHOR, 10))
    (r,) = _run(
        tmp_path, _vcf(tmp_path, [(ANCHOR + 1, ref[ANCHOR], alt)]), bam, _fasta(tmp_path, ref)
    )
    assert not any(f.startswith("CLIP_CANDIDATES") for f in _diag_flags(r)), r["gbcms_diagnostic"]


def test_clips_beyond_duplication_reach_are_not_candidates(tmp_path):
    """Guard: clips 150bp away from the insertion are unrelated."""
    ref = _itd_ref()
    far = []
    for i in range(6):
        s = ANCHOR + 100 + i
        far.append(make_read(f"far{i}", ref[s : s + READ_LEN], s, ((0, 70), (4, 30))))
    bam = _bam(tmp_path, ref, far + _wt(ref, ANCHOR, 10))
    (r,) = _run(tmp_path, _vcf(tmp_path, _itd_rows(ref)), bam, _fasta(tmp_path, ref))
    assert not any(f.startswith("CLIP_CANDIDATES") for f in _diag_flags(r)), r["gbcms_diagnostic"]


# ═════════════ T4: --umi-tag never seen ═════════════════════════════════
def _umi_setup(tmp_path, with_tag):
    rng = random.Random(9)
    ref = "".join(rng.choice("ACGT") for _ in range(600))
    reads = _wt(ref, 300, 8)
    if with_tag:
        for i, a in enumerate(reads):
            a.set_tag("RX", f"ACGT{i:04d}", value_type="Z")
    bam = _bam(tmp_path, ref, reads)
    v = _rs.Variant("1", 300, ref[300], "A" if ref[300] != "A" else "C", "SNP")
    return bam, v


def test_umi_tag_never_seen_warns_once(tmp_path, caplog):
    bam, v = _umi_setup(tmp_path, with_tag=False)
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _binned(bam, [v], umi_tag="RX")
    warns = [r for r in caplog.records if r.levelno == logging.WARNING and "RX" in r.message]
    assert len(warns) == 1, [r.message for r in caplog.records]


def test_umi_tag_present_does_not_warn(tmp_path, caplog):
    bam, v = _umi_setup(tmp_path, with_tag=True)
    with caplog.at_level(logging.WARNING, logger=ENGINE_LOGGER):
        _binned(bam, [v], umi_tag="RX")
    assert not [r for r in caplog.records if r.levelno == logging.WARNING and "RX" in r.message]


def test_missing_umi_tag_leaves_counts_unchanged(tmp_path):
    """Guard: the fallback to QNAME grouping is behavior-preserving."""
    bam, v = _umi_setup(tmp_path, with_tag=False)
    (with_flag,) = _binned(bam, [v], umi_tag="RX")
    (without,) = _binned(bam, [v])
    for f in ("dp", "rd", "ad", "dpf", "rdf", "adf"):
        assert getattr(with_flag, f) == getattr(without, f), f
