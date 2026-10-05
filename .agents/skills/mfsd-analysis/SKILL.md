---
name: mfsd-analysis
description: Reference for gbcms mFSD (fragment size distribution) — the opt-in --mfsd flag, MAF/VCF column gating, KS test + LLR statistics, and Rust-native ZSTD Parquet output. Use when touching mfsd.rs, parquet_writer.rs, the mFSD report, or column-count tests.
---

# mFSD Analysis Patterns

## Overview

mFSD (modified Fragment Size Distribution) analysis is **opt-in** via `--mfsd`.

## Column Gating

- Without `--mfsd`: 145 MAF columns (0 `mfsd_*`)
- With `--mfsd`: 186 MAF columns (+41 `mfsd_*`)
- VCF INFO: exactly 13 `MFSD_*` fields added only when `--mfsd` is set
- Columns are **absent** when off, not NA-filled

## Parquet Output

- `--mfsd-parquet` requires `--mfsd` (validated at both CLI and Pydantic)
- Written by `write_fsd_parquet()` in Rust via `arrow`/`parquet` crates
- ZSTD(1) compression: `parquet = { default-features = false, features = ["arrow", "zstd"] }`
- `ref_sizes`/`alt_sizes` are Rust-internal — no `#[pyo3(get)]`
- Called from `pipeline.py` after `count_bam_binned()`, not inside the counting engine

## Statistics (watch-outs)

- KS p-value is exact up to 10⁷ lattice cells (shares of in-band paths, no
  overflow; every realistic class pair), Stephens-corrected asymptotic above. The
  uncorrected series overstated p 1.7–45x for few ALT fragments vs a deep REF.
  Exact p is conservative with tied integer sizes. At small ALT counts the test has
  little power (about 8% at 5 fragments), so a non-significant result is not
  evidence of "no shift" — and never evidence for CH.
- LLR is a fragment-size Gaussian log-ratio (distinct from the PairHMM LLR),
  reported as the mean per fragment (n = `mfsd_*_count`); NaN for an empty class,
  as are empty classes' mean sizes. Closed-form log-ratio, so no ±∞.
- The report grades evidence (LEANS-SOMATIC / NO-SIZE-EVIDENCE / INSUFFICIENT);
  nothing leans CH and gene membership is a note, not a gate. The CH-vs-tumor
  prediction lives outside gbcms (a separate model consuming `fsd.parquet`).

## Key Files
- `rust/src/counting/mfsd.rs`: KS test, LLR, pairwise comparisons
- `rust/src/counting/parquet_writer.rs`: Arrow/ZSTD native Parquet
- `src/gbcms/io/output.py`: MafWriter column gating
- `src/gbcms/report/mfsd_report.py`: HTML report generator
- `tests/test_mfsd_flag.py`: column count assertions
