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
  (`count_bam` is the feature-gated legacy parity oracle, not the production path)

## Statistics (watch-outs)

- KS p-value uses an asymptotic approximation — weak at the small ALT-fragment
  counts (n≈5–20) typical of low-input cfDNA.
- LLR is a fragment-size Gaussian log-ratio (distinct from the PairHMM LLR);
  guard against ±∞ from tail fragments by using the closed-form log-ratio.

## Key Files
- `rust/src/counting/mfsd.rs`: KS test, LLR, pairwise comparisons
- `rust/src/counting/parquet_writer.rs`: Arrow/ZSTD native Parquet
- `src/gbcms/io/output.py`: MafWriter column gating
- `src/gbcms/report/mfsd_report.py`: HTML report generator
- `tests/test_mfsd_flag.py`: column count assertions
