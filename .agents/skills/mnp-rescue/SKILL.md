---
name: mnp-rescue
description: Reference for gbcms MNP rescue and diagnostic flags — rescue intent and candidate gate, row adoption and the gbcms_rescue audit format, --rescue-mnp-threshold semantics, and the strongly-typed gbcms_diagnostic / gbcms_rescue fields. Use when touching the MNP rescue pass, diagnostic-flag computation, or their tests/stubs.
---

# MNP Rescue & Diagnostic Flags

## Intent
For annotated MNPs whose carriers hold only a *component* of the haplotype (SNVs on
different molecules). The engine correctly puts those carriers in `partial_alt`; rescue
reports the best-supported component instead. Opt-in (`--rescue-mnp`) — the flag's
contract is component counting.

## Candidate gate (`pipeline._rescue_mnp_pass`)
PASS + `MNP_RESCUE_ELIGIBLE` + `partial_alt > ad` + not in a `multi_allelic_group`.
- Not `ad == 0`: masked per-position evaluation counts a component carrier whose other
  discriminating base is low-BQ as full ALT; one such read must not block rescue.
- Grouped rows → `outcome=skipped_grouped` (exclusive assignment owns their reads).

## Adoption
Best component (highest `ad`, leftmost on ties) must beat the MNP's `ad`; its whole
`BaseCounts` replaces the row's (never graft single fields — mixed classifications break
`dpf >= rdf + adf`), and diagnostics are recomputed via `_diagnostic_flags`. Both count
calls share `_engine_kwargs()` so components classify reads exactly like the main count.

## `gbcms_rescue` (built only by `_format_rescue_audit`)
`method=decomposed;outcome=<rescued|skipped_grouped|no_improvement|ref_validation_failed>;original_ref=R;original_alt=A;original_partial=P[;adopted=chr:pos(R>A)][;positions=chr:pos(R>A):<ad|ref_fail>,...]`
- Reset for every sample (the prepared list is shared across BAMs).
- `no_improvement` is unreachable for consistent counts (best component ≥
  (partial_alt + ad)/2) → logged WARNING; so is `ref_validation_failed`.

## Threshold (`--rescue-mnp-threshold`, 0.0–1.0)
`MNP_RESCUE_ELIGIBLE` when disc/len ≤ threshold. 1.0 (default) = all MNPs; 0.5 = sparse
only; 0.0 = none (ratio > 0 for valid MNPs). `MNP_DISC_RATIO(n/m)` is always emitted.

## Diagnostic flags
`gbcms_diagnostic`: `;`-separated flags, e.g. `PARTIAL_DOMINANT;MNP_DISC_RATIO(2/5);MNP_RESCUE_ELIGIBLE`.
Both fields are `#[pyo3(get, set)]` on `PreparedVariant` (`rust/src/normalize/types.rs`),
stubbed in `src/gbcms/_rs.pyi`.

## Tests
`tests/test_rescue_mnp.py` (end-to-end battery asserting invariants on written rows +
outcome/audit unit tests), `tests/test_diagnostic_flags.py`,
`tests/test_mnp_concordance.py::TestONPCarrierShapes` (which read shapes are ALT vs partial).
