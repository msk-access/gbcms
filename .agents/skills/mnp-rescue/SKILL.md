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
PASS + `MNP_RESCUE_ELIGIBLE` + `partial_alt > ad` + not in a `multi_allelic_group` +
haplotype not shown by any read: `mnp_confirmed_alt == 0`.
- Not `ad == 0`: masked per-position evaluation counts a component carrier whose other
  discriminating base is low-BQ as full ALT; one such read must not block rescue.
- Grouped rows → `outcome=skipped_grouped` (exclusive assignment owns their reads).
- `mnp_confirmed_alt` (engine, internal): MNP ALT reads with every discriminating base read
  (none masked, none N); set only via check_mnp's `Alt(.., confirmed)`. Any such read means
  the BAM shows the annotated allele → `outcome=haplotype_confirmed`, no re-count. (An error
  allowance was tried and rejected — see REJECTED.md: it adopted a germline SNP on ACCESS.)
  Guards against adopting a germline SNP merged into a somatic MNP; cannot help where the
  MNP is absent (fillout of other timepoints/normals).

## Adoption
Best component (highest `ad`, leftmost on ties) must beat the MNP's `ad`; its whole
`BaseCounts` replaces the row's (never graft single fields — mixed classifications break
`dpf >= rdf + adf`), and diagnostics are recomputed via `_diagnostic_flags`. Both count
calls share `_engine_kwargs()` so components classify reads exactly like the main count.

## `gbcms_rescue` (built only by `_format_rescue_audit`)
`method=decomposed;outcome=<rescued|skipped_grouped|haplotype_confirmed|no_improvement|ref_validation_failed>;original_ref=R;original_alt=A;original_partial=P;original_confirmed=C[;adopted=chr:pos(R>A)][;positions=chr:pos(R>A):<ad|ref_fail>+...]`
- Reset for every sample (the prepared list is shared across BAMs). Positions join with `+`,
  never `,` — VCF emits it as the Number=1 `GR` INFO and parsers split commas.
- Labels (flag, audit, logs) use the contig as the output row writes it (`_output_contig`:
  an input MAF's own `Chromosome`, else the internal name).
- `gbcms merge` warns per row when duplex/simplex rescue outcomes differ (counts unchanged).
- Rescued rows: `RESCUED_COMPONENT(chrom:pos:REF>ALT)` appended to `gbcms_diagnostic`, a
  WARNING per row, and a WARNING when `--rescue-mnp` is enabled (counts are replaced — opt-in).
- `no_improvement` is legitimate: partial evidence from indel-disrupted reads
  (complex path: REF + nearby-indel evidence) that no single-base count calls ALT.
  `ref_validation_failed` is an anomaly → WARNING. `ref_fail` marks a failed position.
- `--observations-parquet` keeps the MNP evaluation and `--mfsd-parquet` keeps the MNP's
  coordinates for rescued rows (both logged per sample).

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
