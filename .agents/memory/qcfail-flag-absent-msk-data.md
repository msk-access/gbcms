---
name: qcfail-flag-absent-msk-data
description: "No MSK pipeline stage sets the BAM QC-fail flag (0x200), so gbcms's filter-qc-failed is inert on real data; the filter itself is verified correct."
metadata: 
  node_type: memory
  type: project
  originSessionId: afbf4a49-f216-4b9f-aa60-421bb8c1073c
  modified: 2026-08-25T17:15:05.865Z
---

Flagstat sweeps (2026-08-12) across ACCESS/CMO pipeline stages — uncollapsed_FM,
collapsed_grouped, standard `_cl_aln_srt_MD_IR_FX_BR`, collapsed, duplex, simplex
(~1.2B reads, XS1/XS2 samples) and IMPACT standard BAMs — found **zero** reads with
FLAG 0x200. `--filter-qc-failed` / `--no-filter-qc-failed` therefore changes nothing
on real MSK data today.

The filter itself is correct end-to-end (CLI → FFI → `rust/src/shared/filters.rs`
gated check; no unconditional skip like bug #82 had). Pinned by
`test_qc_failed_filter_contract` in tests/test_filters.py (PR #86): flagged reads
dropped by default, `--no-filter-qc-failed` restores byte-identical counts.

**Why:** avoids re-investigating "does this filter do anything?" and explains why
toggling it never changes MSK results.
**How to apply:** if a future upstream (Illumina chastity passthrough, fgbio
consensus marking) starts setting 0x200, the filter becomes load-bearing — the
contract test already covers it. See [[nextflow-cli-default-divergence]].
