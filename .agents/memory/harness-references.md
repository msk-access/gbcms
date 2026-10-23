---
name: harness-references
description: "Which reference FASTA each local harness arm needs — the IMPACT/ACCESS panel slices are b37, HG002 and FORTE are GRCh38; brief agents with the right one."
metadata:
  node_type: memory
  type: reference
  originSessionId: a8529157-ff34-4400-8321-c28f0357971d
  modified: 2026-10-17T15:52:58.774Z
---

The plan670 harness mixes two builds; a measurement run against the wrong FASTA fails
REF checks or silently skips rows.

- **b37** — the gate panel's IMPACT/ACCESS slices (`hpc_kit_670/slices_indels`,
  `slices_all`, `slices_std`), normals/ACCESS (`across/`), WES: the reference is
  `~/Library/CloudStorage/OneDrive-MemorialSloanKetteringCancerCenter/Downloads/03Jan2023/nucleo_qc_generation/test_data/Homo_sapiens_assembly19.fasta`
  (contigs `1`, `2`, …). The arm scripts (`gate_pi/run_gate_pi.py` via `obs_ab.FA`) carry it.
- **GRCh38** — HG002 (`giab/`) and FORTE RNA (`hpc_kit_670/slices_forte`, with the
  Ensembl 111 GTF): `~/Downloads/genome.fa` (contigs without `chr`; HG002 rows say
  `chr10`, gbcms normalizes).

**How to apply:** when briefing a measurement agent on a panel arm, name the b37 FASTA
explicitly (an agent defaulted to genome.fa for #280 and had to switch). Traces: see
[[trace-output-wraps]] (COLUMNS=60000). Related: [[other-tools-harness]], [[giab-truth-arm]].
