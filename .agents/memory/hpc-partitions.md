---
name: hpc-partitions
description: MSK HPC (SLURM) rules for gbcms runs — never run work on the head node; cmobic_short (≤3 h, priority QOS) for short interactive work, cmobic_cpu (7-day limit) for anything longer.
metadata:
  type: reference
---

The operator's SLURM cluster:
- **Never run work on the head node.** Setup and analysis scripts run in an
  interactive allocation: `srun --pty bash` (or `salloc`, then check `hostname`, since
  a plain `salloc` can leave the shell on the head node). Submitting with `sbatch` from
  the head node is fine; it only queues the job.
- **`cmobic_short`**: up to 3 h, used with `--qos=priority` for interactive work
  (the operator's form: `--mem=40G --time=02:59:00 --partition=cmobic_short --qos=priority`).
- **`cmobic_cpu`**: 7-day limit, for anything longer, interactive or batch. The
  regression panel's `submit_panel.sbatch` and the Nextflow retry tier use it.

Whether compute nodes reach GitHub/PyPI is not yet known (2026-10-07); an offline
wheelhouse is the fallback for installs. Related: [[long-runs-visible]].
