#!/bin/bash
#SBATCH --job-name=gbcms-panel
#SBATCH --partition=cmobic_cpu
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --array=0-19
#SBATCH --output=panel_%A_%a.log
# A regression-panel build as a SLURM array: each task runs every 20th row of the run
# list (docs/development/regression-panel.md). Edit the CONFIG block, then:
#   sbatch submit_panel.sbatch
# Re-submitting skips finished runs.

# ---- CONFIG -----------------------------------------------------------------
SCRIPTS=/path/to/gbcms/scripts/regression_panel   # SLURM runs a spooled copy of this file
TIER=tier_660                     # the tier folder carried from the selection host
RUNS=$TIER/runs.tsv               # or an attribution run list: attrib/runs.tsv
BUILD=6.6.0                       # output label: 6.5.0, 6.6.0, cp_<sha>
GBCMS="singularity exec docker://ghcr.io/msk-access/gbcms:6.6.0 gbcms"
DMP_ROOT=/path/to/dmp_share       # BAMs of runs with root "dmp" (b37)
FORTE_ROOT=/path/to/forte         # BAMs of runs with root "forte" (hg38)
FASTA_B37=/path/to/Homo_sapiens_assembly19.fasta
FASTA_HG38=/path/to/GRCh38.fasta
GTF=/path/to/Homo_sapiens.GRCh38.111.gtf
OUT=panel_out
# -----------------------------------------------------------------------------

export SHARD=${SLURM_ARRAY_TASK_ID:-0} NSHARDS=${SLURM_ARRAY_TASK_COUNT:-1}
bash "$SCRIPTS/run_panel.sh" "$TIER" "$RUNS" "$BUILD" "$GBCMS" "$DMP_ROOT" "$FORTE_ROOT" \
    "$FASTA_B37" "$FASTA_HG38" "$GTF" "$OUT"
