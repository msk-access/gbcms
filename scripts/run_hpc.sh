#!/bin/bash
#SBATCH --job-name=gbcms
#SBATCH --partition=cmobic_cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=gbcms_%j.log
#SBATCH --error=gbcms_%j.err

# ============================================================================
# gbcms Nextflow Pipeline — HPC (SLURM + Singularity) submission template
# ============================================================================
# Usage:
#   Edit the CONFIG block below (or export the vars), then:
#     sbatch scripts/run_hpc.sh
#   Resume a failed/interrupted run:
#     sbatch scripts/run_hpc.sh --resume
#
# The head job (this script) is a lightweight Nextflow orchestrator on
# cmobic_cpu (24h). Nextflow submits the per-sample child jobs (GBCMS_DNA/RNA,
# merge, …) to cmobic_short → cmobic_cpu (retries) per nextflow.config.
# The container is pinned per-module to gbcms:<version> and pulled via
# Singularity by the `iris` profile.
#
# INPUTS
#   BAM list  — EITHER a Nextflow samplesheet CSV (header: sample,bam[,bai])
#               OR a 2-column TSV "sample<TAB>bam" (e.g. bam_list.tsv). A TSV is
#               auto-converted to a CSV samplesheet here; BAI/CRAI indexes are
#               auto-discovered by the pipeline.
#   VARIANTS  — VCF or MAF variant list to genotype.
#   FASTA     — reference FASTA; a matching <fasta>.fai MUST sit beside it.
#   MODE      — dna | rna.
# ============================================================================

set -euo pipefail

# ── CONFIG (edit here, or override via environment) ──────────────────────────
# Pipeline entry defaults to this repo's own nextflow/main.nf (script lives in
# scripts/), so it needs no editing when run from a clone.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GBCMS_REPO="${GBCMS_REPO:-${SCRIPT_DIR}/../nextflow/main.nf}"

INPUT="${GBCMS_INPUT:-/path/to/bam_list.tsv}"        # CSV samplesheet OR sample<TAB>bam TSV
VARIANTS="${GBCMS_VARIANTS:-/path/to/input_mutations.maf}"   # VCF or MAF
FASTA="${GBCMS_FASTA:-/path/to/Homo_sapiens_assembly19.fasta}"  # needs <fasta>.fai beside it
MODE="${GBCMS_MODE:-dna}"                            # dna | rna
FORMAT="${GBCMS_FORMAT:-maf}"                        # vcf | maf
PROFILE="${GBCMS_PROFILE:-iris}"                     # iris | slurm | singularity
PRESERVE_BARCODE="${GBCMS_PRESERVE_BARCODE:-true}"   # keep original Tumor_Sample_Barcode (MAF→MAF)
OUTDIR_BASE="${GBCMS_OUTDIR:-$PWD}"                  # results parent dir
# RNA-only (used when MODE=rna):
GTF="${GBCMS_GTF:-}"                                 # optional Ensembl/GENCODE GTF
# ─────────────────────────────────────────────────────────────────────────────

GBCMS_CONFIG="$(dirname "${GBCMS_REPO}")/nextflow.config"
[[ -f "${GBCMS_CONFIG}" ]] || { echo "ERROR: pipeline config not found at ${GBCMS_CONFIG} — set GBCMS_REPO to <clone>/nextflow/main.nf" >&2; exit 1; }
# Parse the manifest version only (anchored ^version to skip custom_config_version).
GBCMS_VERSION=$(grep -m1 -E "^[[:space:]]*version[[:space:]]*=" "${GBCMS_CONFIG}" | sed "s/.*= *['\"]//;s/['\"].*//" | tr -d '[:space:]')
echo ">>> gbcms pipeline version: v${GBCMS_VERSION:-unknown}"

# ── Preflight (fail fast, before any child job is submitted) ──────────────────
for f in "${INPUT}" "${VARIANTS}" "${FASTA}"; do
    [[ -f "${f}" ]] || { echo "ERROR: input not found: ${f}" >&2; exit 1; }
done
[[ -f "${FASTA}.fai" ]] || { echo "ERROR: FASTA index missing: ${FASTA}.fai  (run: samtools faidx ${FASTA})" >&2; exit 1; }
[[ "${MODE}" == "dna" || "${MODE}" == "rna" ]] || { echo "ERROR: MODE must be dna|rna, got '${MODE}'" >&2; exit 1; }

OUTDIR="${OUTDIR_BASE}/gbcms_v${GBCMS_VERSION}_${MODE}"
mkdir -p "${OUTDIR}"

# ── Normalise the BAM list into a Nextflow samplesheet (sample,bam) ──────────
# Accepts an existing CSV samplesheet (starts with "sample,") as-is; otherwise
# treats INPUT as a 2-column "sample<TAB>bam" TSV and builds the CSV.
SAMPLESHEET="${OUTDIR}/samplesheet.csv"
if head -1 "${INPUT}" | grep -qiE '^sample,'; then
    cp "${INPUT}" "${SAMPLESHEET}"
    echo ">>> Using CSV samplesheet as provided ($(($(wc -l < "${SAMPLESHEET}")-1)) samples)"
else
    { echo "sample,bam"; awk -F'\t' 'NF>=2 && $1 !~ /^#/ {print $1","$2}' "${INPUT}"; } > "${SAMPLESHEET}"
    echo ">>> Converted TSV → CSV samplesheet: ${SAMPLESHEET} ($(($(wc -l < "${SAMPLESHEET}")-1)) samples)"
fi

# ── Nextflow environment + Singularity cache on fast/local storage ───────────
eval "$(micromamba shell hook --shell bash)"
micromamba activate nf-env
export SINGULARITY_CACHEDIR="${SINGULARITY_CACHEDIR:-$PWD/.singularity_cache}"
export SINGULARITY_TMPDIR="${SINGULARITY_TMPDIR:-$PWD/.singularity_tmp}"
mkdir -p "${SINGULARITY_CACHEDIR}" "${SINGULARITY_TMPDIR}"

# ── Optional flags ────────────────────────────────────────────────────────────
RESUME_FLAG=""; [[ "${1:-}" == "--resume" ]] && { RESUME_FLAG="-resume"; echo ">>> Resuming previous run..."; }
PB_FLAG="";     [[ "${PRESERVE_BARCODE}" == "true" ]] && PB_FLAG="--preserve_barcode"
GTF_FLAG="";    [[ "${MODE}" == "rna" && -n "${GTF}" ]] && GTF_FLAG="--gtf ${GTF}"

echo ">>> Starting gbcms v${GBCMS_VERSION} (${MODE}, profile=${PROFILE}) at $(date)"
echo ">>> Variants: ${VARIANTS}"
echo ">>> Output:   ${OUTDIR}"

nextflow run "${GBCMS_REPO}" \
  --input     "${SAMPLESHEET}" \
  --variants  "${VARIANTS}" \
  --fasta     "${FASTA}" \
  --mode      "${MODE}" \
  --format    "${FORMAT}" \
  --outdir    "${OUTDIR}" \
  ${PB_FLAG} \
  ${GTF_FLAG} \
  -profile    "${PROFILE}" \
  ${RESUME_FLAG}

echo ">>> gbcms v${GBCMS_VERSION} completed at $(date)"
echo ">>> Results: ${OUTDIR}/gbcms/"
