#!/usr/bin/env bash
# Run one gbcms build over a regression-panel run list (docs/development/regression-panel.md).
#
# Every row of the run list is one gbcms call. A run that exits 0 leaves OUT/BUILD/<run>.ok;
# re-running skips those and retries every other run (a run that failed after writing
# its MAF is retried, not taken as finished). Each attempt is a row of
# OUT/BUILD/times.<shard>.tsv: exit code, wall time and peak memory (GNU time; wall
# time only elsewhere) and its end time, so the latest attempt of a run is its outcome.
# A SLURM array runs the shards side by side: SHARD=$SLURM_ARRAY_TASK_ID (0-based),
# NSHARDS=$SLURM_ARRAY_TASK_COUNT (see submit_panel.sbatch).
#
# usage: run_panel.sh TIER_DIR RUNS_TSV BUILD "GBCMS_CMD" DMP_ROOT FORTE_ROOT \
#                     FASTA_B37 FASTA_HG38 GTF OUT_ROOT [ARM ...]
#   TIER_DIR   the tier folder (runs.tsv, mafs/, forte/); variant paths are relative to it
#   RUNS_TSV   TIER_DIR/runs.tsv, or an attribution run list (attribute.py prepare)
#   BUILD      a label for the output folder: 6.5.0, 6.6.0, cp_<sha>
#   GBCMS_CMD  e.g. "singularity exec docker://ghcr.io/msk-access/gbcms:6.5.0 gbcms"
#              or "/path/to/cp_<sha>/bin/gbcms"
#   DMP_ROOT   holds the BAMs of runs whose root is "dmp" (b37)
#   FORTE_ROOT holds the BAMs of runs whose root is "forte" (hg38, RNA mode with GTF)
#   ARM ...    optional: run only these arms (dna mfsd mq0 sw vcf normal cohort rna)
set -uo pipefail
if [ $# -lt 10 ]; then sed -n '2,24p' "$0"; exit 2; fi
TIER=$1 RUNS=$2 BUILD=$3 GB=$4 DMP=$5 FORTE=$6 FA37=$7 FA38=$8 GTF=$9 OUT=${10}
shift 10
ARMS=" $* "
SHARD=${SHARD:-0} NSHARDS=${NSHARDS:-1}
if ! [ "$SHARD" -ge 0 ] 2>/dev/null || ! [ "$SHARD" -lt "$NSHARDS" ]; then
  echo "SHARD must be 0..NSHARDS-1 (got SHARD=$SHARD NSHARDS=$NSHARDS): a SLURM array is 0-based" >&2
  exit 2
fi
mkdir -p "$OUT/$BUILD"
TIMES="$OUT/$BUILD/times.$SHARD.tsv"
[ -f "$TIMES" ] || printf 'run_id\tarm\texit\twall_s\tmax_rss_kb\tend_epoch\n' > "$TIMES"
TIMER=""
if /usr/bin/time -f 'GBCMS_TIME %e %M' true 2>/dev/null; then TIMER="gnu"; fi

n=-1
tail -n +2 "$RUNS" | while IFS=$'\t' read -r rid tag arm mode root rel variants extra; do
  n=$((n + 1))
  [ $((n % NSHARDS)) -eq "$SHARD" ] || continue
  if [ "$ARMS" != "  " ] && [[ "$ARMS" != *" $arm "* ]]; then continue; fi
  o="$OUT/$BUILD/$rid"
  [ -f "$o.ok" ] && continue
  if [ "$root" = "forte" ]; then bam="$FORTE/$rel" fa="$FA38" ref_args=(--gtf "$GTF")
  else bam="$DMP/$rel" fa="$FA37" ref_args=()
  fi
  fmt=maf
  if [[ " $extra " == *" --format vcf "* ]]; then fmt=vcf; extra=${extra/--format vcf/}; fi
  # shellcheck disable=SC2086  # extra_args is a word list by design
  cmd=($GB "$mode" -v "$TIER/$variants" -b "S:$bam" -f "$fa" -o "$o" --format "$fmt" \
       ${ref_args[@]+"${ref_args[@]}"} $extra)
  start=$(date +%s)
  wall=NA rss=NA
  if [ "$TIMER" = "gnu" ]; then
    /usr/bin/time -f 'GBCMS_TIME %e %M' -o "$o.time" "${cmd[@]}" < /dev/null > "$o.log" 2>&1; rc=$?
    # On failure GNU time writes a status line first: read only the tagged line.
    line=$(grep '^GBCMS_TIME ' "$o.time" 2>/dev/null | tail -n 1)
    [ -n "$line" ] && read -r _ wall rss <<< "$line"
  else
    "${cmd[@]}" < /dev/null > "$o.log" 2>&1; rc=$?
    wall=$(( $(date +%s) - start ))
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$rid" "$arm" "$rc" "$wall" "$rss" "$(date +%s)" >> "$TIMES"
  if [ "$rc" -eq 0 ]; then touch "$o.ok"; fi
done
echo "shard $SHARD/$NSHARDS done: $(find "$OUT/$BUILD" -maxdepth 1 -name '*.ok' | wc -l) runs finished in this build"
