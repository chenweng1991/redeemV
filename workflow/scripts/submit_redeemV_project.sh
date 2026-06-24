#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=scripts/common.sh
source "$SCRIPT_DIR/common.sh"

usage() {
  cat >&2 <<'USAGE'
Usage:
  submit_redeemV_project.sh PROJECT_ENV SAMPLES_TSV [SAMPLE ...]

If no SAMPLE names are supplied, all samples in SAMPLES_TSV are submitted.

Set DRY_RUN=1 to print the sbatch commands without submitting.
USAGE
}

[ "${1:-}" != "-h" ] || { usage; exit 0; }
[ "$#" -ge 2 ] || { usage; exit 1; }

PROJECT_ENV=$(realpath "$1")
SAMPLES_TSV=$(realpath "$2")
shift 2

source_project_env "$PROJECT_ENV"
require_file "$SAMPLES_TSV"

ACCOUNT="${ACCOUNT:-your_slurm_account}"
PARTITION="${PARTITION:-20}"
MEM_PER_CPU_STEP1="${MEM_PER_CPU_STEP1:-8gb}"
MEM_PER_CPU_STEP2="${MEM_PER_CPU_STEP2:-8gb}"
MEM_PER_CPU_STEP3="${MEM_PER_CPU_STEP3:-16G}"
MEM_PER_CPU_STEP4="${MEM_PER_CPU_STEP4:-8gb}"
THREADS_STEP1="${THREADS_STEP1:-24}"
PREPROCESS_CORES="${PREPROCESS_CORES:-48}"
FILES_PER_CORE="${FILES_PER_CORE:-1}"
CONSENSUS_ARRAY_TASKS=$(consensus_array_tasks)
DRY_RUN="${DRY_RUN:-0}"

run_sbatch() {
  if [ "$DRY_RUN" = "1" ]; then
    printf 'DRY_RUN sbatch' >&2
    printf ' %q' "$@" >&2
    printf '\n' >&2
    echo "DRYRUN_JOB"
  else
    sbatch --parsable "$@"
  fi
}

if [ "$#" -gt 0 ]; then
  samples=("$@")
else
  mapfile -t samples < <(awk -F '\t' 'NR > 1 && $1 != "" {print $1}' "$SAMPLES_TSV")
fi

[ "${#samples[@]}" -gt 0 ] || die "No samples to submit"

for sample in "${samples[@]}"; do
  load_sample_row "$SAMPLES_TSV" "$sample"
  out_dir=$(sample_output_dir "$sample")
  mkdir -p "$out_dir/Log/Array.log"

  export_arg="ALL,PROJECT_ENV=$PROJECT_ENV,SAMPLES_TSV=$SAMPLES_TSV,SAMPLE=$sample"

  jid1=$(run_sbatch \
    --job-name="rd_${sample}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$THREADS_STEP1" \
    --mem-per-cpu="$MEM_PER_CPU_STEP1" \
    --partition="$PARTITION" \
    --account="$ACCOUNT" \
    --output="$out_dir/Log/%x_%j.log" \
    --export="$export_arg" \
    "$SCRIPT_DIR/redeem_step1_map.sh")

  jid2=$(run_sbatch \
    --dependency="afterok:$jid1" \
    --job-name="prp_${sample}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$PREPROCESS_CORES" \
    --mem-per-cpu="$MEM_PER_CPU_STEP2" \
    --partition="$PARTITION" \
    --account="$ACCOUNT" \
    --output="$out_dir/Log/%x_%j.log" \
    --export="$export_arg" \
    "$SCRIPT_DIR/redeem_step2_preprocess.sh")

  jid3=$(run_sbatch \
    --dependency="afterok:$jid2" \
    --array="1-${CONSENSUS_ARRAY_TASKS}" \
    --job-name="S3_${sample}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task=1 \
    --mem-per-cpu="$MEM_PER_CPU_STEP3" \
    --partition="$PARTITION" \
    --account="$ACCOUNT" \
    --output="$out_dir/Log/Array.log/ARRAY-%A-%a.out" \
    --export="$export_arg" \
    "$SCRIPT_DIR/redeem_step3_array_consensus.sh")

  jid4=$(run_sbatch \
    --dependency="afterok:$jid3" \
    --job-name="fin_${sample}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="$THREADS_STEP1" \
    --mem-per-cpu="$MEM_PER_CPU_STEP4" \
    --partition="$PARTITION" \
    --account="$ACCOUNT" \
    --output="$out_dir/Log/%x_%j.log" \
    --export="$export_arg" \
    "$SCRIPT_DIR/redeem_step4_finalize.sh")

  printf '%s\tstep1=%s\tstep2=%s\tstep3=%s\tstep4=%s\n' "$sample" "$jid1" "$jid2" "$jid3" "$jid4"
done
