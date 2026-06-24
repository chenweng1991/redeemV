#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=scripts/common.sh
source "$SCRIPT_DIR/common.sh"

[ -n "${PROJECT_ENV:-}" ] || die "PROJECT_ENV is not set"
[ -n "${SAMPLES_TSV:-}" ] || die "SAMPLES_TSV is not set"
[ -n "${SAMPLE:-}" ] || die "SAMPLE is not set"
[ -n "${SLURM_ARRAY_TASK_ID:-}" ] || die "SLURM_ARRAY_TASK_ID is not set"

source_project_env "$PROJECT_ENV"
activate_runtime_env
load_sample_row "$SAMPLES_TSV" "$SAMPLE"

OUT_DIR=$(sample_output_dir "$SAMPLE")
WD="$OUT_DIR/Out_mitoConsensus"
mkdir -p "$WD"
cd "$OUT_DIR"

python "$MITOCONSENSUS/mitoConsensus.py" \
  "barcodes.${SLURM_ARRAY_TASK_ID}" \
  "$WD" \
  "$BARCODE_TAG" \
  "$MIN_BASE_QUALITY"
