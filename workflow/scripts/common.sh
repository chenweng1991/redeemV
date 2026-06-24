#!/bin/bash
set -euo pipefail

die() {
  echo "ERROR: $*" >&2
  exit 1
}

require_file() {
  local path="$1"
  [ -f "$path" ] || die "Missing file: $path"
}

source_project_env() {
  local env_file="$1"
  require_file "$env_file"
  # shellcheck disable=SC1090
  source "$env_file"

  REDEEM_V="${REDEEM_V:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
  MYMULTIOME="${MYMULTIOME:-$REDEEM_V/MyMultiome}"
  MITOCONSENSUS="${MITOCONSENSUS:-$REDEEM_V/mitoConsensus}"
  CUT="${CUT:-0}"
  THREADS_STEP1="${THREADS_STEP1:-24}"
  PREPROCESS_CORES="${PREPROCESS_CORES:-48}"
  FILES_PER_CORE="${FILES_PER_CORE:-1}"
  MIN_BASE_QUALITY="${MIN_BASE_QUALITY:-30}"
  MITO_GENOME="${MITO_GENOME:-rCRS}"
  BARCODE_TAG="${BARCODE_TAG:-BC}"
  REDEEM_QUICK="${REDEEM_QUICK:-0}"
  CONDA_SH="${CONDA_SH:-}"
  CONDA_ENV="${CONDA_ENV:-redeemv}"

  [ -n "${OUT_ROOT:-}" ] || die "OUT_ROOT is not set in $env_file"
  [ -n "${GENOME_PREFIX:-}" ] || die "GENOME_PREFIX is not set in $env_file"
}

activate_runtime_env() {
  if [ -n "$CONDA_SH" ]; then
    require_file "$CONDA_SH"
    # shellcheck disable=SC1090
    source "$CONDA_SH"
  fi
  if command -v conda >/dev/null 2>&1; then
    conda activate "$CONDA_ENV"
  else
    die "conda is not available; activate the REDEEM-V environment before submitting jobs or set CONDA_SH"
  fi
}

load_sample_row() {
  local samples_tsv="$1"
  local sample="$2"
  require_file "$samples_tsv"

  local row
  row=$(awk -F '\t' -v s="$sample" 'NR > 1 && $1 == s {print; found=1; exit} END {if (!found) exit 1}' "$samples_tsv") \
    || die "Sample '$sample' not found in $samples_tsv"

  IFS=$'\t' read -r SAMPLE R1 R2_BARCODE R3 CELLRANGER_BARCODES <<< "$row"
  [ "$SAMPLE" = "$sample" ] || die "Loaded wrong sample row for $sample"
  require_file "$R1"
  require_file "$R2_BARCODE"
  require_file "$R3"
  require_file "$CELLRANGER_BARCODES"
}

sample_output_dir() {
  local sample="$1"
  echo "${OUT_ROOT}/${sample}"
}

consensus_array_tasks() {
  echo $((PREPROCESS_CORES * FILES_PER_CORE))
}
