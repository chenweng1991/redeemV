#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
# shellcheck source=scripts/common.sh
source "$SCRIPT_DIR/common.sh"

[ -n "${PROJECT_ENV:-}" ] || die "PROJECT_ENV is not set"
[ -n "${SAMPLES_TSV:-}" ] || die "SAMPLES_TSV is not set"
[ -n "${SAMPLE:-}" ] || die "SAMPLE is not set"

source_project_env "$PROJECT_ENV"
activate_runtime_env
load_sample_row "$SAMPLES_TSV" "$SAMPLE"

OUT_DIR=$(sample_output_dir "$SAMPLE")
mkdir -p "$OUT_DIR/Log"
cd "$OUT_DIR"

ATAC_BARCODES="${SAMPLE}_atac.barcodes.tsv"
if [ ! -f "$ATAC_BARCODES" ]; then
  Rscript "$MYMULTIOME/Helpers/RNAbc2ATAC.R" "$REDEEM_V" "$CELLRANGER_BARCODES" "$ATAC_BARCODES"
fi

quick_args=()
if [ "$REDEEM_QUICK" = "1" ]; then
  quick_args=(-q)
fi

"$MYMULTIOME/MultiomeATAC_mito.sh" \
  -n "$SAMPLE" \
  -1 "$R1" \
  -2 "$R3" \
  -i "$R2_BARCODE" \
  -c "$CUT" \
  -t "$THREADS_STEP1" \
  -m "$MYMULTIOME" \
  -b "$GENOME_PREFIX" \
  -w "$ATAC_BARCODES" \
  "${quick_args[@]}"
