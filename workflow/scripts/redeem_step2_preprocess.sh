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

if [ ! -f "${SAMPLE}.uniqmapped.mito.bam.bai" ]; then
  samtools index "${SAMPLE}.uniqmapped.mito.bam"
fi

python "$MITOCONSENSUS/Preprocess.py" \
  -i "${SAMPLE}.uniqmapped.mito.bam" \
  -c "$PREPROCESS_CORES" \
  -f "$FILES_PER_CORE" \
  -b "$ATAC_BARCODES" \
  -o ./Out_mitoConsensus/ \
  -g "$MITO_GENOME" \
  -bt "$BARCODE_TAG" \
  -sd "$MITOCONSENSUS"
