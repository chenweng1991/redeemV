# Shared Slurm Workflow

This document describes the current maintained way to run REDEEM-V for multi-sample projects on the Weissman lab cluster.

## Locations

```text
Workflow root: /lab/solexa_weissman/cweng/workflows/redeemV
Submitter:     /lab/solexa_weissman/cweng/workflows/redeemV/scripts/submit_redeemV_project.sh
Package:       /lab/solexa_weissman/cweng/Packages/REDEEM-V
```

Project-specific files, logs, and outputs should stay in the project directory. Do not copy and edit the shared workflow scripts for each project.

## Project Layout

```text
PROJECT/
  Analysis_redeemV/
    project.env
    samples.tsv
    SAMPLE_A/
    SAMPLE_B/
```

Create the config files from the shared templates:

```bash
mkdir -p PROJECT/Analysis_redeemV
cp /lab/solexa_weissman/cweng/workflows/redeemV/config/example.project.env \
  PROJECT/Analysis_redeemV/project.env
cp /lab/solexa_weissman/cweng/workflows/redeemV/templates/sample_sheet.tsv \
  PROJECT/Analysis_redeemV/samples.tsv
```

## `project.env`

Required settings:

| Variable | Purpose |
| --- | --- |
| `PROJECT_DIR` | Project root. |
| `OUT_ROOT` | Output root, usually `${PROJECT_DIR}/Analysis_redeemV`. |
| `REDEEM_V` | Path to this REDEEM-V package. |
| `GENOME_PREFIX` | Bowtie2 genome index prefix for the mitochondrial-mask genome. |

Common tuning settings:

| Variable | Purpose |
| --- | --- |
| `THREADS_STEP1` | Threads for mapping and QC. |
| `PREPROCESS_CORES` | Parallel workers used by `Preprocess.py`. |
| `FILES_PER_CORE` | Number of barcode chunk files per preprocessing core. |
| `MIN_BASE_QUALITY` | Minimum base quality for consensus calling. |
| `MITO_GENOME` | Mitochondrial reference, usually `rCRS`. |
| `BARCODE_TAG` | BAM barcode tag, usually `BC`. |
| `REDEEM_QUICK` | Set to `1` to skip the mapping QC section. |

The step 3 Slurm array size is `PREPROCESS_CORES * FILES_PER_CORE`.

## `samples.tsv`

The file is tab-separated and must contain this header:

```text
sample	r1	r2_barcode	r3	cellranger_barcodes
```

Column meanings:

| Column | Required path |
| --- | --- |
| `sample` | Sample name used for output folder and Slurm job names. |
| `r1` | Read 1 FASTQ. Passed to `MultiomeATAC_mito.sh -1`. |
| `r2_barcode` | Barcode read FASTQ. Passed to `MultiomeATAC_mito.sh -i`. |
| `r3` | Read 2 FASTQ. Passed to `MultiomeATAC_mito.sh -2`. |
| `cellranger_barcodes` | Matching `filtered_feature_bc_matrix/barcodes.tsv.gz`. |

Illumina `I1` FASTQs are not used by this workflow.

## Validation Checklist

Before submitting jobs, confirm that:

- every FASTQ path in `samples.tsv` exists and is readable
- every `cellranger_barcodes` path exists and points to the matching sample
- `R1`, barcode read, and `R3` are assigned to the correct columns
- sample names are unique
- `OUT_ROOT` has enough storage for BAMs and consensus outputs
- `GENOME_PREFIX` points to a valid Bowtie2 index prefix

## Dry Run

Always dry-run before submission:

```bash
DRY_RUN=1 /lab/solexa_weissman/cweng/workflows/redeemV/scripts/submit_redeemV_project.sh \
  PROJECT/Analysis_redeemV/project.env \
  PROJECT/Analysis_redeemV/samples.tsv
```

The dry run prints the `sbatch` commands and validates required input files without submitting jobs.

## Submit

Submit all samples:

```bash
/lab/solexa_weissman/cweng/workflows/redeemV/scripts/submit_redeemV_project.sh \
  PROJECT/Analysis_redeemV/project.env \
  PROJECT/Analysis_redeemV/samples.tsv
```

Submit selected samples:

```bash
/lab/solexa_weissman/cweng/workflows/redeemV/scripts/submit_redeemV_project.sh \
  PROJECT/Analysis_redeemV/project.env \
  PROJECT/Analysis_redeemV/samples.tsv \
  SAMPLE_A SAMPLE_B
```

For each sample, the submitter schedules:

1. `redeem_step1_map.sh`: trim, barcode, map, extract mitochondrial BAM, and QC.
2. `redeem_step2_preprocess.sh`: convert Cell Ranger RNA barcodes to ATAC barcodes and split the mitochondrial BAM.
3. `redeem_step3_array_consensus.sh`: run `mitoConsensus.py` over barcode chunks as a Slurm array.
4. `redeem_step4_finalize.sh`: run `Finalize.sh` to concatenate, annotate, filter, and summarize final outputs.

Dependencies are per sample. If one sample fails, the other samples continue independently.

## Outputs

Sample outputs are written under:

```text
${OUT_ROOT}/${sample}/
```

Final consensus outputs are written under:

```text
${OUT_ROOT}/${sample}/Out_mitoConsensus/final/
```

Main downstream inputs for REDEEM-R:

```text
QualifiedTotalCts
RawGenotypes.Total.StrandBalance
RawGenotypes.VerySensitive.StrandBalance
RawGenotypes.Sensitive.StrandBalance
RawGenotypes.Specific.StrandBalance
```
