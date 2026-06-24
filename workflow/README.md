# redeemV Workflow

Reusable Slurm wrapper for the REDEEM-V multiome mitochondrial variant workflow.

The shared workflow lives here:

```text
workflow
```

Each project should keep only its own config, sample sheet, logs, and outputs. Do
not edit the shared scripts for each project unless the workflow itself changes.

## Project Layout

Recommended per-project layout:

```text
PROJECT/
  Analysis_redeemV/
    project.env
    samples.tsv
    SAMPLE_1/
    SAMPLE_2/
```

## Sample Sheet

`samples.tsv` is tab-separated and must contain:

```text
sample  r1  r2_barcode  r3  cellranger_barcodes
```

For this workflow:

- `r1` maps to `MultiomeATAC_mito.sh -1`
- `r3` maps to `MultiomeATAC_mito.sh -2`
- `r2_barcode` maps to `MultiomeATAC_mito.sh -i`
- Illumina `I1` FASTQs are not used by this REDEEM-V step

## Config

Copy `config/example.project.env` into the project analysis folder and edit paths.

Key settings:

- `PROJECT_DIR`: project root
- `OUT_ROOT`: where sample output folders are written
- `REDEEM_V`: REDEEM-V package path
- `GENOME_PREFIX`: Bowtie2 mitochondrial-mask genome prefix
- `PREPROCESS_CORES` and `FILES_PER_CORE`: determine the number of consensus array tasks

## Submit

From any directory:

```bash
workflow/scripts/submit_redeemV_project.sh \
  /path/to/project.env \
  /path/to/samples.tsv
```

To submit selected samples:

```bash
workflow/scripts/submit_redeemV_project.sh \
  /path/to/project.env \
  /path/to/samples.tsv \
  NK_D6S NK_D6D
```

For each sample the submitter schedules:

1. map and QC ATAC mitochondrial reads
2. convert CellRanger RNA barcodes to ATAC barcodes and split mito BAMs
3. run `mitoConsensus.py` as a Slurm array
4. finalize genotype outputs

Outputs are written under `${OUT_ROOT}/${sample}`.
