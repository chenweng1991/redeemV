# REDEEM-V

REDEEM-V is the mitochondrial variant-calling component of the REDEEM framework. It starts from multiome mitochondrial FASTQs plus Cell Ranger ARC cell barcodes and produces per-cell mitochondrial coverage and consensus genotype tables for downstream lineage tracing with [REDEEM-R](https://github.com/chenweng1991/REDEEM-R).

This repository contains the core mapping, barcode handling, QC, and consensus calling code. On the Weissman lab cluster, the recommended production entry point is the shared Slurm workflow:

```text
/lab/solexa_weissman/cweng/workflows/redeemV
```

The workflow wraps this package into one dependency chain per sample:

```text
map/QC -> preprocess -> consensus array -> finalize
```

## Installation

Clone the repository and create the conda environment:

```bash
git clone https://github.com/chenweng1991/REDEEM-V.git
cd REDEEM-V
conda env create -f environment.yml
conda activate redeemv
```

On the Weissman lab cluster, production runs usually use the shared conda environment configured in `project.env`. The `environment.yml` file documents the package requirements for reproducible development and local testing.

## Quickstart

A small workflow-style example is available in [examples/quickstart](examples/quickstart). From the repository root, validate the bundled example data with:

```bash
DRY_RUN=1 /lab/solexa_weissman/cweng/workflows/redeemV/scripts/submit_redeemV_project.sh \
  examples/quickstart/Analysis_redeemV/project.env \
  examples/quickstart/Analysis_redeemV/samples.tsv
```

The dry run checks input files and prints the Slurm commands without submitting jobs.

## Recommended Workflow

See [docs/slurm-workflow.md](docs/slurm-workflow.md) for the full project setup and submission checklist, and [docs/io.md](docs/io.md) for detailed input and output definitions.

For each project, create:

```text
PROJECT/
  Analysis_redeemV/
    project.env
    samples.tsv
```

Start from the maintained templates:

```bash
mkdir -p PROJECT/Analysis_redeemV
cp /lab/solexa_weissman/cweng/workflows/redeemV/config/example.project.env \
  PROJECT/Analysis_redeemV/project.env
cp /lab/solexa_weissman/cweng/workflows/redeemV/templates/sample_sheet.tsv \
  PROJECT/Analysis_redeemV/samples.tsv
```

Edit `project.env` for the project root, output root, genome index, Slurm resources, and package path. Edit `samples.tsv` with one row per sample:

```text
sample	r1	r2_barcode	r3	cellranger_barcodes
```

FASTQ assignment for the current workflow:

| `samples.tsv` column | REDEEM-V option | Meaning |
| --- | --- | --- |
| `r1` | `MultiomeATAC_mito.sh -1` | Read 1 FASTQ |
| `r2_barcode` | `MultiomeATAC_mito.sh -i` | Barcode read FASTQ |
| `r3` | `MultiomeATAC_mito.sh -2` | Read 2 FASTQ |
| `cellranger_barcodes` | `RNAbc2ATAC.R` input | `filtered_feature_bc_matrix/barcodes.tsv.gz` |

Illumina `I1` FASTQs are not used by this REDEEM-V workflow.

Always dry-run first:

```bash
DRY_RUN=1 /lab/solexa_weissman/cweng/workflows/redeemV/scripts/submit_redeemV_project.sh \
  PROJECT/Analysis_redeemV/project.env \
  PROJECT/Analysis_redeemV/samples.tsv
```

Submit all samples after validation:

```bash
/lab/solexa_weissman/cweng/workflows/redeemV/scripts/submit_redeemV_project.sh \
  PROJECT/Analysis_redeemV/project.env \
  PROJECT/Analysis_redeemV/samples.tsv
```

Submit selected samples by appending sample names:

```bash
/lab/solexa_weissman/cweng/workflows/redeemV/scripts/submit_redeemV_project.sh \
  PROJECT/Analysis_redeemV/project.env \
  PROJECT/Analysis_redeemV/samples.tsv \
  SAMPLE1 SAMPLE2
```

Outputs are written under `${OUT_ROOT}/${sample}`. Final consensus outputs are in:

```text
${OUT_ROOT}/${sample}/Out_mitoConsensus/final/
```

Key files:

- `QualifiedTotalCts`
- `RawGenotypes.Total.StrandBalance`
- `RawGenotypes.VerySensitive.StrandBalance`
- `RawGenotypes.Sensitive.StrandBalance`
- `RawGenotypes.Specific.StrandBalance`

## Core Package Entry Points

- `MyMultiome/MultiomeATAC_mito.sh`: trim, attach barcodes, map with Bowtie2, extract unique mitochondrial alignments, and generate QC.
- `MyMultiome/Helpers/RNAbc2ATAC.R`: convert Cell Ranger RNA barcodes to ATAC barcodes for multiome mitochondrial reads.
- `mitoConsensus/Preprocess.py`: split mitochondrial BAMs by barcode chunks for parallel consensus calling.
- `mitoConsensus/mitoConsensus.py`: call mitochondrial consensus variants for one barcode chunk.
- `mitoConsensus/Finalize.sh`: concatenate chunk outputs, add coverage depth, remove strand-biased calls, and write final tables.

## Dependencies

See [environment.yml](environment.yml) for the reproducible development environment. Runtime tools include:

- Python with `click`, `numpy`, `pandas`, `pysam`, and `progress`
- R with `ggplot2`, `dplyr`, `gridExtra`, `plyr`, and `labeling`
- `bowtie2`, `cutadapt`, `samtools`, and `bedtools`

The shared Slurm workflow activates the configured conda environment from `project.env`.

## Manual Tutorial

[`Tutorial_20221025.md`](Tutorial_20221025.md) is kept as a historical manual walkthrough and small example. For new multi-sample cluster projects, use the shared Slurm workflow above.
