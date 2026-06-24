# Inputs and Outputs

This document defines the files REDEEM-V expects and the files it produces in the current shared Slurm workflow.

## Required Inputs

### Per-sample FASTQs

Each row in `samples.tsv` describes one sample:

```text
sample	r1	r2_barcode	r3	cellranger_barcodes
```

| Column | Passed to | Description |
| --- | --- | --- |
| `sample` | Slurm job names and output folder | Unique sample identifier. |
| `r1` | `MultiomeATAC_mito.sh -1` | Read 1 FASTQ. |
| `r2_barcode` | `MultiomeATAC_mito.sh -i` | Barcode read FASTQ used to attach cell barcodes to read names. |
| `r3` | `MultiomeATAC_mito.sh -2` | Read 2 FASTQ. |
| `cellranger_barcodes` | `RNAbc2ATAC.R` | Matching Cell Ranger ARC `filtered_feature_bc_matrix/barcodes.tsv.gz`. |

Illumina `I1` FASTQs are not used by the current REDEEM-V workflow.

### Cell Ranger Barcodes

`cellranger_barcodes` must point to the matching sample's filtered Cell Ranger ARC barcodes:

```text
filtered_feature_bc_matrix/barcodes.tsv.gz
```

These are RNA barcodes. The workflow converts them to ATAC barcodes with `MyMultiome/Helpers/RNAbc2ATAC.R` before mitochondrial read filtering.

### Genome Index

`GENOME_PREFIX` in `project.env` must point to a Bowtie2 genome index prefix for the mitochondrial-mask genome, for example:

```text
/lab/solexa_weissman/cweng/Genomes/GRCH38/GRCH38_Bowtie2_MitoMask/hg38.mitoMask
```

The prefix is passed to `bowtie2 -x` by `MyMultiome/MultiomeATAC_mito.sh`.

## Project Configuration

`project.env` controls paths, runtime settings, and Slurm resources.

Required variables:

| Variable | Description |
| --- | --- |
| `PROJECT_DIR` | Project root. |
| `OUT_ROOT` | Output root, usually `${PROJECT_DIR}/Analysis_redeemV`. |
| `REDEEM_V` | Path to this REDEEM-V package. |
| `GENOME_PREFIX` | Bowtie2 mitochondrial-mask genome index prefix. |

Common analysis variables:

| Variable | Default | Description |
| --- | --- | --- |
| `CUT` | `0` | Unique fragment cutoff used during mapping QC deduplication. |
| `THREADS_STEP1` | `24` | Threads for trimming, mapping, BAM processing, and QC. |
| `PREPROCESS_CORES` | `48` | Parallel workers used by `Preprocess.py`. |
| `FILES_PER_CORE` | `1` | Barcode chunk files produced per preprocessing core. |
| `MIN_BASE_QUALITY` | `30` | Minimum base quality used by `mitoConsensus.py`. |
| `MITO_GENOME` | `rCRS` | Mitochondrial genome configuration or FASTA. |
| `BARCODE_TAG` | `BC` | BAM tag containing cell barcodes. |
| `REDEEM_QUICK` | `0` | Set to `1` to skip mapping QC plots after mitochondrial BAM extraction. |

The consensus Slurm array size is:

```text
PREPROCESS_CORES * FILES_PER_CORE
```

## Intermediate Outputs

For each sample, outputs are written under:

```text
${OUT_ROOT}/${sample}/
```

Important intermediate files include:

| File or directory | Producer | Description |
| --- | --- | --- |
| `${sample}_atac.barcodes.tsv` | `RNAbc2ATAC.R` | ATAC-space whitelist converted from Cell Ranger RNA barcodes. |
| `${sample}.bam` | `MultiomeATAC_mito.sh` | Sorted mapped BAM. |
| `${sample}.tagged.bam` | `AddBC2BAM.py` | BAM with barcode tags added. |
| `${sample}.uniqmapped.bam` | `samtools view -bf 2 -q30` | Properly paired, uniquely mapped reads. |
| `${sample}.uniqmapped.mito.bam` | `samtools view chrM` | Unique mitochondrial alignments used for consensus calling. |
| `${sample}.QCplot.png` | `MultiATAC_mito.QC_v2.R` | Mapping and library-complexity QC plot, unless `REDEEM_QUICK=1`. |
| `Out_mitoConsensus/temp/barcoded_bams/` | `Preprocess.py` | Barcode-chunk BAMs for array consensus calling. |
| `Out_mitoConsensus/temp/barcode_files/` | `Preprocess.py` | Barcode chunk lists matching the BAM chunks. |
| `Out_mitoConsensus/temp/sparse_matrices2.0/` | `mitoConsensus.py` | Per-chunk genotype and coverage tables before finalization. |

## Final Outputs

Final consensus files are written under:

```text
${OUT_ROOT}/${sample}/Out_mitoConsensus/final/
```

Main files:

| File | Description |
| --- | --- |
| `QualifiedTotalCts` | Per-cell, per-position mitochondrial coverage table. |
| `RawGenotypes.Total` | Least stringent raw consensus variant calls. |
| `RawGenotypes.VerySensitive` | Less stringent raw consensus variant calls. |
| `RawGenotypes.Sensitive` | Stringent raw consensus variant calls. |
| `RawGenotypes.Specific` | Most stringent raw consensus variant calls. |
| `RawGenotypes.Total.StrandBalance` | Strand-bias-filtered total calls. |
| `RawGenotypes.VerySensitive.StrandBalance` | Strand-bias-filtered very-sensitive calls. |
| `RawGenotypes.Sensitive.StrandBalance` | Strand-bias-filtered sensitive calls. |
| `RawGenotypes.Specific.StrandBalance` | Strand-bias-filtered specific calls. |
| `TotalRawBamRows.Mito` | Count of mitochondrial BAM rows included across barcode chunks. |

The primary downstream inputs for REDEEM-R are:

```text
QualifiedTotalCts
RawGenotypes.Total.StrandBalance
RawGenotypes.VerySensitive.StrandBalance
RawGenotypes.Sensitive.StrandBalance
RawGenotypes.Specific.StrandBalance
```

## Output Table Schemas

### `QualifiedTotalCts`

This table records mitochondrial coverage per cell and position. The expected columns are:

| Column | Description |
| --- | --- |
| cell barcode | Cell barcode. |
| mitochondrial coordinate | Position on the mitochondrial genome. |
| total unique fragments | Least stringent coverage. |
| less-stringent unique fragments | Very-sensitive coverage. |
| stringent unique fragments | Sensitive coverage. |
| very-stringent unique fragments | Specific coverage. |

### `RawGenotypes.*`

Each row represents one molecule with a candidate mitochondrial variant.

| Column | Description |
| --- | --- |
| `MoleculeID` | Cell barcode plus molecule start/end identifier. |
| `CellBC` | Cell barcode. |
| `Pos` | Mitochondrial coordinate. |
| `Variant` | Variant label. |
| `V` | Variant base. |
| `Ref` | Reference base. |
| `FamSize` | Consensus family size. |
| `V-counts` | Reads supporting the variant. |
| `CSS` | Consensus score. |
| `DB_Cts` | Double-covered copies. |
| `SG_Cts` | Single-covered copies. |
| `Is+` | Variant observed on plus strand. |
| `Is-` | Variant observed on minus strand. |
| `TotalDepth` | Total unique-fragment depth at this cell and position. |
