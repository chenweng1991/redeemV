# REDEEM-V Quickstart

This quickstart uses the small example files bundled in `source/` and the current shared Slurm workflow layout.

## Files

```text
examples/quickstart/
  Analysis_redeemV/
    project.env
    samples.tsv
```

The sample sheet maps the bundled example files as:

| Column | Example file |
| --- | --- |
| `r1` | `source/Example.R1.fq.gz` |
| `r2_barcode` | `source/Example.i5.fq.gz` |
| `r3` | `source/Example.R2.fq.gz` |
| `cellranger_barcodes` | `source/barcodes.tsv.gz` |

`Example.i5.fq.gz` is the barcode read for this bundled tutorial dataset. Illumina `I1` FASTQs are not used by the current project workflow.

## Dry Run

From the repository root:

```bash
DRY_RUN=1 /lab/solexa_weissman/cweng/workflows/redeemV/scripts/submit_redeemV_project.sh \
  examples/quickstart/Analysis_redeemV/project.env \
  examples/quickstart/Analysis_redeemV/samples.tsv
```

The dry run validates that the FASTQs and barcode file exist and prints the Slurm commands without submitting jobs.

## Run

After the dry run succeeds, submit the quickstart sample with:

```bash
/lab/solexa_weissman/cweng/workflows/redeemV/scripts/submit_redeemV_project.sh \
  examples/quickstart/Analysis_redeemV/project.env \
  examples/quickstart/Analysis_redeemV/samples.tsv
```

Outputs are written under:

```text
examples/quickstart/Analysis_redeemV/quickstart/
```

Final consensus outputs are written under:

```text
examples/quickstart/Analysis_redeemV/quickstart/Out_mitoConsensus/final/
```
