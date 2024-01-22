# redeemV
## Introduction
ReDeeM: single-cell **Re**gulatory multi-omics with **Dee**p **M**itochondrial mutation profiling. ReDeeM is a single-cell multiomics platform featuring ultra-sensitive mitochondrial DNA (mtDNA) variant calling and joint RNA+ATAC profiling. ReDeeM enables fine-scale lineage tracing at single cell level, allowing for subclonal and phylogenetic analyses, with simultaneous integrative analyses of cell-state and gene regulatory circuits.</br> 

The analytical pipelines for ReDeeM analysis includes two parts:
- [redeemV](https://github.com/chenweng1991/redeemV) is set of in-house Bash pipeline and python scripts for mapping and deep mitochondrial mutation calling. (**This page**, Input from rawdata) 
- [redeemR](https://github.com/chenweng1991/redeemR) is an in-house R package for downstream lineage tracing and single cell integrative analysis. (Input from redeemV)

**redeemV** is a streamlined pipeline taking advantage of endogenous unique molecular identifier (eUMI) for consensus-based error correction in single-cell mitochondrial DNA mutation detection. 
![Github fig variant calling strategy](https://github.com/chenweng1991/redeemV/assets/43254272/e20c55b3-056a-4e1f-bacf-73c9b4a5ae63)
## Installation and usage
redeemV includes a set of ready-to-use Bash pipeline and Python scripts
```
git clone https://github.com/chenweng1991/redeemV.git
```
Please check the [tutorial](https://github.com/chenweng1991/REDEEM-V/blob/master/Tutorial_20221025.md)
(A small set of example fastq files are included)

## Additional documentation
- [ReDeeM data structure](https://github.com/chenweng1991/redeemV/wiki/Organize-ReDeeM-full-data)
- [Run Cellranger-arc](https://github.com/chenweng1991/redeemV/wiki/Run-cellranger%E2%80%90arc)
- [Cell hashing processing](https://github.com/chenweng1991/redeemV/wiki/Cell-Hashing-Demultiplexing)

## Citation
Please check out our study of human hematopoiesis using ReDeeM [Deciphering cell states and genealogies of human hematopoiesis](https://doi.org/10.1038/s41586-024-07066-z)

## Contact
If you have any questions or suggestions, please feel free to contact us. Feedbacks are very welcome! (Chen Weng, cweng@broadinstitute.org)






