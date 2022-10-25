l<-'/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal'
.libPaths(c(.libPaths(),'/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal'))
library(Rcpp,lib=l)
library(fastmatch,lib=l)
library(TreeTools)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)


#' Wrap Seurat RNA clustering 
#'
#' This function allows you to 
#' @param mtx sparse Matrix of class "dgCMatrix", each row is a gene, each column is a cell, 
#' @param exp The name of this sample/experiment
#' @export ob Standard Seurat object
#' @examples
#' bmmc.data=Read10X(data.dir = "/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_BMMC_1/CellRanger/Donor01_BMMC_1/outs/filtered_feature_bc_matrix")
#' docluster_GEM(mtx=bmmc.data$`Gene Expression`,exp="DN1_BMMC1")
docluster_GEM<-function(mtx=bmmc.data$`Gene Expression`,exp="DN1_BMMC1"){
require(Seurat)
ob <- CreateSeuratObject(counts = mtx, project = exp, min.cells = 3, min.features = 200)
ob <- NormalizeData(ob, normalization.method = "LogNormalize", scale.factor = 10000)
ob <- FindVariableFeatures(ob, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ob)
ob <- ScaleData(ob, features = all.genes)
ob <- RunPCA(ob, features = VariableFeatures(object = ob))
ob <- FindNeighbors(ob, dims = 1:10)
ob <- FindClusters(ob, resolution = 0.5)
ob <- RunUMAP(ob, dims = 1:10)
return(ob)
}

#' 
#'
#' This function allows you to express your love of cats.
#' @param love Do you love cats? Defaults to TRUE.
#' @keywords cats
#' @export
#' @examples
#' cat_function()
MultiWrapper<-function(path="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_CD34_1_Multiomekit/CellRanger/Donor01_CD34_1/outs"){
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
inputdata.10x <- Read10X_h5(paste(path,"/raw_feature_bc_matrix.h5",sep=""))
per_barcode_metrics<-read.csv(paste(path,"/per_barcode_metrics.csv",sep=""))
CellID<-subset(per_barcode_metrics,is_cell==1)$barcode
# Extract rna and atac counts
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
# Only use standard chromasome for atac counts
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
# Filter atac and rna counts
rna_counts.filtered<-rna_counts[,CellID]
atac_counts.filtered<-atac_counts[,CellID]
# Use RNA to create the default object
ob<-CreateSeuratObject(counts = rna_counts.filtered)    
# Create chrome_assay
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations)<-"UCSC" 
genome(annotations) <- "hg38"
frag.file <- paste(path,"/atac_fragments.tsv.gz",sep="")
chrom_assay <- CreateChromatinAssay(
   counts = atac_counts.filtered,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
)
# Add chrom_assay 
ob[["ATAC"]]<-chrom_assay
## Further filter the object
ob <- subset(
  x = ob,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 1e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000
)
# RNA analysis
DefaultAssay(ob) <- "RNA"
ob <- SCTransform(ob, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(ob) <- "ATAC"
ob <- RunTFIDF(ob)
ob <- FindTopFeatures(ob, min.cutoff = 'q0')
ob <- RunSVD(ob)
ob <- RunUMAP(ob, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
ob <- FindMultiModalNeighbors(ob, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
ob <- RunUMAP(ob, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
ob <- FindClusters(ob, graph.name = "wsnn", algorithm = 3, verbose = FALSE) 
return(list(seurat=ob,metric=per_barcode_metrics))
}    