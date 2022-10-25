.libPaths(c("/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor4Donor9/renv/library/R-4.1/x86_64-pc-linux-gnu",.libPaths()))
library(scMitoTracing)
library(ggplot2)
library(dplyr)
library(BuenColors)
library(gridExtra)
args = commandArgs(trailingOnly=TRUE)
WD<-paste(args[1],"/final",sep="")
#WD<-paste(getwd(),"/DN4_mitoV/final",sep="")
#args[1] ## mitoV/final, eg WD<-"/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor4Donor9/Donor4/DN4_BMMC/MTenrichCombine/mitoV/final"



## Processing
depth<-DepthSummary(WD,CellSubset = NA,cellSubSetName = NA)
VariantsGTSummary<-CW_mgatk.read(WD,Processed =F)
sink("README_VariantFiltering.txt")
Variants.feature.lst<-Vfilter_v3(InputSummary=VariantsGTSummary,depth=depth)
cat("##-------##\n")
cat("##--Read in processed files by below--##\n")
cat('result_depth<-readRDS("',WD,'/result_depth.rds")\n',sep="")
cat('result_VariantsGTSummary<-readRDS("',WD,'/result_VariantsGTSummary.rds")\n',sep="")
cat('result_Variants.feature<-readRDS("',WD,'/result_Variants.feature.lst.rds")\n',sep="")
sink()

## Plotting- depth plot
png(paste(WD,"/results_depth.png",sep=""),width = 3000, height = 1600, res=300)
plot_depth(depth$Total,"Total")
dev.off()

## Plotting- mutation signature plot
png(paste(WD,"/results_mutSig.png",sep=""),width = 3000, height = 1600, res=300)
ps<-list()
for(name in names(Variants.feature.lst)){
p<-MutationProfile.bulk(Variants.feature.lst[[name]]$Variants)+ggtitle(name)+theme(title =element_text(size=20))
ps<-c(ps,list(p))
}
grid.arrange(grobs=ps)

pdf("results_VariantMetrics.pdf",width=14)
if(dim(Variants.feature.lst[[1]])[1]!=0){
plot_variant(VariantsGTSummary,Variants.feature.lst,depth=depth,cat=c("Total"),p4xlim = 30)
}else{
    ggplot()
}

saveRDS(depth,paste(WD,"/result_depth.rds",sep=""))
saveRDS(VariantsGTSummary,paste(WD,"/result_VariantsGTSummary.rds",sep=""))
saveRDS(Variants.feature.lst,paste(WD,"/result_Variants.feature.lst.rds",sep=""))
