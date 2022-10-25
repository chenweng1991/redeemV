# .libPaths(c('/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal',.libPaths()))
# library(Rcpp,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(dplyr,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(Matrix,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(Matrix.utils,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(ggplot2,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(MASS,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(gridExtra,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(Matrix.utils,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(BuenColors,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(EZsinglecell,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(reshape2,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(ape,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(treeio,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(ggtree,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(ggnewscale,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(ggtreeExtra,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(phangorn,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
# library(gsubfn,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
l<-'/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal'
.libPaths(c(.libPaths(),'/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal'))
library(Rcpp,lib=l)
library(fastmatch,lib=l)
library(TreeTools)
library(dplyr)
library(Matrix)
library(Matrix.utils)
library(ggplot2)
library(MASS)
library(gridExtra)
library(Matrix.utils)
library(BuenColors)
library(EZsinglecell)
library(reshape2)
library(ape)
library(treeio)
library(ggtree)
library(ggnewscale)
library(ggtreeExtra)
library(ggExtra)
library(phangorn)
library(gsubfn)
##Prepare Mito sequence context
mitoref<-read.table("/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/Run210706L2/COMBINE_CALLV/CW_mgatk_test/final/chrM_refAllele.txt")
mitobasebias<-table(toupper(mitoref$V2)) %>% as.data.frame()
# Important! Make Context dictionary
# library("insect",lib="/home/cweng/R/x86_64-pc-linux-gnu-library/4.1-focal")
library(stringi)

rc<-function(input){
o<-stri_reverse(gsubfn(".", list("A" = "T", "T" = "A","C"="G", "G"="C"), input))
return(o)    
}

Contexts<-c()
for(i in 2:(nrow(mitoref)-1)){
    Context<-paste(mitoref$V2[i-1],mitoref$V2[i],mitoref$V2[i+1],sep="")
    Contexts<-c(Contexts,toupper(Context))
}
ContextsDic<-c("GGA",Contexts,"TGG")
names(ContextsDic)<-as.character(1:nrow(mitoref))
head(ContextsDic)
idx<-which(substr(ContextsDic,2,2) %in% c("G","A"))
ContextsDic[idx]<-stri_reverse(rc(ContextsDic[idx]))
#-------------------------------------------------------------------------


#Helper functions---------------------------------------------------------


## A function to summarize the depth (Total that passed Q30)
DepthSummary<-function(path,Processed=T){
setwd(path)
if(Processed){
    depth<-readRDS("depth.RDS")
}else{
    QualifiedTotalCts<-read.table("QualifiedTotalCts")
    Pos.MeanCov<-QualifiedTotalCts %>% group_by(V2) %>% dplyr::summarise(meanCov=mean(V3)) 
    Cell.MeanCov<-QualifiedTotalCts %>% group_by(V1) %>% dplyr::summarise(meanCov=mean(V3))
    depth_Total<-list(Pos.MeanCov,Cell.MeanCov)
    
    Pos.MeanCov<-QualifiedTotalCts %>% group_by(V2) %>% dplyr::summarise(meanCov=mean(V4)) 
    Cell.MeanCov<-QualifiedTotalCts %>% group_by(V1) %>% dplyr::summarise(meanCov=mean(V4))
    depth_VerySensitive<-list(Pos.MeanCov,Cell.MeanCov)
    
    Pos.MeanCov<-QualifiedTotalCts %>% group_by(V2) %>% dplyr::summarise(meanCov=mean(V5)) 
    Cell.MeanCov<-QualifiedTotalCts %>% group_by(V1) %>% dplyr::summarise(meanCov=mean(V5))
    depth_Sensitive<-list(Pos.MeanCov,Cell.MeanCov)
    
    Pos.MeanCov<-QualifiedTotalCts %>% group_by(V2) %>% dplyr::summarise(meanCov=mean(V6)) 
    Cell.MeanCov<-QualifiedTotalCts %>% group_by(V1) %>% dplyr::summarise(meanCov=mean(V6))
    depth_Specific<-list(Pos.MeanCov,Cell.MeanCov)

    depth<-list(Total=depth_Total,VerySensitive=depth_VerySensitive,Sensitive=depth_Sensitive,Specific=depth_Specific)
    saveRDS(depth,"depth.RDS")    
return(depth)
}
}




## Function to generate GTS summary
GTSummary<-function(RawGenotypes,filterN=T){ ## At this moment, the context with N is probably prone to error due to mapping, in the future should work on realignment
# Make Depth dictionary
Depth<-unique(RawGenotypes[,c("Cell","Pos","Depth")])
Depthdic<-Depth$Depth
names(Depthdic)<-paste(Depth$Cell, Depth$Pos,sep="")
# Summarise    
Genotypes.summary<-table(paste(RawGenotypes$Cell,RawGenotypes$Variants,sep="_")) %>% as.data.frame()
Genotypes.summary$Cell<-strsplit(as.character(Genotypes.summary$Var1),"_") %>% sapply(.,function(x){x[1]})
Genotypes.summary$Variants<-strsplit(as.character(Genotypes.summary$Var1),"_") %>% sapply(.,function(x){paste(x[2:4],collapse="_")})
Genotypes.summary$cellPos<-strsplit(as.character(Genotypes.summary$Var1),"_") %>% sapply(.,function(x){paste(x[1:2],collapse="")}) 
Genotypes.summary$depth<-Depthdic[Genotypes.summary$cellPos]
Genotypes.summary<-Genotypes.summary[,c("Var1","Cell","Variants","Freq","depth")] 
Genotypes.summary$Type<-strsplit(Genotypes.summary$Variants,"_") %>% sapply(.,function(x){paste(x[2],x[3],sep="_")})
Genotypes.summary$Context<-ContextsDic[strsplit(Genotypes.summary$Variants,"_") %>% sapply(.,function(x){x[1]})]
if(filterN){
    Genotypes.summary<-subset(Genotypes.summary,!grepl("N",Genotypes.summary$Context))
}
return(Genotypes.summary)
}


CW_mgatk.read<-function(path,Processed=F){
setwd(path)
if(Processed){
    VariantsGTSummary<-readRDS("VariantsGTSummary.RDS")
}else{
    RawGenotypes.Sensitive.StrandBalance<-read.table("RawGenotypes.Sensitive.StrandBalance")
    RawGenotypes.VerySensitive.StrandBalance<-read.table("RawGenotypes.VerySensitive.StrandBalance")
    RawGenotypes.Specific.StrandBalance<-read.table("RawGenotypes.Specific.StrandBalance")
    RawGenotypes.Total.StrandBalance<-read.table("RawGenotypes.Total.StrandBalance")
    GiveName<-c("UMI","Cell","Pos","Variants","Call","Ref","FamSize","GT_Cts","CSS","DB_Cts","SG_Cts","Plus","Minus","Depth")
    colnames(RawGenotypes.Sensitive.StrandBalance)<-GiveName
    colnames(RawGenotypes.VerySensitive.StrandBalance)<-GiveName
    colnames(RawGenotypes.Specific.StrandBalance)<-GiveName
    colnames(RawGenotypes.Total.StrandBalance)<-GiveName
    Specific.GTSummary<-GTSummary(RawGenotypes.Specific.StrandBalance)
    Sensitive.GTSummary<-GTSummary(RawGenotypes.Sensitive.StrandBalance)
    VerySensitive.GTSummary<-GTSummary(RawGenotypes.VerySensitive.StrandBalance)
    Total.GTSummary<-GTSummary(RawGenotypes.Total.StrandBalance)
    ##Calculate heteroplasmy
    Specific.GTSummary$hetero<-with(Specific.GTSummary,Freq/depth)
    Sensitive.GTSummary$hetero<-with(Sensitive.GTSummary,Freq/depth)
    VerySensitive.GTSummary$hetero<-with(VerySensitive.GTSummary,Freq/depth)
    Total.GTSummary$hetero<-with(Total.GTSummary,Freq/depth)
    VariantsGTSummary<-list(Total=Total.GTSummary,VerySensitive=VerySensitive.GTSummary,Sensitive=Sensitive.GTSummary,Specific=Specific.GTSummary)
    saveRDS(VariantsGTSummary,"VariantsGTSummary.RDS")
}
return(VariantsGTSummary)
}

plot_depth<-function(depth=DN1CD34_1.depth,name=""){
d1<-depth[[1]] 
d2<-depth[[2]] 
names(d1)<-c("pos","meanCov")
names(d2)<-c("cell","meanCov")
options(repr.plot.width=10, repr.plot.height=3)
p1<-ggplot(d1)+aes(pos,meanCov)+geom_point()+theme_bw()
p2<-ggplot(d1)+aes("cell",meanCov)+geom_violin()+geom_boxplot()+theme_bw()
grid.arrange(p1,p2,layout_matrix=rbind(c(1,1,1,1,1,1,2)),top=name)
}


Vfilter_v3<-function(InputSummary,depth,Rmvhomo=F,Min_Cells=2, Max_Count_perCell=2,QualifyCellCut=10){
    CV<-function(x){
    var(x)/mean(x)
    }
Names<-names(InputSummary)
feature.list<-list()    
for(i in Names){    
VariantFeature0<- InputSummary[[i]] %>% group_by(Variants) %>% dplyr::summarise(CellN=n(),PositiveMean=mean(hetero),maxcts=max(Freq),CV=CV(hetero),TotalVcount=sum(Freq))
VariantFeature0$pos<-strsplit(VariantFeature0$Variants,"_") %>% sapply(.,function(x){x[1]}) %>% as.numeric
VariantFeature0<-merge(VariantFeature0,depth[[i]][[1]],by.x="pos",by.y="V2")   ## This generate different meanCov for each threahold
VariantFeature0$TotalCov<-length(unique(InputSummary[[i]]$Cell))*VariantFeature0$meanCov
VariantFeature0$VAF<-VariantFeature0$TotalVcount/VariantFeature0$TotalCov

qualifiedCell<-subset(depth[["Total"]][[2]],meanCov>=QualifyCellCut)[,1,drop=T]  ## Filter Qualified cell based on total depth
InputSummary.qualified<-subset(InputSummary[[i]],Cell %in% qualifiedCell)
VariantFeature<- InputSummary.qualified %>% group_by(Variants) %>% dplyr::summarise(CellN=n(),PositiveMean=mean(hetero),maxcts=max(Freq),CV=CV(hetero),TotalVcount=sum(Freq))
print(paste(i,":\n",nrow(VariantFeature0),"variants to start"))
print(paste(nrow(VariantFeature),"variants after remove low quality cells"))

VariantFeature$CellNPCT<-VariantFeature$CellN/length(unique(InputSummary.qualified$Cell))
VariantFeature<-merge(VariantFeature[,c("Variants","CellN","PositiveMean","maxcts","CellNPCT")],VariantFeature0[,c("Variants","TotalVcount","TotalCov","VAF","CV")],by="Variants")
#HomoVariants<-subset(VariantFeature,VAF>0.9 & CV<0.01)$Variants
HomoVariants<-subset(VariantFeature,CellNPCT>0.75 & PositiveMean>0.75 & CV<0.01)$Variants
VariantFeature$HomoTag<-ifelse(VariantFeature$Variants %in% HomoVariants,"Homo","Hetero")        
if (Rmvhomo){
    VariantFeature<-subset(VariantFeature,!Variants %in% HomoVariants)
    print(paste(length(HomoVariants),"Homoplasmy variants to remove"))
    print(HomoVariants)
}else{
    print(paste("Tag Homoplasmy:",HomoVariants))
}
out<-subset(VariantFeature,CellN>=Min_Cells & maxcts>=Max_Count_perCell)    
print(paste("After filtering,",nrow(out), "Variants left"))
print("\n\n")
feature.list<-c(feature.list,list(out))
} 
    names(feature.list)<-Names
    return(feature.list) 
}   


BinaryDist<-function(M,method="jaccard"){
print("This function compute pairwise distance(row-row) for binary matrix, input sparse matrix(Each row is cell, each column is variant)")
print("Available method:")
print(c("Jaccard","Dice","3WJaccard","Simpson","Kulczynski2","Ochiai","Hamming"))  
Total<-rowSums(M) # Compute total variant number for each cell
a<-M %*% t(M)   ## Compute the overlaped variants across any two cells
b<-Total-a  ## Compute the variant only for the give row but not for the given column
c<-t(b)  ## Compute the variant only for the give column but not for the given row
if(method=="Jaccard"){
disimilarity<-1-a/(a+b+c)
distance<-as.dist(disimilarity)
}else if(method=="Dice"){
disimilarity<-1-2*a/(2*a+b+c)
distance<-as.dist(disimilarity)
}else if(method=="Simpson"){
bcmin<-pmin(b,c)
disimilarity<-1-a/(a+bcmin)
distance<-as.dist(disimilarity)    
}else if(method=="Kulczynski2"){
pr1<-a/(a+b)
pr2<-a/(a+c)  
disimilarity<-1-(pr1+pr2)/2
distance<-as.dist(disimilarity)
}else if (method=="Ochiai"){
pr1<-a/(a+b)
pr2<-a/(a+c)  
disimilarity<-1-sqrt(pr1*pr2)
distance<-as.dist(disimilarity)  
}else if(method=="Hamming"){
disimilarity<-b+c
distance<-as.dist(disimilarity)   
}else if(method=="3WJaccard"){
disimilarity<-1-3*a/(3*a+b+c)
distance<-as.dist(disimilarity)
}    
return(distance)
}


#########----------------------------------------------------------------------------------
## Mutation signature functions

## Make ref_all_long for bulk level assesment (from Caleb mgatk), Juts run it, will be used
library(data.table)
library(dplyr)
# Simple reverse complement function
reverse_complement <- function(s){
  chartr("ATGC","TACG",s)
}

# Process 3 digit signature based on letters
ref_all <- fread("/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/Run210706L2/COMBINE_CALLV/CW_mgatk_test/final/chrM_refAllele.txt")
colnames(ref_all) <- c("pos", "ref")
ref_all$ref <- toupper(ref_all$ref)
l <- as.character(ref_all$ref)

# Gs happen to be at the first and last position
ref_all$three <- paste0(c("G", l[-length(l)]), l, c(l[-1], "G"))

# Remove Ns
ref_all <- ref_all[!grepl("N", ref_all$three),]

# Make every possible mutation
ref_all_long <- rbind(ref_all,ref_all, ref_all,ref_all)
ref_all_long$alt <- rep(c("A", "C", "G", "T"), each = dim(ref_all)[1])
ref_all_long <- ref_all_long[ref_all_long$ref != ref_all_long$alt,]

# add some meta data
ref_all_long$variant <- paste0(as.character(ref_all_long$pos), ref_all_long$ref, ">", ref_all_long$alt)
ref_all_long$change <- paste0(ref_all_long$ref, ref_all_long$alt)
ref_all_long$change_rc <- reverse_complement(paste0(ref_all_long$ref, ref_all_long$alt))

# A/G rich strand is "heavy" -- https://en.wikipedia.org/wiki/Heavy_strand
table(ref_all$ref) # so the reference strand is light (more C/T)
ref_all_long$strand <- ifelse(ref_all_long$ref %in% c("C","T"), "L", "H")

# Change to C/T as ref allele
ref_all_long$rc3 <- reverse_complement(ref_all_long$three)
ref_all_long$three_plot <- ifelse(ref_all_long$strand == "L", ref_all_long$three, ref_all_long$rc3)
ref_all_long$group_change <- ifelse(ref_all_long$strand == "L", ref_all_long$change, ref_all_long$change_rc)


##Below Define the function to plot bulk level mutation signatures
MutationProfile.bulk<-function(cell_variants){  ## cell_variants look like c('93_A_G''103_G_A''146_T_C''150_C_T''152_T_C''182_C_T')
# Annotate with called variants
called_variants <- strsplit(cell_variants,"_") %>% sapply(.,function(x){paste(x[1],x[2],">",x[3],sep="")})
ref_all_long$called <- ref_all_long$variant %in% called_variants
# Compute changes in expected/observed
total <- dim(ref_all_long)[1]
total_called <- sum(ref_all_long$called)
prop_df <- ref_all_long %>% group_by(three_plot, group_change, strand) %>%
  dplyr::summarize(observed_prop_called = sum(called)/total_called, expected_prop = n()/total, n = n()) %>%
  mutate(fc_called = observed_prop_called/expected_prop)
prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)
# Visualize
p1 <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
  geom_bar(stat = "identity", position = "dodge") + pretty_plot(fontsize = 8) + L_border() + 
  theme(axis.title.x=element_blank(),
        axis.text.x =element_blank())+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y = "Substitution Rate")+
   facet_grid(.~group_change,scales = "free",space="free") 
return(p1)   
}



CountVperCell<-function(x,name,CellN){
s<-c(CellN-length(x),length(which(x==1)),length(which(x>=2 &x<=5)),length(which(x>=6 &x<=10)),length(which(x>10)))
names(s)<-c("0","1","2-5","6-10",">10")
paste(paste(s,"Cells have",names(s),"Variants"),collapse="\n")
}

plot_variant<-function(GTSummary,feature.list,depth,cat=c("Total","VerySensitive","Sensitive","Specific"),p4xlim=50,QualifyCellCut=10){
options(repr.plot.width=20, repr.plot.height=6)
require(ggExtra)
qualifiedCell<-subset(depth[["Total"]][[2]],meanCov>=QualifyCellCut)[,1,drop=T]  ## Filter Qualified cell based on total depth
for (c in cat){
#p1
p1<-ggplot(feature.list[[c]])+aes(log10(VAF),log10(CV),color=HomoTag)+geom_point(size=0.2)+theme_bw()+scale_color_brewer(palette = "Set1")+theme(axis.text=element_text(size=20,color="black"))
p1<-ggMarginal(p1, type = "histogram",)
#p2
p2<-subset(GTSummary[[c]],depth>20) %>% .$hetero %>% data.frame(HeteroPlasmy=.) %>% ggplot()+aes(HeteroPlasmy)+geom_histogram(binwidth = 0.05)+theme_bw()+theme(axis.text=element_text(size=20,color="black"))
#p3
QualifiedV<-subset(feature.list[[c]],maxcts>=2 & CellN>=2 & HomoTag=="Hetero")
p3title=paste(c,":",nrow(QualifiedV),"Qualified Hetero Variants\nMedian Cell # per V:",median(QualifiedV$CellN),"\nVariants # maxcts>=3:",length(which(QualifiedV$maxcts>=3)))
p3<-ggplot(feature.list[[c]])+aes(log2(CellN),log2(maxcts))+geom_jitter(color="grey80")+geom_point(data=subset(feature.list[[c]],maxcts>=2 & CellN>=2 & HomoTag=="Hetero"),color="black",size=1)+theme_classic()+ggtitle(p3title)
p3<-ggMarginal(p3, type = "histogram",)
#p4

CellVar.Sum<-subset(GTSummary[[c]],Variants %in% QualifiedV$Variants & Cell %in% qualifiedCell) %>% group_by(Cell) %>% dplyr::summarise(VN=n(), maxcts=max(Freq),mediancts=median(Freq)) 
p4title<-paste("Qualified Cell number:",length(qualifiedCell),"\nMedian V number is",median(CellVar.Sum$VN),"\n",CountVperCell(CellVar.Sum$VN,c,CellN=nrow(CellVar.Sum)))
p4<-ggplot(CellVar.Sum)+aes(VN)+geom_histogram(binwidth = 1,color="black",fill="white")+xlim(0,p4xlim)+ggtitle(p4title)+theme(axis.text=element_text(size=20))+geom_vline(xintercept = median(CellVar.Sum$VN),linetype=2)    
grid.arrange(p1,p2,p3,p4,ncol=4,top=c)
}
}

SeuratLSIClustering<-function(ob,res=0.3){
require(Seurat)
require(Signac)
Cell_Variant.seurat<-CreateSeuratObject(counts = t(as.matrix(ob)), assay = "mitoV")
# Cell_Variant.seurat <- FindVariableFeatures(Cell_Variant.seurat)
# Cell_Variant.seurat <- NormalizeData(Cell_Variant.seurat)
# Cell_Variant.seurat <- ScaleData(Cell_Variant.seurat)
VariableFeatures(Cell_Variant.seurat) <- row.names(Cell_Variant.seurat) #names(which(Matrix::rowSums(Cell_Variant.seurat) > 100))
#Cell_Variant.seurat <- RunPCA(Cell_Variant.seurat, npcs = 10)
Cell_Variant.seurat <- RunTFIDF(Cell_Variant.seurat, n = 50)
Cell_Variant.seurat<- FindTopFeatures(Cell_Variant.seurat, min.cutoff = 'q0')
Cell_Variant.seurat <- RunSVD(Cell_Variant.seurat, n = 50)
Cell_Variant.seurat <- RunUMAP(Cell_Variant.seurat, reduction = "lsi", dims = 1:20)
Cell_Variant.seurat <- FindNeighbors(Cell_Variant.seurat,reduction ="lsi"  ,dims = 1:20)
# Cell_Variant.seurat <- FindClusters(Cell_Variant.seurat, resolution = 0.05)
Cell_Variant.seurat <- FindClusters(Cell_Variant.seurat, resolution = res)
}

Translate_RNA2ATAC<-function(meta=bmmc.filtered@meta.data,celltepe="CellType",RNAclusterPost="-1"){
ATACWhite<-read.table("/lab/solexa_weissman/cweng/Genomes/10X/WhiteList_10X_Multiome.ATAC")
RNAWhite<-read.table("/lab/solexa_weissman/cweng/Genomes/10X/WhiteList_10X_Multiome.RNA")
Dic2<-ATACWhite$V1
names(Dic2)<-as.character(RNAWhite$V1)
meta$ATACName<-Dic2[gsub(RNAclusterPost,"",row.names(meta))]
meta<-meta[,c("ATACName",celltepe)]
return(meta)
}


distance_jaccard<-function(bi){
Binary.jaccard<-Signac::Jaccard(x = bi, y = bi)
Binary.jaccard.dissimilarity<-1-Binary.jaccard
Binary.jaccard.distance<-as.dist(Binary.jaccard.dissimilarity) 
return(Binary.jaccard.distance)    
} 


#-------------------Functions for Cell Type analysis

    
CellTypeRandomTest<-function(TreeCellType,CellTypeFocus="MEP",n=100){
p<-data.frame(y=sample(1:nrow(TreeCellType),nrow(subset(TreeCellType,CellType==CellTypeFocus)))) %>% ggplot()+aes(y)+geom_histogram(binwidth=25)+theme(axis.text = element_text(size=25))
randomData<-ggplot_build(p)$data[[1]][,c("x","y")] 
for (i in 1:n){
p<-data.frame(y=sample(1:nrow(TreeCellType),nrow(subset(TreeCellType,CellType==CellTypeFocus)))) %>% ggplot()+aes(y)+geom_histogram(binwidth=25)+theme(axis.text = element_text(size=25))
randomData<-merge(randomData,ggplot_build(p)$data[[1]][,c("x","y")],by="x") 
}
    return(randomData)
}
    
    
GenerateHistData<-function(treeplot,pData=p.CD4,Stat=CD4Stat,background=ExpectC["CD4"],k=25){
expendhist<-function(x,k){
    out<-data.frame(pos=seq(as.integer(x[4]),as.integer(x[5])-1,1),count=as.integer(x[2]),pvalues=as.double(x[19]))
    return(out)
}
ps<-ggplot_build(pData)$data[[1]][,c("x","count")] %>% merge(.,Stat,by="x",all.x = T) %>% apply(.,1,function(x){1-length(which(x[2] >x[3:length(x)]))/(length(x)-2)}) 

HistData<-cbind(ggplot_build(pData)$data[[1]],pvalues=ps) %>% apply(.,1,expendhist,k) %>% do.call(rbind,.) %>% subset(.,pos>0 & pos<length(unique(treeplot$data$label))) %>%
merge(.,subset(treeplot$data[,c("label","y")],!is.na(label)),by.x="pos",by.y="y")
HistData$Fold<-HistData$count/background
HistData$log10p<--log10(HistData$pvalues+1e-16)    
return(HistData)
}
    
    
    
#--------------------mitoTracing 
library(phytools)
library(ape)
library(phangorn)
library(treeio)
library(ggtree)
library(Signac)
library(Seurat)

Datatoplots<-setClass(
    "Datatoplots",
    slots=c(
    clustering="data.frame"
    )
)

DistObjects<-setClass(
    "DistObjects",
    slots=c(
    jaccard="dist",
    Dice="dist",
    jaccard3W="dist"    
    )
)


TREE<-setClass(
    "TREE",
    slots=c(
    phylo="phylo",
    treedata="treedata",
    records="character"
    )
)

mitoTracing<-setClass(
    "mitoTracing",
    slots=c(GTsummary.filtered="data.frame",
            CellMeta="data.frame",
            V.fitered.list="list",
            UniqueV="character",
            Cts.Mtx="dgCMatrix",
            Cts.Mtx.bi="dgCMatrix",
            para="character",
            Seurat="Seurat",
            DataToplotList="Datatoplots",
            DistObjects="DistObjects",
            TREE="TREE"
           )
)


setGeneric(name="Make_matrix", def=function(object)standardGeneric("Make_matrix"))
setGeneric(name="SeuratLSIClustering", def=function(object,...)standardGeneric("SeuratLSIClustering"))
setGeneric(name="AddDatatoplot_clustering", def=function(object,...)standardGeneric("AddDatatoplot_clustering"))
setGeneric(name="AddDist_jaccard", def=function(object,...)standardGeneric("AddDist_jaccard"))
setGeneric(name="AddDist", def=function(object,...) standardGeneric("AddDist"))     
setGeneric(name="Make_tree", def=function(object,d="jaccard", algorithm="upgma",onlyreturntree=F,...) standardGeneric("Make_tree"))         
setGeneric(name="AddTree", def=function(object,phylo,...) standardGeneric("AddTree"))     
  
setMethod(f="show",
          signature="mitoTracing",
          definition=function(object){
              print(object@para)
              print(paste("Total Cell number:",nrow(object@CellMeta)))
              print(table(object@CellMeta$Label))
              print(paste("Total Variant number:",length(object@UniqueV)))
              print(paste("Slot:",slotNames(object)))
          })

setMethod(f="Make_matrix",
          signature="mitoTracing",
          definition=function(object){
                require(dplyr)
                require(Matrix.utils)
                Cts.Mtx<-dMcast(object@GTsummary.filtered,Cell~Variants,value.var = "Freq")
                colnames(Cts.Mtx)<-strsplit(as.character(colnames(Cts.Mtx)),"_") %>% sapply(.,function(x){paste(x[1],x[2],x[3],sep="")})
                Cts.Mtx.bi<-Cts.Mtx
                Cts.Mtx.bi[Cts.Mtx.bi>=1]<-1
                object@Cts.Mtx.bi<-Cts.Mtx.bi
                object@Cts.Mtx<-Cts.Mtx
#                 validObject(object)
                return(object)
})


setMethod(f="Make_tree",
          signature="mitoTracing",
          definition=function(object,d,algorithm,onlyreturntree=F){
          dist<-slot(object@DistObjects,d)
          if(algorithm=="nj"){
          phylo<-nj(dist)
          }else if (algorithm=="upgma"){
          phylo<-upgma(dist)
          }
          treedata<-as.treedata(phylo)
          TREEobject<-new("TREE",phylo=phylo,treedata=treedata,records=paste(d,algorithm,sep="-"))       
          if(onlyreturntree){    
          return(TREEobject)
          }else{
          object@TREE<-TREEobject
          return(object)
          }             
})           
 
           
           
setMethod(f="SeuratLSIClustering",
          signature="mitoTracing",
          definition=function(object,binary=T,res=0.3){
          require(Signac)
          if(binary){    
              Cell_Variant.seurat<-CreateSeuratObject(counts = t(as.matrix(object@Cts.Mtx.bi)), assay = "mitoV")
          }else{
              Cell_Variant.seurat<-CreateSeuratObject(counts = t(as.matrix(object@Cts.Mtx)), assay = "mitoV")
          }
          VariableFeatures(Cell_Variant.seurat) <- row.names(Cell_Variant.seurat) #names(which(Matrix::rowSums(Cell_Variant.seurat) > 100))
          Cell_Variant.seurat <- RunTFIDF(Cell_Variant.seurat, n = 50)
          Cell_Variant.seurat<- FindTopFeatures(Cell_Variant.seurat, min.cutoff = 'q0')
          Cell_Variant.seurat <- RunSVD(Cell_Variant.seurat, n = 50)
          Cell_Variant.seurat <- RunUMAP(Cell_Variant.seurat, reduction = "lsi", dims = 1:20)
          Cell_Variant.seurat <- FindNeighbors(Cell_Variant.seurat,reduction ="lsi"  ,dims = 1:20)
          Cell_Variant.seurat <- FindClusters(Cell_Variant.seurat, resolution = 0.3)
          object@Seurat<-Cell_Variant.seurat
          return(object)
})


setMethod(f="AddDatatoplot_clustering",
          signature="mitoTracing",
          definition=function(object){
          row.names(object@CellMeta)<-object@CellMeta$Cell
          datatoplot<-Tomerge_v2(object@Seurat@meta.data,object@Seurat@reductions$umap@cell.embeddings) %>% Tomerge_v2(.,object@Seurat@reductions$lsi@cell.embeddings[,1:6]) %>% Tomerge_v2(.,object@CellMeta)
          object@DataToplotList<-Datatoplots()
          object@DataToplotList@clustering<-datatoplot
          return(object)    
})   
           

# Legacy function
# setMethod(f="AddDist_jaccard",
#           signature="mitoTracing",
#           definition=function(object){
#           d<-distance_jaccard(object@Cts.Mtx.bi)
#           object@DistObjects<-new("DistObjects",jaccard=d)
#           return(object)    
#           })

setMethod(f="AddDist",
          signature="mitoTracing",
          definition=function(object){
          d.Jaccard<-BinaryDist(object@Cts.Mtx.bi,method="Jaccard")
          d.Dice<-BinaryDist(object@Cts.Mtx.bi,method="Dice")
          d.3WJaccard<-BinaryDist(object@Cts.Mtx.bi,method="3WJaccard")    
          object@DistObjects<-new("DistObjects",jaccard=d.Jaccard, Dice=d.Dice,jaccard3W=d.3WJaccard)
          return(object)    
          })

setMethod(f="AddTree",
          signature="mitoTracing",
          definition=function(object,phylo,record=""){
          TREEobject<-new("TREE",phylo=phylo,treedata=as.treedata(phylo),records=record)
          object@TREE<-TREEobject
          return(object)    
          })
           
# Advanced version to create mitoTracing object 
Create_mitoTracing<-function(GTsummary_list,depth_list,feature.list_list,meta_list,labels,thr="VerySensitive",qualifiedCellCut=10,OnlyHetero=T,VAFcut=1,Cellcut=2,maxctscut=2,sampleCell=F){
CellMeta.all<-c()
GTsummary.all<-c()
V.union<-c()
V.fitered.list<-list()
len<-length(GTsummary_list)
for(i in 1:len){
CellMeta<-subset(depth_list[[i]]$Total[[2]],meanCov>=qualifiedCellCut)
names(CellMeta)[1]<-"Cell"
CellMeta<-merge(CellMeta,meta_list[[i]],by.x="Cell",by.y="ATACName")
CellMeta$Cell<-paste(CellMeta$Cell,i,sep="_")
CellMeta$Label<-labels[i]
GTsummary<-GTsummary_list[[i]][[thr]]
GTsummary$Cell<-paste(GTsummary$Cell,i,sep="_")
V.filtered<-subset(feature.list_list[[i]][[thr]],VAF<=VAFcut & CellN>=Cellcut & maxcts>=maxctscut) 
if(OnlyHetero){
    V.filtered<-subset(V.filtered,HomoTag=="Hetero") 
}
CellMeta.all<-rbind(CellMeta.all,CellMeta)
GTsummary.all<-rbind(GTsummary.all,GTsummary)
V.union<-c(V.union,as.character(V.filtered$Variants))
V.fitered.list<-c(V.fitered.list,list(V.filtered))
}
V.union.unique<-unique(V.union)
names(V.fitered.list)<-labels
GTsummary.all.filtered<-subset(GTsummary.all,Variants %in% V.union.unique & Cell %in% CellMeta.all$Cell)
ob<-mitoTracing()
ob@GTsummary.filtered<-GTsummary.all.filtered
ob@CellMeta<-CellMeta.all
ob@V.fitered.list=V.fitered.list
ob@UniqueV<-V.union.unique
ob@para<-c(Threhold=thr,qualifiedCellCut=qualifiedCellCut,OnlyHetero=OnlyHetero,VAFcut=VAFcut,Cellcut=Cellcut,maxctscut=maxctscut,sampleCell=sampleCell)
return(ob)
}





