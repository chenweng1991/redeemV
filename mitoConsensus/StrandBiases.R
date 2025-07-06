#.libPaths(c('/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal',.libPaths()))
library(plyr)
library(dplyr)
library(ggplot2)
library(qvalue)
args = commandArgs(trailingOnly=TRUE)
WD<-args[1]  # "RawGenotypes.Total"
SimBoundry<-as.integer(args[2])
qlimit<-as.double(args[3])

# Define file paths for both regular and combined
RawGenotypesFile<-paste(WD,"/final/RawGenotypes.Total",sep="")
RawGenotypesFile_combined<-paste(WD,"/final/RawGenotypes.Total.combined",sep="")

Blackout<-paste(WD,"/final/StrandBiaseBlackList",sep="")
Blackout_combined<-paste(WD,"/final/StrandBiaseBlackList.combined",sep="")

img<-paste(WD,"/final/StrandBiase.png",sep="")
img_combined<-paste(WD,"/final/StrandBiase.combined.png",sep="")

# QualifiedDepth<-args[2]
log2fold=1

if(is.na(SimBoundry)){
  SimBoundry<-200
}

print("Processing regular files:")
print(RawGenotypesFile)
print(Blackout)
print("Processing combined files:")
print(RawGenotypesFile_combined)
print(Blackout_combined)
print(SimBoundry)

# Function to process strand bias analysis
process_strand_bias <- function(genotypes_file, blacklist_file, image_file, file_type) {
  print(paste("Processing", file_type, "file..."))
  
  # Check if file exists
  if (!file.exists(genotypes_file)) {
    print(paste("Warning:", genotypes_file, "not found, skipping..."))
    return()
  }
  
  RawGenotypes <- read.table(genotypes_file, header = FALSE, skip = 1)
  print(paste("RawGenotypes", file_type, "In"))
  
  x<-ddply(RawGenotypes,.(V4),summarise,plus=sum(V12))
  print(paste("x", file_type, "In"))
  
  y<-ddply(RawGenotypes,.(V4),summarise,minus=sum(V13))
  print(paste("y", file_type, "In"))
  
  datatoplot<-data.frame(Variant=x[,1],plus=x[,2],minus=y[,2])
  datatoplot$N<-datatoplot$plus+datatoplot$minus
  print(paste("datatoplot", file_type, "In"))
  
  # sim<-matrix(seq)
  plus.sims<-c()
  minus.sims<-c()
  ps<-c()
  for (plus.sim in 0:200){
    for (minus.sim in 0:200){
      if((plus.sim+minus.sim)>0){
        #print(plus.sim)
        mod<-binom.test(plus.sim,plus.sim+minus.sim,0.5)
        plus.sims<-c(plus.sims,plus.sim)
        minus.sims<-c(minus.sims,minus.sim)
        ps<-c(ps,mod$p.value)
      }
    }
  }
  simdatatoplot<-data.frame(plus.sims,minus.sims,ps)
  
  p<-ggplot(data=simdatatoplot,aes(plus.sims,minus.sims,color=-log10(ps)))+geom_point()+xlim(0,200)+ylim(0,200)+scale_color_gradient(low="white",high="red",limits=c(2, 5), oob=scales::squish)+geom_point(data=datatoplot,aes(plus,minus),color="black",size=0.1)+theme_classic()
  
  ## Make strand biase blacklist
  pvalues<-as.matrix(datatoplot[,2:4]) %>% apply(.,1,function(x){binom.test(x[1],x[3],p=0.5)$p.value})
  datatoplot$pvalues<-pvalues
  qvalues<-qvalue(datatoplot$pvalues)$qvalues
  datatoplot$odds<-abs(log2(datatoplot$plus/datatoplot$minus))
  BlackListTable<-subset(datatoplot,qvalues<qlimit & odds>log2fold)
  #print(BlackListTable)
  write.table(BlackListTable[,1],blacklist_file,quote=F,row.names=F,col.names=F)
  p<-p+ggtitle(paste(nrow(BlackListTable),"suspected strandbiased variants",file_type,",\n when qcut=",qlimit,"fold>",2^log2fold))
  
  png(image_file,width = 1800, height = 1600, res=300)
  print(p)
  dev.off()
  
  print(paste("Completed processing", file_type, "- found", nrow(BlackListTable), "strand-biased variants"))
}

# Process regular files
process_strand_bias(RawGenotypesFile, Blackout, img, "regular")

# Process combined files
process_strand_bias(RawGenotypesFile_combined, Blackout_combined, img_combined, "combined")