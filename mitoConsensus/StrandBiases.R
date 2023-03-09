#.libPaths(c('/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal',.libPaths()))
library(plyr)
library(dplyr)
library(ggplot2)
library(qvalue)
args = commandArgs(trailingOnly=TRUE)
WD<-args[1]  # "RawGenotypes.Total"
SimBoundry<-as.integer(args[2])
qlimit<-as.double(args[3])
RawGenotypesFile<-paste(WD,"/final/RawGenotypes.Total",sep="")
Blackout<-paste(WD,"/final/StrandBiaseBlackList",sep="")
img<-paste(WD,"/final/StrandBiase.png",sep="")
# QualifiedDepth<-args[2]
log2fold=1

if(is.na(SimBoundry)){
  SimBoundry<-200
}
print(RawGenotypesFile)
print(Blackout)
print(SimBoundry)

RawGenotypes<-read.table(RawGenotypesFile)
print("RawGenotypes In")
x<-ddply(RawGenotypes,.(V4),summarise,plus=sum(V12))
print("x In")
y<-ddply(RawGenotypes,.(V4),summarise,minus=sum(V13))
print("y In")
datatoplot<-data.frame(Variant=x[,1],plus=x[,2],minus=y[,2])
datatoplot$N<-datatoplot$plus+datatoplot$minus
print("datatoplot In")
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
write.table(BlackListTable[,1],paste(Blackout),quote=F,row.names=F,col.names=F)
p<-p+ggtitle(paste(nrow(BlackListTable),"suspected strandbiased variants,\n when qcut=",qlimit,"fold>",2^log2fold))

png(img,width = 1800, height = 1600, res=300)
print(p)
dev.off()
