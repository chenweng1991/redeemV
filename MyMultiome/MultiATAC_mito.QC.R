#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(gridExtra)
library(plyr)
args = commandArgs(trailingOnly=TRUE)
name=args[1]
Cut=args[2]

ReadsCounts<-read.table(paste(name,"ReadsCounts",sep="."))
CellReads<-read.table(paste(name,"uniqmapped.fragment",Cut,"cut.summary",sep="."))
Cell.Mito.Reads<-read.table(paste(name,"uniqmapped.fragment",Cut,"cut.mito.summary",sep="."))
# FragmentsOnPeak<-read.table(paste(name,"uniqmapped.Peak_bc_sparse_mtx",sep="."))
tsvOnCell<-read.table(paste(name,"uniqmapped.fragment",Cut,"cut.tsv",sep="."))
tsvOnCell.mito<-read.table(paste(name,"uniqmapped.fragment",Cut,"cut.mito.tsv",sep="."))
##Plot Reads processing summary
data<-ReadsCounts
data$V1<-factor(data$V1,levels=c("TotalReads","TotalMitoReads","UniqProperpairReads","MitoReadsPair"))
datatoplot<-data.frame(Name=c("Raw","TotalMito","Mapped","MappedMito","Read_Cell","Frag_Cell"),N=c(ReadsCounts$V2[1],ReadsCounts$V2[2],ReadsCounts$V2[3],ReadsCounts$V2[4],sum(CellReads$V3),sum(CellReads$V2)))
datatoplot$Name<-factor(datatoplot$Name,levels=c("Raw","TotalMito","Mapped","MappedMito","Read_Cell","Frag_Cell"))
p1<-subset(datatoplot,Name %in%c("Raw","TotalMito","Mapped","MappedMito")) %>% ggplot(.)+aes(Name,N,fill=Name)+geom_bar(color="black",stat="identity")+ggtitle(paste("Reads processing:\nTotalmitoFraction: ",100*round(data$V2[2]/data$V2[1],2),"%\nPropermappedMitoFraction: ",100*round(data$V2[4]/data$V2[3],2),"%",sep=""))+theme_classic()+theme(axis.text=element_text(size=6,color="black"))


OnCellRatio<-paste(100*round(datatoplot[5,2]/datatoplot[3,2],2),"%",sep="")
Complexity<-paste(100*round(median(CellReads$V2/CellReads$V3),2),"%",sep="")
ComplexityMito<-paste(100*round(median(Cell.Mito.Reads$V2/Cell.Mito.Reads$V3),2),"%",sep="")
tt<-paste("OnCellRatio:",OnCellRatio,"\nComplexity:",Complexity,"\nComplexityMito:",ComplexityMito)

p2<-subset(datatoplot,Name %in%c("Raw","Mapped","Read_Cell","Frag_Cell")) %>% ggplot(.)+aes(Name,N,fill=Name)+geom_bar(color="black",stat="identity")+ggtitle(tt)+theme_classic()+theme(axis.text=element_text(size=6,color="black"))
# ##Plot single cell level QC
# FragmentsOnPeak.summary<-ddply(FragmentsOnPeak,.(V2),summarise,FragOnPK=sum(V3))
# FragmentsOnPeak.summary<-merge(CellReads,FragmentsOnPeak.summary,by.x="V1",by.y="V2")
# FragmentsOnPeak.summary$OnPeakRatio<-FragmentsOnPeak.summary$FragOnPK/FragmentsOnPeak.summary$V2
# FragmentsOnPeak.summary$rank<-rank(-FragmentsOnPeak.summary$V2)
# p3<-ggplot(FragmentsOnPeak.summary)+aes(rank,V2)+geom_point()+theme_classic()+ggtitle("Knee plot, rank vs uniq fragment #")
# p4<-ggplot(FragmentsOnPeak.summary)+aes("UniqFrag",V2)+geom_violin()+geom_boxplot(width=0.2)+theme_classic()+ggtitle("uniq fragment violin")
# p5<-ggplot(FragmentsOnPeak.summary)+aes(V2,OnPeakRatio)+geom_point()+theme_classic()+ggtitle(paste("uniq fragment # vs on peak ratio\nUsing cutoff  uniq frag>1000 \nOnPeakRatio>0.15\nIdentify Cell #: ", nrow(subset(FragmentsOnPeak.summary,V2>1000 & OnPeakRatio>0.15))))+geom_hline(yintercept=0.15,linetype=2)+geom_vline(xintercept=1000,linetype=2)



##Plot Mitofragments per cell and
MitoFraction<-merge(CellReads[,1:2],Cell.Mito.Reads[,1:2],by="V1")
names(MitoFraction)<-c("Cell","Fragments","MitoFragments")
MitoFraction$mitoFraction<-with(MitoFraction,MitoFragments/Fragments)

p4<-ggplot(MitoFraction)+aes("UniqFrag.all",log10(Fragments))+geom_violin()+geom_boxplot(width=0.2)+theme_classic()
p5<-ggplot(MitoFraction)+aes("UniqFrag.mito",log10(MitoFragments))+geom_violin()+geom_boxplot(width=0.2)+theme_classic()
p4_5<-ggplot(MitoFraction)+aes(log10(Fragments),log10(MitoFragments))+geom_point(size=0.2)
p6<-ggplot(MitoFraction)+aes("MitoFraction",mitoFraction)+geom_boxplot()+theme_classic()

##plot fragment size distribution
tsvOnCell$size<-tsvOnCell$V3-tsvOnCell$V2
tsvOnCell.mito$size<-tsvOnCell.mito$V3-tsvOnCell.mito$V2
p7<-ggplot(tsvOnCell)+aes(size)+geom_density()+theme_classic()+ggtitle("ALL.insert.distribution")
p7.mito<-ggplot(tsvOnCell.mito)+aes(size)+geom_density()+theme_classic()+ggtitle("Mito.insert.distribution")


##complexities
Data4Pie<-as.data.frame(table(tsvOnCell$V5))
if(nrow(Data4Pie)>=8){
  Data4Pie<-rbind(Data4Pie[1:8,],data.frame(Var1=">=8",Freq=sum(Data4Pie[9:nrow(Data4Pie),2])))
}
ttPie<-paste(paste("1-copy:",round((Data4Pie$Freq/sum(Data4Pie$Freq))[1],2)),"\n",
paste("2-copy:",round((Data4Pie$Freq/sum(Data4Pie$Freq))[2],2)),"\n",
paste("3-copy:",round((Data4Pie$Freq/sum(Data4Pie$Freq))[3],2)))


Data4Pie.mito<-as.data.frame(table(tsvOnCell.mito$V5))
if(nrow(Data4Pie.mito)>=8){
  Data4Pie.mito<-rbind(Data4Pie.mito[1:8,],data.frame(Var1=">=8",Freq=sum(Data4Pie.mito[9:nrow(Data4Pie.mito),2])))
}

ttPie.mito<-paste(paste("1-copy:",round((Data4Pie.mito$Freq/sum(Data4Pie.mito$Freq))[1],2)),"\n",
paste("2-copy:",round((Data4Pie.mito$Freq/sum(Data4Pie.mito$Freq))[2],2)),"\n",
paste("3-copy:",round((Data4Pie.mito$Freq/sum(Data4Pie.mito$Freq))[3],2)))

p8<-ggplot(Data4Pie)+aes("",Freq,fill=Var1)+geom_bar(stat="identity",width=1)+coord_polar("y", start=0)+scale_fill_brewer(palette="Set1")+theme_classic()+ggtitle(paste("Complexity\n",ttPie))
p9<-ggplot(Data4Pie.mito)+aes("",Freq,fill=Var1)+geom_bar(stat="identity",width=1)+coord_polar("y", start=0)+scale_fill_brewer(palette="Set1")+theme_classic()+ggtitle(paste("ComplexityMito\n",ttPie.mito))



##Print out the plots
pdf(paste(name,"QCplot.pdf",sep="."),width=9)
grid.arrange(p1,p2,ncol=2,nrow=2)
grid.arrange(p4,p5,p4_5,p6,ncol=2,nrow=2)
grid.arrange(p7,p7.mito,ncol=2,nrow=2)
grid.arrange(p8,p9,ncol=2,nrow=2)
dev.off()
