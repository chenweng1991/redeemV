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
FragmentsOnPeak<-read.table(paste(name,"uniqmapped.Peak_bc_sparse_mtx",sep="."))
tsvOnCell<-read.table(paste(name,"uniqmapped.fragment",Cut,"cut.tsv",sep="."))
tsvOnCell.mito<-read.table(paste(name,"uniqmapped.fragment",Cut,"cut.mito.tsv",sep="."))
##Plot Reads processing summary
data<-ReadsCounts
data$V1<-factor(data$V1,levels=c("TotalReads","UniqProperpairReads","MitoReadsPair"))
datatoplot<-data.frame(Name=c("Raw","Uniq","UniqMito","Read_Cell","Frag_Cell"),N=c(ReadsCounts$V2[1],ReadsCounts$V2[2],ReadsCounts$V2[3],sum(CellReads$V3),sum(CellReads$V2)))
datatoplot$Name<-factor(datatoplot$Name,levels=c("Raw","Uniq","UniqMito","Read_Cell","Frag_Cell"))
p1<-subset(datatoplot,Name %in%c("Raw","Uniq","UniqMito")) %>% ggplot(.)+aes(Name,N,fill=Name)+geom_bar(color="black",stat="identity")+ggtitle("Reads processing")+theme_classic()+theme(axis.text=element_text(size=15,color="black",angle=45))
p2<-subset(datatoplot,Name %in%c("Raw","Uniq","Read_Cell","Frag_Cell")) %>% ggplot(.)+aes(Name,N,fill=Name)+geom_bar(color="black",stat="identity")+ggtitle("Reads processing")+theme_classic()+theme(axis.text=element_text(size=15,color="black",angle=45))
##Plot single cell level QC
FragmentsOnPeak.summary<-ddply(FragmentsOnPeak,.(V2),summarise,FragOnPK=sum(V3))
FragmentsOnPeak.summary<-merge(CellReads,FragmentsOnPeak.summary,by.x="V1",by.y="V2")
FragmentsOnPeak.summary$OnPeakRatio<-FragmentsOnPeak.summary$FragOnPK/FragmentsOnPeak.summary$V2
FragmentsOnPeak.summary$rank<-rank(-FragmentsOnPeak.summary$V2)
p3<-ggplot(FragmentsOnPeak.summary)+aes(rank,V2)+geom_point()+theme_classic()+ggtitle("Knee plot, rank vs uniq fragment #")
p4<-ggplot(FragmentsOnPeak.summary)+aes("UniqFrag",V2)+geom_violin()+geom_boxplot(width=0.2)+theme_classic()+ggtitle("uniq fragment violin")
p4.mito<-ggplot(Cell.Mito.Reads)+aes("UniqFrag.Mito",V2)+geom_violin()+geom_boxplot(width=0.2)+theme_classic()+ggtitle("uniq Mito fragment violin")
p5<-ggplot(FragmentsOnPeak.summary)+aes(V2,OnPeakRatio)+geom_point()+theme_classic()+ggtitle(paste("uniq fragment # vs on peak ratio\nUsing cutoff  uniq frag>1000 \nOnPeakRatio>0.15\nIdentify Cell #: ", nrow(subset(FragmentsOnPeak.summary,V2>1000 & OnPeakRatio>0.15))))+geom_hline(yintercept=0.15,linetype=2)+geom_vline(xintercept=1000,linetype=2)

##plot fragment size distribution
tsvOnCell$size<-tsvOnCell$V3-tsvOnCell$V2
tsvOnCell.mito$size<-tsvOnCell.mito$V3-tsvOnCell.mito$V2
p6<-ggplot(tsvOnCell)+aes(size)+geom_density()+theme_classic()+ggtitle("ALL.insert.distribution")
p6.mito<-ggplot(tsvOnCell.mito)+aes(size)+geom_density()+theme_classic()+ggtitle("Mito.insert.distribution")

##Plot Mitofragments per cell and
MitoFraction<-merge(CellReads[,1:2],Cell.Mito.Reads[,1:2],by="V1")
MitoFraction$mitoFraction<-with(MitoFraction,V2.y/V2.x)
p7<-ggplot(MitoFraction)+aes("MitoFraction",mitoFraction)+geom_boxplot()+theme_classic()

##complexities
p8<-as.data.frame(table(tsvOnCell$V5))%>% ggplot()+aes("",Freq,fill=Var1)+geom_bar(stat="identity",width=1)+coord_polar("y", start=0)+scale_fill_brewer(palette="Set1")+theme_classic()+ggtitle("Complexity on all uniq fragments")
p9<-as.data.frame(table(tsvOnCell.mito$V5))%>% ggplot()+aes("",Freq,fill=Var1)+geom_bar(stat="identity",width=1)+coord_polar("y", start=0)+scale_fill_brewer(palette="Set1")+theme_classic()+ggtitle("Complexity on mito uniq fragments")



##Print out the plots
pdf(paste(name,"QCplot.pdf",sep="."))
grid.arrange(p1,p2,ncol=2)
grid.arrange(p3,p4,p4.mito,p5,ncol=2,nrow=2)
grid.arrange(p7,p8,p9,ncol=2,nrow=2)
grid.arrange(p6,p6.mito,ncol=2,nrow=2)
dev.off()
