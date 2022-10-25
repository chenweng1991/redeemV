#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(gridExtra)
library(plyr)
args = commandArgs(trailingOnly=TRUE)
name=args[1]


ReadsCounts<-read.table(paste(name,"ReadsCounts",sep="."))
##Plot Reads processing summary
names(ReadsCounts)<-c("Name","N")
ReadsCounts$Name<-factor(c("Raw","RawMito","Mapped","MappedMito"),levels=c("Raw","RawMito","Mapped","MappedMito"))

p1<-ggplot(ReadsCounts)+aes(Name,N,fill=Name)+geom_bar(color="black",stat="identity")+ggtitle(paste("Reads processing:\nTotalmitoFraction: ",100*round(ReadsCounts$N[2]/ReadsCounts$N[1],2),"%\nPropermappedMitoFraction: ",100*round(ReadsCounts$N[4]/ReadsCounts$N[3],2),"%",sep=""))+theme_classic()+theme(axis.text=element_text(size=15,color="black"))



##Print out the plots
pdf(paste(name,"QCplot.pdf",sep="."),width=5,height=4)
print(p1)
dev.off()
