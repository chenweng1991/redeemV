#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(gridExtra)
library(plyr)
args = commandArgs(trailingOnly=TRUE)
f=args[1]
name=args[2]
Cell.Reads<-read.table(f)


Cell.Reads.order<-Cell.Reads[order(Cell.Reads$V2,decreasing=T),]
Cell.Reads.order$cumFragRatio<-cumsum(Cell.Reads.order$V2)/sum(Cell.Reads.order$V2)
Cell.Reads.order$rk<-rank(-Cell.Reads.order$V2)
p<-ggplot(Cell.Reads.order[1:40000,])+aes(rk,cumFragRatio)+geom_point()+theme_bw()+geom_vline(xintercept=10000,linetype=2)+ggtitle(name)+ylim(0,1)

png(paste(name,"QCplot.png",sep="."),width=500,height=500,res=150)
print(p)
dev.off()
