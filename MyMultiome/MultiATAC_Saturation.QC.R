#!/usr/bin/env Rscript
library(ggplot2)
library(dplyr)
library(gridExtra)
library(plyr)
args = commandArgs(trailingOnly=TRUE)
name=args[1]
WhiteListFile<-ead.table(args[2]
CellWhiteList=
MitoSummary<-
for (f in c(0.002,0.005,0.01,0.0125,0.02,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)){
paste("tmp_",f,".summary",sep="")
}
