# Simple code to translate cell ranger generated RNA barcode into ATAC barcode
# The input is the barcodes.tsv.gz from Cellranger  "barcodes.tsv.gz"
args = commandArgs(trailingOnly=TRUE)

# Make the dictionary
ATACWhite<-read.table("/lab/solexa_weissman/cweng/Genomes/10X/WhiteList_10X_Multiome.ATAC")
RNAWhite<-read.table("/lab/solexa_weissman/cweng/Genomes/10X/WhiteList_10X_Multiome.RNA")
Dic2<-ATACWhite$V1
names(Dic2)<-as.character(RNAWhite$V1)


# Read in the barcodes.tsv.gz
RNAbc<-read.table(args[1],col.names="RNABarcode")

# Translate
Translated<-data.frame(ATACBArcode=Dic2[gsub("-1","",RNAbc$RNABarcode)])
write.table(Translated,args[2],row.names=F,col.names=F, quote=F)
