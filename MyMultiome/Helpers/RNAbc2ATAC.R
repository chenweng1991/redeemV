# Simple code to translate cell ranger generated RNA barcode into ATAC barcode
# The input is the barcodes.tsv.gz from Cellranger  "barcodes.tsv.gz"
args = commandArgs(trailingOnly=TRUE)

# Make the dictionary
REDEEM_V<-args[1]
ATACWhite<-read.table(paste(REDEEM_V,"/source/WhiteList_10X_Multiome.ATAC",sep=""))
RNAWhite<-read.table(paste(REDEEM_V,"/source/WhiteList_10X_Multiome.RNA",sep=""))
Dic2<-ATACWhite$V1
names(Dic2)<-as.character(RNAWhite$V1)


# Read in the barcodes.tsv.gz
RNAbc<-read.table(args[2],col.names="RNABarcode")

# Translate
Translated<-data.frame(ATACBArcode=Dic2[gsub("-1","",RNAbc$RNABarcode)])
write.table(Translated,args[3],row.names=F,col.names=F, quote=F)
