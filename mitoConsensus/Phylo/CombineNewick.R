library(Rcpp)
library(Matrix)
library(ape)
library(phangorn)
args = commandArgs(trailingOnly=TRUE)
path=args[1]
N=args[2]
prefix=args[3]
print("Input: (1)bootstrap folder    (2)The number of trees to combine    (3) prefix(Combine.XXXX is suggested)")

all.nws<-list.files(path)[grep("Boot.",list.files(path))]
print(length(all.nws))

if(N>length(all.nws)){
  stop(paste("N should be less than the file numebr is in",path))
}else{
  print(paste("Choose",N,"trees out of",length(all.nws),"Generated trees"))
  chosen.nw<-sample(all.nws,N)
}


Combine.trees<-read.tree(paste(path,chosen.nw[1],sep=""))
for(i in 2:length(chosen.nw)){
  Combine.trees<-c(Combine.trees,read.tree(paste(path,chosen.nw[i],sep="")))
}
## Write out the combined tree
write.tree(Combine.trees,file=paste(path,prefix,"_",N,"trees.nw",sep=""))
print(paste("Output",paste(path,prefix,"_",N,"trees.nw",sep="")))
