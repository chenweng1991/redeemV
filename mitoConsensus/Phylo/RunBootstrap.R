.libPaths(c('/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal',"/home/cweng/R/x86_64-pc-linux-gnu-library/4.1-focal"))
# .libPaths(c("/home/cweng/R/x86_64-pc-linux-gnu-library/4.1-focal","/nfs/apps/lib/R/4.1-focal/site-library.2021q2","/opt/R/4.1.0/lib/R/library","/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal"))
library(Rcpp)
library(Matrix)
library(ape)
library(phangorn) #,lib='/home/cweng/R/x86_64-pc-linux-gnu-library/4.1-focal'
# library(Signac,lib="/lab/solexa_weissman/cweng/Packages/R/x86_64-pc-linux-gnu-library/4.1-focal")
args = commandArgs(trailingOnly=TRUE)
#print(c("Available distance methods:",c("Jaccard","Dice","Simpson","Kulczynski2","Ochiai","Hamming")))
if(length(args)<5){
  stop("Please enter: matrix_file BootStrap(T/F) distance()   al(nj/upgma) out(output prefix)")
}
matrix_file=args[1]
BootStrap=args[2]  ## T(Do bootstrap) or F(Reference tree)
d=args[3] # Choose from jaccard,
al=args[4]  ## choose from nj, upgma, and others
out=args[5]

print(paste("Matix_file:",matrix_file,"\nBootStrap:",BootStrap,"\nd:",d,"\nal:",al,"\nout:",out))

## Define distance function
BinaryDist<-function(M,method="jaccard"){
print("This function compute pairwise distance(row-row) for binary matrix, input sparse matrix(Each row is cell, each column is variant)")
print("Available method:")
print(c("Jaccard","Dice","3WJaccard","Simpson","Kulczynski2","Ochiai","Hamming"))
Total<-rowSums(M) # Compute total variant number for each cell
a<-M %*% t(M)   ## Compute the overlaped variants across any two cells
b<-Total-a  ## Compute the variant only for the give row but not for the given column
c<-t(b)  ## Compute the variant only for the give column but not for the given row
if(method=="Jaccard"){
disimilarity<-1-a/(a+b+c)
distance<-as.dist(disimilarity)
}else if(method=="Dice"){
disimilarity<-1-2*a/(2*a+b+c)
distance<-as.dist(disimilarity)
}else if(method=="Simpson"){
bcmin<-pmin(b,c)
disimilarity<-1-a/(a+bcmin)
distance<-as.dist(disimilarity)
}else if(method=="Kulczynski2"){
pr1<-a/(a+b)
pr2<-a/(a+c)
disimilarity<-1-(pr1+pr2)/2
distance<-as.dist(disimilarity)
}else if (method=="Ochiai"){
pr1<-a/(a+b)
pr2<-a/(a+c)
disimilarity<-1-sqrt(pr1*pr2)
distance<-as.dist(disimilarity)
}else if(method=="Hamming"){
disimilarity<-b+c
distance<-as.dist(disimilarity)
}else if(method=="3WJaccard"){
disimilarity<-1-3*a/(3*a+b+c)
distance<-as.dist(disimilarity)
}
return(distance)
}


## If necessary, use the in-house writeTree,  borrow from http://blog.phytools.org/search?q=write
writeTree<-function(tree){
  tree<-reorder.phylo(tree,"cladewise")
  n<-length(tree$tip)
  string<-vector(); string[1]<-"("; j<-2
  for(i in 1:nrow(tree$edge)){
    if(tree$edge[i,2]<=n){
      string[j]<-tree$tip.label[tree$edge[i,2]]; j<-j+1
      if(!is.null(tree$edge.length)){
        string[j]<-paste(c(":",round(tree$edge.length[i],10)), collapse="")
        j<-j+1
      }
      v<-which(tree$edge[,1]==tree$edge[i,1]); k<-i
      while(length(v)>0&&k==v[length(v)]){
        string[j]<-")"; j<-j+1
        w<-which(tree$edge[,2]==tree$edge[k,1])
        if(!is.null(tree$edge.length)){
          string[j]<-paste(c(":",round(tree$edge.length[w],10)), collapse="")
          j<-j+1
        }
        v<-which(tree$edge[,1]==tree$edge[w,1]); k<-w
      }
      string[j]<-","; j<-j+1
    } else if(tree$edge[i,2]>=n){
      string[j]<-"("; j<-j+1
    }
  }
  if(is.null(tree$edge.length)) string<-c(string[1:(length(string)-1)], ";")
  else string<-c(string[1:(length(string)-2)],";")
  string<-paste(string,collapse="")
  return(string)
}



# matrix_file<-"CD34_BMMC_Cts.Mtx.bi.csv"
Cts.Mtx.bi<-Matrix(as.matrix(read.csv(matrix_file,row.names=1)))
## If This is for Bootstrap
if (BootStrap){
print("Bootstrapping...")
Cts.Mtx.bi<-Cts.Mtx.bi[,sample(1:ncol(Cts.Mtx.bi),replace = T)]
}else{
print("Compute for Ref tree")
}



## Compute distance
dist<-BinaryDist(Cts.Mtx.bi,method=d)
if(any(is.na(dist))){
dist[is.na(dist)]<-1
}

## Reconstruct
print("Reconstructing...")
if (al=="nj"){
tree<-nj(dist)
}else if(al=="upgma"){
tree<-upgma(dist)
}


res=try(write.tree(tree,file=paste(out,"nw",sep=".")), silent=TRUE)
if(class(res)=="try-error") {
  print("There is an C stackq error, use in-house writeTree instead")
  string<-writeTree(tree)
  conn <- file(paste(out,"nw",sep="."), "w")
  writeLines(string,conn, sep = "\n")
  close(conn)
}
