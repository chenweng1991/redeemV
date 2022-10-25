#!/bin/bash
## Note:  this version is for miseq or novaseq,  nextseq is slightly different
name=$1
Cut=$2 # Minimum uniq fragment per cell to be considered
MyMultiome=/lab/solexa_weissman/cweng/Packages/MyMultiome
lib=/lab/solexa_weissman/cweng/Packages/cxw486/scATAClib/


##Step11 ReadsCount
echo "Step11 ReadsCount"
$MyMultiome/Counts.2.sh $name

##Simple Plot QC
echo "Step12 Plot QC"
$MyMultiome/MultiATAC_mito.QC.R $name $Cut

##cleanup
rm -rf tmp
