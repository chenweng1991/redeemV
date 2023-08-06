#!/bin/bash
echo "This is not excutable, but to find a good exampl, go to"
echo "/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/CoverageAnalysis/DN1_BMMC2_multikit"

die()
## Note:  this should be run in a filder that has already run the MultiomeATAC_mito_Miseq.sh or
## Note this is only a template, the DownSampleFrac needs to be customized aacordingly
name=$1
DownSampleFrac=(0.1 0.2 0.3) ## As an example it can be like this  (0.1 0.2 0.3)
ReadBarcode=$4
Cut=$5 # Minimum uniq fragment per cell to be considered
MyMultiome=/lab/solexa_weissman/cweng/Packages/MyMultiome
lib=/lab/solexa_weissman/cweng/Packages/cxw486/scATAClib/



##Count the totalDepth for $name.uniqmapped.RawBed.Sort.Tag
Depth=`cat $name.uniqmapped.RawBed.Sort.Tag | wc -l`
for f in ${DownSampleFrac[@]}
do
echo $f
cat $name.uniqmapped.RawBed.Sort.Tag | awk -v myf=$f 'BEGIN {srand()} !/^$/ { if (rand() <= myf) print $0}' > tmp"_"$f.RawBed &
done


#Step7 deduplicate at single cell monoclonal tsv
wait
echo "Deduplicate at single cell monoclonal tsv"
for f in ${DownSampleFrac[@]}
do
echo $f
python3 $MyMultiome/DeduplicateRawBed.10X.py tmp"_"$f.RawBed 0 & #  This will generate $name.uniqmapped.fragment.tsv an $name.uniqmapped.fragment.1000cut.tsv
wait
cat tmp"_"$f.fragment.tsv | python3 $MyMultiome/SubsetTableByDictionary.py CellWhiteList >tmp"_"$f.fragment.InWL.tsv
done

##Step8 Extract Mito monoclonal tsv
wait
echo "Extract Mito monoclonal tsv"
for f in ${DownSampleFrac[@]}
do
echo $f
grep chrM tmp"_"$f.fragment.InWL.tsv > tmp"_"$f.fragment.InWL.mito.tsv &
done

wait
echo "Summarize"
for f in ${DownSampleFrac[@]}
do
echo $f
cat tmp"_"$f.fragment.InWL.tsv | python3 $MyMultiome/Summarize.TagDedup.10X.py > tmp"_"$f.summary &
cat tmp"_"$f.fragment.InWL.mito.tsv | python3 $MyMultiome/Summarize.TagDedup.10X.py > tmp"_"$f.mito.summary &
done
