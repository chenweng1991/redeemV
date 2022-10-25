#!/bin/bash
RawBedSortTag=$1

DownSampleFrac=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99)
MyMultiome=/lab/solexa_weissman/cweng/Packages/MyMultiome

mkdir DownsampleTest

echo "Step1, DownSampling the Rawbed"
for f in ${DownSampleFrac[@]}
do
echo $f
cat $RawBedSortTag | awk -v myf=$f 'BEGIN {srand()} !/^$/ { if (rand() <= myf) print $0}' > DownsampleTest/tmp"_"$f.RawBed.Sort.Tag &
done


wait
cd DownsampleTest
echo "Step2, Deduplicate"
for F in `ls tmp*`
do
  python3 $MyMultiome/DeduplicateRawBed.10X.py $F 0 &
done

wait
echo "Step4 Extract Mito monoclonal tsv"
for F in `ls tmp*0.cut.tsv`
do
grep chrM $F > ${F/tsv/mito.tsv} &
done

wait
echo "Step5 Summarize"
for F in `ls tmp*tsv`
do
cat $F | python3 $MyMultiome/Summarize.TagDedup.10X.py > ${F/tsv/summary} &
done
