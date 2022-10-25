#!/bin/bash
GexBam=$1

DownSampleFrac=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 0.99)
MyMultiome=/lab/solexa_weissman/cweng/Packages/MyMultiome


python3 $MyMultiome/Helpers/Gex2RawBed.py $GexBam > Gex.RawBed
mkdir DownsampleTest.Gex


echo "Step1, DownSampling the Rawbed"
for f in ${DownSampleFrac[@]}
do
echo $f
cat Gex.RawBed | awk -v myf=$f 'BEGIN {srand()} !/^$/ { if (rand() <= myf) print $0}' > DownsampleTest.Gex/tmp"_"$f.RawBed.Gex &
done


wait
cd DownsampleTest.Gex
echo "Step2, Count UMI"
for F in `ls tmp*RawBed.Gex`
do
  cut -f1,3 $F| sort -u | cut -f1 | uniq -c | awk '{print $2,$1}' OFS="\t" | sort -rgk 2 > ${F/RawBed.Gex/UniqUMI} &
done

wait

echo "Step3, Count Gene"
for F in `ls tmp*RawBed.Gex`
do
  cut -f1,2 $F| sort -u | cut -f1 | uniq -c | awk '{print $2,$1}' OFS="\t" |sort -rgk 2 > ${F/RawBed.Gex/UniqGene} &
done

echo "Step4, Count TotalReads"
for F in `ls tmp*RawBed.Gex`
do
  cut -f1 $F| sort | uniq -c | awk '{print $2,$1}' OFS="\t" |sort -rgk 2 > ${F/RawBed.Gex/TotalReads} &
done
