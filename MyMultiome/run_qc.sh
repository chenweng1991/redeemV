#!/bin/bash
source /broad/software/scripts/useuse
use .cellranger-arc-2.0.2
use Anaconda3
conda init --all
source activate redeemV

mitobam=$1
ncores=$2
Cut=0
redeemV=/broad/sankaranlab/tgao/redeemV
MyMultiome=$redeemV/MyMultiome

outdir=$(dirname $mitobam)
cd $outdir
name=$(basename "${mitobam%.uniqmapped.mito.bam}")

# stop if running into error
set -e

echo "Starting QC steps"

##Step7 Get raw bed file and Add cell barcode to the fragment bed files
if [ ! -s "$name.uniqmapped.RawBed.mito.Tag" ]; then
  echo "Running step 7 Get raw bed file and Add cell barcode to the fragment bed files..."
  samtools sort -@ $ncores -n $mitobam -o $name.uniqmapped.mito.namesorted.bam
  bedtools bamtobed -bedpe -i $name.uniqmapped.mito.namesorted.bam | awk -v OFS='\t' '{split($7,name,"|"); print name[2],$1,$2,$6,$7}' > $name.uniqmapped.RawBed.mito.Tag
else
  echo "$name.uniqmapped.RawBed.mito.Tag. Skip Step 7."
fi

##Step7.5 
if [ ! -s "$name.uniqmapped.RawBed.mito.Sort.Tag" ]; then
  echo "Running step 7.5"
  sort --parallel=$ncores -S 4G -T $outdir -k1,1 -k2,2 -k3,3n -k4,4n -k5,5 $name.uniqmapped.RawBed.mito.Tag > $name.uniqmapped.RawBed.mito.Sort.Tag
else
  echo "$name.uniqmapped.RawBed.mito.Sort.Tag. Skip Step 7."
fi


#Step8 deduplicate at single cell monoclonal tsv
if [ ! -s "$name.uniqmapped.fragment.$Cut.cut.mito.tsv" ]; then
  echo "Running step7 deduplicate at single cell monoclonal tsv..."
  python3 $MyMultiome/DeduplicateRawBed.10X.py $name.uniqmapped.RawBed.mito.Sort.Tag $name.uniqmapped.fragment.$Cut.cut.mito.tsv --cutoff $Cut 
else
  echo "$name.uniqmapped.fragment.$Cut.cut.mito.tsv eise. Skip Step 8."
fi

##Step9 Summarize
if [ ! -s "$name.uniqmapped.fragment.$Cut.cut.mito.summary" ]; then
  echo "Running step9 Summarize..."
  cat $name.uniqmapped.fragment.$Cut.cut.mito.tsv | python3 $MyMultiome/Summarize.TagDedup.10X.py > $name.uniqmapped.fragment.$Cut.cut.mito.summary
else
  echo "$name.uniqmapped.fragment.$Cut.cut.mito.summary exist. Skip step 9."
fi

##Step10 ReadsCount
if [ ! -s "$name.ReadsCounts" ]; then
  echo "Running step10 ReadsCount..."
  $MyMultiome/Counts.2.sh $name.tagged.bam $name.ReadsCounts $ncores
else
  echo "$name.ReadsCounts exist. Skip step 10."
fi

##Step 11 Plot QC
if [ ! -s "$name.QCplot.png" ]; then
  echo "Running step 11 Plot QC..."
  $MyMultiome/MultiATAC_mito.QC_v2.R $name $name.ReadsCounts $name.uniqmapped.fragment.$Cut.cut.mito.summary $name.uniqmapped.fragment.$Cut.cut.mito.tsv
else
  echo "$name.QCplot.png exist. Skip Step 11"
fi