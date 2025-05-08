#!/bin/bash
## Simplified script for mitochondrial DNA QC
## Based on original by Chen Weng
## This version starts from uniqmapped.mito.bam file
## Updated for use with Nextflow pipeline

Help()
{
  # Display Help
  echo "This script performs QC on mitochondrial DNA data starting from uniqmapped.mito.bam file"
  echo
  echo "MitoQC.sh -h for this page"
  echo "Syntax: MitoQC.sh -n -b -c -t -m"
  echo "Options"
  echo "-n name: The prefix of all analyzed files"
  echo "-b bam: Path to the input BAM file"
  echo "-c Cut: the cutoff the uniq fragment per cell"
  echo "-t CORE: The number of cores to use"
  echo "-m MyMultiome: The path to the folder of MyMultiome"
}

while getopts "hn:b:c:t:m:" option; do
  case $option in
    h) # display help
        Help
        exit;;
    n) # The prefix of all analyzed files
        name=$OPTARG;;
    b) # Path to the BAM file
        bam_file=$OPTARG;;
    c) # The cutoff the uniq fragment per cell
        Cut=$OPTARG;;
    t) # The number of cores to use
        CORE=$OPTARG;;
    m) # The path to the folder of MyMultiome
        MyMultiome=$OPTARG;;
   \?) # Invalid option
        echo "Error: Invalid option"
        exit;;
  esac
done

## Exit if any of the necessary input is empty
echo "Running MitoQC script - QC only version"
echo "Please run MitoQC.sh -h for help info"
set -u
: "$name$bam_file$Cut$CORE$MyMultiome"

# Check if input file exists
if [ ! -f "$bam_file" ]; then
  echo "Error: Input file $bam_file does not exist!"
  exit 1
else
  echo "Found input file $bam_file"
fi

##Step 7 Get raw bed file and Add cell barcode to the fragment bed files
if [ ! -f "$name.uniqmapped.RawBed.mito.Sort.Tag" ]; then
  echo "Running step 1 (original step 7): Get raw bed file and add cell barcode to the fragment bed files..."
  samtools sort -@ $CORE -n $bam_file | bedtools bamtobed -bedpe -i stdin | awk -v OFS='\t' '{split($7,name,"|"); print name[2],$1,$2,$6,$7}' | sort -k1,1 -k2,2 -k3,3n -k4,4n -k5,5 > $name.uniqmapped.RawBed.mito.Sort.Tag
else
  echo "$name.uniqmapped.RawBed.mito.Sort.Tag exists. Skip Step 1."
fi

##Step 2 deduplicate at single cell monoclonal tsv
if [ ! -f "$name.uniqmapped.fragment.$Cut.cut.mito.tsv" ]; then
  echo "Running step 2 (original step 8): Deduplicate at single cell monoclonal tsv..."
  python3 $MyMultiome/DeduplicateRawBed.10X.py $name.uniqmapped.RawBed.mito.Sort.Tag $name.uniqmapped.fragment.$Cut.cut.mito.tsv --cutoff $Cut 
else
  echo "$name.uniqmapped.fragment.$Cut.cut.mito.tsv exists. Skip Step 2."
fi

##Step 3 Summarize
if [ ! -f "$name.uniqmapped.fragment.$Cut.cut.mito.summary" ]; then
  echo "Running step 3 (original step 9): Summarize..."
  cat $name.uniqmapped.fragment.$Cut.cut.mito.tsv | python3 $MyMultiome/Summarize.TagDedup.10X.py > $name.uniqmapped.fragment.$Cut.cut.mito.summary
else
  echo "$name.uniqmapped.fragment.$Cut.cut.mito.summary exists. Skip step 3."
fi

##Step 4 ReadsCount
if [ ! -f "$name.ReadsCounts" ]; then
  echo "Running step 4 (original step 10): ReadsCount..."
  $MyMultiome/Counts.2.sh $bam_file $name.ReadsCounts $CORE
else
  echo "$name.ReadsCounts exists. Skip step 4."
fi

##Step 5 Plot QC
if [ ! -f "$name.QCplot.png" ]; then
  echo "Running step 5 (original step 11): Plot QC..."
  $MyMultiome/MultiATAC_mito.QC_v2.R $name $name.ReadsCounts $name.uniqmapped.fragment.$Cut.cut.mito.summary $name.uniqmapped.fragment.$Cut.cut.mito.tsv
else
  echo "$name.QCplot.png exists. Skip Step 5"
fi

##Final step cleanup
if [ ! -f "$name.QCplot.png" ]; then
  echo "Script running with error. Please check Log."
else
  echo "Script run complete successfully and got QC plot."
  echo "Cleaning up..."
  rm -rf tmp
  rm -rf *uniqmapped.RawBed
fi