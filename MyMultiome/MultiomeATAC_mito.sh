#!/bin/bash
## Note:  this version is for illumina reverse complement workflow, including nextseq, Novaseq 6000 with v1.5 reagent kits (i.e., Reverse complement chemistry)
## Author: Chen Weng
## Date: 2022-10-27
## Updated: 2023-08-04
Help()
{
  # Display Help
  echo "This script trim and map the mtDNA fastqs, generating uniqmapped bam and QC plots "
  echo
  echo "MultiomeATAC_mito.sh -h for this page"
  echo "Syntax: MultiomeATAC_mito.sh -n -1 -2 -i -c -t -m -b"
  echo "Options"
  echo "-n name: The prefix of all analyzed files"
  echo "-1 Read1: Read1 of fastq file (150nt is recommended)"
  echo "-2 Read2: Read2 of fastq file (150nt is recommended)"
  echo "-i ReadBarcode: i5 of fastq file (24nt)"
  echo "-c Cut: the cutoff the uniq fragment per cell"
  echo "-t CORE: The number of cores to use"
  echo "-m MyMultiome:The path to the folder of MyMultiome"
  echo "-b bowtie2Index: the bowtie2 index path/prefix"
  echo "-o output directory"
  echo "-q quick, defult is false, if true then exit after uniqmapped.mito.bam, skipping QC step"
  echo "-p premap, default is false, if true then exit after mapping."
}

quick=0
premap=0
while getopts "hn:1:2:i:c:t:m:b:o:q:p" option; do
  case $option in
    h) # display help
        Help
        exit;;
    n) # The prefix of all analyzed files
        name=$OPTARG;;
    1) # Read1 of fastq file
        Read1=$OPTARG;;
    2) # Read2 of fastq file
        Read2=$OPTARG;;
    i) # index5 fastq file
        ReadBarcode=$OPTARG;;
    c) # The cutoff the uniq fragment per cell
        Cut=$OPTARG;;
    t) # The number of cores to use
        CORE=$OPTARG;;
    m) # The path to the folder of MyMultiome
        MyMultiome=$OPTARG;;
    b) # The bowtie2 index path/prefix
        bowtie2Index=$OPTARG;;
    o) # output folder
        outdir=$OPTARG;;
    q) # if use this option then exit after uniqmapped.mito.bam, skipping QC step
        quick=1;;
    p) # If use this option then exit after mapping, skipping the rest. 
        premap=1;;
   \?) # Invalid option
        echo "Error: Invalid option"
        exit;;
  esac
done

## Exit if any of the necessary input is empty

echo "This is Memory saving version updated 2023-8-4"
echo "Please run MultiomeATAC_mito.sh -h for help info"
echo "R1: $Read1"
echo "R2: $Read2"
echo "outdir: $outdir"
set -u
: "$name$Read1$Read2$ReadBarcode$Cut$CORE$MyMultiome$bowtie2Index$quick$outdir" 

# change to output dir
if [ ! -d "$outdir" ]; then
  mkdir -p "$outdir"
  echo "Directory $outdir created."
else
  echo "Directory $outdir already exists."
fi

cd $outdir

# stop if running into error
set -e

Read1T=$(basename "${Read1%.fastq.gz}").trim.fastq.gz
Read2T=$(basename "${Read2%.fastq.gz}").trim.fastq.gz

##Step 1 trim adaptor (Important)
if [ ! -s "$Read1T" ]; then
  echo "Running step 1 trim adaptor..."
  cutadapt --cores=$CORE -a CTGTCTCTTATA -A CTGTCTCTTATA -o $Read1T -p $Read2T $Read1 $Read2
else
  echo "$Read1T and $Read2T exist. Skip Step 1"
fi

Read1TB=$(basename "${Read1T%.fastq.gz}").bc.fastq.gz
Read2TB=$(basename "${Read2T%.fastq.gz}").bc.fastq.gz

##Step 2 Add cell barcode to the readname
if [ ! -s "$Read1TB" ]; then
  echo "Running step 2 Add cell barcode to the readname..."
  # python3 $MyMultiome/AddBC2Fastq.py $Read1T $Read2T $ReadBarcode $Read1TB $Read2TB
  for readfile in $Read1T $Read2T; do 
      outfile=$(basename "${readfile%.fastq.gz}").bc.fastq.gz
      paste <(zcat "$ReadBarcode" | awk 'NR%4==2' | cut -c9-24 | rev | tr 'ATGC' 'TACG') \
          <(zcat "$readfile" | awk 'NR%4==1 {split($0, a, " "); header=a[1]} NR%4==2 {seq=$0} NR%4==3 {plus=$0} NR%4==0 {qual=$0; print header, seq, plus, qual}') | \
      awk '{print $2 "|" $1 "\n" $3 "\n" $4 "\n" $5}' | gzip > "$outfile" &
  done
  wait;
else
  echo "$Read1TB and $Read2TB exist. Skip Step 2."
fi

##Step3 Mapping Sorting and Indexing; needs at least 5G of RAM
if [ ! -s "$name.bam" ]; then
  echo "Running step3 Mapping Sorting and Indexing..."
  bowtie2 -X 1200  --very-sensitive -p $CORE -x $bowtie2Index -1 $Read1TB -2 $Read2TB | samtools sort -@ $CORE -m 4G > $name.bam
  samtools index -@ $CORE $name.bam
else 
  echo "$name.bam exist. Skip Step 3."
fi

# Exit here if -p option is enabled
if [[ premap -eq 1 ]]
  then
    echo "Mapping has completed, Only mapping, exit"
    exit
  else
    echo "Mapping has completed, Next------"
fi

#Step4 Extract cell barcode
if [ ! -s "$name.tagged.bam" ]; then
  echo "Running step 4 Extract cell barcode..."
  python3 $MyMultiome/AddBC2BAM.py $name.bam $name.tagged.bam
else
  echo "$name.tagged.bam exist. Skip Step 4."
fi

##Step5 Get uniq mapped bam 
if [ ! -s "$name.uniqmapped.bam" ]; then
  echo "Runing step 5 Get uniq mapped bam..."
  samtools view -@ $CORE -bf 2 -q30 $name.tagged.bam > $name.uniqmapped.bam
else
  echo "$name.uniqmapped.bam exist. Skip Step 5."
fi

##Step6 Get Mito uniqmapped.bam
if [ ! -s "$name.uniqmapped.mito.bam" ]; then
  echo "Running step 6 Get Mito uniqmapped.bam..."
  samtools index -@ $CORE $name.uniqmapped.bam
  samtools view -@ $CORE -b $name.uniqmapped.bam chrM > $name.uniqmapped.mito.bam
else
  echo "$name.uniqmapped.mito.bam. Skip Step 6."
fi