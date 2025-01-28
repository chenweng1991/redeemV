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
  echo "-q quick, defult is false, if true then exit after uniqmapped.mito.bam, skipping QC step"
  echo "-p premap, default is false, if true then exit after mapping."
}

quick=0
premap=0
while getopts "hn:1:2:i:c:t:m:b:q:p" option; do
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
set -u
: "$name$Read1$Read2$ReadBarcode$Cut$CORE$MyMultiome$bowtie2Index$quick" 


##Step 1 trim adaptor (Important)
if [ ! -f "$Read1.trim" ]; then
  echo "Running step 1 trim adaptor..."
  cutadapt --cores=$CORE -a CTGTCTCTTATA -A CTGTCTCTTATA -o $Read1.trim -p $Read2.trim $Read1 $Read2
else
  echo "$Read1.trim and $Read2.trim exist. Skip Step 1"
fi

##Step 2 Add cell barcode to the readname
if [ ! -f "$Read1.trim.BC" ]; then
  echo "Running step 2 Add cell barcode to the readname..."
  python3 $MyMultiome/AddBC2Fastq.py $Read1.trim $Read2.trim $ReadBarcode $Read1.trim.BC $Read2.trim.BC
else
  echo "$Read1.trim.BC and $Read2.trim.BC exist. Skip Step 2."
fi

##Step3 Mapping Sorting and Indexing
#bowtie2Index=/lab/solexa_weissman/cweng/Genomes/GRCH38/GRCH38_Bowtie2_MitoMask/hg38.mitoMask
if [ ! -f "$name.bam" ]; then
  echo "Running step3 Mapping Sorting and Indexing..."
  bowtie2 -X 1200  --very-sensitive -p $CORE -x $bowtie2Index -1 $Read1.trim.BC  -2 $Read2.trim.BC | samtools view -@ $CORE -bS - > $name.tmp.bam
  samtools sort -@ $CORE $name.tmp.bam > $name.bam
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
if [ ! -f "$name.bam" ]; then
  echo "Running step 4 Extract cell barcode..."
  python3 $MyMultiome/AddBC2BAM.py $name.bam $name.tagged.bam
else
  echo "$name.tagged.bam exist. Skip Step 4."
fi

##Step5 Get uniq mapped bam 
if [ ! -f "$name.uniqmapped.bam" ]; then
  echo "Runing step 5 Get uniq mapped bam..."
  samtools view -@ $CORE -bf 2 -q30 $name.tagged.bam > $name.uniqmapped.bam 
else
  echo "$name.uniqmapped.bam exist. Skip Step 5."
fi

##Step6 Get Mito uniqmapped.bam
if [ ! -f "$name.uniqmapped.mito.bam" ]; then
  echo "Running step 6 Get Mito uniqmapped.bam..."
  samtools index -@ $CORE $name.uniqmapped.bam
  samtools view -@ $CORE -b $name.uniqmapped.bam chrM > $name.uniqmapped.mito.bam
else
  echo "$name.uniqmapped.mito.bam. Skip Step 6."
fi

# Exit here if -q option is enabled
if [[ quick -eq 1 ]]
  then
    echo "Skipping QC steps, exit"
    rm -rf *tmp.bam
    rm *trim
    exit
  else
    echo "Starting QC steps"
fi

##Step7 Get raw bed file and Add cell barcode to the fragment bed files
if [ ! -f "$name.uniqmapped.RawBed.mito.Sort.Tag" ]; then
  echo "Running step 7 Get raw bed file and Add cell barcode to the fragment bed files..."
  samtools sort -@ $CORE -n $name.uniqmapped.mito.bam | bedtools bamtobed -bedpe -i stdin | awk -v OFS='\t' '{split($7,name,"|"); print name[2],$1,$2,$6,$7}' | sort -k1,1 -k2,2 -k3,3n -k4,4n -k5,5 > $name.uniqmapped.RawBed.mito.Sort.Tag
else
  echo "$name.uniqmapped.RawBed.mito.Sort.Tag. Skip Step 7."
fi

#Step8 deduplicate at single cell monoclonal tsv
if [ ! -f "$name.uniqmapped.fragment.$Cut.cut.mito.tsv" ]; then
  echo "Running step7 deduplicate at single cell monoclonal tsv..."
  python3 $MyMultiome/DeduplicateRawBed.10X.py $name.uniqmapped.RawBed.mito.Sort.Tag $name.uniqmapped.fragment.$Cut.cut.mito.tsv --cutoff $Cut 
else
  echo "$name.uniqmapped.fragment.$Cut.cut.mito.tsv eise. Skip Step 8."
fi

##Step9 Summarize
if [ ! -f "$name.uniqmapped.fragment.$Cut.cut.mito.summary" ]; then
  echo "Running step9 Summarize..."
  cat $name.uniqmapped.fragment.$Cut.cut.mito.tsv | python3 $MyMultiome/Summarize.TagDedup.10X.py > $name.uniqmapped.fragment.$Cut.cut.mito.summary
else
  echo "$name.uniqmapped.fragment.$Cut.cut.mito.summary exist. Skip step 9."
fi

##Step10 ReadsCount
if [ ! -f "$name.ReadsCounts" ]; then
  echo "Running step10 ReadsCount..."
  $MyMultiome/Counts.2.sh $name
else
  echo "$name.ReadsCounts exist. Skip step 10."
fi

##Step 11 Plot QC
if [ ! -f "$name.QCplot.png" ]; then
  echo "Running step 11 Plot QC..."
  $MyMultiome/MultiATAC_mito.QC_v2.R $name $name.ReadsCounts $name.uniqmapped.fragment.$Cut.cut.mito.summary $name.uniqmapped.fragment.$Cut.cut.mito.tsv
else
  echo "$name.QCplot.png exist. Skip Step 11"
fi



##Final step cleanup
if [ ! -f "$name.QCplot.png" ]; then
  echo "Script running with error. Please check Log."
else
  echo "Script run complete successfully and got QC plot."
  echo "cleaning up."
  rm -rf tmp
  rm *trim
  rm -rf *tmp.bam
  rm -rf *uniqmapped.RawBed
  rm -rf *uniqmapped.fragment.0.cut.tsv
  rm -rf *uniqmapped.fragment.0.cut.summary
fi