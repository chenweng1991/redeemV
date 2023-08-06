#!/bin/bash
## Note:  This script is to QC ReDeeM mito reads. 
## It works with Multiome_ATAC_mito.sh, only if -q is enabled
## It basically only run from step 6 from Multiome_ATAC_mito.sh

Help()
{
  # Display Help
  echo "This script QC mito reads and generate QC plots "
  echo
  echo "MultiomeATAC_mitoQC.sh -h for this page"
  echo "Syntax: MultiomeATAC_mito.sh -n -c -m"
  echo "Options"
  echo "-n name: The prefix of all analyzed files"
  echo "-c Cut: the cutoff the uniq fragment per cell"
  echo "-t CORE: The number of cores to use"
  echo "-m MyMultiome:The path to the folder of MyMultiome"
}

while getopts "hn:1:2:i:c:t:m:b:q:p" option; do
  case $option in
    h) # display help
        Help
        exit;;
    n) # The prefix of all analyzed files
        name=$OPTARG;;
    c) # The cutoff the uniq fragment per cell
        Cut=$OPTARG;;
    t) # The number of cores to use
        CORE=$OPTARG;;
    m) # The path to the folder of MyMultiome
        MyMultiome=$OPTARG;;
    b) # The bowtie2 index path/prefix
        bowtie2Index=$OPTARG;;
   \?) # Invalid option
        echo "Error: Invalid option"
        exit;;
  esac
done

##Step6 Get raw bed file and Add cell barcode to the fragment bed files
samtools sort -@ $CORE -n $name.uniqmapped.mito.bam | bedtools bamtobed -bedpe -i stdin | awk -v OFS='\t' '{split($7,name,"|"); print name[2],$1,$2,$6,$7}' | sort -k1,1 -k2,2 -k3,3n -k4,4n -k5,5 > $name.uniqmapped.RawBed.mito.Sort.Tag


#Step7 deduplicate at single cell monoclonal tsv
echo "##Step7 deduplicate at single cell monoclonal tsv"
python3 $MyMultiome/DeduplicateRawBed.10X.py $name.uniqmapped.RawBed.mito.Sort.Tag $name.uniqmapped.fragment.$Cut.cut.mito.tsv --cutoff $Cut 


##Step10 Summarize
echo "Step10 Summarize"
cat $name.uniqmapped.fragment.$Cut.cut.mito.tsv | python3 $MyMultiome/Summarize.TagDedup.10X.py > $name.uniqmapped.fragment.$Cut.cut.mito.summary


##Step11 ReadsCount
echo "Step11 ReadsCount"
$MyMultiome/Counts.2.sh $name

##Simple Plot QC
echo "Step12 Plot QC"
$MyMultiome/MultiATAC_mito.QC_v2.R $name $name.ReadsCounts $name.uniqmapped.fragment.$Cut.cut.mito.summary $name.uniqmapped.fragment.$Cut.cut.mito.tsv

##cleanup
rm -rf tmp
rm *trim
rm -rf *tmp.bam
rm -rf *uniqmapped.RawBed
rm -rf *uniqmapped.fragment.0.cut.tsv
rm -rf *uniqmapped.fragment.0.cut.summary