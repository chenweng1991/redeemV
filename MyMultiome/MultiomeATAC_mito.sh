#!/bin/bash
## Note:  this version is for miseq or novaseq,  nextseq is slightly different

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
echo "MultiomeATAC_mito.sh -h for help info"
set -u
: "$name$Read1$Read2$ReadBarcode$Cut$CORE$MyMultiome$bowtie2Index$quick" 


! test -e $name.bam && {
##Step 1 trim adaptor (Important)
cutadapt --cores=$CORE -a CTGTCTCTTATA -A CTGTCTCTTATA -o $Read1.trim -p $Read2.trim $Read1 $Read2

##Step2 Mapping
#bowtie2Index=/lab/solexa_weissman/cweng/Genomes/GRCH38/GRCH38_Bowtie2_MitoMask/hg38.mitoMask
bowtie2 -X 1200  --very-sensitive -p $CORE -x $bowtie2Index -1 $Read1.trim  -2 $Read2.trim | samtools view -@ $CORE -bS - > $name.tmp.bam
samtools sort -@ $CORE $name.tmp.bam > $name.bam
samtools index -@ $CORE $name.bam
}

# Exit here if -p option is enabled
if [[ premap -eq 1 ]]
  then
    echo "Mapping has completed, Only mapping, exit"
    exit
  else
    echo "Mapping has completed, Next------"
fi

#rm -rf $CORE $name.tmp.bam
#Step3 Extract cell barcode
echo "##Step3 Extract cell barcode"
python3 $MyMultiome/MakeBC.ATAC.I2.10X.Nextseq_Nova1.5.py ${ReadBarcode}  > $name.BC.dic
python3 $MyMultiome/MergeBC2Bam.py $name.BC.dic $name.bam  $name


##Step4 Get uniq mapped bam and make bulk bigwig and call peaks
echo "##Step4 Get uniq mapped bam and make bulk bigwig and call peaks"
samtools view -@ $CORE -bf 2 -q30 $name.tagged.bam > $name.uniqmapped.bam   #309450630/=82% uniq and properly paired
# $MyMultiome/Bam2bw.sh $name.uniqmapped   ## No need for mito, this is for making ATAC bigwig for visulization

##Step5 Get Mito uniqmapped.bam
echo "##Step5 Get uniq mapped bam and make bulk bigwig and call peaks"
samtools index -@ $CORE $name.uniqmapped.bam
samtools view -@ $CORE -b $name.uniqmapped.bam chrM > $name.uniqmapped.mito.bam

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

##Step6 Get raw bed file and Add cell barcode to the fragment bed files
echo "##Step6 Get raw bed file and Add cell barcode to the fragment bed files"
samtools sort -@ $CORE -n $name.uniqmapped.bam| bedtools bamtobed -bedpe -i stdin | awk -v OFS='\t' '{print $1,$2,$6,$7}' >$name.uniqmapped.RawBed
python3 $MyMultiome/MergeBC2Bed.10X.py $name.BC.dic  $name.uniqmapped.RawBed  | sort -k1,1 -k2,2 -k3,3n -k4,4n -k5,5 >$name.uniqmapped.RawBed.Sort.Tag

#Step7 deduplicate at single cell monoclonal tsv
echo "##Step7 deduplicate at single cell monoclonal tsv"
python3 $MyMultiome/DeduplicateRawBed.10X.py $name.uniqmapped.RawBed.Sort.Tag $Cut #  This will generate $name.uniqmapped.fragment.tsv and $name.uniqmapped.fragment.1000cut.tsv

##Step8 Extract Mito monoclonal tsv
echo "Step8 Extract Mito monoclonal tsv"
grep chrM $name.uniqmapped.fragment.$Cut.cut.tsv > $name.uniqmapped.fragment.$Cut.cut.mito.tsv

##Step9 Make Peak VS Cell sparse matrix
echo "Skip Step9 Make Peak VS Cell sparse matrix"
# mkdir tmp
# cat $name.uniqmapped.fragment.$Cut.cut.tsv | awk '{print > "tmp/"$4}'
# for f in `ls tmp/`
# do
#   bedtools intersect -a $name.uniqmapped_peaks.narrowPeak -b tmp/$f -c | awk -v name=$f '{if($11>0){print $1"_"$2"_"$3,name,$11}}' OFS="\t" >> $name.uniqmapped.Peak_bc_sparse_mtx
# done

##Step10 Summarize
echo "Step10 Summarize"
cat $name.uniqmapped.fragment.$Cut.cut.tsv | python3 $MyMultiome/Summarize.TagDedup.10X.py > $name.uniqmapped.fragment.$Cut.cut.summary
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
