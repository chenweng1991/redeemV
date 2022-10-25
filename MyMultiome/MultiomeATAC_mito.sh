#!/bin/bash
## Note:  this version is for miseq or novaseq,  nextseq is slightly different

echo "7 inputs: $1:Prefix, $2:Read1.fq, $3:Read2.fq, $4:i5.fq, $5:0; $6:threads, $7:MyMultiome path"
name=$1
Read1=$2
Read2=$3
ReadBarcode=$4
Cut=$5 # Minimum uniq fragment per cell to be considered
CORE=$6
MyMultiome=$7


##Step 1 trim adaptor (Important)
cutadapt --cores=$CORE -a CTGTCTCTTATA -A CTGTCTCTTATA -o $Read1.trim -p $Read2.trim $Read1 $Read2

##Step2 Mapping
#bowtie2Index=/lab/solexa_weissman/cweng/Genomes/GRCH38/GRCH38_Bowtie2_MitoMask/hg38.mitoMask
bowtie2 -X 1200  --very-sensitive -p $CORE -x $bowtie2Index -1 $Read1.trim  -2 $Read2.trim | samtools view -@ $CORE -bS - > $name.tmp.bam
#samtools sort -@ $CORE $name.tmp.bam > $name.bam
#samtools index -@ $CORE $name.bam

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
