#!/bin/bash
## Note:  this version is for miseq or novaseq,  nextseq is slightly different
name=$1
Read1=$2
Read2=$3
MyMultiome=/lab/solexa_weissman/cweng/Packages/MyMultiome
lib=/lab/solexa_weissman/cweng/Packages/cxw486/scATAClib/


##Step 1 trim adaptor (Important)
cutadapt --cores=8 -a CTGTCTCTTATA -A CTGTCTCTTATA -o ${Read1/fastq.gz/trim.fastq.gz} -p ${Read2/fastq.gz/trim.fastq.gz} $Read1 $Read2

##Step2 Mapping
bowtie2Index=/lab/solexa_weissman/cweng/Genomes/GRCH38/GRCH38_Bowtie2_MitoMask/hg38.mitoMask
bowtie2 -X 1200  --very-sensitive -p 16 -x $bowtie2Index -1 ${Read1/fastq.gz/trim.fastq.gz}  -2 ${Read2/fastq.gz/trim.fastq.gz}  |  samtools view -u -  |  samtools sort -  >$name.bam
samtools index $name.bam

##Step4 Get uniq mapped bam and make bulk bigwig and call peaks
echo "##Step4 Get uniq mapped bam"
samtools view -bf 2 -q30 $name.bam > $name.uniqmapped.bam
$lib/Bam2bw.sh $name.uniqmapped

##Step5 Get Mito uniqmapped.bam
echo "##Step5 Get uniq mapped bam and make bulk bigwig and call peaks"
samtools index $name.uniqmapped.bam
samtools view -b $name.uniqmapped.bam chrM > $name.uniqmapped.mito.bam

##Step11 ReadsCount
echo "Step6 ReadsCount"
$MyMultiome/Counts.2.sh $name

##Simple Plot QC
echo "Step12 Plot QC"
$MyMultiome/BulkATAC_Map_SummaryPlot.R $name
