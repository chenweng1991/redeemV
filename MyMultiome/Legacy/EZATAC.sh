#!/bin/bash
## Note:  this version is for miseq or novaseq,  nextseq is slightly different
name=$1
MyMultiome=/lab/solexa_weissman/cweng/Packages/MyMultiome
lib=/lab/solexa_weissman/cweng/Packages/cxw486/scATAClib/



##Step2 Mapping
bowtie2Index=/lab/solexa_weissman/cweng/Genomes/GRCH38/GRCH38_Bowtie2_MitoMask/hg38.mitoMask
bowtie2 -X 1200  --very-sensitive -p 6 -x $bowtie2Index -1 $name"_"R1_001.fastq.gz -2 $name"_"R2_001.fastq.gz  |  samtools view -u -  |  samtools sort -  >$name.bam
samtools index $name.bam



##Step4 Get uniq mapped bam and make bulk bigwig and call peaks
echo "##Step4 Get uniq mapped bam and make bulk bigwig and call peaks"
samtools view -bf 2 -q30 $name.bam > $name.uniqmapped.bam   #309450630/=82% uniq and properly paired
$lib/Bam2bw.sh $name.uniqmapped
