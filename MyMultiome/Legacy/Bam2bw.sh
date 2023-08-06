#! /bin/bash
picard=/lab/solexa_weissman/cweng/Packages/Picard/picard.jar
dropseq=/lab/solexa_weissman/cweng/Packages/cxw486/Drop-seq_tools-1.0
name=$1
scATAClib=/lab/solexa_weissman/cweng/Packages/cxw486/scATAClib
Genome=/lab/solexa_weissman/cweng/Genomes/GRCH38
##################################################### Make aggregate track
java -jar -Xmx4g $picard MarkDuplicates   -I $name.bam -O $name.dupMark.bam -M $name.matrix -REMOVE_DUPLICATES true
samtools sort -n $name.dupMark.bam |samtools view -bf 2  | bedtools bamtobed -bedpe -i stdin | awk -v OFS='\t' '{print $1,$2,$6,".",1,$9}' |sort -k1,1 -k2,2n -k3,3n>$name.monoclonal.bed
bedtools genomecov -bg -i $name.monoclonal.bed -g  $Genome/hg38.chrom.sizes > $name.monoclonal.bg
bedGraphToBigWig $name.monoclonal.bg $Genome/hg38.chrom.sizes $name.monoclonal.bw
macs2 callpeak -f BAMPE -t $name.dupMark.bam -g hs  -n $name
# less $name.$genome'_'peaks.narrowPeak | cut -f1-4 > $name.$genome.narrowpeak.bed
# $scATAClib/bedToBigBed $name.$genome.narrowpeak.bed ~/Genome/$genome/$genome.chrom.sizes $name.$genome.narroepeak.bb
