#! /bin/bash
name=$1
Genome=/lab/solexa_weissman/cweng/Genomes/GRCH38
##################################################### Make aggregate track
bedtools genomecov -bg -i $name -g  $Genome/hg38.chrom.sizes > ${name/.bed/.bg}
bedGraphToBigWig ${name/.bed/.bg} $Genome/hg38.chrom.sizes ${name/.bed/.bw}
macs2 callpeak -f BEDPE -t $name -g hs  -n ${name/.monoclonal.bed/}
# less $name.$genome'_'peaks.narrowPeak | cut -f1-4 > $name.$genome.narrowpeak.bed
# $scATAClib/bedToBigBed $name.$genome.narrowpeak.bed ~/Genome/$genome/$genome.chrom.sizes $name.$genome.narroepeak.bb