name=$1
TotalReadsPair=$(expr `samtools view $name.bam | wc -l` / 2)
TotalMitoReadsPair=$(expr `samtools view $name.bam |grep chrM | wc -l` / 2)
UniqProperpairReadsPair=$(expr `samtools view $name.uniqmapped.bam | wc -l` / 2)
MitoReadsPair=$(expr `samtools view $name.uniqmapped.mito.bam| wc -l` / 2)




echo -ne "TotalReads\t$TotalReadsPair\n" >$name.ReadsCounts
echo -ne "TotalMitoReads\t$TotalMitoReadsPair\n" >>$name.ReadsCounts
echo -ne "UniqProperpairReads\t$UniqProperpairReadsPair\n" >>$name.ReadsCounts
echo -ne "MitoReadsPair\t$MitoReadsPair\n" >>$name.ReadsCounts
