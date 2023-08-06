name=$1
TotalReadsPair=$(expr `samtools view $name.bam | wc -l` / 2)
UniqProperpairReadsPair=`less $name.RawBed | wc -l`
MitoReadsPair=`less $name.RawBed | grep chrM | wc -l`

echo -ne "TotalReads\t$TotalReadsPair\n" >ReadsCounts
echo -ne "UniqProperpairReads\t$UniqProperpairReadsPair\n" >>ReadsCounts
echo -ne "MitoReadsPair\t$MitoReadsPair\n" >>ReadsCounts
