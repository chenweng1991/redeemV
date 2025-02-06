inbam=$1
outfile=$2
ncores=$3

if [ ! -s "$inbam.bai" ]; then
    samtools index $inbam
fi

TotalReadsPair=$(expr `samtools view -@ $ncores -c $inbam` / 2)
TotalMitoReadsPair=$(expr `samtools view -@ $ncores -c $inbam chrM` / 2)
UniqProperpairReadsPair=$(expr `samtools view -@ $ncores -c $inbam -bf 2 -q30` / 2)
MitoReadsPair=$(expr `samtools view -@ $ncores -c $inbam -bf 2 -q30 chrM` / 2)

echo -ne "TotalReads\t$TotalReadsPair\n" >$outfile
echo -ne "TotalMitoReads\t$TotalMitoReadsPair\n" >>$outfile
echo -ne "UniqProperpairReads\t$UniqProperpairReadsPair\n" >>$outfile
echo -ne "MitoReadsPair\t$MitoReadsPair\n" >>$outfile
