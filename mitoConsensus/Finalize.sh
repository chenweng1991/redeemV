#!/bin/bash
WD=$1
Threads=$2
mitoConsensus=$3

if [ -z "$2" ]
then
	echo "Number of splited bams are not supplied, please provide, or count number of bams in ./mitoV/temp/barcoded_bams/"
	exit 1
fi


cat $WD/temp/sparse_matrices2.0/*QualifiedTotalCts >$WD/final/QualifiedTotalCts
cat $WD/temp/sparse_matrices2.0/*RawGenotypes.Total >$WD/final/RawGenotypes.Total
cat $WD/temp/sparse_matrices2.0/*RawGenotypes.VerySensitive >$WD/final/RawGenotypes.VerySensitive
cat $WD/temp/sparse_matrices2.0/*RawGenotypes.Sensitive >$WD/final/RawGenotypes.Sensitive
cat $WD/temp/sparse_matrices2.0/*RawGenotypes.Specific >$WD/final/RawGenotypes.Specific

python3 $mitoConsensus/AddQualifiedTotalCts.py $WD
Rscript $mitoConsensus/StrandBiases.R $WD 200 0.01
python $mitoConsensus/RemoveStrandBias.py $WD

##Final counting (Now it is included in the Finalize.sh script)
let N=0
for i in $(seq 1 $Threads)
        do
        echo $i
        let N=N+`samtools view $WD/temp/barcoded_bams/barcodes.$i.bam | wc -l`
        done
echo -e $N > $WD/final/TotalRawBamRows.Mito

