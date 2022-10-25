#!/bin/bash
# Example MyMultiome=/lab/solexa_weissman/cweng/Packages/MyMultiome
# Example cd /lab/solexa_weissman/cweng/Projects/Collaborator/Florian_ArmstrongLab
# Exampole $MyMultiome/RunCWmgatk.sh patient_2_relapse.bam barcodes_patient2_relapse.tsv /lab/solexa_weissman/cweng/Projects/Collaborator/Florian_ArmstrongLab/patient2_relapse_CWmgatk2.0.Test.Wrap 48
## Step 1 will be the same as CW_mgatk
CW_mgatk=/lab/solexa_weissman/cweng/Packages/CW_mgatk/

##STEP1 Preprocess, this will create a folder CW_mgatk
bamfile=$1
Celltb=$2
WD=$3 ## The WD or working directory is a name of new folder to be created in the same folder of the bam
core=$4

bsub python $CW_mgatk/Preprocess.py -i $bamfile -c $core -b $Celltb -o $WD -g rCRS

## Step 2 Run CW_sumstatsBP2.0.py --> for each barcode group, 4 RawGenotypes files and 1 QualifiedTotalCts will be generated.
## Columns of
## WD=/lab/solexa_weissman/cweng/Packages/CW_mgatk/TestingData/CW_mgatk_test/  Working directory is created by $CW_mgatk/Preprocess.py
wait
for i in {1..$core}
do
bsub python $CW_mgatk/CW_sumstatsBP2.0.py barcodes.$i $WD BC 30
done
## Step 3 concat together
wait
cat $WD/temp/sparse_matrices2.0/*QualifiedTotalCts >$WD/final/QualifiedTotalCts
cat $WD/temp/sparse_matrices2.0/*RawGenotypes.Total >$WD/final/RawGenotypes.Total
cat $WD/temp/sparse_matrices2.0/*RawGenotypes.VerySensitive >$WD/final/RawGenotypes.VerySensitive
cat $WD/temp/sparse_matrices2.0/*RawGenotypes.Sensitive >$WD/final/RawGenotypes.Sensitive
cat $WD/temp/sparse_matrices2.0/*RawGenotypes.Specific >$WD/final/RawGenotypes.Specific
## Step4 remove strandbias variants by biomial modeling
python3 $CW_mgatk/AddQualifiedTotalCts.py $WD
Rscript $CW_mgatk/StrandBiases.R $WD 200 0.01
python $CW_mgatk/RemoveStrandBias.py $WD
