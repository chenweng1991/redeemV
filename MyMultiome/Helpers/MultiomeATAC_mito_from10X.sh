#!/bin/bash
## Note:  this version is starting from the 10x genomiocs atac bam file, and generate the fragment summary files needed for QC
bam=$1
CORE=24
# dic=$2
#Cut=0 # Minimum uniq fragment per cell to be considered
MyMultiome=/lab/solexa_weissman/cweng/Packages/MyMultiome
lib=/lab/solexa_weissman/cweng/Packages/cxw486/scATAClib/

# Step 1 extract dic table
#samtools view -bS $bam chrM > ${bam/.bam/.mito.bam}
#samtools index ${bam/.bam/.mito.bam} 
#python3 $MyMultiome/Helpers/ExtractDic.py ${bam/.bam/.mito.bam} CB > ${bam/.bam/.mito.dic}
python3 $MyMultiome/MergeBC2Bam.py  ${bam/.bam/.mito.dic} ${bam/.bam/.mito.bam} ${bam/.bam/.mito.bam}
samtools view -bf 2 -q30 ${bam/.bam/.mito.tagged.bam} > ${bam/.bam/.mito.tagged.unimapped.bam}


# Below comment out the order version for Step 1
python3 $MyMultiome/Helpers/ExtractDic.py $bam CB > ${bam/.bam/.dic}
python3 $MyMultiome/MergeBC2Bam.py  ${bam/.bam/.dic} $bam $bam
samtools index ${bam/.bam/.tagged.bam}
samtools view -bS ${bam/.bam/.tagged.bam} chrM > ${bam/.bam/.tagged.mito.bam}
samtools view -@ 24 -bf 2 -q30 ${bam/.bam/.tagged.mito.bam} > ${bam/.bam/.tagged.unimapped.mito.bam}  ## Added later, have not tested this yet.
rm -rf ${bam/.bam/.tagged.bam}
rm -rf ${bam/.bam/.tagged.bam.bai}

#Step2 Get raw bed file and Add cell barcode to the fragment bed files (Comment out)
echo "##Step2 Get raw bed file and Add cell barcode to the fragment bed files"
samtools sort -@ $CORE -n $bam| bedtools bamtobed -bedpe -i stdin | awk -v OFS='\t' '{print $1,$2,$6,$7}' > ${bam/.bam/.RawBed}
python3 $MyMultiome/MergeBC2Bed.10X.py ${bam/.bam/.dic}  ${bam/.bam/.RawBed}  | sort -k1,1 -k2,2 -k3,3n -k4,4n -k5,5 | awk '{if($3>0){print $0}}' >${bam/.bam/.RawBed.Sort.Tag}

Step3 deduplicate at single cell monoclonal tsv  (Comment out)
echo "##Step3 deduplicate at single cell monoclonal tsv"
python3 $MyMultiome/DeduplicateRawBed.10X.py ${bam/.bam/.RawBed.Sort.Tag} $Cut #  This will generate $name.uniqmapped.fragment.tsv and $name.uniqmapped.fragment.1000cut.tsv

#Step4 Extract Mito monoclonal tsv  (Comment out)
echo "Step4 Extract Mito monoclonal tsv"
grep chrM ${bam/.bam/.RawBed}.fragment.$Cut.cut.tsv > ${bam/.bam/.RawBed}.fragment.$Cut.cut.mito.tsv

#Step10 Summarize  (Comment out)
echo "Step5 Summarize"
cat ${bam/.bam/.RawBed}.fragment.$Cut.cut.tsv | python3 $MyMultiome/Summarize.TagDedup.10X.py > ${bam/.bam/.RawBed}.fragment.$Cut.cut.summary
cat ${bam/.bam/.RawBed}.fragment.$Cut.cut.mito.tsv | python3 $MyMultiome/Summarize.TagDedup.10X.py > ${bam/.bam/.RawBed}.fragment.$Cut.cut.mito.summary

#Step6 ReadsCount
echo "Step6 ReadsCount"
$MyMultiome/Counts.2.sh $name

##Simple Plot QC
echo "Step12 Plot QC"
$MyMultiome/MultiATAC_mito.QC.R $name $Cut

##cleanup
rm -rf tmp
