#!/bin/bash
## Note:  this version is for illumina reverse complement workflow, including nextseq, Novaseq 6000 with v1.5 reagent kits (i.e., Reverse complement chemistry)
## Author: Chen Weng
## Date: 2022-10-27
## Updated: 2023-08-04
Help()
{
 # Display Help
 echo "This script trim and map the mtDNA fastqs, generating uniqmapped bam and QC plots "
 echo
 echo "MultiomeATAC_mito.sh -h for this page"
 echo "Syntax: MultiomeATAC_mito.sh -h -n -1 -2 -i -c -t -m -b -x -o -q -p -a -M"
 echo "Options"
 echo "-h : Display this help page"
 echo "-n name: The prefix of all analyzed files"
 echo "-1 Read1: Read1 of fastq file (150nt is recommended)"
 echo "-2 Read2: Read2 of fastq file (150nt is recommended)"
 echo "-i ReadBarcode: i5 of fastq file (24nt)"
 echo "-c Cut: the cutoff the uniq fragment per cell"
 echo "-t CORE: The number of cores to use"
 echo "-m MyMultiome: The path to the folder of MyMultiome"
 echo "-b regMitoIndex: the bowtie2 or bwa-mem index path/prefix"
 echo "-x shiftedMitoIndex: the  shifted mito index path/prefix"
 echo "-M MitoIndexOnly: the mito index prefix"
 echo "-o output directory"
 echo "-q quick, default is false, if true then exit after uniqmapped.mito.bam, skipping QC step"
 echo "-p premap, default is false, if true then exit after mapping."
 echo "-a aligner: Choose aligner (bowtie2 or bwa-mem2) (we've mostly used bowtie2 in the past)"

}

quick=0
premap=0
while getopts "hn:1:2:i:c:t:m:o:qpa:b:x:" option; do
    case $option in
    h) # display help
        Help
        exit;;
    n) # The prefix of all analyzed files
        name=$OPTARG;;
    1) # Read1 of fastq file
        Read1=$OPTARG;;
    2) # Read2 of fastq file
        Read2=$OPTARG;;
    i) # index5 fastq file
        ReadBarcode=$OPTARG;;
    c) # The cutoff the uniq fragment per cell
        Cut=$OPTARG;;
    t) # The number of cores to use
        CORE=$OPTARG;;
    m) # The path to the folder of MyMultiome
        MyMultiome=$OPTARG;;
    o) # output folder
        outdir=$OPTARG;;
    q) # if use this option then exit after uniqmapped.mito.bam, skipping QC step
        quick=1;;
    p) # If use this option then exit after mapping, skipping the rest. 
        premap=1;;
    a) # Validate aligner choice immediately
      if [ "$OPTARG" != "bwamem2" ] && [ "$OPTARG" != "bowtie2" ]; then
          echo "Error: aligner (-a) must be either 'bwamem2' or 'bowtie2'"
          exit 1
      fi
      aligner=$OPTARG;;
    b) # The bowtie2 index path/prefix
        regMitoIndex=$OPTARG;;
    x) # the bowtie2 shifted mito ref index path 
        shiftedMitoIndex=$OPTARG;;
   \?) # Invalid option
        echo "Error: Invalid option"
        exit;;
  esac
done

## Exit if any of the necessary input is empty

echo "This is Memory saving version updated 2023-8-4"
echo "Please run MultiomeATAC_mito.sh -h for help info"
echo "R1: $Read1"
echo "R2: $Read2"
echo "outdir: $outdir"
set -u
echo "aligner: $aligner" 
if [ "$aligner" = "bowtie2" ]; then
   : "$name$Read1$Read2$ReadBarcode$Cut$CORE$MyMultiome$regMitoIndex$quick$outdir"
elif [ "$aligner" = "bwamem2" ]; then
   : "$name$Read1$Read2$ReadBarcode$Cut$CORE$MyMultiome$regMitoIndex$quick$outdir"
fi

if [ ! -d "$outdir" ]; then
  mkdir -p "$outdir"
  echo "Directory $outdir created."
else
  echo "Directory $outdir already exists."
fi

echo "reg Mito index: $regMitoIndex"

cd $outdir

# stop if running into error
set -e

Read1T=$(basename "${Read1%.fastq.gz}").trim.fastq.gz
Read2T=$(basename "${Read2%.fastq.gz}").trim.fastq.gz

##Step 1 trim adaptor (Important)
if [ ! -s "$Read1T" ]; then
  echo "Running step 1 trim adaptor..."
  cutadapt --cores=$CORE -a CTGTCTCTTATA -A CTGTCTCTTATA -o $Read1T -p $Read2T $Read1 $Read2
else
  echo "$Read1T and $Read2T exist. Skip Step 1"
fi

Read1TB=$(basename "${Read1T%.fastq.gz}").bc.fastq.gz
Read2TB=$(basename "${Read2T%.fastq.gz}").bc.fastq.gz

##Step 2 Add cell barcode to the readname
if [ ! -s "$Read1TB" ]; then
  echo "Running step 2 Add cell barcode to the readname..."
  # python3 $MyMultiome/AddBC2Fastq.py $Read1T $Read2T $ReadBarcode $Read1TB $Read2TB
  for readfile in $Read1T $Read2T; do 
      outfile=$(basename "${readfile%.fastq.gz}").bc.fastq.gz
      paste <(zcat "$ReadBarcode" | awk 'NR%4==2' | cut -c9-24 | rev | tr 'ATGC' 'TACG') \
          <(zcat "$readfile" | awk 'NR%4==1 {split($0, a, " "); header=a[1]} NR%4==2 {seq=$0} NR%4==3 {plus=$0} NR%4==0 {qual=$0; print header, seq, plus, qual}') | \
      awk '{print $2 "|" $1 "\n" $3 "\n" $4 "\n" $5}' | gzip > "$outfile" &
  done
  wait;
else
  echo "$Read1TB and $Read2TB exist. Skip Step 2."
fi

##Step3 Mapping Sorting and Indexing; needs at least 5G of RAM
## uses bowtie2 or bwamem2 based on indicated input
if [ ! -s "${name}.bam" ]; then
 echo "Running step3 Mapping Sorting and Indexing using $aligner..."
 if [ "$aligner" = "bowtie2" ]; then
   echo "Using bowtie2 for alignment..."
     bowtie2 -X 1200 --very-sensitive -p $CORE -x $regMitoIndex \
     -1 $Read1TB -2 $Read2TB | \
     samtools sort -@ $CORE -m 4G > ${name}.bam
 elif [ "$aligner" = "bwamem2" ]; then
  echo "Using bwa-mem2 for alignment..."
  
  # Create an intermediate SAM file instead of direct piping
  # Write SAM output to a file first
  bwa-mem2 mem -t $CORE $regMitoIndex $Read1TB $Read2TB > $name.sam

  # Convert SAM to BAM
  samtools view -b $name.sam | samtools sort -@ $CORE -m 4G > $name.bam
  

    
 
 else
   echo "Error: Invalid aligner specified: $aligner"
   exit 1
 fi
 samtools index -@ $CORE $name.bam
else 
 echo "$name.bam exist. Skip Step 3."
fi

# Exit here if -p option is enabled
if [[ premap -eq 1 ]]
 then
   echo "Mapping has completed, Only mapping, exit"
   exit
 else
   echo "Mapping has completed, Next------"
fi
#Step4 Extract cell barcode
if [ ! -s "$name.tagged.bam" ]; then
  echo "Running step 4 Extract cell barcode..."
  python3 $MyMultiome/AddBC2BAM.py $name.bam $name.tagged.bam
else
  echo "$name.tagged.bam exist. Skip Step 4."
fi

##Step5 Get uniq mapped bam - updated to save the non mapped lines for re-mapping
if [ ! -s "$name.uniqmapped.bam" ]; then
  echo "Runing step 5 Get uniq mapped bam..."
  samtools view -@ $CORE -bf 2 -q30 -U $name.NOT_uniqmapped.bam $name.tagged.bam > $name.uniqmapped.bam
else
  echo "$name.uniqmapped.bam exist. Skip Step 5."
fi

##Step6 Get Mito uniqmapped.bam
if [ ! -s "$name.uniqmapped.mito.bam" ]; then
  echo "Running step 6 Get Mito uniqmapped.bam..."
  samtools index -@ $CORE $name.uniqmapped.bam
  samtools view -@ $CORE -b $name.uniqmapped.bam chrM > $name.uniqmapped.mito.bam
else
  echo "$name.uniqmapped.mito.bam. Skip Step 6."
fi



## Step 7: Remap unmapped reads to shifted mito reference

### Step 7A: Convert unmapped BAM to FASTQ pairs ###
if [ ! -s "${name}_remap_R1.fastq.gz" ]; then
 echo "Running step 7A: Converting unmapped BAM to FASTQ..."
 samtools fastq -@ $CORE -1 ${name}_remap_R1.fastq.gz -2 ${name}_remap_R2.fastq.gz $name.NOT_uniqmapped.bam
else
 echo "${name}_remap_R1.fastq.gz exists. Skip Step 7A."
fi

## Step 7B: Remap to shifted mitochondrial reference ###
if [ ! -s "${name}_remapped.bam" ]; then
 echo "Running step 7B: Remapping to shifted mito reference using $aligner..."
 if [ "$aligner" = "bowtie2" ]; then
   echo "Using bowtie2 for shifted mito alignment..."
   bowtie2 -X 1200 --very-sensitive -p $CORE -x $shiftedMitoIndex \
     -1 ${name}_remap_R1.fastq.gz -2 ${name}_remap_R2.fastq.gz | \
     samtools sort -@ $CORE -m 4G > ${name}_remapped.bam
 elif [ "$aligner" = "bwamem2" ]; then
   echo "Using bwa-mem2 for shifted mito alignment..."
  bwa-mem2 mem -t $CORE $shiftedMitoIndex $Read1TB $Read2TB > ${name}_remapped.sam
  samtools view -b ${name}_remapped.sam | samtools sort -@ $CORE -m 4G > ${name}_remapped.bam
 else
   echo "Error: Invalid aligner specified: $aligner"
   exit 1
 fi
else
 echo "${name}_remapped.bam exists. Skip Step 7B."
fi


### Step 7C: Filter remapped reads, shift genome back, and cleanup ###
if [ ! -s "$name.reshifted.alt_mito_mapped_uniqmapped.bam" ]; then
 echo "Running step 7C: Filtering remapped reads..."
samtools view -@ $CORE -bf 2 -q30 \
  -U ${name}.alt_mito_filtered_OUT.bam \
  ${name}_remapped.bam > ${name}.alt_mito_mapped_uniqmapped.bam

 samtools index -@ $CORE ${name}.alt_mito_mapped_uniqmapped.bam

 python3 $MyMultiome/revert_mito_shift.py ${name}.alt_mito_mapped_uniqmapped.bam ${name}.reshifted.alt_mito_mapped_uniqmapped.bam



 # Clean up intermediate files after successful filtering
 echo "Cleaning up intermediate remapping files..."
  # rm -f ${name}_remap_R1.fastq.gz ${name}_remap_R2.fastq.gz ${name}_remapped.bam ${name}.alt_mito_mapped_uniqmapped.bam
else
 echo "$name.alt_mito_mapped_uniqmapped.bam exists. Skip Step 7C."
fi

#Step7D Extract cell barcode for alt 
if [ ! -s "$name.alt_mito_mapped_uniqmapped.tagged.bam" ]; then
  echo "Running step 7D Extract cell barcode..."
  python3 $MyMultiome/AddBC2BAM.py ${name}.reshifted.alt_mito_mapped_uniqmapped.bam ${name}.reshifted.alt_mito_mapped_uniqmapped.tagged.bam
else
  echo "${name}.reshifted.alt_mito_mapped_uniqmapped.tagged.bam exist. Skip Step 7D."
fi

#Step8 make output version that is combined 
if [ ! -s "${name}.merged_with_shifted_mito.sorted.bam" ]; then
  echo "Running step 8 combine bams..."

    # For shifted BAM
  samtools addreplacerg \
    -r "ID:shifted" \
    -o "${name}.shifted.tmp.bam" \
    "${name}.reshifted.alt_mito_mapped_uniqmapped.tagged.bam"
  mv "${name}.shifted.tmp.bam" "${name}.reshifted.alt_mito_mapped_uniqmapped.tagged.bam"



  # Merge
  samtools merge -f "${name}.merged_with_shifted_mito.bam" \
    "${name}.reshifted.alt_mito_mapped_uniqmapped.tagged.bam" \
    "${name}.uniqmapped.mito.bam"

  # Sort
  samtools sort -o "${name}.merged_with_shifted_mito.sorted.bam" \
    "${name}.merged_with_shifted_mito.bam"

  # Index
  samtools index "${name}.merged_with_shifted_mito.sorted.bam"




else
  echo "${name}.reshifted.alt_mito_mapped_uniqmapped.tagged.bam Skip step 8."

fi

# ### Step 8: alternate mito mapping to only mitochondrial reference
# if [ ! -s "$name.mito_ONLY_mapped.bam" ]; then
#  echo "Running step 8  Mapping sorting and indexing using only mito reference..."
#  if [ "$aligner" = "bowtie2" ]; then
#    echo "Using bowtie2 for alignment..."
#    bowtie2 -X 1200 --very-sensitive -p $CORE -x $MitoIndexOnly \
#      -1 $Read1TB -2 $Read2TB | \
#      samtools sort -@ $CORE -m 4G | samtools view -@ $CORE -bf 2 -q30 > $name.mito_ONLY_mapped.bam
#      samtools index -@ $CORE $name.mito_ONLY_mapped.bam
#  else
#    echo "Skipping Step 8: This step requires bowtie2 aligner, but $aligner was specified."
#  fi
 
# else 
#  echo "$name.alt_mito_ONLY_mapped.bam exist. Skip Step 8."
# fi


# # Step 8B: Extract cell barcode for alt
# if [ "$aligner" = "bowtie2" ]; then
#   if [ ! -s "$name.mito_ONLY_mapped.tagged.bam" ]; then
#     echo "Running step 4 Extract cell barcode..."
#     python3 $MyMultiome/AddBC2BAM.py $name.mito_ONLY_mapped.bam $name.mito_ONLY_mapped.tagged.bam
#   else
#     echo "$name.mito_ONLY_mapped.tagged.bam exists. Skip Step 4."
#   fi
#   echo "aligner is not bowtie2, so this file is not needed and we skip this step." 
# fi
