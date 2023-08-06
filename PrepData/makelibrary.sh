#!/usr/bin/env bash
RNAPATH=$(realpath "`pwd`/../FASTQ/RNA/")
ATACPATH=$(realpath "`pwd`/../FASTQ/ATAC/")
Prefix=`dirname $(pwd) | xargs basename`
echo -e "fastqs,sample,library_type\n$RNAPATH,$Prefix,Gene Expression\n$ATACPATH,$Prefix,Chromatin Accessibility" > Libraries
