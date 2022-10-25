#!/bin/bash

## Note: I think the following 4 files can be removed without any influence
# $name.bam   # We are fine as long as we have the $name.tagged.bam
# $name.bam.bai
# $name.uniqmapped.mito.bam   ## This is can be easily generated from $name.uniqmapped.bam
# $name.uniqmapped.RawBed
name=$1

echo $name.bam   
echo $name.bam.bai
echo $name.uniqmapped.mito.bam
echo $name.uniqmapped.RawBed

rm -rf $name.bam   
rm -rf $name.bam.bai
rm -rf $name.uniqmapped.mito.bam
rm -rf $name.uniqmapped.RawBed