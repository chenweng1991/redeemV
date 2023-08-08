#!/usr/bin/python
"""
chunk_barcoded_bam.py

Author: Chen Weng
Created Date: 2023-8-4
Last Updated: 2023-8-4
Description:
This script extracts and writes reads with specific barcodes to a new BAM file.
"""

import sys
import pysam
import os

def subset_bam(bamfile, outfolder, barcodeTag, bcfile, mtchr):
    """Extracts and writes reads with specific barcodes to a new BAM file."""
    
    basename = os.path.basename(os.path.splitext(bcfile)[0])
    outname = os.path.join(outfolder, f"{basename}.bam")

    # Read in the barcodes
    with open(bcfile, 'r') as barcode_file_handle:
        bc = {x.strip() for x in barcode_file_handle.readlines()}  # Using a set for O(1) lookup
    
    with pysam.AlignmentFile(bamfile, "rb") as bam, pysam.AlignmentFile(outname, "wb", header=bam.header) as out:
        for read in bam.fetch(str(mtchr),multiple_iterators=False):
            try:
                barcode_id = read.get_tag(barcodeTag)
                if barcode_id in bc:
                    out.write(read)
            except KeyError:
                pass  # Handle cases where the barcode tag doesn't exist for a read
            
    pysam.index(outname)

if __name__ == "__main__":
    bamfile = sys.argv[1]
    outfolder = sys.argv[2]
    barcodeTag = sys.argv[3]
    bcfile = sys.argv[4]
    mtchr = sys.argv[5]
    
    subset_bam(bamfile, outfolder, barcodeTag, bcfile, mtchr)
