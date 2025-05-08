#!/usr/bin/python
"""
chunk_barcoded_bam.py

Author: Chen Weng
Created Date: 2023-8-4
Last Updated: 2025-3-10
Description:
This script extracts and writes reads with specific barcodes to a new BAM file.
optimized to run on a single cluster core, not loading all reads into memory
"""

import sys
import pysam
import os
import time
import tempfile
import subprocess

def read_barcodes(bcfile):
    """Read barcodes from file into a set for efficient lookup."""
    print(f"Reading barcodes from {bcfile}")
    try:
        with open(bcfile, 'r') as barcode_file_handle:
            bc = {x.strip() for x in barcode_file_handle.readlines() if x.strip()}
        print(f"Loaded {len(bc)} barcodes")
        return bc
    except Exception as e:
        print(f"Error reading barcode file: {e}")
        raise

def subset_bam(bamfile, outfolder, barcodeTag, bcfile, mtchr, prefix="", chunk_size=3000):
    """
    Extracts and writes reads with specific barcodes to a new BAM file.
    Uses chunking for better memory efficiency and ensures proper sorting.
    """
    start_time = time.time()
    
    # Create output directory if it doesn't exist
    if not os.path.exists(outfolder):
        os.makedirs(outfolder)
        print(f"Created output directory: {outfolder}")
    
    basename = os.path.basename(os.path.splitext(bcfile)[0])
    outname = os.path.join(outfolder, f"{prefix}{basename}.bam")
    
    # Create a temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Read barcodes
        bc = read_barcodes(bcfile)
        if not bc:
            print("No valid barcodes found. Exiting.")
            return

        print(f"Processing BAM file: {bamfile}")
        print(f"Output file will be: {outname}")
        
        # Open input BAM file
        with pysam.AlignmentFile(bamfile, "rb") as bam:
            try:
                chrom_length = bam.get_reference_length(str(mtchr))
                print(f"Chromosome {mtchr} length: {chrom_length}")
            except ValueError:
                print(f"Chromosome {mtchr} not found in BAM file")
                return
                
            # Define chunks
            chunks = [(i, min(i + chunk_size, chrom_length)) 
                    for i in range(0, chrom_length, chunk_size)]
            print(f"Processing in {len(chunks)} chunks")
            
            # Process each chunk and write matching reads to separate temp files
            temp_files = []
            total_matching_reads = 0
            
            for i, (chunk_start, chunk_end) in enumerate(chunks):
                chunk_matching_reads = 0
                temp_file = os.path.join(temp_dir, f"chunk_{i}.bam")
                temp_files.append(temp_file)
                
                # Process this chunk
                print(f"Processing chunk {i+1}/{len(chunks)}: {chunk_start}-{chunk_end}")
                
                with pysam.AlignmentFile(temp_file, "wb", header=bam.header) as chunk_out:
                    for read in bam.fetch(str(mtchr), start=chunk_start, end=chunk_end, multiple_iterators=True):
                        try:
                            barcode_id = read.get_tag(barcodeTag)
                            if barcode_id in bc:
                                chunk_out.write(read)
                                chunk_matching_reads += 1
                        except KeyError:
                            continue  # Skip reads without the barcode tag
                
                total_matching_reads += chunk_matching_reads
                print(f"Chunk {i+1}/{len(chunks)} complete, found {chunk_matching_reads} matching reads")
            
            print(f"Total matching reads: {total_matching_reads}")
        
        # Sort and merge all temp files into the final output
        if temp_files:
            print("Merging and sorting chunks...")
            
            # Two approaches to merge:
            # 1. Using pysam.merge (clean but might have issues with very large datasets)
            try:
                pysam.merge("-f", outname, *temp_files)
                print("Successfully merged chunks using pysam.merge")
            except Exception as e:
                print(f"Error with pysam.merge: {e}")
                print("Falling back to samtools sort...")
                
                # 2. Alternative: Using samtools sort directly
                try:
                    # First concatenate all BAMs
                    concat_bam = os.path.join(temp_dir, "concat.bam")
                    with pysam.AlignmentFile(concat_bam, "wb", header=pysam.AlignmentFile(temp_files[0]).header) as out:
                        for temp_file in temp_files:
                            with pysam.AlignmentFile(temp_file) as infile:
                                for read in infile:
                                    out.write(read)
                    
                    # Then sort the concatenated BAM
                    subprocess.run(["samtools", "sort", "-o", outname, concat_bam], check=True)
                    print("Successfully merged and sorted using samtools sort")
                except Exception as e:
                    print(f"Error with samtools sort fallback: {e}")
                    return
        else:
            print("No matching reads found. No output file created.")
            return
    
    # Index the output BAM file
    print(f"Indexing {outname}")
    try:
        pysam.index(outname)
        print("Indexing complete")
    except Exception as e:
        print(f"Error indexing BAM file: {e}")
        print("Trying to sort the BAM file again before indexing...")
        try:
            # Create a temporary sorted BAM
            sorted_bam = outname + ".sorted.bam"
            pysam.sort("-o", sorted_bam, outname)
            # Replace the original with the sorted one
            os.replace(sorted_bam, outname)
            # Try indexing again
            pysam.index(outname)
            print("Indexing complete after re-sorting")
        except Exception as e2:
            print(f"Final error indexing BAM file: {e2}")
    
    elapsed_time = time.time() - start_time
    print(f"Processing complete in {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python chunk_barcoded_bam.py bamfile outfolder barcodeTag bcfile mtchr [prefix]")
        sys.exit(1)
    bamfile = sys.argv[1]
    outfolder = sys.argv[2]
    barcodeTag = sys.argv[3]
    bcfile = sys.argv[4]
    mtchr = sys.argv[5]
    
    # Handle optional prefix parameter
    prefix = ""
    if len(sys.argv) >= 7:
        prefix = sys.argv[6]
    
    subset_bam(bamfile, outfolder, barcodeTag, bcfile, mtchr, prefix)