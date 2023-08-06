'''
File: AddBC2BAM.py
Author: Chen Weng
Date: 2023-08-04
Updated: 2023-08-04

Description:
    This script is designed to extract cell barcodes from read names in a BAM file and 
    add them as a new tag ("BC"). It reads an input BAM file, extracts the cell barcode 
    from each read name, adds the cell barcode as a new tag to the read, and writes the 
    read to an output BAM file. 

Usage:
    python3 AddBC2BAM.py <input_bam> <output_bam>
'''

import pysam
import argparse

def add_cell_barcode_to_bam(input_bam_file, output_bam_file):
    # Open input and output BAM files
    with pysam.AlignmentFile(input_bam_file, "rb") as infile, \
        pysam.AlignmentFile(output_bam_file, "wb", header=infile.header) as outfile:

        # Iterate over each read in the BAM file
        for read in infile:
            # Split the read name by "|", and take the first and second elements
            read_name_parts = read.query_name.split("|")
            #new_read_name = read_name_parts[0]
            cell_barcode = read_name_parts[1]

            # Add the Cell Barcode as a new tag ("BC") to the read
            #read.query_name = new_read_name
            read.set_tag("BC", cell_barcode)

            # Write the read to the output BAM file
            outfile.write(read)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Add Cell barcode as a new tag to a BAM file.')
    parser.add_argument('input_bam', help='Path to the input BAM file.')
    parser.add_argument('output_bam', help='Path to the output BAM file.')

    args = parser.parse_args()

    add_cell_barcode_to_bam(args.input_bam, args.output_bam)