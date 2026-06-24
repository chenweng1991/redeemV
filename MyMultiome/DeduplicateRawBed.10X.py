'''
File: DeduplicateRawBed.10X.py
Author: Chen Weng
Date: 2023-08-04
Updated: 2023-08-04

Description:
    This script is specifically designed for the deduplication of 10X raw BED files. It processes
    the BED files by grouping entries based on unique fragment numbers for each cell. After deduplication,
    only cells with fragment counts surpassing a specified cutoff (Default is 0) are retained and written to the output.

Usage:
    python3 DeduplicateRawBed.10X.py <input_bed> <output_bed> --cutoff <cutoff_value>
'''

import argparse

def process_bed_file(bed_filename, out_filename, Cutoff):
    def write_cell(container, cell_bc, handle):
        if cell_bc != "NA" and len(container) > Cutoff:
            for fragment in container:
                handle.write("\t".join([str(i) for i in fragment]) + "\n")

    with open(bed_filename) as f, open(out_filename, "w") as f_cut_file:
        first = f.readline().strip().split()
        if not first:
            return

        cell_bc, chrom, start, end = first[:4]
        container = [[chrom, start, end, cell_bc, 1]]

        for line in f:
            cur_cell_bc, cur_chrom, cur_start, cur_end = line.strip().split()[:4]
            if cur_cell_bc != cell_bc:
                write_cell(container, cell_bc, f_cut_file)
                cell_bc, chrom, start, end = cur_cell_bc, cur_chrom, cur_start, cur_end
                container = [[cur_chrom, cur_start, cur_end, cur_cell_bc, 1]]
                continue

            if cur_chrom == chrom and cur_start == start and cur_end == end:
                container[-1][4] += 1
            else:
                chrom, start, end = cur_chrom, cur_start, cur_end
                container.append([cur_chrom, cur_start, cur_end, cur_cell_bc, 1])

        write_cell(container, cell_bc, f_cut_file)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a BED file based on unique fragment number cutoff for each cell.')
    parser.add_argument('bed_filename', type=str, help='Input BED filename')
    parser.add_argument('out_filename', type=str, help='Output filename')
    parser.add_argument('--cutoff', type=int, default=0, help='Number of unique fragment number for each cell, as cutoff. Default is 0.')

    args = parser.parse_args()

    process_bed_file(args.bed_filename, args.out_filename, args.cutoff)