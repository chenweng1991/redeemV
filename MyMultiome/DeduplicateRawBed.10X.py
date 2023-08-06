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
    with open(bed_filename) as f:
        line1 = f.readline().strip().split()
        CellBC, Chr, Start, End = line1[:4]
        PCR = 1
        
        with open(out_filename, "w") as f_cut_file:
            Container=[[Chr,Start,End,CellBC,PCR]]
            for line in f:
                Cur_CellBC, Cur_Chr, Cur_Start, Cur_End, =line.strip().split()[:4]
                Cur_PCR=1
                if Cur_CellBC==CellBC:
                    if Cur_Chr==Chr and Cur_Start==Start and Cur_End==End:
                        Container[-1][4]=Container[-1][4]+1
                    else:
                        Container.append([Cur_Chr,Cur_Start,Cur_End,Cur_CellBC,Cur_PCR])
                        CellBC=Cur_CellBC
                        Chr=Cur_Chr
                        Start=Cur_Start
                        End=Cur_End
                        PCR=1
                else:
                    if(len(Container)>Cutoff):
                        if not CellBC=="NA":
                            for fragment in Container:
                                f_cut_file.write("\t".join([str(i) for i in fragment])+"\n")
                    CellBC=Cur_CellBC
                    Container=[[Cur_Chr,Cur_Start,Cur_End,Cur_CellBC,Cur_PCR]]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a BED file based on unique fragment number cutoff for each cell.')
    parser.add_argument('bed_filename', type=str, help='Input BED filename')
    parser.add_argument('out_filename', type=str, help='Output filename')
    parser.add_argument('--cutoff', type=int, default=0, help='Number of unique fragment number for each cell, as cutoff. Default is 0.')

    args = parser.parse_args()

    process_bed_file(args.bed_filename, args.out_filename, args.cutoff)