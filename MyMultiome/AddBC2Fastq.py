'''
File: AddBC2Fastq.py
Author: Chen Weng
Date: 2023-07-30
Updated: 2023-07-30

Description:
    This script is designed to append barcodes from an i5 barcode fastq file to the names of sequences in paired-end FASTQ reads (R1 and R2).
    It accepts two input FASTQ files (or gzipped FASTQ files), an i5 barcode FASTQ file, and two output FASTQ files as command-line arguments.
    By extracting and transforming the barcodes, it creates two output files with the barcodes appended to the sequence identifiers.

Usage:
    python3 AddBC2Fastq.py <input1> <input2> <barcodes> <output1> <output2>
'''


import gzip
import argparse

def open_file(filepath, mode, buffering=-1):
    if filepath.endswith('.gz'):
        return gzip.open(filepath, mode)
    else:
        return open(filepath, mode, buffering=buffering)

TRANSLATION_TABLE = str.maketrans('ATCG', 'TAGC')
def reverse_complement(seq):
    return seq[::-1].translate(TRANSLATION_TABLE)

def write_fastq_file(in_file, out_file, barcode):
    read_line = next(in_file).strip().split()[0]
    out_file.write((read_line + '|' + barcode + '\n'))
    out_file.write(next(in_file))  # sequence
    out_file.write(next(in_file))  # separator
    out_file.write(next(in_file))  # quality

def append_barcodes(input1, input2, barcodes, output1, output2,buffering):
    with open_file(input1, 'rt',buffering = buffering) as in1, \
         open_file(input2, 'rt',buffering = buffering) as in2, \
         open_file(barcodes, 'rt',buffering = buffering) as bc_file, \
         open_file(output1, 'wt',buffering = buffering) as out1, \
         open_file(output2, 'wt',buffering = buffering) as out2:

        while True:
            # Read 4 lines from the barcode file
            L1 = bc_file.readline().strip()
            L2 = bc_file.readline().strip()
            L3 = bc_file.readline().strip()
            L4 = bc_file.readline().strip()

            if L1 == "":
                break

            ReadName = L1.split()[0]
            Sequence = L2
            CellBC = reverse_complement(Sequence[8:24])

            # Process input1 and input2
            write_fastq_file(in1, out1, CellBC)
            write_fastq_file(in2, out2, CellBC)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get Cell barcode from i5 and Append barcodes to FASTQ file names.')
    parser.add_argument('input1', help='Path to the R1 input FASTQ (or gzipped FASTQ) file.')
    parser.add_argument('input2', help='Path to the R2 input FASTQ (or gzipped FASTQ) file.')
    parser.add_argument('barcodes', help='Path to the i5 barcode FASTQ (or gzipped FASTQ) file.')
    parser.add_argument('output1', help='Path to the R1 output FASTQ file.')
    parser.add_argument('output2', help='Path to the R2 output FASTQ file.')
    parser.add_argument('--buffering', type=int, default=-1, help='Buffer size for reading/writing files (in bytes). Default is -1.')

    args = parser.parse_args()

    append_barcodes(args.input1, args.input2, args.barcodes, args.output1, args.output2,buffering=args.buffering)
