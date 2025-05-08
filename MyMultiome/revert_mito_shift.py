import pysam
import argparse
import os 
def convert_position(pos, shift, genome_size): 
    orig_pos = (pos + shift) % genome_size
    if orig_pos == 0:
        orig_pos = genome_size 
    return orig_pos

def convert_bam_position(input_bam, output_bam, shift=8000, genome_size=16569):
    in_bam = pysam.AlignmentFile(input_bam, "rb")
    out_bam = pysam.AlignmentFile(output_bam, "wb", header=in_bam.header)

    read_count = 0
    for read in in_bam:
        read_count += 1

        if read.is_unmapped or read.mate_is_unmapped:
            continue

        # Shift read position
        read_start_1b = read.reference_start + 1
        shifted_start_1b = convert_position(read_start_1b, shift, genome_size)
        read.reference_start = shifted_start_1b - 1

        # Shift mate position
        mate_start_1b = read.next_reference_start + 1
        shifted_mate_1b = convert_position(mate_start_1b, shift, genome_size)
        read.next_reference_start = shifted_mate_1b - 1

        # Compute TLEN 
        if read.reference_start <= read.next_reference_start:
            tlen = (read.next_reference_start + read.query_length) - read.reference_start
        else:
            tlen = -((read.reference_start + read.query_length) - read.next_reference_start)

        read.template_length = tlen

        if read_count <= 10 or read_count % 10000 == 0:
            print(f"Read {read.query_name}: {read_start_1b} -> {shifted_start_1b}")

        out_bam.write(read)

    in_bam.close()
    out_bam.close()

    print(f"Processed {read_count} reads")
    print(f"Indexing {output_bam}...")
    tmp_sorted_bam = output_bam + ".sorted"
    pysam.sort("-o", tmp_sorted_bam, output_bam)
    os.replace(tmp_sorted_bam, output_bam)
    pysam.index(output_bam)
    print(f"Results written to {output_bam}")
    print("Warning: This is a simplified conversion. For proper BAM conversion,")
    print("consider additional factors like read pairs, CIGAR operations, etc.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Shift read positions in a BAM file for circular genomes.")
    parser.add_argument("input_bam", help="Input BAM file")
    parser.add_argument("output_bam", help="Output BAM file name")
    parser.add_argument("--shift", type=int, default=8000, help="Shift amount (default: 8000)")
    parser.add_argument("--genome_size", type=int, default=16569, help="Genome size (default: 16569)")
    args = parser.parse_args()

    convert_bam_position(args.input_bam, args.output_bam, shift=args.shift, genome_size=args.genome_size)
