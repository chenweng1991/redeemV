import pysam
import sys
import sqlite3
import argparse
import subprocess
import os
from concurrent.futures import ProcessPoolExecutor


def process_chromosome(chromosome, Bamfilename, DB_path, outSamFileName):
    Bam = pysam.AlignmentFile(Bamfilename, "rb")
    # Connect to SQLite database
    conn = sqlite3.connect(DB_path)
    cursor = conn.cursor()

    # Open the output SAM file
    outSam = pysam.AlignmentFile(outSamFileName, "wb", template=Bam)

    results = []
    i=0
    for read in Bam.fetch(contig=chromosome):
        Qname = "@"+read.query_name

        # Debug
        i=i+1
        print(i)

        # Query the database for the value associated with Qname
        # cursor.execute("SELECT value FROM kv WHERE key = ?", (Qname,))
        # result = cursor.fetchone()

        # if result:
        #     print(result[0])
        #     read.set_tag('BC', result[0])
        #     outSam.write(read)  # Write the modified read directly to the output file

    # Close database connection
    conn.close()
    Bam.close()

    return results

def merge_and_remove_temp_files(temp_bam_files, merged_bam_file, num_cores):
    # Merge the temporary BAM files using samtools
    command = ['samtools', 'merge', '-@', str(num_cores), merged_bam_file] + temp_bam_files
    subprocess.run(command, check=True)

    # Sort the merged BAM file
    sorted_bam_file = merged_bam_file.replace('.bam', '.sorted.bam')
    command = ['samtools', 'sort', '-o', sorted_bam_file, '-@', str(num_cores), merged_bam_file]
    subprocess.run(command, check=True)

    # Remove the temporary BAM files
    for temp_file in temp_bam_files:
        os.remove(temp_file)

    # Optional: Remove the unsorted merged BAM file if not needed
    os.remove(merged_bam_file)
    



def main(DB_path, Bamfilename, OutPre, num_cores):
    TPLT = pysam.AlignmentFile(Bamfilename, "rb")
    Bam = pysam.AlignmentFile(Bamfilename, "rb")
    chromosomes = Bam.references

    temp_bam_files = []
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = []
        for chromosome in chromosomes:
            out_filename = f"{OutPre}_{chromosome}.tmp.bam"
            temp_bam_files.append(out_filename)
            futures.append(executor.submit(process_chromosome, chromosome, Bamfilename, DB_path, out_filename))

        # Wait for all futures to complete
        for future in futures:
            future.result()

    # Merge and remove temporary BAM files
    merge_and_remove_temp_files(temp_bam_files, OutPre + ".tagged.bam")

    TPLT.close()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge Barcodes into BAM file using SQLite database.")
    parser.add_argument("DB_path", help="Path to the SQLite database.")
    parser.add_argument("Bam", help="Path to the BAM file.")
    parser.add_argument("output_prefix", help="Output prefix for the tagged BAM file.")
    parser.add_argument("--num_cores", type=int, default=16, help="Number of cores to use for processing. Default is 16")
    
    args = parser.parse_args()
    main(args.DB_path, args.Bam, args.output_prefix, args.num_cores)