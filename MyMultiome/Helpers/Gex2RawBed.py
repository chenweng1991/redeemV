import sys
import pysam

bam_file=sys.argv[1]

barcode_tag="CB"
UMI_tag="UB"
Gene_tag="GN"

bam_input = pysam.AlignmentFile(bam_file, "rb")

for read in bam_input:
    if read.has_tag(barcode_tag) and read.has_tag(Gene_tag) and read.has_tag(UMI_tag):
        print(read.reference_name+"\t"+str(read.reference_start)+"\t"+str(read.reference_end)+"\t"+read.get_tag(barcode_tag)+"\t"+read.get_tag(UMI_tag)+"\t"+read.get_tag(Gene_tag))
        ##print(read.get_tag(barcode_tag)+"\t"+read.get_tag(Gene_tag)+"\t"+read.get_tag(UMI_tag)+"\t"+read.reference_name+"\t"+str(read.reference_start)+"\t"+str(read.reference_end))
