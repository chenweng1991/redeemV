import sys
import pysam
from collections import defaultdict


bam_file=sys.argv[1]
barcode_tag=sys.argv[2]
bam_input = pysam.AlignmentFile(bam_file, "rb")

ReadPairDict=defaultdict(list)
for read in bam_input:
    ReadPairDict[read.query_name].append(read)

print("Dictionary in...")

for read_name in ReadPairDict:
    # disregard singlets and multiplets
    if len(ReadPairDict[read_name]) != 2:
        continue

    # identify fwd and rev in a pair
    read0, read1 = ReadPairDict[read_name]
    if read0.is_reverse and not read1.is_reverse:
            fwd_read, rev_read = read1, read0
    elif not read0.is_reverse and read1.is_reverse:
            fwd_read, rev_read = read0, read1
    else:
            # disregard a pair if both are the same strand
        continue
    if fwd_read.has_tag(barcode_tag):
        CellBC=fwd_read.get_tag(barcode_tag)
        fwd_chr=fwd_read.reference_name
        rev_chr=rev_read.reference_name
        FragSize=abs(fwd_read.tlen)
        if fwd_chr==rev_chr and FragSize<5000:
            if(fwd_read.tlen>0):
                print(CellBC+"\t"+fwd_chr+"\t"+str(fwd_read.pos)+"\t"+str(fwd_read.pos+FragSize)+"\t"+read_name)
            else:
                print(CellBC+"\t"+fwd_chr+"\t"+str(fwd_read.pos+FragSize)+"\t"+str(fwd_read.pos)+read_name)
