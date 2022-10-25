import sys
import pysam
import gzip
from collections import defaultdict


bam_file=sys.argv[1]
barcode_tag=sys.argv[2]
# CellBC_file=sys.argv[3]
bam_input = pysam.AlignmentFile(bam_file, "rb")

### Make dictionary for real cells
# CellBC_dic=[]
# with gzip.open(CellBC_file) as f:
#     while True:
#         L=f.readline().strip()
#         CellBC_dic.append(L.decode())
#         if L==b"":
#             break


### tuple the cells
# CellBC_dic=tuple(CellBC_dic)

##
ReadPairDict=defaultdict(list)
for read in bam_input:
    if read.has_tag(barcode_tag):
        ReadPairDict[read.query_name].append(read)
        

##
for read_name in ReadPairDict:
    # disregard singlets and multiplets
    if len(ReadPairDict[read_name]) != 2:
        continue
    else:
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
            if(fwd_read.tlen>0):
                print(CellBC+'\t'+fwd_read.reference_name+"\t"+str(fwd_read.pos)+'\t'+str(fwd_read.pos+fwd_read.tlen))  #fwd_read.qname
            else:
                print(CellBC+'\t'+fwd_read.reference_name+"\t"+str(fwd_read.pos-fwd_read.tlen)+"\t"+str(fwd_read.pos))

                
                