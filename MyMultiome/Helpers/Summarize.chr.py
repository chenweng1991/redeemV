'''Usage
This is a simple Python for chrY loss analysis
It take a monoclonal bed-like files, sorted by both chrmosome(column1) and cell name(column4)
It count the uniq reads on each chromosome for each cell and out put a long matrix

Note, to prepare the monoclonal bed-like files, simply do bash "cat atac_fragments.tsv.gz | grep -v '^#' | sort -k1,1 -k4,4 > atac_fragments.tsv.sorted "

Example: 
path=/lab/solexa_weissman/cweng/Packages/MyMultiome/Helpers/
python $path/Summarize.chr.py atac_fragments.tsv.sorted filtered_feature_bc_matrix/barcodes.tsv.gz > atac_chr_summarise

'''

import sys
import gzip
Monoclonal=sys.argv[1]
ValidCellBC=sys.argv[2]

ValidCellTupple=[]
with gzip.open(ValidCellBC,mode="rb") as f:
    while True:
        BC=f.readline().strip()
        if BC==b"":
            break
        else:
            ValidCellTupple.append(BC.decode())

ValidCellTupple=tuple(ValidCellTupple)


with open(Monoclonal) as f:
    Line=f.readline().strip()
    LineArray=Line.split()
    Group_chr=LineArray[0]
    Group_Cell=LineArray[3]
    Count=1
    while True:
        # print(i)
        # i=i+1
        Line=f.readline().strip()
        if Line=="":
            break
        else:
            LineArray=Line.split()
            chr=LineArray[0]
            Cell=LineArray[3]
            if chr==Group_chr and Cell==Group_Cell:
                Count=Count+1
            else:
                if Group_Cell in ValidCellTupple:                    
                    print(Group_Cell+"\t"+Group_chr+"\t"+str(Count)+"\t")
                Count=1
                Group_chr=LineArray[0]
                Group_Cell=LineArray[3]


