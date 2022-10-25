'''
This small program read in the monoclonal fragment bed.gz file with column 4 cell barcode, this is file1 
File 2 is a simple two column table that includes the cell barcodes and cell types

The outputs are seperated bed files 
'''
import sys
import gzip

Monoclonal=sys.argv[1]
CellTypeFile=sys.argv[2]

'''Define RegChr'''
RegChr=("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")

'''Make CellTypeDict, a dictionary from cell name to cell type'''
CellTypeDict={}
with open(CellTypeFile) as f:
    for Line in f:
        LineArray=Line.strip().split()
        CellTypeDict[LineArray[0]]=LineArray[1]
f.close()

'''Make CellType2File, a dictionary from cell type to file handler'''
UniqValues=set(CellTypeDict.values())
CellType2File={}
for V in UniqValues:
    CellType2File[V]=open(V+".monoclonal.bed","w+")

'''Make CellName2File, a dictionary from cellname to file handler'''
CellName2File={}
with open(CellTypeFile) as f:
    for Line in f:
        LineArray=Line.strip().split()
        CellName2File[LineArray[0]]=CellType2File[LineArray[1]]
f.close()

print("Dictionarys are in ...")
print("Start spliting...")

with gzip.open(Monoclonal, mode="rt") as f:
    for Line in f:
        LineArray=Line.strip().split()
        chr=LineArray[0]
        start=LineArray[1]
        end=LineArray[2]
        cell=LineArray[3]
        Count=LineArray[4]
        if cell in CellName2File.keys():
            if chr in RegChr:
            #print(cell+"is"+CellTypeDict[cell])
                CellName2File[cell].write(chr+"\t"+start+"\t"+end+"\t"+cell+"\t"+Count+"\n")
        #else:
            #print(cell+" not in dictionary")

for V in CellType2File.keys():
    CellType2File[V].close()