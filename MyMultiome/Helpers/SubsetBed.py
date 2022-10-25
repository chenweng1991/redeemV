'''
This small program read in the bed.gz stype file with column 4 cell barcode, this is file1
File 2 is a simple one column table that includes the cell barcodes of interest to extract from file1

The output is the extracted bed file
'''
import sys
import gzip

Monoclonal=sys.argv[1]
CellChoose=sys.argv[2]

CellNameList=list()
with open(CellChoose) as f:
    while True:
        Line=f.readline().strip()
        if Line=="":
            break
        else:
           CellNameList.append(Line)

CellNameList=tuple(CellNameList)


with gzip.open(Monoclonal) as f:
    while True:
        Line=f.readline().strip()
        if Line==b"":
            break
        else:
            LineArray=Line.split()
            chr=LineArray[0]
            start=LineArray[1]
            end=LineArray[2]
            cell=LineArray[3]
            Count=LineArray[4]
            if cell.decode() in CellNameList:
                print(chr.decode()+"\t"+start.decode()+"\t"+end.decode()+"\t"+cell.decode()+"\t"+Count.decode())