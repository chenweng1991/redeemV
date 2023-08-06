import sys

f=sys.stdin

FirstLine=f.readline().strip().split()
CellBC=FirstLine[3]
Counts=0
TotalReads=0
Counts+=1
TotalReads=TotalReads+int(FirstLine[4])
for line in f:
    Cur_line=line.strip().split()
    Cur_CellBC=Cur_line[3]
    if CellBC==Cur_CellBC:
        TotalReads=TotalReads+int(Cur_line[4])
        Counts+=1
    else:
        sys.stdout.write(CellBC+"\t"+str(Counts)+"\t"+str(TotalReads)+"\n")
        CellBC=Cur_CellBC
        Counts=0
        TotalReads=0
        Counts+=1
        TotalReads=TotalReads+int(Cur_line[4])
