import sys
import os
f=sys.stdin

FirstlineAy=f.readline().strip().split()
CellBC=FirstlineAy[3]
SingleCellBed=[]
SingleCellBed.append("\t".join(FirstlineAy[0:3]))
out=sys.argv[1]
print(out)
i=0
for line in f:
    lineAy=line.strip().split()
    Cur_CellBC=lineAy[3]
    if Cur_CellBC==CellBC:
        SingleCellBed.append("\t".join(lineAy[0:3]))
    else:
        i+=1
        print(i)
        fout=open("tmp.bed","w")
        fout.writelines("%s\n" % l for l in SingleCellBed)
        fout.close()
        CellBCs=[CellBC]*16569
        fnameout=open("tmp.name.bed","w")
        fnameout.writelines("%s\n" % l for l in CellBCs)
        fnameout.close()
        os.system("bedtools coverage -a /lab/solexa_weissman/cweng/Genomes/mitoGenome.bed -b tmp.bed -d | cut -f4,5 | paste - tmp.name.bed  >>"+out)
        CellBC=Cur_CellBC
        SingleCellBed=[]
        SingleCellBed.append("\t".join(lineAy[0:3]))
