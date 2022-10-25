import sys


Wlist=sys.argv[1]
WlistDic=[]

with open(Wlist) as WL:
    for line in WL:
        WlistDic.append(line.strip())

WlistDic=tuple(WlistDic)

f=sys.stdin
for line in f:
    Cur_line=line.strip().split()
    if (Cur_line[3] in WlistDic):
        print(Cur_line[0]+'\t'+Cur_line[1]+'\t'+Cur_line[2]+'\t'+Cur_line[3]+'\t'+Cur_line[4])
