''' This small program take sorted two column data in, if duplicated, only leave one and add a third column indicating the numbers'''

import sys
file=sys.argv[1]


f=open(file)

FirstLineArray=f.readline().strip().split()
PreviousCell=FirstLineArray[0]
PreviousPeak=FirstLineArray[1]
Count=1
for Line in f:
    #LineArray=f.readline().strip().split()
    LineArray=Line.strip().split()
    Cur_Cell=LineArray[0]
    Cur_Peak=LineArray[1]
    if (Cur_Cell==PreviousCell and Cur_Peak==PreviousPeak):
        Count=Count+1
    else:
        print(PreviousCell+"\t"+PreviousPeak+"\t"+str(Count))
        PreviousCell=Cur_Cell
        PreviousPeak=Cur_Peak
        Count=1
