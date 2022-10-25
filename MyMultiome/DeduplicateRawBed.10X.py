import sys

Bedfilename=sys.argv[1]
Cutoff=int(sys.argv[2])  # # of unique fragment number for each cell. 1000 as cutoff
f=open(Bedfilename)
Line1=f.readline().strip().split()
CellBC=Line1[0]
Chr=Line1[1]
Start=Line1[2]
End=Line1[3]
PCR=1
fcutfile=open(Bedfilename.split(".")[0]+"."+Bedfilename.split(".")[1]+".fragment"+"."+str(Cutoff)+".cut.tsv","w")
#fallfile=open(Bedfilename.split(".")[0]+"."+Bedfilename.split(".")[1]+".fragment.tsv","w")

Container=[]
Container.append([Chr,Start,End,CellBC,PCR])
for line in f:
    Cur_CellBC=line.strip().split()[0]
    Cur_Chr=line.strip().split()[1]
    Cur_Start=line.strip().split()[2]
    Cur_End=line.strip().split()[3]
    Cur_PCR=1
    if Cur_CellBC==CellBC:
        if Cur_Chr==Chr and Cur_Start==Start and Cur_End==End:
            Container[-1][4]=Container[-1][4]+1
        else:
            Container.append([Cur_Chr,Cur_Start,Cur_End,Cur_CellBC,Cur_PCR])
            CellBC=Cur_CellBC
            Chr=Cur_Chr
            Start=Cur_Start
            End=Cur_End
            PCR=1
    else:
        #for fragment in Container:
            #fallfile.write("\t".join([str(i) for i in fragment])+"\n")
        if(len(Container)>Cutoff):
            if not CellBC=="NA":
                for fragment in Container:
                    fcutfile.write("\t".join([str(i) for i in fragment])+"\n")
        CellBC=Cur_CellBC
        Container=[]
        Container.append([Cur_Chr,Cur_Start,Cur_End,Cur_CellBC,Cur_PCR])

fcutfile.close()
#fallfile.close()
