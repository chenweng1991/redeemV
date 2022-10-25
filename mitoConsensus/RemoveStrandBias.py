import sys
out_dir=sys.argv[1]

StrandBiaseBlackList_file=out_dir+"/final/StrandBiaseBlackList"
Total=out_dir+"/final/RawGenotypes.Total"
VerySensitive=out_dir+"/final/RawGenotypes.VerySensitive"
Sensitive=out_dir+"/final/RawGenotypes.Sensitive"
Specific=out_dir+"/final/RawGenotypes.Specific"

Total_o_file=Total+".StrandBalance"
VerySensitive_o_file=VerySensitive+".StrandBalance"
Sensitive_o_file=Sensitive+".StrandBalance"
Specific_o_file=Specific+".StrandBalance"

StrandBiasDic=[]
with open(StrandBiaseBlackList_file) as f:
    for line in f:
        StrandBiasDic.append(line.strip())

StrandBiasDic=tuple(StrandBiasDic)

##Open output files
Total_o=open(Total_o_file,"w")
VerySensitive_o=open(VerySensitive_o_file,"w")
Sensitive_o=open(Sensitive_o_file,"w")
Specific_o=open(Specific_o_file,"w")

##
with open(Total) as f:
    for line in f:
        content=line.strip().split()
        V=content[3]
        DB=float(content[9])
        # print(V)
        # print(DB)
        if V in StrandBiasDic and DB==0:
            continue
        else:
            Total_o.write(line)
Total_o.close()

##
with open(VerySensitive) as f:
    for line in f:
        content=line.strip().split()
        V=content[3]
        DB=float(content[9])
        # print(V)
        # print(DB)
        if V in StrandBiasDic and DB==0:
            continue
        else:
            VerySensitive_o.write(line)
Total_o.close()


##
with open(Sensitive) as f:
    for line in f:
        content=line.strip().split()
        V=content[3]
        DB=float(content[9])
        # print(V)
        # print(DB)
        if V in StrandBiasDic and DB==0:
            continue
        else:
            Sensitive_o.write(line)
Total_o.close()


##
with open(Specific) as f:
    for line in f:
        content=line.strip().split()
        V=content[3]
        DB=float(content[9])
        # print(V)
        # print(DB)
        if V in StrandBiasDic and DB==0:
            continue
        else:
            Specific_o.write(line)
Total_o.close()
