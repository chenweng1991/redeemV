import sys
out_dir=sys.argv[1]

# out_dir="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor01_BMMC2_ATACkit/MTenrich/CW_mgatk/"

Total={}
VerySensitive={}
Sensitive={}
Specific={}

## Read in the genotypes as
print("Start...")
with open(out_dir+"/final/RawGenotypes.Total") as f:
    for line in f:
        content=line.strip().split()
        CellPos=content[1]+content[2]
        if CellPos in Total.keys():
            Total[CellPos].append(line.strip())
        else:
            Total[CellPos]=[line.strip()]

print("Dic1 In")

with open(out_dir+"/final/RawGenotypes.VerySensitive") as f:
    for line in f:
        content=line.strip().split()
        CellPos=content[1]+content[2]
        if CellPos in VerySensitive.keys():
            VerySensitive[CellPos].append(line.strip())
        else:
            VerySensitive[CellPos]=[line.strip()]


print("Dic2 In")

with open(out_dir+"/final/RawGenotypes.Sensitive") as f:
    for line in f:
        content=line.strip().split()
        CellPos=content[1]+content[2]
        if CellPos in Sensitive.keys():
            Sensitive[CellPos].append(line.strip())
        else:
            Sensitive[CellPos]=[line.strip()]


print("Dic3 In")

with open(out_dir+"/final/RawGenotypes.Specific") as f:
    for line in f:
        content=line.strip().split()
        CellPos=content[1]+content[2]
        if CellPos in Specific.keys():
            Specific[CellPos].append(line.strip())
        else:
            Specific[CellPos]=[line.strip()]


print("Dic4 In")

o_Total=open(out_dir+"/final/RawGenotypes.Total","w")
o_VerySensitive=open(out_dir+"/final/RawGenotypes.VerySensitive","w")
o_Sensitive=open(out_dir+"/final/RawGenotypes.Sensitive","w")
o_Specific=open(out_dir+"/final/RawGenotypes.Specific","w")

with open(out_dir+"/final/QualifiedTotalCts") as f:
    for line in f:
        content=line.strip().split()
        CellPos=content[0]+content[1]
        if CellPos in Total.keys():
            for molecule in Total[CellPos]:
                o_Total.write(molecule+"\t"+content[2]+"\n")
            if CellPos in VerySensitive.keys():
                for molecule in VerySensitive[CellPos]:
                    o_VerySensitive.write(molecule+"\t"+content[3]+"\n")
                if CellPos in Sensitive.keys():
                    for molecule in Sensitive[CellPos]:
                        o_Sensitive.write(molecule+"\t"+content[4]+"\n")
                    if CellPos in Specific.keys():
                        for molecule in Specific[CellPos]:
                            o_Specific.write(molecule+"\t"+content[5]+"\n")


o_Total.close()
o_VerySensitive.close()
o_Sensitive.close()
o_Specific.close()
