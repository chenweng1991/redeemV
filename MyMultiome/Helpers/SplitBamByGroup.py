
#############
################   DO NOT USE, THERE ARE SOME bugs not fixed yet, but should be fixed

##########################
import pysam
import sys
import numpy as np
import pandas as pd

HashDic_file=sys.argv[1]
Hashcall_file=sys.argv[2]
# LibraryType="RNA"   ## The CB tag in both ATAC and RNA bams are RNA barcode
BamInput=sys.argv[3]
tag="CB"
SigCut=0.6


HashDic={}
f=open(HashDic_file)
for line in f:
    L=line.strip().split(",")
    HashDic[L[0]]=L[1]
f.close()

## Make a dictionary to convert from RNA names to ATAC names(Not in use in this script, but leave it here for other cases)
ATACWhite=pd.read_table("/lab/solexa_weissman/cweng/Genomes/10X/WhiteList_10X_Multiome.ATAC",names=["ATACName"])
RNAWhite=pd.read_table("/lab/solexa_weissman/cweng/Genomes/10X/WhiteList_10X_Multiome.RNA",names=["RNAName"])
RNAWhite["ATAC_Name"]=ATACWhite.ATACName
RNA2ATACDic=RNAWhite.set_index("RNAName")["ATAC_Name"].to_dict()
# ATAC2RNADic=RNAWhite.set_index("ATAC_Name")["RNA_Name"].to_dict()

Hashcall=pd.read_table(Hashcall_file)
Hashcall["ATACName"]=[RNA2ATACDic[n] for n in Hashcall.cell]

## Make HashCall_RNA_dict convert RNA cell names into group bam file;
## Make HashCall_ATAC_dict convert ATAC cell names into group bam file;
Hashcall_sig=Hashcall[Hashcall["sig"]>SigCut]
Hashcall_sig["group"]=[HashDic[x] for x in Hashcall_sig["call"]]
Hashcall_sig.cell=Hashcall_sig.cell+"-1"
Hashcall_sig.ATACName=Hashcall_sig.ATACName+"-1"


# Define Bam files
# BamInput="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor04_BMMC1_ConditionTest/CellRanger/Donor04_BMMC1_GDN/outs/atac_possorted_bam.bam"
BamInput="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/Donor04_BMMC1_ConditionTest/CellRanger/Donor04_BMMC1_GDN/outs/gex_possorted_bam.bam"

out_handles_Dic={}
TPLT= pysam.AlignmentFile(BamInput, "rb")

for group in set(Hashcall_sig.group):
    Bam=pysam.AlignmentFile("Post10x."+group+".bam","wb",template=TPLT)
    out_handles_Dic[group]=Bam

Hashcall_sig["output"]=[out_handles_Dic[n] for n in Hashcall_sig.group]
HashCall_RNA2out_dict=Hashcall_sig[["cell","output"]].set_index("cell")["output"].to_dict()
HashCall_RNA2group_dict=Hashcall_sig[["cell","group"]].set_index("cell")["group"].to_dict()
# HashCall_ATAC_dict=Hashcall_sig[["ATACName","output"]].set_index("ATACName")["output"].to_dict()
## The HashCall_ATAC_dict is not in use becasue the CB tag in ATAC bam is also gex(rna) barcode


## ReadIn the bam

Bam=pysam.AlignmentFile(BamInput,"rb")
for read in Bam.fetch(until_eof=True):
    try:
        CB=read.get_tag('CB')
    except KeyError:
        continue
    if CB in HashCall_RNA2out_dict.keys():
        # print(CB)
        read.set_tag('GP',HashCall_RNA2group_dict[CB])
        HashCall_RNA2out_dict[CB].write(read)


print("BamSpliting successfully completed!")

with open("SplitLog","w") as f:
    f.write("BamSpliting successfully completed!\n")

f.close()
