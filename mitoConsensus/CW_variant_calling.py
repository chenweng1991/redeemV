import click
import sys
import glob
import pickle
import pandas as pd
import numpy as np

##Function define
##Simple function to change to upper case for a list
def upper_list(lst):
    New=[]
    for x in lst:
        if x.isupper():
            New.append(x)
        else:
            New.append(x.upper())
    return New


## organize the variant into dataframe for any given cell
def scOrganizeVariant(Genotypemtx, cellname, ref,refidx):
    dna_letters = ['A','C','G','T']
    variantslist=[]
    for idx, ACGT in enumerate(Genotypemtx):
        total=ACGT.sum()
        for baseidx in range(0,4):
            if not ACGT[baseidx]==0:
                variant=[cellname,str(idx+1)+ref[idx]+">"+dna_letters[baseidx],not refidx[idx]==baseidx,ACGT[baseidx],total]
                variantslist.append(variant)
    return variantslist

def CellFilteredDepth(Genotypemtx, cellname):
    df=pd.DataFrame({'Cell':[cellname]*16569, 'FilteredDepth':np.sum(Genotypemtx,axis=1)})
    return(df)
###############################################################################################

filename=sys.argv[1]

with open(filename,"rb") as fl:
    mtx=pickle.load(fl)
# Prepare mito reference
mito_ref_file ="CW_mgatk/final/chrM_refAllele.txt"
mito_ref=list(pd.read_table(mito_ref_file,names=["pos","base"])["base"])
mito_ref[mito_ref.index('N')]="A"   ## This is an ok solution for the 3106-N, maybe we can try better later
dna_letters = ['A','C','G','T']
mito_ref_idx=[dna_letters.index(base) for base in upper_list(mito_ref)]

 #"temp/sparse_matrices/barcodes.1.supermtx.data"




variants_df_refrm_header=pd.DataFrame(columns=["CellName","Variant","IsMutant","Counts","Total"])
FilteredDepthSummary_header=pd.DataFrame(columns=["Position","CellName","FilteredDepth"])
variants_df_refrm_header.to_csv(filename.replace(".supermtx.data",".Variants.csv"),index=False,header=True)
FilteredDepthSummary_header.to_csv(filename.replace(".supermtx.data",".FilteredDepthSummary.csv"),index=False,header=True)
for cell in mtx.keys():
    print(cell)
    variants_df=pd.DataFrame(scOrganizeVariant(mtx[cell], cell, mito_ref,mito_ref_idx),columns=["CellName","Variant","IsMutant","Counts","Total"])
    variants_df_refrm=variants_df[variants_df["IsMutant"]]
    FilteredDepthSummary=CellFilteredDepth(mtx[cell], cell)
    FilteredDepthSummary.insert(0,"Position",FilteredDepthSummary.index.values+1)
    variants_df_refrm.to_csv(filename.replace(".supermtx.data",".Variants.csv"),index=False,mode='a', header=False)
    FilteredDepthSummary.to_csv(filename.replace(".supermtx.data",".FilteredDepthSummary.csv"),index=False,mode='a', header=False)
