import click
import glob
import pickle
import pandas as pd
from progress.bar import Bar

@click.command()
@click.version_option()
@click.option('--out-dir', '-o', default = ".", required=True, help='The folder of the whole analysis. REQUIRED.')

def main(out_dir):
    # Prepare mito reference
    mito_ref_file =out_dir + "/final/chrM_refAllele.txt"
    mito_ref=list(pd.read_table(mito_ref_file,names=["pos","base"])["base"])
    mito_ref[mito_ref.index('N')]="A"   ## This is an ok solution for the 3106-N, maybe we can try better later
    dna_letters = ['A','C','G','T']
    mito_ref_idx=[dna_letters.index(base) for base in upper_list(mito_ref)]

    ## Read in the intermediate data processed by CW_sumstatsBP.by
    mtx_dir=out_dir + "/temp/sparse_matrices/"
    supermtx_files=glob.glob(mtx_dir+"*[0-9].supermtx.data")   # we can also read in the qc file by qcmtx_files=glob.glob(mtx_dir+"*[0-9].QC.data")
    supermtx_full=dict()
    for filename in supermtx_files:
        with open(filename,"rb") as fl:
            mtx=pickle.load(fl)
            supermtx_full.update(mtx)

    ##
    Variants_allcell=[]
    bar = Bar('Genotyping', max=len(supermtx_full))
    for cell in supermtx_full.keys():
        Variants_allcell=Variants_allcell+scOrganizeVariant(supermtx_full[cell], cell, mito_ref,mito_ref_idx)
        bar.next()

    bar.finish()

    variants_df=pd.DataFrame(Variants_allcell,columns=["CellName","Variant","IsMutant","Counts","Total"])
    variants_df_refrm=variants_df[variants_df["IsMutant"]]  # Remove the WT genotypes
    # optional, calculate heteroplasmy variants_df_sparse["Heteroplasmy"]=variants_df_refrm["Counts"]/variants_df_refrm["Total"] # Add a colum of heteroplasmy
    # Write to the final folder
    variants_df_refrm.to_csv(out_dir+'/final/Variants.csv')


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
    bar = Bar('Process', max=len(Genotypemtx))
    for idx, ACGT in enumerate(Genotypemtx):
        total=ACGT.sum()
        for baseidx in range(0,4):
            if not ACGT[baseidx]==0:
                variant=[cellname,str(idx+1)+ref[idx]+">"+dna_letters[baseidx],not refidx[idx]==baseidx,ACGT[baseidx],total]
                variantslist.append(variant)
    return variantslist


if __name__ == '__main__':
    main()
