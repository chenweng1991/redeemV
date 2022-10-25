import pysam
import click
import os
import sys
import pickle
import numpy as np
import pandas as pd
from progress.bar import Bar
from collections import defaultdict


########## Section 1, define argument and constant and read files
##########

dna_letters = ['A','C','G','T','N']
##Take argument from command line
barcodes_file=sys.argv[1]
bam_file=sys.argv[2]
sample=os.path.basename(barcodes_file).replace(".txt","")  #The prefix, should be barcodes.# check the prefix of /temp/barcoded_bams/ for example. REQUIRED
out_dir=sys.argv[3] #The folder of the whole analysis. It is usually called CW_mgatk, under which you could find temp, final, etc REQUIRED.
barcode_tag=sys.argv[4]   ## It is usually "BC"
BaseQ_thld_hi=int(sys.argv[5]) ## It is usually set as 30


# barcodes_file =out_dir + "/temp/barcode_files/"+sample+".txt"  ## For example "temp/barcode_files/barcodes.1.txt"
# bam_file =out_dir + "/temp/barcoded_bams/"+sample+".bam"  ## For example "temp/barcoded_bams/barcodes.1.bam"
mito_ref_file =out_dir + "/final/chrM_refAllele.txt"
out_genotypeTotal_file=out_dir + "/temp/sparse_matrices2.0/"+sample+".RawGenotypes.Total"
out_genotypeVerySensitive_file=out_dir + "/temp/sparse_matrices2.0/"+sample+".RawGenotypes.VerySensitive"
out_genotypeSensitive_file=out_dir + "/temp/sparse_matrices2.0/"+sample+".RawGenotypes.Sensitive"
out_genotypeSpecific_file=out_dir + "/temp/sparse_matrices2.0/"+sample+".RawGenotypes.Specific"
out_totalCts_file=out_dir + "/temp/sparse_matrices2.0/"+sample+".QualifiedTotalCts"

if not 'sparse_matrices2.0' in os.listdir(out_dir +"/temp"):
    os.mkdir(out_dir+"/temp/sparse_matrices2.0")


#Read in the mito reference
mito_ref=pd.read_table(mito_ref_file,names=["pos","base"])
max_bp=mito_ref.shape[0]
# Import all cell barcodes
with open(barcodes_file) as barcode_file_handle:
    content = barcode_file_handle.readlines()


bcs = [x.strip() for x in content]


########## Section 2 Read bam file make a dictionary from Rname->read1,read2   Basically stick each Read pair, and search by ReadNames
##########
bam_input = pysam.AlignmentFile(bam_file, "rb")
MoleculeDict=defaultdict(list)
ReadPairDict=defaultdict(list)
for read in bam_input:
    ReadPairDict[read.query_name].append(read)

# Make dictionary from Molecule(Cell_Start_End)-->Rname
for read_name in ReadPairDict:
    # disregard singlets and multiplets
    if len(ReadPairDict[read_name]) != 2:
        continue
    # identify fwd and rev in a pair
    read0, read1 = ReadPairDict[read_name]
    if read0.is_reverse and not read1.is_reverse:
            fwd_read, rev_read = read1, read0
    elif not read0.is_reverse and read1.is_reverse:
            fwd_read, rev_read = read0, read1
    else:
            # disregard a pair if both are the same strand
        continue
    if fwd_read.has_tag(barcode_tag):
        CellBC=fwd_read.get_tag(barcode_tag)
        if(fwd_read.tlen>0):
            Molecule=CellBC+'_'+str(fwd_read.pos)+'_'+str(fwd_read.pos+fwd_read.tlen)
        else:
            Molecule=CellBC+'_'+str(fwd_read.pos-fwd_read.tlen)+"_"+str(fwd_read.pos)
        MoleculeDict[Molecule].append(read_name)

print("Dictionaries are in...")
########## Section 3 Main function, loop through each molecule for genotyping
########## Four .RawGenotypes outputs + one .QualifiedTotalCts with 4 columns output
########## Total: Without any consensus level filtering
########## Very_sensitive a=2 b=1 c= 0.75
########## Sensitive a=3 b=2 c= 0.75
########## Specific a=4 b=3 c= 0.9
TotalMoleculeCtsMatrix={key:np.zeros((max_bp,4)) for key in bcs} # Matrix to collect 1,Total;Very_sensitive;Sensitive;Specific

## Open the 4 files for write
out_genotypeTotal=open(out_genotypeTotal_file,"w")
out_genotypeVerySensitive=open(out_genotypeVerySensitive_file,"w")
out_genotypeSensitive=open(out_genotypeSensitive_file,"w")
out_genotypeSpecific=open(out_genotypeSpecific_file,"w")

MolecularN=len(MoleculeDict)
bar = Bar('Genotyping', max=MolecularN)
for m in MoleculeDict:
    # print(m)
    CellBC=m.split('_')[0]
    ## Create predefined arrary to collect single and double stranded genotype, both of which are an array, each coloum is a base, A, C, G, T, N; each row is a position
    SG_Genotypes=np.zeros((max_bp,5))
    DB_Genotypes=np.zeros((max_bp,5))
    Strand_mtx=np.zeros((max_bp,2)) ## 0 will be + or forward,  1 will be - or reverse
    for read_pair in MoleculeDict[m]:
        seq_0=ReadPairDict[read_pair][0].seq
        seq_1=ReadPairDict[read_pair][1].seq
        quality_0=ReadPairDict[read_pair][0].query_qualities
        quality_1=ReadPairDict[read_pair][1].query_qualities
        pos_array_0=np.asarray(ReadPairDict[read_pair][0].get_aligned_pairs(matches_only=True))
        pos_array_1=np.asarray(ReadPairDict[read_pair][1].get_aligned_pairs(matches_only=True))
        pos_array_overlap=np.intersect1d(pos_array_0[:,1],pos_array_1[:,1])
        pos_array_specific_0=pos_array_0[~np.isin(pos_array_0[:,1],pos_array_overlap)]
        pos_array_specific_1=pos_array_1[~np.isin(pos_array_1[:,1],pos_array_overlap)]
        pos_array_overlap_0=pos_array_0[np.isin(pos_array_0[:,1],pos_array_overlap)]
        pos_array_overlap_1=pos_array_1[np.isin(pos_array_1[:,1],pos_array_overlap)]
        ## Collect genotype for the specific_0, or the non-overlapped left part
        if len(pos_array_specific_0)>0:
            for base_0 in pos_array_specific_0:
                if quality_0[base_0[0]]>BaseQ_thld_hi:
                    SG_Genotypes[base_0[1],dna_letters.index(seq_0[base_0[0]])]+=1
                    Strand_mtx[base_0[1],int(ReadPairDict[read_pair][0].is_reverse)]+=1
                else:
                    SG_Genotypes[base_0[1],4]+=1
        ## Collect genotype for the overlap part
        if len(pos_array_overlap)>0:
            for base_0,base_1 in zip(pos_array_overlap_0,pos_array_overlap_1):
                if (seq_0[base_0[0]]==seq_1[base_1[0]]):
                    if(quality_0[base_0[0]]>BaseQ_thld_hi or quality_1[base_1[0]]>BaseQ_thld_hi):
                        DB_Genotypes[base_0[1],dna_letters.index(seq_0[base_0[0]])]+=1
                        Strand_mtx[base_0[1],0]+=1
                        Strand_mtx[base_0[1],1]+=1
                    else:
                        DB_Genotypes[base_0[1],4]+=1
                # else:  # I figured this step is too confusing and not making sense the DB with conflict is worse then even the SG, not woth it to be included into the DB_Genotypes, maybe we can include those into SG_Genotypes, but here wanted to keep simple
                #     ConflictBases=(seq_0[base_0[0]],seq_1[base_1[0]])
                #     ConflictQ=(quality_0[base_0[0]],quality_1[base_1[0]])
                #     if max(ConflictQ)> BaseQ_thld_hi and min(ConflictQ)< BaseQ_thld_mid:
                #         DB_Genotypes[base_0[1],dna_letters.index(ConflictBases[ConflictQ.index(max(ConflictQ))])]+=1
                #         Strand_mtx[base_0[1],0]+=1
                #         Strand_mtx[base_0[1],1]+=1
                else:
                    DB_Genotypes[base_0[1],4]+=1
        ## Collect genotype for the specific_1, or the non-overlapped right part
        if len(pos_array_specific_1)>0:
            for base_1 in pos_array_specific_1:
                if quality_1[base_1[0]]>BaseQ_thld_hi:
                    SG_Genotypes[base_1[1],dna_letters.index(seq_1[base_1[0]])]+=1
                    Strand_mtx[base_1[1],int(ReadPairDict[read_pair][1].is_reverse)]+=1
                else:
                    SG_Genotypes[base_1[1],4]+=1
    for i in np.where(np.sum((SG_Genotypes+DB_Genotypes),axis=1)>0)[0]:
        Cur_Genotype_array=(SG_Genotypes+DB_Genotypes)[i][0:4]
        FamSize=sum(Cur_Genotype_array)
        if FamSize>0:
            CallIndex=Cur_Genotype_array.tolist().index(max(Cur_Genotype_array))
            Call=dna_letters[CallIndex]
            Ref=mito_ref["base"][i].upper()
            Variant=str(i+1)+"_"+Ref+"_"+Call
            GT_Cts=Cur_Genotype_array[CallIndex]
            SG_Cts=SG_Genotypes[i][CallIndex]
            DB_Cts=DB_Genotypes[i][CallIndex]
            CSS=GT_Cts/FamSize
            Strand=((Strand_mtx>0).astype(int))[i]
            TotalMoleculeCtsMatrix[CellBC][i][0]+=1
            OUT=m+"\t"+m.split("_")[0]+"\t"+str(i+1)+"\t"+Variant+"\t"+Call+"\t"+Ref+"\t"+str(FamSize)+"\t"+str(GT_Cts)+"\t"+str(CSS)+"\t"+str(DB_Cts)+"\t"+str(SG_Cts)+"\t"+str(Strand[0])+"\t"+str(Strand[1])+"\n" ## Note those matrix are 0-16568, therefore position ot i need to add 1 to match the 1-based corrdinates
            if not Call==Ref:
                out_genotypeTotal.write(OUT)
            if DB_Cts==0:
                if CSS>0.75 and FamSize>=2:  ##VerySensitive
                    TotalMoleculeCtsMatrix[CellBC][i][1]+=1
                    if not Call==Ref:
                        out_genotypeVerySensitive.write(OUT)
                if CSS>0.75 and FamSize>=3:  ##Sensitive
                    TotalMoleculeCtsMatrix[CellBC][i][2]+=1
                    if not Call==Ref:
                        out_genotypeSensitive.write(OUT)
                if CSS>0.9 and FamSize>=4:   ##Specific
                    TotalMoleculeCtsMatrix[CellBC][i][3]+=1
                    if not Call==Ref:
                        out_genotypeSpecific.write(OUT)
            else:
                if CSS>0.75 and FamSize>=1:  ##VerySensitive
                    TotalMoleculeCtsMatrix[CellBC][i][1]+=1
                    if not Call==Ref:
                        out_genotypeVerySensitive.write(OUT)
                if CSS>0.75 and FamSize>=2:  ##Sensitive
                    TotalMoleculeCtsMatrix[CellBC][i][2]+=1
                    if not Call==Ref:
                        out_genotypeSensitive.write(OUT)
                if CSS>0.9 and FamSize>=3:   ##Specific
                    TotalMoleculeCtsMatrix[CellBC][i][3]+=1
                    if not Call==Ref:
                        out_genotypeSpecific.write(OUT)
    bar.next()

bar.finish()

out_genotypeTotal.close()
out_genotypeVerySensitive.close()
out_genotypeSensitive.close()
out_genotypeSpecific.close()

######### Print out the qualified total counts, aka qualified depth
print("Start to print Qualified Total Molecule Counts")
with open(out_totalCts_file,"w") as out_totalCts:
    for Cell in TotalMoleculeCtsMatrix.keys():
        for pos in range(0,len(TotalMoleculeCtsMatrix[Cell])):
            out_totalCts.write(Cell+"\t"+str(pos+1)+"\t"+str(TotalMoleculeCtsMatrix[Cell][pos][0])+"\t"+str(TotalMoleculeCtsMatrix[Cell][pos][1])+"\t"+str(TotalMoleculeCtsMatrix[Cell][pos][2])+"\t"+str(TotalMoleculeCtsMatrix[Cell][pos][3])+"\n")

########## Below is a session of records that might be helpful in refreshing memory and debugging

#To check Genotype array
#An example how the Genotype array look like,  this is an 16569 row * 5 columns array
#         A,  C,  G,  T,  N
# array([[0., 1., 0., 0., 0.],
#        [1., 0., 0., 0., 0.],
#        [0., 1., 0., 0., 0.],
#        [1., 0., 0., 0., 0.],
#        [0., 1., 0., 0., 0.]])
# SG_Genotypes[SG_Genotypes.sum(axis=1)>0]  ## Show SG_Genotypes non zero chunk
# DB_Genotypes[DB_Genotypes.sum(axis=1)>0]  ## Show DB_Genotypes non zero chunk
# Strand[Strand.sum(axis=1)>0]

## A chunk quickly examine the Start_0 End_0, Start_1, End_1.
# for read_pair in MoleculeDict[m]:
#     seq_0=ReadPairDict[read_pair][0].seq
#     seq_1=ReadPairDict[read_pair][1].seq
#     quality_0=ReadPairDict[read_pair][0].query_qualities
#     quality_1=ReadPairDict[read_pair][1].query_qualities
#     align_qual_0 = ReadPairDict[read_pair][0].mapping_quality
#     align_qual_1 = ReadPairDict[read_pair][1].mapping_quality
#     # Get the overlaped position array for both strand
#     pos_array_0=np.asarray(ReadPairDict[read_pair][0].get_aligned_pairs(matches_only=True, with_seq=True))
#     pos_array_1=np.asarray(ReadPairDict[read_pair][1].get_aligned_pairs(matches_only=True, with_seq=True))
#     pos_array_overlap=np.intersect1d(pos_array_0[:,1],pos_array_1[:,1])
#     pos_array_0[~np.isin(pos_array_0[:,1],pos_array_overlap)]
#     pos_array_overlap_0=pos_array_0[np.isin(pos_array_0[:,1],pos_array_overlap)]
#     pos_array_overlap_1=pos_array_1[np.isin(pos_array_1[:,1],pos_array_overlap)]
#     Start_0=pos_array_0[0][1]
#     End_0=pos_array_0[len(pos_array_0)-1][1]
#     Start_1=pos_array_1[0][1]
#     End_1=pos_array_1[len(pos_array_1)-1][1]
#     CurrentStich=Start_0+'.'+End_0+'.'+Start_1+'.'+End_1+'.'+str(len(pos_array_overlap))
#     print(CurrentStich)
