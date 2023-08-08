import os
import pysam
import numpy as np
import pandas as pd
import concurrent.futures
import argparse
from collections import defaultdict
from tqdm import tqdm
import warnings

def build_molecule_dict(bam_file, barcode_tag):
    """
    Build a dictionary from Rname->read1,read2 and Molecule(Cell_Start_End)-->Rname.

    Parameters:
    - bam_file: path to the BAM file.
    - barcode_tag: the tag used to retrieve barcodes.

    Returns:
    - local_read_pair_dict: dictionary containing read pairs.
    - local_molecule_dict: dictionary containing molecular information.
    """
    
    # Open the BAM file for reading
    bam_input = pysam.AlignmentFile(bam_file, "rb")

    # Create a dictionary to hold the read pairs
    local_read_pair_dict = defaultdict(list)
    # Create a dictionary to hold the molecule data
    local_molecule_dict = defaultdict(list)
    
    # Populate local_read_pair_dict with read query names and their corresponding reads
    for read in bam_input:
        local_read_pair_dict[read.query_name].append(read)

    # Populate local_molecule_dict based on the molecule information
    # (cell barcode and position) and the read query names
    for read_name, reads in local_read_pair_dict.items():
        # disregard singlets and multiplets
        if len(reads) != 2:
            continue
        read0, read1 = reads
        # identify fwd and rev in a pair
        if read0.is_reverse and not read1.is_reverse:
            fwd_read, rev_read = read1, read0
        elif not read0.is_reverse and read1.is_reverse:
            fwd_read, rev_read = read0, read1
        else:
            # disregard a pair if both are the same strand
            continue
        if fwd_read.has_tag(barcode_tag):
            cell_bc = fwd_read.get_tag(barcode_tag)
            molecule = f"{cell_bc}_{fwd_read.pos}_{fwd_read.pos + abs(fwd_read.tlen)}"
            local_molecule_dict[molecule].append(read_name)

    return local_read_pair_dict, local_molecule_dict


def generate_genotype_matrices(molecule_dict, read_pair_dict, bcs, mito_ref, BaseQ_thld_hi, out_genotypeTotal_file, out_genotypeVerySensitive_file,out_genotypeSensitive_file, out_genotypeSpecific_file, out_totalCts_file):
    """
    Section 3 Main Function: Genotype each molecule from the given BAM file and barcode set.

    Outputs:
    - Four .RawGenotypes outputs:
        1. Total: Without any consensus level filtering.
        2. Very Sensitive (a=2, b=1, c=0.75)
        3. Sensitive (a=3, b=2, c=0.75)
        4. Specific (a=4, b=3, c=0.9)
    - One .QualifiedTotalCts output with 4 columns.

    :param bam_file: Path to the BAM file.
    :param bcs: Set of barcodes.
    :param mito_ref: Reference mitogenome.
    :param dna_letters: DNA letters for indexing (e.g., ['A', 'C', 'G', 'T']).
    :param BaseQ_thld_hi: Quality threshold.
    :param max_bp: Maximum base pairs.
    :param out_files: List of output file paths (Total, VerySensitive, Sensitive, Specific).
    :return: A dictionary containing genotype matrices.
    """
    dna_letters = ['A','C','G','T','N']
    max_bp=mito_ref.shape[0]
    TotalMoleculeCtsMatrix={key:np.zeros((max_bp,4)) for key in bcs} # Matrix to collect 1,Total;Very_sensitive;Sensitive;Specific

    ## Open the 4 files for write
    out_genotypeTotal=open(out_genotypeTotal_file,"w")
    out_genotypeVerySensitive=open(out_genotypeVerySensitive_file,"w")
    out_genotypeSensitive=open(out_genotypeSensitive_file,"w")
    out_genotypeSpecific=open(out_genotypeSpecific_file,"w")

    for m in molecule_dict:
        # print(m)
        CellBC=m.split('_')[0]
        ## Create predefined arrary to collect single and double stranded genotype, both of which are an array, each coloum is a base, A, C, G, T, N; each row is a position
        SG_Genotypes=np.zeros((max_bp,5))
        DB_Genotypes=np.zeros((max_bp,5))
        Strand_mtx=np.zeros((max_bp,2)) ## 0 will be + or forward,  1 will be - or reverse
        for read_pair in molecule_dict[m]:
            seq_0=read_pair_dict[read_pair][0].seq
            seq_1=read_pair_dict[read_pair][1].seq
            quality_0=read_pair_dict[read_pair][0].query_qualities
            quality_1=read_pair_dict[read_pair][1].query_qualities
            pos_array_0=np.asarray(read_pair_dict[read_pair][0].get_aligned_pairs(matches_only=True))
            pos_array_1=np.asarray(read_pair_dict[read_pair][1].get_aligned_pairs(matches_only=True))
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
                        Strand_mtx[base_0[1],int(read_pair_dict[read_pair][0].is_reverse)]+=1
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
                    else:
                        DB_Genotypes[base_0[1],4]+=1
            ## Collect genotype for the specific_1, or the non-overlapped right part
            if len(pos_array_specific_1)>0:
                for base_1 in pos_array_specific_1:
                    if quality_1[base_1[0]]>BaseQ_thld_hi:
                        SG_Genotypes[base_1[1],dna_letters.index(seq_1[base_1[0]])]+=1
                        Strand_mtx[base_1[1],int(read_pair_dict[read_pair][1].is_reverse)]+=1
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

    out_genotypeTotal.close()
    out_genotypeVerySensitive.close()
    out_genotypeSensitive.close()
    out_genotypeSpecific.close()
    ######### Print out the qualified total counts, aka qualified depth
    with open(out_totalCts_file,"w") as out_totalCts:
        for Cell in TotalMoleculeCtsMatrix.keys():
            for pos in range(0,len(TotalMoleculeCtsMatrix[Cell])):
                out_totalCts.write(Cell+"\t"+str(pos+1)+"\t"+str(TotalMoleculeCtsMatrix[Cell][pos][0])+"\t"+str(TotalMoleculeCtsMatrix[Cell][pos][1])+"\t"+str(TotalMoleculeCtsMatrix[Cell][pos][2])+"\t"+str(TotalMoleculeCtsMatrix[Cell][pos][3])+"\n")



def generate_filenames(sample, out_dir):
    """
    Generate filenames based on the given parameters.

    Parameters:
    - sample: The prefix. e.g. barcodes. Should check the prefix of /temp/barcoded_bams/.
    - out_dir: The root directory for the analysis. e.g. Out_mitoConsensus (check --output in Preprocess.py)

    Returns:
    A list containing paths to various files based on the parameters.
    """
    barcodes_file = f"{out_dir}/temp/barcode_files/{sample}.txt"
    bam_file = f"{out_dir}/temp/barcoded_bams/{sample}.bam"
    mito_ref_file = f"{out_dir}/final/chrM_refAllele.txt"
    out_genotypeTotal_file = f"{out_dir}/temp/sparse_matrices2.0/{sample}.RawGenotypes.Total"
    out_genotypeVerySensitive_file = f"{out_dir}/temp/sparse_matrices2.0/{sample}.RawGenotypes.VerySensitive"
    out_genotypeSensitive_file = f"{out_dir}/temp/sparse_matrices2.0/{sample}.RawGenotypes.Sensitive"
    out_genotypeSpecific_file = f"{out_dir}/temp/sparse_matrices2.0/{sample}.RawGenotypes.Specific"
    out_totalCts_file = f"{out_dir}/temp/sparse_matrices2.0/{sample}.QualifiedTotalCts"

    return [
        barcodes_file,
        bam_file,
        mito_ref_file,
        out_genotypeTotal_file,
        out_genotypeVerySensitive_file,
        out_genotypeSensitive_file,
        out_genotypeSpecific_file,
        out_totalCts_file
    ]









def run_mito_consensus(prefix, work_dir, BaseQ_thld_hi=30):
    """
    Process the data based on the given prefix and working directory.

    Parameters:
    - prefix: The prefix for the filenames, e.g., "barcodes.1".
    - work_dir: The root directory where the files are located.
    - BaseQ_thld_hi: Base Quality threshold, default is 30.

    Returns:
    None. Operations are applied in-place and output is saved to the respective files.
    """

    # Generate filenames
    file_list = generate_filenames(sample=prefix, out_dir=work_dir)

    # Import a subset of cell barcodes
    with open(file_list[0], 'r') as barcode_file_handle:
        cell_barcodes = [line.strip() for line in barcode_file_handle]

    # Build molecule dictionary
    subset_read_pair_dict, subset_molecule_dict = build_molecule_dict(file_list[1], "BC")

    # Read mito reference
    mito_ref_pd = pd.read_table(file_list[2], names=["pos", "base"])

    # Generate genotype matrices
    generate_genotype_matrices(
        molecule_dict=subset_molecule_dict, 
        read_pair_dict=subset_read_pair_dict, 
        bcs=cell_barcodes, 
        mito_ref=mito_ref_pd, 
        BaseQ_thld_hi=BaseQ_thld_hi, 
        out_genotypeTotal_file=file_list[3], 
        out_genotypeVerySensitive_file=file_list[4],
        out_genotypeSensitive_file=file_list[5], 
        out_genotypeSpecific_file=file_list[6], 
        out_totalCts_file=file_list[7]
    )









# work_dir="/lab/solexa_weissman/cweng/Projects/Collaborator/Roman_NK/Data_220907_Pilot3/REDEEM-V/source/Test/Out_mitoConsensus/"
# num_cores=12

if __name__=="__main__":
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Process barcodes in the given directory.")
    parser.add_argument('work_dir', type=str, help='Path to the working directory.')
    parser.add_argument('--num_cores', type=int, default=12, help='Number of cores to use. Default is 12')
    args = parser.parse_args()
    work_dir=args.work_dir
    num_cores=args.num_cores

    # Check if the directory exists and create if not
    if 'sparse_matrices2.0' not in os.listdir(f"{work_dir}/temp"):
        os.mkdir(f"{work_dir}/temp/sparse_matrices2.0")
    
    
    # Define the directory where the barcode files are located
    barcode_dir = os.path.join(work_dir, "temp/barcode_files")
    prefix_list = [filename[:-4] for filename in os.listdir(barcode_dir)]

    # This will run the function for each prefix in parallel
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = {executor.submit(run_mito_consensus, prefix, work_dir): prefix for prefix in prefix_list}
        
        # Using tqdm with as_completed to display progress
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(prefix_list), desc="Processing", unit="prefix"):
            result = future.result()  # if your function returns something, capture it here
        
        # Final check the completeness of the processing    
        if len(os.listdir(os.path.join(work_dir, "temp/sparse_matrices2.0/")))==5*len(os.listdir(barcode_dir)):
                print("ReDeeM Consensus Variant Calling Completed")
        else:
                warnings.warn("Warnning: There are some processing incomplete!")

