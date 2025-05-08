import os
import sys
import pysam
import numpy as np
import pandas as pd
import concurrent.futures
import argparse
from collections import defaultdict
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
    with pysam.AlignmentFile(bam_file, "rb") as bam_input:
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

import numpy as np
import logging
from typing import Dict, List, Set, Tuple, Any

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("genotype_matrices.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

def generate_genotype_matrices(
    molecule_dict: Dict, 
    read_pair_dict: Dict, 
    bcs: Set[str], 
    mito_ref: Any, 
    BaseQ_thld_hi: int, 
    out_genotypeTotal_file: str, 
    out_genotypeVerySensitive_file: str,
    out_genotypeSensitive_file: str, 
    out_genotypeSpecific_file: str, 
    out_totalCts_file: str,
    debug: bool = False
) -> Dict:
    """
    Section 3 Main Function: Genotype each molecule from the given data and barcode set.

    Outputs:
    - Four .RawGenotypes outputs:
        1. Total: Without any consensus level filtering.
        2. Very Sensitive (a=2, b=1, c=0.75)
        3. Sensitive (a=3, b=2, c=0.75)
        4. Specific (a=4, b=3, c=0.9)
    - One .QualifiedTotalCts output with 6 columns (Cell, Position, Total, VerySensitive, Sensitive, Specific).

    Parameters:
    -----------
    molecule_dict : Dict
        Dictionary containing molecule information
    read_pair_dict : Dict
        Dictionary containing read pair information
    bcs : Set[str]
        Set of cell barcodes
    mito_ref : Any
        Reference mitogenome
    BaseQ_thld_hi : int
        Base quality threshold (high)
    out_genotypeTotal_file : str
        Output file path for total genotypes
    out_genotypeVerySensitive_file : str
        Output file path for very sensitive genotypes
    out_genotypeSensitive_file : str
        Output file path for sensitive genotypes
    out_genotypeSpecific_file : str
        Output file path for specific genotypes
    out_totalCts_file : str
        Output file path for total counts
    debug : bool, optional
        Enable detailed debug output, by default False

    Returns:
    --------
    Dict
        Dictionary containing genotype matrices for each barcode
    """
    try:
        dna_letters = ['A', 'C', 'G', 'T', 'N']
        max_bp = mito_ref.shape[0]
        logger.info(f"Processing genotype matrices for {len(bcs)} barcodes with max_bp={max_bp}")
        
        # Matrix to collect Total, Very_sensitive, Sensitive, Specific counts
        TotalMoleculeCtsMatrix = {key: np.zeros((max_bp, 4)) for key in bcs}
        
        # Track statistics for reporting
        stats = {
            "total_molecules": 0,
            "processed_molecules": 0,
            "total_variants": 0,
            "very_sensitive_variants": 0,
            "sensitive_variants": 0,
            "specific_variants": 0
        }

        ## Open the 4 files for write
        with open(out_genotypeTotal_file, "w") as out_genotypeTotal, \
             open(out_genotypeVerySensitive_file, "w") as out_genotypeVerySensitive, \
             open(out_genotypeSensitive_file, "w") as out_genotypeSensitive, \
             open(out_genotypeSpecific_file, "w") as out_genotypeSpecific:

            # Write headers to all files
            header = "MoleculeID\tCellBC\tPosition\tVariant\tCall\tRef\tFamSize\tGT_Cts\tCSS\tDB_Cts\tSG_Cts\tForwardStrand\tReverseStrand\n"
            out_genotypeTotal.write(header)
            out_genotypeVerySensitive.write(header)
            out_genotypeSensitive.write(header)
            out_genotypeSpecific.write(header)
            
            # Track progress
            total_molecules = len(molecule_dict)
            stats["total_molecules"] = total_molecules
            
            for idx, m in enumerate(molecule_dict):
                if idx % 100 == 0:
                    logger.info(f"Processing molecule {idx+1}/{total_molecules} ({(idx+1)/total_molecules*100:.1f}%)")
                
                if debug:
                    logger.debug(f"Processing molecule: {m}")
                
                try:
                    CellBC = m.split('_')[0]
                    
                    # Create predefined arrays to collect single and double stranded genotype
                    # Each column is a base (A, C, G, T, N); each row is a position
                    SG_Genotypes = np.zeros((max_bp, 5))
                    DB_Genotypes = np.zeros((max_bp, 5))
                    Strand_mtx = np.zeros((max_bp, 2))  # 0=forward, 1=reverse
                    
                    for read_pair in molecule_dict[m]:
                        try:
                            # Extract sequence data
                            seq_0 = read_pair_dict[read_pair][0].seq
                            seq_1 = read_pair_dict[read_pair][1].seq
                            quality_0 = read_pair_dict[read_pair][0].query_qualities
                            quality_1 = read_pair_dict[read_pair][1].query_qualities
                            
                            # Get aligned positions
                            pos_array_0 = np.asarray(read_pair_dict[read_pair][0].get_aligned_pairs(matches_only=True))
                            pos_array_1 = np.asarray(read_pair_dict[read_pair][1].get_aligned_pairs(matches_only=True))
                            
                            # Handle empty arrays
                            if len(pos_array_0) == 0 or len(pos_array_1) == 0:
                                logger.warning(f"Empty position array for read pair {read_pair} in molecule {m}")
                                continue
                                
                            pos_array_overlap = np.intersect1d(pos_array_0[:, 1], pos_array_1[:, 1])
                            pos_array_specific_0 = pos_array_0[~np.isin(pos_array_0[:, 1], pos_array_overlap)]
                            pos_array_specific_1 = pos_array_1[~np.isin(pos_array_1[:, 1], pos_array_overlap)]
                            pos_array_overlap_0 = pos_array_0[np.isin(pos_array_0[:, 1], pos_array_overlap)]
                            pos_array_overlap_1 = pos_array_1[np.isin(pos_array_1[:, 1], pos_array_overlap)]
                            
                            # Collect genotype for non-overlapped left part
                            if len(pos_array_specific_0) > 0:
                                for base_0 in pos_array_specific_0:
                                    if quality_0[base_0[0]] > BaseQ_thld_hi:
                                        SG_Genotypes[base_0[1], dna_letters.index(seq_0[base_0[0]])] += 1
                                        Strand_mtx[base_0[1], int(read_pair_dict[read_pair][0].is_reverse)] += 1
                                    else:
                                        SG_Genotypes[base_0[1], 4] += 1
                                        
                            # Collect genotype for the overlap part
                            if len(pos_array_overlap) > 0:
                                for base_0, base_1 in zip(pos_array_overlap_0, pos_array_overlap_1):
                                    if seq_0[base_0[0]] == seq_1[base_1[0]]:
                                        if quality_0[base_0[0]] > BaseQ_thld_hi or quality_1[base_1[0]] > BaseQ_thld_hi:
                                            DB_Genotypes[base_0[1], dna_letters.index(seq_0[base_0[0]])] += 1
                                            Strand_mtx[base_0[1], 0] += 1
                                            Strand_mtx[base_0[1], 1] += 1
                                        else:
                                            DB_Genotypes[base_0[1], 4] += 1
                                    else:
                                        DB_Genotypes[base_0[1], 4] += 1
                                        
                            # Collect genotype for non-overlapped right part
                            if len(pos_array_specific_1) > 0:
                                for base_1 in pos_array_specific_1:
                                    if quality_1[base_1[0]] > BaseQ_thld_hi:
                                        SG_Genotypes[base_1[1], dna_letters.index(seq_1[base_1[0]])] += 1
                                        Strand_mtx[base_1[1], int(read_pair_dict[read_pair][1].is_reverse)] += 1
                                    else:
                                        SG_Genotypes[base_1[1], 4] += 1
                                        
                        except Exception as e:
                            logger.error(f"Error processing read pair {read_pair} in molecule {m}: {str(e)}")
                            continue
                    
                    # Buffer output data for batch writing
                    buffer_total = []
                    buffer_very_sensitive = []
                    buffer_sensitive = []
                    buffer_specific = []
                    
                    # Find positions with data
                    positions_with_data = np.where(np.sum((SG_Genotypes + DB_Genotypes), axis=1) > 0)[0]
                    
                    for i in positions_with_data:
                        Cur_Genotype_array = (SG_Genotypes + DB_Genotypes)[i][0:4]
                        FamSize = sum(Cur_Genotype_array)
                        
                        if FamSize > 0:
                            CallIndex = Cur_Genotype_array.tolist().index(max(Cur_Genotype_array))
                            Call = dna_letters[CallIndex]
                            
                            if debug:
                                logger.debug(f"Position {i+1}, Call: {Call}, FamSize: {FamSize}")
                                
                            Ref = mito_ref["base"][i].upper()
                            Variant = str(i + 1) + "_" + Ref + "_" + Call
                            GT_Cts = Cur_Genotype_array[CallIndex]
                            SG_Cts = SG_Genotypes[i][CallIndex]
                            DB_Cts = DB_Genotypes[i][CallIndex]
                            CSS = GT_Cts / FamSize
                            Strand = ((Strand_mtx > 0).astype(int))[i]
                            
                            TotalMoleculeCtsMatrix[CellBC][i][0] += 1
                            
                            OUT = (
                                f"{m}\t{CellBC}\t{i+1}\t{Variant}\t{Call}\t{Ref}\t{FamSize}\t"
                                f"{GT_Cts}\t{CSS:.4f}\t{DB_Cts}\t{SG_Cts}\t{Strand[0]}\t{Strand[1]}\n"
                            )

                            if Call != Ref:
                                buffer_total.append(OUT)
                                stats["total_variants"] += 1

                            # Apply filtering criteria
                            if DB_Cts == 0:  # Single strand
                                if CSS > 0.75 and FamSize >= 2:
                                    TotalMoleculeCtsMatrix[CellBC][i][1] += 1
                                    if Call != Ref:
                                        buffer_very_sensitive.append(OUT)
                                        stats["very_sensitive_variants"] += 1
                                if CSS > 0.75 and FamSize >= 3:
                                    TotalMoleculeCtsMatrix[CellBC][i][2] += 1
                                    if Call != Ref:
                                        buffer_sensitive.append(OUT)
                                        stats["sensitive_variants"] += 1
                                if CSS > 0.9 and FamSize >= 4:
                                    TotalMoleculeCtsMatrix[CellBC][i][3] += 1
                                    if Call != Ref:
                                        buffer_specific.append(OUT)
                                        stats["specific_variants"] += 1
                            else:  # Double strand
                                if CSS > 0.75 and FamSize >= 1:
                                    TotalMoleculeCtsMatrix[CellBC][i][1] += 1
                                    if Call != Ref:
                                        buffer_very_sensitive.append(OUT)
                                        stats["very_sensitive_variants"] += 1
                                if CSS > 0.75 and FamSize >= 2:
                                    TotalMoleculeCtsMatrix[CellBC][i][2] += 1
                                    if Call != Ref:
                                        buffer_sensitive.append(OUT)
                                        stats["sensitive_variants"] += 1
                                if CSS > 0.9 and FamSize >= 3:
                                    TotalMoleculeCtsMatrix[CellBC][i][3] += 1
                                    if Call != Ref:
                                        buffer_specific.append(OUT)
                                        stats["specific_variants"] += 1
                    
                    # Write buffered data for this molecule
                    if buffer_total:
                        out_genotypeTotal.writelines(buffer_total)
                    if buffer_very_sensitive:
                        out_genotypeVerySensitive.writelines(buffer_very_sensitive)
                    if buffer_sensitive:
                        out_genotypeSensitive.writelines(buffer_sensitive)
                    if buffer_specific:
                        out_genotypeSpecific.writelines(buffer_specific)
                    
                    stats["processed_molecules"] += 1
                    
                except Exception as e:
                    logger.error(f"Error processing molecule {m}: {str(e)}")
            
            logger.info(f"Molecule processing complete. Processed {stats['processed_molecules']} of {stats['total_molecules']} molecules.")
        
        # Print out the qualified total counts (qualified depth)
        with open(out_totalCts_file, "w") as out_totalCts:
            # Write header
            out_totalCts.write("CellBC\tPosition\tTotal\tVerySensitive\tSensitive\tSpecific\n")
            
            for Cell in TotalMoleculeCtsMatrix.keys():
                for pos in range(len(TotalMoleculeCtsMatrix[Cell])):
                    out_totalCts.write(
                        f"{Cell}\t{pos+1}\t"
                        f"{TotalMoleculeCtsMatrix[Cell][pos][0]}\t"
                        f"{TotalMoleculeCtsMatrix[Cell][pos][1]}\t"
                        f"{TotalMoleculeCtsMatrix[Cell][pos][2]}\t"
                        f"{TotalMoleculeCtsMatrix[Cell][pos][3]}\n"
                    )
                    
        logger.info(f"Results summary:")
        logger.info(f"  - Total molecules: {stats['total_molecules']}")
        logger.info(f"  - Processed molecules: {stats['processed_molecules']}")
        logger.info(f"  - Total variants: {stats['total_variants']}")
        logger.info(f"  - Very sensitive variants: {stats['very_sensitive_variants']}")
        logger.info(f"  - Sensitive variants: {stats['sensitive_variants']}")
        logger.info(f"  - Specific variants: {stats['specific_variants']}")
        
        return TotalMoleculeCtsMatrix
    
    except Exception as e:
        logger.error(f"Fatal error in generate_genotype_matrices: {str(e)}")
        raise
# def generate_genotype_matrices(molecule_dict, read_pair_dict, bcs, mito_ref, BaseQ_thld_hi, out_genotypeTotal_file, out_genotypeVerySensitive_file,out_genotypeSensitive_file, out_genotypeSpecific_file, out_totalCts_file):
#     """
#     Section 3 Main Function: Genotype each molecule from the given BAM file and barcode set.

#     Outputs:
#     - Four .RawGenotypes outputs:
#         1. Total: Without any consensus level filtering.
#         2. Very Sensitive (a=2, b=1, c=0.75)
#         3. Sensitive (a=3, b=2, c=0.75)
#         4. Specific (a=4, b=3, c=0.9)
#     - One .QualifiedTotalCts output with 4 columns.

#     :param bam_file: Path to the BAM file.
#     :param bcs: Set of barcodes.
#     :param mito_ref: Reference mitogenome.
#     :param dna_letters: DNA letters for indexing (e.g., ['A', 'C', 'G', 'T']).
#     :param BaseQ_thld_hi: Quality threshold.
#     :param max_bp: Maximum base pairs.
#     :param out_files: List of output file paths (Total, VerySensitive, Sensitive, Specific).
#     :returen: A dictionary containing genotype matrices.
#     """
#     dna_letters = ['A','C','G','T','N']
#     max_bp=mito_ref.shape[0]
#     TotalMoleculeCtsMatrix={key:np.zeros((max_bp,4)) for key in bcs} # Matrix to collect 1,Total;Very_sensitive;Sensitive;Specific

#     ## Open the 4 files for write
#     with open(out_genotypeTotal_file, "w") as out_genotypeTotal, \
#          open(out_genotypeVerySensitive_file, "w") as out_genotypeVerySensitive, \
#          open(out_genotypeSensitive_file, "w") as out_genotypeSensitive, \
#          open(out_genotypeSpecific_file, "w") as out_genotypeSpecific:

#         for m in molecule_dict:
#             print("molecule: " + m)
#             CellBC=m.split('_')[0]
#             ## Create predefined arrary to collect single and double stranded genotype, both of which are an array, each coloum is a base, A, C, G, T, N; each row is a position
#             SG_Genotypes=np.zeros((max_bp,5))
#             DB_Genotypes=np.zeros((max_bp,5))
#             Strand_mtx=np.zeros((max_bp,2)) ## 0 will be + or forward,  1 will be - or reverse
#             for read_pair in molecule_dict[m]:
#                 seq_0=read_pair_dict[read_pair][0].seq
#                 seq_1=read_pair_dict[read_pair][1].seq
#                 quality_0=read_pair_dict[read_pair][0].query_qualities
#                 quality_1=read_pair_dict[read_pair][1].query_qualities
#                 pos_array_0=np.asarray(read_pair_dict[read_pair][0].get_aligned_pairs(matches_only=True))
#                 pos_array_1=np.asarray(read_pair_dict[read_pair][1].get_aligned_pairs(matches_only=True))
#                 pos_array_overlap=np.intersect1d(pos_array_0[:,1],pos_array_1[:,1])
#                 pos_array_specific_0=pos_array_0[~np.isin(pos_array_0[:,1],pos_array_overlap)]
#                 pos_array_specific_1=pos_array_1[~np.isin(pos_array_1[:,1],pos_array_overlap)]
#                 pos_array_overlap_0=pos_array_0[np.isin(pos_array_0[:,1],pos_array_overlap)]
#                 pos_array_overlap_1=pos_array_1[np.isin(pos_array_1[:,1],pos_array_overlap)]
#                 ## Collect genotype for the specific_0, or the non-overlapped left part
#                 if len(pos_array_specific_0)>0:
#                     for base_0 in pos_array_specific_0:
#                         if quality_0[base_0[0]]>BaseQ_thld_hi:
#                             SG_Genotypes[base_0[1],dna_letters.index(seq_0[base_0[0]])]+=1
#                             Strand_mtx[base_0[1],int(read_pair_dict[read_pair][0].is_reverse)]+=1
#                         else:
#                             SG_Genotypes[base_0[1],4]+=1
#                 ## Collect genotype for the overlap part
#                 if len(pos_array_overlap)>0:
#                     for base_0,base_1 in zip(pos_array_overlap_0,pos_array_overlap_1):
#                         if (seq_0[base_0[0]]==seq_1[base_1[0]]):
#                             if(quality_0[base_0[0]]>BaseQ_thld_hi or quality_1[base_1[0]]>BaseQ_thld_hi):
#                                 DB_Genotypes[base_0[1],dna_letters.index(seq_0[base_0[0]])]+=1
#                                 Strand_mtx[base_0[1],0]+=1
#                                 Strand_mtx[base_0[1],1]+=1
#                             else:
#                                 DB_Genotypes[base_0[1],4]+=1
#                         else:
#                             DB_Genotypes[base_0[1],4]+=1
#                 ## Collect genotype for the specific_1, or the non-overlapped right part
#                 if len(pos_array_specific_1)>0:
#                     for base_1 in pos_array_specific_1:
#                         if quality_1[base_1[0]]>BaseQ_thld_hi:
#                             SG_Genotypes[base_1[1],dna_letters.index(seq_1[base_1[0]])]+=1
#                             Strand_mtx[base_1[1],int(read_pair_dict[read_pair][1].is_reverse)]+=1
#                         else:
#                             SG_Genotypes[base_1[1],4]+=1
#             buffer_total = []
#             buffer_very_sensitive = []
#             buffer_sensitive = []
#             buffer_specific = []
#             for i in np.where(np.sum((SG_Genotypes + DB_Genotypes), axis=1) > 0)[0]:
#                 Cur_Genotype_array = (SG_Genotypes + DB_Genotypes)[i][0:4]
#                 FamSize = sum(Cur_Genotype_array)
#                 if FamSize > 0:
#                     CallIndex = Cur_Genotype_array.tolist().index(max(Cur_Genotype_array))
#                     Call = dna_letters[CallIndex]
#                     print(f"call: {Call}")
#                     Ref = mito_ref["base"][i].upper()
#                     Variant = str(i + 1) + "_" + Ref + "_" + Call
#                     GT_Cts = Cur_Genotype_array[CallIndex]
#                     SG_Cts = SG_Genotypes[i][CallIndex]
#                     DB_Cts = DB_Genotypes[i][CallIndex]
#                     CSS = GT_Cts / FamSize
#                     Strand = ((Strand_mtx > 0).astype(int))[i]
#                     TotalMoleculeCtsMatrix[CellBC][i][0] += 1
#                     OUT = (
#                         f"{m}\t{m.split('_')[0]}\t{i+1}\t{Variant}\t{Call}\t{Ref}\t{FamSize}\t"
#                         f"{GT_Cts}\t{CSS}\t{DB_Cts}\t{SG_Cts}\t{Strand[0]}\t{Strand[1]}\n"
#                     )

#                     if not Call == Ref:
#                         buffer_total.append(OUT)

#                     if DB_Cts == 0:
#                         if CSS > 0.75 and FamSize >= 2:
#                             TotalMoleculeCtsMatrix[CellBC][i][1] += 1
#                             if not Call == Ref:
#                                 buffer_very_sensitive.append(OUT)
#                         if CSS > 0.75 and FamSize >= 3:
#                             TotalMoleculeCtsMatrix[CellBC][i][2] += 1
#                             if not Call == Ref:
#                                 buffer_sensitive.append(OUT)
#                         if CSS > 0.9 and FamSize >= 4:
#                             TotalMoleculeCtsMatrix[CellBC][i][3] += 1
#                             if not Call == Ref:
#                                 buffer_specific.append(OUT)
#                     else:
#                         if CSS > 0.75 and FamSize >= 1:
#                             TotalMoleculeCtsMatrix[CellBC][i][1] += 1
#                             if not Call == Ref:
#                                 buffer_very_sensitive.append(OUT)
#                         if CSS > 0.75 and FamSize >= 2:
#                             TotalMoleculeCtsMatrix[CellBC][i][2] += 1
#                             if not Call == Ref:
#                                 buffer_sensitive.append(OUT)
#                         if CSS > 0.9 and FamSize >= 3:
#                             TotalMoleculeCtsMatrix[CellBC][i][3] += 1
#                             if not Call == Ref:
#                                 buffer_specific.append(OUT)

#         # Write all buffered data at once
#         out_genotypeTotal.writelines(buffer_total)
#         out_genotypeVerySensitive.writelines(buffer_very_sensitive)
#         out_genotypeSensitive.writelines(buffer_sensitive)
#         out_genotypeSpecific.writelines(buffer_specific)

       
    
#     ######### Print out the qualified total counts, aka qualified depth
#     with open(out_totalCts_file,"w") as out_totalCts:
#         for Cell in TotalMoleculeCtsMatrix.keys():
#             for pos in range(0,len(TotalMoleculeCtsMatrix[Cell])):
#                 out_totalCts.write(Cell+"\t"+str(pos+1)+"\t"+str(TotalMoleculeCtsMatrix[Cell][pos][0])+"\t"+str(TotalMoleculeCtsMatrix[Cell][pos][1])+"\t"+str(TotalMoleculeCtsMatrix[Cell][pos][2])+"\t"+str(TotalMoleculeCtsMatrix[Cell][pos][3])+"\n")

def generate_filenames(sample, bam_path=None, barcode_path=None):
    """
    Generate filenames based on the given parameters.

    Parameters:
    - sample: The prefix. e.g. barcodes. Should check the prefix of /temp/barcoded_bams/.
    - out_dir: The root directory for the analysis. e.g. Out_mitoConsensus (check --output in Preprocess.py)
    - bam_path: Optional. Direct path to the BAM file. If provided, this will be used instead of constructing the path.
    - barcode_path: Optional. Direct path to the barcode file. If provided, this will be used instead of constructing the path.

    Returns:
    A list containing paths to various files based on the parameters.
    """

    if not bam_path:
        sys.exit("Error: bam_path not provided.")

    if not barcode_path: 
        sys.exit("Error: barcode_path not provided.")
    
    

    
    # Determine the directory where temp files are located, use bam_path directory for outputs

    
    
   # Print the current working directory for debugging
    current_working_dir = os.getcwd()
    print(f"Current working directory: {current_working_dir}")

    base_dir = os.path.dirname(os.path.dirname(bam_path))  # Go up one level from barcoded_bams
    sparse_dir = os.path.join(current_working_dir, "sparse_matrices2.0")  # Use absolute path to current working directory

    # These paths don't typically change, so we keep the same construction
    mito_ref_file = f"{base_dir}/final/chrM_refAllele.txt"
    
    # Ensure the sparse matrices directory exists
    os.makedirs(sparse_dir, exist_ok=True)
    
    # Output files
    out_genotypeTotal_file = f"{sparse_dir}/{sample}.RawGenotypes.Total"
    out_genotypeVerySensitive_file = f"{sparse_dir}/{sample}.RawGenotypes.VerySensitive"
    out_genotypeSensitive_file = f"{sparse_dir}/{sample}.RawGenotypes.Sensitive"
    out_genotypeSpecific_file = f"{sparse_dir}/{sample}.RawGenotypes.Specific"
    out_totalCts_file = f"{sparse_dir}/{sample}.QualifiedTotalCts"

    return [
        barcode_path,
        bam_path,
        mito_ref_file,
        out_genotypeTotal_file,
        out_genotypeVerySensitive_file,
        out_genotypeSensitive_file,
        out_genotypeSpecific_file,
        out_totalCts_file
    ]








def run_mito_consensus(bam_file, barcode_file, chrM_ref, BaseQ_thld_hi=30):
    """
    Process the data based on the given BAM file, barcode file, and output directory.

    Parameters:
    - bam_file: The BAM file.
    - barcode_file: The barcode file.
    - output_dir: The directory where the output files will be saved.
    - BaseQ_thld_hi: Base Quality threshold, default is 30.

    Returns:
    None. Operations are applied in-place and output is saved to the respective files.
    """

    # Generate filenames based on the BAM file name (extracted from bam_file)
    prefix = os.path.splitext(os.path.basename(bam_file))[0]
    file_list = generate_filenames(sample=prefix, bam_path = bam_file, barcode_path =barcode_file)
    print(file_list)
    # Import a subset of cell barcodes
    with open(barcode_file, 'r') as barcode_file_handle:
        cell_barcodes = [line.strip() for line in barcode_file_handle]

    # Build molecule dictionary
    subset_read_pair_dict, subset_molecule_dict = build_molecule_dict(file_list[1], "BC")

    # Read mito reference
    mito_ref_pd = pd.read_table(chrM_ref, names=["pos", "base"])


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







if __name__ == "__main__":
    # Argument parser setup for Nextflow parameters
    parser = argparse.ArgumentParser(description="Process barcodes and bam files for mito consensus.")
    parser.add_argument('--bam', type=str, required=True, help='Path to the BAM file.')
    parser.add_argument('--barcode', type=str, required=True, help='Path to the barcode file.')
    parser.add_argument('--chrM_ref', type=str, required=True, help='Path to the chrm_ref.')
    parser.add_argument('--BaseQ_thld_hi', type=int, default=30, help='Base Quality threshold (default is 30).')
    args = parser.parse_args()

    # Run mito consensus
    run_mito_consensus(args.bam, args.barcode, args.chrM_ref, args.BaseQ_thld_hi)
