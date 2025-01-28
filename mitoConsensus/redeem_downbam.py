import argparse
import pysam
import random
import os
from collections import defaultdict

os.getcwd()
bam_file = "/lab/solexa_weissman/cweng/Packages/REDEEM-V/example_data/barcoded_bams/barcodes.10.bam"
barcode_tag='BC'
fraction=0.5
def build_dictionaries(
    bam_file, 
    barcode_tag="BC"
):
    """
    Reads a BAM file and builds:
      1) ReadPairDict: read_name -> list of pysam.AlignedSegment
      2) MoleculeDict: molecule_id -> list of read_names
      3) CellDict: cell_barcode -> list of molecule_ids

    :param bam_file: Path to the input BAM file
    :param barcode_tag: The BAM tag where the cell barcode is stored
    :return: (ReadPairDict, MoleculeDict, CellDict, bam_header)
    """
    
    # Open the BAM
    bam_input = pysam.AlignmentFile(bam_file, "rb")
    bam_header = bam_input.header
    
    # Prepare dictionaries
    ReadPairDict = defaultdict(list)
    MoleculeDict = defaultdict(list)
    CellDict = defaultdict(list)
    
    # 1) Read all alignments and group by read_name
    for read in bam_input:
        ReadPairDict[read.query_name].append(read)
    bam_input.close()  # Done reading
    
    # 2) Build MoleculeDict and CellDict
    for read_name, reads in ReadPairDict.items():
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
        
        # extract cell barcode, build molecule ID
        if fwd_read.has_tag(barcode_tag):
            CellBC = fwd_read.get_tag(barcode_tag)
            if fwd_read.tlen > 0:
                start = fwd_read.pos
                end   = fwd_read.pos + fwd_read.tlen
            else:
                start = fwd_read.pos - fwd_read.tlen
                end   = fwd_read.pos

            Molecule = f"{CellBC}_{start}_{end}"
            
            # add to dictionaries
            MoleculeDict[Molecule].append(read_name)
            CellDict[CellBC].append(Molecule)
    
    return ReadPairDict, MoleculeDict, CellDict, bam_header



def downsample_by_molecule_in_memory(
    ReadPairDict,
    MoleculeDict,
    CellDict,
    fraction,
    out_bam,
    bam_header
):
    """
    Downsample reads by randomly keeping a fraction of molecules for each cell
    *in memory*, then write out only the chosen reads to `out_bam`.

    :param ReadPairDict: dict -> read_name -> list of pysam.AlignedSegment
    :param MoleculeDict: dict -> molecule_id -> list of read_names
    :param CellDict: dict -> cell_barcode -> list of molecule_ids
    :param fraction: float (e.g. 0.5 for 50%)
    :param out_bam: str, path to output BAM
    :param bam_header: dict, header object (e.g. from pysam.AlignmentFile(...).header)
    """

    # 1) Decide which molecules to keep per cell
    molecules_to_keep = set()

    for cell_bc, mol_list in CellDict.items():
        unique_mols = list(set(mol_list))  # remove duplicates if any
        num_mols = len(unique_mols)

        # Number of molecules to keep for this cell
        keep_count = int(num_mols * fraction)
        if keep_count < 1:
            # If fraction is small and num_mols is low, keep_count could be zero.
            # You can decide if you want to force at least 1.
            keep_count = 0
        
        # Randomly choose that many molecules
        chosen_molecules = random.sample(unique_mols, keep_count)
        molecules_to_keep.update(chosen_molecules)

    # 2) Gather all read names we want to keep
    read_names_to_keep = set()
    for mol_id in molecules_to_keep:
        read_names_to_keep.update(MoleculeDict[mol_id])

    # 3) Write only the chosen reads to the output BAM
    with pysam.AlignmentFile(out_bam, "wb", header=bam_header) as bam_out:
        for rname in read_names_to_keep:
            if rname in ReadPairDict:
                for read_obj in ReadPairDict[rname]:
                    bam_out.write(read_obj)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Downsample ReDeeM mito bam file at molecular level')
    parser.add_argument('-i', '--bam-file', type=str, help='input bam file')
    parser.add_argument('-b','--barcode-tag', default="BC", type=str, help='BAM tag storing the cell barcode (default: BC)')
    parser.add_argument('-f','--fraction', type=float, help='Fraction of molecule to take for each cell. Input a number between 0 and 1')
    parser.add_argument('-d','--directory', type=str, default="./", help='The directory where the bam file will be write out. Default is ./')

    args = parser.parse_args()

    bam_file = args.bam_file
    barcode_tag = args.barcode_tag
    fraction = args.fraction
    directory = args.directory
    if not (0.0 < fraction <= 1.0):
        raise ValueError("Fraction must be in the range (0, 1].")
    
    bam_basename = os.path.basename(bam_file)            # e.g. "barcodes.10.bam"
    bam_stem = os.path.splitext(bam_basename)[0]         # e.g. "barcodes.10"
    out_bam = f"{directory}{bam_stem}.ds{fraction}.bam"             # e.g. "barcodes.10.ds0.5.bam"

    # 1) Build dictionaries
    ReadPairDict, MoleculeDict, CellDict, bam_header = build_dictionaries(
        bam_file=bam_file,
        barcode_tag=barcode_tag
    )

    # 2) Downsample and write out
    downsample_by_molecule_in_memory(
        ReadPairDict=ReadPairDict,
        MoleculeDict=MoleculeDict,
        CellDict=CellDict,
        fraction=fraction,
        out_bam=out_bam,
        bam_header=bam_header
    )
    
print(out_bam)
directory = '/lab/solexa_weissman/cweng/Packages/REDEEM-V/example_data/test/'