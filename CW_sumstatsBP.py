import pysam
import click
import sys
import pickle
import numpy as np
import pandas as pd
from progress.bar import Bar
from collections import defaultdict


@click.command()
@click.version_option()
@click.option('--out-dir', '-o', default = ".", required=True, help='The folder of the whole analysis. REQUIRED.')
@click.option('--sample', '-s', default = ".", required=True, help='The prefix, should be barcodes.# check the prefix of /temp/barcoded_bams/ for example. REQUIRED.')
@click.option('--barcode-tag', '-bt', default = "BC", help='The tag shown in the bam file used as cell barcode; Default is CB')
@click.option('--base-qual', '-bq', default = 10, help='The threshold of minimum base quality that can be considered; Default is 10')
@click.option('--pcr-qual', '-pq', default = 0.8, help='The percentage of PCR copies that support a mutation, Default is 80%, aka. for one molecule 9 out of 10 PCR copies showing the same mutation')

def main(out_dir, sample, barcode_tag, base_qual, pcr_qual):
    barcodes_file =out_dir + "/temp/barcode_files/"+sample+".txt"
    bam_file =out_dir + "/temp/barcoded_bams/"+sample+".bam"
    mito_ref_file =out_dir + "/final/chrM_refAllele.txt"
    out_pre =out_dir + "/temp/sparse_matrices/"+sample
    dna_letters = ['A','C','G','T']
    # Import all cell barcodes
    with open(barcodes_file) as barcode_file_handle:
        content = barcode_file_handle.readlines()

    bcs = [x.strip() for x in content]

    # Read bam file make a dictionary a dictionary from Rname->read1,read2
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
        CellBC=fwd_read.get_tag(barcode_tag)
        Molecule=CellBC+'_'+str(fwd_read.pos)+'_'+str(fwd_read.pos+fwd_read.tlen)
        MoleculeDict[Molecule].append(read_name)

    ## Predefine the SuperMatrix
    mito_ref=pd.read_table(mito_ref_file,names=["pos","base"])
    max_bp=mito_ref.shape[0]
    SuperMatrix={key:np.zeros((max_bp,4)) for key in bcs}
    Quality_report={key:list() for key in bcs}

    ## Loop through each molecule for genotyping
    MolecularN=len(MoleculeDict)
    bar = Bar('Genotyping', max=MolecularN)
    for index, m in enumerate(MoleculeDict):
        CellBC=m.split('_')[0]
        FragmentSize=int(m.split('_')[2])-int(m.split('_')[1])
        MoleculeGenotypes=np.zeros((max_bp,4)) # count 'A' 'C' 'G' 'T' in order
        overlap=np.empty(0,dtype=int) # to summzrize the overlaped sizes for each copy of a molecule
        overlap_bases=0
        Conflict_strand_bases=0
        Conflict_pcr_bases=0
        LowQ_bases=0
        for read_pair in MoleculeDict[m]:
            seq_0=ReadPairDict[read_pair][0].seq
            seq_1=ReadPairDict[read_pair][1].seq
            quality_0=ReadPairDict[read_pair][0].query_qualities
            quality_1=ReadPairDict[read_pair][1].query_qualities
            align_qual_0 = ReadPairDict[read_pair][0].mapping_quality
            align_qual_1 = ReadPairDict[read_pair][1].mapping_quality
            # Get the overlaped position array for both strand
            pos_array_0=np.asarray(ReadPairDict[read_pair][0].get_aligned_pairs(True))
            pos_array_1=np.asarray(ReadPairDict[read_pair][1].get_aligned_pairs(True))
            pos_array_overlap=np.intersect1d(pos_array_0[:,1],pos_array_1[:,1])
            pos_array_overlap_0=pos_array_0[np.isin(pos_array_0[:,1],pos_array_overlap)]
            pos_array_overlap_1=pos_array_1[np.isin(pos_array_1[:,1],pos_array_overlap)]
            overlap=np.append(overlap,len(pos_array_overlap))
            ## Identify reliable genotype that is concensus between the overlaped R1 and R2
            for base_0,base_1 in zip(pos_array_overlap_0,pos_array_overlap_1):
                overlap_bases+=1
                if(quality_0[base_0[0]]>base_qual and quality_1[base_1[0]]>base_qual):
                    if (seq_0[base_0[0]]==seq_1[base_1[0]]):
                        MoleculeGenotypes[base_0[1],dna_letters.index(seq_0[base_0[0]])]+=1
                    else:
                        Conflict_strand_bases+=1
                else:
                    LowQ_bases+=1
        ## Summarize the concensus genotype for the current molecule
        Molecule_con=MoleculeGenotypes.transpose()/MoleculeGenotypes.sum(axis=1)
        Molecule_con=Molecule_con.transpose()
        Molecule_con[np.isnan(Molecule_con)]=0
        Conflict_pcr_bases=np.count_nonzero(np.logical_and(np.amax(Molecule_con,axis=1)>0,np.amax(Molecule_con,axis=1)<pcr_qual))
        Molecule_con[Molecule_con<pcr_qual]=0
        Molecule_con[Molecule_con>pcr_qual]=1
        SuperMatrix[CellBC]=SuperMatrix[CellBC]+Molecule_con
        ## Summarize the molecule report
        if (np.mean(overlap)==0):
            Molecule_report=[m,len(MoleculeDict[m]),FragmentSize,0,0,0,0]
        else:
            Molecule_report=[m,len(MoleculeDict[m]),FragmentSize,np.mean(overlap),LowQ_bases/overlap_bases,Conflict_strand_bases/overlap_bases,Conflict_pcr_bases]
            Quality_report[CellBC].append(Molecule_report)
        bar.next()


    bar.finish()
    # Write out the quality report and the SuperMatrix
    out_file_mtx = out_pre + "."+"supermtx.data"
    out_file_qc = out_pre + "."+"QC.data"

    with open(out_file_mtx,"wb") as file_handle_mtx:
        pickle.dump(SuperMatrix, file_handle_mtx)

    with open(out_file_qc,"wb") as file_handle_qc:
        pickle.dump(Quality_report, file_handle_qc)


if __name__ == '__main__':
    main()
