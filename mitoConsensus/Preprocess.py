import click
from multiprocessing import Pool
from mgatkHelp import *
from itertools import repeat
""" Preprocess.py
This script split the mtDNA bam files for in-paralell variant calling

Syntax: python Preprocess.py -i -c -b -o -g -bt -sd
Options: --input or -i: uniqmapped.bam files generated by MultiomeATAC_mito.sh, REQUIRED.
	 --ncores or -c: Number of cores to run the main job in parallel. REQUIRED.
	 --barcodes or -b: The csv file indicating the qualified cell names. barcodes.tsv generated by RNAbc2ATAC.R, required
  	 --output or -o: The name of folder to be created that will contain the analyzed results, default=Out_mitoConsensus
	 --mito-genome or -gmitochondrial genome configuration. Choose hg19, hg38, mm10, (etc.) or a custom .fasta file,  Default = rCRS
     --barcode-tag or -bt  tag name in bam file, 10X is BC. Default = BC
     --script-dir or -sd  default = "/lab/solexa_weissman/cweng/Packages/REDEEM-V/mitoConsensus/ The path to REDEEM-V/mitoConsensus/

Acknowledgment: Inherited some of the codes from https://github.com/caleblareau/mgatk by Caleb Lareau et al. 
"""

@click.command()
@click.version_option()
@click.option('--input', '-i', default = ".", required=True, help='Uniqmapped.bam files generated by MultiomeATAC_mito.sh, REQUIRED')
@click.option('--ncores', '-c', default = 24, help='Number of cores to run the main job in parallel. default')
@click.option('--barcodes', '-b', default = ".",  help='The csv file indicating the qualified cell names. barcodes.tsv generated by RNAbc2ATAC.R, required')
@click.option('--output', '-o', default="Out_mitoConsensus", help='The name of folder to be created that will that will contail the analyzed results, default=Out_mitoConsensus')
@click.option('--mito-genome', '-g', default = "rCRS", required=True, help='mitochondrial genome configuration. Choose hg19, hg38, mm10, (etc.) or a custom .fasta file,  Default = rCRS.')
@click.option('--barcode-tag', '-bt', default = "BC", required=True, help='tag name in bam file, 10X is BC. Default = BC.')
@click.option('--script-dir', '-sd', default = "/lab/solexa_weissman/cweng/Packages/REDEEM-V/mitoConsensus/", required=True, help='The path to REDEEM-V/mitoConsensus/')


def Preprocess(input,ncores,barcodes,output,mito_genome,barcode_tag,script_dir):
    """
    This script split the mtDNA bam files for in-paralell variant calling
    
    Syntax: python Preprocess.py -i -c -b -o -g -bt -sd

    Acknowledgment: Inherited some of the codes from https://github.com/caleblareau/mgatk by Caleb Lareau et al. 
    """
    umi_barcode="XX"
    of = output; tf = of + "/temp"; bcbd = tf + "/barcoded_bams" # bcdb = barcoded bam directory
    folders = [of, tf, bcbd, of + "/final"] #,tf+"/sparse_matrices"
    mkfolderout = [make_folder(x) for x in folders]

    # Handle fasta requirements
    rawsg = os.popen('ls ' + script_dir + "/bin/anno/fasta/*.fasta").read().strip().split("\n")
    supported_genomes = [x.replace(script_dir + "/bin/anno/fasta/", "").replace(".fasta", "") for x in rawsg]
    fastaf, mito_chr, mito_length = handle_fasta_inference(mito_genome, supported_genomes, script_dir, "tenx", of)
    idxs = pysam.idxstats(input).split("\n")

    # Split Cell barcode files based on the cores to use
    barcode_files = split_barcodes_file(barcodes, math.ceil(file_len(barcodes)/int(ncores)), output)
    samples = [os.path.basename(os.path.splitext(sample)[0]) for sample in barcode_files]
    samplebams = [of + "/temp/barcoded_bams/" + sample + ".bam" for sample in samples]

    # Split bam files in a parallel manner
    pool = Pool(processes=int(ncores))
    pmblah = pool.starmap(split_chunk_file, zip(barcode_files, repeat(script_dir), repeat(input), repeat(bcbd), repeat(barcode_tag), repeat(mito_chr),repeat(umi_barcode)))
    pool.close()
    print("Splitted bam into "+str(ncores)+" bam files")


    # Print the help document
    click.echo()


if __name__ == '__main__':
    Preprocess()
