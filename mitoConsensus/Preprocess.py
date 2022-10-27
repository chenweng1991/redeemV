import click
from multiprocessing import Pool
from mgatkHelp import *
from itertools import repeat

# ncores=12
# input="Merge_Mito.sub.bam"
# barcodes="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/Run210706L1/Qualified_CellBC_ATAC_2k0.2.tb"
output="TestFolder"
mito_genome="rCRS"
#script_dir="/lab/solexa_weissman/cweng/Packages/CW_mgatk"

@click.command()
@click.version_option()
@click.option('--input', '-i', default = ".", required=True, help='Input; either directory of singular .bam file; see documentation. REQUIRED.')
@click.option('--ncores', '-c', default = "detect", help='Number of cores to run the main job in parallel.')
@click.option('--barcodes', '-b', default = "",  help='File path to barcodes that will be extracted; useful only in `bcall` mode. If none supplied, mgatk will learn abundant barcodes from the bam file (threshold defined by the -mb tag).')
@click.option('--output', '-o', default="mgatk_out", help='Output directory for analysis required for `call` and `bcall`. Default = mgatk_out')
@click.option('--mito-genome', '-g', default = "rCRS", required=True, help='mitochondrial genome configuration. Choose hg19, hg38, mm10, (etc.) or a custom .fasta file; see documentation. Default = rCRS.')
@click.option('--barcode-tag', '-bt', default = "BC", required=True, help='tag name see documentation. Default = BC.')
@click.option('--script-dir', '-sd', default = "/lab/solexa_weissman/cweng/Packages/REDEEM-V/mitoConsensus/", required=True, help='')
# For testing purpose


def main(input,ncores,barcodes,output,mito_genome,barcode_tag,script_dir):
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


if __name__ == '__main__':
    main()
