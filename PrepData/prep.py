"""
This script is used prepare data from ReDeeM.
According the configuration file and the Data.summary in the fastq folder(s)
https://github.com/chenweng1991/redeemV/wiki/Organize-ReDeeM-full-data

Usage:
    python prep.py ini

Arguments:
    ini: This is the only argument. It is the path to the configuration file

Examples:
    python REDEEM-V/PrepData/prep.py prepdata.ini > PrepData.log

Author:
    Chen Weng <cweng@wi.mit.edu>
"""

import pandas as pd
import configparser
import argparse
import sys
import os

# Parse argument
parser = argparse.ArgumentParser(description='This script is used prepare data from ReDeeM.According the configuration file and the Data.summary in the fastq folder(s).')
parser.add_argument('ini', help='It is the path to the configuration file')

args = parser.parse_args()

config = configparser.ConfigParser()
config.read(args.ini)

# Access input folders from the config file
fq_folders = config.get('Input', 'fq_folders').split(', ')
dfs=[]
for folder in fq_folders:
    summary=pd.read_csv(folder+'/Data.summary',sep=',',header=None)
    summary.iloc[:,0]=folder+"/"+summary.iloc[:,0]
    dfs.append(summary)
summary=pd.concat(dfs,axis=0)

# Define samples and assays
ExpectMito = config.getboolean('Parameters', 'mitofq')
samples=summary.iloc[:,1].unique()
assays=summary.iloc[:,2].unique()
MitoAssay='Mito' in assays
if ExpectMito and MitoAssay:
    print("This analysis will have a mtDNA specific set of fastqs!")
elif not ExpectMito and not MitoAssay:
    print("This analysis do not have a mtDNA specific set of fastqs, use ATAC fastqs!")
else:
    print("Error: do you have a mtDNA specific set of fastqs? Clarify in the ini file")

# Create folders  
Create=sys.argv[0].replace("prep.py","Create.sh")  
out = config.get('output', 'out')
for sample in samples:
    shell_script=Create+' '+out+"/"+sample
    os.system(shell_script)
    
# softlink or perform fastx_trimmer for each 
parallel= " &" if config.getboolean('Parameters', 'parallel')  else " "
## Write ATAC, RNA and Mito (if Mito is in summary table, and mitofq=True)
for sample in samples:
    for assay in assays:
        summary_sample_assay=summary[(summary.iloc[:,1]==sample) & (summary.iloc[:,2]==assay)]
        infiles=summary_sample_assay.iloc[:,0].str.split("/").str[-1]
        summary_sample_assay.loc[:,"outfile"]=sample+"_"+infiles.str.extract(r'(S\d.*)').squeeze()
        outfolder=out+'/'+sample+'/FASTQ/'+assay+"/"
        for index, row in summary_sample_assay.iterrows():
            if row[3] =="notrim":
                shell_script="ln -s "+row[0]+" "+outfolder+row["outfile"]
            else:
                shell_script="zcat " +row[0]+" | fastx_trimmer -f 1 -l " +row[3]+" -z -o "+outfolder+row["outfile"]+parallel
            print(shell_script)
            os.system(shell_script)

## Write Mito if Mito not in summary table
if not 'Mito' in assays:
    for sample in samples:
        summary_sample_assay=summary[(summary.iloc[:,1]==sample) & (summary.iloc[:,2]=="ATAC")]
        infiles=summary_sample_assay.iloc[:,0].str.split("/").str[-1]
        summary_sample_assay.loc[:,"outfile"]=sample+"_"+infiles.str.extract(r'(S\d.*)').squeeze()
        outfolder=out+'/'+sample+'/FASTQ/Mito/'
        for index, row in summary_sample_assay.iterrows():
            shell_script="ln -s "+row[0]+" "+outfolder+row["outfile"]
            print(shell_script)
            os.system(shell_script)


