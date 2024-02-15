import sys
args = sys.argv
REDEEM_V=args[1]
sys.path.append(REDEEM_V+"/JohnnyCellHash/")
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from collections import Counter
from scEasyMode import pymulti
import os
import pickle



sys.setrecursionlimit(100000)
cwd = os.getcwd()
Name=args[2]
lib10x = Name
libmulti = Name

fq1=args[3]
fq2=args[4]

#fqpath="/lab/solexa_weissman/cweng/Projects/MitoTracing_Velocity/SecondaryAnalysis/ConditionTest2_Donor04_BMMC2/FASTQ/Hash/"
####define files
R1 = cwd+"/"+fq1
R2 = cwd+"/"+fq2
bc10xfile = "barcodes.tsv.16bp"
bcfile = 'HashDic.csv'
####define metadata
v10x = 'v3.1'
expname = Name
sampname = Name
####define length of sequences
len_10x = 16
len_umi = 12
len_multi = 15
####define multiseq barcodes file
bcsmulti = pd.read_csv(bcfile,sep=',',index_col=1,header=None)
bcsmulti.columns = ['multi']
bcsmulti = bcsmulti['multi'].tolist()
####define 10x barcodes whitelist
bcs10x = pd.read_csv(bc10xfile,sep='\t',header=None)
bcs10x = bcs10x[0].tolist()
####Run pymult
fig=pymulti.pymulti(R1,R2,bcsmulti,bcs10x,len_10x=len_10x,len_multi=len_multi,len_umi=len_umi,split=True,hamming=True,med_factor=1,median_only=True,sampname=lib10x)
fig.savefig("plot.pdf")

#### Write out umi count and umicount histogram
readtable = pd.DataFrame(pickle.load(open(cwd+'/pymulti/'+Name+"_reads.p", "rb" ) ))
readtable.columns = ['cell','umi','multi']
readtable.set_index('cell')
filtd = pymulti.filter_readtable(readtable,bcsmulti,bcs10x)
umicount=filtd.groupby("cell")["multi"].value_counts().unstack()
umicount["total"]=umicount.sum(axis=1)
umicount.to_csv("pymulti/umicount.csv")
plot_umicount=sns.histplot(np.log10(umicount['total']), kde=False, bins=20)
plot_umicount.figure.savefig("pymulti/umicount_hist.pdf")