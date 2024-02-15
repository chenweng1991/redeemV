#!/usr/bin/env python3
"""
Module sceasy
"""

__author__ = "Johnny Yu"
__version__ = "0.1.0"
__license__ = "MIT"


#####################
#####################
import pandas as pd
import copy
import os
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import scanpy as sc
#####################
#####################

#######################################################
####################################################### IO functions
#######################################################

def read_species(human=True):
    if human == True:
        return(sc.read_h5ad("DataBySpecies/human.anndata.h5ad"))
    else:
        return(sc.read_h5ad("DataBySpecies/mouse.anndata.h5ad"))

def save(adata,filename):
    adata.write(filename+".h5ad")

def read(filename):
    adata = sc.read_h5ad(filename+".h5ad")
    return(adata)
#######################################################
####################################################### Preliminary metadata functions
#######################################################
    
def overlay_meta(adata,LMOfile):
    ###get barcodes
    adata.obs['barcode'] = adata.obs.index.str[:-2]
    ###get dictionaries of metadata
    LMOdict = get_LMOfile(LMOfile)
    sig,call = get_pymulti()
    ###check indices
    adata = adata[adata.obs.barcode.isin(sig.keys())]
    adata = adata[adata.obs.barcode.isin(call.keys())]
    ###write into .obs
    adata.obs['sig'] = adata.obs.apply(lambda row: sig[row.barcode],axis=1) 
    adata.obs['call'] = adata.obs.apply(lambda row: call[row.barcode],axis=1)
    ###write sample name into .obs
    adata.obs['sample'] = adata.obs.apply(lambda row: LMOdict[row.call],axis=1)
    ###return
    return(adata)
    
def get_LMOfile(LMOfile):
    """ read in LMOfile and turn into dictionary. 
        using the multiseq barcode sequence as the keys. """
    bcsmulti = pd.read_csv(LMOfile,sep=',',index_col=0,header=None)
    bcsmulti.columns = ['multi']
    return(bcsmulti.to_dict()['multi'])

def get_pymulti():
    """ read in python multi seq file and convert to two dictionaries, one for significance and one for call.
        using the cell barcodes as the keys. """
    bcs = pd.read_csv('pymulti/pymulti__calls.tsv',sep='\t',index_col=0)
    bcs = bcs[['sig','call']]
    sig = bcs.sig.to_dict()
    call = bcs.call.to_dict()
    return(sig,call)

#######################################################
####################################################### Basic QC
#######################################################

def qcstats(adata,MT_prefix):
    ###stat
    adata.obs['n_counts'] = adata.X.sum(1)
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    adata.obs['n_genes'] = (adata.X > 0).sum(1)
    ###calculate # mitochondrial
    mito_genes = adata.var.index.str.startswith(MT_prefix)
    adata.obs['mt_frac'] = np.sum(
        adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
    ###return
    return(adata)

def qcplots(adata):
    #Sample quality plots
    t1 = sc.pl.violin(adata, 'n_genes', groupby='sample', size=2, log=True, cut=0)
    t1 = sc.pl.violin(adata, 'n_counts', groupby='sample', size=2, log=True, cut=0)
    t2 = sc.pl.violin(adata, 'mt_frac', groupby='sample')
    #Data quality summary plots
    colors = ['mt_frac','sample']
    for color in colors:
        print('plotting '+color+' below.')
        sc.pl.scatter(adata, 'n_counts', 'n_genes', color=color)
        sc.pl.scatter(adata[adata.obs['n_counts']<10000], 'n_counts', 'n_genes', color=color)
        sc.pl.scatter(adata[adata.obs['n_counts']<5000], 'n_counts', 'n_genes', color=color)
        sc.pl.scatter(adata[adata.obs['n_counts']<2500], 'n_counts', 'n_genes', color=color)
        
def qccheck(adata,MT_prefix='MT-'):
    adata = qcstats(adata,MT_prefix)
    qcplots(adata)
    return(adata)

def qccounts(adata,threshold=4000):
    #Thresholding decision: counts
    sns.distplot(adata.obs['n_counts'])
    plt.figure()
    sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']<threshold],bins=100)
    plt.figure()
    sns.distplot(adata.obs['n_counts'][adata.obs['n_counts']>threshold],bins=100)
    plt.figure()

def qcgenes(adata,threshold=1000):
    #Thresholding decision: genes
    sns.distplot(adata.obs['n_genes'], bins=100)
    plt.figure()
    sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']<threshold], bins=100)
    plt.figure()
    sns.distplot(adata.obs['n_genes'][adata.obs['n_genes']>threshold], bins=100)
    plt.figure()

def cellcycle(adata):
    adata.var_names_make_unique()
    #Score cell cycle and visualize the effect:
    cell_cycle_genes = [x.strip() for x in open('scEasyMode/regev_lab_cell_cycle_genes.txt')]
    s_genes = cell_cycle_genes[:43]
    g2m_genes = cell_cycle_genes[43:]
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    return(adata)

def qc_all(adata):
    ###qc
    adata = qccheck(adata)
    adata = cellcycle(adata)
    qccounts(adata)
    qcgenes(adata)
    return(adata)

#######################################################
####################################################### Basic Filtering
#######################################################

def annotate_mito(adata,mt_thresh):
    adata.obs['dead'] = adata.obs.apply(lambda row: row.mt_frac>mt_thresh,axis=1)
    return(adata)

def filter_cells(adata,mt_thresh,min_counts=0,max_counts=10000000,min_genes=0):
    ###Filter cells according to identified QC thresholds:
    print('Total number of cells: {:d}'.format(adata.n_obs))
    sc.pp.filter_cells(adata, min_counts = min_counts)
    print('Number of cells after min count filter: {:d}'.format(adata.n_obs))
    sc.pp.filter_cells(adata, max_counts = max_counts)
    print('Number of cells after max count filter: {:d}'.format(adata.n_obs))
    adata = adata[adata.obs['mt_frac'] < mt_thresh]
    print('Number of cells after MT filter: {:d}'.format(adata.n_obs))
    sc.pp.filter_cells(adata, min_genes = min_genes)
    print('Number of cells after gene filter: {:d}'.format(adata.n_obs))
    return(adata)

def filter_genes(adata,min_cells=20):
    print('Total number of genes: {:d}'.format(adata.n_vars))
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print('Number of genes after cell filter: {:d}'.format(adata.n_vars))
    return(adata)

def normalize(adata):
    ###Keep the count data in a counts layer
    adata.layers["counts"] = adata.X.copy()
    ###normalize
    sc.pp.normalize_total(adata,exclude_highly_expressed=True)
    sc.pp.log1p(adata)
    ##Store the full data set in 'raw' as log-normalised data for statistical testing
    adata.raw = adata
    ##return
    return(adata)

def filters(adata,mt_thresh,min_cells,sig_pct,min_counts,min_genes):
    ###annotate dead cells in full data set
    adata = annotate_mito(adata,mt_thresh)
    ###filter by sig from multiseq calls
    adata = adata[adata.obs.sig>np.percentile(adata.obs.sig,sig_pct)]
    ###filter mitochondrial high genes out and create clean dataset
    clean = filter_cells(adata,mt_thresh,min_counts=min_counts,min_genes=min_genes)
    ###filter genes
    adata = filter_genes(adata,min_cells)
    clean = filter_genes(clean,min_cells)
    ###normalize
    adata = normalize(adata)
    clean = normalize(clean)
    ###get hvgs
    adata = define_hvgs(adata)
    clean = define_hvgs(clean)
    ###get shapes
    print(adata.shape)
    print(clean.shape)
    ###return
    return(adata,clean)

def define_hvgs(adata,n_genes=3000):
    sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=n_genes)
    print('\n','Number of highly variable genes: {:d}'.format(np.sum(adata.var['highly_variable'])))
    sc.pl.highly_variable_genes(adata)
    return(adata)

#######################################################
####################################################### Clustering and visualization
#######################################################

def visualize(adata,covariates=['n_counts','n_genes','mt_frac','phase','sample','louvain','dead','sig']):
    ###Calculate the visualizations
    sc.pp.pca(adata, n_comps=50, use_highly_variable=True, svd_solver='arpack')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata,random_state=10)
    sc.tl.diffmap(adata)
    sc.tl.draw_graph(adata)
    ###cluster
    adata = cluster(adata)
    ###plot
    for covariate in covariates:
        sc.pl.pca_scatter(adata, color=covariate)
        sc.pl.umap(adata, color=covariate)
        sc.pl.diffmap(adata, color=covariate, components=['1,2','1,3'])
        sc.pl.draw_graph(adata, color=covariate)
    ###return
    return(adata)

def regress(adata,factors):
    sc.pp.regress_out(adata, factors)
    sc.pp.scale(adata)
    return(adata)

def cluster(adata):
    sc.tl.louvain(adata, resolution=0.5, key_added='louvain', random_state=10)
    return(adata)

def densitymap(adata,sample_key='sample'):
    sc.tl.embedding_density(adata, basis='umap', groupby=sample_key)
    sc.pl.embedding_density(adata, basis='umap', key='umap_density_'+sample_key)

        