#!/usr/bin/env python3
"""
Module MouseFilter
"""

__author__ = "Johnny Yu"
__version__ = "0.1.0"
__license__ = "MIT"


#####################
#####################
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import igraph
import seaborn as sns
import seaborn as sb
import os
#####################
#####################

def main(mousefile,humanfile=None,cutoff_top=0.5,cutoff_bottom=0.5,objects=False):
    """ Returns truemouse and truehuman Anndata objects filtered by cutoff of % mouse genes 
        mousefile = the raw anndata object that has been aligned to mm10hg19 (or something equivalent)
        humanfile = the raw anndata object that has been aligned to human only
        cutoff = the threshold of % mouse genes that is used to group cells as mouse or human
    """
    if objects == True:
        mouse = mousefile
        human = humanfile
    else:
        ###make sure inputs are valid h5 files
        check_inputs(mousefile,humanfile,cutoff_top,cutoff_bottom)
        ###read files into anndata objects
        mouse,human = read_files(mousefile,humanfile)
    ###format mouse dataframe
    process_mouse(mouse)
    ###plot mouse reads
    plot_mouse(mouse)
    ###split dataset into mouse and human
    truemouse,truehuman = filter_mouse(mouse,cutoff_top,cutoff_bottom)
    truehuman = subset_human(human,truehuman)
    ###return anndata objects for each species
    write_anndata(truemouse,truehuman)
    return(truemouse,truehuman)
    
def check_inputs(mousefile,humanfile,cutoff_top,cutoff_bottom):
    """ Make sure you have two matrices: one from mouse, one from human. mouse should be aligned to mm10hg19 and human should be hg38 """
    if cutoff_top > 1 or cutoff_top < 0:
        raise NameError('cutoff has to be within 0 < cutoff < 1')
    if cutoff_bottom > 1 or cutoff_bottom < 0:
        raise NameError('cutoff has to be within 0 < cutoff < 1')
    try:
        sc.read_10x_h5(mousefile)
        if humanfile != None:
            sc.read_10x_h5(humanfile)
    except:
        print('the files are not readable as h5.')
    
def read_files(mousefile,humanfile):
    """ Read into anndata objects and return """
    mouse = sc.read_10x_h5(mousefile)
    if humanfile != None:
        human = sc.read_10x_h5(humanfile)
    else:
        human = None
    return(mouse,human)

def process_mouse(mouse):
    """ Format and calculate the % mouse cells with mouse reads for the mouse anndata object prior to filtering """
    mouse.obs['barcode'] = mouse.obs.index.str[:-2]
    mouse.var_names_make_unique()
    mouse_genes = [name for name in mouse.var_names if name.startswith('mm10')]
    ### for each cell compute fraction of counts in mouse genes vs. all genes
    mouse.obs['percent_mouse'] = np.sum(
        mouse[:, mouse_genes].X, axis=1).A1 / np.sum(mouse.X, axis=1).A1
    ### add the total counts per cell as observations-annotation 
    mouse.obs['n_counts'] = np.sum(mouse.X, axis=1).A1

def plot_mouse(mouse):
    """ Plot the distribution of cells by % mouse genes """
    plt.hist(mouse.obs.percent_mouse,bins=100)
    plt.title('distribution of cells by % mouse genes')
    plt.xlim(-0.1,1.1)
    plt.show()

def filter_mouse(mouse,cutoff_top,cutoff_bottom):
    """ Saves two lists, one of mouse cells and one of human cells with the cell barcodes """
    ####filter by cutoff
    truemouse = mouse[mouse.obs['percent_mouse']>cutoff_top]
    truehuman = mouse[mouse.obs['percent_mouse']<cutoff_bottom]
    ####format genes
    keep_genes = list(truemouse.var[truemouse.var.genome == 'mm10'].index)
    truemouse = truemouse[:,keep_genes]
    ##
    keep_genes = list(truehuman.var[truehuman.var.genome == 'hg19'].index)
    truehuman = truehuman[:,keep_genes]
    ####get cell names
    mousecells = truemouse.obs.index.tolist()
    humancells = truehuman.obs.index.tolist()
    ####write cell names
    os.system('mkdir CellsBySpecies')
    with open('CellsBySpecies/mousecells.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % cell for cell in mousecells)
    with open('CellsBySpecies/humancells.txt', 'w') as filehandle:
        filehandle.writelines("%s\n" % cell for cell in humancells)      
    ####return anndata objects
    return(truemouse,truehuman)

def subset_human(human,truehuman):
    if human == None:
        return(truehuman)
    else:
        human = human[human.obs.index.isin(truehuman.obs.index)]
        return(human)

def write_anndata(truemouse,truehuman):
    os.system('mkdir DataBySpecies')
    truemouse.write("DataBySpecies/mouse.anndata.h5ad")
    truehuman.write("DataBySpecies/human.anndata.h5ad")
    
if __name__ == "__main__":
    """ This is executed when run from the command line """
    main(mousefile,humanfile,cutoff)
    
