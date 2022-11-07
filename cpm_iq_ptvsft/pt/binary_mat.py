import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
from scipy.stats import kstest
import statistics
from matplotlib import pyplot as plt
from cpmFunctions import *


all_masks = np.load('WISC_FULL_CS_8yo_pearson_fold_80_pval_0.01_all_masks.npy', allow_pickle=True).item()


thresh = 1

for tail in all_masks.keys():
    edge_frac = (all_masks[tail].sum(axis=0))/(all_masks[tail].shape[0])
    
    # Convert a vector-form distance vector to a square-form distance matrix
    edge_frac_square = sp.spatial.distance.squareform(edge_frac)
    
    # mask
    # selected edges as 1
    # masked edges as 0
    mat_mask = edge_frac_square >= thresh
    edge_frac_square[mat_mask] = 1
    edge_frac_square[~mat_mask] = 0
    
    # number of selected edges
    n_edges = edge_frac_square.sum()
    print('{} edges were selected for {} network'.format(n_edges, tail))

    # export
    fn = tail + 'WISC_FULL_CS_8yo_pearson_fold_80_pval_0.01'
    np.savetxt(fn + '.csv', edge_frac_square, delimiter=',', fmt='%d') 
