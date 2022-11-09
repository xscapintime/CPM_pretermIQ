import os, glob
import numpy as np
import pandas as pd
import scipy as sp
from matplotlib import pyplot as plt

# files
# predicted table
files = sorted(glob.glob(os.path.join('./','*_all_masks.npy')))


# mask and convert to symmetrical matrices 
thresh = 1 # for all the folds

for f in files:
    all_masks = np.load(f, allow_pickle=True).item()

    for tail in all_masks.keys():
        edge_frac = (all_masks[tail].sum(axis=0))/(all_masks[tail].shape[0])
        print("For the {} tail, {} edges were selected in at least {}% of folds".format(tail, (edge_frac>=thresh).sum(), thresh*100))
        # Convert a vector-form distance vector to a square-form distance matrix
        edge_frac_square = sp.spatial.distance.squareform(edge_frac)
    
        # mask
        # selected edges as 1
        # masked edges as 0
        mat_mask = edge_frac_square >= thresh
        edge_frac_square[mat_mask] = 1
        edge_frac_square[~mat_mask] = 0
    
        # number of selected edges
        n_edges = edge_frac_square.sum()/2 # devide 2 is becuz now it's a symmetrical matrix
        print(f)
        print('{} edges were selected for {} network'.format(n_edges, tail))

        # export
        fn = tail + '-'+ os.path.split(f)[1].replace('_all_masks.npy', '') 
        np.savetxt(fn + '.csv', edge_frac_square, delimiter=',', fmt='%d') 
