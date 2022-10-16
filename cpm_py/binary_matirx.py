# %%
import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
from Funcs import *
# from cpmFunctions import *

# %%
## FC files
path = '../data/fc_individual/'
file = glob.glob(os.path.join(path,'*.txt'))

# get subject id
subj_list = [ os.path.split(f)[1].replace('.fc.txt', '') for f in file ] # to get the subj id from file name
subj_list = np.array(subj_list, dtype=str) # data structure conversion, list to arrary
subj_list

# %%
# read in FC matrices
path = '../data/fc_individual/'
all_fc_data = read_in_matrices(subj_list, file_suffix=None, data_dir=Path(path))

all_fc_data.shape
# 82 subj, 69751 edges
# edge number = n_nodes*(n_nodes-1)/2, 69751 in this case

# %%
all_fc_data.head()


# %%
# your_asymm_matrix = ... from some transformation


# %%
edge_frac = (your_asymm_matrix.sum(axis=0))/(your_asymm_matrix.shape[0])
# all_masks['pos'].sum(axis=0), column sum
# all_masks['pos'].shape[0], row number

edge_frac_square = sp.spatial.distance.squareform(edge_frac) # Convert a vector-form distance vector to a square-form distance matrix

# %%
## to mask
thresh = 1
mat_mask = edge_frac_square >= thresh

edge_frac_square[mat_mask] = 1
edge_frac_square[~mat_mask] = 0

edge_frac_square.shape



# %%
## save it
np.savetxt('binary_mat.csv', edge_frac_square, delimiter=',', fmt='%d') 

# %% [markdown]
# ### Get binary matrix without cpm
# 
# We only need `all_masks` to get the binary matrices. So how CPM works doesn't matter, we just need to create a object has similar data structure with `all masks`.
# 
# - load FC matrices with `read_in_matrices`, get `fc_df`
#     - dataframe with n rows and m columns (n=subj, m=edges)
# - whatever transformation
# - **KEY STEP**: get symmetric matrix
#     ```python
#         edge_frac = (fc_df.sum(axis=0))/(fc_df.shape[0]) # maybe need to convert to array
#         edge_frac_square = sp.spatial.distance.squareform(edge_frac) # Convert a vector-form distance vector to a square-form distance matrix
#     ```
# - mask: make it binary
#     - some as 0, some as 1


