"""
Test the optimal threshold regarding r value, p value and percentile
Find the fastest way
"""

import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
import statistics
from matplotlib import pyplot as plt
from cpmFinctions import *
import time

## check FC file names
path = '../data/fc_individual/'
file = glob.glob(os.path.join(path,'*.txt'))
file

# get subject id
subj_list = [ os.path.split(f)[1].replace('.fc.txt', '') for f in file ]
subj_list = np.array(subj_list, dtype=str)
subj_list


## beahve data
all_behav_data = pd.read_csv('../data/id_vars_fin.csv', dtype={'Eprime_ID': str})
all_behav_data.dtypes

# filter by AP ID
all_behav_data = all_behav_data[all_behav_data['AP_ID'].isin(subj_list)]
all_behav_data = all_behav_data.set_index('AP_ID')

# for modeling, no NA allowed
behav_data = all_behav_data.iloc[:,6:].dropna()
print(behav_data.shape)

## PCs
all_pca_data = pd.read_csv('../data/var19_pca.csv', index_col=0)
all_pca_data.dtypes

# filter
all_pca_data = all_pca_data[all_pca_data.index.isin(subj_list)]
print(all_pca_data.shape)


## CPM
# read in FC matrices
condition = 'fc' ## actually no need, should be like REST or EMO
all_fc_data = read_in_matrices(subj_list, file_suffix=condition, data_dir=Path(path))
fc_data = all_fc_data.loc[behav_data.index]

## test time
# r val, no intermediate pandas dataframe, might be faster
cpm_kwargs = {'r_thresh': 0.3, 'corr_type': 'pearson', 'verbose': False} ## use spearman if the distribution is skewed

start_time = time.time()
for behav in all_pca_data.columns[0:1]: ## PC1 as the example
    print(behav)
    behav_obs_pred, all_masks, corr = cpm_wrapper(fc_data, all_pca_data, behav=behav, **cpm_kwargs)
print("--- %s seconds, select by rval---" %(time.time() - start_time))
# --- 3.880763292312622 seconds, select by rval---


# p val
cpm_kwargs = {'p_thresh': 0.01, 'corr_type': 'pearson', 'verbose': False} ## use spearman if the distribution is skewed

start_time = time.time()
for behav in all_pca_data.columns[0:1]: ## PC1 as the example
    print(behav)
    behav_obs_pred, all_masks, corr, pval = cpm_wrapper_pval(fc_data, all_pca_data, behav=behav, **cpm_kwargs)
print("--- %s seconds, select by pval ---" %(time.time() - start_time))
# --- 191.23768043518066 seconds, select by pval ---
# pandas dataframe made it slow?

# change pval & corr dataframe to dict
# not helping
# --- 193.67261004447937 seconds, select by pval ---

#### loop to iterate the columns made it slow ####



# top x% strong edges
cpm_kwargs = {'perc_thresh': 1, 'corr_type': 'pearson', 'verbose': False} ## use spearman if the distribution is skewed

start_time = time.time()
for behav in all_pca_data.columns[0:1]: ## PC1 as the example
    print(behav)
    behav_obs_pred, all_masks, corr, pval = cpm_wrapper_alt(fc_data, all_pca_data, behav=behav, **cpm_kwargs)
print("--- %s seconds, select by top x ---" %(time.time() - start_time))
# the time will not be different with selecting by p val


### so the original one without doing loops between columns is the fastest one ###



# test equivalent to ...
# use the second result


# corr dist
sum(abs(np.array(corr)) > 0.3) # number of delected edges
# 2538, too many
g = sns.displot(np.array(corr)[abs(np.array(corr)) > 0.3])
plt.savefig(os.path.join('dist', 'r0.3edges_dist.png'))
plt.close()

sum(abs(np.array(corr)) > 0.38) # number of delected edges
# 452, equivalent to p val < 0.01
g = sns.displot(np.array(corr)[abs(np.array(corr)) > 0.38])
plt.savefig(os.path.join('dist', 'r0.38edges_dist.png'))
plt.close()


# pval as threshold, cor dist
sum(np.array(pval) < 0.01) # number of delected edges
# 492
g = sns.displot(np.array(corr)[np.array(pval) < 0.01])
plt.savefig(os.path.join('dist', 'p0.01edges_dist.png'))
plt.close()


# percentile as thredhold, cor dist
sum(abs(np.array(corr)) > np.percentile(abs(np.array(corr)), 100-0.7))
# 489
g = sns.displot(np.array(corr)[abs(np.array(corr)) > np.percentile(abs(np.array(corr)), 100-0.7)])
plt.savefig(os.path.join('dist', 'top0.7percedges_dist.png'))
plt.close()


# for glasser plot, would be good to have ~ 50 edges with fold threshold 0.8
# for one edge, sum the value of it in each fold
# If a edge was selected every time (10 folds here), it would be 10/10 
edge_frac = all_masks['pos'].sum(axis=0)/10
n_passed = (edge_frac >= 0.9).sum()
print(n_passed)
# 89

plot_consistent_edges(all_masks, "pos", thresh = 0.9, color = 'red', node_coords = coords)
# 89
plot_consistent_edges(all_masks, "neg", thresh = 0.9, color = 'blue', node_coords = coords)
# 48
