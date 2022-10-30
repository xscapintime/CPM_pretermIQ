import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
from scipy.stats import kstest
import statistics
from matplotlib import pyplot as plt
from cpmFunctions import *


## check FC file names
path = '../data/fc_individual/'
file = glob.glob(os.path.join(path,'*.txt'))
# file

# get subject id
subj_list = [ os.path.split(f)[1].replace('.fc.txt', '') for f in file ]
subj_list = np.array(subj_list, dtype=str)
subj_list
# remove 7 for external validation
mask = np.ones(len(subj_list), dtype=bool)
mask[[75,76,77,78,79,80,81]] = False
subj_list = subj_list[mask]
len(subj_list)


## beahve data
# all_behav_data = pd.read_csv('../data/id_vars_fin.csv', dtype={'Eprime_ID': str})
# all_behav_data.dtypes

# # filter by AP ID
# all_behav_data = all_behav_data[all_behav_data['AP_ID'].isin(subj_list)]
# all_behav_data = all_behav_data.set_index('AP_ID')


# # for modeling, no NA allowed
# behav_data = all_behav_data.iloc[:,6:].dropna()
# print(behav_data.shape)


## PCs
# 8 yo vars only
all_pca_data = pd.read_csv('../data/var6_8yo_sbj82_imp_pca.csv', index_col=0)
all_pca_data.dtypes

# filter
pca_data = all_pca_data[all_pca_data.index.isin(subj_list)]
print(pca_data.shape)


## CPM
# read in FC matrices
condition = 'fc' ## actually no need, should be something like REST or EMO to distinguish different matrices
all_fc_data = read_in_matrices(subj_list, file_suffix=condition, data_dir=Path(path))
print(all_fc_data.shape)
fc_data = all_fc_data.loc[pca_data.index]
print(fc_data.shape) # edge number = n_nodes*(n_nodes-1)/2, 69751 in this case



## load covariates data
cova_data = pd.read_csv('../data/pt_sbj82withfc_newvars_imp_fd_pca_merged.csv', index_col=1)
cova_data = cova_data.loc[pca_data.index]
cova_data = cova_data[['sex', 'age8', 'Neonatal Sickness', 'IMD Score', 'fd.m']]#.dropna() ##
cova_data.shape



## viz edges
coords = pd.read_csv("../data/coords/MMP_MNI_374_UCCHILD_coords.txt", index_col=None, header=None, sep=" ")
print(coords.shape)

# paras for CPM
cpm_kwargs = {'r_thresh': 0.35, 'corr_type': 'spearman', 'verbose': False} ## use spearman if the distribution is skewed
covar = ['sex', 'age8', 'Neonatal Sickness', 'IMD Score', 'fd.m']


# plots, would be one of the random cases
# set k as the number of subjects, Leave-One-Out Cross-Validation (LOOCV)


### k as number of subjects, LOOCV
k = fc_data.shape[0]
corr_type = 'spearman'

for behav in pca_data.columns[:2]:
    print(behav)

    behav_obs_pred, all_masks, corr = cpm_wrapper_seed_part(fc_data, pca_data, behav=behav, cova_data=cova_data, covar=covar, k=k, seed=202209, **cpm_kwargs)
    ## count selected edges
    print('{:.2f} pos edges passed the threshold at all folds: '.format(((all_masks['pos'].sum(axis=0)/k) >= 1).sum()))
    print('{:.2f} neg edges passed the threshold at all folds: '.format(((all_masks['neg'].sum(axis=0)/k) >= 1).sum()))
    
    ## export pred table
    fn = behav + '_8yo_' + corr_type + '_fold_' + str(k) + '_loocvpearson_pred.csv'
    behav_obs_pred.to_csv(fn)

    ## plot pred vs obs
    g = plot_predictions(behav_obs_pred)
    g.figure.tight_layout()
    plt.savefig(os.path.join('models', behav + '_8yo_' + corr_type + '_fold_' + str(k) + '_loocvpearson_model.png'))
    plt.close()

    ## save all_mask
    np.save(fn + '_all_masks.npy',  all_masks)    

    ## edges
    for tail in all_masks.keys():
        edfn = tail + '_' + behav + '_8yo_' + corr_type + '_fold_' + str(k)
        
        edge_frac = (all_masks[tail].sum(axis=0))/(all_masks[tail].shape[0])
        edge_frac_square = sp.spatial.distance.squareform(edge_frac)
        
        # export pos and neg matrices
        # edmat = edfn + '_loocvpearson_bn.csv'
        # np.savetxt(edmat, edge_frac_square, delimiter=',', fmt='%1.0f')
    
        # plot edges
        # all_masks(pos/neg matrices) is a binary array, k fold * n total edges, in this case is 10,69751
        if tail == 'pos':
            color = 'red'
        else:
            color = 'blue'
        
        g = plot_consistent_edges(all_masks, tail, thresh = 0.8, color = color, node_coords = coords)
        plt.savefig(os.path.join('edges', edfn + '_loocvpearson_edges.png'))
        plt.close()
    
