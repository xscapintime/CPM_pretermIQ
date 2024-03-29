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
path = '../../data/fc_individual_pt/'
file = sorted(glob.glob(os.path.join(path,'*.txt')))
print(file[0:9])
print(len(file))

# get subject id
subj_list = [ os.path.split(f)[1].replace('.fc.txt', '') for f in file ]
subj_list = np.array(subj_list, dtype=str)
print(subj_list[0:9])

####### use all subject, 8 is not enough for external validation #######
# mask = np.ones(len(subj_list), dtype=bool)
# mask[[80,81,82,83,84,85,86,87]] = False    # remove 8 for external validation
# subj_list = subj_list[mask]
# print(len(subj_list))




## beahve data
all_behav_data = pd.read_csv('../../data/pt_vars_fd_imp.txt', sep=' ', index_col=0) #use imputed data
all_behav_data.dtypes
print(all_behav_data.shape)

# filter by AP ID
all_behav_data = all_behav_data[all_behav_data.index.isin(subj_list)]
print(all_behav_data.shape)


# for modeling, no NA allowed
behav_data = all_behav_data.iloc[:,2:12].drop(columns='group').dropna()
print(behav_data.shape)



## CPM
# read in FC matrices
condition = 'fc' ## actually no need, should be something like REST or EMO in the file name, just for string matching
all_fc_data = read_in_matrices(subj_list, file_suffix=condition, data_dir=Path(path))
print(all_fc_data.shape) # edge number = n_nodes*(n_nodes-1)/2, 69751 in this case
# (88, 69751)
# fc_data = all_fc_data.copy()

## viz edges
coords = pd.read_csv("../data/coords/MMP_MNI_374_UCCHILD_coords.txt", index_col=None, header=None, sep=" ")
print(coords.shape)

# paras for CPM
cpm_kwargs = {'p_thresh': 0.01, 'corr_type': 'spearman', 'verbose': False} ## use spearman if the distribution is skewed

### k as number of subjects, LOOCV
k = all_fc_data.shape[0]
corr_type = 'partial_' + cpm_kwargs['corr_type']
# corr_type = 'partial_spearman'
selby = '_pval_' + str(cpm_kwargs['p_thresh'])
covar = behav_data[['fd.m']]

for behav in behav_data.columns[:1]:
    print(behav)

    behav_obs_pred, all_masks, corr = cpm_wrapper_seed_part_pval(all_fc_data, behav_data, behav=behav, covar=covar, k=k, seed=202210, **cpm_kwargs)
    ## count selected edges
    print('{:.2f} pos edges passed the threshold at all folds: '.format(((all_masks['pos'].sum(axis=0)/k) >= 1).sum()))
    print('{:.2f} neg edges passed the threshold at all folds: '.format(((all_masks['neg'].sum(axis=0)/k) >= 1).sum()))
    
    ## export pred table
    fn = behav + '_8yo_' + corr_type + '_fold_' + str(k) + selby
    behav_obs_pred.to_csv(fn + '_pred.csv')

    ## plot pred vs obs
    # g = plot_predictions(behav_obs_pred)
    # g.figure.tight_layout()
    # plt.savefig(os.path.join('models', fn + '_model.png'))
    # plt.close()

    ## save all_mask
    np.save(fn + '_all_masks.npy',  all_masks) 

    ## edges
    for tail in all_masks.keys():
        edfn = tail + '_' + fn
        
        edge_frac = (all_masks[tail].sum(axis=0))/(all_masks[tail].shape[0])
        edge_frac_square = sp.spatial.distance.squareform(edge_frac)
        
        # export pos and neg matrices
        # edmat = edfn + '_pearson_bn.csv'
        # np.savetxt(edmat, edge_frac_square, delimiter=',', fmt='%1.0f')
    
        # plot edges
        # all_masks(pos/neg matrices) is a binary array, k fold * n total edges, in this case is 10,69751
        if tail == 'pos':
            color = 'red'
        else:
            color = 'blue'
        
        g = plot_consistent_edges(all_masks, tail, thresh = 0.8, color = color, node_coords = coords)
        plt.savefig(os.path.join('edges', edfn + '_edges.png'))
        plt.close()