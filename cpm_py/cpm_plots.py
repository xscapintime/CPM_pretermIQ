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
all_behav_data = pd.read_csv('../data/vars_8yo_82subj_imputed.csv', index_col=0)
all_behav_data.dtypes

# filter by AP ID
all_behav_data = all_behav_data[all_behav_data.index.isin(subj_list)]
# all_behav_data = all_behav_data.set_index('AP_ID')


# for modeling, no NA allowed
# only keep 8yo
# behav_data = all_behav_data.iloc[:,6:11]#.dropna()
behav_data = all_behav_data
print(behav_data.shape)

## PCs
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


# heatmap, check FC of every subject, not necessary
plt.figure(figsize = (6,5))
for s in range(0,all_fc_data.shape[0]):
    g = sns.heatmap(sp.spatial.distance.squareform(all_fc_data.iloc[s,:]), square=True, xticklabels=20, yticklabels=20)
    # g.figure.tight_layout()
    plt.title(subj_list[s])
    # plt.show()
    plt.savefig(os.path.join('heatmaps', subj_list[s] + '_fc.png'))
    plt.close()


# compare all the vars and PCs
# plt.figure(figsize = (12,12))
# p = sns.pairplot(behav_data, kind ='scatter')
# plt.savefig(os.path.join('dist', 'var17_pairdist.pdf'))
# plt.close()

# plt.figure(figsize = (5,5))
# p = sns.pairplot(pca_data.iloc[:,:2], kind ='scatter')
# plt.savefig(os.path.join('dist', 'pcs2_pairdist.pdf'))
# plt.close()

# plt.figure(figsize = (12,12))
# p = sns.pairplot(pca_data, kind ='scatter')
# plt.savefig(os.path.join('dist', 'pca_pairdist.pdf'))
# plt.close()

plt.figure(figsize = (6,6))
p = sns.pairplot(behav_data, kind ='scatter')
plt.savefig(os.path.join('dist', 'var6_8yo_pairdist.pdf'))
plt.close()

# plt.figure(figsize = (5,5))
# p = sns.pairplot(pca_data.iloc[:,:2], kind ='scatter')
# plt.savefig(os.path.join('dist', 'pcs2_pairdist.pdf'))
# plt.close()

plt.figure(figsize = (12,12))
p = sns.pairplot(pca_data, kind ='scatter')
plt.savefig(os.path.join('dist', 'pca_pairdist.pdf'))
plt.close()



## Check correlation between behavioral stat vs head motion
fd = pd.read_csv('../data/fd_mean_max.txt', sep=' ', header=0)
fd = fd[fd.index.isin(pca_data.index)]

sp.stats.pearsonr(fd['fd.m'], pca_data['PC1'])
# (0.3105841565689169, 0.006688262332681814)
sp.stats.pearsonr(fd['fd.m'], pca_data['PC2'])
# (0.05783213135402286, 0.622124464908075)

for pc in ['PC1', 'PC2'] :
    
    r = sp.stats.pearsonr(pca_data[pc], fd['fd.m'])[0]
    p = sp.stats.pearsonr(pca_data[pc], fd['fd.m'])[1]

    plt.figure(figsize = (5,5))
    g = sns.regplot(x=pca_data[pc], y=fd['fd.m'], color='blue')
    g.figure.tight_layout()
    g.annotate('r = {0:.2f}'.format(r), xy = (0.7, 0.1), xycoords = 'axes fraction')
    g.annotate('p = {0:.2f}'.format(p), xy = (0.7, 0.05), xycoords = 'axes fraction') ##
    plt.savefig(os.path.join('dist', pc + '_vs_fd.png'))
    plt.close()



## check normality of PCs
kstest(pca_data['PC1'], 'norm')
# KstestResult(statistic=0.1660644636021128, pvalue=0.02820521696212752)
# pval < 0.05, rejct H0, that PC1 is not normally distributed ###

#### pval > 0.05, reject H1, PC1 normal

kstest(pca_data['PC2'], 'norm')
# KstestResult(statistic=0.1054905957612926, pvalue=0.3496145154997725)


# ## load covariates data
# cova_data = pd.read_csv('../data/pt_sbj75withfc_newvars_imp_fd_pca_merged.csv', index_col=1)
# cova_data = cova_data[['sex', 'age8', 'Neonatal Sickness', 'IMD Score']]#.dropna() ##
# cova_data.shape


# ##### for NA in covariates ####
# ### might remove this step ###
# fc_data = all_fc_data.loc[cova_data.index]
# print(fc_data.shape)

# pca_data = pca_data.loc[cova_data.index]
# print(pca_data.shape)



# ## viz edges
# coords = pd.read_csv("../data/coords/MMP_MNI_374_UCCHILD_coords.txt", index_col=None, header=None, sep=" ")
# print(coords.shape)

# # paras for CPM
# cpm_kwargs = {'r_thresh': 0.35, 'corr_type': 'spearman', 'verbose': False} ## use spearman if the distribution is skewed
# # perc_thresh=1 for top edges selection
# # cpm_kwargs = {'p_thresh': 0.01, 'corr_type': 'pearson', 'verbose': False} ## use spearman if the distribution is skewed
# covar = ['sex', 'age8', 'Neonatal Sickness', 'IMD Score']


# # plots, would be one of the random cases
# # set k as the number of subjects, Leave-One-Out Cross-Validation (LOOCV)

# ## make a test fc_data of 50 nodes, 68*67/2
# fc_data_t = fc_data.loc[:,:2277]

# ### k as number of subjects, LOOCV
# # k = fc_data_t.shape[0]
# k = 10
# corr_type = 'spearman'

# for behav in pca_data.columns[0:1]: #PC1
#     print(behav)

#     behav_obs_pred, all_masks, corr = cpm_wrapper_seed_part(fc_data, pca_data, behav=behav, cova_data=cova_data, covar=covar, k=k, seed=202209, **cpm_kwargs)
#     ## count selected edges
#     print('{:.2f} pos edges passed all folds: '.format(((all_masks['pos'].sum(axis=0)/10) >= 1).sum()))
#     print('{:.2f} neg edges passed all folds: '.format(((all_masks['neg'].sum(axis=0)/10) >= 1).sum()))
    
#     ## export pred table
#     fn = behav + '_' + corr_type + '_fold_' + str(k) + '_partcorr_pred.csv'
#     behav_obs_pred.to_csv(fn)

#     ## plot pred vs obs
#     g = plot_predictions(behav_obs_pred)
#     g.figure.tight_layout()
#     plt.savefig(os.path.join('models', behav + '_' + corr_type + '_fold_' + str(k) + '_partcorr_model.png'))
#     plt.close()

    
#     ## edges
#     for tail in all_masks.keys():
#         edfn = tail + '_' + behav + '_' + corr_type + '_fold_' + str(k)
        
#         edge_frac = (all_masks[tail].sum(axis=0))/(all_masks[tail].shape[0])
#         edge_frac_square = sp.spatial.distance.squareform(edge_frac)
        
#         # export pos and neg matrices
#         edmat = edfn + '_partcorr_bn.csv'
#         np.savetxt(edmat, edge_frac_square, delimiter=',', fmt='%1.0f')
    
#         # plot edges
#         # all_masks(pos/neg matrices) is a binary array, k fold * n total edges, in this case is 10,69751
#         if tail == 'pos':
#             color = 'red'
#         else:
#             color = 'blue'
        
#         g = plot_consistent_edges(all_masks, tail, thresh = 1, color = color, node_coords = coords)
#         plt.savefig(os.path.join('edges', edfn + '_partcorr_edges.png'))
#         plt.close()
    
