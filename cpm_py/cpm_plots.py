import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
from scipy.stats import kstest
import statistics
from matplotlib import pyplot as plt
from cpmFinctions import *


## check FC file names
path = '../data/fc_individual/'
file = glob.glob(os.path.join(path,'*.txt'))
# file

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
all_pca_data = pd.read_csv('../data/var17_imp_pca.csv', index_col=0)
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
plt.figure(figsize = (12,12))
p = sns.pairplot(behav_data, kind ='scatter')
plt.savefig(os.path.join('dist', 'var17_pairdist.pdf'))
plt.close()

plt.figure(figsize = (5,5))
p = sns.pairplot(pca_data.iloc[:,:2], kind ='scatter')
plt.savefig(os.path.join('dist', 'pcs2_pairdist.pdf'))
plt.close()

plt.figure(figsize = (12,12))
p = sns.pairplot(pca_data, kind ='scatter')
plt.savefig(os.path.join('dist', 'pca_pairdist.pdf'))
plt.close()



## Check correlation between behavioral stat vs head motion
fd = pd.read_csv('../data/fd_mean_max.txt', sep=' ', header=0)
fd = fd[fd.index.isin(pca_data.index)]

sp.stats.pearsonr(fd['fd.m'], pca_data['PC1'])
# (0.1783240798203215, 0.1258458793145855)
sp.stats.pearsonr(fd['fd.m'], pca_data['PC2'])
# (-0.2585155832221153, 0.02512842358871094)

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
# KstestResult(statistic=0.22012991583245278, pvalue=0.0011335296969197106)
# pval < 0.05, rejct H0, that PC1 is not normally distributed

kstest(pca_data['PC2'], 'norm')
# KstestResult(statistic=0.20643213219397566, pvalue=0.002792582138901656)


## load covariates data
cova_data = pd.read_csv('../data/pt_sbj75withfc_newvars_fd_pca_merged.csv', index_col=1)
cova_data = cova_data[['sex', 'age22', 'age4', 'age8', 'Neonatal Sickness', 'IMD Score']].dropna() ##
cova_data.shape


##### for NA in covariates ####
### might remove this step ###
fc_data = all_fc_data.loc[cova_data.index]
print(fc_data.shape)

pca_data = pca_data.loc[cova_data.index]
print(pca_data.shape)



## viz edges
coords = pd.read_csv("../data/coords/MMP_MNI_374_UCCHILD_coords.txt", index_col=None, header=None, sep=" ")
print(coords.shape)

# paras for CPM
cpm_kwargs = {'r_thresh': 0.38, 'corr_type': 'spearman', 'verbose': False} ## use spearman if the distribution is skewed
# perc_thresh=1 for top edges selection
# cpm_kwargs = {'p_thresh': 0.01, 'corr_type': 'pearson', 'verbose': False} ## use spearman if the distribution is skewed
covar = ['sex', 'age22', 'age4', 'age8', 'Neonatal Sickness', 'IMD Score']


# plots, would be one of the random cases
# set k as the number of subjects, Leave-One-Out Cross-Validation (LOOCV)

## make a test fc_data of 50 nodes, 68*67/2
fc_data_t = fc_data.loc[:,:2277]

### k as number of subjects, LOOCV
for behav in pca_data.columns[:2]:
    print(behav)

    behav_obs_pred, all_masks, corr = cpm_wrapper_seed_part(fc_data_t, pca_data, behav=behav, cova_data=cova_data, covar=covar, k=fc_data_t.shape[0], seed=202209, **cpm_kwargs)
    ## count selected edges
    print('{:.2f} pos edges passed all folds: '.format(((all_masks['pos'].sum(axis=0)/10) >= 1).sum()))
    print('{:.2f} neg edges passed all folds: '.format(((all_masks['neg'].sum(axis=0)/10) >= 1).sum()))
    
    ## plot pred vs obs
    g = plot_predictions(behav_obs_pred)
    g.figure.tight_layout()
    plt.savefig(os.path.join('models', behav + k + '_partcorr_model.png'))
    plt.close()
    
    ## plot edges
    # all_masks(pos/neg matrices) is a binary array, k fold * n total edges, in this case is 10,69751
    # g = plot_consistent_edges(all_masks, "pos", thresh = 1, color = 'red', node_coords = coords)
    # plt.savefig(os.path.join('edges', behav + k +'_partcorr_pos_edges.png'))
    # plt.close()
    # g = plot_consistent_edges(all_masks, "neg", thresh = 1, color = 'blue', node_coords = coords)
    # plt.savefig(os.path.join('edges', behav + k + '_partcorr_neg_edges.png'))
    # plt.close()


#### non-correlated and non-sig, few subjects = 51 ####
