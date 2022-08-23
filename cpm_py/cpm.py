import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from cpmFinctions import *


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
plt.savefig(os.path.join('dist', 'var19_pairdist.pdf'))
plt.close()

plt.figure(figsize = (5,5))
p = sns.pairplot(all_pca_data.iloc[:,:2], kind ='scatter')
plt.savefig(os.path.join('dist', 'pcs2_pairdist.pdf'))
plt.close()


# paras for CPM
cpm_kwargs = {'r_thresh': 0.3, 'corr_type': 'pearson', 'verbose': False} ## use spearman if the distribution is skewed

# pred vs obs scatter
# all vars, don't need this, use PCs
# for behav in behav_data.columns[6:]:
#     print(behav)
#     behav_obs_pred, all_masks = cpm_wrapper(fc_data, behav_data, behav=behav, **cpm_kwargs)
#     g = plot_predictions(behav_obs_pred) 
#     g.figure.tight_layout()
#     plt.savefig(os.path.join('models', behav + '_model.png'))
#     plt.close()

for behav in all_pca_data.columns[:2]:
    print(behav)
    behav_obs_pred, all_masks, mask_dict, corr = cpm_wrapper(fc_data, all_pca_data, behav=behav, **cpm_kwargs)
    g = plot_predictions(behav_obs_pred) 
    g.figure.tight_layout()
    plt.savefig(os.path.join('models', behav + '_model.png'))
    plt.close()

#### non-correlated and non-sig, few subjects = 51 ####


## viz edges
coords = pd.read_csv("../data/coords/MMP_MNI_374_UCCHILD_coords.txt", index_col=None, header=None, sep=" ")
print(coords.shape)

plot_consistent_edges(all_masks, "pos", thresh = 0.8, color = 'red', node_coords = coords)
plot_consistent_edges(all_masks, "neg", thresh = 0.8, color = 'blue', node_coords = coords)
