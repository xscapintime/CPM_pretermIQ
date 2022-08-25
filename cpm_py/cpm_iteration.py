import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
import statistics
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



# paras for CPM
# cpm_kwargs = {'r_thresh': 0.3, 'corr_type': 'pearson', 'verbose': False} ## use spearman if the distribution is skewed
# np.random.seed(202208)
cpm_kwargs = {'p_thresh': 0.01, 'corr_type': 'pearson', 'verbose': False} ## use spearman if the distribution is skewed

## viz edges
coords = pd.read_csv("../data/coords/MMP_MNI_374_UCCHILD_coords.txt", index_col=None, header=None, sep=" ")
print(coords.shape)


# pred vs obs scatter
# all vars, don't need this, use PCs
# for behav in behav_data.columns[6:]:
#     print(behav)
#     behav_obs_pred, all_masks = cpm_wrapper(fc_data, behav_data, behav=behav, **cpm_kwargs)
#     g = plot_predictions(behav_obs_pred) 
#     g.figure.tight_layout()
#     plt.savefig(os.path.join('models', behav + '_model.png'))
#     plt.close()


# iterations for null distribution
iters=10
for behav in all_pca_data.columns[:2]:
    print(behav)

    rval = []
    for n in range(iters):
        print('Iteration: ' + '{:.0f}'.format(n))
        behav_obs_pred, all_masks, mask_dict, corr, pval = cpm_wrapper(fc_data, all_pca_data, behav=behav, **cpm_kwargs)

        x = np.squeeze(behav_obs_pred.filter(regex=("obs")).astype(float))
        y = np.squeeze(behav_obs_pred.filter(regex=("glm")).astype(float))
        rval.append(sp.stats.pearsonr(x,y)[0])

        ## count selected edges
        print('Selected pos edges: %' + '{:.2f}'.format(sum(mask_dict['pos'])/len(mask_dict['pos']) * 100))
        print('Selected neg edges: %' + '{:.2f}'.format(sum(mask_dict['neg'])/len(mask_dict['neg']) * 100))


    print('Mode of rvals: ' + '{:.2f}'.format(statistics.mode(rval)))

    ## histogram
    # g = sns.displot(rval, kde=True)
    # plt.title('r-vals of Pred vs Obs of edges correlated to ' + behav + '\n  Iter=' + '{:.0f}'.format(iters) + ', mode=' + '{:.2f}'.format(statistics.mode(rval)))
    # plt.savefig(os.path.join('dist', behav + '_{:.0f}'.format(iters) + '_rvals_dist.png'))
    # plt.close()

    ## export r-vals, sort p-vals later # list to txt, line by line
    fn = behav + '_predvsobs_rvals' + '_iter{:.0f}'.format(iters) + '.txt'
    with open(fn, 'w') as file:
        for v in rval:
            file.write('%s\n' % v)
    

