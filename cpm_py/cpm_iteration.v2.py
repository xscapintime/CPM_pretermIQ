import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
from scipy.stats import kstest
import statistics
import time
from matplotlib import pyplot as plt
from cpmFinctions import *
from itertools import product



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
all_pca_data = pd.read_csv('../data/var19_pca.csv', index_col=0)
all_pca_data.dtypes

# filter
pca_data = all_pca_data[all_pca_data.index.isin(subj_list)]
print(pca_data.shape)


## Check correlation between behavioral stat vs head motion
fd = pd.read_csv('../data/fd_mean_max.txt', sep=' ', header=0)
fd = fd[fd.index.isin(pca_data.index)]

sp.stats.pearsonr(fd['fd.m'], pca_data['PC1'])
# (-0.2773401388237849, 0.048798998416344996)

sp.stats.pearsonr(fd['fd.m'], pca_data['PC2'])
# (-0.214343195323698, 0.13094555084819162)


## check normality of PCs
kstest(pca_data['PC1'], 'norm')
# KstestResult(statistic=0.31659944232769943, pvalue=4.812505160566234e-05)

kstest(pca_data['PC2'], 'norm')
# KstestResult(statistic=0.18070356799566734, pvalue=0.06273635080068618)


##### ============================= #####
## So for PC1 


## CPM
# read in FC matrices
condition = 'fc' ## actually no need, should be like REST or EMO
fc_data = read_in_matrices(pca_data.index, file_suffix=condition, data_dir=Path(path))
# fc_data = all_fc_data.loc[behav_data.index]


## pearson
# paras for CPM
# cpm_kwargs = {'r_thresh': 0.38, 'corr_type': 'pearson', 'verbose': False} ## use spearman if the distribution is skewed
# cpm_kwargs = {'p_thresh': 0.01, 'corr_type': 'pearson', 'verbose': False} ## use spearman if the distribution is skewed

cpm_kwargs = {'corr_type': 'pearson', 'verbose': False} ## use spearman if the distribution is skewed


k = [2, 5, 10]
rthresh = [0.1, 0.2, 0.38]

keys = list(product(k, rthresh))

for behav in pca_data.columns[:2]:
    print(behav)
    
    for key in keys: 
        
        start_time = time.time()

        # k only affect mask as it's a k * edges array
        # behav_obs_pred and corr are actually the output of last fold
        behav_obs_pred, all_masks, corr = cpm_wrapper_seed(fc_data, pca_data, behav=behav, k=key[0], seed = 202208, r_thresh = key[1], **cpm_kwargs)
    
        print("--- %s seconds for 1 PC, select by rval---" %(time.time() - start_time))

        ## count selected edges
        print('{:.2f} pos edges passed all folds at {}-fold & {} r threshold: '.format(((all_masks['pos'].sum(axis=0)/10) >= 1).sum(), key[0], key[1]))
        print('{:.2f} neg edges passed all folds at {}-fold & {} r threshold: '.format(((all_masks['neg'].sum(axis=0)/10) >= 1).sum(), key[0], key[1]))
        

        ## plot pred vs obs
        g = plot_predictions(behav_obs_pred)
        g.figure.tight_layout()
        plt.savefig(os.path.join('models', behav + '_{}fold_{}rthresh'.format(key[0], key[1]) + '_model.v2.png'))
        plt.close()




## spearman
# paras for CPM
cpm_kwargs = {'r_thresh': 0.38, 'corr_type': 'spearman', 'verbose': False} ## use spearman if the distribution is skewed
# np.random.seed(202208)
# cpm_kwargs = {'p_thresh': 0.01, 'corr_type': 'pearson', 'verbose': False} ## use spearman if the distribution is skewed

start_time = time.time()

for behav in pca_data.columns[:2]:
    print(behav)
    
    behav_obs_pred, all_masks, corr = cpm_wrapper_seed(fc_data, pca_data, behav=behav, k=10, seed = 202208, **cpm_kwargs)
    
    ## count selected edges
    print('{:.2f} pos edges passed all folds: '.format(((all_masks['pos'].sum(axis=0)/10) >= 1).sum()))
    print('{:.2f} neg edges passed all folds: '.format(((all_masks['neg'].sum(axis=0)/10) >= 1).sum()))

    ## plot pred vs obs
    g = plot_predictions(behav_obs_pred)
    g.figure.tight_layout()
    plt.savefig(os.path.join('models', behav + '_model.spearman.v2.png'))
    plt.close()


print("--- %s seconds for 2 PCs, select by Spearman rval---" %(time.time() - start_time))
# --- 2589.731003522873 seconds for 2 PCs, select by Spearman rval---



## Shuffle the results
# iterations for null distribution

iters=1000
for n in range(iters):
    print('Iteration: ' + '{:.0f}'.format(n))
    rval = []

    rval.append(sp.stats.pearsonr(x,y)[0])
    print('Mode of rvals: ' + '{:.2f}'.format(statistics.mode(rval)))

    # histogram
    plt.figure(figsize = (5,6))
    g = sns.displot(rval, kde=True)
    plt.subplots_adjust(top=0.85)
    plt.title('r-vals of Pred vs Obs of edges correlated to ' + behav + '\n  Iter={:.0f}, mode={:.2f}'.format(iters, statistics.mode(rval)))
    plt.savefig(os.path.join('dist', behav + '_{:.0f}_rvals_dist.png'.format(iters)))
    plt.close()


    ## export r-vals, sort p-vals later # list to txt, line by line
    fn = behav + '_predvsobs_rvals_iter{:.0f}.txt'.format(iters)
    with open(fn, 'w') as file:
        for v in rval:
            file.write('%s\n' % v)





## viz edges
coords = pd.read_csv("../data/coords/MMP_MNI_374_UCCHILD_coords.txt", index_col=None, header=None, sep=" ")
print(coords.shape)




## test spearman speed

# scipy
start_time = time.time()

corr = []
for edge in fc_data.columns:
    r_val = sp.stats.spearmanr(fc_data.loc[:,edge], pca_data['PC1'])[0].item()
    corr.append(r_val)
    # corr = np.array(corr) #

print("--- %s seconds for 1 PC, Spearman from Scipy ---" %(time.time() - start_time))
# --- 131.64980578422546 seconds for 1 PC, Spearman from Scipy ---


# pandas
# len_df1 = df1.shape[0]
df2_index = pca_data.index.values.tolist()

df = pd.concat([pca_data[['PC1']], fc_data], axis=1).reset_index(drop=True)#.T


start_time = time.time()

corr_pd = [ df[['PC1',i]].corr('spearman')['PC1'][i] for i in range(fc_data.shape[1]) ]
# 3m 38.4a

print("--- %s seconds for 1 PC, Spearman from Pandas ---" %(time.time() - start_time))



# by hand
df = pd.concat([pca_data[['PC1']], fc_data], axis=1).reset_index(drop=True)#.T

