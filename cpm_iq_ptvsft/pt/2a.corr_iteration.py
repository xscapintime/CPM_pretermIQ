"""
1,000 times permutation of k-fold cross-validation
pearson r value threshold: 0.25
"""
import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
from scipy.stats import kstest
import statistics
from matplotlib import pyplot as plt
from cpmFunctions import *
import time


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



##### pearson, select by r-val ######
# paras for CPM
cpm_kwargs = {'r_thresh': 0.25, 'corr_type': 'pearson', 'verbose': False} ## use spearman if the distribution is skewed

### k as number of subjects
k = 5
corr_type = cpm_kwargs['corr_type']
selby = '_rval_' + str(cpm_kwargs['r_thresh'])

# iterations for null distribution
start_time = time.time()

rvals = {'glm':[], 'pos':[], 'neg':[]}

iters=1000
for behav in behav_data.columns[:1]:
    print(behav)

    for n in range(iters):
        print('Iteration: ' + '{:.0f}'.format(n))

        # shuffle behavioral index
        # to make one participant have other's behavioral score    
        np.random.seed(n)
        behav_shuff = behav_data
        behav_shuff.index = np.random.permutation(behav_data.index)
        
        # reorder fc data
        all_fc_data = all_fc_data.loc[behav_shuff.index,]

        behav_obs_pred, all_masks, corr = cpm_wrapper_seed(all_fc_data, behav_shuff, behav=behav, k=k, seed=202211, **cpm_kwargs)
        
        for tail in ['glm', 'pos', 'neg']:
            x = np.squeeze(behav_obs_pred.filter(regex=("obs")).astype(float))
            y = np.squeeze(behav_obs_pred.filter(regex=(tail)).astype(float))
        
            rhos = sp.stats.spearmanr(x,y)[0]
            rvals[tail].append(rhos)

    ## export r-vals, sort p-vals later # list to txt, line by line
    fn = behav + '_8yo_' + corr_type + '_fold_' + str(k) + selby + '_nullrvals_iter{:.0f}'.format(iters)
    # with open(fn, 'w') as file:
    #     for v in rval:
    #         file.write('%s\n' % v)
    pd.DataFrame(rvals).to_csv(fn + '.csv', index=False)
    

print("--- %s seconds for 1 behav, select by pearson rval---" %(time.time() - start_time))
