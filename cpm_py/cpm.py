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



## CPM
# read in FC matrices
condition = 'fc' ## actually no need
all_fc_data = read_in_matrices(subj_list, file_suffix=condition, data_dir=Path(path))

# heatmap
plt.figure(figsize = (6,6))
for s in range(0,all_fc_data.shape[0]):
    g = sns.heatmap(sp.spatial.distance.squareform(all_fc_data.iloc[s,:]), square=True)
    g.figure.tight_layout()
    plt.suptitle(subj_list[s])
    # plt.show()
    plt.savefig(os.path.join('heatmaps', subj_list[s] + '_fc.png'))
    plt.close()



# paras for CPM
cpm_kwargs = {'r_thresh': 0.2, 'corr_type': 'pearson', 'verbose': False}


# for modeling, no NA allowed
behav_data = all_behav_data.dropna()
fc_data = all_fc_data.loc[behav_data.index]


for behav in behav_data.columns[6:]:
    # behav = 'ListSort_Unadj'
    print(behav)
    behav_obs_pred, all_masks = cpm_wrapper(fc_data, behav_data, behav=behav, **cpm_kwargs)
    g = plot_predictions(behav_obs_pred) 
    g.figure.tight_layout()
    plt.savefig(os.path.join('models', behav + '_model.png'))
    plt.close()

