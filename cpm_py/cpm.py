import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
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
all_fc_data = read_in_matrices(subj_list, file_suffix=condition, data_dir=Path(path))

# heatmap
s = 0

sns.heatmap(sp.spatial.distance.squareform(all_fc_data.iloc[s,:]))
plt.show()




# paras for CPM
condition = 'fc' ## actually no need
cpm_kwargs = {'r_thresh': 0.2, 'corr_type': 'pearson', 'verbose': False}



for behav in all_behav_data.columns:
    # behav = 'ListSort_Unadj'
    print(behav)



behav_obs_pred, all_masks = cpm_wrapper(all_fc_data, all_behav_data, behav=behav, **cpm_kwargs)
g = plot_predictions(behav_obs_pred)
g.set_title(condition)
plt.show()
