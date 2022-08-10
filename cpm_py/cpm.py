import os, glob
import numpy as np
import pandas as pd
import cpmFinctions

## check FC file names
path = '../data/fc_individual/'
file = glob.glob(os.path.join(path,'*.txt'))
file

# get subject id
subj_list = [ os.path.split(f)[1].replace('.fc.txt', '') for f in file ]
subj_list = np.array(subj_list, dtype=str)


## beahve data
all_behav_data = pd.read_csv('../data/unrestricted_behav_data_n337.csv', dtype={'Subject': str})


condition = 'fc' ## actually no need
cpm_kwargs = {'r_thresh': 0.2, 'corr_type': 'pearson', 'verbose': False}


behav = 'ListSort_Unadj'



all_fc_data = read_in_matrices(subj_list, file_suffix=condition)
behav_obs_pred, all_masks = cpm_wrapper(all_fc_data, all_behav_data, behav=behav, **cpm_kwargs)
g = plot_predictions(behav_obs_pred)
g.set_title(condition)
plt.show()
