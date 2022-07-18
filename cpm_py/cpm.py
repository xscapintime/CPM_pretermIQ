import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
%matplotlib inline
import pandas as pd
import seaborn as sns
from pathlib import Path
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)



condition = 'WM'
behav = 'ListSort_Unadj'
cpm_kwargs = {'r_thresh': 0.2, 'corr_type': 'pearson', 'verbose': False}


all_fc_data = read_in_matrices(subj_list, file_suffix=condition)
behav_obs_pred, all_masks = cpm_wrapper(all_fc_data, all_behav_data, behav=behav, **cpm_kwargs)
g = plot_predictions(behav_obs_pred)
g.set_title(condition)
plt.show()
