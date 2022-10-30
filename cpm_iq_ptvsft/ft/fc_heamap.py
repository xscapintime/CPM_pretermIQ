import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
from scipy.stats import kstest
from scipy.stats import shapiro
# import statistics
from matplotlib import pyplot as plt
from cpmFunctions import *


## check FC file names
path = '../../data/fc_individual_ft/'
file = glob.glob(os.path.join(path,'*.txt'))
print(file[0:9])
print(len(file))

# get subject id
subj_list = [ os.path.split(f)[1].replace('.fc.txt', '') for f in file ]
subj_list = np.array(subj_list, dtype=str)
print(subj_list[0:9])



## CPM
# read in FC matrices
condition = 'fc' ## actually no need, should be something like REST or EMO to distinguish different matrices
all_fc_data = read_in_matrices(subj_list, file_suffix=condition, data_dir=Path(path))
print(all_fc_data.shape) # edge number = n_nodes*(n_nodes-1)/2, 69751 in this case
# (45, 69751)


# heatmaps, check FC of every subject, not necessary
plt.figure(figsize = (6,5))
for s in range(0,all_fc_data.shape[0]):
    g = sns.heatmap(sp.spatial.distance.squareform(all_fc_data.iloc[s,:]), square=True, xticklabels=20, yticklabels=20)
    # g.figure.tight_layout()
    plt.title(subj_list[s])
    # plt.show()
    plt.savefig(os.path.join('heatmaps', subj_list[s] + '_fc.png'))
    plt.close()


# check ouliers
np.histogram(all_fc_data)
