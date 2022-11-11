"""
Combine r null distribution, manually parallel processed
"""

import os, glob
import csv
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt



path = 'permutations/'
files = sorted(glob.glob(os.path.join(path,'*_nullrvals_iter*.csv')))


pearson_perms = []
spearman_perms = []
parsp_perms = []

for f in files:
    if '_8yo_pearson_' in f:
        perms = pd.read_csv(f)
        pearson_perms.append(perms)
        pearson_df = pd.concat(pearson_perms, axis=0)
    
    elif '_8yo_spearman_' in f:
        perms = pd.read_csv(f)
        spearman_perms.append(perms)
        spearman_df = pd.concat(spearman_perms, axis=0)

    elif '_8yo_partial_spearman_' in f:
        perms = pd.read_csv(f)
        parsp_perms.append(perms)
        parsp_df = pd.concat(parsp_perms, axis=0)


# export
pr_fn = [ os.path.split(f)[1].split('_iter')[0] for f in files if '_8yo_pearson_' in f ][0]
sp_fn = [ os.path.split(f)[1].split('_iter')[0] for f in files if '_8yo_spearman_' in f ][0]
pa_fn = [ os.path.split(f)[1].split('_iter')[0] for f in files if '_8yo_partial_spearman_' in f ][0]


pearson_df.to_csv(pr_fn +'_iter1000.csv', index=False)
spearman_df.to_csv(sp_fn +'_iter1000.csv', index=False)
parsp_df.to_csv(pa_fn +'_iter1000.csv', index=False)