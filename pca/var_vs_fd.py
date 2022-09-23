import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
from scipy.stats import kstest
import statistics
from matplotlib import pyplot as plt


# load impted vars table
vars = pd.read_csv('../data/vars_75subj_imputed.csv', index_col=0)

# laod fd table
fd = pd.read_csv('../data/fd_mean_max.txt', index_col=0, sep=' ')
fd = fd[['fd.m', 'fd.max']]


# check correlation between fd mean and all vars

for var in vars.columns :
    
    r = sp.stats.pearsonr(vars[var], fd['fd.m'])[0]
    p = sp.stats.pearsonr(vars[var], fd['fd.m'])[1]

    plt.figure(figsize = (5,5))
    g = sns.regplot(x=vars[var], y=fd['fd.m'], color='blue')
    g.figure.tight_layout()
    g.annotate('r = {0:.2f}'.format(r), xy = (0.7, 0.1), xycoords = 'axes fraction')
    g.annotate('p = {0:.2f}'.format(p), xy = (0.7, 0.05), xycoords = 'axes fraction') ##
    plt.savefig(os.path.join('.', var + '_vs_fd.png'))
    plt.close()
