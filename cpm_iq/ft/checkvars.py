import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
from scipy.stats import kstest
from scipy.stats import shapiro
# import statistics
from matplotlib import pyplot as plt
# from cpmFunctions import *


## check FC file names
path = '../../data/fc_individual_ft/'
file = glob.glob(os.path.join(path,'*.txt'))
print(file[0:9])
print(len(file))

# get subject id
subj_list = [ os.path.split(f)[1].replace('.fc.txt', '') for f in file ]
subj_list = np.array(subj_list, dtype=str)
print(subj_list[0:9])

# # remove 5 for external validation
# mask = np.ones(len(subj_list), dtype=bool)
# mask[[40,41,42,43,44]] = False
# subj_list = subj_list[mask]
# print(len(subj_list))


## beahve data
# all_behav_data = pd.read_csv('../data/id_vars_fin.csv', dtype={'Eprime_ID': str})
all_behav_data = pd.read_csv('../../data/ft_vars_fd.txt', sep=' ', index_col=0)
all_behav_data.dtypes
print(all_behav_data.shape)

# filter by AP ID
all_behav_data = all_behav_data[all_behav_data.index.isin(subj_list)]
print(all_behav_data.shape)


# for modeling, no NA allowed
behav_data = all_behav_data.iloc[:,2:12].drop(columns='group').dropna()
print(behav_data.shape)

##### something wrong with AP095 #####
print(behav_data.loc['AP095','WISC_PS_CS'])
# 1001
behav_data.loc['AP095','WISC_PS_CS'] = 101


## check distributions
# compare all the vars and FD
plt.figure(figsize = (8,8))
p = sns.pairplot(behav_data, kind ='scatter')
plt.savefig(os.path.join('dist', 'ft_8yo_wisc_srs_fd_pairdist.pdf'))
plt.close()

# wisc only
plt.figure(figsize = (5,5))
p = sns.pairplot(behav_data.iloc[:,:5], kind ='scatter')
plt.savefig(os.path.join('dist', 'ft_8yo_wisc_pairdist.pdf'))
plt.close()

# wisc full vs fd mean
plt.figure(figsize = (3,3))
p = sns.pairplot(behav_data.loc[:,['WISC_FULL_CS', 'fd.m']], kind ='scatter')
plt.savefig(os.path.join('dist', 'ft_8yo_wiscfull_fd_pairdist.pdf'))
plt.close()




## Check correlation between behavioral stat vs head motion
# pearson
for var in behav_data.columns[:7]:
    
    r = sp.stats.pearsonr(behav_data[var], behav_data['fd.m'])[0]
    p = sp.stats.pearsonr(behav_data[var], behav_data['fd.m'])[1]

    plt.figure(figsize = (5,5))
    g = sns.regplot(x=behav_data[var], y=behav_data['fd.m'], color='blue')
    g.figure.tight_layout()
    g.annotate('r = {0:.2f}'.format(r), xy = (0.7, 0.1), xycoords = 'axes fraction')
    g.annotate('p = {0:.2f}'.format(p), xy = (0.7, 0.05), xycoords = 'axes fraction') ##
    plt.savefig(os.path.join('dist', var + '_vs_fd.pdf'))
    plt.close()

# spearman
for var in behav_data.columns[:7]:
    
    r = sp.stats.spearmanr(behav_data[var], behav_data['fd.m'])[0]
    p = sp.stats.spearmanr(behav_data[var], behav_data['fd.m'])[1]

    plt.figure(figsize = (5,5))
    g = sns.regplot(x=behav_data[var], y=behav_data['fd.m'], color='blue')
    g.figure.tight_layout()
    g.annotate('r = {0:.2f}'.format(r), xy = (0.7, 0.1), xycoords = 'axes fraction')
    g.annotate('p = {0:.2f}'.format(p), xy = (0.7, 0.05), xycoords = 'axes fraction') ##
    plt.savefig(os.path.join('dist', var + '_vs_fd_spear.pdf'))
    plt.close()

##### not correlated with motion much #####



## check normality of vars
# Kolmogorov-Smirnov test  
for var in behav_data.columns[:7]:
    ktest = kstest(behav_data[var], 'norm')
    if ktest[1] > 0.05:
        print(var + ': normally distributed.')
    else:
        print(var + ': NOT normally distributed. pval=' + str(ktest[1]))


# Shapiro-Wilk test
for var in behav_data.columns[:7]:
    sptest = shapiro(behav_data[var])
    if ktest[1] > 0.05:
        print(var + ': normally distributed.')
    else:
        print(var + ': NOT normally distributed. pval=' + str(sptest[1]))

