"""
Check variables distribution and their correlation with motion
PT subjects
"""
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
path = '../../data/fc_individual_pt/'
file = glob.glob(os.path.join(path,'*.txt'))
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
all_behav_data = pd.read_csv('../../data/pt_vars_fd_imp.txt', sep=' ', index_col=0) # use imputed data
all_behav_data.dtypes
print(all_behav_data.shape)

# filter by AP ID
all_behav_data = all_behav_data[all_behav_data.index.isin(subj_list)]
print(all_behav_data.shape)


# for modeling, no NA allowed
# already imputed
behav_data = all_behav_data.iloc[:,2:12].drop(columns='group').dropna()
print(behav_data.shape)


## check distributions
# compare all the vars and FD
plt.figure(figsize = (8,8))
p = sns.pairplot(behav_data, kind ='scatter')
plt.savefig(os.path.join('dist', 'pt_8yo_wisc_srs_fd_pairdist.png'))
plt.close()

# wisc only
plt.figure(figsize = (5,5))
p = sns.pairplot(behav_data.iloc[:,:5], kind ='scatter')
plt.savefig(os.path.join('dist', 'pt_8yo_wisc_pairdist.png'))
plt.close()

# wisc full vs fd mean
plt.figure(figsize = (3,3))
p = sns.pairplot(behav_data.loc[:,['WISC_FULL_CS', 'fd.m']], kind ='scatter')
plt.savefig(os.path.join('dist', 'pt_8yo_wiscfull_fd_pairdist.png'))
plt.close()


## check normality of vars
# Kolmogorov-Smirnov test  
for var in behav_data.columns[:7]:
    ktest = kstest(behav_data[var], 'norm')
    if ktest[1] > 0.05:
        print(var + ': normally distributed.' + str(ktest[1]))
    else:
        print(var + ': NOT normally distributed. pval=' + str(ktest[1]))

# WISC_FULL_CS: NOT normally distributed. pval=0.0
# WISC_VCI_CS: NOT normally distributed. pval=0.0
# WISC_PR_CS: NOT normally distributed. pval=0.0
# WISC_WM_CS: NOT normally distributed. pval=0.0
# WISC_PS_CS: NOT normally distributed. pval=0.0
# srs.rrb: NOT normally distributed. pval=0.0
# srs.sci: NOT normally distributed. pval=0.0


# Shapiro-Wilk test
for var in behav_data.columns[:7]:
    sptest = shapiro(behav_data[var])
    if sptest[1] > 0.05:
        print(var + ': normally distributed. pval=' + str(sptest[1]))
    else:
        print(var + ': NOT normally distributed. pval=' + str(sptest[1]))


# WISC_FULL_CS: normally distributed. pval=0.23627786338329315
# WISC_VCI_CS: normally distributed. pval=0.4142579734325409
# WISC_PR_CS: normally distributed. pval=0.06336785852909088
# WISC_WM_CS: normally distributed. pval=0.18117009103298187
# WISC_PS_CS: normally distributed. pval=0.1424320489168167
# srs.rrb: NOT normally distributed. pval=7.27403359768175e-10
# srs.sci: NOT normally distributed. pval=5.657304882333847e-06



## Check correlation between behavioral stat vs head motion
# pearson and spearman
for var in behav_data.columns[:7]:

    print(var)
    print(sp.stats.pearsonr(behav_data[var], behav_data['fd.m']))
    print(sp.stats.spearmanr(behav_data[var], behav_data['fd.m']))


    # r = sp.stats.pearsonr(behav_data[var], behav_data['fd.m'])[0]
    # p = sp.stats.pearsonr(behav_data[var], behav_data['fd.m'])[1]

    # plt.figure(figsize = (5,5))
    # g = sns.regplot(x=behav_data[var], y=behav_data['fd.m'], color='blue')
    # g.figure.tight_layout()
    # g.annotate('r = {0:.2f}'.format(r), xy = (0.7, 0.1), xycoords = 'axes fraction')
    # g.annotate('p = {0:.2f}'.format(p), xy = (0.7, 0.05), xycoords = 'axes fraction') ##
    # plt.savefig(os.path.join('dist', var + '_vs_fd.png'))
    # plt.close()

#####  #####


# WISC_FULL_CS
# (-0.25987592561851885, 0.014478648887248726)
# SpearmanrResult(correlation=-0.2494823894163985, pvalue=0.019073827165565438)
# WISC_VCI_CS
# (-0.20963316396793588, 0.049965954830521665)
# SpearmanrResult(correlation=-0.15596824177624358, pvalue=0.1467552258002097)
# WISC_PR_CS
# (-0.30249991649284874, 0.004175851587835593)
# SpearmanrResult(correlation=-0.29852429114100243, pvalue=0.004726778680059869)
# WISC_WM_CS
# (-0.08811584081307375, 0.4142882303431613)
# SpearmanrResult(correlation=-0.026869294467390475, pvalue=0.8037497277584218)
# WISC_PS_CS
# (-0.14146817251752836, 0.18859828875867224)
# SpearmanrResult(correlation=-0.16383118923994244, pvalue=0.12720331013623595)
# srs.rrb
# (0.12227547338772883, 0.2564132326712549)
# SpearmanrResult(correlation=0.08355653534931147, pvalue=0.43894426367078)
# srs.sci
# (0.14377381617291793, 0.18142238384489476)
# SpearmanrResult(correlation=0.09670161417181553, pvalue=0.3701064839327348)
