import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
from matplotlib import pyplot as plt

# files
# predicted table
files = sorted(glob.glob(os.path.join('./','*_pred.csv')))

# model performance
for f in files:
    pred = pd.read_csv(f, index_col=0)
    x = np.squeeze(pred.filter(regex=("obs")).astype(float))

    for tail in ['glm', 'pos', 'neg']:
        y = np.squeeze(pred.filter(regex=(tail)).astype(float))

        color = 'green' if tail == 'glm' else 'red' if tail == 'pos' else 'blue'
        print(tail + '_' + color)
        
        g = sns.regplot(x=x.T.squeeze(), y=y.T.squeeze(), color=color)
        ax_min = min(min(g.get_xlim()), min(g.get_ylim()))
        ax_max = max(max(g.get_xlim()), max(g.get_ylim()))
        g.set_xlim(ax_min, ax_max)
        g.set_ylim(ax_min, ax_max)
        g.set_aspect('equal', adjustable='box')

        r = sp.stats.spearmanr(x,y)[0]
        g.annotate('r = {0:.2f}'.format(r), xy = (0.7, 0.1), xycoords = 'axes fraction')
        
        fn = os.path.split(f)[1].replace('.csv', '') 
        plt.savefig(os.path.join('models', tail + '_' + fn + '.png'))
        plt.close()


# add empirical p-value
## p-val
rvals = pd.read_csv('PC1_predvsobs_rvals_iter1000.txt', header=None)[0]
emp_pval = (sum(rvals >= rvals.mean())+1)/(1000+1)

rvals = pd.read_csv('PC2_predvsobs_rvals_iter1000.txt', header=None)[0]
emp_pval = (sum(rvals >= rvals.mean())+1)/(1000+1)


# predicted table
rfiles = sorted(glob.glob(os.path.join('./','*_iter1000.csv')))

for rf in rfiles:
    rv = pd.read_csv(rf, header=0)
