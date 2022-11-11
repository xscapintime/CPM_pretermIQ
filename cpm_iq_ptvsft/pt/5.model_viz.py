import os, glob
import numpy as np
import pandas as pd
import seaborn as sns
import scipy as sp
from matplotlib import pyplot as plt
from matplotlib import rcParams
plt.rcParams['figure.figsize'] = (6, 6)


# for empirical p-value
rfiles = sorted(glob.glob(os.path.join('./','*_fold_5_*_iter1000.csv')))


# files
# predicted table
files = sorted(glob.glob(os.path.join('./','*_pred.csv')))

# model performance
for rf,f in zip(rfiles, files):
    pred = pd.read_csv(f, index_col=0)
    null = pd.read_csv(rf, header=0)
    
    x = np.squeeze(pred.filter(regex=("obs")).astype(float))

    for tail in ['glm', 'pos', 'neg']:
        y = np.squeeze(pred.filter(regex=(tail)).astype(float))

        color = '#458B00' if tail == 'glm' else '#B22222' if tail == 'pos' else '#1874CD'
        print(tail + '_' + color)
        
        # scatter of obs vs pred
        g = sns.regplot(x=x.T.squeeze(), y=y.T.squeeze(), color=color, scatter_kws={'alpha':0.7})
        ax_min = min(min(g.get_xlim()), min(g.get_ylim()))
        ax_max = max(max(g.get_xlim()), max(g.get_ylim()))
        g.set_xlim(ax_min, ax_max)
        g.set_ylim(ax_min, ax_max)
        g.set_aspect('equal', adjustable='box')
        plt.xlabel("Obs. Full-scale IQ")
        plt.ylabel("Pred. Full-scale IQ")

        r = sp.stats.spearmanr(x,y)[0]
        g.annotate('r = {0:.3f}'.format(r), xy = (0.7, 0.1), xycoords = 'axes fraction')
        
        # empirical p value
        emp_pval = (sum(null[tail] >= r)+1)/(null.shape[0]+1)
        g.annotate('p = {0:.3f}'.format(emp_pval), xy = (0.7, 0.05), xycoords = 'axes fraction')

        fn = os.path.split(f)[1].replace('.csv', '') 
        plt.savefig(os.path.join('models', tail + '_' + fn + '.png'), dpi=300)
        plt.close()


        # null r value distribution
        g = sns.histplot(data=null, x=tail, bins=30, color=color, alpha=0.5)
        plt.xlabel("Correlation coefficient by permutation testing")
        plt.ylabel("Frequency")
        plt.axvline(r, color=color)

        plt.savefig(os.path.join('models', tail + '_' + fn + 'nullhist' + '.png'), dpi=300)
        plt.close()

