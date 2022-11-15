"""
Functions for connectome-based prodctive modeling.
Adaptation based on https://github.com/esfinn/cpm_tutorial and @shenUsingConnectomebasedPredictive2017
Added seed setting sampling.
Added robust correlation in feature selection.
Added partial correlaton in feature selection.
Replaced for-loop with apply-lambda method, maybe faster.
"""

import numpy as np
import scipy as sp
from matplotlib import pyplot as plt
#%matplotlib inline
import pandas as pd
import seaborn as sns
from pathlib import Path
import os
import warnings
from nilearn.plotting import plot_connectome
import pingouin as pg
import statsmodels.api as sm

# from sklearn.preprocessing import PolynomialFeatures
# from sklearn.linear_model import LinearRegression
# from sklearn.pipeline import make_pipeline
# from sklearn.metrics import mean_squared_error
warnings.simplefilter(action='ignore', category=FutureWarning)
os.environ['OUTDATED_IGNORE'] = '1'


# defualt path
top_dir = Path("./")
data_dir = top_dir/"fc_data/"

# Read in individual-subject connectivity matrices
def read_in_matrices(subj_list, file_suffix=None, data_dir=data_dir, zscore=False):
    """
    Reads in a set of individual-subject connectivity matrices stored in data_dir,
    
    Returns a dataframe that is subjects x edges (by vectorizing the upper triangle of each FC matrix).
    
    Assumes:
    - each matrix is stored in a separate file beginning with the subject ID, and
    - matrices are symmetric (squareform); i.e., for a parcellation with 268 nodes, matrices should be 268 x 268
    """
    
    all_fc_data = {}
            
    for subj in subj_list:
        # try to find this subject's matrix
        if file_suffix:
            file = [f for f in os.listdir(data_dir) if subj in f and file_suffix in f]
        else:
            file = [f for f in os.listdir(data_dir) if subj in f]
            
        # make sure there is one and only one file    
        if len(file) ==0:
            raise ValueError("No data found for subject {}".format(subj))
        if len(file) >1:
            raise ValueError("More than one matrix found for subject {}! Specify a suffix?".format(subj))
        
        # read it in and make sure it's symmetric and has reasonable dimensions
        tmp = np.loadtxt(data_dir / file[0])
        assert tmp.shape[0]==tmp.shape[1]>1, "Matrix seems to have incorrect dimensions: {}".format(tmp.shape)
        
        # take just the upper triangle and store it in a dictionary
        if ~zscore:
            all_fc_data[subj] = tmp[np.triu_indices_from(tmp, k=1)]
        if zscore:
            all_fc_data[subj] = sp.stats.zscore(tmp[np.triu_indices_from(tmp, k=1)])
        
    # Convert dictionary into dataframe
    all_fc_data = pd.DataFrame.from_dict(all_fc_data, orient='index')
    
    return all_fc_data



# CPM functions

# sampling, about 9:1 train:test
def mk_kfold_indices(subj_list, k = 10):
    """
    Splits list of subjects into k folds for cross-validation.
    """
    
    n_subs = len(subj_list)
    n_subs_per_fold = n_subs//k # floor integer for n_subs_per_fold

    indices = [[fold_no]*n_subs_per_fold for fold_no in range(k)] # generate repmat list of indices
    remainder = n_subs % k # figure out how many subs are left over
    remainder_inds = list(range(remainder))
    indices = [item for sublist in indices for item in sublist]    ## list comprehension
    [indices.append(ind) for ind in remainder_inds] # add indices for remainder subs

    assert len(indices)==n_subs, "Length of indices list does not equal number of subjects, something went wrong"

    np.random.shuffle(indices) # generator, shuffles in place

    return np.array(indices)


# fix the sampling by using seed
def mk_kfold_indices_seed(subj_list, k = 10, seed = 202208):
    """
    Splits list of subjects into k folds for cross-validation.
    """
    
    n_subs = len(subj_list)
    n_subs_per_fold = n_subs//k # floor integer for n_subs_per_fold

    indices = [[fold_no]*n_subs_per_fold for fold_no in range(k)] # generate repmat list of indices
    remainder = n_subs % k # figure out how many subs are left over
    remainder_inds = list(range(remainder))
    indices = [item for sublist in indices for item in sublist]    ## list comprehension
    [indices.append(ind) for ind in remainder_inds] # add indices for remainder subs

    assert len(indices)==n_subs, "Length of indices list does not equal number of subjects, something went wrong"

    np.random.seed(seed)
    indices = np.random.permutation(indices)
    # indices = np.random.shuffle(indices) # shuffles in place

    return np.array(indices)



def split_train_test(subj_list, indices, test_fold):
    """
    For a subj list, k-fold indices, and given fold, returns lists of train_subs and test_subs
    """

    train_inds = np.where(indices!=test_fold)
    test_inds = np.where(indices==test_fold)

    train_subs = []
    for sub in subj_list[train_inds]:
        train_subs.append(sub)

    test_subs = []
    for sub in subj_list[test_inds]:
        test_subs.append(sub)

    return (train_subs, test_subs)



def get_train_test_data(all_fc_data, train_subs, test_subs, behav_data, behav):
    
    """
    Extracts requested FC and behavioral data for a list of train_subs and test_subs
    """

    train_vcts = all_fc_data.loc[train_subs, :]
    test_vcts = all_fc_data.loc[test_subs, :]

    train_behav = behav_data.loc[train_subs, behav]

    return (train_vcts, train_behav, test_vcts)



# def get_train_test_cova(train_subs, test_subs, cova_data):

#     """
#     Extracts covariates data for a list of train_subs and test_subs
#     """
#     train_cova = cova_data.loc[train_subs, :]
#     # test_cova = cova_data.loc[test_subs, :]

#     return (train_cova)




################################# feature selection by coeffient #################################
#################################             PEARSON            #################################
#################################             SPEARMAN           #################################
#################################        ROBUST REGRESSION       #################################
def select_features(train_vcts, train_behav, r_thresh=0.2, corr_type='pearson', verbose=False):
    
    """
    Runs the CPM feature selection step: 
    - select edges by correlation coeffient
    - options: pearson, spearman, robust regression
    - correlates each edge with behavior, and returns a mask of edges that are correlated above some threshold, one for each tail (positive and negative)
    """

    assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

    # Correlate all edges with behav vector
    if corr_type =='pearson':
        cov = np.dot(train_behav.T - train_behav.mean(), train_vcts - train_vcts.mean(axis=0)) / (train_behav.shape[0]-1)
        corr = cov / np.sqrt(np.var(train_behav, ddof=1) * np.var(train_vcts, axis=0, ddof=1))
    
    elif corr_type =='spearman':
        ########################## original code  #######################
        corr = []
        for edge in train_vcts.columns:
            r_val = sp.stats.spearmanr(train_vcts.loc[:,edge], train_behav)[0]
            corr.append(r_val)
        #########################################################################
        # corr = train_vcts.apply(lambda x : sp.stats.spearmanr(x, train_behav)[0])
        corr = np.array(corr) # list bug
    
    elif corr_type == 'robust':
        # corr = train_vcts.apply(lambda x : sm.RLM(x, train_behav).fit().params[0])
        corr = []
        for edge in train_vcts.columns:
            r_val = sm.RLM(train_vcts.loc[:,edge], train_behav).fit().params[0]
            corr.append(r_val)
        corr = np.array(corr) # list bug


    # Define positive and negative masks
    mask_dict = {}
    mask_dict["pos"] = corr > r_thresh
    mask_dict["neg"] = corr < -r_thresh
    
    if verbose:
        print("Found ({}/{}) edges positively/negatively correlated with behavior in the training set".format(mask_dict["pos"].sum(), mask_dict["neg"].sum())) # for debugging

    return mask_dict, corr


#################################    feature selection by p-val  #################################
#################################             PEARSON            #################################
#################################             SPEARMAN           #################################
#################################        ROBUST REGRESSION       #################################
def select_features_pval(train_vcts, train_behav, p_thresh=0.05, corr_type='pearson', verbose=False):
    
    """
    Runs the CPM feature selection step: 
    - select edges by p-value
    - options: pearson, spearman, robust regression
    - correlates each edge with behavior, and returns a mask of edges that are correlated above some threshold, one for each tail (positive and negative)
    """

    assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

    # Correlate all edges with behav vector
    if corr_type =='pearson':
        # corr = train_vcts.apply(lambda x : sp.stats.pearsonr(x, train_behav)[0])
        # pval = train_vcts.apply(lambda x : sp.stats.pearsonr(x, train_behav)[1])
        corr = []
        pval = []
        for edge in train_vcts.columns:
            r_val = sp.stats.pearsonr(train_vcts.loc[:,edge], train_behav)[0]
            p_val = sp.stats.pearsonr(train_vcts.loc[:,edge], train_behav)[1]
            corr.append(r_val)
            pval.append(p_val)
        stat = {'corr':np.array(corr), 'pval' : np.array(pval)}


    elif corr_type =='spearman':
        # corr = train_vcts.apply(lambda x : sp.stats.spearmanr(x, train_behav)[0])
        # pval = train_vcts.apply(lambda x : sp.stats.spearmanr(x, train_behav)[1])
        corr = []
        pval = []
        for edge in train_vcts.columns:
            r_val = sp.stats.spearmanr(train_vcts.loc[:,edge], train_behav)[0]
            p_val = sp.stats.spearmanr(train_vcts.loc[:,edge], train_behav)[1]
            corr.append(r_val)
            pval.append(p_val)
        
        stat = {'corr':np.array(corr), 'pval' : np.array(pval)}


    elif corr_type == 'robust':
        # corr = train_vcts.apply(lambda x : sm.RLM(x, train_behav).fit().params[0])
        # pval = train_vcts.apply(lambda x : sm.RLM(x, train_behav).fit().bse[0])
        corr = []
        pval = []
        for edge in train_vcts.columns:
            r_val = sm.RLM(train_vcts.loc[:,edge], train_behav).fit().params[0]
            p_val = sm.RLM(train_vcts.loc[:,edge], train_behav).fit().bse[0]
            corr.append(r_val)
            pval.append(p_val)

        stat = {'corr':np.array(corr), 'pval' : np.array(pval)}


    # Define positive and negative masks
    mask_dict = {}

    mask_dict["pos"] = (stat['pval'] < p_thresh) & (stat['corr'] > 0)
    mask_dict["neg"] = (stat['pval'] < p_thresh) & (stat['corr'] < 0)


    if verbose:
        print("Found ({}/{}) edges positively/negatively correlated with behavior in the training set".format(mask_dict["pos"].sum(), mask_dict["neg"].sum())) # for debugging

    return mask_dict, corr #, pval




##################################################################################################
#################################        PARTIAL CORRELATION     #################################
##################################################################################################

##### to make it apply-able
def partcor_cpm(xedge, ybehav, cova_df, pcorr_type):
    """
    xedge: pandas series
    ybehav: pandas series
    cova_df: pandas dataframe
    """
    tmp = pd.concat([xedge, ybehav, cova_df.loc[xedge.index]], axis=1)
    pg_parcorr = pg.partial_corr(tmp, x=tmp.columns[0], y=tmp.columns[1], covar=cova_df.columns, method=pcorr_type)
    return pg_parcorr['r'][0], pg_parcorr['p-val'][0] 


################################# feature selection by coeffient #################################
def select_features_part(train_vcts, train_behav, covar, r_thresh=0.2, corr_type='pearson', verbose=False):
    
    """
    Runs the CPM feature selection step: 
    - select edges by correlation coeffient
    - partial correlation of each edge with behavior, and returns a mask of edges that are correlated above some threshold, one for each tail (positive and negative)
    - covar: pandas dataframe
    """

    assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

    # Correlate all edges with behav vector
    
    # parstats = train_vcts.apply(lambda x : partcor_cpm(xedge=x, ybehav=train_behav, cova_df=covar, pcorr_type=corr_type))
    # corr = parstats.loc[0,:]
    # pval = parstats.loc[1,:]

    corr = []
    # pval = []
    for edge in train_vcts.columns:
        r_val = partcor_cpm(xedge=train_vcts.loc[:,edge], ybehav=train_behav, cova_df=covar, pcorr_type=corr_type)[0]
        # p_val = partcor_cpm(xedge=train_vcts.loc[:,edge], ybehav=train_behav, cova_df=covar, pcorr_type=corr_type)[1]
        corr.append(r_val)
        # pval.append(p_val)

    corr = np.array(corr) # list bug
    
    ###############################################################################################
    # tmp_df = pd.concat([train_vcts, train_behav, train_cova], axis=1)
    # # tmp_df.head().pcorr() ## too big

    # corr = []
    # for edge in train_vcts.columns:
    #     r_val = pg.partial_corr(tmp_df, x=edge, y=behav, covar=covar, method=pcorr_type)['r'][0]
    #     corr.append(r_val)
    # corr = pd.Series(corr) # list bug
    ################################################################################################


    # Define positive and negative masks
    mask_dict = {}
    mask_dict["pos"] = corr > r_thresh
    mask_dict["neg"] = corr < -r_thresh
    
    if verbose:
        print("Found ({}/{}) edges positively/negatively correlated with behavior in the training set".format(mask_dict["pos"].sum(), mask_dict["neg"].sum())) # for debugging

    return mask_dict, corr


#################################    feature selection by p-val  #################################
def select_features_part_pval(train_vcts, train_behav, covar, p_thresh=0.05, corr_type='pearson', verbose=False):
    
    """
    Runs the CPM feature selection step: 
    - select edges by correlation coeffient
    - partial correlation of each edge with behavior, and returns a mask of edges that are correlated above some threshold, one for each tail (positive and negative)
    - covar: pandas dataframe
    """

    assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

    # Correlate all edges with behav vector
    
    # parstats = train_vcts.apply(lambda x : partcor_cpm(xedge=x, ybehav=train_behav, cova_df=covar, pcorr_type=corr_type))
    # corr = parstats.loc[0,:]
    # pval = parstats.loc[1,:]
    # stat = {'corr':np.array(corr), 'pval' : np.array(pval)}

    corr = []
    pval = []
    for edge in train_vcts.columns:
        r_val = partcor_cpm(xedge=train_vcts.loc[:,edge], ybehav=train_behav, cova_df=covar, pcorr_type=corr_type)[0]
        p_val = partcor_cpm(xedge=train_vcts.loc[:,edge], ybehav=train_behav, cova_df=covar, pcorr_type=corr_type)[1]
        corr.append(r_val)
        pval.append(p_val)
    stat = {'corr':np.array(corr), 'pval' : np.array(pval)}

    ##########################################################################################################
    # tmp_df = pd.concat([train_vcts, train_behav, train_cova], axis=1)
    # # tmp_df.head().pcorr() ## too big

    # corr = []
    # for edge in train_vcts.columns:
    #     r_val = pg.partial_corr(tmp_df, x=edge, y=behav, covar=covar, method=pcorr_type)['r'][0]
    #     corr.append(r_val)
    # corr = pd.Series(corr) # list bug
    ##########################################################################################################

    # Define positive and negative masks
    mask_dict = {}
    mask_dict["pos"] = (stat['pval'] < p_thresh) & (stat['corr'] > 0)
    mask_dict["neg"] = (stat['pval'] < p_thresh) & (stat['corr'] < 0)
    
    if verbose:
        print("Found ({}/{}) edges positively/negatively correlated with behavior in the training set".format(mask_dict["pos"].sum(), mask_dict["neg"].sum())) # for debugging

    return mask_dict, corr






### feature selection func, top edges ###
## resultes are not that different from selecting by p-val ##
## UNREFINED ##
def select_features_top(train_vcts, train_behav, perc_thresh=3, corr_type='pearson', verbose=False):
    
    """
    Runs the CPM feature selection step: 
    - correlates each edge with behavior, and returns a mask of edges that are correlated above some threshold, one for each tail (positive and negative)
    """

    assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

    # Correlate all edges with behav vector
    if corr_type =='pearson':
        # cov = np.dot(train_behav.T - train_behav.mean(), train_vcts - train_vcts.mean(axis=0)) / (train_behav.shape[0]-1)
        # corr = cov / np.sqrt(np.var(train_behav, ddof=1) * np.var(train_vcts, axis=0, ddof=1))
        corr = []
        pval = []
        for edge in train_vcts.columns:
            r_val = sp.stats.pearsonr(train_vcts.loc[:,edge], train_behav)[0]
            p_val = sp.stats.pearsonr(train_vcts.loc[:,edge], train_behav)[1]
            corr.append(r_val)
            pval.append(p_val)
        # stat = pd.DataFrame({'corr': corr, 'pval': pval})
        stat = {'corr':np.array(corr), 'pval' : np.array(pval)}


    elif corr_type =='spearman':
        corr = []
        pval = []
        for edge in train_vcts.columns:
            r_val = sp.stats.spearmanr(train_vcts.loc[:,edge], train_behav)[0]
            p_val = sp.stats.spearmanr(train_vcts.loc[:,edge], train_behav)[1]
            corr.append(r_val)
            pval.append(p_val)
        # stat = pd.DataFrame({'corr': corr, 'pval': pval})
        stat = {'corr':np.array(corr), 'pval' : np.array(pval)}


    # Define positive and negative masks
    mask_dict = {}
    # mask_dict["pos"] = corr > r_thresh
    # mask_dict["neg"] = corr < -r_thresh

    mask_dict["pos"] = (abs(stat['corr']) > np.percentile(abs(stat['corr']), 100-perc_thresh)) & (stat['corr'] > 0)
    mask_dict["neg"] = (abs(stat['corr']) > np.percentile(abs(stat['corr']), 100-perc_thresh)) & (stat['corr'] < 0)

    

    if verbose:
        print("Found ({}/{}) edges positively/negatively correlated with behavior in the training set".format(mask_dict["pos"].sum(), mask_dict["neg"].sum())) # for debugging

    return mask_dict, corr, pval





##### build model #####
def build_model(train_vcts, mask_dict, train_behav):
    """
    Builds a CPM model:
    - takes a feature mask, sums all edges in the mask for each subject, and uses simple linear regression to relate summed network strength to behavior
    """

    assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

    model_dict = {}

    # Loop through pos and neg tails
    X_glm = np.zeros((train_vcts.shape[0], len(mask_dict.items())))

    t = 0
    for tail, mask in mask_dict.items():
        X = train_vcts.values[:, mask].sum(axis=1)
        X_glm[:, t] = X
        y = train_behav
        (slope, intercept) = np.polyfit(X, y, 1)
        model_dict[tail] = (slope, intercept)
        t+=1

    X_glm = np.c_[X_glm, np.ones(X_glm.shape[0])]
    model_dict["glm"] = tuple(np.linalg.lstsq(X_glm, y, rcond=None)[0])

    return model_dict



def apply_model(test_vcts, mask_dict, model_dict):
    """
    Applies a previously trained linear regression model to a test set to generate predictions of behavior.
    """

    behav_pred = {}


    X_glm = np.zeros((test_vcts.shape[0], len(mask_dict.items())))

    # Loop through pos and neg tails
    t = 0
    for tail, mask in mask_dict.items():
        X = test_vcts.loc[:, mask].sum(axis=1)
        X_glm[:, t] = X

        slope, intercept = model_dict[tail]
        behav_pred[tail] = slope*X + intercept
        t+=1

    X_glm = np.c_[X_glm, np.ones(X_glm.shape[0])]
    behav_pred["glm"] = np.dot(X_glm, model_dict["glm"])

    return behav_pred




###### multivariable polynomial ######
#### unfinished ####
#### seems can't change the linear regression part as the relationship should be liner ####
# def build_model_cova(train_vcts, mask_dict, train_behav, train_cova, degree=2):
#     """
#     Builds a CPM model:
#     - takes a feature mask, sums all edges in the mask for each subject, and uses simple linear regression to relate summed network strength to behavior
#     """

#     assert train_vcts.index.equals(train_behav.index), "Row indices of FC vcts and behavior don't match!"

#     model_dict = {}

#     # Loop through pos and neg tails
#     X_glm = np.zeros((train_vcts.shape[0], len(mask_dict.items())))

#     t = 0
#     for tail, mask in mask_dict.items():
#         X = train_vcts.values[:, mask].sum(axis=1)
#         X_cova = np.hstack((X.reshape(train_vcts.shape[0],1), train_cova[['sex', 'age22']].values))
#         X_glm[:, t] = X
#         y = train_behav

#         poly_model = PolynomialFeatures(degree=degree)
#         poly_model = make_pipeline(PolynomialFeatures(), LinearRegression())
#         poly_model.fit(X_cova, y)


#         # (slope, intercept) = np.polyfit(X, y, 1)
#         # model_dict[tail] = (slope, intercept)
#         t+=1

#     X_glm = np.c_[X_glm, np.ones(X_glm.shape[0])]
#     model_dict["glm"] = tuple(np.linalg.lstsq(X_glm, y, rcond=None)[0])

#     return model_dict








# use r_thresh
def cpm_wrapper(all_fc_data, all_behav_data, behav, k=10, **cpm_kwargs):
    
    assert all_fc_data.index.equals(all_behav_data.index), "Row (subject) indices of FC vcts and behavior don't match!"

    subj_list = all_fc_data.index # get subj_list from df index
    
    indices = mk_kfold_indices(subj_list, k=k)
    
    # Initialize df for storing observed and predicted behavior
    col_list = []
    for tail in ["pos", "neg", "glm"]:
        col_list.append(behav + " predicted (" + tail + ")")
    col_list.append(behav + " observed")
    behav_obs_pred = pd.DataFrame(index=subj_list, columns = col_list)
    
    # Initialize array for storing feature masks
    n_edges = all_fc_data.shape[1]
    all_masks = {}
    all_masks["pos"] = np.zeros((k, n_edges))
    all_masks["neg"] = np.zeros((k, n_edges))
    
    for fold in range(k):
        print("doing fold {}".format(fold))
        train_subs, test_subs = split_train_test(subj_list, indices, test_fold=fold)
        train_vcts, train_behav, test_vcts = get_train_test_data(all_fc_data, train_subs, test_subs, all_behav_data, behav=behav)
        mask_dict, corr = select_features(train_vcts, train_behav, **cpm_kwargs)
        all_masks["pos"][fold,:] = mask_dict["pos"]
        all_masks["neg"][fold,:] = mask_dict["neg"]
        model_dict = build_model(train_vcts, mask_dict, train_behav)
        behav_pred = apply_model(test_vcts, mask_dict, model_dict)
        for tail, predictions in behav_pred.items():
            behav_obs_pred.loc[test_subs, behav + " predicted (" + tail + ")"] = predictions
            
    behav_obs_pred.loc[subj_list, behav + " observed"] = all_behav_data[behav]
    
    return behav_obs_pred, all_masks, corr


# use p_thresh
def cpm_wrapper_pval(all_fc_data, all_behav_data, behav, k=10, **cpm_kwargs):
    
    assert all_fc_data.index.equals(all_behav_data.index), "Row (subject) indices of FC vcts and behavior don't match!"

    subj_list = all_fc_data.index # get subj_list from df index
    
    indices = mk_kfold_indices(subj_list, k=k)
    
    # Initialize df for storing observed and predicted behavior
    col_list = []
    for tail in ["pos", "neg", "glm"]:
        col_list.append(behav + " predicted (" + tail + ")")
    col_list.append(behav + " observed")
    behav_obs_pred = pd.DataFrame(index=subj_list, columns = col_list)
    
    # Initialize array for storing feature masks
    n_edges = all_fc_data.shape[1]
    all_masks = {}
    all_masks["pos"] = np.zeros((k, n_edges))
    all_masks["neg"] = np.zeros((k, n_edges))
    
    for fold in range(k):
        print("doing fold {}".format(fold))
        train_subs, test_subs = split_train_test(subj_list, indices, test_fold=fold)
        train_vcts, train_behav, test_vcts = get_train_test_data(all_fc_data, train_subs, test_subs, all_behav_data, behav=behav)
        mask_dict, corr = select_features_pval(train_vcts, train_behav, **cpm_kwargs)
        all_masks["pos"][fold,:] = mask_dict["pos"]
        all_masks["neg"][fold,:] = mask_dict["neg"]
        model_dict = build_model(train_vcts, mask_dict, train_behav)
        behav_pred = apply_model(test_vcts, mask_dict, model_dict)
        for tail, predictions in behav_pred.items():
            behav_obs_pred.loc[test_subs, behav + " predicted (" + tail + ")"] = predictions
            
    behav_obs_pred.loc[subj_list, behav + " observed"] = all_behav_data[behav]
    
    return behav_obs_pred, all_masks, corr # to see the selected features




## samping with seed
## use r_thresh
def cpm_wrapper_seed(all_fc_data, all_behav_data, behav, k=10, seed=202208, **cpm_kwargs):
    
    assert all_fc_data.index.equals(all_behav_data.index), "Row (subject) indices of FC vcts and behavior don't match!"

    subj_list = all_fc_data.index # get subj_list from df index
    
    indices = mk_kfold_indices_seed(subj_list, k=k, seed=seed)
    
    # Initialize df for storing observed and predicted behavior
    col_list = []
    for tail in ["pos", "neg", "glm"]:
        col_list.append(behav + " predicted (" + tail + ")")
    col_list.append(behav + " observed")
    behav_obs_pred = pd.DataFrame(index=subj_list, columns = col_list)
    
    # Initialize array for storing feature masks
    n_edges = all_fc_data.shape[1]
    all_masks = {}
    all_masks["pos"] = np.zeros((k, n_edges))
    all_masks["neg"] = np.zeros((k, n_edges))
    
    for fold in range(k):
        print("doing fold {}".format(fold))
        train_subs, test_subs = split_train_test(subj_list, indices, test_fold=fold)
        train_vcts, train_behav, test_vcts = get_train_test_data(all_fc_data, train_subs, test_subs, all_behav_data, behav=behav)
        mask_dict, corr = select_features(train_vcts, train_behav, **cpm_kwargs)
        all_masks["pos"][fold,:] = mask_dict["pos"]
        all_masks["neg"][fold,:] = mask_dict["neg"]
        model_dict = build_model(train_vcts, mask_dict, train_behav)
        behav_pred = apply_model(test_vcts, mask_dict, model_dict)
        for tail, predictions in behav_pred.items():
            behav_obs_pred.loc[test_subs, behav + " predicted (" + tail + ")"] = predictions
            
    behav_obs_pred.loc[subj_list, behav + " observed"] = all_behav_data[behav]
    
    return behav_obs_pred, all_masks, corr


## samping with seed
## use p_thresh
def cpm_wrapper_seed_pval(all_fc_data, all_behav_data, behav, k=10, seed=202208, **cpm_kwargs):
    
    assert all_fc_data.index.equals(all_behav_data.index), "Row (subject) indices of FC vcts and behavior don't match!"

    subj_list = all_fc_data.index # get subj_list from df index
    
    indices = mk_kfold_indices_seed(subj_list, k=k, seed=seed)
    
    # Initialize df for storing observed and predicted behavior
    col_list = []
    for tail in ["pos", "neg", "glm"]:
        col_list.append(behav + " predicted (" + tail + ")")
    col_list.append(behav + " observed")
    behav_obs_pred = pd.DataFrame(index=subj_list, columns = col_list)
    
    # Initialize array for storing feature masks
    n_edges = all_fc_data.shape[1]
    all_masks = {}
    all_masks["pos"] = np.zeros((k, n_edges))
    all_masks["neg"] = np.zeros((k, n_edges))
    
    for fold in range(k):
        print("doing fold {}".format(fold))
        train_subs, test_subs = split_train_test(subj_list, indices, test_fold=fold)
        train_vcts, train_behav, test_vcts = get_train_test_data(all_fc_data, train_subs, test_subs, all_behav_data, behav=behav)
        mask_dict, corr = select_features_pval(train_vcts, train_behav, **cpm_kwargs)
        all_masks["pos"][fold,:] = mask_dict["pos"]
        all_masks["neg"][fold,:] = mask_dict["neg"]
        model_dict = build_model(train_vcts, mask_dict, train_behav)
        behav_pred = apply_model(test_vcts, mask_dict, model_dict)
        for tail, predictions in behav_pred.items():
            behav_obs_pred.loc[test_subs, behav + " predicted (" + tail + ")"] = predictions
            
    behav_obs_pred.loc[subj_list, behav + " observed"] = all_behav_data[behav]
    
    return behav_obs_pred, all_masks, corr



## samping with seed
## use partial correlation r_thresh
def cpm_wrapper_seed_part(all_fc_data, all_behav_data, behav, covar, k=10, seed=202208, **cpm_kwargs):
    
    assert all_fc_data.index.equals(all_behav_data.index), "Row (subject) indices of FC vcts and behavior don't match!"

    subj_list = all_fc_data.index # get subj_list from df index
    
    indices = mk_kfold_indices_seed(subj_list, k=k, seed=seed)
    
    # Initialize df for storing observed and predicted behavior
    col_list = []
    for tail in ["pos", "neg", "glm"]:
        col_list.append(behav + " predicted (" + tail + ")")
    col_list.append(behav + " observed")
    behav_obs_pred = pd.DataFrame(index=subj_list, columns = col_list)
    
    # Initialize array for storing feature masks
    n_edges = all_fc_data.shape[1]
    all_masks = {}
    all_masks["pos"] = np.zeros((k, n_edges))
    all_masks["neg"] = np.zeros((k, n_edges))
    
    for fold in range(k):
        print("doing fold {}".format(fold))
        train_subs, test_subs = split_train_test(subj_list, indices, test_fold=fold)
        train_vcts, train_behav, test_vcts = get_train_test_data(all_fc_data, train_subs, test_subs, all_behav_data, behav=behav)
        # train_cova = get_train_test_cova(train_subs, test_subs, cova_data)
        mask_dict, corr = select_features_part(train_vcts, train_behav, covar, **cpm_kwargs)
        all_masks["pos"][fold,:] = mask_dict["pos"]
        all_masks["neg"][fold,:] = mask_dict["neg"]

        model_dict = build_model(train_vcts, mask_dict, train_behav) ##
        behav_pred = apply_model(test_vcts, mask_dict, model_dict) ##
        
        for tail, predictions in behav_pred.items():
            behav_obs_pred.loc[test_subs, behav + " predicted (" + tail + ")"] = predictions
            
    behav_obs_pred.loc[subj_list, behav + " observed"] = all_behav_data[behav]
    
    return behav_obs_pred, all_masks, corr


## samping with seed
## use partial correlation p_thresh
def cpm_wrapper_seed_part_pval(all_fc_data, all_behav_data, behav, covar, k=10, seed=202208, **cpm_kwargs):
    
    assert all_fc_data.index.equals(all_behav_data.index), "Row (subject) indices of FC vcts and behavior don't match!"

    subj_list = all_fc_data.index # get subj_list from df index
    
    indices = mk_kfold_indices_seed(subj_list, k=k, seed=seed)
    
    # Initialize df for storing observed and predicted behavior
    col_list = []
    for tail in ["pos", "neg", "glm"]:
        col_list.append(behav + " predicted (" + tail + ")")
    col_list.append(behav + " observed")
    behav_obs_pred = pd.DataFrame(index=subj_list, columns = col_list)
    
    # Initialize array for storing feature masks
    n_edges = all_fc_data.shape[1]
    all_masks = {}
    all_masks["pos"] = np.zeros((k, n_edges))
    all_masks["neg"] = np.zeros((k, n_edges))
    
    for fold in range(k):
        print("doing fold {}".format(fold))
        train_subs, test_subs = split_train_test(subj_list, indices, test_fold=fold)
        train_vcts, train_behav, test_vcts = get_train_test_data(all_fc_data, train_subs, test_subs, all_behav_data, behav=behav)
        # train_cova = get_train_test_cova(train_subs, test_subs, cova_data)
        mask_dict, corr = select_features_part_pval(train_vcts, train_behav, covar, **cpm_kwargs)
        all_masks["pos"][fold,:] = mask_dict["pos"]
        all_masks["neg"][fold,:] = mask_dict["neg"]

        model_dict = build_model(train_vcts, mask_dict, train_behav) ##
        behav_pred = apply_model(test_vcts, mask_dict, model_dict) ##
        
        for tail, predictions in behav_pred.items():
            behav_obs_pred.loc[test_subs, behav + " predicted (" + tail + ")"] = predictions
            
    behav_obs_pred.loc[subj_list, behav + " observed"] = all_behav_data[behav]
    
    return behav_obs_pred, all_masks, corr




### alternative feature selection. seledt top x% edges, unrefined function
def cpm_wrapper_alt(all_fc_data, all_behav_data, behav, k=10, **cpm_kwargs):
    
    assert all_fc_data.index.equals(all_behav_data.index), "Row (subject) indices of FC vcts and behavior don't match!"

    subj_list = all_fc_data.index # get subj_list from df index
    
    indices = mk_kfold_indices(subj_list, k=k)
    
    # Initialize df for storing observed and predicted behavior
    col_list = []
    for tail in ["pos", "neg", "glm"]:
        col_list.append(behav + " predicted (" + tail + ")")
    col_list.append(behav + " observed")
    behav_obs_pred = pd.DataFrame(index=subj_list, columns = col_list)
    
    # Initialize array for storing feature masks
    n_edges = all_fc_data.shape[1]
    all_masks = {}
    all_masks["pos"] = np.zeros((k, n_edges))
    all_masks["neg"] = np.zeros((k, n_edges))
    
    for fold in range(k):
        print("doing fold {}".format(fold))
        train_subs, test_subs = split_train_test(subj_list, indices, test_fold=fold)
        train_vcts, train_behav, test_vcts = get_train_test_data(all_fc_data, train_subs, test_subs, all_behav_data, behav=behav)
        mask_dict, corr, pval = select_features_top(train_vcts, train_behav, **cpm_kwargs)
        all_masks["pos"][fold,:] = mask_dict["pos"]
        all_masks["neg"][fold,:] = mask_dict["neg"]
        model_dict = build_model(train_vcts, mask_dict, train_behav)
        behav_pred = apply_model(test_vcts, mask_dict, model_dict)
        for tail, predictions in behav_pred.items():
            behav_obs_pred.loc[test_subs, behav + " predicted (" + tail + ")"] = predictions
            
    behav_obs_pred.loc[subj_list, behav + " observed"] = all_behav_data[behav]
    
    return behav_obs_pred, all_masks, corr, pval # to see the selected features


##### plot #####
def plot_predictions(behav_obs_pred, tail="glm"):
    # x = behav_obs_pred.filter(regex=("obs")).astype(float)
    # y = behav_obs_pred.filter(regex=(tail)).astype(float)

    ## dimension bug
    ## https://stackoverflow.com/a/71546253/14498100
    x = np.squeeze(behav_obs_pred.filter(regex=("obs")).astype(float))
    y = np.squeeze(behav_obs_pred.filter(regex=(tail)).astype(float))


    g = sns.regplot(x=x.T.squeeze(), y=y.T.squeeze(), color='gray')
    ax_min = min(min(g.get_xlim()), min(g.get_ylim()))
    ax_max = max(max(g.get_xlim()), max(g.get_ylim()))
    g.set_xlim(ax_min, ax_max)
    g.set_ylim(ax_min, ax_max)
    g.set_aspect('equal', adjustable='box')
    
    r = sp.stats.pearsonr(x,y)[0]
    # p = sp.stats.pearsonr(x,y)[1] ## parametric p-value, not correct to use
    g.annotate('r = {0:.2f}'.format(r), xy = (0.7, 0.1), xycoords = 'axes fraction')
    # g.annotate('p = {0:.2f}'.format(p), xy = (0.7, 0.05), xycoords = 'axes fraction') ##
    
    return g


## visualize edges
# default shen 268, change to 374
# shen268_coords = pd.read_csv("../data/coords/shen268_coords.csv", index_col="NodeNo")
coords = pd.read_csv("../data/coords/MMP_MNI_374_UCCHILD_coords.txt", index_col=None, header=None, sep=" ")

def plot_consistent_edges(all_masks, tail, thresh = 1., color='gray', node_coords=coords):
    
    edge_frac = (all_masks[tail].sum(axis=0))/(all_masks[tail].shape[0])
    print("For the {} tail, {} edges were selected in at least {}% of folds".format(tail, (edge_frac>=thresh).sum(), thresh*100))
    edge_frac_square = sp.spatial.distance.squareform(edge_frac)

    node_mask = np.amax(edge_frac_square, axis=0) >= thresh # find nodes that have at least one edge that passes the threshold
    node_size = edge_frac_square.sum(axis=0)*node_mask*20 # size nodes based on how many suprathreshold edges they have

    plot_connectome(adjacency_matrix=edge_frac_square, edge_threshold=thresh,
                    node_color = color,
                    node_coords=node_coords, node_size=node_size, ## need define coords first
                    display_mode= 'lzry',
                    edge_kwargs={"linewidth": 1, 'color': color})
