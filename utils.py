import pandas as pd
import numpy as np

############################## Loading data from Google Drive #############################

import requests
import openpyxl
import os
from io import BytesIO

def read_gsheet(sheet_id, tabname):
    url = "https://docs.google.com/spreadsheets/export?exportFormat=xlsx&id=" + sheet_id
    res = requests.get(url)
    data = BytesIO(res.content)
    xlsx = openpyxl.load_workbook(filename=data)
    return(pd.read_excel(data, sheet_name=tabname))

def read_dataset_log(tabname='Current Dataset Paths'):
    return(read_gsheet('1gJMdNBd7qkn_peRQEuz2ZHQaCO2EBMPCJnJZrjjzkBU', tabname))

def read_dataset_log_full(tabname='Dataset Paths', sheetid='1VyQURqVDzQQdqJfwn4A6Zg9Hwmpjy8iNfMnH4ltoJdc'):
    return(read_gsheet(sheetid, tabname))

def relpath_to_parent(start: str | None = None, target: str = "BCAT") -> str:
    """
    Return the relative path (string) that reaches the first ancestor directory
    named *target* (default 'BCAT') from *start* (default: current working dir).

    Raises FileNotFoundError if the target directory is not found.
    """
    # Absolute path of the starting directory
    start_path = os.path.abspath(start or os.getcwd())
    current = start_path
    depth = 0                         # how many levels we walk upward

    while True:
        if os.path.basename(current) == target:
            # Build ".." segments equal to the number of levels climbed
            if depth == 0:
                return os.curdir      # we're already inside BCAT
            return os.path.join(*[".."] * depth)

        parent = os.path.dirname(current)
        if parent == current:         # reached filesystem root
            raise FileNotFoundError(
                f"No parent directory named {target!r} above {start_path}"
            )

        current = parent
        depth += 1

############################## Plotting #############################

from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.gridspec import GridSpec
from datashader.mpl_ext import dsshow
import datashader as ds

def compute_smooth_scatter_color(x, y, z, n_neighbors=5):
    '''Compute local averages of z averaging over KNN to make smoothed scatter plot'''
    nbrs = NearestNeighbors(n_neighbors=5, algorithm='auto').fit(np.column_stack([x, y]))
    distances, indices = nbrs.kneighbors(np.column_stack([x, y]))
    z_avg = np.mean(z[indices], axis=1)
    return z_avg


def gate_biaxial(data, g1, g2, ind=None, vertical_gate=None, horizontal_gate=None, quadrant_gate=None,
                 labfontsize=9, plot_labeled=False, upper_only=False, ax=None, xlim=None, ylim=None,
                 xlabel=None, ylabel=None):
    '''Make a biaxial density plot of the columns g1 and g2 in the pandas.DataFrame data. Optionally
    gates the data points based on information in vertical_gate, horizontal_gate, or quadrant_gate'''
    cmap = plt.cm.hsv
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cmaplist = cmaplist[0:round(len(cmaplist)*0.7)]
    cmaplist.reverse()
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list('', cmaplist, cmap.N)
    
    if ax is None:
        fig = plt.figure(constrained_layout=True, figsize=(2., 2.), dpi=200)
        gs = GridSpec(1, 1, figure=fig, left=0.2, bottom=.2, right=.95, top=.92)
        ax = fig.add_subplot(gs[0])

    if ind  is None:
        ind = pd.Series(np.array([True]*data.shape[0]), data.index)
    
    
    x = data.loc[ind, g1]
    y = data.loc[ind, g2]
    
    dsshow(pd.DataFrame({'x':x, 'y':y}), ds.Point('x', 'y'), norm='eq_hist', cmap=cmap, ax = ax, aspect = 'auto')

    if xlim is None: xlim = ax.get_xlim()
    if ylim is None: ylim = ax.get_ylim()
    if xlabel is None: xlabel=g1
    if ylabel is None: ylabel=g2
    ax.set_xlabel(xlabel, fontsize=labfontsize)
    ax.set_ylabel(ylabel, fontsize=labfontsize)
    
    if vertical_gate is not None:
        ax.hlines(y=vertical_gate['vthresh'], xmin=xlim[0], xmax=xlim[1], linestyle='--', color='k', linewidth=1)
        ax.set_xlim(xlim)        
        initial_res = (data.loc[ind, g2]>vertical_gate['vthresh']).replace({True:vertical_gate['above_name'], False:vertical_gate['below_name']})
        final_res = ind.copy()
        final_res.loc[ind] = initial_res
        final_res.loc[~ind] = np.nan
    elif horizontal_gate is not None:
        ax.vlines(x=horizontal_gate['hthresh'], ymin=ylim[0], ymax=ylim[1], linestyle='--', color='k', linewidth=1)
        ax.set_ylim(ylim)
        
        initial_res = (data.loc[ind, g1]>horizontal_gate['hthresh']).replace({False:horizontal_gate['left_name'], True:horizontal_gate['right_name']})
        final_res = ind.copy()
        final_res.loc[ind] = initial_res
        final_res.loc[~ind] = np.nan
        
        
    elif quadrant_gate is not None:
        if not upper_only:
            ax.hlines(y=quadrant_gate['vthresh'], xmin=xlim[0]-1, xmax=xlim[1]+1, linestyle='--', color='k', linewidth=1)
            ax.set_xlim(xlim)
            ax.vlines(x=quadrant_gate['hthresh'], ymin=ylim[0]-1, ymax=ylim[1]+1, linestyle='--', color='k', linewidth=1)
            ax.set_ylim(ylim)
        else:
            ax.hlines(y=quadrant_gate['vthresh'], xmin=quadrant_gate['hthresh'], xmax=xlim[1]+1, linestyle='--', color='k', linewidth=1)
            ax.set_xlim(xlim)
            ax.vlines(x=quadrant_gate['hthresh'], ymin=quadrant_gate['vthresh'], ymax=ylim[1]+1, linestyle='--', color='k', linewidth=1)
            ax.set_ylim(ylim)            
            

        indh = data.loc[ind, g1]> quadrant_gate['hthresh']
        indv = data.loc[ind, g2]> quadrant_gate['vthresh']
        final_res = ind.replace(False,np.nan)
        for hval,vval, lab in [[False, False, quadrant_gate['ll']], [True, False, quadrant_gate['lr']], [False, True, quadrant_gate['ul']], [True, True, quadrant_gate['ur']]]:
            tolab = (indh==hval) & (indv==vval)
            final_res.loc[tolab.index[tolab]] = lab
    else:
        final_res = None
            
    if plot_labeled:
        fig = plt.figure(constrained_layout=True, figsize=(2.4, 2.), dpi=200)
        gs = gridspec.GridSpec(1, 1, figure=fig, left=0.2,
                       bottom=.2, right=.8, top=.92)
        ax2 = fig.add_subplot(gs[0])
        dat = pd.concat([x,y,final_res], axis=1)
        dat.columns = [g1, g2, 'label']
        for (k,g) in dat.groupby('label'):
            ax2.scatter(g[g1], g[g2], label=k, s=.5, edgecolor='None')
        ax2.legend(bbox_to_anchor=(1,1), markerscale=3, fontsize=6)

    return(final_res, ax)

############################## Matching GEPs between datasets #############################
def df_col_corr(X, Y):
    '''
    Compute pairwise Pearson correlation matrix of columns of X and Y. Returns
    R which is n_X_cols X n_Y_cols
    '''

    if X.isnull().any().any() or Y.isnull().any().any():
        raise Exception("Does not support NaNs currently")

    X_norm = X.subtract(X.mean(axis=0), axis=1)
    X_norm= X_norm.divide(X_norm.std(axis=0), axis=1)
    Y_norm = Y.subtract(Y.mean(axis=0), axis=1)
    Y_norm= Y_norm.divide(Y_norm.std(axis=0), axis=1)
    
    X_mask = np.ma.array(X_norm.values, mask=np.isnan(X_norm.values))
    Y_mask = np.ma.array(Y_norm.values, mask=np.isnan(Y_norm.values))
    R = np.ma.dot(X_mask.T, Y_mask) / (X.shape[0]-1)
    
    R = pd.DataFrame(R, index=X.columns, columns=Y.columns)
    return(R)


def match_columns(X, Y):
    '''
    Assign each column in Y to a distinct column in X via greedy search of max Pearson correlation.
    The # columns in X must be >= the # columns in Y
    '''

    if X.shape[1] < Y.shape[1]:
        raise Exception('X must have more columns than Y')
    
    R = df_col_corr(X, Y)
    mapping = R.unstack().reset_index().sort_values(by=0, ascending=False)
    mapping.columns = ['Y_columns', 'X_columns', 'R']

    used_X = []
    used_Y = []
    todrop = []
    for i in mapping.index:
        if mapping.at[i, 'X_columns'] in used_X:
            todrop.append(i)
        elif mapping.at[i, 'Y_columns'] in used_Y:
            todrop.append(i)
        else:
            used_X.append(mapping.at[i, 'X_columns'])
            used_Y.append(mapping.at[i, 'Y_columns'])  
            
    mapping = mapping.drop(todrop)
    unmatched = list(set(X.columns) - set(used_X))
    return(mapping, unmatched, R)


############################## Hypothesis testing #############################

from scipy.stats import ttest_rel

def ttest_paired_allcols(X, Y):
    '''
    Takes 2 pandas.DataFrames X and Y with the same number of rows and columns.
    Performs a paired t-test for each column of X and Y. 
    
    Returns
    ---------------------------------
    Ts - pandas.Series of T statistics
    Ps - pandas.Series of P-values
    '''
    
    Ts = []
    Ps = []
    for g in X.columns:
        T, P = ttest_rel(X[g], Y[g])
        Ts.append(T)
        Ps.append(P)

    Ts = pd.Series(Ts, index=X.columns)
    Ps = pd.Series(Ps, index=X.columns)
    return(Ts, Ps)


def permute_within_group(df, group_col, vartopermute):
    '''
    Performs permutation of the column vartopermute in the pandas.DataFrame df 
    but only permutes within the column group_col (which is assumed to contain discrete
    categories)

    Returns:
    permuted - a copy of df containing the additional column _permuted which can be used for
    downstream hypothesis testing
    '''
    permuted = df.groupby(group_col, group_keys=True).apply(lambda x: x.assign(_permuted=np.random.permutation(x[vartopermute]))).reset_index(drop=True)
    return(permuted)
