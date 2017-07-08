"""
Functions that can be called to qa a dataframe. These are generally called from
the apply_qa_flags function in the flag module. Functions must return a 
dataframe (often the same as the input), a boolean array mask indicating which
dataframe values the qa flag points to, and a boolean value indicating whether
the flagged data should be masked (True = flag for removal).
"""

import pandas as pd
import numpy as np
import pdb

def dtrange_rm_all(df, idxrange, colrange):
    """
    Mask all matching idxrange and colrange
    """
    mask = pd.DataFrame(False, index=df.index, columns=df.columns)
    mask.loc[idxrange, colrange] = True
    return [df, mask, True]

def dtrange_rm_threshold(df, idxrange, colrange, direction, threshval):
    """
    Mask values in matching idxrange and colrange AND colrange variables
    are above/below threshval 
    """
    mask = pd.DataFrame(False, index=df.index, columns=df.columns)
    for c in colrange:
        if direction=='above':
            idxrange_th = np.logical_and(idxrange, df[c] > threshval)
        elif direction=='below':
            idxrange_th = np.logical_and(idxrange, df[c] < threshval)
        else:
            raise ValueError('Incorrect direction')
        mask.loc[idxrange_th, c] = True
    return [df, mask, True]

def dtrange_rm_var_threshold(df, idxrange, colrange, threshvar,
        direction, threshval):
    """
    Mask values in matching idxrange and colrange AND where another variable
    (threshvar) is above/below threshval 
    """
    mask = pd.DataFrame(False, index=df.index, columns=df.columns)
    if direction=='above':
        idxrange_thv = np.logical_and(idxrange, df[threshvar] > threshval)
    elif direction=='below':
        idxrange_thv = np.logical_and(idxrange, df[threshvar] < threshval)
    else:
        raise ValueError('Incorrect direction')
    mask.loc[idxrange_thv, colrange] = True
    return [df, mask, True]

