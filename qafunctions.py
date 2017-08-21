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

def mask_by_datetime(df, idxrange, colrange):
    """
    Mask all matching idxrange and colrange
    """
    mask = pd.DataFrame(False, index=df.index, columns=df.columns)
    mask.loc[idxrange, colrange] = True
    return [df, mask, True]

def mask_by_comparison(df, idxrange, colrange, direction, cval):
    """
    Mask values in matching idxrange and colrange AND colrange variables
    are above/below cval 
    """
    mask = pd.DataFrame(False, index=df.index, columns=df.columns)
    for c in colrange:
        if direction=='above':
            idxrange_th = np.logical_and(idxrange, df[c] > cval)
        elif direction=='below':
            idxrange_th = np.logical_and(idxrange, df[c] < cval)
        elif direction=='equals':
            idxrange_th = np.logical_and(idxrange, df[c] == cval)
        else:
            raise ValueError('Invalid comparison (above, below, equals)')
        mask.loc[idxrange_th, c] = True
    return [df, mask, True]

def dtrange_rm_var_threshold(df, idxrange, colrange, indvar,
        direction, cval):
    """
    Mask values in matching idxrange and colrange AND where an independent
    variable (indvar) is above/below cval 
    """
    mask = pd.DataFrame(False, index=df.index, columns=df.columns)
    if direction=='above':
        idxrange_thv = np.logical_and(idxrange, df[indvar] > cval)
    elif direction=='below':
        idxrange_thv = np.logical_and(idxrange, df[indvar] < cval)
    elif direction=='equals':
        idxrange_thv = np.logical_and(idxrange, df[indvar] == cval)
    else:
        raise ValueError('Invalid comparison (above, below, equals)')
    mask.loc[idxrange_thv, colrange] = True
    return [df, mask, True]

