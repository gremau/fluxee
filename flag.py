"""
Functions for flagging and masking data. 

QA functions that can be applied are called from the qafunctions module.

"""

import pandas as pd
import numpy as np
from datetime import datetime
import imp
import qafunctions as qaf
imp.reload(qaf)
import pdb


def get_qafunction(flag):
    if 'qa_function' in flag:
        outfunc = getattr(qaf, flag['qa_function'])
        if 'qa_args' in flag:
            outargs = flag['qa_args']
        else:
            outargs = ''
    else:
        outfunc = getattr(qaf, 'dtrange_rm_all')
        outargs = ''
    return [outfunc, outargs]

def apply_qa_flags(df, flags):
    """
    Apply qa flags to a dataframe. There are two types of operations that can
    be performed on the input dataframe, depending on the QA flag:

    1. Transform data values based on qa flag input
    2. Mask data based on qa flag input
    
    These changes are logged in the flag array (df_flag) with a number
    corresponding to the qa flag in the site configuration files.

    Args:
        df      : input dataframe
        flags   : qa_flag dictionary from the site's ecoflux configuration dir
    Returns:
        Three pandas dataframes with identical dimensions to the input df
        df_new  : original data with any qa transformations applied
        df_mask : Mask dataframe containing boolean values (True = remove)
        df_flag : Flag dataframe with values corresponding to qa flags

    TODO - may want to add a flag for data already missing
    """
    # Make a copy to be a qa'd dataframe, aboolean array, and a flag dataframe
    df_new = df.copy()
    df_mask = pd.DataFrame(False, index=df.index, columns=df.columns)
    df_flag = pd.DataFrame(0, index=df.index, columns=df.columns)
    for i in flags.keys():
        if i == 0:
            raise ValueError('QA flag key cannot be zero (0)!')
        flag_cols_in = flags[i]['columns']
        st = flags[i]['start']
        en = flags[i]['end']
        if en is None:
            en = datetime.now()
        qafunc, qaargs = get_qafunction(flags[i])
        print('Apply QA flag {0}, using {1}.'.format(i, qafunc))
        if flag_cols_in=='all':
            # If "all" columns to be flagged select all
            colrange = df.columns
        else:
            # Or find dataframe columns matching those in qa_flags
            test = [any(s in var for s in flag_cols_in) for var in df.columns] 
            colrange = df.columns[test]
        # Get the index range to be flagged
        idxrange = np.logical_and(df.index >= st, df.index <= en)
        # Get the mask for flag i and set appropriate flag
        df_new, mask_i, rm = qafunc(df_new, idxrange, colrange, *qaargs)
        # Add mask_i to df_flag and to df_mask if data are to be masked
        df_flag[mask_i] = i
        if rm:
            df_mask = np.logical_or(df_mask, mask_i)

    # Rewrite df_flag column names
    df_flag.columns = df_flag.columns + '_flag'
    return df_new, df_mask, df_flag # df_new[df_mask]=np.nan will apply mask

def qa_dataframes(df, flags):
    """
    Get qa dataframes with flags appended and values masked

    Args:
        df: input dataframe
        flags: qa_flag dictionary from the site's ecoflux configuration dir
    Returns:
        df_qa       : QA'd dataframe with flags appended
        df_qa_masked: QA'd dataframe with flags appended and mask applied
    """
    df_qa, df_mask, df_flag = apply_qa_flags(df, flags)
    #df_qa_fl = pd.concat([df_qa, df_flag], axis=1)
    df_qa_masked = df_qa.copy()
    df_qa_masked[df_mask] = np.nan
    return df_qa, df_qa_masked, df_flag 
