import pandas as pd
import numpy as np
import pdb

def flag_dataframe(df, flags):
    """
    Flag dataframe using a dictionary of flags (ecoflux format).

    Args:
        df      : input dataframe
        flags   : qa_flag dictionary from the site's ecoflux configuration dir
    Returns:
        df_bool : pandas dataframe containing boolean values (True = flagged)
        df_flag : pandas dataframe containing flag values (all 1 currently)
    """
    # Make a copy to be a boolean array and a flag dataframe
    df_bool = pd.DataFrame(False, index=df.index, columns=df.columns)
    df_flag = pd.DataFrame(0, index=df.index, columns=df.columns+'_flag')
    for i in flags.keys():
        flag_cols = flags[i]['columns']
        st = flags[i]['start']
        en = flags[i]['end']
        if flag_cols=='all':
            # If "all" columns to be flagged select all
            bool_cols = df_bool.columns
            flag_cols = df_flag.columns
        else:
            # Or find dataframe columns matching those in qa_flags
            test = [any(s in var for s in flag_cols) for var in df.columns] 
            bool_cols = df.columns[test]
            flag_cols = bool_cols + '_flag'
        # Get the index range to be flagged and set boolean and flag columns
        # to the appropriate value
        idxrange = np.logical_and(df.index >= st, df.index <= en)
        df_bool.loc[idxrange, bool_cols] = True
        df_flag.loc[idxrange, flag_cols] = 1
    return df_bool, df_flag # now df[!df_bool] should remove values

def qa_dataframe(df, flags):
    """
    Get a qa dataframe with bad data removed and flags appended

    Args:
        df: input dataframe
        flags: qa_flag dictionary from the site's ecoflux configuration dir
    Returns:
        df_qa: fully qa'd dataframe with bad data removed and flags appended
    """
    df_qa = df.copy()
    df_bool, df_flag = flag_dataframe(df, flags)
    df_qa[df_bool] = np.nan
    df_qa = pd.concat([df_qa, df_flag], axis=1)
    return df_qa
