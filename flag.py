import pandas as pd
import numpy as np
import pdb

def flag_dataframe(df, flags):
    """
    """
    # Make a copy to be a boolean array and a flag dataframe
    df_bool = pd.DataFrame(False, index=df.index, columns=df.columns)
    cols_fl = df.columns + '_flag'
    df_fl = pd.DataFrame(False, index=df.index, columns=cols_fl)
    for i in flags.keys():
        flagcols = flags[i]['columns']
        st = flags[i]['start']
        en = flags[i]['end']
        if flagcols=='all':
            bool_cols = df_bool.columns
            flagcols=cols_fl
        else:
            bool_cols = flagcols
            flagcols = [s + '_flag' for s in flagcols]
        flagrange = np.logical_and(df.index >= st, df.index <= en)
        df_bool.loc[flagrange, bool_cols] = True
        df_fl.loc[flagrange, flagcols] = True
    return df_bool, df_fl # now df[!df_bool] should remove values

def qa_dataframe(df, flags):
    df_qa = df.copy()
    df_bool, df_flags = flag_dataframe(df, flags)
    df_qa[df_bool] = np.nan
    df_qa = pd.concat([df_qa, df_flags], axis=1)
    return df_qa
