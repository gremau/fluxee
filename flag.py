import yaml
import pandas as pd
import numpy as np
import pdb

def readflags(yamlfile):
    stream = open(yamlfile, 'r')
    flags = yaml.load(stream)
    return flags

def add_flag_cols(df):
    cols = df.columns
    cols_fl = cols + '_flag'
    df_fl = pd.concat([df, pd.DataFrame(0, index=df.index, columns=cols_fl)],
            axis=1)
    return df_fl

def set_flags(flags, df):
    cols = df.columns
    cols_fl = cols + '_flag'
    df_fl = pd.DataFrame(False, index=df.index, columns=cols_fl)
    for i in flags.keys():
        flagcols = flags[i]['columns']
        st = flags[i]['start']
        en = flags[i]['end']
        if flagcols=='all':
            flagcols=cols_fl
        else:
            flagcols = flagcols + '_flag'
        flagrange = np.logical_and(df.index >= st, df.index <= en)
        df_fl.loc[flagrange, flagcols] = True
    return df_fl # now df.iloc[df_fl] should remove values
