import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def gapfill_series(s_withgaps, s_gapfiller, makeplots=False):
    '''
    Fill nans in one series with data from another

    IN:
        s_withgaps  : (pandas series) series containing gaps
        s_gapfiller : (pandas series) series (continuous) used to fill gaps 
    
    OUT:
        s_filled    : pandas series) series with gaps filled
    '''

    s_filled = pd.DataFrame(s_withgaps.copy())

    if s_withgaps.index.equals(s_gapfiller.index):
        gapfill = np.isnan(s_withgaps)
        s_filled.columns = [s_filled.columns[0] + '_gf']
        s_filled.loc[gapfill, s_filled.columns] = s_gapfiller[gapfill]
        s_filled[s_filled.columns[0] + 'FLAG'] = gapfill
    else:
        print('Error - indices are not the same')
    if makeplots:
        f, ax = plt.subplots(1)
        ax.set_title(s_withgaps.name)
        ax.plot(s_filled.iloc[:,0], 'og', mfc='w')
        ax.plot(s_gapfiller, '.r', ms=3)
        ax.plot(s_withgaps, '.b', ms=3)
        ax.set_ylabel('Y')
        ax.legend(['filled series', 'gapfill data', 'original series'], ncol=3)

    return s_filled

def resample_dataframe_by_col( df, freq='1D', avg_cols=[ 'TA_F'],
        min_cols=[ 'TA_F', 'VPD_F' ], max_cols=['LE_F', 'H_F'],
        sum_cols=[ 'P_F' ]):
    """
    Resample a dataframe, specifying resample statistic for given columns.

    Args:
        df          : pandas DataFrame (usually derived from datalog file)
        freq        : frequency to resample to (default daily)
        avg_cols    : list of header names (strings) to average
        minmax_cols : list of header names (strings) to convert to min/max
        int_cols    : list of header names (strings) to integrate (*1800)
        sum_cols    : list of header names (strings) to sum

    Return:
        df_resamp   : pandas dataframe with data at new frequency
    """

    # Subset site data into summable, averagable, etc data
    df_sum = df[ sum_cols ]
    #df_int = df[ int_cols ]*1800
    df_avg = df[ avg_cols ]
    df_min = df[ min_cols ]
    df_max = df[ max_cols ]
    
    # Resample to daily using sum or mean
    sums_resamp = df_sum.resample( freq ).sum()
    # Sometimes only C fluxes are provided, handle exceptions
    try: 
        avg_resamp = df_avg.resample( freq ).mean()
        min_resamp = df_min.resample( freq ).min()
        max_resamp = df_max.resample( freq ).max()
    
        # Rename the int columns
        #for i in int_cols:
        #    int_resamp.rename(columns={ i:i + '_int'}, inplace=True)

        # Rename the avg columns
        for i in avg_cols:
            avg_resamp.rename(columns={ i:i + '_avg'}, inplace=True)
        
        # Rename the min/max columns
        for i in min_cols:
            min_resamp.rename(columns={ i:i + '_min'}, inplace=True)
        for i in max_cols:
            max_resamp.rename(columns={ i:i + '_max'}, inplace=True)
    except:
        #int_resamp = pd.DataFrame(index=sums_resamp.index)
        avg_resamp = pd.DataFrame(index=sums_resamp.index)
        min_resamp = pd.DataFrame(index=sums_resamp.index)
        max_resamp = pd.DataFrame(index=sums_resamp.index)

    # Rename the sum columns
    for i in sum_cols:
        sums_resamp.rename(columns={ i:i + '_sum'}, inplace=True)


    # Put to dataframes back together
    df_resamp = pd.concat( [ sums_resamp, avg_resamp,
        min_resamp, max_resamp ], axis=1 )

    return df_resamp


    
