import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def measurement_h_v_dict(cols, meas, str_exclude=None):
    '''
    Extract measurement horizontal and vertical location configuration. This
    relies on column names to follow the 'MEASTYPE_H_V_R' convention
    (Horiz, Vert, Rep)

    IN:
        cols: (string list) column names for each measurement
        meas: (string) the measurement type as represented in column names
        str_exclude: (string) exclude columns containing this string
    OUT:
        hv_dict: (dict) dict with vertical location list for each meas_horiz key
    '''
    # Count underscores in meas variable
    meas_uscores = meas.count('_')
    # Match column names with meas variable and split into H and V
    if str_exclude is not None:
        meas_cols = [c for c in cols if meas + '_' in c and 
                str_exclude not in c]
    else:
        meas_cols = [c for c in cols if meas + '_' in c]
    horiz = [n.split('_')[1 + meas_uscores] for n in meas_cols]
    vert = [n.split('_')[2 + meas_uscores] for n in meas_cols]
    # Create dictionary - meas_H = keys, V = values
    hv_dict = {meas + '_' + p:[] for p in set(horiz)}
    for i, pnum in enumerate(horiz):
        hv_dict[meas + '_' + pnum].append(vert[i])

    return hv_dict

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
    
