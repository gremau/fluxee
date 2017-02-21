import pdb
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
