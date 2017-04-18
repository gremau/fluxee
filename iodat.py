""" 
Functions for loading various environmental monitoring data files

Greg Maurer
"""

import numpy as np
import datetime as dt
import pandas as pd
import os

now = dt.datetime.now()

def calculate_freq(idx):
    cfreq = (idx.max()-idx.min())/(len(idx)-1)
    cfreq = cfreq.seconds/60
    print("Calculated frequency is " + "{:5.3f}".format(cfreq) + " minutes")
    print("Rounding to " + str(round(cfreq)) + 'min')
    return str(round(cfreq)) + "min"

def load_toa5( fpathname, reindex=False ) :
    """
    Load a specified TOA5 datalogger file (a Campbell standard output format)
    and return a pandas DataFrame object. DataFrame has a datetime index and
    user is warned if any measurement periods appear to be missing. Dataframe
    can be reindexed to fill in missing periods with NAN.

    Args:
        fpathname (str) : path and filename of desired AF file
        efreq (str)     : expected frequency of data file, used to reindex
    Return:
        parsed_df   : pandas DataFrame    
    """

    print('Parsing ' + fpathname)

    # Parse using Campbell timestamp
    parsed_df = pd.read_csv(fpathname, skiprows=( 0,2,3 ), header=0,
            parse_dates = { 'Date': [0]}, index_col='Date',
            na_values=['NaN', 'NAN', 'INF', '-INF'])
    
    cfreq = calculate_freq(parsed_df.index)
    # Create an index that includes every period between the first and
    # last datetimes in the file
    startd = parsed_df.index.min()
    endd = parsed_df.index.max()
    fullidx = pd.date_range( startd, endd, freq=cfreq)
    # Warn if observations are missing
    if len( parsed_df.index ) < len( fullidx ):
        print("WARNING: index frequency is less than expected!")
        print("Reindexing will introduce NaN values")
    elif len( parsed_df.index ) > len( fullidx ):
        print("WARNING: index frequency is greater than expected!")
        print("Reindexing may remove valid values")
    if reindex:
        print("Reindexing dataframe...")
        parsed_df_ret =  parsed_df.reindex(fullidx)
    else:
        parsed_df_ret = parsed_df
    return parsed_df_ret

def load_century_lis( fpathname ) :
    """
    Load a specified century output file and return a pandas DataFrame object.
    DataFrame has a datetime index and has been reindexed to include all
    one-month periods. 

    Args:
        fpathname (str) : path and filename of desired century (.lis) file
    Return:
        parsed_df   : pandas DataFrame    
    """

    print('Parsing ' + fpathname)

    # Parse using Campbell timestamp
    parsed_df = pd.read_csv(fpathname, delim_whitespace=True, skiprows=( 1, ),
            header=0, parse_dates = { 'Date': [0]}, index_col='Date',
            na_values=['NaN', 'NAN', 'INF', '-INF'])
    # This is tricky around year zero because Century has a wierd gap
#    year = np.floor(parsed_df.index))
#    months = np.round((parsed_df.index - year)*12)
#    index = dt.datetime(
#    # Warn if observations are missing
#    if len( parsed_df.index ) < len( full_idx ):
#        print( "WARNING: some observations may be missing!" )
        
    return parsed_df

def site_datafile_concat(sitename, datapath, func=load_toa5) :
    """
    Load a list of datalogger files, append them, and then return a pandas
    DataFrame object. Also returns a list of collection dates. This function
    expects to find a directory full of files with the format 
    "prefix_<sitename>_<Y>_<m>_<d>_<H>_<M>_optionalsuffix.dat". For example:

    MNPclimoseq_Creosote_2017_03_17_11_55_00.dat

    Only files matching the 'sitename' argument are concatenated into the
    resulting file.
    
    Note that files don't load in chronological order, so the resulting 
    dataframe is reindexed based on the min/max dates in the indices. This 
    will fill in any missing values with NAN.
    
    Warns if observations are missing. Fails if indices of concatenated files
    have different frequencies.
    
    Args:
        sitename    : Site name
        datapath    : Path to directory of data files
        func        : function used to load each file
    Return:
        sitedf      : pandas DataFrame containing concatenated raw data
                      from one site
        collect_dt  : list of datetime objects parsed from the filenames 
                      (data collection date)
    """
            
    # Get a list of filenames in the raw directory for the site
    files = os.listdir( datapath )
    # Select desired files from file_list (by site)
    site_files = [ f for f in files if sitename in f ]

    # Get collection date for each file. The collection date is in the 
    # filename with fields delimited by '_'
    collect_dt = []
    for i in site_files:
        tokens = i.split('_')
        collect_dt.append(dt.datetime.strptime('-'.join(tokens[-6:-1]),
            '%Y-%m-%d-%H-%M'))

    # Initialize DataFrame
    sitedf = pd.DataFrame()
    # Loop through each year and fill the dataframe
    for j in site_files:
        # Call load_toa5_file
        filedf = func(datapath + j)
        # And append to site_df, 'verify_integrity' warns if there are
        # duplicate indices
        sitedf = sitedf.append(filedf, verify_integrity=True)
    # Calculate frequency
    cfreq = calculate_freq(sitedf.index)
    # Create index spanning all days from min to max date
    fullidx = pd.date_range(sitedf.index.min(), sitedf.index.max(),
            freq = cfreq)
    # Warn if observations are missing
    if len( sitedf.index ) < len( fullidx ):
        print("WARNING: index frequency is less than expected!")
        print("Reindexing will introduce NaN values")
    elif len( sitedf.index ) > len( fullidx ):
        print("WARNING: index frequency is greater than expected!")
        print("Reindexing may remove valid values")
    # Now reindex the dataframe
    print("Reindexing dataframe...")
    sitedf = sitedf.reindex( fullidx )
    return sitedf, collect_dt
