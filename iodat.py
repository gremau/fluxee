""" 
Functions for loading various environmental monitoring data files

Greg Maurer
"""

import numpy as np
import datetime as dt
import pandas as pd
import subprocess as sp
import yaml
import os
import pdb

def get_file_collection(sitename, datpath, ext='.dat'):
    """
    Read a list of files from a data directory, match against desired site,
    and return the list and the collection dates of each file. This function
    expects to find a directory full of files with the format 
    "prefix_<sitename>_<Y>_<m>_<d>_<H>_<M>_optionalsuffix.dat". For example:

    MNPclimoseq_Creosote_2017_03_17_11_55_00.dat

    Only returns files matching the sitename and extension (.dat by 
    default) strings.
    """
    # Get a list of filenames in provided data directory
    files = os.listdir(datpath)
    # Select desired files from the list (by site)
    site_files = [ f for f in files if sitename in f and ext in f ]
    # Get collection date for each file. The collection date is in the 
    # filename with fields delimited by '_'
    collect_dt = []
    for i in site_files:
        tokens = i.split('_')
        collect_dt.append(dt.datetime.strptime('-'.join(tokens[-6:-1]),
            '%Y-%m-%d-%H-%M'))
    return site_files, collect_dt


def read_yaml_conf(sitename, yamltype, confdir='ecoflux_config'):
    """
    Read a specified YAML configuration file from a given site's ecoflux
    configuration directory. Checks the YAML file meta dictionary to ensure
    configuration is for the correct site and type.

    Args:
        sitename (string): site identifier
        yamltype (string): identifies the YAML filename and configuration type
        confdir (string): directory to look for YAML configuration files

    Returns:
        yamlf (dict): Returns the "items" dictionary of configuration values
                      from the YAML file

    Raises:
        ValueError: Raises ValueError if the sitename or yamltype parameters
                    do not match those found in the yaml file
    """
    yamlfile = os.path.join(confdir, sitename, yamltype + '.yaml')
    stream = open(yamlfile, 'r')
    yamlf = yaml.load(stream)
    ysite = yamlf['meta']['site']==sitename
    ytype = yamlf['meta']['conftype']==yamltype
    if not(ysite) or not(ytype):
        raise ValueError('YAML file site/type mismatch.')
    else:
        return yamlf['items']

def calculate_freq(idx):
    """
    """
    cfreq = (idx.max()-idx.min())/(len(idx)-1)
    cfreq = cfreq.seconds/60
    print("Calculated frequency is " + "{:5.3f}".format(cfreq) + " minutes")
    print("Rounding to " + str(round(cfreq)) + 'min')
    return str(round(cfreq)) + "min"

def load_toa5(fpathname, reindex=False) :
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

def site_datafile_concat(sitename, datpath, iofunc=load_toa5):
    """
    Load a list of datalogger files, append them, and then return a pandas
    DataFrame object. Also returns a list of collection dates. 

    Note that files don't load in chronological order, so the resulting 
    dataframe is reindexed based on the min/max dates in the indices. This 
    will fill in any missing values with NAN.
    
    Warns if observations are missing. Fails if indices of concatenated files
    have different frequencies.
    
    Args:
        sitename    : Site name
        datpath     : Path to directory of data files
        iofunc      : function used to load each file
    Return:
        sitedf      : pandas DataFrame containing concatenated raw data
                      from one site
        collect_dt  : list of datetime objects parsed from the filenames 
                      (data collection date)
    """
            
    # Get list of files and thier collection dates from the provided directory
    files, collect_dt = get_file_collection(sitename, datpath)
    # Initialize DataFrame
    sitedf = pd.DataFrame()
    # Loop through each year and fill the dataframe
    for i in files:
        # Call load_toa5_file
        filedf = iofunc(datapath + i)
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


def rename_raw_variables(sitename, rawpath, rnpath, confdir='ecoflux_config'):
    """
    Rename raw datalogger variables according to YAML configuration file

    Args:
        sitename (string): site identifier
        rawpath (string): raw datalogger file path
        rnpath (string): path to write renamed data files
        confdir (string): directory to look for YAML configuration files

    Returns:
        Does not return anything but writes new files in rnpath
    """
    # Get var_rename configuration file for site
    yamlf = read_yaml_conf(sitename, 'var_rename', confdir=confdir)
    # Get list of files and their collection dates from the raw directory
    files, collected = get_file_collection(sitename, rawpath)
    # Loop through each rename event and change the file
    for i, filename in enumerate(files):
        findvars, newvars = ([],[])
        # For each rename event, add variable changes to findvar/repvar
        # if the changes occurred after the file collection date
        for j, key in enumerate(sorted(yamlf.keys())):
            rn = yamlf[key] # Get rename event j
            if rn["first_changed_dt"] > collected[i]:
                findvars = findvars + rn["from"]
                newvars = newvars + rn["to"]

        # Now read in file, replace target strings, write a new file
        filepath = os.path.join(rawpath, filename)
        with open(filepath, 'r') as file:
            filedata = file.read()
        for k, findvar in enumerate(findvars):
            newvar = newvars[k]
            filedata = filedata.replace(findvar, newvar)
        rn_filepath = os.path.join(rnpath, filename)
        with open(rn_filepath, 'w') as file:
            file.write(filedata)


def ecoflux_out(df,foutpath, sitename='Undefined', cscript='Undefined'):
    """
    """
    git_sha = sp.check_output(
            ['git', 'rev-parse', 'HEAD']).decode('ascii').strip()
    meta_data = pd.Series([('location: {0}'.format(sitename)),
        ('date generated: {0}'.format(str(dt.datetime.now()))),
        ('writer: ecoflux.iodat.ecoflux_out'),
        ('git HEAD SHA: {0}'.format(git_sha)),
        ('calling script: {0}'.format(cscript)),
        ('--------')])
    with open(foutpath, 'w') as fout:
        fout.write('---file metadata---\n')
        meta_data.to_csv(fout, index=False)
        df.to_csv(fout, na_rep='NA')
