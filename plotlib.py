"""
Plot functions for checking data from a MojaveCarbon datalogger
"""

import sys
sys.path.append('/home/greg/data/current/ecoflux_tools/')

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import dattool as dtool

def mc_met1_tsplot(df, sitename, colldates):
    """
    Air temp, humidity, PAR, precip
    """
    fig, ax = plt.subplots(4, figsize=(11.5, 8), sharex=True)
    fig.canvas.set_window_title(sitename + ' Met data 1') 
    # Temperature
    ax[0].plot(df.index, df.PTemp_C_Avg, color='k', ls='--', lw=1.25)
    ax[0].plot(df.index, df.AirTC_Avg, color='b', lw=1.25)
    ax[0].plot(df.index, df.AirTC_Min, color='b',ls=':', lw=1.25)
    ax[0].plot(df.index, df.AirTC_Max, color='b',ls=':', lw=1.25)
    ax[0].set_ylabel('Tair ($^\circ$C)')
    ax[0].legend(fontsize=10, loc='best')
    # RH
    ax[1].plot(df.index, df.RH, color='b', lw=1.25)
    ax[1].set_ylabel('RH (%)')
    # PAR
    ax[2].plot(df.index, df.PPFD_Up_umol_Avg, color='r', lw=1.25)
    ax[2].set_ylabel('PPFD ($\mu mol/m^2/sec$)')
    # Precip
    ax[3].plot(df.index, df.Rain_mm_Tot, color='b', marker='.', ms=7, ls='None')
    ax[3].set_ylabel('Precip. (mm)')
    return fig

def mc_met2_tsplot(df, sitename):
    """
    Wind and Pressure
    
    There may be a way to get windbarbs on the plot,
    see this code: https://code.google.com/archive/p/windbarb/downloads
    and the barbs() and quiver() functions in matplotlib
    """
    fig, ax = plt.subplots(3, figsize=(11.5, 8), sharex=True)
    fig.canvas.set_window_title(sitename + ' Met data 2') 
    # Wind speed
    ax[0].plot( df.index, df.WS_ms_Avg, color='g', lw=1.25 )
    ax[0].plot( df.index, df.WS_ms_Max, color='g', ls=':', lw=1.25 )
    #plt.plot( df.index, df.WS_ms_MEAN1_WVT,'--m', lw=1.25 )
    ax[0].set_ylabel('Wind speed (m/s)')
    # Wind direction
    ax[1].plot( df.index, df.WindDir_MEAN1_WVT, color='b', lw=1.25 )
    ax[1].plot( df.index, df.WindDir_MEAN1_WVT + df.WindDir_SD1_WVT,
            color='b', ls=':', lw=1.25 )
    ax[1].plot( df.index, df.WindDir_MEAN1_WVT - df.WindDir_SD1_WVT,
            color='b', ls=':', lw=1.25 )
    ax[1].set_ylabel('Wind Dir ($^\circ$)')
    # Atmospheric pressure
    ax[2].plot( df.index, df.BP_hPa, color='k', lw=1.25 )
    ax[2].set_ylabel('Pressure (hPa)')
    return fig

def mc_power_tsplot(df, sitename):
    """
    Battery/station diagnostics
    
    Plot power diagnostics
    """
    hasppfd = True if 'PPFD_Up_umol_Avg' in df.columns else False
    fig, ax = plt.subplots(3, figsize=(11.5, 8), sharex=True)
    fig.canvas.set_window_title(sitename + ' Power diagnostics')
    if hasppfd:
        # PAR
        ax[0].plot( df.index, df.PPFD_Up_umol_Avg, color='r', lw=1.25 )
    ax[0].set_ylabel('PPFD (umol/m2/sec)')
    # Panel T
    ax[1].plot( df.index, df.PTemp_C_Avg, color='k', lw=1.25 )
    ax[1].set_ylabel('Panel Temp (C)')
    # Battery voltage
    ax[2].plot( df.index, df.BattV, color='r', lw=1.25)
    ax[2].set_ylabel('Battery (V)')
    return fig


def meas_profile_tsplot(df, sitename, var, ylabel, ylimit=None):
    """
    Make a time series plot for sensors in a measurement profile
    """
    # Get measurement dictionary
    measdict = dtool.measurement_h_v_dict(df.columns, var)
    nplots = len(measdict.keys())
    # Set up plot
    fig, ax = plt.subplots(nplots, figsize=(11.5, 8), sharex=True)
    if nplots==1: ax = [ax]
    fig.canvas.set_window_title(sitename + ' ' + var + ' timeseries') 
    # Loop through each profile and depth and plot
    for i, pnum in enumerate(sorted(measdict.keys())):
        for d in measdict[pnum]:
            colname = pnum + '_' + d + '_Avg'
            ax[i].plot( df.index, df[colname], lw=1.25, label=str(d)+'cm' )
        ax[i].legend(loc='upper left', bbox_to_anchor=(0, 1.05),
                ncol=4, fontsize=10)
        if ylimit is not None:
            ax[i].set_ylim(ylimit)
        ax[i].set_title('Profile ' + pnum)
        ax[i].set_ylabel(ylabel)
    return fig

def meas_profile_scatter(df, sitename, var, ylabel, ylimit=[-155,0]):
    """
    Make a scatterplot for sensors in a measurement profile
    """
    # Get measurement dictionary
    measdict = dtool.measurement_h_v_dict(df.columns, var)
    nplots = len(measdict.keys())
    # Set up plot
    fig, ax = plt.subplots(1, nplots, figsize=(7, 5), sharey=True)
    if nplots==1: ax = [ax]
    fig.canvas.set_window_title(sitename + ' ' + var + ' profile') 
    # Loop through each profile and depth and plot againt depth
    for i, pnum in enumerate(sorted(measdict.keys())):
        for d in measdict[pnum]:
            colname = pnum + '_' + d + '_Avg'
            ax[i].plot(df[colname],np.tile(-int(d), [len(df), 1]),
                    marker='o', ls='None', label=str(d)+'cm' )
        ax[i].set_title('Profile ' + pnum)
        ax[i].set_ylim(ylimit)
        ax[i].set_xlabel(ylabel)
        if i==0:
            ax[i].set_ylabel('Depth (cm)')
    return fig

def tsplot_add_colldates(fig, colldates):
    for a in fig.axes:
        ymin, ymax = a.get_ylim()
        a.vlines(colldates, ymin, ymax, linestyles='dotted',lw=0.5)

