import numpy as np
import matplotlib.pyplot as plt
import pdb

def fluxprod_3(CO2,Ts,SVWC,P,poros,S,z_vals=[.05, .10, .30], makeplots=False):
    '''
    This code is adapted from MATLAB code written by: Rodrigo Vargas,
    University of Delaware (rvargas@udel.edu)

    When using this function please cite the following papers:

    Vargas R. and M. F. Allen, 2008. Dynamics of fine root, fungal rhizomorphs,
    and soil respiration in a mixed temperate forest: Integrating sensors and
    observations. Vadose Zone Journal 7:1055-1064.

    Vargas et al. 2010. Looking deeper into the soil: biophysical controls and
    seasonal lags of soil CO2 production and efflux. Ecol. Appl. 20:1569-1582. 

    This function computes the soil surface CO2 flux using the Gradient method
    
    Inputs :
        CO2  : (Nx3 ndarray) soil CO2 concentrations [PPM] at depth X
        Ts   : (Nx3 ndarray) soil temperature [Celsiuis] at depth X
        SVWC : (Nx3 ndarray) volumetric soil water content at depth X
        P    : (ndarray) atmospheric pressure [Kpa]
        S:   : (float) silt + sand content
        poros: (float) porosity (see Moldrup model and Vargas papers, but
                  this may need to be an array since it changes with VWC)
        z_vals: (optional list) Soil depths in m for d1, d2, and d3

    Outputs (ndarrays):
        flux_arr   : surface, d1, d2, d3 CO2 flux  [mumol/m2/s]
        prod_arr   : surface-d1, d1-d2, d2-d3, d1-d3 CO2 production [mumol/m3/s]
        diffity_arr: d1-d2, d1-d3, d2-d3 CO2 diffusivity [?]

    Usage: [F0 Pco2]=fluxprod_3(CO2,Ts,SVWC,P,poros,S,z_vals=[.05, .10, .30])
    '''
    # constants and parameters
    beta = 2.9 # Constant for Moldrup model
    R = 8.3144 # Universal gas constant
    # Distances
    z2 = np.array([np.nan, np.nan, np.nan])
    z2[0]=(z_vals[1]+z_vals[0])/2
    z2[1]=(z_vals[2]+z_vals[0])/2
    z2[2]=(z_vals[2]+z_vals[1])/2
    
    # Not sure why this is here
    #z3[0]=(z2[1]+z2[0])/2
    #z3[1]=(z2[2]+z2[0])/2
    #z3[2]=(z2[2]+z2[1])/2
    #Transform T soil in Kelvin
    TsK = Ts + 273.15
    Pa = P * 1000  # Transform air pressure from kPa to Pascals (Pa)
    
    #Transform from ppm(umol mol-1) to mole concentration (umol m-3)
    # Note: Air pressure should be in Pascals
    CO2mc = (CO2*Pa[:,np.newaxis])/(R*TsK)
    #C_d1 = (CO2[:,0]*Pa) /(R*TsK[:,0])
    #C_d2 = (CO2[:,1]*Pa) /(R*TsK[:,1])
    #C_d3 = (CO2[:,2]*Pa)/(R*TsK[:,2])
    # /////////////////////////////////////////////////////////////
    # calculation of CO2 DIFFUSION in the soil
    # The CO2 diffusion coefficient in the soil, Ds, is equal to the product of
    # CO2 diffusion coefficient in the free air, Da, and the gas tortuosity 
    # factor .... Temperature(oK) Pressure (kPa)

    #Da for layers Da_d1_d2 (2-8 cm) Da16 (8-16 cm), Da2_16 (2-16 cm)
    Da_d1_d2 = (1.47e-5*((((TsK[:,0] + TsK[:,1])/2)/293.15)**1.75)*
            (1.013e+5/(Pa)))
    Da_d2_d3 = (1.47e-5*((((TsK[:,2] + TsK[:,1])/2)/293.15)**1.75)*
            (1.013e+5/(Pa)))
    Da_d1_d3 = (1.47e-5*((((TsK[:,2] + TsK[:,0])/2)/293.15)**1.75)*
            (1.013e+5/(Pa)))

    # Diffusivity calculation using Moldrup et al. model 1999    
    # Calculate the air-filled porosity (using mean SWC on sensor range)
    Ds_Da_d1_d2 = (poros**2)*(((poros-SVWC[:,(0,1)].mean(axis=1))/poros)
            **(beta*S))* Da_d1_d2
    Ds_Da_d2_d3 = (poros**2)*(((poros-SVWC[:,(1,2)].mean(axis=1))/poros)
            **(beta*S))* Da_d2_d3
    Ds_Da_d1_d3 = (poros**2)*(((poros-SVWC[:,1])/poros)
            **(beta*S))* Da_d1_d3
    
    #=== Calculation of the Flux at each layer ===
    
    # calculation of CO2 flux from Fd2_d1,  Fd3_d2, Fd3_d1......
    # Note: depth should be in meters
    # z1 = 0.06..... or distance in between 0.08 and 0.02 m
    # z2 = 0.09..... or distance in between 0.16 and 0.02 m
    # z2 = 0.12..... or distance in between 0.16 and 0.08 m
    F1 = Ds_Da_d1_d2*((CO2mc[:,1]-CO2mc[:,0])/(z_vals[1]-z_vals[0]));
    F2 = Ds_Da_d1_d3*((CO2mc[:,2]-CO2mc[:,0])/(z_vals[2]-z_vals[0]));
    F3 = Ds_Da_d2_d3*((CO2mc[:,2]-CO2mc[:,1])/(z_vals[2]-z_vals[1]));
    #=== Calculation of the Flux at soil surface ===
    F0 = np.zeros(F1.shape)
    for i in range(0, len(F1)): #Calculation of means for Node 1
        #length(F1)-i
        try:
            fluxes = np.array([F1[i], F2[i], F3[i]])
            idx = np.isfinite(z2) & np.isfinite(fluxes)
            p = np.polyfit(z2[idx], fluxes[idx], 1)
            F0[i]=p[1]
        except:
            F0[i]=np.nan

    #%%%%%Calculation of the soil CO2 production%%%%
    # p1 = 0.07..... or distance in between 0.05 and 0.09 m
    # p2 = 0.085..... or distance in between 0.05 and 0.12 m
    # p3 = 0.105..... or distance in between 0.09 and 0.12 m
    PR_0=(F0 - F1)/(z2[0]-0);
    PR_1=(F1 - F2)/(z2[1]-z2[0]);
    PR_2=(F2 - F3)/(z2[2]-z2[1]);
    PR_T=(F1 - F3)/(z2[2]-z2[0]);

    flux_arr = np.array([F0,F1,F2,F3]).T
    prod_arr = np.array([PR_0,PR_1,PR_2,PR_T]).T
    diffity_arr = np.array([Ds_Da_d1_d2,Ds_Da_d1_d3,Ds_Da_d2_d3]).T

    if makeplots:
        f, axarr = plt.subplots(3, sharex=True)
        axarr[0].plot(flux_arr)
        axarr[0].set_ylabel('CO2 flux')
        axarr[0].legend(['surface', 'd1', 'd2', 'd3'], ncol=4)
        axarr[1].plot(prod_arr)
        axarr[1].set_ylabel('CO2 production')
        axarr[1].legend(['surface-d1', 'd1-d2', 'd2-d3', 'd1-d3'], ncol=4)
        axarr[2].plot(diffity_arr, label=['d1', 'd2', 'd3'])
        axarr[2].set_ylabel('Diffusivity')
        axarr[2].legend(['d1', 'd2', 'd3'], ncol=4)
    return flux_arr, prod_arr, diffity_arr


