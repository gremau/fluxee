import numpy as np
import matplotlib.pyplot as plt
import pdb

def get_adjusted_Da(TsK_layer, Pa, coeff = 1.47e-5):
    """
    Get gas diffusion coefficient in the free air, Da (adjusted for T and P). 
    Default is for CO2 (1.47e-5 from Jones 1992 ???)
    """
    Da_adj = (coeff*((TsK_layer/293.15)**1.75)*(1.013e+5/(Pa)))
    return Da_adj

def get_airfilled_poros(SWC_layer, poros):
    """
    Get air-filled soil porosity given total porosity and soil water content
    of a layer
    """
    return poros - SWC_layer

def soil_diff_moldrup_1999(Da_d1_d2, af_poros, poros, S):
    # Calculation of CO2 DIFFUSIVITY in the soil
    # The CO2 diffusion coefficient in the soil, Ds, is equal to the product of
    # CO2 diffusion coefficient in the free air, Da (adjusted for T and P), 
    # and the gas tortuosity factor in soil. This calculation is based on the 
    # Moldrup 1999 model that adjusts air-filled porosity using soil water.
    
    beta = 2.9 # b Constant (empirically determined -  see Moldrup 99 paper)
    # Soil diffusivity calculation using Moldrup et al. model 1999 (Eqn7)    
    # Calculate the air-filled porosity using mean SWC for depth interval
    Ds_Da_d1_d2 = (poros**2.)*((af_poros/poros)**(beta*S))* Da_d1_d2

    return Ds_Da_d1_d2

def soil_diff_millington_1959(Da_d1_d2, af_poros, poros, S):
    # Calculation of CO2 DIFFUSIVITY in the soil
    # Millington 1959 model
    # The CO2 diffusion coefficient in the soil, Ds, is equal to the product of
    # CO2 diffusion coefficient in the free air, Da (adjusted for T and P), 
    # and the gas tortuosity factor in soil. This calculation is based on the 
    # Moldrup 1999 model that adjusts air-filled porosity using soil water.
    
    Ds_Da_d1_d2 =  (af_poros**(4/3)) * Da_d1_d2

    return Ds_Da_d1_d2

def soil_diff_millington_1961(Da_d1_d2, af_poros, poros, S):
    # Calculation of CO2 DIFFUSIVITY in the soil
    # Millington 1959 model
    # The CO2 diffusion coefficient in the soil, Ds, is equal to the product of
    # CO2 diffusion coefficient in the free air, Da (adjusted for T and P), 
    # and the gas tortuosity factor in soil. This calculation is based on the 
    # Moldrup 1999 model that adjusts air-filled porosity using soil water.
    
    Ds_Da_d1_d2 =  ((af_poros**(10/3))/(poros**2)) * Da_d1_d2

    return Ds_Da_d1_d2

def soil_diff_penman_1940(Da_d1_d2, af_poros, poros, S):
    # Calculation of CO2 DIFFUSIVITY in the soil
    # Penman 1940
    # The CO2 diffusion coefficient in the soil, Ds, is equal to the product of
    # CO2 diffusion coefficient in the free air, Da (adjusted for T and P), 
    # and an empirical relationship to soil porosity.

    Ds_Da_d1_d2 =  (0.66 * af_poros) * Da_d1_d2

    return Ds_Da_d1_d2


def gradient_flux_prod(CO2,Ts,SVWC,P,poros,S,z_vals=[.05, .10, .30],
        Ds_func=soil_diff_moldrup_1999,makeplots=False):
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
        CO2  : (ndarray) soil CO2 concentrations [PPM] at depth X
        Ts   : (ndarray) soil temperature [Celsiuis] at depth X
        SVWC : (ndarray) volumetric soil water content at depth (fraction)X
        P    : (ndarray) atmospheric pressure [Kpa]
        S:   : (float) silt + sand content
        poros: (float) porosity (see Moldrup model and Vargas papers, but
                  this may need to be a list since it changes with depth)
        z_vals: (optional list) Soil depths in m for each depth

    Outputs (ndarrays):
        flux_arr   : surface, d1, d2, d3 CO2 flux  [mumol/m2/s]
        prod_arr   : surface-d1, d1-d2, d2-d3, d1-d3 CO2 production [mumol/m3/s]
        diffity_arr: d1-d2, d1-d3, d2-d3 CO2 diffusivity [?]

    Usage: [F0 Pco2]=fluxprod_3(CO2,Ts,SVWC,P,poros,S,z_vals=[.05, .10, .30])
    '''

    # Check inputs
    n_CO2depth = CO2.shape[1]
    n_TSdepth = Ts.shape[1]
    n_VWCdepth = SVWC.shape[1]
    n_depths = len(z_vals)
    try:
        if all(x == n_CO2depth for x in (n_TSdepth, n_VWCdepth, n_depths)):
            print("Calculating production for " + str(n_CO2depth) + " depths.")
        else:
            raise ValueError("Different # of CO2, TS, VWC, and z_vals given!")
    except ValueError as err:
        print(err.args)
        raise

    # constants and parameters
    R = 8.3144 # Universal gas constant
    # conversions of T and P
    TsK = Ts + 273.15
    Pa = P * 1000  # Transform air pressure from kPa to Pascals (Pa)
    
    #Transform CO2 from ppm(umol mol-1) to mole concentration (umol m-3)
    # Note: Air pressure should be in Pascals
    CO2mc = (CO2*Pa[:,np.newaxis])/(R*TsK)

    # Depth intervals (in original vargas code, but not sure how these are
    # supposed to work) They give the midpoint depths between each CO2 sensor.
    z2 = np.tile(np.nan,(n_CO2depth))
    z2[0]=(z_vals[1]+z_vals[0])/2
    z2[1]=(z_vals[2]+z_vals[1])/2
    z2[2]=(z_vals[2]+z_vals[0])/2

    F_out = np.zeros((CO2mc.shape[0], n_CO2depth))
    Ds_out = np.zeros((CO2mc.shape[0], n_CO2depth))
    P_out = np.zeros((CO2mc.shape[0], n_CO2depth+1))
    for i in range(n_CO2depth):
        if i==(n_CO2depth - 1):
            # Use T and SWC mean of entire soil profile
            Da_di_dj = get_adjusted_Da(TsK.mean(axis=1), Pa)
            afporos_di_dj = poros - SVWC.mean(axis=1)
            Ds_Da = Ds_func(Da_di_dj, afporos_di_dj, poros, S)
            F = Ds_Da*((CO2mc[:,i]-CO2mc[:,0])/(z_vals[i]-z_vals[0]))
        else:
            Da_di_dj = get_adjusted_Da((TsK[:,i] + TsK[:,i+1])/2, Pa)
            afporos_di_dj = poros - (SVWC[:,i] + SVWC[:,i+1])/2
            Ds_Da = Ds_func(Da_di_dj, afporos_di_dj, poros, S)
            F = Ds_Da*((CO2mc[:,i+1]-CO2mc[:,i])/(z_vals[i+1]-z_vals[i]))

        F_out[:,i] = F
        Ds_out[:,i] = Ds_Da

    #=== Calculation of the Flux at soil surface ===
    # Do this only if there is no surface measurement
    if z_vals[0] > 0:
        alen = F_out.shape[0]
        F0 = np.zeros((alen, 1))
        for i in range(alen): #Calculation of means for Node 1
            try:
                if n_CO2depth > 3:
                    fluxes = np.array([F_out[i,0], F_out[i,1], F_out[i,2]])
                    z_ind = z2[0:3]
                else:
                    # Reorder to use the middle (full profile) flux value in
                    # the fit
                    fluxes = np.array([F_out[i,0], F_out[i,2], F_out[i,1]])
                    z_ind = np.array([z2[0], z2[2], z2[1]])
                idx = np.isfinite(z_ind) & np.isfinite(fluxes)
                p = np.polyfit(z_ind[idx], fluxes[idx], 1)
                F0[i]=p[1]
            except:
                F0[i]=np.nan
        F_out = np.concatenate((F0, F_out), axis=1)

    #%%%%%Calculation of the soil CO2 production%%%%
    # Still not sure I am calculating production right
    # not sure about how these depths are being used
    for i in range(n_CO2depth+1):
        if i==0:
            P = (F_out[:,0] - F_out[:,1])/(z2[0]-0)
        elif i==n_CO2depth:
            P = (F_out[:,1] - F_out[:,i])/(z2[i-1]-z2[0])
        else:
            P = (F_out[:,i] - F_out[:,i+1])/(z2[i]-z2[i-1])
        P_out[:,i] = P

    if makeplots:
        f, axarr = plt.subplots(3, sharex=True)
        axarr[0].plot(F_out)
        axarr[0].set_ylabel('CO2 flux')
        axarr[0].legend(['surface', 'd1', 'd2', 'dT'], ncol=4)
        axarr[1].plot(P_out)
        axarr[1].set_ylabel('CO2 production')
        axarr[1].legend(['surface-d1', 'd1-d2', 'd1-d3'], ncol=4)
        axarr[2].plot(Ds_out, label=['d1', 'd2', 'dT'])
        axarr[2].set_ylabel('Diffusivity')
        axarr[2].legend(['d1', 'd2', 'dT'], ncol=4)
    return F_out, P_out, Ds_out
