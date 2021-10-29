import numpy as np
import matplotlib.pyplot as plt
from IPython.core.debugger import set_trace

def get_adjusted_Da(TsK_layer, Pa, coeff=1.47e-5):
    """
    Get gas diffusion coefficient in the free air in m^2 s^-1, Da (adjusted
    for T and P). Default is for CO2 (1.47e-5 from Jones 1992 ???). But also
    see Massman 1998 that recommends 1.38e-5
    """
    return coeff*((TsK_layer/293.15)**1.75)*(1.013e+5/(Pa))

def get_airfilled_poros(SWC_layer, poros):
    """
    Get air-filled soil porosity given total porosity and soil water content
    of a layer
    """
    return poros - SWC_layer

def soil_diff_moldrup_1999(Da_d1_d2, af_poros, poros, S):
    """
    Calculation of gas DIFFUSIVITY in the soil based on the Moldrup 1999 model.
    Requires total soil porosity, air-filled porosity value, sand+silt
    percentange and a free-air gas diffusivity value (which may be T and P
    corrected)
    """
    beta = 2.9 # b Constant (empirically determined)
    # From Moldrup et al. model 1999 (Eqn7)    
    Ds_Da_d1_d2 = (poros**2.)*((af_poros/poros)**(beta*S))* Da_d1_d2

    return Ds_Da_d1_d2

def soil_diff_millington_1959(Da_d1_d2, af_poros, poros=None, S=None):
    """
    Calculation of gas DIFFUSIVITY in the soil based on the Millington 1959
    model. Requires an air-filled porosity value and a free-air gas diffusivity
    value (which may be T and P corrected)
    """
    Ds_Da_d1_d2 =  (af_poros**(4/3)) * Da_d1_d2

    return Ds_Da_d1_d2

def soil_diff_millington_1961(Da_d1_d2, af_poros, poros, S=None):
    """
    Calculation of gas DIFFUSIVITY in the soil based on the Millington and
    Quirk 1961 model. Requires total soil porosity, air-filled porosity, and a
    free-air gas diffusivity value (which may be T and P corrected)
    """
    Ds_Da_d1_d2 =  ((af_poros**(10/3))/(poros**2)) * Da_d1_d2

    return Ds_Da_d1_d2

def soil_diff_penman_1940(Da_d1_d2, af_poros, poros, S):
    """
    Calculation of gas DIFFUSIVITY in the soil based on the Penman 1940 model.
    Requires an air-filled porosity value and a free-air gas diffusivity value
    (which may be T and P corrected)
    """
    return (0.66 * af_poros) * Da_d1_d2

def gradient_flux_prod(CO2,Ts,SVWC,P,poros,S,z_vals=[.05, .10, .30],
        Ds_func=soil_diff_moldrup_1999,adjust_Da=True,makeplots=False):
    '''
    This function computes the soil surface CO2 flux using the gradient method.
    To some degree it is adapted from the method (and some MATLAB code) used in
    these papers:

    Vargas R. and M. F. Allen, 2008. Dynamics of fine root, fungal rhizomorphs,
    and soil respiration in a mixed temperate forest: Integrating sensors and
    observations. Vadose Zone Journal 7:1055-1064.

    Vargas et al. 2010. Looking deeper into the soil: biophysical controls and
    seasonal lags of soil CO2 production and efflux. Ecol. Appl. 20:1569-1582. 

    Args :
        CO2      : (ndarray) soil CO2 concentrations [PPM] at depth X
        Ts       : (ndarray) soil temperature [Celsius] at depth X
        SVWC     : (ndarray) volumetric soil water content at depth X (fraction)
        P        : (ndarray) atmospheric pressure in hPa, gets converted
        poros    : (list) porosity for each depth interval) 
        S:       : (list) silt + sand content for each depth interval
        z_vals   : (optional list) Soil depths in m for each depth
        Ds_func  : function for calculating gas diffusivity
        adjust_Da: (bool) flag for to allow T and P correction of free air gas
                   diffusivity

    Returns:
        Three ndarrays: 
        F_out   : surface, d1, d2, d3 CO2 flux  [mumol/m2/s]
        P_out   : surface-d1, d1-d2, d2-d3, d1-d3 CO2 production [mumol/m3/s]
        Ds_out  : d1-d2, d1-d3, d2-d3 CO2 diffusivity [?]
    '''
    import warnings
    warnings.warn("Warning: This function is deprecated!!! See profile_flux()")

    # Check inputs
    n_CO2depth = CO2.shape[1]
    n_TSdepth = Ts.shape[1]
    n_VWCdepth = SVWC.shape[1]
    n_depths = len(z_vals)
    n_S = len(S)
    n_poros = len(poros)
    try:
        if all(x == n_CO2depth for x in (n_TSdepth, n_VWCdepth, n_depths,
            n_S, n_poros)):
            print("Calculating production for " + str(n_CO2depth) + " depths.")
            surfCO2 = None
        elif all(x == n_CO2depth-1 for x in (n_TSdepth, n_VWCdepth, n_depths,
            n_S, n_poros)):
            surfCO2 = CO2[:,0]
            CO2 = CO2[:,1::]
            n_CO2depth = CO2.shape[1]
        else:
            raise ValueError("Different # of CO2, TS, VWC, S, poros, "
                    "and z_vals given!")
    except ValueError as err:
        print(err.args)
        raise
    # constants and parameters
    R = 8.3144 # Universal gas constant
    # conversions of T and P
    TsK = Ts + 273.15
    Pa = P * 100  # Transform air pressure from hPa to Pascals (Pa)
    
    #Transform CO2 from ppm(umol mol-1) to mole concentration (umol m-3)
    # Note: Air pressure should be in Pascals
    CO2mc = (CO2*Pa[:,np.newaxis])/(R*TsK)
    if surfCO2 is not None:
        surfCO2 = (surfCO2*Pa)/(R*TsK[:,0])

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
        # If we are at the last depth calculate flux for entire profile
        if i==(n_CO2depth - 1):
            # Use T, SWC, poros and S mean
            S_i = sum(S)/len(S)
            poros_i = sum(poros)/len(poros)
            if adjust_Da:
                Da_di_dj = get_adjusted_Da(TsK.mean(axis=1), Pa)
            else:
                Da_di_dj = 1.47e-5
            afporos_di_dj = poros_i - SVWC.mean(axis=1)
            Ds_Da = Ds_func(Da_di_dj, afporos_di_dj, poros_i, S_i)
            F = Ds_Da*((CO2mc[:,i]-CO2mc[:,0])/(z_vals[i]-z_vals[0]))
        # Calculate porosities and flux between sensor i and the one below
        else:
            S_i = (S[i] + S[i+1])/2
            poros_i = (poros[i] + poros[i+1])/2
            if adjust_Da:
                Da_di_dj = get_adjusted_Da((TsK[:,i] + TsK[:,i+1])/2, Pa)
            else:
                Da_di_dj = 1.47e-5
            afporos_di_dj = poros_i - (SVWC[:,i] + SVWC[:,i+1])/2
            Ds_Da = Ds_func(Da_di_dj, afporos_di_dj, poros_i, S_i)
            F = Ds_Da*((CO2mc[:,i+1]-CO2mc[:,i])/(z_vals[i+1]-z_vals[i]))

        F_out[:,i] = F
        Ds_out[:,i] = Ds_Da

    #=== Calculation of the Flux at soil surface ===
    # If there is no surface measurement fit a curve to the fluxes below and
    # use the intercept value as the surface flux
    alen = F_out.shape[0]
    F0 = np.zeros((alen, 1))
    if surfCO2 is None:
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
    # Otherwise calculate using method above with surface T/VWC/poros/S values
    else:
        if adjust_Da:
            Da_di_dj = get_adjusted_Da(TsK[:,0], Pa)
        else:
            Da_di_dj = 1.47e-5
        afporos_di_dj = poros[0] - SVWC[:,0]
        Ds_Da = Ds_func(Da_di_dj, afporos_di_dj, poros[0], S[0])
        F0[:,0] = Ds_Da*((CO2mc[:,0]-surfCO2)/(z_vals[0]))

    F_out = np.concatenate((F0, F_out), axis=1)
    #%%%%%Calculation of the soil CO2 production%%%%
    # According to Davidson and Trumbore 1995 and Davidson et. al 2006 (GCB)
    # CO2 production of a layer can be calculated by subtracting the flux at 
    # one depth from the flux of depth below (Fi - Fi+1). Since we have
    # Still not sure I am calculating production right
    # not sure about how these depths are being used
    # Seems to give a production value for the between sensor
    # depth in micromoles per m^3
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


def gradient_flux_layer_old(gas_mf,Ts,SVWC,P,poros,S,Dcoeff=1.47e-5,
        z_vals=[.05, .10, .30],Ds_func=soil_diff_moldrup_1999,
        adjust_Da=True,makeplots=False):
    '''
    This function computes the soil gas flux using the gradient method,
    i.e Fick's first law, with soil gas diffusivity computed using one of
    several available methods.
    
    It is a refined version of gradient_flux_prod generalized for any gas
    diffusing through porous media (soil). This process is most similar to the
    methods used by Tang et al 2005 and Vargas et al 2008

    Returns a surface flux estimate and a layer-based flux profile. The
    surface estimate comes from either a linear interpolation from
    below-soil fluxes, or, if a boundary condition is provided, the flux
    between the topmost sensor to the boundary. Within-soil fluxes are
    calculated from the gradient between belowground sensors, giving a result
    that is the average flux between measurement depth i and i+1.

    For an integrated profile approach see the inverse model method (though
    results from both should be comparable).
    
    Args :
        gas_mf   : (ndarray) soil gas concentrations [PPM] at depth X
        Ts       : (ndarray) soil temperature [Celsius] at depth X
        SVWC     : (ndarray) volumetric soil water content at depth X (fraction)
        P        : (ndarray) atmospheric pressure in hPa, gets converted
        poros    : (list) porosity for each depth interval) 
        S:       : (list) silt + sand content for each depth interval
        Dcoeff   : (numeric) free-air diffusion coefficient for gas species
        z_vals   : (optional list) Soil depths in m for each depth
        Ds_func  : function for calculating gas diffusivity
        adjust_Da: (bool) flag for to allow T and P correction of free air gas
                   diffusivity

    Returns:
        Three ndarrays: 
        F_out   : flux [umol/m2/s] at depth intervals d0d1, d1d2, d2d3, d3d4...
        Ds_out  : gas diffusivity values [m2/s2] at d0d1, d1d2, d2d3, d3d4...
    '''

    # Check inputs
    n_gasdepth = gas_mf.shape[1]
    n_TSdepth = Ts.shape[1]
    n_VWCdepth = SVWC.shape[1]
    n_depths = len(z_vals)
    n_S = len(S)
    n_poros = len(poros)
    try:
        if all(x == n_gasdepth for x in (n_TSdepth, n_VWCdepth, n_depths,
            n_S, n_poros)):
            print("Calculating production for " + str(n_gasdepth) + " depths.")
            surfgas = None
        elif all(x == n_gasdepth-1 for x in (n_TSdepth, n_VWCdepth, n_depths,
            n_S, n_poros)):
            surfgas = gas_mf[:,0]
            gas_mf = gas_mf[:,1::]
            n_gasdepth = gas_mf.shape[1]
        else:
            raise ValueError("Different # of gas, TS, VWC, S, poros, "
                    "and z_vals given!")
    except ValueError as err:
        print(err.args)
        raise
    # constants and parameters
    R = 8.3144 # Universal gas constant
    # conversions of T and P
    TsK = Ts + 273.15
    Pa = P * 100  # Transform air pressure from hPa to Pascals (Pa)
    
    #Transform gas from ppm(umol mol-1) to mole concentration (umol m-3)
    # Note: Air pressure should be in Pascals
    gas_mc = (gas_mf*Pa[:,np.newaxis])/(R*TsK)
    if surfgas is not None:
        surfgas = (surfgas*Pa)/(R*TsK[:,0])
    
    # Depth intervals (in original vargas code, but not sure how these are
    # supposed to work) They give the midpoint depths between each CO2 sensor.
    z2 = np.tile(np.nan,(n_gasdepth))
    z2[0]=(z_vals[1]+z_vals[0])/2
    z2[1]=(z_vals[2]+z_vals[1])/2
    if len(z_vals) > 3:
        z2[2]=(z_vals[3]+z_vals[2])/2
    else:
        z2[2]=(z_vals[2]+z_vals[0])/2

    F_out = np.zeros((gas_mc.shape[0], n_gasdepth))
    Ds_out = np.zeros((gas_mc.shape[0], n_gasdepth))
    for i in range(n_gasdepth):
        # If we are at the last depth calculate flux for entire profile
        if i==(n_gasdepth - 1):
            # Use T, SWC, poros and S mean
            S_i = sum(S)/len(S)
            poros_i = sum(poros)/len(poros)
            if adjust_Da:
                Da_di_dj = get_adjusted_Da(TsK.mean(axis=1), Pa, coeff=Dcoeff)
            else:
                Da_di_dj = Dcoeff
            afporos_di_dj = poros_i - SVWC.mean(axis=1)
            Ds_Da = Ds_func(Da_di_dj, afporos_di_dj, poros_i, S_i)
            if surfgas is None:
                # If no boundary condition flux calculated with deepest
                # and shallowest values
                F = Ds_Da*((gas_mc[:,i]-gas_mc[:,0])/(z_vals[i]-z_vals[0]))
            else:
                # Otherwise use the surface and lowest value
                F = Ds_Da*((gas_mc[:,i]-surfgas)/(z_vals[i]))

        # Calculate porosities and flux between sensor i and i+1
        else:
            S_i = (S[i] + S[i+1])/2
            poros_i = (poros[i] + poros[i+1])/2
            if adjust_Da:
                Da_di_dj = get_adjusted_Da((TsK[:,i] + TsK[:,i+1])/2, Pa,
                        coeff=Dcoeff)
            else:
                Da_di_dj = Dcoeff
            afporos_di_dj = poros_i - (SVWC[:,i] + SVWC[:,i+1])/2
            Ds_Da = Ds_func(Da_di_dj, afporos_di_dj, poros_i, S_i)
            F = Ds_Da*((gas_mc[:,i+1]-gas_mc[:,i])/(z_vals[i+1]-z_vals[i]))

        F_out[:,i] = F
        Ds_out[:,i] = Ds_Da

    #=== Calculation of the Flux at soil surface ===
    # If there is no surface measurement fit a curve to the fluxes below and
    # use the intercept value as the surface flux
    alen = F_out.shape[0]
    F0 = np.zeros((alen, 1))
    DsDa0 = np.zeros((alen, 1))
    if surfgas is None:
        for i in range(alen): #Calculation of means for Node 1
            try:
                if n_gasdepth > 3:
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
    # Otherwise calculate using method above with surface T/VWC/poros/S values
    else:
        if adjust_Da:
            Da_di_dj = get_adjusted_Da(TsK[:,0], Pa, coeff=Dcoeff)
        else:
            Da_di_dj = Dcoeff
        afporos_di_dj = poros[0] - SVWC[:,0]
        Ds_Da = Ds_func(Da_di_dj, afporos_di_dj, poros[0], S[0])
        DsDa0[:,0] = Ds_Da
        # I have occasionally added a boundary layer to surface,
        # which lowers the flux
        F0[:,0] = Ds_Da*((gas_mc[:,0]-surfgas)/(z_vals[0]))

    F_out = np.concatenate((F0, F_out), axis=1)
    Ds_out = np.concatenate((DsDa0, Ds_out), axis=1)
    
    if makeplots:
        f, axarr = plt.subplots(2, sharex=True)
        axarr[0].plot(F_out)
        axarr[0].set_ylabel('Gas flux')
        axarr[0].legend(['surface', 'd1', 'd2', 'dT'], ncol=4)
        axarr[1].plot(Ds_out, label=['d1', 'd2', 'dT'])
        axarr[1].set_ylabel('Diffusivity')
        axarr[1].legend(['d1', 'd2', 'dT'], ncol=4)
    return F_out, Ds_out

def diff_profile_to_flux(gasmf, Ts, VWC, Poros, SS, Zvals, Pa,
        boundary=None, Dcoeff=1.47e-5, method='layer',
        Ds_model=soil_diff_moldrup_1999, adjust_Da=True, makeplots=False):
    '''
    This function assembles a dataset (python dict) from a diffusive soil
    trace gas profile and related measurements (soil T, VWC, pressure, soil 
    texture, etc) and calls a requested flux measurement method.

    Args :
        gasmf    : (list/ndarray) soil gas concentrations [PPM] at depth X
        Ts       : (list/ndarray) soil temperature [Celsius] at depth X
        VWC      : (list/ndarray) vol. soil water content at depth X (fraction)
        Poros    : (list/ndarray) porosity for each depth interval) 
        SS       : (list/ndarray) silt + sand content for each depth interval
        Zvals    : (list) Soil depths in m for each depth
        Pa       : (float/list/ndarray) atm pressure in hPa, gets converted
        boundary : (float/list/ndarray)
        Dcoeff   : (numeric) free-air diffusion coefficient for gas species
        method   : (string) method to calculate
        Ds_model : function for calculating gas diffusivity
        adjust_Da: (bool) flag for to allow T and P correction of free air gas
                   diffusivity
        makeplots: (bool)

    Returns:
        2 ndarrays: 
        flux   : flux [umol/m2/s] at depth intervals d0d1, d1d2, d2d3, d3d4...
        Ds  : gas diffusivity values [m2/s2] at d0d1, d1d2, d2d3, d3d4...
    '''
    # Convert inputs to arrays and put in dict
    keys = ['gasmf','soilt','vwc','poros','ss','z','patm', 'bcond']
    invars = [gasmf, Ts, VWC, Poros, SS, Zvals, Pa, boundary]
    prof = {k: np.array(v, ndmin=1) for k, v in zip(keys, invars)}
    
    # Check for 1d arrays and convert to array(1,d) where d = # of depths
    for k in keys[0:-2]:
        if len(prof[k].shape) < 2:
            prof[k] = prof[k][np.newaxis, :]
    # patm and bcond must be array(n,1) where n = # of profile observations
    for k in keys[-2:]:
        if len(prof[k].shape) < 2:
            prof[k] = prof[k][:, np.newaxis]
        if prof[k].shape[0] < prof['gasmf'].shape[0]:
            prof[k] = np.repeat(prof[k], prof['gasmf'].shape[0], axis=0)
    # If present, concatenate the boundary condition to the gas profile,
    # otherwise the surface flux will be interpolated
    if prof['bcond'][0,0] is not None:
        prof['gasmf'] = np.c_[prof['bcond'], prof['gasmf']]
        interpsurf = False
    else:
        interpsurf = True

    # Get number of depths for each incoming data stream
    prof_d = {k:prof[k].shape[1] for k in keys[:-2]}
    
    # Calculate flux with equal gas and non-gas observation depths
    if all(prof_d[k]==prof_d['gasmf'] for k in keys[1:-2]):
        print("Calculating flux for " + str(prof_d['gasmf']) + " depths.")
    # If upper non-gas observations (boundary) are missing fill in using
    # uppermost profile values
    elif all(prof_d[k]==prof_d['gasmf']-1 for k in keys[1:-2]):
        print("Calculating flux for boundary and " +
                str(prof_d['gasmf']-1) + " depths.")
        print("Boundary VWC, Tsoil, and soil texture derived from shallowest"
                " measurements")
        for k in ['soilt', 'vwc', 'poros', 'ss']:
            prof[k] = np.c_[prof[k][:,0], prof[k]]
        prof['z'] = np.insert(prof['z'], 0, 0, axis=1)
    # Mismatched observation depths raise an error
    else:
        raise ValueError("Different # of gas, TS, VWC, S, poros, "
                "and z given!")
    # Now call the requested flux calculation
    if method is 'layer':
        flux, Ds = profile_flux_layer(prof, surfinterp=interpsurf,
                Ds_model=Ds_model, adjust_Da=adjust_Da,
                Dcoeff=Dcoeff, makeplots=makeplots)
    elif method is 'inverse':
        print('In development!!!')
        flux = np.nan
        Ds = np.nan
        #flux, Ds = profile_flux_inverse_model(prof, surfinterp=interpsurf,
        #        Ds_model=soil_diff_moldrup_1999,adjust_Da=adjust_Da,
        #        Dcoeff=Dcoeff, makeplots=makeplots)
    else:
        raise ValueError('No method named {0}'.format(method))

    return flux, Ds


def profile_flux_layer(profdat,surfinterp=False,
        Ds_model=soil_diff_moldrup_1999, Dcoeff=1.47e-5,
        adjust_Da=True,makeplots=False):
    '''
    This function computes the soil gas flux using the gradient method,
    i.e Fick's first law, with soil gas diffusivity computed using one of
    several available methods.
    
    It is a refined version of gradient_flux_prod generalized for any gas
    diffusing through porous media (soil). This process is most similar to the
    methods used by Tang et al 2005 and Vargas et al 2008

    Returns a layer-based flux profile and soil gas diffusivity profile. 
    Within-soil fluxes are the product of the gradient and calculated Ds
    between belowground sensors, giving a result that is the average flux
    between measurement depth i and i+1. If a boundary condition is 

    For an integrated profile approach see the inverse model method (though
    results from both should be comparable).
    
    Args :
        profdat     : (dict) 
        Dcoeff      : (numeric) free-air diffusion coefficient for gas species
        z_vals      : (optional list) Soil depths in m for each depth
        Ds_func     : function for calculating gas diffusivity
        adjust_Da   : (bool) flag to T and P correct free air gas diffusivity

    Returns:
        Two ndarrays: 
        flux    : flux [umol/m2/s] at depth intervals d0d1, d1d2, d2d3, d3d4...
        Ds      : gas diffusivity values [m2/s2] at d0d1, d1d2, d2d3, d3d4...
    '''
    # constants and parameters
    R = 8.3144 # Universal gas constant
    # conversions of T and P
    TsK = profdat['soilt'] + 273.15
    SVWC = profdat['vwc']
    Pa = profdat['patm'] * 100  # Transform from hPa to Pascals (Pa)
    #Transform gas from ppm(umol mol-1) to mole concentration (umol m-3)
    gas_mc = (profdat['gasmf']*Pa)/(R*TsK)
    SVWC = profdat['vwc']
    zvals = profdat['z']
    poros = profdat['poros']
    S = profdat['ss']
    # Compute dCdz (gradient) and diffusivity (Ds) for profile, then flux
    dCdz = gas_profile_to_dCdz(gas_mc, zvals)
    Ds = soil_profile_to_Ds(TsK, SVWC, poros, S, Pa, Ds_method=Ds_model)
    flux = Ds*dCdz
    return flux, Ds



def profile_flux_inverse_model(profdat,surfinterp=False,
        Ds_model=soil_diff_moldrup_1999, Dcoeff=1.47e-5,
        adjust_Da=True,makeplots=False):
    '''
    This function computes the soil gas flux using the gradient method,
    i.e Fick's first law, with soil gas diffusivity computed using one of
    several available methods.
    
    The gradient is solved for at any point in the soil allowing flux
    estimation at for any depth.

    This work is based on the Inverse Model in th origin_files directory
    by Laura Lammers (translated from MATLAB), and it has close analogs in 
    the literature in the work of Davidson and Trumbore 1995, Gaudinski et
    al. 2000, and Davidson et al 2006.

    
    Args :
        gas_mf   : (ndarray) soil gas concentrations [PPM] at depth X
        Ts       : (ndarray) soil temperature [Celsius] at depth X
        SVWC     : (ndarray) volumetric soil water content at depth X (fraction)
        P        : (float or ndarray) atm pressure in hPa, gets converted
        poros    : (list) porosity for each depth interval) 
        S:       : (list) silt + sand content for each depth interval
        Dcoeff   : (numeric) free-air diffusion coefficient for gas species
        z_vals   : (optional list) Soil depths in m for each depth
        Ds_func  : function for calculating gas diffusivity
        adjust_Da: (bool) flag for to allow T and P correction of free air gas
                   diffusivity

    Returns:
        Three ndarrays: 
        F_out   : flux [umol/m2/s] at depth intervals d0d1, d1d2, d2d3, d3d4...
        Ds_out  : gas diffusivity values [m2/s2] at d0d1, d1d2, d2d3, d3d4...
    '''
    # constants and parameters
    R = 8.3144 # Universal gas constant
    # conversions of T and P
    TsK = profdat['soilt'] + 273.15
    SVWC = profdat['vwc']
    Pa = profdat['patm'] * 100  # Transform from hPa to Pascals (Pa)
    #Transform gas from ppm(umol mol-1) to mole concentration (umol m-3)
    gas_mc = (profdat['gasmf']*Pa)/(R*TsK)
    SVWC = profdat['vwc']
    zvals = profdat['z']
    poros = profdat['poros']
    S = profdat['ss']
    # Compute dCdz (gradient) and diffusivity for the profile, then flux
    dCdz = gas_profile_to_dCdz(gas_mc, zvals)
    Ds = soil_profile_to_Ds(TsK, SVWC, poros, S, Pa,
            Ds_method=Ds_model)
    flux = Ds*dCdz
    return flux, Ds


def gas_profile_to_dCdz(C_array, z_obs, full_gradient=True):
    """
    Compute gas concentration gradient (dCdz) from an array of concentration
    observations and observation depth information.

    Input rows are observations (date) and input columns are depths

    Args:
        C_array : (ndarray) array of concentration/depth observations in ppm
        z_obs   : (ndarray) array of observation depths in meters
    Returns:
        dCdz    : (ndarray) conc. gradient for each observation/depth interval
        
    """
    dC = np.diff(C_array, axis=1)
    dz = np.diff(z_obs)
    dCdz = dC/dz
    if full_gradient:
        fulldC = C_array[:,-1] - C_array[:,0]
        fulldz = dz.sum() 
        fdCdz = fulldC/fulldz
        dCdz = np.concatenate((dCdz, fdCdz[:,np.newaxis]), axis=1)
    return dCdz

def soil_profile_to_Ds(Ts_array, VWC_array, Poros, SS, Pa,
        Dcoeff=1.47e-5, Ds_method=soil_diff_moldrup_1999, full_profile=True):
    """
    Compute trace gas diffusivity through soil (Ds) given soil temperature,
    water content, and texture observations along a depth profile. Several
    methods are available for this calculation.

    Input rows are observations (date) and input columns are depths

    If full_profile is set (bool), the last column is the total profile mean

    Args:
        Ts_array    : (ndarray) soil T observations in celsius
        VWC_array   : (ndarray) soil water content in m3/m3
        Poros       : (ndarray) soil porosity observations
        Ts_array    : (ndarray) soil sand + silt fraction
        Zvals       : (ndarray) observation depths in meters
        Pa          : (ndarray) atmospheric pressure in Pascals
    Returns:
        Ds    : (ndarray) soil gas diffusion coeffient for each
                depth interval
    """

    n_dates = Ts_array.shape[0]
    # Get interval means for soil profile measurements
    poros = profile_z_interval_mean(Poros,
            rep_dates=n_dates, full_profile=full_profile)
    ss = profile_z_interval_mean(SS,
            rep_dates=n_dates, full_profile=full_profile)
    soilt = profile_z_interval_mean(Ts_array, full_profile=full_profile)
    vwc = profile_z_interval_mean(VWC_array, full_profile=full_profile)
    # vvv Why do I need this vvv ?
    #z_obs = z_interval_size(Zvals, rep_dates=n_dates,
    #        full_profile=full_profile)
    af_poros = get_airfilled_poros(vwc, poros)
    Da_adj = get_adjusted_Da(soilt, Pa, coeff=Dcoeff)
    Ds = Ds_method(Da_adj, af_poros, poros, ss)
    return Ds

def profile_z_interval_mean(profile, rep_dates=None, full_profile=False):
    """
    Convert a profile of measurements to the mean of each depth interval.

    If rep_dates is set (int), the result is replicated for the specified
    number of rows (each representing a date)

    If full_profile is set (bool), the last column is the total profile mean
    """
    if len(profile.shape) < 2:
        profile = profile[np.newaxis, :]
    # Mean of each interval along profile
    zint_m = (profile[:,1:] + profile[:,:-1]) / 2.
    # If requested get the mean between the shallowest and deepest value
    # Not sure of best way to do this yet - depth-weighted mean may be better
    if full_profile:
        #fullprofmean = (zint_m[:,0] + zint_m[:,-1]) / 2.
        fullprofmean = profile[:,1::].sum(axis=1)/(profile.shape[1]-1)
        zint_m =  np.c_[zint_m, fullprofmean]
    if rep_dates is not None:
        zint_m = np.tile(zint_m, (rep_dates, 1))
    return zint_m

def z_interval_size(z_obs, rep_dates=None, full_profile=False):
    """
    Convert a profile of discrete depths to the size of each interval,
    assuming the surface is 0cm.

    If rep_dates is set (int), the result is replicated for the specified
    number of rows (each representing a date)

    If full_profile is set (bool), the last column is the total profile size
    """
    if len(z_obs.shape) < 2:
        z_obs = z_obs[np.newaxis, :]
    # Distance between observations
    zint_s = np.diff(z_obs)
    # If requested get the distance between the shallowest and deepest depth
    if full_profile:
        zint_s =  np.c_[zint_s, zint_s.sum()]
    if rep_dates is not None:
        zint_s = np.tile(zint_s, (rep_dates, 1))
    return zint_s

def production_from_flux_profile(df):
    ## Calculation of the soil gas production by layer
    # According to Davidson and Trumbore 1995 and Davidson et. al 2006 (GCB)
    # CO2 production of a layer can be calculated by subtracting the flux at 
    # one depth from the flux of depth below (Fi - Fi+1). 
    # Seems to give a production value for the between sensor
    # depth in micromoles per m^2
    # Could also be possible to give flux per unit volume by dividing by
    # depth between sensors (in m)
    P_out = df.copy()
    n_gasdepth = df.shape[1]
    for i in range(n_gasdepth):
        P = df.iloc[:,i] - df.iloc[:,i+1]
        P_out.iloc[:,i] = P

def interpolate_surface_flux(fprofile, z_obs, full_profile=False):
    """
    Calculate the surface flux from a flux profile (or a timeseries of 
    profiles) by linearly interpolating to the top of the profile (0m)
    """
    nobs = fprofile.shape[0]
    ndep = len(z_obs)
    fluxdeps = profile_z_interval_mean(z_obs, full_profile=full_profile)
    F0 = np.zeros((nobs, 1)) * np.nan
    # Loop through observations and fit surface flux to flux profile
    for i in range(nobs):
        if ndep > 3:
            fluxes = fprofile[i,0:3]
            fluxes_z = fluxdeps[0,0:3]
        elif ndep==3 and full_profile:
            fluxes = np.array([fprofile[i,0],fprofile[i,2],fprofile[i,1]])
            fluxes_z = np.array([fluxdeps[0,0],fluxdeps[0,2],fluxdeps[0,1]])
        elif ndep==3 and not full_profile:
            fluxes = fprofile[i,:]
            fluxes_z = fluxdeps[0,:]
        else:
            # Reorder to use the middle (full profile) flux value in
            raise ValueError('Not enough values to interpolate')
        # Check for missing values
        idx = np.isfinite(fluxes_z) & np.isfinite(fluxes)
        # Linear fit, intercept is the surface flux
        try:
            p = np.polyfit(fluxes_z[idx], fluxes[idx], 1)
            F0[i]=p[1]
        except:
            F0[i]=np.nan
    return np.c_[F0, fprofile]
