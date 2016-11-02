import numpy as np

def fluxprod_3(CO2_d1,CO2_d2,CO2_d3,Ts_d1,Ts_d2,Tsd3,SVWC,P,porosity,S)
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
        CO2_dX  : (ndarray) soil CO2 concentrations [PPM] at depth X
        Ts_dX   : (ndarray) soil temperature [Celsiuis] at depth X
        SVWC_dX : (ndarray) volumetric soil water content at depth X
        P       : (ndarray) atmospheric pressure [Kpa]
        S:      : (float) silt + sand content
        porosity: (float) porosity (see Moldrup model and Vargas papers, but
                  this may need to be an array since it changes with VWC)

    Outputs (numpy arrays):
        F0   : (ndarray) surface CO2 flux  [mumol/m2/s]
        Pco2 : (ndarray) production of CO2 between 2 sensors [mumol/m3/s]
    
    Usage: [F0 Pco2]=fluxprod_3(CO2_d1,CO2_d2,CO2_d3,Ts_d1,Ts_d2,Tsd3,
                                SVWC,P,porosity,S)
    '''
    import numpy

    # constants and parameters
    beta = 2.9 # Constant for Moldrup model
    R = 8.3144 # Universal gas constant
    z = [.05 .10 .20] #soil depth at where CO2 sensors have installed.

    z2(1)=(z(2)+z(1))/2
    z2(2)=(z(3)+z(1))/2
    z2(3)=(z(3)+z(2))/2
    
    z3(1)=(z2(2)+z2(1))/2
    z3(2)=(z2(3)+z2(1))/2
    z3(3)=(z2(3)+z2(2))/2
    
    #Transform T soil in Kelvin
    Ts_d1k = Ts_d1 + 273.15
    Ts_d2k = Ts_d2 + 273.15
    Ts_d3k = Ts_d3 + 273.15
    Pa = P* 1000  # Transform air pressure from kPa to Pascals (Pa)
    
    #Transform from ppm(umol mol-1) to mole concentration (umol m-3)
    # Note: Air pressure should be in Pascals
    C_d1 = (CO2_d1*Pa) /(R*Ts_d1k)
    C_d2 = (CO2_d2*Pa) /(R*Ts_d2k)
    C_13 = (CO2_d3*Pa)/(R*Ts_d3k)

    # /////////////////////////////////////////////////////////////
    # calculation of CO2 DIFFUSION in the soil
    # The CO2 diffusion coefficient in the soil, Ds, is equal to the product of
    # CO2 diffusion coefficient in the free air, Da, and the gas tortuosity 
    # factor .... Temperature(oK) Pressure (kPa)

    #Da for layers Da_d1_d2 (2-8 cm) Da16 (8-16 cm), Da2_16 (2-16 cm)
    Da_d1_d2 = (1.47e-5*((((Ts_d1k + Ts_d2k)/2)/293.15)**1.75)*
            (1.013e+5/(Pa)))
    Da_d2_d3 = (1.47e-5*((((Ts_d3k +Ts_d2k)/2)/293.15)**1.75)*
            (1.013e+5/(Pa)))
    Da_d1_d3 = (1.47e-5*((((Ts_d3k +Ts_d1k)/2)/293.15)**1.75)*
            (1.013e+5/(Pa)))
    
    # Diffusivity calculation using Moldrup et al. model 1999
    
    # Calculate the air-filled porosity
    Ds_Da_d1_d2 =  (porosity**2)*(((porosity-((SVWC)))/porosity)**
            (beta*S))* Da_d1_d2;
    Ds_Da_d2_d3 = (porosity**2)*(((porosity-((SVWC)))/porosity)**
            (beta*S))* Da_d2_d3;
    Ds_Da_d1_d3 = (porosity**2)*(((porosity-((SVWC)))/porosity)**
            (beta*S))* Da_d1_d3;
    
    #=== Calculation of the Flux at each layer ===
    
    # calculation of CO2 flux from Fd2_d1,  Fd3_d2, Fd3_d1......
    # Note: depth should be in meters
    # z1 = 0.06..... or distance in between 0.08 and 0.02 m
    # z2 = 0.09..... or distance in between 0.16 and 0.02 m
    # z2 = 0.12..... or distance in between 0.16 and 0.08 m

    F1 = Ds_Da_d1_d2*((C_8-C_2)/(z(2)-z(1)));
    F2 = Ds_Da_d1_d3*((C_16-C_2)/(z(3)-z(1)));
    F3 = Ds_Da_d2_d3*((C_16-C_8)/(z(3)-z(2)));


    #=== Calculation of the Flux at soil surface ===
    
    for i=1:length(F1) #Calculation of means for Node 1
        #length(F1)-i
        p = polyfit(z2,[F1(i) F2(i) F3(i)], 1)
        F0(i)=p(2)

    #F0=F0'

    #%%%%%Calculation of the soil CO2 production%%%%
    # p1 = 0.07..... or distance in between 0.05 and 0.09 m
    # p2 = 0.085..... or distance in between 0.05 and 0.12 m
    # p3 = 0.105..... or distance in between 0.09 and 0.12 m

    PR=  (F1 - F3)/(z2(3)-z2(1));

    return([F0, PR])


