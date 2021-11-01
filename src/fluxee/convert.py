import pandas as pd

def co2_mol_to_C_mass_flux( df, n_seconds ) :
    """
    Convert molar CO2 flux to mass C flux ( umolCO2/m^2/s to gC/m^2 )
    and sum (integrate) for each period in the timeseries for each column
    in data frame
    """
    # Initialize returned variables
    df_int = df.copy()
    export_cols = []
    # For each input column create a new header and convert values to mass flux
    for cname in df.columns :
        export_cname = cname + '_g_int'
        export_cols.append( export_cname )
        # Integrate based on number of seconds
        df_int[ export_cname ] = df_int[ cname ] * ( 12.011/1e+06 ) * n_seconds

    return df_int[ export_cols ]

def mol_to_mass_flux( df, molarmass, n_seconds, mol_denom=1e+06 ) :
    """
    Similar to above, but abstracted for other mass fluxes. Expects incoming
    flux to be in micromoles (1e06), but this can be changed

    Convert molar flux to mass flux ( umolX/m^2/s to gX/m^2/n_seconds )
    """
    # Initialize returned variables
    df_int = df.copy()
    export_cols = []
    # For each input column create a new header and convert values to mass flux
    for cname in df.columns :
        export_cname = cname + '_g_int'
        export_cols.append( export_cname )
        # Integrate based on number of seconds
        df_int[ export_cname ] = (df_int[ cname ] * ( molarmass/mol_denom ) * 
                n_seconds)

    return df_int[ export_cols ]

def umol_m2_s_to_kg_ha_yr(df, molarmass, n_seconds=1):
    # First convert to g/m2 and then
    # Initialize returned variables
    df_mass = mol_to_mass_flux(df, molarmass, 1)
    df_int = df_mass.copy()
    export_cols = []
    # For each input column create a new header and convert values to mass flux
    for i, cname in enumerate(df_mass.columns) :
        export_cname = df.columns[i] + '_kgHaYr'
        export_cols.append( export_cname )
        # Integrate based on number of seconds
        df_int[ export_cname ] = ((df_mass[ cname ] * 10000 * 
                ((60*60*24*365)/n_seconds))/1000)
    return df_int[ export_cols ]


