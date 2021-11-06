# fluxee 

[WIP] This is a collection of tools written in python that can be used for 
calculating ecosystem fluxes (CO2, H2O, etc), unit and measurement conversions,
statistical analysis (uncertainty, etc), and plotting with environmental data.
Methods have been collected and adapted from various sources and all operations
should have appropriate references (in docstrings).

## What is in here?

so far....

### Conversions and scaling of fluxes

Convert a flux unit like umol/m^2/s to g/m^2/day and other fun stuff.

### Gas fluxes in soils (or other porous media)

When you have a set of gas concentration measurements that spans the depth
profile of a soil or other porous media, you can calculate diffusive fluxes
and production or consumption of gases throughout the profile. The
`soil_gas_profile` helps with this.

### Assorted gapfilling and plots

Probably some of this will get axed...

## Install

Fluxee requires `pandas`, `matplotlib`, and `ruamel.yaml` to work properly. You
can install with pip.

    pip install git+git://github.com/gremau/fluxee@main
