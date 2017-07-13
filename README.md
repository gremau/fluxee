# ecoflux_tools

WIP!

This is a collection of tools written in python and R that can be used for calculating ecosystem fluxes (CO2, H2O, etc), unit and measurement conversions, and various statistics (uncertainty, etc) with environmental data. These have been collected and adapted from various sources and all operations should have appropriate references.

## Classes:

### Projects

The top-level class. Configured with `proj.yml`. Each project may consist of many individual research sites, each of which may have many measurements from dataloggers or other sources.

### Site

Class referring to a research location. Configured with `sitename.yaml`. The location has a bunch of attributes and a set of measurement locations. It may contain zero or more dataloggers, and zero or more data tables.

Within each site, measuremnets are identified with H_V_R (horizontal, vertical, replicate) notation.

### Datalogger

### Datatable


