## TSMP-PDAF Environment variables ##

A number of environment variables are set during the build process of
TSMP-PDAF in order to make the code behave in specific ways.

TSMP-PDAF environemnt variables are set in the machine-specific
TSMP-PDAF build-scripts of the form:

```
bldsva/intf_DA/pdaf1_1/arch/<machine>/build_interface_pdaf1_1_<machine>.ksh
```

where `<machine>` could be f.e. `JURECA` or `JUWELS`.

### CLMSA ###

Environment variable `CLMSA` is true if CLM-standalone is used in
TSMP-PDAF, i.e. no coupling to other component models (ParFlow,
atmospheric model).

`CLMSA` is used in many places in TSMP-PDAF, where CLM-standalone
specific code is introduced. This includes
- observation reading
- setting observation vector
- setting state vector
- communicator handling
- localized filters
- TSMP-PDAF-wrapper routines

### PARFLOW_STAND_ALONE ###

Environment variable `PARFLOW_STAND_ALONE` is true if
ParFlow-standalone is used in TSMP-PDAF, i.e. no coupling to other
component models (CLM, atmospheric model).

It is used less frequently than [CLMSA](#clmsa), only at code places
where the behavior of ParFlow-CLM-PDAF and ParFlow-PDAF should differ.

### OBS_ONLY_PARFLOW ###

Environment variable `OBS_ONLY_PARFLOW` is true if observations in
TSMP-PDAF are of ParFlow-type.

This will remove unnecessary code during observation reading, when
ParFlow-CLM-PDAF is built, but no CLM-type observations are included.

### OBS_ONLY_CLM ###

Environment variable `OBS_ONLY_CLM` is true if observations in
TSMP-PDAF are of CLM-type.

This will remove unnecessary code during observation reading, when
ParFlow-CLM-PDAF is built, but no ParFlow-type observations are
included.

