# TSMP-PDAF Environment variables #

A number of environment variables are set during the build process of
TSMP-PDAF in order to make the code behave in specific ways.

TSMP-PDAF environemnt variables are set in the machine-specific
TSMP-PDAF build-script:

```
bldsva/intf_DA/pdaf/arch/build_interface_pdaf.ksh
```

## CLMSA ##

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

## PARFLOW_STAND_ALONE ##

Environment variable `PARFLOW_STAND_ALONE` is true if
ParFlow-standalone is used in TSMP-PDAF, i.e. no coupling to other
component models (CLM, atmospheric model).

It is used less frequently than [CLMSA](#clmsa), only at code places
where the behavior of ParFlow-CLM-PDAF and ParFlow-PDAF should differ.

## OBS_ONLY_PARFLOW ##

Environment variable `OBS_ONLY_PARFLOW` is true if observations in
TSMP-PDAF are of ParFlow-type.

This will remove unnecessary code during observation reading, when
ParFlow-CLM-PDAF is built, but no CLM-type observations are included.

## OBS_ONLY_CLM ##

Environment variable `OBS_ONLY_CLM` is true if observations in
TSMP-PDAF are of CLM-type.

This will remove unnecessary code during observation reading, when
ParFlow-CLM-PDAF is built, but no ParFlow-type observations are
included.

## WATSAT3D ##

Environment variable `WATSAT3D` is set in `common_build_interface.ksh`
in function `c_configure_clm`.

If it is turned on, the possibility of a read-in porosity is
implemented in CLM's `iniTimeConst.F90`.

Usage: Enforce consistent porosity in ParFlow and CLM.

## CLMFIVE ##

Currently only in feature branch `TSMP_pdaf-clm5`.

If defined, CLM5.0 is used.

If undefined, CLM3.5 is used.

This distinction is important in many parts of the TSMP-PDAF-wrapper
source code, when, f.e., function calls have changed from version 3 to
version 5 of CLM.

## OLD_TRUNCATE_SAT ##

If `OLD_TRUNCATE_SAT` is defined, the saturation truncation in
function `SaturationToPressure` (file
`problem_saturationtopressure.c`) is changed to an outdated default
(before October 2023).

Reference for the current default saturation truncation in Hung et al
(<https://doi.org/10.1029/2021WR031549>, Equation 11)

``` c++
	if(psdat[ips] <= (s_res + 0.003) ) psdat[ips] = s_res + 0.003;
```

`OLD_TRUNCATE_SAT` will change this to the outdated default:

``` c++
	if(psdat[ips] <= s_res) psdat[ips] = s_res + 0.01;
```
