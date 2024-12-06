# Control file `enkfpf.par` #

An additional file named `enkfpf.par` needs to be present in the
TSMP-PDAF run directory.

`enkfpf.par` is read by the routine `read_enkfpar` in
`model/common/read_enkfpar.c`.

`enkfpf.par` contains information about the different model components
(ParFlow, CLM and data assimilation) like, e.g., the timing
information of the models, the prefixes for the ParFlow/ CLM input
files, etc.

This input is structured in four sections: 
- [`[PF]`](#pf) which holds information about ParFlow, 
- [`[CLM]`](#clm) which holds information about CLM 
- [`[COSMO]`](#cosmo) which holds information about COSMO
- [`[DA]`](#da) which holds information about the data assimilation process

An example for `enkfpf.par` is given below.  Note that the sequence of
the entries within these four categories can be changed if desired.

``` text
[PF]
problemname        = ""
nprocs             =
starttime          =
dt                 =
endtime            =
simtime            =
updateflag         =
gwmasking          =
paramupdate        =
paramupdate_frequency =
dampingfactor_param =
dampingfactor_state =
aniso_perm_y       =
aniso_perm_z       =
aniso_use_parflow  =
printensemble      =
printstat          =
paramprintensemble =
paramprintstat     =
olfmasking         =

[CLM]
problemname = ""
nprocs      =
update_swc  =
update_texture  =
print_swc   =
print_et   =

[COSMO]
nprocs      =
dtmult      =

[DA]
outdir = ""
nreal =
startreal =
da_interval       =
stat_dumpoffset   =
point_obs =
```

In the following the individual entries of `enkfpf.par` are described:

## [PF] ##

### PF:problemname ###

`PF:problemname`: (string) Problem prefix for ParFlow.

### PF:nprocs ###

`PF:nprocs`: (integer) Number of processors per ParFlow instance. Must
match with the specifications in the `*.pfidb` input.

This number of processors specifies a subset of the processors
available in a single `COMM_model`. Each `COMM_model` contains - if
there are no remainders - the following number of processes:
`npes_world / nreal`.

### PF:starttime ###

`PF:starttime`: (real) ParFlow start time. Must match with the
specifications in the `*.pfidb` input (`TimingInfo.StartTime`).

### PF:dt ###

`PF:dt`: (real) Length of ParFlow time step. Must match with the
specifications in the `*.pfidb` input. 

It is implicitly assumed that ParFlow and CLM calculate the same
amount of time steps (i.e., have the same time step
length). (Johannes: This may not be true, CLM time steps are computed)

### PF:endtime (deprecated) ###

Deprecated. Use `PF:simtime` instead.

`PF:endtime`: (real) Total simulation time (in terms of ParFlow
timing). Must match with the specifications in the `*.pfidb` input.

### PF:simtime ###

`PF:simtime`: (real) Total simulation time (in terms of ParFlow
timing). 

Must match with the specifications in the `*.pfidb` input.

`PF:simtime` must correspond to `TimingInfo.StopTime` MINUS
`TimingInfo.StartTime`!

### PF:updateflag ###

`PF:updateflag`: (integer) Type of state vector update in ParFlow.

-   1: Assimilation of pressure data. State vector consists of
    pressure values and is directly updated with pressure
    observations.

-   2: Assimilation of soil moisture data. State vector consists of
    soil moisture content values and is updated with soil moisture
    observations. The updated soil moisture content is transformed
    back to pressure via the inverse van Genuchten functions.

-   3: Assimilation of soil moisture data. State vector consists of
    soil moisture content and pressure values. Soil moisture data are
    used to update pressure indirectly.

### PF:gwmasking ###

`PF:gwmasking`: (integer) Groundwater masking for assimilation of
pressure data (updateflag=1) in ParFlow.

-   0: No groundwater masking.

-   1: Groundwater masking using saturated cells only.

-   2: Groundwater masking using mixed state vector.

### PF:paramupdate ###

`PF:paramupdate`: (integer) Flag for parameter update

-   0: No parameter update

-   1: Update of saturated hydraulic conductivity

-   2: Update of Mannings coefficient

- 3: Update porosity

- 4: Update van Genuchten parameters

- 5: Update hydraulic conductivity and porosity

- 6: Update hydraulic conductivity and van Genuchten parameters

- 7: Update porosity and van Genuchten parameters

- 8 Update hydraulic conductivity, porosity and van Genuchten parameters


### PF:paramupdate_frequency ###
`PF:paramupdate_frequency`: (integer) Frequency of parameter
updates. Default: `1`

For each assimilation cycle it is checked whether `tcycle mod
paramupdate_frequency == 0`. When this happens, a parameter update is
applied. For the default of `1`, each assimilation cycle contains a
parameter update.

### PF:dampingfactor_param ###
`PF:dampingfactor_param`: (real) Damping factor for parameter
updates. Default `1.0`.

The damping factor should be chosen between `0.0` and `1.0`, where the
inputs yield the following behavior:
- `0.0`: No parameter update.
- `1.0`: Parameter update without damping (default)

General equation for the damped, updated parameter vector

\begin{align*}
p_{update, damped} &= p + \mathtt{dampingfactor\_param} \cdot ( p_{update} - p )
\end{align*}

where $p$ is the parameter vector that is part of the state vector
before assimilation and $p_{update}$ is $p$ after the assimilation.

### PF:dampingfactor_state ###
`PF:dampingfactor_state`: (real) Damping factor for state
updates. Default `1.0`.

Remark: Currently implemented for `pf_updateflag==1`.

The damping factor should be chosen between `0.0` and `1.0`, where the
inputs yield the following behavior:
- `0.0`: No state update.
- `1.0`: State update without damping (default)

General equation for the damped, updated state vector

\begin{align*}
x_{update, damped} &= x + \mathtt{dampingfactor\_state} \cdot ( x_{update} - x )
\end{align*}

where $x$ is the state vector (without parameters) before assimilation
and $x_{update}$ is the state vector after the assimilation.


### PF:aniso_perm_y ###

`PF:aniso_perm_y`: (real) Anisotropy factor of saturated hydraulic
conductivity in y-direction. Only used when hydraulic conductivity is
updated (`[PF]paramupdate = 1`). Important: Compare this anisotropy
factor with the corresponding anisotropy factor from
ParFlow-input. Currently, we believe that the ParFlow-factor is used
in the beginning of the simulation, while the factor set here, is used
after each parameter update.

### PF:aniso_perm_z ###

`PF:aniso_perm_z`: (real) Anisotropy factor of saturated hydraulic
conductivity in z-direction. Only used when hydraulic conductivity is
updated (`[PF]paramupdate = 1`). Important: Compare this anisotropy
factor with the corresponding anisotropy factor from
ParFlow-input. Currently, we believe that the ParFlow-factor is used
in the beginning of the simulation, while the factor set here, is used
after each parameter update.

### PF:aniso_use_parflow ###

`PF:aniso_use_parflow`: (integer) If 1, heterogeneous anisotropy
factors (computed from ParFlow permeability arrays before the PDAF
update) are used to compute the permeability in y and z direction
based on the permeability in x direction that has been updated by
PDAF.

If 0, two homogeneous anisotropy factors are used to compute the
updated permeabilities in y and z direction. These factors are read
from the input variables `PF:aniso_perm_y` and `PF:aniso_perm_z`.

### PF:printensemble ###

`PF:printensemble`: (integer) If set to `1`, the updated state
variables for all ensemble members is printed out as `pfb` files after
each assimilation cycle.  Files are printed to `[DA]outdir` and follow
the file naming convention of ParFlow. They include the specifier
`update` in the file name.

### PF:printstat ###

`PF:printstat`: (integer) If set to `1` (default) the ensemble
statistics (mean and standard deviation) of the forecasted state
variable are calculated and printed out as `pfb` files after each
assimilation cycle. Files are printed to `DA:outdir` and follow the
file naming convention of ParFlow. The output files include the
specifiers `press.mean`/ `press.sd` in case of assimilated pressure
update `PF:updateflag = 1` and the specifiers `swc.mean`/ `swc.sd` in
case of assimilated soil moisture data `PF:updateflag = 2` or
`PF:updateflag = 3`.

### PF:paramprintensemble ###

`PF:paramprintensemble`: (integer) Only used in case of parameter
update. If set to `1`, the updated parameters are printed to `pfb`
files (similar to `[PF]printensemble`). Output files include a
specifier depending on `PF:paramupdate`, f.e. `update.param.ksat`,
`update.param.ksat`.

### PF:paramprintstat ###

`PF:paramprintstat`: (integer) Only used in case of parameter
update. If set to `1` statistics on the updated parameters are printed
to `pfb` files (similar to `[PF]printstat`).

The output files include a specifier depending on `PF:paramupdate`,
f.e.  `param.ksat`, `param.mannings` or `param.poro`.

### PF:olfmasking ###

`PF:olfmasking`: (integer) Only used in case you do not want to update
the state on certain grid-cells during DA with pdaf. eg. not update
the cell which is saturated. Option \"1\" means that all saturated
cells at surface are not used for an update. Option \"2\" reads a pfb
for masking the stream.

## [CLM] ##

### CLM:problemname ###

`CLM:problemname`: (string) Problem prefix for CLM

### CLM:nprocs ###

`CLM:nprocs`: (integer) Number of processors for each CLM instance.

This number of processors specifies a subset of the processors
available in a single `COMM_model`. Each `COMM_model` contains - if
there are no remainders - the following number of processes:
`npes_world / nreal`.

### CLM:update_swc ###

`CLM:update_swc`: (integer) Flag for update of soil moisture content
in CLM (standalone only).

-  0: No update of soil moisture content

-  1: Update of soil moisture content

- (only branch `TSMP_pdaf-crns`) 2: Update and average of soil
  moisture content for use with Cosmic-Ray data (Hui Pung
  implementation).

- (only CLM5.0) 3: Update and average of soil moisture content for use
  with Cosmic-Ray data (Strebel implementation of Schrön2017,
  <https://research-information.bris.ac.uk/en/publications/improving-calibration-and-validation-of-cosmic-ray-neutron-sensor>).

### CLM:update_texture ###

`CLM:update_texture`: (integer) Flag for update of soil parameter in
CLM (standalone only).

-  0: No update of soil parameter

-  1: Update of soil parameter (clay or sand)

-  (only CLM5.0) 2: Update of soil parameter (clay, sand or organic matter)

-  (only CLM5.0) 3: Update of hydraulic conductivity log-transformed

-  (only CLM5.0) 4: Update of hydraulic conductivity, porosity,
   suction, hydraulic conductivity exponent B (see CLM5.0 technical
   manual
   <https://escomp.github.io/ctsm-docs/versions/release-clm5.0/html/tech_note/index.html>)

### CLM:print_swc ###

`CLM:print_swc`: (integer) If set to `1`, the updated soil moisture
content in CLM for all ensemble members is printed out as `netcdf`
files (one file per realisation for the whole simulation
period). Files are printed to the run directory and are named
according to the CLM problem prefix of the corresponding realisation
and include the specifier `update` in the file name.

### CLM:print_et ###

`CLM:print_et`: (integer) Invoke function `write_clm_statistics`. For
further information, see source code.

## [COSMO] ##

### COSMO:nprocs ###

`COSMO:nprocs`: (integer) Number of processors for each COSMO
instance.

Currently, `COSMO:nprocs` is NOT USED by TSMP-PDAF. Instead, COSMO
will get processors if there are processes left after giving
`PF:nprocs + CLM:nprocs` processes to ParFlow and CLM.

### COSMO:dtmult ###

`COSMO:dtmult`: (integer) Number of COSMO time steps within one
ParFlow time step.

## [DA] ##

### DA:outdir ###

`DA:outdir`: (string) Directory where assimilation results should be
written.

### DA:nreal ###

`DA:nreal`: (integer) Number of realisations used in the
simulation. `DA:nreal` Must be equal to command line input
`n_modeltasks`.

### DA:startreal ###

`DA:startreal`: (integer) Added to suffix-numbers for input file
creation.

### DA:da_interval ###

`DA:da_interval`: (double) Time interval (units of ParFlow timing,
usually hours, check ParFlow input `TimingInfo.BaseUnit`,
<https://parflow.readthedocs.io/en/latest/keys.html#timing-information>).

After the time interval `da_interval` the forward simulation
(integration) of all component models is stopped in the main loop of
TSMP-PDAF (`pdaf_terrsysmp.f90`). Afterwards PDAF is called and it
depends on the command line option `delt_obs` ([Command line
options](./input_cmd.md#command-line-options)) whether data
assimilation is performed during this iteration. `delt_obs` will
specify the number of steps (loop iterations) that PDAF will iterate
back to the forward simulation, before performing the actual data
assimilation.

One exception is that the observation file is empty at a data
assimilation step. Then, no data assimilation takes place and the next
set forward integration steps is executed.

Formula for the simulation time between data assimilation steps (given
non-empty observation files)

\begin{gather*} 
	t_{\mathtt{betweenDA}} =  \mathtt{da\_interval} \cdot \mathtt{delt\_obs} \cdot \mathtt{BaseUnit}
\end{gather*}


**Example 1 (FallSchoolCase):** `da_interval=1.0`, `delt_obs 12`,
`TimingInfo.BaseUnit=1.0`. In this case, component models will be
simulated for 12 1-hour-steps between data assimilation times. So an
assimilation will be applied each 12 hours.

**Example 2:** `da_interval=24.0`, `delt_obs 1`,
`TimingInfo.BaseUnit=1.0`. In this case, component models will be
simulated for 1 24-hour-step between data assimilation times. So an
assimilation will be applied each 24 hours.

### DA:stat_dumpoffset ###

`DA:stat_dumpoffset`: File number offset for the data assimilation
output files for ParFlow.  Can be used when an assimilation run is
restarted.

### DA:screen_wrapper ###

`DA:screen_wrapper` (added 02/2022): Variable for controlling
output. Analogous to `screen` variable for PDAF. Values: (0) no
outputs, (1) medium outputs (default), (2) debugging output.

### DA:point_obs ###
`DA:point_obs`: (integer) Only used in case of multi-scale data
assimilation. Set to value 0 for using multi-scale data assimilation
(eg. using SMAP satellite data over a large area which is not point
data). If not specified its default value is set to 1, which is for
using point observation for data assimilation run.

### DA:obs_interp_switch ###

`DA:obs_interp_switch`: (integer) Switch for using an interpolation of
simulated measurements from the closest grid cells to the observation
location.

Default: `0` (no interpolation)

Effect of `obs_interp_switch=1`:
- For [CLM-type
  observations](./input_obs.md#clm-observation-file-variables) the
  longitude and latitude inputs in the observation file are compared
  to the four neihboring grid points and the simulated measurements
  from these grid points are averaged with the distances to the grid
  points as weights.
- For [ParFlow-type
  observations](./input_obs.md#parflow-observation-file-variables) the
  observation files need two more inputs: `ix_interp_d` and
  `iy_interp_d`. These index distances are used to compute an average
  of the simulated measurements from the value at the surrounding grid
  points in the x-y-plane.

## Parameter Summary ##

 | section   | parameter               | default value |
 |:---------:|:-----------------------:|:-------------:|
 | `[PF]`    |                         |               |
 |           | `problemname`           | \-            |
 |           | `nprocs`                | 0             |
 |           | `starttime`             | 0.0           |
 |           | `endtime`               | 0             |
 |           | `simtime`               | 0             |
 |           | `dt`                    | 0.0           |
 |           | `updateflag`            | 1             |
 |           | `paramupdate`           | 0             |
 |           | `paramupdate_frequency` | 1             |
 |           | `dampingfactor_param`   | 1.0           |
 |           | `dampingfactor_state`   | 1.0           |
 |           | `aniso_perm_y`          | 1.0           |
 |           | `aniso_perm_z`          | 1.0           |
 |           | `aniso_use_parflow`     | 0             |
 |           |                         |               |
 |           | `printensemble`         | 1             |
 |           | `printstat`             | 1             |
 |           | `paramprintensemble`    | 1             |
 |           | `paramprintstat`        | 1             |
 |           | `olfmasking`            | 0             |
 | `[CLM]`   |                         |               |
 |           | `problemname`           | \-            |
 |           | `nprocs`                | 0             |
 |           | `update_swc`            | 1             |
 |           | `print_swc`             | 0             |
 | `[COSMO]` |                         |               |
 |           | `nprocs`                | 0             |
 |           | `dtmult`                | 0             |
 | `[DA]`    |                         |               |
 |           | `nreal`                 | 0             |
 |           | `outdir`                | \-            |
 |           | `da_interval`           | 1             |
 |           | `stat_dumpoffset`       | 0             |
 |           | `screen_wrapper`        | 1             |
 |           | `point_obs`             | 1             |
 |           | `obs_interp_switch`     | 0             |

Default values for parameter file `enkfpf.par`.

