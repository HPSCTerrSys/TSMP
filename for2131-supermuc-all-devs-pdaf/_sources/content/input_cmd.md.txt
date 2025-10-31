## Command line options ##

The following command line options must be specified when TSMP-PDAF is
executed:

### n_modeltasks ###

`n_modeltasks` (integer) Number of realisations. Must be consistent
with `[DA] nreal`.

### filtertype ###

`filtertype` (integer) Type of filter used for data assimilation. For
more details see the PDAF documentation. Currently, the following
filter options are available

-   (Ensemble Kalman Filter(Enkf)) is supported.

-   (Ensemble Transform Kalman Filter(Etkf)) is supported.

-   (Local Ensemble Transform Kalman Filter(Letkf)) is supported.

-   (Error Subspace Transform Kalman Filter(Estkf)) is supported.

-   (Local Error Subspace Transform Kalman Filter(Lestkf)) is supported.

-   (Local Ensemble Kalman Filter(LEnkf)) is supported.

To see all available command line options for the different filter
algorithms please refer to PDAF wiki link

<http://pdaf.awi.de/trac/wiki/AvailableOptionsforInitPDAF>

### subtype ###

`subtype` (integer) Parameter subtype, different options for each
filter. See [Command Line Examples](#command-line-examples).

### obs_filename ###

`obs_filename` (string) Prefix for observation files.

### rms_obs ###

`rms_obs` (real) Measurement error.

### delt_obs ###

`delt_obs` (integer) Number of data assimilation intervals (see
[`[DA]da_interval`](./input_enkfpf.md#dada_interval)) that are forward
computed, before assimilation is actually performed.

### screen ###

`screen` (integer) Control verbosity of PDAF
- 0: no outputs
- 1: basic output (default)
- 2: 1 plus timing output
- 3: 2 plus debug output

## Command Line Examples ##

### Command Line Example: EnKF ###

`subtype` (integer) Parameter subtype for Ensemble Kalman Filter (EnKF).

-   Full ensemble integration; analysis for $2\cdot \mathtt{dim\_obs} > \mathtt{dim\_ens}$

-   Full ensemble integration; analysis for $2\cdot \mathtt{dim\_obs} \leq \mathtt{dim\_ens}$

-   Offline mode

Example of the model execution using EnKF analysis:

        mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 2 -subytype 1 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run with
64 realisations (8 processors for each realisation) and EnKF analysis
will be performed every assimilation cycle with an observation error of
0.1 and observations read from the files `myobs.xxxxx`.

### Command Line Example: ETKF ###

`subtype` (integer) Parameter subtype for Ensemble Transform Kalman
Filter (ETKF).

-   full ensemble integration; apply T-matrix analogously to SEIK

-   full ensemble integration; formulation without T matrix

-   Fixed error space basis; analysis with T-matrix

-   Fixed state covariance matrix; analysis with T-matrix

-   Offline mode; analysis with T-matrix

Example of the model execution using EtKF analysis:

        mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 4 -subytype 1 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run with
64 realisations (8 processors for each realisation) and EtKF analysis
will be performed every assimilation cycle with an observation error of
0.1 and observations read from the files `myobs.xxxxx`.

### Command Line Example: ESTKF ###

`subtype` (integer) Parameter subtype for Error Subspace Transform Kalman Filter
(ESTKF).

-   Standard implementation with ensemble integration

-   Fixed error space basis

-   Fixed state covariance matrix

-   Offline mode

Example of the model execution using EstKF analysis:

        mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 6 -subytype 1 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run with
64 realisations (8 processors for each realisation) and EstKF analysis
will be performed every assimilation cycle with an observation error of
0.1 and observations read from the files `myobs.xxxxx`.

### Command Line Example: LETKF ###

`local_range` (real) set localization radius, 0.0 by default, any
positive value should work.

`locweight` (integer) set weight function for localization, default=0
for constant weight of 1; possible are integer values 0 to 4.

`subtype` (integer) Parameter subtype for Local Ensemble Transform
Kalman Filter (Letkf).

-   full ensemble integration; apply T-matrix analogously to SEIK

-   Fixed error space basis; analysis with T-matrix

-   Fixed state covariance matrix; analysis with T-matrix

-   Offline mode; analysis with T-matrix.

Example of the model execution using Letkf analysis:

        mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 5 -local_range 3 -locweight 0 -subtype 0 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run with
64 realisations (8 processors for each realisation) and Letkf analysis
will be performed every assimilation cycle with an observation error of
0.1 and observations read from the files `myobs.xxxxx`.

### Command Line Example: LESTKF ###

`local_range` (real) set localization radius, 0.0 by default, any
positive value should work.

`locweight` (integer) set weight function for localization, default=0
for constant weight of 1; possible are integer values 0 to 4.

`subtype` (integer) Parameter subtype for Local Error Subspace Transform Kalman
Filter (LESTKF).

-   Standard implementation with ensemble integration

-   Fixed error space basis

-   Fixed state covariance matrix

-   Offline mode

Example of the model execution using Lestkf analysis:

        mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 7 -local_range 3 -locweight 0 -subtype 0 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run with
64 realisations (8 processors for each realisation) and Lestkf analysis
will be performed every assimilation cycle with an observation error of
0.1 and observations read from the files `myobs.xxxxx`.

### Command Line Example: LENKF ###

`local_range` (real) set localization radius (cut-off-radius)
- 0.0 by default
- any positive value should work.

`srange` (integer): the support radius of the localization
- the localization weight radius; support range for 5th-order
  polynomial
- by default set to `local_range`

`locweight` (integer): set weight function for localization,
- default=0
- for constant weight of 1
- possible are integer values 0 to 4.(see `init_pdaf`)

``` text
locweight == 0    ! Uniform (unit) weighting
locweight == 1    ! Exponential weighting
locweight == 2    ! 5th-order polynomial (Gaspari&Cohn, 1999)
```

`subtype` (integer) Parameter subtype for Localized Ensemble Kalman
Filter (Lenkf).
-   Full ensemble integration; analysis with covariance localization
-   Offline mode

Example of the model execution using Lenkf analysis:

``` bash
mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 8 -local_range 3 -locweight 0 -subtype 0 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs
```

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run
with 64 realisations (8 processors for each realisation) and Lenkf
analysis will be performed every assimilation cycle with an
observation error of 0.1 and observations read from the files
`myobs.xxxxx`.
