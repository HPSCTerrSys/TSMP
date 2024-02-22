# Command line options #

The following command line options must be specified when TSMP-PDAF is
executed:

## n_modeltasks ##

`n_modeltasks` (integer) Number of realisations. Must be equal to
[`[DA] nreal`](./input_enkfpf.md#danreal).

Note that both `n_modeltasks` (or `nreal`) and the total number of MPI
tasks specified, f.e. with `srun -n` or `mpiexec -np` (internally:
`npes_world`) have to match with the total number of processes
specified for the component models:
[`PF:nprocs`](./input_enkfpf.md#pfnprocs),
[`CLM:nprocs`](./input_enkfpf.md#clmnprocs),
[`COSMO:nprocs`](./input_enkfpf.md#cosmonprocs).

## filtertype ##

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

## subtype ##

`subtype` (integer) Parameter subtype, different options for each
filter. See [Command Line Examples](#command-line-examples).

## obs_filename ##

`obs_filename` (string) Prefix for observation files.

## rms_obs ##

`rms_obs` (real) Measurement error. This is the standard deviation of
the Gaussian probability distribution of the measurements.

## delt_obs ##

`delt_obs` (integer) Number of data assimilation intervals (see
[`[DA]da_interval`](./input_enkfpf.md#dada_interval)) that are forward
computed, before assimilation is actually performed.

In general, `delt_obs` should be as small as possible, in order to
avoid performance loss. See remark in
[`[DA]da_interval`](./input_enkfpf.md#dada_interval).

## screen ##

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

`cradius` (deprecated: `local_range`) (real) set cut-off radius, 0.0
by default, any positive value should work.

`locweight` (integer) set weight function for localization, default=0
for constant weight of 1; possible are integer values 0 to 4.

`subtype` (integer) Parameter subtype for Local Ensemble Transform
Kalman Filter (Letkf).

-   full ensemble integration; apply T-matrix analogously to SEIK

-   Fixed error space basis; analysis with T-matrix

-   Fixed state covariance matrix; analysis with T-matrix

-   Offline mode; analysis with T-matrix.

Example of the model execution using Letkf analysis:

        mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 5 -cradius 3 -locweight 0 -subtype 0 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run with
64 realisations (8 processors for each realisation) and Letkf analysis
will be performed every assimilation cycle with an observation error of
0.1 and observations read from the files `myobs.xxxxx`.

### Command Line Example: LESTKF ###

`cradius` (deprecated: `local_range`) (real) set cut-off radius, 0.0 by default, any
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

        mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 7 -cradius 3 -locweight 0 -subtype 0 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run with
64 realisations (8 processors for each realisation) and Lestkf analysis
will be performed every assimilation cycle with an observation error of
0.1 and observations read from the files `myobs.xxxxx`.

### Command Line Example: LENKF ###

`cradius` (deprecated: `local_range`) (real) set cut-off radius
- 0.0 by default
- any positive value should work.

`sradius` (deprecated `srange`) (real): the support radius of the
localization, default: equal to `cradius`. Usage in
`PDAF_local_weight.F90`
(<https://pdaf.awi.de/trac/wiki/PDAF_local_weight>)
- `locweight == 0`: `sradius` is not used
- `locweight == 1`: `sradius` is distance with weight `1/e`, for
  avoiding sharp cut-off, `cradius` should be significantly larger
  than `sradius`
- `locweight == 2`: `sradius` is the support radius for 5th-order
  polynomial, can safely be chosen the same as `cradius`,
  parametrization change of the 5th-order polynomial is at `sradius/2`,
  see
  <https://github.com/PDAF/PDAF/blob/ed631034956dece8e91e8b588c4cf3aaa7916f49/src/PDAF_local_weight.F90#L147-L176>.

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
mpiexec -np 512 ./tsmp-pdaf -n_modeltasks 64 -filtertype 8 -cradius 3 -locweight 0 -subtype 0 -delt_obs 1 -rms_obs 0.1 -obs_filename myobs
```

With this command, the TSMP-PDAF executable `tsmp-pdaf` will be run
with 64 realisations (8 processors for each realisation) and Lenkf
analysis will be performed every assimilation cycle with an
observation error of 0.1 and observations read from the files
`myobs.xxxxx`.
