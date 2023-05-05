# Misc #

Various topics related to TSMP-PDAF.

## Gitlab: TSMP-PDAF group ##

The main repository of TSMP-PDAF is found on GitHub
([https://github.com/HPSCTerrSys/TSMP](https://github.com/HPSCTerrSys/TSMP)).

Additional development of TSMP-PDAF is organized on Gitlab in the
group TSMP-PDAF:

[https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf](https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf)

### Repositories ###

The Gitlab-group TSMP-PDAF currently contains the following
repositories

- `TSMP`: main repository, clone of GitHub TSMP from
  [https://github.com/HPSCTerrSys/TSMP](https://github.com/HPSCTerrSys/TSMP)
- `TSMP-PDAF Manual`: Repository for this manual.
- `TSMP-PDAF Build Scripts`: Scripts for (remote) building of TSMP
- `TSMP-PDAF Component Models`: Component Models used for building TSMP-PDAF
- `TSMP-PDAF Testcases`: Testcases for TSMP-PDAF
- `TSMP-PDAF Presentations`: Presentation slides for TSMP-PDAF

## Open loop compilation of (TSMP-)PDAF ##

For compiling TSMP-PDAF with PDAF but without data assimilation, you
can set the flag `PDAF_NO-UPDATE` in Makefiles under
`./TSMP/bldsva/intf_DA/pdaf1_1/arch/<machine>/config/<compiler>.h.`
and in the line, where `CPP_DEFS` are set. Usually only `-DUSE_PDAF`
is set, here you can add `-DPDAF_NO_UPDATE`.

If you compile once with this flag and once without, you should obtain
two build-versions of TSMP-PDAF, one for open-loop runs and one for
data assimilation runs.

## Scripts and other resources ##

- Script for calculating discharge for ParFlow along the river network
by Carina Furusho. Discharge can be calculated for each grid cell and
each time step:
<https://icg4geo.icg.kfa-juelich.de/ModelSystems/realpep_p4/blob/master/tools_postprocessing/functions_output_UNIVERSAL.py>


# CLM5-PDAF Starting instructions

## 1.  Install components

The following instructions assume you are working on the JSC machine JURECA. If you are working on another JSC machine switch the name appropriately and the instructions should still work. Running CLM5-PDAF on non JSC machines is not covered in these instructions but should be similar if TSMP-PDAF, PDAF, and CLM5 can be installed on the machine. 

### 1.1 TSMP-PDAF

Clone the TSMP-PDAF repository to your work folder and switch to the CLM5-PDAF branch like this:

```bash
 git clone https://github.com/HPSCTerrSys/TSMP.git tsmp

 cd tsmp

 git checkout TSMP_pdaf-clm5
```


### 1.2 PDAF
Get PDAF from [http://pdaf.awi.de/](http://pdaf.awi.de/) at least
version 1.13 (suggested version 2.0) and unpack it in the tsmp
directory and re-name the folder to `pdaf1_1`. 

Users with internal IBG-3 access may also look here:
[https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/pdaf2.0_fresh](https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/pdaf2.0_fresh)

Example: On your local machine:

```bash
 scp your_local_download_folder/PDAF-D_V1.16.tar.gz your_username@jureca.fz-juelich.de:/full_path_to_your_tsmp_folder
```

On the JSC machine in the tsmp folder:

```bash
 tar -xzf PDAF-D_V1.16.tar.gz
 mv PDAF-D_V1.16 pdaf1_1
```

### 1.3 CLM5

Clone CLM5 and its dependencies from the official repo to the tsmp folder like this:

```bash
git clone -b release-clm5.0 https://github.com/ESCOMP/CTSM.git clm5_0

cd clm5_0

source ../bldsva/machines/JURECA/loadenvs.Intel

./manage_externals/checkout_externals
```

- `source ...`: Loads the software environment on JURECA that is later
  used for build with Intel
- `./manage_externals/checkout_externals`: Loads dependencies of CLM5,
  mostly a collection of `git clone` commands
### 1.4 Compile the models together

The TSMP framework works with its own build system. With the models in
the folders as described before the compilation of both CLM5 and PDAF
can be done with a single build command after the following
preparations.

If not already configured for CLM5 the following environmental
variables are needed to compile CLM5 and can be for example be set
like this:

User defined data paths

Path to the root directory for all CESM / CLM input files
```bash
export CESMDATAROOT=$SCRATCH/cesm
```

This is the default value set by the build process. However, you can
export your own `CESMDATAROOT`.

Path to the CSM specific input files (usually subfolder of CESMDATAROOT)

```bash
export CSMDATA=$CESMDATAROOT/inputdata
```

In the `tsmp/bldsva` directory:

```bash
./build_tsmp.ksh -v 4.4.0MCTPDAF -c clm -m JURECA -O Intel
```

Potential problems that can happen during this step:

  * If you worked with CLM5 standalone and have created a `.cime`
    folder in your home directory it can cause a conflict during
    compilation. The error file will contain a line that says `ERROR:
    No machine jureca found` for example. The easiest solution to this
    is to move the `~/.cime` folder to `~/.cime_deactivated` while
    working with CLM5-PDAF.
  *  Using the newest software stages from JSC there is an import
    conflict with bigint for PERL.  See
    [https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf/tsmp/-/issues/67](https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf/tsmp/-/issues/67)
    for details. The error message will contain a line like this:
    `err=Can't locate bigint.pm in @INC`. The solution for now is to
    comment all lines with `bigint` and `bignum` in the script
    `tsmp/clm5_0/bld/config_files/clm_phys_vers.pm` since this is only
    needed to check some version numbers the PERL modules is not
    strictly needed.

Once successfully compiled the executable can be found in path
`tsmp/bin/JURECA_4.4.0MCTPDAF_clm/tsmp-pdaf`.

## 2. File preparation

For example files see [Example Case](#5-example-case).

### 2.1 Ensemble generation

For data assimilation with the ensemble Kalman filter we usually
create an ensemble by perturbing the atmospheric forcing files as well
as the sand and clay fractions in the surface file.

For CLM5 how these files are named and numbered can be customized, as
shown in the next section. But a recommended approach is to have a
suffix with zero filled numbering like
`surfdata_DE-Wue_hist_78pfts_CMIP6_simyr2000_c191001.nc_00001.nc`.

The exact details of perturbation can vary from study to study, but a
usable script to perturb the forcings and surface files can be found
at
`/p/largedata/jicg41/strebel2/CLM5PDAF_Example/perturb_forcings_and_soil_properties_in_range.py`.

In this script the following variables should be modified or at least considered for modification before running it for your case.

    * years (range of years for atm. forcing perturb)
    * num_ensembles (how many ensemble members should be created)
    * sname (path and naming for the output surface files)
    * sorig (path for the original surface file)
    * the attribute "perturbed by" (to indicate for future use who created the files)
    * the values of the random.uniform function (range of perturbation for soil characteristics)
    * sd, mean, and correl (the statistical characteristics for the atm. forcing perturbations)
    * fname (path and naming for the original atm. files)
    * outname (path and naming for the output atm. files)

### 2.2 Case generation

These instructions assume you already have a running CLM5 standalone
case before moving to CLM5PDAF. If not instructions of how to generate
the domain file, surface file, and download the default CLM5 inputfile
can be found on the GitLab.

We setup a new case in a new directory. In this directory we will need a few subfolders:

```bash
 mkdir logs

 mkdir timing/checkpoints -p
```
CLM5-PDAF does not use the same case setup as CLM5 standalone, instead we skip to the final namelist files that are needed to run CLM5. The minimal set of namelists to run CLM5 is:

    * atm_modelio.nml
    * datm_in
    * drv_flds_in
    * drv_in
    * esp_modelio.nml
    * glc_modelio.nml
    * ice_modelio.nml
    * lnd_in
    * lnd_modelio.nml
    * mosart_in
    * ocn_modelio.nml
    * rof_modelio.nml
    * wav_modelio.nml

Since these files are needed for each ensemble member we create them with a script. An example of  such a script can be found at `/p/largedata/jicg41/strebel2/CLM5PDAF_Example/create_ensemble_namelists.py`

In this script the following variables should be modified or at least considered for modification before running it for your case.

    *  domain file (appears in write_datm_in(), as fatmlndfrc in write_lnd_in(), in write_stream_files() with path and name separated)
    *  streams (in write_datm_in() modify to name of your case)
    *  finidat (in write_lnd_in() your initial condition file from spinup)
    *  fsurdat (in write_lnd_in() your surface files including ensemble numbering)
    *  hist_mfilt, hist_nhtfrq (customize your history file settings)
    *  filePath, fileNames (in write_stream_files customize the name and list of years for forcings)
    *  name of the case (at the end of write_stream_files() needs to be consistent with the ones given in write_datm_in())
    *  year_start, year_end (in main() give start and end year)
    *  prefix (in main() for naming of output files)
    *  num_ensemble (in main() number of ensemble members)
    *  b_dir (in main() path to build directory of clm within tsmp)
    *  r_dir (in main() path to the directory where the case will be run from i.e. the case directory)
    *  Every other variable uses the default from a standalone CLM5 case, but may vary in your case, check the namelists in the run directory of your case against the variables listed in this script if you have customized the CLM5 case.

The script at the moment can either be called without arguments if all
variables are modified within the script or with year_start and prefix
as arguments.

To run PDAF an additional namelist is needed, enkfpf.par for details
on the contents of this namelist refer to the TSMP-PDAF manual.

Lastly, we create a symbolic link to the executable in the case
directory. This way we do not create unnecessary copies of the
executable for each case, we make sure that all cases are run with the
same executable, and modification and re-compilations of the
executable will not need any changes in the case directories. The
symbolic link can be created like this:

```bash
 ln -s PATH_TO_THE_EXECUTABLE/tsmp-pdaf .
```

Similarly, we can use a symbolic link to the loadenvs file from TSMP to make sure we use the same module versions during runtime as during compilation.

```bash
 ln -s PATH_TO_THE_TSMP_FOLDER/bldsva/machines/JURECA/loadenvs.Intel loadenvs
```

## 3. Running the simulation

To run the simulation we use can use a simple jobscript, an example can be found at `/p/largedata/jicg41/strebel2/CLM5PDAF_Example/jobscript_da.slurm`

The arguments for the tsmp-pdaf executable are explained in the TSMP-PDAF manual.

Important parameters are `-n_modeltasks` is the number of ensemble members and `nodes` times `ntasks-per-node` should be at least the same or a multiple of `n_modeltasks`. 

The argument `-delt_obs` can be used to switch between open loop and DA simulation by increasing the value to be larger than the number of simulation steps set in `enkfpar.par` no assimilation will take place and the simulation will be an open loop ensemble simulation. 

Additionally, `-obs_filename` has the path to the observations and the prefix of the observation files without the numbering.


## 4. Post-processing the output

CLM5-PDAF will create history files for each ensemble members
according to the settings defined in the script
`create_ensemble_namelists.py` with the hist variables. For large
ensembles this will create a lot of files and should be processed and
moved / archived away from `$SCRATCH` to avoid filling the number of
files limit. One way to process the output is to collect the ensemble
statistics of the variables in one file. An example script to perform
this post-processing can be found at
`/p/largedata/jicg41/strebel2/CLM5PDAF_Example/collect_some_ensemble_stats.py`

The script should be modified before use, by changing the
`collect_vars` variable to a list with the CLM5 variable names that
should be collected. The output netcdf file will be called
`Ensemble_Collection_` and contain the listed variables with their
ensemble minimum, maximum, mean, and standard deviation.

At the moment the script assumes that each year is in a separate
history file, but could be customized to work on multi-year or less
than a year files.  Currently, the script takes the arguments: year,
OL or DA, num_ensemble, case name. Where case name refers to the
previously mentioned prefix that can be used to differentiate between
different simulations within the same case.

## 5. Example case:

All the input files necessary to create a 10 ensemble case for the
year 2009 are located at
`/p/largedata/jicg41/strebel2/CLM5PDAF_Example/Example_Case` and can
be used to test the setup before working on your own case.




