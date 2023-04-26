# Getting Started

**Throughout this quick start guide a distinction is made between HPSC-TerrSys users who use specific HPC systems at HPC centres in Germany and external users. Through centralized machinefiles and compute environment initialisation scripts, TSMP can be very quickly adjusted to any other HPC site.**

For a general documentation on how to use TSMP, please refer to the more detailed [Building TSMP](./build_tsmp/build_examples.md#building-tsmp) section.

## General concept

TSMP uses a multiple program multiple data (MPMD) execution model where states (e.g., near surface air temperature) and fluxes (e.g., latent and sensible heat fluxes) are passed between the atmospheric, land surface and hydrologic/groundwater component models along compartmental interfaces through the OASIS3 or the parallel OASIS3-MCT coupling interface.

TSMP essentially consists of the interface that combines the specific versions of different component models in a modular (various combinations of component models) and flexible (various spatial and temporal resolutions) way. During the install process, before the compilation of the component models, the source code of the component models is patched by the setup shell scripts, implementing the necessary coupling functionalities.

Please note: For licensing and maintenance reasons, we do not provide any source code of any of the component models or the coupler library. The users should download the component models from their respective websites. TSMP has various setup options: each component model standalone (ParFlow, CLM, COSMO), only land surface and hydrology / groundwater (CLM+ParFlow), only land surface and atmosphere (CLM+COSMO), or fully coupled (ParFlow, CLM, COSMO). Since the patches are designed for specific versions of the models, it is required to get the exact original, unaltered version as specified.

During the build and compile step ([step 5 below](#step-5-build-tsmp-interface-and-component-models)) multiple predefined coupled model combinations (of component models and their versions) as well as predefined HPC systems and operating environments can be chosen. Via so-called machinefiles the built and execution environment can be centrally adjusted to the respective HPC systems. A generic machinefile for x86-64 GNU/Linux systems with a GNU Compiler Collection and a standard scientific software stack is also provided. Many predefined built environments are provided for the machines which are available to the HPSC-TerrSys users.

During the setup and configuration step ([step 6 below](#step-6-setup-and-configuration-of-the-respective-usage-and-test-case)), TSMP is configured for different, predefined test cases or numerical experiments, for which we also provide input data via a dedicated data server. Here we document one specific experiment. The [TSMP user guide](https://hpscterrsys.github.io/TSMP/index.html) will contain a complete list.

The two most important scripts to build and setup TSMP are:

* `build_tsmp.ksh` to built and compile TSMP (one component model after the other)
* `setup_tsmp.ksh` to setup one of the predefined experiments

To get a man-page for the usage of this scripts, do:
```shell
./build_tsmp.ksh --man
./setup_tsmp.ksh --man
```

## The fully coupled pan-European EURO-CORDEX evaluation experiment with TSMP

This test case uses the current TSMP release version v1.4.0 with OASIS3-MCT, COSMO v5.01, CLM v3.5 and ParFlow 3.12.  A short 3hr simulation in a climate-mode configuration over Europe is set up, driven by ERA-Interim reanalysis, following the [EURO-CORDEX project](https://euro-cordex.net/) experiment guidelines. Simulated time span: 2016-05-01_12:00:00 to 2016-05-01_15:00:00.

### Step 1: Dependencies
For the users who use Jülich Supercomputing Centre facilities JUWELS and JURECA-DC, all necessary software modules are loaded automatically through a "loadenv" file located in directory JUWELS or JURECA-DC in machines directory. The users of other HPC systems should provide an appropriate "loadenv" files for loading the modules and locate it in `machines/<machine_name>`, similar to JURECA-DC and JUWELS. For the users who want to port TSMP on GENERIC_X86 Linux platform, a script is provided by TSMP team which installs the following libraries automatically and create a "loadenv" file in the directory `machines/GENERIC_X86`. For more information on using this script please see the README in branch **TSMP_x86_64**.

* gfortran
* gcc
* g++
* ksh
* bash
* zlib
* curl
* make
* python
* OpenMPI
* netCDF
* HDF5
* GRIBAPI
* TCL
* Hypre
* Silo
* Lapack

A short guide on how the TSMP built system can be expanded to account for your local situation (software modules, compiler, MPI wrapper) will be provided shortly. The best starting point is to use the generic GNU compiler built. A hint on how to proceed is given in [step 5 below](#step-5-build-tsmp-interface-and-component-models).
The following libraries are required to run the EURO-CORDEX experiment:

#### Porting TSMP on GENERIC_X86 Linux platform

The users who want to port TSMP on GENERIC_X86 Linux, the TSMP team provided a script to install all the necessary libraries as mentioned in [step 1 above](#step-1-dependencies) automatically in TSMP root directory. Please run the script "lib_install.sh" located in `bldsva` directory to install the libraries. Note that if you exported already one of the libraries Netcdf, HDF5, GRIBAPI, Silo, Hypre and TCL in the .baschrc or .profile, you need to comment them out in order to not mess up the installation via the script lib_install.sh.

### Step 2: Get the TSMP interface

Go to your preferred root directory (e.g., on a scratch file system) for the TSMP installation and get the `HEAD` of the `master` branch from `github` and set the environment variable `TSMP_DIR` to TSMP installation directory; for bash:

```shell
git clone https://github.com/HPSCTerrSys/TSMP.git
export TSMP_DIR=$(realpath TSMP)
```

### Step 3: Get the component models for this experiment

```shell
cd $TSMP_DIR
```

#### HPSC-TerrSys users

Authenticate with your GitLab web GUI user name and password and clone
the repositories (instead of "fresh", also "legacy" repositories with
specific code modifications may be retrieved):

```shell
git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/cosmo5.01_fresh.git  cosmo5_1
git clone -b v3.12.0 https://github.com/parflow/parflow.git                              parflow
git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/clm3.5_fresh.git     clm3_5
git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/oasis3-mct.git       oasis3-mct
```

#### External users

Make sure directories are renamed as shown above. The component models can be downloaded from the respective websites as indicated below. Due to licensing reasons HPSC-TerrSys cannot provide these source codes to external users.

##### COSMO v5.01

Available from http://www.cosmo-model.org. A license agreement is needed.

##### CLM v3.5

Available from https://github.com/HPSCTerrSys/CLM3.5/tree/clm3.5_rel3.code.c070524.

##### ParFlow v3.12
ParFlow 3.12  is available from
```shell
git clone -b v3.12.0 https://github.com/parflow/parflow.git
```

##### OASIS3-MCT v2.0

Available from https://oasis.cerfacs.fr/en/downloads/.

### Step 4: Retrieving the test case input data (EURO-CORDEX)

For each officially documented and supported test case experiment, HPSC-TerrSys provides all required, pre-processed input data (as well as reference simulation results in the future).

The input files which are necessary for running the `EURO-CORDEX evaluation run` test case:

```shell
cd $TSMP_DIR/bldsva
./download_data_for_test_cases.ksh cordex
```

### Step 5: Build TSMP, interface and component models

Before building TerrSysMP, first check what build options are there

```shell
cd $TSMP_DIR/bldsva
./build_tsmp.ksh -a
```

Building the fully coupled TSMP with ParFlow (pfl), the Community Land Model (clm) and the COSMO NWP and regional climate model (cos); this is a built on the JUWELS HPC system of Jülich Supercomputing Centre using Intel compilers and ParaStation MPI:

```shell
cd $TSMP_DIR/bldsva
./build_tsmp.ksh -c clm3-cos5-pfl -m JUWELS -O Intel
```

A note to external users:

The path to the modules, compiler and MPI wrapper can be set in:  
`$TSMP_DIR/bldsva/machines/config_<your_machine_name>.ksh`

An example for the JUWELS HPC system at the Jülich Supercomputing Centre is here: `$TSMP_DIR/bldsva/machines/config_JUWELS.ksh` \
and an example to load the environment modules or the needed software compatible with Intel compiler version 2023 is in `$TSMP_DIR/bldsva/machines/JUWELS/loadenvs.Intel`.

#### Build TSMP on GENERIC_X86 Linux platform

After successful installation of all the necessary libraries as mentioned in [step 1 above](#step-1-dependencies), TSMP can be built using: 

```shell
cd $TSMP_DIR/bldsva
./build_tsmp.ksh -c clm3-cos5-pfl -m GENERIC_X86 -O Gnu
```

### Step 6: Setup and configuration of the respective usage and test case

For HPSC-TerrSys users:

To configure TSMP for the [EURO-CORDEX test case experiment](#the-fully-coupled-pan-european-euro-cordex-evaluation-experiment-with-tsmp) on JUWELS machine (on JURECA-DC just change -m JUWELS to -m JURECA):

```shell
cd $TSMP_DIR/bldsva
./setup_tsmp.ksh -c clm3-cos5-pfl  -V cordex -m JUWELS -I _cordex -O Intel
```

This includes the creation of a run directory, the copying of namelists, the provisioning of run control scripts for the job scheduler, incl. mapping files which pin the MPI tasks of the component model to specific CPU cores, as well as copying and linking of forcing data.

#### Setup and configuration of the test case on GENERIC_X86 Linux platform

Now that the all the necessary libraries are installed as mentioned in [step 1 above](#step-1-dependencies) and TSMP is also build ([step 5 above](#step-5-build-tsmp-interface-and-component-models)), we can configure TSMP on GENERIC_X86 Linux for the the [EURO-CORDEX test case experiment](#the-fully-coupled-pan-european-euro-cordex-evaluation-experiment-with-tsmp) as follows: 

```shell
cd $TSMP_DIR/bldsva
./setup_tsmp.ksh -c clm3-cos5-pfl -V cordex -m GENERIC_X86 -I _cordex -O Gnu
```

### Step 7: Run the test case

For HPSC-TerrSys users:

Change into the run directory:

```shell
cd $TSMP_DIR/run/JUWELS_clm3-cos5-pfl _cordex_cordex
```

Edit the scheduler settings (`#SBATCH` lines):

```shell
vim tsmp_slm_run.bsh
```

Submit the job:

```shell
sbatch tsmp_slm_run.bsh
```

Monitor the running job by checking some COSMO and ParFlow logs:

```shell
tail -f YUPRMASS
tail -f cordex0.11.out.kinsol.log
```

### Step 8: Simulation results

After a successful model run you will have model outputs form ParFlow, CLM and COSMO in the run directory:

```shell
cd $TSMP_DIR/run/JUWELS_clm3-cos5-pfl cordex_cordex
ls -l *satur*.pfb *press*.pfb
ls -l clmoas*.nc
ls -l cosmo_out/*.nc
```

