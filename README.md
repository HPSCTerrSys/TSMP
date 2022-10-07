# Table of contents
1. [Introduction](#introduction)
	1. [TSMP](#TSMP)
	2. [Citing TSMP](#Citing)
	3. [Quick Start on Linux](#Quick-Start)
	4. [General concept](#General-concept)
	5. [TSMP version history](#ver_his)
2. [The fully coupled pan-European EURO-CORDEX evaluation experiment with TSMP](#ref_exp)
    1. [Step 1: Dependencies](#ref_step1)
    2. [Step 2: Get the TSMP interface](#ref_step2)
    3. [Step 3: Get the component models for this experiment](#ref_step3)
		1. [HPSC-TerrSys users](#HPSC-TerrSys-users)
		2. [External users](#External-users)
    4. [Step 4: Retrieving the test case input data (NRW and EURO-CORDEX)](#ref_step4)
    5. [Step 5: Build TSMP, interface and component models](#ref_step5)
    6. [Step 6: Setup and configuration of the respective usage and test case](#ref_step6)
    7. [Step 7: Run the test case](#ref_step7)
    8. [Step 8: Simulation results](#ref_step8)
3. [For ICON Users (HPSC-TerrSys users)](#foricon)	
4. [Heterogeneous Job using TSMP](#hetero_job)
5. [NRW Test case](#nrw_test)
6. [IdealRTD (Idealized) Test case](#idealrtd_test)
7. [Patching the orginal source code](#patching)
8. [Long time climate simulation](#climate_sim)
9. [Restart functionality for EURO-CORDEX test case](#restart)
10. [Automatic Porting of TSMP on x86 machines](#TSMP_x86)
11. [To come](#To-come)
12. [Documentation](#ref_doc)

# Introduction <a name="introduction"></a>

## TSMP <a name="TSMP"></a>

The Terrestrial System Modeling Platform (TSMP or TerrSysMP, https://www.terrsysmp.org) is an open source scale-consistent, highly modular, massively parallel regional Earth system model. TSMP essentially consists of an interface which couples dedicated versions of the Consortium for Small-scale Modeling (COSMO, http://www.cosmo-model.org) atmospheric model in NWP or climate mode, the Community Land Model (CLM, http://www.cesm.ucar.edu/models/clm/), and the hydrologic model ParFlow (https://www.parflow.org) through the OASIS3-MCT coupler (https://portal.enes.org/oasis, https://www.mcs.anl.gov/research/projects/mct/).

TSMP allows for a physically-based representation of transport processes of mass, energy and momentum and interactions between the different compartments of the geo-ecosystem across scales, explicitly reproducing feedbacks in the hydrological cycle from the groundwater into the atmosphere.

TSMP is extensively used for idealized and real data process and sensitivity studies in water cycle research, for climate change simulations, data assimilation studies including reanalyses, as well as experimental real time forecasting and monitoring simulations, ranging from individual catchments to continental model domains. TSMP runs on notebooks as well on latest supercomputers using a range of compilers.

TSMP development has been driven by groups within the [Center for High-Performance Scientific Computing in Terrestrial Systems](http://www.hpsc-terrsys.de) (HPSC-TerrSys), as part of the [Geoverbund ABC/J](http://www.geoverbund-abcj.de/geoverbund/EN/Home/home_node.html), the geoscientific network of the University of Cologne, Bonn University, RWTH Aachen University, and the Research Centre Jülich. The current team is anchored in Jülich and Bonn in Germany.

**Visit**

**https://www.terrsysmp.org**

**for information on the features of TSMP, ongoing developments, citation, usage examples, links to documentation, the team, contact information and publications.**

## Citing TSMP <a name="TSMP"></a>

If you use TSMP in a publication, please cite the these papers that describe the model's basic functionalities:

* Shrestha, P., Sulis, M., Masbou, M., Kollet, S., and Simmer, C. (2014). A Scale-Consistent Terrestrial Systems Modeling Platform Based on COSMO, CLM, and ParFlow. Monthly Weather Review, 142(9), 3466–3483. doi:[10.1175/MWR-D-14-00029.1](https://dx.doi.org/10.1175/MWR-D-14-00029.1).
* Gasper, F., Goergen, K., Kollet, S., Shrestha, P., Sulis, M., Rihani, J., and Geimer, M. (2014). Implementation and scaling of the fully coupled Terrestrial Systems Modeling Platform (TerrSysMP) in a massively parallel supercomputing environment &ndash; a case study on JUQUEEN (IBM Blue Gene/Q). Geoscientific Model Development, 7(5), 2531-2543. doi:[10.5194/gmd-7-2531-2014](https://dx.doi.org/10.5194/gmd-7-2531-2014).

## Quick Start on Linux <a name="Quick-Start"></a>

This very short user guide only covers one variant on how the model can be setup and configured for *one* specific experiment, which we use as one of the default test cases. To get an overview on possible TSMP applications refer to the [TSMP website](https://www.terrsysmp.org) and the [user guide](#ref_doc).

## General concept <a name="General-concept"></a>

TSMP uses a multiple program multiple data (MPMD) execution model where states (e.g., near surface air temperature) and fluxes (e.g., latent and sensible heat fluxes) are passed between the atmospheric, land surface and hydrologic/groundwater component models along compartmental interfaces through the OASIS3 or the parallel OASIS3-MCT coupling interface.

TSMP essentially consists of the interface that combines the specific versions of different component models in a modular (various combinations of component models) and flexible (various spatial and temporal resolutions) way. During the install process, before the compilation of the component models, the source code of the component models is patched by the setup shell scripts, implementing the necessary coupling functionalities.

Please note: For licensing and maintenance reasons, we do not provide any source code of any of the component models or the coupler library. The users should download the component models from their respective websites. TSMP has various setup options: each component model standalone (ParFlow, CLM, COSMO), only land surface and hydrology / groundwater (CLM+ParFlow), only land surface and atmosphere (CLM+COSMO), or fully coupled (ParFlow, CLM, COSMO). Since the patches are designed for specific versions of the models, it is required to get the exact original, unaltered version as specified.

During the build and compile step ([step 5 below](#ref_step5)) multiple predefined coupled model combinations (of component models and their versions) as well as predefined HPC systems and operating environments can be chosen. Via so-called machinefiles the built and execution environment can be centrally adjusted to the respective HPC systems. A generic machinefile for x86-64 GNU/Linux systems with a GNU Compiler Collection and a standard scientific software stack is also provided. Many predefined built environments are provided for the machines which are available to the HPSC-TerrSys users.

During the setup and configuration step ([step 6 below](#ref_step6)), TSMP is configured for different, predefined test cases or numerical experiments, for which we also provide input data via a dedicated data server. Here we document one specific experiment. The [TSMP user guide](#ref_doc) will contain a complete list.

**Throughout this quick start guide a destinction is made between HPSC-TerrSys users who use specific HPC systems at HPC centres in Germany and external users. Through centralized machinefiles and compute environment initialisation scripts, TSMP can be very quickly adjusted to any other HPC site.**

The two most important scripts to build and setup TSMP are:

* `build_tsmp.ksh` to built and compile TSMP (one component model after the other)
* `setup_tsmp.ksh` to setup one of the predefined experiments

To get a man-page for the usage of this scripts, do:
```shell
   ./build_tsmp.ksh --man
   ./setup_tsmp.ksh --man
```
## TSMP version history <a name="ver_his"></a> 
The model components used in TSMP are OASIS3-MCT v2, COSMO v5.01, CLM v3.5, ParFlow 3.2 for TSMP versions v1.2.1, v1.2.2 and v1.2.3, and ParFlow 3.9 for version v1.3.3. TSMP supports ParFlow 3.7 onwards from version v1.3.3 onward.

Those who need to work with ParFlow 3.2, should use the branch `TSMP_pdaf`.

# The fully coupled pan-European EURO-CORDEX evaluation experiment with TSMP <a name="ref_exp"></a>

This test case uses the current TSMP release version v1.2.1 with OASIS3-MCT, COSMO v5.01, CLM v3.5 and ParFlow 3.2 (ParFlow 3.9 from TSMP version v1.3.3 onward).  A short 3hr simulation in a climate-mode configuration over Europe is set up, driven by ERA-Interim reanalysis, following the [EURO-CORDEX project](https://euro-cordex.net/) experiment guidelines. Simulated time span: 2016-05-01_12:00:00 to 2016-05-01_15:00:00.

### Step 1: Dependencies <a name="ref_step1"></a>
For the users who use Jülich Supercomputing Centre facilities JUWELS and JURECA, all necessary software modules are loaded automatically through a "loadenv" file located in directory JUWELS or JURECA in machines directory. The users of other HPC systems should provide an appropriate "loadenv" files for loading the modules and locate it in `machines/<machine_name>`, similar to JURECA and JUWELS. For the users who want to port TSMP on GENERIC_X86 Linux platform, a script is provided by TSMP team which installs the following libraries automatically and create a "loadenv" file in the directory `machines/GENERIC_X86`. For more information on using this script please see the README in branch **TSMP_x86_64**.

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

A short guide on how the TSMP built system can be expanded to account for your local situation (software modules, compiler, MPI wrapper) will be provided shortly. The best starting point is to use the generic GNU compiler built. A hint on how to proceed is given in [step 5 below](#ref_step5).
The following libraries are required to run the EURO-CORDEX experiment:

#### <a name="ref_step11"></a> Porting TSMP on GENERIC_X86 Linux platform

The users who want to port TSMP on GENERIC_X86 Linux, the TSMP team provided a script to install all the necessary libraries as mentioned in [step 1 above](#ref_step1) automatically in TSMP root directory. Please run the script "lib_install.sh" located in `bldsva` directory to install the libraries. Note that if you exported already one of the libraries Netcdf, HDF5, GRIBAPI, Silo, Hypre and TCL in the .baschrc or .profile, you need to comment them out in order to not mess up the installation via the script lib_install.sh.

### Step 2: Get the TSMP interface <a name="ref_step2"></a>

Go to your preferred root directory (e.g., on a scratch file system) for the TSMP installation and get the `HEAD` of the `master` branch from `github` and set the environment variable `TSMP_DIR` to TSMP installation directory; for bash:

```shell
   git clone https://github.com/HPSCTerrSys/TSMP.git
   cd TSMP
   export TSMP_DIR=$(pwd)
```

### Step 3: Get the component models for this experiment <a name="ref_step3"></a>

```shell
   cd $TSMP_DIR
```

#### HPSC-TerrSys users <a name="HPSC-TerrSys-users"></a>

Authenticate with your GitLab web GUI user name and password and clone
the repositories (instead of "fresh", also "legacy" repositories with
specific code modifications may be retrieved):

```shell
   git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/cosmo5.01_fresh.git  cosmo5_1
   git clone -b v3.9.0 https://github.com/parflow/parflow.git                              parflow
   git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/clm3.5_fresh.git     clm3_5
   git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/oasis3-mct.git       oasis3-mct
```

#### External users <a name="External-users"></a>

Make sure directories are renamed as shown above. The component models can be downloaded from the respective websites as indicated below. Due to licensing reasons HPSC-TerrSys cannot provide these source codes to external users.

##### COSMO v5.01

Available from http://www.cosmo-model.org. A license agreement is needed.

##### CLM v3.5

Available from http://www.cgd.ucar.edu/tss/clm/distribution/clm3.5/.

##### ParFlow v3.2

Available from https://github.com/parflow/.

```shell
   git clone --branch v3.2.0 https://github.com/parflow/parflow.git parflow3_2
```
##### ParFlow v3.9
ParFlow 3.9  is available from
```shell
git clone -b v3.9.0 https://github.com/parflow/parflow.git
```

##### OASIS3-MCT v2.0

Available from https://verc.enes.org/oasis/.

### Step 4: Retrieving the test case input data (NRW and EURO-CORDEX) <a name="ref_step4"></a>

For each officially documented and supported test case experiment, HPSC-TerrSys provides all required, pre-processed input data (as well as reference simulation results in the future).

The input files which are necessary for running the `EURO-CORDEX evaluation run` test case:

```shell
   cd $TSMP_DIR/bldsva
   ./download_data_for_test_cases.ksh cordex
```

and for `NRW` test case:

```shell
   cd $TSMP_DIR/bldsva
   ./download_data_for_test_cases.ksh nrw
```

### Step 5: Build TSMP, interface and component models <a name="ref_step5"></a>

Before building TerrSysMP, first check what build options are there

```shell
   cd $TSMP_DIR/bldsva
   ./build_tsmp.ksh -a
```

Building the fully coupled TSMP with ParFlow (pfl), the Community Land Model (clm) and the COSMO NWP and regional climate model (cos); this is a built on the JURECA HPC system of Jülich Supercomputing Centre using Intel compilers and ParaStation MPI:

```shell
   cd $TSMP_DIR/bldsva
   ./build_tsmp.ksh -v 3.1.0MCT -c clm-cos-pfl -m JUWELS -O Intel
```
For building ParFlow 3.9 on GPU:

```shell
   cd $TSMP_DIR/bldsva
   ./build_tsmp.ksh -v 3.1.0MCT -c clm-cos-pfl -m JUWELS -O Intel -A GPU
```


A note to external users:

The path to the modules, compiler and MPI wrapper can be set in:  
`$TSMP_DIR/bldsva/machines/<your_machine_name>/build_interface_<your_machine_name>.ksh`

An example for the JUWELS HPC system at the Jülich Supercomputing Centre is here: `$TSMP_DIR/bldsva/machines/JUWELS/build_interface_JURECA.ksh` \
and an example to load the environment modules or the needed software compatible with Intel compiler version 2020 is in `$TSMP_DIR/bldsva/machines/JUWELS/loadenvs.Intel`.

#### Build TSMP on GENERIC_X86 Linux platform

After successful installation of all the necessary libraries as mentioned in [step 1 above](#ref_step1), TSMP can be built using: 

```shell
   cd $TSMP_DIR/bldsva
   ./build_tsmp.ksh -v 3.1.0MCT -c clm-cos-pfl -m GENERIC_X86 -O Gnu
```

### Step 6: Setup and configuration of the respective usage and test case <a name="ref_step6"></a>

For HPSC-TerrSys users:

To configure TSMP for the [EURO-CORDEX test case experiment](#ref_exp) on JUWELS machine (on JURECA just change -m JUWELS to -m JURECA):

```shell
   cd $TSMP_DIR/bldsva
   ./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m JUWELS -I _cordex -O Intel
```

For configuring TSMP for a heterogeneous job:

```shell
   cd $TSMP_DIR/bldsva
   ./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m JUWELS -I _cordex -O Intel -A GPU
```
In this heterogeneous job ParFlow 3.9 will run on GPU while Cosmo and CLM on CPU.

This includes the creation of a run directory, the copying of namelists, the provisioning of run control scripts for the job scheduler, incl. mapping files which pin the MPI tasks of the component model to specific CPU cores, as well as copying and linking of forcing data.

#### Setup and configuration of the test case on GENERIC_X86 Linux platform

Now that the all the necessary libraries are installed as mentioned in [step 1 above](#ref_step1) and TSMP is also build ([step 5 above](#ref_step4)), we can configure TSMP on GENERIC_X86 Linux for the the [EURO-CORDEX test case experiment](#ref_exp) as follows: 

```shell
   cd $TSMP_DIR/bldsva
   ./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m GENERIC_X86 -I _cordex -O Gnu
```

### Step 7: Run the test case <a name="ref_step7"></a>

For HPSC-TerrSys users:

Change into the run directory:

```shell
   cd $TSMP_DIR/run/JURECA_3.1.0MCT_clm-cos-pfl_cordex_cordex
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

### Step 8: Simulation results <a name="ref_step8"></a>

After a successful model run you will have model outputs form ParFlow, CLM and COSMO in the run directory:

```shell
   cd $TSMP_DIR/run/JURECA_3.1.0MCT_clm-cos-pfl_cordex_cordex
   ls -l *satur*.pfb *press*.pfb
   ls -l clmoas*.nc
   ls -l cosmo_out/*.nc
```
# For ICON Users (HPSC-TerrSys users) <a name="foricon"></a>

### Step 0: JUDAC Storage Managment

Please inform yourself about the storage managment. If you are new to the Jülich supercomputers, it is recomannded to visit "Introduction to the usage and programming of supercomputer resources in Jülich" course. The next date of this ocurse can be found here: https://www.fz-juelich.de/ias/jsc/EN/Expertise/Workshops/Courses/courses_node.html . Please, make yourself familiar with the different partitions of the JUDAC system as well:  
> https://www.fz-juelich.de/ias/jsc/EN/Expertise/Datamanagement/JUDAC/FAQ/judac-FAQ_node.html

It is recomannded to store no data at $HOME, your model data at $PROJECT and your experiment data at $SCRATCH.

### Step 1: Dependencies
See [step 1 above](#ref_step1)

### Step 2: Get the TSMP interface

Go to your preferred root directory (e.g., your $PROJECT directory) for the TSMP installation and get the `TSMP_iconcoup` branch from `gitlab` and set the environment variable `TSMP_DIR` to TSMP installation directory; for bash:

```shell
   git clone https://icg4geo.icg.kfa-juelich.de/spoll/tsmp.git
   cd tsmp
   export TSMP_DIR=$(pwd)
   git checkout TSMP_iconcoup
```

### Step 3: Get the component models for this experiment

Authenticate with your GitLab web GUI user name and password and clone the repositories (instead of "fresh", also "legacy" repositories with specific code modifications may be retrieved). Choose the model components needed from your experiment, in case of any coupled model one need the external coupler "oasis3-mct".

```shell
   git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/icon2.1_legacy.git   icon2-1
   git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/parflow3.2_fresh.git parflow3_2
   git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/clm3.5_fresh.git     clm3_5
   git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/oasis3-mct.git       oasis3-mct
```

### Step 4: Build TSMP, interface and component models

Before building TerrSysMP, first check what build options are there

```shell
   cd $TSMP_DIR/bldsva
   ./build_tsmp.ksh -a
```

Building the fully coupled TSMP with ParFlow (pfl), the Community Land Model (clm) and the ICON NWP and regional climate model (icon); this is a built on the JUWELS HPC system of Jülich Supercomputing Centre using Intel compilers and ParaStation MPI:

```shell
   cd $TSMP_DIR/bldsva
   ./build_tsmp.ksh -v 4.1.0MCT -c clm-icon-pfl -m JUWELS -O Intel
```

For ICON standalone version use:
```shell
   cd $TSMP_DIR/bldsva
   ./build_tsmp.ksh -v 4.1.0MCT -c icon -m JUWELS -O Intel
```

### Step 5: Setup and configuration of the respective usage and test case

It is recommended to store your model data at scratch (please notice that data will be automatically deleted if there are not touched for a while). Please, edit *YOUR_PROJECT* and *YOUR_DIRECTORY* correspondingly.

```shell
   export EXP=/p/scratch/YOUR_PROJECT/YOUR_DIRECTORY
```
To configure TSMP for the Germany test case for ICON stand-alone on JUWELS machine:

```shell
   cd $TSMP_DIR/bldsva
   ./setup_tsmp.ksh -v 4.1.0MCT -V germany -m JUWELS -c icon -r $EXP/JUWELS_icon_4.1.0MCT_germany -I _experiment -O Intel
```

This includes the creation of a run directory, the copying of namelists, the provisioning of run control scripts for the job scheduler, incl. mapping files which pin the MPI tasks of the component model to specific CPU cores, as well as copying and linking of forcing data.

### Step 6: Run the test case

Change into the run directory:

```shell
   cd $EXP/JUWELS_icon_4.1.0MCT_germany_experiment
```

Edit the scheduler settings (`#SBATCH` lines) accordingly to *YOUR_PROJECT* and run time of your experiment. Please do not change *nodes*, *ntasks* and *ntasks-per-node*:

```shell
   vi tsmp_slm_run.bsh
```

Proof if all necessary files are in your experiment directory:

```shell
   ls -lrth
```

Submit the job:

```shell
   sbatch tsmp_slm_run.bsh
```

Monitor your ICON simulation by

```shell
    tail -f mpiMPMD-err.*
```

### Step 7: Create own case (optional)

> This step is for advanced users.

In case you want create your own setup you can create your own direcory (*YOUR_CASE*) in the `TSMP_DIR`/setups and add the name of your case in the `setupsAvail` of your maschine (e.g. JUWELS) in the `supported_versions.ksh` file.

```shell
    mkdir $TSMP_DIR/bldsva/setups/YOUR_CASE
    vi $TSMP_DIR/bldsva/supported_versions.ksh
```
Each setup needs a setup ksh-script with the naming convention *'YOUR_CASE'_'YOUR_MACHINE'_setup.ksh* (e.g. `germany_JUWELS_setup.ksh` for the germany test case on the JUWELS machine). Please, adjust the setup carefully, you can use the setups scripts from various test cases as blueprint. You can add the namelist of your component models in the direcotry of your case (e.g. icon_master.namelist and NAMELIST_icon for ICON).  

# Heterogeneous Job using TSMP <a name="hetero_job"></a>

TSMP has the possibility of submitting heterogeneous job for EURO-CORDEX test experiment, meaning that ParFlow3.9  will run on GPU while Cosmo5.1 and CLM3.5 on CPU.
After locating the model components in the TSMP root as mentioned in [Step 3](#ref_step3), the following command should be executed in order to build TSMP and to create the run directory for the [EURO-CORDEX test case experiment](#ref_exp) on JUWELS machine (on JURECA just change -m JUWELS to -m JURECA):

Building TSMP for HPSC-TerrSys users:
```shell
   cd $TSMP_DIR/bldsva
   ./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m JUWELS -I _cordex -O Intel -A GPU
```

Creating run directory for HPSC-TerrSys users:

```shell
   cd $TSMP_DIR/bldsva
   ./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m JUWELS -I _cordex -O Intel -A GPU
```

Not that heterogeneous job will be possible only for TSMP with ParFlow3.9 which is available from version v1.3.3.

# NRW Test case <a name="nrw_test"></a>
NRW Test case covers a geographical domain of 150 km x 150 km encompassing the North Rhine-Westphalia region, located in western Germany, Belgium, the Netherlands, and Luxembourg. The experiment is carried out in a clear sky day condition (08 May 2008) using the fully coupled (COSMO5.01-CLM3.5-
ParFlow) configuration of TSMP. The atmospheric component uses a constant lateral spatial resolution of about 1 km and a variable vertical discretization into 50 levels gradually coarsening from the bottom (20 m) to the top (22000 m). Initial and lateral boundary conditions for the atmospheric model are obtained from the operational weather forecast model
COSMO-DE of the German Weather Service (DWD).For downloading the necessary INPUT data for the NRW test case please see [Step 3](#ref_step3) and for building the TSMP (with -v 3.1.0MCT) refer to [Step 5](#ref_step5). For more information about NRW test case please refer to https://doi.org/10.3390/w10111697. \
To configure TSMP for the NRW test case on JUWELS machine (on JURECA just change -m JUWELS to -m JURECA):

```shell
   cd $TSMP_DIR/bldsva
 ./setup_tsmp.ksh -v 3.1.0MCT -V nrw -m JUWELS -I _ClearSkyDay  -s 2008-05-08_00 -S 2008-05-08_00 -T 24 -O Intel
```

# IdealRTD (Idealized) Test case <a name="idealrtd_test"></a>
We will use a horizontally homogeneous domain (Figure 1) to simulate the diurnal
cycle of the atmospheric boundary layer with radiative forcing. A periodic boundary
condition will be used in X and Y direction for COSMO, while the X- and Y-slopes are
set to zero for ParFlow, so there is no lateral flow. The domain setup is summarized in the Table below.


|        |  NX    |  NY  | NZ | dt(s) |
| -----  | :----- | ----:|----:|----:|
| cosmo  | 60    | 30  | 50 | 18 |
| CLM    | 54    | 24  | 10 | 18 |
| ParFlow| 54    | 24  | 30 | 18 |

NX_clm = NX_pfl =NX_cos -2 *nboundlines

The initial soil moisture content for the different soil textures (e.g., sandy loam, clay
loam) is varied based on the ($\psi - \theta_vol$) van Genuchten relationship, to
simulate the physical processes from water stressed environment to an atmospheric
controlled environment. The groundwater table is specified at a depth of 5 m below
the surface, with a spatially homogeneous and constant unsaturated zone above. The atmosphere is initialized with a semi-idealized sounding data obtained from
Stuttgart (observations) at 00Z 07 Aug 2015. Ground and vegetation temperature are
initialized horizontally homogeneous with a value of 287 K. The simulation will be
integrated for 24 hours starting at midnight with an hourly output frequency.

## Preprocessing of Input Data <a name="Preprocessing_idealrtd"></a>
Here the user can generate the case-specific input data for ParFlow
and CLM. Note that the COSMO model uses the same setup for all the ensemble
members, that is, no need for ad-hoc editing. 
Please download the Input data and preproccesing sript using:

```shell
   cd $TSMP_DIR/bldsva
   ./download_data_for_test_cases.ksh idealrtd
```

then 

```shell
   cd $TSMP_DIR/tsmp_idealrtd/pre-processing
   source loadenvs.Gnu_2020
   cd ../external
   chmod 755 compile
   ./compile
   cd ../pre-processing-tools
   ncl sva_surfdata_clm.ncl
   ncl sva_iniPress_pfl.ncl
   mv *.nc  ../input/clm
   mv *.pfb ../input/parflow
```

Note that you have the possibilty to generate different Input data for CLM and parFlow by in the scripts `sva_iniPress_pfl.ncl` and `sva_surfdata_clm.ncl`.
To configure TSMP for the IdealRTD test case on JUWELS machine (on JURECA just change -m JUWELS to -m JURECA):

```shell
   cd $TSMP_DIR/bldsva
   ./setup_tsmp.ksh -v 3.1.0MCT -V idealRTD -m JUWELS -O Intel
 ```


# Patching the orginal source code  <a name="patching "></a>

In order to prepare the original Cosmo5_1 source code for coupling, the Cosmo source code is patched using the diff files which are located in `$TSMP_DIR/bldsva/intf_oas3/cosmo5_1/pfile`.
The diff files are generated using the script `fpatch.sh` which is located in the same directory. Note that the patching is done only for cosmo5_1 source code.
If the user aims to modify the coupling source files, should modify the files located in `$TSMP_DIR/bldsva/intf_oas3/cosmo5_1/tsmp` and run again the script `fpatch.sh` in order to generate the the new diff files accordingly.
The necessary changes for making the ParFlow3.9 ready for coupling are already included in the official ParFlow3.9 release.
The changed source files for CLM3.5 are located in  `$TSMP_DIR/bldsva/intf_oas3/clm3_5/tsmp` and will be copied by TSMP scripts to the user's source code directory of CLM3.5 '(`bld/usr.src`)

# Long time climate simulation  <a name="climate_sim"></a>

## Implementation of non-const CO2 in TSMP <a name="CO2_Implementation"></a>

CO2 is different for different RCP scenarios, and it can impact largely the results of TSMP simulations. In order to change CO2 value in CLM according to the RCP scenario you need to change the value of the variable `co2_s` which is sent from COSMO to CLM and is hard-coded in `$TSMP_DIR/ bldsva/intf_oas3/cosmo5_1/oas3/send_fld_2clm.F90`.
Therefor the variable `co2_s` should be changed in order to vary CO2 in CLM depending on RCP scenarios. This means that if you want to do simulation for different RCP scenarios, you need to recompile the TSMP for each of them. Note that in order to change the CO2 value in COSMO according to RCP scenarios you need to change the COSMO namelist parameter `ico2_rad` in `INPUT_PHY`. For more

# Restart functionality for EURO-CORDEX test case  <a name="restart"></a>

TSMP generates automatically a script which does the simulation in two cycles.If interested please submit the job using the following:

```shell
   sbatch tsmp_restart.sh
```
This script can be a hint for you when you want to do long time climate simulation. For more info please contact "Niklas Wagner; n.wagner@fz-juelich.de"


# Automatic Porting of TSMP on x86 machines  <a name="TSMP_x86 "></a>

The users who want to port TSMP on GENERIC_X86 Linux, the TSMP team provided a script to install all the necessary libraries (Netcdf, GRIBAPI, OpenMPI, HDF5, TCL, Hypre and Silo) automatically in TSMP root directory. Please run the script "lib_install.sh" located in bldsva directory to install the libraries. Note that if you exported already one of the libraries Netcdf, HDF5, GRIBAPI, Silo, Hypre and TCL in the .baschrc or .profile, you need to comment them out in order to not mess up the installation via the script lib_install.sh.

The instruction on how to build and how to configure TSMP for the cordex test case ([EURO-CORDEX test case experiment](#ref_exp)) are given in [step 5 ](#ref_step4) and [step 6 ](#ref_step5), respectively.

# To come <a name="To-come"></a>

HPSC-TerrSys will also stage run control scripts for different usage scenarios such as long climate runs, and a variety of pre- and postprocessing tools we developed specifically for TSMP. Documented test case experiments also encompass convection permitting simulations at 1km for COSMO and 0.5km for CLM and ParFlow, idealized experiments, a data assimilation experiment based on the TSMP-PDAF version. We will also update the generic machinefiles to make an adjustment of the built system more straightforward for external users.

# Documentation <a name="ref_doc"></a>

An updated, extensive user guide, which describes the TSMP features, the technical implementation of the coupling, the different configurations, adjustments to different HPC systems as well as various test cases incl. all needed forcing datasets is currently under revision and will be available soon with the model system. See the TSMP website for a list of the individual websites and manuals of the component mdoels. To operate TSMP, the functioning and configuration of the component models has to be understood.
