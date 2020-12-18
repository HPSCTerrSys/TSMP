# TSMP

The Terrestrial System Modeling Platform (TSMP or TerrSysMP, https://www.terrsysmp.org) is an open source scale-consistent, highly modular, massively parallel regional Earth system model. TSMP essentially consists of an interface which couples dedicated versions of the Consortium for Small-scale Modeling (COSMO, http://www.cosmo-model.org) atmospheric model in NWP or climate mode, the Community Land Model (CLM, http://www.cesm.ucar.edu/models/clm/), and the hydrologic model ParFlow (https://www.parflow.org) through the OASIS3-MCT coupler (https://portal.enes.org/oasis, https://www.mcs.anl.gov/research/projects/mct/).

TSMP allows for a physically-based representation of transport processes of mass, energy and momentum and interactions between the different compartments of the geo-ecosystem across scales, explicitly reproducing feedbacks in the hydrological cycle from the groundwater into the atmosphere.

TSMP is extensively used for idealized and real data process and sensitivity studies in water cycle research, for climate change simulations, data assimilation studies including reanalyses, as well as experimental real time forecasting and monitoring simulations, ranging from individual catchments to continental model domains. TSMP runs on notebooks as well on latest supercomputers using a range of compilers.

TSMP development has been driven by groups within the [Centre for High-Performance Scientific Computing in Terrestrial Systems](http://www.hpsc-terrsys.de) (HPSC-TerrSys), as part of the [Geoverbund ABC/J](http://www.geoverbund-abcj.de/geoverbund/EN/Home/home_node.html), the geoscientific network of the University of Cologne, Bonn University, RWTH Aachen University, and the Research Centre Jülich. The current team is anchored in Jülich and Bonn in Germany.

**Visit**

**https://www.terrsysmp.org**

**for information on the features of TSMP, ongoing developments, citation, usage examples, links to documentation, the team, contact information and publications.**

## Citing TSMP

If you use TSMP in a publication, please cite the these papers that describe the model's basic functionalities:

* Shrestha, P., Sulis, M., Masbou, M., Kollet, S., and Simmer, C. (2014). A Scale-Consistent Terrestrial Systems Modeling Platform Based on COSMO, CLM, and ParFlow. Monthly Weather Review, 142(9), 3466–3483. doi:[10.1175/MWR-D-14-00029.1](https://dx.doi.org/10.1175/MWR-D-14-00029.1).
* Gasper, F., Goergen, K., Kollet, S., Shrestha, P., Sulis, M., Rihani, J., and Geimer, M. (2014). Implementation and scaling of the fully coupled Terrestrial Systems Modeling Platform (TerrSysMP) in a massively parallel supercomputing environment &ndash; a case study on JUQUEEN (IBM Blue Gene/Q). Geoscientific Model Development, 7(5), 2531-2543. doi:[10.5194/gmd-7-2531-2014](https://dx.doi.org/10.5194/gmd-7-2531-2014).

## Quick Start on Linux

This very short user guide only covers one variant on how the model can be setup and configured for *one* specific experiment, which we use as one of the default test cases. To get an overview on possible TSMP applications refer to the [TSMP website](https://www.terrsysmp.org) and the [user guide](#ref_doc).

### General concept

TSMP uses a multiple program multiple data (MPMD) execution model where states (e.g., near surface air temperature) and fluxes (e.g., latent and sensible heat fluxes) are passed between the atmospheric, land surface and hydrologic/groundwater component models along compartmental interfaces through the OASIS3 or the parallel OASIS3-MCT coupling interface.

TSMP essentially consists of the interface that combines the specific versions of different component models in a modular (various combinations of component models) and flexible (various spatial and temporal resolutions) way. During the install process, before the compilation of the component models, the source code of the component models is patched by the setup shell scripts, implementing the necessary coupling functionalities.

Please note: For licensing and maintenance reasons, we do not provide any source code of any of the component models or the coupler library. The users should download the component models from their respective websites. TSMP has various setup options: each component model standalone (ParFlow, CLM, COSMO), only land surface and hydrology / groundwater (CLM+ParFlow), only land surface and atmosphere (CLM+COSMO), or fully coupled (ParFlow, CLM, COSMO). Since the patches are designed for specific versions of the models, it is required to get the exact original, unaltered version as specified.

During the build and compile step ([step 5 below](#ref_step4)) multiple predefined coupled model combinations (of component models and their versions) as well as predefined HPC systems and operating environments can be chosen. Via so-called machinefiles the built and execution environment can be centrally adjusted to the respective HPC systems. A generic machinefile for x86-64 GNU/Linux systems with a GNU Compiler Collection and a standard scientific software stack is also provided. Many predefined built environments are provided for the machines which are available to the HPSC-TerrSys users.

During the setup and configuration step ([step 6 below](#ref_step5)), TSMP is configured for different, predefined test cases or numerical experiments, for which we also provide input data via a dedicated data server. Here we document one specific experiment. The [TSMP user guide](#ref_doc) will contain a complete list.

**Throughout this quick start guide a destinction is made between HPSC-TerrSys users who use specific HPC systems at HPC centres in Germany and external users. Through centralized machinefiles and compute environment initialisation scripts, TSMP can be very quickly adjusted to any other HPC site.**

The two most important scripts to build and setup TSMP are:

* `build_tsmp.ksh` to built and compile TSMP (one component model after the other)
* `setup_tsmp.ksh` to setup one of the predefined experiments

To get a man-page for the usage of this scripts, do:
```shell
   ./build_tsmp.ksh --man
   ./setup_tsmp.ksh --man
```
### <a name="ver_his"></a> TSMP version history
The model components used in TSMP are OASIS3-MCT v2, COSMO v5.01, CLM v3.5 and ParFlow 3.2 (for versions v1.2.1, v1.2.2 and v1.2.3) and ParFlow 3.7 for version v1.3.3. TSMP supports ParFlow 3.7 from version v1.3.3 onward.
Note that **ParFlow 3.7 is cloned automatically via TSMP from the GitHub repository** https://github.com/hokkanen/parflow.git.

### <a name="ref_exp"></a> The fully coupled pan-European EURO-CORDEX evaluation experiment with TSMP  

This test case uses the current TSMP release version v1.2.1 with OASIS3-MCT, COSMO v5.01, CLM v3.5 and ParFlow 3.2 (ParFlow 3.7 from TSMP version v1.3.3 onward).  A short 3hr simulation in a climate-mode configuration over Europe is set up, driven by ERA-Interim reanalysis, following the [EURO-CORDEX project](https://euro-cordex.net/) experiment guidelines. Simulated time span: 2016-05-01_12:00:00 to 2016-05-01_15:00:00.

### Step 1: Dependencies
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

A short guide on how the TSMP built system can be expanded to account for your local situation (software modules, compiler, MPI wrapper) will be provided shortly. The best starting point is to use the generic GNU compiler built. A hint on how to proceed is given in [step 5 below](#ref_step4).
The following libraries are required to run the EURO-CORDEX experiment:

### Step 2: Get the TSMP interface

Go to your preferred root directory (e.g., on a scratch file system) for the TSMP installation and get the `HEAD` of the `master` branch from `github` and set the environment variable `TSMP_DIR` to TSMP installation directory; for bash:

```shell
   git clone https://github.com/HPSCTerrSys/TSMP.git
   cd TSMP
   export TSMP_DIR=$(pwd)
```

### Step 3: Get the component models for this experiment

```shell
   cd $TSMP_DIR
```

#### HPSC-TerrSys users

Authenticate with your GitLab web GUI user name and password and clone the repositories (instead of "fresh", also "legacy" repositories with specific code modifications may be retrieved):

```shell
   git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/cosmo5.01_fresh.git
   git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/parflow3.2_fresh.git
   git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/clm3.5_fresh.git
   git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/oasis3-mct.git
```
It should be noted that ParFlow 3.7 is cloned automatically via TSMP v1.3.3 from the GitHub repository https://github.com/hokkanen/parflow.git. \
Rename the component model directories:

```shell
   mv cosmo5.01_fresh cosmo5_1
   mv clm3.5_fresh clm3_5
   mv parflow3.2_fresh parflow3_2
```

#### External users

Make sure directories are renamed as shown above. The component models can be downloaded from the respective websites as indicated below. Due to licensing reasons HPSC-TerrSys cannot provide these source codes to external users.

##### COSMO v5.01

Available from http://www.cosmo-model.org. A license agreement is needed.

##### CLM v3.5

Available from http://www.cgd.ucar.edu/tss/clm/distribution/clm3.5/.

##### ParFlow v3.2

Available from https://github.com/parflow/.

```shell
   git clone https://github.com/parflow/parflow.git
   mv parflow parflow3_2
   cd parflow
   git checkout v3.2.0
```
Alternatively:

```shell
   git clone --branch v3.2.0 https://github.com/parflow/parflow.git
   mv parflow parflow3_2
```
##### ParFlow v3.7
As indicated in [TSMP version history](#ver_his), this version is automatically cloned  and configured via TSMP interface. It is available from https://github.com/hokkanen/parflow.git under branch oas-gpu. 


##### OASIS3-MCT v2.0

Available from https://verc.enes.org/oasis/.

### Step 4: Retrieving the test case input data

For each officially documented and supported test case experiment, HPSC-TerrSys provides all required, pre-processed input data (as well as reference simulation results in the future).

The input files which are necessary for running the `EURO-CORDEX evaluation run` test case:

```shell
   cd $TSMP_DIR/bldsva
   ./download_data_for_test_cases.ksh
```

### <a name="ref_step4"></a> Step 5: Build TSMP, interface and component models

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
For building ParFlow 3.7 on GPU:

```shell
   cd $TSMP_DIR/bldsva
   ./build_tsmp.ksh -v 3.1.0MCT -c clm-cos-pfl -m JUWELS -O Intel -A GPU
```


A note to external users:

The path to the modules, compiler and MPI wrapper can be set in:  
`$TSMP_DIR/bldsva/machines/<your_machine_name>/build_interface_<your_machine_name>.ksh`

An example for the JUWELS HPC system at the Jülich Supercomputing Centre is here: `$TSMP_DIR/bldsva/machines/JUWELS/build_interface_JURECA.ksh` \
and an example to load the environment modules or the needed software compatible with Intel compiler version 2020 is in `$TSMP_DIR/bldsva/machines/JUWELS/loadenvs.Intel`.

### <a name="ref_step5"></a> Step 6: Setup and configuration of the respective usage and test case

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
In this heterogeneous job ParFlow 3.7 will run on GPU while Cosmo and CLM on CPU.

This includes the creation of a run directory, the copying of namelists, the provisioning of run control scripts for the job scheduler, incl. mapping files which pin the MPI tasks of the component model to specific CPU cores, as well as copying and linking of forcing data.

### Step 7: Run the test case

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

### Step 8: Simulation results

After a successful model run you will have model outputs form ParFlow, CLM and COSMO in the run directory:

```shell
   cd $TSMP_DIR/run/JURECA_3.1.0MCT_clm-cos-pfl_cordex_cordex
   ls -l *satur*.pfb *press*.pfb
   ls -l clmoas*.nc
   ls -l cosmo_out/*.nc
```

## To come

HPSC-TerrSys will also stage run control scripts for different usage scenarios such as long climate runs, and a variety of pre- and postprocessing tools we developed specifically for TSMP. Documented test case experiments also encompass convection permitting simulations at 1km for COSMO and 0.5km for CLM and ParFlow, idealized experiments, a data assimilation experiment based on the TSMP-PDAF version. We will also update the generic machinefiles to make an adjustment of the built system more straightforward for external users.

## <a name="ref_doc"></a> Documentation

An updated, extensive user guide, which describes the TSMP features, the technical implementation of the coupling, the different configurations, adjustments to different HPC systems as well as various test cases incl. all needed forcing datasets is currently under revision and will be available soon with the model system. See the TSMP website for a list of the individual websites and manuals of the component mdoels. To operate TSMP, the functioning and configuration of the component models has to be understood.
