# Setup Examples #

This Section need to be revised! Please take a look into the
[FallSchool
setup/configuration](https://gitlab.jsc.fz-juelich.de/sdlts/FallSchool_HPSC_TerrSys)
to test your model.

Active setups.

## Pan-European EURO-CORDEX

Please find the instructions for this case in the [getting started section](./../gettingstarted.md/#the-fully-coupled-pan-european-euro-cordex-evaluation-experiment-with-tsmp).

## Heterogeneous Job using TSMP

TSMP has the possibility of submitting heterogeneous job for EURO-CORDEX test experiment, meaning that ParFlow3.12  will run on GPU while Cosmo5.1 and CLM3.5 on CPU.
After locating the model components in the TSMP root as mentioned in [Step 3](./../gettingstarted.md/#step-3-get-the-component-models-for-this-experiment), the following command should be executed in order to build TSMP and to create the run directory for the [EURO-CORDEX test case experiment](./../gettingstarted.md/#the-fully-coupled-pan-european-euro-cordex-evaluation-experiment-with-tsmp) on JUWELS machine (on JURECA just change -m JUWELS to -m JURECA):

Building TSMP for HPSC-TerrSys users:
```shell
cd $TSMP_DIR/bldsva
./setup_tsmp.ksh -c clm3-cos5-pfl -V cordex -m JUWELS -I _cordex -O Intel -A GPU
```

Creating run directory for HPSC-TerrSys users:

```shell
cd $TSMP_DIR/bldsva
./setup_tsmp.ksh -c clm3-cos5-pfl -V cordex -m JUWELS -I _cordex -O Intel -A GPU
```

Not that heterogeneous job will be possible only for TSMP with ParFlow3.12 which is available from version v1.3.3.

## NRW Test case

### Retrieving the test case input data (NRW)

For `NRW` test case:

```shell
cd $TSMP_DIR/bldsva
./download_data_for_test_cases.ksh nrw
```

### Setup and configuration of the NRW test case

NRW Test case covers a geographical domain of 150 km x 150 km encompassing the North Rhine-Westphalia region, located in western Germany, Belgium, the Netherlands, and Luxembourg. The experiment is carried out in a clear sky day condition (08 May 2008) using the fully coupled (COSMO5.01-CLM3.5-
ParFlow) configuration of TSMP. The atmospheric component uses a constant lateral spatial resolution of about 1 km and a variable vertical discretization into 50 levels gradually coarsening from the bottom (20 m) to the top (22000 m). Initial and lateral boundary conditions for the atmospheric model are obtained from the operational weather forecast model
COSMO-DE of the German Weather Service (DWD).For downloading the necessary INPUT data for the NRW test case please see [Step 3](./../gettingstarted.md/#step-3-get-the-component-models-for-this-experiment) and for building the TSMP (with -c clm3-cos5-pfl) refer to [Step 5](./../gettingstarted.md/#step-5-build-tsmp-interface-and-component-models). For more information about NRW test case please refer to https://doi.org/10.3390/w10111697. \
To configure TSMP for the NRW test case on JUWELS machine (on JURECA just change -m JUWELS to -m JURECA):

```shell
cd $TSMP_DIR/bldsva
./setup_tsmp.ksh -c clm3-cos5-pfl -V nrw -m JUWELS -I _ClearSkyDay  -s 2008-05-08_00 -S 2008-05-08_00 -T 24 -O Intel
```

## IdealRTD (Idealized) Test case
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

### Preprocessing of Input Data
Here the user can generate the case-specific input data for ParFlow
and CLM. Note that the COSMO model uses the same setup for all the ensemble
members, that is, no need for ad-hoc editing. 
Please download the Input data and preproccesing script using:

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
./setup_tsmp.ksh -c clm3-cos5-pfl -V idealRTD -m JUWELS -O Intel
 ```

## Long time climate simulation

### Implementation of non-const CO2 in TSMP

CO2 is different for different RCP scenarios, and it can impact largely the results of TSMP simulations. In order to change CO2 value in CLM according to the RCP scenario you need to change the value of the variable `co2_s` which is sent from COSMO to CLM and is hard-coded in `$TSMP_DIR/ bldsva/intf_oas3/cosmo5_1/oas3/send_fld_2clm.F90`.
Therefor the variable `co2_s` should be changed in order to vary CO2 in CLM depending on RCP scenarios. This means that if you want to do simulation for different RCP scenarios, you need to recompile the TSMP for each of them. Note that in order to change the CO2 value in COSMO according to RCP scenarios you need to change the COSMO namelist parameter `ico2_rad` in `INPUT_PHY`. For more

## CLM5 in TSMP setup: nrw_5x

By Lukas Strebel

Cloning `TSMP` and `CLM5`
``` bash
	git clone https://icg4geo.icg.kfa-juelich.de/ExternalRepos/tsmp-pdaf/tsmp.git TSMP
	cd TSMP
	git checkout clm5-coupling

	git clone -b release-clm5.0.29 https://github.com/ESCOMP/ctsm.git clm5_0
	cd clm5_0
	./manage_externals/checkout_externals
	cd ..
```

Build and setup commands
``` bash
	cd bldsva
	./build_tsmp.ksh -v 4.4.0MCT -c clm -m JUWELS -O Intel
	# wait until build is finished - can take some time
	./setup_tsmp.ksh -v 4.4.0MCT -c clm -V nrw_5x -m JUWELS -O Intel
	cd #Rundir given as last output line

	# Change number of tasks / cores:
	# change all *_ntasks = x in drv_in
	# change nodes and task_per_node (and account for 1 node maybe increase time
	# limit to 2h) in tsmp_clm_run.bsh and then finally:

	sbatch tsmp_clm_run.bsh
```
If you want to submit multiple experiments at the same time, copy the full
Rundir and do changes in individual folders. Otherwise one at a time and save
timing folder between submits (is overwritten each time).

## For ICON Users (HPSC-TerrSys users)

### Step 0: JUDAC Storage Managment

Please inform yourself about the storage managment. If you are new to the Jülich supercomputers, it is recomannded to visit "Introduction to the usage and programming of supercomputer resources in Jülich" course. The next date of this ocurse can be found here: https://www.fz-juelich.de/ias/jsc/EN/Expertise/Workshops/Courses/courses_node.html . Please, make yourself familiar with the different partitions of the JUDAC system as well:  
> https://www.fz-juelich.de/ias/jsc/EN/Expertise/Datamanagement/JUDAC/FAQ/judac-FAQ_node.html

It is recomannded to store no data at $HOME, your model data at $PROJECT and your experiment data at $SCRATCH.

### Step 1: Dependencies
See [step 1 above](./../gettingstarted.md/#step-1-dependencies)

### Step 2: Get the TSMP interface

Go to your preferred root directory (e.g., your $PROJECT directory) for the TSMP installation and get the `TSMP_iconcoup` branch from `gitlab` and set the environment variable `TSMP_DIR` to TSMP installation directory; for bash:

```shell
git clone https://icg4geo.icg.kfa-juelich.de/spoll/tsmp.git
export TSMP_DIR=$(realpath tsmp)
git checkout TSMP_iconcoup
```

### Step 3: Get the component models for this experiment

Authenticate with your GitLab web GUI user name and password and clone the repositories (instead of "fresh", also "legacy" repositories with specific code modifications may be retrieved). Choose the model components needed from your experiment, in case of any coupled model one need the external coupler "oasis3-mct".

```shell
git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/icon2.1_legacy.git   icon2-1
git clone -b v3.12.0 https://github.com/parflow/parflow.git 							parflow
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
./build_tsmp.ksh -c clm-icon21-pfl -m JUWELS -O Intel
```

For ICON standalone version use:
```shell
cd $TSMP_DIR/bldsva
./build_tsmp.ksh -c icon21 -m JUWELS -O Intel
```

### Step 5: Setup and configuration of the respective usage and test case

It is recommended to store your model data at scratch (please notice that data will be automatically deleted if there are not touched for a while). Please, edit *YOUR_PROJECT* and *YOUR_DIRECTORY* correspondingly.

```shell
export EXP=/p/scratch/YOUR_PROJECT/YOUR_DIRECTORY
```
To configure TSMP for the Germany test case for ICON stand-alone on JUWELS machine:

```shell
cd $TSMP_DIR/bldsva
./setup_tsmp.ksh -V germany -m JUWELS -c icon21 -r $EXP/JUWELS_icon_4.1.0MCT_germany -I _experiment -O Intel
```

This includes the creation of a run directory, the copying of namelists, the provisioning of run control scripts for the job scheduler, incl. mapping files which pin the MPI tasks of the component model to specific CPU cores, as well as copying and linking of forcing data.

### Step 6: Run the test case

Change into the run directory:

```shell
cd $EXP/JUWELS_icon21_germany_experiment
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



## Setup TSMP-PDAF ##

For starting our data assimilation experiments with TSMP-PDAF, we will
consider a down- scaled version of the test case described in Kurtz et
al. (2016) which deals with a relatively simple land
surface-subsurface problem. With this forward model we will perform
experiments for the assimilation of soil moisture data including
parameter estimation. This testcase is also available in the TSMP
virtual machine or can be downloaded from gitlab in IBG-3
institute. For external users you can contact Prof. Dr. Harrie-Jan
Hendricks-Franssen at IBG-3 institute in Juelich Research Center.

## Day4: Testcase FallSchool 2022

Input files for the Testcase are found at
[https://gitlab.jsc.fz-juelich.de/sdlts/FallSchool_HPSC_TerrSys/tsmp_pdaf_exercise_data](https://gitlab.jsc.fz-juelich.de/sdlts/FallSchool_HPSC_TerrSys/tsmp_pdaf_exercise_data)

A detailed step-by-step explanation of how to build and run this setup
is given at
[https://gitlab.jsc.fz-juelich.de/sdlts/FallSchool_HPSC_TerrSys/fall-school-tutorials/-/wikis/Ensemble-Data-Assimilation](https://gitlab.jsc.fz-juelich.de/sdlts/FallSchool_HPSC_TerrSys/fall-school-tutorials/-/wikis/Ensemble-Data-Assimilation)

### Build

The common build for the Testcase FallSchool 2019 is [Compile
Parflow + CLM with
PDAF](./../build_tsmp/build_examples_tsmppdaf.md#compile-parflow-and-clm-with-pdaf). Choose the
correct command for your machine and, if possible, the newest version.


## Setups ###

### Setup CORDEX test case for fully coupled TSMP on JUWELS with Oasis3-MCT ####

      ./setup_tsmp.ksh -m JUWELS -c clm3-cos4-pfl -V cordex


### Setup NRW test case for ParFlow standalone on JURECA ####

      ./setup_tsmp.ksh -m JURECA -c pfl -V nrw

### Setup NRW test case for new CLM4.0 + Cosmo5.1 on JURECA ####

      ./setup_tsmp.ksh -m JURECA -c clm4-cos5 -C false

The new models are currently only supported for clm-cos and with
alternative coupling-scheme. 

### Setup NRW test case with 24h runtime, 3h wallclock time, latest(lets say 01.07.2016 00:00) Cosmo forcing, ParFlow and CLM restart, fully coupled on JURECA on `$WORK`

      ./setup_tsmp.ksh -m JURECA -c clm3-cos4-pfl -V nrw -r "$WORK/tsmp/nrw20160701/run" -Q 3 -T 24 -s 2016-07-01_00 -F "$WORK/tsmp/nrw20160701/cosforcing" -j "$WORK/tsmp/nrw20160630/run/clmoas.clm2.r.2016-06-30-00000.nc" -l "$WORK/tsmp/nrw20160630/run/rurlaf.out.press.00024.pfb"

### Setup NRW test case with 20 member ensemble (parflow) on JURECA, but only use members 10-15. perturbed namelists reside in a special folder ####

      ./setup_tsmp.ksh -m JURECA -c clm3-cos4-pfl  -V nrw -N 5 -n 10 -g "$HOME/ensemble/namelists/coup\_oas.tcl"

Note, that even though you give the base namelist-name for the `-g`
flag the actual namelists in this folder must be named like:
`coup_oas.tcl_0`, `coup_oas.tcl_1`, `coup_oas.tcl_2`, etc.

### Setup NRW spinup with init date (lets say 01.01.2005 00:00) 1 month runtime, 24h wallclock time, restart from arbitrary date (lets say 01.03.2005 00:00) in indexed (3rd) rundir on JURECA ####

      ./setup_tsmp.ksh -m JURECA -c clm3-cos4-pfl -V nrw -r "$WORK/tsmp/nrwSpinup/run" -I 3 -Q 24 -T 744 -s 2005-03-01_00 -S 2005-01-01_00 -j "$WORK/tsmp/nrwSpinup/run2/clmoas.clm2.r.2005-02-28-00000.nc" -k "$WORK/tsmp/nrwSpinup/run2/cosmo_out/lfff59000000" -l "$WORK/tsmp/nrwSpinup/run2/rurlaf.out.press.01416.pfb"
