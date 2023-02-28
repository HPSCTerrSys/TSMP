# Setup Examples #

This Section need to be revised! Please take a look into the
[FallSchool
setup/configuration](https://gitlab.jsc.fz-juelich.de/sdlts/FallSchool_HPSC_TerrSys)
to test your model.

Active setups. Information on Backups by former colleagues
[here](#backups)

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
PDAF](./build_examples.md#compile-parflow-and-clm-with-pdaf). Choose the
correct command for your machine and, if possible, the newest version.

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

## CORDEX on master branch

TSMP with Stage2020 (with Intel and ParaStationMPI) for JUWELS is
ready to clone. Please do note forget to add the option `-O Intel` for
building and setup.

- to clone:

``` bash
git clone https://github.com/HPSCTerrSys/TSMP.git
```

-to build:

``` bash
./build_tsmp.ksh -v 3.1.0MCT -c clm-cos-pfl -m JUWELS -O Intel
```

- to create Cordex test case:

``` bash
./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m JUWELS -I _cordex -O Intel
```

## CORDEX with v1.3.3 with ParFlow3.7

TSMP version v1.3.3 is now upgraded to use ParFlow3.7 with a
possibility of submitting heterogeneous job (Cosmo and CLM on CPU and
ParFlow on GPU). TSMP will only support ParFlow 3.7 from version
v1.3.3 onward. Please read README in GitHub page for further
information.  Those who want to do a quick test, submitting a
heterogeneous job, please follow the instruction below:

1-

``` bash
git clone https://github.com/HPSCTerrSys/TSMP.git
```

2- download the model components (except ParFlow3.7) and place them in
TSMP root directory as written in "*Step 3*" of README".

3-

``` bash
cd bldsva
```

4-

``` bash
./build_tsmp.ksh -v 3.1.0MCT -c clm-cos-pfl -m JUWELS -O Intel *-A GPU*
```

5-

``` bash
./download_data_for_test_cases.ksh
```

6-

``` bash
./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m JUWELS -I _cordex -O Intel*-A GPU
```

*7- cd to run directory 8- modify the script `tsmp_slm_run.bsh`
according to your project account and submit a job using `sbatch
tsmp_slm_run.bsh`

If you want to submit TSMP job on CPU (Cosmo, CLM and ParFlow run on
CPU ), just drop `*-A GPU*` from build and setup commands.

## Setups ###

### Setup CORDEX test case for fully coupled TSMP on JUQUEEN with Oasis3-MCT ####

      ./setup_tsmp.ksh -m JUQUEEN -c clm-cos-pfl -v 1.1.0MCT -V cordex

The `-c` and `-v` flag are optional in this case because they default.

### Setup NRW test case for ParFlow standalone on CLUMA2 ####

      ./setup_tsmp.ksh -m CLUMA2 -c pfl -v 1.1.0MCT -V nrw

The `-m` , `-v` and `-V` flag are optional in this case because they are
default.

### Setup NRW test case for new CLM4.0 + Cosmo5.1 on JURECA ####

      ./setup_tsmp.ksh -m JURECA -c clm-cos -v 2.1.0MCT -C false

The new models are currently only supported for clm-cos and with
alternative coupling-scheme. 7.2.4) Setup NRW test case with 24h
runtime, 3h wallclock time, latest(lets say 01.07.2016 00:00) Cosmo
forcing, ParFlow and CLM restart, fully coupled on JURECA on `$WORK`

      ./setup_tsmp.ksh -m JURECA -c clm-cos-pfl -v 1.1.0MCT -V nrw -r "$WORK/tsmp/nrw20160701/run" -Q 3 -T 24 -s 2016-07-01_00 -F "$WORK/tsmp/nrw20160701/cosforcing" -j "$WORK/tsmp/nrw20160630/run/clmoas.clm2.r.2016-06-30-00000.nc" -l "$WORK/tsmp/nrw20160630/run/rurlaf.out.press.00024.pfb"

### Setup NRW test case with 20 member ensemble (parflow) on JURECA, but only use members 10-15. perturbed namelists reside in a special folder ####

      ./setup_tsmp.ksh -m JURECA -c clm-cos-pfl -v 1.1.0MCT -V nrw -N 5 -n 10 -g "$HOME/ensemble/namelists/coup\_oas.tcl"

Note, that even though you give the base namelist-name for the `-g`
flag the actual namelists in this folder must be named like:
`coup_oas.tcl_0`, `coup_oas.tcl_1`, `coup_oas.tcl_2`, etc.

### Setup NRW spinup with init date (lets say 01.01.2005 00:00) 1 month runtime, 24h wallclock time, restart from arbitrary date (lets say 01.03.2005 00:00) in indexed (3rd) rundir on JURECA ####

      ./setup_tsmp.ksh -m JURECA -c clm-cos-pfl -v 1.1.0MCT -V nrw -r "$WORK/tsmp/nrwSpinup/run" -I 3 -Q 24 -T 744 -s 2005-03-01_00 -S 2005-01-01_00 -j "$WORK/tsmp/nrwSpinup/run2/clmoas.clm2.r.2005-02-28-00000.nc" -k "$WORK/tsmp/nrwSpinup/run2/cosmo_out/lfff59000000" -l "$WORK/tsmp/nrwSpinup/run2/rurlaf.out.press.01416.pfb"



## Backups ##

- Zhenlei Yang: `icg4lts` under `DATAASSIMILATION/Zhenlei\ Yang`
