# Getting Started

This small Getting Started section is mainly for users with access to JSC 
machines. For a general documentation on how to use TSMP, please refer to the 
more detailed [Building TSMP](./build_examples.md#building-tsmp) section.

## Step 1: Get the TSMP interface

Go to your preferred root directory (e.g., on a scratch file system) for the 
TSMP installation and get the `HEAD` of the `master` branch from `github` and 
set the environment variable `TSMP_DIR` to TSMP installation directory

```shell
git clone https://github.com/HPSCTerrSys/TSMP.git
cd TSMP
export TSMP_DIR=$(pwd)
```

## Step 2: Get the component models

Authenticate with your GitLab web GUI user name and password and clone
the repositories:

```shell
cd $TSMP_DIR
git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/cosmo5.01_fresh.git cosmo5_1
git clone -b v3.9.0 https://github.com/parflow/parflow.git parflow
git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/clm3.5_fresh.git clm3_5
git clone https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src/oasis3-mct.git oasis3-mct
```

## Step 3: Build TSMP, interface and component models

Building the fully coupled TSMP with ParFlow (pfl), the Community Land Model 
(clm) and the COSMO model (cos). This example is built on the JUWELS HPC system 
of Jülich Supercomputing Centre (JSC) using Intel compilers and ParaStation MPI:

```shell
cd $TSMP_DIR/bldsva
./build_tsmp.ksh --readclm=true -v 3.1.0MCT -c clm-cos-pfl -m JURECA -O Intel
```

## Step 4: Setup and prepare a test case

To configure TSMP for the EURO-CORDEX test case on JUWELS machine:
```shell
cd $TSMP_DIR/bldsva
./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m JUWELS -I _cordex -O Intel
```

Download needed input files for the EURO-CORDEX test case:
```shell
cd $TSMP_DIR/bldsva
./download_data_for_test_cases.ksh cordex
```

## Step 5: Run the test case

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
