# TSMP
This readme contains a short instruction for external and internal users (within Forschungszentrum Jülich) on how to prepare TSMP package and how to provide the model components. 
As an example to test the TSMP system, you will find in the step 6, an instruction to build, setup and run Euro-cordex.11 test case.

## Table of contents
* [Introduction]
* [Downloading TSMP Interface]
* [Downloading Model Components]
* [Creating and running the Euro-cordex standard test case]
* [Contact](#contact)

## Introduction
The TAMP package contains only a coupling interface for different versions of model components.
IGB3 institute provides the model components for internal users, who are working in Forschungszentrum Jülich, as explained below. 
External users should provide the model components by thier own. 
Through this readme $basedir is the address of terrsysmp directory in your machine (cd to terrsysmp directory and then export basedir=`pwd`).  

* Internal users
The model components cane be downloaded from IGB3-Gitlab repository (https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src) .
Since the authentification to IGB3-Gitlab is based on a SSH-key, then internal users are requested to copy thier public SSH-key to IGB3-Gitlab repository.
The model components available in IGB3-Gitlab are: Cosmo5.01, Parflow3.2, CLM3.5, oasis3-mct.

* External users
The model components should be downloaded from the respective websites as indicated in the following section.
External users are also requested to take care of the software modules as well as compiler and MPI wrapper that needed to be loaded for building TSMP interface. 
The path to the modules,compiler and MPI wrapper should be set in $basedir/bldsva/machines/<machin name>/build_interface_<machin name>.ksh. 
See for example  $basedir/bldsva/machines/JUWELS/build_interface_JUWELS.ksh.

Note also that model components should be copied inside the $basedir directory and be renamed as suggested in the follwoing.

## Downloading TSMP Interface
In your preferred directory and for current release version, type

git clone https://......GitHub/terrsysmp
cd terrsysmp
export basedir=`pwd`

## Downloading Model Components
* Atmospheric Model (COSMO)
COSMO is free for academic R&D work.

- COSMO5_1
External users need to get a license (https://git2.meteo.uni-bonn.de/projects/terrsysmp/wiki/License) and download the model from DWD or http://www.cosmo-model.org .
Internal users:
cd to $basedir
git clone git@icg4geo.icg.kfa-juelich.de:ModelSystems/tsmp_src/cosmo5.01.git
Rename the folder to cosmo5_1:
mv cosmo5.01 COSMO5_1

* Land Surface Model (CLM)
The repositories for CLMX.X are different as CLM3.5 is available as an offline model but the newer versions are part of the CESM repository

- CLM3.5
External users need to get CLM3.5 from http://www.cgd.ucar.edu/tss/clm/distribution/clm3.5/
Internal users:
cd to $basedir and download the model:
git clone git@icg4geo.icg.kfa-juelich.de:ModelSystems/tsmp_src/clm3.5.git
Rename the folder to clm3_5:
mv clm3.5 clm3_5


* Groundwater Model (ParFlow)
External users need to get parflow3.2 from https://github.com/parflow/parflow.
Information on releases available at https://github.com/parflow/parflow/releases.
cd to $basedir and download the model and Download the model:
git clone --branch v3.2.0 https://github.com/parflow/parflow.git

Internal users:
cd to $basedir and download the model:
git clone git@icg4geo.icg.kfa-juelich.de:ModelSystems/tsmp_src/parflow3.2.git
git clone --branch v3.2.0 https://github.com/parflow/parflow.git

Rename the folder to parflow3_2:
mv parflow3.2 parflow3_2

* External Coupler (OASIS3-MCT)
OASIS3-MCT-V2.0 is used to coupl the model components, and is available in TSMP-GitHub repository under the CERFACS license (https://verc.enes.org/oasis).
Both Internal and External users can download the model from TSMP-GitHub.
----

Internal users can also downlaod it from IGB3-Gitlab.
cd to $basedir and download it via:
git clone git@icg4geo.icg.kfa-juelich.de:ModelSystems/tsmp_src/oasis3-mct.git

## Creating and running the Euro-cordex standard test case.
After locating the model components in $basedir directory, 
follwoing the steps below, you can build TSMP system and setup the cordex.11 test case.

1- Downloading the input files necessary for running Euro-cordex.11 test case:
cd to $basedir directory and run the script "download_input_cordex11.ksh"

2- cd to $basedir/bldsva directory and build the coupled "clm-cos-pfl" system using Intel compiler.

./build_tsmp.ksh -v 3.1.0MCT -c clm-cos-pfl -m <machin name> -O Intel

3- (Optional, better not to mess up the root directory and run directory) export your work directory using
export WORK=

4- If you export your work directory then run:
./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m <machin name> -I _cordex -r $WORK/run -O Intel
if NOT
./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m <machin name> -I _cordex  -O Intel

5- cd to the run directory and change "#SBATCH " lines in "tsmp_slm_run.bsh" script accordingly.

6- Submit the job using "sbatch tsmp_slm_run.bsh".

## Contact
Dr. Abouzar Ghasemi
Jülich Supercomputing Centre (JSC)
Wilhelm-Johnen-Straße
52425 Jülich, Germany
Phone: +49 2461 61-2312
Fax: +49 2461 61-6656
email: a.ghasemi@fz-juelich.de
