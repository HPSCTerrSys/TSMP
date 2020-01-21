# TSMP Interface
  This readme contains a short instruction for external and internal users (within Forschungszentrum Jülich) on how to prepare TSMP package and how to provide the model components.  
  As an example to test the TSMP system, you will find in the step 6, an instruction to build, setup and run Euro-cordex.11 test case.

## Table of contents
  * [Introduction]()
  * [Downloading TSMP Interface]()
  * [Downloading Model Components]()
  * [Creating and running the Euro-cordex standard test case]()
  * [Contact]()

## Introduction

  The TAMP package contains only a coupling interface for different versions of model components.
  IGB3 institute provides the model components for the internal users, who are working in Forschungszentrum Jülich, as explained below.
  External users should provide the model components by their own.
  Through this readme $basedir is the address of terrsysmp directory in your machine.
  Please cd to terrsysmp directory and then run "export basedir=`pwd`".

  * Internal users:  
    The model components cane be downloaded from IGB3-Gitlab repository [`https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src`]() .  
    Since the authentication to IGB3-Gitlab is based on a SSH-key, then internal users are requested to copy their public SSH-key to IGB3-Gitlab repository.  
    The model components available in IGB3-Gitlab are:  
    Cosmo5.01, Parflow3.2 and CLM3.5   
    Note that `_legacy` suffix means that the source codes are partially modified by IGB3,  
    and `_fresh` suffix means the fresh and original source codes without any manipulation.

  * External users:  
    The model components should be downloaded from the respective websites as indicated in the following section.  
    External users are also requested to take care of the software modules as well as compiler and MPI wrapper that needed to be loaded for building TSMP interface.    
    The path to the modules, compiler and MPI wrapper should be set in:  
    `$basedir/bldsva/machines/<machin name>/build_interface_<machin name>.ksh`.<br/>
    See for example  `$basedir/bldsva/machines/JUWELS/build_interface_JUWELS.ksh`.

  * Renaming the folders containing the model components:  
  The model components should be copied inside the `$basedir` directory and be renamed as:
  `mv Cosmo5.01_<  > cosmo5_1`  
  `mv CLM3.5_<  > clm3_5`  
  `mv parflow3.2_<  > parflow3_2`  

  * For further Information the user can read the manuals provided in:  
    [`https://www.terrsysmp.org/`]()  
    [`https://git.meteo.uni-bonn.de/projects/terrsysmp/wiki/wiki`]()  

## Downloading TSMP Interface
  In your preferred directory and for current release version, type  
  `git clone https://......GitHub/terrsysmp`  
  `cd terrsysmp`  
  `export basedir=`pwd` `

## Downloading Model Components
  * Atmospheric Model (COSMO)
    COSMO is free for academic R&D work.

    - COSMO5_1  
      - External users:  
        First get a license from [`https://git2.meteo.uni-bonn.de/projects/terrsysmp/wiki/License`]() and then download the model from DWD or http://www.cosmo-model.org.<br/>

      - Internal users:  
        `cd $basedir`  
        Download the legacy source code:  
        `git clone git@icg4geo.icg.kfa-juelich.de:ModelSystems/tsmp_src/cosmo5.01_legacy.git`  
        Download the fresh source code:  
        `git clone git@icg4geo.icg.kfa-juelich.de:ModelSystems/tsmp_src/cosmo5.01_fresh.git`

      - Rename the folder to cosmo5_1:
        `mv Cosmo5.01_<  > cosmo5_1`

  * Land Surface Model (CLM)<br/>

    - CLM3.5  
      - External users:  
        Can get CLM3.5 from [`http://www.cgd.ucar.edu/tss/clm/distribution/clm3.5/`]()

      - Internal users:  
        `cd $basedir`  
        Downlaod the legacy source code:
        `git clone  git@icg4geo.icg.kfa-juelich.de:ModelSystems/tsmp_src/clm3.5_legacy.git`  
        Download the fresh source code:
        `git clone git@icg4geo.icg.kfa-juelich.de:ModelSystems/tsmp_src/clm3.5_fresh.git`

      - Rename the folder to clm3_5:  
      `mv CLM3.5_<  > clm3_5`

  * Groundwater Model (ParFlow)  
    - Parflow3.2  
      - External users:
        can be retrieved from [`https://github.com/parflow/parflow`]().<br/>
        Information on releases are available at [`https://github.com/parflow/parflow/releases`](). <br/>
        `cd $basedir`  
        Download the model:  
        `git clone --branch v3.2.0 https://github.com/parflow/parflow.git`

      - Internal users:  
        `cd $basedir`  
        Download the legacy source code:  
        `git clone git@icg4geo.icg.kfa-juelich.de:ModelSystems/tsmp_src/parflow3.2_legacy.git`  
        Download the fresh source code:  
        `git clone git@icg4geo.icg.kfa-juelich.de:ModelSystems/tsmp_src/parflow3.2_fresh.git`  

      - Rename the folder to parflow3_2:
        `mv parflow3.2_<  > parflow3_2`

  * External Coupler (OASIS3-MCT)  
    OASIS3-MCT-V2.0 is used to couple the model components, and is available in TSMP-GitHub repository under the CERFACS license [`https://verc.enes.org/oasis`](https://verc.enes.org/oasis).<br/>
    Both Internal and External users can download the model from TSMP-GitHub.
    - External users:  
      Download it via:    
      `github link`

    - Internal users:  
      Can also download it from IGB3-Gitlab.  
      `cd $basedir`  
      Download it via:  
      `git clone git@icg4geo.icg.kfa-juelich.de:ModelSystems/tsmp_src/oasis3-mct.git`

## Creating and running the Euro-cordex standard test case.
  After locating the model components in `$basedir` directory,
  following the steps below, you can build TSMP system, setup and run the Euro-cordex test case.

  - Downloading the input files which are necessary for running Euro-cordex.11 test case:  
    `cd $basedir/bldsva` and then run the script "download_input_cordex11.ksh"  

  - Building the coupled "clm-cos-pfl" system.  
    `cd $basedir/bldsva`  
    Compile the system using Intel compiler.  
    `./build_tsmp.ksh -v 3.1.0MCT -c clm-cos-pfl -m <machin name> -O Intel`

  - This step is optional. Better not to mess up the base directory and run directory.   
    Export your work directory using:  
    `export WORK=  `

  - If you export your work directory then run:  
  `./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m <machin name> -I _cordex -r $WORK/run -O Intel`  
  If NOT:  
  `./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m <machin name> -I _cordex  -O Intel`

  - Go to the run directory and change `#SBATCH ` lines in `tsmp_slm_run.bsh` script accordingly.  

  - Submit the job using `sbatch tsmp_slm_run.bsh`.

## Contact
Abouzar Ghasemi  
Jülich Supercomputing Centre (JSC)  
E-mail: a.ghasemi@fz-juelich.de  
Klaus Görgen  
E-Mail: k.goergen@fz-juelich.de
