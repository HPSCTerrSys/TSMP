
#********* History********************* #

15.01.2020 written by Abouzar Ghasesmi

#********* History********************* #

What comes in the follwoing is a short readme for Terrsysmp package for external and internal users (within Forschungszentrum Jülich). 
As an example to test the Terrsysmp system you will find in the step 6 an instruction to build, setup and run Euro-cordex.11 test case.

1- Cloning the terrsysmp package which is available in TSMP-GitHub repository using:
git clone https://......GitHub/terrsysmp
The terrsysmp package contains only a coupling interface for different versions of model components.
IGB3 institue provides the model components for internal users (within Forschungszentrum Jülich) as explained in the step 1-2, while external users 
should provide them by thier own as explained in the step 1-1.

1-1 Externla users should download the model component from the respective websites as indicated in the steps 2,3,4 and 5.
External users are also requested to take care of the software modules as well as compiler and MPI wrapper that needed to be loaded for building Terrsysmp system. 
The path to the modules,compiler and MPI wrapper should be set in .../terrsysmp/bldsva/machines/<machin name>/build_interface_<machin name>.ksh. 
See for example  ../terrsysmp/bldsva/machines/JUWELS/build_interface_JUWELS.ksh.

1-2 Internal users can download the model component from IGB3-gitlab repository (https://icg4geo.icg.kfa-juelich.de/ModelSystems/tsmp_src)
Since the authentification to IGB3-gitlab is based on a ssh key, then internal users are requested to copy public thier SSH-key to IGB3-gitlab repository.
The model components available in IGB3-gitlab are: Cosmo5.01, Cosmo4_21, Parflow3.2, Parflow, CLM3.5, oasis3-mct, pdaf1_1.

1-3 Note also that model components should be copied inside the terrsysmp directory and be renamed as suggested in the follwoing.

2- Atmospheric Model (COSMO/ICON)
COSMO/ICON itself is free for academic R&D work.

- COSMO5_1
Get license (https://git2.meteo.uni-bonn.de/projects/terrsysmp/wiki/License) and download from DWD.
Rename the folder to cosmo5_1

- ICON
The ICON LEM source code can be obtained via GIT. To access the GIT repository you have to be registered. The authentification to GIT is based on a ssh key.
Create an ssh-key (see below) and send public key to Monika Esch (MPI-M) at monika.esch@mpimet.mpg.de and describe why you need ICON for.
If your access to ICON git repository is approved, download ICON as follows:
1) git clone --recursive git@git.mpimet.mpg.de:icon-lem.git
2) cd icon-lem
3) git checkout -b icon-lem-2.1.00-alpha remotes/origin/icon-lem-2.1.00-alpha
4) git submodule update
Rename the folder to icon-lem

3- Land Surface Model (CLM)
The repository for CLMX.X are different as CLM3.5 is available as an offline model but the newer versions are part of the CESM repository

- CLM3.5
Available from http://www.cgd.ucar.edu/tss/clm/distribution/clm3.5/
Rename the folder to clm3_5

- CLM4.0
Users will need to register in the CESM portal: http://www2.cesm.ucar.edu/models/cesm1.2/
svn co -r 76850 --username guestuser --password [password] https://svn-ccsm-models.cgd.ucar.edu/cesm1/release_tags/cesm1_2_1 cesm1_2_1
OR
svn co -r 76850 https://svn-ccsm-models.cgd.ucar.edu/cesm1/release_tags/cesm1_2_1 cesm1_2_1
Rename the folder to clm4_0 and copy it to ../terrsysmp directory
- CLM5.0
UPDATE PENDING ...

3- Groundwater Model (ParFlow)
Available from https://github.com/parflow/parflow
Information on releases available at https://github.com/parflow/parflow/releases
git clone --branch v3.2.0 https://github.com/parflow/parflow.git

4- External Coupler (OASIS3-MCT)
OASIS3-MCT-V2.0 is used to coupl the model components, and is available in Terrsysmp GitHub repository under the CERFACS license (https://verc.enes.org/oasis).
Rename the folder to oasis3-mct

5- Data Assimilation Framework (PDAF)
wikipage for Parallel Data Assimilation Framework (PDAF) http://pdaf.awi.de/trac/wiki
To download new releases for PDAF, you need to register at http://pdaf.awi.de/register/index.php
Currently the newest available PDAF Version 1.13.2, has been successfully tested with Terrsysmp Interface. Release notes for PDAF's different version can be found at http://pdaf.awi.de/trac/wiki/ReleaseNotes
When you download the PDAF after registering, you need to extract or untar the PDAF files and folders. you need to replace the older pdaf1_1 in terrsysmp/pdaf1_1 with newer one. So if you are downloading PDAF Version 1.13.2. After extracting the tar files you get a folder PDAF-D_V1.13.2. Just rename this folder from PDAF-D_V1.13.2 to pdaf1_1 and place it in terrsysmp folder.
Please note pdaf1_1 in path terrsysmp/bldsva/intf_DA/pdaf1_1 is the terrsysmp interface to terrsysmp/pdaf1_1, it doesnt have the data assimilation algorithms.

6- Building terrsysmp system and setup the Cordex.11 test case.
After locating the model components in terrsysmp directory, 
follwoing the steps below, you can build terrsysmp system and setup the cordex.11 test case.

6-1 Downloading the input files necessary for running Euro-cordex.11 test case:
cd to the root directory .../terrsysmp and run the script "download_input_cordex11.ksh"

6-2 cd to .../terrsysmp/bldsva directory and build the coupled "clm-cos-pfl" system using Intel compiler.

./build_tsmp.ksh -v 3.1.0MCT -c clm-cos-pfl -m <machin name> -O Intel

6-3 (Optional, better not to mess up the root directory and run directory) export your work directory using
export WORK=

6-4 If you export your work directory then run:
./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m <machin name> -I _cordex -r $WORK/run -O Intel
if NOT
./setup_tsmp.ksh -v 3.1.0MCT -V cordex -m <machin name> -I _cordex  -O Intel

6-5 cd to the run directory and change "#SBATCH " lines in "tsmp_slm_run.bsh" script accordingly.

6-6 Submit the job using "sbatch tsmp_slm_run.bsh".






