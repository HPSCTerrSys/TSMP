#!/bin/csh
# Always make clean build for the first time
# For COSMO-CLM runs, we always update the CLM source codes for each soil moisture sensitivity run
# For COSMO runs, we only update the source code for one time for each soil moisture sensitivity run 
# So, when changing from RWU to CONSTANT, again a fresh compile is necessary **** WARNING
#
#1 time build is enough for COSMO runs (cos)
#The system has to be first build for one time for  COSMO-CLM runs (clm-cos) using "fresh"
#Then it will exit and from now on submit the run with "skip" but compile option set as 1
#due to change in source code for initializing soil moisture but COSMO and OASIS3 build can be
#skipped after the first compile, because changes will be only made in CLM source code later
# Only submit jobs for one vegetation type but different soil moisture at a time for COSMO-CLM runs
#
set cosSND = "raso_IdealSnd_0000LT_U01.dat"
set soiMOIS = "RWUMOD"   # OPTIONS, CONSTANT, RWU, RWUMOD
echo " "
echo "----------------------------------------"
echo "Usage: ./submit_smresponse.csh MODEL LANDCOVER SOILMOISTURE COMPILE  COPTION"
echo "MODEL : cos ,clm-cos"
echo "LANDCOVER: veg, bs"
echo "SOILMOUSTURE: 02,03,04,05,06"
echo "COMPILE: 0, 1"
echo "COPTION: fresh, skip"
echo "Examples ..."
echo "Usage: ./submit_smresponse.csh cos vg 03 0 fresh"
echo "Usage: ./submit_smresponse.csh clm-cos vg 03 1 skip"
echo "COSMO sounding file used: " ${cosSND}
echo "Soil Moisture Initialization: " ${soiMOIS}
echo "----------------------------------------"
echo " "

set expID  = "smresponse"
set mach   = "JURECA" 
set modID  = $1
set vegID  = $2
set swID   = $3
set cmpID  = $4
set bldmod = $5  #optional for clm-cos only
set smID   = ${vegID}${swID}
set nameID = ${modID}_${smID}

#Build Model

if ($modID == "cos" && $cmpID == 1) then
  ./build_tsmp.ksh -v 1.2.0MCT -c ${modID} -m $mach
  #Select which modifed source codes in /setup/smresponse folder to use for COSMO
  if (${soiMOIS} == "RWU") then
    cp $HOME/terrsysmp/bldsva/setups/$expID/cosmoSrc/src_artifdata.f90 $HOME/terrsysmp/cosmo5_1_${mach}_1.2.0MCT_${modID}/src/
    cd $HOME/terrsysmp/bldsva
    ./build_tsmp.ksh -v 1.2.0MCT -c ${modID} -m $mach -Y 'make' 
  endif
  exit 0
endif

#Select which modifed source codes in /setup/smresponse folder to use for CLM
if (${soiMOIS} == "RWU") then
  set clmSrc = "clmSrcRWU"
else if (${soiMOIS} == "RWUMOD") then
  set clmSrc = "clmSrcRWUMOD"
else if (${soiMOIS} == "CONSTANT") then
  set clmSrc = "clmSrc"
else
  echo "Not Defined ... " ${clmSrc}
  exit 1
endif

if ($modID == "clm-cos" && $cmpID == 1) then
  cd $HOME/terrsysmp/bldsva
  if ($bldmod == "fresh") then
   ./build_tsmp.ksh -v 1.2.0MCT -c ${modID} -m $mach -r 'true' 
   exit 0
  else
   cp $HOME/terrsysmp/bldsva/setups/$expID/${clmSrc}/mkarbinitMod_${swID}.F90 $HOME/terrsysmp/clm3_5_${mach}_1.2.0MCT_${modID}/bld/usr.src/mkarbinitMod.F90
   cp $HOME/terrsysmp/bldsva/setups/$expID/${clmSrc}/SurfaceAlbedoMod.F90 $HOME/terrsysmp/clm3_5_${mach}_1.2.0MCT_${modID}/bld/usr.src/
   #dummy copy for other except RWUMOD
   cp $HOME/terrsysmp/bldsva/setups/$expID/${clmSrc}/iniTimeConst.F90 $HOME/terrsysmp/clm3_5_${mach}_1.2.0MCT_${modID}/bld/usr.src/
   ./build_tsmp.ksh -v 1.2.0MCT -c ${modID} -m $mach -r 'true' -X 'make' -Y ${bldmod} -W ${bldmod}
  endif
endif

#Setup Experiment
#for cosmo
cp $HOME/terrsysmp/bldsva/setups/$expID/landcover/lmrun_uc5_1_${smID} $HOME/terrsysmp/bldsva/setups/$expID/lmrun_uc5_1
cp $HOME/database/smresponse/cosmo/${cosSND} $HOME/database/smresponse/cosmo/raso_IdealSnd_0000LT.dat 
#for clm
if ($vegID == "bs") then
  cp $HOME/database/smresponse/clm/bare-soil/*.nc $HOME/database/smresponse/clm/
else if ($vegID == "vg") then
  cp $HOME/database/smresponse/clm/vegetated/*.nc $HOME/database/smresponse/clm/
endif

cd $HOME/terrsysmp/bldsva

./setup_tsmp.ksh -Q "00:10:00" -v 1.2.0MCT -V $expID -c ${modID} -m $mach -r $WORK/run -I "_"$nameID

cd $WORK/run_${nameID}
sbatch tsmp_slm_run.bsh 
exit 0
endif
