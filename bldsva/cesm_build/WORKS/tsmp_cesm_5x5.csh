#! /bin/csh -f 
set cesm_rut = $HOME/cesm1_2_1
set cesm_dir = $HOME/test1

#userdefine input data parameters based on:
#models/lnd/clm/bld/namelist_files/namelist_defaults_usr_files.xml
#queryDefaultNamelist.pl -usrname "1x1_boulderCO" -options mask=navy,sim_year=1850,sim_year_range="constant" -csmdata /daten01/z4/database/cesm/inputdata
set GRIDNAME = 0300x0300NRW
set MASK = navy
set ATMDOM = domain.lnd.${GRIDNAME}_$MASK.nc
set yyyy_start = 2009
set yyyy_end = 2009

#cpl7 parameters
set cesm_res = 5x5_amazon
#set cesm_res = CLM_USRDAT
#set cesm_res = f45_g37 
set cesm_cst = I
#set cesm_cst = X
set cesm_mch = cluma2

#cluma2 machine configuration
cp config_compilers.xml $cesm_rut/scripts/ccsm_utils/Machines/
cp config_machines.xml  $cesm_rut/scripts/ccsm_utils/Machines/
cp env_mach_specific.cluma2 $cesm_rut/scripts/ccsm_utils/Machines/

cd
if (-d $cesm_dir) then 
  echo "Directory already exists!"
  exit 0
else
  rm -rf scratch
  echo "Running cesm setup"
  cd
  cd $cesm_rut/scripts
  ./create_newcase -v -case $cesm_dir -res $cesm_res -compset $cesm_cst -mach $cesm_mch 
  mkdir $cesm_dir/inputdata
  #CPS not needed ./link_dirtree $cesm_dat $cesm_dir/inputdata/
  cd $cesm_dir
  #
  echo "Updatexml for userspecific run"
  #env_run.xml
  #./xmlchange CLM_FORCE_COLDSTART="on"  
  #./xmlchange CLM_CO2_TYPE="constant"
  #./xmlchange CLM_NAMELIST_OPTS="hist_dov2xy=.true."
  #./xmlchange CLM_NAMELIST_OPTS="hist_avgflag_pertape=&apos;A&apos;"
  #./xmlchange CLM_USRDAT_NAME=$GRIDNAME
  #./xmlchange ATM_DOMAIN_PATH=$GENDOM_PATH,LND_DOMAIN_PATH=$GENDOM_PATH
  #./xmlchange ATM_DOMAIN_FILE=$ATMDOM,LND_DOMAIN_FILE=$ATMDOM
  #./xmlchange DATM_MODE=CLM1PT
  #./xmlchange DATM_CLMNCEP_YR_START=$yyyy_start,DATM_CLMNCEP_YR_END=$yyyy_end
  ./xmlchange DATM_CLMNCEP_YR_END=1972
  #
  #
  ./cesm_setup
  exit 0
endif

