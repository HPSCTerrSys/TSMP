#! /bin/csh -f 
#Author P. Shrestha
#email: pshrestha@uni-bonn.de
#28-11-2015 Porting of CESM to cluma2 for offline CLM RUNS
#           Users need to edit the following files for their runs or use the default NRW case
#           namelist_defaults_datm.xml  -- for domain directory and domain file of atm. forcing
#           user_nl_clm -- This is similar to lnd.stdin of CLM3.5, add your own namelists
#           user_nl_datm -- For datm specific parameters
#           $cesm_bld/SourceMods/ -- Include update CLM source codes          

#------------------------------------------------------------------------------------------------
# USER SETTINGS
#------------------------------------------------------------------------------------------------
set cesm_rut = $SVAROOT/cesm1_2_1
#set cesm_rut = $HOME/post4.5crop_slevis
set cesm_cas = tempbld_clm
#set cesm_cas = test_slevis
set cesm_dir = $SVAROOT/$cesm_cas
set cesm_bld = $SVAROOT/bldsva/build_cesm
#userdefine input data parameters based on:
#models/lnd/clm/bld/namelist_files/namelist_defaults_usr_files.xml
#queryDefaultNamelist.pl -usrname "1x1_boulderCO" -options mask=navy,sim_year=1850,sim_year_range="constant" -csmdata /daten01/z4/database/cesm/inputdata
set GRIDNAME = 300x300pt_NRW 
set LMASK = navy
set ATMDOM = domain.lnd.${GRIDNAME}_$LMASK.nc
set yyyy_start = $yyyy 
set mm_start   = $mm 
set dd_start   = $dd 
set yyyy_end = $yyyy
set stop_option = nhours
set stop_optionVal = $runhours 
set nATM_CPL = 24
#Customizing PE Layout
@ NPROC_LND = ( $nprocx_clm * $nprocy_clm )
set NTHRD_LND = 1
set NROOT_LND = 0 
#cpl7 parameters
#set cesm_res = 5x5_amazon
set cesm_res = CLM_USRDAT
#set cesm_res = f45_g37 
set cesm_cst = I
set clm_ver  = clm4_0
#set cesm_cst = X
set cesm_mch = cluma2

#------------------------------------------------------------------------------------------------
# USER SETTINGS END, DO NOT CHANGE BELOW
#------------------------------------------------------------------------------------------------

#Update Machine Configuration
cp $cesm_bld/Machines/$cesm_mch/config_compilers.xml $cesm_rut/scripts/ccsm_utils/Machines/
cp $cesm_bld/Machines/$cesm_mch/config_machines.xml  $cesm_rut/scripts/ccsm_utils/Machines/
cp $cesm_bld/Machines/$cesm_mch/env_mach_specific.cluma2 $cesm_rut/scripts/ccsm_utils/Machines/
cp $cesm_bld/Machines/$cesm_mch/mkbatch.cluma2 $cesm_rut/scripts/ccsm_utils/Machines/
#Update domain directory and domain file for DATM OFFLINE MODE
cp namelist_defaults_datm.xml $cesm_rut/models/atm/datm/bld/namelist_files/
cd
if (-d $cesm_dir) then 
  echo "Directory already exists!"
  exit 0
else
  #CPS rm -rf scratch
  echo "Running cesm setup"
  cd
  cd $cesm_rut/scripts
  ./create_newcase -v -case $cesm_dir -res $cesm_res -compset $cesm_cst -mach $cesm_mch 
  #CPS not needed ./link_dirtree $cesm_dat $cesm_dir/inputdata/
  cd $cesm_dir
  #
  echo "Updatexml for userspecific run"
  #env_machine_pes.xml
  ./xmlchange NTASKS_LND=$NPROC_LND
  ./xmlchange NTHRDS_LND=$NTHRD_LND
  ./xmlchange ROOTPE_LND=$NROOT_LND
  #env_build.xml
  ./xmlchange MASK_GRID=$LMASK
  #env_run.xml
  ./xmlchange CLM_BLDNML_OPTS="-mask $LMASK" 
  ./xmlchange CLM_FORCE_COLDSTART="off"
  ./xmlchange CLM_CO2_TYPE="constant"
  #./xmlchange CLM_NAMELIST_OPTS="hist_dov2xy=.true.;"
  #./xmlchange CLM_NAMELIST_OPTS="hist_avgflag_pertape=&apos;A&apos;"
  ./xmlchange CLM_USRDAT_NAME=$GRIDNAME
  ./xmlchange ATM_DOMAIN_FILE=$ATMDOM,LND_DOMAIN_FILE=$ATMDOM
  ./xmlchange DATM_MODE=CLM1PT
  ./xmlchange DATM_CLMNCEP_YR_START=$yyyy_start,DATM_CLMNCEP_YR_END=$yyyy_end
  ./xmlchange DATM_CLMNCEP_YR_ALIGN=$yyyy_start
  ./xmlchange CLM_BLDNML_OPTS="-sim_year 2000"
  ./xmlchange RUN_STARTDATE=$yyyy_start"-"$mm_start"-"$dd_start
  ./xmlchange RUN_REFDATE=$yyyy_start"-"$mm_start"-"$dd_start
  ./xmlchange STOP_OPTION=$stop_option
  ./xmlchange STOP_N=$stop_optionVal
  ./xmlchange ATM_NCPL=$nATM_CPL
  #
  ./cesm_setup
  #Change Source Files
  cp $cesm_bld/SourceMods/models/lnd/clm/src/$clm_ver/* $cesm_dir/SourceMods/src.clm/ 
  cp $cesm_bld/user_nl_clm . 
  cp $cesm_bld/user_nl_datm .
  ./$cesm_cas.build
#  ./$cesm_cas.submit  
  exit 0
endif

