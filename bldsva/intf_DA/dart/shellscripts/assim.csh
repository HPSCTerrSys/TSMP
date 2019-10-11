#!/bin/csh
# Extract Synthetic Observations in five steps (Step 1 to 5)
# "Usage: ./assim.csh model" 
# P. Shrestha
# User Settings
# model ="clm", "cosmo", "parflow"
# inst -> randomly selected member number as PM
# st_day -> start day of experiment
# en_day -> end day of experiment - 1
set MODEL_PATH=$WORK/OL/rundart
set inst = 19      
set st_day = 30 
set en_day = 89 
# Do not change below
#

# SELECT COMPONENT MODEL FOR ASSIMILATION
if ($1 == "cosmo") then
  set DART_DIR="$TRAINHOME/DART/lanai/models/terrsysmp/cosmo/work"
else if ($1 == "clm") then
  set DART_DIR="$TRAINHOME/DART/lanai/models/terrsysmp/clm/work"
else if ($1 == "parflow") then
  set DART_DIR="$TRAINHOME/DART/lanai/models/terrsysmp/parflow/work"
endif

## SELECT COMPILATION
#if ($1 == 0) then
#echo "-------------------------------------------------------------------------------"
#echo "Update DART Source Code..."
#echo "-------------------------------------------------------------------------------"
#./git_to_dart.csh
#endif
#cd $DART_DIR 

#if ($1 == 0) then
#  echo "-------------------------------------------------------------------------------"
#  echo "Compiling..."
#  echo "-------------------------------------------------------------------------------"
#  ./quickbuild.csh -mpi
#  rm input.nml.*
#  echo "-------------------------------------------------------------------------------"
#  echo "Creating Perfect Obs..."
#  echo "-------------------------------------------------------------------------------"
#else if ($1 == 1) then
#  echo "-------------------------------------------------------------------------------"
#  echo "Creating Perfect Obs..."
#  echo "-------------------------------------------------------------------------------"
#else
#  echo "-------------------------------------------------------------------------------"
#  echo "Code not writtten for ..." $1
#  echo "Usage: ./assim.csh N c", N=0 or 1, c="clm", "cosmo", "parflow"
#  echo "-------------------------------------------------------------------------------"
#  exit 1
#endif
#
cd $DART_DIR 
#
foreach irun (`seq $st_day $en_day`)

# counter for column_rand outputs which start from 1 always
set ctr = `echo "($irun - $st_day + 1)" | bc`

if ($1 == "parflow") then
set model_dir = `printf $MODEL_PATH%02d $irun`/tsmp_instance_$inst
if ($irun == $st_day) then
#  set model_dir0 = $model_dir
   set model_dir0 = `printf $MODEL_PATH%02d 1`/tsmp_instance_$inst
endif
set prspfb = `ls -1 $model_dir/rurlaf.out.press*.pfb | tail -n -1`
set satpfb = `ls -1 $model_dir/rurlaf.out.satur*.pfb | tail -n -1`
set clmrst = `ls -1 $model_dir/clmoas.clm2.r.*.nc | tail -n -2 | head -n 1`
set clmrst0 = `ls -1 $model_dir0/clmoas.clm2.r.*.nc | tail -n -2 | head -n 1`
set pfbsoil = `ls -1 $TSMP_DATA/idealRTD/parflow/pfb_sID*.nc`
echo $clmrst
set pflgrd = `ls $model_dir/grids.nc`
 
set colum_exp = `printf column_export_out_%02d $ctr`
echo "------------------------------------------------------------------------    -------"
echo $colum_exp "    " $irun
echo " "
echo "------------------------------------------------------------------------    -------"
rm dart_prior pflgrid.nc clm_restart.nc 
rm perfect_restart
rm obs_seq.in dart_log.*
ln -sf $prspfb pfl_press.pfb 
ln -sf $satpfb pfl_satur.pfb 
ln -sf $pflgrd pflgrid.nc
ln -sf $clmrst clm_restart.nc
ln -sf $clmrst0 clm_restart_s.nc
ln -sf $pfbsoil pfb_soil.nc
#&model_nml has the above names for restart and netcdf file
endif  #$1 == "parflow"

if ($1 == "clm") then
set model_dir = `printf $MODEL_PATH%02d $irun`/tsmp_instance_$inst
if ($irun == $st_day) then
#  set model_dir0 = $model_dir
  set model_dir0 = `printf $MODEL_PATH%02d 1`/tsmp_instance_$inst
endif
set clmrst = `ls -1 $model_dir/clmoas.clm2.r.*.nc | tail -n -2 | head -n 1`
set clmrst0 = `ls -1 $model_dir0/clmoas.clm2.r.*.nc | tail -n -2 | head -n 1`
set clmout = `ls -1 $model_dir/clmoas.clm2.h0*.nc | tail -n -1`

set colum_exp = `printf column_export_out_%02d $ctr`
echo "-------------------------------------------------------------------------------"
echo `pwd`
echo " "
echo "-------------------------------------------------------------------------------"
rm clm_restart.nc clm_restart_s.nc clm_history.nc 
rm perfect_restart
rm obs_seq.in dart_log.*
ln -sf $clmout clm_history.nc
ln -sf $clmrst clm_restart.nc 
ln -sf $clmrst0 clm_restart_s.nc
#&model_nml has the above names for restart and netcdf file
endif  #$1 == "clm"

if ($1 == "cosmo") then
set model_dir = `printf $MODEL_PATH%02d $irun`/tsmp_instance_$inst
set cosrst = `ls -1 $model_dir/cosrst/lrff* | tail -n -1`
set cosout = `ls -1 $model_dir/cosout/lfff* | tail -n -1`
 
set colum_exp = `printf column_export_out_%02d $ctr`
 
echo "-------------------------------------------------------------------------------"
echo `pwd`
echo " "
echo "-------------------------------------------------------------------------------"
rm cosmo.nc cosmo_prior
rm perfect_restart
rm obs_seq.in dart_log.*
ln -sf $cosout cosmo.nc
ln -sf $cosrst cosmo_prior
endif #$1 == "cosmo"

echo $DART_DIR
if ($irun == $st_day) then
#Generates output for create_obs_sequence, hardwired for Radiosonde Tempearture
echo "-------------------------------------------------------------------------------"
echo "Step 1: Run column_rand ..."
echo "-------------------------------------------------------------------------------"
./column_rand
#Generate set_def.out
echo "-------------------------------------------------------------------------------"
echo "Step 2: Run create_obs_sequence"
echo "-------------------------------------------------------------------------------"
./create_obs_sequence < column_rand.out
endif

#creates obs_seq.in, do it for one restart time 
set defaultInitDate  = `grep defaultInitDate $model_dir/cosmo_prior_time.txt`
echo $defaultInitDate[2]
sleep 3 
echo " "
echo "-------------------------------------------------------------------------------"
echo "Step 3: Run create_fixed_network_seq..."
echo "-------------------------------------------------------------------------------"
./create_fixed_network_seq <$colum_exp
wait

echo "-------------------------------------------------------------------------------"
echo " Step 4: Run "${1}"_to_dart"
echo "-------------------------------------------------------------------------------"
./${1}_to_dart

# creates 3 files including obs_seq.perfect 
echo "-------------------------------------------------------------------------------"
echo "Step 5: Harvest model data: perfect_model_obs"
echo "-------------------------------------------------------------------------------"
./perfect_model_obs
mv obs_seq.perfect obs_seq.$defaultInitDate[2]
mv True_State.nc True_State.$defaultInitDate[2].nc
echo "-------------------------------------------------------------------------------"
echo "Exiting ..."
echo "-------------------------------------------------------------------------------"

end
exit 0
