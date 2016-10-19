#!/bin/csh
#
set DART_DIR="/home/pshrestha/DART/lanai/models/terrsysmp/cosmo/work"
set MODEL_PATH="/daten01/z4/database/dart/ideal_RTD_perfect/CLUMA2_1.2.0MCT_clm-cos-pfl_idealRTD_perfectModel"
set num_days = 13
#
if ($1 == 0) then
echo "-------------------------------------------------------------------------------"
echo "Update DART Source Code..."
echo "-------------------------------------------------------------------------------"
  ./git_to_dart_cosmo.csh
endif

cd $DART_DIR 
if ($1 == 0) then
  echo "-------------------------------------------------------------------------------"
  echo "Compiling..."
  echo "-------------------------------------------------------------------------------"
  ./quickbuild.csh -mpi
  rm input.nml.*
  echo "-------------------------------------------------------------------------------"
  echo "Creating Perfect Obs..."
  echo "-------------------------------------------------------------------------------"
else if ($1 == 1) then
  echo "-------------------------------------------------------------------------------"
  echo "Creating Perfect Obs..."
  echo "-------------------------------------------------------------------------------"
else
  echo "-------------------------------------------------------------------------------"
  echo "Code not writtten for ..." $1
  echo "Usage: ./assim.csh N", N=0 or 1
  echo "-------------------------------------------------------------------------------"
  exit 1
endif
#
foreach irun (`seq 1 $num_days`)

set model_dir = `printf $MODEL_PATH%02d $irun`
set cosrst = `ls -1 $model_dir/cosrst/lrff* | tail -n -1`
set cosout = `ls -1 $model_dir/cosout/lfff* | tail -n -1`

set colum_exp = `printf column_export_out_%02d $irun`

echo "-------------------------------------------------------------------------------"
echo `pwd`
echo " "
echo "-------------------------------------------------------------------------------"
rm cosmo.nc cosmo_prior
rm True_State.nc perfect_restart
rm obs_seq.in dart_log.*
ln -sf $cosout cosmo.nc
ln -sf $cosrst cosmo_prior
#&model_nml has the above names for restart and netcdf file

if ($irun == 1) then
#Generates output for create_obs_sequence, hardwired for Radiosonde Tempearture
echo "-------------------------------------------------------------------------------"
echo "Run column_rand ..."
echo "-------------------------------------------------------------------------------"
./column_rand
#Generate set_def.out
echo "-------------------------------------------------------------------------------"
echo "Run create_obs_sequence"
echo "-------------------------------------------------------------------------------"
./create_obs_sequence < column_rand.out
endif

#creates obs_seq.in, do it for one restart time 
set defaultInitDate  = `grep defaultInitDate $model_dir/cosmo_prior_time.txt`
echo $defaultInitDate[2]
sleep 3 
echo " "
echo "-------------------------------------------------------------------------------"
echo "Run create_fixed_network_seq..."
echo "-------------------------------------------------------------------------------"
./create_fixed_network_seq <$colum_exp
wait

#creates dart_prior
echo "-------------------------------------------------------------------------------"
echo " Run cosmo_to_dart"
echo "-------------------------------------------------------------------------------"
./cosmo_to_dart

# creates 3 files including obs_seq.perfect 
echo "-------------------------------------------------------------------------------"
echo "Harvest model data: perfect_model_obs"
echo "-------------------------------------------------------------------------------"
./perfect_model_obs
mv obs_seq.perfect obs_seq.$defaultInitDate[2]

echo "-------------------------------------------------------------------------------"
echo "Exiting ..."
echo "-------------------------------------------------------------------------------"

end
exit 0
