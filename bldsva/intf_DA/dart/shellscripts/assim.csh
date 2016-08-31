#!/bin/csh
#
set DART_DIR="/home/pshrestha/DART/lanai/models/terrsysmp/cosmo/work"
set MODEL_PATH="/home/pshrestha/terrsysmp/run/CLUMA2_1.2.0MCT_clm-cos-pfl_idealRTD_perfectModel"
set num_days = 14 
#
cd $DART_DIR 
if ($1 == 0) then
  echo "Compiling..."
  ./quickbuild.csh -mpi
  rm input.nml.*
  echo "Creating Perfect Obs..."
else if ($1 == 1) then
  echo "Creating Perfect Obs..."
else
  echo "Code not writtten for ..." $1
  exit 1
endif
#
foreach irun (`seq 1 $num_days`)

set model_dir = `printf $MODEL_PATH%02d $irun`
set cosrst = `ls -1 $model_dir/cosrst/lrff* | tail -n -1`
set cosout = `ls -1 $model_dir/cosout/lfff* | tail -n -1`

echo `pwd`
echo " "
echo "-----------------------------------------------------------------------------"
rm cosmo.nc cosmo_prior
rm True_State.nc perfect_restart
rm obs_seq.in dart_log.*
ln -sf $cosout cosmo.nc
ln -sf $cosrst cosmo_prior
#&model_nml has the above names for restart and netcdf file

#./create_obs_sequence
#CPS Already generated set_def.out
#5.6 E, 49.8976 N, creates set_def.out


#creates obs_seq.in, do it for one restart time 
set defaultInitDate  = `grep defaultInitDate $model_dir/cosmo_prior_time.txt`
echo $defaultInitDate[2]
sleep 3 
echo " "
echo "-----------------------------------------------------------------------------"
./create_fixed_network_seq
wait

#creates dart_prior
./cosmo_to_dart

# creates 3 files including obs_seq.perfect 
./perfect_model_obs
mv obs_seq.perfect obs_seq.$defaultInitDate[2]
end
exit 0
