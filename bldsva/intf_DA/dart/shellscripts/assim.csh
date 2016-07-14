#echo "`date` -- Start initial run"
#qsub tsmp_pbs_run.ksh
#while (1)
#  if (-f "ready.txt") then
#    echo "`date` -- END initial run"
#    echo " "
#    break
#  else
#    echo "Job is running..."
#    sleep 10 
#  endif
#end

./quickbuild.csh -mpi
rm input.nml.*
#
ln -sf /home/pshrestha/terrsysmp/run/CLUMA2_1.2.0MCT_clm-cos-pfl_idealRTD_dart01/tsmp_instance_0/cosout/lfff01000000.nc cosmo.nc
ln -sf /home/pshrestha/terrsysmp/run/CLUMA2_1.2.0MCT_clm-cos-pfl_idealRTD_dart01/tsmp_instance_0/cosrst/lrff01000000o cosmo_prior
#&model_nml has the above names for restart and netcdf file
./create_obs_sequence
#5.61481 E, 49.9066 N, creates set_def.out
cat set_def.out
/create_fixed_network_seq
#creates obs_seq.in, not sure about input period of days in sequence
cat obs_seq.in
./cosmo_to_dart
#creates dart_prior
# change namelist 
./perfect_model_obs
# creates 3 files including obs_seq.perfect 

