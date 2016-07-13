#!/bin/ksh

#Job Submission to Cluma
#PBS -N TerrSysMP_run
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=64
#PBS -V 
#PBS -u pshrestha
#PBS -q batch

export LOGNAME="/home/pshrestha/terrsysmp/run/CLUMA2_1.2.0MCT_clm-cos-pfl_idealRTD_dart01"
export DART_DIR="/home/pshrestha/DART/lanai/models/terrsysmp/cosmo/work"
cd /home/pshrestha/terrsysmp/run/CLUMA2_1.2.0MCT_clm-cos-pfl_idealRTD_dart01

cp $DART_DIR/input.nml .
cp $DART_DIR/filter .
cp $DART_DIR/dart_to_cosmo .
cp $DART_DIR/cosmo_to_dart .
cp $DART_DIR/obs_seq.perfect obs_seq.out

numInst=2

for instance in {0..$(($numInst-1))}
do
  cd tsmp_instance_$instance
  ln -s ../input.nml .
  ln -s cosrst/lrff00030000o cosmo_prior
  ln -s cosout/lfff00030000.nc cosmo.nc
  ../cosmo_to_dart
  dartinstance=$(( $instance + 1 ))
  mv dart_prior ../filter_ics.000$dartinstance
  cd ..
done

export LD_LIBRARY_PATH=/daten01/z4/netcdf4.1.3_gfortran_4.2.1/lib/
ln -s tsmp_instance_0/cosrst/lrff00030000o cosmo_prior
ln -s tsmp_instance_0/cosout/lfff00030000.nc cosmo.nc

date
mpirun  -np 64  ./filter   >> log_file 2>> err_file
date

for instance in {0..$(($numInst-1))}
do
  cd tsmp_instance_$instance
  dartinstance=$(( $instance + 1 ))
  ln -s ../filter_restart.000$dartinstance dart_posterior 
  ../dart_to_cosmo
  cd ..
done

echo "ready" > ready.txt
exit 0

