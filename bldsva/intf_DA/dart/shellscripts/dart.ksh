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

rm filter_ics.*
rm P*Diag.nc
rm filter_res*
rm dart_log.*
rm obs_seq.out
rm log_file
rm err_file

cp $DART_DIR/input.nml .
cp $DART_DIR/filter .
cp $DART_DIR/dart_to_cosmo .
cp $DART_DIR/cosmo_to_dart .
cp $DART_DIR/obs_seq.perfect obs_seq.out


numInst=16

for instance in {0..$(($numInst-1))}
do
  cd tsmp_instance_$instance
  rm dart_log.*
  rm dart_posterior
  rm cosmo.nc
  rm cosmo_prior

  ln -s ../input.nml .
  ln -s cosrst/lrff01000000o cosmo_prior
  ln -s cosout/lfff01000000.nc cosmo.nc
  ../cosmo_to_dart || exit 1
  dartinstance=$(( $instance + 1 ))
  filterName=`printf ../filter_ics.%04d $dartinstance`
  mv dart_prior $filterName 
  cd ..
done

export LD_LIBRARY_PATH=/daten01/z4/netcdf4.1.3_gfortran_4.2.1/lib/
ln -s tsmp_instance_0/cosrst/lrff01000000o cosmo_prior
ln -s tsmp_instance_0/cosout/lfff01000000.nc cosmo.nc

date
mpirun  -np 64  ./filter || exit 1  >> log_file 2>> err_file
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

