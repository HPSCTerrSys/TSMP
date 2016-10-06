#!/bin/ksh

#Job Submission to Cluma
#PBS -N TerrSysMP-DART
#PBS -l walltime=00:90:00
#PBS -l nodes=1:ppn=64
#PBS -V 
#PBS -u pshrestha
#PBS -q batch
#PBS -o tsmp.log
#PBS -e tsmp.err


export LOGNAME="/home/pshrestha/terrsysmp/run/CLUMA2_1.2.0MCT_clm-cos-pfl_idealRTD_dart14"
export DART_DIR="/home/pshrestha/DART/lanai/models/terrsysmp/cosmo/work"
export LD_LIBRARY_PATH=/daten01/z4/netcdf4.1.3_gfortran_4.2.1/lib/

cd $LOGNAME

# Cleanup---------------
rm filter_ics.*
rm P*Diag.nc
rm filter_res*
rm dart_log.*
rm obs_seq.*
rm log_file
rm err_file
rm ready.txt
rm cosmo.nc
rm cosmo_prior

# Copy namelist and executable to rundirectory--------#TODO-
cp $DART_DIR/input.nml .
cp $DART_DIR/filter .
cp $DART_DIR/dart_to_cosmo .
cp $DART_DIR/cosmo_to_dart .

numInst=48

for instance in {0..$(($numInst-1))}
do
  cd tsmp_instance_$instance
  rm dart_log.*
  rm dart_posterior
  rm cosmo.nc
  rm cosmo_prior

  ln -s ../input.nml .
  cosrst=`ls -1 cosrst/lrff* | tail -n -1`
  cosout=`ls -1 cosout/lfff* | tail -n -1`
  ln -s $cosrst cosmo_prior
  ln -s $cosout cosmo.nc
  ../cosmo_to_dart || exit 1

  dartinstance=$(( $instance + 1 ))
  filterName=`printf ../filter_ics.%04d $dartinstance`
  mv dart_prior $filterName 
  cd ..
done

wait

# Grab the date/time of the cosmo state so we know which 
# observation sequence file to use.
timefile_path="tsmp_instance_0"
cosrbin=`grep cosrbin $timefile_path/cosmo_prior_time.txt | cut -d' ' -f2`
coshist=`grep coshist $timefile_path/cosmo_prior_time.txt | cut -d' ' -f2`
clmext=`grep clmext $timefile_path/cosmo_prior_time.txt | cut -d' ' -f2`
dIndat=`grep defaultInitDate $timefile_path/cosmo_prior_time.txt | cut -d' ' -f2`

echo $clmext $cosrbin $coshist

ln -s $DART_DIR/obs_seq.$dIndat obs_seq.out || exit 2

ln -s tsmp_instance_0/cosrst/$cosrbin cosmo_prior
ln -s tsmp_instance_0/cosout/$coshist cosmo.nc

echo "CPS 2"
date
mpirun  -np 64  ./filter || exit 3  >> log_file 2>> err_file
date

for instance in {0..$(($numInst-1))}
do
  cd tsmp_instance_$instance
  dartinstance=$(( $instance + 1 ))
  filterRestartName=`printf ../filter_restart.%04d $dartinstance`
  ln -s $filterRestartName dart_posterior 
  ../dart_to_cosmo || exit 4
  cd ..
done

wait

echo "ready" > ready.txt
exit 0

