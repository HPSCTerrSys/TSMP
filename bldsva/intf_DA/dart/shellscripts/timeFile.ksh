#!/bin/ksh

export LOGNAME="/home/pshrestha/terrsysmp/run/CLUMA2_1.2.0MCT_clm-cos-pfl_idealRTD_dart13"
export DART_DIR="/home/pshrestha/DART/lanai/models/terrsysmp/cosmo/work"

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

numInst=1

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

echo "ready" > ready.txt
exit 0

