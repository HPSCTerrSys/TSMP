#!/bin/ksh

#SBATCH --job-name="dartCLM"
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --output=mpiMPMD-out.%j
#SBATCH --error=mpiMPMD-err.%j
#SBATCH --time=00:10:00
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
 
export LOGNAME="$WORK/$1"
export LOGNAME_S="$WORK/rundart01"
export DART_DIR="$HOME/DART/lanai/models/terrsysmp/clm/work"
export LD_LIBRARY_PATH="$EBROOTNETCDFMINFORTRAN/lib/":$LD_LIBRARY_PATH
cd $LOGNAME
source $LOGNNAME/loadenvs

# Cleanup---------------
rm input.nml
rm clm_history.nc
rm clm_restart.nc
rm filter_ics.*
rm P*Diag.nc
rm filter_res*
rm dart_log.*
rm obs_seq.*
rm log_file
rm err_file
rm ready.txt
rm clm_restart.nc
rm clm_restart_s.nc

# Copy namelist and executable to rundirectory--------#TODO-
cp $DART_DIR/input.nml .
#cp $DART_DIR/pfidb_dz .
cp $DART_DIR/filter .
cp $DART_DIR/dart_to_clm .
cp $DART_DIR/clm_to_dart .

numInst=48

for instance in {0..$(($numInst-1))}
do
  cd tsmp_instance_$instance

  rm input.nml
  rm clm_restart.nc
  rm clm_history.nc 
  rm clm_restart_s.nc
  rm dart_log.*
  rm dart_posterior
  rm clm_prior_time.txt 
  rm dart_posterior_times.txt

  ln -s ../input.nml .
  clmhst=`ls -1 clmoas.clm2.h0.*.nc | tail -n -1`
  clmrst=`ls -1 clmoas.clm2.r.*.nc | tail -n -2 | head -n 1`
  clmrst_s=`ls -1 $LOGNAME_S/tsmp_instance_${instance}/clmoas.clm2.r.*.nc | tail -n -2 | head -n 1`
  #
  ln -s $clmhst clm_history.nc
  ln -s $clmrst clm_restart.nc
  ln -s $clmrst_s clm_restart_s.nc
  ../clm_to_dart || exit 1

  dartinstance=$(( $instance + 1 ))
  filterName=`printf ../filter_ics.%04d $dartinstance`
  mv dart_prior $filterName 
  cd ..
done

wait

# GRAB PARFLOW TIME FOR ASSIMILATION
timefile_path="tsmp_instance_0"
dIndat=`grep defaultInitDate $timefile_path/clm_prior_time.txt | cut -d' ' -f2`

ln -s $DART_DIR/obs_seq.$dIndat obs_seq.out || exit 2

#for filter?
ln -s tsmp_instance_0/$clmhst clm_history.nc
ln -s tsmp_instance_0/$clmrst clm_restart.nc
ln -s tsmp_instance_0/$clmrst clm_restart_s.nc
echo "CPS 2"
#for filter?

date
srun  -n 48  ./filter || exit 3  >> log_file 2>> err_file
date

for instance in {0..$(($numInst-1))}
do
  cd tsmp_instance_$instance
  dartinstance=$(( $instance + 1 ))
  filterRestartName=`printf ../filter_restart.%04d $dartinstance`
  ln -s $filterRestartName dart_posterior 
  ../dart_to_clm || exit 4
  cd ..
done

wait

echo "ready" > ready.txt
exit 0


