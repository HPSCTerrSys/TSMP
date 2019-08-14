#!/bin/ksh

#SBATCH --ignore-pbs
#SBATCH --job-name="dartDATE"
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --output=mpiMPMD-out.%j
#SBATCH --error=mpiMPMD-err.%j
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --mail-type=ALL
#SBATCH --account="hbn33"

#PBS -N dartCOS 
#PBS -l walltime=2:30:00
#PBS -l nodes=1:ppn=64
#PBS -V 
#PBS -u pshrestha 
#PBS -q batch
#PBS -e err.txt
#PBS -o out.txt
 
export LOGNAME="$WORK/$1"
export DART_DIR="$HOME/DART/lanai/models/terrsysmp/cosmo/work"
#Machine Specific
export LD_LIBRARY_PATH="$EBROOTNETCDFMINFORTRAN/lib/":$LD_LIBRARY_PATH
cd $LOGNAME
#Machine Specific
source $LOGNNAME/loadenvs

export numInst=$2
# Cleanup---------------
rm input.nml
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
cp $DART_DIR/cosmo_to_dart .


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

echo "ready" > ready.txt
exit 0

