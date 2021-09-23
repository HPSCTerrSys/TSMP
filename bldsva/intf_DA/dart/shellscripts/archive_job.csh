#!/bin/csh

#SBATCH --job-name="dart_archive"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=mpiMPMD-out.%j
#SBATCH --error=mpiMPMD-err.%j
#SBATCH --time=01:00:00
#SBATCH --partition=esm
#SBATCH --mail-type=ALL
#SBATCH --account="hbn33"

foreach instance (`seq 50 90`)
  set model_dir = `printf $WORK/rundart%02d $instance`
 ./archive_hist.csh $model_dir 49
end
#set instance = $1
#set model_dir = `printf $WORK/rundart%02d $instance`
#echo Archiving $model_dir 
#./archive_hist.csh $model_dir 49
#exit 0
