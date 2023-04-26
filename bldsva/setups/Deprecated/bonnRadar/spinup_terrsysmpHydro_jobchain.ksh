#!/bin/ksh

#SBATCH --job-name="TerrSysMP"
#SBATCH --nodes=1
#SBATCH --account=hbn33
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=mpiMPMD-out.%j
#SBATCH --error=mpiMPMD-err.%j
#SBATCH --time=00:10:00
#SBATCH --partition=batch
#SBATCH --mail-type=NONE

JOBSCRIPT0="spinup_terrsysmpHydro_jobchain.ksh"
JOBSCRIPT1="tsmp_slm_run.bsh"
#JOBCHAIN NAMELIST
#------------

# ARGUMENTS
STEP=$1
module load Tcl/.8.6.8
cd /p/scratch/chbn33/hbn331//runHET

#------------
if [[ $STEP == "run" ]] then
  
  #Clear directory
  ./remove.ksh

  clmrst=`ls -1 clmoas.clm2.r.*.nc | tail -n -1`
  pflrst=`ls -1 rurlaf.out.press.*.pfb | tail -n -1`

  echo $clmrst
  #Parse date from CLM restart file
  dd=$(echo "$clmrst" | cut -d"-" -f3)
  mm=$(echo "$clmrst" | cut -d"-" -f2)
  prefx=$(echo "$clmrst" | cut -d"-" -f1)
  yyyy=$(echo "$prefx" | cut -d"." -f4)
  yyyymmdd=$yyyy$mm$dd
  echo $yyyymmdd

  echo $pflrst
  #Parse StartCount from ParFlow restart file
  sc=$(echo "$pflrst" | cut -d"." -f4)
  typeset -i sc 
  echo $sc

  #Update namelist
  cp namelist/* .
  sed "s,__yyyymmdd__,${yyyymmdd}," -i lnd.stdin
  sed "s,__restart__,${clmrst}," -i lnd.stdin
  #
  sed "s,__sc__,${sc}," -i coup_oas.tcl
  sed "s,__restart__,${pflrst}," -i coup_oas.tcl
  /usr/bin/tclsh coup_oas.tcl

  echo " TerrSysMP-HYDRO run ...."
  echo " "

  echo "sbatch $JOBSCRIPT1 "
  JOBID=$(sbatch $JOBSCRIPT1 2>&1 | awk '{print $(NF)}')
  echo $JOBID

  echo "sbatch --dependency=afterok:${JOBID} $JOBSCRIPT0 run "
  JOBID=$(sbatch --dependency=afterok:${JOBID} $JOBSCRIPT0 "run" 2>&1 | awk '{print $(NF)}')
fi

exit 0
