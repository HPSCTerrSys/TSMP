#! /bin/ksh

#SBATCH --job-name="TSMPLoopCntr"
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --output=TSMPLoopCntr-out.%j
#SBATCH --error=TSMPLoopCntr-err.%j
#SBATCH --time=00:05:00
#SBATCH --partition=batch
#SBATCH --mail-type=ALL


#JOBCHAIN NAMELIST
SPATH="$HOME/terrsysmp/bldsva/intf_DA/dart/shellscripts/"
MACHINE="JUWELS"        #(which machine are your running on)
NUMCYCLE=120             #(number of days to run , number of JOBS = 2*$numCycle - 1)
NRST=1                  #, 1 or 2 or 3  (Which component to assimilate, 0: no assimilation, 1 cos, 2: clm, 3: parflow)
MAP_FN=9                #Mapping matrix (selected ensemble member for perfect model run)
RUNSFX="perfectModel"
ASSIMC=""
JOBSCRIPT0="jobchain_perfect.ksh "
JOBSCRIPT1="tsmp_slm_run.bsh"
#JOBCHAIN NAMELIST
#------------

if [[ $NRST == 1 ]] ; then ; ASSIMC="cosmo" ; fi
if [[ $NRST == 2 ]] ; then ; ASSIMC="clm" ; fi
if [[ $NRST == 3 ]] ; then ; ASSIMC="parflow" ; fi

JOBSCRIPT2="dart_"${ASSIMC}".ksh"
#------------
# ARGUMENTS
ICYCLE=$1
STEP=$2

#------------
cd $SPATH


if [[ $STEP == "run" ]] then
    echo " TerrSysMP run ...."
    echo " "
    # 1: (Script to create the rundirectory with multiple instances e.g. rundar01)
    if [[ $ICYCLE = "1" ]] then ; NRST=0 ; fi
    ./perfect_tsmp_setup.csh $ICYCLE $NRST $MAP_FN $MACHINE

    # 2: Submit jobscript for terrsysmp run
    RUNNAME=`printf $RUNSFX%02d ${ICYCLE}`
    echo "sbatch $WORK/$RUNNAME/$JOBSCRIPT1"
    JOBID=$(sbatch $WORK/$RUNNAME/$JOBSCRIPT1 2>&1 | awk '{print $(NF)}')
    echo $JOBID
    if [[  $(($NUMCYCLE - $ICYCLE)) == 0 ]] ; then ; exit 0 ; fi
    #
    echo "sbatch --dependency=afterok:${JOBID} $JOBSCRIPT0 $((ICYCLE+1)) run"
    JOBID=$(sbatch --dependency=afterok:${JOBID} $JOBSCRIPT0 $((ICYCLE+1)) "run" 2>&1 | awk '{print $(NF)}')
fi

if [[ $STEP == "filter" ]] then

    echo "Nothing here ..."
fi

exit 0
