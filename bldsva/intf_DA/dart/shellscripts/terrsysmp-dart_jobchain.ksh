#!/bin/ksh

#SBATCH --job-name="TSMPLoopCntr"
#SBATCH --nodes=1
#SBATCH --ntasks=48
#SBATCH --ntasks-per-node=48
#SBATCH --output=TSMPLoopCntr-out.%j
#SBATCH --error=TSMPLoopCntr-err.%j
#SBATCH --time=00:30:00
#SBATCH --account="hbn33"
#SBATCH --partition=batch
#SBATCH --mail-type=ALL


#JOBCHAIN NAMELIST
SPATH="$TRAINHOME/terrsysmp/bldsva/intf_DA/dart/shellscripts/"
DPATH="$TRAINHOME/DART/lanai/models/terrsysmp/"
MACHINE="JUWELS"        #(which machine are your running on)
COMPILER="Intel"        #which compiler
NUMCYCLE=30             #(number of days to run , number of JOBS = 2*$numCycle - 1)
NRST=0                  #, 1 or 2 or 3  (Which component to assimilate, 0: no assimilation, 1 cos, 2: clm, 3: parflow)
NENS=48                #Ensemble Size (max 49 based on spinup)
MAP_FN="$SPATH/map_fn.txt"  #Mapping matrix for ensemble runs
RUNSFX="rundart"
ASSIMC=""
JOBSCRIPT0="terrsysmp-dart_jobchain.ksh "
JOBSCRIPT1="tsmp_slm_run.bsh"
#JOBCHAIN NAMELIST
#------------

# Clean up space 
rm mpiMPMD*
rm TSMPLoopCntr*
 
# Need to run dart to get time output files
if [[ $NRST == 0 ]] ; then ; ASSIMC="date" ; fi
if [[ $NRST == 1 ]] ; then ; ASSIMC="cosmo" ; fi
if [[ $NRST == 2 ]] ; then ; ASSIMC="clm" ; fi
if [[ $NRST == 3 ]] ; then ; ASSIMC="parflow" ; fi

JOBSCRIPT2="dart_"${ASSIMC}".ksh"
#------------
# ARGUMENTS
ICYCLE=$1
STEP=$2

#------------
if [[ $ICYCLE == 1 ]] then
    echo " Updating DART namelist ...."
    echo " "
    ./dart_setup.csh $NENS $DPATH
fi

cd $SPATH

if [[ $STEP == "run" ]] then
    echo " TerrSysMP run ...."
    echo " "
    # 1: (Script to create the rundirectory with multiple instances e.g. rundar01)
    if [[ $ICYCLE = "1" ]] then ; NRST=0 ; fi
    ./tsmp_setup.csh $ICYCLE $NRST $NENS $MAP_FN $MACHINE $COMPILER

    # 2: Submit jobscript for terrsysmp run
    RUNNAME=`printf $RUNSFX%02d ${ICYCLE}`
    echo "sbatch $WORK/$RUNNAME/$JOBSCRIPT1"
    JOBID=$(sbatch $WORK/$RUNNAME/$JOBSCRIPT1 2>&1 | awk '{print $(NF)}')
    echo $JOBID
    if [[  $(($NUMCYCLE - $ICYCLE)) == 0 ]] ; then ; exit 0 ; fi

    echo "sbatch --dependency=afterok:${JOBID} $JOBSCRIPT0 $ICYCLE filter"
    JOBID=$(sbatch --dependency=afterok:${JOBID} $JOBSCRIPT0 $ICYCLE "filter" 2>&1 | awk '{print $(NF)}')
fi



if [[ $STEP == "filter" ]] then

    echo " Running DART filter ..."
    echo " "
    RUNNAME=`printf $RUNSFX%02d ${ICYCLE}`

    # 4: Submit jobscript for dart run
    echo "sbatch $JOBSCRIPT2 $RUNNAME $NENS"
    JOBID=$(sbatch $JOBSCRIPT2 $RUNNAME $NENS 2>&1 | awk '{print $(NF)}')
    #
    echo "sbatch --dependency=afterok:${JOBID} $JOBSCRIPT0 $(($ICYCLE+1)) run"
    JOBID=$(sbatch --dependency=afterok:${JOBID} $JOBSCRIPT0 $(($ICYCLE+1)) "run" 2>&1 | awk '{print $(NF)}')
fi

exit 0
