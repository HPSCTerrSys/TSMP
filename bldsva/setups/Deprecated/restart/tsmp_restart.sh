#!/bin/bash

#SBATCH --job-name="TSMP"
#SBATCH --nodes=6
#SBATCH --ntasks=288
#SBATCH --ntasks-per-node=48
#SBATCH --output=mpiMPMD-out.%j
#SBATCH --error=mpiMPMD-err.%j
#SBATCH --time=00:20:00
#SBATCH --partition=devel
#SBATCH --mail-type=NONE
#SBATCH --account=slts


startDate=__startDate_bldsva__
initDate=__initDate_bldsva__
dt_clm=__dt_clm_bldsva__
dt_cosmo=__dt_cos_bldsva__
start_h=`echo $startDate|cut -b 12-`
start_sec=$(($start_h*3600))

hstart=0
Tcycle=2
numHours=12
cycle_h=$(($numHours/$Tcycle))

rundir=`pwd`
source $rundir/loadenvs
export PARFLOW_DIR=__PARFLOW_DIR__
cd $rundir

for cycle in $(seq 1 1 $Tcycle) ; # CHANGE LINES 53++ -> dependent jobs
do

cur_cycle=${cycle}
##############################################################
# Modifying COSMO namelists
##############################################################

nstop=$(($cycle_h*$cur_cycle*3600/$dt_cosmo))

sed -i "s/hstart = .*,/hstart = ${hstart},/g" $rundir/INPUT_ORG
sed -i "s/nstop = .*,/nstop = ${nstop},/g" $rundir/INPUT_ORG
echo "hstart:: " $hstart
echo "nstop:: " $nstop
cp $rundir/INPUT_ORG $rundir/INPUT_ORG_cycle-${cur_cycle}
##############################################################
# Modifying CLM namelists
##############################################################
nelapse=$(($cycle_h*3600/$dt_clm +1 ))
sed -i "/start_tod/cstart_tod     = ${start_sec}" lnd.stdin
sed -i "/nelapse/cnelapse     = ${nelapse}" lnd.stdin

if [ "$cycle" -gt 1 ]
then

   sed -i "/nsrest/cnsrest     = 1" lnd.stdin
   startday=`echo $startDate | cut -c1-10`
   repl=" finidat      ='${rundir}/clmoas.clm2.r.${startday}-${start_sec}.nc'"
   sed -i "/finidat.*/c \\$repl" lnd.stdin

fi
cp lnd.stdin lnd.stdin_cycle-${cur_cycle}
##############################################################
# Modifying ParFlow TCL flags
##############################################################
pfl_stop=$(($cycle_h*$cycle))
sed -i "/pfset TimingInfo.StopTime/cpfset TimingInfo.StopTime ${pfl_stop}.0025" coup_oas.tcl

if [ "$cycle" -gt 1 ]
then
    ic_pressure=`ls -1rt cordex0.11.out.press.*.pfb | tail -1`
#    sed -i "/pfset Geom.domain.ICPressure.FileName/a pfset Geom.domain.ICPressure.FileName    ${ic_pressure}/" coup_oas.tcl
    sed -i "/pfset Geom.domain.ICPressure.Value.*/a pfset Geom.domain.ICPressure.FileName    ${ic_pressure}" coup_oas.tcl
    sed -i "/pfdist/cpfdist ${ic_pressure}" coup_oas.tcl
fi

tclsh coup_oas.tcl
cp coup_oas.tcl coup_oas.tcl_cycle-${cur_cycle}

cd $rundir

date
echo "started cycle = $cycle at $(date '+%Y-%m-%d %H:%M:%S')" >> started.txt
rm -rf YU*

srun --multi-prog slm_multiprog_mapping.conf
date
echo "ready cycle = $cycle at $(date '+%Y-%m-%d %H:%M:%S')" >> ready.txt

sim_hour=$(($cycle_h*$cur_cycle))
hstart=${sim_hour}
start_sec=$(($start_sec+$cycle_h*3600))

done
