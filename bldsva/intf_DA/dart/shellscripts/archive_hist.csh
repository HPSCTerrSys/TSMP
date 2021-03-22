#!/bin/csh

# usage: ./archive_hist.csh $rundir $ensemble_size
#
set cmd = cp 
echo "-------------------------------------------------------------------"
echo "History files will be archived in the "$WORK" directory"
echo "usage: ./archive_hist.csh rundir ensemble_size"
echo "-------------------------------------------------------------------"
echo " "

set rundir = $1
set ensemble_size = $2
set outputdir = $rundir:t
set numInst = `echo "($ensemble_size - 1)" | bc`

echo $outputdir
cd $WORK
if ( -d ./TEMP ) then
  echo TEMP exists
else
  mkdir TEMP
endif
cd TEMP
mkdir $outputdir
cd $outputdir

foreach instance (`seq 0 $numInst`)
  set dir = tsmp_instance_$instance
  if ( -d ./$dir ) then
    echo Directory already exists! $dir
    exit 1
  else
    mkdir $dir
    cd $dir
    mkdir clmout
    mkdir pflout
    ${cmd} -rf $rundir/$dir/cosout .
    cd cosout
    mkdir ivr
    #always mv here
    mv lfff00000000c* ivr/
    cd ..
    cd clmout
    ${cmd} $rundir/$dir/clmoas.clm2.h0* .
    cd ..
    cd pflout
    mkdir satur
    mkdir press
    cd satur
    ${cmd} $rundir/$dir/rurlaf.out.satur* .
    rm *.dist
    cd ..
    cd press
    ${cmd} $rundir/$dir/rurlaf.out.press* .
    rm *.dist
    cd ..
    cd ..
    ${cmd} $rundir/$dir/timing_all .
    ${cmd} $rundir/$dir/YUTIMING .
    ${cmd} $rundir/$dir/YUPRMASS .
    ${cmd} $rundir/$dir/lnd.stdin .
    ${cmd} $rundir/$dir/lmrun_uc .
    ${cmd} $rundir/$dir/coup_oas.tcl .
    ${cmd} $rundir/$dir/rurlaf.out.log .
    cd ..
    chmod -R a=+rwx $dir
  endif
  echo "Processed output " $dir
end
exit 0
