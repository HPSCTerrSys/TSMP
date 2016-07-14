#!/bin/csh

# usage: ./archive_hist.csh $rundir $ensemble_size
#
echo "-------------------------------------------------------------------"
echo "Block 6: History files will be archived in the "$HOME" directory"
echo "-------------------------------------------------------------------"
echo " "

set rundir = $1
set ensemble_size = $2

set numInst = `echo "($ensemble_size - 1)" | bc`

foreach instance (`seq 0 $numInst`)
  set dir = tsmp_instance_$instance
  if ( -d ./$dir ) then
    echo Directory already exists! $dir
    echo Moving $dir to archive...
  else
    mkdir $dir
    cd $dir
    mkdir clmout
    mkdir pflout
    mv $rundir/cosout .
    cd cosout
    mkdir ivr
    mv lfff00000000c* ivr/
    cd ..
    cd clmout
    mv $rundir/clmoas.clm2.h0* .
    cd ..
    cd pflout
    mkdir satur
    mkdir press
    cd satur
    mv $rundir/rurlaf.out.satur* .
    rm *.dist
    cd ..
    cd press
    mv $rundir/rurlaf.out.press* .
    rm *.dist
    cd ..
    cd ..
    mv $rundir/timing_all .
    mv $rundir/YUTIMING .
    mv $rundir/YUPRMASS .
    mv $rundir/lnd.stdin .
    mv $rundir/lmrun_uc .
    mv $rundir/coup_oas.tcl .
    mv $rundir/rurlaf.out.log .
    cd ..
    chmod -R a=+rwx $dir
  endif
  echo "Processed output " $dir
end
exit 0
