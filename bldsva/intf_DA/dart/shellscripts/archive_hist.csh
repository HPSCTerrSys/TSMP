#!/bin/csh

# usage: ./archive_hist.csh $rundir $ensemble_size
#
echo "-------------------------------------------------------------------"
echo "History files will be archived in the "$HOME" directory"
echo "usage: ./archive_hist.csh rundir ensemble_size"
echo "-------------------------------------------------------------------"
echo " "

cd $HOME

set rundir = $1
set ensemble_size = $2
set outputdir = $rundir:t
set numInst = `echo "($ensemble_size - 1)" | bc`

echo $outputdir
cd $HOME
mkdir $outputdir
cd $outputdir
exit 0

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
    mv $rundir/$dir/cosout .
    cd cosout
    mkdir ivr
    mv lfff00000000c* ivr/
    cd ..
    cd clmout
    mv $rundir/$dir/clmoas.clm2.h0* .
    cd ..
    cd pflout
    mkdir satur
    mkdir press
    cd satur
    mv $rundir/$dir/rurlaf.out.satur* .
    rm *.dist
    cd ..
    cd press
    mv $rundir/$dir/rurlaf.out.press* .
    rm *.dist
    cd ..
    cd ..
    mv $rundir/$dir/timing_all .
    mv $rundir/$dir/YUTIMING .
    mv $rundir/$dir/YUPRMASS .
    mv $rundir/$dir/lnd.stdin .
    mv $rundir/$dir/lmrun_uc .
    mv $rundir/$dir/coup_oas.tcl .
    mv $rundir/$dir/rurlaf.out.log .
    cd ..
    chmod -R a=+rwx $dir
  endif
  echo "Processed output " $dir
end
exit 0
