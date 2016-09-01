#!/bin/csh

# usage: ./copy.csh $rundir 
#
echo "-------------------------------------------------------------------"
echo "History files will be archived in the "$HOME" directory"
echo "usage: ./copy.csh rundir"
echo "-------------------------------------------------------------------"
echo " "

cd $HOME

set rundir = $1
set outputdir = $rundir:t

echo $outputdir
cd $HOME

  if ( -d ./$outputdir ) then
    echo Directory already exists! $outputdir
    exit 1
  else
    mkdir $outputdir
    cd $outputdir
    mkdir clmout
    mkdir pflout
    cp -r $rundir/cosout .
    cd cosout
    mkdir ivr
    cp lfff00000000c* ivr/
    cd ..
    cd clmout
    cp $rundir/clmoas.clm2.h0* .
    cd ..
    cd pflout
    mkdir satur
    mkdir press
    cd satur
    cp $rundir/rurlaf.out.satur* .
    rm *.dist
    cd ..
    cd press
    cp $rundir/rurlaf.out.press* .
    rm *.dist
    cd ..
    cd ..
    cp $rundir/timing_all .
    cp $rundir/YUTIMING .
    cp $rundir/YUPRMASS .
    cp $rundir/lnd.stdin .
    cp $rundir/lmrun_uc .
    cp $rundir/coup_oas.tcl .
    cp $rundir/rurlaf.out.log .
    cd ..
    chmod -R a=+rwx $outputdir
  endif
  echo "Processed output " $outputdir
exit 0
