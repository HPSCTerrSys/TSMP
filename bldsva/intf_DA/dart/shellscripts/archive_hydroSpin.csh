#!/bin/csh

# usage: ./archive_hydroSpin.csh $rundir $ensemble_size
# ensemble_size = ensemble_size + 1
# Change map_fn to the above size
echo "-------------------------------------------------------------------"
echo "History files will be archived in the "$HOME" directory"
echo "usage: ./archive_hydroSpin.csh rundir ensemble_size"
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
    #mv $rundir/$dir/cosout .
    #cd cosout
    #mkdir ivr
    #mv lfff00000000c* ivr/
    #cd ..
    cd clmout
    mv $rundir/$dir/clmoas.clm2.h0* .
    mkdir restart
    set clmrestart = `ls -1 $rundir/$dir/clmoas.clm2.r.* | tail -n -1`
    mv $clmrestart restart/
    cd ..
    cd pflout
    mkdir restart
    mkdir satur
    mkdir press
    cd satur
    mv $rundir/$dir/rurlaf.out.satur* .
    rm *.dist
    cd ..
    set pflrestart = `ls -1 $rundir/$dir/rurlaf.out.press*.pfb | tail -n -1`
    cp $pflrestart restart/
    mv $rundir/$dir/rurlaf.out.perm_* restart/
    cd press
    mv $rundir/$dir/rurlaf.out.press* .
    rm *.dist
    cd ..
    cd ..
    mv $rundir/$dir/timing_all .
    #mv $rundir/$dir/YUTIMING .
    #mv $rundir/$dir/YUPRMASS .
    mv $rundir/$dir/lnd.stdin .
    #mv $rundir/$dir/lmrun_uc .
    mv $rundir/$dir/coup_oas.tcl .
    mv $rundir/$dir/rurlaf.out.log .
    cd ..
    chmod -R a=+rwx $dir
  endif
  echo "Processed output " $dir
end
exit 0
