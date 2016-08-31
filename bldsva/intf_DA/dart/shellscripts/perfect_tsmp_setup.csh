#!/bin/csh
# Script to setup perfect terrsymp run with cycle
# usage: ./perfect_tsmp_setup.csh cycle nrst
#  cycle  = 1 for initial run
#  cycle  > 1 for restart run
#  nrst   = 0 for normal restart
#  nrst   = 1 for restart with assimilation
# Input needed are 
# --- defaultStartDate = "2008-05-08 00"
# --- defaultInitDate  = "2008-05-09 00"
# --- clmext  = "2008-05-09-00000"
# --- cosrbin  = "lrff01000000o"
# --- pflhist  = "00024" 
#
#-------------------------------------------------------------------------
# User Settings
# 
#-------------------------------------------------------------------------
# experiment setup flags
set tsmpver       = "1.2.0MCT"
set refsetup      = "idealRTD"
set machine       = "CLUMA2"
# paths 
set tsmpdir       = $HOME/terrsysmp
set archivedir    = "tsmp"
set shellpath     = $HOME/terrsysmp/bldsva/intf_DA/dart/shellscripts
set DART_DIR      = "/home/pshrestha/DART/lanai/models/terrsysmp/cosmo/work"
#
#User Settings End Here

echo " "
echo "`date` -- BEGIN setup of terrsysmp perfect model"
echo " "

set icycle = $1
#-------------------------------------------------------------------------
# Block 1, 
# Setup the terrsysmp version, reference setup to use, and ensemble numbers
# and the machine to use. Run setup for the initial run, and restart runs
#-------------------------------------------------------------------------

#foreach icycle (`seq 1 $runcycle`)


  set sdate = `printf perfectModel%02d $icycle`

  if ($icycle == 1) then
    #
    echo "-------------------------------------------------------------------"
    echo "Block 1:  Setting up the inital terrsysmp"
    echo "-------------------------------------------------------------------"
    echo " "
    #
    cd $tsmpdir/bldsva
    ./setup_tsmp.ksh -v $tsmpver -V $refsetup -m $machine -I $sdate
    set temp_dir      = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$sdate
    set rundir        = $tsmpdir"/run/"$temp_dir
    #
    rm $rundir/cosmo_in/raso_IdealSnd_0000LT_*
  else if ($icycle > 1) then
    #
    echo "-------------------------------------------------------------------"
    echo "Block 2: Make the restart runs ....TODO dates !!!!"
    echo "-------------------------------------------------------------------"
    echo " "
    #
    #Use the last rundir
    set oldicycle = `echo "($icycle - 1)" | bc` 
    set oldsdate  = `printf perfectModel%02d $oldicycle`
    set temp_dir  = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$oldsdate
    set oldrundir = $tsmpdir"/run/"$temp_dir
    #
    # Date manager
    cd $oldrundir
    cp $DART_DIR/input.nml .
    cp $DART_DIR/filter .
    cp $DART_DIR/cosmo_to_dart .

    set cosrst = `ls -1 cosrst/lrff* | tail -n -1`
    set cosout = `ls -1 cosout/lfff* | tail -n -1`
    ln -s $cosrst cosmo_prior
    ln -s $cosout cosmo.nc
    ./cosmo_to_dart || exit 1
    # Cleanup folder
    rm cosmo_prior cosmo.nc input.nml filter cosmo_to_dart dart_prior dart_log.out dart_log.nml
    cd
    #
    set timefile_path    = $oldrundir
    set defaultStartDate = `grep defaultStartDate $timefile_path/cosmo_prior_time.txt`
    set defaultInitDate  = `grep defaultInitDate $timefile_path/cosmo_prior_time.txt`
    set clmext = `grep clmext $timefile_path/cosmo_prior_time.txt`
    set cosrbin = `grep cosrbin $timefile_path/cosmo_prior_time.txt`
    set pflhist = `grep pflhist $timefile_path/cosmo_prior_time.txt`
    #
    #
    set clmrstfil = "$oldrundir/clmoas.clm2.r.$clmext[2].nc"
    set cosrstfil = "$oldrundir/cosrst/$cosrbin[2]"
    set pflrstfil = "$oldrundir/rurlaf.out.press.$pflhist[2].pfb"

    cd $tsmpdir/bldsva

    ./setup_tsmp.ksh -v $tsmpver -V $refsetup -m $machine -I $sdate -s "$defaultInitDate[2]" -S "$defaultStartDate[2]"  -j "$clmrstfil" -k "$cosrstfil" -l "$pflrstfil"

    # Update the new rundir
    set temp_dir      = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$sdate
    #
    set rundir = $tsmpdir"/run/"$temp_dir
    rm $rundir/cosmo_in/raso_IdealSnd_0000LT_*   

  else
    echo "ERROR : icycle < 1"
    exit 1
  endif
#end

exit 0


