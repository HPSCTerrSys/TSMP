#!/bin/csh
# Script to setup terrsymp runs with cycle
# usage: ./tsmp_setup.csh cycle
#  cycle  = 1 for initial run
#  cycle  > 1 for restart run
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
set ensemble_size = 16 
set machine       = "CLUMA2"
# paths 
set tsmpdir       = $HOME/terrsysmp
set archivedir    = "tsmp"
set shellpath     = $HOME/terrsysmp/bldsva/intf_DA/dart/shellscripts
#
#User Settings End Here

echo " "
echo "`date` -- BEGIN setup of terrsysmp for assimilation with dart"
echo " "

set icycle = $1
#-------------------------------------------------------------------------
# Block 1, 
# Setup the terrsysmp version, reference setup to use, and ensemble numbers
# and the machine to use. Run setup for the initial run, and restart runs
#-------------------------------------------------------------------------

#foreach icycle (`seq 1 $runcycle`)


  set sdate = `printf dart_cycle_%02d $icycle`

  if ($icycle == 1) then
    #
    echo "-------------------------------------------------------------------"
    echo "Block 1:  Setting up the inital terrsysmp"
    echo "-------------------------------------------------------------------"
    echo " "
    #
    cd $tsmpdir/bldsva
    ./setup_tsmp.ksh -v $tsmpver -V $refsetup -m $machine -I $sdate -N $ensemble_size
    set temp_dir      = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$sdate
    set rundir        = $tsmpdir"/run/"$temp_dir
    #
    # Perturb
    $shellpath/perturb.csh $rundir $ensemble_size
  else if ($icycle > 1) then
    #
    echo "-------------------------------------------------------------------"
    echo "Block 2: Make the restart runs ....TODO dates !!!!"
    echo "-------------------------------------------------------------------"
    echo " "
    # Date manager
    set timefile_path = "./"
    set defaultStartDate = `grep defaultStartDate $timefile_path/cosmo_prior_time.txt`
    set defaultInitDate = `grep defaultInitDate $timefile_path/cosmo_prior_time.txt`
    set clmext = `grep clmext $timefile_path/cosmo_prior_time.txt`
    set cosrbin = `grep cosrbin $timefile_path/cosmo_prior_time.txt`
    set pflhist = `grep pflhist $timefile_path/cosmo_prior_time.txt`

#TODO

    #
    #Use the last rundir
    set oldicycle = `echo "($icycle - 1)" | bc` 
    set oldsdate  = `printf dart_cycle_%02d $oldicycle`
    set temp_dir  = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$oldsdate
    set oldrundir = $tsmpdir"/run/"$temp_dir
    #
    set clmrstfil = "$oldrundir/tsmp_instance_X/clmoas.clm2.r.$clmext[2].nc"
    set cosrstfil = "$oldrundir/tsmp_instance_X/cosrst/$cosrbin[2]"
    set pflrstfil = "$oldrundir/tsmp_instance_X/rurlaf.out.press.$pflhist[2].pfb"

    cd $tsmpdir/bldsva

    ./setup_tsmp.ksh -v $tsmpver -V $refsetup -m $machine -I $sdate -N $ensemble_size -s "$defaultInitDate[2]" -S "$defaultStartDate[2]"  -j "$clmrstfil" -k "$cosrstfil" -l "$pflrstfil"

    echo " "
    echo " Update the restart files for the corresponding instances ..."
    echo " "
    #
    set temp_dir      = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$sdate
    #
    # Update the new rundir
    set rundir = $tsmpdir"/run/"$temp_dir
    cd $rundir
    #
    echo "Moving into " $tsmpdir"/run/"$temp_dir/

    set numInst = `echo "($ensemble_size - 1)" | bc`
    #
    foreach instance (`seq 0 $numInst`)
      cd  "tsmp_instance_"$instance
      #COSMO
      ln -sf $oldrundir/"tsmp_instance_"$instance/cosrst/$cosrbin[2] ./cosmo_in/$cosrbin[2]
      #CLM
      sed "s,tsmp_instance_X,tsmp_instance_$instance," -i lnd.stdin
      #ParFlow
      sed "s,tsmp_instance_X,tsmp_instance_$instance," -i coup_oas.tcl
      tclsh coup_oas.tcl
      cd ..
    end

    # Perturb
    $shellpath/perturb.csh $rundir $ensemble_size
  else
    echo "ERROR : icycle < 1"
    exit 1
  endif
#end

exit 0


