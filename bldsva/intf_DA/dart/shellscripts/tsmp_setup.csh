#!/bin/csh
# Script to setup terrsymp runs with cycle
# Usage: ./tsmp_setup.csh cycle
#  cycle  = 1 for initial run
#  cycle  > 1 for restart run
# Input needed are 
# --- defaultStartDate = "2008-05-08 00"
# --- defaultInitDate  = "2008-05-09 00"
# --- clmext  = "2008-05-09-00000"
# --- cosext  = "01000000o"
# --- pflext  = "00024" 
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

#TODO needed for restart
set defaultStartDate = "2008-05-08 00"
set defaultInitDate  = "2008-05-09 00"

set clmext  = "2008-05-09-00000"
set cosext  = "01000000o"
set pflext  = "00024"
#TODO

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
    #
    #Use the last rundir
    set oldicycle = `echo "($icycle - 1)" | bc` 
    set oldsdate  = `printf dart_cycle_%02d $oldicycle`
    set temp_dir  = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$oldsdate
    set oldrundir = $tsmpdir"/run/"$temp_dir
    #
    set cosmo_rfile    = "lrff$cosext" 
    set clmrstfil = "$oldrundir/tsmp_instance_X/clmoas.clm2.r.$clmext.nc"
    set cosrstfil = "$oldrundir/tsmp_instance_X/cosrst/$cosmo_rfile"
    set pflrstfil = "$oldrundir/tsmp_instance_X/rurlaf.out.press.$pflext.pfb"

    cd $tsmpdir/bldsva

    ./setup_tsmp.ksh -v $tsmpver -V $refsetup -m $machine -I $sdate -N $ensemble_size -s "$defaultInitDate" -S "$defaultStartDate"  -j "$clmrstfil" -k "$cosrstfil" -l "$pflrstfil"

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
      ln -sf $oldrundir/"tsmp_instance_"$instance/cosrst/$cosmo_rfile ./cosmo_in/$cosmo_rfile
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


