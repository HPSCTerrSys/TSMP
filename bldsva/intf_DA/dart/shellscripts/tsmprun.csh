#!/bin/csh
# Script to setup terrsymp runs with cycle
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
set runcycle      = 2 
set year          = 2008
set month         = 5
set day           = 8
set hour          = 0

set defaultStartDate = "2008-05-08 00"
#User Settings End Here

echo " "
echo "`date` -- BEGIN setup of terrsysmp for assimilation with dart"
echo " "

#-------------------------------------------------------------------------
# Block 1, 
# Setup the terrsysmp version, reference setup to use, and ensemble numbers
# and the machine to use. Run setup for the initial run, and restart runs
#-------------------------------------------------------------------------

#set defaultStartDate = `"2008-05-08 00"
foreach icycle (`seq 1 $runcycle`)


  set sdate = `printf dart_cycle_%02d $icycle`
  #TODO needed for restart
  set defaultStartDate = "2008-05-08 00"
  set defaultInitDate  = "2008-05-09 00"

  set clmext  = "2008-05-09-00000"
  set cosext  = "01000000o"
  set pflext  = "00024"
  #TODO

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
  else
    #
    echo "-------------------------------------------------------------------"
    echo "Block 2: Make the restart runs ....TODO dates !!!!"
    echo "-------------------------------------------------------------------"
    echo " "
    #
    #Use the last rundir
    set cosmo_rfile    = "lrff$cosext" 
    set clmrstfil = "$rundir/tsmp_instance_X/clmoas.clm2.r.$clmext.nc"
    set cosrstfil = "$rundir/tsmp_instance_X/cosrst/$cosmo_rfile"
    set pflrstfil = "$rundir/tsmp_instance_X/rurlaf.out.press.$pflext.pfb"

    cd $tsmpdir/bldsva

    ./setup_tsmp.ksh -v $tsmpver -V $refsetup -m $machine -I $sdate -N $ensemble_size -s "$defaultInitDate" -S "$defaultStartDate"  -j "$clmrstfil" -k "$cosrstfil" -l "$pflrstfil"

    echo " "
    echo " Update the restart files for the corresponding instances ..."
    echo " "
    #
    set temp_dir      = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$sdate
    #
    cd $tsmpdir"/run/"$temp_dir/
    #
    echo "Moving into " $tsmpdir"/run/"$temp_dir/

    set numInst = `echo "($ensemble_size - 1)" | bc`
    #
    foreach instance (`seq 0 $numInst`)
      cd  "tsmp_instance_"$instance
      #COSMO
      ln -sf $rundir/"tsmp_instance_"$instance/cosrst/$cosmo_rfile ./cosmo_in/$cosmo_rfile
      #CLM
      sed "s,tsmp_instance_X,tsmp_instance_$instance," -i lnd.stdin
      #ParFlow
      sed "s,tsmp_instance_X,tsmp_instance_$instance," -i coup_oas.tcl
      tclsh coup_oas.tcl
      cd ..
    end

    #Update run dir
    set rundir = $tsmpdir"/run/"$temp_dir

    # Perturb
    $shellpath/perturb.csh $rundir $ensemble_size
  endif
end

exit 0

#-------------------------------------------------------------------------
# Block 3 
# Archive the model and dart outputs ... 
#-------------------------------------------------------------------------
./archive_hist.csh $rundir $ensemble_size

exit 0


