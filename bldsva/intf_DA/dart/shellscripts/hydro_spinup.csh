#!/bin/csh
# Script to create perturbed spinup for the ensemles
# usage: ./hydro_spinup.csh machine
# machine = JURECA, CLUMA2
#
#-------------------------------------------------------------------------
# User Settings
# 
#-------------------------------------------------------------------------
# experiment setup flags
set tsmpver       = "1.2.0MCT"
set refsetup      = "idealRTD"
# number of ensemble + 1
set ensemble_size = 49
# wall clock
set wc            = "02:00:00"
# paths 
set tsmpdir       = $HOME/terrsysmp
set archivedir    = "tsmp"
set shellpath     = $HOME/terrsysmp/bldsva/intf_DA/dart/shellscripts
set sdate         = "TBspinup"
#
#User Settings End Here

set machine       = $1

echo " "
echo "`date` -- BEGIN setup of spinup for .."
echo " "

#-------------------------------------------------------------------------
# Block 1, 
# Setup the terrsysmp version, 
#-------------------------------------------------------------------------

echo "-------------------------------------------------------------------"
echo "Block 1:  Setting up the inital terrsysmp"
echo "-------------------------------------------------------------------"
echo " "
#
set defDir = $tsmpdir/bldsva/setups/$refsetup
cp $defDir/def/idealRTD_JURECA_setup_TB.ksh $defDir/idealRTD_JURECA_setup.ksh
cp $defDir/def/coup_oas_TB.tcl $defDir/coup_oas.tcl
cp $defDir/def/lnd.stdin_TB $defDir/lnd.stdin

cd $tsmpdir/bldsva
./setup_tsmp.ksh -v $tsmpver -c "clm-pfl" -V $refsetup -m $machine -I $sdate -Q $wc -N $ensemble_size -r $WORK/run
set rundir         = $WORK/run$sdate

echo "-------------------------------------------------------------------"
echo "Block 2:  Perturb seed number and leafcn"
echo "-------------------------------------------------------------------"
echo " "

$shellpath/perturb_model_param.csh $rundir $ensemble_size

echo "-------------------------------------------------------------------"
echo "Block 3:  Updating the setup directory with default script"
echo "-------------------------------------------------------------------"
echo " "

# Update the setup directory again
cp $defDir/def/idealRTD_JURECA_setup_default.ksh $defDir/idealRTD_JURECA_setup.ksh
cp $defDir/def/coup_oas_default.tcl $defDir/coup_oas.tcl
cp $defDir/def/lnd.stdin_default $defDir/lnd.stdin

exit 0


