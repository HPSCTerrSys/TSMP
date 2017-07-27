#!/bin/csh
# Script to setup perfect terrsymp run with cycle
# usage: ./perfect_tsmp_setup.csh cycle nrst map_fn machine
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
# paths 
set tsmpdir       = $HOME/terrsysmp
set archivedir    = "tsmp"
set shellpath     = $HOME/terrsysmp/bldsva/intf_DA/dart/shellscripts
set DART_DIR      = "$HOME/DART/lanai/models/terrsysmp/cosmo/work"
set clm_forcing_dir = "$HOME/database/idealRTD/clm/"
#
#User Settings End Here

set icycle = $1
set inrst  = $2
set map_fn = $3
set machine = $4

echo " "
echo "`date` -- BEGIN setup of terrsysmp perfect model"
echo " icycle: " $icycle
echo " "

#-------------------------------------------------------------------------
# Block 1, 
# Setup the terrsysmp version, reference setup to use, and ensemble numbers
# and the machine to use. Run setup for the initial run, and restart runs
#-------------------------------------------------------------------------

#
set defDir = $tsmpdir/bldsva/setups/$refsetup
cp $defDir/def/idealRTD_JURECA_setup_DA.ksh $defDir/idealRTD_JURECA_setup.ksh
cp $defDir/def/coup_oas_DA.tcl $defDir/coup_oas.tcl
cp $defDir/def/lnd.stdin_DA $defDir/lnd.stdin
cp $defDir/def/lmrun_uc5_1_DA $defDir/lmrun_uc5_1
#
set sdate = `printf perfectModel%02d $icycle`

if ($icycle == 1) then
    #
    echo "-------------------------------------------------------------------"
    echo "Block 1:  Setting up the inital terrsysmp"
    echo "-------------------------------------------------------------------"
    echo " "
    #
    cd $tsmpdir/bldsva
    ./setup_tsmp.ksh -v $tsmpver -V $refsetup -m $machine -I $sdate -r $WORK/
    set temp_dir      = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$sdate
    set rundir        = $WORK/$sdate
    #
    cd $rundir

     # Update Model States
     #ParFlow
     echo " Using the spinup parflow states ...." $map_fn
     cp $HOME/database/idealRTD/restart/tsmp_instance_${map_fn}/rurlaf.out.press.00096.pfb ./rur_ic_press.pfb
     tclsh ascii2pfb.tcl

     #cosmo
     set rasonum = `printf raso_IdealSnd_0000LT_%02d $map_fn`
     sed "s,raso_IdealSnd_0000LT.dat,$rasonum.dat," -i lmrun_uc
     rm cosmo_in/raso_IdealSnd_0000LT_*
     cp $HOME/database/idealRTD/cosmo/$rasonum.dat cosmo_in/
     ./lmrun_uc execluma

     #clm
     echo " Using the spinup clm states ...."
     cp $HOME/database/idealRTD/restart/tsmp_instance_${map_fn}/clmoas.clm2.r.2008-05-08-00000.nc ./clm_restart.nc

else if ($icycle > 1) then
    #
    echo "-------------------------------------------------------------------"
    echo "Block 2: Make the restart runs ...."
    echo "-------------------------------------------------------------------"
    echo " "
    #
    #Use the last rundir
    set oldicycle = `echo "($icycle - 1)" | bc` 
    set oldsdate  = `printf perfectModel%02d $oldicycle`
    set oldrundir = $WORK/$oldsdate
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

    ./setup_tsmp.ksh -v $tsmpver -V $refsetup -m $machine -I $sdate -r $WORK/ -s "$defaultInitDate[2]" -S "$defaultStartDate[2]"  -j "$clmrstfil" -k "$cosrstfil" -l "$pflrstfil"

    # Update the new rundir
    #
    set rundir = $WORK/$sdate
    rm $rundir/cosmo_in/raso_IdealSnd_0000LT_*   

else
    echo "ERROR : icycle < 1"
    exit 1
endif

#Perturb model Parameters

 cd $rundir
 #parflow, manual recommends large odd number, not sure why?
 #without root distribution perturbation, simulation fails with inst 41


 set seedno = `echo "($map_fn*1000+13111)" | bc`
 if ($map_fn == 41) then
   echo "Changing seed no. for 41"
   set seedno = `echo "($map_fn*1000+11311)" | bc`
 endif  

 sed "s,__seedno__,${seedno}," -i coup_oas.tcl
 tclsh coup_oas.tcl

 #clm
 #Perturb leaf c:n and root distribution
 cp ${clm_forcing_dir}/inputdata/lnd/clm2/pftdata/pft-physiology.c070207 .
 cp ${clm_forcing_dir}/perturb_surf/surfdata_${map_fn}_0014x0024.nc ./surfdata_0014x0024.nc

 set leafcn_def = `echo "(24.1+49.*0.125)" | bc`
 set leafcn = `echo "($leafcn_def-$map_fn*0.25)" | bc`
 set roota  = `echo "(1.+$map_fn*0.22)" | bc`
 set rootb  = `echo "(1.+$map_fn*0.05)" | bc`

 echo " " 
 sed "s,0.00000 25.0 0.100,0.00000 $leafcn 0.100," -i  pft-physiology.c070207
 #sed "s,80 -0.30  6.0 3.0 0.05000,80 -0.30 $roota $rootb 0.05000," -i  pft-physiology.c070207 

 #cosmo
 set turlength = `echo "(200.-$map_fn*2.5)" | bc`
 if ( -f lmrun_uc ) then
   sed "s,__turlen__,${turlength}," -i lmrun_uc
   ./lmrun_uc execluma
 endif
 cd ..

echo "-------------------------------------------------------------------"
echo "Block 4:  Updating the setup directory with default script"
echo "-------------------------------------------------------------------"
echo " "

cp $defDir/def/idealRTD_JURECA_setup_default.ksh $defDir/idealRTD_JURECA_setup.ksh
cp $defDir/def/coup_oas_default.tcl $defDir/coup_oas.tcl
cp $defDir/def/lnd.stdin_default $defDir/lnd.stdin
cp $defDir/def/lmrun_uc5_1_default $defDir/lmrun_uc5_1
 
exit 0
