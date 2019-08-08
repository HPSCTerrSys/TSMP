#!/bin/csh
# Script to setup terrsymp runs with cycle
# usage: ./tsmp_setup.csh cycle nrst ens map_fn machine
#  cycle  = 1 for initial run
#  cycle  > 1 for restart run
#  nrst   = 0 for normal restart
#  nrst   = 1,2,3 for restart with cosmo,clm,parflow assimilation
#  ens    = ensemble size
#  map_fn = mapping of ensemble members
# machine = JURECA, CLUMA2
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
set testcasedir   = $TSMP_DATA/idealRTD 

#User Settings End Here

echo " "
echo "`date` -- BEGIN setup of terrsysmp for assimilation with dart"
echo " "

set icycle  = $1
set inrst   = $2
set ensemble_size = $3
set map_fn = $4
set machine = $5
if ($inrst <= 1) then
  set assimC = 'cosmo'
else if ($inrst == 2) then
  set assimC = 'clm'
else if ($inrst == 3) then
  set assimC = 'parflow'
endif
#-------------------------------------------------------------------------
# tsmp_setup: Block 1, 
# Setup the terrsysmp version, reference setup to use, and ensemble numbers
# and the machine to use. Run setup for the initial run, and restart runs
#-------------------------------------------------------------------------

#
set defDir = $tsmpdir/bldsva/setups/$refsetup
cp $defDir/def/idealRTD_${machine}_setup_DA.ksh $defDir/idealRTD_${machine}_setup.ksh
cp $defDir/def/coup_oas_DA.tcl $defDir/coup_oas.tcl
cp $defDir/def/lnd.stdin_DA $defDir/lnd.stdin
cp $defDir/def/lmrun_uc5_1_DA $defDir/lmrun_uc5_1

#foreach icycle (`seq 1 $runcycle`)


  set sdate = `printf dart%02d $icycle`

  if ($icycle == 1) then
    #
    echo "-------------------------------------------------------------------"
    echo "tsmp_setup: Block 1:  Setting up the inital terrsysmp"
    echo "-------------------------------------------------------------------"
    echo " "
    #
    cd $tsmpdir/bldsva
    # Added clmrstfil to run with CLM spinup states
    # Parflow spinup states is used directly as initial PFB file, no need to specify here
    set clmrstfil = "./clm_restart.nc"
    ./setup_tsmp.ksh -v $tsmpver -V $refsetup -m $machine -j "$clmrstfil" -I $sdate -N $ensemble_size -r $WORK/run --wtime=02:00:00
    set rundir         = $WORK/run$sdate
    #set temp_dir      = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$sdate
    #set rundir        = $tsmpdir"/run/"$temp_dir
    #
    # Perturb the model state
    $shellpath/perturb_model_state.csh $rundir $testcasedir $ensemble_size $map_fn
  else if ($icycle > 1) then
    #
    echo "-------------------------------------------------------------------"
    echo "tsmp_setup: Block 2 Make the restart runs ....TODO dates "
    echo "-------------------------------------------------------------------"
    echo " "
    #
    #Use the last rundir
    set oldicycle = `echo "($icycle - 1)" | bc` 
    set oldsdate  = `printf dart%02d $oldicycle`
    set oldrundir = $WORK/run$oldsdate
    #set temp_dir  = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$oldsdate
    #set oldrundir = $tsmpdir"/run/"$temp_dir
    #
    # Date manager, but this file is obtained only with dart job submission, so you have to submit
    # dart_....ksh anyway for open loop runs also
    set timefile_path = $oldrundir"/tsmp_instance_0"
    set defaultStartDate = `grep defaultStartDate $timefile_path/${assimC}_prior_time.txt`
    set defaultInitDate = `grep defaultInitDate $timefile_path/${assimC}_prior_time.txt`
    set clmext = `grep clmext $timefile_path/${assimC}_prior_time.txt`
    set cosrbin = `grep cosrbin $timefile_path/${assimC}_prior_time.txt`
    set pflhist = `grep pflhist $timefile_path/${assimC}_prior_time.txt`
    #
    #
    set cosrstfil = "$oldrundir/tsmp_instance_X/cosrst/$cosrbin[2]"

    if ($inrst == 2) then
      #Assimilation overwrites the same restart file ..
      set clmrstfil = "$oldrundir/tsmp_instance_X/clmoas.clm2.r.$clmext[2].nc" 
    else
      set clmrstfil = "$oldrundir/tsmp_instance_X/clmoas.clm2.r.$clmext[2].nc"
    endif

    if ($inrst == 3) then
      # for restart run with parflow assimilation
      set pflrstfil = "$oldrundir/tsmp_instance_X/parflow_restart.pfb"
    else
      #TODO for normal restart run without assimilation
      set pflrstfil = "$oldrundir/tsmp_instance_X/rurlaf.out.press.$pflhist[2].pfb"
    endif

    cd $tsmpdir/bldsva

    ./setup_tsmp.ksh -v $tsmpver -V $refsetup -m $machine -I $sdate -N $ensemble_size -s "$defaultInitDate[2]" -S "$defaultStartDate[2]"  -j "$clmrstfil" -k "$cosrstfil" -l "$pflrstfil" -r $WORK/run

    echo " "
    echo " Update the restart files for the corresponding instances ..."
    echo " "
    #
    #set temp_dir      = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$sdate
    #
    # Update the new rundir
    set rundir = $WORK/"run"$sdate
    cd $rundir
    #
    echo "Moving into " $WORK/"run"$sdate 

    set numInst = `echo "($ensemble_size - 1)" | bc`
    #
    foreach instance (`seq 0 $numInst`)
      cd  "tsmp_instance_"$instance
      #COSMO
      rm cosmo_in/raso_IdealSnd_0000LT_*
      if ($inrst == 1) then 
        # for restart run with cosmo assimilation
        rm cosmo_in/$cosrbin[2]
        ln -sf $oldrundir/"tsmp_instance_"$instance/cosmo_restart ./cosmo_in/$cosrbin[2]
      else
        #TODO for normal restart run without assimilation
        ln -sf $oldrundir/"tsmp_instance_"$instance/cosrst/$cosrbin[2] ./cosmo_in/$cosrbin[2]
      endif
      #CLM
      sed "s,tsmp_instance_X,tsmp_instance_$instance," -i lnd.stdin
      #ParFlow
      sed "s,tsmp_instance_X,tsmp_instance_$instance," -i coup_oas.tcl
      tclsh coup_oas.tcl
      if ($inrst == 3) then
        #Distribute assimilated parflow restart files
        sed "s,rur_ic_press.pfb,$oldrundir/tsmp_instance_$instance/parflow_restart.pfb," -i ascii2pfb.tcl
        tclsh ascii2pfb.tcl
      endif
      cd ..
    end
    #Move the prior and posterior dart debug files
    if ($inrst > 0) then
      cp $oldrundir/Posterior_Diag.nc $rundir/
      cp $oldrundir/Prior_Diag.nc $rundir/
    endif
  else
    echo "ERROR : icycle < 1"
    exit 1
  endif
#end

echo "-------------------------------------------------------------------"
echo "tsmp_setup: Block 2:  Perturb seed number and leafcn"
echo "-------------------------------------------------------------------"
echo " "

$shellpath/perturb_model_param.csh $rundir $ensemble_size $map_fn

echo "-------------------------------------------------------------------"
echo "tsmp_setup: Block 3:  Updating the setup directory with default script"
echo "-------------------------------------------------------------------"
echo " "

cp $defDir/def/idealRTD_${machine}_setup_default.ksh $defDir/idealRTD_${machine}_setup.ksh
cp $defDir/def/coup_oas_default.tcl $defDir/coup_oas.tcl
cp $defDir/def/lnd.stdin_default $defDir/lnd.stdin
cp $defDir/def/lmrun_uc5_1_default $defDir/lmrun_uc5_1

exit 0


