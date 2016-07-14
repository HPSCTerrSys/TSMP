#!/bin/csh
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

echo " "
echo "`date` -- BEGIN terrsysmp assimilation with dart"
echo " "

#User Settings End Here
#-------------------------------------------------------------------------
# Block 1, 
# Setup the terrsysmp version, reference setup to use, and ensemble numbers
# and the machine to use. Run setup for the initial run
#-------------------------------------------------------------------------
echo "-------------------------------------------------------------------"
echo "Block 1:  Setting up the terrsysmp"
echo "-------------------------------------------------------------------"
echo " "
#
set sdate         = "dart01"      #Later will be adapted for cycling
#
cd $tsmpdir/bldsva
./setup_tsmp.ksh -v $tsmpver -V $refsetup -m $machine -I $sdate -N $ensemble_size
set temp_dir      = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$sdate
set rundir        = $tsmpdir"/run/"$temp_dir

#
# Perturb
$shellpath/perturb.csh $rundir $ensemble_size

#-------------------------------------------------------------------------
# Block 3, 
# Integrate the ensemble terrsysmp run
#-------------------------------------------------------------------------
echo "-------------------------------------------------------------------"
echo "Block 3: Integrating terrsysmp"
echo "-------------------------------------------------------------------"
echo " "
#
#cd $rundir

#echo "`date` -- Start initial run"
#qsub tsmp_pbs_run.ksh
#while (1)
#  if (-f "ready.txt") then
#    echo "`date` -- END initial run"
#    echo " "
#    break
#  else
#    echo "Job is running..."
#    sleep 10 
#  endif
#end

#-------------------------------------------------------------------------
# Block 4, 
# dart assimilation 
#-------------------------------------------------------------------------
echo "-------------------------------------------------------------------"
echo "Block 4:  Code not written data assimilation ..."
echo " So we move ahead with the restart runs"
echo "-------------------------------------------------------------------"
echo " "

#-------------------------------------------------------------------------
# Block 5, 
# Restart runs 
#-------------------------------------------------------------------------

echo "-------------------------------------------------------------------"
echo "Block 5: Make the restart runs ....TODO dates !!!!"
echo "-------------------------------------------------------------------"
echo " "
#
#TODO
set defaultStartDate = "2008-05-08 00"
set defaultInitDate  = "2008-05-09 00"
set sdate = "dart02"

#Use the last rundir
#TODO
set clmrstfil = "$rundir/tsmp_instance_X/clmoas.clm2.r.2008-05-09-00000.nc"
set cosrstfil = "$rundir/tsmp_instance_X/cosrst/lrff01000000o"
set pflrstfil = "$rundir/tsmp_instance_X/rurlaf.out.press.00024.pfb"

cd $tsmpdir/bldsva

./setup_tsmp.ksh -v $tsmpver -V $refsetup -m $machine -I $sdate -N $ensemble_size -s "$defaultInitDate" -S "$defaultStartDate"  -j "$clmrstfil" -k "$cosrstfil" -l "$pflrstfil"

echo "./setup_tsmp.ksh -v 1.2.0MCT -V "idealRTD" -s "$defaultInitDate" -S "$defaultStartDate"  -j "$clmrstfil" -k "$cosrstfil" -l "$pflrstfil" "
echo " "
echo " Update the restart files for the corresponding instances ..."
echo " "
#
#CPS
set temp_dir      = $machine"_"$tsmpver"_clm-cos-pfl_"$refsetup"_"$sdate
set cosmo_rfile    = "lrff01000000o"   #TODO 
#
cd $tsmpdir"/run/"$temp_dir/
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
#cd $rundir
#echo "`date` -- Start restart run"
#qsub tsmp_pbs_run.ksh
#while (1)
#  if (-f "ready.txt") then
#    echo "`date` -- END restart run"
#    echo " "
#    break
#  else
#    echo "Job is running..."
#    sleep 10
#  endif
#end


exit 0
#-------------------------------------------------------------------------
# Block 6, 
# Archive the model and dart outputs ... 
#-------------------------------------------------------------------------
echo "-------------------------------------------------------------------"
echo "Block 6: History files will be archived in the "$HOME" directory"
echo "-------------------------------------------------------------------"
echo " "
set numInst = `echo "($ensemble_size - 1)" | bc`

foreach instance (`seq 0 $numInst`)
  set dir = tsmp_instance_$instance
  if ( -d ./$dir ) then
    echo Directory already exists! $dir
    echo Moving $dir to archive... 
  else
    mkdir $dir
    cd $dir
    mkdir clmout
    mkdir pflout
    mv $rundir/cosout .
    cd cosout
    mkdir ivr
    mv lfff00000000c* ivr/
    cd ..
    cd clmout
    mv $rundir/clmoas.clm2.h0* .
    cd ..
    cd pflout
    mkdir satur
    mkdir press
    cd satur
    mv $rundir/rurlaf.out.satur* . 
    rm *.dist
    cd ..
    cd press
    mv $rundir/rurlaf.out.press* .
    rm *.dist
    cd ..
    cd ..
    mv $rundir/timing_all .
    mv $rundir/YUTIMING .
    mv $rundir/YUPRMASS .
    mv $rundir/lnd.stdin .
    mv $rundir/lmrun_uc .
    mv $rundir/coup_oas.tcl .
    mv $rundir/rurlaf.out.log .
    cd ..
    chmod -R a=+rwx $dir
  endif
  echo "Processed output " $dir 
end

exit 0


