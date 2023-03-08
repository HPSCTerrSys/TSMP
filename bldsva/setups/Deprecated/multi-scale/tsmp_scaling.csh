#!/bin/csh
# Script for running multiple scale simulation ...
echo " "
echo "----------------------------------------"
echo "Usage: ./tsmp_scaling dxxx machine N    "
echo "----------------------------------------"
echo " "

set expID  = "multi-scale"
set nameID = $1
set mach   = $2
set rstID  = $3

cp $HOME/terrsysmp/bldsva/setups/$expID/$nameID/multi-scale_{$mach}_setup.ksh $HOME/terrsysmp/bldsva/setups/$expID/

cd $HOME/terrsysmp/bldsva

if ($rstID == 1) then
  ./setup_tsmp.ksh -v 1.2.0MCT -V $expID -c "clm-pfl" -m $mach -I "_"$nameID"_"$rstID
  exit 0
else
  set iniID = `echo "($rstID - 1)" | bc`
  set rstDate   = "2009-01-01 00"
  set clmfname  = "clmoas.clm2.r.2010-01-01-00000.nc"
  set pflfname  = "rurlaf.out.press.00073.pfb"
  #set rstDate   = "2009-10-08 00"
  #set clmfname  = "clmoas.clm2.r.2009-10-08-00000.nc"
  #set pflfname  = "rurlaf.out.press.00056.pfb"
  set clmrstfil = "/home/pshrestha/terrsysmp/run/"$mach"_1.2.0MCT_clm-pfl_multi-scale_"$nameID"_"$iniID/$clmfname
  set pflrstfil = "/home/pshrestha/terrsysmp/run/"$mach"_1.2.0MCT_clm-pfl_multi-scale_"$nameID"_"$iniID/$pflfname

  ./setup_tsmp.ksh -v 1.2.0MCT -V $expID -c "clm-pfl" -S "2009-01-01 00" -I "_"$nameID"_"$rstID -s "$rstDate"  -j "$clmrstfil" -l "$pflrstfil" -m $mach 
  exit 0
endif
