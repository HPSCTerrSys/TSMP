#!/bin/csh
#
# Usage: ./perturb.csh rundir ensemble_szie
#-------------------------------------------------------------------------
# Block 2, 
# Create perturbation for the ensemble runs 
#-------------------------------------------------------------------------
echo "-------------------------------------------------------------------"
echo "Block 2: Perturbing the parameters of clm and cosmo ...."
echo "-------------------------------------------------------------------"
echo " "
#
set rundir = $1
set ensemble_size = $2

cd $rundir

set numInst = `echo "($ensemble_size - 1)" | bc`

foreach instance (`seq 0 $numInst`)
 echo " "
 echo "tsmp_instance_"$instance
 echo " "
 #parflow
 cd tsmp_instance_$instance
 cp /daten01/z4/database/ParFlow/idealRTD/rur_ic_press_Sv06.pfb ./rur_ic_press.pfb
 tclsh ascii2pfb.tcl
 cd ..

 #clm
 cd tsmp_instance_$instance
 cp /daten01/z4/database/clm3.5/inputdata/lnd/clm2/pftdata/pft-physiology.c070207 ./pft-physiology_$instance
 set leafcn = `echo "(40.-$instance*2.)" | bc`

 echo "TODO, I ADD DECIMAL HERE >>>"
 echo " "
 sed "s,0.00000 25.0 0.100,0.00000 $leafcn.0 0.100," -i  pft-physiology_$instance
 sed "s,fpftcon        = '/daten01/z4/database/clm3.5/inputdata/lnd/clm2/pftdata/pft-physiology.c070207',fpftcon        = './pft-physiology_$instance'," -i lnd.stdin
 cd ..

 #cosmo
 set turlength = `echo "(400.-$instance*20.)" | bc`
 cd tsmp_instance_$instance
 sed "s,tur_len=150,tur_len=$turlength," -i lmrun_uc
 ./lmrun_uc execluma
 cd ..
end

