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
set clm_forcing_dir = "/homea/hbn33/hbn331/database/idealRTD/clm/"

cd $rundir

set numInst = `echo "($ensemble_size - 1)" | bc`

foreach instance (`seq 0 $numInst`)
 echo " "
 echo "tsmp_instance_"$instance
 echo " "
 #parflow
 set seedno = `echo "($instance*100)" | bc`
 cd tsmp_instance_$instance
 sed "s,__seedno__,${seedno}," -i coup_oas.tcl
 tclsh coup_oas.tcl
 cd ..

 #clm
 cd tsmp_instance_${instance}
 cp ${clm_forcing_dir}/inputdata/lnd/clm2/pftdata/pft-physiology.c070207 .
 cp ${clm_forcing_dir}/perturb_surf/surfdata_${instance}_0014x0024.nc ./surfdata_0014x0024.nc

 set leafcn_def = `echo "(24.1+$ensemble_size*0.125)" | bc`
 set leafcn = `echo "($leafcn_def-$instance*0.25)" | bc`

 echo "TODO, I ADD DECIMAL HERE >>>"
 echo " "
 sed "s,0.00000 25.0 0.100,0.00000 $leafcn 0.100," -i  pft-physiology.c070207
 cd ..

 #cosmo
 #set turlength = `echo "(400.-$instance*20.)" | bc`
 #cd tsmp_instance_$instance
 #sed "s,tur_len=150,tur_len=$turlength," -i lmrun_uc
 #./lmrun_uc execluma
 #cd ..
end

