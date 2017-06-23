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
 #parflow, manual recommends large odd number, not sure why?
 #without root distribution perturbation, simulation fails with inst 41
 
 set seedno = `echo "($instance*1000+13111)" | bc`
 if ($instance == 41) then
   echo "Changing seed no. for 41"
   set seedno = `echo "($instance*1000+11311)" | bc`
 endif  
 cd tsmp_instance_$instance
 sed "s,__seedno__,${seedno}," -i coup_oas.tcl
 tclsh coup_oas.tcl
 cd ..

 #clm
 #Perturb leaf c:n and root distribution
 cd tsmp_instance_${instance}
 cp ${clm_forcing_dir}/inputdata/lnd/clm2/pftdata/pft-physiology.c070207 .
 cp ${clm_forcing_dir}/perturb_surf/surfdata_${instance}_0014x0024.nc ./surfdata_0014x0024.nc

 set leafcn_def = `echo "(24.1+$ensemble_size*0.125)" | bc`
 set leafcn = `echo "($leafcn_def-$instance*0.25)" | bc`
 set roota  = `echo "(1.+$instance*0.22)" | bc`
 set rootb  = `echo "(1.+$instance*0.05)" | bc`

 echo " "
 sed "s,0.00000 25.0 0.100,0.00000 $leafcn 0.100," -i  pft-physiology.c070207
 #sed "s,80 -0.30  6.0 3.0 0.05000,80 -0.30 $roota $rootb 0.05000," -i  pft-physiology.c070207 
 cd ..

 #cosmo
 set turlength = `echo "(200.-$instance*2.5)" | bc`
 cd tsmp_instance_${instance}
 if ( -f lmrun_uc ) then
   sed "s,__turlen__,${turlength}," -i lmrun_uc
   ./lmrun_uc execluma
 endif
 cd ..
end

