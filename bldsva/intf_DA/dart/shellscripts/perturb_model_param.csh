#!/bin/csh
#
# Usage: ./perturb_model_param rundir ensemble_szie map_fn
#-------------------------------------------------------------------------
# Create perturbation for the ensemble runs 
#-------------------------------------------------------------------------
echo "-------------------------------------------------------------------"
echo "perturb_model_param: Perturbing the parameters of terrsysmp ...."
echo "-------------------------------------------------------------------"
echo " "
#
set rundir = $1
set ensemble_size = $2
set map_fn_file = $3
set testcasedir = "$TSMP_DATA/idealRTD/"
# set clm_forcing_dir = "/daten01/z4/database/TestCases/idealRTD/clm/"
#Read the map_fn_file
set map_fn =
foreach line (`cat ${map_fn_file}`)
  set map_fn = ( $map_fn $line )
end

cd $rundir

set numInst = `echo "($ensemble_size - 1)" | bc`

foreach instance (`seq 0 $numInst`)
 set sinst = `echo "($instance + 1)" | bc`
 set minstance = $map_fn[$sinst]
 echo " "
 echo "tsmp_instance_"$instance
 echo "Mapping..."$instance"->"$minstance
 echo " "

 cd tsmp_instance_$instance
 #
 #IF ANY CHANGES MADE HERE< MAKE ALSO CHANGE IN perfect_model_setup###
 #-------------------------------------------------------------------#

 #parflow, manual recommends large odd number, not sure why?
 #without root distribution perturbation, simulation fails with inst 41

 set seedno = `echo "($minstance*1000+13111)" | bc`
 if ($minstance == 41) then
   echo "Changing seed no. for 41"
   set seedno = `echo "($minstance*1000+11311)" | bc`
 endif  

 sed "s,__seedno__,${seedno}," -i coup_oas.tcl
 tclsh coup_oas.tcl

 #clm
 #Perturb leaf c:n and root distribution
 cp ${testcasedir}/clm/inputdata/lnd/clm2/pftdata/pft-physiology.c070207 .
 cp ${testcasedir}/clm/perturb_surf/surfdata_${minstance}_0014x0024.nc ./surfdata_0014x0024.nc

 set leafcn_def = `echo "(24.1+49.*0.125)" | bc`
 set leafcn = `echo "($leafcn_def-$minstance*0.25)" | bc`
 set roota  = `echo "(1.+$minstance*0.22)" | bc`
 set rootb  = `echo "(1.+$minstance*0.05)" | bc`

 echo " "
 sed "s,0.00000 25.0 0.100,0.00000 $leafcn 0.100," -i  pft-physiology.c070207
 #sed "s,80 -0.30  6.0 3.0 0.05000,80 -0.30 $roota $rootb 0.05000," -i  pft-physiology.c070207 

 #cosmo
 echo "no perturbation of turbulennt length scale"
# set turlength = `echo "(200.-$minstance*2.5)" | bc`
 set turlength = 100.
 if ( -f lmrun_uc ) then
   sed "s,__turlen__,${turlength}," -i lmrun_uc
   ./lmrun_uc execluma
 else
   echo " No cosmo ...."
 endif
 #-------------------------------------------------------------------#
 #
 cd ..

end

