#!/bin/csh
#
# Usage: ./perturb_model_state.csh rundir testcasedir ensemble_size map_fn
#-------------------------------------------------------------------------
# Block 2, 
# Create perturbation for the ensemble runs 
#-------------------------------------------------------------------------
echo "-------------------------------------------------------------------"
echo "Block 1: Perturbing the model states...."
echo "-------------------------------------------------------------------"
echo " "
#
set rundir = $1
set testcasedir = $2
set ensemble_size = $3
set map_fn_file = $4

#Read the map_fn_file
set map_fn =
foreach line (`cat ${map_fn_file}`)
  set map_fn = ( $map_fn $line )
end
#
cd $rundir

set numInst = `echo "($ensemble_size - 1)" | bc`

foreach instance (`seq 0 $numInst`)
 set sinst = `echo "($instance + 1)" | bc`
 set minstance = $map_fn[$sinst]
 echo " "
 echo "tsmp_instance_"$instance
 echo "Mapping..."$instance"->"$minstance
 echo " "
 #
 cd tsmp_instance_${instance}

 #IF ANY CHANGES MADE HERE< MAKE ALSO CHANGE IN perfect_model_setup###
 #-------------------------------------------------------------------#
 #ParFlow
 echo " Using the spinup parflow states ...." $instance " .." $minstance
 cp ${testcasedir}/restart/tsmp_instance_${minstance}/rurlaf.out.press.00096.pfb ./rur_ic_press.pfb
 tclsh ascii2pfb.tcl

 #cosmo
 set rasonum = `printf raso_IdealSnd_0000LT_%02d $minstance`
 sed "s,raso_IdealSnd_0000LT.dat,$rasonum.dat," -i lmrun_uc
 rm cosmo_in/raso_IdealSnd_0000LT_* 
 cp ${testcasedir}/cosmo/$rasonum.dat cosmo_in/
 ./lmrun_uc execluma

 #clm
 echo " Using the spinup clm states ...."
 cp ${testcasedir}/restart/tsmp_instance_${minstance}/clmoas.clm2.r.2008-05-08-00000.nc ./clm_restart.nc
 #-------------------------------------------------------------------#

 cd ..
end

