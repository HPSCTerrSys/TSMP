#!/bin/csh
#
# Usage: ./perturb_model_state.csh rundir ensemble_szie
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
set ensemble_size = $2

cd $rundir

set numInst = `echo "($ensemble_size - 1)" | bc`

foreach instance (`seq 0 $numInst`)
 echo " "
 echo "tsmp_instance_"$instance
 echo " "
 cd tsmp_instance_${instance}
# cp $HOME/database/idealRTD/parflow/rur_ic_press_Sv06.pfb ./rur_ic_press.pfb
 echo " Using the spinup parflow states ...."
 cp $HOME/database/idealRTD/restart/tsmp_instance_${instance}/rurlaf.out.press.00096.pfb ./rur_ic_press.pfb
 tclsh ascii2pfb.tcl
 cd ..

 #cosmo
 cd tsmp_instance_${instance}
 set rasonum = `printf raso_IdealSnd_0000LT_%02d $instance`
 sed "s,raso_IdealSnd_0000LT.dat,$rasonum.dat," -i lmrun_uc
 rm cosmo_in/raso_IdealSnd_0000LT_* 
 cp $HOME/database/idealRTD/cosmo/$rasonum.dat cosmo_in/
 ./lmrun_uc execluma
 cd ..

 #clm
 cd tsmp_instance_${instance}
 echo " Using the spinup clm states ...."
 cp $HOME/database/idealRTD/restart/tsmp_instance_${instance}/clmoas.clm2.r.2008-05-08-00000.nc ./clm_restart.nc
 cd ..
end

