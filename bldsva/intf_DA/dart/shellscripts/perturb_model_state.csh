#!/bin/csh
#
# Usage: ./perturb_model_state.csh rundir ensemble_szie
#-------------------------------------------------------------------------
# Block 2, 
# Create perturbation for the ensemble runs 
#-------------------------------------------------------------------------
echo "-------------------------------------------------------------------"
echo "Block 2: Perturbing the cosmo model states...."
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
 cd tsmp_instance_$instance
 cp $HOME/database/idealRTD/parflow/rur_ic_press_Sv06.pfb ./rur_ic_press.pfb
 tclsh ascii2pfb.tcl
 cd ..

 #cosmo
 cd tsmp_instance_$instance
 set rasonum = `printf raso_IdealSnd_0000LT_%02d $instance`
 sed "s,raso_IdealSnd_0000LT.dat,$rasonum.dat," -i lmrun_uc
 rm cosmo_in/raso_IdealSnd_0000LT_* 
 cp $HOME/database/idealRTD/cosmo/$rasonum.dat cosmo_in/
 ./lmrun_uc execluma
 cd ..
end

