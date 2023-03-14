#! /bin/ksh


getMachineDefaults(){
route "${cyellow}>> getMachineDefaults${cnormal}"
  defaultMpiPath="$rootdir/lib/openmpi"
  defaultNcdfPath="$rootdir/lib/netcdf"
  defaultGribPath="$rootdir/lib/gribapi"
  defaultGribapiPath="$rootdir/lib/gribapi"
  defaultJasperPath=""
#  defaultTclPath="/usr/lib/x86_64-linux-gnu"
  defaultTclPath="$rootdir/lib/tcl"
  defaultHyprePath="$rootdir/lib/hypre"
  defaultSiloPath="$rootdir/lib/silo"
  hdf5path="$rootdir/lib/hdf5"
  defaultLapackPath=""
  defaultPncdfPath=""
  export PATH="$defaultTclPath/bin:$PATH"
  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$defaultSiloPath/lib"
  #echo "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
  #echo "tclsh is loaded "

  # Default Compiler/Linker optimization
  defaultOptC="-O2"

  profilingImpl=" no scalasca "
  profComp=""
  profRun=""
  profVar=""

  # Default Processor settings
  defaultwtime=""
  defaultQ=""

route "${cyellow}<< getMachineDefaults${cnormal}"
}

finalizeMachine(){
route "${cyellow}>> finalizeMachine${cnormal}"
route "${cyellow}<< finalizeMachine${cnormal}"
}


createRunscript(){
route "${cyellow}>> createRunscript${cnormal}"

comment "   copy x86 module load script into rundirectory"
  cp $rootdir/bldsva/machines/loadenv_x86 $rundir
check



mpitasks=$((numInst * ($nproc_cos + $nproc_clm + $nproc_pfl + $nproc_oas)))
nnodes=`echo "scale = 2; $mpitasks / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`

#DA
if [[ $withPDAF == "true" ]] ; then
  srun="mpiexec -n $mpitasks ./tsmp-pdaf -n_modeltasks $(($numInst-$startInst)) -filtertype 2 -subtype 1 -delt_obs $delta_obs -rms_obs 0.03 -obs_filename swc_crp"
else
  srun="mpiexec -np $nproc_cos ./lmparbin_pur : -np $nproc_pfl ./parflow $pflrunname : -np $nproc_clm ./clm"
fi

cat << EOF >> $rundir/tsmp_slm_run.bsh
#!/bin/bash

cd $rundir
source $rundir/loadenv_x86 
date
echo "started" > started.txt
rm -rf YU*
export $profVar
$srun
date
echo "ready" > ready.txt
exit 0

EOF





counter=0
#for mapfile
start_oas=$counter
end_oas=$(($start_oas+$nproc_oas-1))

start_cos=$(($nproc_oas+$counter))
end_cos=$(($start_cos+($numInst*$nproc_cos)-1))

start_pfl=$(($numInst*$nproc_cos+$nproc_oas+$counter))
end_pfl=$(($start_pfl+($numInst*$nproc_pfl)-1))

start_clm=$((($numInst*$nproc_cos)+$nproc_oas+($numInst*$nproc_pfl)+$counter))
end_clm=$(($start_clm+($numInst*$nproc_clm)-1))


if [[ $numInst > 1 &&  $withOASMCT == "true" ]] then
 for instance in {$startInst..$(($startInst+$numInst-1))}
 do
  for iter in {1..$nproc_cos}
  do
    if [[ $withCOS == "true" ]] then ; echo $instance >>  $rundir/instanceMap.txt ;fi
  done
 done
 for instance in {$startInst..$(($startInst+$numInst-1))}
 do
  for iter in {1..$nproc_pfl}
  do
    if [[ $withPFL == "true" ]] then ; echo $instance >>  $rundir/instanceMap.txt ;fi
  done
 done
 for instance in {$startInst..$(($startInst+$numInst-1))}
 do
  for iter in {1..$nproc_clm}
  do
    if [[ $withCLM == "true" ]] then ; echo $instance >>  $rundir/instanceMap.txt ;fi
  done
 done
fi

comment "   change permission of runscript and mapfile"
chmod 755 $rundir/tsmp_slm_run.bsh >> $log_file 2>> $err_file
check
route "${cyellow}<< createRunscript${cnormal}"
}

Npp=4

PFLProcXg=1
PFLProcYg=1
CLMProcXg=1
CLMProcYg=1
COSProcXg=1
COSProcYg=1
if [[ $refSetup == "cordex" ]] then
	PFLProcX=1
	PFLProcY=1
	CLMProcX=1
	CLMProcY=1
	COSProcX=1
	COSProcY=1
	elif [[ $refSetup == "nrw" ]] then
	PFLProcX=1
	PFLProcY=1
	CLMProcX=1
	CLMProcY=1
	COSProcX=1
	COSProcY=1
	elif [[ $refSetup == "idealRTD" ]] then
	PFLProcX=1
	PFLProcY=1
	CLMProcX=1
	CLMProcY=1
	COSProcX=1
	COSProcY=1
	else 
	PFLProcX=1
	PFLProcY=1
	CLMProcX=1
	CLMProcY=1
	COSProcX=1
	COSProcY=1
fi

export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$rootdir/lib/silo/lib"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$rootdir/lib/hdf5/lib"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$rootdir/lib/hypre/lib"
