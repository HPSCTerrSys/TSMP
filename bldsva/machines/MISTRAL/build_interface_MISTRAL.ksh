#! /bin/ksh


getMachineDefaults(){

  comment "   init lmod functionality"
  source /sw/rhel6-x64/etc/profile.mistral
  check
  comment "   source and load Modules on MISTRAl compatibale with Intel_2018 ADN bull OPENMPI, NETCDF loaded locallyF"
  . $rootdir/bldsva/machines/$platform/loadenv.Intel >> $log_file 2>> $err_file
  check

 # Default library paths
  defaultMpiPath="$OMPI_HOME"
  defaultNcdfPath="/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-intel14" # netcdf fortran
  defaultGribPath="/pf/b/b380190/lib/grib_api_dir13"
  defaultGribapiPath="/pf/b/b380190/lib/grib_api_dir13"
  defaultTclPath="/usr" # need to be installed in MISTRAL
  efaultJasperPath="/usr" # need to be installed in MISTRAL
  defaultHyprePath="/pf/b/b380190/opt/hypre-2.11.2/src/hypre" #??
  defaultSiloPath="/pf/b/b380190/opt/silo-4.10.2" # this is compiled with gcc (gcc/6.4.0, openmpi/2.0.2p2_hpcx-gcc64 ) 
  defaultPncdfPath="/sw/rhel6-x64/netcdf/netcdf_fortran-4.4.3-parallel-openmpi2-intel14" #netcdf parallel fortran
  defaultLapackPath="/pf/b/b380190/lib/lapack-3.3.0"

  export LD_LIBRARY_PATH=$defaultTclPath/lib64/:/home/f.gasper/local/libs/netcdf/lib:$defaultGribapiPath/lib # need to be corrected for first two
# better to export grib deffinition files!
 export GRIB_DEFINITION_PATH=/pf/b/b380190/lib/grib_api_dir13/share/grib_api/definitions/
 export GRIB_SAMPLES_PATH=/pf/b/b380190/lib/grib_api_dir13/share/grib_api/samples/

  # Default Compiler/Linker optimization
  defaultOptC="-O2"

  profilingImpl=" no scalasca "
  if [[ $profiling == "scalasca" ]] ; then ; profComp="" ; profRun="scalasca -analyse" ; profVar=""  ;fi
  # Default Processor settings
  defaultwtime="001:00:00"
  defaultQ="bb0722"

route "${cblue}<< getMachineDefaults${cnormal}"

# sourcing and loading module can be added later on here


route "${cblue}<< getMachineDefaults${cnormal}"
}

finalizeMachine(){
route "${cblue}>> finalizeMachine${cnormal}"
route "${cblue}<< finalizeMachine${cnormal}"
}


createRunscript(){
route "${cblue}>> createRunscript${cnormal}"
comment "not that the module are loaded manully"
#comment "   copy MISTRAL module load script into rundirectory"
#  cp $rootdir/bldsva/machines/$platform/loadenvs.$compiler $rundir/loadenvs
check

mpitasks=$((numInst * ($nproc_icon + $nproc_cos + $nproc_clm + $nproc_pfl + $nproc_oas)))
nnodes=`echo "scale = 2; $mpitasks / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`

#DA
if [[ $withPDAF == "true" ]] ; then
  srun="srun -n $mpitasks ./tsmp-pdaf -n_modeltasks $(($numInst-$startInst)) -filtertype 2 -subtype 1 -delt_obs $delta_obs -rms_obs 0.03 -obs_filename swc_crp"
else
  srun="srun --multi-prog slm_multiprog_mapping.conf"
fi

cat << EOF >> $rundir/tsmp_slm_run.bsh
#!/bin/bash

#SBATCH --job-name="TerrSysMP"
#SBATCH --nodes=$nnodes
#SBATCH --ntasks=$mpitasks
##SBATCH --ntasks-per-node=$nppn
#SBATCH --output=mpiMPMD-out.%j
#SBATCH --error=mpiMPMD-err.%j
#SBATCH --time=$wtime
#SBATCH --partition=compute
#SBATCH --mail-type=NONE
#SBATCH --account=bb0722 

export PSP_RENDEZVOUS_OPENIB=-1

cd $rundir
source $rundir/loadenvs
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

start_icon=$(($nproc_oas+$counter))
end_icon=$(($start_icon+($numInst*$nproc_icon)-1))

start_pfl=$(($numInst*($nproc_cos+$nproc_icon)+$nproc_oas+$counter))
end_pfl=$(($start_pfl+($numInst*$nproc_pfl)-1))

start_clm=$((($numInst*($nproc_cos+$nproc_icon))+$nproc_oas+($numInst*$nproc_pfl)+$counter))
end_clm=$(($start_clm+($numInst*$nproc_clm)-1))

if [[ $numInst > 1 &&  $withOASMCT == "true" ]] then
 for instance in {$startInst..$(($startInst+$numInst-1))}
 do
  for iter in {1..$nproc_cos}
  do
    if [[ $withCOS == "true" ]] then ; echo $instance >>  $rundir/instanceMap.txt ;fi
    if [[ $withICON == "true" ]] then ; echo $instance >>  $rundir/instanceMap.txt ;fi
  done
 done
fi

if [[ $numInst > 1 &&  $withOASMCT == "true" ]] then
 for instance in {$startInst..$(($startInst+$numInst-1))}
 do
  for iter in {1..$nproc_icon}
  do
    if [[ $withICON == "true" ]] then ; echo $instance >>  $rundir/instanceMap.txt ;fi
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


if [[ $withCOS == "true" ]]; then
cat << EOF >> $rundir/slm_multiprog_mapping.conf
__oas__
__cos__
__pfl__
__clm__

EOF
else
cat << EOF >> $rundir/slm_multiprog_mapping.conf
__oas__
__icon__
__pfl__
__clm__

EOF
fi

comment "   sed executables and processors into mapping file"
	if [[ $withOAS == "false" ||  $withOASMCT == "true" ]] then ; sed "s/__oas__//" -i $rundir/slm_multiprog_mapping.conf  >> $log_file 2>> $err_file; check; fi
	if [[ $withCOS == "false" ]] then ; sed "s/__cos__//" -i $rundir/slm_multiprog_mapping.conf  >> $log_file 2>> $err_file; check; fi
	if [[ $withICON == "false" ]] then ; sed "s/__icon__//" -i $rundir/slm_multiprog_mapping.conf  >> $log_file 2>> $err_file; check; fi
        if [[ $withPFL == "false" ]] then ; sed "s/__pfl__//" -i $rundir/slm_multiprog_mapping.conf  >> $log_file 2>> $err_file; check; fi
        if [[ $withCLM == "false" ]] then ; sed "s/__clm__//" -i $rundir/slm_multiprog_mapping.conf  >> $log_file 2>> $err_file; check; fi


sed "s/__oas__/$start_oas  $profRun \.\/oasis3.MPI1.x/" -i $rundir/slm_multiprog_mapping.conf >> $log_file 2>> $err_file
check
sed "s/__cos__/$start_cos-$end_cos  $profRun \.\/lmparbin_pur/" -i $rundir/slm_multiprog_mapping.conf >> $log_file 2>> $err_file
check
sed "s/__icon__/$start_icon-$end_icon  $profRun \.\/icon/" -i $rundir/slm_multiprog_mapping.conf >> $log_file 2>> $err_file
check
sed "s/__pfl__/$start_pfl-$end_pfl  $profRun \.\/parflow $pflrunname/" -i $rundir/slm_multiprog_mapping.conf >> $log_file 2>> $err_file
check
sed "s/__clm__/$start_clm-$end_clm  $profRun \.\/clm/" -i $rundir/slm_multiprog_mapping.conf >> $log_file 2>> $err_file
check



comment "   change permission of runscript and mapfile"
chmod 755 $rundir/tsmp_slm_run.bsh >> $log_file 2>> $err_file
check
chmod 755 $rundir/slm_multiprog_mapping.conf >> $log_file 2>> $err_file
check
route "${cblue}<< createRunscript${cnormal}"
}

