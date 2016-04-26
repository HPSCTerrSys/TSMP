#! /bin/ksh


getMachineDefaults(){
route "${cblue}>> getMachineDefaults${cnormal}"

  # Default library paths
  defaultMpiPath="/opt/ibg3/openmpi-1.8.4"
  defaultNcdfPath="/home/w.kurtz/libs/netcdf_openmpi-1.8.4"
  defaultGrib1Path="/home/w.kurtz/libs/DWD-libgrib1_061107"
  defaultTclPath="/usr"
  defaultHyprePath="/home/w.kurtz/libs/hypre_openmpi-1.8.4"
  defaultSiloPath="/home/w.kurtz/libs/silo-4.8_openmpi-1.8.4"
	
  export LD_LIBRARY_PATH=$defaultTclPath/lib64/:/home/w.kurtz/libs/netcdf_openmpi-1.8.4/lib


  # Default Compiler/Linker optimization
  defaultOptC="-O2"


  # Default Processor settings
  defaultwtime="00:30:00"
  defaultQ="batch"

route "${cblue}<< getMachineDefaults${cnormal}"
}



finalizeMachine(){
route "${cblue}>> finalizeMachine${cnormal}"

  comment "   init lmod functionality"
    . /usr/share/Modules/init/ksh  >> $log_file 2>> $err_file
  check
  comment "   source and load Modules on $platform"
    module load gcc-ibg3  >> $log_file 2>> $err_file
  check

route "${cblue}<< finalizeMachine${cnormal}"
}



createRunscript(){
route "${cblue}>> createRunscript${cnormal}"

mpitasks=$((numInst*($nproc_cos + $nproc_clm + $nproc_pfl + $nproc_oas)))
nnodes=`echo "scale = 2; $mpitasks / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`


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


exel="mpirun "


if [[ $withOAS == "true" && $withOASMCT == "false" ]] ; then ; exel=$exel" -np $nproc_oas  ./oasis3.MPI1.x :" ; fi
if [[ $withCOS == "true" ]] ; then ; exel=$exel" -np $(($numInst*$nproc_cos))  ./lmparbin_pur :" ; fi  
if [[ $withPFL == "true" ]] ; then ; exel=$exel" -np $(($numInst*$nproc_pfl))  ./parflow $pflrunname  :" ; fi
if [[ $withCLM == "true" ]] ; then ; exel=$exel" -np $(($numInst*$nproc_clm))  ./clm  :" ; fi

exel=${exel%?} #remove trailing ":"


cat << EOF >> $rundir/tsmp_pbs_run.ksh

#Job Submission to Agrocluster
#PBS -S /bin/csh
#PBS -N TerrSysMP_run
#PBS -l walltime=$wtime
#PBS -l nodes=$nnodes:ppn=$nppn
#PBS -V 
#PBS -u $USER
#PBS -q $queue


cd $rundir

rm -rf  YU*


date
$exel 
date

exit 0

EOF


comment "   change permission of runscript"
chmod 755 $rundir/tsmp_pbs_run.ksh >> $log_file 2>> $err_file
check
route "${cblue}<< createRunscript${cnormal}"
}

