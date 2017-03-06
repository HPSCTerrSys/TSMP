#! /bin/ksh


getMachineDefaults(){
route "${cblue}>> getMachineDefaults${cnormal}"

  # Default library paths
  defaultMpiPath="/usr"
  defaultNcdfPath="/home/f.gasper/local/libs/netcdf"
  defaultGrib1Path="/home/w.kurtz/libs/DWD-libgrib1_061107"
  defaultTclPath="/usr"
  defaultHyprePath="/home/w.kurtz/libs/terrsysmp"
  defaultSiloPath="/home/w.kurtz/libs/terrsysmp"
  defaultPncdfPath="/usr"
  defaultLapackPath="/usr"
	
  export LD_LIBRARY_PATH=$defaultTclPath/lib64/:/home/f.gasper/local/libs/netcdf/lib


  # Default Compiler/Linker optimization
  defaultOptC="-O2"

  profilingImpl=" no "

  # Default Processor settings
  defaultwtime="00:30:00"
  defaultQ="batch"

route "${cblue}<< getMachineDefaults${cnormal}"
}



finalizeMachine(){
route "${cblue}>> finalizeMachine${cnormal}"
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
#PBS -S /bin/ksh
#PBS -N TerrSysMP_run
#PBS -l walltime=$wtime
#PBS -l nodes=$nnodes:ppn=$nppn
#PBS -V 
#PBS -u $USER
#PBS -q $queue
#PBS -l pvmem=2800mb

cd $rundir

rm -rf  YU*
export LD_LIBRARY_PATH=$defaultTclPath/lib64/:$defaultNcdfPath/lib
echo "started" > started.txt

date
$exel 
date

echo "ready" > ready.txt

exit 0

EOF


comment "   change permission of runscript"
chmod 755 $rundir/tsmp_pbs_run.ksh >> $log_file 2>> $err_file
check
route "${cblue}<< createRunscript${cnormal}"
}

