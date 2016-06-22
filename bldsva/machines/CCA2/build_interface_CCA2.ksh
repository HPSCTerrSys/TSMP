#! /bin/ksh


getMachineDefaults(){
route "${cblue}>> getMachineDefaults${cnormal}"

  comment "   init lmod functionality"
    . /opt/modules/default/init/ksh  >> $log_file 2>> $err_file
  check
  comment "   source and load Modules on $platform"
    . $rootdir/bldsva/machines/$platform/loadenvs  >> $log_file 2>> $err_file
  check

  # Default library paths
  defaultMpiPath="$CRAY_MPICH_DIR"
  defaultNcdfPath="$NETCDF_DIR"
  defaultGrib1Path="/home/ms/de/dwd/libdwd/1.1.6/GNU/lib"
  defaultTclPath="/perm/ms/spde/de5f/libs/tcl"
  defaultHyprePath="/perm/ms/spde/de5f/libs/hypre"
  defaultSiloPath="/perm/ms/spde/de5f/libs/silo"
  defaultPncdfPath=""
  defaultLapackPath=""

	
  export LD_LIBRARY_PATH=$defaultNcdfPath/lib


  # Default Compiler/Linker optimization
  defaultOptC="-O2"


  # Default Processor settings
  defaultwtime="00:30:00"
  defaultQ="np"

route "${cblue}<< getMachineDefaults${cnormal}"
}

finalizeMachine(){
route "${cblue}>> finalizeMachine${cnormal}"
route "${cblue}<< finalizeMachine${cnormal}"
}



createRunscript(){
route "${cblue}>> createRunscript${cnormal}"

comment "   copy CCA module load script into rundirectory"
  cp $rootdir/bldsva/machines/$platform/loadenvs $rundir
check


mpitasks=`expr $nproc_cos + $nproc_clm + $nproc_pfl + $nproc_oas`
nnodes1=`echo "scale = 2; $numInst * $nproc_cos / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`
nnodes2=`echo "scale = 2; $numInst * $nproc_clm / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`
nnodes3=`echo "scale = 2; $numInst * $nproc_pfl / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`
nnodes4=`echo "scale = 2; $nproc_oas / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`
nnodes=$(($nnodes1+$nnodes2+$nnodes3+$nnodes4))


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


exel="aprun "


if [[ $withOAS == "true" && $withOASMCT == "false" ]] ; then ; exel=$exel" -n $nproc_oas  ./oasis3.MPI1.x :" ; fi
if [[ $withCOS == "true" ]] ; then ; exel=$exel" -n $(($numInst*$nproc_cos))  ./lmparbin_pur :" ; fi  
if [[ $withPFL == "true" ]] ; then ; exel=$exel" -n $(($numInst*$nproc_pfl))  ./parflow $pflrunname  :" ; fi
if [[ $withCLM == "true" ]] ; then ; exel=$exel" -n $(($numInst*$nproc_clm))  ./clm  :" ; fi

exel=${exel%?} #remove trailing ":"


cat << EOF >> $rundir/tsmp_pbs_run.ksh
#!/bin/ksh

#Job Submission to cca
#PBS -N TerrSysMP_run
#PBS -l walltime=$wtime
#PBS -l EC_total_tasks=$mpitasks
#PBS -l EC_tasks_per_node=$nppn
#PBS -l EC_threads_per_task=1
#PBS -l EC_nodes=$nnodes
#PBS -V 
#PBS -q $queue

cd $rundir
. $rundir/loadenvs

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

