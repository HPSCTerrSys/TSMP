#! /bin/ksh


getMachineDefaults(){
route "${cblue}>> getMachineDefaults${cnormal}"
  # Default library paths
  defaultMpiPath="/daten01/z4/openmpi-1.4.3_gnu"
  defaultNcdfPath="/daten01/z4/netcdf4.1.3_gfortran_4.2.1"
  defaultGrib1Path="/daten01/z4/DWD_libgrib/libgrib1_061107gnu"
  defaultTclPath="/daten01/z4/tcl8.5.13_gnu"
  defaultHyprePath="/daten01/z4/hypre2.9_gnu"
  defaultSiloPath="/daten01/z4/silo4.8_gnu"
  defaultPncdfPath="/daten01/z4/p-netcdf"
  defaultLapackPath="/daten01/z4/lapack-3.6.0"

  export LD_LIBRARY_PATH=$defaultTclPath/lib/


  # Default Compiler/Linker optimization
  defaultOptC="-O2"

  profilingImpl=" no "

  # Default Processor settings
  defaultwtime="00:300:00"
  defaultQ="batch"

route "${cblue}<< getMachineDefaults${cnormal}"
}

finalizeMachine(){
route "${cblue}>> finalizeMachine${cnormal}"
route "${cblue}<< finalizeMachine${cnormal}"
}


createRunscript(){
route "${cblue}>> createRunscript${cnormal}"

mpitasks=$(($numInst*($nproc_cos + $nproc_clm + $nproc_pfl + $nproc_oas)))
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
#!/bin/ksh

#Job Submission to Cluma
#PBS -N TerrSysMP_run
#PBS -l walltime=$wtime
#PBS -l nodes=$nnodes:ppn=$nppn
#PBS -V 
#PBS -u $USER
#PBS -q $queue

export LOGNAME="$rundir"

cd $rundir

rm -rf  YU*

export LD_LIBRARY_PATH=$defaultNcdfPath/lib/
echo "started" > started.txt
date
$exel >> log_file 2>> err_file
date

echo "ready" > ready.txt
exit 0

EOF


comment "   change permission of runscript"
chmod 755 $rundir/tsmp_pbs_run.ksh >> $log_file 2>> $err_file
check
route "${cblue}<< createRunscript${cnormal}"
}

