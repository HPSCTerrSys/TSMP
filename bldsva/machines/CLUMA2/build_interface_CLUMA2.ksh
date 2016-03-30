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

  export LD_LIBRARY_PATH=$defaultTclPath/lib/


  # Default Compiler/Linker optimization
  defaultOptC="-O0"


  # Default Processor settings
  defaultwtime="00:30:00"
  defaultQ="batch"

route "${cblue}<< getMachineDefaults${cnormal}"
}

createRunscript(){
route "${cblue}>> createRunscript${cnormal}"

mpitasks=`expr $nproc_cos + $nproc_clm + $nproc_pfl + $nproc_oas`
nnodes=`echo "scale = 2; $mpitasks / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`



exel=""
if [[ $withCOS == "true" && $withOAS == "false" ]] ; then ; exel="mpirun -np $nproc_cos  ./lmparbin_pur" ; fi 
if [[ $withPFL == "true" && $withOAS == "false" ]] ; then ; exel="mpirun -np $nproc_pfl  ./parflow $pflrunname" ; fi
if [[ $withCLM == "true" && $withOAS == "false" ]] ; then ; exel="mpirun -np $nproc_clm  ./clm" ; fi

if [[ $withCLM == "true" && $withCOS == "true"  && $withOASMCT == "false" ]] ; then ; exel="mpirun -np $nproc_oas ./oasis3.MPI1.x : -np $nproc_clm  ./clm : -np $nproc_cos ./lmparbin_pur" ; fi
if [[ $withCLM == "true" && $withPFL == "true"  && $withOASMCT == "false" ]] ; then ; exel="mpirun -np $nproc_oas ./oasis3.MPI1.x : -np $nproc_clm  ./clm : -np $nproc_pfl ./parflow $pflrunname" ; fi
if [[ $withCLM == "true" && $withPFL == "true" && $withCOS == "true"  && $withOASMCT == "false" ]] ; then ; exel="mpirun -np $nproc_oas ./oasis3.MPI1.x : -np $nproc_clm  ./clm : -np $nproc_cos ./lmparbin_pur : -np $nproc_pfl ./parflow $pflrunname" ; fi


if [[ $withCLM == "true" && $withCOS == "true"  && $withOASMCT == "true" ]] ; then ; exel="mpirun -np $nproc_clm  ./clm : -np $nproc_cos ./lmparbin_pur" ; fi
if [[ $withCLM == "true" && $withPFL == "true"  && $withOASMCT == "true" ]] ; then ; exel="mpirun -np $nproc_clm  ./clm : -np $nproc_pfl ./parflow $pflrunname" ; fi
if [[ $withCLM == "true" && $withPFL == "true" && $withCOS == "true"  && $withOASMCT == "true" ]] ; then ; exel="mpirun -np $nproc_clm  ./clm : -np $nproc_cos ./lmparbin_pur : -np $nproc_pfl ./parflow $pflrunname" ; fi



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

