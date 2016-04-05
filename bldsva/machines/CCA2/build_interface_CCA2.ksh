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
  defaultMpiPath="/opt/cray/mpt/6.3.1/gni/mpich2-gnu/48"
  defaultNcdfPath="/usr/local/apps/netcdf/3.6.3/GNU/48"
  defaultGrib1Path="/home/ms/de/dwd/libdwd/1.1.6/GNU/lib"
  defaultTclPath="/perm/ms/spde/de5f/libs/tcl"
  defaultHyprePath="/perm/ms/spde/de5f/libs/hypre"
  defaultSiloPath="/perm/ms/spde/de5f/libs/silo"
	
  export LD_LIBRARY_PATH=$defaultNcdfPath/lib


  # Default Compiler/Linker optimization
  defaultOptC="-O2"


  # Default Processor settings
  defaultwtime="00:30:00"
  defaultQ="np"

route "${cblue}<< getMachineDefaults${cnormal}"
}

createRunscript(){
route "${cblue}>> createRunscript${cnormal}"

comment "   copy CCA module load script into rundirectory"
  cp $rootdir/bldsva/machines/$platform/loadenvs $rundir
check


mpitasks=`expr $nproc_cos + $nproc_clm + $nproc_pfl + $nproc_oas`
nnodes1=`echo "scale = 2; $nproc_cos / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`
nnodes2=`echo "scale = 2; $nproc_clm / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`
nnodes3=`echo "scale = 2; $nproc_pfl / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`
nnodes4=`echo "scale = 2; $nproc_oas / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`
nnodes=$(($nnodes1+$nnodes2+$nnodes3+$nnodes4))


exel=""
if [[ $withCOS == "true" && $withOAS == "false" ]] ; then ; exel="aprun -n $nproc_cos  -N $nppn  ./lmparbin_pur" ; fi 
if [[ $withPFL == "true" && $withOAS == "false" ]] ; then ; exel="aprun -n $nproc_pfl  -N $nppn ./parflow $pflrunname" ; fi
if [[ $withCLM == "true" && $withOAS == "false" ]] ; then ; exel="aprun -n $nproc_clm  -N $nppn ./clm" ; fi

if [[ $withCLM == "true" && $withCOS == "true"  && $withOASMCT == "false" ]] ; then ; exel="aprun -n $nproc_oas ./oasis3.MPI1.x : -n $nproc_cos ./lmparbin_pur : -n $nproc_clm  ./clm" ; fi
if [[ $withCLM == "true" && $withPFL == "true"  && $withOASMCT == "false" ]] ; then ; exel="aprun -n $nproc_oas ./oasis3.MPI1.x : -n $nproc_pfl  ./parflow $pflrunname : -n $nproc_clm  ./clm" ; fi
if [[ $withCLM == "true" && $withPFL == "true" && $withCOS == "true"  && $withOASMCT == "false" ]] ; then ; exel="aprun -n $nproc_oas  ./oasis3.MPI1.x : -n $nproc_cos ./lmparbin_pur : -n $nproc_pfl ./parflow $pflrunname : -n $nproc_clm  ./clm" ; fi


if [[ $withCLM == "true" && $withCOS == "true"  && $withOASMCT == "true" ]] ; then ; exel="aprun -n $nproc_cos ./lmparbin_pur : -n $nproc_clm  ./clm" ; fi
if [[ $withCLM == "true" && $withPFL == "true"  && $withOASMCT == "true" ]] ; then ; exel="aprun -n $nproc_pfl ./parflow $pflrunname : -n $nproc_clm  ./clm" ; fi
if [[ $withCLM == "true" && $withPFL == "true" && $withCOS == "true"  && $withOASMCT == "true" ]] ; then ; exel="aprun -n $nproc_cos ./lmparbin_pur : -n $nproc_pfl ./parflow $pflrunname : -n $nproc_clm ./clm" ; fi






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

