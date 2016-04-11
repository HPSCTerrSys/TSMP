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
  defaultOptC="-O2"


  # Default Processor settings
  defaultwtime="00:30:00"
  defaultQ="batch"

route "${cblue}<< getMachineDefaults${cnormal}"
}

createRunscript(){
route "${cblue}>> createRunscript${cnormal}"

mpitasks=$(($numInst*($nproc_cos + $nproc_clm + $nproc_pfl + $nproc_oas)))
nnodes=`echo "scale = 2; $mpitasks / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`



exel="mpirun "

for instance in {0..$(($numInst-1))}
do

if [[ $withCOS == "true" && $withOAS == "false" ]] ; then ; exel=$exel" -np $nproc_cos  ./cos_starter.ksh $instance :" ; fi 
if [[ $withPFL == "true" && $withOAS == "false" ]] ; then ; exel=$exel" -np $nproc_pfl  ./pfl_starter.ksh $instance :" ; fi
if [[ $withCLM == "true" && $withOAS == "false" ]] ; then ; exel=$exel" -np $nproc_clm  ./clm_starter.ksh $instance :" ; fi

if [[ $withCLM == "true" && $withCOS == "true"  && $withPFL == "false" && $withOASMCT == "false" ]] ; then ; exel=$exel" -np $nproc_oas ./oasis3.MPI1.x : -np $nproc_cos ./cos_starter.ksh $instance : -np $nproc_clm  ./clm_starter.ksh $instance :" ; fi
if [[ $withCLM == "true" && $withCOS == "false" && $withPFL == "true"  && $withOASMCT == "false" ]] ; then ; exel=$exel" -np $nproc_oas ./oasis3.MPI1.x : -np $nproc_pfl ./pfl_starter.ksh $instance : -np $nproc_clm  ./clm_starter.ksh $instance :" ; fi
if [[ $withCLM == "true" && $withPFL == "true" && $withCOS == "true"  && $withOASMCT == "false" ]] ; then ; exel=$exel" -np $nproc_oas ./oasis3.MPI1.x : -np $nproc_cos ./cos_starter.ksh $instance : -np $nproc_pfl ./pfl_starter.ksh $instance : -np $nproc_clm  ./clm_starter.ksh $instance :" ; fi


if [[ $withCLM == "true" && $withCOS == "true" && $withPFL == "false" && $withOASMCT == "true" ]] ; then ; exel=$exel" -np $nproc_cos ./cos_starter.ksh $instance : -np $nproc_clm  ./clm_starter.ksh $instance :" ; fi
if [[ $withCLM == "true" && $withCOS == "false" && $withPFL == "true"  && $withOASMCT == "true" ]] ; then ; exel=$exel" -np $nproc_pfl ./pfl_starter.ksh $instance : -np $nproc_clm  ./clm_starter.ksh $instance :" ; fi
if [[ $withCLM == "true" && $withPFL == "true" && $withCOS == "true"  && $withOASMCT == "true" ]] ; then ; exel=$exel" -np $nproc_cos ./cos_starter.ksh $instance : -np $nproc_pfl ./pfl_starter.ksh $instance : -np $nproc_clm  ./clm_starter.ksh $instance :" ; fi

done

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

date
$exel >> log_file 2>> err_file
date

echo "ready" > ready.txt
exit 0

EOF

if [[ $numInst > 1 && ( $withOASMCT == "true" || $withOAS == "false"   ) ]] ; then

cat << EOF >> $rundir/cos_starter.ksh
#!/bin/ksh
cd tsmp_instance_\$1
./lmparbin_pur
EOF

cat << EOF >> $rundir/clm_starter.ksh
#!/bin/ksh
cd tsmp_instance_\$1
./clm
EOF

cat << EOF >> $rundir/pfl_starter.ksh
#!/bin/ksh
cd tsmp_instance_\$1
./parflow $pflrunname
EOF

else

cat << EOF >> $rundir/cos_starter.ksh
#!/bin/ksh
./lmparbin_pur
EOF

cat << EOF >> $rundir/clm_starter.ksh
#!/bin/ksh
./clm
EOF

cat << EOF >> $rundir/pfl_starter.ksh
#!/bin/ksh
./parflow $pflrunname
EOF

fi


comment "   change permission of module starter scripts"
chmod 755 $rundir/cos_starter.ksh >> $log_file 2>> $err_file
check
chmod 755 $rundir/clm_starter.ksh >> $log_file 2>> $err_file
check
chmod 755 $rundir/pfl_starter.ksh >> $log_file 2>> $err_file
check



comment "   change permission of runscript"
chmod 755 $rundir/tsmp_pbs_run.ksh >> $log_file 2>> $err_file
check
route "${cblue}<< createRunscript${cnormal}"
}

