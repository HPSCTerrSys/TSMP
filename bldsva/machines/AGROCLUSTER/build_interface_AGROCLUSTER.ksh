#! /bin/ksh


getMachineDefaults(){
route "${cblue}>> getMachineDefaults${cnormal}"

  comment "   init lmod functionality"
    . /usr/share/Modules/init/ksh  >> $log_file 2>> $err_file
  check
  comment "   source and load Modules on $platform"
    module load gcc-ibg3  >> $log_file 2>> $err_file
  check
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

createRunscript(){
route "${cblue}>> createRunscript${cnormal}"

mpitasks=$((numInst*($nproc_cos + $nproc_clm + $nproc_pfl + $nproc_oas)))
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

if [[ $numInst > 1 && $withOASMCT == "true"   ]] ; then

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

