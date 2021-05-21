#! /bin/ksh


getMachineDefaults(){
route "${cyellow}>> getMachineDefaults${cnormal}"

  # Default library paths
  comment "   init lmod functionality"
    eval `tclsh /usr/local/module/modulecmd.tcl sh autoinit`
  check
  comment "   load modules"
    module purge >> $log_file 2>> $err_file
  check
    module load lapack >> $log_file 2>> $err_file
  check
    module load parallel-netcdf >> $log_file 2>> $err_file 
  check

  defaultMpiPath="/homea/slts/slts00/local/juqueen/mpi"
  defaultNcdfPath="/homea/slts/slts00/local/juqueen/netcdf"
  defaultGrib1Path="/homea/slts/slts00/local/juqueen/grib1"
  defaultTclPath="/homea/slts/slts00/local/juqueen/tcl"
  defaultHyprePath="/homea/slts/slts00/local/juqueen/hypre"
  defaultSiloPath="/homea/slts/slts00/local/juqueen/silo"
  defaultLapackPath="$LAPACK_DIR"
  defaultPncdfPath="$PNETCDF_DIR"


  # Default Compiler/Linker optimization
  defaultOptC="-O2 -qhot -qarch=qp -qtune=qp"

  profilingImpl=" no scalasca "
  if [[ $profiling == "scalasca" ]] ; then ;
    comment "   load scalasca module" 
      module load UNITE scalasca >> $log_file 2>> $err_file
    check
    profComp="scalasca -instrument" ; profRun="scalasca -analyze" ; profVar="ESD_BUFFER_SIZE=300000"  
  fi

  # Default Processor settings
  defaultwtime="00:30:00"
  defaultQ="queue"

route "${cyellow}<< getMachineDefaults${cnormal}"
}


finalizeMachine(){
route "${cyellow}>> finalizeMachine${cnormal}"
route "${cyellow}<< finalizeMachine${cnormal}"
}


createRunscript(){
route "${cyellow}>> createRunscript${cnormal}"

mpitasks=$(($numInst*($nproc_cos + $nproc_clm + $nproc_pfl + $nproc_oas)))
nnodes=`echo "scale = 2; $mpitasks / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`
bgs_oas=`echo "scale = 2; $nproc_oas / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'` 
bgs_cos=`echo "scale = 2; $numInst * $nproc_cos / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'` 
bgs_clm=`echo "scale = 2; $numInst * $nproc_clm / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'` 
bgs_pfl=`echo "scale = 2; $numInst * $nproc_pfl / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'` 



errorname='$(job_name).$(jobid).err'
outputname='$(job_name).$(jobid).out'
runflags="--envs LOGNAME=$rundir --envs $profVar "

cat << EOF >> $rundir/tsmp_ll_run.ksh
#!/bin/ksh

#@ job_name= TerrSysMP_run
#@ error=$rundir/$errorname
#@ output=$rundir/$outputname
#@ environment = COPY_ALL
#@ notification = never
#@ wall_clock_limit = $wtime
#@ bg_size = $nnodes
#@ job_type = bluegene
#@ bg_connectivity = TORUS
#@ $queue

cd $rundir

eval \`tclsh /usr/local/module/modulecmd.tcl sh autoinit\`
module load UNITE scalasca

sed -i -r '/^[ ]/d' ll_multiprog_mapping.conf
rm -rf partinfos.txt base YU*

runjob -p $nppn --verbose 4 : /bgsys/local/samples/personality/personality.elf > partinfos.txt
grep -i "task" partinfos.txt | sort | cut -d " " -f 17 | sed 's/,/ /g' | sed 's/(/ /;s/)//' > base

EOF





counter=0
counter2=1

#for mapfile
start_oas=$counter
end_oas=$(($start_oas+$nproc_oas-1))

start_cos=$(($nproc_oas+$counter))
end_cos=$(($start_cos+($numInst*$nproc_cos)-1))

start_pfl=$(($numInst*$nproc_cos+$nproc_oas+$counter))
end_pfl=$(($start_pfl+($numInst*$nproc_pfl)-1))

start_clm=$((($numInst*$nproc_cos)+$nproc_oas+($numInst*$nproc_pfl)+$counter))
end_clm=$(($start_clm+($numInst*$nproc_clm)-1))


#for cutting from personality
start_oasM=$counter2
end_oasM=$(($start_oasM + $nproc_oas - 1))

start_cosM=$(($bgs_oas * $nppn + $counter2))
end_cosM=$(($start_cosM + ($numInst*$nproc_cos) - 1))

start_pflM=$((($bgs_cos  + $bgs_oas) * $nppn + $counter2))
end_pflM=$(($start_pflM + ($numInst*$nproc_pfl) - 1))

start_clmM=$((($bgs_cos + $bgs_oas + $bgs_pfl ) * $nppn + $counter2))
end_clmM=$(($start_clmM + ($numInst*$nproc_clm) - 1))

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



cat << EOF >> $rundir/tsmp_ll_run.ksh
__oas__
__cos__
__pfl__
__clm__

EOF


comment "   sed personality creation into run script for instance $instance"
  if [[ $withOAS == "false" ||  $withOASMCT == "true" ]] then ; sed "s/__oas__//" -i $rundir/tsmp_ll_run.ksh  >> $log_file 2>> $err_file; check; fi
  if [[ $withCOS == "false" ]] then ; sed "s/__cos__//" -i $rundir/tsmp_ll_run.ksh  >> $log_file 2>> $err_file; check; fi
  if [[ $withCLM == "false" ]] then ; sed "s/__clm__//" -i $rundir/tsmp_ll_run.ksh  >> $log_file 2>> $err_file; check; fi
  if [[ $withPFL == "false" ]] then ; sed "s/__pfl__//" -i $rundir/tsmp_ll_run.ksh  >> $log_file 2>> $err_file; check; fi

sed "s@__oas__@sed -n \"${start_oasM},${end_oasM}p\" base >> $rundir/ll_multiprog_mapping.conf@" -i $rundir/tsmp_ll_run.ksh >> $log_file 2>> $err_file
check
sed "s@__cos__@sed -n \"${start_cosM},${end_cosM}p\" base >> $rundir/ll_multiprog_mapping.conf@" -i $rundir/tsmp_ll_run.ksh >> $log_file 2>> $err_file
check
sed "s@__pfl__@sed -n \"${start_pflM},${end_pflM}p\" base >> $rundir/ll_multiprog_mapping.conf@" -i $rundir/tsmp_ll_run.ksh >> $log_file 2>> $err_file
check
sed "s@__clm__@sed -n \"${start_clmM},${end_clmM}p\" base >> $rundir/ll_multiprog_mapping.conf@" -i $rundir/tsmp_ll_run.ksh >> $log_file 2>> $err_file
check




cat << EOF >> $rundir/ll_multiprog_mapping.conf
__oas1__
__oas2__
__oas3__
__cos1__
__cos2__
__cos3__
__pfl1__
__pfl2__
__pfl3__
__clm1__
__clm2__
__clm3__

EOF

if [[ $withCLM == "true" ]] then ; exe=clm ; fi
if [[ $withPFL == "true" ]] then ; exe=parflow ; fi
if [[ $withCOS == "true" ]] then ; exe=lmparbin_pur ; fi

comment "   sed executables and processors into mapping file for instance $instance"
  if [[ $withOAS == "false" ||  $withOASMCT == "true" ]] then ; sed "s/__oas.__//" -i $rundir/ll_multiprog_mapping.conf  >> $log_file 2>> $err_file; check; fi
  if [[ $withCOS == "false" ]] then ; sed "s/__cos.__//" -i $rundir/ll_multiprog_mapping.conf  >> $log_file 2>> $err_file; check; fi
  if [[ $withCLM == "false" ]] then ; sed "s/__clm.__//" -i $rundir/ll_multiprog_mapping.conf  >> $log_file 2>> $err_file; check; fi
  if [[ $withPFL == "false" ]] then ; sed "s/__pfl.__//" -i $rundir/ll_multiprog_mapping.conf  >> $log_file 2>> $err_file; check; fi

  sed "s/__cos1__/#mpmdbegin $start_cos-$end_cos/" -i $rundir/ll_multiprog_mapping.conf >> $log_file 2>> $err_file
  check
  sed "s/__cos2__/#mpmdcmd lmparbin_pur/" -i $rundir/ll_multiprog_mapping.conf >> $log_file 2>> $err_file
  check
  sed "s/__cos3__/#mpmdend/" -i $rundir/ll_multiprog_mapping.conf >> $log_file 2>> $err_file
  check

  sed "s/__clm1__/#mpmdbegin $start_clm-$end_clm/" -i $rundir/ll_multiprog_mapping.conf >> $log_file 2>> $err_file
  check
  sed "s/__clm2__/#mpmdcmd clm/" -i $rundir/ll_multiprog_mapping.conf >> $log_file 2>> $err_file
  check
  sed "s/__clm3__/#mpmdend/" -i $rundir/ll_multiprog_mapping.conf >> $log_file 2>> $err_file
  check

  sed "s/__pfl1__/#mpmdbegin $start_pfl-$end_pfl/" -i $rundir/ll_multiprog_mapping.conf >> $log_file 2>> $err_file
  check
  sed "s/__pfl2__/#mpmdcmd parflow $pflrunname/" -i $rundir/ll_multiprog_mapping.conf >> $log_file 2>> $err_file
  check
  sed "s/__pfl3__/#mpmdend/" -i $rundir/ll_multiprog_mapping.conf >> $log_file 2>> $err_file
  check

  sed "s/__oas1__/#mpmdbegin $start_oas/" -i $rundir/ll_multiprog_mapping.conf >> $log_file 2>> $err_file
  check
  sed "s/__oas2__/#mpmdcmd oasis3.MPI1.x/" -i $rundir/ll_multiprog_mapping.conf >> $log_file 2>> $err_file
  check
  sed "s/__oas3__/#mpmdend/" -i $rundir/ll_multiprog_mapping.conf >> $log_file 2>> $err_file
  check

#DA
if [[ $withPDAF == "true" ]] ; then
  runjob="runjob -p $nppn -n $mpitasks $runflags : ./tsmp-pdaf -n_modeltasks $(($numInst-$startInst)) -filtertype 2 -delt_obs $delta_obs -rms_obs 0 -obs_filename noname"
else
  runjob="runjob -p $nppn -n $mpitasks --mapping ll_multiprog_mapping.conf $runflags : ./$exe"
fi



cat << EOF >> $rundir/tsmp_ll_run.ksh
echo "started" > started.txt
date
$runjob
date

echo "ready" > ready.txt
exit 0

EOF



comment "   change permission of runscript and mapfile"
chmod 755 $rundir/tsmp_ll_run.ksh >> $log_file 2>> $err_file
check
chmod 755 $rundir/ll_multiprog_mapping.conf >> $log_file 2>> $err_file
check
route "${cyellow}<< createRunscript${cnormal}"
}

