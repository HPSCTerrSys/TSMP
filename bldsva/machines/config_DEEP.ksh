#! /bin/ksh


getMachineDefaults(){
route "${cyellow}>> getMachineDefaults${cnormal}"
  comment "   init lmod functionality"
  . /beegfs/usr_local/software/lmod/lmod/init/ksh >> $log_file 2>> $err_file
  check
  comment "   source and load Modules on DEEP: loadenvs.$compiler"
  . $rootdir/bldsva/machines/loadenvs.$compiler >> $log_file 2>> $err_file
  check

  defaultMpiPath="$EBROOTPSMPI"
  defaultNcdfPath="$EBROOTNETCDFMINFORTRAN"
  defaultGrib1Path=""
  defaultGribPath="$EBROOTECCODES"
  defaultGribapiPath="$EBROOTECCODES"
  defaultJasperPath="$EBROOTJASPER"
  defaultTclPath="$EBROOTTCL"
  defaultHyprePath="$EBROOTHYPRE"
  defaultSiloPath="$EBROOTSILO"
  defaultLapackPath="$EBROOTIMKL"
  defaultPncdfPath="$EBROOTPARALLELMINNETCDF"
  # Additional option for GPU compilation  
  gpuMpiSettings=
  cuda_architectures="-DCMAKE_CUDA_ARCHITECTURES=70"
#
  # Default Compiler/Linker optimization
  if [[ $compiler == "Gnu" ]] ; then
      defaultOptC="-O2" # Gnu
  elif [[ $compiler == "Intel" ]] ; then
      defaultOptC="-O2 -xHost" # Intel
  else
      defaultOptC="-O2" # Default
  fi

  profilingImpl=" no scalasca "  
  if [[ $profiling == "scalasca" ]] ; then ; profComp="" ; profRun="scalasca -analyse" ; profVar=""  ;fi

  # Default Processor settings
  defaultwtime="01:00:00"
  defaultQ="dp-cn"

  # DEEP only (because of problems with the perl path in CLM 3.5)
  if [[ $platform == "DEEP" ]]; then
      echo
      echo "Platform: $platform"
      echo
      
      sed -i -e 's+/usr/bin/env perl+/usr/bin/perl+' ${rootdir}/clm3_5/bld/mkDepends
      sed -i -e 's+/usr/bin/env perl+/usr/bin/perl+' ${rootdir}/clm3_5/bld/configure
      sed -i -e 's+/usr/bin/env perl+/usr/bin/perl+' ${rootdir}/clm3_5/bld/mkSrcfiles
      sed -i -e 's+/usr/bin/env perl+/usr/bin/perl+' ${rootdir}/clm3_5/bld/queryDefaultNamelist.pl
      sed -i -e 's+/usr/bin/env perl+/usr/bin/perl+' ${rootdir}/bldsva/intf_oas3/clm3_5/arch/config/configure
  fi
  
  route "${cyellow}<< getMachineDefaults${cnormal}"
}

finalizeMachine(){
route "${cyellow}>> finalizeMachine${cnormal}"
route "${cyellow}<< finalizeMachine${cnormal}"
}

# computes nodes based on number of processors and resources
computeNodes(){
processes=$1
resources=$2
echo $((processes%resources?processes/resources+1:processes/resources))
}

createRunscript(){
route "${cyellow}>> createRunscript${cnormal}"
comment "   copy DEEP module load script into rundirectory"
  cp $rootdir/bldsva/machines/loadenvs.$compiler $rundir/loadenvs
check

mpitasks=$((numInst * ($nproc_icon + $nproc_cos + $nproc_clm + $nproc_pfl + $nproc_oas)))
nnodes=`echo "scale = 2; $mpitasks / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`

#DA
if [[ $withPDAF == "true" ]] ; then
  srun="srun -n $mpitasks ./tsmp-pdaf -n_modeltasks $(($numInst-$startInst)) -filtertype 2 -subtype 1 -delt_obs $delta_obs -rms_obs 0.03 -obs_filename swc_crp"
else
  srun="srun --multi-prog slm_multiprog_mapping.conf"
fi
## Heterogeneous and modular jobs

if [[ $processor == "GPU" || $processor == "MSA" ]]; then
nnode_cos=$((($nproc_cos)/$nppn))
nnode_clm=$((($nproc_clm)/$nppn))
nnode_pfl=$((($nproc_pfl)/$ngpn))

nnode_cos=$(computeNodes $nproc_cos $nppn)
nnode_clm=$(computeNodes $nproc_clm $nppn)
nnode_pfl=$(computeNodes $nproc_pfl $ngpn)

route "${cyellow}<< setting up heterogeneous/modular job${cnormal}"
comment "  nppn=$nppn\t\t ngpn=$ngpn\n"
comment "  Nproc: COSMO=$nproc_cos\tCLM=$nproc_clm\tPFL=$nproc_pfl\n"
comment "  Nnode: COSMO=$nnode_cos\tCLM=$nnode_clm\tPFL=$nnode_pfl\n"

cat << EOF >> $rundir/tsmp_slm_run.bsh
#!/bin/bash
#SBATCH --account=deepsea
#SBATCH --job-name="TSMP_Hetero"
#SBATCH --output=hetro_job-out.%j
#SBATCH --error=hetro_job-err.%j
#SBATCH --time=$wtime
#SBATCH -N $nnode_cos --ntasks-per-node=$nppn -p dp-cn
#SBATCH hetjob
#SBATCH -N $nnode_clm --ntasks-per-node=$nppn  -p dp-cn
#SBATCH hetjob
#SBATCH -N $nnode_pfl --ntasks-per-node=$ngpn --gres=gpu:$ngpn -p dp-esb

cd $rundir
source $rundir/loadenvs
export LD_LIBRARY_PATH="$rootdir/${mList[3]}_${platform}_${combination}/rmm/lib:\$LD_LIBRARY_PATH"
date
echo "started" > started.txt
rm -rf YU*

module unload nvidia-driver/.default
export CUDA_VISIBLE_DEVICES=0

srun --het-group=0 ./lmparbin_pur :\\
     --het-group=1 ./clm :\\
     --het-group=2 ./parflow $pflrunname
date
echo "ready" > ready.txt
exit 0
EOF

else

cat << EOF >> $rundir/tsmp_slm_run.bsh
#!/bin/bash

#SBATCH --job-name="TSMP"
#SBATCH --nodes=$nnodes
#SBATCH --ntasks=$mpitasks
#SBATCH --ntasks-per-node=$nppn
#SBATCH --output=mpiMPMD-out.%j
#SBATCH --error=mpiMPMD-err.%j
#SBATCH --time=$wtime
#SBATCH --partition=$queue
#SBATCH --mail-type=NONE
#SBATCH --account=deepsea 

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
fi


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
route "${cyellow}<< createRunscript${cnormal}"
}

Npp=24
Ngp=1

PFLProcXg=1
PFLProcYg=1
CLMProcXg=3
CLMProcYg=8
COSProcXg=12
COSProcYg=16
if [[ $refSetup == "cordex" ]] then
	PFLProcX=9
	PFLProcY=8
	CLMProcX=3
	CLMProcY=8
	COSProcX=12
	COSProcY=16
	elif [[ $refSetup == "nrw" ]] then
	PFLProcX=4
	PFLProcY=6
	CLMProcX=8
	CLMProcY=3
	COSProcX=12
	COSProcY=16
	elif [[ $refSetup == "idealRTD" ]] then
	PFLProcX=2
	PFLProcY=2
	CLMProcX=2
	CLMProcY=6
	COSProcX=6
	COSProcY=5
	else 
	PFLProcX=4
	PFLProcY=4
	CLMProcX=4
	CLMProcY=2
	COSProcX=8
	COSProcY=8
fi

