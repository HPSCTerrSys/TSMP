#! /bin/ksh


getMachineDefaults(){
print "${cblue}>> getMachineDefaults${cnormal}"
  # Default library paths
  print -n "   init lmod functionality"
  . /usr/local/software/lmod/lmod/init/ksh >> $log_file 2>> $err_file
  check
  print -n "   source and load Modules on JURECA"
  . $rootdir/bldsva/machines/JURECA/loadenvs >> $log_file 2>> $err_file
  check
  defaultMpiPath="$EBROOTPSMPI"
  defaultNcdfPath="$EBROOTNETCDFMINFORTRAN"
  defaultGrib1Path="/homea/slts/slts00/sandbox/grib1-DWD20061107.jureca_tc2015.07_psintel_opt_KGo2/lib"
  defaultGribapiPath="$EBROOTGRIB_API"
  defaultJasperPath="$EBROOTJASPER"
  defaultTclPath="$EBROOTTCL"
  defaultHyprePath="$EBROOTHYPRE"
  defaultSiloPath="$EBROOTSILO"

  # Default Compiler/Linker optimization
  defaultOptC="-O0"

  # Default Processor settings
  defaultNppn=48
  defaultwtime="01:00:00"
  defaultQ="batch"
  defaultCLMProcX=4
  defaultCLMProcY=2
  defaultCOSProcX=8
  defaultCOSProcY=8
  defaultPFLProcX=4
  defaultPFLProcY=4


  # Default setups
  setups="nrw ideal300150 ideal600300 ideal1200600 ideal24001200"
  defaultRefSetup="nrw"
print "${cblue}<< getMachineDefaults${cnormal}"
}


createRunscript(){

mpitasks=`expr $nproc_cos + $nproc_clm + $nproc_pfl + $nproc_oas`
nnodes=`echo "scale = 2; $mpitasks / $nppn" | bc | perl -nl -MPOSIX -e 'print ceil($_);'`
start_oas=0
end_oas=$(($start_oas+$nproc_oas-1))
start_cos=$nproc_oas
end_cos=$(($start_cos+$nproc_cos-1))

start_pfl=$(($nproc_cos+$nproc_oas))
end_pfl=$(($start_pfl+$nproc_pfl-1))

start_clm=$(($nproc_cos+$nproc_oas+$nproc_pfl))
end_clm=$(($start_clm+$nproc_clm-1))

cat << EOF >> $rundir/tsmp_slm_run.bsh
#!/bin/bash

USAGE="sbatch <scriptname>"

# JUROPATEST module load intel-para/2014.11
# all in same directory, copy beforehand manually to $WORK

#SBATCH --job-name="TerrSysMP"
#SBATCH --nodes=$nnodes
#SBATCH --ntasks=$mpitasks
#SBATCH --ntasks-per-node=$nppn
#SBATCH --output=mpiMPMD-out.%j
#SBATCH --error=mpiMPMD-err.%j
#SBATCH --time=$wtime
#SBATCH --partition=$queue
#SBATCH --mail-type=ALL

source $rundir/loadenvs
date
rm -rf YU*
srun --multi-prog slm_multiprog_mapping.conf
date
echo "ready" > ready.txt
exit 0

EOF


cat << EOF >> $rundir/slm_multiprog_mapping.conf
__oas__
__cos__
__pfl__
__clm__
EOF

if [[ $withOAS == "false" ||  $withOASMCT == "true" ]] then ; sed "s/__oas__//" -i $rundir/slm_multiprog_mapping.conf ; fi

if [[ $withCOS == "true" && $withCLM == "false" ]] then
sed "s/__cos__/$start_cos-$end_cos  \.\/lmparbin_pur/" -i $rundir/slm_multiprog_mapping.conf
fi
if [[ $withCLM == "true" && $withCOS == "false" && $withPFL == "false" ]] then
sed "s/__clm__/$start_clm-$end_clm  \.\/clm/" -i $rundir/slm_multiprog_mapping.conf
fi
if [[ $withPFL == "true" && $withCLM == "false" ]] then
sed "s/__pfl__/$start_pfl-$end_pfl  \.\/parflow rurlaf/" -i $rundir/slm_multiprog_mapping.conf
fi
if [[ $withCOS == "true" && $withCLM == "true" && $withPFL == "false" ]] then
sed "s/__oas__/$start_oas  \.\/oasis3.MPI1.x/" -i $rundir/slm_multiprog_mapping.conf
sed "s/__cos__/$start_cos-$end_cos  \.\/lmparbin_pur/" -i $rundir/slm_multiprog_mapping.conf
sed "s/__clm__/$start_clm-$end_clm  \.\/clm/" -i $rundir/slm_multiprog_mapping.conf
fi
if [[ $withPFL == "true" && $withCLM == "true" && $withCOS == "false" ]] then
sed "s/__oas__/$start_oas  \.\/oasis3.MPI1.x/" -i $rundir/slm_multiprog_mapping.conf
sed "s/__clm__/$start_clm-$end_clm  \.\/clm/" -i $rundir/slm_multiprog_mapping.conf
sed "s/__pfl__/$start_pfl-$end_pfl  \.\/parflow rurlaf/" -i $rundir/slm_multiprog_mapping.conf
fi
if [[ $withCOS == "true" && $withCLM == "true" && $withPFL == "true" ]] then
sed "s/__oas__/$start_oas  \.\/oasis3.MPI1.x/" -i $rundir/slm_multiprog_mapping.conf
sed "s/__clm__/$start_clm-$end_clm  \.\/clm/" -i $rundir/slm_multiprog_mapping.conf
sed "s/__pfl__/$start_pfl-$end_pfl  \.\/parflow rurlaf/" -i $rundir/slm_multiprog_mapping.conf
sed "s/__cos__/$start_cos-$end_cos  \.\/lmparbin_pur/" -i $rundir/slm_multiprog_mapping.conf
fi

sed "s/__oas__//" -i $rundir/slm_multiprog_mapping.conf
sed "s/__cos__//" -i $rundir/slm_multiprog_mapping.conf
sed "s/__pfl__//" -i $rundir/slm_multiprog_mapping.conf
sed "s/__clm__//" -i $rundir/slm_multiprog_mapping.conf
sed '/^$/d' -i $rundir/slm_multiprog_mapping.conf

chmod 755 $rundir/tsmp_slm_run.bsh
chmod 755 $rundir/slm_multiprog_mapping.conf

}

