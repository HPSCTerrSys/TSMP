#! /bin/ksh

initSetup(){
  defaultFDCLM="/p/scratch/cslts/shared_data/rcmod_TSMP-ref_SLTS/TestCases/nrw_5x/clm"
  defaultFDCOS="/p/scratch/cslts/shared_data/rcmod_TSMP-ref_SLTS/TestCases/nrw_5x/cosmo"
  defaultFDOAS="/p/scratch/cslts/shared_data/rcmod_TSMP-ref_SLTS/TestCases/nrw_5x/oasis3"
  defaultFDPFL="/p/scratch/cslts/shared_data/rcmod_TSMP-ref_SLTS/TestCases/nrw_5x/parflow"

  defaultNLCLM=$rootdir/bldsva/setups/nrw_5x/lnd.stdin 
  defaultNLCOS=$rootdir/bldsva/setups/nrw_5x/lmrun_uc 
  defaultNLPFL=$rootdir/bldsva/setups/nrw_5x/coup_oas.tcl
  defaultNLDA=$rootdir/bldsva/setups/nrw_5x/enkfpf.par

  defaultNppn=24
  defaultCLMProcX=4
  defaultCLMProcY=24
  defaultCOSProcX=1
  defaultCOSProcY=1
  defaultPFLProcX=1
  defaultPFLProcY=1
  
  defaultStartDate="2016-01-01 00"
  defaultInitDate="2016-01-02 00"
  defaultRunhours=24

  defaultDumpCLM=1
  defaultDumpCOS=1
  defaultDumpPFL=1	

  gx_clm=300
  gy_clm=300
  dt_clm=1800
  res="0300x0300"

  gx_cos=150
  gy_cos=150
  dt_cos=10
  nbndlines=4

  gx_pfl=300
  gy_pfl=300
  dt_pfl=0.25
  pflrunname="rurlaf"
  base_pfl=0.0025

  cplfreq1=900
  cplfreq2=900

  delta_obs=1
 
#  if [[ $withPFL == "false" && $withCOS == "true" ]]; then
#    if [[ $cplscheme == "false" ]]; then
#      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm_a1
#    else
#      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm
#    fi
#  fi
#  if [[ $withPFL == "true" && $withCOS == "false" ]]; then
#    defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_pfl_clm
#  fi
#  if [[ $withPFL == "true" && $withCOS == "true" ]]; then
#    if [[ $cplscheme == "false" ]]; then
#      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm_pfl_a1
#    else
#      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm_pfl
#    fi
#  fi


}

finalizeSetup(){
route "${cblue}>> finalizeSetup${cnormal}"

if [[ $withPDAF == "true" ]]; then

comment "copy PDAF files into rundir"
cp $rootdir/bldsva/setups/nrw_5x/create_ensemble_namelists.py $rundir
check
  
comment "create CLM5+PDAF run script"
check

cat << EOF >> $rundir/clm5pdaf_run.bsh
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=02:00:00
#SBATCH --partition=batch
#SBATCH --job-name=C5P
#SBATCH --output=out-clm5pdaf.%j
#SBATCH --error=err-clm5pdaf.%j
#SBATCH --account=

nensemble=2

source $rundir/loadenvs
srun ./tsmp-pdaf -n_modeltasks \$nensemble -screen 3 -filtertype 2 -subtype 1 -delt_obs 1 -rms_obs 0.02 -obs_filename /p/scratch/cslts/shared_data/rcmod_TSMP-ref_SLTS/TestCases/nrw_5x/pdaf/obs/swc_obs

EOF
check

chmod +x $rundir/clm5pdaf_run.bsh >> $log_file 2>> $err_file
check

mkdir -p $rundir/timing/checkpoints
check
mkdir  $rundir/logs
check

else
comment " create CLM5 specific run script "

cat << EOF >> $rundir/tsmp_clm_run.bsh
#!/bin/bash

#SBATCH --job-name="TerrSysMP"
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --output=out-clm5_std.%j
#SBATCH --error=err-clm5_std.%j
#SBATCH --time=01:00:00
#SBATCH --partition=batch
#SBATCH --mail-type=NONE
#SBATCH --account=slts

cd $rundir 
source $rundir/loadenvs
date
echo "started" > started.txt
datefmt=\`date +%d%m%y-%H%M%S\`
for fname in atm cpl esp glc ice rof ocn wav; do
sed -i 's/.*logfile.*/  logfile = "'\${fname}'.log.'\${datefmt}'"/' \${fname}_modelio.nml ;
done
srun clm >> cesm.log.\$datefmt 2>&1
date
echo "ready" > ready.txt
exit 0
EOF

check
chmod +x $rundir/tsmp_clm_run.bsh >> $log_file 2>> $err_file  
check

mkdir -p $rundir/timing/checkpoints
check

mkdir $rundir/logs

fi
route "${cblue}<< finalizeSetup${cnormal}"
}
