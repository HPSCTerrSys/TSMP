#! /bin/ksh

initSetup(){
  defaultFDCLM="$rootdir/tsmp_idealscal/input_1200600/clm"
  defaultFDCOS=""
  defaultFDOAS="$rootdir/tsmp_idealscal/input_1200600/oasis3"
  defaultFDPFL=""


  defaultNLCLM=$rootdir/bldsva/setups/ideal1200600/lnd.stdin 
  defaultNLCOS=$rootdir/bldsva/setups/ideal1200600/lmrun_uc 
  defaultNLPFL=$rootdir/bldsva/setups/ideal1200600/coup_oas.tcl 


  defaultNppn=128
  defaultCLMProcX=16
  defaultCLMProcY=8
  defaultCOSProcX=32
  defaultCOSProcY=16
  defaultPFLProcX=16
  defaultPFLProcY=16

  defaultStartDate="2008-05-08 00"
  defaultInitDate="2008-05-08 00"
  defaultRunhours=3

  defaultDumpCLM=1
  defaultDumpCOS=1
  defaultDumpPFL=1

  gx_clm=1200
  gy_clm=1200
  dt_clm=900
  res="1200x1200"

  gx_cos=600
  gy_cos=600
  dt_cos=10
  nbndlines=4

  gx_pfl=1200
  gy_pfl=1200
  dt_pfl=0.25
  pflrunname="rurlaf"
  base_pfl=0.0025

  cplfreq1=900
  cplfreq2=900


  if [[ $withPFL == "false" && $withCOS == "true" ]]; then
    if [[ $cplscheme == "false" ]]; then
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm_a1
    else
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm
    fi
  fi
  if [[ $withPFL == "true" && $withCOS == "false" ]]; then
    defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_pfl_clm
  fi
  if [[ $withPFL == "true" && $withCOS == "true" ]]; then
    if [[ $cplscheme == "false" ]]; then
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm_pfl_a1
    else
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm_pfl
    fi
  fi

}

finalizeSetup(){
route "${cyellow}>> finalizeSetup${cnormal}"
  if [[ $withOAS == "true" ]] then
    comment "   copy clmgrid into rundir"
      cp $forcingdir_clm/grid* $rundir/clmgrid.nc >> $log_file 2>> $err_file
    check
    comment "   copy oasis remappingfiles into rundir"
      cp $forcingdir_oas/* $rundir >> $log_file 2>> $err_file
    check
  fi  

route "${cyellow}<< finalizeSetup${cnormal}"
}
