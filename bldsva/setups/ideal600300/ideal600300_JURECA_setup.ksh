#! /bin/ksh

initSetup(){
  defaultFDCLM="/work/slts/slts00/tsmp/TerrSysMPdb/idealized/600"
  defaultFDCOS=""
  defaultFDOAS="/work/slts/slts00/tsmp/TerrSysMPdb/idealized/600/oasis3"
  defaultFDPFL=""


  defaultNLCLM=$rootdir/bldsva/setups/ideal600300/lnd.stdin 
  defaultNLCOS=$rootdir/bldsva/setups/ideal600300/lmrun_uc 
  defaultNLPFL=$rootdir/bldsva/setups/ideal600300/coup_oas.tcl 


  defaultNppn=48
  defaultCLMProcX=4
  defaultCLMProcY=2
  defaultCOSProcX=8
  defaultCOSProcY=8
  defaultPFLProcX=4
  defaultPFLProcY=4

  defaultStartDate="2008-05-08 00"
  defaultInitDate="2008-05-08 00"
  defaultRunhours=3

  defaultDumpCLM=1
  defaultDumpCOS=1
  defaultDumpPFL=1

  gx_clm=600
  gy_clm=600
  dt_clm=900
  res="0600x0600"

  gx_cos=300
  gy_cos=300
  dt_cos=10
  nbndlines=4

  gx_pfl=600
  gy_pfl=600
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
route "${cblue}>> finalizeSetup${cnormal}"
  if [[ $withOAS == "true" ]] then
    comment "   copy clmgrid into rundir"
      cp $forcingdir_clm/clm3.5/ideal/grid* $rundir/clmgrid.nc >> $log_file 2>> $err_file
    check

    comment "   copy oasis remappingfiles into rundir"
      cp $forcingdir_oas/* $rundir >> $log_file 2>> $err_file
    check
    if [[ $withOASMCT == "true" ]] then
      for x in $rundir/*BILINEA* ;do 
        comment "   rename oasis3 remapping files" 
          mv $x $(echo $x | sed "s/BILINEA/BILINEAR/") >> $log_file 2>> $err_file
        check 
      done
    fi  
  fi  

route "${cblue}<< finalizeSetup${cnormal}"
}
