#! /bin/ksh

initSetup(){
  defaultFDCLM="/homea/hbn33/hbn331/database/smresponse/clm"
  defaultFDCOS="/homea/hbn33/hbn331/database/smresponse/cosmo"
  defaultFDOAS="/homea/hbn33/hbn331/database/smresponse/oasis3"
  defaultFDPFL="/homea/hbn33/hbn331/database/smresponse/parflow"

  defaultNLCLM=$rootdir/bldsva/setups/smresponse/lnd.stdin 
  defaultNLCOS=$rootdir/bldsva/setups/smresponse/lmrun_uc 
  defaultNLPFL=$rootdir/bldsva/setups/smresponse/coup_oas.tcl


  defaultNppn=48
  defaultCLMProcX=2
  defaultCLMProcY=2
  defaultCOSProcX=2
  defaultCOSProcY=2

  defaultStartDate="2015-08-07 00"
  defaultInitDate="2015-08-07 00"
  defaultRunhours=24

  defaultDumpCLM=1
  defaultDumpCOS=1

  gx_clm=14
  gy_clm=14
  dt_clm=10
  res="0014x0014"

  gx_cos=20
  gy_cos=20
  dt_cos=10
  nbndlines=3

  cplfreq1=10

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
      cp $forcingdir_clm/grid* $rundir/clmgrid.nc >> $log_file 2>> $err_file
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
