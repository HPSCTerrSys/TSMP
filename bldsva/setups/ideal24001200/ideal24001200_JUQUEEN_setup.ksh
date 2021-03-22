#! /bin/ksh

initSetup(){
  defaultFDCLM="/work/slts/slts00/tsmp/TerrSysMPdb/idealized/2400"
  defaultFDCOS=""
  defaultFDOAS="/work/slts/slts00/tsmp/TerrSysMPdb/idealized/2400/oasis3"
  defaultFDPFL=""


  defaultNLCLM=$rootdir/bldsva/setups/ideal24001200/lnd.stdin 
  defaultNLCOS=$rootdir/bldsva/setups/ideal24001200/lmrun_uc 
  defaultNLPFL=$rootdir/bldsva/setups/ideal24001200/coup_oas.tcl 


  defaultNppn=16
  defaultCLMProcX=8
  defaultCLMProcY=8
  defaultCOSProcX=23
  defaultCOSProcY=16
  defaultPFLProcX=8
  defaultPFLProcY=8

  defaultStartDate="2008-05-08 00"
  defaultInitDate="2008-05-08 00"
  defaultRunhours=3

  defaultDumpCLM=1
  defaultDumpCOS=1
  defaultDumpPFL=1

  gx_clm=2400
  gy_clm=2400
  dt_clm=900
  res="2400x2400"

  gx_cos=1200
  gy_cos=1200
  dt_cos=10
  nbndlines=4

  gx_pfl=2400
  gy_pfl=2400
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

route "${cyellow}<< finalizeSetup${cnormal}"
}
