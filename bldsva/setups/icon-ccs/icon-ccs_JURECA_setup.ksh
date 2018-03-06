#! /bin/ksh

initSetup(){
  defaultFDCLM="/work/slts/slts00/tsmp/TestCases/ideal300150/clm"
  defaultFDCOS=""
  defaultFDOAS="/work/slts/slts00/tsmp/TestCases/ideal300150/oasis3"
  defaultFDPFL=""


  defaultNLCLM=$rootdir/bldsva/setups/ideal300150/lnd.stdin 
  defaultNLICON="$rootdir/bldsva/setups/icon-ccs/icon_master.namelist $rootdir/bldsva/setups/icon-ccs/NAMELIST_ccs"
  #defaultNLPFL=$rootdir/bldsva/setups/ideal300150/coup_oas.tcl 


  defaultNppn=2
  defaultCLMProcX=1
  defaultCLMProcY=1
  #defaultCOSProcX=8
  #defaultCOSProcY=8
  defaultICONProc=1
  #defaultPFLProcX=4
  #defaultPFLProcY=4

  defaultStartDate="2008-05-08 00"
  defaultInitDate="2008-05-08 00"
  defaultRunhours=1

  defaultDumpCLM=1
  #defaultDumpCOS=1
  defaultDumpICON=1
  #defaultDumpPFL=1

  gx_clm=300
  gy_clm=300
  dt_clm=1
  res="0300x0300"

  #gx_cos=150
  #gy_cos=150
  #dt_cos=10
  #nbndlines=4

  #gx_pfl=300
  #gy_pfl=300
  #dt_pfl=0.25
  #pflrunname="rurlaf"
  #base_pfl=0.0025

  cplfreq1=1
  #cplfreq2=900


  if [[ $withPFL == "false" && $withICON == "true" ]]; then
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_icon_clm
  fi
  if [[ $withPFL == "true" && $withICON == "false" ]]; then
    defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_pfl_clm
  fi
  if [[ $withPFL == "true" && $withICON == "true" ]]; then
      defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_icon_clm_pfl
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
