#! /bin/ksh

initSetup(){
  defaultFDCLM="/work/slts/slts00/tsmp/TerrSysMPdb/idealized/1200"
  defaultFDCOS=""
  defaultFDOAS="/work/slts/slts00/tsmp/TerrSysMPdb/idealized/1200/oasis3"
  defaultFDPFL=""


  defaultNLCLM=$rootdir/bldsva/setups/ideal1200600/lnd.stdin 
  defaultNLCOS=$rootdir/bldsva/setups/ideal1200600/lmrun_uc 
  defaultNLPFL=$rootdir/bldsva/setups/ideal1200600/coup_oas.tcl 


  defaultNppn=16
  defaultCLMProcX=8
  defaultCLMProcY=8
  defaultCOSProcX=23
  defaultCOSProcY=16
  defaultPFLProcX=8
  defaultPFLProcY=8

  defaultStartDate="2008-05-08 00"
  defaultRestDate=""
  defaultRunhours=3

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

  cplfreq1=900
  cplfreq2=900


  if [[ $withPFL == "false" && $withCOS == "true" ]]; then
    defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm
  fi	
  if [[ $withPFL == "true" && $withCOS == "false" ]]; then
    defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_pfl_clm
  fi
  if [[ $withPFL == "true" && $withCOS == "true" ]]; then
    defaultNLOAS=$rootdir/bldsva/data_oas3/namcouple_cos_clm_pfl
  fi

  fn_finidat="$WORK/tsmp/TSMPForecastNRW$restDate-00/run/clmoas.clm2.r.${yyyy}-${mm}-${dd}-00000.nc"
  pfbfilename="/work/slts/slts06/tsmp/TSMPForecastNRW$restDate-00/run/rurlaf.out.press.00024.pfb"

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
