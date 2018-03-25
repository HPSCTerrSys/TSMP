#! /bin/ksh

initSetup(){
  defaultFDOAS="/homea/slts/slts23/data/tsmp/setups/icon-ccs/oasis"
  defaultFDCLM="/homea/slts/slts23/data/tsmp/setups/icon-ccs/clm"
  defaultFDPFL="/homea/slts/slts23/data/tsmp/setups/icon-ccs/pfl"

  defaultNLICON="$rootdir/bldsva/setups/icon-ccs/icon_master.namelist $rootdir/bldsva/setups/icon-ccs/NAMELIST_ccs"
  defaultNLCLM=$rootdir/bldsva/setups/icon-ccs/lnd.stdin 
  defaultNLPFL=$rootdir/bldsva/setups/icon-ccs/coup_oas.tcl 


  defaultNppn=24
  defaultICONProc=52
  defaultCLMProcX=2
  defaultCLMProcY=2
  defaultPFLProcX=4
  defaultPFLProcY=4

  defaultStartDate="2008-05-08 00"
  defaultInitDate="2008-05-08 00"
  defaultRunhours=0

  defaultDumpICON=1
  defaultDumpCLM=1
  defaultDumpPFL=1

  gx_icon=73728

  gx_clm=192
  gy_clm=192
  dt_clm=1
  res="0192x0192"

  gx_pfl=192
  gy_pfl=192
  dt_pfl=1.0
  pflrunname="rurlaf"
  base_pfl=0.0025

  cplfreq1=10
  cplfreq2=10

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
      comment "   rename oasis3 remapping files (SBr: skip for now!)"
      #for x in $rundir/*BILINEA* ;do 
      #  comment "   rename oasis3 remapping files" 
      #    mv $x $(echo $x | sed "s/BILINEA/BILINEAR/") >> $log_file 2>> $err_file
      #  check 
      #done
    fi  
  fi  

  comment "  copy grid for icon"
    ln -s /homea/slts/slts23/data/tsmp/setups/icon-ccs/icon/torus_grid_x192_y192_e70m.nc $rundir/ >> $log_file 2>> $err_file
  check

route "${cblue}<< finalizeSetup${cnormal}"
}
