#! /bin/ksh

initSetup(){
  defaultFDCLM="/daten01/z4/database/TestCases/idealRTD/clm"
  defaultFDCOS="/daten01/z4/database/TestCases/idealRTD/cosmo"
  defaultFDOAS="/daten01/z4/database/TestCases/idealRTD/oasis3"
  defaultFDPFL="/daten01/z4/database/TestCases/idealRTD/parflow"

  defaultNLCLM=$rootdir/bldsva/setups/idealRTD/lnd.stdin 
  defaultNLCOS=$rootdir/bldsva/setups/idealRTD/lmrun_uc 
  defaultNLPFL=$rootdir/bldsva/setups/idealRTD/coup_oas.tcl


  defaultNppn=64
  defaultCLMProcX=2
  defaultCLMProcY=2
  defaultCOSProcX=3
  defaultCOSProcY=2
  defaultPFLProcX=2
  defaultPFLProcY=2

  defaultStartDate="2015-04-01 00"
  defaultInitDate="2015-04-01 00"
  defaultRunhours=24

  defaultDumpCLM=3
  defaultDumpCOS=3
  defaultDumpPFL=3


  gx_clm=24
  gy_clm=14
  dt_clm=18
  res="0014x0024"

  gx_cos=30
  gy_cos=20
  dt_cos=18
  nbndlines=3

  gx_pfl=24
  gy_pfl=14
  dt_pfl=0.005
  pflrunname="rurlaf"
  base_pfl=0.0025

  cplfreq1=18
  cplfreq2=18


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

  if [[ $withPFL == "true" ]] then


        comment "   cd to rundir"
          cd $rundir >> $log_file 2>> $err_file
        check
        comment "   copy parflow soilID and geodata into rundir"
          cp $forcingdir_pfl/pfb*.nc $rundir/ >> $log_file 2>> $err_file
        check
        comment "   copy initial pressure and script into rundir"
          cp $forcingdir_pfl/ascii2pfb.tcl $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
        # Specific perturbed surface data and spinup pressure data are copied using DART scripts

        #  cp $forcingdir_pfl/rur_ic_press.pfb $rundir >> $log_file 2>> $err_file
        #check
        #  chmod u+w $rundir/rur_ic_press.pfb  $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        chmod u+w $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
        comment "   sed procs into pfbscript"
          sed "s,lappend auto_path.*,lappend auto_path $bindir/bin," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
          sed "s,pfset Process\.Topology\.P.*,pfset Process\.Topology\.P $px_pfl," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
          sed "s,pfset Process\.Topology\.Q.*,pfset Process\.Topology\.Q $py_pfl," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
        #comment "   create sloap pfb with tclsh"
        #  tclsh ./ascii2pfb.tcl >> $log_file 2>> $err_file
        #check
              

  fi 
route "${cblue}<< finalizeSetup${cnormal}"
}
