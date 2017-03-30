#! /bin/ksh

initSetup(){
  defaultFDCLM="/daten01/z4/database/TestCases/idealLES/clm"
  defaultFDCOS="/daten01/z4/database/TestCases/idealLES/cosmo"
  defaultFDOAS="/daten01/z4/database/TestCases/idealLES/oasis3"
  defaultFDPFL="/daten01/z4/database/TestCases/idealLES/parflow"

  defaultNLCLM=$rootdir/bldsva/setups/idealLES/lnd.stdin 
  defaultNLCOS=$rootdir/bldsva/setups/idealLES/lmrun_uc 
  defaultNLPFL=$rootdir/bldsva/setups/idealLES/coup_oas.tcl


  defaultNppn=64
  defaultCLMProcX=1
  defaultCLMProcY=1
  defaultCOSProcX=8
  defaultCOSProcY=8
  defaultPFLProcX=1
  defaultPFLProcY=1

  defaultStartDate="2008-05-08 00"
  defaultInitDate="2008-05-08 00"
  defaultRunhours=2

  defaultDumpCLM=1
  defaultDumpCOS=1
  defaultDumpPFL=1

  gx_clm=116
  gy_clm=116
  dt_clm=6
  res="0116x0116"

  gx_cos=120
  gy_cos=120
  dt_cos=6
  nbndlines=3
  dump_cos=0.0833333333333333333333333333

  gx_pfl=116
  gy_pfl=116
  dt_pfl=0.00166667
  pflrunname="rurlaf"
  base_pfl=0.00166667
  dump_pfl=3.0

  cplfreq1=6
  cplfreq2=6


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
      cp $forcingdir_clm/clm3.5/idealLES/grid* $rundir/clmgrid.nc >> $log_file 2>> $err_file
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

        comment "   copy initial pressure and script into rundir"
          cp $forcingdir_pfl/ascii2pfb.tcl $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
          cp $forcingdir_pfl/rur_ic_press.pfb $rundir >> $log_file 2>> $err_file
        check
          chmod u+w $rundir/rur_ic_press.pfb  $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
        comment "   sed procs into pfbscript"
          sed "s,lappend auto_path.*,lappend auto_path $bindir/bin," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
          sed "s,pfset Process\.Topology\.P.*,pfset Process\.Topology\.P $px_pfl," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
          sed "s,pfset Process\.Topology\.Q.*,pfset Process\.Topology\.Q $py_pfl," -i $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
        check
        comment "   create slope pfb with tclsh"
          tclsh ./ascii2pfb.tcl >> $log_file 2>> $err_file
        check
              

  fi 
route "${cblue}<< finalizeSetup${cnormal}"
}
