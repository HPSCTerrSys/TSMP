#! /bin/ksh

initSetup(){
  defaultFDCLM="/daten01/z4/database/multi-scale/D0480subcatchment"
  defaultFDOAS="/daten01/z4/database/multi-scale/D0480subcatchment"
  defaultFDPFL="/daten01/z4/database/multi-scale/D0480subcatchment"


  defaultNLCLM=$rootdir/bldsva/setups/multi-scale/d480/lnd.stdin 
  defaultNLPFL=$rootdir/bldsva/setups/multi-scale/d480/coup_oas.tcl


  defaultNppn=64
  defaultCLMProcX=2
  defaultCLMProcY=7
  defaultPFLProcX=10
  defaultPFLProcY=10

  defaultStartDate="2009-01-01 00"
  defaultInitDate="2009-01-01 00"
  defaultRunhours=8760

  gx_clm=50
  gy_clm=70
  dt_clm=3600
  res="0070x0050"

  gx_pfl=50
  gy_pfl=70
  dt_pfl=1.0
  pflrunname="rurlaf"
  base_pfl=1.0
  dump_pfl=120.0

  cplfreq1=3600
  cplfreq2=3600


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

        comment "   copy slopes and slope script into rundir"
          cp $forcingdir_pfl/ascii2pfb_slope.tcl $rundir/ascii2pfb_slope.tcl >> $log_file 2>> $err_file
        check
          cp $forcingdir_pfl/xslope*.pfb $rundir >> $log_file 2>> $err_file
          cp $forcingdir_pfl/yslope*.pfb $rundir >> $log_file 2>> $err_file
        check
          chmod u+w $rundir/xslope*.pfb  $rundir/ascii2pfb_slope.tcl >> $log_file 2>> $err_file
          chmod u+w $rundir/yslope*.pfb  $rundir/ascii2pfb_slope.tcl >> $log_file 2>> $err_file
        check
        comment "   sed procs into slopescript"
          sed "s,lappend auto_path.*,lappend auto_path $bindir/bin," -i $rundir/ascii2pfb_slope.tcl >> $log_file 2>> $err_file
        check
          sed "s,pfset Process\.Topology\.P.*,pfset Process\.Topology\.P $px_pfl," -i $rundir/ascii2pfb_slope.tcl >> $log_file 2>> $err_file
        check
          sed "s,pfset Process\.Topology\.Q.*,pfset Process\.Topology\.Q $py_pfl," -i $rundir/ascii2pfb_slope.tcl >> $log_file 2>> $err_file
        check
        comment "   create sloap pfb with tclsh"
          tclsh ./ascii2pfb_slope.tcl >> $log_file 2>> $err_file
        check

        comment "   copy initial pressure and script into rundir"
          cp $forcingdir_pfl/ascii2pfb_psi.tcl $rundir/ascii2pfb_psi.tcl >> $log_file 2>> $err_file
        check
          cp $forcingdir_pfl/ptb*.pfb $rundir >> $log_file 2>> $err_file
        check
          chmod u+w $rundir/ptb*.pfb  $rundir/ascii2pfb_psi.tcl >> $log_file 2>> $err_file
        check
        comment "   sed procs into pfbscript"
          sed "s,lappend auto_path.*,lappend auto_path $bindir/bin," -i $rundir/ascii2pfb_psi.tcl >> $log_file 2>> $err_file
        check
          sed "s,pfset Process\.Topology\.P.*,pfset Process\.Topology\.P $px_pfl," -i $rundir/ascii2pfb_psi.tcl >> $log_file 2>> $err_file
        check
          sed "s,pfset Process\.Topology\.Q.*,pfset Process\.Topology\.Q $py_pfl," -i $rundir/ascii2pfb_psi.tcl >> $log_file 2>> $err_file
        check
        comment "   create slope pfb with tclsh"
          tclsh ./ascii2pfb_psi.tcl >> $log_file 2>> $err_file
        check
              

  fi 
route "${cblue}<< finalizeSetup${cnormal}"
}
