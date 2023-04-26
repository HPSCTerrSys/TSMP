#! /bin/ksh

initSetup(){
  defaultFDCLM="/p/project/chbn33/hbn331/database/rur/clm"
  defaultFDOAS="/p/project/chbn33/hbn331/database/rur/oasis3"
  defaultFDPFL="/p/project/chbn33/hbn331/database/rur/parflow"


  defaultNLCLM=$rootdir/bldsva/setups/rur/lnd.stdin 
  defaultNLPFL=$rootdir/bldsva/setups/rur/coup_oas.tcl


  defaultNppn=48
  defaultCLMProcX=4
  defaultCLMProcY=3
  defaultPFLProcX=12
  defaultPFLProcY=7

  defaultStartDate="2008-01-01 00"
  defaultInitDate="2008-01-01 00"
  defaultRunhours=8760

  gx_clm=360
  gy_clm=376
  dt_clm=3600
  res="0376x0360"

  gx_pfl=360
  gy_pfl=376
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
route "${cyellow}>> finalizeSetup${cnormal}"
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
          cp $forcingdir_pfl/*.pfsol $rundir >> $log_file 2>> $err_file
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

        comment "   copy initial pressure, ksat, porosity and script into rundir"
          cp $forcingdir_pfl/ascii2pfb_psi.tcl $rundir/ascii2pfb_psi.tcl >> $log_file 2>> $err_file
        check
          cp $forcingdir_pfl/Rur*.pfb $rundir >> $log_file 2>> $err_file
        check
          chmod u+w $rundir/Rur*.pfb  $rundir/ascii2pfb_psi.tcl >> $log_file 2>> $err_file
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
route "${cyellow}<< finalizeSetup${cnormal}"
}
