#! /bin/ksh

initSetup(){
  defaultFDCLM="/perm/ms/spde/de5j/inputdata/cordex11.testcase/clm"
  defaultFDCOS="/scratch/ms/spde/de5j/bound4eur11/20020828"
  defaultFDOAS="/scratch/ms/spde/de5j/juq2ecmwf/rmp"
  defaultFDPFL="/perm/ms/spde/de5j/inputdata/cordex11.testcase/parflow"

  defaultNLCLM=$rootdir/bldsva/setups/cordex/lnd.stdin_jk 
  defaultNLCOS=$rootdir/bldsva/setups/cordex/lmrun_uc_jk 
  defaultNLPFL=$rootdir/bldsva/setups/cordex/coup_oas_jk.tcl 

  defaultNppn=36
  defaultCLMProcX=2
  defaultCLMProcY=2
  defaultCOSProcX=9
  defaultCOSProcY=4
  defaultPFLProcX=8
  defaultPFLProcY=4

  defaultStartDate="2002-08-28 00"
  defaultInitDate="2002-08-28 00"
  defaultRunhours=3

  defaultDumpCLM=1
  defaultDumpCOS=1
  defaultDumpPFL=1

  gx_clm=436
  gy_clm=424
  dt_clm=3600
  res="436x424"

  gx_cos=444
  gy_cos=432
  dt_cos=60
  nbndlines=4

  gx_pfl=436
  gy_pfl=424
  dt_pfl=1.0
  pflrunname="cordex0.11"
  base_pfl=0.001

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
      #cp $forcingdir_clm/clm3.5/grid* $rundir/clmgrid.nc >> $log_file 2>> $err_file
    check

    comment "   copy oasis remappingfiles into rundir"
      cp $forcingdir_oas/* $rundir >> $log_file 2>> $err_file
    check
    #if [[ $withOASMCT == "true" ]] then
      #for x in $rundir/*BILINEA* ;do 
      #  comment "   rename oasis3 remapping files" 
      #    mv $x $(echo $x | sed "s/BILINEA/BILINEAR/") >> $log_file 2>> $err_file
      #  check 
      #done
    #fi  
  fi  

  if [[ $withPFL == "true" ]] then
        comment "   cd to rundir"
          cd $rundir >> $log_file 2>> $err_file
        check

        comment "   copy slopes and slope script into rundir"
          #cp $forcingdir_pfl/slopes/ascii2pfb_slopes.tcl.template $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
          cp $forcingdir_pfl/slopes/alps.new/ascii2pfb_slopes.tcl.template $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check
          #cp $forcingdir_pfl/slopes/slope*.sa $rundir >> $log_file 2>> $err_file
          cp $forcingdir_pfl/slopes/alps.new/slope*.sa $rundir >> $log_file 2>> $err_file
	check
          chmod u+w $rundir/slope*.sa  $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
        check
	comment "   sed procs into slopescript"
          sed "s,lappend auto_path.*,lappend auto_path $bindir/bin," -i $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocx_pfl__,$px_pfl," -i $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocy_pfl__,$py_pfl," -i $rundir/ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check
	comment "   create sloap pfb with tclsh"
          tclsh ./ascii2pfb_slopes.tcl >> $log_file 2>> $err_file
	check		
	
	comment "   copy soilind and soilind script into rundir"
          #cp $forcingdir_pfl/soilInd/ascii2pfb_ind.tcl.template $rundir/ascii2pfb.tcl >> $log_file 2>> $err_file
          cp $forcingdir_pfl/indi/ascii2pfb_ind.tcl.template $rundir/ascii2pfb_ind.tcl >> $log_file 2>> $err_file
	check
        indipath=parflow_436x424x15_cosmomask_indicator_FAOonly.sa
	comment "   copy $indipath  and soilind script into rundir"
          cp $forcingdir_pfl/indi/$indipath $rundir/pfl_ind.sa >> $log_file 2>> $err_file
	check
          chmod u+w $rundir/pfl_ind.sa $rundir/ascii2pfb_ind.tcl >> $log_file 2>> $err_file
        check
	comment "   sed procs into soilindscript"
          sed "s,lappend auto_path.*,lappend auto_path $bindir/bin," -i $rundir/ascii2pfb_ind.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocx_pfl__,$px_pfl," -i $rundir/ascii2pfb_ind.tcl >> $log_file 2>> $err_file
	check
          sed "s,__nprocy_pfl__,$py_pfl," -i $rundir/ascii2pfb_ind.tcl >> $log_file 2>> $err_file
	check
	comment "   create soilInd pfb with tclsh"
        tclsh ./ascii2pfb_ind.tcl >> $log_file 2>> $err_file
	check
  fi 
route "${cyellow}<< finalizeSetup${cnormal}"
}
